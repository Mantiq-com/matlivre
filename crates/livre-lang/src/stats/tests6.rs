//! `hougen`, `ecdfhist`, `testcholdout`, `invpred` and `partialcorri`.

use super::anova::{grp2idx, group_labels, Label};
use super::dist::{binocdf, chi2cdf, normcdf_std, tcdf, tinv};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["hougen", "ecdfhist", "testcholdout", "invpred", "partialcorri"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "hougen" => hougen(a),
        "ecdfhist" => ecdfhist(it, a, nargout),
        "testcholdout" => testcholdout(it, a, nargout),
        "invpred" => invpred(a, nargout),
        "partialcorri" => partialcorri(a, nargout),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn ulp(x: f64) -> f64 {
    let x = x.abs();
    if x == 0.0 {
        return f64::from_bits(1);
    }
    f64::from_bits(x.to_bits() + 1) - x
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// A name/value parser with MATLAB's partial matching; `names` lists the canonical names.
fn parse_pairs<'a>(rest: &'a [Value], names: &[&str]) -> Result<Vec<Option<&'a Value>>, crate::error::LangError> {
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let mut out: Vec<Option<&Value>> = vec![None; names.len()];
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match pick {
            Some(i) => out[i] = Some(&pair[1]),
            None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            None => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    Ok(out)
}

/// `validatestring`: the full name matching `v` (case-insensitive prefix).
fn validate_string(v: &Value, choices: &[&str], func: &str, var: &str) -> Result<usize, crate::error::LangError> {
    let list: Vec<String> = choices.iter().map(|c| format!("'{c}'")).collect();
    let bad = |given: &str| rt::<usize>(format!("Expected {var} to match one of these values:\n\n{}\n\nThe input, '{given}', did not match any of the valid values.", list.join(", ")));
    let _ = func;
    if !v.is_text() {
        return rt(format!("Expected {var} to be text."));
    }
    let given = v.text().unwrap_or_default();
    let t = given.to_lowercase();
    if let Some(i) = choices.iter().position(|c| c.to_lowercase() == t) {
        return Ok(i);
    }
    let hits: Vec<usize> = choices.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.to_lowercase().starts_with(&t)).map(|(i, _)| i).collect();
    match hits.len() {
        1 => Ok(hits[0]),
        0 => bad(&given),
        _ => rt(format!("Expected {var} to match one of these values:\n\n{}\n\nThe input, '{given}', matched more than one valid value.", list.join(", "))),
    }
}

// ---------------------------------------------------------------------------------------------
// hougen

fn hougen(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let (Value::Num(beta), Value::Num(x)) = (&a[0], &a[1]) else { return rt("Input must be numeric.") };
    let nb = beta.data.len();
    if nb < 5 {
        return rt(format!("Index exceeds the number of array elements. Index must not exceed {nb}."));
    }
    let (r, c) = (x.dims[0], x.dims[1..].iter().product::<usize>());
    if c < 3 {
        if c == 0 {
            return rt("Index in position 2 exceeds array bounds.");
        }
        return rt(format!("Index in position 2 exceeds array bounds. Index must not exceed {c}."));
    }
    let (b1, b2, b3, b4, b5) = (beta.data[0], beta.data[1], beta.data[2], beta.data[3], beta.data[4]);
    let single = is_single(beta) || is_single(x);
    let col = |j: usize| &x.data[j * r..(j + 1) * r];
    let (x1, x2, x3) = (col(0), col(1), col(2));
    let y: Vec<f64> = (0..r).map(|i| (b1 * x2[i] - x3[i] / b5) / (1.0 + b2 * x1[i] + b3 * x2[i] + b4 * x3[i])).collect();
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![r, 1], y), single))]))
}

// ---------------------------------------------------------------------------------------------
// ecdfhist

fn ecdfhist(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let (fv, xv) = (&a[0], &a[1]);
    if fv.is_text() || xv.is_text() {
        return rt("Input arguments must be numeric.");
    }
    let (Value::Num(fm), Value::Num(xm)) = (fv, xv) else { return rt("Input arguments must be numeric.") };
    let vecp = |m: &Matrix| m.data.len() == m.dims.iter().copied().max().unwrap_or(0) || m.data.len() <= 1;
    if !vecp(fm) || !vecp(xm) {
        return rt("Both F and X must be vectors produced by the ecdf function.");
    }
    let mut f = fm.data.clone();
    let x0 = xm.data.clone();
    if f.len() < 2 || x0.len() < 2 {
        return rt("Requires both F and X.");
    }
    if x0.windows(2).any(|w| w[1] - w[0] < 0.0) {
        return rt("X must be non-decreasing.");
    }
    let wassurv = f[0] == 1.0;
    if wassurv {
        f = f.iter().map(|v| 1.0 - v).collect();
    }
    let dif: Vec<f64> = f.windows(2).map(|w| w[1] - w[0]).collect();
    if dif.iter().any(|&d| d < 0.0) {
        return if wassurv { rt("The survivor function F must be non-increasing") } else { rt("The cdf F must be non-decreasing") };
    }
    let x: Vec<f64> = x0[1..].to_vec();
    if dif.len() != x.len() {
        return rt("Second input must be a vector with one element for each row in the first input, or a scalar.");
    }
    // bin specification
    let rest = &a[2..];
    let (mut nbins, mut ctrs, mut edges): (f64, Vec<f64>, Vec<f64>) = (10.0, vec![], vec![]);
    let mut mode = "nbins";
    let spec_vec = |v: &Value| -> Vec<f64> {
        match v {
            Value::Num(m) => m.data.clone(),
            _ => vec![],
        }
    };
    if let Some(b) = rest.first() {
        if b.is_text() {
            match b.text().unwrap_or_default().as_str() {
                "ctrs" => {
                    mode = "ctrs";
                    ctrs = rest.get(1).map(spec_vec).unwrap_or_default();
                }
                "edges" => {
                    mode = "edges";
                    edges = rest.get(1).map(spec_vec).unwrap_or_default();
                }
                "nbins" => {
                    nbins = rest.get(1).map(spec_vec).and_then(|v| v.first().copied()).unwrap_or(10.0);
                }
                _ => return rt("You must specify 'nbins', 'ctrs', or 'edges'."),
            }
        } else {
            let v = spec_vec(b);
            if v.len() == 1 {
                nbins = v[0];
            } else {
                mode = "ctrs";
                ctrs = v;
            }
        }
    }
    let mut lo = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let mut hi = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let use_centers = mode != "edges";
    let mut wng = false;
    let mut cent: Vec<f64>;
    match mode {
        "nbins" => {
            if lo == hi {
                lo = lo - (nbins / 2.0).floor() - 0.5;
                hi = hi + (nbins / 2.0).ceil() - 0.5;
            }
            let bw = (hi - lo) / nbins;
            let nb = nbins as usize;
            edges = (0..=nb).map(|k| lo + bw * k as f64).collect();
            *edges.last_mut().unwrap() = hi;
            cent = (0..nb).map(|k| edges[k] + bw / 2.0).collect();
        }
        "ctrs" => {
            cent = ctrs.clone();
            let mut bw: Vec<f64> = cent.windows(2).map(|w| w[1] - w[0]).collect();
            if bw.iter().any(|&b| b <= 0.0) {
                return rt("The vector C of bin centers must be increasing.");
            }
            bw.push(*bw.last().unwrap_or(&1.0));
            edges = std::iter::once(cent[0] - bw[0] / 2.0).chain(cent.iter().zip(&bw).map(|(c, b)| c + b / 2.0)).collect();
            edges[0] = edges[0].min(lo);
            let l = edges.len() - 1;
            edges[l] = edges[l].max(hi);
        }
        _ => {
            let bw: Vec<f64> = edges.windows(2).map(|w| w[1] - w[0]).collect();
            if bw.iter().any(|&b| b <= 0.0) {
                return rt("The vector E of bin edges must be increasing.");
            }
            cent = vec![0.0; bw.len()];
            if nargout > 1 || true {
                cent[0] = edges[0] + bw[0] / 2.0;
                for j in 1..cent.len() {
                    cent[j] = 2.0 * edges[j] - cent[j - 1];
                }
                let bad = (0..cent.len()).any(|j| cent[j] <= edges[j]) || (cent[cent.len() - 1] - (edges[edges.len() - 1] - bw[bw.len() - 1] / 2.0)).abs() > 1000.0 * ulp(bw[bw.len() - 1]);
                if bad {
                    wng = true;
                    cent = (0..bw.len()).map(|j| edges[j] + bw[j] / 2.0).collect();
                }
            }
            if lo < edges[0] || edges[edges.len() - 1] < hi {
                it.emit_warning("The vector E of bin edges does not span the range of X.".into(), "stats:ecdfhist:DataNotSpanned".into());
            }
        }
    }
    let nbn = edges.len() - 1;
    let mut counts = vec![0.0; nbn];
    let es: Vec<f64> = if use_centers { edges.iter().map(|&e| e + ulp(e)).collect() } else { edges.clone() };
    for (i, &v) in x.iter().enumerate() {
        let bin: Option<usize> = if use_centers {
            // intervals [-Inf, e2), [e2, e3), ..., [e_n, Inf)
            let inner = &es[1..es.len() - 1];
            let k = inner.iter().take_while(|&&e| e <= v).count();
            Some(k)
        } else if v < es[0] || v > es[nbn] || v.is_nan() {
            None
        } else if v == es[nbn] {
            Some(nbn - 1)
        } else {
            Some(es.iter().take_while(|&&e| e <= v).count() - 1)
        };
        if let Some(b) = bin {
            if b < nbn {
                counts[b] += dif[i];
            }
        }
    }
    let widths: Vec<f64> = edges.windows(2).map(|w| w[1] - w[0]).collect();
    let nout: Vec<f64> = counts.iter().zip(&widths).map(|(c, w)| c / w).collect();
    if nargout == 0 {
        return Ok(Some(vec![]));
    }
    let single_f = is_single(fm);
    let mut out = vec![Value::Num(finish(Matrix::from_dims(vec![1, nbn], nout), single_f))];
    if nargout >= 2 {
        out.push(Value::Num(Matrix::from_dims(vec![1, cent.len()], cent)));
        if wng {
            it.emit_warning("Cannot compute centers that are consistent with EDGES.".into(), "stats:ecdfhist:InconsistentEdges".into());
        }
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// testcholdout

fn labels_of(v: &Value) -> Result<Vec<Label>, crate::error::LangError> {
    match v {
        Value::Num(_) | Value::Str(_) => group_labels(v),
        Value::Cell(c) if c.data.iter().all(|e| e.is_text()) => group_labels(v),
        _ => rt("Classification labels must be a character array, a string array, a cell array of character vectors, an array of logicals, a numeric array, a categorical array, or an object of class ClassLabel."),
    }
}

fn is_missing(l: &Label) -> bool {
    match l {
        Label::Num(v) => v.is_nan(),
        Label::Str(s) => s.is_empty(),
    }
}

fn testcholdout(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let _ = it;
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    let y = labels_of(&a[2])?;
    let y1 = labels_of(&a[0])?;
    let y2 = labels_of(&a[1])?;
    let n1 = y1.len();
    if y.len() != n1 {
        return rt("The number of elements in Yhat1 does not match the number of elements in Y.");
    }
    if y2.len() != n1 {
        return rt("The number of elements in Yhat2 does not match the number of elements in Y.");
    }
    let opts = parse_pairs(&a[3..], &["alpha", "alternative", "test", "classnames", "cost", "costtest"])?;
    let alpha = match opts[0] {
        None => 0.05,
        Some(Value::Num(m)) if m.kind != Kind::Char && m.kind != Kind::Logical && m.data.len() == 1 && m.im.is_none() && m.data[0] > 0.0 && m.data[0] < 1.0 => m.data[0],
        Some(_) => return rt("You must pass 'Alpha' as a scalar between 0 and 1."),
    };
    let cost_v = opts[4].filter(|v| !is_empty_val(v));
    let fitcost = cost_v.is_some();
    let alt_v = opts[1].cloned().unwrap_or_else(|| Value::chars("unequal"));
    let alternative = ["unequal", "less", "greater"][validate_string(&alt_v, &["unequal", "less", "greater"], "testcholdout", "Alternative")?];
    let mode_v = match opts[2] {
        Some(v) if !is_empty_val(v) && !(v.is_text() && v.text().unwrap_or_default().is_empty()) => v.clone(),
        _ => Value::chars(if fitcost { "asymptotic" } else { "midp" }),
    };
    let mode = ["asymptotic", "exact", "midp"][validate_string(&mode_v, &["asymptotic", "exact", "midp"], "testcholdout", "Test")?];
    if !(mode == "asymptotic" && alternative == "unequal") && fitcost {
        return rt("If you pass misclassification costs, you must set 'Alternative' to 'unequal' and 'Test' to 'asymptotic'.");
    }
    let ct_v = opts[5].cloned().unwrap_or_else(|| Value::chars("likelihood"));
    let costtest = ["likelihood", "chisquare"][validate_string(&ct_v, &["likelihood", "chisquare"], "testcholdout", "CostTest")?];
    // observations with a missing true label drop out
    let keep: Vec<bool> = y.iter().map(|l| !is_missing(l)).collect();
    if keep.iter().all(|k| !k) {
        return rt("There are no valid class labels in Y.");
    }
    let sel = |v: &[Label]| -> Vec<Label> { v.iter().zip(&keep).filter(|(_, &k)| k).map(|(l, _)| l.clone()).collect() };
    let (mut y, mut y1, mut y2) = (sel(&y), sel(&y1), sel(&y2));
    // the classes of the true labels, in sorted order
    let mut classes: Vec<Label> = Vec::new();
    for l in &y {
        if !classes.contains(l) {
            classes.push(l.clone());
        }
    }
    classes.sort_by(|p, q| match (p, q) {
        (Label::Num(a), Label::Num(b)) => a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal),
        (Label::Str(a), Label::Str(b)) => a.cmp(b),
        (Label::Num(_), Label::Str(_)) => std::cmp::Ordering::Less,
        _ => std::cmp::Ordering::Greater,
    });
    let _ = grp2idx;
    if let Some(cn) = opts[3] {
        if !is_empty_val(cn) {
            let user = labels_of(cn)?;
            if user.iter().all(|u| !classes.contains(u)) {
                return rt("Requested class names are not found among passed class labels.");
            }
            // observations of the classes that were not asked for drop out
            let dropc: Vec<Label> = classes.iter().filter(|c| !user.contains(c)).cloned().collect();
            if !dropc.is_empty() {
                let keep2: Vec<bool> = y.iter().map(|l| !dropc.contains(l)).collect();
                let f = |v: &Vec<Label>| -> Vec<Label> { v.iter().zip(&keep2).filter(|(_, &k)| k).map(|(l, _)| l.clone()).collect() };
                y = f(&y);
                y1 = f(&y1);
                y2 = f(&y2);
                classes.retain(|c| !dropc.contains(c));
            }
        }
    }
    let n = y.len();
    let k = classes.len();
    let (h, p, err1, err2): (bool, f64, f64, f64);
    let class_idx = |l: &Label| classes.iter().position(|c| c == l);
    if fitcost {
        let cm = match cost_v.unwrap() {
            Value::Num(m) if m.dims.len() == 2 && m.dims[0] == k && m.dims[1] == k && m.kind != Kind::Char => m.clone(),
            _ => return rt(format!("Misclassification cost matrix must be {k}-by-{k}.")),
        };
        let cost = |i: usize, j: usize| cm.data[j * k + i];
        if (0..k).any(|i| cost(i, i) != 0.0) {
            return rt("Misclassification cost matrix must have zeros on the diagonal.");
        }
        if cm.data.iter().any(|&c| c < 0.0 || c.is_nan()) {
            return rt("Misclassification cost matrix must contain non-negative values.");
        }
        let mut e1 = 0.0;
        let mut e2 = 0.0;
        let mut mobs = vec![0.0; k * k * k];
        for i in 0..n {
            let t = class_idx(&y[i]);
            let a1 = class_idx(&y1[i]);
            let a2 = class_idx(&y2[i]);
            if let (Some(t), Some(a1)) = (t, a1) {
                e1 += cost(t, a1);
            }
            if let (Some(t), Some(a2)) = (t, a2) {
                e2 += cost(t, a2);
            }
            if let (Some(t), Some(a1), Some(a2)) = (t, a1, a2) {
                mobs[t * k * k + a2 * k + a1] += 1.0;
            }
        }
        err1 = e1 / n as f64;
        err2 = e2 / n as f64;
        // fit of the cell counts
        let kc = k * k * k;
        if cm.data.iter().all(|&c| c == 0.0) {
            p = 1.0;
        } else {
            let mut dcost = vec![0.0; kc];
            for kk in 0..k {
                for j in 0..k {
                    for i in 0..k {
                        dcost[kk * k * k + j * k + i] = cost(kk, i) - cost(kk, j);
                    }
                }
            }
            let mobs1: Vec<f64> = mobs.iter().map(|v| v + 1.0).collect();
            if costtest == "chisquare" {
                return rt("You need an Optimization Toolbox license to set the 'Cost' parameter to 'chisquare'.");
            }
            let nn: f64 = mobs1.iter().sum();
            let maxd = dcost.iter().fold(0.0f64, |m, d| m.max(d.abs()));
            let mut maxl = nn / maxd;
            maxl *= 1.0 - f64::EPSILON;
            let f = |lam: f64| -> f64 { (0..kc).map(|i| mobs1[i] * dcost[i] / (nn + lam * dcost[i])).sum() };
            let (mut lo, mut hi) = (-maxl, maxl);
            let (flo, fhi) = (f(lo), f(hi));
            let fzero_fail = |why: &str| -> Ret { rt(format!("Unable to solve for cell probabilities under the null hypothesis. FZERO fails with the following message: {why}. Try using a cost matrix with off-diagonal elements closer to 1 or set 'CostTest' to 'chisquare'.")) };
            if !flo.is_finite() || !fhi.is_finite() {
                return fzero_fail("Function values at interval endpoints must be finite and real.");
            }
            if flo * fhi > 0.0 {
                return fzero_fail("Function values at interval endpoints must differ in sign.");
            }
            let mut flo = flo;
            for _ in 0..300 {
                let mid = 0.5 * (lo + hi);
                let fm = f(mid);
                if fm == 0.0 {
                    lo = mid;
                    hi = mid;
                    break;
                }
                if (fm > 0.0) == (flo > 0.0) {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
                if (hi - lo).abs() <= 2.0 * f64::EPSILON * lo.abs().max(hi.abs()) {
                    break;
                }
            }
            let lambda = 0.5 * (lo + hi);
            let chisq: f64 = (0..kc).map(|i| 2.0 * mobs1[i] * (lambda * dcost[i] / nn).ln_1p()).sum();
            p = chi2cdf(chisq, 1.0, true);
        }
        h = p < alpha;
    } else {
        let good1: Vec<bool> = (0..n).map(|i| y1[i] == y[i]).collect();
        let good2: Vec<bool> = (0..n).map(|i| y2[i] == y[i]).collect();
        err1 = good1.iter().filter(|&&g| !g).count() as f64 / n as f64;
        err2 = good2.iter().filter(|&&g| !g).count() as f64 / n as f64;
        let n01 = (0..n).filter(|&i| !good1[i] && good2[i]).count() as f64;
        let n10 = (0..n).filter(|&i| good1[i] && !good2[i]).count() as f64;
        let early = (n01 == 0.0 && n10 == 0.0) || (n01 == n10 && alternative == "unequal");
        if early {
            let mut o = vec![Value::bool(false), Value::num(1.0)];
            if nargout > 2 {
                o.push(Value::num(err1));
            }
            if nargout > 3 {
                o.push(Value::num(err2));
            }
            return Ok(Some(o));
        }
        let m = n01 + n10;
        let mut binopdf = |kk: f64| -> f64 {
            let r = it.call_function("binopdf", vec![Value::num(kk), Value::num(m), Value::num(0.5)], 1);
            match r {
                Ok(v) => v[0].scalar().unwrap_or(f64::NAN),
                Err(_) => f64::NAN,
            }
        };
        let cdf = |kk: f64| -> f64 { if kk < 0.0 { 0.0 } else { binocdf(kk, m, 0.5) } };
        p = match (mode, alternative) {
            ("asymptotic", "unequal") => chi2cdf((n01 - n10) * (n01 - n10) / m, 1.0, true),
            ("asymptotic", "less") => normcdf_std((n10 - n01) / m.sqrt()),
            ("asymptotic", _) => normcdf_std((n01 - n10) / m.sqrt()),
            ("exact", "unequal") => 2.0 * cdf(n01.min(n10)),
            ("exact", "less") => cdf(n10),
            ("exact", _) => cdf(n01),
            (_, "unequal") => 2.0 * cdf(n01.min(n10) - 1.0) + binopdf(n01.min(n10)),
            (_, "less") => cdf(n10 - 1.0) + 0.5 * binopdf(n10),
            _ => cdf(n01 - 1.0) + 0.5 * binopdf(n01),
        };
        h = p < alpha;
    }
    let mut out = vec![Value::bool(h), Value::num(p)];
    if nargout > 2 {
        out.push(Value::num(err1));
    }
    if nargout > 3 {
        out.push(Value::num(err2));
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// invpred

fn invpred(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    let vecof = |v: &Value| -> Option<Vec<f64>> {
        match v {
            Value::Num(m) if m.im.is_none() && m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1) => Some(m.data.clone()),
            _ => None,
        }
    };
    let Some(x0) = vecof(&a[0]) else { return rt("X must be a vector of real values.") };
    let Some(y0v) = vecof(&a[1]) else { return rt("Y must be a vector of real values.") };
    if x0.len() != y0v.len() {
        return rt("X and Y must have the same length.");
    }
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for (xi, yi) in x0.iter().zip(&y0v) {
        if !xi.is_nan() && !yi.is_nan() {
            x.push(*xi);
            y.push(*yi);
        }
    }
    let Value::Num(y0m) = &a[2] else { return rt("Y0 must be an array of real values.") };
    if y0m.im.is_some() {
        return rt("Y0 must be an array of real values.");
    }
    let opts = parse_pairs(&a[3..], &["alpha", "predopt"])?;
    let alpha = match opts[0] {
        None => 0.05,
        Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() == 1 && m.im.is_none() && m.data[0] > 0.0 && m.data[0] < 1.0 => m.data[0],
        Some(_) => return rt("ALPHA must be a scalar strictly between 0 and 1."),
    };
    let doobs = match opts[1] {
        None => true,
        Some(v) if v.is_text() => {
            let t = v.text().unwrap_or_default().to_lowercase();
            let choices = ["curve", "observation"];
            let hits: Vec<usize> = choices.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.starts_with(&t)).map(|(i, _)| i).collect();
            if hits.len() != 1 {
                return rt("'PredOpt' value must be 'curve' or 'observation'.");
            }
            hits[0] == 1
        }
        Some(_) => return rt("'PredOpt' value must be 'curve' or 'observation'."),
    };
    let nn = x.len() as f64;
    let mx = x.iter().sum::<f64>() / nn;
    let xc: Vec<f64> = x.iter().map(|v| v - mx).collect();
    let normx = dot(&xc, &xc);
    if normx == 0.0 {
        return rt("Cannot compute inverse predictions if X is constant.");
    }
    let my = y.iter().sum::<f64>() / nn;
    let slope = (0..x.len()).map(|i| xc[i] * y[i]).sum::<f64>() / normx;
    let icpt = my;
    let normr = ((0..x.len()).map(|i| (y[i] - (slope * xc[i] + icpt)).powi(2)).sum::<f64>()).sqrt();
    let single_out = matches!(&a[0], Value::Num(m) if is_single(m)) || matches!(&a[1], Value::Num(m) if is_single(m)) || is_single(y0m);
    let n0 = y0m.data.len();
    let x0v: Vec<f64> = y0m.data.iter().map(|&yy| if slope == 0.0 { f64::NAN } else { (yy - icpt) / slope }).collect();
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(y0m.dims.clone(), x0v.iter().map(|v| v + mx).collect()), single_out))];
    if nargout > 1 {
        let df = nn - 2.0;
        let tq = -tinv(alpha / 2.0, df);
        let s = if df > 0.0 { normr / df.sqrt() } else { f64::NAN };
        let varfact = if doobs { 1.0 + 1.0 / nn } else { 1.0 / nn };
        let mut lo = vec![0.0; n0];
        let mut up = vec![0.0; n0];
        for i in 0..n0 {
            let yy = y0m.data[i];
            let aa = slope * slope - tq * tq * s * s / normx;
            let b = -2.0 * slope * (yy - my);
            let c = (yy - my).powi(2) - tq * tq * s * s * varfact;
            let bsq = b * b;
            let mut discr = bsq - 4.0 * aa * c;
            if discr < 0.0 && discr > -100.0 * ulp(bsq) {
                discr = 0.0;
            }
            let x0i = x0v[i];
            let complex = discr < 0.0;
            let sgn = if b > 0.0 { 1.0 } else if b < 0.0 { -1.0 } else { 0.0 };
            let (mut d1, mut d2) = (f64::NAN, f64::NAN);
            if !complex {
                let q = -0.5 * (b + sgn * discr.sqrt());
                d1 = q / aa;
                d2 = c / q;
                if (d1 - d2).abs() <= ulp(normx) {
                    d1 = x0i;
                    d2 = x0i;
                }
                if d1 > d2 {
                    std::mem::swap(&mut d1, &mut d2);
                }
            }
            let (mut xlo, mut xup) = (d1, d2);
            if complex {
                xlo = f64::NEG_INFINITY;
                xup = f64::INFINITY;
            } else {
                if x0i < d1 {
                    xup = d1;
                    xlo = f64::NEG_INFINITY;
                }
                if x0i > d2 {
                    xlo = d2;
                    xup = f64::INFINITY;
                }
            }
            lo[i] = x0i - xlo;
            up[i] = xup - x0i;
        }
        outs.push(Value::Num(finish(Matrix::from_dims(y0m.dims.clone(), lo), single_out)));
        if nargout > 2 {
            outs.push(Value::Num(finish(Matrix::from_dims(y0m.dims.clone(), up), single_out)));
        }
    }
    Ok(Some(outs))
}

// ---------------------------------------------------------------------------------------------
// partialcorri

/// An orthonormal basis of the column space (pivoted Gram-Schmidt), and its rank.
fn orth_basis(cols: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let mut a: Vec<Vec<f64>> = cols.to_vec();
    let maxnorm = a.iter().map(|c| dot(c, c).sqrt()).fold(0.0, f64::max);
    let tol = n.max(cols.len()) as f64 * f64::EPSILON * maxnorm;
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut used = vec![false; a.len()];
    loop {
        let mut best = None;
        let mut bn = tol;
        for (j, c) in a.iter().enumerate() {
            if used[j] {
                continue;
            }
            let nn = dot(c, c).sqrt();
            if nn > bn {
                bn = nn;
                best = Some(j);
            }
        }
        let Some(j) = best else { break };
        used[j] = true;
        let qj: Vec<f64> = a[j].iter().map(|v| v / bn).collect();
        for (k, c) in a.iter_mut().enumerate() {
            if used[k] {
                continue;
            }
            for _ in 0..2 {
                let d = dot(&qj, c);
                for i in 0..n {
                    c[i] -= d * qj[i];
                }
            }
        }
        q.push(qj);
    }
    q
}

/// The part of `y` outside the span of `q`; a residual smaller than the rounding level of the
/// fit counts as zero.
fn residual(y: &[f64], q: &[Vec<f64>], rank: usize) -> Vec<f64> {
    let mut r = y.to_vec();
    for _ in 0..2 {
        for qj in q {
            let d = dot(qj, &r);
            for i in 0..r.len() {
                r[i] -= d * qj[i];
            }
        }
    }
    let tol = y.len().max(rank) as f64 * f64::EPSILON * dot(y, y).sqrt();
    if dot(&r, &r).sqrt() < tol {
        r.iter_mut().for_each(|v| *v = 0.0);
    }
    r
}

fn tied_ranks(c: &[f64]) -> Vec<f64> {
    let n = c.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| c[i].partial_cmp(&c[j]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && c[idx[j + 1]] == c[idx[i]] {
            j += 1;
        }
        let rank = (i + j) as f64 / 2.0 + 1.0;
        for k in i..=j {
            r[idx[k]] = rank;
        }
        i = j + 1;
    }
    r
}

fn mat_cols(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1]);
    (0..c).map(|j| m.data[j * r..(j + 1) * r].to_vec()).collect()
}

fn partialcorri(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let (Value::Num(ym), Value::Num(xm)) = (&a[0], &a[1]) else { return rt("All input arrays must be matrices.") };
    if ym.dims.len() != 2 || xm.dims.len() != 2 || ym.kind == Kind::Char || xm.kind == Kind::Char {
        return rt("All input arrays must be matrices.");
    }
    let nobs = ym.dims[0];
    if xm.dims[0] != nobs {
        return rt("All input arrays must have the same number of rows.");
    }
    let mut rest = &a[2..];
    let mut zcols: Option<Vec<Vec<f64>>> = None;
    if let Some(Value::Num(zm)) = rest.first() {
        if zm.kind != Kind::Char {
            if zm.dims.len() != 2 {
                return rt("All input arrays must be matrices.");
            }
            if zm.dims[0] != nobs {
                return rt("All input arrays must have the same number of rows.");
            }
            zcols = Some(mat_cols(zm));
            rest = &rest[1..];
        }
    }
    if ym.is_empty() {
        return rt("Input matrices cannot be empty.");
    }
    let opts = parse_pairs(rest, &["type", "rows", "tail"])?;
    let rows = match opts[1] {
        None => 'a',
        Some(v) => ['a', 'c', 'p'][super::tests3::param_val(v, &["all", "complete", "pairwise"], "'Rows'", true)?],
    };
    let ty = match opts[0] {
        None => 'p',
        Some(v) => {
            let t = if v.is_text() { v.text().unwrap_or_default().to_lowercase() } else { String::new() };
            let choices = ["pearson", "kendall", "spearman"];
            let hits: Vec<usize> = choices.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.starts_with(&t)).map(|(i, _)| i).collect();
            if hits.len() != 1 {
                return rt("The 'type' parameter value must be 'Pearson' or 'Spearman'.");
            }
            ['p', 'k', 's'][hits[0]]
        }
    };
    if ty == 'k' {
        return rt("Cannot compute Kendall's partial rank correlation.");
    }
    let tail = match opts[2] {
        None => 'b',
        Some(v) if v.is_text() => {
            let t = v.text().unwrap_or_default().to_lowercase();
            let long = ["left", "both", "right"];
            let short = ["lt", "ne", "gt"];
            let mut hits: Vec<usize> = long.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.starts_with(&t)).map(|(i, _)| i).collect();
            if hits.is_empty() {
                hits = short.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.starts_with(&t)).map(|(i, _)| i).collect();
            }
            if hits.len() != 1 {
                return rt("Tail must be 'both', 'right', or 'left'.");
            }
            long[hits[0]].chars().next().unwrap()
        }
        Some(_) => return rt("Tail must be 'both', 'right', or 'left'."),
    };
    let two_stage = zcols.is_some();
    let mut ycols = mat_cols(ym);
    let mut xcols = mat_cols(xm);
    let mut zc = zcols.unwrap_or_default();
    let mut n = nobs;
    let (p, q) = (ycols.len(), xcols.len());
    if rows != 'a' {
        let mut missing: Vec<bool> = (0..n).map(|i| xcols.iter().any(|c| c[i].is_nan())).collect();
        if two_stage {
            for i in 0..n {
                missing[i] = missing[i] || zc.iter().any(|c| c[i].is_nan());
            }
        }
        let keepc = |cols: &mut Vec<Vec<f64>>, miss: &[bool]| {
            for c in cols.iter_mut() {
                *c = c.iter().zip(miss).filter(|(_, &m)| !m).map(|(v, _)| *v).collect();
            }
        };
        keepc(&mut xcols, &missing);
        keepc(&mut ycols, &missing);
        if two_stage {
            keepc(&mut zc, &missing);
        }
        n = missing.iter().filter(|&&m| !m).count();
        if rows == 'c' {
            let miss2: Vec<bool> = (0..n).map(|i| ycols.iter().any(|c| c[i].is_nan())).collect();
            keepc(&mut ycols, &miss2);
            keepc(&mut xcols, &miss2);
            if two_stage {
                keepc(&mut zc, &miss2);
            }
            n = miss2.iter().filter(|&&m| !m).count();
        }
    }
    let pairwise = rows == 'p';
    if ty == 's' {
        for c in xcols.iter_mut().chain(ycols.iter_mut()).chain(zc.iter_mut()) {
            *c = tied_ranks(c);
        }
    }
    let nan_out = |p: usize, q: usize| -> Value { Value::Num(Matrix::from_dims(vec![p, q], vec![f64::NAN; p * q])) };
    let anynan = |cols: &[Vec<f64>]| cols.iter().any(|c| c.iter().any(|v| v.is_nan()));
    if anynan(&xcols) || anynan(&zc) {
        let mut o = vec![nan_out(p, q)];
        if nargout > 1 {
            o.push(nan_out(p, q));
        }
        return Ok(Some(o));
    }
    if n == 0 || ycols.iter().all(|c| c.is_empty()) {
        return rt("Warning-level condition: no rows are left after removing missing values.");
    }
    let mut coef = vec![f64::NAN; p * q];
    let mut pval = vec![f64::NAN; p * q];
    let groups: Vec<(Vec<usize>, Vec<bool>)> = if !(pairwise && anynan(&ycols)) {
        vec![((0..p).collect(), vec![true; n])]
    } else {
        (0..p).map(|i| (vec![i], ycols[i].iter().map(|v| !v.is_nan()).collect())).collect()
    };
    for (ycs, mask) in groups {
        let sub = |c: &Vec<f64>| -> Vec<f64> { c.iter().zip(&mask).filter(|(_, &m)| m).map(|(v, _)| *v).collect() };
        let nsub = mask.iter().filter(|&&m| m).count();
        if nsub == 0 {
            continue;
        }
        let yc: Vec<Vec<f64>> = ycs.iter().map(|&i| sub(&ycols[i])).collect();
        let mut xs: Vec<Vec<f64>> = xcols.iter().map(&sub).collect();
        let zs: Vec<Vec<f64>> = zc.iter().map(&sub).collect();
        let mut ys = yc;
        let demean = |c: &mut Vec<f64>| {
            let m = c.iter().sum::<f64>() / c.len() as f64;
            for v in c.iter_mut() {
                *v -= m;
            }
        };
        for c in xs.iter_mut().chain(ys.iter_mut()) {
            demean(c);
        }
        let mut zs = zs;
        for c in zs.iter_mut() {
            demean(c);
        }
        let mut rank_z = 0;
        if two_stage {
            let qz = orth_basis(&zs, nsub);
            rank_z = qz.len();
            for c in xs.iter_mut().chain(ys.iter_mut()) {
                *c = residual(c, &qz, rank_z);
                demean(c);
            }
        }
        let rank_x = orth_basis(&xs, nsub).len();
        for (yi, ycol) in ys.iter().enumerate() {
            for j in 0..q {
                let others: Vec<Vec<f64>> = (0..q).filter(|&k| k != j).map(|k| xs[k].clone()).collect();
                let qo = orth_basis(&others, nsub);
                let rk = rank_x;
                let ry = residual(ycol, &qo, rk);
                let rx = residual(&xs[j], &qo, rk);
                let nrx = dot(&rx, &rx).sqrt();
                let c0 = dot(&ry, &rx) / nrx;
                let c = c0 / dot(&ry, &ry).sqrt();
                let row = ycs[yi];
                coef[j * p + row] = c;
                if nargout > 1 {
                    let df = (nsub as f64 - 1.0 - rank_x as f64 - rank_z as f64).max(0.0);
                    let mut t = if c > 0.0 { f64::INFINITY } else if c < 0.0 { f64::NEG_INFINITY } else { 0.0 * f64::INFINITY };
                    if c.abs() < 1.0 {
                        t = c / (1.0 - c * c).sqrt();
                    }
                    t *= df.sqrt();
                    pval[j * p + row] = match tail {
                        'b' => 2.0 * tcdf(-t.abs(), df),
                        'r' => tcdf(-t, df),
                        _ => tcdf(t, df),
                    };
                }
            }
        }
    }
    let mut o = vec![Value::Num(Matrix::from_dims(vec![p, q], coef))];
    if nargout > 1 {
        o.push(Value::Num(Matrix::from_dims(vec![p, q], pval)));
    }
    Ok(Some(o))
}
