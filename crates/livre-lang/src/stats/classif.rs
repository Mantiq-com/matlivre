//! `classify` (discriminant analysis), `perfcurve` (ROC and other performance curves) and
//! `robustfit` (iteratively reweighted robust regression).

use super::dist_matrix::chol_upper;
use super::regress::{qr_pivot, solve_upper};
use super::tests2::{labels_of, Labels};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["classify", "perfcurve", "robustfit"];

fn matches_prefix(given: &str, full: &str) -> bool {
    !given.is_empty() && full.to_lowercase().starts_with(&given.to_lowercase())
}

fn rows_of(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1]);
    (0..r).map(|i| (0..c).map(|j| m.data[j * r + i]).collect()).collect()
}

fn col(v: Vec<f64>, single: bool) -> Value {
    let n = v.len();
    Value::Num(finish(Matrix::from_dims(vec![n, 1], v), single))
}

fn mat(rows: usize, cols: usize, d: Vec<f64>, single: bool) -> Value {
    Value::Num(finish(Matrix::from_dims(vec![rows, cols], d), single))
}

// ---- classify ------------------------------------------------------------------------------------

/// `z` solved against the transposed upper Cholesky factor: the vector `a` with `|a|^2 = z S^-1 z'`.
fn whiten(r: &[Vec<f64>], z: &[f64]) -> Vec<f64> {
    let d = z.len();
    let mut y = vec![0.0; d];
    for j in 0..d {
        let mut s = z[j];
        for i in 0..j {
            s -= y[i] * r[i][j];
        }
        y[j] = s / r[j][j];
    }
    y
}

fn classify(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt("Requires at least three arguments.");
    }
    let (Value::Num(sample0), Value::Num(training0)) = (&a[0], &a[1]) else { return rt("SAMPLE and TRAINING must be numeric matrices.") };
    let single = is_single(sample0) || is_single(training0);
    let gr = it.call_function("grp2idx", vec![a[2].clone()], 3)?;
    let Value::Num(gidx0) = &gr[0] else { return rt("grp2idx failed") };
    let gnames_n = match &gr[1] {
        Value::Cell(c) => c.data.len(),
        _ => 0,
    };
    let ngroups = gnames_n;
    let train_rows = rows_of(training0);
    let (mut tr, mut gindex): (Vec<Vec<f64>>, Vec<usize>) = (vec![], vec![]);
    if gidx0.numel() != train_rows.len() {
        return rt("The length of GROUP must equal the number of rows in TRAINING.");
    }
    for (i, &g) in gidx0.data.iter().enumerate() {
        if !g.is_nan() {
            tr.push(train_rows[i].clone());
            gindex.push(g as usize - 1);
        }
    }
    let n = tr.len();
    let d = training0.dims[1];
    let sample_rows: Vec<Vec<f64>> = if sample0.is_empty() {
        vec![]
    } else if sample0.dims[1] != d {
        return rt("SAMPLE and TRAINING must have the same number of columns.");
    } else {
        rows_of(sample0)
    };
    let m = sample_rows.len();
    let mut gsize = vec![0usize; ngroups];
    for &g in &gindex {
        gsize[g] += 1;
    }
    let nonempty: Vec<usize> = (0..ngroups).filter(|&k| gsize[k] > 0).collect();
    let nused = nonempty.len();
    // type
    let types = ["linear", "quadratic", "diaglinear", "diagquadratic", "mahalanobis"];
    let ty = match a.get(3) {
        None => "linear",
        Some(v) if matches!(v, Value::Num(mm) if mm.is_empty()) => "linear",
        Some(v) if v.is_text() => {
            let t = v.text()?;
            let hits: Vec<&&str> = types.iter().filter(|n| matches_prefix(&t, n)).collect();
            if hits.len() == 1 {
                hits[0]
            } else {
                return rt(format!("'{t}' is not a valid value for the TYPE argument.  Valid values are: 'linear', 'quadratic', 'diaglinear', 'diagquadratic', 'mahalanobis'."));
            }
        }
        _ => return rt("TYPE must be a character vector or string scalar."),
    };
    // prior
    let mut prior: Vec<f64> = vec![0.0; ngroups];
    match a.get(4) {
        None => nonempty.iter().for_each(|&k| prior[k] = 1.0 / nused as f64),
        Some(v) if matches!(v, Value::Num(mm) if mm.is_empty()) => nonempty.iter().for_each(|&k| prior[k] = 1.0 / nused as f64),
        Some(v) if v.is_text() && matches_prefix(&v.text()?, "empirical") => {
            for k in 0..ngroups {
                prior[k] = gsize[k] as f64 / n as f64;
            }
        }
        Some(Value::Num(p)) if p.kind != Kind::Char => {
            if p.dims.iter().min().copied().unwrap_or(0) != 1 || p.numel() != ngroups {
                return rt("PRIOR must be a vector one element for each group.");
            }
            if p.data.iter().any(|&v| v < 0.0) {
                return rt("PRIOR cannot contain negative values.");
            }
            let mut pv = p.data.clone();
            for k in 0..ngroups {
                if gsize[k] == 0 {
                    pv[k] = 0.0;
                }
            }
            let s: f64 = pv.iter().sum();
            prior = pv.iter().map(|v| v / s).collect();
        }
        _ => return rt("Prior must be a vector, a structure, or the character vector or string scalar 'empirical'."),
    }
    let all_rows: Vec<Vec<f64>> = if nargout > 1 { sample_rows.iter().chain(tr.iter()).cloned().collect() } else { sample_rows.clone() };
    let mm = all_rows.len();
    let mut means = vec![vec![f64::NAN; d]; ngroups];
    for &k in &nonempty {
        for j in 0..d {
            means[k][j] = (0..n).filter(|&i| gindex[i] == k).map(|i| tr[i][j]).sum::<f64>() / gsize[k] as f64;
        }
    }
    let mut dmat = vec![vec![f64::NAN; ngroups]; mm];
    let cov_chol = |rows: &[Vec<f64>], cen: &dyn Fn(usize) -> Vec<f64>, denom: f64| -> Option<Vec<Vec<f64>>> {
        let mut c = vec![vec![0.0; d]; d];
        for (i, _) in rows.iter().enumerate() {
            let z = cen(i);
            for p in 0..d {
                for q in 0..d {
                    c[p][q] += z[p] * z[q];
                }
            }
        }
        for p in 0..d {
            for q in 0..d {
                c[p][q] /= denom;
            }
        }
        chol_upper(&c)
    };
    let mut logdet_k = vec![0.0; ngroups];
    let mut chol_k: Vec<Vec<Vec<f64>>> = vec![vec![]; ngroups];
    let mut sd_k: Vec<Vec<f64>> = vec![vec![]; ngroups];
    match ty {
        "linear" => {
            if n <= nused {
                return rt("TRAINING must have more observations than the number of groups.");
            }
            let Some(r) = cov_chol(&tr, &|i| (0..d).map(|j| tr[i][j] - means[gindex[i]][j]).collect(), (n - nused) as f64) else { return rt("The pooled covariance matrix of TRAINING must be positive definite.") };
            let logdet = 2.0 * (0..d).map(|i| r[i][i].ln()).sum::<f64>();
            for &k in &nonempty {
                for (i, row) in all_rows.iter().enumerate() {
                    let z: Vec<f64> = (0..d).map(|j| row[j] - means[k][j]).collect();
                    let av = whiten(&r, &z);
                    dmat[i][k] = prior[k].ln() - 0.5 * (av.iter().map(|v| v * v).sum::<f64>() + logdet);
                }
            }
        }
        "diaglinear" => {
            if n <= nused {
                return rt("TRAINING must have more observations than the number of groups.");
            }
            let s: Vec<f64> = (0..d).map(|j| ((0..n).map(|i| (tr[i][j] - means[gindex[i]][j]).powi(2)).sum::<f64>() / (n - nused) as f64).sqrt()).collect();
            let smax = s.iter().cloned().fold(0.0, f64::max);
            if s.iter().any(|&v| v <= n as f64 * crate::stats::util::eps_of(smax)) {
                return rt("The pooled variance of each column of TRAINING must be positive.");
            }
            let logdet = 2.0 * s.iter().map(|v| v.ln()).sum::<f64>();
            for &k in &nonempty {
                for (i, row) in all_rows.iter().enumerate() {
                    let q: f64 = (0..d).map(|j| ((row[j] - means[k][j]) / s[j]).powi(2)).sum();
                    dmat[i][k] = prior[k].ln() - 0.5 * (q + logdet);
                }
            }
        }
        "quadratic" | "mahalanobis" => {
            if gsize.iter().any(|&s| s == 1) {
                return rt("Each group in TRAINING must have at least two observations.");
            }
            for &k in &nonempty {
                let idx: Vec<usize> = (0..n).filter(|&i| gindex[i] == k).collect();
                let Some(r) = cov_chol(&vec![vec![]; idx.len()], &|t| (0..d).map(|j| tr[idx[t]][j] - means[k][j]).collect(), (gsize[k] - 1) as f64) else { return rt("The covariance matrix of each group of TRAINING must be positive definite.") };
                logdet_k[k] = 2.0 * (0..d).map(|i| r[i][i].ln()).sum::<f64>();
                for (i, row) in all_rows.iter().enumerate() {
                    let z: Vec<f64> = (0..d).map(|j| row[j] - means[k][j]).collect();
                    let q: f64 = whiten(&r, &z).iter().map(|v| v * v).sum();
                    dmat[i][k] = if ty == "quadratic" { prior[k].ln() - 0.5 * (q + logdet_k[k]) } else { -q };
                }
                chol_k[k] = r;
            }
        }
        _ => {
            if gsize.iter().any(|&s| s == 1) {
                return rt("Each group in TRAINING must have at least two observations.");
            }
            for &k in &nonempty {
                let s: Vec<f64> = (0..d).map(|j| ((0..n).filter(|&i| gindex[i] == k).map(|i| (tr[i][j] - means[k][j]).powi(2)).sum::<f64>() / (gsize[k] - 1) as f64).sqrt()).collect();
                let smax = s.iter().cloned().fold(0.0, f64::max);
                if s.iter().any(|&v| v <= gsize[k] as f64 * crate::stats::util::eps_of(smax)) {
                    return rt("The variance of each column within each group of TRAINING must be positive.");
                }
                logdet_k[k] = 2.0 * s.iter().map(|v| v.ln()).sum::<f64>();
                for (i, row) in all_rows.iter().enumerate() {
                    let q: f64 = (0..d).map(|j| ((row[j] - means[k][j]) / s[j]).powi(2)).sum();
                    dmat[i][k] = prior[k].ln() - 0.5 * (q + logdet_k[k]);
                }
                sd_k[k] = s;
            }
        }
    }
    let _ = (&chol_k, &sd_k);
    let argmax = |row: &Vec<f64>| -> (usize, f64) {
        let mut best = (0usize, f64::NEG_INFINITY);
        let mut found = false;
        for (k, &v) in row.iter().enumerate() {
            if !v.is_nan() && (!found || v > best.1) {
                best = (k, v);
                found = true;
            }
        }
        if !found { (0, f64::NAN) } else { best }
    };
    let picks: Vec<(usize, f64)> = dmat.iter().map(argmax).collect();
    let glevels = |idx: &[usize]| -> Res<Value> {
        let src = it_levels(&gr[2], idx)?;
        Ok(src)
    };
    let out_idx: Vec<usize> = picks[..m].iter().map(|p| p.0).collect();
    let outclass = glevels(&out_idx)?;
    let mut outs = vec![outclass];
    if nargout > 1 {
        let mut e = vec![0.0; ngroups];
        for &k in &nonempty {
            let miss = (0..n).filter(|&i| gindex[i] == k && picks[m + i].0 != k).count();
            e[k] = miss as f64 / gsize[k] as f64;
        }
        let err: f64 = (0..ngroups).map(|k| prior[k] * e[k]).sum();
        outs.push(Value::Num(finish(Matrix::scalar(err), single)));
    }
    if nargout > 2 {
        if ty == "mahalanobis" {
            outs.push(Value::empty());
            outs.push(Value::empty());
        } else {
            let mut post = vec![0.0; m * ngroups];
            let mut logp = vec![0.0; m];
            for i in 0..m {
                let maxd = picks[i].1;
                let p: Vec<f64> = (0..ngroups).map(|k| (dmat[i][k] - maxd).exp()).collect();
                let sump: f64 = p.iter().filter(|v| !v.is_nan()).sum();
                for k in 0..ngroups {
                    post[k * m + i] = p[k] / sump;
                }
                logp[i] = sump.ln() + maxd - 0.5 * d as f64 * (2.0 * std::f64::consts::PI).ln();
            }
            outs.push(mat(m, ngroups, post, single));
            outs.push(col(logp, single));
        }
    }
    if nargout > 4 {
        outs.push(classify_coeffs(it, ty, ngroups, &nonempty, &means, &prior, &gr[2], &logdet_k, &tr, &gindex, &gsize, d)?);
    }
    Ok(Some(outs))
}

fn it_levels(levels: &Value, idx: &[usize]) -> Res<Value> {
    match levels {
        Value::Num(m) => {
            let mut out = Matrix::from_dims(vec![idx.len(), 1], idx.iter().map(|&i| m.data[i]).collect());
            out.kind = m.kind;
            Ok(Value::Num(out))
        }
        Value::Cell(c) => Ok(Value::Cell(CellArray::from_dims(vec![idx.len(), 1], idx.iter().map(|&i| c.data[i].clone()).collect()))),
        Value::Str(s) => Ok(Value::Str(StrArray::col(idx.iter().map(|&i| s.data[i].clone().unwrap_or_default()).collect()))),
        other => Ok(other.clone()),
    }
}

#[allow(clippy::too_many_arguments)]
fn classify_coeffs(it: &mut Interpreter, ty: &str, ngroups: usize, nonempty: &[usize], means: &[Vec<f64>], prior: &[f64], levels: &Value, _logdet: &[f64], tr: &[Vec<f64>], gindex: &[usize], gsize: &[usize], d: usize) -> Res<Value> {
    let n = tr.len();
    let nused = nonempty.len();
    let level_at = |i: usize| -> Res<Value> { it_levels(levels, &[i]) };
    let mut st = StructArray { dims: vec![ngroups, ngroups], fields: vec!["type".into(), "name1".into(), "name2".into(), "const".into(), "linear".into()], elems: vec![] };
    for _ in 0..ngroups * ngroups {
        st.elems.push(vec![Value::chars(ty), Value::empty(), Value::empty(), Value::empty(), Value::empty()]);
    }
    if ty == "quadratic" || ty == "mahalanobis" || ty == "diagquadratic" {
        st.add_field("quadratic");
    }
    let inv = |it: &mut Interpreter, c: Vec<f64>, k: usize| -> Res<Vec<f64>> {
        let r = it.call_function("inv", vec![Value::Num(Matrix::from_dims(vec![k, k], c))], 1)?;
        let Value::Num(m) = &r[0] else { return rt("inv failed") };
        Ok(m.data.clone())
    };
    let mut pairs = Vec::new();
    for (ai, &i1) in nonempty.iter().enumerate() {
        for &i2 in &nonempty[ai + 1..] {
            pairs.push((i1, i2));
        }
    }
    // covariance pieces
    let pooled: Vec<f64> = if ty == "linear" || ty == "diaglinear" {
        let mut c = vec![0.0; d * d];
        for i in 0..n {
            for p in 0..d {
                for q in 0..d {
                    c[q * d + p] += (tr[i][p] - means[gindex[i]][p]) * (tr[i][q] - means[gindex[i]][q]);
                }
            }
        }
        c.iter().map(|v| v / (n - nused) as f64).collect()
    } else {
        vec![]
    };
    let group_cov = |k: usize| -> Vec<f64> {
        let idx: Vec<usize> = (0..n).filter(|&i| gindex[i] == k).collect();
        let mut c = vec![0.0; d * d];
        for &i in &idx {
            for p in 0..d {
                for q in 0..d {
                    c[q * d + p] += (tr[i][p] - means[k][p]) * (tr[i][q] - means[k][q]);
                }
            }
        }
        c.iter().map(|v| v / (gsize[k] - 1) as f64).collect()
    };
    let at = |i: usize, j: usize| i + ngroups * j;
    for &(i1, i2) in &pairs {
        let (mu1, mu2) = (&means[i1], &means[i2]);
        let (mut kk, lin, quad): (f64, Vec<f64>, Option<Vec<f64>>) = match ty {
            "linear" => {
                let si = inv(it, pooled.clone(), d)?;
                let diff: Vec<f64> = (0..d).map(|j| mu1[j] - mu2[j]).collect();
                let b: Vec<f64> = (0..d).map(|i| (0..d).map(|j| si[j * d + i] * diff[j]).sum()).collect();
                let kk = 0.5 * (0..d).map(|j| (mu1[j] + mu2[j]) * b[j]).sum::<f64>();
                (kk, b, None)
            }
            "diaglinear" => {
                let b: Vec<f64> = (0..d).map(|j| (mu1[j] - mu2[j]) / pooled[j * d + j]).collect();
                let kk = 0.5 * (0..d).map(|j| (mu1[j] + mu2[j]) * b[j]).sum::<f64>();
                (kk, b, None)
            }
            _ => {
                let (c1, c2) = (group_cov(i1), group_cov(i2));
                let (s1, s2) = if ty == "diagquadratic" {
                    ((0..d).map(|j| 1.0 / c1[j * d + j]).collect::<Vec<f64>>(), (0..d).map(|j| 1.0 / c2[j * d + j]).collect::<Vec<f64>>())
                } else {
                    (inv(it, c1.clone(), d)?, inv(it, c2.clone(), d)?)
                };
                let (ld1, ld2) = (_logdet[i1], _logdet[i2]);
                if ty == "diagquadratic" {
                    let k0 = 0.5 * ((0..d).map(|j| s1[j] * mu1[j] * mu1[j]).sum::<f64>() - (0..d).map(|j| s2[j] * mu2[j] * mu2[j]).sum::<f64>()) + 0.5 * (ld1 - ld2);
                    let l: Vec<f64> = (0..d).map(|j| s1[j] * mu1[j] - s2[j] * mu2[j]).collect();
                    let q: Vec<f64> = (0..d * d).map(|idx| if idx % d == idx / d { -0.5 * (s1[idx % d] - s2[idx % d]) } else { 0.0 }).collect();
                    (k0, l, Some(q))
                } else {
                    let quad_form = |s: &Vec<f64>, mu: &Vec<f64>| -> f64 { (0..d).map(|i| (0..d).map(|j| mu[i] * s[j * d + i] * mu[j]).sum::<f64>()).sum() };
                    let mut k0 = 0.5 * (quad_form(&s1, mu1) - quad_form(&s2, mu2));
                    if ty != "mahalanobis" {
                        k0 += 0.5 * (ld1 - ld2);
                    }
                    let l: Vec<f64> = (0..d).map(|i| (0..d).map(|j| s1[j * d + i] * mu1[j] - s2[j * d + i] * mu2[j]).sum()).collect();
                    let q: Vec<f64> = (0..d * d).map(|idx| -0.5 * (s1[idx] - s2[idx])).collect();
                    (k0, l, Some(q))
                }
            }
        };
        if ty != "mahalanobis" {
            kk = kk - prior[i1].ln() + prior[i2].ln();
        }
        let lv = |v: &Vec<f64>| Value::Num(Matrix::from_dims(vec![d, 1], v.clone()));
        let neg = |v: &Vec<f64>| Value::Num(Matrix::from_dims(vec![d, 1], v.iter().map(|x| -x).collect()));
        let k_ij = at(i1, i2);
        let k_ji = at(i2, i1);
        st.set_at(k_ij, "name1", level_at(i1)?);
        st.set_at(k_ij, "name2", level_at(i2)?);
        st.set_at(k_ij, "const", Value::num(-kk));
        st.set_at(k_ij, "linear", lv(&lin));
        st.set_at(k_ji, "name1", level_at(i2)?);
        st.set_at(k_ji, "name2", level_at(i1)?);
        st.set_at(k_ji, "const", Value::num(kk));
        st.set_at(k_ji, "linear", neg(&lin));
        if let Some(q) = quad {
            st.set_at(k_ij, "quadratic", Value::Num(Matrix::from_dims(vec![d, d], q.clone())));
            st.set_at(k_ji, "quadratic", Value::Num(Matrix::from_dims(vec![d, d], q.iter().map(|x| -x).collect())));
        }
    }
    Ok(Value::Struct(st))
}

// ---- perfcurve ----------------------------------------------------------------------------------

#[derive(Clone, Copy)]
enum Crit {
    Tp,
    Fn,
    Fp,
    Tn,
    TpFp,
    Rpp,
    Rnp,
    Accu,
    Tpr,
    Fnr,
    Fpr,
    Tnr,
    Ppv,
    Npv,
    Ecost,
}

fn crit_of(name: &str) -> Option<Crit> {
    Some(match name.to_lowercase().as_str() {
        "tp" | "truepositives" => Crit::Tp,
        "fn" | "falsenegatives" => Crit::Fn,
        "fp" | "falsepositives" => Crit::Fp,
        "tn" | "truenegatives" => Crit::Tn,
        "tp+fp" | "sumoftrueandfalsepositives" => Crit::TpFp,
        "rpp" | "rateofpositivepredictions" => Crit::Rpp,
        "rnp" | "rateofnegativepredictions" => Crit::Rnp,
        "accu" | "accuracy" => Crit::Accu,
        "tpr" | "sens" | "reca" | "truepositiverate" => Crit::Tpr,
        "fnr" | "miss" | "falsenegativerate" => Crit::Fnr,
        "fpr" | "fall" | "falsepositiverate" => Crit::Fpr,
        "tnr" | "spec" | "truenegativerate" => Crit::Tnr,
        "ppv" | "prec" | "positivepredictivevalue" => Crit::Ppv,
        "npv" | "negativepredictivevalue" => Crit::Npv,
        "ecost" | "expectedcost" => Crit::Ecost,
        _ => return None,
    })
}

fn crit_value(c: Crit, tp: f64, fn_: f64, fp: f64, tn: f64, s: [f64; 2], cost: &[f64; 4]) -> f64 {
    let denom = s[0] * (tp + fn_) + s[1] * (fp + tn);
    match c {
        Crit::Tp => s[0] * tp,
        Crit::Fn => s[0] * fn_,
        Crit::Fp => s[1] * fp,
        Crit::Tn => s[1] * tn,
        Crit::TpFp => s[0] * tp + s[1] * fp,
        Crit::Rpp => (s[0] * tp + s[1] * fp) / denom,
        Crit::Rnp => (s[0] * fn_ + s[1] * tn) / denom,
        Crit::Accu => (s[0] * tp + s[1] * tn) / denom,
        Crit::Tpr => tp / (tp + fn_),
        Crit::Fnr => fn_ / (tp + fn_),
        Crit::Fpr => fp / (fp + tn),
        Crit::Tnr => tn / (fp + tn),
        Crit::Ppv => s[0] * tp / (s[0] * tp + s[1] * fp),
        Crit::Npv => s[1] * tn / (s[0] * fn_ + s[1] * tn),
        // cost is [c11 c12; c21 c22] stored column-major
        Crit::Ecost => (s[0] * (tp * cost[0] + fn_ * cost[2]) + s[1] * (fp * cost[1] + tn * cost[3])) / denom,
    }
}

fn perfcurve(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt(format!("Invalid argument list. Function requires {} more input(s).", 3 - a.len()));
    }
    let names = ["NegClass", "XCrit", "YCrit", "XVals", "TVals", "UseNearest", "ProcessNaN", "Prior", "Cost", "Weights", "NBoot", "BootType", "BootArg", "Alpha", "Options"];
    let rest = &a[3..];
    if rest.len() % 2 != 0 {
        return rt("Invalid argument list. Name-value arguments must come in pairs.");
    }
    let (mut neg, mut xc, mut yc): (Option<Value>, String, String) = (None, "fpr".into(), "tpr".into());
    let (mut xvals, mut tvals): (Option<Value>, Option<Value>) = (None, None);
    let (mut nearest, mut add_false) = (true, false);
    let mut prior: Option<Value> = None;
    let mut cost = [0.0, 1.0, 1.0, 0.0];
    let mut weights: Option<Vec<f64>> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        let Some(nm) = names.iter().find(|n| n.eq_ignore_ascii_case(&key)) else {
            return rt(format!("Invalid argument name '{key}'. Name must be one of these:     {}", names.iter().map(|n| format!("'{n}'")).collect::<Vec<_>>().join("     ") + " "));
        };
        match *nm {
            "NegClass" => neg = Some(pair[1].clone()),
            "XCrit" => xc = pair[1].text().unwrap_or_default(),
            "YCrit" => yc = pair[1].text().unwrap_or_default(),
            "XVals" => xvals = Some(pair[1].clone()),
            "TVals" => tvals = Some(pair[1].clone()),
            "UseNearest" => {
                nearest = match &pair[1] {
                    v if v.is_text() => v.text()?.eq_ignore_ascii_case("on"),
                    v => v.scalar()? != 0.0,
                }
            }
            "ProcessNaN" => add_false = matches_prefix(&pair[1].text().unwrap_or_default(), "addtofalse"),
            "Prior" => prior = Some(pair[1].clone()),
            "Cost" => {
                let m = pair[1].as_matrix()?;
                if m.dims != [2, 2] {
                    return rt("COST must be a 2-by-2 matrix.");
                }
                cost = [m.data[0], m.data[1], m.data[2], m.data[3]];
            }
            "Weights" => weights = Some(pair[1].as_matrix()?.data.clone()),
            _ => {}
        }
    }
    let (Some(fx), Some(fy)) = (crit_of(&xc), crit_of(&yc)) else { return rt("You specified an unknown criterion for X or Y.") };
    let Value::Num(sc) = &a[1] else { return rt("SCORES must be a numeric vector.") };
    let single = is_single(sc);
    let scores = sc.data.clone();
    let labs = labels_of(&a[0])?;
    let nobs = match &labs {
        Labels::Num(v, _) => v.len(),
        Labels::Text(v, _) => v.len(),
    };
    if scores.len() != nobs {
        return rt("The size of scores does not match the size of labels.");
    }
    let w: Vec<f64> = match &weights {
        Some(wv) => {
            if wv.iter().any(|&v| v < 0.0) || wv.iter().all(|&v| v == 0.0) {
                return rt("WEIGHTS must be nonnegative with at least one positive value.");
            }
            if wv.len() != nobs {
                return rt("The size of weights does not match the size of scores.");
            }
            wv.clone()
        }
        None => vec![1.0; nobs],
    };
    // class membership: 0 = none, 1 = positive, 2.. = negative classes
    let same = |i: usize, which: &Value| -> bool {
        match (&labs, which) {
            (Labels::Num(v, _), Value::Num(m)) if m.kind != Kind::Char => m.data.contains(&v[i]),
            (Labels::Text(v, _), x) if x.is_text() => v[i].as_deref() == x.text().ok().as_deref(),
            (Labels::Text(v, _), Value::Cell(c)) => c.data.iter().any(|e| v[i].as_deref() == e.text().ok().as_deref()),
            _ => false,
        }
    };
    let distinct: Vec<usize> = {
        let mut seen: Vec<usize> = Vec::new();
        for i in 0..nobs {
            let dup = seen.iter().any(|&j| match &labs {
                Labels::Num(v, _) => v[j] == v[i],
                Labels::Text(v, _) => v[j] == v[i],
            });
            if !dup {
                seen.push(i);
            }
        }
        seen
    };
    let pos = &a[2];
    let pos_present = (0..nobs).any(|i| same(i, pos));
    if distinct.len() < 2 {
        return rt("Less than two classes are found in the array of true class labels.");
    }
    if !pos_present {
        return rt("Positive class is not found in the input data.");
    }
    let is_pos: Vec<bool> = (0..nobs).map(|i| same(i, pos)).collect();
    let is_neg: Vec<bool> = (0..nobs)
        .map(|i| match &neg {
            None => !is_pos[i],
            Some(nv) => same(i, nv),
        })
        .collect();
    // sort by descending score (stable), drop NaN scores and non-members
    let mut order: Vec<usize> = (0..nobs).collect();
    order.sort_by(|&i, &j| scores[j].partial_cmp(&scores[i]).unwrap_or(std::cmp::Ordering::Equal));
    let (mut nan_p, mut nan_n) = (0.0, 0.0);
    let mut srt: Vec<(f64, f64, f64)> = Vec::new(); // score, wpos, wneg
    for &i in &order {
        if !(is_pos[i] || is_neg[i]) || w[i] == 0.0 {
            continue;
        }
        if scores[i].is_nan() {
            if add_false {
                if is_pos[i] {
                    nan_p += w[i];
                } else {
                    nan_n += w[i];
                }
            }
            continue;
        }
        srt.push((scores[i], if is_pos[i] { w[i] } else { 0.0 }, if is_neg[i] { w[i] } else { 0.0 }));
    }
    // cumulative table: row 0 = FP of NaN scores, rows 1.. = cumulative weights, last row adds NaN positives
    let mut wcum: Vec<(f64, f64)> = Vec::new();
    let mut cp = 0.0;
    let mut cn = nan_n;
    wcum.push((0.0, cn));
    let mut sscores: Vec<f64> = Vec::new();
    for (s, wp, wn) in &srt {
        cp += wp;
        cn += wn;
        wcum.push((cp, cn));
        sscores.push(*s);
    }
    wcum.push((cp + nan_p, cn));
    // merge ties
    let eps = |v: f64| crate::stats::util::eps_of(v);
    let mut drop_rows: Vec<usize> = Vec::new();
    for i in 0..sscores.len().saturating_sub(1) {
        let (s1, s2) = (sscores[i], sscores[i + 1]);
        if s1 < s2 + eps(s1).max(eps(s2)) {
            drop_rows.push(i + 1);
        }
    }
    let keep_row = |r: usize| !drop_rows.contains(&r);
    let wcum: Vec<(f64, f64)> = wcum.iter().enumerate().filter(|(r, _)| keep_row(*r)).map(|(_, v)| *v).collect();
    let sscores: Vec<f64> = sscores.iter().enumerate().filter(|(i, _)| !drop_rows.contains(&(i + 1))).map(|(_, v)| *v).collect();
    let (wp_tot, wn_tot) = (wcum.last().unwrap().0, wcum.last().unwrap().1);
    // scale
    let scale: [f64; 2] = match &prior {
        None => [1.0, 1.0],
        Some(Value::Num(p)) if p.kind == Kind::Char => {
            let t = p.to_string_lossy();
            if matches_prefix(&t, "empirical") {
                [1.0, 1.0]
            } else if matches_prefix(&t, "uniform") {
                let (s0, s1) = (wn_tot, wp_tot);
                [s0 / (s0 + s1), s1 / (s0 + s1)]
            } else {
                return rt("PRIOR must be 'empirical', 'uniform', or a vector of two positive values.");
            }
        }
        Some(Value::Num(p)) if p.numel() == 2 => {
            if p.data.iter().any(|&v| v <= 0.0) {
                return rt("Priors must be positive.");
            }
            let (s0, s1) = (p.data[0] * wn_tot, p.data[1] * wp_tot);
            [s0 / (s0 + s1), s1 / (s0 + s1)]
        }
        _ => return rt("PRIOR must be 'empirical', 'uniform', or a vector of two positive values."),
    };
    let ndiv = wcum.len() - 1;
    let fx_all: Vec<(f64, f64, f64)> = (0..ndiv).map(|r| (wcum[r].0, wcum[r].1, 0.0)).collect();
    let eval = |c: Crit, r: usize| -> f64 {
        let (tp, fp) = (fx_all[r].0, fx_all[r].1);
        crit_value(c, tp, wp_tot - tp, fp, wn_tot - fp, scale, &cost)
    };
    // select divisions
    let mut div: Vec<usize> = (0..ndiv).collect();
    let mut tgiven: Option<Vec<f64>> = None;
    let xall = matches!(&xvals, Some(v) if v.is_text());
    let tall = matches!(&tvals, Some(v) if v.is_text());
    let mut xrange: Option<(f64, f64)> = None;
    let mut trange: Option<(f64, f64)> = None;
    if xvals.as_ref().is_some_and(|v| !matches!(v, Value::Num(m) if m.is_empty())) && tvals.as_ref().is_some_and(|v| !matches!(v, Value::Num(m) if m.is_empty())) {
        return rt("You cannot specify both XVals and TVals.");
    }
    if let (Some(Value::Num(tv)), false) = (&tvals, tall) {
        if !tv.is_empty() {
            let mut t = tv.data.clone();
            t.sort_by(|a, b| b.partial_cmp(a).unwrap());
            let mut d: Vec<usize> = Vec::new();
            let mut iall = 0usize;
            let mut cur = 1usize;
            for (idx, &tvv) in t.iter().enumerate() {
                let _ = idx;
                let found = (iall..sscores.len()).rev().find(|&j| sscores[j] >= tvv);
                if let Some(j) = found {
                    iall = j;
                    cur = j + 1;
                } else {
                    cur = 1.max(cur.min(1));
                }
                let mut dv = cur + 1;
                if sscores.is_empty() || tvv > sscores[0] {
                    dv = 1;
                }
                d.push(dv - 1);
            }
            // divisions are 0-based rows of the cumulative table
            if nearest {
                d.sort();
                d.dedup();
            }
            div = d;
            tgiven = Some(t);
            trange = Some((tv.data[0].min(tv.data[tv.numel() - 1]), tv.data[0].max(tv.data[tv.numel() - 1])));
        }
    } else if let (Some(Value::Num(xv)), false) = (&xvals, xall) {
        if !xv.is_empty() {
            let allv: Vec<f64> = (0..ndiv).map(|r| eval(fx, r)).collect();
            let finite: Vec<f64> = allv.iter().copied().filter(|v| !v.is_nan()).collect();
            let inc = !finite.windows(2).any(|w| w[0] > w[1]);
            let inc_s = if inc { 1.0 } else { -1.0 };
            let mut xs: Vec<f64> = xv.data.clone();
            xs.sort_by(|a, b| (inc_s * a).partial_cmp(&(inc_s * b)).unwrap());
            let mut d = vec![0usize; xs.len()];
            let mut iall = 0usize;
            for (i, &xvv) in xs.iter().enumerate() {
                if let Some(j) = (iall..ndiv).rev().find(|&j| inc_s * allv[j] <= inc_s * xvv) {
                    iall = j;
                    for t in d.iter_mut().skip(i) {
                        *t = iall;
                    }
                }
            }
            if nearest {
                if inc_s * allv[0] >= inc_s * xs[0] {
                    d.insert(0, 0);
                }
                d.sort();
                d.dedup();
            }
            div = d;
            xrange = Some((xv.data[0].min(xv.data[xv.numel() - 1]), xv.data[0].max(xv.data[xv.numel() - 1])));
        }
    }
    let xs: Vec<f64> = div.iter().map(|&r| eval(fx, r)).collect();
    let ys: Vec<f64> = div.iter().map(|&r| eval(fy, r)).collect();
    let xout: Vec<f64> = match (&xvals, nearest) {
        (Some(Value::Num(xv)), false) if !xv.is_empty() => xv.data.clone(),
        _ => xs.clone(),
    };
    let th = |r: usize| -> f64 {
        if sscores.is_empty() {
            f64::NAN
        } else if r == 0 {
            sscores[0]
        } else {
            sscores[r - 1]
        }
    };
    let tout: Vec<f64> = match (&tgiven, nearest) {
        (Some(t), false) => t.clone(),
        _ => div.iter().map(|&r| th(r)).collect(),
    };
    let mut outs = vec![col(xout, single), col(ys.clone(), single), col(tout, single)];
    if nargout > 3 {
        // area under the curve over all divisions (restricted to the X or T range)
        let xa: Vec<f64> = (0..ndiv).map(|r| eval(fx, r)).collect();
        let ya: Vec<f64> = (0..ndiv).map(|r| eval(fy, r)).collect();
        let ta: Vec<f64> = (0..ndiv).map(th_closure(&sscores)).collect();
        let inr: Vec<bool> = (0..ndiv)
            .map(|i| match (xrange, trange) {
                (Some((lo, hi)), _) => xa[i] >= lo && xa[i] <= hi,
                (_, Some((lo, hi))) => ta[i] >= lo && ta[i] <= hi,
                _ => true,
            })
            .collect();
        let mut xa2: Vec<f64> = (0..ndiv).filter(|&i| inr[i]).map(|i| xa[i]).collect();
        let mut ya2: Vec<f64> = (0..ndiv).filter(|&i| inr[i]).map(|i| ya[i]).collect();
        let (trim_first, trim_last) = (inr.first().copied().unwrap_or(false), inr.last().copied().unwrap_or(false));
        if trim_last && !xa2.is_empty() && (xa2[xa2.len() - 1].is_nan() || ya2[ya2.len() - 1].is_nan()) {
            xa2.pop();
            ya2.pop();
        }
        if trim_first && !xa2.is_empty() && (xa2[0].is_nan() || ya2[0].is_nan()) {
            xa2.remove(0);
            ya2.remove(0);
        }
        let auc = if xa2.len() < 2 { 0.0 } else { (0..xa2.len() - 1).map(|i| (xa2[i + 1] - xa2[i]) * (ya2[i + 1] + ya2[i])).sum::<f64>().mul_add(0.5, 0.0).abs() };
        outs.push(Value::Num(finish(Matrix::scalar(auc), single)));
    }
    if nargout > 4 {
        let isroc = matches!(fx, Crit::Fpr) && matches!(fy, Crit::Tpr);
        if isroc {
            let (wp, wn) = (scale[0] * wp_tot, scale[1] * wn_tot);
            let m = (cost[1] - cost[3]) / (cost[2] - cost[0]) * wn / wp;
            let mut best = 0usize;
            for i in 0..xs.len() {
                if xs[i] - ys[i] / m < xs[best] - ys[best] / m {
                    best = i;
                }
            }
            outs.push(Value::Num(finish(Matrix::from_dims(vec![1, 2], vec![xs[best], ys[best]]), single)));
        } else {
            outs.push(Value::Num(finish(Matrix::from_dims(vec![1, 2], vec![f64::NAN, f64::NAN]), single)));
        }
    }
    let _ = it;
    Ok(Some(outs))
}

fn th_closure(sscores: &[f64]) -> impl Fn(usize) -> f64 + '_ {
    move |r| {
        if sscores.is_empty() {
            f64::NAN
        } else if r == 0 {
            sscores[0]
        } else {
            sscores[r - 1]
        }
    }
}

// ---- robustfit ----------------------------------------------------------------------------------

#[derive(Clone)]
enum WFun {
    Builtin(&'static str),
    User(Value),
}

fn weight_default(name: &str) -> Option<(&'static str, f64)> {
    Some(match name.to_lowercase().as_str() {
        "andrews" => ("andrews", 1.339),
        "bisquare" => ("bisquare", 4.685),
        "cauchy" => ("cauchy", 2.385),
        "fair" => ("fair", 1.400),
        "huber" => ("huber", 1.345),
        "logistic" => ("logistic", 1.205),
        "ols" => ("ols", 1.0),
        "talwar" => ("talwar", 2.795),
        "welsch" => ("welsch", 2.985),
        _ => return None,
    })
}

fn apply_w(it: &mut Interpreter, f: &WFun, r: &[f64]) -> Res<Vec<f64>> {
    Ok(match f {
        WFun::Builtin(n) => r
            .iter()
            .map(|&v| match *n {
                "andrews" => {
                    let a = v.abs().max(f64::EPSILON.sqrt());
                    if a < std::f64::consts::PI { a.sin() / a } else { 0.0 }
                }
                "bisquare" => if v.abs() < 1.0 { (1.0 - v * v).powi(2) } else { 0.0 },
                "cauchy" => 1.0 / (1.0 + v * v),
                "fair" => 1.0 / (1.0 + v.abs()),
                "huber" => 1.0 / v.abs().max(1.0),
                "logistic" => {
                    let a = v.abs().max(f64::EPSILON.sqrt());
                    a.tanh() / a
                }
                "ols" => 1.0,
                "talwar" => if v.abs() < 1.0 { 1.0 } else { 0.0 },
                _ => (-(v * v)).exp(),
            })
            .collect(),
        WFun::User(Value::Func(fv)) => {
            let out = it.call_value(fv, vec![col(r.to_vec(), false)], 1)?;
            let Some(Value::Num(m)) = out.first() else { return rt("The weight function must return a numeric vector.") };
            m.data.clone()
        }
        WFun::User(other) => {
            let name = other.text()?;
            let out = it.call_function(&name, vec![col(r.to_vec(), false)], 1)?;
            let Some(Value::Num(m)) = out.first() else { return rt("The weight function must return a numeric vector.") };
            m.data.clone()
        }
    })
}

fn median_sorted(v: &[f64]) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 }
}

fn mad_sigma(r: &[f64], p: usize) -> f64 {
    let mut rs: Vec<f64> = r.iter().map(|v| v.abs()).collect();
    rs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let start = p.max(1) - 1;
    median_sorted(&rs[start.min(rs.len())..]) / 0.6745
}

fn std_of(v: &[f64]) -> f64 {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
}

fn weighted_fit(y: &[f64], x: &[Vec<f64>], w: &[f64]) -> (Vec<f64>, usize) {
    let n = y.len();
    let p = x[0].len();
    let sw: Vec<f64> = w.iter().map(|v| v.sqrt()).collect();
    let xw: Vec<Vec<f64>> = (0..n).map(|i| x[i].iter().map(|v| v * sw[i]).collect()).collect();
    let (q, r, perm) = qr_pivot(&xw, n, p);
    let tol = r[0][0].abs() * (n.max(p) as f64) * f64::EPSILON;
    let rank = (0..p).filter(|&j| r[j][j].abs() > tol).count();
    let rq: Vec<Vec<f64>> = r.iter().take(rank).map(|row| row[..rank].to_vec()).collect();
    let qty: Vec<f64> = (0..rank).map(|c| (0..n).map(|i| q[i][c] * y[i] * sw[i]).sum()).collect();
    let bp = solve_upper(&rq, &qty, rank);
    let mut b = vec![0.0; p];
    for (j, &pc) in perm.iter().take(rank).enumerate() {
        b[pc] = bp[j];
    }
    (b, rank)
}

fn robustfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("ROBUSTFIT requires at least two input arguments.");
    }
    let (Value::Num(xm), Value::Num(ym)) = (&a[0], &a[1]) else { return rt("X and Y must be numeric.") };
    let single = is_single(xm) || is_single(ym);
    // weight function and tuning constant
    let wval = match a.get(2) {
        Some(v) if !matches!(v, Value::Num(m) if m.is_empty() && m.kind != Kind::Char) => v.clone(),
        _ => Value::chars("bisquare"),
    };
    let mut t_default = None;
    let wfun = if wval.is_text() {
        let n = wval.text()?;
        match weight_default(&n) {
            Some((b, t)) => {
                t_default = Some(t);
                WFun::Builtin(b)
            }
            None => WFun::User(wval.clone()),
        }
    } else {
        WFun::User(wval.clone())
    };
    let tune = match a.get(3) {
        Some(Value::Num(m)) if !m.is_empty() => {
            if m.data[0] <= 0.0 {
                return rt("Tuning constant must be positive.");
            }
            m.data[0]
        }
        _ => t_default.unwrap_or(1.0),
    };
    let doconst = match a.get(4) {
        None => true,
        Some(v) if v.is_text() => match v.text()?.as_str() {
            "on" => true,
            "off" => false,
            _ => return rt("CONST must be 'on' or 'off'."),
        },
        Some(Value::Num(m)) if m.numel() == 1 && (m.data[0] == 1.0 || m.data[0] == 0.0) => m.data[0] == 1.0,
        _ => return rt("CONST must be 'on' or 'off'."),
    };
    let n0 = ym.numel();
    let xrows0 = if xm.dims[0] == 1 && xm.dims[1] > 1 && xm.numel() == n0 { (0..n0).map(|i| vec![xm.data[i]]).collect::<Vec<_>>() } else { rows_of(xm) };
    if xrows0.len() != n0 {
        return rt("Lengths of X and Y must match.");
    }
    let priorw0: Vec<f64> = match a.get(5) {
        Some(Value::Num(m)) if !m.is_empty() => m.data.clone(),
        _ => vec![1.0; n0],
    };
    if priorw0.len() != n0 {
        return rt("Lengths of X and Y must match.");
    }
    // the weight function must be usable
    if let WFun::User(v) = &wfun {
        let ok = match v {
            Value::Func(_) => true,
            t if t.is_text() => {
                let nm = t.text()?;
                it.call_function("exist", vec![Value::chars(&nm)], 1)?[0].scalar().map(|e| e != 0.0).unwrap_or(false)
            }
            _ => false,
        };
        if !ok {
            return rt("Weight function is not valid.");
        }
    }
    let keep: Vec<usize> = (0..n0).filter(|&i| !ym.data[i].is_nan() && !priorw0[i].is_nan() && xrows0[i].iter().all(|v| !v.is_nan())).collect();
    let n = keep.len();
    let mut xr: Vec<Vec<f64>> = keep.iter().map(|&i| xrows0[i].clone()).collect();
    let mut y: Vec<f64> = keep.iter().map(|&i| ym.data[i]).collect();
    let pw: Vec<f64> = keep.iter().map(|&i| priorw0[i]).collect();
    if doconst {
        for r in xr.iter_mut() {
            r.insert(0, 1.0);
        }
    }
    let p = xr.first().map_or(0, |r| r.len());
    if n <= p {
        return rt("There must be more observations than predictor variables (including the constant).");
    }
    let sw: Vec<f64> = if pw.iter().all(|&v| v == 1.0) { vec![1.0; n] } else { pw.iter().map(|v| v.sqrt()).collect() };
    if !pw.iter().all(|&v| v == 1.0) {
        for i in 0..n {
            for v in xr[i].iter_mut() {
                *v *= sw[i];
            }
            y[i] *= sw[i];
        }
    }
    let (q, r, perm) = qr_pivot(&xr, n, p);
    let tol = r[0][0].abs() * (n.max(p) as f64) * f64::EPSILON;
    let xrank = (0..p).filter(|&j| r[j][j].abs() > tol).count();
    let rq: Vec<Vec<f64>> = r.iter().take(xrank).map(|row| row[..xrank].to_vec()).collect();
    let qty: Vec<f64> = (0..xrank).map(|c| (0..n).map(|i| q[i][c] * y[i]).sum()).collect();
    let bp = solve_upper(&rq, &qty, xrank);
    let mut b = vec![0.0; p];
    for (j, &pc) in perm.iter().take(xrank).enumerate() {
        b[pc] = bp[j];
    }
    let perm_r: Vec<usize> = perm.iter().take(xrank).copied().collect();
    // leverage
    let mut h = vec![0.0; n];
    for i in 0..n {
        // E = X(:,perm)/R: solve R' e = x'
        let mut e = vec![0.0; xrank];
        for j in 0..xrank {
            let mut s = xr[i][perm_r[j]];
            for k in 0..j {
                s -= e[k] * rq[k][j];
            }
            e[j] = s / rq[j][j];
        }
        h[i] = e.iter().map(|v| v * v).sum::<f64>().min(0.9999);
    }
    let adj: Vec<f64> = h.iter().map(|v| 1.0 / (1.0 - v).sqrt()).collect();
    let dfe = n - xrank;
    let fitted = |b: &[f64]| -> Vec<f64> { (0..n).map(|i| (0..p).map(|j| xr[i][j] * b[j]).sum()).collect() };
    let ols_s = {
        let f = fitted(&b);
        ((0..n).map(|i| ((y[i] - f[i]) / sw[i]).powi(2)).sum::<f64>()).sqrt() / (dfe as f64).sqrt()
    };
    let mut tiny_s = 1e-6 * std_of(&y);
    if tiny_s == 0.0 {
        tiny_s = 1.0;
    }
    let dtol = f64::EPSILON.sqrt();
    let mut iter = 0;
    let mut wxrank = xrank;
    let mut b0 = vec![0.0; p];
    let mut w = vec![1.0; n];
    while iter == 0 || (0..p).any(|j| (b[j] - b0[j]).abs() > dtol * b[j].abs().max(b0[j].abs())) {
        iter += 1;
        if iter > 50 {
            break;
        }
        let f = fitted(&b);
        let radj: Vec<f64> = (0..n).map(|i| (y[i] - f[i]) * adj[i] / sw[i]).collect();
        let s = mad_sigma(&radj, wxrank);
        let scaled: Vec<f64> = radj.iter().map(|v| v / (s.max(tiny_s) * tune)).collect();
        w = apply_w(it, &wfun, &scaled)?;
        b0 = b.clone();
        let xp: Vec<Vec<f64>> = (0..n).map(|i| perm_r.iter().map(|&c| xr[i][c]).collect()).collect();
        let (bn, rk) = weighted_fit(&y, &xp, &w);
        for (j, &pc) in perm_r.iter().enumerate() {
            b[pc] = bn[j];
        }
        wxrank = rk;
    }
    let bv = col(b.clone(), single);
    if nargout < 2 {
        return Ok(Some(vec![bv]));
    }
    let f = fitted(&b);
    let resid: Vec<f64> = (0..n).map(|i| y[i] - f[i]).collect();
    let radj: Vec<f64> = (0..n).map(|i| resid[i] * adj[i] / sw[i]).collect();
    let mad_s = mad_sigma(&radj, xrank);
    let robust_s = if w.iter().all(|&v| v < dtol || v > 1.0 - dtol) {
        let inc: Vec<usize> = (0..n).filter(|&i| w[i] > 1.0 - dtol).collect();
        (inc.iter().map(|&i| resid[i] * resid[i]).sum::<f64>()).sqrt() / ((inc.len() as f64) - xrank as f64).sqrt()
    } else {
        let st = mad_s.max(tiny_s) * tune;
        let u: Vec<f64> = radj.iter().map(|v| v / st).collect();
        let wu = apply_w(it, &wfun, &u)?;
        let phi: Vec<f64> = u.iter().zip(&wu).map(|(a, b)| a * b).collect();
        let delta = 0.0001;
        let u0: Vec<f64> = u.iter().map(|v| v - delta).collect();
        let u1: Vec<f64> = u.iter().map(|v| v + delta).collect();
        let (w0, w1) = (apply_w(it, &wfun, &u0)?, apply_w(it, &wfun, &u1)?);
        let dphi: Vec<f64> = (0..n).map(|i| (u1[i] * w1[i] - u0[i] * w0[i]) / (2.0 * delta)).collect();
        let m1 = dphi.iter().sum::<f64>() / n as f64;
        let m2 = (0..n).map(|i| (1.0 - h[i]) * phi[i] * phi[i]).sum::<f64>() / (n as f64 - xrank as f64);
        let kk = 1.0 + (xrank as f64 / n as f64) * (1.0 - m1) / m1;
        kk * m2.sqrt() * st / m1
    };
    let nn = n as f64;
    let xr2 = (xrank * xrank) as f64;
    let sigma = robust_s.max(((ols_s * ols_s * xr2 + robust_s * robust_s * nn) / (xr2 + nn)).sqrt());
    // (R'R)^-1 sigma^2
    let mut rinv = vec![vec![0.0; xrank]; xrank];
    for j in 0..xrank {
        let mut e = vec![0.0; xrank];
        e[j] = 1.0;
        let c = solve_upper(&rq, &e, xrank);
        for i in 0..xrank {
            rinv[i][j] = c[i];
        }
    }
    let tc: Vec<Vec<f64>> = (0..xrank).map(|i| (0..xrank).map(|j| (0..xrank).map(|m| rinv[i][m] * rinv[j][m]).sum::<f64>() * sigma * sigma).collect()).collect();
    let tse: Vec<f64> = (0..xrank).map(|i| tc[i][i].max(f64::EPSILON).sqrt()).collect();
    let mut covb = vec![0.0; p * p];
    let mut cc = vec![f64::NAN; p * p];
    let mut se = vec![0.0; p];
    for i in 0..xrank {
        se[perm_r[i]] = tse[i];
        for j in 0..xrank {
            covb[perm_r[j] * p + perm_r[i]] = tc[i][j];
            cc[perm_r[j] * p + perm_r[i]] = tc[i][j] / (tse[i] * tse[j]);
        }
    }
    let t: Vec<f64> = (0..p).map(|i| if se[i] > 0.0 { b[i] / se[i] } else { f64::NAN }).collect();
    let negabs = Matrix::from_dims(vec![p, 1], t.iter().map(|v| -v.abs()).collect());
    let pv = it.call_function("tcdf", vec![Value::Num(negabs), Value::num(dfe as f64)], 1)?;
    let Value::Num(pm) = &pv[0] else { return rt("tcdf failed") };
    let pvals: Vec<f64> = pm.data.iter().map(|v| 2.0 * v).collect();
    let mut rr = vec![0.0; p * p];
    for i in 0..xrank {
        for j in 0..xrank {
            rr[perm_r[j] * p + perm_r[i]] = rq[i][j];
        }
    }
    let mut qy = vec![0.0; p];
    for i in 0..xrank {
        qy[perm_r[i]] = qty[i];
    }
    let reinsert = |v: &[f64]| -> Vec<f64> {
        let mut out = vec![f64::NAN; n0];
        for (r, &i) in keep.iter().enumerate() {
            out[i] = v[r];
        }
        out
    };
    let rstud: Vec<f64> = (0..n).map(|i| resid[i] * adj[i] / sigma).collect();
    let mut st = StructArray::scalar();
    st.set("ols_s", Value::num(ols_s));
    st.set("robust_s", Value::num(robust_s));
    st.set("mad_s", Value::num(mad_s));
    st.set("s", Value::num(sigma));
    st.set("resid", col(reinsert(&resid), single));
    let rst_full = reinsert(&rstud);
    st.set("rstud", col(rst_full, single));
    st.set("se", col(se, single));
    st.set("covb", mat(p, p, covb, single));
    st.set("coeffcorr", mat(p, p, cc, single));
    st.set("t", col(t, single));
    st.set("p", col(pvals, single));
    st.set("w", col(reinsert(&w), single));
    st.set("Qy", col(qy, single));
    st.set("R", mat(p, p, rr, single));
    st.set("dfe", Value::num(dfe as f64));
    st.set("h", col(reinsert(&h), single));
    st.set("Rtol", Value::num(tol));
    Ok(Some(vec![bv, Value::Struct(st)]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "classify" => classify(it, a, nargout),
        "perfcurve" => perfcurve(it, a, nargout),
        _ => robustfit(it, a, nargout),
    }
}
