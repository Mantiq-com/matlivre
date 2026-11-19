//! `capability` and `hist3`.

use super::dist::normcdf_std;
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["capability", "hist3", "controlrules", "fracfact", "ccdesign"];

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn row(v: Vec<f64>) -> Value {
    let n = v.len();
    Value::Num(Matrix::from_dims(vec![1, n], v))
}

fn capability(a: &[Value]) -> Ret {
    if a.len() > 2 {
        return rt("Too many input arguments.");
    }
    if a.len() != 2 {
        return rt("CAPABILITY must have exactly two input arguments, DATA and SPECS.");
    }
    let data = match &a[0] {
        Value::Num(m) if m.dims.len() == 2 && m.kind != Kind::Char => m,
        _ => return rt("DATA must be a vector or matrix."),
    };
    let (n, m) = if data.dims[0] == 1 { (data.dims[1], 1) } else { (data.dims[0], data.dims[1]) };
    let col = |j: usize| -> Vec<f64> { data.data[j * n..(j + 1) * n].to_vec() };
    let specs_in = match &a[1] {
        Value::Num(s) if s.kind != Kind::Char => Some(s),
        _ => None,
    };
    let mut specs: Vec<[f64; 2]> = Vec::new();
    if m == 1 {
        match specs_in {
            Some(s) if s.data.len() == 2 => specs.push([s.data[0], s.data[1]]),
            _ => return rt("SPECS must be a vector of lower and upper specification limits."),
        }
    } else {
        match specs_in {
            Some(s) if s.data.len() == 2 => specs = vec![[s.data[0], s.data[1]]; m],
            Some(s) if s.dims == [2, m] => specs = (0..m).map(|j| [s.data[2 * j], s.data[2 * j + 1]]).collect(),
            _ => return rt(format!("SPECS must be a 2-element vector or 2-by-{m} matrix of specification limits.")),
        }
    }
    let mut lb = vec![0.0; m];
    let mut ub = vec![0.0; m];
    for j in 0..m {
        let [mut l, mut u] = specs[j];
        if u < l {
            std::mem::swap(&mut l, &mut u);
        }
        lb[j] = if l.is_nan() { f64::NEG_INFINITY } else { l };
        ub[j] = if u.is_nan() { f64::INFINITY } else { u };
    }
    if (0..m).any(|j| lb[j].is_infinite() && ub[j].is_infinite()) {
        return rt("The SPECS argument cannot have both lower and upper bounds infinite.");
    }
    if (0..m).any(|j| lb[j] == ub[j]) {
        return rt("The SPECS argument cannot have equal lower and upper bounds.");
    }
    let (mut mu, mut sigma, mut pp, mut pl, mut pu) = (vec![], vec![], vec![], vec![], vec![]);
    let (mut cp, mut cpl, mut cpu, mut cpk) = (vec![], vec![], vec![], vec![]);
    for j in 0..m {
        let v: Vec<f64> = col(j).into_iter().filter(|x| !x.is_nan()).collect();
        let cnt = v.len() as f64;
        let mean = if v.is_empty() { f64::NAN } else { v.iter().sum::<f64>() / cnt };
        let sd = if v.len() > 1 { (v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (cnt - 1.0)).sqrt() } else if v.len() == 1 { 0.0 } else { f64::NAN };
        let (mean, sd) = if data.kind == Kind::Single { (mean as f32 as f64, sd as f32 as f64) } else { (mean, sd) };
        let cdf = |x: f64, mu: f64| normcdf_std((x - mu) / sd);
        let l = cdf(lb[j], mean);
        let u = cdf(-ub[j], -mean);
        mu.push(mean);
        sigma.push(sd);
        pl.push(l);
        pu.push(u);
        pp.push(1.0 - (l + u));
        let fin = |x: f64| x.is_finite();
        let c_p = if fin(ub[j]) && fin(lb[j]) { (ub[j] - lb[j]) / (6.0 * sd) } else { f64::NAN };
        let c_l = if fin(lb[j]) { (mean - lb[j]) / (3.0 * sd) } else { f64::NAN };
        let c_u = if fin(ub[j]) { (ub[j] - mean) / (3.0 * sd) } else { f64::NAN };
        cp.push(c_p);
        cpl.push(c_l);
        cpu.push(c_u);
        cpk.push(if c_l.is_nan() { c_u } else if c_u.is_nan() { c_l } else { c_l.min(c_u) });
    }
    let single = data.kind == Kind::Single;
    let rw = |v: Vec<f64>| {
        let n = v.len();
        Value::Num(finish(Matrix::from_dims(vec![1, n], v), single))
    };
    let mut st = StructArray::scalar();
    st.set("mu", rw(mu));
    st.set("sigma", rw(sigma));
    st.set("P", rw(pp));
    st.set("Pl", rw(pl));
    st.set("Pu", rw(pu));
    st.set("Cp", row(cp));
    st.set("Cpl", row(cpl));
    st.set("Cpu", row(cpu));
    st.set("Cpk", row(cpk));
    Ok(Some(vec![Value::Struct(st)]))
}

/// The bin of `v` among `edges` as `histc` numbers it (0 when outside; the last edge itself is bin n).
fn histc_bin(edges: &[f64], v: f64) -> usize {
    let n = edges.len();
    if v.is_nan() || n == 0 || v < edges[0] || v > edges[n - 1] {
        return 0;
    }
    if v == edges[n - 1] {
        return n;
    }
    edges.iter().filter(|&&e| e <= v).count()
}

fn cell_pair(v: &Value) -> Option<[Vec<f64>; 2]> {
    let Value::Cell(c) = v else { return None };
    if c.data.len() != 2 {
        return None;
    }
    let get = |x: &Value| match x {
        Value::Num(m) if m.kind != Kind::Char => Some(m.data.clone()),
        _ => None,
    };
    Some([get(&c.data[0])?, get(&c.data[1])?])
}

fn hist3(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Requires X.");
    }
    let x = a[0].clone();
    let (spec, rest): (Option<Value>, &[Value]) = if a.len() > 1 && !a[1].is_text() { (Some(a[1].clone()), &a[2..]) } else { (None, &a[1..]) };
    let mut nbins_v: Option<Value> = None;
    let mut ctrs_v: Option<Value> = None;
    let mut edges_v: Option<Value> = None;
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        let names = ["nbins", "ctrs", "edges"];
        if let Some(i) = names.iter().position(|n| !key.is_empty() && n.starts_with(&key)) {
            let v = if is_empty_val(&pair[1]) { None } else { Some(pair[1].clone()) };
            match i {
                0 => nbins_v = v,
                1 => ctrs_v = v,
                _ => edges_v = v,
            }
        }
    }
    let nspec = [nbins_v.is_none(), ctrs_v.is_none(), edges_v.is_none(), spec.as_ref().map(is_empty_val).unwrap_or(true)];
    if nspec.iter().filter(|b| **b).count() < 3 {
        return rt("Ambiguous specification of bins.");
    }
    if let Some(s) = spec.filter(|s| !is_empty_val(s)) {
        if matches!(s, Value::Cell(_)) {
            ctrs_v = Some(s);
        } else {
            nbins_v = Some(s);
        }
    }
    #[derive(PartialEq)]
    enum Mode {
        Auto,
        Ctrs,
        Edges,
    }
    let mut nbins = [10usize, 10usize];
    let mut nbins_nonint = false;
    let mut mode = Mode::Auto;
    let mut hist_behavior = true;
    let mut ctrs_in: [Vec<f64>; 2] = [vec![], vec![]];
    let mut edges_in: [Vec<f64>; 2] = [vec![], vec![]];
    if let Some(nv) = &nbins_v {
        match nv {
            Value::Num(m) if m.data.len() == 2 && m.kind != Kind::Char => {
                nbins_nonint = m.data.iter().any(|v| v.fract() != 0.0);
                nbins = [m.data[0].max(0.0) as usize, m.data[1].max(0.0) as usize]
            }
            _ => return rt("The number of bins must be specified with a 2-element numeric vector."),
        }
    } else if let Some(cv) = &ctrs_v {
        let Some(c) = cell_pair(cv) else { return rt("Bin centers must be specified with a cell array containing two numeric vectors.") };
        nbins = [c[0].len(), c[1].len()];
        ctrs_in = c;
        mode = Mode::Ctrs;
    } else if let Some(ev) = &edges_v {
        let Some(e) = cell_pair(ev) else { return rt("Bin edges must be specified with a cell array containing two numeric vectors.") };
        nbins = [e[0].len(), e[1].len()];
        edges_in = e;
        mode = Mode::Edges;
        hist_behavior = false;
    }
    let Value::Num(xm) = &x else { return rt("Invalid data type. First argument must be numeric or logical.") };
    if xm.kind == Kind::Char {
        return rt("Invalid data type. First argument must be numeric or logical.");
    }
    if mode == Mode::Auto && nbins_nonint {
        return rt("Size inputs must be integers.");
    }
    if xm.dims.len() != 2 || xm.dims[1] != 2 {
        return rt("X must be a matrix with two columns.");
    }
    let nrows = xm.dims[0];
    let mut ctrs_out: [Vec<f64>; 2] = [vec![], vec![]];
    let mut counts = vec![0.0; nbins[0] * nbins[1]];
    if nrows == 0 {
        if mode == Mode::Auto {
            ctrs_out = [(1..=nbins[0]).map(|v| v as f64).collect(), (1..=nbins[1]).map(|v| v as f64).collect()];
        } else if mode == Mode::Ctrs {
            ctrs_out = ctrs_in.clone();
        }
    } else {
        let mut bins = vec![[0usize; 2]; nrows];
        for i in 0..2 {
            let xs = &xm.data[i * nrows..(i + 1) * nrows];
            let (mn, mx) = xs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
            let hist_edges: Vec<f64>;
            match mode {
                Mode::Auto => {
                    let (mut lo, mut hi) = (mn, mx);
                    if lo.is_infinite() || hi.is_infinite() {
                        return rt("Bin centers or edges must be specified when data contain infinite values.");
                    }
                    let nb = nbins[i] as f64;
                    if lo == hi {
                        lo = lo - (nb / 2.0).floor() - 0.5;
                        hi = hi + (nb / 2.0).ceil() - 0.5;
                    }
                    let bw = (hi - lo) / nb;
                    let edges: Vec<f64> = (0..=nbins[i]).map(|k| lo + bw * k as f64).collect();
                    ctrs_out[i] = (0..nbins[i]).map(|k| edges[k] + bw / 2.0).collect();
                    let mut h = vec![f64::NEG_INFINITY];
                    h.extend_from_slice(&edges[1..edges.len() - 1]);
                    h.push(f64::INFINITY);
                    hist_edges = h;
                }
                Mode::Ctrs => {
                    let c = &ctrs_in[i];
                    ctrs_out[i] = c.clone();
                    let dc: Vec<f64> = c.windows(2).map(|w| w[1] - w[0]).collect();
                    if dc.is_empty() {
                        return rt("Index exceeds the number of array elements. Index must not exceed 0.");
                    }
                    let dc0 = *dc.first().unwrap_or(&f64::NAN);
                    let dcl = *dc.last().unwrap_or(&f64::NAN);
                    let mut left = vec![c[0]];
                    left.extend_from_slice(c);
                    let mut offs = vec![-dc0];
                    offs.extend_from_slice(&dc);
                    offs.push(dcl);
                    let edges: Vec<f64> = left.iter().zip(&offs).map(|(e, d)| e + d / 2.0).collect();
                    let mut h = vec![f64::NEG_INFINITY];
                    if edges.len() >= 2 {
                        h.extend_from_slice(&edges[1..edges.len() - 1]);
                    }
                    h.push(f64::INFINITY);
                    if h.windows(2).any(|w| w[1] < w[0]) {
                        return rt("Edge vector must be monotonically non-decreasing.");
                    }
                    hist_edges = h;
                }
                Mode::Edges => {
                    let e = &edges_in[i];
                    if e.windows(2).any(|w| w[1] < w[0]) {
                        return rt("Edge vector must be monotonically non-decreasing.");
                    }
                    if nargout > 1 && e.len() >= 2 {
                        let de: Vec<f64> = e.windows(2).map(|w| w[1] - w[0]).collect();
                        let mut c = vec![0.0; de.len()];
                        c[0] = e[0] + de[0] / 2.0;
                        for j in 1..c.len() {
                            c[j] = 2.0 * e[j] - c[j - 1];
                        }
                        let dl = *de.last().unwrap();
                        let bad = (0..c.len()).any(|j| c[j] <= e[j]) || (c[c.len() - 1] - (e[e.len() - 1] - dl / 2.0)).abs() > 1000.0 * eps_of(dl);
                        if bad {
                            it.emit_warning("Cannot compute centers that are consistent with EDGES.".into(), "stats:hist3:InconsistentEdges".into());
                            c = (0..de.len()).map(|j| e[j] + de[j] / 2.0).collect();
                        }
                        c.push(e[e.len() - 1] + dl / 2.0);
                        ctrs_out[i] = c;
                    }
                    hist_edges = e.clone();
                }
            }
            for r in 0..nrows {
                bins[r][i] = histc_bin(&hist_edges, xs[r]).min(nbins[i]);
            }
        }
        for b in &bins {
            if b[0] > 0 && b[1] > 0 {
                counts[(b[1] - 1) * nbins[0] + (b[0] - 1)] += 1.0;
            }
        }
    }
    let _ = hist_behavior;
    if nargout == 0 {
        return Ok(Some(vec![]));
    }
    let mut out = vec![Value::Num(Matrix::from_dims(vec![nbins[0], nbins[1]], counts))];
    if nargout > 1 {
        let single = xm.kind == Kind::Single;
        let rw = |v: Vec<f64>| {
            let n = v.len();
            Value::Num(finish(Matrix::from_dims(vec![1, n], v), single))
        };
        out.push(Value::Cell(CellArray::from_dims(vec![1, 2], vec![rw(ctrs_out[0].clone()), rw(ctrs_out[1].clone())])));
    }
    Ok(Some(out))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "capability" => capability(a),
        "controlrules" => controlrules(a),
        "fracfact" => fracfact(a, nargout),
        "ccdesign" => ccdesign(it, a, nargout),
        _ => hist3(it, a, nargout),
    }
}

// ---------------------------------------------------------------------------------------------
// controlrules

/// Counts, for each position, how many of the last `w` flags (including the current one) are set.
fn window_count(flags: &[bool], w: usize) -> Vec<usize> {
    (0..flags.len()).map(|i| flags[i.saturating_sub(w - 1)..=i].iter().filter(|&&f| f).count()).collect()
}

fn rule_names() -> Vec<String> {
    let mut v: Vec<String> = vec!["we".into(), "n".into()];
    v.extend((1..=10).map(|k| format!("we{k}")));
    v.extend((1..=8).map(|k| format!("n{k}")));
    v
}

fn controlrules(a: &[Value]) -> Ret {
    if a.len() < 4 {
        return rt("Not enough input arguments.");
    }
    if a.len() > 4 {
        return rt("Too many input arguments.");
    }
    let mut rules: Vec<String> = Vec::new();
    let given: Vec<Value> = match &a[0] {
        Value::Cell(c) => c.data.clone(),
        v if v.is_text() => vec![v.clone()],
        v if is_empty_val(v) => vec![],
        _ => return rt("Invalid RULES argument, must be a character array or cell array of character vectors"),
    };
    let all = rule_names();
    let names: Vec<&str> = all.iter().map(|s| s.as_str()).collect();
    for g in &given {
        if !g.is_text() {
            return rt("Invalid RULES argument, must be a character array or cell array of character vectors");
        }
        let t = g.text().unwrap_or_default();
        let low = t.to_lowercase();
        let idx = if let Some(i) = names.iter().position(|n| *n == low) {
            i
        } else {
            let hits: Vec<usize> = names.iter().enumerate().filter(|(_, n)| !low.is_empty() && n.starts_with(&low)).map(|(i, _)| i).collect();
            match hits.len() {
                1 => hits[0],
                0 => return rt(format!("'{t}' is not a valid value for the RULES argument.")),
                _ => return rt(format!("The value '{t}' is a partial match for more than one valid value for the RULES argument.")),
            }
        };
        let nm = names[idx];
        if nm == "we" {
            rules.extend((1..=10).map(|k| format!("we{k}")));
        } else if nm == "n" {
            rules.extend((1..=8).map(|k| format!("n{k}")));
        } else {
            rules.push(nm.to_string());
        }
    }
    let num = |v: &Value| -> Option<Vec<f64>> {
        match v {
            Value::Num(m) if m.kind != Kind::Char && m.kind != Kind::Logical => Some(m.data.clone()),
            Value::Num(m) if m.kind == Kind::Logical => Some(m.data.clone()),
            _ => None,
        }
    };
    let (Some(x0), Some(cl0), Some(se0)) = (num(&a[1]), num(&a[2]), num(&a[3])) else { return rt("X, CL, and SE must be numeric vectors.") };
    let n = x0.len();
    let cl: Vec<f64> = if cl0.len() == 1 {
        vec![cl0[0]; n]
    } else if cl0.len() == n {
        cl0
    } else {
        return rt("CL must be a scalar or must have the same length as X.");
    };
    let se: Vec<f64> = if se0.len() == 1 {
        vec![se0[0]; n]
    } else if se0.len() == n {
        se0
    } else {
        return rt("SE must be a scalar or must have the same length as X.");
    };
    let notnan: Vec<bool> = (0..n).map(|i| !(x0[i].is_nan() || cl[i].is_nan() || se[i].is_nan())).collect();
    let keep: Vec<usize> = (0..n).filter(|&i| notnan[i]).collect();
    let x: Vec<f64> = keep.iter().map(|&i| x0[i]).collect();
    let cl: Vec<f64> = keep.iter().map(|&i| cl[i]).collect();
    let se: Vec<f64> = keep.iter().map(|&i| se[i]).collect();
    let m = x.len();
    let idx: Vec<usize> = (0..m).collect();
    let gt = |off: f64, k: f64| -> Vec<bool> { idx.iter().map(|&i| x[i] > cl[i] + off * se[i] * k).collect() };
    let _ = gt;
    let above = |k: f64| -> Vec<bool> { idx.iter().map(|&i| x[i] > cl[i] + k * se[i]).collect() };
    let below = |k: f64| -> Vec<bool> { idx.iter().map(|&i| x[i] < cl[i] - k * se[i]).collect() };
    let mut cols: Vec<Vec<bool>> = Vec::new();
    for r in &rules {
        let c: Vec<bool> = match r.as_str() {
            "we1" => above(3.0),
            "we2" => {
                let hi = above(2.0);
                let cnt = window_count(&hi, 3);
                (0..m).map(|i| cnt[i] >= 2 && hi[i]).collect()
            }
            "we3" => {
                let hi = above(1.0);
                let cnt = window_count(&hi, 5);
                (0..m).map(|i| cnt[i] >= 4 && hi[i]).collect()
            }
            "we4" => {
                let f: Vec<bool> = idx.iter().map(|&i| x[i] > cl[i]).collect();
                window_count(&f, 8).iter().map(|&c| c >= 8).collect()
            }
            "we5" => below(3.0),
            "we6" => {
                let lo = below(2.0);
                let cnt = window_count(&lo, 3);
                (0..m).map(|i| cnt[i] >= 2 && lo[i]).collect()
            }
            "we7" => {
                let lo = below(1.0);
                let cnt = window_count(&lo, 5);
                (0..m).map(|i| cnt[i] >= 4 && lo[i]).collect()
            }
            "we8" => {
                let f: Vec<bool> = idx.iter().map(|&i| x[i] < cl[i]).collect();
                window_count(&f, 8).iter().map(|&c| c >= 8).collect()
            }
            "we9" | "n7" => {
                let f: Vec<bool> = idx.iter().map(|&i| x[i] > cl[i] - se[i] && x[i] < cl[i] + se[i]).collect();
                window_count(&f, 15).iter().map(|&c| c >= 15).collect()
            }
            "we10" | "n8" => {
                let f: Vec<bool> = idx.iter().map(|&i| x[i] < cl[i] - se[i] || x[i] > cl[i] + se[i]).collect();
                window_count(&f, 8).iter().map(|&c| c >= 8).collect()
            }
            "n1" => idx.iter().map(|&i| x[i] > cl[i] + 3.0 * se[i] || x[i] < cl[i] - 3.0 * se[i]).collect(),
            "n2" => {
                let up: Vec<bool> = idx.iter().map(|&i| x[i] > cl[i]).collect();
                let dn: Vec<bool> = idx.iter().map(|&i| x[i] < cl[i]).collect();
                let (c1, c2) = (window_count(&up, 9), window_count(&dn, 9));
                (0..m).map(|i| c1[i] >= 9 || c2[i] >= 9).collect()
            }
            "n3" => {
                let d: Vec<f64> = (0..m).map(|i| if i == 0 { 0.0 } else { x[i] - x[i - 1] }).collect();
                let up: Vec<bool> = d.iter().map(|&v| v > 0.0).collect();
                let dn: Vec<bool> = d.iter().map(|&v| v < 0.0).collect();
                let (c1, c2) = (window_count(&up, 5), window_count(&dn, 5));
                (0..m).map(|i| c1[i] >= 5 || c2[i] >= 5).collect()
            }
            "n4" => {
                let d: Vec<f64> = (0..m)
                    .map(|i| {
                        let sg = if i == 0 { 0.0 } else { (x[i] - x[i - 1]).signum() * if x[i] == x[i - 1] { 0.0 } else { 1.0 } };
                        sg * if i % 2 == 1 { -1.0 } else { 1.0 }
                    })
                    .collect();
                (0..m).map(|i| d[i.saturating_sub(12)..=i].iter().sum::<f64>()).map(|s| s <= -13.0 || s >= 13.0).collect()
            }
            "n5" => {
                let (hi, lo) = (above(2.0), below(2.0));
                let (c1, c2) = (window_count(&hi, 3), window_count(&lo, 3));
                (0..m).map(|i| (c1[i] >= 2 && hi[i]) || (c2[i] >= 2 && lo[i])).collect()
            }
            _ => {
                // n6
                let (hi, lo) = (above(1.0), below(1.0));
                let (c1, c2) = (window_count(&hi, 5), window_count(&lo, 5));
                (0..m).map(|i| (c1[i] >= 4 && hi[i]) || (c2[i] >= 4 && lo[i])).collect()
            }
        };
        cols.push(c);
    }
    let rows = n;
    let nr = rules.len();
    let mut out = vec![0.0; rows * nr];
    for (j, c) in cols.iter().enumerate() {
        for (k, &i) in keep.iter().enumerate() {
            out[j * rows + i] = if c[k] { 1.0 } else { 0.0 };
        }
    }
    let rm = Matrix::from_dims(vec![rows, nr], out).cast_to(Kind::Logical);
    let rcell = Value::Cell(CellArray::from_dims(vec![1, nr], rules.iter().map(|s| Value::chars(s)).collect()));
    Ok(Some(vec![Value::Num(rm), rcell]))
}

// ---------------------------------------------------------------------------------------------
// fracfact

/// All k-element subsets of 1..=n in lexicographic order.
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if k == 0 || k > n {
        return out;
    }
    let mut idx: Vec<usize> = (1..=k).collect();
    loop {
        out.push(idx.clone());
        let Some(i) = (0..k).rev().find(|&i| idx[i] < n - (k - 1 - i)) else { return out };
        idx[i] += 1;
        for j in i + 1..k {
            idx[j] = idx[j - 1] + 1;
        }
    }
}

fn fracfact(a: &[Value], nargout: usize) -> Ret {
    const BAD: &str = "Invalid generator character vector or string scalar.";
    let Some(g0) = a.first() else { return rt("Not enough input arguments.") };
    let gen_text: String = match g0 {
        Value::Cell(c) => {
            let mut s = String::new();
            for v in &c.data {
                if !v.is_text() {
                    return rt(BAD);
                }
                s.push_str(&v.text().unwrap_or_default());
                s.push(' ');
            }
            s
        }
        v if v.is_text() => {
            if let Value::Num(m) = v {
                if m.dims[0] > 1 {
                    return rt(BAD);
                }
            }
            v.text().unwrap_or_default()
        }
        _ => return rt(BAD),
    };
    if !gen_text.chars().all(|c| c.is_ascii_alphabetic() || c.is_whitespace() || c == '-') {
        return rt(BAD);
    }
    let mut maxint_v: Option<Value> = None;
    let mut names_v: Option<Value> = None;
    let rest = &a[1..];
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        if !key.is_empty() && "maxint".starts_with(&key) {
            maxint_v = Some(pair[1].clone());
        } else if !key.is_empty() && "factornames".starts_with(&key) {
            names_v = Some(pair[1].clone());
        } else {
            return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default()));
        }
    }
    let trimmed = gen_text.trim();
    if trimmed.is_empty() {
        return rt("First input can not be empty. Expected a non-empty character vector or a valid file-id.");
    }
    let mut words: Vec<String> = trimmed.split_whitespace().map(|s| s.to_string()).collect();
    if !trimmed.chars().any(|c| c.is_ascii_lowercase()) {
        for w in words.iter_mut() {
            *w = w.to_lowercase();
        }
    }
    let seq: Vec<char> = ('a'..='z').chain('A'..='Z').collect();
    let letters: Vec<char> = seq.into_iter().filter(|c| words.iter().any(|w| w.contains(*c))).collect();
    let k = letters.len();
    if k == 0 {
        return rt(BAD);
    }
    let mut alpha: Vec<char> = vec!['-'];
    alpha.extend(letters.iter().copied());
    let n = 1usize << k;
    let p = words.len();
    let full = |row: usize, col: usize| -> f64 {
        if col == 0 {
            -1.0
        } else {
            let bit = (row >> (k - col)) & 1;
            -1.0 + 2.0 * bit as f64
        }
    };
    let mut x = vec![0.0; n * p];
    let mut bitmat = vec![vec![false; k + 1]; p];
    for (j, w) in words.iter().enumerate() {
        let mut cols: Vec<usize> = w.chars().filter_map(|c| alpha.iter().position(|&a| a == c)).collect();
        cols.sort();
        for r in 0..n {
            x[j * n + r] = cols.iter().map(|&c| full(r, c)).product();
        }
        for &c in &cols {
            bitmat[j][c] = !bitmat[j][c];
        }
    }
    let mut out = vec![Value::Num(Matrix::from_dims(vec![n, p], x))];
    if nargout > 1 {
        let maxint = match &maxint_v {
            None => 2.0,
            Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] > 0.0 && m.data[0] == m.data[0].floor() => m.data[0],
            Some(_) => return rt("MaxInt parameter must be a positive integer."),
        } as usize;
        let names: Vec<String> = match &names_v {
            None => (1..=p).map(|i| format!("X{i}")).collect(),
            Some(Value::Cell(c)) if c.data.len() == p && c.data.iter().all(|v| v.is_text()) => c.data.iter().map(|v| v.text().unwrap_or_default()).collect(),
            Some(_) => return rt("FactorNames value must be a cell array with one name for each factor."),
        };
        let klen = alpha.len();
        let gen_of = |mask: &[bool]| -> String {
            let s: String = (0..klen).filter(|&i| mask[i]).map(|i| alpha[i]).collect();
            if s.is_empty() { "1".to_string() } else { s }
        };
        let mut rows: Vec<(String, String, Option<String>)> = vec![("Term".into(), "Generator".into(), Some("Confounding".into()))];
        for j in 0..p {
            rows.push((names[j].clone(), gen_of(&bitmat[j]), None));
        }
        for j in 2..=maxint.min(klen) {
            if p == 1 {
                return rt("K must be an integer between 0 and N.");
            }
            for g in combinations(p, j) {
                let name = g.iter().map(|&i| rows[i].0.clone()).collect::<Vec<_>>().join("*");
                let mask: Vec<bool> = (0..=k).map(|c| g.iter().map(|&i| bitmat[i - 1][c]).filter(|&b| b).count() % 2 == 1).collect();
                rows.push((name, gen_of(&mask), None));
            }
        }
        for j in 1..rows.len() {
            if rows[j].2.is_some() {
                continue;
            }
            let gen = rows[j].1.clone();
            let other: String = if gen.is_empty() {
                "- -".to_string()
            } else if gen.starts_with('-') {
                gen[1..].to_string()
            } else {
                format!("-{gen}")
            };
            let mut sel: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].1 == gen || rows[i].1 == other).collect();
            sel.sort();
            let rows1: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].1 == gen).collect();
            let mut t = rows[sel[0]].0.clone();
            if rows[sel[0]].1 != gen {
                t = format!("-{t}");
            }
            for &i in &sel[1..] {
                if rows[i].1 == gen {
                    t = format!("{t} + {}", rows[i].0);
                } else {
                    t = format!("{t} - {}", rows[i].0);
                }
            }
            for i in rows1 {
                rows[i].2 = Some(t.clone());
            }
        }
        let nr = rows.len();
        let mut cells: Vec<Value> = Vec::with_capacity(nr * 3);
        cells.extend(rows.iter().map(|r| Value::chars(&r.0)));
        cells.extend(rows.iter().map(|r| Value::chars(&r.1)));
        cells.extend(rows.iter().map(|r| Value::chars(r.2.as_deref().unwrap_or(""))));
        out.push(Value::Cell(CellArray::from_dims(vec![nr, 3], cells)));
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// ccdesign

/// The generator words for the cube portion (and the block generators) of a central composite design.
fn cc_getgen(nfactors: usize, frac: Option<usize>, blocked: bool) -> (String, usize, String) {
    let frac = frac.unwrap_or_else(|| {
        let best = [0usize, 0, 0, 0, 1, 1, 1, 2, 2, 3, 4];
        if nfactors <= best.len() { best[nfactors - 1] } else { 1 }
    });
    let mut g = String::new();
    let mut blk = String::new();
    let nletters = nfactors.saturating_sub(frac);
    let alf: Vec<char> = (0..nletters).map(|i| (b'a' + i as u8) as char).collect();
    match nfactors {
        3 => blk = "abc".into(),
        4 => blk = "abcd".into(),
        5 => {
            if frac == 0 {
                blk = "abc bcde".into();
            } else if frac == 1 {
                g = "abcd".into();
            }
        }
        6 => {
            if frac == 0 {
                blk = "abcd cdef".into();
            } else if frac == 1 {
                g = "abcde".into();
                blk = "abc".into();
            }
        }
        7 => {
            if frac == 1 {
                g = "abcdef".into();
                blk = "abcd abef bdf".into();
            }
        }
        8 => {
            if frac == 1 {
                g = "abcdefg".into();
            } else if frac == 2 {
                g = "abcd abef".into();
                blk = "ace bdf".into();
            }
        }
        9 => {
            if frac == 1 {
                g = "abcdefgh".into();
            } else if frac == 2 {
                g = "abcde cdefg".into();
            }
        }
        10 => {
            if frac == 1 {
                g = "abcdefghi".into();
            } else if frac == 2 {
                g = "abcde defgh".into();
            } else if frac == 3 {
                g = "abcd bceg acfg".into();
            }
        }
        11 => {
            if frac == 1 {
                g = "abcdefghij".into();
            } else if frac == 2 {
                g = "abcdef defghi".into();
            } else if frac == 3 {
                g = "abcde bcdgh acefg".into();
            } else if frac == 4 {
                g = "abcd bceg acfg abcdefg".into();
            }
        }
        _ => {}
    }
    if nfactors > 11 && frac == 1 {
        g = alf.iter().collect();
    }
    if frac == 0 || !g.is_empty() {
        let prefix: String = alf.iter().map(|c| format!("{c} ")).collect();
        g = format!("{prefix}{g}");
    }
    if !blocked {
        blk.clear();
    }
    (g, frac, blk)
}

fn cc_ncenter_uniform(k: usize, frac: usize) -> Option<f64> {
    if (2..=8).contains(&k) && frac == 0 {
        Some([f64::NAN, 5.0, 6.0, 7.0, 10.0, 15.0, 21.0, 28.0][k - 1])
    } else if (5..=8).contains(&k) && frac == 1 {
        Some([f64::NAN, f64::NAN, f64::NAN, f64::NAN, 6.0, 9.0, 14.0, 20.0][k - 1])
    } else if k == 8 && frac == 2 {
        Some(13.0)
    } else {
        None
    }
}

fn cc_ctrblks(ncenter: f64, blocksize: f64, nb: usize, cblocksize: f64, sblocksize: f64) -> Vec<f64> {
    let nbf = nb as f64;
    let nper = (ncenter / nbf).ceil();
    let maxncs = blocksize - sblocksize;
    let maxncc = blocksize - cblocksize;
    let (ncs, ncc);
    if sblocksize + nper > blocksize {
        ncs = blocksize - sblocksize;
        ncc = (maxncc.min((ncenter - ncs) / (nbf - 1.0))).ceil();
    } else if cblocksize + nper > blocksize {
        ncc = blocksize - cblocksize;
        ncs = maxncs.min(ncenter - ncc * (nbf - 1.0));
    } else {
        ncs = nper;
        ncc = nper;
    }
    let ncenter = ncenter.min(ncs + (nbf - 1.0) * ncc);
    let mut base: Vec<f64> = Vec::new();
    for _ in 0..(ncc.max(0.0) as usize) {
        base.extend((1..nb).map(|v| v as f64));
    }
    let take = (ncenter - ncs).max(0.0) as usize;
    let mut out: Vec<f64> = base.into_iter().take(take).collect();
    out.extend(std::iter::repeat(nbf).take(ncs.max(0.0) as usize));
    out
}

fn ccdesign(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let nf = match a.first() {
        Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] >= 2.0 && m.data[0] <= 26.0 && m.data[0] == m.data[0].round() => m.data[0] as usize,
        Some(_) => return rt("Number of factors must be an integer between 2 and 26."),
        None => return rt("Not enough input arguments."),
    };
    let rest = &a[1..];
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let mut center_v: Option<Value> = None;
    let mut frac_v: Option<Value> = None;
    let mut type_v: Option<Value> = None;
    let mut blocksize_v: Option<Value> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        let names = ["center", "fraction", "type", "blocksize"];
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, n)| !key.is_empty() && n.starts_with(&key)).map(|(i, _)| i).collect();
        if hits.len() != 1 {
            return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default()));
        }
        let v = Some(pair[1].clone());
        match hits[0] {
            0 => center_v = v,
            1 => frac_v = v,
            2 => type_v = v,
            _ => blocksize_v = v,
        }
    }
    let oktypes = ["circumscribed", "inscribed", "faced"];
    let startype = match &type_v {
        None => "circumscribed",
        Some(v) if v.is_text() && !v.text().unwrap_or_default().is_empty() => {
            let t = v.text().unwrap_or_default().to_lowercase();
            match oktypes.iter().find(|o| o.starts_with(&t)) {
                Some(o) => o,
                None => return rt("Valid types are 'inscribed', 'circumscribed', and 'faced'."),
            }
        }
        _ => return rt("Valid types are 'inscribed', 'circumscribed', and 'faced'."),
    };
    let blocksize = match &blocksize_v {
        None => f64::INFINITY,
        Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() == 1 && (m.data[0].is_infinite() && m.data[0] > 0.0 || (m.data[0] >= 1.0 && m.data[0] == m.data[0].floor())) => m.data[0],
        Some(_) => return rt("Value of 'blocksize' parameter must be a positive integer."),
    };
    let blocked = nargout >= 2 || blocksize.is_finite();
    let frac_opt: Option<usize> = match &frac_v {
        Some(Value::Num(m)) if !m.is_empty() && m.data.len() == 1 && m.data[0] >= 0.0 => Some(m.data[0] as usize),
        _ => None,
    };
    let (g, frac, bg) = cc_getgen(nf, frac_opt, blocked);
    if g.is_empty() {
        return rt(format!("Cannot find design for {nf} factors with a 1/2^{frac} fractional cube."));
    }
    let nc = 1usize << (nf - frac);
    let ns = 2 * nf;
    let alpha = if startype == "faced" { 1.0 } else { (nc as f64).powf(0.25) };
    let mut useorth = false;
    let mut ncenter_default = false;
    let mut ncenter: f64 = 0.0;
    match &center_v {
        None => ncenter_default = true,
        Some(v) if is_empty_val(v) => ncenter_default = true,
        Some(v) if v.is_text() => {
            let t = v.text().unwrap_or_default().to_lowercase();
            if !t.is_empty() && "uniform".starts_with(&t) {
                match cc_ncenter_uniform(nf, frac) {
                    Some(n0) => ncenter = n0,
                    None => {
                        useorth = true;
                        it.emit_warning("Number of center points for uniform precision unknown for\nthis design, using center points for orthogonality instead.".into(), "stats:ccdesign:NoUniformPrecision".into());
                    }
                }
            } else if !t.is_empty() && "orthogonal".starts_with(&t) {
                ncenter_default = true;
            } else {
                return rt("Bad value for 'center' parameter.");
            }
        }
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] >= 0.0 && m.data[0] == m.data[0].floor() => ncenter = m.data[0],
        Some(_) => return rt("Bad value for 'center' parameter."),
    }
    if ncenter_default || useorth {
        ncenter = (4.0 * (nc as f64).sqrt() + 4.0 - 2.0 * nf as f64).ceil().max(1.0);
    }
    let full = it.call_function("fracfact", vec![Value::chars(&format!("{g} {bg}"))], 1)?;
    let Value::Num(dcm) = &full[0] else { return rt("fracfact") };
    let ncols = dcm.dims[1];
    let mut dc: Vec<Vec<f64>> = (0..nc).map(|r| (0..nf).map(|c| dcm.data[c * nc + r]).collect()).collect();
    let mut blk: Vec<f64> = vec![];
    let mut nblocks = 0usize;
    let mut ctrblks: Vec<f64> = vec![];
    if blocked {
        blk = vec![1.0; nc];
        let numblkgen = if !bg.is_empty() { ncols - nf } else { 0 };
        let blkcols: Vec<Vec<f64>> = (0..nc).map(|r| (nf..ncols).map(|c| dcm.data[c * nc + r]).collect()).collect();
        let minblk = 1.0 + (ns as f64).max(nc as f64 / 2f64.powi(numblkgen as i32));
        if minblk > blocksize {
            it.emit_warning(format!("Cannot find a design with fewer than {} runs per block.", minblk), "stats:ccdesign:BlockSizeExceeded".into());
        }
        let mut cblocksize = nc as f64;
        if numblkgen > 0 {
            let b = (nc as f64 / (blocksize - 1.0)).max(1.0).log2().ceil();
            if b > 0.0 {
                let b = b as usize;
                let used: Vec<Vec<f64>> = blkcols.iter().map(|r| r[..b.min(r.len())].to_vec()).collect();
                let mut uniq = used.clone();
                uniq.sort_by(|x, y| x.partial_cmp(y).unwrap());
                uniq.dedup();
                blk = used.iter().map(|r| (uniq.iter().position(|u| u == r).unwrap() + 1) as f64).collect();
                cblocksize = blk.iter().filter(|&&v| v == 1.0).count() as f64;
            }
        }
        nblocks = blk.iter().cloned().fold(0.0, f64::max) as usize + 1;
        ctrblks = cc_ctrblks(ncenter, blocksize, nblocks, cblocksize, ns as f64);
        if ncenter > ctrblks.len() as f64 {
            it.emit_warning(format!("Block size limit allows only {} center points.", ctrblks.len()), "stats:ccdesign:CenterPointLimit".into());
            ncenter = ctrblks.len() as f64;
        } else if ncenter < nblocks as f64 {
            it.emit_warning(format!("Number of center points ({}) less than number of blocks ({}).", ncenter, nblocks), "stats:ccdesign:CenterPointDeficit".into());
        }
    }
    let ncen = ncenter as usize;
    let nruns = ns + nc + ncen;
    let mut d = vec![0.0; nruns * nf];
    for r in 0..nc {
        for c in 0..nf {
            d[c * nruns + r] = dc[r][c];
        }
    }
    for i in 0..nf {
        d[i * nruns + nc + 2 * i] = -alpha;
        d[i * nruns + nc + 2 * i + 1] = alpha;
    }
    dc.clear();
    if startype == "inscribed" {
        for v in d.iter_mut() {
            *v /= alpha;
        }
    }
    let mut out = vec![Value::Num(Matrix::from_dims(vec![nruns, nf], d))];
    if nargout > 1 {
        let mut bv = vec![0.0; nruns];
        bv[..nc].copy_from_slice(&blk[..nc]);
        for k in nc..nc + ns {
            bv[k] = nblocks as f64;
        }
        for k in 0..ncen {
            bv[nc + ns + k] = ctrblks[k];
        }
        out.push(Value::Num(Matrix::from_dims(vec![nruns, 1], bv)));
    }
    Ok(Some(out))
}
