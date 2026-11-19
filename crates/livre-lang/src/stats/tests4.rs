//! `capability` and `hist3`.

use super::dist::normcdf_std;
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["capability", "hist3"];

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
        _ => hist3(it, a, nargout),
    }
}
