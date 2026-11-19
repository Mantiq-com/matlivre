//! `knnsearch`, `rangesearch` (exhaustive search with the toolbox's distances), `crosstab`, `polyconf`.

use super::tests2::{fmt_num, labels_of, Labels};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["knnsearch", "rangesearch", "crosstab", "polyconf"];

#[derive(Clone)]
pub(crate) enum Metric {
    Euclidean,
    SEuclidean(Vec<f64>),
    CityBlock,
    Chebychev,
    Minkowski(f64),
    Mahalanobis(Vec<Vec<f64>>),
    Cosine,
    Correlation,
    Spearman,
    Hamming,
    Jaccard,
    Custom(FuncValue),
}

struct Opts {
    k: usize,
    metric: Metric,
    include_ties: bool,
    sort: bool,
}

pub(crate) fn cols_stats(rows: &[Vec<f64>], p: usize) -> (Vec<f64>, Vec<f64>) {
    let ok: Vec<&Vec<f64>> = rows.iter().filter(|r| r.iter().all(|v| !v.is_nan())).collect();
    let n = ok.len() as f64;
    let mean: Vec<f64> = (0..p).map(|j| ok.iter().map(|r| r[j]).sum::<f64>() / n).collect();
    let sd: Vec<f64> = (0..p).map(|j| (ok.iter().map(|r| (r[j] - mean[j]).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()).collect();
    (mean, sd)
}

pub(crate) fn rows_of(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1]);
    (0..r).map(|i| (0..c).map(|j| m.data[j * r + i]).collect()).collect()
}

fn ranks(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        for t in i..=j {
            r[idx[t]] = (i + j) as f64 / 2.0 + 1.0;
        }
        i = j + 1;
    }
    r
}

fn corr_dist(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let za: Vec<f64> = a.iter().map(|v| v - ma).collect();
    let zb: Vec<f64> = b.iter().map(|v| v - mb).collect();
    let dot: f64 = za.iter().zip(&zb).map(|(x, y)| x * y).sum();
    1.0 - dot / (za.iter().map(|x| x * x).sum::<f64>().sqrt() * zb.iter().map(|x| x * x).sum::<f64>().sqrt())
}

pub(crate) fn dist_to_all(it: &mut Interpreter, m: &Metric, y: &[f64], xs: &[Vec<f64>], xm: &Matrix) -> Res<Vec<f64>> {
    let p = y.len();
    Ok(match m {
        Metric::Custom(f) => {
            let ym = Value::Num(Matrix::from_dims(vec![1, p], y.to_vec()));
            let r = it.call_value(f, vec![ym, Value::Num(xm.clone())], 1)?;
            let Some(Value::Num(d)) = r.first() else { return rt("The distance function must return a numeric column vector.") };
            if d.numel() != xs.len() {
                return rt(format!("The distance function must return a column vector with {} rows.", xs.len()));
            }
            d.data.clone()
        }
        _ => xs
            .iter()
            .map(|x| match m {
                Metric::Euclidean => x.iter().zip(y).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt(),
                Metric::SEuclidean(s) => (0..p).map(|j| ((x[j] - y[j]) / s[j]).powi(2)).sum::<f64>().sqrt(),
                Metric::CityBlock => x.iter().zip(y).map(|(a, b)| (a - b).abs()).sum(),
                Metric::Chebychev => x.iter().zip(y).map(|(a, b)| (a - b).abs()).fold(0.0, |m, v| if v.is_nan() || v > m { v } else { m }),
                Metric::Minkowski(q) => x.iter().zip(y).map(|(a, b)| (a - b).abs().powf(*q)).sum::<f64>().powf(1.0 / q),
                Metric::Mahalanobis(ci) => {
                    let d: Vec<f64> = (0..p).map(|j| x[j] - y[j]).collect();
                    let mut s = 0.0;
                    for i in 0..p {
                        for j in 0..p {
                            s += d[i] * ci[i][j] * d[j];
                        }
                    }
                    s.sqrt()
                }
                Metric::Cosine => {
                    let dot: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
                    1.0 - dot / (x.iter().map(|a| a * a).sum::<f64>().sqrt() * y.iter().map(|a| a * a).sum::<f64>().sqrt())
                }
                Metric::Correlation => corr_dist(x, y),
                Metric::Spearman => corr_dist(&ranks(x), &ranks(y)),
                Metric::Hamming => x.iter().zip(y).filter(|(a, b)| a != b).count() as f64 / p as f64,
                Metric::Jaccard => {
                    let nz = x.iter().zip(y).filter(|(a, b)| **a != 0.0 || **b != 0.0).count() as f64;
                    if nz == 0.0 {
                        0.0
                    } else {
                        x.iter().zip(y).filter(|(a, b)| (**a != 0.0 || **b != 0.0) && a != b).count() as f64 / nz
                    }
                }
                Metric::Custom(_) => unreachable!(),
            })
            .collect(),
    })
}

pub(crate) fn make_metric(it: &mut Interpreter, dist: Option<Value>, pexp: f64, cov: Option<Matrix>, scale: Option<Matrix>, xs: &[Vec<f64>], p: usize) -> Res<Metric> {
    Ok(match dist {
        None => Metric::Euclidean,
        Some(Value::Func(f)) => Metric::Custom(f),
        Some(v) => {
            let name = v.text().unwrap_or_default();
            match name.to_lowercase().as_str() {
                "euclidean" => Metric::Euclidean,
                "seuclidean" => Metric::SEuclidean(match scale {
                    Some(s) => s.data.clone(),
                    None => cols_stats(xs, p).1,
                }),
                "cityblock" => Metric::CityBlock,
                "chebychev" => Metric::Chebychev,
                "minkowski" => Metric::Minkowski(pexp),
                "mahalanobis" => {
                    let c: Vec<f64> = match cov {
                        Some(c) => c.data.clone(),
                        None => {
                            let ok: Vec<&Vec<f64>> = xs.iter().filter(|r| r.iter().all(|v| !v.is_nan())).collect();
                            let n = ok.len() as f64;
                            let mean: Vec<f64> = (0..p).map(|j| ok.iter().map(|r| r[j]).sum::<f64>() / n).collect();
                            (0..p * p).map(|idx| { let (i, j) = (idx % p, idx / p); ok.iter().map(|r| (r[i] - mean[i]) * (r[j] - mean[j])).sum::<f64>() / (n - 1.0) }).collect()
                        }
                    };
                    let inv = it.call_function("inv", vec![Value::Num(Matrix::from_dims(vec![p, p], c))], 1)?;
                    let Value::Num(im) = &inv[0] else { return rt("inv failed") };
                    Metric::Mahalanobis((0..p).map(|i| (0..p).map(|j| im.data[j * p + i]).collect()).collect())
                }
                "cosine" => Metric::Cosine,
                "correlation" => Metric::Correlation,
                "spearman" => Metric::Spearman,
                "hamming" => Metric::Hamming,
                "jaccard" => Metric::Jaccard,
                _ => return rt(format!("'{name}' is not a valid value for the Distance argument.")),
            }
        }
    })
}

fn parse(it: &mut Interpreter, a: &[Value], who: &str, xs: &[Vec<f64>], p: usize) -> Res<Opts> {
    let mut o = Opts { k: 1, metric: Metric::Euclidean, include_ties: false, sort: true };
    let rest = if who == "knnsearch" { &a[2..] } else { &a[3..] };
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let mut dist: Option<Value> = None;
    let (mut pexp, mut cov, mut scale): (f64, Option<Matrix>, Option<Matrix>) = (2.0, None, None);
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        match key.to_lowercase().as_str() {
            "k" if who == "knnsearch" => {
                let v = pair[1].scalar().unwrap_or(f64::NAN);
                if !(v >= 1.0 && v.fract() == 0.0 && v.is_finite()) {
                    return rt("K must be a positive integer specifying the number of neighbors.");
                }
                o.k = v as usize;
            }
            "distance" => dist = Some(pair[1].clone()),
            "p" => pexp = pair[1].scalar()?,
            "cov" => cov = Some(pair[1].as_matrix()?.clone()),
            "scale" => scale = Some(pair[1].as_matrix()?.clone()),
            "includeties" if who == "knnsearch" => o.include_ties = pair[1].scalar()? != 0.0,
            "sortindices" => o.sort = pair[1].scalar()? != 0.0,
            "nsmethod" | "bucketsize" => {}
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    o.metric = make_metric(it, dist, pexp, cov, scale, xs, p)?;
    Ok(o)
}

fn xy_args<'a>(a: &'a [Value], who: &str) -> Res<(&'a Matrix, &'a Matrix)> {
    let need = if who == "knnsearch" { 2 } else { 3 };
    if a.len() < need {
        return rt("Not enough input arguments.");
    }
    let (Value::Num(x), Value::Num(y)) = (&a[0], &a[1]) else { return rt("X and Y must be numeric matrices.") };
    if x.dims.len() != 2 || y.dims.len() != 2 {
        return rt("X and Y must be matrices.");
    }
    if y.dims[1] != x.dims[1] {
        return rt(format!("Y must be a matrix with {} columns.", x.dims[1]));
    }
    Ok((x, y))
}

fn row_val(v: Vec<f64>, single: bool) -> Value {
    let n = v.len();
    Value::Num(finish(Matrix::from_dims(vec![1, n], v), single))
}

fn knnsearch(it: &mut Interpreter, a: &[Value]) -> Ret {
    let (x, y) = xy_args(a, "knnsearch")?;
    let single = is_single(x) || is_single(y);
    let (xs, ys) = (rows_of(x), rows_of(y));
    let o = parse(it, a, "knnsearch", &xs, x.dims[1])?;
    let n = xs.len();
    let k = o.k.min(n);
    let nq = ys.len();
    let mut idx_rows: Vec<Vec<f64>> = Vec::with_capacity(nq);
    let mut d_rows: Vec<Vec<f64>> = Vec::with_capacity(nq);
    for q in &ys {
        let d = dist_to_all(it, &o.metric, q, &xs, x)?;
        let mut cand: Vec<(f64, usize)> = d.iter().copied().enumerate().filter(|(_, v)| !v.is_nan()).map(|(i, v)| (v, i)).collect();
        cand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        if cand.is_empty() || q.iter().any(|v| v.is_nan()) {
            idx_rows.push((1..=k).map(|v| v as f64).collect());
            d_rows.push(vec![f64::NAN; k]);
            continue;
        }
        let mut take = k.min(cand.len());
        if o.include_ties && take > 0 {
            let kth = cand[take - 1].0;
            while take < cand.len() && cand[take].0 <= kth {
                take += 1;
            }
        }
        let sel = &cand[..take];
        idx_rows.push(sel.iter().map(|c| c.1 as f64 + 1.0).collect());
        d_rows.push(sel.iter().map(|c| c.0).collect());
    }
    if o.include_ties {
        let mk = |rows: Vec<Vec<f64>>, sg: bool, dbl: bool| {
            let n = rows.len();
            let _ = dbl;
            Value::Cell(CellArray::from_dims(vec![n, 1], rows.into_iter().map(|r| row_val(r, sg)).collect()))
        };
        return Ok(Some(vec![mk(idx_rows, false, true), mk(d_rows, single, false)]));
    }
    let cols = idx_rows.first().map_or(k, |r| r.len());
    let (mut idat, mut ddat) = (vec![0.0; nq * cols], vec![0.0; nq * cols]);
    for i in 0..nq {
        for j in 0..cols {
            idat[j * nq + i] = idx_rows[i][j];
            ddat[j * nq + i] = d_rows[i][j];
        }
    }
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![nq, cols], idat)), Value::Num(finish(Matrix::from_dims(vec![nq, cols], ddat), single))]))
}

fn rangesearch(it: &mut Interpreter, a: &[Value]) -> Ret {
    let (x, y) = xy_args(a, "rangesearch")?;
    let single = is_single(x) || is_single(y);
    let r = match &a[2] {
        Value::Num(m) if m.numel() == 1 && m.kind != Kind::Char && m.data[0] >= 0.0 => m.data[0],
        _ => return rt("RADIUS must be a non-negative numeric."),
    };
    let (xs, ys) = (rows_of(x), rows_of(y));
    let o = parse(it, a, "rangesearch", &xs, x.dims[1])?;
    let mut idx_cells = Vec::new();
    let mut d_cells = Vec::new();
    for q in &ys {
        let d = if q.iter().any(|v| v.is_nan()) { vec![f64::NAN; xs.len()] } else { dist_to_all(it, &o.metric, q, &xs, x)? };
        let mut cand: Vec<(f64, usize)> = d.iter().copied().enumerate().filter(|(_, v)| !v.is_nan() && *v <= r).map(|(i, v)| (v, i)).collect();
        if o.sort {
            cand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        }
        idx_cells.push(row_val(cand.iter().map(|c| c.1 as f64 + 1.0).collect(), false));
        d_cells.push(row_val(cand.iter().map(|c| c.0).collect(), single));
    }
    let n = idx_cells.len();
    Ok(Some(vec![Value::Cell(CellArray::from_dims(vec![n, 1], idx_cells)), Value::Cell(CellArray::from_dims(vec![n, 1], d_cells))]))
}

fn crosstab(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("First input must be a real, full, numeric matrix with at least one column.");
    }
    let mut groups: Vec<(Vec<Option<usize>>, Vec<String>)> = Vec::new();
    let mut len = None;
    for v in a {
        let lab = labels_of(v)?;
        let (g, names): (Vec<Option<usize>>, Vec<String>) = match lab {
            Labels::Num(vals, kind) => {
                let mut u: Vec<f64> = vals.iter().copied().filter(|t| !t.is_nan()).collect();
                u.sort_by(|a, b| a.partial_cmp(b).unwrap());
                u.dedup();
                let _ = kind;
                (vals.iter().map(|t| if t.is_nan() { None } else { u.iter().position(|y| y == t) }).collect(), u.iter().map(|&t| fmt_num(t)).collect())
            }
            Labels::Text(items, _) => {
                let mut u: Vec<String> = Vec::new();
                let g = items
                    .iter()
                    .map(|t| {
                        t.as_ref().map(|s| match u.iter().position(|y| y == s) {
                            Some(i) => i,
                            None => {
                                u.push(s.clone());
                                u.len() - 1
                            }
                        })
                    })
                    .collect();
                (g, u)
            }
        };
        if let Some(l) = len {
            if l != g.len() {
                return rt("All arguments must have the same length");
            }
        }
        len = Some(g.len());
        groups.push((g, names));
    }
    let n = len.unwrap_or(0);
    let k = groups.len();
    let dims: Vec<usize> = groups.iter().map(|g| g.1.len()).collect();
    let total: usize = dims.iter().product();
    let mut table = vec![0.0; total];
    let mut used = 0.0;
    for i in 0..n {
        if groups.iter().any(|g| g.0[i].is_none()) {
            continue;
        }
        let mut off = 0;
        let mut stride = 1;
        for (d, g) in groups.iter().enumerate() {
            off += g.0[i].unwrap() * stride;
            stride *= dims[d];
        }
        table[off] += 1.0;
        used += 1.0;
    }
    // marginals
    let mut marg: Vec<Vec<f64>> = dims.iter().map(|&d| vec![0.0; d]).collect();
    for (off, &c) in table.iter().enumerate() {
        let mut rem = off;
        for d in 0..k {
            marg[d][rem % dims[d]] += c;
            rem /= dims[d];
        }
    }
    let (mut chi, mut pv) = (f64::NAN, f64::NAN);
    let df = total as f64 - 1.0 - dims.iter().map(|&d| d as f64 - 1.0).sum::<f64>();
    if k > 1 && df >= 1.0 {
        let mut s = 0.0;
        for (off, &c) in table.iter().enumerate() {
            let mut rem = off;
            let mut e = used;
            for d in 0..k {
                e *= marg[d][rem % dims[d]] / used;
                rem /= dims[d];
            }
            if e > 0.0 {
                s += (c - e) * (c - e) / e;
            }
        }
        chi = s;
        pv = it.call_function("chi2cdf", vec![Value::num(s), Value::num(df), Value::chars("upper")], 1)?[0].scalar()?;
    }
    let tdims = if k == 1 { vec![dims[0], 1] } else { dims.clone() };
    let maxd = dims.iter().copied().max().unwrap_or(0);
    let mut labels = Vec::with_capacity(maxd * k);
    for g in &groups {
        for i in 0..maxd {
            labels.push(match g.1.get(i) {
                Some(s) => Value::chars(s),
                None => Value::empty(),
            });
        }
    }
    Ok(Some(vec![Value::Num(Matrix::from_dims(tdims, table)), Value::num(chi), Value::num(pv), Value::Cell(CellArray::from_dims(vec![maxd, k], labels))]))
}

fn polyconf(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let (p, x) = (a[0].clone(), a[1].clone());
    let s = a.get(2).cloned();
    let mut alpha: Option<Value> = None;
    let mut mu: Option<Value> = None;
    let (mut doobs, mut dosim) = (true, false);
    if a.len() > 3 {
        if a[3].is_text() {
            let rest = &a[3..];
            if rest.len() % 2 != 0 {
                return rt("Invalid parameter name/value pair.");
            }
            for pair in rest.chunks(2) {
                let key = pair[0].text().unwrap_or_default();
                match key.to_lowercase().as_str() {
                    "alpha" => alpha = Some(pair[1].clone()),
                    "mu" => mu = Some(pair[1].clone()),
                    "predopt" => {
                        let t = pair[1].text().unwrap_or_default().to_lowercase();
                        doobs = match ["curve", "observation"].iter().position(|n| !t.is_empty() && n.starts_with(&t)) {
                            Some(i) => i == 1,
                            None => return rt("'PredOpt' value must be 'curve' or 'observation'."),
                        };
                    }
                    "simopt" => {
                        let t = pair[1].text().unwrap_or_default().to_lowercase();
                        dosim = match ["on", "off"].iter().position(|n| !t.is_empty() && n.starts_with(&t)) {
                            Some(i) => i == 0,
                            None => return rt("'SimOpt' value must be 'on' or 'off'."),
                        };
                    }
                    _ => return rt(format!("Invalid parameter name: {key}.")),
                }
            }
        } else {
            alpha = Some(a[3].clone());
            mu = a.get(4).cloned();
        }
    }
    if nargout <= 1 {
        let args = match &mu {
            Some(m) if !matches!(m, Value::Num(mm) if mm.is_empty()) => vec![p, x, Value::empty(), m.clone()],
            _ => vec![p, x],
        };
        return Ok(Some(it.call_function("polyval", args, 1)?.into_iter().take(1).collect()));
    }
    let Some(sv) = s else { return rt("S is required to compute error estimates.") };
    let al = match &alpha {
        None => 0.05,
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        Some(Value::Num(m)) if m.is_empty() => 0.05,
        _ => return rt("ALPHA must be a scalar between 0 and 1."),
    };
    if let Some(Value::Num(m)) = &alpha {
        if !m.is_empty() && (m.numel() != 1 || !(m.data[0] > 0.0 && m.data[0] < 1.0)) {
            return rt("ALPHA must be a scalar between 0 and 1.");
        }
    }
    if !matches!(&sv, Value::Struct(_)) {
        return rt("S matrix must be n+2-by-n where n = length(p).");
    }
    let pv = match &mu {
        Some(m) if !matches!(m, Value::Num(mm) if mm.is_empty()) => it.call_function("polyval", vec![p.clone(), x, sv.clone(), m.clone()], 2)?,
        _ => it.call_function("polyval", vec![p.clone(), x, sv.clone()], 2)?,
    };
    let (y, delta) = (pv[0].clone(), pv[1].clone());
    let Value::Struct(ss) = &sv else { return rt("S matrix must be n+2-by-n where n = length(p).") };
    let df = ss.get("df").map(|v| v.scalar()).transpose()?.unwrap_or(f64::NAN);
    let normr = ss.get("normr").map(|v| v.scalar()).transpose()?.unwrap_or(f64::NAN);
    let Value::Num(dm) = &delta else { return rt("polyval failed") };
    let predvar: Vec<f64> = if doobs {
        dm.data.clone()
    } else {
        let s = normr / df.sqrt();
        dm.data.iter().map(|d| { let dd = d / s; s * (dd * dd - 1.0).sqrt() }).collect()
    };
    let k = p.as_matrix()?.numel() as f64;
    let crit = if dosim {
        let f = it.call_function("finv", vec![Value::num(1.0 - al), Value::num(k), Value::num(df)], 1)?[0].scalar()?;
        (k * f).sqrt()
    } else {
        it.call_function("tinv", vec![Value::num(1.0 - al / 2.0), Value::num(df)], 1)?[0].scalar()?
    };
    let mut out = dm.clone();
    out.data = predvar.iter().map(|v| crit * v).collect();
    Ok(Some(vec![y, Value::Num(out)]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "knnsearch" => knnsearch(it, a),
        "rangesearch" => rangesearch(it, a),
        "crosstab" => crosstab(it, a),
        _ => polyconf(it, a, nargout),
    }
}

/// `polyval(p,x,S,mu)` with the error estimate `delta` (the error-estimate structure of `polyfit`).
pub fn polyval_full(a: &[Value], nargout: usize) -> Ret {
    let (Some(Value::Num(p)), Some(Value::Num(x))) = (a.first(), a.get(1)) else { return rt("Not enough input arguments.") };
    let mu: Option<&Matrix> = match a.get(3) {
        Some(Value::Num(m)) if !m.is_empty() => Some(m),
        _ => None,
    };
    let xs: Vec<f64> = match mu {
        Some(m) if m.numel() >= 2 => x.data.iter().map(|v| (v - m.data[0]) / m.data[1]).collect(),
        Some(_) => return rt("MU must be a two-element vector."),
        None => x.data.clone(),
    };
    let single = is_single(p) || is_single(x);
    let y: Vec<f64> = xs.iter().map(|&xv| p.data.iter().fold(0.0, |acc, &c| acc * xv + c)).collect();
    let yv = Value::Num(finish(Matrix::from_dims(x.dims.clone(), y), single));
    if nargout < 2 {
        return Ok(Some(vec![yv]));
    }
    let n = p.numel();
    let s = a.get(2);
    let Some(Value::Struct(ss)) = s else {
        return match s {
            None | Some(Value::Num(_)) if s.map_or(true, |v| matches!(v, Value::Num(m) if m.is_empty())) => rt("S is required to compute error estimates."),
            _ => rt("S matrix must be n+2-by-n where n = length(p)."),
        };
    };
    let (Some(Value::Num(r)), Some(df), Some(nr)) = (ss.get("R"), ss.get("df"), ss.get("normr")) else { return rt("S must contain the fields R, df and normr.") };
    if r.dims != [n, n] {
        return rt("S matrix must be n+2-by-n where n = length(p).");
    }
    let (df, normr) = (df.scalar()?, nr.scalar()?);
    let scale = normr / df.sqrt();
    let delta: Vec<f64> = xs
        .iter()
        .map(|&xv| {
            // row of the Vandermonde matrix times inv(R): solve R' w = v
            let v: Vec<f64> = (0..n).map(|j| xv.powi((n - 1 - j) as i32)).collect();
            let mut w = vec![0.0; n];
            for j in 0..n {
                let mut sacc = v[j];
                for i in 0..j {
                    sacc -= w[i] * r.data[j * n + i];
                }
                w[j] = sacc / r.data[j * n + j];
            }
            (1.0 + w.iter().map(|t| t * t).sum::<f64>()).sqrt() * scale
        })
        .collect();
    Ok(Some(vec![yv, Value::Num(finish(Matrix::from_dims(x.dims.clone(), delta), single))]))
}
