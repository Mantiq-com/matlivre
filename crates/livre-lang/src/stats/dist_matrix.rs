//! Pairwise distances: `pdist`, `pdist2`, `squareform`, `mahal`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["pdist", "pdist2", "squareform", "mahal"];

const METHODS: &[&str] = &[
    "euclidean", "seuclidean", "cityblock", "chebychev", "mahalanobis", "minkowski", "cosine", "correlation", "spearman", "hamming", "jaccard",
    "squaredEuclidean", "fasteuclidean", "fastseuclidean", "fastsquaredeuclidean",
];

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Metric {
    Euc,
    Sqe,
    Seu,
    Cit,
    Che,
    Mah,
    Min,
    Cos,
    Cor,
    Spe,
    Ham,
    Jac,
}

/// Distance method of a name: Ok(None) when the name is not a built-in one (a function name).
pub fn parse_metric(name: &str, who: &str) -> Res<Option<Metric>> {
    let low = name.to_lowercase();
    let hits: Vec<usize> = METHODS.iter().enumerate().filter(|(_, m)| m.to_lowercase().starts_with(&low) && !low.is_empty()).map(|(i, _)| i).collect();
    if hits.len() > 1 {
        let _ = who;
        return rt(format!("Ambiguous 'DISTANCE' argument:  {name}."));
    }
    if hits.is_empty() {
        return Ok(None);
    }
    Ok(Some(match hits[0] {
        0 | 12 => Metric::Euc,
        1 | 13 => Metric::Seu,
        2 => Metric::Cit,
        3 => Metric::Che,
        4 => Metric::Mah,
        5 => Metric::Min,
        6 => Metric::Cos,
        7 => Metric::Cor,
        8 => Metric::Spe,
        9 => Metric::Ham,
        10 => Metric::Jac,
        _ => Metric::Sqe,
    }))
}

fn rows_of(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1]);
    (0..r).map(|i| (0..c).map(|j| m.data[j * r + i]).collect()).collect()
}

pub fn var_cols(rows: &[Vec<f64>]) -> Vec<f64> {
    let n = rows.len();
    let p = rows.first().map_or(0, |r| r.len());
    (0..p)
        .map(|j| {
            let vals: Vec<f64> = rows.iter().map(|r| r[j]).filter(|v| !v.is_nan()).collect();
            let k = vals.len() as f64;
            if vals.len() < 2 {
                return if vals.is_empty() { f64::NAN } else { 0.0 };
            }
            let m = vals.iter().sum::<f64>() / k;
            let _ = n;
            vals.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (k - 1.0)
        })
        .collect()
}

/// Covariance of the rows (observations), omitting rows with NaN.
pub fn cov_rows(rows: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let ok: Vec<&Vec<f64>> = rows.iter().filter(|r| r.iter().all(|v| !v.is_nan())).collect();
    let n = ok.len() as f64;
    let p = rows.first().map_or(0, |r| r.len());
    let mean: Vec<f64> = (0..p).map(|j| ok.iter().map(|r| r[j]).sum::<f64>() / n).collect();
    let mut c = vec![vec![0.0; p]; p];
    for a in 0..p {
        for b in 0..p {
            c[a][b] = ok.iter().map(|r| (r[a] - mean[a]) * (r[b] - mean[b])).sum::<f64>() / (n - 1.0);
        }
    }
    c
}

/// Upper Cholesky factor R (C = R'R), or None when the matrix is not positive definite.
pub fn chol_upper(c: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let p = c.len();
    let mut r = vec![vec![0.0; p]; p];
    for j in 0..p {
        for i in 0..=j {
            let mut s = c[i][j];
            for k in 0..i {
                s -= r[k][i] * r[k][j];
            }
            if i == j {
                if !(s > 0.0) {
                    return None;
                }
                r[i][i] = s.sqrt();
            } else {
                r[i][j] = s / r[i][i];
            }
        }
    }
    Some(r)
}

/// Inverse of an upper triangular matrix.
fn inv_upper(t: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let p = t.len();
    let mut inv = vec![vec![0.0; p]; p];
    for j in 0..p {
        inv[j][j] = 1.0 / t[j][j];
        for i in (0..j).rev() {
            let mut s = 0.0;
            for k in i + 1..=j {
                s += t[i][k] * inv[k][j];
            }
            inv[i][j] = -s / t[i][i];
        }
    }
    inv
}

fn normalize_rows(rows: &mut [Vec<f64>]) -> bool {
    let mut flag = false;
    for r in rows.iter_mut() {
        let nrm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if nrm == 0.0 {
            flag = true;
            r.iter_mut().for_each(|v| *v = f64::NAN);
        } else {
            r.iter_mut().for_each(|v| *v /= nrm);
        }
    }
    flag
}

/// Everything a metric needs besides the two observations.
pub struct Prep {
    pub metric: Metric,
    pub w: Vec<f64>,
    pub tinv: Vec<Vec<f64>>,
    pub expo: f64,
}

pub fn dist(pr: &Prep, x: &[f64], y: &[f64]) -> f64 {
    match pr.metric {
        Metric::Euc => x.iter().zip(y).map(|(a, b)| (a - b) * (a - b)).sum::<f64>().sqrt(),
        Metric::Sqe => x.iter().zip(y).map(|(a, b)| (a - b) * (a - b)).sum::<f64>(),
        Metric::Seu => x.iter().zip(y).zip(&pr.w).map(|((a, b), w)| (a - b) * (a - b) * w).sum::<f64>().sqrt(),
        Metric::Cit => x.iter().zip(y).map(|(a, b)| (a - b).abs()).sum(),
        Metric::Che => {
            if x.iter().zip(y).any(|(a, b)| (a - b).is_nan()) {
                return f64::NAN;
            }
            x.iter().zip(y).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
        }
        Metric::Mah => {
            let d: Vec<f64> = x.iter().zip(y).map(|(a, b)| a - b).collect();
            let p = d.len();
            let mut s = 0.0;
            for j in 0..p {
                let t: f64 = (0..p).map(|i| d[i] * pr.tinv[i][j]).sum();
                s += t * t;
            }
            s.sqrt()
        }
        Metric::Min => x.iter().zip(y).map(|(a, b)| (a - b).abs().powf(pr.expo)).sum::<f64>().powf(1.0 / pr.expo),
        Metric::Cos | Metric::Cor | Metric::Spe => {
            let d: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
            1.0 - d
        }
        Metric::Ham => x.iter().zip(y).filter(|(a, b)| a != b).count() as f64 / x.len() as f64,
        Metric::Jac => {
            let nz = x.iter().zip(y).filter(|(a, b)| **a != 0.0 || **b != 0.0).count() as f64;
            let diff = x.iter().zip(y).filter(|(a, b)| a != b).count() as f64;
            if nz == 0.0 { 0.0 } else { diff / nz }
        }
    }
}

/// Parameters of the metric for `rows` (the data that fixes the defaults). `extra` is the
/// optional numeric argument (weights, covariance, exponent).
pub fn prepare(metric: Metric, rows: &mut Vec<Vec<f64>>, extra: Option<&Matrix>, who: &str) -> Res<Prep> {
    let p = rows.first().map_or(0, |r| r.len());
    let mut pr = Prep { metric, w: vec![], tinv: vec![], expo: 2.0 };
    match metric {
        Metric::Seu => match extra {
            None => pr.w = var_cols(rows).iter().map(|v| 1.0 / v).collect(),
            Some(e) => {
                let ok = (e.dims.len() == 2 && (e.dims[0] == 1 || e.dims[1] == 1)) && e.numel() == p && e.data.iter().all(|&v| v >= 0.0);
                if !ok {
                    return rt("The weights for the standardized Euclidean metric must be a vector of non-negative values with one element for each column of X.");
                }
                pr.w = e.data.iter().map(|v| 1.0 / (v * v)).collect();
            }
        },
        Metric::Mah => {
            let c = match extra {
                None => cov_rows(rows),
                Some(e) => {
                    if e.dims != [p, p] {
                        return rt("The covariance matrix for the Mahalanobis metric must be a square matrix with the same number of columns as X.");
                    }
                    (0..p).map(|i| (0..p).map(|j| e.data[j * p + i]).collect()).collect()
                }
            };
            let t = match chol_upper(&c) {
                Some(t) if extra.map_or(true, |_| is_symmetric(&c)) => t,
                _ => {
                    return rt(if extra.is_some() {
                        "The covariance matrix for the Mahalanobis metric must be symmetric and positive definite."
                    } else {
                        "The covariance matrix of X must be positive definite for the Mahalanobis metric; the data may be singular."
                    })
                }
            };
            pr.tinv = inv_upper(&t);
        }
        Metric::Min => {
            if let Some(e) = extra {
                if !(e.numel() == 1 && e.data[0] > 0.0) {
                    return rt("The exponent for the Minkowski metric must be a positive scalar.");
                }
                pr.expo = e.data[0];
            }
            if pr.expo.is_infinite() {
                pr.metric = Metric::Che;
            } else if pr.expo == 2.0 {
                pr.metric = Metric::Euc;
            } else if pr.expo == 1.0 {
                pr.metric = Metric::Cit;
            }
        }
        Metric::Cos => {
            normalize_rows(rows);
        }
        Metric::Cor => {
            for r in rows.iter_mut() {
                let m = r.iter().sum::<f64>() / r.len() as f64;
                r.iter_mut().for_each(|v| *v -= m);
            }
            normalize_rows(rows);
        }
        Metric::Spe => {
            for r in rows.iter_mut() {
                let (rk, _) = super::descr::tiedrank_vec(r, false, false);
                let m = (r.len() as f64 + 1.0) / 2.0;
                *r = rk.iter().map(|v| v - m).collect();
            }
            normalize_rows(rows);
        }
        _ => {}
    }
    let _ = who;
    Ok(pr)
}

fn is_symmetric(c: &[Vec<f64>]) -> bool {
    let p = c.len();
    (0..p).all(|i| (0..p).all(|j| c[i][j] == c[j][i]))
}

/// A cell of data for `pdist`: metric, optional numeric argument, or a distance function.
enum Spec {
    Builtin(Metric, Option<Matrix>),
    User(FuncValue, Vec<Value>),
}

fn parse_dist_args(it: &mut Interpreter, rest: &[Value], who: &str) -> Res<Spec> {
    let _ = it;
    match rest.first() {
        None => Ok(Spec::Builtin(Metric::Euc, None)),
        Some(Value::Func(f)) => Ok(Spec::User(f.clone(), rest[1..].to_vec())),
        Some(v) if v.is_text() => {
            let name = v.text()?;
            match parse_metric(&name, who)? {
                Some(m) => {
                    let extra = match rest.get(1) {
                        Some(Value::Num(e)) if e.kind != Kind::Char && matches!(m, Metric::Seu | Metric::Mah | Metric::Min) => Some(e.clone()),
                        _ => None,
                    };
                    Ok(Spec::Builtin(m, extra))
                }
                None => rt(format!("The distance function '{name}' was not found.")),
            }
        }
        _ => rt("The distance must be a character vector, string scalar, or function handle."),
    }
}

fn real_matrix<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(_) => rt("Invalid data type."),
        None => rt("Not enough input arguments."),
    }
}

pub fn pdist_core(it: &mut Interpreter, x: &Matrix, spec_args: &[Value]) -> Res<(Vec<f64>, bool)> {
    let single = is_single(x);
    let spec = parse_dist_args(it, spec_args, "pdist")?;
    if matches!(spec, Spec::Builtin(..)) && x.is_complex() {
        return rt("PDIST does not accept complex data for built-in distances.");
    }
    if x.dims.len() > 2 {
        return rt("N-D arrays are not supported.");
    }
    let n = x.dims[0];
    if n < 2 {
        return Ok((vec![], single));
    }
    let mut rows = rows_of(x);
    let out = match spec {
        Spec::Builtin(m, extra) => {
            let pr = prepare(m, &mut rows, extra.as_ref(), "pdist")?;
            let mut out = Vec::with_capacity(n * (n - 1) / 2);
            for i in 0..n - 1 {
                for j in i + 1..n {
                    out.push(dist(&pr, &rows[i], &rows[j]));
                }
            }
            out
        }
        Spec::User(f, extra) => {
            let p = x.dims[1];
            let mut out = Vec::with_capacity(n * (n - 1) / 2);
            for i in 0..n - 1 {
                let m = n - i - 1;
                let xi = Matrix::from_dims(vec![1, p], rows[i].clone());
                let mut data = vec![0.0; m * p];
                for (jj, r) in rows[i + 1..].iter().enumerate() {
                    for c in 0..p {
                        data[c * m + jj] = r[c];
                    }
                }
                let xj = Matrix::from_dims(vec![m, p], data);
                let mut args = vec![Value::Num(xi), Value::Num(xj)];
                args.extend(extra.iter().cloned());
                let res = it.call_value(&f, args, 1)?;
                let Some(Value::Num(d)) = res.first() else { return rt("The distance function must return a numeric column vector.") };
                if d.numel() != m {
                    return rt(format!("The distance function must return a column vector with {m} rows."));
                }
                out.extend(d.data.iter().copied());
            }
            out
        }
    };
    Ok((out, single))
}

fn pdist(it: &mut Interpreter, a: &[Value]) -> Ret {
    let x = real_matrix(a, 0)?;
    let (out, single) = pdist_core(it, x, &a[1.min(a.len())..])?;
    let n = out.len();
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, n], out), single))]))
}

fn squareform(a: &[Value]) -> Ret {
    let y = match a.first() {
        Some(Value::Num(m)) if m.kind != Kind::Char && m.dims.len() <= 2 => m,
        Some(_) => return rt("Y must be a numeric or logical vector or matrix."),
        None => return rt("Not enough input arguments."),
    };
    let (m, n) = (y.dims[0], y.dims[1]);
    let is_vec = m == 1 || n == 1;
    let dir = match a.get(1) {
        Some(v) if !matches!(v, Value::Num(x) if x.is_empty()) => {
            let t = v.text().unwrap_or_default().to_lowercase();
            if !t.is_empty() && "tovector".starts_with(&t) {
                1
            } else if !t.is_empty() && "tomatrix".starts_with(&t) {
                0
            } else {
                return rt(format!("'{}' is not a valid value for the DIRECTION argument.  Valid values are: 'tovector', 'tomatrix'.", v.text().unwrap_or_default()));
            }
        }
        _ => if is_vec { 0 } else { 1 },
    };
    if dir == 0 {
        if !is_vec {
            return rt("Y must be a vector to convert it to a distance matrix");
        }
        let len = y.numel();
        let k = ((2.0 * len as f64).sqrt()).ceil() as usize;
        if k * (k.saturating_sub(1)) / 2 != len {
            return rt("The size of the vector Y is not correct");
        }
        let mut z = vec![0.0; k * k];
        let mut idx = 0;
        for j in 0..k {
            for i in j + 1..k {
                z[j * k + i] = y.data[idx];
                z[i * k + j] = y.data[idx];
                idx += 1;
            }
        }
        let mut out = Matrix::from_dims(vec![k, k], z);
        out.kind = y.kind;
        Ok(Some(vec![Value::Num(out)]))
    } else {
        if m != n || (0..n).any(|i| y.data[i * n + i] != 0.0) {
            return rt("The distance matrix Z must be square with 0 along the diagonal.");
        }
        let mut v = Vec::new();
        for j in 0..n {
            for i in j + 1..n {
                v.push(y.data[j * n + i]);
            }
        }
        let len = v.len();
        let mut out = Matrix::from_dims(vec![1, len], v);
        out.kind = y.kind;
        Ok(Some(vec![Value::Num(out)]))
    }
}

fn pdist2(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let x = real_matrix(a, 0)?;
    let y = real_matrix(a, 1)?;
    if x.dims[1] != y.dims[1] {
        return rt("X and Y must have the same number of columns.");
    }
    let single = is_single(x) || is_single(y);
    let mut rest = &a[2..];
    let mut dist_args: Vec<Value> = Vec::new();
    // distance, optional argument
    if let Some(v) = rest.first() {
        if v.is_text() {
            let t = v.text()?.to_lowercase();
            let is_opt = ["smallest", "largest"].iter().any(|o| o.starts_with(&t)) && !t.is_empty() && parse_metric(&t, "pdist2").ok().flatten().is_none();
            if !is_opt {
                dist_args.push(v.clone());
                rest = &rest[1..];
                if let Some(Value::Num(e)) = rest.first() {
                    if e.kind != Kind::Char {
                        dist_args.push(Value::Num(e.clone()));
                        rest = &rest[1..];
                    }
                }
            }
        } else if let Value::Func(_) = v {
            dist_args.push(v.clone());
            rest = &rest[1..];
        }
    }
    let mut smallest: Option<usize> = None;
    let mut largest: Option<usize> = None;
    let mut i = 0;
    while i < rest.len() {
        let key = if rest[i].is_text() { rest[i].text()?.to_lowercase() } else { String::new() };
        let val = rest.get(i + 1);
        let get_k = |v: Option<&Value>, name: &str| -> Res<usize> {
            match v {
                Some(Value::Num(m)) if m.numel() == 1 && m.data[0] >= 1.0 && m.data[0].fract() == 0.0 => Ok(m.data[0] as usize),
                _ => rt(format!("{name} must be a positive integer.")),
            }
        };
        if !key.is_empty() && "smallest".starts_with(&key) {
            smallest = Some(get_k(val, "SMALLEST")?);
        } else if !key.is_empty() && "largest".starts_with(&key) {
            largest = Some(get_k(val, "LARGEST")?);
        } else {
            return rt(format!("Invalid parameter name: {}.", rest[i].text().unwrap_or_default()));
        }
        i += 2;
    }
    let (nx, ny) = (x.dims[0], y.dims[0]);
    let xr = rows_of(x);
    let yr = rows_of(y);
    let mut d = vec![0.0; nx * ny];
    match parse_dist_args(it, &dist_args, "pdist2") {
        Err(e) if e.message.starts_with("The distance function") => {
            let name = dist_args[0].text().unwrap_or_default();
            return rt(format!("Unrecognized 'DISTANCE' argument: {name}."));
        }
        Err(e) => return Err(e),
        Ok(Spec::Builtin(m, extra)) => {
            // Defaults (variance, covariance) come from X; cosine and the like normalise both.
            let mut xs = xr.clone();
            let mut ys = yr.clone();
            let pr = prepare(m, &mut xs, extra.as_ref(), "pdist2")?;
            match m {
                Metric::Cos | Metric::Cor | Metric::Spe => {
                    let _ = prepare(m, &mut ys, None, "pdist2")?;
                }
                _ => {}
            }
            for j in 0..ny {
                for i in 0..nx {
                    d[j * nx + i] = dist(&pr, &xs[i], &ys[j]);
                }
            }
        }
        Ok(Spec::User(f, extra)) => {
            let p = x.dims[1];
            for j in 0..ny {
                let mut args = vec![Value::Num(Matrix::from_dims(vec![1, p], yr[j].clone())), Value::Num(x.clone())];
                let _ = &mut args;
                // distfun(ZI, ZJ): ZI is one row of X, ZJ the rows of Y
                let mut out = vec![0.0; nx];
                for i in 0..nx {
                    let zi = Matrix::from_dims(vec![1, p], xr[i].clone());
                    let mut data = vec![0.0; ny * p];
                    for (jj, r) in yr.iter().enumerate() {
                        for c in 0..p {
                            data[c * ny + jj] = r[c];
                        }
                    }
                    let zj = Matrix::from_dims(vec![ny, p], data);
                    let mut av = vec![Value::Num(zi), Value::Num(zj)];
                    av.extend(extra.iter().cloned());
                    let res = it.call_value(&f, av, 1)?;
                    let Some(Value::Num(r)) = res.first() else { return rt("The distance function must return a numeric vector.") };
                    if r.numel() != ny {
                        return rt(format!("The distance function must return a vector with {ny} elements."));
                    }
                    out[i] = r.data[j];
                }
                d[j * nx..(j + 1) * nx].copy_from_slice(&out);
            }
        }
    }
    let k = smallest.or(largest);
    let mut outs = Vec::new();
    match k {
        None => {
            outs.push(Value::Num(finish(Matrix::from_dims(vec![nx, ny], d), single)));
            if nargout > 1 {
                outs.push(Value::Num(Matrix::from_dims(vec![nx, ny], (0..nx * ny).map(|t| (t % nx + 1) as f64).collect())));
            }
        }
        Some(k) => {
            let kk = k.min(nx);
            let mut dv = Vec::with_capacity(kk * ny);
            let mut iv = Vec::with_capacity(kk * ny);
            for j in 0..ny {
                let mut idx: Vec<usize> = (0..nx).collect();
                idx.sort_by(|&p, &q| {
                    let (a, b) = (d[j * nx + p], d[j * nx + q]);
                    let o = match (a.is_nan(), b.is_nan()) {
                        (true, true) => std::cmp::Ordering::Equal,
                        (true, false) => std::cmp::Ordering::Greater,
                        (false, true) => std::cmp::Ordering::Less,
                        _ => if largest.is_some() { b.partial_cmp(&a).unwrap() } else { a.partial_cmp(&b).unwrap() },
                    };
                    o
                });
                for &t in idx.iter().take(kk) {
                    dv.push(d[j * nx + t]);
                    iv.push((t + 1) as f64);
                }
            }
            outs.push(Value::Num(finish(Matrix::from_dims(vec![kk, ny], dv), single)));
            if nargout > 1 {
                outs.push(Value::Num(Matrix::from_dims(vec![kk, ny], iv)));
            }
        }
    }
    Ok(Some(outs))
}

fn mahal(a: &[Value]) -> Ret {
    let y = real_matrix(a, 0)?;
    let x = real_matrix(a, 1)?;
    let (rx, cx) = (x.dims[0], x.dims[1]);
    let (ry, cy) = (y.dims[0], y.dims[1]);
    if cx != cy {
        return rt("Requires the inputs to have the same number of columns.");
    }
    if rx < cx {
        return rt("The number of rows of X must exceed the number of columns.");
    }
    let single = is_single(x) || is_single(y);
    let xr = rows_of(x);
    let yr = rows_of(y);
    let c = cov_rows(&xr);
    let mean: Vec<f64> = (0..cx).map(|j| xr.iter().map(|r| r[j]).sum::<f64>() / rx as f64).collect();
    let t = match chol_upper(&c) {
        Some(t) => t,
        None => {
            // Rank-deficient data: the QR-based formula divides by zero.
            return Ok(Some(vec![Value::Num(finish(Matrix::filled(ry, 1, f64::NAN), single))]));
        }
    };
    let ti = inv_upper(&t);
    let out: Vec<f64> = yr
        .iter()
        .map(|r| {
            let d: Vec<f64> = r.iter().zip(&mean).map(|(a, b)| a - b).collect();
            let mut s = 0.0;
            for j in 0..cx {
                let u: f64 = (0..cx).map(|i| d[i] * ti[i][j]).sum();
                s += u * u;
            }
            s
        })
        .collect();
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![ry, 1], out), single))]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "pdist" => pdist(it, a),
        "pdist2" => pdist2(it, a, nargout),
        "squareform" => squareform(a),
        _ => mahal(a),
    }
}
