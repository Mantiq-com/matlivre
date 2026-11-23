//! Small closed-form functions: `combnk`, `ff2n`, `fullfact`, `mnpdf`, `mvnpdf`, `mvtpdf`, `leverage`,
//! `ridge`, `nancov`, `tabulate` (numeric), `ecdf`.

use super::dist_matrix::chol_upper;
use super::dist::normcdf_std;
use super::dist::norminv_std;
use super::util::*;
use crate::builtins::Ret;
use crate::builtins_specfun::gammaln;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["combnk", "ff2n", "fullfact", "mnpdf", "mvnpdf", "mvtpdf", "leverage", "ridge", "nancov", "tabulate", "ecdf"];

fn num<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(_) => rt("Input must be numeric."),
        None => rt("Not enough input arguments."),
    }
}

fn combnk(a: &[Value]) -> Ret {
    let v = num(a, 0)?;
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let kf = a[1].scalar()?;
    if kf.fract() != 0.0 || kf < 0.0 {
        return rt("Size arguments must be real integers.");
    }
    let k = kf as usize;
    if v.dims.len() > 2 || (v.dims[0] != 1 && v.dims[1] != 1) {
        return rt("First argument has to be a vector.");
    }
    let vals = v.data.clone();
    let n = vals.len();
    let col_in = v.dims[1] == 1 && v.dims[0] > 1;
    let mk = |rows: Vec<Vec<f64>>, kk: usize| -> Matrix {
        let r = rows.len();
        let mut d = vec![0.0; r * kk];
        for (i, row) in rows.iter().enumerate() {
            for j in 0..kk {
                d[j * r + i] = row[j];
            }
        }
        let mut m = Matrix::from_dims(vec![r, kk], d);
        m.kind = v.kind;
        m
    };
    let out: Matrix = if n == k {
        let mut m = Matrix::from_dims(vec![1, n], vals.clone());
        m.kind = v.kind;
        m
    } else if n == k + 1 {
        let rows: Vec<Vec<f64>> = (0..n).map(|i| (0..n).filter(|&j| j != n - 1 - i).map(|j| vals[j]).collect()).collect();
        mk(rows, n - 1)
    } else if k == 1 {
        let mut m = Matrix::from_dims(vec![n, 1], vals.clone());
        m.kind = v.kind;
        m
    } else if n < 17 && (k > 3 || (n as i64 - k as i64) < 4) {
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for code in 0u32..(1u32 << n) {
            if code.count_ones() as usize == k {
                let row: Vec<f64> = (0..n).filter(|&c| (code >> (n - 1 - c)) & 1 == 1).map(|c| vals[c]).collect();
                rows.push(row);
            }
        }
        mk(rows, k)
    } else {
        // lexicographic recursion
        fn rec(vals: &[f64], k: usize) -> Vec<Vec<f64>> {
            if k == 0 {
                return vec![vec![]];
            }
            let n = vals.len();
            if k > n {
                return vec![];
            }
            let mut out = Vec::new();
            for i in 0..=n - k {
                for tail in rec(&vals[i + 1..], k - 1) {
                    let mut r = vec![vals[i]];
                    r.extend(tail);
                    out.push(r);
                }
            }
            out
        }
        if k == 0 {
            let mut m = Matrix::from_dims(vec![1, 0], vec![]);
            m.kind = v.kind;
            m
        } else if k > n {
            Matrix::from_dims(vec![0, k], vec![])
        } else {
            mk(rec(&vals, k), k)
        }
    };
    let _ = col_in;
    Ok(Some(vec![Value::Num(out)]))
}

fn ff2n(a: &[Value]) -> Ret {
    let m = num(a, 0)?;
    if m.numel() != 1 {
        return rt("Size inputs must be scalar.");
    }
    let nf = m.data[0];
    if nf < 0.0 || nf.fract() != 0.0 {
        return rt("Size inputs must be integers.");
    }
    let n = nf as usize;
    let rows = 1usize << n;
    let mut d = vec![0.0; rows * n];
    for r in 0..rows {
        for c in 0..n {
            d[c * rows + r] = ((r >> (n - 1 - c)) & 1) as f64;
        }
    }
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![rows, n], d))]))
}

fn fullfact(a: &[Value]) -> Ret {
    let m = num(a, 0)?;
    if m.dims.len() > 2 || (m.dims[0] > 1 && m.dims[1] > 1) {
        return rt("Requires a vector input.");
    }
    if m.data.iter().any(|&v| v.fract() != 0.0 || v < 1.0) {
        return rt("The input values must be positive integers.");
    }
    let lv: Vec<usize> = m.data.iter().map(|&v| v as usize).collect();
    let total: usize = lv.iter().product();
    let cols = lv.len();
    let mut d = vec![0.0; total * cols];
    let mut ncycles = total;
    for (k, &l) in lv.iter().enumerate() {
        let nreps = total / ncycles;
        ncycles /= l;
        for r in 0..total {
            d[k * total + r] = ((r / nreps) % l + 1) as f64;
        }
    }
    let mut out = Matrix::from_dims(vec![total, cols], d);
    if m.kind == Kind::Single {
        out = out.cast_to(Kind::Single);
    }
    Ok(Some(vec![Value::Num(out)]))
}

fn mnpdf(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let (mut x, mut p) = (num(a, 0)?.clone(), num(a, 1)?.clone());
    if p.dims[1] == 1 && p.dims[0] > 1 && (p.data.iter().sum::<f64>() - 1.0).abs() <= p.dims[0] as f64 * f64::EPSILON {
        p = p.transpose();
        if x.dims[1] == 1 && x.dims[0] > 1 {
            x = x.transpose();
        }
    }
    let single = is_single(&x) || is_single(&p);
    let (mut m, k) = (x.dims[0], x.dims[1]);
    if k < 1 {
        return rt("X must have at least one column.");
    }
    let (mm, kk) = (p.dims[0], p.dims[1]);
    if kk != k {
        return rt("P and X must have the same number of columns.");
    }
    let xr = |i: usize, j: usize| x.data[j * x.dims[0] + i];
    let pr = |i: usize, j: usize| p.data[j * mm + i];
    let (xi, pi): (Box<dyn Fn(usize) -> usize>, Box<dyn Fn(usize) -> usize>);
    if mm == 1 {
        pi = Box::new(|_| 0);
        xi = Box::new(|i| i);
    } else if m == 1 {
        m = mm;
        xi = Box::new(|_| 0);
        pi = Box::new(|i| i);
    } else if mm != m {
        return rt("P and X must have the same number of rows, or either can be a row vector.");
    } else {
        xi = Box::new(|i| i);
        pi = Box::new(|i| i);
    }
    let mut out = vec![0.0; m];
    for i in 0..m {
        let (r, s) = (xi(i), pi(i));
        let n: f64 = (0..k).map(|j| xr(r, j)).sum();
        let xbad = (0..k).any(|j| xr(r, j) < 0.0 || xr(r, j) != xr(r, j).round());
        let pbad = (0..k).any(|j| pr(s, j) < 0.0 || pr(s, j) > 1.0) || ((0..k).map(|j| pr(s, j)).sum::<f64>() - 1.0).abs() > k as f64 * f64::EPSILON;
        let nbad = n < 0.0 || n.round() != n;
        if pbad || nbad {
            out[i] = f64::NAN;
        } else if xbad {
            out[i] = 0.0;
        } else {
            let xlogp: f64 = (0..k).map(|j| if xr(r, j) > 0.0 { xr(r, j) * pr(s, j).ln() } else { 0.0 }).sum();
            out[i] = (gammaln(n + 1.0) - (0..k).map(|j| gammaln(xr(r, j) + 1.0)).sum::<f64>() + xlogp).exp();
        }
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![m, 1], out), single))]))
}

fn rows_of(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1]);
    (0..r).map(|i| (0..c).map(|j| m.data[j * r + i]).collect()).collect()
}

fn mvnpdf(a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Requires at least one input.");
    }
    let x = num(a, 0)?;
    if x.dims.len() != 2 {
        return rt("X must be a matrix.");
    }
    let (mut n, mut d) = (x.dims[0], x.dims[1]);
    let single = a.iter().any(|v| matches!(v, Value::Num(m) if is_single(m)));
    let mut xr = rows_of(x);
    match a.get(1) {
        Some(Value::Num(mu)) if !mu.is_empty() => {
            if mu.numel() == 1 {
                xr.iter_mut().for_each(|r| r.iter_mut().for_each(|v| *v -= mu.data[0]));
            } else {
                let (n2, d2) = (mu.dims[0], mu.dims[1]);
                if d2 != d {
                    return rt("X and MU must have the same number of columns.");
                }
                let mr = rows_of(mu);
                if n2 == 1 {
                    xr.iter_mut().for_each(|r| r.iter_mut().zip(&mr[0]).for_each(|(v, m)| *v -= m));
                } else if n2 == n {
                    for (r, m) in xr.iter_mut().zip(&mr) {
                        r.iter_mut().zip(m).for_each(|(v, mv)| *v -= mv);
                    }
                } else if n == 1 {
                    n = n2;
                    xr = mr.iter().map(|m| xr[0].iter().zip(m).map(|(a, b)| a - b).collect()).collect();
                } else {
                    return rt("X or MU must be a row vector, or X and MU must have the same number of rows.");
                }
            }
        }
        _ => {}
    }
    let (xr_inv, log_sqrt_det): (Vec<Vec<f64>>, Vec<f64>);
    match a.get(2) {
        Some(Value::Num(s)) if !s.is_empty() => {
            if s.dims.len() == 2 {
                let (mut s0, s1) = (s.dims[0], s.dims[1]);
                let diag = s0 == 1 && s1 > 1;
                if diag {
                    s0 = s1;
                }
                if d == 1 && x.numel() > 1 && s0 == n {
                    // a column of scalars: treat as rows of dimension n
                    let t: Vec<Vec<f64>> = vec![xr.iter().map(|r| r[0]).collect()];
                    xr = t;
                    d = n;
                    n = 1;
                }
                if s0 != s1 && !diag {
                    return rt("SIGMA must be a square covariance matrix, the diagonal vector of a covariance matrix, or a 3-D array of covariance matrices or covariance diagonal vectors.");
                }
                if s0 != d {
                    return rt("SIGMA must be a square matrix with size equal to the number of columns in X, or a row vector with length equal to the number of columns in X.");
                }
                if diag {
                    if s.data.iter().any(|&v| v <= 0.0) {
                        return rt("All elements of diagonal SIGMA must be positive.");
                    }
                    let r: Vec<f64> = s.data.iter().map(|v| v.sqrt()).collect();
                    xr_inv = xr.iter().map(|row| row.iter().zip(&r).map(|(a, b)| a / b).collect()).collect();
                    log_sqrt_det = vec![r.iter().map(|v| v.ln()).sum()];
                } else {
                    let sm = rows_of(s);
                    let sym = (0..d).all(|i| (0..d).all(|j| sm[i][j] == sm[j][i]));
                    let Some(r) = chol_upper(&sm).filter(|_| sym) else { return rt("SIGMA must be a square, symmetric, positive definite matrix.") };
                    // x / R: solve y R = x, i.e. R' y' = x'
                    let solve = |row: &Vec<f64>| -> Vec<f64> {
                        let mut y = vec![0.0; d];
                        for j in 0..d {
                            let mut sacc = row[j];
                            for i in 0..j {
                                sacc -= y[i] * r[i][j];
                            }
                            y[j] = sacc / r[j][j];
                        }
                        y
                    };
                    xr_inv = xr.iter().map(solve).collect();
                    log_sqrt_det = vec![(0..d).map(|i| r[i][i].ln()).sum()];
                }
            } else if s.dims.len() == 3 && s.dims[0] == d && s.dims[1] == d && s.dims[2] == n {
                let mut res = Vec::new();
                for pg in 0..n {
                    let sm: Vec<Vec<f64>> = (0..d).map(|i| (0..d).map(|j| s.data[pg * d * d + j * d + i]).collect()).collect();
                    let Some(r) = chol_upper(&sm) else { return rt("SIGMA must be a square, symmetric, positive definite matrix.") };
                    let mut y = vec![0.0; d];
                    for j in 0..d {
                        let mut sacc = xr[pg][j];
                        for i in 0..j {
                            sacc -= y[i] * r[i][j];
                        }
                        y[j] = sacc / r[j][j];
                    }
                    let lsd: f64 = (0..d).map(|i| r[i][i].ln()).sum();
                    res.push((-0.5 * y.iter().map(|v| v * v).sum::<f64>() - lsd - d as f64 * (2.0 * std::f64::consts::PI).ln() / 2.0).exp());
                }
                let len = res.len();
                return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![len, 1], res), single))]));
            } else {
                return rt("SIGMA must be a square covariance matrix, the diagonal vector of a covariance matrix, or a 3-D array of covariance matrices or covariance diagonal vectors.");
            }
        }
        _ => {
            if d == 1 && x.numel() > 1 {
                xr = vec![xr.iter().map(|r| r[0]).collect()];
                d = xr[0].len();
                n = 1;
            }
            xr_inv = xr.clone();
            log_sqrt_det = vec![0.0];
        }
    }
    let _ = n;
    let out: Vec<f64> = xr_inv.iter().map(|row| (-0.5 * row.iter().map(|v| v * v).sum::<f64>() - log_sqrt_det[0] - d as f64 * (2.0 * std::f64::consts::PI).ln() / 2.0).exp()).collect();
    let len = out.len();
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![len, 1], out), single))]))
}

fn mvtpdf(a: &[Value]) -> Ret {
    if a.len() < 3 {
        return rt("Requires three inputs.");
    }
    let x = num(a, 0)?;
    let c = num(a, 1)?;
    let df = num(a, 2)?;
    if x.dims.len() != 2 {
        return rt("X must be a matrix.");
    }
    let single = is_single(x) || is_single(c) || is_single(df);
    let (mut n, mut d) = (x.dims[0], x.dims[1]);
    let mut xr = rows_of(x);
    let is_vec = n == 1 || d == 1;
    if is_vec && c.dims[0] == n {
        xr = vec![xr.iter().map(|r| r[0]).collect()];
        n = 1;
        d = xr[0].len();
    }
    if c.dims[0] != c.dims[1] {
        return rt("C must be a square matrix.");
    }
    if c.dims != [d, d] {
        return rt("C must be a square matrix with size equal to the number of columns in X.");
    }
    let cm = rows_of(c);
    let s: Vec<f64> = (0..d).map(|i| cm[i][i].sqrt()).collect();
    let cn: Vec<Vec<f64>> = (0..d).map(|i| (0..d).map(|j| if s.iter().any(|&v| v != 1.0) { cm[i][j] / (s[i] * s[j]) } else { cm[i][j] }).collect()).collect();
    let sym = (0..d).all(|i| (0..d).all(|j| cn[i][j] == cn[j][i]));
    let Some(r) = chol_upper(&cn).filter(|_| sym) else { return rt("C must be symmetric and positive definite.") };
    if !(df.numel() == 1 || (df.dims.len() == 2 && (df.dims[0] == 1 || df.dims[1] == 1) && df.numel() == n)) {
        return rt("DF must be a scalar or a vector with one element for each row in X.");
    }
    if df.data.iter().any(|&v| v <= 0.0) {
        return rt("DF must be positive.");
    }
    let lsd: f64 = (0..d).map(|i| r[i][i].ln()).sum();
    let df_at = |i: usize| if df.numel() == 1 { df.data[0] } else { df.data[i] };
    let out: Vec<f64> = (0..n)
        .map(|i| {
            let mut y = vec![0.0; d];
            for j in 0..d {
                let mut sacc = xr[i][j];
                for t in 0..j {
                    sacc -= y[t] * r[t][j];
                }
                y[j] = sacc / r[j][j];
            }
            let ss: f64 = y.iter().map(|v| v * v).sum();
            let v = df_at(i);
            let log_numer = -((v + d as f64) / 2.0) * (1.0 + ss / v).ln();
            let log_denom = lsd + (d as f64 / 2.0) * (v * std::f64::consts::PI).ln();
            (gammaln((v + d as f64) / 2.0) - gammaln(v / 2.0) + log_numer - log_denom).exp()
        })
        .collect();
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![n, 1], out), single))]))
}

fn leverage(it: &mut Interpreter, a: &[Value]) -> Ret {
    let data = num(a, 0)?;
    let model = a.get(1).cloned().unwrap_or_else(|| Value::chars("linear"));
    let single = is_single(data);
    let (n, p) = (data.dims[0], data.dims[1]);
    let wasnan: Vec<bool> = (0..n).map(|i| (0..p).any(|j| data.data[j * n + i].is_nan())).collect();
    let keep: Vec<usize> = (0..n).filter(|&i| !wasnan[i]).collect();
    let nk = keep.len();
    let mut d = vec![0.0; nk * p];
    for (r, &i) in keep.iter().enumerate() {
        for j in 0..p {
            d[j * nk + r] = data.data[j * n + i];
        }
    }
    let x = it.call_function("x2fx", vec![Value::Num(Matrix::from_dims(vec![nk, p], d)), model], 1)?;
    let qr = it.call_function("qr", vec![x[0].clone(), Value::num(0.0)], 2)?;
    let Value::Num(q) = &qr[0] else { return rt("qr failed") };
    let k = q.dims[1];
    let h: Vec<f64> = (0..nk).map(|i| (0..k).map(|j| q.data[j * nk + i].powi(2)).sum()).collect();
    let mut out = vec![f64::NAN; n];
    for (r, &i) in keep.iter().enumerate() {
        out[i] = h[r];
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![n, 1], out), single))]))
}

pub(crate) fn lstsq(a: &[Vec<f64>], b: &[f64], p: usize) -> Vec<f64> {
    let n = a.len();
    let (q, r, perm) = super::regress::qr_pivot(&a.to_vec(), n, p);
    let rq: Vec<Vec<f64>> = r.iter().take(p).map(|row| row[..p].to_vec()).collect();
    let qty: Vec<f64> = (0..p).map(|c| (0..n).map(|i| q[i][c] * b[i]).sum()).collect();
    let bp = super::regress::solve_upper(&rq, &qty, p);
    let mut out = vec![0.0; p];
    for (j, &pc) in perm.iter().enumerate() {
        out[pc] = bp[j];
    }
    out
}

fn ridge(a: &[Value]) -> Ret {
    if a.len() < 3 {
        return rt("Requires at least three input arguments.");
    }
    let (y, x, k) = (num(a, 0)?, num(a, 1)?, num(a, 2)?);
    let unscale = match a.get(3) {
        None => false,
        Some(Value::Num(m)) if m.is_empty() => false,
        Some(Value::Num(m)) if m.numel() == 1 && m.data[0] == 1.0 => false,
        Some(Value::Num(m)) if m.numel() == 1 && m.data[0] == 0.0 => true,
        _ => return rt("The scaling flag must be 0 or 1."),
    };
    let single = is_single(x) || is_single(y);
    let (n, p) = (x.dims[0], x.dims[1]);
    if y.dims[0] != n {
        return rt("The number of rows in Y must equal the number of rows in X.");
    }
    if y.dims[1] != 1 {
        return rt("Y must be a column vector.");
    }
    let keep: Vec<usize> = (0..n).filter(|&i| !y.data[i].is_nan() && (0..p).all(|j| !x.data[j * n + i].is_nan())).collect();
    let nn = keep.len();
    let xr: Vec<Vec<f64>> = keep.iter().map(|&i| (0..p).map(|j| x.data[j * n + i]).collect()).collect();
    let yv: Vec<f64> = keep.iter().map(|&i| y.data[i]).collect();
    let mx: Vec<f64> = (0..p).map(|j| xr.iter().map(|r| r[j]).sum::<f64>() / nn as f64).collect();
    let mut stdx: Vec<f64> = (0..p).map(|j| (xr.iter().map(|r| (r[j] - mx[j]).powi(2)).sum::<f64>() / (nn as f64 - 1.0)).sqrt()).collect();
    let idx: Vec<usize> = (0..p).filter(|&j| stdx[j].abs() < f64::EPSILON.sqrt()).collect();
    for &j in &idx {
        stdx[j] = 1.0;
    }
    let mut z: Vec<Vec<f64>> = xr.iter().map(|r| (0..p).map(|j| (r[j] - mx[j]) / stdx[j]).collect()).collect();
    for &j in &idx {
        for r in z.iter_mut() {
            r[j] = 1.0;
        }
    }
    let nk = k.numel();
    let mut cols: Vec<Vec<f64>> = Vec::new();
    for t in 0..nk {
        let mut zp = z.clone();
        for i in 0..p {
            let mut row = vec![0.0; p];
            row[i] = k.data[t].sqrt();
            zp.push(row);
        }
        let mut yp = yv.clone();
        yp.extend(vec![0.0; p]);
        cols.push(lstsq(&zp, &yp, p));
    }
    let ymean = yv.iter().sum::<f64>() / nn as f64;
    let mut out: Vec<f64> = Vec::new();
    if unscale {
        for c in cols.iter_mut() {
            for j in 0..p {
                c[j] /= stdx[j];
            }
        }
        for c in &cols {
            out.push(ymean - (0..p).map(|j| mx[j] * c[j]).sum::<f64>());
            out.extend(c.iter().copied());
        }
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![p + 1, nk], out), single))]));
    }
    for c in &cols {
        out.extend(c.iter().copied());
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![p, nk], out), single))]))
}

fn local_cov(rows: &[Vec<f64>], n: usize, mle: bool) -> Vec<f64> {
    let m = rows.len();
    let mean: Vec<f64> = (0..n).map(|j| rows.iter().map(|r| r[j]).sum::<f64>() / m as f64).collect();
    let denom = if mle || m == 1 { m as f64 } else { m as f64 - 1.0 };
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            c[j * n + i] = rows.iter().map(|r| (r[i] - mean[i]) * (r[j] - mean[j])).sum::<f64>() / denom;
        }
    }
    c
}

fn nancov(a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    let x0 = num(a, 0)?;
    let single = is_single(x0);
    let mut rest: Vec<&Value> = a[1..].iter().collect();
    let mut pairwise = false;
    if let Some(l) = rest.last() {
        if l.is_text() {
            let t = l.text()?.to_lowercase();
            match ["pairwise", "complete"].iter().position(|c| !t.is_empty() && c.starts_with(&t)) {
                Some(j) => pairwise = j == 0,
                None => return rt(format!("Invalid argument '{}'.", l.text()?)),
            }
            rest.pop();
        }
    }
    let mut mle = false;
    if let Some(Value::Num(m)) = rest.last() {
        if m.numel() == 1 && (m.data[0] == 0.0 || m.data[0] == 1.0) {
            mle = m.data[0] == 1.0;
            rest.pop();
        }
    }
    let mut x = x0.clone();
    let mut scalarxy = false;
    if rest.len() > 1 {
        return rt("Too many arguments.");
    }
    if let Some(Value::Num(y)) = rest.first() {
        if x.numel() != y.numel() {
            return rt("The lengths of X and Y must match.");
        }
        scalarxy = x.numel() == 1 && y.numel() == 1;
        let mut d = x.data.clone();
        d.extend(y.data.iter().copied());
        x = Matrix::from_dims(vec![x.numel(), 2], d);
    } else if x.dims.len() > 2 {
        return rt("Inputs must be 2-D.");
    }
    if (x.dims[0] == 1 || x.dims[1] == 1) && !scalarxy && x.dims.len() == 2 {
        x = Matrix::from_dims(vec![x.numel(), 1], x.data.clone());
    }
    let (m, n) = (x.dims[0], x.dims[1]);
    if x.is_empty() {
        let out = if m == 0 && n == 0 { Matrix::scalar(f64::NAN) } else { Matrix::filled(n, n, f64::NAN) };
        return Ok(Some(vec![Value::Num(finish(out, single))]));
    }
    let rows = rows_of(&x);
    let has_nan = rows.iter().any(|r| r.iter().any(|v| v.is_nan()));
    let c: Vec<f64> = if !pairwise || !has_nan {
        let ok: Vec<Vec<f64>> = rows.into_iter().filter(|r| r.iter().all(|v| !v.is_nan())).collect();
        if ok.is_empty() { vec![f64::NAN; n * n] } else { local_cov(&ok, n, mle) }
    } else {
        let mut c = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let ok: Vec<(f64, f64)> = rows.iter().filter(|r| !r[i].is_nan() && !r[j].is_nan()).map(|r| (r[i], r[j])).collect();
                let k = ok.len() as f64;
                if k == 0.0 {
                    c[j * n + i] = f64::NAN;
                    continue;
                }
                let (mi, mj) = (ok.iter().map(|p| p.0).sum::<f64>() / k, ok.iter().map(|p| p.1).sum::<f64>() / k);
                let denom = if mle || k == 1.0 { k } else { k - 1.0 };
                c[j * n + i] = ok.iter().map(|p| (p.0 - mi) * (p.1 - mj)).sum::<f64>() / denom;
            }
        }
        c
    };
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![n, n], c), single))]))
}

fn tabulate(a: &[Value]) -> Ret {
    let x = match a.first() {
        Some(Value::Num(m)) if m.kind != Kind::Logical && m.kind != Kind::Char => m,
        Some(_) => return rt("This form of tabulate is not implemented."),
        None => return rt("Not enough input arguments."),
    };
    if x.dims.len() > 2 || (x.dims[0] > 1 && x.dims[1] > 1) {
        return rt("Requires a vector input.");
    }
    if x.is_empty() {
        return Ok(Some(vec![Value::Num(Matrix::from_dims(vec![0, 1], vec![]))]));
    }
    let single = is_single(x);
    let vals: Vec<f64> = x.data.iter().copied().filter(|v| !v.is_nan()).collect();
    let lo = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let k = (hi - lo) as usize + 1;
    let mut counts = vec![0.0; k];
    for v in &vals {
        counts[(v - lo) as usize] += 1.0;
    }
    let total = vals.len() as f64;
    let mut d = Vec::with_capacity(3 * k);
    d.extend((0..k).map(|i| lo + i as f64));
    d.extend(counts.iter().copied());
    d.extend(counts.iter().map(|c| 100.0 * c / total));
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![k, 3], d), single))]))
}

fn ecdf(a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(xm)) = a.first() else { return rt("Not enough input arguments.") };
    let single = is_single(xm);
    let rest = &a[1..];
    if rest.len() % 2 == 1 {
        return rt("Invalid parameter name or missing value.");
    }
    if xm.dims.len() > 2 || xm.is_empty() || (xm.dims[0] > 1 && xm.dims[1] > 1) {
        return rt("Sample data must be a vector, or a two-column matrix.");
    }
    let mut alpha = 0.05;
    let mut function = "cdf".to_string();
    let mut cens: Option<Vec<f64>> = None;
    let mut freq: Option<Vec<f64>> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        match ["alpha", "function", "bounds", "censoring", "frequency", "npoints"].iter().position(|n| !key.is_empty() && n.starts_with(&key)) {
            Some(0) => alpha = pair[1].scalar().unwrap_or(0.05),
            Some(1) => {
                let t = pair[1].text().unwrap_or_default().to_lowercase();
                match ["cdf", "survivor", "cumulative hazard"].iter().position(|n| !t.is_empty() && n.starts_with(&t)) {
                    Some(i) => function = ["cdf", "survivor", "cumhaz"][i].into(),
                    None => return rt(format!("'{}' is not a valid value for the Function argument.", pair[1].text().unwrap_or_default())),
                }
            }
            Some(2) | Some(5) => {}
            Some(3) => cens = pair[1].as_matrix().ok().map(|m| m.data.clone()),
            Some(4) => freq = pair[1].as_matrix().ok().map(|m| m.data.clone()),
            _ => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
        }
    }
    let n0 = xm.numel();
    let cn: Vec<f64> = cens.unwrap_or_else(|| vec![0.0; n0]);
    let fq: Vec<f64> = freq.unwrap_or_else(|| vec![1.0; n0]);
    let keep: Vec<usize> = (0..n0).filter(|&i| !xm.data[i].is_nan()).collect();
    if keep.is_empty() {
        return rt("Sample data contains only missing data.");
    }
    let mut order = keep.clone();
    order.sort_by(|&i, &j| xm.data[i].partial_cmp(&xm.data[j]).unwrap());
    let total: f64 = order.iter().map(|&i| fq[i]).sum();
    let mut times: Vec<f64> = Vec::new();
    let mut surv = vec![1.0];
    let mut cum_h = vec![0.0];
    let mut gw = vec![0.0];
    let mut hv = 0.0;
    let mut hvar = vec![0.0];
    let mut atrisk = total;
    let mut s = 1.0;
    let mut h = 0.0;
    let mut g = 0.0;
    let mut i = 0;
    while i < order.len() {
        let t = xm.data[order[i]];
        let (mut d, mut c) = (0.0, 0.0);
        let mut j = i;
        while j < order.len() && xm.data[order[j]] == t {
            if cn[order[j]] != 0.0 {
                c += fq[order[j]];
            } else {
                d += fq[order[j]];
            }
            j += 1;
        }
        if d > 0.0 {
            s *= 1.0 - d / atrisk;
            h += d / atrisk;
            g += d / (atrisk * (atrisk - d));
            hv += d / (atrisk * atrisk);
            hvar.push(hv);
            times.push(t);
            surv.push(s);
            cum_h.push(h);
            gw.push(g);
        }
        atrisk -= d + c;
        i = j;
    }
    let first = xm.data[order[0]];
    let mut xs = vec![first];
    xs.extend(times.iter().copied());
    let z = norminv_std(1.0 - alpha / 2.0);
    let m = xs.len();
    let (f, lo, up): (Vec<f64>, Vec<f64>, Vec<f64>) = match function.as_str() {
        "cdf" => {
            let f: Vec<f64> = surv.iter().map(|v| 1.0 - v).collect();
            let se: Vec<f64> = (0..m).map(|k| surv[k] * gw[k].sqrt()).collect();
            let lo: Vec<f64> = (0..m).map(|k| if k == 0 || k == m - 1 && f[k] >= 1.0 { f64::NAN } else { (f[k] - z * se[k]).max(0.0) }).collect();
            let up: Vec<f64> = (0..m).map(|k| if k == 0 || k == m - 1 && f[k] >= 1.0 { f64::NAN } else { (f[k] + z * se[k]).min(1.0) }).collect();
            (f, lo, up)
        }
        "survivor" => {
            let se: Vec<f64> = (0..m).map(|k| surv[k] * gw[k].sqrt()).collect();
            let lo: Vec<f64> = (0..m).map(|k| if k == 0 || surv[k] <= 0.0 { f64::NAN } else { (surv[k] - z * se[k]).max(0.0) }).collect();
            let up: Vec<f64> = (0..m).map(|k| if k == 0 || surv[k] <= 0.0 { f64::NAN } else { (surv[k] + z * se[k]).min(1.0) }).collect();
            (surv.clone(), lo, up)
        }
        _ => {
            let lo: Vec<f64> = (0..m).map(|k| if k == 0 { f64::NAN } else { (cum_h[k] - z * hvar[k].sqrt()).max(0.0) }).collect();
            let up: Vec<f64> = (0..m).map(|k| if k == 0 { f64::NAN } else { cum_h[k] + z * hvar[k].sqrt() }).collect();
            (cum_h.clone(), lo, up)
        }
    };
    let col = |v: Vec<f64>| Value::Num(finish(Matrix::from_dims(vec![v.len(), 1], v), single));
    let mut outs = vec![col(f), col(xs)];
    if nargout > 2 {
        outs.push(col(lo));
        outs.push(col(up));
    }
    let _ = normcdf_std;
    Ok(Some(outs))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "combnk" => combnk(a),
        "ff2n" => ff2n(a),
        "fullfact" => fullfact(a),
        "mnpdf" => mnpdf(a),
        "mvnpdf" => mvnpdf(a),
        "mvtpdf" => mvtpdf(a),
        "leverage" => leverage(it, a),
        "ridge" => ridge(a),
        "nancov" => nancov(a),
        "tabulate" => tabulate(a),
        _ => ecdf(a, nargout),
    }
}
