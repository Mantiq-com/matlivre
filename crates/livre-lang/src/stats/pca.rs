//! Principal components and related decompositions: `pca`, `pcacov`, `pcares`, `cmdscale`,
//! `corrcov`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["pca", "pcacov", "pcares", "cmdscale", "corrcov"];

type Mat = Vec<Vec<f64>>; // row-major rows

fn to_rows(m: &Matrix) -> Mat {
    let (r, c) = (m.dims[0], m.dims[1]);
    (0..r).map(|i| (0..c).map(|j| m.data[j * r + i]).collect()).collect()
}

fn from_rows(rows: &Mat, ncols: usize) -> Matrix {
    let r = rows.len();
    let mut data = vec![0.0; r * ncols];
    for (i, row) in rows.iter().enumerate() {
        for j in 0..ncols {
            data[j * r + i] = row[j];
        }
    }
    Matrix::from_dims(vec![r, ncols], data)
}

fn first_matrix(v: Vec<Value>, k: usize) -> Res<Matrix> {
    match v.into_iter().nth(k) {
        Some(Value::Num(m)) => Ok(m),
        _ => rt("internal error: expected a matrix"),
    }
}

/// Economy SVD of the rows: (U rows, singular values, V rows)
fn svd_econ(it: &mut Interpreter, x: &Mat, ncols: usize, econ: bool) -> Res<(Mat, Vec<f64>, Mat)> {
    let m = from_rows(x, ncols);
    let args = if econ { vec![Value::Num(m), Value::chars("econ")] } else { vec![Value::Num(m)] };
    let out = it.call_function("svd", args, 3)?;
    let u = first_matrix(out.clone(), 0)?;
    let s = first_matrix(out.clone(), 1)?;
    let v = first_matrix(out, 2)?;
    let k = s.dims[0].min(s.dims[1]);
    let sv: Vec<f64> = (0..k).map(|i| s.data[i * s.dims[0] + i]).collect();
    Ok((to_rows(&u), sv, to_rows(&v)))
}

/// Eigen-decomposition of a symmetric matrix, eigenvalues descending: (values, vectors as columns in rows)
fn eig_sym(it: &mut Interpreter, c: &Mat) -> Res<(Vec<f64>, Mat)> {
    let p = c.len();
    let m = from_rows(c, p);
    let out = it.call_function("eig", vec![Value::Num(m)], 2)?;
    let v = to_rows(&first_matrix(out.clone(), 0)?);
    let d = first_matrix(out, 1)?;
    let mut vals: Vec<(f64, usize)> = (0..p).map(|i| (d.data[i * p + i], i)).collect();
    vals.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let sorted: Vec<f64> = vals.iter().map(|x| x.0).collect();
    let vec_rows: Mat = (0..p).map(|r| vals.iter().map(|(_, ci)| v[r][*ci]).collect()).collect();
    Ok((sorted, vec_rows))
}

fn on_off(v: &Value, pname: &str) -> Res<bool> {
    match v {
        Value::Num(m) if m.numel() == 1 && m.kind == Kind::Logical => Ok(m.data[0] != 0.0),
        Value::Num(m) if m.numel() == 1 && m.kind != Kind::Char => Ok(m.data[0] != 0.0),
        _ if v.is_text() => {
            let t = v.text()?.to_lowercase();
            if !t.is_empty() && "on".starts_with(&t) && t.len() >= 2 || t == "true" {
                Ok(true)
            } else if !t.is_empty() && "off".starts_with(&t) && t.len() >= 2 || t == "false" {
                Ok(false)
            } else {
                rt(format!("'{}' is not a valid value for the '{pname}' argument.  Valid values are: 'on', 'off'.", v.text()?))
            }
        }
        _ => rt(format!("'{pname}' must be logical or 'on' or 'off'.")),
    }
}

fn choose(v: &Value, valid: &[&str], pname: &str) -> Res<usize> {
    let t = if v.is_text() { v.text()? } else { String::new() };
    let low = t.to_lowercase();
    let hits: Vec<usize> = valid.iter().enumerate().filter(|(_, c)| !low.is_empty() && c.starts_with(&low)).map(|(i, _)| i).collect();
    if hits.len() == 1 {
        return Ok(hits[0]);
    }
    let list = valid.iter().map(|c| format!("'{c}'")).collect::<Vec<_>>().join(", ");
    rt(format!("'{t}' is not a valid value for the '{pname}' argument.  Valid values are: {list}."))
}

pub struct PcaOut {
    pub coeff: Matrix,
    pub score: Matrix,
    pub latent: Matrix,
    pub tsq: Matrix,
    pub explained: Matrix,
    pub mu: Matrix,
}

fn mk_col(v: Vec<f64>) -> Matrix {
    let n = v.len();
    Matrix::from_dims(vec![n, 1], v)
}

fn tsquared(score: &Mat, latent: &[f64], dof: f64, p: usize, ncols_score: usize) -> Vec<f64> {
    let r = dof.min(p as f64);
    let q = if dof > 1.0 {
        let thr = dof.max(p as f64) * eps_of(latent[0]);
        let q = latent.iter().filter(|&&l| l > thr).count();
        let _ = r;
        q
    } else {
        0
    };
    let q = q.min(ncols_score);
    score.iter().map(|row| (0..q).map(|j| (row[j] / latent[j].sqrt()).powi(2)).sum()).collect()
}

pub fn pca_core(it: &mut Interpreter, a: &[Value]) -> Res<PcaOut> {
    let x = match a.first() {
        Some(Value::Num(m)) if m.kind != Kind::Char && !matches!(m.kind, Kind::Int(_) | Kind::Logical) => m,
        Some(Value::Num(m)) if m.kind == Kind::Char => return rt("Value of X must be a numeric array."),
        Some(Value::Num(_)) => return rt("Value of X must be a floating point array."),
        Some(_) => return rt("Value of X must be a numeric array."),
        None => return rt("Not enough input arguments."),
    };
    if x.dims.len() > 2 {
        return rt("Value of X must be a matrix.");
    }
    let (n, p) = (x.dims[0], x.dims[1]);
    let single = is_single(x);
    let rest = &a[1..];
    if rest.len() % 2 == 1 {
        return rt("Wrong number of arguments.");
    }
    let (mut algo, mut centered, mut economy, mut ncomp, mut rows) = (0usize, true, true, p as f64, 0usize);
    let (mut weights, mut vweights): (Vec<f64>, Option<Vec<f64>>) = (vec![1.0; n], None);
    let mut vw_variance = false;
    let (mut algo_set, mut rows_set) = (false, false);
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        let names = ["algorithm", "centered", "economy", "numcomponents", "rows", "weights", "variableweights", "coeff0", "score0", "options"];
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !key.is_empty() && nm.starts_with(&key)).map(|(i, _)| i).collect();
        let idx = names.iter().position(|nm| *nm == key).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match idx {
            Some(0) => {
                algo = choose(&pair[1], &["svd", "eig", "als"], "Algorithm")?;
                algo_set = true;
            }
            Some(1) => centered = on_off(&pair[1], "Centered")?,
            Some(2) => economy = on_off(&pair[1], "Economy")?,
            Some(3) => {
                ncomp = match &pair[1] {
                    Value::Num(m) if m.numel() == 1 => m.data[0],
                    _ => f64::NAN,
                }
            }
            Some(4) => {
                rows = choose(&pair[1], &["complete", "pairwise", "all"], "Rows")?;
                rows_set = true;
            }
            Some(5) => {
                match &pair[1] {
                    Value::Num(m) if m.numel() == n && (m.dims[0] == 1 || m.dims[1] == 1) => weights = m.data.clone(),
                    _ => return rt(format!("Observation weights must be a vector of length {n}.")),
                };
            }
            Some(6) => {
                if pair[1].is_text() {
                    choose(&pair[1], &["variance"], "VariableWeights")?;
                    vw_variance = true;
                } else if let Value::Num(m) = &pair[1] {
                    if m.numel() == p && (m.dims[0] == 1 || m.dims[1] == 1) {
                        vweights = Some(m.data.clone());
                    } else {
                        return rt(format!("Variable weights must be the character vector or string scalar 'variance', or a vector of length {p}."));
                    }
                } else {
                    return rt(format!("Variable weights must be the character vector or string scalar 'variance', or a vector of length {p}."));
                }
            }
            Some(_) => {}
            None => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
        }
    }
    let _ = (algo_set, rows_set);
    if !x.is_empty() && !(ncomp.fract() == 0.0 && ncomp >= 1.0 && ncomp <= p as f64) {
        return rt(format!("The value of the 'NumComponents' must be a positive scalar integer and no greater than {p}."));
    }
    if algo == 2 {
        return rt("The 'als' algorithm is not implemented.");
    }
    if algo == 0 && rows == 1 {
        algo = 1;
    }
    let xr = to_rows(x);
    if vw_variance {
        // 1 ./ weighted variance of each column (NaN omitted)
        let mut v = Vec::new();
        for j in 0..p {
            let (mut sw, mut sx) = (0.0, 0.0);
            for i in 0..n {
                if !xr[i][j].is_nan() {
                    sw += weights[i];
                    sx += weights[i] * xr[i][j];
                }
            }
            let m = sx / sw;
            let mut ss = 0.0;
            let mut sw2 = 0.0;
            for i in 0..n {
                if !xr[i][j].is_nan() {
                    ss += weights[i] * (xr[i][j] - m).powi(2);
                    sw2 += weights[i] * weights[i];
                }
            }
            // weighted variance with the reliability-weights correction
            let denom = sw - sw2 / sw;
            v.push(1.0 / (ss / denom));
        }
        vweights = Some(v);
    }
    let vw = vweights.unwrap_or_else(|| vec![1.0; p]);
    if weights.iter().any(|&w| w <= 0.0) || vw.iter().any(|&w| w <= 0.0) {
        return rt("Weights must be positive.");
    }
    let mut zero = |dims: Vec<usize>| Matrix::from_dims(dims.clone(), vec![0.0; dims.iter().product()]);
    if x.is_empty() {
        let p_or_zero = if economy { 0 } else { p };
        let mut c = zero(vec![p, p_or_zero]);
        for i in 0..p.min(p_or_zero) {
            c.data[i * p + i] = 1.0;
        }
        return Ok(PcaOut { coeff: c, score: zero(vec![n, p_or_zero]), latent: zero(vec![p_or_zero, 1]), tsq: zero(vec![n, 1]), explained: zero(vec![0, 0]), mu: zero(vec![0, 0]) });
    }
    let nan_all = x.data.iter().all(|v| v.is_nan());
    if nan_all {
        let s = Matrix::scalar(f64::NAN);
        return Ok(PcaOut { coeff: s.clone(), score: s.clone(), latent: s.clone(), tsq: s.clone(), explained: s.clone(), mu: s });
    }
    if n == 1 && p == 1 {
        let v = x.data[0];
        let nc = if centered { 0.0 } else { 1.0 };
        return Ok(PcaOut {
            coeff: Matrix::scalar(1.0),
            score: Matrix::scalar(nc * v),
            latent: Matrix::scalar(nc * v * v),
            tsq: Matrix::scalar(nc),
            explained: Matrix::scalar(100.0),
            mu: Matrix::scalar(if centered { v } else { 0.0 }),
        });
    }
    let was_nan: Vec<bool> = xr.iter().map(|r| r.iter().any(|v| v.is_nan())).collect();
    let n_nan = was_nan.iter().filter(|&&b| b).count();
    if rows == 2 && algo != 2 {
        if n_nan > 0 {
            return rt("Raw data contains NaN missing value while 'Rows' option is set to 'all'. Consider using 'complete' or 'pairwise' option, or using 'als' algorithm instead.");
        }
        rows = 0;
    }
    let cen = if centered { 1.0 } else { 0.0 };
    let dof = if rows == 0 {
        (n as f64 - cen - n_nan as f64).max(0.0)
    } else {
        let mut best: f64 = 0.0;
        for i in 0..p {
            for j in 0..p {
                if i != j {
                    let c = (0..n).filter(|&r| !xr[r][i].is_nan() && !xr[r][j].is_nan()).count() as f64;
                    best = best.max(c);
                }
            }
        }
        best - cen
    };
    // weighted column means (NaN omitted)
    let mu: Vec<f64> = if centered {
        (0..p)
            .map(|j| {
                let (mut sw, mut sx) = (0.0, 0.0);
                for i in 0..n {
                    if !xr[i][j].is_nan() {
                        sw += weights[i];
                        sx += weights[i] * xr[i][j];
                    }
                }
                sx / sw
            })
            .collect()
    } else {
        vec![0.0; p]
    };
    let mut xc: Mat = xr.iter().map(|r| r.iter().zip(&mu).map(|(a, m)| a - m).collect()).collect();
    let (coeff_m, score_m, latent, tsq): (Mat, Mat, Vec<f64>, Vec<f64>);
    let mut ncols_c;
    if algo == 1 {
        // eigenvalue algorithm on the weighted covariance
        let sw: Vec<f64> = weights.iter().map(|w| w.sqrt()).collect();
        let sv: Vec<f64> = vw.iter().map(|w| w.sqrt()).collect();
        let xw: Mat = xc.iter().enumerate().map(|(i, r)| r.iter().enumerate().map(|(j, v)| v * sw[i] * sv[j]).collect()).collect();
        let d = cen;
        let mut c = vec![vec![0.0; p]; p];
        if rows == 1 {
            for i in 0..p {
                for j in 0..=i {
                    let ok: Vec<usize> = (0..n).filter(|&r| !xw[r][i].is_nan() && !xw[r][j].is_nan()).collect();
                    let denom = (ok.len() as f64 - d).max(0.0);
                    let v = ok.iter().map(|&r| xw[r][i] * xw[r][j]).sum::<f64>() / denom;
                    c[i][j] = v;
                    c[j][i] = v;
                }
            }
        } else {
            let okrows: Vec<usize> = (0..n).filter(|&r| !was_nan[r]).collect();
            let denom = (okrows.len() as f64 - d).max(0.0);
            for i in 0..p {
                for j in 0..p {
                    c[i][j] = okrows.iter().map(|&r| xw[r][i] * xw[r][j]).sum::<f64>() / denom;
                }
            }
        }
        let (mut ev, vecs) = eig_sym(it, &c)?;
        // coeff = coeff ./ PhiSqrt' : divide each row i by sv[i]
        let mut co2: Mat = vecs.clone();
        for (i, r) in co2.iter_mut().enumerate() {
            r.iter_mut().for_each(|v| *v /= sv[i]);
        }
        let mut keep = p;
        if dof < p as f64 {
            let d0 = dof.max(0.0) as usize;
            if economy {
                keep = d0;
            } else {
                for e in ev.iter_mut().skip(d0) {
                    *e = 0.0;
                }
            }
        }
        let thr = eps_of(*ev.first().unwrap()) * ev.len() as f64;
        for e in ev.iter_mut() {
            if *e < 0.0 && e.abs() < thr {
                *e = 0.0;
            }
        }
        if ev.iter().any(|&e| e < 0.0) {
            return rt("The covariance matrix is not positive semi-definite.");
        }
        ncols_c = keep;
        // score = x / coeff' (x already centered, original weights not applied)
        let ck: Mat = co2.iter().map(|r| r[..keep].to_vec()).collect();
        // solve score * ck' = xc  =>  score = xc * pinv(ck')
        let xcm = from_rows(&xc.iter().map(|r| r.iter().map(|v| if v.is_nan() { 0.0 } else { *v }).collect()).collect(), p);
        let ckm = from_rows(&ck, keep);
        let out = it.call_function("mrdivide", vec![Value::Num(xcm), Value::Num(ckm.transpose())], 1)?;
        let sc = to_rows(&first_matrix(out, 0)?);
        let lat: Vec<f64> = ev[..keep.max(0)].to_vec();
        let mut scm = sc;
        for (i, row) in scm.iter_mut().enumerate() {
            if was_nan[i] {
                row.iter_mut().for_each(|v| *v = f64::NAN);
            }
        }
        let ts = tsquared(&scm, &lat, n as f64, p, keep);
        coeff_m = ck;
        score_m = scm;
        latent = lat;
        tsq = ts;
    } else {
        // SVD on the complete rows
        let okrows: Vec<usize> = (0..n).filter(|&r| !was_nan[r]).collect();
        let w_ok: Vec<f64> = okrows.iter().map(|&r| weights[r].sqrt()).collect();
        let sv: Vec<f64> = vw.iter().map(|w| w.sqrt()).collect();
        let xw: Mat = okrows.iter().enumerate().map(|(k, &r)| xc[r].iter().enumerate().map(|(j, v)| v * w_ok[k] * sv[j]).collect()).collect();
        let (u, s, v) = svd_econ(it, &xw, p, economy)?;
        let nr = okrows.len();
        let kk = s.len();
        let _ = nr;
        let lat_all: Vec<f64> = s.iter().map(|x| x * x / dof).collect();
        let mut sc_ok: Mat = u.iter().enumerate().map(|(k, row)| row.iter().take(kk).enumerate().map(|(j, val)| val / w_ok[k] * s[j]).collect()).collect();
        let mut co: Mat = v.iter().enumerate().map(|(i, row)| row.iter().map(|val| val / sv[i]).collect()).collect();
        let mut lat = lat_all;
        let tsq_ok = {
            let ts = tsquared(&sc_ok, &lat, dof, p, kk);
            ts
        };
        let mut ts_full = vec![f64::NAN; n];
        let mut sc_full: Mat = vec![vec![f64::NAN; kk]; n];
        for (k, &r) in okrows.iter().enumerate() {
            ts_full[r] = tsq_ok[k];
            sc_full[r] = sc_ok[k].clone();
        }
        sc_ok = sc_full;
        let mut kcols = co.first().map_or(0, |r| r.len());
        if dof < p as f64 {
            let d0 = dof.max(0.0) as usize;
            if economy {
                for r in co.iter_mut() {
                    r.truncate(d0);
                }
                for r in sc_ok.iter_mut() {
                    r.truncate(d0);
                }
                lat.truncate(d0);
                kcols = d0;
            } else {
                for r in sc_ok.iter_mut() {
                    for v in r.iter_mut().skip(d0) {
                        *v = 0.0;
                    }
                }
                for v in lat.iter_mut().skip(d0) {
                    *v = 0.0;
                }
            }
        }
        ncols_c = kcols;
        coeff_m = co;
        score_m = sc_ok;
        latent = lat;
        tsq = ts_full;
    }
    let _ = &mut xc;
    // NumComponents
    let (mut coeff_m, mut score_m) = (coeff_m, score_m);
    let explained: Vec<f64> = {
        let tot: f64 = latent.iter().sum();
        latent.iter().map(|l| 100.0 * l / tot).collect()
    };
    if ncomp < dof {
        let k = ncomp as usize;
        for r in coeff_m.iter_mut() {
            r.truncate(k);
        }
        for r in score_m.iter_mut() {
            r.truncate(k);
        }
        ncols_c = ncols_c.min(k);
    }
    // sign convention: the largest element of each coefficient column is positive
    let nc = coeff_m.first().map_or(0, |r| r.len());
    for j in 0..nc {
        let mut best = 0usize;
        for i in 0..p {
            if coeff_m[i][j].abs() > coeff_m[best][j].abs() {
                best = i;
            }
        }
        let sgn = if coeff_m[best][j] > 0.0 { 1.0 } else if coeff_m[best][j] < 0.0 { -1.0 } else { 0.0 };
        for r in coeff_m.iter_mut() {
            r[j] *= sgn;
        }
        for r in score_m.iter_mut() {
            if j < r.len() {
                r[j] *= sgn;
            }
        }
    }
    let _ = ncols_c;
    let nl = latent.len();
    let ne = explained.len();
    let sc_cols = score_m.first().map_or(0, |r| r.len());
    let fin = |m: Matrix| finish(m, single);
    Ok(PcaOut {
        coeff: fin(from_rows(&coeff_m, nc)),
        score: fin(from_rows(&score_m, sc_cols)),
        latent: fin(Matrix::from_dims(vec![nl, 1], latent)),
        tsq: fin(mk_col(tsq)),
        explained: fin(Matrix::from_dims(vec![ne, 1], explained)),
        mu: fin(Matrix::from_dims(vec![1, p], mu)),
    })
}

fn pca(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let o = pca_core(it, a)?;
    let all = vec![o.coeff, o.score, o.latent, o.tsq, o.explained, o.mu];
    Ok(Some(all.into_iter().take(nargout.max(1)).map(Value::Num).collect()))
}

fn check_cov(it: &mut Interpreter, v: &Matrix, who: &str) -> Res<()> {
    let bad = |m: &str| -> Res<()> { rt(m.to_string()) };
    let msg = if who == "pcacov" { "Covariance matrix must be square, symmetric, and positive semi-definite." } else { "C must be a square, symmetric, positive semidefinite matrix." };
    if v.dims.len() != 2 || v.dims[0] != v.dims[1] || v.is_empty() {
        return bad(msg);
    }
    let p = v.dims[0];
    let c = to_rows(v);
    let dmax = (0..p).map(|i| c[i][i].abs()).fold(0.0, f64::max);
    let tol = eps_of(dmax) * p as f64;
    for i in 0..p {
        for j in 0..i {
            if (c[i][j] - c[j][i]).abs() > 10.0 * tol.max(f64::EPSILON * c[i][j].abs().max(c[j][i].abs())) {
                return bad(msg);
            }
        }
    }
    let sym: Mat = (0..p).map(|i| (0..p).map(|j| (c[i][j] + c[j][i]) / 2.0).collect()).collect();
    let (ev, _) = eig_sym(it, &sym)?;
    if ev.iter().any(|&e| e < -tol) {
        return bad(msg);
    }
    Ok(())
}

fn pcacov(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(v)) = a.first() else { return rt("Not enough input arguments.") };
    check_cov(it, v, "pcacov")?;
    let single = is_single(v);
    let p = v.dims[0];
    let rows = to_rows(v);
    let (_, s, vv) = svd_econ(it, &rows, p, false)?;
    let tot: f64 = s.iter().sum();
    let expl: Vec<f64> = s.iter().map(|x| 100.0 * x / tot).collect();
    let mut co = vv;
    for j in 0..p {
        let mut best = 0;
        for i in 0..p {
            if co[i][j].abs() > co[best][j].abs() {
                best = i;
            }
        }
        let sg = co[best][j].signum();
        for r in co.iter_mut() {
            r[j] *= sg;
        }
    }
    let outs = vec![
        Value::Num(finish(from_rows(&co, p), single)),
        Value::Num(finish(Matrix::from_dims(vec![p, 1], s), single)),
        Value::Num(finish(Matrix::from_dims(vec![p, 1], expl), single)),
    ];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn pcares(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(x)) = a.first() else { return rt("Not enough input arguments.") };
    let ndim_v = match a.get(1) {
        Some(Value::Num(m)) => m,
        _ => return rt("Not enough input arguments."),
    };
    let (n, p) = (x.dims[0], x.dims[1]);
    let r = (n as i64 - 1).min(p as i64).max(0) as usize;
    if ndim_v.numel() != 1 || ndim_v.data[0] > p as f64 {
        return rt("NDIM must be a scalar less than or equal to the number of columns in X.");
    }
    let ndim = (ndim_v.data[0].min(r as f64)).max(0.0) as usize;
    let o = pca_core(it, &[Value::Num(x.clone())])?;
    let xr = to_rows(x);
    let mean: Vec<f64> = (0..p).map(|j| xr.iter().map(|r| r[j]).sum::<f64>() / n as f64).collect();
    let (cr, sr) = (to_rows(&o.coeff), to_rows(&o.score));
    let mut rec = vec![vec![0.0; p]; n];
    for i in 0..n {
        for j in 0..p {
            let mut s = mean[j];
            for k in 0..ndim {
                s += sr[i][k] * cr[j][k];
            }
            rec[i][j] = s;
        }
    }
    let res: Mat = (0..n).map(|i| (0..p).map(|j| xr[i][j] - rec[i][j]).collect()).collect();
    let single = is_single(x);
    let outs = vec![Value::Num(finish(from_rows(&res, p), single)), Value::Num(finish(from_rows(&rec, p), single))];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn cmdscale(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(d0)) = a.first() else { return rt("Not enough input arguments.") };
    let single = is_single(d0);
    let (mut n, m) = (d0.dims[0], d0.dims[1]);
    let del = 10.0 * if single { f32::EPSILON as f64 } else { f64::EPSILON };
    let mut d: Mat;
    if n == 1 {
        n = ((2.0 * m as f64).sqrt()).ceil() as usize;
        if n * (n.saturating_sub(1)) / 2 == m && d0.data.iter().all(|&v| v >= 0.0) {
            d = vec![vec![0.0; n]; n];
            let mut k = 0;
            for j in 0..n {
                for i in j + 1..n {
                    d[i][j] = d0.data[k];
                    d[j][i] = d0.data[k];
                    k += 1;
                }
            }
        } else {
            return rt("Not a valid dissimilarity or distance matrix.");
        }
    } else if n == m && {
        let dm = to_rows(d0);
        let mx = dm.iter().flatten().cloned().fold(0.0, f64::max);
        dm.iter().all(|r| r.iter().all(|&v| v >= 0.0)) && (0..n).all(|i| (0..n).all(|j| (dm[i][j] - dm[j][i]).abs() <= del * mx))
    } {
        d = to_rows(d0);
        if (0..n).all(|i| d[i][i] < del) {
        } else if (0..n).all(|i| (d[i][i] - 1.0).abs() < del) && d.iter().flatten().all(|&v| v < 1.0 + del) {
            d = d.iter().map(|r| r.iter().map(|v| (1.0 - v).sqrt()).collect()).collect();
        } else {
            return rt("Not a valid dissimilarity or distance matrix.");
        }
    } else {
        return rt("Not a valid dissimilarity or distance matrix.");
    }
    let p = match a.get(1) {
        None => n,
        Some(Value::Num(pm)) if pm.numel() == 1 && pm.data[0] >= 1.0 && pm.data[0] <= n as f64 && pm.data[0].fract() == 0.0 => pm.data[0] as usize,
        _ => return rt(format!("p must be an integer between 1 and {n}.")),
    };
    let nf = n as f64;
    let sq: Mat = d.iter().map(|r| r.iter().map(|v| v * v).collect()).collect();
    let colsum: Vec<f64> = (0..n).map(|j| (0..n).map(|i| sq[i][j]).sum()).collect();
    let rowsum: Vec<f64> = sq.iter().map(|r| r.iter().sum()).collect();
    let tot: f64 = sq.iter().flatten().sum();
    let b: Mat = (0..n).map(|i| (0..n).map(|j| (((sq[i][j] - colsum[j] / nf) - rowsum[i] / nf) + tot / (nf * nf)) * -0.5).collect()).collect();
    let bs: Mat = (0..n).map(|i| (0..n).map(|j| (b[i][j] + b[j][i]) / 2.0).collect()).collect();
    let (ev, vecs) = eig_sym(it, &bs)?;
    let evp: Vec<f64> = ev.iter().take(p).cloned().collect();
    let emax = evp.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let thr = emax * (if single { f32::EPSILON as f64 } else { f64::EPSILON }).powf(0.75);
    let keep: Vec<usize> = (0..evp.len()).filter(|&k| evp[k] > thr).collect();
    let yv: Matrix;
    if keep.is_empty() {
        yv = Matrix::from_dims(vec![n, 1], vec![0.0; n]);
    } else {
        let mut cols: Vec<Vec<f64>> = keep.iter().map(|&k| (0..n).map(|i| vecs[i][k] * evp[k].sqrt()).collect()).collect();
        for c in cols.iter_mut() {
            let mut best = 0;
            for i in 0..n {
                if c[i].abs() > c[best].abs() {
                    best = i;
                }
            }
            let sg = c[best].signum();
            c.iter_mut().for_each(|v| *v *= sg);
        }
        yv = Matrix::from_dims(vec![n, cols.len()], cols.into_iter().flatten().collect());
    }
    let e = Matrix::from_dims(vec![ev.len(), 1], ev);
    let outs = vec![Value::Num(finish(yv, single)), Value::Num(finish(if p < n { Matrix::from_dims(vec![evp.len(), 1], evp) } else { e }, single))];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn corrcov(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(c)) = a.first() else { return rt("Not enough input arguments.") };
    if a.len() < 2 {
        check_cov(it, c, "corrcov")?;
    }
    let single = is_single(c);
    let (m, n) = (c.dims[0], c.dims[1]);
    let _ = n;
    let sigma: Vec<f64> = (0..m).map(|i| c.data[i * m + i].sqrt()).collect();
    let mut r = vec![0.0; m * m];
    for j in 0..m {
        for i in 0..m {
            let v = c.data[j * m + i] / sigma[i] / sigma[j];
            r[j * m + i] = if v.abs() > 1.0 { v / v.abs() } else { v };
        }
    }
    for i in 0..m {
        let d = r[i * m + i];
        r[i * m + i] = if d.is_nan() { f64::NAN } else { d.signum() };
    }
    let outs = vec![Value::Num(finish(Matrix::from_dims(vec![m, m], r), single)), Value::Num(finish(mk_col(sigma), single))];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "pca" => pca(it, a, nargout),
        "pcacov" => pcacov(it, a, nargout),
        "pcares" => pcares(it, a, nargout),
        "cmdscale" => cmdscale(it, a, nargout),
        _ => corrcov(it, a, nargout),
    }
}
