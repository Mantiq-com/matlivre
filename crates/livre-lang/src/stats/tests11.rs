//! `bootci` and `lhsdesign`.

use super::dist::{norminv_std, normcdf_std};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["bootci", "lhsdesign"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "bootci" => bootci(it, a, nargout),
        "lhsdesign" => lhsdesign(it, a),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn eval_fun(it: &mut Interpreter, f: &Value, args: Vec<Value>) -> Result<Vec<Value>, crate::error::LangError> {
    if f.is_text() {
        let name = f.text().unwrap_or_default();
        it.call_function(&name, args, 1)
    } else if let Value::Func(fv) = f {
        it.call_value(fv, args, 1)
    } else {
        rt("The bootstrap function must be a function handle or a function name.")
    }
}

fn num_of(v: &Value) -> Option<Matrix> {
    match v {
        Value::Num(m) => Some(m.clone()),
        _ => None,
    }
}

fn prctile1(it: &mut Interpreter, col: &[f64], pct: f64) -> Result<f64, crate::error::LangError> {
    let r = it.call_function("prctile", vec![Value::Num(Matrix::col(col.to_vec())), Value::num(pct)], 1)?;
    r[0].scalar()
}

fn cols_of(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1..].iter().product::<usize>());
    (0..c).map(|j| m.data[j * r..(j + 1) * r].to_vec()).collect()
}

fn std_of(c: &[f64]) -> f64 {
    let n = c.len() as f64;
    let m = c.iter().sum::<f64>() / n;
    (c.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
}

fn bootci(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("BOOTCI requires at least two arguments.");
    }
    let nboot = match &a[0] {
        Value::Num(m) if m.data.len() == 1 && m.kind != Kind::Char => m.data[0],
        _ => f64::NAN,
    };
    if !(nboot > 0.0) || nboot != nboot.round() {
        return rt("NBOOT must be a positive integer.");
    }
    let (fun, data, opts): (Value, Vec<Value>, Vec<Value>);
    let mut ty = "bca".to_string();
    let mut alpha = 0.05;
    let mut stderrfun: Option<Value> = None;
    let mut nbootstd = 100.0;
    let mut weights = Value::empty();
    match &a[1] {
        Value::Cell(c) => {
            if c.data.is_empty() {
                return rt("The first cell of BOOTFUN must hold the bootstrap function.");
            }
            fun = c.data[0].clone();
            data = c.data[1..].to_vec();
            opts = a[2..].to_vec();
            if opts.len() % 2 != 0 {
                return rt("Wrong number of arguments.");
            }
            let names = ["type", "alpha", "stderr", "nbootstd", "weights", "options"];
            for pair in opts.chunks(2) {
                let key = match pair[0].text() {
                    Ok(k) if pair[0].is_text() => k,
                    _ => return rt("Parameter name must be text."),
                };
                let kl = key.to_lowercase();
                let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
                let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
                match pick {
                    Some(0) => ty = pair[1].text().unwrap_or_default(),
                    Some(1) => alpha = match &pair[1] {
                        Value::Num(m) if m.data.len() == 1 => m.data[0],
                        _ => f64::NAN,
                    },
                    Some(2) => stderrfun = Some(pair[1].clone()).filter(|v| !is_empty_val(v)),
                    Some(3) => nbootstd = match &pair[1] {
                        Value::Num(m) if m.data.len() == 1 => m.data[0],
                        _ => f64::NAN,
                    },
                    Some(4) => weights = pair[1].clone(),
                    Some(5) => {}
                    None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
                    _ => return rt(format!("Invalid parameter name: {key}.")),
                }
            }
        }
        f => {
            fun = f.clone();
            data = a[2..].to_vec();
        }
    }
    // the statistic on the observed data
    let f2s = |it: &mut Interpreter, f: &Value| -> String {
        match it.call_function("func2str", vec![f.clone()], 1) {
            Ok(v) => v[0].text().unwrap_or_default(),
            Err(_) => f.text().unwrap_or_default(),
        }
    };
    let obs = match eval_fun(it, &fun, data.clone()) {
        Ok(v) => v,
        Err(_) => {
            let nm = f2s(it, &fun);
            return rt(format!("Error evaluating BOOTFUN function '{nm}'."));
        }
    };
    let Some(Value::Num(om)) = obs.first() else { return rt("The bootstrap function must return numeric values.") };
    if om.data.iter().any(|v| !v.is_finite()) {
        return rt("BOOTFUN returns a NaN or Inf.");
    }
    let sz: Option<Vec<usize>> = if om.dims.len() == 2 && (om.dims[0] == 1 || om.dims[1] == 1) { None } else { Some(om.dims.clone()) };
    let stat: Vec<f64> = om.data.clone();
    let m = stat.len();
    let boot_args = |it: &mut Interpreter, nb: f64, data: &[Value], w: &Value| -> Result<Vec<Value>, crate::error::LangError> {
        let mut args = vec![Value::num(nb), fun.clone()];
        args.extend(data.iter().cloned());
        args.push(Value::chars("weights"));
        args.push(w.clone());
        args.push(Value::chars("Options"));
        args.push(Value::empty());
        it.call_function("bootstrp", args, 2)
    };
    let tl = ty.to_lowercase();
    let lower_type = tl.as_str();
    let za = norminv_std(alpha / 2.0);
    let mut lo = vec![0.0; m];
    let mut hi = vec![0.0; m];
    let bstat_val: Value;
    match lower_type {
        "norm" | "normal" => {
            let r = boot_args(it, nboot, &data, &weights)?;
            bstat_val = r[0].clone();
            let bm = num_of(&bstat_val).unwrap();
            for (i, c) in cols_of(&bm).iter().enumerate() {
                let se = std_of(c);
                let bias = c.iter().map(|v| v - stat[i]).sum::<f64>() / c.len() as f64;
                lo[i] = stat[i] - bias + se * za;
                hi[i] = stat[i] - bias - se * za;
            }
        }
        "per" | "percentile" => {
            let r = boot_args(it, nboot, &data, &weights)?;
            bstat_val = r[0].clone();
            let bm = num_of(&bstat_val).unwrap();
            let (p1, p2) = (100.0 * alpha / 2.0, 100.0 - 100.0 * alpha / 2.0);
            for (i, c) in cols_of(&bm).iter().enumerate() {
                lo[i] = prctile1(it, c, p1)?;
                hi[i] = prctile1(it, c, p2)?;
            }
        }
        "cper" | "corrected percentile" | "bca" => {
            let r = boot_args(it, nboot, &data, &weights)?;
            bstat_val = r[0].clone();
            let bm = num_of(&bstat_val).unwrap();
            let bcols = cols_of(&bm);
            let z0: Vec<f64> = (0..m)
                .map(|i| {
                    let n = bcols[i].len() as f64;
                    let less = bcols[i].iter().filter(|&&v| v < stat[i]).count() as f64 / n;
                    let eq = bcols[i].iter().filter(|&&v| v == stat[i]).count() as f64 / n;
                    norminv_std(less + eq / 2.0)
                })
                .collect();
            if lower_type != "bca" {
                for i in 0..m {
                    let pct1 = 100.0 * normcdf_std(2.0 * z0[i] - za);
                    let pct2 = 100.0 * normcdf_std(2.0 * z0[i] + za);
                    lo[i] = prctile1(it, &bcols[i], pct2)?;
                    hi[i] = prctile1(it, &bcols[i], pct1)?;
                }
            } else {
                let mut jargs = vec![fun.clone()];
                jargs.extend(data.iter().cloned());
                let jr = it.call_function("jackknife", jargs, 1)?;
                let Some(Value::Num(jm)) = jr.first() else { return rt("jackknife returned no value") };
                let jcols = cols_of(jm);
                let nn = jm.dims[0] as f64;
                let w: Vec<f64> = match &weights {
                    Value::Num(wm) if !wm.is_empty() => {
                        let s: f64 = wm.data.iter().sum();
                        wm.data.iter().map(|v| v / s).collect()
                    }
                    _ => vec![1.0 / nn; jm.dims[0]],
                };
                for i in 0..m {
                    let mj: f64 = jcols[i].iter().zip(&w).map(|(v, wt)| v * wt).sum();
                    let score: Vec<f64> = jcols[i].iter().map(|v| mj - v).collect();
                    let iszer = score.iter().all(|&v| v == 0.0);
                    let s3: f64 = score.iter().zip(&w).map(|(v, wt)| v.powi(3) * wt).sum();
                    let s2: f64 = score.iter().zip(&w).map(|(v, wt)| v * v * wt).sum();
                    let mut skew = s3 / s2.powf(1.5) / nn.sqrt();
                    if iszer {
                        skew = 0.0;
                    }
                    let acc = skew / 6.0;
                    let z1 = za;
                    let z2 = -za;
                    let mut pct1 = 100.0 * normcdf_std(z0[i] + (z0[i] + z1) / (1.0 - acc * (z0[i] + z1)));
                    let mut pct2 = 100.0 * normcdf_std(z0[i] + (z0[i] + z2) / (1.0 - acc * (z0[i] + z2)));
                    if z0[i] == f64::INFINITY {
                        pct1 = 100.0;
                        pct2 = 100.0;
                    }
                    if z0[i] == f64::NEG_INFINITY {
                        pct1 = 0.0;
                        pct2 = 0.0;
                    }
                    let l = prctile1(it, &bcols[i], pct2)?;
                    let u = prctile1(it, &bcols[i], pct1)?;
                    lo[i] = l.min(u);
                    hi[i] = l.max(u);
                }
            }
        }
        "stud" | "student" => {
            if stderrfun.is_none() && (!(nbootstd > 0.0) || nbootstd != nbootstd.round()) {
                return rt("NBOOTSTD must be a positive integer.");
            }
            let r = boot_args(it, nboot, &data, &weights)?;
            bstat_val = r[0].clone();
            let bootsam = num_of(&r[1]).unwrap();
            let bm = num_of(&bstat_val).unwrap();
            let bcols = cols_of(&bm);
            let nb = nboot as usize;
            let nn = bootsam.dims[0];
            // the data as columns, scalars kept
            let mut dvals: Vec<Value> = data.clone();
            let mut scalard = vec![false; dvals.len()];
            for (k, d) in dvals.iter_mut().enumerate() {
                if let Value::Num(mm) = d {
                    let (row, col) = (mm.dims[0], mm.dims[1..].iter().product::<usize>());
                    if row.max(col) == 1 {
                        scalard[k] = true;
                    }
                    if row == 1 && col != 1 {
                        let mut c = mm.clone();
                        c.dims = vec![col, 1];
                        *d = Value::Num(c);
                    }
                }
            }
            if let Some(sf) = &stderrfun {
                match eval_fun(it, sf, data.clone()) {
                    Ok(v) => {
                        if let Some(Value::Num(sm)) = v.first() {
                            if sm.data.iter().any(|x| !x.is_finite()) {
                                return rt("STDERR returns a NaN or Inf.");
                            }
                        }
                    }
                    Err(_) => {
                        let nm = f2s(it, sf);
                        return rt(format!("Error evaluating STDERR function '{nm}'."));
                    }
                }
            }
            let mut sd_t = vec![vec![0.0; nb]; m];
            for b in 0..nb {
                let idx: Vec<usize> = (0..nn).map(|i| bootsam.data[b * nn + i] as usize).collect();
                let db: Vec<Value> = dvals
                    .iter()
                    .enumerate()
                    .map(|(k, d)| match d {
                        Value::Num(mm) if !scalard[k] => {
                            let (rr, cc) = (mm.dims[0], mm.dims[1..].iter().product::<usize>());
                            let mut out = Vec::with_capacity(idx.len() * cc);
                            for j in 0..cc {
                                for &i in &idx {
                                    out.push(mm.data[j * rr + i - 1]);
                                }
                            }
                            let mut mo = Matrix::from_dims(vec![idx.len(), cc], out);
                            mo.kind = mm.kind;
                            Value::Num(mo)
                        }
                        other => other.clone(),
                    })
                    .collect();
                if let Some(sf) = &stderrfun {
                    let rv = eval_fun(it, sf, db)?;
                    let Some(Value::Num(sm)) = rv.first() else { return rt("The standard error function must return numeric values.") };
                    if sm.data.len() != m {
                        return rt("Number of elements must not change. Use [] as one of the size inputs to automatically calculate the appropriate size for that dimension.");
                    }
                    for (i, v) in sm.data.iter().enumerate().take(m) {
                        sd_t[i][b] = *v;
                    }
                } else {
                    let w = Value::empty();
                    let r2 = boot_args(it, nbootstd, &db, &w)?;
                    let sm = num_of(&r2[0]).unwrap();
                    for (i, c) in cols_of(&sm).iter().enumerate().take(m) {
                        sd_t[i][b] = std_of(c);
                    }
                }
            }
            let (p1, p2) = (100.0 * alpha / 2.0, 100.0 * (1.0 - alpha / 2.0));
            for i in 0..m {
                let tstat: Vec<f64> = (0..nb).map(|b| if sd_t[i][b] > 0.0 { (bcols[i][b] - stat[i]) / sd_t[i][b] } else { 0.0 }).collect();
                let l = prctile1(it, &tstat, p1)?;
                let u = prctile1(it, &tstat, p2)?;
                let sb = std_of(&bcols[i]);
                lo[i] = l * sb + stat[i];
                hi[i] = u * sb + stat[i];
            }
        }
        _ => return rt("Bad confidence interval type"),
    }
    let mut ci_data = Vec::with_capacity(2 * m);
    for i in 0..m {
        ci_data.push(lo[i]);
        ci_data.push(hi[i]);
    }
    let ci = match sz {
        Some(s) => {
            let mut dims = vec![2];
            dims.extend(s);
            Matrix::from_dims(dims, ci_data)
        }
        None => Matrix::from_dims(vec![2, m], ci_data),
    };
    let mut out = vec![Value::Num(ci)];
    if nargout > 1 {
        out.push(bstat_val);
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// lhsdesign

fn rank_of(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| x[i].partial_cmp(&x[j]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; n];
    for (k, &i) in idx.iter().enumerate() {
        r[i] = (k + 1) as f64;
    }
    r
}

fn rand_mat(it: &mut Interpreter, n: usize, p: usize) -> Result<Vec<f64>, crate::error::LangError> {
    let r = it.call_function("rand", vec![Value::num(n as f64), Value::num(p as f64)], 1)?;
    match r.first() {
        Some(Value::Num(m)) => Ok(m.data.clone()),
        _ => rt("rand returned no value"),
    }
}

fn getsample(it: &mut Interpreter, n: usize, p: usize, smooth: bool) -> Result<Vec<f64>, crate::error::LangError> {
    let mut x = rand_mat(it, n, p)?;
    for j in 0..p {
        let r = rank_of(&x[j * n..(j + 1) * n]);
        x[j * n..(j + 1) * n].copy_from_slice(&r);
    }
    if smooth {
        let u = rand_mat(it, n, p)?;
        for i in 0..n * p {
            x[i] -= u[i];
        }
    } else {
        x.iter_mut().for_each(|v| *v -= 0.5);
    }
    x.iter_mut().for_each(|v| *v /= n as f64);
    Ok(x)
}

fn score(x: &[f64], n: usize, p: usize, crit: &str) -> Result<f64, crate::error::LangError> {
    if n < 2 {
        return Ok(0.0);
    }
    if p == 0 && crit == "maximin" {
        return rt("X and Y can't be empty.");
    }
    if crit == "correlation" {
        let cols: Vec<&[f64]> = (0..p).map(|j| &x[j * n..(j + 1) * n]).collect();
        let mean = |c: &[f64]| c.iter().sum::<f64>() / n as f64;
        let mut s = 0.0;
        for j in 1..p {
            for k in 0..j {
                let (mj, mk) = (mean(cols[j]), mean(cols[k]));
                let cov: f64 = (0..n).map(|i| (cols[j][i] - mj) * (cols[k][i] - mk)).sum();
                let vj: f64 = (0..n).map(|i| (cols[j][i] - mj).powi(2)).sum();
                let vk: f64 = (0..n).map(|i| (cols[k][i] - mk).powi(2)).sum();
                let c = cov / (vj.sqrt() * vk.sqrt());
                s -= c * c;
            }
        }
        Ok(s)
    } else {
        let mut best = f64::INFINITY;
        for i in 0..n {
            let mut d = f64::INFINITY;
            for j in 0..n {
                if i != j {
                    let dd: f64 = (0..p).map(|c| (x[c * n + i] - x[c * n + j]).powi(2)).sum();
                    d = d.min(dd);
                }
            }
            best = best.min(d.sqrt());
        }
        Ok(best)
    }
}

fn lhsdesign(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let dim = |v: &Value| -> Result<usize, crate::error::LangError> {
        match v {
            Value::Num(m) if m.data.len() == 1 && m.data[0] == m.data[0].round() => Ok(m.data[0].max(0.0) as usize),
            Value::Num(m) if m.data.len() == 1 => rt("Size inputs must be integers."),
            _ => rt("Size inputs must be scalar."),
        }
    };
    let (n, p) = (dim(&a[0])?, dim(&a[1])?);
    let rest = &a[2..];
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let names = ["iterations", "criterion", "smooth"];
    let mut maxiter = f64::NAN;
    let mut crit_v: Option<Value> = None;
    let mut smooth_v: Option<Value> = None;
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match pick {
            Some(0) => {
                maxiter = match &pair[1] {
                    v if is_empty_val(v) => f64::NAN,
                    Value::Num(m) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] >= 0.0 => m.data[0],
                    _ => return rt("Value of 'iterations' parameter must be a scalar >= 0."),
                }
            }
            Some(1) => crit_v = Some(pair[1].clone()),
            Some(2) => smooth_v = Some(pair[1].clone()),
            None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let maxiter0 = if maxiter.is_nan() { 5 } else { maxiter as usize };
    let crit = match &crit_v {
        None => "maximin".to_string(),
        Some(v) => {
            let choices = ["none", "maximin", "correlation"];
            let i = super::tests3::param_val(v, &choices, "Criterion", true)?;
            choices[i].to_string()
        }
    };
    let smooth = match &smooth_v {
        None => true,
        Some(v) if is_empty_val(v) => true,
        Some(v) if v.is_text() && v.text().unwrap_or_default() == "on" => true,
        Some(v) if v.is_text() && v.text().unwrap_or_default() == "off" => false,
        _ => return rt("Value of 'smooth' parameter must be 'on' or 'off'."),
    };
    let mut x = getsample(it, n, p, smooth)?;
    let maxiter = if crit == "none" || n < 2 { 0 } else { maxiter0 };
    if crit == "maximin" {
        let mut best = score(&x, n, p, &crit)?;
        for _ in 2..=maxiter {
            let xn = getsample(it, n, p, smooth)?;
            let ns = score(&xn, n, p, &crit)?;
            if ns > best {
                x = xn;
                best = ns;
            }
        }
    } else if crit == "correlation" {
        let mut best = score(&x, n, p, &crit)?;
        let takeout = |xv: &[f64], yv: &[f64]| -> Vec<f64> {
            let mx = xv.iter().sum::<f64>() / n as f64;
            let my = yv.iter().sum::<f64>() / n as f64;
            let xc: Vec<f64> = xv.iter().map(|v| v - mx).collect();
            let mxc = xc.iter().sum::<f64>() / n as f64;
            let xcc: Vec<f64> = xc.iter().map(|v| v - mxc).collect();
            let _ = my;
            let ycc: Vec<f64> = {
                let yc: Vec<f64> = yv.iter().map(|v| v - my).collect();
                let m = yc.iter().sum::<f64>() / n as f64;
                yc.iter().map(|v| v - m).collect()
            };
            let b = xcc.iter().zip(&ycc).map(|(u, v)| u * v).sum::<f64>() / xcc.iter().map(|u| u * u).sum::<f64>();
            yv.iter().zip(&xc).map(|(y, xcv)| y - b * xcv).collect()
        };
        for _ in 2..=maxiter {
            for j in 1..p {
                for k in 0..j {
                    let z = takeout(&x[j * n..(j + 1) * n].to_vec(), &x[k * n..(k + 1) * n].to_vec());
                    let r = rank_of(&z);
                    for i in 0..n {
                        x[k * n + i] = (r[i] - 0.5) / n as f64;
                    }
                }
            }
            for j in (0..p.saturating_sub(1)).rev() {
                for k in ((j + 1)..p).rev() {
                    let z = takeout(&x[j * n..(j + 1) * n].to_vec(), &x[k * n..(k + 1) * n].to_vec());
                    let r = rank_of(&z);
                    for i in 0..n {
                        x[k * n + i] = (r[i] - 0.5) / n as f64;
                    }
                }
            }
            let ns = score(&x, n, p, &crit)?;
            if ns <= best {
                break;
            }
            best = ns;
        }
    }
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![n, p], x))]))
}
