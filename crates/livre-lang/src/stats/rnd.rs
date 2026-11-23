//! Random number generators of the univariate families (`unifrnd`, `exprnd`, `binornd`, ...),
//! drawing from the global stream in the order MATLAB's implementations do.

use super::dist::FAMILIES;
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::nd::Dims;
use crate::random::{global_normal, global_uniform};
use crate::value::*;

pub const NAMES: &[&str] = &[
    "normrnd", "lognrnd", "unifrnd", "exprnd", "wblrnd", "raylrnd", "evrnd", "gevrnd", "gprnd", "geornd", "unidrnd", "binornd", "hygernd",
];

/// (minimum number of parameters, number of parameters)
fn arity(name: &str) -> (usize, usize) {
    match name {
        "normrnd" | "lognrnd" => (0, 2),
        "unifrnd" | "wblrnd" | "evrnd" | "binornd" => (2, 2),
        "gevrnd" | "gprnd" | "hygernd" => (3, 3),
        _ => (1, 1),
    }
}

/// The output size from the parameter sizes and the trailing size arguments.
fn out_dims(params: &[&Matrix], sizes: &[Value]) -> Res<Dims> {
    let mut pdims: Option<&Dims> = None;
    for m in params {
        if m.numel() != 1 {
            match pdims {
                None => pdims = Some(&m.dims),
                Some(d) if d != &m.dims => return rt("Size information is inconsistent."),
                _ => {}
            }
        }
    }
    let size_vals = |v: &Value| -> Res<Vec<f64>> {
        match v {
            Value::Num(m) => Ok(m.data.clone()),
            _ => rt("Size inputs must be integers."),
        }
    };
    let dims: Dims = if sizes.is_empty() {
        pdims.cloned().unwrap_or_else(|| vec![1, 1])
    } else if sizes.len() == 1 {
        let v = size_vals(&sizes[0])?;
        if v.iter().any(|x| x.fract() != 0.0 || !x.is_finite()) {
            return rt("Size inputs must be integers.");
        }
        let v: Vec<usize> = v.iter().map(|&x| x.max(0.0) as usize).collect();
        match v.len() {
            0 => return rt("Size inputs must be integers."),
            1 => vec![v[0], v[0]],
            _ => v,
        }
    } else {
        let mut d = Vec::new();
        for s in sizes {
            let v = size_vals(s)?;
            if v.len() != 1 {
                return rt("Size inputs must be scalar.");
            }
            if v[0].fract() != 0.0 || !v[0].is_finite() {
                return rt("Size inputs must be integers.");
            }
            d.push(v[0].max(0.0) as usize);
        }
        d
    };
    let dims = crate::nd::normalize(dims);
    if let Some(pd) = pdims {
        if *pd != dims {
            return rt("Size information is inconsistent.");
        }
    }
    Ok(dims)
}

fn p_at(m: &Matrix, i: usize) -> f64 {
    if m.numel() == 1 { m.data[0] } else { m.data[i] }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    let (minp, np) = arity(name);
    if a.len() < minp {
        return rt(match minp {
            1 => "Requires at least one input argument.",
            2 => "Requires at least two input arguments.",
            _ => "Requires at least three input arguments.",
        });
    }
    let nparams = a.len().min(np);
    let mut pm: Vec<Matrix> = Vec::new();
    for v in a.iter().take(nparams) {
        match v {
            Value::Num(m) => pm.push(m.clone()),
            other => return rt(format!("Undefined function '{name}' for input arguments of type '{}'.", other.class_name())),
        }
    }
    // Defaults of normrnd and lognrnd.
    while pm.len() < np && matches!(name, "normrnd" | "lognrnd") {
        pm.push(Matrix::scalar(if pm.is_empty() { 0.0 } else { 1.0 }));
    }
    let single = pm.iter().any(is_single);
    let prefs: Vec<&Matrix> = pm.iter().collect();
    let dims = out_dims(&prefs, &a[nparams..])?;
    let n: usize = dims.iter().product();
    let p = |k: usize, i: usize| p_at(&pm[k], i);
    let mut out = vec![0.0; n];
    match name {
        "normrnd" | "lognrnd" => {
            let z = global_normal(it, n);
            for i in 0..n {
                let (mu, mut sg) = (p(0, i), p(1, i));
                if sg < 0.0 {
                    sg = f64::NAN;
                }
                let v = z[i] * sg + mu;
                out[i] = if name == "lognrnd" { v.exp() } else { v };
            }
        }
        "raylrnd" => {
            let z1 = global_normal(it, n);
            let z2 = global_normal(it, n);
            for i in 0..n {
                let b = if p(0, i) < 0.0 { f64::NAN } else { p(0, i) };
                out[i] = (z1[i] * z1[i] + z2[i] * z2[i]).sqrt() * b;
            }
        }
        "binornd" => {
            let (nn, pp) = (&pm[0], &pm[1]);
            if nn.numel() == 1 && pp.numel() == 1 {
                let (nv, pv) = (nn.data[0], pp.data[0]);
                if (0.0..=1.0).contains(&pv) && nv >= 0.0 && nv.round() == nv {
                    let u = global_uniform(it, n * nv as usize);
                    for i in 0..n {
                        out[i] = (0..nv as usize).filter(|&t| u[i + n * t] < pv).count() as f64;
                    }
                } else {
                    out.iter_mut().for_each(|v| *v = f64::NAN);
                }
            } else if nn.numel() == 1 {
                let nv = nn.data[0];
                if nv >= 0.0 && nv.round() == nv {
                    let u = global_uniform(it, n * nv as usize);
                    for i in 0..n {
                        let pv = p(1, i);
                        out[i] = if !(0.0..=1.0).contains(&pv) { f64::NAN } else { (0..nv as usize).filter(|&t| u[i + n * t] < pv).count() as f64 };
                    }
                } else {
                    out.iter_mut().for_each(|v| *v = f64::NAN);
                }
            } else {
                let nmax = nn.data.iter().cloned().fold(0.0, f64::max) as usize;
                for t in 1..=nmax {
                    let ks: Vec<usize> = (0..n).filter(|&i| p(0, i) >= t as f64).collect();
                    let u = global_uniform(it, ks.len());
                    for (j, &k) in ks.iter().enumerate() {
                        if u[j] < p(1, k) {
                            out[k] += 1.0;
                        }
                    }
                }
                for i in 0..n {
                    let (nv, pv) = (p(0, i), p(1, i));
                    if pv < 0.0 || pv > 1.0 || nv < 0.0 || nv.round() != nv || nv.is_nan() || pv.is_nan() {
                        out[i] = f64::NAN;
                    }
                }
            }
        }
        "hygernd" => {
            let valid = |i: usize| {
                let (m, k, nn) = (p(0, i), p(1, i), p(2, i));
                nn <= m && k <= m && nn >= 0.0 && nn.round() == nn && k >= 0.0 && k.round() == k && m >= 0.0 && m.round() == m
            };
            let hf = FAMILIES.iter().find(|f| f.name == "hyge").unwrap();
            let ok: Vec<usize> = (0..n).filter(|&i| valid(i)).collect();
            let scalar_params = pm.iter().all(|m| m.numel() == 1);
            let all_ok = if scalar_params { ok.len() == 1 } else { ok.len() == n };
            out.iter_mut().for_each(|v| *v = f64::NAN);
            let idx: Vec<usize> = if all_ok { (0..n).collect() } else { ok };
            let u = global_uniform(it, idx.len());
            for (j, &i) in idx.iter().enumerate() {
                out[i] = (hf.inv)(u[j], &[p(0, i), p(1, i), p(2, i)]);
            }
        }
        _ => {
            let u = global_uniform(it, n);
            for i in 0..n {
                out[i] = match name {
                    "unifrnd" => {
                        let (a2, b2) = (p(0, i) / 2.0, p(1, i) / 2.0);
                        let (mu, sig) = (a2 + b2, b2 - a2);
                        if p(0, i) > p(1, i) { f64::NAN } else { mu + sig * (2.0 * u[i] - 1.0) }
                    }
                    "exprnd" => {
                        let mu = if p(0, i) < 0.0 { f64::NAN } else { p(0, i) };
                        -mu * u[i].ln()
                    }
                    "wblrnd" => {
                        let (aa, bb) = (if p(0, i) < 0.0 { f64::NAN } else { p(0, i) }, if p(1, i) < 0.0 { f64::NAN } else { p(1, i) });
                        aa * (-u[i].ln()).powf(1.0 / bb)
                    }
                    "evrnd" => {
                        let sg = if p(1, i) < 0.0 { f64::NAN } else { p(1, i) };
                        (-u[i].ln()).ln() * sg + p(0, i)
                    }
                    "gevrnd" => {
                        let (k, sg, mu) = (p(0, i), if p(1, i) < 0.0 { f64::NAN } else { p(1, i) }, p(2, i));
                        let r = if k.abs() < f64::EPSILON { -(-u[i].ln()).ln() } else { (-k * (-u[i].ln()).ln()).exp_m1() / k };
                        mu + sg * r
                    }
                    "gprnd" => {
                        let (k, sg, th) = (p(0, i), if p(1, i) < 0.0 { f64::NAN } else { p(1, i) }, p(2, i));
                        let r = if k.abs() < f64::EPSILON { -u[i].ln() } else { (-k * u[i].ln()).exp_m1() / k };
                        th + sg * r
                    }
                    "geornd" => {
                        let pr = p(0, i);
                        if pr <= 0.0 || pr > 1.0 {
                            f64::NAN
                        } else {
                            let r = ((u[i].ln() / (1.0 - pr).ln()).abs() - 1.0).ceil();
                            if r < 0.0 { 0.0 } else { r }
                        }
                    }
                    _ => {
                        let nn = p(0, i);
                        if nn <= 0.0 || nn.round() != nn { f64::NAN } else { (nn * u[i]).ceil() }
                    }
                };
            }
        }
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(dims, out), single))]))
}
