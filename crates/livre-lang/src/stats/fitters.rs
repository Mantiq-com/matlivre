//! Closed-form distribution fits and negative log-likelihoods: `normfit`, `expfit`, `poissfit`,
//! `binofit`, `lognfit`, `unifit`, `raylfit`, `normlike`, `lognlike`, `explike`.

use super::dist::{chi2inv, norminv_std, normcdf_std, tinv};
use super::util::*;
use crate::builtins::Ret;
use crate::builtins_specfun::{betaincinv, erfc};
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["normfit", "expfit", "poissfit", "binofit", "lognfit", "unifit", "raylfit", "normlike", "lognlike", "explike"];

fn alpha_arg(v: Option<&Value>) -> f64 {
    match v {
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        _ => 0.05,
    }
}

fn is_vec(m: &Matrix) -> bool {
    m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1) && m.numel() >= 1
}

fn arg_matrix<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(_) => rt("Input must be numeric."),
        None => rt("Not enough input arguments."),
    }
}

/// Columns of `m` (a vector is one column).
fn columns(m: &Matrix) -> Vec<Vec<f64>> {
    if is_vec(m) {
        return vec![m.data.clone()];
    }
    let r = m.dims[0];
    let c: usize = m.dims[1..].iter().product();
    (0..c).map(|j| m.data[j * r..(j + 1) * r].to_vec()).collect()
}

fn mat(rows: usize, cols: usize, data: Vec<f64>, single: bool) -> Value {
    Value::Num(finish(Matrix::from_dims(vec![rows, cols], data), single))
}

fn normfit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    let single = is_single(x);
    let vector = is_vec(x);
    let alpha = alpha_arg(a.get(1));
    let has_cens = matches!(a.get(2), Some(Value::Num(m)) if !m.is_empty());
    if !vector && a.len() >= 3 {
        return rt("X must be a vector.");
    }
    if has_cens {
        let Some(Value::Num(c)) = a.get(2) else { unreachable!() };
        if c.dims != x.dims {
            return rt("X and CENSORING must have the same size.");
        }
        if c.data.iter().any(|&v| v != 0.0) {
            return rt("Censored data are not supported by normfit in this release.");
        }
    }
    let mut cols = columns(x);
    let mut freq: Option<Vec<f64>> = None;
    if let Some(Value::Num(f)) = a.get(3) {
        if !f.is_empty() {
            if f.dims != x.dims {
                return rt("X and FREQ must have the same size.");
            }
            freq = Some(f.data.clone());
        }
    }
    let ncols = cols.len();
    let mut n = cols.first().map_or(0, |c| c.len()) as f64;
    if let Some(f) = &freq {
        n = f.iter().sum();
        let keep: Vec<usize> = (0..f.len()).filter(|&i| f[i] != 0.0).collect();
        cols = cols.iter().map(|c| keep.iter().map(|&i| c[i]).collect()).collect();
        freq = Some(keep.iter().map(|&i| f[i]).collect());
    }
    let fq = |i: usize| freq.as_ref().map_or(1.0, |f| f[i]);
    if n == 0.0 {
        let nan = vec![f64::NAN; ncols];
        let mut outs = vec![mat(1, ncols, nan.clone(), single), mat(1, ncols, nan, single)];
        outs.push(mat(2, ncols, vec![f64::NAN; 2 * ncols], single));
        outs.push(mat(2, ncols, vec![f64::NAN; 2 * ncols], single));
        return Ok(Some(outs.into_iter().take(nargout.max(1)).collect()));
    }
    let (mut mu, mut sg, mut muci, mut sgci) = (vec![], vec![], vec![], vec![]);
    for c in &cols {
        let sumx: f64 = c.iter().enumerate().map(|(i, v)| fq(i) * v).sum();
        let m = sumx / n;
        let s = if n > 1.0 { (c.iter().enumerate().map(|(i, v)| (v - m) * (v - m) * fq(i)).sum::<f64>() / (n - 1.0)).sqrt() } else { 0.0 };
        mu.push(m);
        sg.push(s);
        if n > 1.0 {
            let t = tinv(1.0 - alpha / 2.0, n - 1.0) * s / n.sqrt();
            muci.extend([m - t, m + t]);
            let (l, u) = (chi2inv(1.0 - alpha / 2.0, n - 1.0), chi2inv(alpha / 2.0, n - 1.0));
            sgci.extend([s * ((n - 1.0) / l).sqrt(), s * ((n - 1.0) / u).sqrt()]);
        } else {
            muci.extend([f64::NEG_INFINITY, f64::INFINITY]);
            sgci.extend([0.0, f64::INFINITY]);
        }
    }
    let outs = vec![mat(1, ncols, mu, single), mat(1, ncols, sg, single), mat(2, ncols, muci, single), mat(2, ncols, sgci, single)];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn expfit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    let single = is_single(x);
    let vector = is_vec(x);
    if !vector && a.len() >= 3 {
        return rt("X must be a vector.");
    }
    if x.data.iter().any(|&v| v < 0.0) {
        return rt("The data in X must be non-negative.");
    }
    let alpha = alpha_arg(a.get(1));
    let cens: Option<Vec<f64>> = match a.get(2) {
        Some(Value::Num(m)) if !m.is_empty() => {
            if m.dims != x.dims {
                return rt("X and CENSORING must have the same size.");
            }
            Some(m.data.iter().map(|&v| if v != 0.0 { 1.0 } else { 0.0 }).collect())
        }
        _ => None,
    };
    let freq: Option<Vec<f64>> = match a.get(3) {
        Some(Value::Num(m)) if !m.is_empty() => {
            if m.dims != x.dims {
                return rt("X and FREQ must have the same size.");
            }
            Some(m.data.clone())
        }
        _ => None,
    };
    if x.is_empty() {
        let outs = vec![Value::Num(Matrix::scalar(f64::NAN)), Value::Num(Matrix::from_dims(vec![2, 1], vec![f64::NAN, f64::NAN]))];
        return Ok(Some(outs.into_iter().take(nargout.max(1)).collect()));
    }
    let cols = columns(x);
    let r = if vector { x.numel() } else { x.dims[0] };
    let (mut hat, mut ci) = (vec![], vec![]);
    for (j, c) in cols.iter().enumerate() {
        let off = j * r;
        let f = |i: usize| freq.as_ref().map_or(1.0, |f| f[off + i]);
        let cn = |i: usize| cens.as_ref().map_or(0.0, |c| c[off + i]);
        let n: f64 = if freq.is_some() { (0..c.len()).map(f).sum() } else { c.len() as f64 };
        let sumx: f64 = c.iter().enumerate().map(|(i, v)| if v.is_finite() { f(i) * v } else { f64::NAN }).sum();
        let nunc = n - (0..c.len()).map(|i| f(i) * cn(i)).sum::<f64>();
        if nunc > 0.0 {
            let mu = sumx / nunc;
            hat.push(mu);
            ci.extend([2.0 * nunc * mu / chi2inv(1.0 - alpha / 2.0, 2.0 * nunc), 2.0 * nunc * mu / chi2inv(alpha / 2.0, 2.0 * nunc)]);
        } else {
            hat.push(f64::NAN);
            ci.extend([f64::NAN, f64::NAN]);
        }
    }
    let nc = hat.len();
    let outs = vec![mat(1, nc, hat, single), mat(2, nc, ci, single)];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn poissfit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    let single = is_single(x);
    let alpha = alpha_arg(a.get(1));
    let cols = columns(x);
    if x.data.iter().any(|&v| v < 0.0) {
        return rt("X must not be negative.");
    }
    let m = cols.first().map_or(0, |c| c.len()) as f64;
    let (mut hat, mut ci) = (vec![], vec![]);
    for c in &cols {
        let sum: f64 = c.iter().sum();
        hat.push(sum / m);
        let (lo, hi) = if sum < 100.0 {
            ((if sum == 0.0 { 0.0 } else { chi2inv(alpha / 2.0, 2.0 * sum) / 2.0 }), chi2inv(1.0 - alpha / 2.0, 2.0 * (sum + 1.0)) / 2.0)
        } else {
            (sum + sum.sqrt() * norminv_std(alpha / 2.0), sum + sum.sqrt() * norminv_std(1.0 - alpha / 2.0))
        };
        ci.extend([lo / m, hi / m]);
    }
    let nc = hat.len();
    let outs = vec![mat(1, nc, hat, single), mat(2, nc, ci, single)];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn binofit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    let n = arg_matrix(a, 1).map_err(|_| crate::error::LangError::new("Not enough input arguments.", 0))?;
    if !is_vec(x) {
        return rt("First argument must be a vector.");
    }
    if n.numel() != 1 && n.dims != x.dims {
        return rt("The first two inputs must match in size.");
    }
    let single = is_single(x) || is_single(n);
    let alpha = alpha_arg(a.get(2));
    let nn = |i: usize| if n.numel() == 1 { n.data[0] } else { n.data[i] };
    for i in 0..x.numel() {
        let v = nn(i);
        if v < 0.0 || v != v.round() || v.is_infinite() || x.data[i] > v {
            return rt("All N values must be non-negative integers at least as large as X.");
        }
    }
    if x.data.iter().any(|&v| v < 0.0) {
        return rt("X must not be negative.");
    }
    let len = x.numel();
    let phat: Vec<f64> = (0..len).map(|i| x.data[i] / nn(i)).collect();
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(x.dims.clone(), phat), single))];
    if nargout > 1 {
        let mut ci = vec![0.0; 2 * len];
        for i in 0..len {
            let (xi, ni) = (x.data[i], nn(i));
            ci[i] = if xi == 0.0 { 0.0 } else { betaincinv(alpha / 2.0, xi, ni - xi + 1.0, false) };
            ci[len + i] = if xi == ni { 1.0 } else { betaincinv(1.0 - alpha / 2.0, xi + 1.0, ni - xi, false) };
        }
        outs.push(mat(len, 2, ci, single));
    }
    Ok(Some(outs))
}

fn lognfit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    if !is_vec(x) {
        return rt("X must be a vector.");
    }
    if x.data.iter().any(|&v| v <= 0.0) {
        return rt("The data in X must be positive.");
    }
    let mut lx = x.clone();
    for v in lx.data.iter_mut() {
        *v = v.ln();
    }
    let mut args = vec![Value::Num(lx)];
    args.extend(a[1.min(a.len())..].iter().cloned());
    let r = normfit(&args, 4)?.unwrap();
    let get = |i: usize| match &r[i] {
        Value::Num(m) => m.clone(),
        _ => Matrix::empty(),
    };
    let (mu, sg, muci, sgci) = (get(0), get(1), get(2), get(3));
    let single = is_single(&mu);
    let mut outs = vec![mat(1, 2, vec![mu.data[0], sg.data[0]], single)];
    if nargout > 1 {
        outs.push(mat(2, 2, vec![muci.data[0], muci.data[1], sgci.data[0], sgci.data[1]], single));
    }
    Ok(Some(outs))
}

fn unifit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    let single = is_single(x);
    let alpha = alpha_arg(a.get(1));
    let cols = columns(x);
    let mut lo = vec![];
    let mut hi = vec![];
    let ssz = x.numel() as f64 / cols.len() as f64;
    for c in &cols {
        lo.push(c.iter().cloned().fold(f64::INFINITY, f64::min));
        hi.push(c.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    }
    let nc = cols.len();
    let mut aci = vec![];
    let mut bci = vec![];
    for j in 0..nc {
        let tmp = (hi[j] - lo[j]) / alpha.powf(1.0 / ssz);
        aci.extend([hi[j] - tmp, lo[j]]);
        bci.extend([hi[j], lo[j] + tmp]);
    }
    let outs = vec![mat(1, nc, lo, single), mat(1, nc, hi, single), mat(2, nc, aci, single), mat(2, nc, bci, single)];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

fn raylfit(a: &[Value], nargout: usize) -> Ret {
    let x = arg_matrix(a, 0)?;
    if x.dims.len() > 2 || (x.dims[0] > 1 && x.dims[1] > 1) {
        return rt("DATA argument must be a vector.");
    }
    let single = is_single(x);
    let alpha = alpha_arg(a.get(1));
    let m = x.numel() as f64;
    let phat = (0.5 * x.data.iter().map(|v| v * v).sum::<f64>() / m).sqrt();
    let mut outs = vec![mat(1, 1, vec![phat], single)];
    if nargout > 1 {
        let ci = vec![(2.0 * m * phat * phat / chi2inv(1.0 - alpha / 2.0, 2.0 * m)).sqrt(), (2.0 * m * phat * phat / chi2inv(alpha / 2.0, 2.0 * m)).sqrt()];
        outs.push(mat(2, 1, ci, single));
    }
    Ok(Some(outs))
}

struct Like {
    data: Vec<f64>,
    cens: Vec<f64>,
    freq: Vec<f64>,
}

fn like_inputs(a: &[Value], who: &str, nparams: usize) -> Res<(Vec<f64>, Like, bool)> {
    if a.len() < 2 {
        return rt("Requires two input arguments");
    }
    let (Value::Num(p), Value::Num(d)) = (&a[0], &a[1]) else { return rt("Inputs must be numeric.") };
    if d.numel() > d.dims[0].max(d.dims[1]) || d.dims.len() > 2 {
        return rt(if who == "normlike" { "DATA must be a vector." } else { "DATA must be a vector." });
    }
    if p.numel() != nparams {
        return rt(if nparams == 2 { "First input must be a vector with 2 elements." } else { "First input must be a scalar." });
    }
    let n = d.numel();
    let mut cens = vec![0.0; n];
    if let Some(Value::Num(c)) = a.get(2) {
        if !c.is_empty() {
            if c.dims != d.dims {
                return rt("DATA and CENSORING must have the same size.");
            }
            cens = c.data.iter().map(|&v| if v != 0.0 && who != "normlike" { 1.0 } else if who == "normlike" { v } else { 0.0 }).collect();
        }
    }
    let mut data = d.data.clone();
    let mut freq = vec![1.0; n];
    if let Some(Value::Num(f)) = a.get(3) {
        if !f.is_empty() {
            if f.dims != d.dims {
                return rt("DATA and FREQ must have the same size.");
            }
            let keep: Vec<usize> = (0..n).filter(|&i| f.data[i] != 0.0).collect();
            freq = keep.iter().map(|&i| f.data[i]).collect();
            data = keep.iter().map(|&i| d.data[i]).collect();
            cens = keep.iter().map(|&i| cens[i]).collect();
        }
    }
    Ok((p.data.clone(), Like { data, cens, freq }, is_single(p) || is_single(d)))
}

fn normlike_core(params: &[f64], l: &Like, want_avar: bool) -> (f64, Vec<f64>) {
    let (mu, mut sg) = (params[0], params[1]);
    if sg <= 0.0 {
        sg = f64::NAN;
    }
    let n = l.data.len();
    let z: Vec<f64> = l.data.iter().map(|d| (d - mu) / sg).collect();
    let mut ll: Vec<f64> = z.iter().map(|z| -0.5 * z * z - ((2.0 * std::f64::consts::PI).sqrt() * sg).ln()).collect();
    let ncen: f64 = (0..n).map(|i| l.freq[i] * l.cens[i]).sum();
    let mut scen = vec![0.0; n];
    if ncen > 0.0 {
        for i in 0..n {
            if l.cens[i] == 1.0 {
                scen[i] = 0.5 * erfc(z[i] / std::f64::consts::SQRT_2);
                ll[i] = scen[i].ln();
            }
        }
    }
    let nlogl = -(0..n).map(|i| l.freq[i] * ll[i]).sum::<f64>();
    let mut avar = vec![];
    if want_avar {
        let (mut d11, mut d12, mut d22) = (vec![-1.0; n], vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            d12[i] = -2.0 * z[i];
            d22[i] = 1.0 - 3.0 * z[i] * z[i];
            if ncen > 0.0 && l.cens[i] == 1.0 {
                let zc = z[i];
                let dlog = (-0.5 * zc * zc).exp() / ((2.0 * std::f64::consts::PI).sqrt() * scen[i]);
                let d2log = dlog * (dlog - zc);
                d11[i] = -d2log;
                d12[i] = -dlog - zc * d2log;
                d22[i] = -zc * (2.0 * dlog + zc * d2log);
            }
        }
        let nh11: f64 = -(0..n).map(|i| l.freq[i] * d11[i]).sum::<f64>();
        let nh12: f64 = -(0..n).map(|i| l.freq[i] * d12[i]).sum::<f64>();
        let nh22: f64 = -(0..n).map(|i| l.freq[i] * d22[i]).sum::<f64>();
        let det = nh11 * nh22 - nh12 * nh12;
        let s2 = sg * sg;
        avar = vec![s2 * nh22 / det, s2 * -nh12 / det, s2 * -nh12 / det, s2 * nh11 / det];
    }
    (nlogl, avar)
}

fn like(name: &str, a: &[Value], nargout: usize) -> Ret {
    let np = if name == "explike" { 1 } else { 2 };
    let (params, l, single) = like_inputs(a, name, np)?;
    let want = nargout > 1;
    let (nlogl, avar_v, avar_dims): (f64, Vec<f64>, (usize, usize)) = match name {
        "normlike" => {
            let (v, av) = normlike_core(&params, &l, want);
            (v, av, (2, 2))
        }
        "lognlike" => {
            let mut l2 = Like { data: l.data.iter().map(|&d| if d < 0.0 { f64::NAN } else { d.ln() }).collect(), cens: l.cens.clone(), freq: l.freq.clone() };
            let (v, av) = normlike_core(&params, &l2, want);
            let extra: f64 = (0..l2.data.len()).map(|i| l.freq[i] * l2.data[i] * (1.0 - l.cens[i])).sum();
            l2.data.clear();
            (v + extra, av, (2, 2))
        }
        _ => {
            let mut mu = params[0];
            if mu <= 0.0 {
                mu = f64::NAN;
            }
            let n: f64 = l.freq.iter().sum();
            let sumz: f64 = (0..l.data.len()).map(|i| l.freq[i] * if l.data[i] < 0.0 { f64::NAN } else { l.data[i] }).sum::<f64>() / mu;
            let nunc = n - (0..l.data.len()).map(|i| l.freq[i] * l.cens[i]).sum::<f64>();
            let v = nunc * mu.ln() + sumz;
            (v, vec![mu * mu / (2.0 * sumz - nunc)], (1, 1))
        }
    };
    let mut outs = vec![mat(1, 1, vec![nlogl], single)];
    if want {
        outs.push(mat(avar_dims.0, avar_dims.1, avar_v, single));
    }
    Ok(Some(outs))
}

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    let _ = normcdf_std;
    match name {
        "normfit" => normfit(a, nargout),
        "expfit" => expfit(a, nargout),
        "poissfit" => poissfit(a, nargout),
        "binofit" => binofit(a, nargout),
        "lognfit" => lognfit(a, nargout),
        "unifit" => unifit(a, nargout),
        "raylfit" => raylfit(a, nargout),
        _ => like(name, a, nargout),
    }
}
