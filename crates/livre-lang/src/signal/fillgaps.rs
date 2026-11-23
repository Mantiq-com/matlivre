//! `fillgaps`: replace NaN gaps by a weighted average of a forward and a backward autoregressive
//! prediction.
//!
//! Gaps are filled left to right for the forward estimate (each model is a Burg fit, with the mean
//! removed, to all the samples before the gap, including the forward estimates of earlier gaps)
//! and right to left for the backward estimate; sample `k` of a gap of length `L` takes the
//! forward value with weight `1 - k/(L+1)` and the backward one with weight `k/(L+1)`. The order
//! is fixed or the one minimising `N*log(e) + 2*order`.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["fillgaps"];

enum Order {
    Aic,
    Fixed(f64),
}

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() > 3 {
        return rt("Too many input arguments.");
    }
    if nargout > 1 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if m.kind == Kind::Double => m.clone(),
        Value::Num(m) if m.kind == Kind::Logical => return rt("NaN values cannot be converted to logicals."),
        Value::Num(m) => return rt(format!("Expected X to be one of these types:\n\ndouble\n\nInstead its type was {}.", m.kind.name())),
        other => return rt(format!("Expected X to be one of these types:\n\ndouble\n\nInstead its type was {}.", other.class_name())),
    };
    if x.dims.len() != 2 {
        return rt("Expected X to be two-dimensional.");
    }
    if x.im.is_some() {
        return rt("fillgaps: complex input is not supported.");
    }
    if x.data.iter().any(|v| v.is_infinite()) {
        return rt("Expected X to be finite.");
    }
    // MAXLEN.
    let mut maxlen: Option<usize> = None;
    if let Some(v) = a.get(1) {
        match v {
            Value::Num(m) if m.kind == Kind::Char || m.kind == Kind::Logical => {
                return rt("Expected MAXLEN to be one of these types:\n\ndouble, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64");
            }
            Value::Num(m) if m.numel() == 0 => {}
            Value::Num(m) => {
                if m.numel() != 1 {
                    return rt("Expected MAXLEN to be a scalar.");
                }
                let v = m.data[0];
                if !(v > 2.0) {
                    return rt("Expected MAXLEN to be a scalar with value > 2.");
                }
                if v.is_finite() && v.fract() != 0.0 {
                    return rt("Expected MAXLEN to be integer-valued.");
                }
                if v.is_finite() {
                    maxlen = Some(v as usize);
                }
            }
            _ => return rt("Expected MAXLEN to be one of these types:\n\ndouble, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64"),
        }
    }
    // ORDER.
    let mut order = Order::Aic;
    if let Some(v) = a.get(2) {
        match v {
            Value::Num(m) if m.kind == Kind::Char => {
                let t = v.text().unwrap_or_default();
                if t.is_empty() || !"aic".starts_with(&t.to_lowercase()) {
                    return rt(format!("Expected ORDER to match one of these values:\n\n'aic'\n\nThe input, '{t}', did not match any of the valid values."));
                }
            }
            Value::Str(_) => {
                let t = v.text().unwrap_or_default();
                if t.is_empty() || !"aic".starts_with(&t.to_lowercase()) {
                    return rt(format!("Expected ORDER to match one of these values:\n\n'aic'\n\nThe input, '{t}', did not match any of the valid values."));
                }
            }
            Value::Num(m) => {
                if m.numel() != 1 {
                    return rt("Expected ORDER to be a scalar.");
                }
                let o = m.data[0];
                if !(o > 0.0) {
                    return rt("Expected ORDER to be positive.");
                }
                if o.is_finite() && o.fract() != 0.0 {
                    return rt("Expected ORDER to be integer-valued.");
                }
                order = if o.is_infinite() { Order::Aic } else { Order::Fixed(o) };
            }
            _ => return rt("Expected ORDER to be a scalar."),
        }
    }
    let row = x.dims[0] == 1 && x.dims[1] > 1;
    let (len, cols) = if row { (x.dims[1], 1) } else { (x.dims[0], x.dims[1]) };
    let mut res = x.data.clone();
    for c in 0..cols {
        let col: Vec<f64> = x.data[c * len..(c + 1) * len].to_vec();
        let y = fill_channel(it, &col, maxlen, &order)?;
        res[c * len..(c + 1) * len].copy_from_slice(&y);
    }
    let mut m = Matrix::from_dims(x.dims.clone(), res);
    m.kind = Kind::Double;
    out(vec![Value::Num(m)])
}

/// The NaN runs of `v` as `(start, end)` inclusive index pairs.
fn gaps(v: &[f64]) -> Vec<(usize, usize)> {
    let mut g = Vec::new();
    let mut i = 0;
    while i < v.len() {
        if v[i].is_nan() {
            let s = i;
            while i + 1 < v.len() && v[i + 1].is_nan() {
                i += 1;
            }
            g.push((s, i));
        }
        i += 1;
    }
    g
}

/// `l` predictions that continue `data` with an AR model of the chosen order.
fn extend(it: &mut Interpreter, data: &[f64], l: usize, maxlen: Option<usize>, order: &Order) -> Result<Vec<f64>, crate::error::LangError> {
    let d: &[f64] = match maxlen {
        Some(m) if data.len() > m => &data[data.len() - m..],
        _ => data,
    };
    let n = d.len();
    let mean = d.iter().sum::<f64>() / n as f64;
    let dd: Vec<f64> = d.iter().map(|v| v - mean).collect();
    let fit = |it: &mut Interpreter, q: usize| -> Result<(Vec<f64>, f64), crate::error::LangError> {
        let r = it.call_function("arburg", vec![real_col(dd.clone()), scalar_value(q as f64)], 2)?;
        Ok((r[0].to_matrix()?.data, r[1].scalar().unwrap_or(f64::NAN)))
    };
    let q = match order {
        Order::Fixed(o) => {
            if let Some(m) = maxlen {
                if *o > m as f64 {
                    return rt(format!("Filter order {o} too large for specified search length {m}.  For automatic order selection, specify ORDER either as 'aic' or Inf."));
                }
            }
            (*o as usize).min(n.saturating_sub(1))
        }
        Order::Aic => {
            let mut best = (0usize, f64::INFINITY);
            for q in 1..n {
                let (_, e) = fit(it, q)?;
                let aic = n as f64 * e.ln() + 2.0 * q as f64;
                if aic < best.1 {
                    best = (q, aic);
                }
            }
            best.0
        }
    };
    let mut buf = dd.clone();
    let mut out = Vec::with_capacity(l);
    if q == 0 {
        return Ok(vec![mean; l]);
    }
    let (a, _) = fit(it, q)?;
    for _ in 0..l {
        let len = buf.len();
        let v: f64 = -(1..=q).map(|k| a[k] * buf[len - k]).sum::<f64>();
        buf.push(v);
        out.push(v + mean);
    }
    Ok(out)
}

fn fill_channel(it: &mut Interpreter, v: &[f64], maxlen: Option<usize>, order: &Order) -> Result<Vec<f64>, crate::error::LangError> {
    let g = gaps(v);
    if g.is_empty() || g.iter().map(|(s, e)| e - s + 1).sum::<usize>() == v.len() {
        return Ok(v.to_vec());
    }
    let n = v.len();
    // Forward estimates, left to right.
    let mut fwd = v.to_vec();
    for &(s, e) in &g {
        if s == 0 {
            continue;
        }
        let pre: Vec<f64> = fwd[..s].iter().cloned().filter(|x| !x.is_nan()).collect();
        let p = extend(it, &pre, e - s + 1, maxlen, order)?;
        fwd[s..=e].copy_from_slice(&p);
    }
    // Backward estimates, right to left.
    let mut bwd = v.to_vec();
    for &(s, e) in g.iter().rev() {
        if e + 1 >= n {
            continue;
        }
        let post: Vec<f64> = bwd[e + 1..].iter().rev().cloned().filter(|x| !x.is_nan()).collect();
        let mut p = extend(it, &post, e - s + 1, maxlen, order)?;
        p.reverse();
        bwd[s..=e].copy_from_slice(&p);
    }
    let mut y = v.to_vec();
    for &(s, e) in &g {
        let l = e - s + 1;
        for k in 0..l {
            let w = (k + 1) as f64 / (l + 1) as f64;
            let (f, b) = (fwd[s + k], bwd[s + k]);
            y[s + k] = if s == 0 {
                b
            } else if e + 1 >= n {
                f
            } else {
                (1.0 - w) * f + w * b
            };
        }
    }
    Ok(y)
}
