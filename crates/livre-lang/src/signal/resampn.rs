//! `resample` for numeric signals: rational-rate resampling with a Kaiser-windowed
//! lowpass FIR filter (the timeseries `resample` lives in `timeseries.rs`).

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if nargout > 2 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected input number 1, X, to be one of these types:  double, single"),
    };
    if x.numel() == 0 {
        return rt("Index exceeds the number of array elements. Index must not exceed 0.");
    }
    if a.len() < 3 {
        return rt("Downsampling parameter Q must be specified together with upsampling parameter P.");
    }
    if a.len() > 5 {
        return rt("Expected input number 2 to be time instants Tx. If the second input represents upsampling factor P, specify it as a scalar at compile time.");
    }
    let scalar = |v: &Value| -> Option<f64> { v.to_matrix().ok().filter(|m| m.numel() == 1 && m.kind != Kind::Char).map(|m| m.data[0]) };
    let p = scalar(&a[1]).unwrap_or(f64::NAN);
    if !(p > 0.0) {
        return rt("Expected input number 2, P, to be positive.");
    }
    if p.fract() != 0.0 {
        return rt("Expected input number 2, P, to be integer-valued.");
    }
    if let Value::Num(m) = &a[2] {
        if m.kind == Kind::Char {
            let t = a[2].text().unwrap_or_default();
            if t == "spline" {
                return rt("Interpolation method 'spline' is supported when time vector Tx is an input.");
            }
            return rt(format!("Expected input to match one of these values:  'linear', 'pchip', 'spline', 'Dimension'  The input, '{t}', did not match any of the valid values."));
        }
    }
    if let Some(Value::Num(m)) = a.get(3) {
        if m.kind == Kind::Char {
            let t = a[3].text().unwrap_or_default();
            if t == "spline" {
                return rt("Interpolation method 'spline' is supported when time vector Tx is an input.");
            }
            return rt(format!("Expected input to match one of these values:  'linear', 'pchip', 'spline', 'Dimension'  The input, '{t}', did not match any of the valid values."));
        }
    }
    let q = scalar(&a[2]).unwrap_or(f64::NAN);
    if !(q > 0.0) {
        return rt("Expected input number 3, Q, to be positive.");
    }
    if q.fract() != 0.0 {
        return rt("Expected input number 3, Q, to be integer-valued.");
    }
    let n = a.get(3).and_then(|v| scalar(v)).unwrap_or(10.0);
    let beta = a.get(4).and_then(|v| scalar(v)).unwrap_or(5.0);
    if beta < 0.0 {
        return rt("Expected BTA to be nonnegative.");
    }
    let g = gcd(p as u64, q as u64) as f64;
    let (p, q) = (p / g, q / g);
    let (pu, qu) = (p as usize, q as usize);
    let vector = x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1);
    let (len, cols) = if vector { (x.numel(), 1) } else { (x.dims[0], x.numel() / x.dims[0]) };
    let ny = (len * pu + qu - 1) / qu;
    // Filter.
    let mut h: Vec<f64> = vec![];
    if n >= 1.0 {
        let pq = pu.max(qu);
        let l = 2 * (n as usize) * pq + 1;
        let kw = it.call_function("kaiser", vec![scalar_value(l as f64), scalar_value(beta)], 1)?.remove(0).to_matrix()?.data;
        let mid = (l as f64 - 1.0) / 2.0;
        h = (0..l)
            .map(|k| {
                let m = (k as f64 - mid) / pq as f64;
                let s = if m == 0.0 { 1.0 } else { (std::f64::consts::PI * m).sin() / (std::f64::consts::PI * m) };
                s * kw[k]
            })
            .collect();
        let sum: f64 = h.iter().sum();
        h.iter_mut().for_each(|v| *v *= p / sum);
    }
    let mut res = vec![0.0; ny * cols];
    for c in 0..cols {
        let col = &x.data[c * len..(c + 1) * len];
        if h.is_empty() {
            for j in 0..ny {
                let src = ((j as f64 * q / p) + 0.5).floor() as usize;
                res[c * ny + j] = col[src.min(len - 1)];
            }
        } else {
            let l = h.len();
            let d = (l - 1) / 2;
            for j in 0..ny {
                let idx = j * qu + d;
                // x[i] sits at upsampled position i*p; h index = idx - i*p must be in [0, l).
                let lo = if idx + 1 > l { (idx + 1 - l + pu - 1) / pu } else { 0 };
                let hi = (idx / pu).min(len - 1);
                let mut s = 0.0;
                let mut i = lo;
                while i <= hi {
                    s += col[i] * h[idx - i * pu];
                    i += 1;
                }
                res[c * ny + j] = s;
            }
        }
    }
    let dims = if vector { if x.dims[0] == 1 { vec![1, ny] } else { vec![ny, 1] } } else { vec![ny, cols] };
    let mut m = Matrix::from_dims(dims, res);
    if x.kind == Kind::Single {
        m.kind = Kind::Single;
    }
    let mut outv = vec![Value::Num(m)];
    if nargout > 1 {
        outv.push(real_row(h));
    }
    out(outv)
}
