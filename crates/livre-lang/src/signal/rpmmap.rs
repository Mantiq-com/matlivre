//! `rpmfreqmap`: short-time Fourier amplitude map against rotational speed.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["rpmfreqmap"];

pub(crate) fn window(it: &mut Interpreter, name: &str, l: usize, par: f64) -> R<Vec<f64>> {
    let args = match name {
        "kaiser" => vec![scalar_value(l as f64), scalar_value(if par.is_nan() { 0.5 } else { par })],
        "chebwin" => vec![scalar_value(l as f64), scalar_value(if par.is_nan() { 100.0 } else { par })],
        _ => vec![scalar_value(l as f64)],
    };
    Ok(it.call_function(name, args, 1)?.remove(0).to_matrix()?.data)
}

/// Resolution bandwidth actually obtained with a window of length `l`.
pub(crate) fn res_of(w: &[f64], fs: f64) -> f64 {
    let s: f64 = w.iter().sum();
    fs * w.iter().map(|v| v * v).sum::<f64>() / (s * s)
}

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    if nargout > 5 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected X to be one of these types:  single, double"),
    };
    let rp = match &a[2] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected RPM to be one of these types:  single, double"),
    };
    if !(x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1)) {
        return rt("Expected X to be a vector.");
    }
    let fsm = a[1].to_matrix()?;
    if fsm.numel() != 1 {
        return rt("Expected Fs to be a scalar.");
    }
    let fs = fsm.data[0];
    if !(fs > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    let n = x.numel();
    if rp.numel() != n {
        return rt("Expected signal and rpm vectors to be of same length.");
    }
    if rp.data.iter().any(|v| !(*v > 0.0)) {
        return rt("Expected RPM to be positive.");
    }
    let mut i = 3;
    let mut res = fs / 128.0;
    if let Some(Value::Num(m)) = a.get(3) {
        if m.kind != Kind::Char {
            if m.numel() != 1 {
                return rt("Expected resolution to be a scalar.");
            }
            if !(m.data[0] > 0.0) {
                return rt("Expected resolution to be positive.");
            }
            res = m.data[0];
            i = 4;
        }
    }
    let mut amp = "rms".to_string();
    let mut scale = "linear".to_string();
    let mut ov = 50.0;
    let mut win = "hann".to_string();
    while i < a.len() {
        let s = a[i].text().unwrap_or_default();
        let Some(v) = a.get(i + 1) else { return rt(format!("No value was given for '{s}'. Name-value pair arguments require a name followed by a value.")) };
        match s.to_lowercase().as_str() {
            "amplitude" => {
                let t = v.text().unwrap_or_default().to_lowercase();
                if !matches!(t.as_str(), "rms" | "peak" | "power") {
                    return rt(format!("Expected Amplitude to match one of these values:  'rms', 'peak', 'power'  The input, '{}', did not match any of the valid values.", v.text().unwrap_or_default()));
                }
                amp = t;
            }
            "scale" => {
                let t = v.text().unwrap_or_default();
                match t.to_lowercase().as_str() {
                    "linear" => scale = "linear".into(),
                    "db" => scale = "db".into(),
                    _ => return rt(format!("Expected Scale to match one of these values:  'linear', 'dB'  The input, '{t}', did not match any of the valid values.")),
                }
            }
            "overlappercent" => {
                let m = v.to_matrix()?;
                if m.numel() != 1 {
                    return rt("Expected OP to be a scalar.");
                }
                ov = m.data[0];
                if ov < 0.0 {
                    return rt("Expected OP to be nonnegative.");
                }
                if ov > 100.0 {
                    return rt("Expected OP to be a scalar with value <= 100.");
                }
            }
            "window" => {
                let t = v.text().unwrap_or_default();
                match t.to_lowercase().as_str() {
                    "hann" | "hamming" | "flattopwin" | "kaiser" | "rectwin" | "chebwin" => win = t.to_lowercase(),
                    _ => return rt(format!("Expected Window to match one of these values:  'hann', 'hamming', 'flattopwin', 'kaiser', 'rectwin', 'chebwin'  The input, '{t}', did not match any of the valid values.")),
                }
            }
            _ => return rt(format!("'{s}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    // Shortest window whose resolution bandwidth does not exceed `res`.
    let mut l = {
        let w = window(it, &win, 129, f64::NAN)?;
        let enbw = res_of(&w, fs) / fs * 129.0;
        ((enbw * fs / res) as usize).max(3)
    };
    let mut guard = 0;
    loop {
        let r = res_of(&window(it, &win, l, f64::NAN)?, fs);
        if r > res * (1.0 + 1e-12) {
            l += 1;
        } else if l > 3 && res_of(&window(it, &win, l - 1, f64::NAN)?, fs) <= res * (1.0 + 1e-12) {
            l -= 1;
        } else {
            break;
        }
        guard += 1;
        if guard > 2000 {
            break;
        }
    }
    if l > n {
        l = n;
    }
    let w = window(it, &win, l, f64::NAN)?;
    let res_out = res_of(&w, fs);
    let hop = ((l as f64 * (1.0 - ov / 100.0)).floor() as usize).max(1);
    let nseg = (n - l) / hop + 1;
    let nfft = l.max(256);
    let nbin = nfft / 2 + 1;
    let sw: f64 = w.iter().sum();
    let mut map = vec![0.0; nbin * nseg];
    let mut rpms = vec![0.0; nseg];
    let mut tms = vec![0.0; nseg];
    for k in 0..nseg {
        let st = k * hop;
        let mut seg: Vec<f64> = (0..l).map(|j| x.data[st + j] * w[j]).collect();
        seg.resize(nfft, 0.0);
        let f = it.call_function("fft", vec![real_col(seg)], 1)?.remove(0).to_matrix()?;
        let fi = f.im.clone().unwrap_or_else(|| vec![0.0; nfft]);
        for b in 0..nbin {
            let mag = f.data[b].hypot(fi[b]) / sw;
            let edge = b == 0 || (nfft % 2 == 0 && b == nfft / 2);
            let v = match amp.as_str() {
                "rms" => if edge { mag } else { mag * 2f64.sqrt() },
                "peak" => if b == 0 { mag } else if edge { 2f64.sqrt() * mag } else { 2.0 * mag },
                _ => if edge { mag * mag } else { 2.0 * mag * mag },
            };
            map[k * nbin + b] = if scale == "db" {
                if amp == "power" { 10.0 * v.log10() } else { 20.0 * v.log10() }
            } else {
                v
            };
        }
        let pos = st as f64 + (l as f64 - 1.0) / 2.0 + 0.5;
        let i0 = (pos.floor() as usize).min(n - 1);
        let i1 = (i0 + 1).min(n - 1);
        rpms[k] = rp.data[i0] + (pos - i0 as f64) * (rp.data[i1] - rp.data[i0]);
        tms[k] = pos / fs;
    }
    let freq: Vec<f64> = (0..nbin).map(|b| fs * b as f64 / nfft as f64).collect();
    let single = x.kind == Kind::Single;
    let mut mm = Matrix::from_dims(vec![nbin, nseg], map);
    let mut fm = Matrix::from_dims(vec![nbin, 1], freq);
    let mut rm = Matrix::from_dims(vec![nseg, 1], rpms);
    let mut tm = Matrix::from_dims(vec![nseg, 1], tms);
    let mut sm = Matrix::from_dims(vec![1, 1], vec![res_out]);
    if single {
        for m in [&mut mm, &mut fm, &mut rm, &mut tm] {
            m.kind = Kind::Single;
        }
    }
    out(vec![Value::Num(mm), Value::Num(fm), Value::Num(rm), Value::Num(tm), Value::Num(sm)].into_iter().take(nargout.max(1)).collect())
}
