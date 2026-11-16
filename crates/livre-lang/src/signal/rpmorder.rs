//! `rpmordermap`: order-RPM map from a signal resampled to constant samples per revolution.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::signal::rpmmap::{res_of, window};
use crate::value::*;

pub const NAMES: &[&str] = &["rpmordermap"];

/// Shortest window length whose resolution bandwidth (at sample rate `fs`) does not exceed `res`.
fn solve_len(it: &mut Interpreter, win: &str, par: f64, fs: f64, res: f64) -> R<usize> {
    let w = window(it, win, 129, par)?;
    let enbw = res_of(&w, fs) / fs * 129.0;
    let mut l = ((enbw * fs / res) as usize).max(3);
    for _ in 0..4000 {
        let r = res_of(&window(it, win, l, par)?, fs);
        if r > res * (1.0 + 1e-12) {
            l += 1;
        } else if l > 3 && res_of(&window(it, win, l - 1, par)?, fs) <= res * (1.0 + 1e-12) {
            l -= 1;
        } else {
            break;
        }
    }
    Ok(l)
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
    let mut res_in: Option<f64> = None;
    if let Some(Value::Num(m)) = a.get(3) {
        if m.kind != Kind::Char {
            if m.numel() != 1 {
                return rt("Expected resolution to be a scalar.");
            }
            if !(m.data[0] > 0.0) {
                return rt("Expected resolution to be positive.");
            }
            res_in = Some(m.data[0]);
            i = 4;
        }
    }
    let mut amp = "rms".to_string();
    let mut scale = "linear".to_string();
    let mut ov = 50.0;
    let mut win = "flattopwin".to_string();
    let mut par = f64::NAN;
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
                let (t, p) = match v {
                    Value::Cell(c) => {
                        let first = c.data.first().and_then(|e| e.text().ok()).unwrap_or_default();
                        let p = c.data.get(1).and_then(|e| e.scalar().ok()).unwrap_or(f64::NAN);
                        (first, p)
                    }
                    _ => (v.text().unwrap_or_default(), f64::NAN),
                };
                match t.to_lowercase().as_str() {
                    "hann" | "hamming" | "flattopwin" | "kaiser" | "rectwin" | "chebwin" => {
                        win = t.to_lowercase();
                        if win == "chebwin" && !p.is_nan() && p < 45.0 {
                            return rt("Expected chebwin sidelobe attenuation parameter to be a scalar with value >= 45.");
                        }
                        if win == "kaiser" && !p.is_nan() && !(p > 0.0) {
                            return rt("Expected Kaiser beta parameter to be positive.");
                        }
                        par = p;
                    }
                    _ => return rt(format!("Expected Window to match one of these values:  'hann', 'hamming', 'flattopwin', 'kaiser', 'rectwin', 'chebwin'  The input, '{t}', did not match any of the valid values.")),
                }
            }
            _ => return rt(format!("'{s}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    let ts: Vec<f64> = (0..n).map(|k| k as f64 / fs).collect();
    // Phase in revolutions.
    let up = it.call_function("resample", vec![real_col(x.data.clone()), scalar_value(15.0), scalar_value(1.0)], 1)?.remove(0).to_matrix()?;
    let nu = up.numel();
    let tf: Vec<f64> = (0..nu).map(|k| k as f64 / (15.0 * fs)).collect();
    let interp = |it: &mut Interpreter, xs: &[f64], ys: &[f64], q: Vec<f64>| -> R<Vec<f64>> {
        Ok(it.call_function("interp1", vec![real_col(xs.to_vec()), real_col(ys.to_vec()), real_col(q), Value::string("linear"), Value::string("extrap")], 1)?.remove(0).to_matrix()?.data)
    };
    let rpm_f = interp(it, &ts, &rp.data, tf.clone())?;
    let mut phase_f = vec![0.0; nu];
    for k in 1..nu {
        phase_f[k] = phase_f[k - 1] + (rpm_f[k] + rpm_f[k - 1]) / 120.0 / (15.0 * fs);
    }
    let rmax = rp.data.iter().cloned().fold(0.0, f64::max);
    let omax = fs / (2.0 * rmax / 60.0);
    let fp = 8.0 * omax;
    let m = (phase_f[nu - 1] * fp).floor() as usize + 1;
    let grid: Vec<f64> = (0..m).map(|k| k as f64 / fp).collect();
    let xp = interp(it, &phase_f, &up.data, grid)?;
    let res = res_in.unwrap_or(fp / 256.0);
    let mut l = solve_len(it, &win, par, fp, res)?;
    if l > m {
        l = m;
    }
    let w = window(it, &win, l, par)?;
    let res_out = res_of(&w, fp);
    let nov = ((ov / 100.0 * l as f64).ceil() as usize).min(l - 1);
    let hop = l - nov;
    let nseg = (m - l) / hop + 1;
    let nbin = (l / 8 + 2).min(l / 2 + 1);
    let sw: f64 = w.iter().sum();
    let mut map = vec![0.0; nbin * nseg];
    let mut centers = vec![0.0; nseg];
    for k in 0..nseg {
        let st = k * hop;
        let seg: Vec<f64> = (0..l).map(|j| xp[st + j] * w[j]).collect();
        let f = it.call_function("fft", vec![real_col(seg)], 1)?.remove(0).to_matrix()?;
        let fi = f.im.clone().unwrap_or_else(|| vec![0.0; l]);
        for b in 0..nbin {
            let mag = f.data[b].hypot(fi[b]) / sw;
            let edge = b == 0 || (l % 2 == 0 && b == l / 2);
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
        centers[k] = (st as f64 + l as f64 / 2.0) / fp;
    }
    let tcol = interp(it, &phase_f, &tf, centers)?;
    let rcol = interp(it, &ts, &rp.data, tcol.clone())?;
    let orders: Vec<f64> = (0..nbin).map(|b| fp * b as f64 / l as f64).collect();
    let mut mm = Matrix::from_dims(vec![nbin, nseg], map);
    let mut om = Matrix::from_dims(vec![nbin, 1], orders);
    let mut rm = Matrix::from_dims(vec![nseg, 1], rcol);
    let mut tm = Matrix::from_dims(vec![nseg, 1], tcol);
    let sm = Matrix::from_dims(vec![1, 1], vec![res_out]);
    if x.kind == Kind::Single {
        for q in [&mut mm, &mut om, &mut rm, &mut tm] {
            q.kind = Kind::Single;
        }
    }
    out(vec![Value::Num(mm), Value::Num(om), Value::Num(rm), Value::Num(tm), Value::Num(sm)].into_iter().take(nargout.max(1)).collect())
}
