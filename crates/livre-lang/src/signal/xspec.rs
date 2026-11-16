//! `xspectrogram`: cross-spectrogram built from the two short-time spectra.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;
use std::f64::consts::PI;

pub const NAMES: &[&str] = &["xspectrogram"];

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    if nargout > 4 {
        return rt("Too many output arguments.");
    }
    {
        let (mx, my) = (a[0].to_matrix()?, a[1].to_matrix()?);
        let isvec = |m: &Matrix| m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1);
        if !isvec(&mx) {
            return rt("Expected input number 1, X, to be a vector.");
        }
        if !isvec(&my) {
            return rt("Expected input number 2, Y, to be a vector.");
        }
        if isvec(&mx) && isvec(&my) && mx.numel() != my.numel() {
            return rt("The number of rows of the two inputs must be equal to calculate the cross spectral density.");
        }
    }
    // Positional inputs: window, noverlap, nfft (or f), Fs; then options.
    let mut pos: Vec<Value> = vec![];
    let mut opts: Vec<Value> = vec![];
    let mut minthr = f64::NEG_INFINITY;
    let mut downrows = false;
    let mut i = 2;
    while i < a.len() {
        let is_text = matches!(&a[i], Value::Num(m) if m.kind == Kind::Char);
        if !is_text {
            if opts.is_empty() {
                pos.push(a[i].clone());
                i += 1;
                continue;
            }
            return rt("Invalid input argument.");
        }
        let t = a[i].text().unwrap_or_default();
        match t.to_lowercase().as_str() {
            "minthreshold" => {
                if let Some(v) = a.get(i + 1) {
                    match v {
                        Value::Num(m) if m.kind != Kind::Char && m.numel() == 1 => minthr = m.data[0],
                        _ => {
                            return rt("Expected input to match one of these values:  'half', 'onesided', 'whole', 'twosided', 'adapt', 'unity', 'eigen', 'corr', 'ev', 'centered', 'power', 'psd', 'ms', 'reassigned', 'confidencelevel', 'mean', 'maxhold', 'minhold', 'mimo'  The input, 'minthreshold', did not match any of the valid values.");
                        }
                    }
                }
                i += 2;
            }
            "outputtimedimension" => {
                let v = a.get(i + 1).and_then(|v| v.text().ok()).unwrap_or_default().to_lowercase();
                downrows = v == "downrows";
                i += 2;
            }
            "xaxis" | "yaxis" => i += 1,
            _ => {
                opts.push(a[i].clone());
                i += 1;
            }
        }
    }
    let has_fs = pos.len() >= 4 && !matches!(&pos[3], Value::Num(m) if m.numel() == 0);
    let empty = || Value::Num(Matrix::from_dims(vec![0, 0], vec![]));
    while pos.len() < 3 {
        pos.push(empty());
    }
    let build = |sig: &Value| -> Vec<Value> {
        let mut v = vec![sig.clone(), pos[0].clone(), pos[1].clone(), pos[2].clone(), if has_fs { pos[3].clone() } else { scalar_value(2.0 * PI) }];
        v.extend(opts.iter().cloned());
        v
    };
    if matches!(&pos[2], Value::Num(m) if m.numel() > 1) {
        return vector_freq(it, a, &pos, has_fs, minthr, nargout);
    }
    // The spectrogram call needs a window argument; an empty one selects the default of eight segments.
    let rx = it.call_function("spectrogram", build(&a[0]), 4)?;
    let ry = it.call_function("spectrogram", build(&a[1]), 4)?;
    let sx = rx[0].to_matrix()?;
    let sy = ry[0].to_matrix()?;
    let px = rx[3].to_matrix()?;
    let (nr, nc) = (sx.dims[0], sx.dims[1]);
    let sxi = sx.im.clone().unwrap_or_else(|| vec![0.0; sx.numel()]);
    let syi = sy.im.clone().unwrap_or_else(|| vec![0.0; sy.numel()]);
    // Per-row scaling from the power spectrogram of x.
    let mut scale = vec![0.0; nr];
    for r in 0..nr {
        let mut best = 0.0;
        for c in 0..nc {
            let k = c * nr + r;
            let m2 = sx.data[k] * sx.data[k] + sxi[k] * sxi[k];
            if m2 > best && k < px.data.len() {
                best = m2;
                scale[r] = px.data[k] / m2;
            }
        }
    }
    let mut cre = vec![0.0; nr * nc];
    let mut cim = vec![0.0; nr * nc];
    let mut s = vec![0.0; nr * nc];
    for c in 0..nc {
        for r in 0..nr {
            let k = c * nr + r;
            let (xr, xi, yr, yi) = (sx.data[k], sxi[k], sy.data[k], syi[k]);
            // x * conj(y)
            let (re, im) = (xr * yr + xi * yi, xi * yr - xr * yi);
            cre[k] = re * scale[r];
            cim[k] = im * scale[r];
            s[k] = cre[k].hypot(cim[k]);
            if minthr.is_finite() && 10.0 * s[k].log10() <= minthr {
                s[k] = 0.0;
            }
        }
    }
    let centered = opts.iter().any(|o| o.text().map(|t| t.to_lowercase() == "centered").unwrap_or(false));
    if centered && a[0].to_matrix()?.im.is_none() && a[1].to_matrix()?.im.is_none() {
        for c in 0..nc {
            for r in 0..nr / 2 {
                cim[c * nr + r] = -cim[c * nr + r];
            }
        }
    }
    let single = a[0].to_matrix()?.kind == Kind::Single || a[1].to_matrix()?.kind == Kind::Single;
    let mut sm = Matrix::from_dims(vec![nr, nc], s);
    let mut cm = Matrix::complex(vec![nr, nc], cre, cim);
    let mut fm = rx[1].to_matrix()?;
    let mut tm = rx[2].to_matrix()?;
    if downrows {
        let tr = |m: &Matrix| -> Matrix {
            let mut d = vec![0.0; nr * nc];
            for c in 0..nc {
                for r in 0..nr {
                    d[r * nc + c] = m.data[c * nr + r];
                }
            }
            Matrix::from_dims(vec![nc, nr], d)
        };
        sm = tr(&sm);
        let ci = cm.im.clone().unwrap();
        let (a1, a2) = (tr(&Matrix::from_dims(vec![nr, nc], cm.data.clone())), tr(&Matrix::from_dims(vec![nr, nc], ci)));
        cm = Matrix::complex(vec![nc, nr], a1.data, a2.data);
        tm = Matrix::from_dims(vec![nc, 1], tm.data.clone());
    }
    if single {
        sm.kind = Kind::Single;
        cm.kind = Kind::Single;
        fm.kind = Kind::Single;
        tm.kind = Kind::Single;
    }
    out(vec![Value::Num(sm), Value::Num(fm), Value::Num(tm), Value::Num(cm)].into_iter().take(nargout.max(1)).collect())
}

/// Cross spectrogram at an explicit frequency vector (direct DFT of every segment).
fn vector_freq(it: &mut Interpreter, a: &[Value], pos: &[Value], has_fs: bool, minthr: f64, nargout: usize) -> Ret {
    let x = a[0].to_matrix()?;
    let y = a[1].to_matrix()?;
    let n = x.numel();
    if y.numel() != n {
        return rt("The number of rows of the two inputs must be equal to calculate the cross spectral density.");
    }
    let fs = if has_fs { pos[3].to_matrix()?.data[0] } else { 2.0 * PI };
    let fv = pos[2].to_matrix()?.data;
    let wm = pos[0].to_matrix()?;
    let (w, l): (Vec<f64>, usize) = if wm.numel() == 1 {
        let l = wm.data[0] as usize;
        (it.call_function("hamming", vec![scalar_value(l as f64)], 1)?.remove(0).to_matrix()?.data, l)
    } else if wm.numel() == 0 {
        let l = (n as f64 / 4.5).floor() as usize;
        let l = ((n - (l as f64 * 0.0) as usize) as f64 / 4.5).floor() as usize;
        let _ = l;
        return rt("xspectrogram: an empty window with a frequency vector is not implemented.");
    } else {
        (wm.data.clone(), wm.numel())
    };
    let nm = pos[1].to_matrix()?;
    let nov = if nm.numel() == 0 { l / 2 } else { nm.data[0] as usize };
    if nov >= l {
        return rt("The number of samples to overlap must be less than the length of the segments.");
    }
    let nseg = (n - nov) / (l - nov);
    let sw2: f64 = w.iter().map(|v| v * v).sum();
    let nf = fv.len();
    let mut cre = vec![0.0; nf * nseg];
    let mut cim = vec![0.0; nf * nseg];
    let mut s = vec![0.0; nf * nseg];
    let dft = |sig: &Matrix, st: usize, f: f64| -> (f64, f64) {
        let (mut re, mut im) = (0.0, 0.0);
        for j in 0..l {
            let ph = -2.0 * PI * f / fs * j as f64;
            let (sn, co) = ph.sin_cos();
            let v = sig.data[st + j] * w[j];
            re += v * co;
            im += v * sn;
        }
        (re, im)
    };
    for sg in 0..nseg {
        let st = sg * (l - nov);
        for (r, &f) in fv.iter().enumerate() {
            let (xr, xi) = dft(&x, st, f);
            let (yr, yi) = dft(&y, st, f);
            let k = sg * nf + r;
            let sc = 1.0 / (fs * sw2);
            cre[k] = (xr * yr + xi * yi) * sc;
            cim[k] = (xi * yr - xr * yi) * sc;
            s[k] = cre[k].hypot(cim[k]);
            if minthr.is_finite() && 10.0 * s[k].log10() <= minthr {
                s[k] = 0.0;
                cre[k] = 0.0;
                cim[k] = 0.0;
            }
        }
    }
    let t: Vec<f64> = (0..nseg).map(|sg| ((sg * (l - nov)) as f64 + l as f64 / 2.0) / fs).collect();
    out(vec![
        Value::Num(Matrix::from_dims(vec![nf, nseg], s)),
        real_col(fv),
        Value::Num(Matrix::from_dims(vec![1, nseg], t)),
        Value::Num(Matrix::complex(vec![nf, nseg], cre, cim)),
    ]
    .into_iter()
    .take(nargout.max(1))
    .collect())
}
