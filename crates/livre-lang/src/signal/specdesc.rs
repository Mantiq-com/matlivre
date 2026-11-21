//! Spectral descriptors of audio signals and spectrograms: `spectralCrest`, `spectralEntropy`,
//! `spectralFlatness`, `spectralKurtosis` and `spectralSkewness`.
//!
//! A time-domain input is cut into frames (window `rectwin(round(0.03*fs))`, overlap
//! `round(0.02*fs)`), each frame is windowed and transformed, and the one-sided spectrum (power or
//! magnitude, interior bins doubled) over the chosen range gives the descriptor. A frequency
//! domain input (`x` with a frequency vector `f`) is used as the spectrum itself.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["spectralCrest", "spectralEntropy", "spectralFlatness", "spectralKurtosis", "spectralSkewness"];

const WIN_HINT: &str = "The default window length depends on the specified sample rate: round(f*0.03).";

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let max_out = if name == "spectralCrest" { 3 } else { 1 };
    if nargout > max_out {
        return rt("One or more output arguments not assigned during call to \"varargout\".");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        Value::Num(m) => return rt(format!("Expected x to be one of these types:\n\nsingle, double\n\nInstead its type was {}.", m.kind.name())),
        other => return rt(format!("Expected x to be one of these types:\n\nsingle, double\n\nInstead its type was {}.", other.class_name())),
    };
    if x.im.is_some() {
        return rt("Expected x to be real.");
    }
    let fm = match &a[1] {
        Value::Num(m) if m.kind != Kind::Char => m.clone(),
        _ => return rt("Expected f to be numeric."),
    };
    let single = x.kind == Kind::Single;
    // Options.
    let mut window: Option<Vec<f64>> = None;
    let mut overlap: Option<f64> = None;
    let mut nfft: Option<f64> = None;
    let mut range: Option<Vec<f64>> = None;
    let mut power = true;
    let mut i = 2;
    if (a.len() - 2) % 2 != 0 {
        return rt("Invalid property/value pair arguments.");
    }
    while i + 1 < a.len() {
        let key = a[i].text().unwrap_or_default().to_lowercase();
        let v = &a[i + 1];
        match key.as_str() {
            "window" => window = Some(v.to_matrix()?.data),
            "overlaplength" => overlap = Some(v.scalar()?),
            "fftlength" => nfft = Some(v.scalar()?),
            "range" => range = Some(v.to_matrix()?.data),
            "spectrumtype" => {
                let t = v.text().unwrap_or_default().to_lowercase();
                power = match t.as_str() {
                    "power" => true,
                    "magnitude" => false,
                    _ => return rt("Error setting property 'SpectrumType' of class 'audio.internal.spectralDescriptors.SpectralDescriptorValidator'. Value must be a member of this set:\n    'power'\n    'magnitude'\n"),
                };
            }
            _ => return rt(format!("Unrecognized property '{}' for class 'audio.internal.spectralDescriptors.SpectralDescriptorValidator'.", a[i].text().unwrap_or_default())),
        }
        i += 2;
    }
    // Spectra: one column per (frame, channel), frequencies and the output shape.
    let (spectra, freqs, oshape): (Vec<Vec<f64>>, Vec<f64>, Vec<usize>);
    if fm.numel() > 1 {
        // Frequency-domain input: rows follow f.
        let rows = x.dims[0];
        if rows != fm.numel() {
            return rt("Invalid frequency vector. For frequency domain input, the number of rows must match the number of elements in the frequency vector:\nsize(x,1)==numel(f)");
        }
        let cols = x.numel() / rows.max(1);
        spectra = (0..cols).map(|c| x.data[c * rows..(c + 1) * rows].to_vec()).collect();
        freqs = fm.data.clone();
        oshape = if x.dims.len() >= 3 { vec![x.dims[1], x.dims[2]] } else { vec![x.dims.get(1).copied().unwrap_or(1), 1] };
    } else {
        let fs = fm.data[0];
        if fs < 0.0 {
            return rt("Expected f to be nonnegative.");
        }
        if x.dims.len() > 2 {
            return rt("Expected x to be a vector or a matrix.");
        }
        let (n, ch) = if x.dims[0] == 1 && x.dims[1] > 1 { (x.dims[1], 1) } else { (x.dims[0], x.dims[1]) };
        let win: Vec<f64> = match window {
            Some(w) => w,
            None => vec![1.0; (fs * 0.03).round() as usize],
        };
        let wl = win.len();
        if wl < 1 || wl > n || x.dims[0] == 1 && x.dims[1] > 1 {
            return rt(format!("Invalid window. The number of elements in the window must be in the range [1,size(x,1)], where x is the audio input.\n{WIN_HINT}"));
        }
        let ov = overlap.unwrap_or((fs * 0.02).round());
        if ov < 0.0 {
            return rt("Expected OVERLAPLENGTH to be nonnegative.");
        }
        if ov >= wl as f64 {
            return rt(format!("Invalid overlap length. OVERLAPLENGTH must be less than size(WINDOW,1).\n{WIN_HINT}"));
        }
        let nf = nfft.unwrap_or(wl as f64) as usize;
        if nf < wl {
            return rt(format!("The FFT length must be greater than or equal to the window length, {wl}."));
        }
        let hop = wl - ov as usize;
        let frames = (n - wl) / hop + 1;
        let half = nf / 2 + 1;
        let wsum: f64 = win.iter().sum();
        let mut sp = Vec::new();
        for c in 0..ch {
            for fr in 0..frames {
                let seg: Vec<f64> = (0..wl).map(|k| x.data[c * n + fr * hop + k] * win[k]).collect();
                let mut segv = real_col(seg);
                if single {
                    segv = match segv {
                        Value::Num(m) => Value::Num(m.cast_to(Kind::Single)),
                        o => o,
                    };
                }
                let r = it.call_function("fft", vec![segv, scalar_value(nf as f64)], 1)?.remove(0).to_matrix()?;
                let mut s: Vec<f64> = (0..half).map(|k| r.data[k].hypot(r.im_at(k)) / wsum).collect();
                if single {
                    s.iter_mut().for_each(|v| *v = *v as f32 as f64);
                }
                if power {
                    s.iter_mut().for_each(|v| *v = *v * *v);
                }
                let last = if nf % 2 == 0 { half - 1 } else { half };
                for v in s.iter_mut().take(last).skip(1) {
                    *v *= 2.0;
                }
                sp.push(s);
            }
        }
        spectra = sp;
        freqs = (0..half).map(|k| k as f64 * fs / nf as f64).collect();
        oshape = vec![frames, ch];
        // Range checks use the sample rate.
        if let Some(r) = &range {
            check_range(r, (fs / 2.0).floor())?;
        }
    }
    let (lo, hi) = match &range {
        Some(r) if fm.numel() == 1 && r.len() == 2 => (r[0], r[1]),
        _ => (f64::NEG_INFINITY, f64::INFINITY),
    };
    let keep: Vec<usize> = (0..freqs.len()).filter(|&k| freqs[k] >= lo && freqs[k] <= hi).collect();
    let eps = if single { f32::EPSILON as f64 } else { f64::EPSILON };
    let bins = |s: &[f64]| -> Vec<f64> { keep.iter().map(|&k| s[k]).collect() };
    let fsel: Vec<f64> = keep.iter().map(|&k| freqs[k]).collect();
    let mut main = Vec::new();
    let (mut aux1, mut aux2) = (Vec::new(), Vec::new());
    for s in &spectra {
        let v = bins(s);
        let k = v.len() as f64;
        let sum: f64 = v.iter().sum();
        match name {
            "spectralCrest" => {
                let peak = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let mean = sum / k;
                main.push(peak / mean);
                aux1.push(peak);
                aux2.push(mean);
            }
            "spectralEntropy" => {
                let e: f64 = if sum == 0.0 { 0.0 } else { v.iter().map(|x| x / sum).filter(|p| *p > 0.0).map(|p| -p * p.log2()).sum::<f64>() };
                main.push(e / k.log2());
            }
            "spectralFlatness" => {
                let g = (v.iter().map(|x| (x.abs() + eps).ln()).sum::<f64>() / k).exp();
                main.push(g / (v.iter().map(|x| x.abs()).sum::<f64>() / k));
            }
            _ => {
                let cen = fsel.iter().zip(&v).map(|(f, s)| f * s).sum::<f64>() / sum;
                let spread = (fsel.iter().zip(&v).map(|(f, s)| (f - cen) * (f - cen) * s).sum::<f64>() / sum).sqrt();
                let p = if name == "spectralSkewness" { 3 } else { 4 };
                main.push(fsel.iter().zip(&v).map(|(f, s)| (f - cen).powi(p) * s).sum::<f64>() / sum / spread.powi(p));
            }
        }
    }
    let shape = oshape.clone();
    let mk = |v: Vec<f64>| -> Value {
        let m = Matrix::from_dims(shape.clone(), v);
        let m = if single { m.cast_to(Kind::Single) } else { m };
        Value::Num(m)
    };
    let empty = |cols: usize| -> Value { Value::Num(Matrix::from_dims(vec![0, cols], vec![])) };
    let mut outs;
    if keep.is_empty() && name == "spectralCrest" {
        let cols = oshape.get(1).copied().unwrap_or(1);
        outs = vec![empty(cols), empty(cols), mk(vec![f64::NAN; spectra.len()])];
    } else if keep.is_empty() {
        let fill = if name == "spectralEntropy" { 0.0 } else { f64::NAN };
        outs = vec![mk(vec![fill; spectra.len()])];
    } else {
        outs = vec![mk(main)];
        if name == "spectralCrest" {
            outs.push(mk(aux1));
            outs.push(mk(aux2));
        }
    }
    out(outs.into_iter().take(nargout.max(1)).collect())
}

fn check_range(r: &[f64], max: f64) -> Result<(), crate::error::LangError> {
    if r.len() != 2 {
        return rt("Expected RANGE to be an array with number of columns equal to 2.");
    }
    if r[0] < 0.0 || r[1] < 0.0 {
        return rt("Expected RANGE to be nonnegative.");
    }
    if !(r[1] > r[0]) {
        return rt("Expected RANGE to be strictly increasing.");
    }
    if r[1] > max {
        return rt(format!("Expected RANGE(2) to be an array with all of the values <= {max}."));
    }
    Ok(())
}
