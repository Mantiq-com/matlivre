//! `instfreq` with the Hilbert method: the phase derivative of the analytic signal.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;
use std::f64::consts::PI;

pub const NAMES: &[&str] = &["instfreq", "instbw"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    let bw = name == "instbw";
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > 2 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        Value::Num(_) => return rt(if matches!(&a[0], Value::Num(m) if m.kind == Kind::Char || m.kind == Kind::Logical) { "Expected Input Signal to be one of these types:  single, double".to_string() } else { format!("Expected Input Signal to be one of these types:  single, double  Instead its type was {}.", a[0].class_name()) }),
        v => return rt(format!("Expected Input Signal to be one of these types:  single, double  Instead its type was {}.", v.class_name())),
    };
    if a.len() >= 3 && matches!(&a[2], Value::Num(m) if m.kind != Kind::Char) && matches!(&a[1], Value::Num(m) if m.kind != Kind::Char) {
        // Time-frequency distribution: tfd, fd, td.
        let fv = a[1].to_matrix()?.data;
        let tv = a[2].to_matrix()?.data;
        let (nf, nt) = (x.dims[0], x.numel() / x.dims[0].max(1));
        if fv.len() != nf || tv.len() != nt {
            return rt("Expected frequency and time vectors to match the dimensions of the time-frequency distribution.");
        }
        let mut data = vec![];
        for k in 0..nt {
            let col = &x.data[k * nf..(k + 1) * nf];
            let trapz = |g: &dyn Fn(usize) -> f64| -> f64 { (0..nf).map(|j| g(j)).sum() };
            let den = trapz(&|j| col[j]);
            let m1 = trapz(&|j| fv[j] * col[j]) / den;
            if bw {
                let m2 = trapz(&|j| (fv[j] - m1) * (fv[j] - m1) * col[j]) / den;
                data.push((4.0 * PI).sqrt() * m2.sqrt());
            } else {
                data.push(m1);
            }
        }
        return out(vec![Value::Num(Matrix::from_dims(vec![nt, 1], data)), real_col(tv)].into_iter().take(nargout.max(1)).collect());
    }
    if a.len() < 2 {
        return rt("Not enough input arguments. Specify a sample rate, a sample time, or a vector of time values.");
    }
    if x.numel() == 0 {
        return rt("Expected Input Signal to be nonempty.");
    }
    let vector = x.dims[0] == 1 || x.dims[1] == 1;
    let (n, ch) = if vector { (x.numel(), 1) } else { (x.dims[0], x.dims[1]) };
    // Time base.
    let tm = match &a[1] {
        Value::Num(m) if m.kind != Kind::Char => m.clone(),
        _ => return rt("Not enough input arguments. Specify a sample rate, a sample time, or a vector of time values."),
    };
    let mut method = "tfmoment".to_string();
    let mut limits = false;
    let mut flim: Option<[f64; 2]> = None;
    let mut scale = (4.0 * PI).sqrt();
    let mut i = 2;
    while i < a.len() {
        let s = a[i].text().unwrap_or_default();
        let low = s.to_lowercase();
        if i + 1 >= a.len() {
            if matches!(low.as_str(), "method" | "frequencylimits") {
                return rt(format!("No value was given for '{s}'. Name-value pair arguments require a name followed by a value."));
            }
            return rt(format!("'{s}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function."));
        }
        match low.as_str() {
            "method" => {
                let v = a[i + 1].text().unwrap_or_default();
                match v.to_lowercase().as_str() {
                    "tfmoment" | "hilbert" if !bw => method = v.to_lowercase(),
                    _ => return rt(format!("Expected Method to match one of these values:  'tfmoment', 'hilbert'  The input, '{v}', did not match any of the valid values.")),
                }
            }
            "frequencylimits" => {
                let m = a[i + 1].to_matrix()?;
                if m.numel() != 2 {
                    return rt("Expected FrequencyLimits to be an array with number of elements equal to 2.");
                }
                if !(m.data[1] > m.data[0]) {
                    return rt("Expected FrequencyLimits to be monotonically increasing.");
                }
                limits = true;
                flim = Some([m.data[0], m.data[1]]);
            }
            "scalefactor" if bw => {
                let m = a[i + 1].to_matrix()?;
                if m.numel() != 1 || !(m.data[0] > 0.0) {
                    return rt("Expected ScaleFactor to be positive.");
                }
                scale = m.data[0];
            }
            _ => return rt(format!("'{s}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    let (fs, times): (f64, Vec<f64>) = if tm.numel() == 1 {
        if !(tm.data[0] > 0.0) {
            return rt("Expected Fs to be positive.");
        }
        (tm.data[0], (0..n).map(|k| k as f64 / tm.data[0]).collect())
    } else {
        if tm.numel() != n {
            return rt("Time values vector length must equal input signal length.");
        }
        if tm.data.windows(2).any(|w| w[1] <= w[0]) {
            return rt("Expected Tv to be strictly increasing.");
        }
        let step = (tm.data[n - 1] - tm.data[0]) / (n as f64 - 1.0);
        if method == "hilbert" && tm.data.windows(2).any(|w| ((w[1] - w[0]) - step).abs() > 1e-9 * step.abs().max(1.0)) {
            return rt("Time values must be uniformly spaced when 'Method' value is 'hilbert'.");
        }
        (1.0 / step, tm.data.clone())
    };
    if method == "hilbert" && limits {
        return rt("'FrequencyLimits' are not supported when 'Method' value is 'hilbert'.");
    }
    if method != "hilbert" || bw {
        let nyq = fs / 2.0;
        if let Some([_, hi]) = flim {
            if hi > nyq {
                return rt(format!("Expected FrequencyLimits to be an array with all of the values <= {}.", if nyq.fract() == 0.0 { format!("{}", nyq as i64) } else { format!("{nyq}") }));
            }
        }
        let mut data: Vec<f64> = vec![];
        let mut tout: Vec<f64> = vec![];
        let mut nt = 0;
        for c in 0..ch {
            let col = Value::Num(match &x.im {
                Some(im) => Matrix::complex(vec![n, 1], x.data[c * n..(c + 1) * n].to_vec(), im[c * n..(c + 1) * n].to_vec()),
                None => Matrix::from_dims(vec![n, 1], x.data[c * n..(c + 1) * n].to_vec()),
            });
            let mut args = vec![col, if tm.numel() == 1 { scalar_value(tm.data[0]) } else { real_col(times.clone()) }, Value::string("spectrogram")];
            let r = it.call_function("pspectrum", args, 3)?;
            let pm = r[0].to_matrix()?;
            let fv = r[1].to_matrix()?.data;
            let tv = r[2].to_matrix()?.data;
            let (nf0, tt) = (pm.dims[0], pm.dims[1]);
            let (klo, khi) = match flim {
                Some([lo, hi]) => ((0..nf0).find(|&k| fv[k] >= lo).unwrap_or(0), (0..nf0).rev().find(|&k| fv[k] <= hi).unwrap_or(nf0 - 1)),
                None => (0, nf0 - 1),
            };
            let nf = khi + 1 - klo;
            let fv: Vec<f64> = fv[klo..=khi].to_vec();
            nt = tt;
            tout = tv;
            for k in 0..tt {
                let col = &pm.data[k * nf0 + klo..k * nf0 + khi + 1];
                let trapz = |g: &dyn Fn(usize) -> f64| -> f64 { (0..nf).map(|j| g(j)).sum() };
                let den = trapz(&|j| col[j]);
                let m1 = trapz(&|j| fv[j] * col[j]) / den;
                if bw {
                    let m2 = trapz(&|j| (fv[j] - m1) * (fv[j] - m1) * col[j]) / den;
                    data.push(scale * m2.sqrt());
                } else {
                    data.push(m1);
                }
            }
        }
        let mut f = Matrix::from_dims(vec![nt, ch], data);
        if x.kind == Kind::Single {
            f.kind = Kind::Single;
        }
        let mut tcol = Matrix::from_dims(vec![nt, 1], tout);
        if x.kind == Kind::Single {
            tcol.kind = Kind::Single;
        }
        return out(vec![Value::Num(f), Value::Num(tcol)].into_iter().take(nargout.max(1)).collect());
    }
    let mut data = vec![0.0; (n - 1) * ch];
    for c in 0..ch {
        let col = Matrix::from_dims(vec![n, 1], x.data[c * n..(c + 1) * n].to_vec());
        let xa = it.call_function("hilbert", vec![Value::Num(col)], 1)?.remove(0).to_matrix()?;
        let im = xa.im.clone().unwrap_or_else(|| vec![0.0; n]);
        let ph: Vec<f64> = (0..n).map(|k| im[k].atan2(xa.data[k])).collect();
        // Unwrap, then difference.
        let mut un = ph.clone();
        for k in 1..n {
            let mut d = ph[k] - ph[k - 1];
            d -= 2.0 * PI * (d / (2.0 * PI)).round();
            un[k] = un[k - 1] + d;
        }
        for k in 0..n - 1 {
            data[c * (n - 1) + k] = (un[k + 1] - un[k]) / (2.0 * PI) * fs;
        }
    }
    let t: Vec<f64> = (0..n - 1).map(|k| (times[k] + times[k + 1]) / 2.0).collect();
    let mut f = Matrix::from_dims(vec![n - 1, ch], data);
    if x.kind == Kind::Single {
        f.kind = Kind::Single;
    }
    let mut tm2 = Matrix::from_dims(vec![n - 1, 1], t);
    if x.kind == Kind::Single {
        tm2.kind = Kind::Single;
    }
    out(vec![Value::Num(f), Value::Num(tm2)].into_iter().take(nargout.max(1)).collect())
}
