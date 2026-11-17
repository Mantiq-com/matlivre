//! `modalfrf`: frequency response functions of a mechanical structure from excitation and
//! response records (H1, H2 and Hv estimators, single input and single output).

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;
use std::f64::consts::PI;

pub const NAMES: &[&str] = &["modalfrf"];

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 4 {
        if a.is_empty() {
            return rt("Index exceeds the number of array elements. Index must not exceed 0.");
        }
        return rt("Not enough input arguments.");
    }
    if nargout > 3 {
        return rt("Too many output arguments.");
    }
    let signal = |v: &Value, name: &str| -> R<Matrix> {
        match v {
            Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => Ok(m.clone()),
            Value::Num(m) if m.kind == Kind::Char || m.kind == Kind::Logical => rt(format!("Expected {name} to be one of these types:  single, double")),
            other => rt(format!("Expected {name} to be one of these types:  single, double  Instead its type was {}.", other.class_name())),
        }
    };
    let x = signal(&a[0], "X")?;
    let y = signal(&a[1], "Y")?;
    let rows = |m: &Matrix| if m.dims[0] == 1 { m.dims[1] } else { m.dims[0] };
    let (nx, ny) = (rows(&x), rows(&y));
    if nx != ny {
        return rt("The number of rows of the two inputs must be equal to calculate the cross spectral density.");
    }
    let fsm = match &a[2] {
        Value::Num(m) if m.kind != Kind::Char => m.clone(),
        _ => return rt("Expected Fs to be positive."),
    };
    if fsm.numel() != 1 {
        return rt("Expected Fs to be a scalar.");
    }
    let fs = fsm.data[0];
    if !(fs > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    let wm = match &a[3] {
        Value::Num(m) if m.kind != Kind::Char => m.clone(),
        _ => return rt("Expected input to be numeric."),
    };
    let wlen = if wm.numel() == 1 { wm.data[0] as usize } else { wm.numel() };
    if wm.numel() == 1 && wm.data[0] < 1.0 {
        return rt("Expected NOVERLAP to be a scalar with value <= -1.");
    }
    if wm.numel() == 1 && wm.data[0].fract() != 0.0 {
        return rt("Window size must be larger than 1 sample.");
    }
    if wlen < 2 {
        return rt("Window size must be larger than 1 sample.");
    }
    // Optional overlap and name-value pairs.
    let mut i = 4;
    let mut noverlap = 0.0;
    if let Some(Value::Num(m)) = a.get(4) {
        if m.kind != Kind::Char {
            noverlap = m.data.first().copied().unwrap_or(0.0);
            if noverlap < 0.0 {
                return rt("Expected NOVERLAP to be nonnegative.");
            }
            if noverlap >= wlen as f64 {
                return rt(format!("Expected NOVERLAP to be a scalar with value <= {}.", wlen as i64 - 1));
            }
            i = 5;
        }
    }
    let mut sensor = "acc".to_string();
    let mut est = "h1".to_string();
    let mut measurement = "fixed".to_string();
    while i < a.len() {
        let key = match &a[i] {
            Value::Num(m) if m.kind == Kind::Char => a[i].text().unwrap_or_default(),
            _ => return rt("Expected a string scalar or character vector for the parameter name."),
        };
        let Some(v) = a.get(i + 1) else { return rt(format!("The argument '{key}' is a string scalar or character vector and does not match any parameter names. It failed validation for the argument 'noverlap'.")) };
        let val = v.text().unwrap_or_default();
        match key.to_lowercase().as_str() {
            "sensor" => match val.to_lowercase().as_str() {
                "dis" | "vel" | "acc" => sensor = val.to_lowercase(),
                _ => return rt(format!("Expected input to match one of these values:  'dis', 'vel', 'acc'  The input, '{val}', did not match any of the valid values.")),
            },
            "estimator" => match val.to_lowercase().as_str() {
                "h1" | "h2" | "hv" | "subspace" => est = val.to_lowercase(),
                _ => return rt(format!("Expected input to match one of these values:  'h1', 'h2', 'hv', 'subspace'  The input, '{val}', did not match any of the valid values.")),
            },
            "measurement" => match val.to_lowercase().as_str() {
                "rovinginput" | "rovingoutput" | "fixed" => measurement = val.to_lowercase(),
                _ => return rt(format!("Expected input to match one of these values:  'rovinginput', 'rovingoutput', 'fixed'  The input, '{val}', did not match any of the valid values.")),
            },
            _ => return rt(format!("The argument '{key}' is a string scalar or character vector and does not match any parameter names. It failed validation for the argument 'noverlap'.")),
        }
        i += 2;
    }
    if x.im.is_some() || y.im.is_some() {
        return rt("Expected X to be real.");
    }
    let cx = if x.dims[0] == 1 { 1 } else { x.dims[1] };
    let cy = if y.dims[0] == 1 { 1 } else { y.dims[1] };
    if est == "hv" && (cx != 1 || cy != 1) {
        return rt("The measurement must be single-input/single-output when using 'HV' for 'Estimator' value.");
    }
    if est == "h2" && cx != cy {
        return rt("X and Y must contain the same number of signals when using the H2 estimator.");
    }
    if wlen > nx {
        return rt("The length of the segments cannot be greater than the length of the input signal.");
    }
    if est == "subspace" {
        return rt("modalfrf: the subspace estimator is not implemented.");
    }
    let win = Value::Num(Matrix::from_dims(vec![wlen, 1], if wm.numel() == 1 { vec![1.0; wlen] } else { wm.data.clone() }));
    let nov = scalar_value(noverlap.floor());
    let nfft = scalar_value(wlen as f64);
    let fsv = scalar_value(fs);
    let nseg = (nx - noverlap.floor() as usize) / (wlen - noverlap.floor() as usize);
    let row_in = x.dims[0] == 1 && x.dims[1] > 1;
    let col = |m: &Matrix, c: usize| real_col(m.data[c * nx..(c + 1) * nx].to_vec());
    // Pairs of (input column, output column) and the shape of the result.
    let (pairs, shape): (Vec<(usize, usize)>, Vec<usize>) = if measurement == "rovinginput" {
        if cx != cy {
            return rt("X and Y must contain the same number of signals when using roving measurements.");
        }
        ((0..cx).map(|k| (k, k)).collect(), vec![0, 1, cx])
    } else if measurement == "rovingoutput" {
        if cx != cy {
            return rt("X and Y must contain the same number of signals when using roving measurements.");
        }
        ((0..cx).map(|k| (k, k)).collect(), vec![0, cx])
    } else if cx == 1 {
        ((0..cy).map(|k| (0, k)).collect(), vec![0, cy])
    } else {
        return rt("modalfrf: multiple fixed inputs are not implemented.");
    };
    let mut hs: Vec<(f64, f64)> = vec![];
    let mut cs: Vec<f64> = vec![];
    let mut f: Vec<f64> = vec![];
    let cplx = |m: &Matrix| -> Vec<(f64, f64)> {
        let im = m.im.clone().unwrap_or_else(|| vec![0.0; m.numel()]);
        (0..m.numel()).map(|k| (m.data[k], im[k])).collect()
    };
    for &(ci, co) in &pairs {
        let (xc, yc) = (col(&x, ci), col(&y, co));
        let args = |p: &Value, q: &Value| vec![p.clone(), q.clone(), win.clone(), nov.clone(), nfft.clone(), fsv.clone()];
        let h1 = it.call_function("tfestimate", args(&xc, &yc), 2)?;
        f = h1[1].to_matrix()?.data;
        let h1c = cplx(&h1[0].to_matrix()?);
        let mut h: Vec<(f64, f64)> = if est == "h1" {
            h1c.clone()
        } else {
            let sxx = it.call_function("cpsd", args(&xc, &xc), 1)?[0].to_matrix()?.data;
            let syy = it.call_function("cpsd", args(&yc, &yc), 1)?[0].to_matrix()?.data;
            let sxy = cplx(&it.call_function("cpsd", args(&xc, &yc), 1)?[0].to_matrix()?);
            (0..syy.len())
                .map(|k| {
                    let d = sxy[k];
                    let den = d.0 * d.0 + d.1 * d.1;
                    if est == "h2" {
                        // Syy / Syx
                        (syy[k] * d.0 / den, -syy[k] * d.1 / den)
                    } else {
                        // Hv = (Syy - Sxx + sqrt((Sxx - Syy)^2 + 4|Sxy|^2)) / (2 Syx)
                        let num = syy[k] - sxx[k] + ((sxx[k] - syy[k]).powi(2) + 4.0 * den).sqrt();
                        (num * d.0 / (2.0 * den), -num * d.1 / (2.0 * den))
                    }
                })
                .collect()
        };
        let coh = if nseg < 2 || row_in { vec![f64::NAN; h.len()] } else { it.call_function("mscohere", args(&xc, &yc), 1)?[0].to_matrix()?.data };
        for k in 0..h.len() {
            let w = 2.0 * PI * if k == 0 && f.len() > 1 { f[1] } else { f[k] };
            h[k] = match sensor.as_str() {
                "acc" => (-h[k].0 / (w * w), -h[k].1 / (w * w)),
                "vel" => (h[k].1 / w, -h[k].0 / w),
                _ => h[k],
            };
        }
        hs.extend(h);
        cs.extend(coh);
    }
    let nf = f.len();
    let np = pairs.len();
    let mut hm = Matrix::complex(vec![nf, np], hs.iter().map(|c| c.0).collect(), hs.iter().map(|c| c.1).collect());
    let mut cm = Matrix::from_dims(vec![nf, np], cs);
    let mut fm = Matrix::from_dims(vec![nf, 1], f);
    if shape.len() == 3 {
        hm.dims = vec![nf, 1, np];
    }
    if x.kind == Kind::Single || y.kind == Kind::Single {
        hm.kind = Kind::Single;
        cm.kind = Kind::Single;
        fm.kind = Kind::Single;
    }
    out(vec![Value::Num(hm), Value::Num(fm), Value::Num(cm)].into_iter().take(nargout.max(1)).collect())
}
