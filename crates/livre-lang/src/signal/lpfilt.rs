//! `lowpass`, `highpass`, `bandpass` and `bandstop`: zero-phase filtering with a minimum-order
//! Kaiser-window FIR filter designed from the passband edge and the steepness. (The IIR branch
//! of MATLAB's automatic choice needs `filtfilt` with second-order sections and is not implemented.)

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["lowpass", "highpass", "bandpass", "bandstop"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > 2 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) && m.im.is_none() => m.clone(),
        _ => return rt("Expected input number 1, X, to be real and one of these types:  double, single"),
    };
    let band = name == "bandpass" || name == "bandstop";
    if a.len() < 2 || matches!(&a[1], Value::Num(m) if m.kind == Kind::Char) {
        return rt("Not enough inputs. Passband frequency input must be specified.");
    }
    let wm = a[1].to_matrix()?;
    let want = if band { 2 } else { 1 };
    if wm.numel() != want {
        return rt(format!("Expected Wpass to be an array with number of elements equal to {want}."));
    }
    if band && !(wm.data[1] > wm.data[0]) {
        return rt("Expected Wpass to be strictly increasing.");
    }
    let mut i = 2;
    let mut fs: Option<f64> = None;
    if let Some(Value::Num(m)) = a.get(2) {
        if m.kind != Kind::Char {
            if m.numel() != 1 || !(m.data[0] > 0.0) {
                return rt("Expected Fs to be positive.");
            }
            fs = Some(m.data[0]);
            i = 3;
        }
    }
    let mut steep = 0.85;
    let mut att = 60.0;
    let mut resp = "auto".to_string();
    while i < a.len() {
        let key = a[i].text().unwrap_or_default();
        let v = a.get(i + 1);
        match (key.to_lowercase().as_str(), v) {
            ("steepness", Some(v)) => {
                let s = v.scalar().unwrap_or(f64::NAN);
                if !(s >= 0.5) {
                    return rt("Expected Steepness to be a scalar with value >= 0.5.");
                }
                if !(s < 1.0) {
                    return rt("Expected Steepness to be a scalar with value < 1.");
                }
                steep = s;
            }
            ("stopbandattenuation", Some(v)) => {
                let s = v.scalar().unwrap_or(f64::NAN);
                if !(s >= 10.0) {
                    return rt("Expected StopbandAttenuation to be a scalar with value >= 10.");
                }
                att = s;
            }
            ("impulseresponse", Some(v)) => {
                let t = v.text().unwrap_or_default().to_lowercase();
                if !["auto", "fir", "iir"].contains(&t.as_str()) {
                    return rt(format!("Expected input to match one of these values:  'auto', 'fir', 'iir'  The input, '{}', did not match any of the valid values.", v.text().unwrap_or_default()));
                }
                resp = t;
            }
            _ => return rt("Invalid input parameter name. Valid input names are: 'Steepness' and 'StopbandAttenuation'."),
        }
        i += 2;
    }
    let nyq = fs.map(|f| f / 2.0).unwrap_or(1.0);
    let w = &wm.data;
    if fs.is_none() && w.iter().any(|v| !(*v > 0.0 && *v < 1.0)) {
        return rt("Normalized passband frequency values must be in the (0,1) interval. Use a value in the (0,1) interval or specify a sample rate.");
    }
    if let Some(f) = fs {
        if w.iter().any(|v| !(*v > 0.0 && *v < f / 2.0)) {
            return rt("Passband frequency values must be in the (0, Fs/2) interval.");
        }
    }
    let c = 0.99 - 0.98 * steep;
    // Design arguments for designfilt.
    let (resp_name, args): (&str, Vec<(&str, f64)>) = match name {
        "lowpass" => ("lowpassfir", vec![("PassbandFrequency", w[0]), ("StopbandFrequency", w[0] + c * (nyq - w[0]))]),
        "highpass" => ("highpassfir", vec![("StopbandFrequency", w[0] - c * w[0]), ("PassbandFrequency", w[0])]),
        "bandpass" => (
            "bandpassfir",
            vec![("StopbandFrequency1", w[0] - c * w[0]), ("PassbandFrequency1", w[0]), ("PassbandFrequency2", w[1]), ("StopbandFrequency2", w[1] + c * (nyq - w[1]))],
        ),
        _ => (
            "bandstopfir",
            vec![("PassbandFrequency1", w[0]), ("StopbandFrequency1", w[0] + c * (w[1] - w[0]) / 2.0), ("StopbandFrequency2", w[1] - c * (w[1] - w[0]) / 2.0), ("PassbandFrequency2", w[1])],
        ),
    };
    let mut dargs: Vec<Value> = vec![Value::chars(resp_name)];
    for (k, v) in &args {
        dargs.push(Value::chars(k));
        dargs.push(scalar_value(*v));
    }
    let ripple_names: Vec<&str> = match name {
        "bandstop" => vec!["PassbandRipple1", "PassbandRipple2"],
        _ => vec!["PassbandRipple"],
    };
    for r in ripple_names {
        dargs.push(Value::chars(r));
        dargs.push(scalar_value(0.1));
    }
    let att_names: Vec<&str> = match name {
        "bandpass" => vec!["StopbandAttenuation1", "StopbandAttenuation2"],
        _ => vec!["StopbandAttenuation"],
    };
    for r in att_names {
        dargs.push(Value::chars(r));
        dargs.push(scalar_value(att));
    }
    dargs.push(Value::chars("DesignMethod"));
    dargs.push(Value::chars("kaiserwin"));
    dargs.push(Value::chars("MinOrder"));
    dargs.push(Value::chars("even"));
    if let Some(f) = fs {
        dargs.push(Value::chars("SampleRate"));
        dargs.push(scalar_value(f));
    }
    // The required order (the Kaiser estimate) decides between FIR and IIR.
    let dp = (10f64.powf(0.1 / 20.0) - 1.0) / (10f64.powf(0.1 / 20.0) + 1.0);
    let ds = 10f64.powf(-att / 20.0);
    let edges: Vec<f64> = args.iter().map(|(_, v)| *v / nyq).collect();
    let (kf, ka, kd): (Vec<f64>, Vec<f64>, Vec<f64>) = match name {
        "lowpass" => (edges.clone(), vec![1.0, 0.0], vec![dp, ds]),
        "highpass" => (edges.clone(), vec![0.0, 1.0], vec![ds, dp]),
        "bandpass" => (edges.clone(), vec![0.0, 1.0, 0.0], vec![ds, dp, ds]),
        _ => (edges.clone(), vec![1.0, 0.0, 1.0], vec![dp, ds, dp]),
    };
    let kr = it.call_function("kaiserord", vec![real_row(kf), real_row(ka), real_row(kd), scalar_value(2.0)], 1)?;
    let req = kr[0].scalar().unwrap_or(0.0);
    let n = if x.dims[0] == 1 || x.dims[1] == 1 { x.numel() } else { x.dims[0] };
    let use_fir = match resp.as_str() {
        "fir" => true,
        "iir" => false,
        _ => (n as f64) >= 2.0 * req,
    };
    if !use_fir {
        return rt("lowpass/highpass/bandpass/bandstop: the IIR branch (it needs filtfilt with second-order sections) is not implemented.");
    }
    if (n as f64) <= 2.0 * req && resp == "fir" {
        return rt(format!("Signal length must be greater than twice the required FIR filter order (RequiredOrder = {req}).\nRelax the steepness and stopband attenuation values to reduce the required filter order or increase the signal length."));
    }
    let d = it.call_function("designfilt", dargs, 1)?.remove(0);
    let coef = match &d {
        Value::Object(o) => o.props.borrow()[0].to_matrix()?,
        _ => return rt("designfilt failed."),
    };
    let order = coef.numel() - 1;
    let half = order / 2;
    let vector = x.dims[0] == 1 || x.dims[1] == 1;
    let (len, cols) = if vector { (x.numel(), 1) } else { (x.dims[0], x.dims[1]) };
    let mut res = vec![0.0; len * cols];
    for c in 0..cols {
        let mut col: Vec<f64> = x.data[c * len..(c + 1) * len].to_vec();
        col.extend(std::iter::repeat(0.0).take(half));
        let y = it.call_function("filter", vec![Value::Num(coef.clone()), scalar_value(1.0), real_col(col)], 1)?.remove(0).to_matrix()?;
        res[c * len..(c + 1) * len].copy_from_slice(&y.data[half..half + len]);
    }
    let mut ym = Matrix::from_dims(x.dims.clone(), res);
    if x.kind == Kind::Single {
        ym.kind = Kind::Single;
    }
    out(vec![Value::Num(ym), d].into_iter().take(nargout.max(1)).collect())
}
