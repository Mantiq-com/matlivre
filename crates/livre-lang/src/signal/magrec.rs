//! `stftmag2sig`: signal reconstruction from an STFT magnitude by the Griffin-Lim algorithm.
//!
//! Each iteration is `x = istft(Sest); S2 = stft(x); Sest = S .* exp(i*angle(S2))`; it stops after
//! `MaxIterations` or at the first iteration whose spectral convergence `||abs(S2) - S|| / ||S||` is
//! below `InconsistencyTolerance`, and the result is `istft(Sest)`.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["stftmag2sig"];

const OPTIONS: &[&str] = &["Window", "OverlapLength", "FrequencyRange", "Method", "InitializePhaseMethod", "InitialPhase", "MaxIterations", "InconsistencyTolerance", "InputTimeDimension"];

fn choice(v: &Value, name: &str, allowed: &[&str], shown: &str) -> Result<String, crate::error::LangError> {
    let t = v.text().unwrap_or_default();
    let low = t.to_lowercase();
    if !low.is_empty() {
        if let Some(m) = allowed.iter().find(|a| a.starts_with(&low)) {
            return Ok(m.to_string());
        }
    }
    rt(format!("Expected {name} to match one of these values:\n\n{shown}\n\nThe input, '{t}', did not match any of the valid values."))
}

fn scalar_of(v: &Value, name: &str) -> Result<f64, crate::error::LangError> {
    match v {
        Value::Num(m) if m.numel() == 1 => Ok(m.data[0]),
        Value::Num(_) => rt(format!("Expected {name} to be a scalar.")),
        _ => rt(format!("Expected {name} to be a scalar.")),
    }
}

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > 3 {
        return rt("Too many output arguments.");
    }
    let s = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        Value::Num(m) => return rt(format!("Expected S to be one of these types:\n\nsingle, double\n\nInstead its type was {}.", m.kind.name())),
        other => return rt(format!("Expected S to be one of these types:\n\nsingle, double\n\nInstead its type was {}.", other.class_name())),
    };
    if s.data.iter().any(|v| !v.is_finite()) {
        return rt("Expected S to be finite.");
    }
    if s.im.as_ref().is_some_and(|im| im.iter().any(|v| *v != 0.0)) || s.data.iter().any(|v| *v < 0.0) {
        return rt("Expected S to be nonnegative.");
    }
    let nfft = match a.get(1) {
        Some(Value::Num(m)) if m.numel() == 1 && m.kind != Kind::Char => m.data[0],
        Some(Value::Num(m)) if m.kind != Kind::Char => return rt("Expected FFTLength to be a scalar."),
        _ => return rt("FFT length must be provided as a numeric, integer, scalar."),
    };
    // Optional sample rate or sample time, then name-value pairs.
    let mut i = 2;
    let mut rate: Option<f64> = None;
    let mut tstep: Option<Value> = None;
    if let Some(v) = a.get(2) {
        match v {
            Value::Num(m) if m.kind != Kind::Char && m.numel() == 1 => {
                if !(m.data[0] > 0.0) {
                    return rt("Expected sample rate to be positive.");
                }
                rate = Some(m.data[0]);
                i = 3;
            }
            Value::Time(_) => {
                tstep = Some(v.clone());
                i = 3;
            }
            _ => {}
        }
    }
    let mut window: Option<Matrix> = None;
    let mut overlap: Option<f64> = None;
    let mut range = "centered".to_string();
    let mut method = "gla".to_string();
    let mut init_method: Option<String> = None;
    let mut init_phase: Option<Matrix> = None;
    let mut maxit = 100.0;
    let mut tol = 1e-4;
    let mut downrows = false;
    if (a.len() - i) % 2 != 0 {
        return rt("Name-value arguments must come in pairs.");
    }
    while i + 1 < a.len() {
        let key = a[i].text().unwrap_or_default();
        let v = &a[i + 1];
        let Some(name) = OPTIONS.iter().find(|o| o.to_lowercase() == key.to_lowercase()) else {
            return rt(format!("'{key}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function."));
        };
        match *name {
            "Window" => {
                let m = v.to_matrix()?;
                if m.dims.len() != 2 || (m.dims[0] != 1 && m.dims[1] != 1) {
                    return rt("Expected Window to be a vector.");
                }
                if m.numel() < 2 {
                    return rt("Expected WindowLength to be a scalar with value > 1.");
                }
                window = Some(m);
            }
            "OverlapLength" => overlap = Some(scalar_of(v, "OverlapLength")?),
            "FrequencyRange" => range = choice(v, "FrequencyRange", &["twosided", "onesided", "centered"], "'twosided', 'onesided', 'centered'")?,
            "Method" => method = choice(v, "Method", &["gla", "legla", "fgla", "gd"], "'gla', 'legla', 'fgla', 'gd'")?,
            "InitializePhaseMethod" => init_method = Some(choice(v, "InitializePhaseMethod", &["zeros", "random"], "'zeros', 'random'")?),
            "InitialPhase" => init_phase = Some(v.to_matrix()?),
            "MaxIterations" => {
                let m = scalar_of(v, "MaxIterations")?;
                if m < 0.0 {
                    return rt("Expected MaxIterations to be nonnegative.");
                }
                if m.fract() != 0.0 {
                    return rt("Expected MaxIterations to be integer-valued.");
                }
                maxit = m;
            }
            "InconsistencyTolerance" => {
                let m = scalar_of(v, "InconsistencyTolerance")?;
                if !(m > 0.0) {
                    return rt("Expected InconsistencyTolerance to be positive.");
                }
                tol = m;
            }
            _ => downrows = choice(v, "InputTimeDimension", &["acrosscolumns", "downrows"], "'acrosscolumns', 'downrows'")? == "downrows",
        }
        i += 2;
    }
    let win = match window {
        Some(w) => w,
        None => it.call_function("hann", vec![scalar_value(128.0), Value::chars("periodic")], 1)?.remove(0).to_matrix()?,
    };
    let wl = win.numel() as f64;
    if nfft.fract() != 0.0 {
        return rt("Expected FFTLength to be integer-valued.");
    }
    if !(nfft >= wl) {
        return rt(format!("Expected FFTLength to be a scalar with value >= {}.", wl));
    }
    let ov = overlap.unwrap_or((0.75 * wl).floor());
    if ov < 0.0 {
        return rt("Expected OverlapLength to be nonnegative.");
    }
    if ov.fract() != 0.0 {
        return rt("Expected OverlapLength to be integer-valued.");
    }
    if !(ov < wl) {
        return rt(format!("Expected OverlapLength to be a scalar with value < {}.", wl));
    }
    if method == "gd" {
        return rt("Undefined function 'dlarray' for input arguments of type 'double'.");
    }
    if method != "gla" || init_method.as_deref() == Some("random") {
        return rt(format!("stftmag2sig: the '{method}' method and random initial phases are not implemented."));
    }
    if init_method.is_some() && init_phase.is_some() {
        return rt("Specify only one of \"InitializePhaseMethod\" or \"InitialPhase\".");
    }
    // Orientation: frequency down the rows.
    let s = if downrows { transpose_first_two(&s) } else { s };
    let nf = nfft as usize;
    let rows_ok = if range == "onesided" { nf / 2 + 1 } else { nf };
    if s.dims[0] != rows_ok {
        return rt(format!("The number of rows in input S must be equal to the number of frequency samples, {rows_ok}."));
    }
    let phase = match init_phase {
        Some(p) => {
            let p = if downrows { transpose_first_two(&p) } else { p };
            if p.dims != s.dims {
                return rt("The input initial phase must have the same size as the input S.");
            }
            if p.data.iter().any(|v| v.abs() > std::f64::consts::PI) {
                return rt("The input initial phase must be within [-pi,pi].");
            }
            Some(p)
        }
        None => None,
    };
    let (f, t) = (s.dims[0], s.dims[1]);
    let chans = if s.dims.len() == 3 { s.dims[2] } else { 1 };
    let single = s.kind == Kind::Single;
    let mut xs: Vec<Vec<f64>> = Vec::new();
    let (mut flags, mut iters, mut incs) = (Vec::new(), Vec::new(), Vec::new());
    let (mut rec_re, mut rec_im) = (Vec::new(), Vec::new());
    for c in 0..chans {
        let lo = c * f * t;
        let mag: Vec<f64> = s.data[lo..lo + f * t].to_vec();
        let ph: Vec<f64> = match &phase {
            Some(p) => p.data[lo..lo + f * t].to_vec(),
            None => vec![0.0; f * t],
        };
        let r = gla(it, &mag, &ph, f, t, &win, ov, nfft, &range, maxit as usize, tol)?;
        xs.push(r.x);
        flags.push(r.flag);
        iters.push(r.iters as f64);
        incs.push(r.inc);
        rec_re.extend(r.re);
        rec_im.extend(r.im);
    }
    let len = xs[0].len();
    let mut data = Vec::with_capacity(len * chans);
    for o in &xs {
        data.extend_from_slice(o);
    }
    let mut m = Matrix::from_dims(vec![len, chans], data);
    if single {
        m.kind = Kind::Single;
        for v in m.data.iter_mut() {
            *v = *v as f32 as f64;
        }
    }
    let mut outs = vec![Value::Num(m)];
    if nargout >= 2 {
        let idx: Vec<f64> = (0..len).map(|k| k as f64).collect();
        let tv = if let Some(r) = rate {
            real_col(idx.iter().map(|k| k / r).collect())
        } else if let Some(ts) = &tstep {
            it.call_function("times", vec![real_col(idx), ts.clone()], 1)?.remove(0)
        } else {
            real_col(idx)
        };
        outs.push(tv);
    }
    if nargout >= 3 {
        let row = |v: Vec<f64>| if v.len() == 1 { scalar_value(v[0]) } else { real_row(v) };
        let sdims = if chans > 1 { vec![f, t, chans] } else { vec![f, t] };
        let phase_vals: Vec<f64> = rec_re.iter().zip(&rec_im).map(|(r, i)| i.atan2(*r)).collect();
        let info = Value::Struct(StructArray {
            dims: vec![1, 1],
            fields: ["ExitFlag", "NumIterations", "Inconsistency", "ReconstructedPhase", "ReconstructedSTFT"].iter().map(|n| n.to_string()).collect(),
            elems: vec![vec![
                row(flags),
                row(iters),
                row(incs),
                Value::Num(Matrix::from_dims(sdims.clone(), phase_vals)),
                Value::Num(Matrix::complex(sdims, rec_re, rec_im)),
            ]],
        });
        outs.push(info);
    }
    out(outs)
}

fn transpose_first_two(m: &Matrix) -> Matrix {
    let (r, c) = (m.dims[0], m.dims[1]);
    let pages = m.numel() / (r * c).max(1);
    let mut d = vec![0.0; m.numel()];
    for p in 0..pages {
        for i in 0..r {
            for j in 0..c {
                d[p * r * c + i * c + j] = m.data[p * r * c + j * r + i];
            }
        }
    }
    let mut dims = m.dims.clone();
    dims[0] = c;
    dims[1] = r;
    let mut out = Matrix::from_dims(dims, d);
    out.kind = m.kind;
    out
}

struct Rec {
    x: Vec<f64>,
    flag: f64,
    iters: usize,
    inc: f64,
    re: Vec<f64>,
    im: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
fn gla(it: &mut Interpreter, mag: &[f64], phase: &[f64], f: usize, t: usize, win: &Matrix, ov: f64, nfft: f64, range: &str, maxit: usize, tol: f64) -> Result<Rec, crate::error::LangError> {
    let sym = range != "onesided";
    let istft = |it: &mut Interpreter, re: &[f64], im: &[f64]| -> Result<Vec<f64>, crate::error::LangError> {
        let mut args = vec![
            Value::Num(Matrix::complex(vec![f, t], re.to_vec(), im.to_vec())),
            Value::chars("Window"),
            Value::Num(win.clone()),
            Value::chars("OverlapLength"),
            scalar_value(ov),
            Value::chars("FFTLength"),
            scalar_value(nfft),
            Value::chars("FrequencyRange"),
            Value::chars(range),
        ];
        if sym {
            args.push(Value::chars("ConjugateSymmetric"));
            args.push(Value::bool(true));
        }
        Ok(it.call_function("istft", args, 1)?.remove(0).to_matrix()?.data)
    };
    let stft = |it: &mut Interpreter, x: Vec<f64>| -> Result<Matrix, crate::error::LangError> {
        let args = vec![real_col(x), Value::chars("Window"), Value::Num(win.clone()), Value::chars("OverlapLength"), scalar_value(ov), Value::chars("FFTLength"), scalar_value(nfft), Value::chars("FrequencyRange"), Value::chars(range)];
        it.call_function("stft", args, 1)?.remove(0).to_matrix()
    };
    let norm_s: f64 = mag.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut re: Vec<f64> = mag.iter().zip(phase).map(|(m, p)| m * p.cos()).collect();
    let mut im: Vec<f64> = mag.iter().zip(phase).map(|(m, p)| m * p.sin()).collect();
    let (mut iters, mut flag, mut inc) = (0usize, 0.0, f64::INFINITY);
    for _ in 0..maxit {
        let x = istft(it, &re, &im)?;
        let s2 = stft(it, x)?;
        let mut conv = 0.0;
        for k in 0..mag.len() {
            let (r, i) = (s2.data[k], s2.im_at(k));
            let a = r.hypot(i);
            conv += (a - mag[k]) * (a - mag[k]);
            let ang = i.atan2(r);
            re[k] = mag[k] * ang.cos();
            im[k] = mag[k] * ang.sin();
        }
        iters += 1;
        inc = conv.sqrt() / norm_s;
        if inc < tol {
            flag = 1.0;
            break;
        }
    }
    let x = istft(it, &re, &im)?;
    Ok(Rec { x, flag, iters, inc, re, im })
}
