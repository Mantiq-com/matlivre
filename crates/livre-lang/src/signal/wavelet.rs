//! Wavelet Toolbox transforms: `modwt` (maximal overlap discrete wavelet transform).
//!
//! Filter tables (`data/wfilters.txt`: `name|Lo_D|Hi_D|Lo_R|Hi_R`) were recorded from MATLAB's
//! `wfilters`.

use crate::builtins::Ret;
use crate::error::rt;
use crate::signal::poly::*;
use crate::value::*;
use std::collections::HashMap;
use std::sync::OnceLock;

pub const NAMES: &[&str] = &["modwt"];

static TABLE_TEXT: &str = include_str!("../../data/wfilters.txt");

type Table = HashMap<&'static str, [Vec<f64>; 4]>;

fn table() -> &'static Table {
    static T: OnceLock<Table> = OnceLock::new();
    T.get_or_init(|| {
        let mut m = HashMap::new();
        for line in TABLE_TEXT.lines() {
            let parts: Vec<&str> = line.split('|').collect();
            if parts.len() != 5 {
                continue;
            }
            let f = |s: &str| -> Vec<f64> { s.split(',').filter(|t| !t.is_empty()).map(|t| t.parse().unwrap_or(0.0)).collect() };
            m.insert(parts[0], [f(parts[1]), f(parts[2]), f(parts[3]), f(parts[4])]);
        }
        m
    })
}

/// Reconstruction filters `(Lo_R, Hi_R)` of an orthogonal wavelet.
fn orth_filters(name: &str) -> Result<(Vec<f64>, Vec<f64>), crate::error::LangError> {
    let key = if name.eq_ignore_ascii_case("haar") { "haar".to_string() } else { name.to_lowercase() };
    match table().get(key.as_str()) {
        Some(t) if !key.starts_with("bior") && !key.starts_with("rbio") => Ok((t[2].clone(), t[3].clone())),
        Some(_) => rt("Wavelet must be orthogonal."),
        None => if name.len() > 2 && name.to_lowercase().starts_with("db") && name[2..].parse::<u32>().is_ok() { rt("Invalid Input Argument.") } else { rt(format!("Invalid wavelet name: {name}.")) },
    }
}

fn is_text(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.kind == Kind::Char) || matches!(v, Value::Str(_))
}

fn wcol(v: &Value) -> Option<Matrix> {
    match v {
        Value::Num(m) if m.kind != Kind::Char => Some(m.clone()),
        _ => None,
    }
}

pub fn call(_it: &mut crate::interp::Interpreter, _name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() > 7 {
        return rt("Too many input arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        Value::Num(m) if m.kind == Kind::Char => return rt("Expected input to be one of these types:\n\ndouble, single"),
        other => return rt(format!("Expected input to be one of these types:\n\ndouble, single\n\nInstead its type was {}.", other.class_name())),
    };
    if x.dims.len() > 2 {
        return rt("Expected input to be two-dimensional.");
    }
    if x.data.iter().chain(x.im.iter().flatten()).any(|v| v.is_nan()) {
        return rt("Expected input to be non-NaN.");
    }
    if x.data.iter().chain(x.im.iter().flatten()).any(|v| v.is_infinite()) {
        return rt("Expected input to be finite.");
    }
    let single = x.kind == Kind::Single;
    let is_vec = x.dims[0] == 1 || x.dims[1] == 1;
    let (n0, ch) = if x.dims[0] == 1 && x.dims[1] >= 1 { (x.dims[1], 1) } else { (x.dims[0], x.dims[1]) };
    if n0 < 2 {
        return rt("Input must contain at least two samples.");
    }
    // Options.
    let mut rest: Vec<&Value> = Vec::new();
    let mut reflect = false;
    let mut align = false;
    let mut i = 1;
    let args = &a[1..];
    let mut idx: Vec<usize> = (0..args.len()).collect();
    // 'reflection'
    idx.retain(|&k| {
        let is = args[k].text().map(|t| t.eq_ignore_ascii_case("reflection")).unwrap_or(false);
        if is {
            reflect = true;
        }
        !is
    });
    // 'TimeAlign', value
    if let Some(p) = idx.iter().position(|&k| args[k].text().map(|t| !t.is_empty() && "timealign".starts_with(&t.to_lowercase())).unwrap_or(false)) {
        if p + 1 >= idx.len() {
            return rt("Index exceeds the number of array elements. Index must not exceed 1.");
        }
        let v = &args[idx[p + 1]];
        let m = match wcol(v) {
            Some(m) if m.numel() == 1 && m.im.is_none() && m.data[0].is_finite() => m,
            Some(_) => return rt("Expected TimeAlign to be a scalar."),
            None => return rt("Expected TimeAlign to be one of these types:

double, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64, logical"),
        };
        align = m.data[0] != 0.0;
        idx.drain(p..p + 2);
    }
    for &k in &idx {
        rest.push(&args[k]);
    }
    let _ = &mut i;
    let mut level = (n0 as f64).log2().floor() as i64;
    let mut wname = "sym4".to_string();
    let mut custom: Option<(Vec<f64>, Vec<f64>)> = None;
    if !rest.is_empty() {
        let nchar = rest.iter().filter(|v| is_text(v)).count();
        if nchar > 1 {
            return rt("Illegal to specify more than one wavelet name.");
        }
        if let Some(v) = rest.iter().find(|v| is_text(v)) {
            wname = v.text().unwrap_or_default();
        }
        let nums: Vec<&&Value> = rest.iter().filter(|v| !is_text(v)).collect();
        let scal: Vec<&&&Value> = nums.iter().filter(|v| matches!(v, Value::Num(m) if m.numel() == 1)).collect();
        if nums.len() > 3 {
            return rt("Expected at most 3 numeric input arguments in addition to the signal.");
        }
        if nums.len() == 3 && scal.is_empty() {
            return rt("A scalar level is required.");
        }
        if nums.len() == 1 && scal.is_empty() {
            return rt("Both Hi and Lo filters must be specified.");
        }
        if nums.len() > 1 {
            let lo = wcol(nums[0]).map(|m| m.data).unwrap_or_default();
            let hi = wcol(nums[1]).map(|m| m.data).unwrap_or_default();
            if lo.len() < 2 || hi.len() < 2 {
                return rt("Filters must contain at least two coefficients.");
            }
            if rest.iter().any(|v| is_text(v)) {
                return rt("You can not specify both a wavelet name and filters. Specify only a wavelet name or a filter pair.");
            }
            custom = Some((lo, hi));
        }
        if let Some(s) = scal.first() {
            let m = wcol(***s).unwrap();
            if m.data[0].fract() != 0.0 || !m.data[0].is_finite() {
                return rt("Expected LEV to be integer-valued.");
            }
            level = m.data[0] as i64;
        }
    }
    let jmax = (n0 as f64).log2().floor() as i64;
    if level <= 0 || level > jmax {
        return rt("Invalid transform level.".to_string());
    }
    let (lo_r, hi_r) = match custom {
        Some(c) => c,
        None => orth_filters(&wname)?,
    };
    let s2 = std::f64::consts::SQRT_2;
    let g: Vec<f64> = lo_r.iter().map(|v| v / s2).collect();
    let h: Vec<f64> = hi_r.iter().map(|v| v / s2).collect();
    let l = g.len();
    let j = level as usize;
    let cplx = x.im.is_some();
    let t_len = if reflect { 2 * n0 } else { n0 };
    let mut out_re = vec![0.0; (j + 1) * t_len * ch];
    let mut out_im = vec![0.0; (j + 1) * t_len * ch];
    for c in 0..ch {
        // Signal (reflected if requested), periodized when shorter than the filter.
        let mut sr: Vec<f64> = (0..n0).map(|k| x.data[c * n0 + k]).collect();
        let mut si: Vec<f64> = match &x.im {
            Some(im) => (0..n0).map(|k| im[c * n0 + k]).collect(),
            None => vec![0.0; n0],
        };
        if reflect {
            let (r2, i2) = (sr.clone(), si.clone());
            sr.extend(r2.iter().rev());
            si.extend(i2.iter().rev());
        }
        let siglen = sr.len();
        if siglen < l {
            let reps = l - siglen;
            let (r2, i2) = (sr.clone(), si.clone());
            for _ in 0..reps {
                sr.extend(&r2);
                si.extend(&i2);
            }
        }
        let nrep = sr.len();
        let (mut vr, mut vi) = (sr, si);
        for lev in 1..=j {
            let up = 1usize << (lev - 1);
            let step = |filt: &[f64], vr: &[f64], vi: &[f64]| -> (Vec<f64>, Vec<f64>) {
                let mut or = vec![0.0; nrep];
                let mut oi = vec![0.0; nrep];
                for t in 0..nrep {
                    let (mut ar, mut ai) = (0.0, 0.0);
                    for (k, f) in filt.iter().enumerate() {
                        let src = (t + nrep * (1 + (up * k) / nrep) - (up * k)) % nrep;
                        ar += f * vr[src];
                        ai += f * vi[src];
                    }
                    or[t] = ar;
                    oi[t] = ai;
                }
                (or, oi)
            };
            let (wr, wi) = step(&h, &vr, &vi);
            let (nvr, nvi) = step(&g, &vr, &vi);
            for t in 0..siglen {
                out_re[(c * t_len + t) * (j + 1) + lev - 1] = wr[t];
                out_im[(c * t_len + t) * (j + 1) + lev - 1] = wi[t];
            }
            vr = nvr;
            vi = nvi;
        }
        for t in 0..siglen {
            out_re[(c * t_len + t) * (j + 1) + j] = vr[t];
            out_im[(c * t_len + t) * (j + 1) + j] = vi[t];
        }
    }
    // out is stored (level fastest, then time, then channel) = MATLAB (J+1) x T x ch order.
    if align {
        let cen = |f: &[f64]| -> f64 {
            let num: f64 = f.iter().enumerate().fold(0.0, |s, (k, v)| s + k as f64 * (v * v));
            let nrm = f.iter().fold(0.0, |s, v| s + v * v).sqrt();
            num / (nrm * nrm)
        };
        let (gr, hr): (Vec<f64>, Vec<f64>) = (g.clone(), h.clone());
        let (es, ew) = (cen(&gr), cen(&hr));
        let mut shifts = Vec::new();
        for jj in 0..j {
            shifts.push(((2f64.powi(jj as i32) * (es + ew) - es).round()) as i64);
        }
        shifts.push((((2f64.powi((j + 1) as i32) - 1.0) * es).round()) as i64);
        let tl = t_len as i64;
        let (src_r, src_i) = (out_re.clone(), out_im.clone());
        for c in 0..ch {
            for t in 0..t_len {
                for (r, sh) in shifts.iter().enumerate() {
                    let from = ((t as i64 + sh) % tl + tl) % tl;
                    out_re[(c * t_len + t) * (j + 1) + r] = src_r[(c * t_len + from as usize) * (j + 1) + r];
                    out_im[(c * t_len + t) * (j + 1) + r] = src_i[(c * t_len + from as usize) * (j + 1) + r];
                }
            }
        }
    }
    // Rearrange to column-major (J+1) x T x ch: element (r, t, c) at r + (J+1)*(t + T*c): already so.
    let dims = if is_vec && ch == 1 { vec![j + 1, t_len] } else { vec![j + 1, t_len, ch] };
    let dims = if ch == 1 { vec![j + 1, t_len] } else { dims };
    let mut m = if cplx { Matrix::complex(dims, out_re, out_im) } else { Matrix::from_dims(dims, out_re) };
    if single {
        m = m.cast_to(Kind::Single);
    }
    let _ = nargout;
    out(vec![Value::Num(m)])
}
