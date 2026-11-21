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

pub const NAMES: &[&str] = &["modwt", "modwpt"];

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

pub fn call(_it: &mut crate::interp::Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if name == "modwpt" {
        return modwpt(a, nargout);
    }
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

/// Circular convolution of `v` with `filt` dilated by `up` (indices taken modulo the length).
fn dil_conv(filt: &[f64], vr: &[f64], vi: &[f64], up: usize) -> (Vec<f64>, Vec<f64>) {
    let n = vr.len();
    let mut or = vec![0.0; n];
    let mut oi = vec![0.0; n];
    for t in 0..n {
        let (mut ar, mut ai) = (0.0, 0.0);
        for (k, f) in filt.iter().enumerate() {
            let back = (up * k) % n;
            let src = (t + n - back) % n;
            ar += f * vr[src];
            ai += f * vi[src];
        }
        or[t] = ar;
        oi[t] = ai;
    }
    (or, oi)
}

fn bit_reversal(level: usize, n: usize) -> usize {
    let mut seq = vec![0usize; level];
    let mut nn = n;
    for j in (1..=level).rev() {
        let r = nn % 4;
        seq[j - 1] = if r == 0 || r == 3 { 0 } else { 1 };
        nn /= 2;
    }
    seq.iter().enumerate().map(|(k, s)| s << k).sum()
}

fn center_of_energy(f: &[f64]) -> f64 {
    let num: f64 = f.iter().enumerate().fold(0.0, |s, (k, v)| s + k as f64 * (v * v));
    let nrm = f.iter().fold(0.0, |s, v| s + v * v).sqrt();
    num / (nrm * nrm)
}

pub fn modwpt(a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() > 8 {
        return rt("Too many input arguments.");
    }
    let args = &a[1..];
    let n_total = match &a[0] {
        Value::Num(m) => m.numel(),
        _ => 1,
    };
    let mut idx: Vec<usize> = (0..args.len()).collect();
    let mut fulltree = false;
    let mut align = false;
    for (pre, which) in [("fu", 0), ("ti", 1)] {
        let found = idx.iter().position(|&k| is_text(&args[k]) && args[k].text().map(|t| t.len() >= 2 && t.to_lowercase().starts_with(pre)).unwrap_or(false));
        if let Some(p) = found {
            if p + 1 >= idx.len() {
                return rt(format!("Index exceeds the number of array elements. Index must not exceed {}.", idx.len()));
            }
            let v = &args[idx[p + 1]];
            let b = match v {
                Value::Num(m) if m.im.is_none() && m.data.iter().all(|d| *d == 0.0 || *d == 1.0) && m.kind != Kind::Char => !m.data.is_empty() && m.data.iter().all(|d| *d == 1.0),
                _ => return rt("Expected true, false, 0, or 1. The input must be convertible to a logical."),
            };
            if which == 0 {
                fulltree = b;
            } else {
                align = b;
            }
            idx.drain(p..p + 2);
        }
    }
    let rest: Vec<&Value> = idx.iter().map(|&k| &args[k]).collect();
    let mut wname = "fk18".to_string();
    let nchar = rest.iter().filter(|v| is_text(v)).count();
    if nchar == 1 {
        wname = rest.iter().find(|v| is_text(v)).unwrap().text().unwrap_or_default();
    } else if nchar > 1 {
        return rt("There can only be one character vector input and that input must be the wavelet name.");
    }
    let nums: Vec<&&Value> = rest.iter().filter(|v| !is_text(v)).collect();
    if nums.len() > 3 {
        return rt("Expected at most 3 numeric input arguments in addition to the signal.");
    }
    let scal: Vec<&&Value> = rest.iter().filter(|v| matches!(v, Value::Num(m) if m.numel() == 1)).collect();
    if nums.len() == 3 && scal.is_empty() {
        return rt("A scalar level is required.");
    }
    let mut level_v: f64 = (n_total as f64).log2().floor().min(4.0);
    let mut custom: Option<(Vec<f64>, Vec<f64>)> = None;
    let (mut lo_m, mut hi_m): (Option<Matrix>, Option<Matrix>) = (None, None);
    let mut filt_single = false;
    if nums.len() > 1 {
        let lo = wcol(nums[0]);
        let hi = wcol(nums[1]);
        filt_single = lo.as_ref().map(|m| m.kind == Kind::Single).unwrap_or(false) || hi.as_ref().map(|m| m.kind == Kind::Single).unwrap_or(false);
        lo_m = lo.clone();
        hi_m = hi.clone();
        custom = Some((lo.map(|m| m.data).unwrap_or_default(), hi.map(|m| m.data).unwrap_or_default()));
        if nchar > 0 {
            return rt("You can not specify both a wavelet name and filters. Specify only a wavelet name or a filter pair.");
        }
    }
    if scal.len() > 1 && !matches!(scal[0], Value::Num(m) if m.kind == Kind::Char) {
        return rt("Too many input arguments.");
    }
    if let Some(s) = scal.first() {
        level_v = match s {
            Value::Num(m) => m.data[0],
            _ => 0.0,
        };
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        Value::Num(m) if m.kind == Kind::Char => return rt("Expected X to be one of these types:\n\ndouble, single"),
        other => return rt(format!("Expected X to be one of these types:\n\ndouble, single\n\nInstead its type was {}.", other.class_name())),
    };
    if x.data.iter().chain(x.im.iter().flatten()).any(|v| !v.is_finite()) {
        return rt("Expected X to be finite.");
    }
    if x.dims.len() > 2 || (x.dims[0] != 1 && x.dims[1] != 1) {
        return rt("Input data must be 1-D.");
    }
    let n = x.numel();
    if n < 2 {
        return rt("Input must contain at least two samples.");
    }
    if level_v.fract() != 0.0 {
        return rt("Expected LEVEL to be integer-valued.");
    }
    if level_v <= 0.0 {
        return rt("Expected LEVEL to be positive.");
    }
    if level_v > (n as f64).log2().floor() {
        return rt("Invalid transform level.");
    }
    let jl = level_v as usize;
    let single = x.kind == Kind::Single || filt_single;
    let s2 = std::f64::consts::SQRT_2;
    let (lo_r, hi_r) = match custom {
        Some(c) => c,
        None => orth_filters(&wname)?,
    };
    if let (Some(lm), Some(hm)) = (lo_m.as_ref(), hi_m.as_ref()) {
        if lm.im.is_some() {
            return rt("Expected Lo to be real.");
        }
        if hm.im.is_some() {
            return rt("Expected Hi to be real.");
        }
    }
    let lo: Vec<f64> = lo_r.iter().map(|v| v / s2).collect();
    let hi: Vec<f64> = hi_r.iter().map(|v| v / s2).collect();
    let cplx = x.im.is_some();
    let l2x: f64 = x.data.iter().map(|v| v * v).sum::<f64>() + x.im.as_ref().map(|im| im.iter().map(|v| v * v).sum::<f64>()).unwrap_or(0.0);
    let mut sr = x.data.clone();
    let mut si = x.im.clone().unwrap_or_else(|| vec![0.0; n]);
    if n < lo.len() {
        let (r2, i2) = (sr.clone(), si.clone());
        for _ in 0..(lo.len() - n) {
            sr.extend(&r2);
            si.extend(&i2);
        }
    }
    let mut rows: Vec<(Vec<f64>, Vec<f64>)> = vec![(Vec::new(), Vec::new()); (1 << (jl + 1)) - 2];
    let mut parents: Vec<(Vec<f64>, Vec<f64>)> = vec![(sr, si)];
    for kk in 1..=jl {
        let up = 1usize << (kk - 1);
        let base = (1usize << kk) - 2;
        let mut next = Vec::with_capacity(parents.len() * 2);
        for (nn, p) in parents.iter().enumerate() {
            let v = dil_conv(&lo, &p.0, &p.1, up);
            let w = dil_conv(&hi, &p.0, &p.1, up);
            let (first, second) = if nn % 2 == 1 { (w, v) } else { (v, w) };
            rows[base + 2 * nn] = first.clone();
            rows[base + 2 * nn + 1] = second.clone();
            next.push(first);
            next.push(second);
        }
        parents = next;
    }
    let total = rows.len();
    let first = if fulltree { 0 } else { total - (1 << jl) };
    let sel = total - first;
    let mut plev = Vec::new();
    let mut freq = Vec::new();
    for k in 1..=jl {
        let df = 1.0 / 2f64.powi(k as i32 + 1);
        let mut acc = 0.0;
        for nn in 0..(1usize << k) {
            let d = if nn == 0 { df * 0.5 } else { df };
            acc += d;
            plev.push(k as f64);
            freq.push(acc);
        }
    }
    let mut vals_r = vec![0.0; sel * n];
    let mut vals_i = vec![0.0; sel * n];
    for r in 0..sel {
        for t in 0..n {
            vals_r[r + sel * t] = rows[first + r].0[t];
            vals_i[r + sel * t] = rows[first + r].1[t];
        }
    }
    let ener: Vec<f64> = (0..sel).map(|r| (0..n).map(|t| vals_r[r + sel * t].powi(2) + vals_i[r + sel * t].powi(2)).sum()).collect();
    if align {
        let (es, ew) = (center_of_energy(&lo), center_of_energy(&hi));
        let mut shifts: Vec<i64> = Vec::new();
        let lv_list: Vec<usize> = if fulltree { (1..=jl).collect() } else { vec![jl] };
        for jv in lv_list {
            for nn in 0..(1usize << jv) {
                let high = bit_reversal(jv, nn) as f64;
                let low = (1usize << jv) as f64 - 1.0 - high;
                shifts.push((low * es + high * ew).round() as i64);
            }
        }
        let (sr_, si_) = (vals_r.clone(), vals_i.clone());
        for r in 0..sel {
            for t in 0..n {
                let from = (((t as i64 + shifts[r]) % n as i64) + n as i64) % n as i64;
                vals_r[r + sel * t] = sr_[r + sel * from as usize];
                vals_i[r + sel * t] = si_[r + sel * from as usize];
            }
        }
    }
    let mk = |dims: Vec<usize>, re: Vec<f64>, im: Option<Vec<f64>>| -> Value {
        let m = match im {
            Some(im) => Matrix::complex(dims, re, im),
            None => Matrix::from_dims(dims, re),
        };
        Value::Num(if single { m.cast_to(Kind::Single) } else { m })
    };
    let w = mk(vec![sel, n], vals_r, if cplx { Some(vals_i) } else { None });
    let pl = mk(vec![sel, 1], plev[first..].to_vec(), None);
    let fv = mk(vec![sel, 1], freq[first..].to_vec(), None);
    let ev = mk(vec![sel, 1], ener.clone(), None);
    let rv = mk(vec![sel, 1], ener.iter().map(|e| e / l2x).collect(), None);
    let mut outs = vec![w, pl, fv, ev, rv];
    outs.truncate(nargout.max(1));
    out(outs)
}
