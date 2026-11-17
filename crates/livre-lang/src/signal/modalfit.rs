//! `modalfit` (least-squares complex exponential fit of a single FRF).

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;
use std::f64::consts::PI;

pub const NAMES: &[&str] = &["modalfit", "modalsd"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if name == "modalsd" {
        return stabilization(it, a, nargout);
    }
    if a.len() < 4 {
        return rt("Not enough input arguments.");
    }
    if nargout > 4 {
        return rt("Too many output arguments.");
    }
    let frf = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected FRF to be one of these types:  double, single"),
    };
    let fm = a[1].to_matrix()?;
    let fsv = a[2].to_matrix()?;
    let mm = a[3].to_matrix()?;
    let p = frf.dims[0];
    let nfrf = frf.numel() / p.max(1);
    if fm.numel() != p {
        return rt(format!("Expected F to be an array with number of elements equal to {p}."));
    }
    if fsv.numel() != 1 || !(fsv.data[0] > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    let fs = fsv.data[0];
    if mm.numel() != 1 || !(mm.data[0] > 0.0) {
        return rt("Expected MNUM to be positive.");
    }
    if mm.data[0].fract() != 0.0 {
        return rt("Expected MNUM to be integer-valued.");
    }
    let mnum = mm.data[0] as usize;
    let mut phys: Option<Vec<f64>> = None;
    let mut i = 4;
    while i < a.len() {
        let key = a[i].text().unwrap_or_default();
        let Some(v) = a.get(i + 1) else { return rt(format!("No value was given for '{key}'. Name-value pair arguments require a name followed by a value.")) };
        match key.to_lowercase().as_str() {
            "fitmethod" => {
                let t = v.text().unwrap_or_default();
                match t.to_lowercase().as_str() {
                    "lsce" => {}
                    "pp" | "lsrf" => return rt(format!("modalfit: FitMethod '{t}' is not implemented.")),
                    _ => return rt(format!("Expected input to match one of these values:  'pp', 'lsce', 'lsrf'  The input, '{t}', did not match any of the valid values.")),
                }
            }
            "freqrange" => {
                let m = v.to_matrix()?;
                if m.numel() != 2 {
                    return rt("Expected FR to be an array with number of elements equal to 2.");
                }
                if !(m.data[1] > m.data[0]) {
                    return rt("Expected FR to be strictly increasing.");
                }
                return rt("modalfit: FreqRange is not implemented.");
            }
            "physfreq" => phys = Some(v.to_matrix()?.data),
            "driveindex" => {
                let m = v.to_matrix()?;
                if m.numel() == 2 && (m.data[0] > frf.dims.get(1).copied().unwrap_or(1) as f64 || m.data[1] > frf.dims.get(2).copied().unwrap_or(1) as f64) {
                    return rt("At least one element of DI is larger than the corresponding dimension of FRF.");
                }
            }
            "feedthrough" => {}
            _ => return rt(format!("'{key}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    if nfrf != 1 {
        return rt("modalfit: only a single FRF is supported.");
    }
    let (hre, him) = (frf.data.clone(), frf.im.clone().unwrap_or_else(|| vec![0.0; p]));
    // Impulse response from the one-sided FRF.
    let nfull = 2 * (p - 1);
    let mut fr = vec![0.0; nfull];
    let mut fi = vec![0.0; nfull];
    for k in 0..p {
        let (re, im) = if k == 0 || k == p - 1 { (hre[k], 0.0) } else { (hre[k], him[k]) };
        fr[k] = re;
        fi[k] = im;
        if k > 0 && k < p - 1 {
            fr[nfull - k] = re;
            fi[nfull - k] = -im;
        }
    }
    let spec = Value::Num(Matrix::complex(vec![nfull, 1], fr, fi));
    let hm = it.call_function("ifft", vec![spec], 1)?.remove(0).to_matrix()?;
    let h = hm.data.clone();
    let order = 2 * mnum;
    let nsamp = (22 * mnum).min(h.len());
    if nsamp <= order {
        return rt("modalfit: not enough samples for the requested number of modes.");
    }
    let rows = nsamp - order;
    let mut mat = vec![0.0; rows * order];
    let mut rhs = vec![0.0; rows];
    for r in 0..rows {
        let n = r + order;
        for c in 0..order {
            mat[c * rows + r] = h[n - 1 - c];
        }
        rhs[r] = h[n];
    }
    let coef = it
        .call_function("mldivide", vec![Value::Num(Matrix::from_dims(vec![rows, order], mat)), Value::Num(Matrix::from_dims(vec![rows, 1], rhs))], 1)?
        .remove(0)
        .to_matrix()?
        .data;
    let mut poly = vec![1.0];
    poly.extend(coef.iter().map(|c| -c));
    let z = it.call_function("roots", vec![real_col(poly)], 1)?.remove(0).to_matrix()?;
    let zi = z.im.clone().unwrap_or_else(|| vec![0.0; z.numel()]);
    let zc: Vec<(f64, f64)> = (0..z.numel()).map(|k| (z.data[k], zi[k])).collect();
    let sc: Vec<(f64, f64)> = zc.iter().map(|&(x, y)| ((x.hypot(y)).ln() * fs, y.atan2(x) * fs)).collect();
    let mut modes: Vec<usize> = (0..zc.len()).filter(|&k| zc[k].1 > 0.0).collect();
    modes.sort_by(|&x, &y| sc[x].1.partial_cmp(&sc[y].1).unwrap_or(std::cmp::Ordering::Equal));
    // Frequency domain residues of the fitted poles.
    let np = sc.len();
    let mut bre = vec![0.0; p * np];
    let mut bim = vec![0.0; p * np];
    let fv = &fm.data;
    for (k, &(sr, si)) in sc.iter().enumerate() {
        for j in 0..p {
            let w = 2.0 * PI * fv[j];
            // 1 / (jw - s)
            let (dr, di) = (-sr, w - si);
            let den = dr * dr + di * di;
            bre[k * p + j] = dr / den;
            bim[k * p + j] = -di / den;
        }
    }
    let bmat = Value::Num(Matrix::complex(vec![p, np], bre.clone(), bim.clone()));
    let hv = Value::Num(Matrix::complex(vec![p, 1], hre.clone(), him.clone()));
    let amp = it.call_function("mldivide", vec![bmat, hv], 1)?.remove(0).to_matrix()?;
    let ai = amp.im.clone().unwrap_or_else(|| vec![0.0; np]);
    let mut fn_: Vec<f64> = vec![];
    let mut dr_: Vec<f64> = vec![];
    let mut ms_re: Vec<f64> = vec![];
    let mut ms_im: Vec<f64> = vec![];
    for &k in &modes {
        let mag = sc[k].0.hypot(sc[k].1);
        fn_.push(mag / (2.0 * PI));
        dr_.push(-sc[k].0 / mag);
        // sqrt of the residue
        let (r, im) = (amp.data[k], ai[k]);
        let m = r.hypot(im).sqrt();
        let th = im.atan2(r) / 2.0;
        ms_re.push(m * th.cos());
        ms_im.push(m * th.sin());
    }
    while fn_.len() < mnum {
        fn_.push(f64::NAN);
        dr_.push(f64::NAN);
        ms_re.push(f64::NAN);
        ms_im.push(f64::NAN);
    }
    // PhysFreq: keep the mode nearest to each physical frequency.
    if let Some(pf) = phys {
        let nearest = |t: f64| -> usize {
            let mut best = 0;
            let mut bd = f64::INFINITY;
            for (k, v) in fn_.iter().enumerate() {
                if v.is_nan() {
                    continue;
                }
                let d = (v - t).abs();
                if d < bd {
                    bd = d;
                    best = k;
                }
            }
            best
        };
        let sel: Vec<usize> = pf.iter().map(|&t| nearest(t)).collect();
        fn_ = sel.iter().map(|&k| fn_[k]).collect();
        dr_ = sel.iter().map(|&k| dr_[k]).collect();
        ms_re = sel.iter().map(|&k| ms_re[k]).collect();
        ms_im = sel.iter().map(|&k| ms_im[k]).collect();
    }
    // Reconstructed FRF from the identified modes (poles with positive frequency and their conjugates).
    let mut ore = vec![0.0; p];
    let mut oim = vec![0.0; p];
    for &k in &modes {
        for j in 0..p {
            let (br, bi) = (bre[k * p + j], bim[k * p + j]);
            let (ar, aim) = (amp.data[k], ai[k]);
            ore[j] += ar * br - aim * bi;
            oim[j] += ar * bi + aim * br;
        }
    }
    for (k, &(x, y)) in zc.iter().enumerate() {
        let _ = (k, x, y);
    }
    let n = fn_.len();
    let single = frf.kind == Kind::Single;
    let mut fnm = Matrix::from_dims(vec![n, 1], fn_);
    let mut drm = Matrix::from_dims(vec![n, 1], dr_);
    let mut msm = Matrix::complex(vec![1, n], ms_re, ms_im);
    let mut om = Matrix::complex(vec![p, 1], ore, oim);
    if single {
        for q in [&mut fnm, &mut drm, &mut msm, &mut om] {
            q.kind = Kind::Single;
        }
    }
    out(vec![Value::Num(fnm), Value::Num(drm), Value::Num(msm), Value::Num(om)].into_iter().take(nargout.max(1)).collect())
}

/// `modalsd`: stable natural frequencies for model orders 1 to MaxModes.
pub fn stabilization(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    if nargout > 1 {
        return rt("Too many output arguments.");
    }
    let frf = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected FRF to be one of these types:  double, single"),
    };
    let fm = a[1].to_matrix()?;
    let fsv = a[2].to_matrix()?;
    let p = frf.dims[0];
    if fm.numel() != p {
        return rt(format!("Expected F to be an array with number of elements equal to {p}."));
    }
    if fsv.numel() != 1 || !(fsv.data[0] > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    if (a.len() - 3) % 2 != 0 {
        return rt("Optional name-value parameters must be in pairs.");
    }
    let mut maxm = 50.0;
    let mut sc = [0.01, 0.05];
    let mut i = 3;
    while i < a.len() {
        let key = a[i].text().unwrap_or_default();
        let v = &a[i + 1];
        match key.to_lowercase().as_str() {
            "fitmethod" => {
                let t = v.text().unwrap_or_default();
                match t.to_lowercase().as_str() {
                    "lsce" => {}
                    "lsrf" => return rt("modalsd: FitMethod 'lsrf' is not implemented."),
                    _ => return rt(format!("Expected input to match one of these values:  'lsce', 'lsrf'  The input, '{t}', did not match any of the valid values.")),
                }
            }
            "maxmodes" => {
                let m = v.to_matrix()?;
                if m.numel() != 1 || !(m.data[0] > 0.0) {
                    return rt("Expected MM to be positive.");
                }
                if m.data[0].fract() != 0.0 {
                    return rt("Expected MM to be integer-valued.");
                }
                maxm = m.data[0];
            }
            "scriteria" => {
                let m = v.to_matrix()?;
                if m.numel() != 2 {
                    return rt("Expected SC to be an array with number of elements equal to 2.");
                }
                sc = [m.data[0], m.data[1]];
            }
            "freqrange" => return rt("modalsd: FreqRange is not implemented."),
            _ => return rt(format!("'{key}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    let cap = (2 * (p - 1)) / 22;
    let mut mm = maxm as usize;
    if mm > cap {
        mm = cap.max(1);
        it.call_function("warning", vec![Value::string(&format!("Using a maximum model order of {mm} modes."))], 0)?;
    }
    let mut fns: Vec<Vec<f64>> = vec![];
    let mut drs: Vec<Vec<f64>> = vec![];
    for m in 1..=mm {
        let r = it.call_function("modalfit", vec![a[0].clone(), a[1].clone(), a[2].clone(), scalar_value(m as f64)], 2)?;
        fns.push(r[0].to_matrix()?.data);
        drs.push(r[1].to_matrix()?.data);
    }
    let mut res = vec![0.0; mm * mm];
    for m in 0..mm {
        for j in 0..=m {
            let f0 = fns[m][j];
            let mut val = f64::NAN;
            if m + 1 < mm && !f0.is_nan() {
                // nearest pole of the next order
                let mut best = usize::MAX;
                let mut bd = f64::INFINITY;
                for (k, &g) in fns[m + 1].iter().enumerate() {
                    if g.is_nan() {
                        continue;
                    }
                    let d = (g - f0).abs();
                    if d < bd {
                        bd = d;
                        best = k;
                    }
                }
                if best != usize::MAX {
                    let g = fns[m + 1][best];
                    let (d0, d1) = (drs[m][j], drs[m + 1][best]);
                    if (f0 - g).abs() / g < sc[0] && (d0 - d1).abs() < sc[1] {
                        val = f0;
                    }
                }
            }
            res[j * mm + m] = val;
        }
    }
    let mut mat = Matrix::from_dims(vec![mm, mm], res);
    if frf.kind == Kind::Single {
        mat.kind = Kind::Single;
    }
    out(vec![Value::Num(mat)])
}
