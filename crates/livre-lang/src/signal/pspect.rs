//! `pspectrum`: power spectrum (modified periodogram or Welch with a Kaiser window).

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;
use std::f64::consts::PI;

pub const NAMES: &[&str] = &["pspectrum"];

fn enbw(w: &[f64]) -> f64 {
    let s: f64 = w.iter().sum();
    w.len() as f64 * w.iter().map(|v| v * v).sum::<f64>() / (s * s)
}

/// Equivalent noise bandwidth (bins) of the Kaiser window for Leakage 0, 0.01, ..., 1 as pspectrum reports it
/// (read from the frequency resolution limits it enforces).
const ENBW_TABLE: [f64; 101] = [3.600207845317627, 3.582460531672408, 3.564623758992623, 3.546696165795149, 3.528676415913742, 3.510563184824802, 3.492355147210171, 3.474050965693428, 3.4556492806862886, 3.437148701281589, 3.418547797129393, 3.3998450912327773, 3.381039053599825, 3.362128095688307, 3.343110565579682, 3.323984743818837, 3.3047488398561877, 3.2854009890286298, 3.265939250015875, 3.246361602708716, 3.226665946425745, 3.206850098415062, 3.186911792577499, 3.166848678347898, 3.146658319670983, 3.126338194008333, 3.105885691313015, 3.085298112908397, 3.064572670207672, 3.0437064832106326, 3.022696578714221, 3.001539888173387, 2.980233245148798, 2.9587733822779168, 2.937156927705988, 2.915380400913468, 2.893440207876425, 2.871332635496439, 2.849053845236555, 2.826599865899778, 2.8039665854867, 2.781149742068737, 2.758144913613544, 2.7349475066991302, 2.7115527440531904, 2.687955650854219, 2.664151039730908, 2.640133494396369, 2.615897351853727, 2.591436683109594, 2.56674527233198, 2.5418165943891484, 2.5166437907059658, 2.491219643374278, 2.4655365474538313, 2.439586481400285, 2.413360975556839, 2.386851078646015, 2.360047322198123, 2.332939682852931, 2.3055175424710987, 2.277769645991879, 2.2496840569736363, 2.221248110753706, 2.1924483651641404, 2.1632705487398463, 2.133699506355681, 2.103719142229012, 2.0733123602242878, 2.042461001396142, 2.0111457787075677, 1.979346208859706, 1.947040541169751, 1.914205683433541, 1.880817124709333, 1.846848854959321, 1.812273281485415, 1.7770611420958131, 1.741181414938905, 1.704601224941042, 1.667285746784687, 1.62919810436351, 1.590299266650918, 1.550547939918597, 1.509900456241564, 1.468310658226277, 1.425729779898334, 1.382106323686288, 1.3378754839965419, 1.292918387254158, 1.247610679000242, 1.202491339672078, 1.1584220491569601, 1.116679714436776, 1.079003082863279, 1.0474064400133969, 1.023771807086518, 1.008965082647308, 1.0020270625689351, 1.000143068739822, 1.00000397846243];

fn enbw_leak(leak: f64) -> f64 {
    let x = (leak * 100.0).clamp(0.0, 100.0);
    let i = (x.floor() as usize).min(99);
    let t = x - i as f64;
    let p = |k: i64| ENBW_TABLE[k.clamp(0, 100) as usize];
    let (p0, p1, p2, p3) = (p(i as i64 - 1), p(i as i64), p(i as i64 + 1), p(i as i64 + 2));
    0.5 * ((2.0 * p1) + (-p0 + p2) * t + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t * t * t)
}

fn kaiser(it: &mut Interpreter, l: usize, beta: f64) -> R<Vec<f64>> {
    Ok(it.call_function("kaiser", vec![scalar_value(l as f64), scalar_value(beta)], 1)?.remove(0).to_matrix()?.data)
}

fn fmt_e(v: f64) -> String {
    // MATLAB's %.15e with width.
    let s = format!("{:.15e}", v);
    let (m, e) = s.split_once('e').unwrap();
    let ev: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if ev < 0 { '-' } else { '+' }, ev.abs())
}

fn fmt_e6(v: f64) -> String {
    let s = format!("{:.6e}", v);
    let (m, e) = s.split_once('e').unwrap();
    let ev: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if ev < 0 { '-' } else { '+' }, ev.abs())
}

fn fmt_f(v: f64) -> String {
    format!("{:.6}", v)
}

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > 3 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected input to be one of these types:  double, single"),
    };
    let mut i = 1;
    let mut fs = 2.0 * PI;
    let mut has_fs = false;
    let mut times: Option<Vec<f64>> = None;
    if let Some(Value::Num(m)) = a.get(1) {
        if m.kind != Kind::Char {
            if m.numel() == 1 {
                fs = m.data[0];
                if !(fs > 0.0) {
                    return rt("Expected sample rate to be positive.");
                }
                has_fs = true;
            } else {
                times = Some(m.data.clone());
            }
            i = 2;
        }
    }
    let n = if x.dims[0] == 1 && x.dims.len() == 2 { x.dims[1] } else { x.dims[0] };
    let cols = if x.dims[0] == 1 && x.dims.len() == 2 { 1 } else { x.numel() / n.max(1) };
    if let Some(t) = &times {
        if t.len() != n {
            return rt("Expected time vector to have the same number of samples as the signal.");
        }
        let dt = (t[n - 1] - t[0]) / (n as f64 - 1.0);
        fs = 1.0 / dt;
        has_fs = true;
    }
    let mut is_spec = false;
    if i < a.len() {
        let first = a[i].text().unwrap_or_default().to_lowercase();
        if matches!(first.as_str(), "power") {
            i += 1;
        } else if first == "spectrogram" {
            is_spec = true;
            i += 1;
        } else if first == "persistence" {
            return rt("pspectrum: type 'persistence' is not implemented.");
        }
    }
    if i < a.len() {
        let t = a[i].text().unwrap_or_default().to_lowercase();
        let known = ["frequencylimits", "leakage", "minthreshold", "frequencyresolution", "reassign", "timeresolution", "overlappercent", "numpowerbins", "twosided"];
        if !known.contains(&t.as_str()) {
            let orig = a[i].text().unwrap_or_default();
            return rt(format!("Expected input to match one of these values:  'power', 'spectrogram', 'persistence', 'FrequencyLimits', 'Leakage', 'MinThreshold', 'FrequencyResolution', 'Reassign', 'TimeResolution', 'OverlapPercent', 'NumPowerBins', 'TwoSided'  The input, '{orig}', did not match any of the valid values."));
        }
    }
    if (a.len() - i) % 2 != 0 {
        return rt("Name-value inputs must come in pairs.");
    }
    let mut leak = 0.5;
    let mut fres: Option<f64> = None;
    let mut flim: Option<[f64; 2]> = None;
    let mut minthr = f64::NEG_INFINITY;
    let mut tres: Option<f64> = None;
    let mut ovp: Option<f64> = None;
    let cplx = x.im.is_some();
    let mut two = cplx;
    while i < a.len() {
        let key = a[i].text().unwrap_or_default();
        let v = &a[i + 1];
        let num = |v: &Value| v.to_matrix().ok();
        if matches!(key.to_lowercase().as_str(), "leakage" | "frequencyresolution" | "frequencylimits" | "minthreshold") {
            if let Value::Num(m) = v {
                if m.kind == Kind::Char {
                    let t = v.text().unwrap_or_default();
                    return rt(format!("Expected input to match one of these values:  'power', 'spectrogram', 'persistence', 'FrequencyLimits', 'Leakage', 'MinThreshold', 'FrequencyResolution', 'Reassign', 'TimeResolution', 'OverlapPercent', 'NumPowerBins', 'TwoSided'  The input, '{t}', did not match any of the valid values."));
                }
            }
        }
        match key.to_lowercase().as_str() {
            "leakage" => {
                let Some(m) = num(v) else { return rt("Expected Leakage to be numeric.") };
                if m.numel() != 1 {
                    return rt("Expected Leakage to be a scalar.");
                }
                if m.data[0] > 1.0 {
                    return rt("Expected Leakage to be a scalar with value <= 1.");
                }
                if m.data[0] < 0.0 {
                    return rt("Expected Leakage to be a scalar with value >= 0.");
                }
                leak = m.data[0];
            }
            "frequencyresolution" => {
                let m = num(v).unwrap_or_else(|| Matrix::from_dims(vec![0, 0], vec![]));
                if m.numel() != 1 {
                    return rt("Expected FrequencyResolution to be a scalar.");
                }
                if !(m.data[0] > 0.0) {
                    return rt("Expected FrequencyResolution to be positive.");
                }
                fres = Some(m.data[0]);
            }
            "frequencylimits" => {
                let m = num(v).unwrap_or_else(|| Matrix::from_dims(vec![0, 0], vec![]));
                if m.numel() != 2 {
                    return rt("Expected FrequencyLimits to be an array with number of elements equal to 2.");
                }
                if !(m.data[1] > m.data[0]) {
                    return rt("Expected FrequencyLimits to be strictly increasing.");
                }
                flim = Some([m.data[0], m.data[1]]);
            }
            "minthreshold" => minthr = num(v).and_then(|m| m.data.first().copied()).unwrap_or(f64::NEG_INFINITY),
            "twosided" => two = num(v).map(|m| m.data.first().copied().unwrap_or(0.0) != 0.0).unwrap_or(false),
            "timeresolution" if is_spec => {
                let m = num(v).unwrap_or_else(|| Matrix::from_dims(vec![0, 0], vec![]));
                if m.numel() != 1 {
                    return rt("Expected TimeResolution to be a scalar.");
                }
                if !(m.data[0] > 0.0) {
                    return rt("Expected TimeResolution to be positive.");
                }
                tres = Some(m.data[0]);
            }
            "overlappercent" if is_spec => {
                let m = num(v).unwrap_or_else(|| Matrix::from_dims(vec![0, 0], vec![]));
                if m.numel() != 1 {
                    return rt("Expected OverlapPercent to be a scalar.");
                }
                if m.data[0] >= 100.0 {
                    return rt("Expected OverlapPercent to be a scalar with value < 100.");
                }
                if m.data[0] < 0.0 {
                    return rt("Expected OverlapPercent to be a scalar with value >= 0.");
                }
                ovp = Some(m.data[0]);
            }
            "overlappercent" | "timeresolution" | "numpowerbins" | "reassign" => {
                return rt(format!("'{key}' is invalid when TYPE is 'power'."));
            }
            _ => {
                return rt(format!("Expected input to match one of these values:  'power', 'spectrogram', 'persistence', 'FrequencyLimits', 'Leakage', 'MinThreshold', 'FrequencyResolution', 'Reassign', 'TimeResolution', 'OverlapPercent', 'NumPowerBins', 'TwoSided'  The input, '{key}', did not match any of the valid values."));
            }
        }
        i += 2;
    }
    if cplx {
        two = true;
    }
    let fnyq = fs / 2.0;
    let (lo, hi) = if two { (-fnyq, fnyq) } else { (0.0, fnyq) };
    let (f1, f2) = match flim {
        Some([p, q]) => {
            if q < lo || p > hi {
                return rt(format!("Frequency limits are outside the Nyquist range. Specify limits within the [  {},   {}] interval.", fmt_e(lo), fmt_e(hi)));
            }
            (p.max(lo), q.min(hi))
        }
        None => (lo, hi),
    };

    if is_spec {
        if cols > 1 {
            return rt("Multi-channel inputs are not supported when TYPE is 'spectrogram' or 'persistence'.");
        }
        if fres.is_some() && tres.is_some() {
            return rt("'FrequencyResolution' and 'TimeResolution' values are specified simultaneously. Specify only one.");
        }
        let beta = 40.0 * (1.0 - leak);
        let en = enbw_leak(leak);
        let (npts, lseg) = match (fres, tres) {
            (Some(r), _) => {
                let rmin = fs * en / n as f64;
                let rmax = fs * en;
                if r < rmin * (1.0 - 1e-12) || r > rmax * (1.0 + 1e-12) {
                    return rt(format!("'FrequencyResolution' value is not achievable for the specified parameters. 'FrequencyResolution' value must be within [{}, {}] (Hz).", fmt_f(rmin), fmt_f(rmax)));
                }
                let npts = (8.0 * fnyq / r).round() as usize + 1;
                let reff = 8.0 * fnyq / (npts as f64 - 1.0);
                (npts, (fs * en / reff).floor() as usize)
            }
            (None, Some(tr)) => {
                let dur = n as f64 / fs;
                if tr >= dur * (1.0 + 1e-12) || tr < (1.0 / fs) * (1.0 - 1e-12) {
                    return rt(format!("'TimeResolution' value must be smaller than the signal duration and greater than one sample duration. Specify a 'TimeResolution' value smaller than {:>14} and greater than {:>14} seconds.", fmt_e6(dur), fmt_e6(1.0 / fs)));
                }
                (1024, (tr * fs).round() as usize)
            }
            _ => {
                let d = if n < 64 { 2 } else if n < 2048 { 8 } else if n < 4096 { 16 } else if n < 8192 { 32 } else if n < 16384 { 64 } else { 128 };
                (1024, (n + d - 1) / d)
            }
        };
        let lseg = lseg.clamp(1, n);
        let w = kaiser(it, lseg, beta)?;
        let ov = ovp.unwrap_or((1.0 - 1.0 / (2.0 * en - 1.0)) * 100.0);
        let hop = (((lseg as f64) * (1.0 - ov / 100.0)) - 1e-9).ceil().max(1.0) as usize;
        let nseg = if n <= lseg { 1 } else { ((n - lseg) as f64 / hop as f64).ceil() as usize + 1 };
        let grid: Vec<f64> = (0..npts).map(|k| f1 + (f2 - f1) * k as f64 / (npts as f64 - 1.0)).collect();
        let full_span = (f1 - lo).abs() < 1e-12 && (f2 - hi).abs() < 1e-12;
        let m_fold = if two { npts - 1 } else { 2 * (npts - 1) };
        let re = &x.data[..n];
        let im: Vec<f64> = match &x.im {
            Some(v) => v[..n].to_vec(),
            None => vec![0.0; n],
        };
        let mut pm = vec![0.0; npts * nseg];
        let mut tv = vec![0.0; nseg];
        for sg in 0..nseg {
            let st = sg * hop;
            let mut xr = vec![0.0; npts];
            let mut xi = vec![0.0; npts];
            let len = lseg.min(n.saturating_sub(st));
            let wl = if len < lseg { kaiser(it, len, beta)? } else { w.clone() };
            let sw: f64 = wl.iter().sum();
            if full_span {
                let mut fr = vec![0.0; m_fold];
                let mut fi = vec![0.0; m_fold];
                for j in 0..len {
                    let ph = -2.0 * PI * f1 * j as f64 / fs;
                    let (sn, co) = ph.sin_cos();
                    let (ar, ai) = (re[st + j] * wl[j], im[st + j] * wl[j]);
                    fr[j % m_fold] += ar * co - ai * sn;
                    fi[j % m_fold] += ar * sn + ai * co;
                }
                let spec = it.call_function("fft", vec![Value::Num(Matrix::complex(vec![m_fold, 1], fr, fi))], 1)?.remove(0).to_matrix()?;
                let si = spec.im.clone().unwrap_or_else(|| vec![0.0; m_fold]);
                for k in 0..npts {
                    xr[k] = spec.data[k % m_fold];
                    xi[k] = si[k % m_fold];
                }
            } else {
                for (k, f) in grid.iter().enumerate() {
                    let (mut sr, mut si) = (0.0, 0.0);
                    let step = -2.0 * PI * f / fs;
                    for j in 0..len {
                        let (sn, co) = (step * j as f64).sin_cos();
                        let (ar, ai) = (re[st + j] * wl[j], im[st + j] * wl[j]);
                        sr += ar * co - ai * sn;
                        si += ar * sn + ai * co;
                    }
                    xr[k] = sr;
                    xi[k] = si;
                }
            }
            for k in 0..npts {
                let mut v = (xr[k] * xr[k] + xi[k] * xi[k]) / (sw * sw);
                if !two {
                    v *= 2.0;
                }
                if minthr.is_finite() && 10.0 * v.log10() <= minthr {
                    v = 0.0;
                }
                pm[sg * npts + k] = v;
            }
            tv[sg] = (st as f64 + lseg as f64 / 2.0) / if has_fs { fs } else { 1.0 };
        }
        let mut pmat = Matrix::from_dims(vec![npts, nseg], pm);
        let mut fmat = Matrix::from_dims(vec![npts, 1], grid);
        let mut tmat = Matrix::from_dims(vec![nseg, 1], tv);
        if x.kind == Kind::Single {
            pmat.kind = Kind::Single;
            fmat.kind = Kind::Single;
            tmat.kind = Kind::Single;
        }
        return out(vec![Value::Num(pmat), Value::Num(fmat), Value::Num(tmat)].into_iter().take(nargout.max(1)).collect());
    }
    let beta = 40.0 * (1.0 - leak);
    let wfull = kaiser(it, n, beta)?;
    let en = enbw_leak(leak);
    let rmin = fs * en / n as f64;
    let rmax = fs * en;
    let span = f2 - f1;
    let (rbw, npts) = match fres {
        Some(r) => {
            if r < rmin * (1.0 - 1e-12) || r > rmax * (1.0 + 1e-12) {
                return rt(format!("'FrequencyResolution' value is not achievable for the specified parameters. 'FrequencyResolution' value must be within [{}, {}] (Hz).", fmt_f(rmin), fmt_f(rmax)));
            }
            (r, (8.0 * fnyq / r).round() as usize + 1)
        }
        None => (rmin.max(fs / 1024.0), 4096),
    };
    let npts = npts.max(2);
    let grid: Vec<f64> = (0..npts).map(|k| f1 + (f2 - f1) * k as f64 / (npts as f64 - 1.0)).collect();
    // Segmentation.
    let welch = rbw > rmin * (1.0 + 1e-9);
    let (l, starts, w) = if welch {
        let rbw_eff = if fres.is_some() { 8.0 * fnyq / (npts as f64 - 1.0) } else { rbw };
        let l = ((fs * en / rbw_eff).round() as usize).clamp(2, n);
        let w = kaiser(it, l, beta)?;
        if l >= n {
            (n, vec![0usize], wfull.clone())
        } else {
            let stride = l as f64 / (2.0 * en - 1.0);
            let k = ((n - l) as f64 / stride).ceil() as usize + 1;
            let st: Vec<usize> = if k <= 1 { vec![0] } else { (0..k).map(|j| ((n - l) as f64 * j as f64 / (k as f64 - 1.0) + 0.5).floor() as usize).collect() };
            (l, st, w)
        }
    } else {
        (n, vec![0usize], wfull.clone())
    };
    let sw: f64 = w.iter().sum();
    let full_span = (f1 - lo).abs() < 1e-12 && (f2 - hi).abs() < 1e-12;
    let m_fold = if two { npts - 1 } else { 2 * (npts - 1) };
    let out_cols = cols;
    let mut pall = vec![0.0; npts * out_cols];
    for c in 0..cols {
        let re = &x.data[c * n..(c + 1) * n];
        let im: Vec<f64> = match &x.im {
            Some(v) => v[c * n..(c + 1) * n].to_vec(),
            None => vec![0.0; n],
        };
        let mut acc = vec![0.0; npts];
        for &st in &starts {
            let mut xr = vec![0.0; npts];
            let mut xi = vec![0.0; npts];
            if full_span {
                let mut fr = vec![0.0; m_fold];
                let mut fi = vec![0.0; m_fold];
                for j in 0..l {
                    // Shift of the grid origin: f1 * j / fs turns.
                    let ph = -2.0 * PI * f1 * j as f64 / fs;
                    let (s, co) = ph.sin_cos();
                    let (ar, ai) = (re[st + j] * w[j], im[st + j] * w[j]);
                    fr[j % m_fold] += ar * co - ai * s;
                    fi[j % m_fold] += ar * s + ai * co;
                }
                let spec = it.call_function("fft", vec![Value::Num(Matrix::complex(vec![m_fold, 1], fr, fi))], 1)?.remove(0).to_matrix()?;
                let si = spec.im.clone().unwrap_or_else(|| vec![0.0; m_fold]);
                for k in 0..npts {
                    xr[k] = spec.data[k % m_fold];
                    xi[k] = si[k % m_fold];
                }
            } else {
                for (k, f) in grid.iter().enumerate() {
                    let (mut sr, mut si) = (0.0, 0.0);
                    let step = -2.0 * PI * f / fs;
                    let (ds, dc) = step.sin_cos();
                    let (mut cr, mut ci) = (1.0, 0.0);
                    for j in 0..l {
                        let (ar, ai) = (re[st + j] * w[j], im[st + j] * w[j]);
                        sr += ar * cr - ai * ci;
                        si += ar * ci + ai * cr;
                        let nr = cr * dc - ci * ds;
                        ci = cr * ds + ci * dc;
                        cr = nr;
                        if j % 256 == 255 {
                            let ph = step * (j + 1) as f64;
                            let (s2, c2) = ph.sin_cos();
                            cr = c2;
                            ci = s2;
                        }
                    }
                    xr[k] = sr;
                    xi[k] = si;
                }
            }
            for k in 0..npts {
                let mut v = (xr[k] * xr[k] + xi[k] * xi[k]) / (sw * sw);
                if !two {
                    v *= 2.0;
                }
                acc[k] += v;
            }
        }
        for k in 0..npts {
            let mut v = acc[k] / starts.len() as f64;
            if minthr.is_finite() && 10.0 * v.log10() <= minthr {
                v = 0.0;
            }
            pall[c * npts + k] = v;
        }
    }
    let fvec: Vec<f64> = if has_fs { grid } else { grid.iter().map(|f| f).cloned().collect() };
    let mut pm = Matrix::from_dims(vec![npts, out_cols], pall);
    let mut fm = Matrix::from_dims(vec![npts, 1], fvec);
    if x.kind == Kind::Single {
        pm.kind = Kind::Single;
        fm.kind = Kind::Single;
    }
    if nargout == 0 {
        let db: Vec<f64> = pm.data.iter().map(|v| 10.0 * v.log10()).collect();
        let dbm = Matrix::from_dims(pm.dims.clone(), db);
        it.call_function("plot", vec![Value::Num(Matrix::from_dims(vec![npts, 1], fm.data.clone())), Value::Num(dbm)], 0)?;
        it.call_function("xlabel", vec![Value::string(if has_fs { "Frequency (Hz)" } else { "Normalized Frequency (rad/sample)" })], 0)?;
        it.call_function("ylabel", vec![Value::string("Power Spectrum (dB)")], 0)?;
        it.call_function("title", vec![Value::string(&format!("Fres = {:.4} {}", rbw, if has_fs { "Hz" } else { "rad/sample" }))], 0)?;
        return out(vec![]);
    }
    let third = Matrix::from_dims(vec![0, 0], vec![]);
    out(vec![Value::Num(pm), Value::Num(fm), Value::Num(third)].into_iter().take(nargout.max(1)).collect())
}
