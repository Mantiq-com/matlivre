//! Harmonic measurements on a real vector of tones: `thd`, `sfdr`, `toi`.
//!
//! The spectrum is a power-scaled periodogram with a Kaiser window (beta 38); the power of a tone
//! is the sum of the periodogram over its mainlobe divided by the window's equivalent noise
//! bandwidth in bins, and its frequency is the power-weighted centroid of that band. Spectra of
//! pure tones are exact; absent harmonics (window sidelobe floors) are not reproduced.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["thd", "sfdr", "toi", "snr", "sinad"];

/// Half-width, in bins, of the band summed around a tone (the Kaiser(38) mainlobe is 12.2 bins).
const BAND: usize = 14;
/// Bins at the low end that belong to the mainlobe of the DC term.
const DC_GUARD: usize = 13;

struct Spec {
    p: Vec<f64>,
    f: Vec<f64>,
    en: f64,
    fs: f64,
}

fn spectrum(it: &mut Interpreter, x: &[f64], fs: f64) -> Result<Spec, crate::error::LangError> {
    let n = x.len();
    let w = it.call_function("kaiser", vec![scalar_value(n as f64), scalar_value(38.0)], 1)?.remove(0).to_matrix()?;
    let sw: f64 = w.data.iter().sum();
    let sw2: f64 = w.data.iter().map(|v| v * v).sum();
    let r = it.call_function("periodogram", vec![real_col(x.to_vec()), Value::Num(w), scalar_value(n as f64), scalar_value(fs), Value::chars("power")], 2)?;
    Ok(Spec { p: r[0].to_matrix()?.data, f: r[1].to_matrix()?.data, en: n as f64 * sw2 / (sw * sw), fs })
}

impl Spec {
    fn peak(&self, lo: usize, hi: usize) -> Option<usize> {
        (lo..=hi.min(self.p.len() - 1)).max_by(|a, b| self.p[*a].partial_cmp(&self.p[*b]).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// The monotone mainlobe around bin `i`: out to where the periodogram stops decreasing.
    fn lobe(&self, i: usize) -> (usize, usize) {
        let (mut lo, mut hi) = (i, i);
        while lo > 0 && self.p[lo - 1] < self.p[lo] {
            lo -= 1;
        }
        while hi + 1 < self.p.len() && self.p[hi + 1] < self.p[hi] {
            hi += 1;
        }
        // The valley bin on the high side belongs to the next lobe.
        if hi > i && hi + 1 < self.p.len() {
            hi -= 1;
        }
        (lo, hi)
    }

    /// Power and centroid frequency of the mainlobe around bin `i`.
    fn tone(&self, i: usize) -> (f64, f64) {
        let (lo, hi) = self.lobe(i);
        let s: f64 = self.p[lo..=hi].iter().sum();
        let c: f64 = (lo..=hi).map(|k| self.f[k] * self.p[k]).sum();
        (s / self.en, if s > 0.0 { c / s } else { self.f[i] })
    }

    fn nearest(&self, f: f64) -> usize {
        let res = self.f[1] - self.f[0];
        ((f / res).round().max(0.0) as usize).min(self.p.len() - 1)
    }

    /// The fundamental: the largest peak above the DC mainlobe.
    fn fundamental(&self) -> Option<usize> {
        if self.p.len() <= DC_GUARD + 1 {
            return None;
        }
        self.peak(DC_GUARD + 1, self.p.len() - 1)
    }
}

fn scalar_arg(a: &[Value], i: usize, label: &str, checks: &[&str]) -> Result<f64, crate::error::LangError> {
    let m = match &a[i] {
        Value::Num(m) => m.clone(),
        other => return rt(format!("Expected input number {}, {label}, to be one of these types:  numeric  Instead its type was {}.", i + 1, other.class_name())),
    };
    for c in checks {
        match *c {
            "scalar" if m.numel() != 1 => return rt(format!("Expected input number {}, {label}, to be a scalar.", i + 1)),
            "positive" if !(m.data[0] > 0.0) => return rt(format!("Expected input number {}, {label}, to be positive.", i + 1)),
            "integer" if m.data[0].fract() != 0.0 => return rt(format!("Expected input number {}, {label}, to be integer-valued.", i + 1)),
            "gt1" if !(m.data[0] > 1.0) => return rt(format!("Expected input number {}, {label}, to be a scalar with value > 1.", i + 1)),
            _ => {}
        }
    }
    Ok(m.data[0])
}

/// The signal vector, validated the way `thd`, `sfdr` and `toi` do; the flag tells single.
fn signal(a: &[Value]) -> Result<(Vec<f64>, bool, Matrix), crate::error::LangError> {
    let m = match &a[0] {
        Value::Num(m) if m.numel() > 1 && m.dims.len() == 2 => m.clone(),
        _ => return rt("The input must be a non-empty vector or two-dimensional matrix."),
    };
    if !matches!(m.kind, Kind::Double | Kind::Single) {
        return rt(format!("Expected input number 1, X, to be one of these types:\n\nsingle, double\n\nInstead its type was {}.", m.kind.name()));
    }
    if m.im.is_some() {
        return rt("Expected input number 1, X, to be real.");
    }
    if m.data.iter().any(|v| !v.is_finite()) {
        return rt("Expected input number 1, X, to be finite.");
    }
    if m.dims[0] != 1 && m.dims[1] != 1 {
        return rt("Expected input number 1, X, to be a vector.");
    }
    Ok((m.data.clone(), m.kind == Kind::Single, m))
}

fn db(v: f64) -> f64 {
    10.0 * v.log10()
}

fn finish(vals: Vec<Value>, single: bool, nargout: usize) -> Ret {
    let v: Vec<Value> = vals.into_iter().take(nargout.max(1)).collect();
    out(if single { v.into_iter().map(crate::signal::to_single).collect() } else { v })
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > match name {
        "toi" => 5,
        "snr" | "sinad" => 2,
        _ => 3,
    } {
        return rt("Too many output arguments.");
    }
    match name {
        "thd" => thd(it, a, nargout),
        "sfdr" => sfdr(it, a, nargout),
        "snr" => snr(it, a, nargout),
        "sinad" => sinad(it, a, nargout),
        _ => toi(it, a, nargout),
    }
}

fn thd(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let (x, single, _) = signal(a)?;
    let mut i = 1;
    let mut fs = 1.0;
    if let Some(Value::Num(m)) = a.get(i) {
        if m.kind != Kind::Char || m.numel() == 1 {
            fs = scalar_arg(a, i, "Fs", &["scalar", "positive"])?;
            i += 1;
        }
    }
    let mut n = 6.0;
    if let Some(Value::Num(m)) = a.get(i) {
        if m.kind != Kind::Char || m.numel() == 1 {
            n = scalar_arg(a, i, "N", &["scalar", "positive", "integer", "gt1"])?;
            i += 1;
        }
    }
    let mut aliased = false;
    for v in &a[i..] {
        let t = v.text().unwrap_or_default();
        if !t.is_empty() && "aliased".starts_with(&t.to_lowercase()) {
            aliased = true;
        } else {
            return rt(format!("Unknown option:  '{t}'"));
        }
    }
    let n = n as usize;
    if x.len() < 4 {
        return finish(vec![scalar_value(f64::NAN), real_col(vec![f64::NAN]), real_col(vec![f64::NAN])], single, nargout);
    }
    let sp = spectrum(it, &x, fs)?;
    let Some(i1) = sp.fundamental().filter(|i| sp.p[*i] > 0.0) else {
        return finish(vec![scalar_value(f64::NAN), real_col(vec![f64::NAN; n]), real_col(vec![f64::NAN; n])], single, nargout);
    };
    let (p1, f1) = sp.tone(i1);
    let mut pw = vec![p1];
    let mut fr = vec![f1];
    for k in 2..=n {
        let mut fe = k as f64 * f1;
        if fe > sp.fs / 2.0 {
            if !aliased {
                pw.push(f64::NAN);
                fr.push(f64::NAN);
                continue;
            }
            let m = fe % sp.fs;
            fe = if m > sp.fs / 2.0 { sp.fs - m } else { m };
        }
        let c = sp.nearest(fe);
        let ik = sp.peak(c.saturating_sub(2), c + 2).unwrap_or(c);
        let (p, f) = sp.tone(ik);
        pw.push(p);
        fr.push(f);
    }
    let harm: f64 = pw[1..].iter().filter(|v| v.is_finite()).sum();
    let r = db(harm / p1);
    finish(vec![scalar_value(r), real_col(pw.iter().map(|v| db(*v)).collect()), real_col(fr)], single, nargout)
}

fn sfdr(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let (x, single, _) = signal(a)?;
    if a.len() > 3 {
        return rt("Too many input arguments.");
    }
    let fs = if a.len() > 1 { scalar_arg(a, 1, "Fs", &["scalar", "positive"])? } else { 1.0 };
    if let Some(Value::Num(m)) = a.get(2) {
        if m.kind == Kind::Char {
            return rt(format!("Unknown option:  '{}'", a[2].text().unwrap_or_default()));
        }
    }
    let msd = if a.len() > 2 { scalar_arg(a, 2, "MSD", &["scalar", "positive"])? } else { 0.0 };
    if x.len() < 4 {
        return finish(vec![scalar_value(f64::NAN); 3], single, nargout);
    }
    let sp = spectrum(it, &x, fs)?;
    let Some(i1) = sp.fundamental().filter(|i| sp.p[*i] > 0.0) else {
        return finish(vec![scalar_value(f64::NAN); 3], single, nargout);
    };
    let (p1, f1) = sp.tone(i1);
    // Spurs: local maxima outside the fundamental's mainlobe and the minimum spacing.
    let res = sp.f[1] - sp.f[0];
    let guard = (BAND as f64).max(msd / res) as usize;
    let mut best: Option<usize> = None;
    for k in (DC_GUARD + 1)..sp.p.len() {
        if (k as i64 - i1 as i64).unsigned_abs() as usize <= guard || (sp.f[k] - f1).abs() < msd {
            continue;
        }
        let left = if k > 0 { sp.p[k - 1] } else { 0.0 };
        let right = if k + 1 < sp.p.len() { sp.p[k + 1] } else { 0.0 };
        if sp.p[k] >= left && sp.p[k] >= right && best.map_or(true, |b| sp.p[k] > sp.p[b]) {
            best = Some(k);
        }
    }
    let Some(ib) = best else {
        return finish(vec![scalar_value(f64::NAN); 3], single, nargout);
    };
    let (p2, f2) = sp.tone(ib);
    finish(vec![scalar_value(db(p1) - db(p2)), scalar_value(db(p2)), scalar_value(f2)], single, nargout)
}

fn toi(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let (x, single, _) = signal(a)?;
    if a.len() > 2 {
        return rt("Too many input arguments.");
    }
    let fs = if a.len() > 1 { scalar_arg(a, 1, "Fs", &["scalar", "positive"])? } else { 1.0 };
    let nan_out = |single: bool| finish(vec![scalar_value(f64::NAN), real_row(vec![f64::NAN; 2]), real_row(vec![f64::NAN; 2]), real_row(vec![f64::NAN; 2]), real_row(vec![f64::NAN; 2])], single, nargout);
    if x.len() < 4 {
        return nan_out(single);
    }
    let sp = spectrum(it, &x, fs)?;
    // The two largest local maxima above the DC mainlobe are the fundamentals.
    let mut peaks: Vec<usize> = ((DC_GUARD + 1)..sp.p.len().saturating_sub(1)).filter(|&k| sp.p[k] > 0.0 && sp.p[k] >= sp.p[k - 1] && sp.p[k] >= sp.p[k + 1]).collect();
    peaks.sort_by(|a, b| sp.p[*b].partial_cmp(&sp.p[*a]).unwrap_or(std::cmp::Ordering::Equal));
    if peaks.len() < 2 || sp.p[peaks[1]] < sp.p[peaks[0]] * 1e-15 {
        return nan_out(single);
    }
    let mut two = [peaks[0], peaks[1]];
    two.sort();
    let (p1, f1) = sp.tone(two[0]);
    let (p2, f2) = sp.tone(two[1]);
    let mut pim = [f64::NAN; 2];
    let mut fim = [f64::NAN; 2];
    for (k, fe) in [2.0 * f1 - f2, 2.0 * f2 - f1].iter().enumerate() {
        if *fe > 0.0 && *fe < sp.fs / 2.0 {
            let c = sp.nearest(*fe);
            let ik = sp.peak(c.saturating_sub(2), c + 2).unwrap_or(c);
            let (p, f) = sp.tone(ik);
            pim[k] = p;
            fim[k] = f;
        }
    }
    let pf = [db(p1), db(p2)];
    let pimd = [db(pim[0]), db(pim[1])];
    let oip3 = ((pf[0] + (pf[0] - pimd[0]) / 2.0) + (pf[1] + (pf[1] - pimd[1]) / 2.0)) / 2.0;
    finish(vec![scalar_value(oip3), real_row(pf.to_vec()), real_row(vec![f1, f2]), real_row(pimd.to_vec()), real_row(fim.to_vec())], single, nargout)
}

fn snr(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    // snr(x,y): the ratio of the summed squared magnitudes.
    if a.len() == 2 && matches!((&a[0], &a[1]), (Value::Num(x), Value::Num(y)) if x.numel() > 1 && y.numel() > 1 && (x.numel() == y.numel() || x.dims != y.dims)) {
        let (x, y) = (a[0].to_matrix()?, a[1].to_matrix()?);
        if x.numel() != y.numel() {
            return rt("Inputs X and Y must have the same dimensions.");
        }
        let sq = |m: &Matrix| -> f64 { (0..m.numel()).map(|i| m.data[i] * m.data[i] + m.im_at(i) * m.im_at(i)).sum() };
        let r = 10.0 * (sq(&x) / sq(&y)).log10();
        let single = x.kind == Kind::Single || y.kind == Kind::Single;
        return finish(vec![scalar_value(r)], single, 1);
    }
    let (x, single, _) = signal(a)?;
    let mut i = 1;
    let mut fs = 1.0;
    if let Some(Value::Num(m)) = a.get(i) {
        if m.kind != Kind::Char || m.numel() == 1 {
            fs = scalar_arg(a, i, "Fs", &["scalar", "positive"])?;
            i += 1;
        }
    }
    let mut n = 6.0;
    if let Some(Value::Num(m)) = a.get(i) {
        if m.kind != Kind::Char || m.numel() == 1 {
            n = scalar_arg(a, i, "N", &["scalar", "positive", "integer", "gt1"])?;
            i += 1;
        }
    }
    let mut aliased = false;
    for v in &a[i..] {
        let t = v.text().unwrap_or_default();
        if !t.is_empty() && "aliased".starts_with(&t.to_lowercase()) {
            aliased = true;
        } else if !t.is_empty() && "omitaliases".starts_with(&t.to_lowercase()) {
            aliased = false;
        } else {
            return rt(format!("Unknown option:  '{t}'"));
        }
    }
    let n = n as usize;
    let nan = |single: bool| finish(vec![scalar_value(f64::NAN), scalar_value(f64::NAN)], single, nargout);
    if x.len() < 4 {
        return nan(single);
    }
    let sp = spectrum(it, &x, fs)?;
    let Some(i1) = sp.fundamental().filter(|i| sp.p[*i] > 0.0) else {
        return nan(single);
    };
    let (lo, hi) = sp.lobe(i1);
    let (p1, f1) = sp.tone(i1);
    let mut mask = vec![false; sp.p.len()];
    for m in mask.iter_mut().take(DC_GUARD + 1) {
        *m = true;
    }
    for m in mask.iter_mut().take(hi + 1).skip(lo) {
        *m = true;
    }
    for k in 2..=n {
        let mut fe = k as f64 * f1;
        if fe > sp.fs / 2.0 {
            if !aliased {
                continue;
            }
            let m = fe % sp.fs;
            fe = if m > sp.fs / 2.0 { sp.fs - m } else { m };
        }
        let c = sp.nearest(fe);
        let ik = sp.peak(c.saturating_sub(2), c + 2).unwrap_or(c);
        let (l, h) = sp.lobe(ik);
        for m in mask.iter_mut().take(h + 1).skip(l) {
            *m = true;
        }
    }
    let noise: f64 = sp.p.iter().zip(&mask).filter(|(_, m)| !**m).map(|(v, _)| *v).sum::<f64>() / sp.en;
    finish(vec![scalar_value(db(p1 / noise)), scalar_value(db(noise))], single, nargout)
}

fn sinad(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let (x, single, _) = signal(a)?;
    if a.len() > 2 {
        return rt("Too many input arguments.");
    }
    let fs = if a.len() > 1 { scalar_arg(a, 1, "Fs", &["scalar", "positive"])? } else { 1.0 };
    if x.len() < 4 {
        return finish(vec![scalar_value(f64::NAN), scalar_value(f64::NAN)], single, nargout);
    }
    let sp = spectrum(it, &x, fs)?;
    let Some(i1) = sp.fundamental().filter(|i| sp.p[*i] > 0.0) else {
        return finish(vec![scalar_value(f64::NAN), scalar_value(f64::NAN)], single, nargout);
    };
    let (lo, hi) = sp.lobe(i1);
    let (p1, _) = sp.tone(i1);
    let rest: f64 = (DC_GUARD + 1..sp.p.len()).filter(|k| *k < lo || *k > hi).map(|k| sp.p[k]).sum::<f64>() / sp.en;
    finish(vec![scalar_value(db(p1 / rest)), scalar_value(db(rest))], single, nargout)
}
