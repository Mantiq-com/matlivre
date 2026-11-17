//! `ifsst`: inverse Fourier synchrosqueezed transform.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["ifsst"];

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > 1 {
        return rt("Too many output arguments.");
    }
    let s = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected SST to be one of these types:  single, double"),
    };
    if s.dims.len() != 2 || s.dims[0] == 1 || s.dims[1] == 1 {
        return rt("The input SST must be a matrix.");
    }
    let (nr, nc) = (s.dims[0], s.dims[1]);
    // Window.
    let mut next = 1;
    let mut win_default = true;
    let mut wvec: Vec<f64> = vec![];
    if let Some(Value::Num(m)) = a.get(1) {
        if m.kind != Kind::Char {
            next = 2;
            if m.numel() == 1 {
                if !(m.data[0] > 0.0) {
                    return rt("Expected WINDOW to be positive.");
                }
                let l = m.data[0] as usize;
                wvec = it.call_function("kaiser", vec![scalar_value(l as f64), scalar_value(10.0)], 1)?.remove(0).to_matrix()?.data;
                win_default = false;
            } else if m.numel() > 1 {
                wvec = m.data.clone();
                win_default = false;
            }
        }
    }
    if win_default {
        wvec = it.call_function("kaiser", vec![scalar_value(256.0), scalar_value(10.0)], 1)?.remove(0).to_matrix()?.data;
    }
    let l = wvec.len();
    let two_sided = nr == l;
    let one_sided = nr == l / 2 + 1;
    if !two_sided && !one_sided {
        if win_default {
            return rt("The default Kaiser window is inconsistent with the input matrix SST. Specify a window with compatible length.");
        }
        return rt("The length of WINDOW is inconsistent with the input matrix SST.");
    }
    // Remaining: [f, freqrange] or iridge, then NumFrequencyBins.
    let mut rest: Vec<&Value> = vec![];
    let mut i = next;
    let mut nbins = 4usize;
    while i < a.len() {
        if let Value::Num(m) = &a[i] {
            if m.kind == Kind::Char {
                let t = a[i].text().unwrap_or_default();
                if t.to_lowercase() != "numfrequencybins" {
                    return rt(format!("Expected input to match one of these values:  'NumFrequencyBins'  The input, '{t}', did not match any of the valid values."));
                }
                let v = a.get(i + 1).map(|v| v.scalar().unwrap_or(f64::NAN)).unwrap_or(f64::NAN);
                if !(v > 0.0) || v.fract() != 0.0 {
                    return rt("Expected NumFrequencyBins to be a positive integer.");
                }
                nbins = v as usize;
                i += 2;
                continue;
            }
        }
        rest.push(&a[i]);
        i += 1;
    }
    let gc = wvec[l / 2];
    let si = s.im.clone().unwrap_or_else(|| vec![0.0; s.numel()]);
    // Rows used for every column.
    let recon = |rows: &dyn Fn(usize) -> bool, col: usize, work: &Matrix, wi: &[f64]| -> (f64, f64) {
        let (mut re, mut im) = (0.0, 0.0);
        let mut dc = 0.0;
        let mut ny = 0.0;
        for r in 0..nr {
            if !rows(r) {
                continue;
            }
            let k = col * nr + r;
            re += work.data[k];
            im += wi[k];
            if one_sided {
                if r == 0 {
                    dc = work.data[k];
                }
                if l % 2 == 0 && r == nr - 1 {
                    ny = work.data[k];
                }
            }
        }
        if one_sided {
            ((2.0 * re - dc - ny) / (l as f64 * gc), 0.0)
        } else {
            (re / (l as f64 * gc), im / (l as f64 * gc))
        }
    };
    let mut outr: Vec<f64> = vec![];
    let mut outi: Vec<f64> = vec![];
    let mut ncols_out = 1;
    if rest.is_empty() {
        for c in 0..nc {
            let (re, im) = recon(&|_| true, c, &s, &si);
            outr.push(re);
            outi.push(im);
        }
    } else if rest.len() >= 2 {
        let fv = rest[0].to_matrix()?.data;
        let fr = rest[1].to_matrix()?.data;
        if fv.len() != nr {
            return rt("The length of F must match the number of rows of TFM.");
        }
        let (fmin, fmax) = (fv.iter().cloned().fold(f64::INFINITY, f64::min), fv.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
        if fr.len() != 2 || !(fr[1] > fr[0]) || fr[0] < fmin || fr[1] > fmax {
            return rt("FREQRANGE must be a two-element numeric vector with strictly increasing values. The values must lie within the range of the values in F.");
        }
        let lo_i = (0..nr).find(|&r| fv[r] >= fr[0]).unwrap_or(0);
        let hi_i = (0..nr).find(|&r| fv[r] >= fr[1]).unwrap_or(nr - 1);
        for c in 0..nc {
            let (re, im) = recon(&|r| r >= lo_i && r <= hi_i, c, &s, &si);
            outr.push(re);
            outi.push(im);
        }
    } else {
        let rm = rest[0].to_matrix()?;
        let nrid = if rm.dims[0] == 1 { rm.numel() } else { rm.dims[0] };
        let ncid = if rm.dims[0] == 1 { 1 } else { rm.numel() / nrid.max(1) };
        if nrid != nc {
            return rt("The number of rows of IRIDGE must match the number of columns of SST.");
        }
        ncols_out = ncid;
        let work = s.clone();
        let wim = si.clone();
        for j in 0..ncid {
            for c in 0..nc {
                let rid = rm.data[j * nrid + c];
                let (lo, hi) = (rid as i64 - nbins as i64, rid as i64 + nbins as i64);
                let rows = |r: usize| {
                    let rr = r as i64 + 1;
                    rr >= lo && rr <= hi
                };
                let (re, im) = recon(&rows, c, &work, &wim);
                outr.push(re);
                outi.push(im);
            }
        }
    }
    let mut res = if two_sided { Matrix::complex(vec![nc, ncols_out], outr, outi) } else { Matrix::from_dims(vec![nc, ncols_out], outr) };
    if s.kind == Kind::Single {
        res.kind = Kind::Single;
    }
    out(vec![Value::Num(res)])
}
