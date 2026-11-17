//! `strips`: a long signal drawn in horizontal strips, the first strip on top.
//!
//! The signal is cut into strips of N samples (the last one padded with NaN), centred on the
//! midrange of all the data, scaled by SCALE and drawn one strip every 1.25 times the range, all in
//! one line per column joined by NaN separators. The tick labels are the strip start positions.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

const NEG_NAN: f64 = f64::from_bits(0xfff8000000000000);

pub const NAMES: &[&str] = &["strips"];

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() > 4 {
        return rt("Too many input arguments.");
    }
    if nargout > 0 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) || m.kind == Kind::Logical || m.kind.is_int() => m.clone(),
        _ => return rt("Invalid data type.  First input argument must be a numeric array."),
    };
    if x.im.is_some() {
        it.emit(crate::interp::OutputKind::Warning("X vector complex.  I'm ignoring the imaginary part.".into()));
    }
    let (len, cols) = if x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1) { (x.numel(), 1) } else { (x.dims[0], x.numel() / x.dims[0].max(1)) };
    let num = |i: usize, label: &str| -> Result<f64, crate::error::LangError> {
        match a.get(i) {
            Some(Value::Num(m)) if m.numel() == 1 => Ok(m.data[0]),
            _ => rt(format!("{label} must be a scalar.")),
        }
    };
    // Strip length in samples, the sampling rate (1 for samples) and the label step.
    let (n, fs, in_time) = match a.len() {
        1 => (250usize, 1.0, false),
        2 => {
            let n = num(1, "N")?;
            if !(n >= 1.0) || n.fract() != 0.0 {
                return rt("N must be a positive integer.");
            }
            (n as usize, 1.0, false)
        }
        _ => {
            let sd = num(1, "SD")?;
            let fs = num(2, "Fs")?;
            if !(sd > 0.0) || !(fs > 0.0) {
                return rt("SD and Fs must be positive.");
            }
            let n = (sd * fs).round();
            if n < 1.0 {
                return rt("The strip duration is shorter than one sample.");
            }
            (n as usize, fs, true)
        }
    };
    let scale = if a.len() == 4 { num(3, "SCALE")? } else { 1.0 };
    let strips = len.div_ceil(n).max(1);
    // Range and midrange over all columns.
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in x.data.iter().take(len * cols).filter(|v| !v.is_nan()) {
        lo = lo.min(*v);
        hi = hi.max(*v);
    }
    let range = hi - lo;
    let mid = (hi + lo) / 2.0;
    let d = if range > 0.0 { 1.25 * range } else { 1.0 };
    let per = n + 1;
    let mut xd = Vec::with_capacity(strips * per);
    for _ in 0..strips {
        for k in 0..n {
            xd.push(k as f64 / fs);
        }
        xd.push(NEG_NAN);
    }
    let mut ys = Vec::with_capacity(strips * per * cols);
    for c in 0..cols {
        for s in 0..strips {
            let off = (strips - 1 - s) as f64 * d;
            for k in 0..n {
                let idx = s * n + k;
                ys.push(if idx < len { x.data[c * len + idx] * scale - mid + off } else { NEG_NAN });
            }
            ys.push(NEG_NAN);
        }
    }
    let xcol = real_col(xd);
    let ymat = Value::Num(Matrix::from_dims(vec![strips * per, cols], ys));
    it.call_function("plot", vec![xcol, ymat], 0)?;
    let ax = it.call_function("gca", vec![], 1)?.remove(0);
    let ticks: Vec<f64> = (0..strips).map(|s| s as f64 * d).collect();
    let starts: Vec<f64> = (0..strips).rev().map(|s| if in_time { s as f64 * n as f64 / fs } else { (s * n) as f64 }).collect();
    let labels = it.call_function("num2str", vec![real_col(starts)], 1)?.remove(0);
    let xmax = n as f64 / fs;
    let ylim = if range > 0.0 { vec![-0.75 * range, strips as f64 * d - 0.4 * d] } else { vec![0.0, strips as f64 * d] };
    it.call_function(
        "set",
        vec![
            ax,
            Value::chars("XLim"),
            real_row(vec![0.0, xmax]),
            Value::chars("YLim"),
            real_row(ylim),
            Value::chars("YTick"),
            real_row(ticks),
            Value::chars("YTickLabel"),
            labels,
        ],
        0,
    )?;
    out(vec![])
}
