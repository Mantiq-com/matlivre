//! `yulewalk`: recursive (IIR) filter design by a least-squares fit to a magnitude response.
//!
//! The method follows the documented algorithm (modified Yule-Walker equations for the
//! denominator, an additive decomposition of the power response, spectral factorization and a
//! least-squares fit for the numerator; Friedlander and Porat, IEEE Trans. AES-20, 1984). The
//! details that the reference leaves open were fixed against MATLAB's output:
//!
//! * the response is sampled on `npt + 1` points from 0 to 1, with `npt` rounded up to a power of
//!   two (default 512); a repeated frequency becomes a linear transition over `lap` grid points
//!   (default `fix(npt/25)`), placed as in `fir2`, including the rounding of half-integer grid
//!   positions when `lap` is odd;
//! * the correlation is the inverse FFT of the squared response, truncated to `nr = 4n` lags
//!   under the half Hamming window `0.54 + 0.46 cos(pi k/(nr-1))`;
//! * the denominator solves the equations for lags `n+1 .. nr-1` by least squares and is then
//!   stabilised (`polystab`);
//! * the causal part of the correlation (lag 0 halved) is fitted as the impulse response of
//!   `Q/A` over `nr` samples; `2 Re(Q/A)` on a `2 npt` FFT grid is factored through the real
//!   cepstrum (Nyquist term dropped) and the numerator is the least-squares fit of `B/A` to the
//!   first `nr` samples of that minimum-phase impulse response.

use crate::builtins::Ret;
use crate::cplx::C;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["yulewalk"];

const OPERANDS: &str = "Operands to the logical AND (&&) and OR (||) operators must be convertible to logical scalar values. Use the ANY or ALL functions to reduce operands to logical scalar values.";

fn is_numeric(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.kind != Kind::Char && m.kind != Kind::Logical)
}

fn float_arg(a: &[Value], i: usize, label: &str) -> R<Matrix> {
    match &a[i] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => Ok(m.clone()),
        other => rt(format!("Input {label} in yulewalk must be double/single, instead it was {}.", other.class_name())),
    }
}

fn numeric_arg(a: &[Value], i: usize, label: &str) -> R<Matrix> {
    if !is_numeric(&a[i]) {
        return rt(format!("Input {label} in yulewalk must be numeric, instead it was {}.", a[i].class_name()));
    }
    a[i].to_matrix()
}

fn is_vector(m: &Matrix) -> bool {
    m.numel() >= 1 && m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1)
}

/// Rounds to single precision when `s` is set: one arithmetic step on a single operand.
fn rs(x: f64, s: bool) -> f64 {
    if s {
        x as f32 as f64
    } else {
        x
    }
}

/// MATLAB's rounding of a non-integer colon value used as an index (half away from zero).
fn round_index(x: f64) -> f64 {
    (x.abs() + 0.5).floor().copysign(x)
}

/// The desired response on the frequency grid (`npt1` points), as `fir2` builds it.
#[allow(clippy::too_many_arguments)]
fn grid(it: &mut Interpreter, f: &[f64], mr: &[f64], mi: &[f64], npt1: f64, lap: f64, fs: bool, ms: bool) -> R<(Vec<f64>, Vec<f64>)> {
    let len = npt1.max(1.0).ceil() as usize;
    let mut hr = vec![0.0; len];
    let mut hi = vec![0.0; len];
    hr[0] = mr[0];
    hi[0] = mi[0];
    let s = fs || ms;
    let mut nb = 1.0;
    for i in 0..f.len() - 1 {
        let ne;
        if f[i + 1] == f[i] {
            nb = rs(nb - lap / 2.0, fs);
            ne = rs(nb + lap, fs);
        } else {
            ne = rs(f[i + 1] * npt1, fs).trunc();
        }
        if nb < 0.0 || ne > npt1 {
            return rt("Too abrupt an amplitude change near end of frequency interval.");
        }
        if !nb.is_finite() || !ne.is_finite() {
            return rt("Nonfinite endpoints or increment for colon operator in index.");
        }
        if nb.fract() != 0.0 || ne.fract() != 0.0 {
            it.emit_warning("Integer operands are required for colon operator when used as index.".into(), "MATLAB:colon:nonIntegerIndex".into());
        }
        let count = if ne >= nb { ((ne - nb) + 1e-10).floor() as usize + 1 } else { 0 };
        let mut idx = Vec::with_capacity(count);
        for t in 0..count {
            let j = rs(nb + t as f64, fs);
            let k = round_index(j);
            if k < 1.0 {
                return rt("Array indices must be positive integers or logical values.");
            }
            idx.push((j, k as usize));
        }
        for (j, k) in idx {
            let inc = if ne == nb { 0.0 } else { rs(rs(j - nb, fs) / rs(ne - nb, fs), fs) };
            let one_minus = rs(1.0 - inc, fs);
            let vr = rs(rs(inc * mr[i + 1], s) + rs(one_minus * mr[i], s), s);
            let vi = rs(rs(inc * mi[i + 1], s) + rs(one_minus * mi[i], s), s);
            if k > hr.len() {
                hr.resize(k, 0.0);
                hi.resize(k, 0.0);
            }
            hr[k - 1] = vr;
            hi[k - 1] = vi;
        }
        nb = rs(ne + 1.0, fs);
    }
    Ok((hr, hi))
}

fn fft_call(it: &mut Interpreter, name: &str, re: Vec<f64>, im: Vec<f64>, n: Option<usize>) -> R<(Vec<f64>, Vec<f64>)> {
    let len = re.len();
    let m = if im.iter().any(|&v| v != 0.0) { Matrix::complex(vec![1, len], re, im) } else { Matrix::from_dims(vec![1, len], re) };
    let mut args = vec![Value::Num(m)];
    if let Some(n) = n {
        args.push(scalar_value(n as f64));
    }
    let r = it.call_function(name, args, 1)?.remove(0).to_matrix()?;
    let k = r.numel();
    let im = (0..k).map(|i| r.im_at(i)).collect();
    Ok((r.data, im))
}

/// Least-squares solution (MATLAB's backslash) of `G x = y`, `G` given by columns.
fn solve(it: &mut Interpreter, cols: &[Vec<f64>], y: &[f64]) -> R<Vec<f64>> {
    let rows = y.len();
    let mut data = Vec::with_capacity(rows * cols.len());
    for c in cols {
        data.extend_from_slice(c);
    }
    let am = Value::Num(Matrix::from_dims(vec![rows, cols.len()], data));
    let ym = Value::Num(Matrix::from_dims(vec![rows, 1], y.to_vec()));
    // linsolve: column-pivoted QR and its rank warning, as MATLAB's backslash for a
    // rectangular system.
    let x = it.call_function("linsolve", vec![am, ym], 1)?.remove(0).to_matrix()?;
    Ok(x.data)
}

/// Columns of the convolution matrix of the impulse response of `1/a` (length `len`), `k` columns.
fn conv_cols(a: &[f64], len: usize, k: usize) -> Vec<Vec<f64>> {
    let mut g = vec![0.0; len];
    for t in 0..len {
        let mut s = if t == 0 { 1.0 } else { 0.0 };
        for (i, &ai) in a.iter().enumerate().skip(1) {
            if t >= i {
                s -= ai * g[t - i];
            }
        }
        g[t] = s / a[0];
    }
    (0..k).map(|c| (0..len).map(|t| if t >= c { g[t - c] } else { 0.0 }).collect()).collect()
}

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() > 5 {
        return rt("Too many input arguments.");
    }
    if a.len() < 3 {
        return rt("Wrong number of input parameters.");
    }
    if nargout > 2 {
        return rt("Too many output arguments.");
    }
    // Order.
    let nm = numeric_arg(a, 0, "N")?;
    if nm.numel() != 1 {
        return rt("Expected N to be a scalar.");
    }
    if nm.im_at(0) != 0.0 {
        return rt("Expected N to be real.");
    }
    let nv = nm.data[0];
    if !nv.is_finite() || nv.fract() != 0.0 {
        return rt("Expected N to be integer-valued.");
    }
    if nv <= 0.0 {
        return rt("Expected N to be positive.");
    }
    let n = nv as usize;
    // Frequencies and magnitudes.
    let fm = float_arg(a, 1, "F")?;
    let mm = float_arg(a, 2, "M")?;
    if !is_vector(&fm) {
        return rt("Expected F to be a vector.");
    }
    if !is_vector(&mm) {
        return rt("Expected M to be a vector.");
    }
    let nptm = if a.len() > 3 { Some(numeric_arg(a, 3, "NPT")?) } else { None };
    let lapm = if a.len() > 4 { Some(numeric_arg(a, 4, "LAP")?) } else { None };
    if fm.dims != mm.dims {
        return rt("You must specify the same number of frequencies and amplitudes.");
    }
    let k = fm.numel();
    let (f0r, f0i, fer, fei) = (fm.data[0], fm.im_at(0), fm.data[k - 1], fm.im_at(k - 1));
    if f0r != 0.0 || f0i != 0.0 || fer != 1.0 || fei != 0.0 {
        return rt("The first frequency must be 0 and the last 1.");
    }
    if fm.data.windows(2).any(|p| p[1] - p[0] < 0.0) {
        return rt("Frequencies must be non-decreasing.");
    }
    // Grid size and transition width.
    let npt = match &nptm {
        None => 512.0,
        Some(m) if m.numel() == 0 => return rt(OPERANDS),
        Some(m) if m.numel() > 1 => return rt("Size inputs must be scalar."),
        Some(m) => {
            if m.im_at(0) != 0.0 || m.data[0] < 0.0 {
                return rt("Size vector must be a row vector with real elements.");
            }
            if !m.data[0].is_finite() {
                return rt("NaN and Inf not allowed.");
            }
            let p = 2f64.powf(m.data[0].log2().ceil());
            if p.is_nan() {
                return rt("NaN and Inf not allowed.");
            }
            p
        }
    };
    let jumps = fm.data.windows(2).any(|p| p[0] == p[1]);
    let lap = match &lapm {
        None => (npt / 25.0).trunc(),
        Some(m) if m.numel() != 1 => {
            if jumps {
                return rt(OPERANDS);
            }
            0.0
        }
        Some(m) => m.data[0],
    };
    let fs = fm.kind == Kind::Single;
    let ms = mm.kind == Kind::Single;
    let mi: Vec<f64> = (0..k).map(|i| mm.im_at(i)).collect();
    let npt1 = npt + 1.0;
    let (hr, hi) = grid(it, &fm.data, &mm.data, &mi, npt1, lap, fs, ms)?;

    // Correlation: inverse transform of the squared response over the whole circle.
    let nfull = if npt >= 1.0 { 2 * npt as usize } else if npt == 0.0 { 1 } else { 2 };
    let nr = 4 * n;
    if nr > nfull {
        return rt(format!("Index exceeds the number of array elements. Index must not exceed {nfull}."));
    }
    let np = npt as usize;
    let sq: Vec<C> = hr.iter().zip(&hi).map(|(&r, &i)| C::new(r, i).mul(C::new(r, i))).collect();
    let mut pr: Vec<f64> = sq.iter().map(|c| c.re).collect();
    let mut pi: Vec<f64> = sq.iter().map(|c| c.im).collect();
    for j in (1..np).rev() {
        pr.push(sq[j].re);
        pi.push(sq[j].im);
    }
    let (rr, _) = fft_call(it, "ifft", pr, pi, None)?;
    let r: Vec<f64> = (0..nr).map(|k| rr[k] * (0.54 + 0.46 * (std::f64::consts::PI * k as f64 / (nr - 1) as f64).cos())).collect();

    // A NaN or Inf in the response makes the denominator NaN, which polystab's roots rejects.
    if r.iter().any(|v| !v.is_finite()) {
        return rt("Input to ROOTS must not contain NaN or Inf.");
    }
    // Denominator: modified Yule-Walker equations for lags n+1 .. nr-1.
    let rows: Vec<usize> = (n + 1..nr).collect();
    let cols: Vec<Vec<f64>> = (1..=n).map(|i| rows.iter().map(|&kk| r[kk.abs_diff(i)]).collect()).collect();
    let rhs: Vec<f64> = rows.iter().map(|&kk| -r[kk]).collect();
    let x = solve(it, &cols, &rhs)?;
    let mut den = vec![1.0];
    den.extend(x);
    if den.iter().any(|v| !v.is_finite()) {
        return rt("Input to ROOTS must not contain NaN or Inf.");
    }
    let den = it.call_function("polystab", vec![real_row(den)], 1)?.remove(0).to_matrix()?.data;

    // Numerator of the additive decomposition: Q/A fitted to the causal correlation.
    let mut rc = r.clone();
    rc[0] /= 2.0;
    let g = conv_cols(&den, nr, n + 1);
    let q = solve(it, &g, &rc)?;
    // Power response 2 Re(Q/A) on the 2 npt grid, then its minimum-phase factor.
    let nf = 2 * np;
    let (qr, qi) = fft_call(it, "fft", q, vec![0.0; n + 1], Some(nf))?;
    let (ar, ai) = fft_call(it, "fft", den.clone(), vec![0.0; n + 1], Some(nf))?;
    let logp: Vec<C> = (0..nf)
        .map(|j| {
            let p = 2.0 * C::new(qr[j], qi[j]).div(C::new(ar[j], ai[j])).re;
            C::new(p, 0.0).ln().mul(C::new(0.5, 0.0))
        })
        .collect();
    let (cr, ci) = fft_call(it, "ifft", logp.iter().map(|c| c.re).collect(), logp.iter().map(|c| c.im).collect(), None)?;
    let mut fr = vec![0.0; nf];
    let mut fi = vec![0.0; nf];
    fr[0] = cr[0];
    fi[0] = ci[0];
    for j in 1..np {
        fr[j] = 2.0 * cr[j];
        fi[j] = 2.0 * ci[j];
    }
    let (er, ei) = fft_call(it, "fft", fr, fi, None)?;
    let (ex_r, ex_i): (Vec<f64>, Vec<f64>) = er.iter().zip(&ei).map(|(&re, &im)| {
        let e = C::new(re, im).exp();
        (e.re, e.im)
    }).unzip();
    let (h, _) = fft_call(it, "ifft", ex_r, ex_i, None)?;
    let b = solve(it, &g, &h[..nr])?;

    let mut outs = vec![real_row(b)];
    if nargout > 1 {
        outs.push(real_row(den));
    }
    out(outs)
}
