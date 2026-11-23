//! Parametric hypothesis tests on means and variances: `ttest`, `ttest2`, `ztest`, `vartest`,
//! `vartest2`. Each returns `[h,p,ci,stats]` computed per column (or along `Dim`).

use super::dist::{chi2cdf, chi2inv, fcdf, finv, norminv_std, normcdf_std, tcdf, tinv};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::nd::Dims;
use crate::value::*;

pub const NAMES: &[&str] = &["ttest", "ttest2", "ztest", "vartest", "vartest2"];

struct Opts {
    alpha: f64,
    /// -1 left, 0 both, 1 right
    tail: i32,
    dim: Option<usize>,
    /// 1 equal, 2 unequal
    vartype: u8,
}

fn spacing(x: f64) -> f64 {
    let a = x.abs();
    if !a.is_finite() {
        return f64::NAN;
    }
    if a < f64::MIN_POSITIVE {
        return 5e-324;
    }
    f64::EPSILON * 2f64.powi(a.log2().floor() as i32)
}

fn matches_name(given: &str, full: &str) -> bool {
    !given.is_empty() && full.starts_with(&given.to_lowercase())
}

/// Chooses among the valid values of a parameter by prefix.
fn choose(v: &Value, valid: &[&str], pname: &str, upper_style: bool) -> Res<usize> {
    if !v.is_text() {
        return rt(format!("Invalid {} data type. Argument must be a character vector or string scalar.", if upper_style { "TAIL".to_string() } else { format!("'{pname}'") }));
    }
    let t = v.text()?;
    let hits: Vec<usize> = valid.iter().enumerate().filter(|(_, c)| matches_name(&t, c)).map(|(i, _)| i).collect();
    if hits.len() == 1 {
        return Ok(hits[0]);
    }
    let list = valid.iter().map(|c| format!("'{c}'")).collect::<Vec<_>>().join(", ");
    let pn = if upper_style { pname.to_uppercase() } else { format!("'{pname}'") };
    rt(format!("'{t}' is not a valid value for the {pn} argument.  Valid values are: {list}."))
}

fn parse_opts(name: &str, rest: &[Value], has_vartype: bool) -> Res<Opts> {
    let mut o = Opts { alpha: 0.05, tail: 0, dim: None, vartype: 1 };
    let upper_style = name == "ztest";
    let bad_alpha = || -> Res<()> { rt("ALPHA must be a scalar between 0 and 1.") };
    let mut alpha_v: Option<&Value> = None;
    let mut tail_v: Option<&Value> = None;
    let mut dim_v: Option<&Value> = None;
    let mut vt_v: Option<&Value> = None;
    let numeric = |v: &Value| matches!(v, Value::Num(m) if m.kind != Kind::Char);
    if !rest.is_empty() {
        if numeric(&rest[0]) {
            alpha_v = Some(&rest[0]);
            if rest.len() >= 2 {
                tail_v = Some(&rest[1]);
            }
            if rest.len() >= 3 {
                if has_vartype {
                    vt_v = Some(&rest[2]);
                    if rest.len() >= 4 {
                        dim_v = Some(&rest[3]);
                    }
                } else {
                    dim_v = Some(&rest[2]);
                }
            }
        } else if rest.len() == 1 {
            bad_alpha()?;
        } else {
            let mut names = vec!["alpha", "tail", "dim"];
            if has_vartype {
                names.push("vartype");
            }
            if rest.len() % 2 == 1 {
                return rt("Invalid parameter name or missing value.");
            }
            for pair in rest.chunks(2) {
                let key = if pair[0].is_text() { pair[0].text()? } else { String::new() };
                let hit = names.iter().position(|n| matches_name(&key, n));
                match hit {
                    Some(0) => alpha_v = Some(&pair[1]),
                    Some(1) => tail_v = Some(&pair[1]),
                    Some(2) => dim_v = Some(&pair[1]),
                    Some(3) => vt_v = Some(&pair[1]),
                    _ => return rt(format!("Invalid parameter name: {key}.")),
                }
            }
        }
    }
    if let Some(v) = alpha_v {
        match v {
            Value::Num(m) if m.is_empty() => {}
            Value::Num(m) if m.numel() == 1 && m.kind != Kind::Char && m.data[0] > 0.0 && m.data[0] < 1.0 => o.alpha = m.data[0],
            _ => bad_alpha()?,
        }
    }
    if let Some(v) = tail_v {
        match v {
            Value::Num(m) if m.is_empty() => {}
            Value::Num(m) if m.numel() == 1 && m.kind != Kind::Char && [-1.0, 0.0, 1.0].contains(&m.data[0]) => o.tail = m.data[0] as i32,
            _ => o.tail = choose(v, &["left", "both", "right"], "tail", upper_style)? as i32 - 1,
        }
    }
    if let Some(v) = vt_v {
        match v {
            Value::Num(m) if m.is_empty() => {}
            Value::Num(m) if m.numel() == 1 && m.kind != Kind::Char && [1.0, 2.0].contains(&m.data[0]) => o.vartype = m.data[0] as u8,
            _ => o.vartype = choose(v, &["equal", "unequal"], "vartype", false)? as u8 + 1,
        }
    }
    if let Some(v) = dim_v {
        match v {
            Value::Num(m) if m.is_empty() => {}
            Value::Num(m) if m.numel() == 1 && m.data[0] >= 1.0 && m.data[0].fract() == 0.0 => o.dim = Some(m.data[0] as usize),
            _ => return rt("Dimension argument must be a positive integer scalar within indexing range."),
        }
    }
    Ok(o)
}

fn default_dim(m: &Matrix) -> usize {
    m.dims.iter().position(|&d| d != 1).map_or(1, |p| p + 1)
}

/// `cat(dim, lo, hi)` for results with the dimension reduced to size 1.
fn cat2(lo: &[f64], hi: &[f64], gdims: &[usize], dim: usize) -> Matrix {
    let mut dims: Vec<usize> = gdims.to_vec();
    while dims.len() < dim {
        dims.push(1);
    }
    let stride: usize = dims[..dim - 1].iter().product();
    let n = lo.len();
    let mut out = vec![0.0; 2 * n];
    for r in 0..n {
        let (a, b) = (r % stride, r / stride);
        out[a + stride * (2 * b)] = lo[r];
        out[a + stride * (2 * b + 1)] = hi[r];
    }
    dims[dim - 1] = 2;
    Matrix::from_dims(crate::nd::normalize(dims), out)
}

fn real_data<'a>(a: &'a [Value], i: usize, name: &str) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => {
            if matches!(m.kind, Kind::Int(_)) {
                return rt("Invalid data type. First input argument must be single or double.");
            }
            Ok(m)
        }
        Some(other) => rt(format!("Undefined function 'isnan' for input arguments of type '{}'.", other.class_name())),
        None => {
            let _ = name;
            rt("Not enough input arguments.")
        }
    }
}

struct Stats {
    p: Vec<f64>,
    lo: Vec<f64>,
    hi: Vec<f64>,
}

fn finish_test(tail: i32, alpha: f64, st: &Stats, gdims: &[usize], dim: usize, single: bool, nargout: usize, extra: Vec<(&str, Matrix)>, ci_name_is_stats: bool) -> Ret {
    let h: Vec<f64> = st.p.iter().map(|&p| if p.is_nan() { f64::NAN } else if p <= alpha { 1.0 } else { 0.0 }).collect();
    let _ = tail;
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(gdims.to_vec(), h), single))];
    if nargout >= 2 || nargout == 0 {
        outs.push(Value::Num(finish(Matrix::from_dims(gdims.to_vec(), st.p.clone()), single)));
    }
    if nargout >= 3 {
        outs.push(Value::Num(finish(cat2(&st.lo, &st.hi, gdims, dim), single)));
    }
    if nargout >= 4 {
        if ci_name_is_stats {
            let mut s = StructArray::scalar();
            for (k, m) in extra {
                s.set(k, Value::Num(finish(m, single)));
            }
            outs.push(Value::Struct(s));
        } else {
            let (_, m) = extra.into_iter().next().unwrap();
            outs.push(Value::Num(finish(m, single)));
        }
    }
    Ok(Some(outs))
}

fn nonan(v: &[f64]) -> Vec<f64> {
    v.iter().copied().filter(|t| !t.is_nan()).collect()
}

fn mean_std(vals: &[f64]) -> (f64, f64) {
    let n = vals.len() as f64;
    if vals.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let m = vals.iter().sum::<f64>() / n;
    let sd = if vals.len() == 1 { 0.0 } else { (vals.iter().map(|t| (t - m) * (t - m)).sum::<f64>() / (n - 1.0)).sqrt() };
    (m, sd)
}

fn ttest(a: &[Value], nargout: usize) -> Ret {
    let mut x = real_data(a, 0, "ttest")?.clone();
    let mut m = 0.0;
    if let Some(Value::Num(mm)) = a.get(1) {
        if mm.numel() == 1 {
            m = mm.data[0];
        } else if !mm.is_empty() {
            if mm.dims != x.dims {
                return rt("The data in a paired t-test must be the same size.");
            }
            for (xv, yv) in x.data.iter_mut().zip(&mm.data) {
                *xv -= yv;
            }
        }
    } else if a.len() > 1 && !matches!(a[1], Value::Num(_)) {
        return rt("CLASSNAME argument must be a class that supports EPS, such as 'double' or 'single'.");
    }
    if matches!(a.get(1), Some(Value::Num(mm)) if mm.kind == Kind::Char) {
        return rt("CLASSNAME argument must be a class that supports EPS, such as 'double' or 'single'.");
    }
    if x.kind == Kind::Char {
        return rt("Argument following 'like' must be a variable with a data type that supports ONES, such as a double or single variable.");
    }
    let o = parse_opts("ttest", &a[2.min(a.len())..], false)?;
    let dim = o.dim.unwrap_or_else(|| default_dim(&x));
    let single = is_single(&x) || matches!(a.get(1), Some(Value::Num(mm)) if is_single(mm));
    let g = groups(&x.dims, &DimSel::Dims(vec![dim]));
    let nout = g.idx.len();
    let (mut tv, mut dfv, mut sdv) = (Vec::new(), Vec::new(), Vec::new());
    let (mut p, mut lo, mut hi) = (Vec::new(), Vec::new(), Vec::new());
    for ix in &g.idx {
        let raw: Vec<f64> = ix.iter().map(|&i| x.data[i]).collect();
        let vals = nonan(&raw);
        let n = vals.len() as f64;
        let df = (n - 1.0).max(0.0);
        let (mut xmean, mut sd) = mean_std(&vals);
        let sqrtn = n.sqrt();
        let mut xdiff = xmean - m;
        if xdiff != 0.0 && xdiff.abs() < 100.0 * sqrtn * spacing(xmean).max(spacing(m)) {
            let cmin = vals.iter().cloned().fold(f64::INFINITY, f64::min);
            if raw.iter().all(|t| t.is_nan() || *t == cmin) {
                xdiff = 0.0;
                sd = 0.0;
                xmean = m;
            }
        }
        let ser = sd / sqrtn;
        let t = xdiff / ser;
        tv.push(t);
        dfv.push(df);
        sdv.push(sd);
        match o.tail {
            0 => {
                p.push(2.0 * tcdf(-t.abs(), df));
                let crit = tinv(1.0 - o.alpha / 2.0, df) * ser;
                lo.push(xmean - crit);
                hi.push(xmean + crit);
            }
            1 => {
                p.push(tcdf(-t, df));
                lo.push(xmean - tinv(1.0 - o.alpha, df) * ser);
                hi.push(f64::INFINITY);
            }
            _ => {
                p.push(tcdf(t, df));
                lo.push(f64::NEG_INFINITY);
                hi.push(xmean + tinv(1.0 - o.alpha, df) * ser);
            }
        }
    }
    let _ = nout;
    let gd = g.dims.clone();
    let mk = |v: Vec<f64>| Matrix::from_dims(gd.clone(), v);
    finish_test(o.tail, o.alpha, &Stats { p, lo, hi }, &g.dims, dim, single, nargout, vec![("tstat", mk(tv)), ("df", mk(dfv)), ("sd", mk(sdv))], true)
}

fn ttest2(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Requires at least two input arguments");
    }
    let mut x = real_data(a, 0, "ttest2")?.clone();
    let mut y = real_data(a, 1, "ttest2")?.clone();
    let o = parse_opts("ttest2", &a[2..], true)?;
    let is_vec = |m: &Matrix| m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1);
    let dim = match o.dim {
        Some(d) => d,
        None => {
            let d = x.dims.iter().position(|&d| d != 1).or_else(|| y.dims.iter().position(|&d| d != 1)).map_or(1, |p| p + 1);
            if is_vec(&x) && is_vec(&y) {
                y.dims = if d == 2 { vec![1, y.numel()] } else { vec![y.numel(), 1] };
            }
            d
        }
    };
    let _ = &mut x;
    let reduce = |m: &Matrix| {
        let mut d = m.dims.clone();
        while d.len() < dim {
            d.push(1);
        }
        d[dim - 1] = 1;
        crate::nd::normalize(d)
    };
    if reduce(&x) != reduce(&y) {
        return rt("The data in a 2-sample t-test must be commensurate.");
    }
    let single = is_single(&x) || is_single(&y);
    let gx = groups(&x.dims, &DimSel::Dims(vec![dim]));
    let gy = groups(&y.dims, &DimSel::Dims(vec![dim]));
    let (mut tv, mut dfv, mut sd1, mut sd2) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut ses, mut diffs): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    let (mut p, mut lo, mut hi) = (Vec::new(), Vec::new(), Vec::new());
    for (ix, iy) in gx.idx.iter().zip(&gy.idx) {
        let rx: Vec<f64> = ix.iter().map(|&i| x.data[i]).collect();
        let ry: Vec<f64> = iy.iter().map(|&i| y.data[i]).collect();
        let (vx, vy) = (nonan(&rx), nonan(&ry));
        let (nx, ny) = (vx.len() as f64, vy.len() as f64);
        let var = |v: &[f64]| {
            let (_, s) = mean_std(v);
            s * s
        };
        let (s2x, s2y) = (var(&vx), var(&vy));
        let (xm, _) = mean_std(&vx);
        let (ym, _) = mean_std(&vy);
        let mut diff = xm - ym;
        let sqrtn = nx.sqrt() + ny.sqrt();
        let mut fix = false;
        if diff != 0.0 && diff.abs() < sqrtn * 100.0 * spacing(xm).max(spacing(ym)) {
            let cmin = rx.iter().cloned().filter(|t| !t.is_nan()).fold(f64::INFINITY, f64::min);
            fix = rx.iter().all(|t| t.is_nan() || *t == cmin) && ry.iter().all(|t| t.is_nan() || *t == cmin);
            if fix {
                diff = 0.0;
            }
        }
        let (ratio, mut dfe, se);
        if o.vartype == 1 {
            dfe = nx + ny - 2.0;
            let mut sp = (((nx - 1.0) * s2x + (ny - 1.0) * s2y) / dfe).sqrt();
            if fix {
                sp = 0.0;
            }
            se = sp * (1.0 / nx + 1.0 / ny).sqrt();
            ratio = diff / se;
            sd1.push(sp);
            sd2.push(f64::NAN);
        } else {
            let (sxb, syb) = (s2x / nx, s2y / ny);
            dfe = (sxb + syb).powi(2) / (sxb * sxb / (nx - 1.0) + syb * syb / (ny - 1.0));
            let mut s = (sxb + syb).sqrt();
            if fix {
                s = 0.0;
            }
            se = s;
            ratio = diff / se;
            sd1.push(s2x.sqrt());
            sd2.push(s2y.sqrt());
        }
        tv.push(ratio);
        dfv.push(dfe);
        ses.push(se);
        diffs.push(diff);
        let _ = &mut dfe;
    }
    let all_zero = !ses.is_empty() && ses.iter().all(|&s| s == 0.0);
    for k in 0..ses.len() {
        let (ratio, se, diff) = (tv[k], ses[k], diffs[k]);
        let dfe = if o.vartype == 2 && all_zero { 1.0 } else { dfv[k] };
        match o.tail {
            0 => {
                p.push(2.0 * tcdf(-ratio.abs(), dfe));
                let spread = tinv(1.0 - o.alpha / 2.0, dfe) * se;
                lo.push(diff - spread);
                hi.push(diff + spread);
            }
            1 => {
                p.push(tcdf(-ratio, dfe));
                lo.push(diff - tinv(1.0 - o.alpha, dfe) * se);
                hi.push(f64::INFINITY);
            }
            _ => {
                p.push(tcdf(ratio, dfe));
                lo.push(f64::NEG_INFINITY);
                hi.push(diff + tinv(1.0 - o.alpha, dfe) * se);
            }
        }
    }
    let gd = gx.dims.clone();
    let mk = |v: Vec<f64>| Matrix::from_dims(gd.clone(), v);
    let sd = if o.vartype == 1 { mk(sd1) } else { cat2(&sd1, &sd2, &gd, dim) };
    finish_test(o.tail, o.alpha, &Stats { p, lo, hi }, &gx.dims, dim, single, nargout, vec![("tstat", mk(tv)), ("df", mk(dfv)), ("sd", sd)], true)
}

fn ztest(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    let x = real_data(a, 0, "ztest")?;
    let m = match &a[1] {
        Value::Num(m) if m.numel() == 1 => m.data[0],
        _ => return rt("M must be a scalar."),
    };
    let sigma = match &a[2] {
        Value::Num(s) if s.numel() == 1 && s.data[0] >= 0.0 => s.data[0],
        _ => return rt("SIGMA must be a positive scalar."),
    };
    let o = parse_opts("ztest", &a[3..], false)?;
    let dim = o.dim.unwrap_or_else(|| default_dim(x));
    let single = is_single(x);
    let g = groups(&x.dims, &DimSel::Dims(vec![dim]));
    let (mut zv, mut p, mut lo, mut hi) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for ix in &g.idx {
        let vals = nonan(&ix.iter().map(|&i| x.data[i]).collect::<Vec<_>>());
        let n = vals.len() as f64;
        let xm = if vals.is_empty() { f64::NAN } else { vals.iter().sum::<f64>() / n };
        let ser = sigma / n.sqrt();
        let z = (xm - m) / ser;
        zv.push(z);
        match o.tail {
            0 => {
                p.push(2.0 * normcdf_std(-z.abs()));
                let crit = norminv_std(1.0 - o.alpha / 2.0) * ser;
                lo.push(xm - crit);
                hi.push(xm + crit);
            }
            1 => {
                p.push(normcdf_std(-z));
                lo.push(xm - norminv_std(1.0 - o.alpha) * ser);
                hi.push(f64::INFINITY);
            }
            _ => {
                p.push(normcdf_std(z));
                lo.push(f64::NEG_INFINITY);
                hi.push(xm + norminv_std(1.0 - o.alpha) * ser);
            }
        }
    }
    let zm = Matrix::from_dims(g.dims.clone(), zv);
    finish_test(o.tail, o.alpha, &Stats { p, lo, hi }, &g.dims, dim, single, nargout, vec![("z", zm)], false)
}

fn check_dim(x: &Matrix, dim: usize) -> Res<()> {
    if dim < 1 || dim > x.dims.len() {
        return rt(format!("DIM must be an integer between 1 and {}.", x.dims.len()));
    }
    Ok(())
}

fn vartest(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let x = real_data(a, 0, "vartest")?;
    let v = match &a[1] {
        Value::Num(m) if m.numel() == 1 && !m.is_complex() && m.data[0] >= 0.0 => m.data[0],
        _ => return rt("V must be a scalar variance value."),
    };
    let o = parse_opts("vartest", &a[2..], false)?;
    let dim = o.dim.unwrap_or_else(|| default_dim(x));
    check_dim(x, dim)?;
    let single = is_single(x);
    let g = groups(&x.dims, &DimSel::Dims(vec![dim]));
    let (mut cs, mut dfv, mut p, mut lo, mut hi) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for ix in &g.idx {
        let vals = nonan(&ix.iter().map(|&i| x.data[i]).collect::<Vec<_>>());
        let n = vals.len() as f64;
        let df = (n - 1.0).max(0.0);
        let m = vals.iter().sum::<f64>() / n.max(1.0);
        let sumsq: f64 = vals.iter().map(|t| (t - m) * (t - m)).sum();
        let chi = if v > 0.0 { sumsq / v } else if sumsq == 0.0 { f64::NAN } else { f64::INFINITY };
        cs.push(chi);
        dfv.push(df);
        match o.tail {
            0 => {
                let pc = chi2cdf(chi, df, false);
                p.push(2.0 * pc.min(1.0 - pc));
                lo.push(sumsq / chi2inv(1.0 - o.alpha / 2.0, df));
                hi.push(sumsq / chi2inv(o.alpha / 2.0, df));
            }
            1 => {
                p.push(chi2cdf(chi, df, true));
                lo.push(sumsq / chi2inv(1.0 - o.alpha, df));
                hi.push(f64::INFINITY);
            }
            _ => {
                p.push(chi2cdf(chi, df, false));
                lo.push(0.0);
                hi.push(sumsq / chi2inv(o.alpha, df));
            }
        }
    }
    let gd = g.dims.clone();
    let mk = |v: Vec<f64>| Matrix::from_dims(gd.clone(), v);
    finish_test(o.tail, o.alpha, &Stats { p, lo, hi }, &g.dims, dim, single, nargout, vec![("chisqstat", mk(cs)), ("df", mk(dfv))], true)
}

fn vartest2(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let x = real_data(a, 0, "vartest2")?.clone();
    let mut y = real_data(a, 1, "vartest2")?.clone();
    let o = parse_opts("vartest2", &a[2..], false)?;
    let dim = match o.dim {
        Some(d) => d,
        None => {
            let d = default_dim(&x);
            let is_vec = |m: &Matrix| m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1);
            if is_vec(&x) && is_vec(&y) {
                y.dims = if x.dims[0] == 1 { vec![1, y.numel()] } else { vec![y.numel(), 1] };
            }
            d
        }
    };
    check_dim(&x, dim)?;
    let reduce = |m: &Matrix| {
        let mut d = m.dims.clone();
        while d.len() < dim {
            d.push(1);
        }
        d[dim - 1] = 1;
        crate::nd::normalize(d)
    };
    if reduce(&x) != reduce(&y) {
        return rt("The data in a 2-sample F-test must be commensurate.");
    }
    let single = is_single(&x) || is_single(&y);
    let gx = groups(&x.dims, &DimSel::Dims(vec![dim]));
    let gy = groups(&y.dims, &DimSel::Dims(vec![dim]));
    let stat = |m: &Matrix, ix: &Vec<usize>| -> (f64, f64) {
        let vals = nonan(&ix.iter().map(|&i| m.data[i]).collect::<Vec<_>>());
        let n = vals.len() as f64;
        let df = (n - 1.0).max(0.0);
        let mean = vals.iter().sum::<f64>() / n.max(1.0);
        let ss: f64 = vals.iter().map(|t| (t - mean) * (t - mean)).sum();
        (df, if df > 0.0 { ss / df } else { f64::NAN })
    };
    let (mut fv, mut d1v, mut d2v, mut p, mut lo, mut hi) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (ix, iy) in gx.idx.iter().zip(&gy.idx) {
        let (df1, vx) = stat(&x, ix);
        let (df2, vy) = stat(&y, iy);
        let f = if vy > 0.0 {
            vx / vy
        } else if vx > 0.0 {
            f64::INFINITY
        } else {
            f64::NAN
        };
        fv.push(f);
        d1v.push(df1);
        d2v.push(df2);
        match o.tail {
            0 => {
                p.push(2.0 * fcdf(f, df1, df2, false).min(fcdf(f, df1, df2, true)));
                lo.push(f * finv(o.alpha / 2.0, df2, df1));
                hi.push(f / finv(o.alpha / 2.0, df1, df2));
            }
            1 => {
                p.push(fcdf(f, df1, df2, true));
                lo.push(f * finv(o.alpha, df2, df1));
                hi.push(f64::INFINITY);
            }
            _ => {
                p.push(fcdf(f, df1, df2, false));
                lo.push(0.0);
                hi.push(f / finv(o.alpha, df1, df2));
            }
        }
    }
    let gd: Dims = gx.dims.clone();
    let mk = |v: Vec<f64>| Matrix::from_dims(gd.clone(), v);
    finish_test(o.tail, o.alpha, &Stats { p, lo, hi }, &gx.dims, dim, single, nargout, vec![("fstat", mk(fv)), ("df1", mk(d1v)), ("df2", mk(d2v))], true)
}

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "ttest" => ttest(a, nargout),
        "ttest2" => ttest2(a, nargout),
        "ztest" => ztest(a, nargout),
        "vartest" => vartest(a, nargout),
        _ => vartest2(a, nargout),
    }
}
