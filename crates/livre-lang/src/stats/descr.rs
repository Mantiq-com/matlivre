//! Descriptive statistics: `geomean`, `harmmean`, `trimmean`, `skewness`, `kurtosis`, `moment`,
//! `range`, `mad`, `zscore`, `tiedrank` and the `nan*` reductions.

use super::util::*;
use crate::builtins::{one, Ret};
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &[
    "geomean", "harmmean", "trimmean", "skewness", "kurtosis", "moment", "range", "mad", "zscore", "tiedrank", "nanmax", "nanmin", "nanmean",
    "nanmedian", "nanstd", "nanvar", "nansum",
];

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "geomean" | "harmmean" => means(name, a),
        "trimmean" => trimmean(a),
        "skewness" | "kurtosis" => shape(name, a),
        "moment" => moment(a),
        "range" => range(a),
        "mad" => mad(a),
        "zscore" => zscore(a, nargout),
        "tiedrank" => tiedrank(a, nargout),
        "nanmax" | "nanmin" => nanextreme(name, a, nargout),
        _ => nanreduce(name, a),
    }
}

/// `'omitnan'` / `'includenan'` (any case, abbreviations accepted).
fn nan_flag(v: &Value) -> Option<bool> {
    if !v.is_text() {
        return None;
    }
    let t = v.text().ok()?.to_lowercase();
    if t.len() >= 2 && "omitnan".starts_with(&t) {
        Some(true)
    } else if t.len() >= 2 && "includenan".starts_with(&t) {
        Some(false)
    } else {
        None
    }
}

const BAD_OPTION: &str = "Invalid option. Option must be 'double', 'native', 'default', 'omitnan', 'includenan', 'omitmissing' or 'includemissing'.";

// ---- geomean, harmmean ------------------------------------------------------------------------

fn means(name: &str, a: &[Value]) -> Ret {
    let x = data_arg(a, 0)?;
    let mut sel = DimSel::Default;
    let mut omit = false;
    for v in &a[1..] {
        if let Some(f) = nan_flag(v) {
            omit = f;
        } else if v.is_text() && !matches!(parse_dim(v), Ok(DimSel::All)) {
            return rt(BAD_OPTION);
        } else {
            sel = parse_dim(v)?;
        }
    }
    let geo = name == "geomean";
    if geo {
        if matches!(x.kind, Kind::Char | Kind::Logical | Kind::Int(_)) {
            return rt(format!("Undefined function 'log' for input arguments of type '{}'.", x.kind.name()));
        }
        if x.is_complex() || x.data.iter().any(|&v| v < 0.0) {
            return rt("X may not contain negative values or complex values.");
        }
    }
    let single = is_single(x);
    if sel == DimSel::Default && is_empty_00(x) {
        return one(Value::Num(finish(Matrix::scalar(f64::NAN), single)));
    }
    let r = apply(x, &sel, |v| {
        let vals: Vec<f64> = if omit { v.iter().copied().filter(|t| !t.is_nan()).collect() } else { v.to_vec() };
        if vals.is_empty() {
            return f64::NAN;
        }
        let n = vals.len() as f64;
        if geo {
            (vals.iter().map(|t| t.ln()).sum::<f64>() / n).exp()
        } else {
            match x.kind {
                Kind::Int(_) => n / vals.iter().map(|t| x.kind.cast(1.0 / t)).sum::<f64>(),
                _ => n / vals.iter().map(|t| 1.0 / t).sum::<f64>(),
            }
        }
    });
    one(Value::Num(finish(r, single)))
}

// ---- trimmean ---------------------------------------------------------------------------------

fn trimmean(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("TRIMMEAN requires two input arguments.");
    }
    let x = data_arg(a, 0)?;
    let p = match &a[1] {
        Value::Num(m) if m.numel() == 1 && !m.data[0].is_nan() && m.data[0] >= 0.0 && m.data[0] < 100.0 => m.data[0],
        _ => return rt("PERCENT must be a scalar between 0 and 100."),
    };
    let mut flag = "round".to_string();
    let mut sel = DimSel::Default;
    let is_dim = |v: &Value| match v {
        Value::Num(m) if m.kind != Kind::Char => true,
        _ => v.is_text() && matches!(parse_dim(v), Ok(DimSel::All)),
    };
    for v in &a[2..a.len().min(4)] {
        if is_dim(v) {
            sel = parse_dim(v)?;
        } else if v.is_text() {
            flag = v.text()?.to_lowercase();
        } else {
            return rt("FLAG must be 'round', 'floor', or 'weighted'.");
        }
    }
    if !["round", "floor", "weighted"].contains(&flag.as_str()) {
        return rt("FLAG must be 'round', 'floor', or 'weighted'.");
    }
    let single = is_single(x);
    if sel == DimSel::Default && is_empty_00(x) {
        return one(Value::Num(finish(Matrix::scalar(f64::NAN), single)));
    }
    let r = apply(x, &sel, |v| {
        let mut s: Vec<f64> = v.iter().copied().filter(|t| !t.is_nan()).collect();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = s.len();
        if n == 0 {
            return f64::NAN;
        }
        let nf = n as f64;
        let k = nf * p / 200.0;
        match flag.as_str() {
            "weighted" => {
                let k0 = k.floor();
                let f = 1.0 + k0 - k;
                if !(k0 < nf / 2.0 || f > 0.0) {
                    return f64::NAN;
                }
                let k0u = k0 as usize;
                let mut sum = 0.0;
                if n >= k0u + 2 {
                    for i in (k0u + 1)..(n - k0u - 1).max(k0u + 1) {
                        sum += s[i];
                    }
                }
                let w = (nf - 2.0 * k0 - 2.0).max(0.0) + 2.0 * f;
                (sum + f * s[k0u] + f * s[n - k0u - 1]) / w
            }
            other => {
                let k0 = if other == "round" { (k - k * f64::EPSILON).round() } else { k.floor() } as usize;
                if !((k0 as f64) < nf / 2.0) {
                    return f64::NAN;
                }
                mean_of(&s[k0..n - k0])
            }
        }
    });
    one(Value::Num(finish(r, single)))
}

// ---- skewness, kurtosis -----------------------------------------------------------------------

fn flag_arg(v: Option<&Value>, default: f64) -> Res<f64> {
    match v {
        None => Ok(default),
        Some(Value::Num(m)) if m.is_empty() => Ok(default),
        Some(Value::Num(m)) if m.numel() == 1 && (m.data[0] == 0.0 || m.data[0] == 1.0) => Ok(m.data[0]),
        _ => rt("Flag must be 0 or 1."),
    }
}

fn shape(name: &str, a: &[Value]) -> Ret {
    let x = data_arg(a, 0)?;
    if matches!(x.kind, Kind::Int(_)) {
        return rt("Invalid data type. First input argument must be single or double.");
    }
    let flag = flag_arg(a.get(1), 1.0)?;
    let sel = match a.get(2) {
        Some(v) => parse_dim(v)?,
        None => DimSel::Default,
    };
    let single = is_single(x);
    if sel == DimSel::Default && is_empty_00(x) {
        return one(Value::Num(finish(Matrix::scalar(f64::NAN), single)));
    }
    let skew = name == "skewness";
    let r = apply(x, &sel, |v| {
        let vals: Vec<f64> = v.iter().copied().filter(|t| !t.is_nan()).collect();
        let n = vals.len() as f64;
        if vals.is_empty() {
            return f64::NAN;
        }
        let m = mean_of(&vals);
        let s2 = vals.iter().map(|t| (t - m).powi(2)).sum::<f64>() / n;
        if skew {
            let m3 = vals.iter().map(|t| (t - m).powi(3)).sum::<f64>() / n;
            let mut s = m3 / s2.powf(1.5);
            if flag == 0.0 {
                s = if n < 3.0 { f64::NAN } else { s * ((n - 1.0) / n).sqrt() * n / (n - 2.0) };
            }
            s
        } else {
            let m4 = vals.iter().map(|t| (t - m).powi(4)).sum::<f64>() / n;
            let mut k = m4 / (s2 * s2);
            if flag == 0.0 {
                k = if n < 4.0 { f64::NAN } else { ((n + 1.0) * k - 3.0 * (n - 1.0)) * (n - 1.0) / ((n - 2.0) * (n - 3.0)) + 3.0 };
            }
            k
        }
    });
    one(Value::Num(finish(r, single)))
}

// ---- moment -----------------------------------------------------------------------------------

fn moment(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Moment requires two inputs.");
    }
    if a.len() > 3 {
        return rt("Too many input arguments.");
    }
    let x = data_arg(a, 0)?;
    let order = match &a[1] {
        Value::Num(m) if m.numel() == 1 => m.data[0],
        _ => return rt("ORDER must be a scalar."),
    };
    let mut sel = DimSel::Default;
    let mut omit = false;
    if let Some(v) = a.get(2) {
        if let Some(f) = nan_flag(v) {
            omit = f;
        } else if v.is_text() {
            return rt(BAD_OPTION);
        } else {
            sel = parse_dim(v)?;
        }
    }
    let single = is_single(x);
    if sel == DimSel::Default && is_empty_00(x) {
        return one(Value::Num(finish(Matrix::scalar(f64::NAN), single)));
    }
    let g = groups(&x.dims, &sel);
    let (mut re, mut im) = (Vec::new(), Vec::new());
    let mut complex = false;
    for ix in &g.idx {
        let vals: Vec<f64> = ix.iter().map(|&i| x.data[i]).filter(|t| !omit || !t.is_nan()).collect();
        if vals.is_empty() {
            re.push(f64::NAN);
            im.push(0.0);
            continue;
        }
        let m = mean_of(&vals);
        let n = vals.len() as f64;
        let (mut sr, mut si) = (0.0, 0.0);
        for t in &vals {
            let d = t - m;
            if d < 0.0 && order.fract() != 0.0 {
                let r = (-d).powf(order);
                sr += r * (std::f64::consts::PI * order).cos();
                si += r * (std::f64::consts::PI * order).sin();
                complex = true;
            } else {
                sr += d.powf(order);
            }
        }
        re.push(sr / n);
        im.push(si / n);
    }
    let r = if complex { Matrix::complex(g.dims, re, im) } else { Matrix::from_dims(g.dims, re) };
    one(Value::Num(finish(r, single)))
}

// ---- range ------------------------------------------------------------------------------------

fn range(a: &[Value]) -> Ret {
    let x = data_arg(a, 0)?;
    let sel = match a.get(1) {
        Some(v) => parse_dim(v)?,
        None => DimSel::Default,
    };
    if is_empty_00(x) && sel == DimSel::Default {
        return one(Value::Num(Matrix::empty()));
    }
    let g = groups(&x.dims, &sel);
    let (mut re, mut im) = (Vec::new(), Vec::new());
    for ix in &g.idx {
        let mut best_max: Option<(f64, f64)> = None;
        let mut best_min: Option<(f64, f64)> = None;
        for &i in ix {
            let (r, c) = (x.data[i], x.im_at(i));
            if r.is_nan() || c.is_nan() {
                continue;
            }
            let key = |r: f64, c: f64| if x.is_complex() { r.hypot(c) } else { r };
            let ab = key(r, c);
            if best_max.map_or(true, |(br, bc)| ab > key(br, bc)) {
                best_max = Some((r, c));
            }
            if best_min.map_or(true, |(br, bc)| ab < key(br, bc)) {
                best_min = Some((r, c));
            }
        }
        match (best_max, best_min) {
            (Some(mx), Some(mn)) => {
                re.push(mx.0 - mn.0);
                im.push(mx.1 - mn.1);
            }
            _ => {
                re.push(f64::NAN);
                im.push(0.0);
            }
        }
    }
    let mut r = if x.is_complex() { Matrix::complex(g.dims, re, im) } else { Matrix::from_dims(g.dims, re) };
    if let Kind::Int(_) | Kind::Single = x.kind {
        r = r.cast_to(x.kind);
    }
    one(Value::Num(r))
}

// ---- mad --------------------------------------------------------------------------------------

fn mad(a: &[Value]) -> Ret {
    if a.len() > 3 {
        return rt("Too many input arguments.");
    }
    let x = data_arg(a, 0)?;
    let flag = match a.get(1) {
        None => 0.0,
        Some(Value::Num(m)) if m.is_empty() => 0.0,
        Some(Value::Num(m)) if m.numel() == 1 && (m.data[0] == 0.0 || m.data[0] == 1.0) => m.data[0],
        _ => return rt("Flag must be 0 or 1."),
    };
    let sel = match a.get(2) {
        Some(v) => parse_dim(v)?,
        None => DimSel::Default,
    };
    let single = is_single(x);
    if sel == DimSel::Default && is_empty_00(x) {
        return one(Value::Num(finish(Matrix::scalar(f64::NAN), single)));
    }
    let r = apply(x, &sel, |v| {
        let mut vals: Vec<f64> = v.iter().copied().filter(|t| !t.is_nan()).collect();
        if vals.is_empty() {
            return f64::NAN;
        }
        if flag == 1.0 {
            let m = median_of(&mut vals.clone());
            let mut dev: Vec<f64> = vals.iter().map(|t| (t - m).abs()).collect();
            median_of(&mut dev)
        } else {
            vals.shrink_to_fit();
            let m = mean_of(&vals);
            mean_of(&vals.iter().map(|t| (t - m).abs()).collect::<Vec<_>>())
        }
    });
    one(Value::Num(finish(r, single)))
}

// ---- zscore -----------------------------------------------------------------------------------

fn zscore(a: &[Value], nargout: usize) -> Ret {
    if a.len() > 3 {
        return rt("Too many input arguments.");
    }
    let x = data_arg(a, 0)?;
    if matches!(x.kind, Kind::Int(_)) {
        return rt("Invalid data type. First input argument must be single or double.");
    }
    let flag = flag_arg(a.get(1), 0.0)?;
    let sel = match a.get(2) {
        Some(v) => parse_dim(v)?,
        None => DimSel::Default,
    };
    let single = is_single(x);
    if x.is_empty() {
        let mut z = Matrix::from_dims(x.dims.clone(), Vec::new());
        z.kind = x.kind;
        return Ok(Some(vec![Value::Num(z), Value::Num(Matrix::scalar(f64::NAN)), Value::Num(Matrix::scalar(f64::NAN))].into_iter().take(nargout.max(1)).collect()));
    }
    let g = groups(&x.dims, &sel);
    let mut mu = Vec::new();
    let mut sigma = Vec::new();
    let mut z = vec![0.0; x.data.len()];
    for ix in &g.idx {
        let vals: Vec<f64> = ix.iter().map(|&i| x.data[i]).collect();
        let n = vals.len() as f64;
        let m = mean_of(&vals);
        let s = if vals.len() == 1 {
            0.0
        } else {
            let ss = vals.iter().map(|t| (t - m) * (t - m)).sum::<f64>();
            (ss / if flag == 1.0 { n } else { n - 1.0 }).sqrt()
        };
        let div = if s == 0.0 { 1.0 } else { s };
        for (&i, &t) in ix.iter().zip(&vals) {
            z[i] = (t - m) / div;
        }
        mu.push(m);
        sigma.push(s);
    }
    let zm = finish(Matrix::from_dims(x.dims.clone(), z), single);
    let mum = finish(Matrix::from_dims(g.dims.clone(), mu), single);
    let sm = finish(Matrix::from_dims(g.dims, sigma), single);
    Ok(Some(vec![Value::Num(zm), Value::Num(mum), Value::Num(sm)].into_iter().take(nargout.max(1)).collect()))
}

// ---- tiedrank ---------------------------------------------------------------------------------

/// Ranks of one vector with ties averaged; the second value is the tie adjustment.
pub fn tiedrank_vec(x: &[f64], tieflag: bool, bidir: bool) -> (Vec<f64>, Vec<f64>) {
    let mut order: Vec<usize> = (0..x.len()).collect();
    order.sort_by(|&i, &j| match (x[i].is_nan(), x[j].is_nan()) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => x[i].partial_cmp(&x[j]).unwrap(),
    });
    let nn = x.iter().filter(|t| t.is_nan()).count();
    let len = x.len() - nn;
    let mut ranks: Vec<f64> = if bidir {
        let half = len.div_ceil(2);
        let mut r: Vec<f64> = (1..=half).map(|v| v as f64).collect();
        let down = len / 2;
        r.extend((1..=down).rev().map(|v| v as f64));
        r
    } else {
        (1..=len).map(|v| v as f64).collect()
    };
    let mut adj = if tieflag { vec![0.0; 3] } else { vec![0.0] };
    let mut i = 0;
    while i < len {
        let mut j = i;
        while j + 1 < len && x[order[j + 1]] == x[order[i]] {
            j += 1;
        }
        let t = (j - i + 1) as f64;
        if j > i {
            if tieflag {
                let n2 = t * (t - 1.0);
                adj[0] += n2 / 2.0;
                adj[1] += n2 * (t - 2.0);
                adj[2] += n2 * (2.0 * t + 5.0);
            } else {
                adj[0] += t * (t - 1.0) * (t + 1.0) / 2.0;
            }
            let avg = ranks[i..=j].iter().sum::<f64>() / t;
            for r in &mut ranks[i..=j] {
                *r = avg;
            }
        }
        i = j + 1;
    }
    let mut out = vec![f64::NAN; x.len()];
    for (k, &o) in order.iter().enumerate() {
        out[o] = if k < len { ranks[k] } else { f64::NAN };
    }
    (out, adj)
}

pub fn tiedrank(a: &[Value], nargout: usize) -> Ret {
    let x = match a.first() {
        Some(Value::Num(m)) if m.kind != Kind::Char => m,
        Some(_) => return rt("Argument following 'like' must be a variable with a data type that supports ZEROS, such as a double or single variable."),
        None => return rt("Not enough input arguments."),
    };
    let truthy = |v: Option<&Value>| matches!(v, Some(Value::Num(m)) if m.data.first().is_some_and(|&t| t != 0.0));
    let tieflag = truthy(a.get(1));
    let bidir = truthy(a.get(2));
    let single = is_single(x);
    let is_vec = x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1);
    let (r, tadj) = if is_vec {
        let (r, adj) = tiedrank_vec(&x.data, tieflag, bidir);
        (Matrix::from_dims(x.dims.clone(), r), Matrix::from_dims(vec![adj.len(), 1], adj))
    } else {
        let rows = x.dims[0];
        let ncols: usize = x.dims[1..].iter().product();
        let mut r = vec![0.0; x.data.len()];
        let mut adj_all = Vec::new();
        for c in 0..ncols {
            let (rv, adj) = if rows == 0 { (vec![], vec![0.0; if tieflag { 3 } else { 1 }]) } else { tiedrank_vec(&x.data[c * rows..(c + 1) * rows], tieflag, bidir) };
            r[c * rows..(c + 1) * rows].copy_from_slice(&rv);
            adj_all.extend(adj);
        }
        let mut ad = vec![if tieflag { 3 } else { 1 }];
        ad.extend_from_slice(&x.dims[1..]);
        (Matrix::from_dims(x.dims.clone(), r), Matrix::from_dims(crate::nd::normalize(ad), adj_all))
    };
    let mut outs = vec![Value::Num(finish(r, single))];
    if nargout > 1 {
        outs.push(Value::Num(finish(tadj, single)));
    }
    Ok(Some(outs))
}

// ---- nan* reductions --------------------------------------------------------------------------

fn nanextreme(name: &str, a: &[Value], nargout: usize) -> Ret {
    let want_max = name == "nanmax";
    let x = data_arg(a, 0)?;
    let better = |t: f64, best: f64| if want_max { t > best } else { t < best };
    let pick = |p: f64, q: f64| match (p.is_nan(), q.is_nan()) {
        (true, _) => q,
        (_, true) => p,
        _ => if want_max { p.max(q) } else { p.min(q) },
    };
    if a.len() >= 2 && !matches!(a.get(1), Some(Value::Num(m)) if m.is_empty()) {
        // Elementwise: nanmax(x, y)
        let y = data_arg(a, 1)?;
        let (dims, pairs) = crate::nd::broadcast_pairs(&x.dims, &y.dims)?;
        let data = pairs.iter().map(|&(i, j)| pick(x.data[i], y.data[j])).collect();
        let single = is_single(x) || is_single(y);
        return one(Value::Num(finish(Matrix::from_dims(dims, data), single)));
    }
    let sel = match a.get(2) {
        Some(v) => parse_dim(v)?,
        None => DimSel::Default,
    };
    let g = groups(&x.dims, &sel);
    let (mut vals, mut idxs) = (Vec::new(), Vec::new());
    for ix in &g.idx {
        let mut best = f64::NAN;
        let mut at = 1usize;
        for (k, &i) in ix.iter().enumerate() {
            let t = x.data[i];
            if t.is_nan() {
                continue;
            }
            if best.is_nan() || better(t, best) {
                best = t;
                at = k + 1;
            }
        }
        vals.push(best);
        idxs.push(if best.is_nan() { 1.0 } else { at as f64 });
    }
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(g.dims.clone(), vals), is_single(x)))];
    if nargout > 1 {
        outs.push(Value::Num(Matrix::from_dims(g.dims, idxs)));
    }
    Ok(Some(outs))
}

fn nanreduce(name: &str, a: &[Value]) -> Ret {
    let x = data_arg(a, 0)?;
    let std_like = matches!(name, "nanstd" | "nanvar");
    let (flag, dim_at) = if std_like { (flag_arg(a.get(1), 0.0)?, 2) } else { (0.0, 1) };
    let sel = match a.get(dim_at) {
        Some(v) => parse_dim(v)?,
        None => DimSel::Default,
    };
    let single = is_single(x);
    let r = apply(x, &sel, |v| {
        let mut vals: Vec<f64> = v.iter().copied().filter(|t| !t.is_nan()).collect();
        let n = vals.len() as f64;
        match name {
            "nansum" => vals.iter().sum(),
            "nanmean" => mean_of(&vals),
            "nanmedian" => median_of(&mut vals),
            _ => {
                if vals.is_empty() {
                    return f64::NAN;
                }
                if vals.len() == 1 {
                    return 0.0;
                }
                let m = mean_of(&vals);
                let var = vals.iter().map(|t| (t - m) * (t - m)).sum::<f64>() / if flag == 1.0 { n } else { n - 1.0 };
                if name == "nanstd" { var.sqrt() } else { var }
            }
        }
    });
    one(Value::Num(finish(r, single)))
}
