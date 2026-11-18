//! Nonparametric tests: `signtest`, `signrank`, `ranksum`.

use super::dist::{binocdf, normcdf_std};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["signtest", "signrank", "ranksum"];

struct Opts {
    alpha: f64,
    method: Option<String>,
    tail: usize,
}

fn choose(v: &Value, valid: &[&str], pname: &str) -> Res<usize> {
    let t = if v.is_text() { v.text()? } else { String::new() };
    let low = t.to_lowercase();
    let hits: Vec<usize> = valid.iter().enumerate().filter(|(_, c)| !low.is_empty() && c.starts_with(&low)).map(|(i, _)| i).collect();
    if hits.len() == 1 {
        return Ok(hits[0]);
    }
    let list = valid.iter().map(|c| format!("'{c}'")).collect::<Vec<_>>().join(", ");
    rt(format!("'{t}' is not a valid value for the '{pname}' argument.  Valid values are: {list}."))
}

fn parse(rest: &[Value], who: &str, allow_old: bool) -> Res<Opts> {
    let mut rest = rest;
    let mut o = Opts { alpha: 0.05, method: None, tail: 0 };
    let bad_alpha = || -> Res<()> {
        if who == "ranksum" { rt("RANKSUM requires a scalar ALPHA value between 0 and 1.") } else { rt("ALPHA must be a scalar value strictly between zero and one.") }
    };
    if let Some(Value::Num(m)) = rest.first() {
        if m.kind != Kind::Char {
            if m.numel() != 1 || m.data[0].is_nan() || m.data[0] <= 0.0 || m.data[0] >= 1.0 {
                bad_alpha()?;
            }
            o.alpha = m.data[0];
            rest = &rest[1..];
        }
    }
    if rest.len() % 2 == 1 {
        return rt("Invalid parameter name or missing value.");
    }
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        let hit = ["alpha", "method", "tail"].iter().position(|n| !key.is_empty() && n.starts_with(&key));
        match hit {
            Some(0) => match &pair[1] {
                Value::Num(m) if m.numel() == 1 && m.data[0] > 0.0 && m.data[0] < 1.0 => o.alpha = m.data[0],
                _ => bad_alpha()?,
            },
            Some(1) => {
                let t = pair[1].text().unwrap_or_default().to_lowercase();
                if allow_old && t == "oldexact" {
                    o.method = Some("oldexact".into());
                } else {
                    let i = choose(&pair[1], &["exact", "approximate"], "method")?;
                    o.method = Some(["exact", "approximate"][i].into());
                }
            }
            Some(2) => o.tail = choose(&pair[1], &["both", "right", "left"], "tail")?,
            _ => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
        }
    }
    Ok(o)
}

fn vector<'a>(v: Option<&'a Value>, who: &str) -> Res<&'a Matrix> {
    match v {
        Some(Value::Num(m)) if m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1) && !m.is_empty() => Ok(m),
        Some(Value::Num(_)) | None => rt(format!("{who} requires vector rather than matrix data.")),
        _ => rt(format!("{who} requires vector rather than matrix data.")),
    }
}

fn paired_inputs(a: &[Value], who: &str) -> Res<(Vec<f64>, Vec<f64>, bool, Vec<f64>)> {
    let x = vector(a.first(), who)?;
    let (y, onesample): (Vec<f64>, bool) = match a.get(1) {
        None => (vec![0.0; x.numel()], true),
        Some(Value::Num(m)) if m.is_empty() => (vec![0.0; x.numel()], true),
        Some(Value::Num(m)) if m.numel() == 1 && m.kind != Kind::Char => (vec![m.data[0]; x.numel()], false),
        Some(Value::Num(m)) => {
            vector(a.get(1), who)?;
            if m.numel() != x.numel() {
                return rt(format!("{who} requires the data vectors to have the same number of elements."));
            }
            (m.data.clone(), false)
        }
        Some(_) => return rt(format!("{who} requires vector rather than matrix data.")),
    };
    let eps: Vec<f64> = if onesample { vec![0.0; x.numel()] } else { x.data.iter().zip(&y).map(|(a, b)| eps_of(*a) + eps_of(*b)).collect() };
    Ok((x.data.clone(), y, onesample, eps))
}

fn signtest(a: &[Value], nargout: usize) -> Ret {
    let (x, y, _, _) = paired_inputs(a, "SIGNTEST")?;
    let single = false;
    let o = parse(&a[2.min(a.len())..], "signtest", false)?;
    let diff: Vec<f64> = x.iter().zip(&y).map(|(p, q)| p - q).filter(|d| *d != 0.0 && !d.is_nan()).collect();
    let n = diff.len();
    let mut st = StructArray::scalar();
    let (p, h);
    if n == 0 {
        p = 1.0;
        h = 0.0;
        st.set("sign", Value::num(0.0));
        st.set("zval", Value::num(f64::NAN));
    } else {
        let method = o.method.clone().unwrap_or_else(|| if n < 100 { "exact".into() } else { "approximate".into() });
        let npos = diff.iter().filter(|&&d| d > 0.0).count() as f64;
        let nneg = n as f64 - npos;
        let nf = n as f64;
        let z;
        if method == "exact" {
            p = match o.tail {
                0 => (2.0 * binocdf(npos.min(nneg), nf, 0.5)).min(1.0),
                1 => binocdf(nneg, nf, 0.5),
                _ => binocdf(npos, nf, 0.5),
            };
            z = f64::NAN;
        } else {
            let d = npos - nneg;
            match o.tail {
                0 => {
                    z = (d - d.signum() * if d == 0.0 { 0.0 } else { 1.0 }) / nf.sqrt();
                    p = 2.0 * normcdf_std(-z.abs());
                }
                1 => {
                    z = (d - 1.0) / nf.sqrt();
                    p = normcdf_std(-z);
                }
                _ => {
                    z = (d + 1.0) / nf.sqrt();
                    p = normcdf_std(z);
                }
            }
        }
        h = if p <= o.alpha { 1.0 } else { 0.0 };
        st.set("zval", Value::num(z));
        st.set("sign", Value::num(npos));
    }
    let mut outs = vec![Value::Num(finish(Matrix::scalar(p), single)), Value::Num(Matrix::logical(1, 1, vec![h]))];
    if n == 0 {
        outs[1] = Value::num(h);
    }
    if nargout > 2 {
        outs.push(Value::Struct(st));
    }
    Ok(Some(outs))
}

/// Tied ranks with a tolerance: values within `eps` of the next one are tied.
fn ranks_with_eps(v: &[f64], eps: &[f64]) -> (Vec<f64>, f64) {
    let n = v.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| v[i].partial_cmp(&v[j]).unwrap());
    let mut ranks = vec![0.0; n];
    let mut adj = 0.0;
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && v[order[j]] + eps[order[j]] >= v[order[j + 1]] - eps[order[j + 1]] {
            j += 1;
        }
        let t = (j - i + 1) as f64;
        if j > i {
            adj += t * (t - 1.0) * (t + 1.0) / 2.0;
        }
        let avg = ((i + 1 + j + 1) as f64) / 2.0;
        for k in i..=j {
            ranks[order[k]] = avg;
        }
        i = j + 1;
    }
    (ranks, adj)
}

/// Distribution of the sum of a random subset (each element in with probability 1/2) of the
/// doubled ranks: probabilities indexed by the doubled sum.
fn subset_sum_dist(r2: &[usize]) -> Vec<f64> {
    let total: usize = r2.iter().sum();
    let mut cnt = vec![0.0; total + 1];
    cnt[0] = 1.0;
    for &r in r2 {
        for s in (r..=total).rev() {
            cnt[s] += cnt[s - r];
        }
    }
    let scale = 0.5f64.powi(r2.len() as i32);
    cnt.iter().map(|c| c * scale).collect()
}

fn signrank(a: &[Value], nargout: usize) -> Ret {
    let (x, y, onesample, mut eps) = paired_inputs(a, "SIGNRANK")?;
    let is32 = matches!(a.first(), Some(Value::Num(m)) if is_single(m)) || matches!(a.get(1), Some(Value::Num(m)) if is_single(m));
    if is32 && !onesample {
        let e32 = |v: f64| -> f64 {
            let a = (v as f32).abs();
            if a == 0.0 { f32::from_bits(1) as f64 } else { (f32::EPSILON * 2f32.powi(a.log2().floor() as i32)) as f64 }
        };
        eps = x.iter().zip(&y).map(|(p, q)| e32(*p) + e32(*q)).collect();
    }
    let single = false;
    let o = parse(&a[2.min(a.len())..], "signrank", true)?;
    let mut diff = Vec::new();
    let mut epsd = Vec::new();
    for ((p, q), e) in x.iter().zip(&y).zip(&eps) {
        let d = if is32 { (*p as f32 - *q as f32) as f64 } else { p - q };
        if !d.is_nan() {
            diff.push(d);
            epsd.push(*e);
        }
    }
    if diff.is_empty() {
        return rt("SIGNRANK requires at least one non-NaN value.");
    }
    let keep: Vec<usize> = (0..diff.len()).filter(|&i| diff[i].abs() > epsd[i]).collect();
    let diff: Vec<f64> = keep.iter().map(|&i| diff[i]).collect();
    let epsd: Vec<f64> = keep.iter().map(|&i| epsd[i]).collect();
    let n = diff.len();
    let _ = onesample;
    let mut st = StructArray::scalar();
    if n == 0 {
        st.set("signedrank", Value::num(0.0));
        let mut outs = vec![Value::num(1.0), Value::num(0.0)];
        if nargout > 2 {
            outs.push(Value::Struct(st));
        }
        return Ok(Some(outs));
    }
    let method = o.method.clone().unwrap_or_else(|| if n <= 15 { "exact".into() } else { "approximate".into() });
    let absd: Vec<f64> = diff.iter().map(|d| d.abs()).collect();
    let (ranks, tieadj) = ranks_with_eps(&absd, &epsd);
    let w: f64 = (0..n).filter(|&i| diff[i] > 0.0).map(|i| ranks[i]).sum();
    let nf = n as f64;
    let w2 = nf * (nf + 1.0) / 2.0 - w;
    let p;
    let mut zval: Option<f64> = None;
    match method.as_str() {
        "approximate" => {
            let sd = ((nf * (nf + 1.0) * (2.0 * nf + 1.0) - tieadj) / 24.0).sqrt();
            let c = w - nf * (nf + 1.0) / 4.0;
            let z = match o.tail {
                0 => c / sd,
                1 => (c - 0.5) / sd,
                _ => (c + 0.5) / sd,
            };
            p = match o.tail {
                0 => 2.0 * normcdf_std(-z.abs()),
                1 => normcdf_std(-z),
                _ => normcdf_std(z),
            };
            zval = Some(z);
        }
        "oldexact" => {
            let r2: Vec<usize> = (1..=n).map(|i| 2 * i).collect();
            let dist = subset_sum_dist(&r2);
            let le = |v: f64| -> f64 { dist.iter().enumerate().filter(|(s, _)| (*s as f64) / 2.0 <= v + 1e-9).map(|(_, p)| *p).sum() };
            let ge = |v: f64| -> f64 { dist.iter().enumerate().filter(|(s, _)| (*s as f64) / 2.0 >= v - 1e-9).map(|(_, p)| *p).sum() };
            p = match o.tail {
                0 => (le(w.min(w2)) + ge(w.max(w2))).min(1.0),
                1 => ge(w),
                _ => le(w),
            };
        }
        _ => {
            let r2: Vec<usize> = ranks.iter().map(|r| (2.0 * r).round() as usize).collect();
            let dist = subset_sum_dist(&r2);
            let w2d = (2.0 * w).round() as usize;
            let w2d_other = (2.0 * w2).round() as usize;
            let le = |s: usize| -> f64 { dist.iter().take(s + 1).sum() };
            let ge = |s: usize| -> f64 { dist.iter().skip(s).sum() };
            let small = w2d.min(w2d_other);
            p = match o.tail {
                0 => (2.0 * le(small)).min(1.0),
                1 => ge(w2d),
                _ => le(w2d),
            };
        }
    }
    let h = if p <= o.alpha { 1.0 } else { 0.0 };
    st.set("signedrank", Value::Num(finish(Matrix::scalar(w), is32)));
    if let Some(z) = zval {
        st.set("zval", Value::num(z));
    }
    let mut outs = vec![Value::Num(finish(Matrix::scalar(p), single)), Value::Num(Matrix::logical(1, 1, vec![h]))];
    if nargout > 2 {
        outs.push(Value::Struct(st));
    }
    Ok(Some(outs))
}

fn ranksum(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let x = vector(a.first(), "RANKSUM")?;
    let y = vector(a.get(1), "RANKSUM")?;
    let single = is_single(x) || is_single(y);
    let xv: Vec<f64> = x.data.iter().copied().filter(|v| !v.is_nan()).collect();
    let yv: Vec<f64> = y.data.iter().copied().filter(|v| !v.is_nan()).collect();
    if xv.is_empty() || yv.is_empty() {
        return rt("SIGNRANK requires at least one non-NaN value in each sample.");
    }
    let o = parse(&a[2..], "ranksum", true)?;
    let (nx, ny) = (xv.len(), yv.len());
    let ns = nx.min(ny);
    let method = match o.method.as_deref() {
        None => if ns < 10 && nx + ny < 20 { "exact" } else { "approximate" },
        Some(m) => if m == "oldexact" { "exact" } else { m },
    };
    let same_order = nx <= ny;
    let (sm, lg) = if same_order { (&xv, &yv) } else { (&yv, &xv) };
    let mut all: Vec<f64> = sm.clone();
    all.extend(lg.iter());
    let (ranks, tieadj) = super::descr::tiedrank_vec(&all, false, false);
    let w: f64 = ranks[..ns].iter().sum();
    let (nxf, nyf) = (nx as f64, ny as f64);
    let mut st = StructArray::scalar();
    let p;
    let mut zval: Option<f64> = None;
    if method == "approximate" {
        let wmean = ns as f64 * (nxf + nyf + 1.0) / 2.0;
        let tiescor = 2.0 * tieadj[0] / ((nxf + nyf) * (nxf + nyf - 1.0));
        let wvar = nxf * nyf * ((nxf + nyf + 1.0) - tiescor) / 12.0;
        let wc = w - wmean;
        let z;
        match o.tail {
            0 => {
                let mut zz = (wc - 0.5 * if wc == 0.0 { 0.0 } else { wc.signum() }) / wvar.sqrt();
                if !same_order {
                    zz = -zz;
                }
                z = zz;
                p = 2.0 * normcdf_std(-z.abs());
            }
            1 => {
                z = if same_order { (wc - 0.5) / wvar.sqrt() } else { -(wc + 0.5) / wvar.sqrt() };
                p = normcdf_std(-z);
            }
            _ => {
                z = if same_order { (wc + 0.5) / wvar.sqrt() } else { -(wc - 0.5) / wvar.sqrt() };
                p = normcdf_std(z);
            }
        }
        zval = Some(z);
    } else {
        // exact distribution of the rank sum of a random subset of size ns
        let r2: Vec<usize> = ranks.iter().map(|r| (2.0 * r).round() as usize).collect();
        let total: usize = r2.iter().sum();
        let big_n = r2.len();
        // cnt[j][s]: number of ways to pick j elements with doubled sum s
        let mut cnt = vec![vec![0.0f64; total + 1]; ns + 1];
        cnt[0][0] = 1.0;
        for &r in &r2 {
            for j in (1..=ns).rev() {
                for s in (r..=total).rev() {
                    let add = cnt[j - 1][s - r];
                    cnt[j][s] += add;
                }
            }
        }
        let denom: f64 = cnt[ns].iter().sum();
        let _ = big_n;
        let wd = (2.0 * w).round() as usize;
        let plo: f64 = cnt[ns][..=wd].iter().sum::<f64>() / denom;
        let phi: f64 = cnt[ns][wd..].iter().sum::<f64>() / denom;
        p = match o.tail {
            0 => (2.0 * plo.min(phi)).min(1.0),
            1 => if same_order { phi } else { plo },
            _ => if same_order { plo } else { phi },
        };
    }
    let h = if p <= o.alpha { 1.0 } else { 0.0 };
    let rs = if same_order { w } else { ranks[ns..].iter().sum() };
    st.set("ranksum", Value::Num(finish(Matrix::scalar(rs), single)));
    if let Some(z) = zval {
        st.set("zval", Value::Num(finish(Matrix::scalar(z), single)));
    }
    let mut outs = vec![Value::Num(finish(Matrix::scalar(p), single)), Value::Num(Matrix::logical(1, 1, vec![h]))];
    if nargout > 2 {
        outs.push(Value::Struct(st));
    }
    Ok(Some(outs))
}

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "signtest" => signtest(a, nargout),
        "signrank" => signrank(a, nargout),
        _ => ranksum(a, nargout),
    }
}
