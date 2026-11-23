//! `corr`: Pearson, Spearman and Kendall correlation of the columns of one or two matrices, with
//! the `Rows` and `Tail` options and p-values (exact permutation distributions for small samples,
//! the AS 89 approximation for Spearman and a normal approximation for Kendall otherwise).

use super::descr::tiedrank_vec;
use super::util::*;
use crate::builtins::Ret;
use crate::builtins_specfun::{betainc_pair, erfc};
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;
use std::collections::HashMap;

pub const NAMES: &[&str] = &["corr"];

#[derive(Clone, Copy, PartialEq)]
enum Kindc {
    Pearson,
    Kendall,
    Spearman,
}

#[derive(Clone, Copy, PartialEq)]
enum Rows {
    All,
    Complete,
    Pairwise,
}

#[derive(Clone, Copy, PartialEq)]
enum Tail {
    Both,
    Right,
    Left,
}

fn norm_cdf(z: f64) -> f64 {
    0.5 * erfc(-z / std::f64::consts::SQRT_2)
}

/// Student t cumulative distribution.
pub fn t_cdf(x: f64, v: f64) -> f64 {
    if x.is_nan() || v.is_nan() || v <= 0.0 {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return 1.0;
    }
    if x == f64::NEG_INFINITY {
        return 0.0;
    }
    if v == f64::INFINITY {
        return super::dist::normcdf_std(x);
    }
    let z = v / (v + x * x);
    let tail = 0.5 * betainc_pair(z, v / 2.0, 0.5).0;
    if x < 0.0 { tail } else { 1.0 - tail }
}

fn prefix_choice(s: &str, choices: &[&str]) -> Option<usize> {
    let s = s.to_lowercase();
    if s.is_empty() {
        return None;
    }
    let hits: Vec<usize> = choices.iter().enumerate().filter(|(_, c)| c.starts_with(&s)).map(|(i, _)| i).collect();
    if hits.len() == 1 { Some(hits[0]) } else { None }
}

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |a, k| a * k as f64)
}

/// Calls `f` with every permutation of `v` (Heap's algorithm).
fn each_permutation(v: &mut Vec<f64>, f: &mut impl FnMut(&[f64])) {
    let n = v.len();
    let mut c = vec![0usize; n];
    f(v);
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                v.swap(0, i);
            } else {
                v.swap(c[i], i);
            }
            f(v);
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
}

/// Cumulative counts (and counts from each bin up) over a frequency table.
fn tails(freq: &[f64], total: f64) -> (Vec<f64>, Vec<f64>) {
    let mut cum = Vec::with_capacity(freq.len());
    let mut s = 0.0;
    for &f in freq {
        s += f;
        cum.push(s);
    }
    let mut rc = vec![total; freq.len()];
    for k in 1..freq.len() {
        rc[k] = total - cum[k - 1];
    }
    (cum, rc)
}

/// AS 89 upper tail of Spearman's statistic for n > 9.
fn as89(d: f64, n: f64, n3const: f64) -> f64 {
    let c = [0.2274, 0.2531, 0.1745, 0.0758, 0.1033, 0.3932, 0.0879, 0.0151, 0.0072, 0.0831, 0.0131, 0.00046];
    let x = (2.0 * (d - 1.0) / n3const - 1.0) * (n - 1.0).sqrt();
    let y = x * x;
    let u = x
        * (c[0]
            + (c[1] + c[2] / n) / n
            + y * (-c[3] + (c[4] + c[5] / n) / n - y * (c[6] + c[7] / n - y * (c[8] - c[9] / n + y * (c[10] - c[11] * y) / n)) / n))
        / n;
    let p = u / (0.5 * y).exp() + 0.5 * erfc(x / std::f64::consts::SQRT_2);
    p.clamp(0.0, 1.0)
}

struct Cache {
    spearman: HashMap<usize, Vec<f64>>,
}

/// Spearman p-value without ties (D is the sum of squared rank differences).
fn spearman_p_noties(tail: Tail, d: f64, n: usize, cache: &mut Cache) -> f64 {
    if d.is_nan() {
        return f64::NAN;
    }
    let nf = n as f64;
    let n3const = (nf * nf * nf - nf) / 3.0;
    if n < 10 {
        let freq = cache.spearman.entry(n).or_insert_with(|| {
            let bins = (n3const * 2.0).round() as usize + 1;
            let mut freq = vec![0.0; bins];
            let mut v: Vec<f64> = (1..=n).map(|k| k as f64).collect();
            each_permutation(&mut v, &mut |p| {
                let dd: f64 = p.iter().enumerate().map(|(i, &q)| (i as f64 + 1.0 - q).powi(2)).sum();
                freq[(2.0 * dd).round() as usize] += 1.0;
            });
            freq
        });
        return exact_tail(tail, freq, (2.0 * d).round() as usize, factorial(n), true);
    }
    match tail {
        Tail::Both => (2.0 * as89(d.max(n3const - d), nf, n3const)).min(1.0),
        Tail::Right => as89(n3const - d, nf, n3const),
        Tail::Left => as89(d, nf, n3const),
    }
}

/// p-value from a frequency table: `at` is the bin of the observed statistic. For Spearman
/// the right tail is the lower cumulative; Kendall's is the other way round.
fn exact_tail(tail: Tail, freq: &[f64], at: usize, total: f64, spearman: bool) -> f64 {
    let (cum, rc) = tails(freq, total);
    let p = match tail {
        Tail::Both => (2.0 * cum[at].min(rc[at]) / total).min(1.0),
        Tail::Right => (if spearman { cum[at] } else { rc[at] }) / total,
        Tail::Left => (if spearman { rc[at] } else { cum[at] }) / total,
    };
    p
}

/// Spearman p-value with ties in the data.
fn spearman_p_ties(tail: Tail, d: f64, mean_d: f64, std_d: f64, xr: &[f64], yr: &[f64]) -> f64 {
    let n = xr.len();
    if std_d == 0.0 || d.is_nan() {
        return f64::NAN;
    }
    if n < 10 {
        let nf = n as f64;
        let n3const = (nf * nf * nf - nf) / 3.0;
        let bins = ((n3const + 0.5) / 0.5).ceil() as usize + 1;
        let mut freq = vec![0.0; bins];
        let mut v = yr.to_vec();
        each_permutation(&mut v, &mut |p| {
            let dd: f64 = xr.iter().zip(p).map(|(a, b)| (a - b).powi(2)).sum();
            let k = ((dd + 0.25) / 0.5).floor() as usize;
            if k < freq.len() {
                freq[k] += 1.0;
            }
        });
        return exact_tail(tail, &freq, (2.0 * d).round() as usize, factorial(n), true);
    }
    let nf = n as f64;
    let r = (mean_d - d) / ((nf - 1.0).sqrt() * std_d);
    let t = if r.abs() < 1.0 { r * ((nf - 2.0) / (1.0 - r * r)).sqrt() } else { f64::INFINITY * r.signum() };
    match tail {
        Tail::Both => 2.0 * t_cdf(-t.abs(), nf - 2.0),
        Tail::Right => t_cdf(-t, nf - 2.0),
        Tail::Left => t_cdf(t, nf - 2.0),
    }
}

/// Kendall p-value; `ranks` holds the tied ranks when there are ties.
fn kendall_p(tail: Tail, k: f64, std_k: f64, n: usize, ranks: Option<(&[f64], &[f64])>) -> f64 {
    let n2 = n * (n - 1) / 2;
    let noties = ranks.is_none();
    if !noties && std_k == 0.0 || k.is_nan() {
        return f64::NAN;
    }
    let exact = if noties { n < 50 } else { n < 10 };
    if exact {
        let mut freq = vec![0.0; 2 * n2 + 1];
        if noties {
            // Distribution of inversions: product of (1 + x + ... + x^(i-1)).
            let mut poly = vec![1.0, 1.0];
            for i in 3..=n {
                let mut next = vec![0.0; poly.len() + i - 1];
                for (a, &pa) in poly.iter().enumerate() {
                    for b in 0..i {
                        next[a + b] += pa;
                    }
                }
                poly = next;
            }
            for (inv, &c) in poly.iter().enumerate() {
                if inv <= n2 {
                    freq[2 * (n2 - inv)] = c;
                }
            }
        } else {
            let (xr, yr) = ranks.unwrap();
            let mut v = yr.to_vec();
            each_permutation(&mut v, &mut |p| {
                let mut kk = 0.0;
                for a in 0..n - 1 {
                    for b in a + 1..n {
                        kk += (xr[a] - xr[b]).signum_zero() * (p[a] - p[b]).signum_zero();
                    }
                }
                freq[(kk + n2 as f64).round() as usize] += 1.0;
            });
        }
        let at = (k + n2 as f64).round() as usize;
        return exact_tail(tail, &freq, at, factorial(n), false);
    }
    match tail {
        Tail::Both => (2.0 * norm_cdf(-(k.abs() - 1.0) / std_k)).min(1.0),
        Tail::Right => norm_cdf(-(k - 1.0) / std_k),
        Tail::Left => norm_cdf((k + 1.0) / std_k),
    }
}

trait SignZero {
    fn signum_zero(self) -> f64;
}
impl SignZero for f64 {
    fn signum_zero(self) -> f64 {
        if self > 0.0 {
            1.0
        } else if self < 0.0 {
            -1.0
        } else if self.is_nan() {
            f64::NAN
        } else {
            0.0
        }
    }
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (a, b) in xs.iter().zip(ys) {
        let (da, db) = (a - mx, b - my);
        sxy += da * db;
        sxx += da * da;
        syy += db * db;
    }
    let r = sxy / sxx.sqrt() / syy.sqrt();
    if r.abs() > 1.0 { r / r.abs() } else { r }
}

fn pearson_p(tail: Tail, r: f64, n: f64) -> f64 {
    let t = r * ((n - 2.0) / (1.0 - r * r)).sqrt();
    match tail {
        Tail::Both => 2.0 * t_cdf(-t.abs(), n - 2.0),
        Tail::Right => t_cdf(-t, n - 2.0),
        Tail::Left => t_cdf(t, n - 2.0),
    }
}

pub fn call(_it: &mut Interpreter, _name: &str, a: &[Value], nargout: usize) -> Ret {
    let x = match a.first() {
        Some(Value::Num(m)) if !m.is_empty() => m,
        _ => return rt("Requires a data matrix X."),
    };
    if x.dims.len() > 2 {
        return rt("TRANSPOSE does not support N-D arrays. Use PAGETRANSPOSE/PAGECTRANSPOSE to transpose pages or PERMUTE to reorder dimensions of N-D arrays.");
    }
    let (n0, p1) = (x.dims[0], x.dims[1]);
    let mut rest = &a[1..];
    let (y, xx) = match rest.first() {
        Some(v) if !v.is_text() => {
            let Value::Num(m) = v else { return rt("Requires a data matrix X.") };
            rest = &rest[1..];
            if m.dims[0] != n0 {
                return rt("X and Y must have the same number of rows.");
            }
            (m, false)
        }
        _ => (x, true),
    };
    let p2 = y.dims[1];
    let (mut ty, mut rows, mut tail) = ("p".to_string(), "a".to_string(), "both".to_string());
    let mut tail_given = false;
    let mut i = 0;
    while i < rest.len() {
        let key = rest[i].text().unwrap_or_default().to_lowercase();
        let val = rest.get(i + 1).ok_or_else(|| crate::error::LangError::new("Name-value pair arguments require a name followed by a value.", 0))?;
        let vt = val.text().unwrap_or_default();
        match prefix_choice(&key, &["type", "rows", "tail"]) {
            Some(0) => ty = vt,
            Some(1) => rows = vt,
            Some(2) => {
                tail = vt;
                tail_given = true;
            }
            _ => return rt(format!("'{}' is not a recognized parameter name.", rest[i].text().unwrap_or_default())),
        }
        i += 2;
    }
    let _ = tail_given;
    let rows_kind = match prefix_choice(&rows, &["all", "complete", "pairwise"]) {
        Some(0) => Rows::All,
        Some(1) => Rows::Complete,
        Some(2) => Rows::Pairwise,
        _ => return rt("The 'rows' parameter value must be 'all', 'complete', or 'pairwise'."),
    };
    let single = is_single(x) || is_single(y);
    let (mut xd, mut yd) = (x.data.clone(), y.data.clone());
    let mut n = n0;
    if rows_kind == Rows::Complete {
        let ok: Vec<bool> = (0..n0)
            .map(|r| (0..p1).all(|c| !x.data[c * n0 + r].is_nan()) && (xx || (0..p2).all(|c| !y.data[c * n0 + r].is_nan())))
            .collect();
        n = ok.iter().filter(|&&b| b).count();
        if n < n0 {
            let take = |m: &Matrix, p: usize| -> Vec<f64> { (0..p).flat_map(|c| (0..n0).filter(|&r| ok[r]).map(move |r| m.data[c * n0 + r])).collect() };
            xd = take(x, p1);
            yd = if xx { xd.clone() } else { take(y, p2) };
        }
    }
    if n < 2 {
        let nan = Matrix::filled(p1, p2, f64::NAN);
        let v = Value::Num(finish(nan, single));
        return Ok(Some(if nargout > 1 { vec![v.clone(), v] } else { vec![v] }));
    }
    let kind = match prefix_choice(&ty, &["pearson", "kendall", "spearman"]) {
        Some(0) => Kindc::Pearson,
        Some(1) => Kindc::Kendall,
        Some(2) => Kindc::Spearman,
        _ => return rt("The 'type' parameter value must be 'Pearson', 'Kendall', or 'Spearman'."),
    };
    let tl = tail.to_lowercase();
    let tail_kind = if tl.is_empty() {
        None
    } else {
        prefix_choice(&tl, &["left", "both", "right"]).or_else(|| prefix_choice(&tl, &["lt", "ne", "gt"]))
    };
    let tail_kind = match tail_kind {
        Some(0) => Tail::Left,
        Some(1) => Tail::Both,
        Some(2) => Tail::Right,
        _ => return rt("Tail value must be 'both', 'right', or 'left'."),
    };
    if x.is_complex() || y.is_complex() {
        if nargout > 1 {
            return rt("P-values are not available for complex data.");
        }
        if kind != Kindc::Pearson {
            return rt("Rank correlations are not available for complex data.");
        }
    }
    if matches!(x.kind, Kind::Int(_)) || matches!(y.kind, Kind::Int(_)) {
        return rt("MTIMES (*) is not fully supported for integer classes. At least one argument must be scalar.");
    }
    let col = |d: &[f64], c: usize| d[c * n..(c + 1) * n].to_vec();
    let xc: Vec<Vec<f64>> = (0..p1).map(|c| col(&xd, c)).collect();
    let yc: Vec<Vec<f64>> = if xx { xc.clone() } else { (0..p2).map(|c| col(&yd, c)).collect() };
    let any_nan = xc.iter().chain(yc.iter()).any(|c| c.iter().any(|v| v.is_nan()));
    let pairwise = rows_kind == Rows::Pairwise && any_nan;
    // Spearman without pairwise deletion: ties anywhere send every pair to the tie branch.
    let global_ties = !pairwise
        && kind == Kindc::Spearman
        && xc.iter().chain(if xx { [].iter() } else { yc.iter() }).any(|c| {
            let clean: Vec<f64> = c.iter().copied().filter(|v| !v.is_nan()).collect();
            tiedrank_vec(&clean, false, false).1[0] > 0.0
        });
    let mut cache = Cache { spearman: HashMap::new() };
    let mut coef = vec![f64::NAN; p1 * p2];
    let mut pval = vec![f64::NAN; p1 * p2];
    for i in 0..p1 {
        let jmax = if xx { i + 1 } else { p2 };
        for j in 0..jmax {
            let (xs, ys): (Vec<f64>, Vec<f64>) = if pairwise {
                xc[i].iter().zip(&yc[j]).filter(|(a, b)| !a.is_nan() && !b.is_nan()).map(|(a, b)| (*a, *b)).unzip()
            } else {
                (xc[i].clone(), yc[j].clone())
            };
            let nij = xs.len();
            let rows_removed = nij < n;
            let has_nan = xs.iter().chain(ys.iter()).any(|v| v.is_nan());
            let (c, p) = if nij < 2 || has_nan {
                (f64::NAN, f64::NAN)
            } else {
                match kind {
                    Kindc::Pearson => {
                        let r = pearson(&xs, &ys);
                        (r, pearson_p(tail_kind, r, nij as f64))
                    }
                    Kindc::Spearman => {
                        let (xr, xa) = tiedrank_vec(&xs, false, false);
                        let (yr, ya) = tiedrank_vec(&ys, false, false);
                        let nf = nij as f64;
                        let n3const = (nf + 1.0) * nf * (nf - 1.0) / 3.0;
                        let r = if pairwise {
                            let d: f64 = xr.iter().zip(&yr).map(|(a, b)| (a - b).powi(2)).sum();
                            let mean_d = (n3const - (xa[0] + ya[0]) / 3.0) / 2.0;
                            let n3c2 = (nf + 1.0) * nf * (nf - 1.0) / 2.0;
                            let std_d = if xa[0] == n3c2 || ya[0] == n3c2 { 0.0 } else { ((n3const / 2.0 - xa[0] / 3.0) * (n3const / 2.0 - ya[0] / 3.0) / (nf - 1.0)).sqrt() };
                            let r = (mean_d - d) / ((nf - 1.0).sqrt() * std_d);
                            if r.abs() > 1.0 { r / r.abs() } else { r }
                        } else {
                            pearson(&xr, &yr)
                        };
                        let ties = if pairwise { rows_removed || xa[0] > 0.0 || ya[0] > 0.0 } else { global_ties };
                        let p = if nargout > 1 {
                            if ties {
                                let d: f64 = xr.iter().zip(&yr).map(|(a, b)| (a - b).powi(2)).sum();
                                let mean_d = (n3const - (xa[0] + ya[0]) / 3.0) / 2.0;
                                let n3c2 = (nf + 1.0) * nf * (nf - 1.0) / 2.0;
                                let std_d = if xa[0] == n3c2 || ya[0] == n3c2 { 0.0 } else { ((n3const / 2.0 - xa[0] / 3.0) * (n3const / 2.0 - ya[0] / 3.0) / (nf - 1.0)).sqrt() };
                                spearman_p_ties(tail_kind, d, mean_d, std_d, &xr, &yr)
                            } else {
                                let mean_d = n3const / 2.0;
                                let std_d = n3const / (2.0 * (nf - 1.0).sqrt());
                                let d = if pairwise {
                                    xr.iter().zip(&yr).map(|(a, b)| (a - b).powi(2)).sum::<f64>()
                                } else {
                                    mean_d - (r * (nf - 1.0).sqrt() * std_d).round()
                                };
                                spearman_p_noties(tail_kind, d, nij, &mut cache)
                            }
                        } else {
                            f64::NAN
                        };
                        (r, p)
                    }
                    Kindc::Kendall => {
                        let (xr, xa) = tiedrank_vec(&xs, true, false);
                        let (yr, ya) = tiedrank_vec(&ys, true, false);
                        let nf = nij as f64;
                        let n2 = nf * (nf - 1.0) / 2.0;
                        let mut k = 0.0;
                        for a in 0..nij - 1 {
                            for b in a + 1..nij {
                                k += (xr[a] - xr[b]).signum_zero() * (yr[a] - yr[b]).signum_zero();
                            }
                        }
                        let r = k / ((n2 - xa[0]) * (n2 - ya[0])).sqrt();
                        let r = if r.abs() > 1.0 { r / r.abs() } else { r };
                        let ties = xa[0] > 0.0 || ya[0] > 0.0;
                        let p = if nargout > 1 {
                            if ties {
                                let std_k = if xa[0] == n2 || ya[0] == n2 {
                                    0.0
                                } else {
                                    (n2 * (2.0 * nf + 5.0) / 9.0 + xa[0] * ya[0] / n2 + xa[1] * ya[1] / (18.0 * n2 * (nf - 2.0)) - (xa[2] + ya[2]) / 18.0).sqrt()
                                };
                                kendall_p(tail_kind, k, std_k, nij, Some((&xr, &yr)))
                            } else {
                                kendall_p(tail_kind, k, (n2 * (2.0 * nf + 5.0) / 9.0).sqrt(), nij, None)
                            }
                        } else {
                            f64::NAN
                        };
                        (r, p)
                    }
                }
            };
            coef[j * p1 + i] = c;
            pval[j * p1 + i] = p;
            if xx && i != j {
                coef[i * p1 + j] = c;
                pval[i * p1 + j] = p;
            }
        }
        if xx {
            let d = coef[i * p1 + i];
            let keep = kind == Kindc::Kendall && tiedrank_vec(&xc[i], true, false).1[0] > 0.0;
            if !keep {
                coef[i * p1 + i] = if d.is_nan() { f64::NAN } else { d.signum() };
            }
            pval[i * p1 + i] = if d.is_nan() { f64::NAN } else { d.signum() };
        }
    }
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(vec![p1, p2], coef), single))];
    if nargout > 1 {
        outs.push(Value::Num(finish(Matrix::from_dims(vec![p1, p2], pval), single)));
    }
    Ok(Some(outs))
}
