//! Kolmogorov-Smirnov tests: `kstest` (one sample) and `kstest2` (two samples).

use super::dist::normcdf_std;
use super::util::*;
use crate::builtins::Ret;
use crate::builtins_specfun::gammaln;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["kstest", "kstest2"];

/// Exact one-sided critical values of the KS statistic (n = 1..20) at the one-sided
/// significance levels 0.005, 0.01, 0.025, 0.05, 0.10.
const EXACT: [[f64; 5]; 20] = [
    [0.99500, 0.99000, 0.97500, 0.95000, 0.90000],
    [0.92929, 0.90000, 0.84189, 0.77639, 0.68377],
    [0.82900, 0.78456, 0.70760, 0.63604, 0.56481],
    [0.73424, 0.68887, 0.62394, 0.56522, 0.49265],
    [0.66853, 0.62718, 0.56328, 0.50945, 0.44698],
    [0.61661, 0.57741, 0.51926, 0.46799, 0.41037],
    [0.57581, 0.53844, 0.48342, 0.43607, 0.38148],
    [0.54179, 0.50654, 0.45427, 0.40962, 0.35831],
    [0.51332, 0.47960, 0.43001, 0.38746, 0.33910],
    [0.48893, 0.45662, 0.40925, 0.36866, 0.32260],
    [0.46770, 0.43670, 0.39122, 0.35242, 0.30829],
    [0.44905, 0.41918, 0.37543, 0.33815, 0.29577],
    [0.43247, 0.40362, 0.36143, 0.32549, 0.28470],
    [0.41762, 0.38970, 0.34890, 0.31417, 0.27481],
    [0.40420, 0.37713, 0.33760, 0.30397, 0.26588],
    [0.39201, 0.36571, 0.32733, 0.29472, 0.25778],
    [0.38086, 0.35528, 0.31796, 0.28627, 0.25039],
    [0.37062, 0.34569, 0.30936, 0.27851, 0.24360],
    [0.36117, 0.33685, 0.30143, 0.27136, 0.23735],
    [0.35241, 0.32866, 0.29408, 0.26473, 0.23156],
];

fn parse_alpha_tail(alpha: Option<&Value>, tail: Option<&Value>, who: &str) -> Res<(f64, i32)> {
    let a = match alpha {
        Some(Value::Num(m)) if !m.is_empty() => {
            if m.numel() != 1 || m.data[0] <= 0.0 || m.data[0] >= 1.0 {
                return rt("Significance level ALPHA must be a scalar between 0 and 1.");
            }
            m.data[0]
        }
        _ => 0.05,
    };
    let bad = || -> Res<(f64, i32)> { rt("Type-of-test indicator TAIL must be 'unequal', 'smaller', or 'larger'.") };
    let _ = who;
    let t = match tail {
        Some(Value::Num(m)) if m.is_empty() => 0,
        Some(v) if v.is_text() => {
            let s = v.text()?.to_lowercase();
            match ["smaller", "unequal", "larger"].iter().position(|c| !s.is_empty() && c.starts_with(&s)) {
                Some(i) => i as i32 - 1,
                None => return bad(),
            }
        }
        Some(Value::Num(m)) if m.numel() == 1 && [-1.0, 0.0, 1.0].contains(&m.data[0]) => m.data[0] as i32,
        Some(_) => return bad(),
        None => 0,
    };
    Ok((a, t))
}

fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..n {
            let v = a[i][k];
            if v != 0.0 {
                for j in 0..n {
                    c[i][j] += v * b[k][j];
                }
            }
        }
    }
    c
}

fn mat_pow(t: &[Vec<f64>], mut e: usize) -> Vec<Vec<f64>> {
    let n = t.len();
    let mut result: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    let mut base = t.to_vec();
    while e > 0 {
        if e & 1 == 1 {
            result = mat_mul(&result, &base);
        }
        base = mat_mul(&base, &base);
        e >>= 1;
    }
    result
}

fn kstest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(x0)) = a.first() else { return rt("Input 'X' must be a vector.") };
    let rest = &a[1..];
    let (mut cdf, mut alpha, mut tail): (Option<&Value>, Option<&Value>, Option<&Value>) = (None, None, None);
    if let Some(first) = rest.first() {
        if matches!(first, Value::Num(m) if m.kind != Kind::Char) {
            cdf = Some(first);
            alpha = rest.get(1);
            tail = rest.get(2);
        } else {
            if rest.len() % 2 == 1 {
                return rt("Invalid parameter name or missing value.");
            }
            for pair in rest.chunks(2) {
                let key = pair[0].text().unwrap_or_default().to_lowercase();
                match ["cdf", "alpha", "tail"].iter().position(|n| !key.is_empty() && n.starts_with(&key)) {
                    Some(0) => cdf = Some(&pair[1]),
                    Some(1) => alpha = Some(&pair[1]),
                    Some(2) => tail = Some(&pair[1]),
                    _ => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
                }
            }
        }
    }
    let (alpha, tail) = parse_alpha_tail(alpha, tail, "kstest")?;
    if x0.kind == Kind::Char {
        return rt("Error calculating CDF of input 'X'.");
    }
    let single = is_single(x0);
    let mut xs: Vec<f64> = x0.data.iter().copied().filter(|v| !v.is_nan()).collect();
    let n = xs.len();
    if n == 0 {
        return rt("Input 'X' must be a vector.");
    }
    xs.sort_by(|p, q| p.partial_cmp(q).unwrap());
    // ecdf: distinct values with cumulative frequencies
    let mut ux: Vec<f64> = Vec::new();
    let mut sample_cdf: Vec<f64> = vec![0.0];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && xs[j + 1] == xs[i] {
            j += 1;
        }
        ux.push(xs[i]);
        sample_cdf.push((j + 1) as f64 / n as f64);
        i = j + 1;
    }
    let (xcdf, ycdf): (Vec<f64>, Vec<f64>) = match cdf {
        Some(Value::Num(c)) if !c.is_empty() => {
            if c.dims[1] != 2 {
                return rt("Hypothesized CDF matrix must have 2 columns.");
            }
            let r = c.dims[0];
            let mut rows: Vec<(f64, f64)> = (0..r).map(|i| (c.data[i], c.data[r + i])).filter(|(p, q)| !(p + q).is_nan()).collect();
            if rows.is_empty() {
                return rt("Hypothesized CDF matrix must have at least 1 valid row.");
            }
            rows.sort_by(|p, q| p.0.partial_cmp(&q.0).unwrap());
            if rows.windows(2).any(|w| w[1].1 < w[0].1) {
                return rt("CDF must define an increasing function of X.");
            }
            let mut dedup: Vec<(f64, f64)> = Vec::new();
            for w in rows.iter() {
                if let Some(last) = dedup.last() {
                    if last.0 == w.0 {
                        if last.1 != w.1 {
                            return rt("CDF must not have duplicate X values.");
                        }
                        continue;
                    }
                }
                dedup.push(*w);
            }
            (dedup.iter().map(|p| p.0).collect(), dedup.iter().map(|p| p.1).collect())
        }
        _ => (ux.clone(), ux.iter().map(|&v| normcdf_std(v)).collect()),
    };
    let null_cdf: Vec<f64> = if xcdf == ux {
        ycdf.clone()
    } else {
        if ux[0] < xcdf[0] || *ux.last().unwrap() > *xcdf.last().unwrap() {
            return rt("Hypothesized CDF matrix must span the observations interval in X.");
        }
        ux.iter()
            .map(|&v| {
                let k = xcdf.partition_point(|&c| c <= v).saturating_sub(1).min(xcdf.len().saturating_sub(2));
                if xcdf.len() == 1 {
                    return ycdf[0];
                }
                let t = (v - xcdf[k]) / (xcdf[k + 1] - xcdf[k]);
                ycdf[k] + t * (ycdf[k + 1] - ycdf[k])
            })
            .collect()
    };
    let m = ux.len();
    let mut deltas: Vec<f64> = Vec::new();
    for i in 0..m {
        let d1 = sample_cdf[i] - null_cdf[i];
        let d2 = sample_cdf[i + 1] - null_cdf[i];
        match tail {
            0 => {
                deltas.push(d1.abs());
                deltas.push(d2.abs());
            }
            -1 => {
                deltas.push(-d1);
                deltas.push(-d2);
            }
            _ => {
                deltas.push(d1);
                deltas.push(d2);
            }
        }
    }
    let ks = deltas.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let nf = n as f64;
    let pval;
    if tail == 0 {
        let s = nf * ks * ks;
        if s > 7.24 || (s > 3.76 && n > 99) {
            pval = 2.0 * (-(2.000071 + 0.331 / nf.sqrt() + 1.409 / nf) * s).exp();
        } else {
            let k = (ks * nf).ceil() as usize;
            let h = k as f64 - ks * nf;
            let mm = 2 * k - 1;
            let fact = |i: usize| (1..=i).fold(1.0, |a, b| a * b as f64);
            let (lmax, tk);
            if mm > 1 {
                let mut t = vec![vec![0.0; mm]; mm];
                for i in 0..mm {
                    for j in 0..mm {
                        t[i][j] = if i >= j { 1.0 / fact(i - j + 1) } else if j == i + 1 { 1.0 } else { 0.0 };
                    }
                }
                for i in 0..mm {
                    t[i][0] -= h.powi(i as i32 + 1) / fact(i + 1);
                }
                let col: Vec<f64> = (0..mm).map(|i| t[i][0]).collect();
                for j in 0..mm {
                    t[mm - 1][j] = col[mm - 1 - j];
                }
                t[mm - 1][0] = (1.0 - 2.0 * h.powi(mm as i32) + (2.0 * h - 1.0).max(0.0).powi(mm as i32)) / fact(mm);
                let flat: Vec<f64> = (0..mm).flat_map(|j| (0..mm).map(|i| t[i][j]).collect::<Vec<_>>()).collect();
                let ev = it.call_function("eig", vec![Value::Num(Matrix::from_dims(vec![mm, mm], flat))], 1)?;
                let Some(Value::Num(e)) = ev.first() else { return rt("eig failed") };
                lmax = e.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let scaled: Vec<Vec<f64>> = t.iter().map(|r| r.iter().map(|v| v / lmax).collect()).collect();
                let tp = mat_pow(&scaled, n);
                tk = tp[k - 1][k - 1];
            } else {
                lmax = 1.0;
                tk = (1.0 - 2.0 * h + (2.0 * h - 1.0).max(0.0)) / fact(1);
            }
            pval = 1.0 - (gammaln(nf + 1.0) + nf * f64::ln(lmax) - nf * nf.ln()).exp() * tk;
        }
    } else {
        let t = nf * ks;
        let kstart = t.ceil() as usize;
        pval = (kstart..=n)
            .map(|k| {
                let kf = k as f64;
                (t.ln() - nf * nf.ln() + gammaln(nf + 1.0) - gammaln(kf + 1.0) - gammaln(nf - kf + 1.0) + kf * (kf - t).ln() + (nf - kf - 1.0) * (t + nf - kf).ln()).exp()
            })
            .sum();
    }
    let h = if pval < alpha { 1.0 } else { 0.0 };
    let mut outs = vec![Value::Num(Matrix::logical(1, 1, vec![h])), Value::Num(finish(Matrix::scalar(pval), single)), Value::Num(finish(Matrix::scalar(ks), single))];
    if nargout > 3 {
        let alpha1 = if tail == 0 { alpha / 2.0 } else { alpha };
        let cv = if (0.005..=0.10).contains(&alpha1) {
            if n <= 20 {
                let a1 = Matrix::from_dims(vec![1, 5], vec![0.005, 0.01, 0.025, 0.05, 0.10]);
                let y = Matrix::from_dims(vec![1, 5], EXACT[n - 1].to_vec());
                let r = it.call_function("spline", vec![Value::Num(a1), Value::Num(y), Value::num(alpha1)], 1)?;
                match r.first() {
                    Some(Value::Num(m)) => m.data[0],
                    _ => f64::NAN,
                }
            } else {
                let l10 = alpha1.log10();
                let aa = 0.09037 * (-l10).powf(1.5) + 0.01515 * l10 * l10 - 0.08467 * alpha1 - 0.11143;
                let asym = (-0.5 * alpha1.ln() / nf).sqrt();
                let c = asym - 0.16693 / nf - aa / nf.powf(1.5);
                c.min(1.0 - alpha1)
            }
        } else {
            f64::NAN
        };
        outs.push(Value::Num(finish(Matrix::scalar(cv), single)));
    }
    Ok(Some(outs))
}

fn kstest2(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("At least 2 inputs are required.");
    }
    let rest = &a[2..];
    let (mut alpha, mut tail): (Option<&Value>, Option<&Value>) = (None, None);
    if let Some(first) = rest.first() {
        if matches!(first, Value::Num(m) if m.kind != Kind::Char) {
            alpha = Some(first);
            tail = rest.get(1);
        } else {
            if rest.len() % 2 == 1 {
                return rt("Invalid parameter name or missing value.");
            }
            for pair in rest.chunks(2) {
                let key = pair[0].text().unwrap_or_default().to_lowercase();
                match ["alpha", "tail"].iter().position(|n| !key.is_empty() && n.starts_with(&key)) {
                    Some(0) => alpha = Some(&pair[1]),
                    Some(1) => tail = Some(&pair[1]),
                    _ => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
                }
            }
        }
    }
    let (Value::Num(m1), Value::Num(m2)) = (&a[0], &a[1]) else { return rt("The samples X1 and X2 must be vectors.") };
    let isv = |m: &Matrix| m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1) && !m.is_empty();
    if !isv(m1) || !isv(m2) {
        return rt("The samples X1 and X2 must be vectors.");
    }
    let single = is_single(m1) || is_single(m2);
    let x1: Vec<f64> = m1.data.iter().copied().filter(|v| !v.is_nan()).collect();
    let x2: Vec<f64> = m2.data.iter().copied().filter(|v| !v.is_nan()).collect();
    if x1.is_empty() {
        return rt("Sample vector X1 contains no data.");
    }
    if x2.is_empty() {
        return rt("Sample vector X2 contains no data.");
    }
    let (alpha, tail) = parse_alpha_tail(alpha, tail, "kstest2")?;
    let mut all: Vec<f64> = x1.iter().chain(x2.iter()).copied().collect();
    all.sort_by(|p, q| p.partial_cmp(q).unwrap());
    let cdf_at = |s: &[f64], v: f64| s.iter().filter(|&&t| t <= v).count() as f64 / s.len() as f64;
    let mut best = f64::NEG_INFINITY;
    for &v in &all {
        let (f1, f2) = (cdf_at(&x1, v), cdf_at(&x2, v));
        let d = match tail {
            0 => (f1 - f2).abs(),
            -1 => f2 - f1,
            _ => f1 - f2,
        };
        if d > best {
            best = d;
        }
    }
    // the bin below the smallest observation contributes zero
    let ks = best.max(0.0);
    let (n1, n2) = (x1.len() as f64, x2.len() as f64);
    let n = n1 * n2 / (n1 + n2);
    let lambda = ((n.sqrt() + 0.12 + 0.11 / n.sqrt()) * ks).max(0.0);
    let p = if tail != 0 {
        (-2.0 * lambda * lambda).exp()
    } else {
        let s: f64 = (1..=101).map(|j| (if (j - 1) % 2 == 0 { 1.0 } else { -1.0 }) * (-2.0 * lambda * lambda * (j * j) as f64).exp()).sum();
        (2.0 * s).clamp(0.0, 1.0)
    };
    let h = if alpha >= p { 1.0 } else { 0.0 };
    let outs = vec![Value::Num(Matrix::logical(1, 1, vec![h])), Value::Num(Matrix::scalar(p)), Value::Num(Matrix::scalar(ks))];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if name == "kstest" { kstest(it, a, nargout) } else { kstest2(a, nargout) }
}
