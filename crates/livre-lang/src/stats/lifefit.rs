//! Likelihoods and fits of the gamma, Weibull, beta, extreme value, generalized extreme value and
//! generalized Pareto distributions (`gamlike`, `wbllike`, `betalike`, `evlike`, `gevlike`, `gplike`,
//! `evfit`, `wblfit`, `gamfit`, `betafit`), plus `statset` and `statget`.

use super::util::*;
use crate::builtins::Ret;
use crate::builtins_specfun::gammaln;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["gamlike", "wbllike", "betalike", "evlike", "gevlike", "gplike", "evfit", "wblfit", "gamfit", "betafit", "statset", "statget"];

fn matrix_arg<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(_) => rt("Input must be numeric."),
        None => rt("Not enough input arguments."),
    }
}

fn row_of(v: Vec<f64>, single: bool) -> Value {
    let n = v.len();
    Value::Num(finish(Matrix::from_dims(vec![1, n], v), single))
}

fn mat2(d: [f64; 4], single: bool) -> Value {
    Value::Num(finish(Matrix::from_dims(vec![2, 2], d.to_vec()), single))
}

fn inv2(h11: f64, h12: f64, h22: f64) -> [f64; 4] {
    let det = h11 * h22 - h12 * h12;
    [h22 / det, -h12 / det, -h12 / det, h11 / det]
}

fn lit(v: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&format!("{:e}", x));
    }
    s.push(']');
    s
}

fn anon(it: &mut Interpreter, src: &str) -> Res<Value> {
    let r = it.call_function("str2func", vec![Value::chars(src)], 1)?;
    Ok(r.into_iter().next().unwrap())
}

fn optimset_tolx(it: &mut Interpreter, tolx: f64, tolfun: Option<f64>) -> Res<Value> {
    let mut args = vec![Value::chars("Display"), Value::chars("off"), Value::chars("TolX"), Value::num(tolx)];
    if let Some(t) = tolfun {
        args.push(Value::chars("TolFun"));
        args.push(Value::num(t));
    }
    Ok(it.call_function("optimset", args, 1)?.into_iter().next().unwrap())
}

fn psi(it: &mut Interpreter, k: f64, x: f64) -> Res<f64> {
    it.call_function("psi", vec![Value::num(k), Value::num(x)], 1)?[0].scalar()
}

/// `(data, censoring, freq)` after the checks shared by the `*like` functions (zero weights removed).
struct Inputs {
    data: Vec<f64>,
    cens: Vec<bool>,
    freq: Vec<f64>,
    single: bool,
    has_cens: bool,
}

fn like_inputs(a: &[Value], nparams: usize, who: &str, min_args_msg: &str, allow_cens: bool) -> Res<(Vec<f64>, Inputs)> {
    if a.len() < 2 {
        return rt(min_args_msg);
    }
    let data = matrix_arg(a, 1)?;
    if data.numel() > data.dims.iter().copied().max().unwrap_or(0) {
        return rt("DATA must be a vector.");
    }
    let params = matrix_arg(a, 0)?;
    if params.numel() != nparams {
        return rt(format!("First input must be a vector with {nparams} elements."));
    }
    let _ = who;
    let mut cens: Vec<bool> = vec![false; data.numel()];
    let mut has_cens = false;
    if allow_cens {
        if let Some(Value::Num(c)) = a.get(2) {
            if !c.is_empty() {
                if c.dims != data.dims {
                    return rt("DATA and CENSORING must have the same size.");
                }
                cens = c.data.iter().map(|&v| v != 0.0).collect();
                has_cens = true;
            }
        }
    }
    let mut freq = vec![1.0; data.numel()];
    if allow_cens {
        if let Some(Value::Num(f)) = a.get(3) {
            if !f.is_empty() {
                if f.dims != data.dims {
                    return rt("DATA and FREQ must have the same size.");
                }
                freq = f.data.clone();
            }
        }
    }
    let single = is_single(data) || is_single(params);
    let mut d = Vec::new();
    let mut c = Vec::new();
    let mut f = Vec::new();
    for i in 0..data.numel() {
        if freq[i] != 0.0 {
            d.push(data.data[i]);
            c.push(cens[i]);
            f.push(freq[i]);
        }
    }
    Ok((params.data.clone(), Inputs { data: d, cens: c, freq: f, single, has_cens }))
}

fn pair(a: f64, b: Value, nargout_hint: usize) -> Ret {
    let _ = nargout_hint;
    Ok(Some(vec![Value::num(a), b]))
}

fn out2(nl: f64, avar: [f64; 4], single: bool) -> Ret {
    Ok(Some(vec![Value::Num(finish(Matrix::scalar(nl), single)), mat2(avar, single)]))
}

fn gamlike(it: &mut Interpreter, a: &[Value]) -> Ret {
    let (p, inp) = like_inputs(a, 2, "gamlike", "Requires two input arguments", true)?;
    let (sh, sc) = (if p[0] <= 0.0 { f64::NAN } else { p[0] }, if p[1] <= 0.0 { f64::NAN } else { p[1] });
    let n = inp.data.len();
    if n == 0 {
        return out2(-0.0, [f64::NEG_INFINITY, -0.0, -0.0, f64::NEG_INFINITY], inp.single);
    }
    let x: Vec<f64> = inp.data.iter().map(|&v| if v <= 0.0 { f64::NAN } else { v }).collect();
    let z: Vec<f64> = x.iter().map(|v| v / sc).collect();
    let mut l: Vec<f64> = z.iter().map(|&zz| (sh - 1.0) * zz.ln() - zz - gammaln(sh) - sc.ln()).collect();
    let ncen: f64 = inp.cens.iter().zip(&inp.freq).map(|(&c, &f)| if c { f } else { 0.0 }).sum();
    if ncen > 0.0 {
        for i in 0..n {
            if inp.cens[i] {
                let s = it.call_function("gammainc", vec![Value::num(z[i]), Value::num(sh), Value::chars("upper")], 1)?[0].scalar()?;
                l[i] = s.ln();
            }
        }
    }
    for i in 0..n {
        if z[i] == f64::INFINITY {
            l[i] = f64::NEG_INFINITY;
        }
    }
    let nl = -l.iter().zip(&inp.freq).map(|(l, f)| f * l).sum::<f64>();
    if ncen > 0.0 {
        return out2(nl, [f64::NAN; 4], inp.single);
    }
    let tri = psi(it, 1.0, sh)?;
    let nh11 = inp.freq.iter().map(|f| f * tri).sum::<f64>();
    let nh12 = inp.freq.iter().map(|f| f / sc).sum::<f64>();
    let nh22 = z.iter().zip(&inp.freq).map(|(zz, f)| f * (2.0 * zz - sh) / (sc * sc)).sum::<f64>();
    let av = if [nh11, nh12, nh22].iter().any(|v| v.is_nan()) { [f64::NAN; 4] } else { inv2(nh11, nh12, nh22) };
    out2(nl, av, inp.single)
}

fn wbllike(a: &[Value]) -> Ret {
    let (p, inp) = like_inputs(a, 2, "wbllike", "Not enough input arguments.", true)?;
    let (aa, bb) = (if p[0] <= 0.0 { f64::NAN } else { p[0] }, if p[1] <= 0.0 { f64::NAN } else { p[1] });
    let n = inp.data.len();
    if n == 0 {
        return out2(-0.0, [f64::NAN; 4], inp.single);
    }
    let (mut nl, mut h11, mut h12, mut h22) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        let x = if inp.data[i] < 0.0 { f64::NAN } else { inp.data[i] };
        let z = x / aa;
        let lz = z.ln();
        let z2b = (bb * lz).exp();
        let unc = if inp.cens[i] { 0.0 } else { 1.0 };
        let mut l = ((bb - 1.0) * lz + (bb / aa).ln()) * unc - z2b;
        if z == f64::INFINITY {
            l = f64::NEG_INFINITY;
        }
        let f = inp.freq[i];
        nl -= f * l;
        h11 += f * (bb * ((1.0 + bb) * z2b - unc));
        h12 += f * ((1.0 + bb * lz) * z2b - unc);
        h22 += f * (lz * lz * z2b + unc / (bb * bb));
    }
    let nh11 = h11 / (aa * aa);
    let nh12 = -h12 / aa;
    let nh22 = h22;
    let det = nh11 * nh22 - nh12 * nh12;
    let av = [nh22 / det, -nh12 / det, -nh12 / det, nh11 / det];
    out2(nl, av, inp.single)
}

fn betaln(a: f64, b: f64) -> f64 {
    gammaln(a) + gammaln(b) - gammaln(a + b)
}

fn betalike(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Requires two input arguments.");
    }
    let data = matrix_arg(a, 1)?;
    if data.numel() > data.dims.iter().copied().max().unwrap_or(0) || (data.dims.len() > 2) {
        return rt("DATA must be a vector.");
    }
    let params = matrix_arg(a, 0)?;
    if params.numel() != 2 {
        return rt("First input must be a vector with 2 elements.");
    }
    let single = is_single(data) || is_single(params);
    let (pa, pb) = (if params.data[0] <= 0.0 { f64::NAN } else { params.data[0] }, if params.data[1] <= 0.0 { f64::NAN } else { params.data[1] });
    let mut x: Vec<f64> = data.data.clone();
    let (xmin, xmax) = (x.iter().cloned().fold(f64::INFINITY, f64::min), x.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    for v in x.iter_mut() {
        if !(0.0 <= *v && *v <= 1.0) {
            *v = f64::NAN;
        }
    }
    let (xl, xu) = (f64::MIN_POSITIVE.sqrt(), 1.0 - f64::EPSILON / 2.0);
    let (mut n0, mut n1) = (0usize, 0usize);
    let x2: Vec<f64> = if xl <= xmin && xmax <= xu {
        x.clone()
    } else {
        n0 = x.iter().filter(|&&v| v < xl).count();
        n1 = x.iter().filter(|&&v| v > xu).count();
        x.iter().copied().filter(|&v| !(v < xl) && !(v > xu)).collect()
    };
    let n2 = x2.len() as f64;
    let l1: Vec<f64> = x2.iter().map(|v| v.ln()).collect();
    let l2: Vec<f64> = x2.iter().map(|v| (-v).ln_1p()).collect();
    let (s1, s2): (f64, f64) = (l1.iter().sum(), l2.iter().sum());
    let mut nl = n2 * betaln(pa, pb) - (pa - 1.0) * s1 - (pb - 1.0) * s2;
    let binc = |it: &mut Interpreter, xx: f64, aa: f64, bb: f64, tail: &str| -> Res<f64> { it.call_function("betainc", vec![Value::num(xx), Value::num(aa), Value::num(bb), Value::chars(tail)], 1)?[0].scalar() };
    if n0 > 0 {
        nl -= n0 as f64 * binc(it, xl, pa, pb, "lower")?.ln();
    }
    if n1 > 0 {
        nl -= n1 as f64 * binc(it, xu, pa, pb, "upper")?.ln();
    }
    let nlv = Value::Num(finish(Matrix::scalar(nl), single));
    if nargout < 2 {
        return Ok(Some(vec![nlv]));
    }
    if data.numel() < 2 {
        return rt("To compute AVAR, DATA must have at least two elements.");
    }
    let (dga, dgb, dgab) = (psi(it, 0.0, pa)?, psi(it, 0.0, pb)?, psi(it, 0.0, pa + pb)?);
    let mut j: Vec<[f64; 2]> = (0..x2.len()).map(|i| [l1[i] + dgab - dga, l2[i] + dgab - dgb]).collect();
    if n0 > 0 || n1 > 0 {
        let delta = f64::EPSILON.sqrt();
        let (aa, bb) = ([pa + pa * delta, pa - pa * delta], [pb + pb * delta, pb - pb * delta]);
        for (cnt, xx, tail) in [(n0, xl, "lower"), (n1, xu, "upper")] {
            if cnt > 0 {
                let da = (binc(it, xx, aa[1], pb, tail)?.ln() - binc(it, xx, aa[0], pb, tail)?.ln()) / (2.0 * pa * delta);
                let db = (binc(it, xx, pa, bb[1], tail)?.ln() - binc(it, xx, pa, bb[0], tail)?.ln()) / (2.0 * pb * delta);
                for _ in 0..cnt {
                    j.push([da, db]);
                }
            }
        }
    }
    let (mut s11, mut s12, mut s22) = (0.0, 0.0, 0.0);
    for r in &j {
        s11 += r[0] * r[0];
        s12 += r[0] * r[1];
        s22 += r[1] * r[1];
    }
    let av = if [s11, s12, s22].iter().any(|v| v.is_nan()) { [f64::NAN; 4] } else { inv2(s11, s12, s22) };
    Ok(Some(vec![nlv, mat2(av, single)]))
}

fn evlike_core(params: &[f64], inp: &Inputs) -> (f64, [f64; 4]) {
    let (mu, sg) = (params[0], if params[1] <= 0.0 { f64::NAN } else { params[1] });
    let (mut nl, mut h11, mut h12, mut h22) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..inp.data.len() {
        let z = (inp.data[i] - mu) / sg;
        let ez = z.exp();
        let unc = if inp.cens[i] { 0.0 } else { 1.0 };
        let mut l = (z - sg.ln()) * unc - ez;
        if z == f64::INFINITY {
            l = f64::NEG_INFINITY;
        }
        let f = inp.freq[i];
        nl -= f * l;
        h11 += f * ez;
        h12 += f * ((z + 1.0) * ez - unc);
        h22 += f * (z * (z + 2.0) * ez - (2.0 * z + 1.0) * unc);
    }
    let det = h11 * h22 - h12 * h12;
    let s2 = sg * sg;
    (nl, [s2 * h22 / det, -s2 * h12 / det, -s2 * h12 / det, s2 * h11 / det])
}

fn evlike(a: &[Value]) -> Ret {
    let (p, inp) = like_inputs(a, 2, "evlike", "Requires two input arguments", true)?;
    let (nl, av) = if inp.data.is_empty() { (-0.0, [f64::NAN; 4]) } else { evlike_core(&p, &inp) };
    out2(nl, av, inp.single)
}

fn gevlike(it: &mut Interpreter, a: &[Value]) -> Ret {
    let (p, inp) = like_inputs(a, 3, "gevlike", "Requires two input arguments", false)?;
    let (k, sigma, mu) = (p[0], p[1], p[2]);
    let n = inp.data.len() as f64;
    let lns = sigma.ln();
    let z: Vec<f64> = inp.data.iter().map(|d| (d - mu) / sigma).collect();
    let sm = |f: &dyn Fn(f64) -> f64| -> f64 { z.iter().map(|&v| f(v)).sum() };
    let (nl, nh): (f64, [[f64; 3]; 3]);
    if k.abs() > f64::EPSILON {
        let u: Vec<f64> = z.iter().map(|v| 1.0 + k * v).collect();
        if u.iter().cloned().fold(f64::INFINITY, f64::min) > 0.0 {
            let lnu: Vec<f64> = z.iter().map(|v| (k * v).ln_1p()).collect();
            let t: Vec<f64> = lnu.iter().map(|l| (-(1.0 / k) * l).exp()).collect();
            nl = n * lns + t.iter().sum::<f64>() + (1.0 + 1.0 / k) * lnu.iter().sum::<f64>();
            let (mut d2k2, mut d2ks, mut d2km, mut d2s2, mut d2sm, mut d2m2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for i in 0..z.len() {
                let (zi, ui, li, ti) = (z[i], u[i], lnu[i], t[i]);
                let s = (-(1.0 / k) * li).exp_m1();
                let r = (s - k) / ui;
                let q = (ti - k) / (ui * ui);
                d2k2 += 2.0 * (zi / ui - li / k) * s + (k + 1.0) * zi * zi * q - 2.0 * li * zi * ti / (ui * k) + li * li * ti / (k * k);
                d2ks += zi * (zi - 1.0 - (k + 1.0) * zi * ti / k + (k * zi + 1.0) * li * ti / (k * k)) / (ui * ui);
                d2km += (zi - 1.0 - (k + 1.0) * zi * ti / k + (k * zi + 1.0) * li * ti / (k * k)) / (ui * ui);
                d2s2 += (k + 1.0) * zi * zi * q - 2.0 * zi * r - 1.0;
                d2sm += zi * ti / (ui * ui) - r / ui;
                d2m2 += (k + 1.0) * q;
            }
            let (d2k2, d2ks, d2km) = (d2k2 / (k * k), d2ks / sigma, d2km / sigma);
            let (d2s2, d2sm, d2m2) = (d2s2 / (sigma * sigma), d2sm / (sigma * sigma), d2m2 / (sigma * sigma));
            nh = [[d2k2, d2ks, d2km], [d2ks, d2s2, d2sm], [d2km, d2sm, d2m2]];
        } else {
            nl = f64::NAN;
            nh = [[f64::NAN; 3]; 3];
        }
    } else {
        let t: Vec<f64> = z.iter().map(|v| (-v).exp()).collect();
        nl = n * lns + (0..z.len()).map(|i| t[i] + z[i]).sum::<f64>();
        let d2k2 = sm(&|v| (-v).exp() * v.powi(4) / 4.0 - 2.0 * (-v).exp_m1() * v.powi(3) - v * v);
        let d2ks = sm(&|v| ((v * v * (-v).exp() / 2.0) - v * (-v).exp_m1() - 1.0) * v) / sigma;
        let d2km = sm(&|v| v * v * (-v).exp() / 2.0 - v * (-v).exp_m1()) / sigma;
        let d2s2 = sm(&|v| v * (v + 2.0) * (-v).exp() - (2.0 * v + 1.0)) / (sigma * sigma);
        let d2sm = sm(&|v| (v + 1.0) * (-v).exp() - 1.0) / (sigma * sigma);
        let d2m2 = sm(&|v| (-v).exp()) / (sigma * sigma);
        nh = [[d2k2, d2ks, d2km], [d2ks, d2s2, d2sm], [d2km, d2sm, d2m2]];
    }
    let flat: Vec<f64> = (0..9).map(|idx| nh[idx % 3][idx / 3]).collect();
    let inv = it.call_function("inv", vec![Value::Num(Matrix::from_dims(vec![3, 3], flat))], 1)?;
    let Value::Num(im) = &inv[0] else { return rt("inv failed") };
    Ok(Some(vec![Value::Num(finish(Matrix::scalar(nl), inp.single)), Value::Num(finish(im.clone(), inp.single))]))
}

fn gplike(a: &[Value]) -> Ret {
    let (p, inp) = like_inputs(a, 2, "gplike", "Requires two input arguments", false)?;
    let (k, sigma) = (p[0], p[1]);
    let n = inp.data.len() as f64;
    let lns = sigma.ln();
    let z: Vec<f64> = inp.data.iter().map(|d| d / sigma).collect();
    let single = inp.single;
    if k.abs() > f64::EPSILON {
        let zmax = z.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if k > 0.0 || zmax < -1.0 / k {
            let u: Vec<f64> = z.iter().map(|v| 1.0 + k * v).collect();
            let sumlnu: f64 = z.iter().map(|v| (k * v).ln_1p()).sum();
            let nl = n * lns + (1.0 + 1.0 / k) * sumlnu;
            let v: Vec<f64> = (0..z.len()).map(|i| z[i] / u[i]).collect();
            let (sumv, sumvsq): (f64, f64) = (v.iter().sum(), v.iter().map(|x| x * x).sum());
            let nh11 = 2.0 * sumlnu / k.powi(3) - 2.0 * sumv / (k * k) - (1.0 + 1.0 / k) * sumvsq;
            let nh12 = (-sumv + (k + 1.0) * sumvsq) / sigma;
            let nh22 = (-n + 2.0 * (k + 1.0) * sumv - k * (k + 1.0) * sumvsq) / (sigma * sigma);
            let det = nh11 * nh22 - nh12 * nh12;
            return out2(nl, [nh22 / det, -nh12 / det, -nh12 / det, nh11 / det], single);
        }
        return out2(f64::INFINITY, [f64::NAN; 4], single);
    }
    let (sumz, sumzsq, sumzcb): (f64, f64, f64) = (z.iter().sum(), z.iter().map(|v| v * v).sum(), z.iter().map(|v| v * v * v).sum());
    let nl = n * lns + sumz;
    let nh11 = (2.0 / 3.0) * sumzcb - sumzsq;
    let nh12 = (-sumz + sumzsq) / sigma;
    let nh22 = (-n + 2.0 * sumz) / (sigma * sigma);
    let det = nh11 * nh22 - nh12 * nh12;
    out2(nl, [nh22 / det, -nh12 / det, -nh12 / det, nh11 / det], single)
}

fn norminv(p: f64, mu: f64, s: f64) -> f64 {
    mu + s * super::dist::norminv_std(p)
}

struct FitIn {
    x: Vec<f64>,
    cens: Vec<bool>,
    freq: Vec<f64>,
    alpha: f64,
    single: bool,
    tolx: f64,
}

fn fit_inputs(a: &[Value], who: &str, default_tolx: f64, allow_options: bool) -> Res<FitIn> {
    let x = matrix_arg(a, 0)?;
    let single = is_single(x);
    let alpha = match a.get(1) {
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        _ => 0.05,
    };
    let mut cens = vec![false; x.numel()];
    if let Some(Value::Num(c)) = a.get(2) {
        if !c.is_empty() {
            if c.dims != x.dims {
                return rt(if who == "gamfit" { "X and CENSORING must have the same size.".to_string() } else { "X and CENSORING must have the same size.".to_string() });
            }
            cens = c.data.iter().map(|&v| v != 0.0).collect();
        }
    }
    let mut freq = vec![1.0; x.numel()];
    if let Some(Value::Num(f)) = a.get(3) {
        if !f.is_empty() {
            if f.dims != x.dims {
                return rt("X and FREQ must have the same size.");
            }
            freq = f.data.clone();
        }
    }
    let mut tolx = default_tolx;
    if allow_options {
        if let Some(Value::Struct(s)) = a.get(4) {
            if let Some(Value::Num(t)) = s.get("TolX") {
                if !t.is_empty() {
                    tolx = t.data[0];
                }
            }
        }
    }
    let mut xs = Vec::new();
    let mut cs = Vec::new();
    let mut fs = Vec::new();
    for i in 0..x.numel() {
        if freq[i] != 0.0 {
            xs.push(x.data[i]);
            cs.push(cens[i]);
            fs.push(freq[i]);
        }
    }
    Ok(FitIn { x: xs, cens: cs, freq: fs, alpha, single, tolx })
}

fn is_vector(m: &Matrix) -> bool {
    m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1)
}

/// `[parmhat, parmci]` of the extreme value fit; `None` for the interval when it is not wanted.
fn evfit_core(it: &mut Interpreter, f: &FitIn, want_ci: bool) -> Res<([f64; 2], [f64; 4])> {
    let nan = ([f64::NAN; 2], [f64::NAN; 4]);
    let n: f64 = f.freq.iter().sum();
    let ncens: f64 = f.cens.iter().zip(&f.freq).map(|(&c, &w)| if c { w } else { 0.0 }).sum();
    let nunc = n - ncens;
    let maxx = f.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let minx = f.x.iter().cloned().fold(f64::INFINITY, f64::min);
    let rangex = maxx - minx;
    if n == 0.0 || nunc == 0.0 || !rangex.is_finite() {
        return Ok(nan);
    }
    let tiny = f64::MIN_POSITIVE;
    let mut range_unc = 0.0;
    if ncens == 0.0 {
        if rangex < tiny {
            let ph = [f.x[0], 0.0];
            let ci = if n == 1.0 { [f64::NEG_INFINITY, f64::INFINITY, 0.0, f64::INFINITY] } else { [ph[0], ph[0], 0.0, 0.0] };
            return Ok((ph, ci));
        }
    } else {
        let xu: Vec<f64> = (0..f.x.len()).filter(|&i| !f.cens[i]).map(|i| f.x[i]).collect();
        range_unc = xu.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - xu.iter().cloned().fold(f64::INFINITY, f64::min);
        if range_unc < tiny && xu[0] >= maxx {
            let ph = [xu[0], 0.0];
            let ci = if nunc == 1.0 { [f64::NEG_INFINITY, f64::INFINITY, 0.0, f64::INFINITY] } else { [ph[0], ph[0], 0.0, 0.0] };
            return Ok((ph, ci));
        }
    }
    let x0: Vec<f64> = f.x.iter().map(|v| (v - maxx) / rangex).collect();
    let (sigma0, wgt);
    let col = |v: &[f64]| Value::Num(Matrix::from_dims(vec![v.len(), 1], v.to_vec()));
    if ncens == 0.0 {
        let m = x0.iter().sum::<f64>() / x0.len() as f64;
        let sd = (x0.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (x0.len() as f64 - 1.0)).sqrt();
        sigma0 = 6.0f64.sqrt() * sd / std::f64::consts::PI;
        wgt = x0.iter().zip(&f.freq).map(|(v, w)| v * w).sum::<f64>() / n;
    } else {
        if range_unc > 0.0 {
            let cv: Vec<f64> = f.cens.iter().map(|&c| if c { 1.0 } else { 0.0 }).collect();
            let r = it.call_function("ecdf", vec![col(&x0), Value::chars("censoring"), col(&cv), Value::chars("frequency"), col(&f.freq)], 2)?;
            let (Value::Num(p), Value::Num(q)) = (&r[0], &r[1]) else { return rt("ecdf failed") };
            let m = p.numel();
            let xs: Vec<f64> = (0..m - 1).map(|i| { let pm = (p.data[i] + p.data[i + 1]) / 2.0; (-(1.0 - pm).ln()).ln() }).collect();
            let ys: Vec<f64> = q.data[1..].to_vec();
            let pf = it.call_function("polyfit", vec![col(&xs), col(&ys), Value::num(1.0)], 1)?;
            let Value::Num(pm) = &pf[0] else { return rt("polyfit failed") };
            sigma0 = pm.data[0];
        } else {
            sigma0 = 1.0;
        }
        wgt = (0..x0.len()).filter(|&i| !f.cens[i]).map(|i| x0[i] * f.freq[i]).sum::<f64>() / nunc;
    }
    let lk = |s: f64| -> f64 {
        let w: Vec<f64> = x0.iter().zip(&f.freq).map(|(x, w)| w * (x / s).exp()).collect();
        s + wgt - x0.iter().zip(&w).map(|(x, w)| x * w).sum::<f64>() / w.iter().sum::<f64>()
    };
    let (mut lower, mut upper);
    if lk(sigma0) > 0.0 {
        upper = sigma0;
        lower = 0.5 * upper;
        while lk(lower) > 0.0 {
            upper = lower;
            lower = 0.5 * upper;
            if lower < tiny {
                return rt("Unable to find a solution: the fit has no positive root.");
            }
        }
    } else {
        lower = sigma0;
        upper = 2.0 * lower;
        while lk(upper) < 0.0 {
            lower = upper;
            upper = 2.0 * lower;
            if upper > f64::MAX {
                return rt("Unable to find a solution: the fit has no finite root.");
            }
        }
    }
    let src = format!("@(s) s + {:e} - sum({} .* ({} .* exp({} ./ s))) / sum({} .* exp({} ./ s))", wgt, lit(&x0), lit(&f.freq), lit(&x0), lit(&f.freq), lit(&x0));
    let fun = anon(it, &src)?;
    let opts = optimset_tolx(it, f.tolx, None)?;
    let r = it.call_function("fzero", vec![fun, Value::Num(Matrix::from_dims(vec![1, 2], vec![lower, upper])), opts], 1)?;
    let sigmahat = r[0].scalar()?;
    let muhat = sigmahat * ((0..x0.len()).map(|i| f.freq[i] * (x0[i] / sigmahat).exp()).sum::<f64>() / nunc).ln();
    let ph = [rangex * muhat + maxx, rangex * sigmahat];
    if !want_ci {
        return Ok((ph, [f64::NAN; 4]));
    }
    let inp = Inputs { data: f.x.clone(), cens: f.cens.clone(), freq: f.freq.clone(), single: false, has_cens: ncens > 0.0 };
    let (_, acov) = evlike_core(&ph, &inp);
    let (se1, se2) = (acov[0].sqrt(), acov[3].sqrt() / ph[1]);
    let (lo, hi) = (f.alpha / 2.0, 1.0 - f.alpha / 2.0);
    let ci = [norminv(lo, ph[0], se1), norminv(hi, ph[0], se1), norminv(lo, ph[1].ln(), se2).exp(), norminv(hi, ph[1].ln(), se2).exp()];
    Ok((ph, ci))
}

fn evfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = matrix_arg(a, 0)?;
    if !is_vector(x) {
        return rt("X must be a vector.");
    }
    let f = fit_inputs(a, "evfit", 1e-6, true)?;
    if f.x.iter().any(|v| v.is_nan()) && false {
        return rt("X must not contain NaN.");
    }
    if f.x.iter().any(|v| v.is_nan()) {
        return rt("Second argument must be finite.");
    }
    let (ph, ci) = evfit_core(it, &f, nargout > 1)?;
    Ok(Some(vec![row_of(ph.to_vec(), f.single), mat2(ci, f.single)]))
}

fn wblfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = matrix_arg(a, 0)?;
    if !is_vector(x) || x.data.iter().any(|&v| v <= 0.0) {
        return rt("X must be a vector containing positive values.");
    }
    let mut f = fit_inputs(a, "wblfit", 1e-6, true)?;
    f.x = f.x.iter().map(|v| v.ln()).collect();
    let (ph, ci) = evfit_core(it, &f, nargout > 1)?;
    let p = [ph[0].exp(), 1.0 / ph[1]];
    // ci rows: [lo; hi] of (mu, sigma): mu interval [ci0, ci1], sigma interval [ci2, ci3]
    let c = [ci[0].exp(), ci[1].exp(), 1.0 / ci[3], 1.0 / ci[2]];
    Ok(Some(vec![row_of(p.to_vec(), f.single), mat2(c, f.single)]))
}

fn gamfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = matrix_arg(a, 0)?;
    if a.len() > 2 {
        if let Value::Num(c) = &a[2] {
            if !c.is_empty() && c.dims != x.dims {
                return rt("X and CENSORING must have the same size.");
            }
        }
    }
    if let Some(Value::Num(fq)) = a.get(3) {
        if !fq.is_empty() && fq.dims != x.dims {
            return rt("X and FREQ must have the same size.");
        }
    }
    if x.dims.len() > 2 || (x.dims[0] > 1 && x.dims[1] > 1) {
        return rt("X must be a vector.");
    }
    let f = fit_inputs(a, "gamfit", 1e-8, true)?;
    let n: f64 = f.freq.iter().sum();
    let ncen: f64 = f.cens.iter().zip(&f.freq).map(|(&c, &w)| if c { w } else { 0.0 }).sum();
    let nunc = n - ncen;
    if ncen == 0.0 && f.x.iter().any(|&v| v < 0.0) {
        return rt("X must be non-negative.");
    }
    if ncen > 0.0 && f.x.iter().any(|&v| v <= 0.0) {
        return rt("X must be positive when there is censoring.");
    }
    let single = f.single;
    let nanr = |s| Ok(Some(vec![row_of(vec![f64::NAN, f64::NAN], s), mat2([f64::NAN; 4], s)]));
    if n == 0.0 || nunc == 0.0 || f.x.iter().any(|v| !v.is_finite()) {
        return nanr(single);
    }
    if ncen > 0.0 {
        return rt("Censored gamma fits are not supported.");
    }
    let scalex = f.x.iter().zip(&f.freq).map(|(x, w)| x * w).sum::<f64>() / n;
    if scalex < f64::MIN_POSITIVE {
        return Ok(Some(vec![row_of(vec![f64::NAN, 0.0], single), mat2([f64::NAN, f64::NAN, 0.0, 0.0], single)]));
    }
    let x: Vec<f64> = f.x.iter().map(|v| v / scalex).collect();
    let s2 = x.iter().zip(&f.freq).map(|(x, w)| w * (x - 1.0).powi(2)).sum::<f64>() / n;
    if s2 <= 100.0 * f64::EPSILON {
        let ci = if nunc > 1.0 { [f64::INFINITY, f64::INFINITY, 0.0, 0.0] } else { [0.0, f64::INFINITY, 0.0, f64::INFINITY] };
        return Ok(Some(vec![row_of(vec![f64::INFINITY, 0.0], single), mat2(ci, single)]));
    }
    let s2 = s2 * n / (n - 1.0);
    let ahat0 = 1.0 / s2;
    let bhat0 = s2;
    if x.iter().any(|&v| v == 0.0) {
        return Ok(Some(vec![row_of(vec![ahat0, bhat0 * scalex], single), mat2([f64::NAN; 4], single)]));
    }
    let sumlogx: f64 = x.iter().zip(&f.freq).map(|(x, w)| w * x.ln()).sum();
    let c = sumlogx / n;
    let mut lkeq = |it: &mut Interpreter, av: f64| -> Res<f64> { Ok(-c - av.ln() + psi(it, 0.0, av)?) };
    let (mut lower, mut upper);
    if lkeq(it, ahat0)? > 0.0 {
        upper = ahat0;
        lower = 0.5 * upper;
        while lkeq(it, lower)? > 0.0 {
            upper = lower;
            lower = 0.5 * upper;
            if lower < f64::MIN_POSITIVE {
                return rt("Unable to find a solution: the fit has no positive root.");
            }
        }
    } else {
        lower = ahat0;
        upper = 2.0 * lower;
        while lkeq(it, upper)? < 0.0 {
            lower = upper;
            upper = 2.0 * lower;
            if upper > f64::MAX {
                return rt("Unable to find a solution: the fit has no finite root.");
            }
        }
    }
    let fun = anon(it, &format!("@(a) -({:e}) - log(a) + psi(a)", c))?;
    let opts = optimset_tolx(it, f.tolx, None)?;
    let r = it.call_function("fzero", vec![fun, Value::Num(Matrix::from_dims(vec![1, 2], vec![lower, upper])), opts], 1)?;
    let ahat = r[0].scalar()?;
    let ph = [ahat, 1.0 / ahat];
    let mut ci = [f64::NAN; 4];
    if nargout > 1 {
        let tri = psi(it, 1.0, ph[0])?;
        let nh11 = f.freq.iter().map(|w| w * tri).sum::<f64>();
        let nh12 = f.freq.iter().map(|w| w / ph[1]).sum::<f64>();
        let nh22 = x.iter().zip(&f.freq).map(|(z0, w)| { let z = z0 / ph[1]; w * (2.0 * z - ph[0]) / (ph[1] * ph[1]) }).sum::<f64>();
        let acov = inv2(nh11, nh12, nh22);
        let (s1, s2v) = (acov[0].sqrt() / ph[0], acov[3].sqrt() / ph[1]);
        let (lo, hi) = (f.alpha / 2.0, 1.0 - f.alpha / 2.0);
        ci = [norminv(lo, ph[0].ln(), s1).exp(), norminv(hi, ph[0].ln(), s1).exp(), norminv(lo, ph[1].ln(), s2v).exp() * scalex, norminv(hi, ph[1].ln(), s2v).exp() * scalex];
    }
    Ok(Some(vec![row_of(vec![ph[0], ph[1] * scalex], single), mat2(ci, single)]))
}

fn betafit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = matrix_arg(a, 0)?;
    let single = is_single(x);
    if x.is_empty() {
        return Ok(Some(vec![row_of(vec![f64::NAN; 2], single), mat2([f64::NAN; 4], single)]));
    }
    if !is_vector(x) {
        return rt("The first input must be a vector.");
    }
    let alpha = match a.get(1) {
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        _ => 0.05,
    };
    let xs: Vec<f64> = x.data.iter().copied().filter(|v| !v.is_nan()).collect();
    if xs.is_empty() {
        return Ok(Some(vec![row_of(vec![f64::NAN; 2], single), mat2([f64::NAN; 4], single)]));
    }
    let xmin = xs.iter().cloned().fold(f64::INFINITY, f64::min);
    let xmax = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if xmin < 0.0 || xmax > 1.0 {
        return rt("All values must be within the closed interval [0,1].");
    }
    if (xmin - xmax).abs() <= 2.0 * crate::stats::util::eps_of(xmax) {
        return rt("Cannot fit a beta distribution if all data values are the same.");
    }
    if single {
        return rt("FMINSEARCH accepts inputs only of data type double.");
    }
    let n = xs.len() as f64;
    let slx: f64 = xs.iter().map(|v| v.ln()).sum();
    let sl1: f64 = xs.iter().map(|v| (-v).ln_1p()).sum();
    let (t1, t2) = ((sl1 / n).exp(), (slx / n).exp());
    let t3 = 1.0 - t1 - t2;
    let (ah, bh) = (0.5 * (1.0 - t1) / t3, 0.5 * (1.0 - t2) / t3);
    let (xl, xu) = (f64::MIN_POSITIVE.sqrt(), 1.0 - f64::EPSILON / 2.0);
    let src = if xl <= xmin && xmax <= xu {
        format!("@(p) {n:e}*betaln(exp(p(1)),exp(p(2))) - (exp(p(1))-1)*({slx:e}) - (exp(p(2))-1)*({sl1:e})")
    } else {
        let x2: Vec<f64> = xs.iter().copied().filter(|&v| !(v < xl) && !(v > xu)).collect();
        let (n0, n1) = (xs.iter().filter(|&&v| v < xl).count() as f64, xs.iter().filter(|&&v| v > xu).count() as f64);
        let (s1, s2): (f64, f64) = (x2.iter().map(|v| v.ln()).sum(), x2.iter().map(|v| (-v).ln_1p()).sum());
        let mut s = format!("@(p) {:e}*betaln(exp(p(1)),exp(p(2))) - (exp(p(1))-1)*({s1:e}) - (exp(p(2))-1)*({s2:e})", x2.len() as f64);
        if n0 > 0.0 {
            s.push_str(&format!(" - {n0:e}*log(betainc({xl:e},exp(p(1)),exp(p(2)),'lower'))"));
        }
        if n1 > 0.0 {
            s.push_str(&format!(" - {n1:e}*log(betainc({xu:e},exp(p(1)),exp(p(2)),'upper'))"));
        }
        s
    };
    let fun = anon(it, &src)?;
    let opts = optimset_tolx(it, 1e-6, Some(1e-6))?;
    let start = Value::Num(Matrix::from_dims(vec![1, 2], vec![ah.ln(), bh.ln()]));
    let r = it.call_function("fminsearch", vec![fun, start, opts], 1)?;
    let Value::Num(ph) = &r[0] else { return rt("fminsearch failed") };
    let phat = [ph.data[0].exp(), ph.data[1].exp()];
    let mut ci = [f64::NAN; 4];
    if nargout > 1 {
        let args = vec![Value::Num(Matrix::from_dims(vec![1, 2], phat.to_vec())), Value::Num(Matrix::from_dims(vec![xs.len(), 1], xs.clone()))];
        let (_, av) = {
            let r = betalike(it, &args, 2)?.unwrap();
            let Value::Num(av) = &r[1] else { return rt("betalike failed") };
            (0, av.data.clone())
        };
        let (s1, s2) = (av[0].sqrt() / phat[0], av[3].sqrt() / phat[1]);
        let (lo, hi) = (alpha / 2.0, 1.0 - alpha / 2.0);
        ci = [norminv(lo, phat[0].ln(), s1).exp(), norminv(hi, phat[0].ln(), s1).exp(), norminv(lo, phat[1].ln(), s2).exp(), norminv(hi, phat[1].ln(), s2).exp()];
    }
    Ok(Some(vec![row_of(phat.to_vec(), false), mat2(ci, false)]))
}

const STAT_FIELDS: [&str; 20] = ["Display", "MaxFunEvals", "MaxIter", "TolBnd", "TolFun", "TolTypeFun", "TolX", "TolTypeX", "GradObj", "Jacobian", "DerivStep", "FunValCheck", "Robust", "RobustWgtFun", "WgtFun", "Tune", "UseParallel", "UseSubstreams", "Streams", "OutputFcn"];

fn stat_defaults(name: &str) -> Option<Vec<(&'static str, Value)>> {
    let off = || Value::chars("off");
    let d = |v: f64| Value::num(v);
    let fam = |tolfun: f64, tolx: f64, mfe: f64, mi: f64| vec![("Display", off()), ("MaxFunEvals", d(mfe)), ("MaxIter", d(mi)), ("TolBnd", d(1e-6)), ("TolFun", d(tolfun)), ("TolX", d(tolx))];
    Some(match name {
        "evfit" | "wblfit" => vec![("Display", off()), ("TolX", d(1e-6))],
        "gamfit" | "normfit" => fam(1e-8, 1e-8, 200.0, 100.0),
        "gevfit" | "gpfit" | "nbinfit" => fam(1e-6, 1e-6, 400.0, 200.0),
        "kmeans" => vec![("Display", off()), ("MaxIter", d(100.0)), ("UseParallel", d(0.0)), ("UseSubstreams", d(0.0))],
        "nlinfit" | "fitnlm" => vec![("Display", off()), ("MaxIter", d(200.0)), ("TolFun", d(1e-8)), ("TolX", d(1e-8)), ("DerivStep", d(6.0554544523933429e-06)), ("FunValCheck", Value::chars("on")), ("Robust", off()), ("WgtFun", Value::chars("bisquare"))],
        "glmfit" => vec![("Display", off()), ("MaxIter", d(100.0)), ("TolX", d(1e-6))],
        "lasso" | "plsregress" => vec![("Display", off()), ("UseParallel", d(0.0)), ("UseSubstreams", d(0.0))],
        "copulafit" => vec![("Display", off()), ("MaxFunEvals", d(200.0)), ("MaxIter", d(100.0)), ("TolBnd", d(1e-6)), ("TolX", d(1e-6))],
        "bootci" | "bootstrp" | "crossval" => vec![("Display", off()), ("UseParallel", d(0.0)), ("UseSubstreams", d(0.0))],
        "jackknife" => vec![("Display", off()), ("UseParallel", d(0.0))],
        "mdscale" => vec![("Display", off()), ("MaxIter", d(200.0)), ("TolFun", d(1e-6)), ("TolX", d(1e-6))],
        "nnmf" => vec![("Display", off()), ("MaxIter", d(100.0)), ("TolFun", d(1e-4)), ("TolX", d(1e-4)), ("UseParallel", d(0.0)), ("UseSubstreams", d(0.0))],
        "tsne" => vec![("Display", off()), ("MaxIter", d(1000.0)), ("TolFun", d(1e-10))],
        "factoran" => vec![("Display", off()), ("MaxFunEvals", d(400.0)), ("MaxIter", d(100.0)), ("TolFun", d(1e-8)), ("TolX", d(1e-8))],
        "pca" | "ppca" => vec![("Display", off()), ("MaxIter", d(1000.0)), ("TolFun", d(1e-6)), ("TolX", d(1e-6))],
        "fitglme" => vec![("Display", off()), ("MaxIter", d(10000.0)), ("TolFun", d(1e-6)), ("TolX", d(1e-12))],
        "fitrgp" => vec![("Display", off()), ("MaxIter", d(10000.0)), ("TolFun", d(1e-6)), ("TolX", d(1e-12)), ("GradObj", Value::chars("on"))],
        "sequentialfs" => vec![("Display", off()), ("MaxIter", d(f64::INFINITY)), ("TolTypeFun", Value::chars("rel")), ("UseParallel", d(0.0)), ("UseSubstreams", d(0.0))],
        _ => return None,
    })
}

fn empty_stat_struct() -> StructArray {
    let mut s = StructArray::scalar();
    for f in STAT_FIELDS {
        s.set(f, if f == "Streams" { Value::Cell(CellArray::from_dims(vec![0, 0], vec![])) } else { Value::empty() });
    }
    s
}

fn field_match(name: &str) -> Res<Vec<&'static str>> {
    let n = name.to_lowercase();
    let exact: Vec<&'static str> = STAT_FIELDS.iter().copied().filter(|f| f.to_lowercase() == n).collect();
    if !exact.is_empty() {
        return Ok(exact);
    }
    Ok(STAT_FIELDS.iter().copied().filter(|f| !n.is_empty() && f.to_lowercase().starts_with(&n)).collect())
}

fn statset(a: &[Value]) -> Ret {
    let mut s = empty_stat_struct();
    if a.is_empty() {
        return Ok(Some(vec![Value::Struct(s)]));
    }
    let mut i = 0;
    if a.len() == 1 && a[0].is_text() {
        let name = a[0].text()?;
        match stat_defaults(&name) {
            Some(d) => {
                for (k, v) in d {
                    let v = if k == "UseParallel" || k == "UseSubstreams" { Value::bool(false) } else { v };
                    s.set(k, v);
                }
                return Ok(Some(vec![Value::Struct(s)]));
            }
            None => return rt(format!("No default options available for the function '{name}'.")),
        }
    }
    while i < a.len() {
        if let Value::Struct(old) = &a[i] {
            for (k, f) in old.fields.iter().enumerate() {
                if let Some(v) = old.elems.first().map(|e| e[k].clone()) {
                    if !v.as_matrix().map(|m| m.is_empty()).unwrap_or(false) {
                        s.set(f, v);
                    }
                }
            }
            i += 1;
        } else {
            break;
        }
    }
    if i == 0 && !a[0].is_text() {
        return rt("Expected argument 1 to be a character vector or string scalar specifying a parameter name, or an options structure created with STATSET.");
    }
    if (a.len() - i) % 2 != 0 {
        return rt("Arguments must occur in name-value pairs.");
    }
    while i < a.len() {
        if !a[i].is_text() {
            return rt("Expected argument 1 to be a character vector or string scalar specifying a parameter name, or an options structure created with STATSET.");
        }
        let nm = a[i].text()?;
        let m = field_match(&nm)?;
        if m.is_empty() {
            return rt(format!("'{nm}' is not a valid value for the parameter name argument."));
        }
        if m.len() > 1 {
            return rt(format!("The value '{nm}' is a partial match for more than one valid value for the parameter name argument."));
        }
        s.set(m[0], a[i + 1].clone());
        i += 2;
    }
    Ok(Some(vec![Value::Struct(s)]))
}

fn statget(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Requires at least two input arguments.");
    }
    let Value::Struct(s) = &a[0] else { return rt("First argument must be an options structure created with STATSET.") };
    let nm = a[1].text().unwrap_or_default();
    let m = field_match(&nm)?;
    if m.is_empty() {
        return rt(format!("'{nm}' is not a valid value for the parameter name argument."));
    }
    if m.len() > 1 {
        return rt(format!("The value '{nm}' is a partial match for more than one valid value for the parameter name argument."));
    }
    let v = s.get(m[0]).cloned().unwrap_or_else(Value::empty);
    let empty = matches!(&v, Value::Num(x) if x.is_empty());
    if empty {
        if let Some(d) = a.get(2) {
            return Ok(Some(vec![d.clone()]));
        }
    }
    Ok(Some(vec![v]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    let _ = pair;
    match name {
        "gamlike" => gamlike(it, a),
        "wbllike" => wbllike(a),
        "betalike" => betalike(it, a, nargout.max(1)),
        "evlike" => evlike(a),
        "gevlike" => gevlike(it, a),
        "gplike" => gplike(a),
        "evfit" => evfit(it, a, nargout),
        "wblfit" => wblfit(it, a, nargout),
        "gamfit" => gamfit(it, a, nargout),
        "betafit" => betafit(it, a, nargout),
        "statset" => statset(a),
        _ => statget(a),
    }
}
