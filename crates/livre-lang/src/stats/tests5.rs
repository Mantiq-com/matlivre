//! `adtest`: the Anderson-Darling test for a composite hypothesis (normal, exponential or
//! extreme value family).

use super::tests3::{param_val, pchip_at};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["adtest"];

const AD_ALPHAS: [f64; 27] = [
    0.0005, 0.001, 0.0015, 0.002, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99,
];
// Critical values by level: the value for an infinite sample, then the corrections in 1/n and 1/n^2
// (normal), 1/n (exponential) and 1/sqrt(n) (extreme value). Recovered from adtest's own output.
const AD_NORM_A: [f64; 27] = [
    1.5649, 1.4407, 1.3699, 1.3187, 1.1556, 1.0339, 0.8733, 0.7519, 0.6308, 0.5598, 0.5092, 0.4694, 0.4366, 0.4084, 0.3835, 0.3611, 0.3405, 0.3212, 0.3029, 0.2852, 0.2679, 0.2506, 0.2330, 0.2144, 0.1935, 0.1673, 0.1296,
];
const AD_NORM_B: [f64; 27] = [
    -0.9362, -0.9029, -0.8906, -0.8865, -0.8375, -0.7835, -0.6746, -0.5835, -0.4775, -0.4094, -0.3679, -0.3327, -0.3099, -0.2969, -0.2795, -0.2623, -0.2464, -0.2325, -0.2164, -0.1994, -0.1784, -0.1569, -0.1377, -0.1201, -0.0989, -0.0800, -0.0598,
];
const AD_NORM_C: [f64; 27] = [
    -8.3249, -6.6022, -5.6461, -4.9685, -3.2208, -2.1647, -1.2460, -0.7803, -0.4627, -0.3672, -0.2833, -0.2349, -0.1442, -0.0229, 0.0377, 0.0817, 0.1150, 0.1583, 0.1801, 0.1887, 0.1695, 0.1513, 0.1533, 0.1724, 0.2027, 0.3158, 0.6431,
];
const AD_EXP_A: [f64; 27] = [
    3.2371, 2.9303, 2.7541, 2.6307, 2.2454, 1.9621, 1.5928, 1.3223, 1.0621, 0.9153, 0.8134, 0.7355, 0.6725, 0.6194, 0.5734, 0.5326, 0.4957, 0.4617, 0.4301, 0.4001, 0.3712, 0.3428, 0.3144, 0.2849, 0.2527, 0.2131, 0.1581,
];
const AD_EXP_B: [f64; 27] = [
    1.6146, 0.8716, 0.4715, 0.2066, -0.4682, -0.7691, -0.7388, -0.5758, -0.4036, -0.3142, -0.2564, -0.2152, -0.1845, -0.1607, -0.1409, -0.1239, -0.1084, -0.0942, -0.0807, -0.0674, -0.0537, -0.0401, -0.0261, -0.0116, 0.0047, 0.0275, 0.0780,
];
const AD_EV_A: [f64; 27] = [
    1.6473, 1.5095, 1.4301, 1.3742, 1.1974, 1.0667, 0.8961, 0.7683, 0.6416, 0.5680, 0.5156, 0.4744, 0.4405, 0.4115, 0.3858, 0.3626, 0.3415, 0.3217, 0.3029, 0.2848, 0.2672, 0.2496, 0.2315, 0.2124, 0.1909, 0.1633, 0.1223,
];
const AD_EV_B: [f64; 27] = [
    -0.7097, -0.5934, -0.5328, -0.4930, -0.3708, -0.2973, -0.2075, -0.1449, -0.0892, -0.0619, -0.0442, -0.0302, -0.0196, -0.0112, -0.0039, 0.0024, 0.0074, 0.0122, 0.0167, 0.0207, 0.0245, 0.0282, 0.0323, 0.0371, 0.0436, 0.0549, 0.0813,
];

fn adtest_vec(v: &Value) -> Option<Vec<f64>> {
    match v {
        Value::Num(m) if m.im.is_none() && m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1) => Some(m.data.clone()),
        _ => None,
    }
}

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn call_cdf(it: &mut Interpreter, f: &str, x: &[f64], params: &[f64]) -> Res<Vec<f64>> {
    let mut args = vec![Value::Num(Matrix::from_dims(vec![x.len(), 1], x.to_vec()))];
    args.extend(params.iter().map(|&p| Value::num(p)));
    match it.call_function(f, args, 1)?.first() {
        Some(Value::Num(m)) => Ok(m.data.clone()),
        _ => rt(format!("{f} returned no value")),
    }
}

fn ad_stat(z: &[f64]) -> f64 {
    let n = z.len();
    let mut s = z.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut acc = 0.0;
    for k in 0..n {
        let w = 2.0 * (k as f64 + 1.0) - 1.0;
        acc += w * (s[k].ln() + (1.0 - s[n - 1 - k]).ln());
    }
    -acc / n as f64 - n as f64
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "adtest" => adtest(it, a, nargout),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn adtest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    let rest = &a[1..];
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let names = ["distribution", "alpha", "mctol", "asymptotic"];
    let mut distr_v = Value::chars("normal");
    let mut alpha_v = Value::num(0.05);
    let mut mctol_v = Value::empty();
    let mut asym_v = Value::bool(false);
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match pick {
            Some(0) => distr_v = pair[1].clone(),
            Some(1) => alpha_v = pair[1].clone(),
            Some(2) => mctol_v = pair[1].clone(),
            Some(3) => asym_v = pair[1].clone(),
            _ if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let asymptotic = match &asym_v {
        Value::Num(m) if m.data.len() == 1 => m.data[0] != 0.0,
        _ => false,
    };
    if asymptotic && !is_empty_val(&mctol_v) {
        // not reached by a text distribution below; the order of the checks follows MATLAB
    }
    if asymptotic {
        return rt("The asymptotic Anderson-Darling test is not valid when testing for a distribution family.");
    }
    let Some(x0) = adtest_vec(&a[0]) else { return rt("Sample data X must be a real vector.") };
    let single_x = matches!(&a[0], Value::Num(m) if is_single(m));
    let single_alpha = matches!(&alpha_v, Value::Num(m) if is_single(m));
    let x: Vec<f64> = x0.into_iter().filter(|v| !v.is_nan()).collect();
    let n = x.len();
    let alpha = match &alpha_v {
        Value::Num(m) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] > 0.0 && m.data[0] < 1.0 => m.data[0],
        _ => return rt("Significance level ALPHA must be a scalar between 0 and 1."),
    };
    if !is_empty_val(&mctol_v) {
        match &mctol_v {
            Value::Num(m) if m.data.len() == 1 && m.data[0] > 0.0 => {}
            _ => return rt("Monte-Carlo standard error tolerance MCTOL must be a positive scalar value."),
        }
        return rt("The 'MCTol' option (Monte Carlo simulation) is not implemented.");
    }
    if !distr_v.is_text() {
        return rt("'Distribution' value must be a character vector or string scalar representing the distribution family for the composite hypothesis test, or a probability distribution object for the simple hypothesis test.");
    }
    let dnames = ["normal", "exponential", "ev", "lognormal", "weibull", "extreme value"];
    let mut distr = dnames[param_val(&distr_v, &dnames, "distribution name", false)?];
    if n < 4 {
        return rt("Not enough data for the composite test: you must enter at least 4 sample points.");
    }
    let mut xs = x.clone();
    let mut real_ok = true;
    if distr == "lognormal" || distr == "weibull" {
        if xs.iter().any(|&v| v < 0.0) {
            real_ok = false;
        }
        xs = xs.iter().map(|v| v.ln()).collect();
        distr = if distr == "lognormal" { "normal" } else { "ev" };
    }
    let nf = n as f64;
    let adstat = match distr {
        "normal" => {
            if !real_ok {
                it.emit_warning("Data is not compatible with this distribution test. Setting p-value to zero.".into(), "stats:adtest:BadDataForDistribution".into());
                f64::NAN
            } else {
                let mean = xs.iter().sum::<f64>() / nf;
                let sd = (xs.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (nf - 1.0)).sqrt();
                ad_stat(&call_cdf(it, "normcdf", &xs, &[mean, sd])?)
            }
        }
        "exponential" => {
            let mean = xs.iter().sum::<f64>() / nf;
            ad_stat(&call_cdf(it, "expcdf", &xs, &[mean])?)
        }
        _ => {
            if !real_ok {
                it.emit_warning("Data is not compatible with this distribution test. Setting p-value to zero.".into(), "stats:adtest:BadDataForDistribution".into());
                f64::NAN
            } else {
                let params = it.call_function("evfit", vec![Value::Num(Matrix::from_dims(vec![n, 1], xs.clone()))], 1)?;
                let Some(Value::Num(pm)) = params.first() else { return rt("evfit returned no value") };
                ad_stat(&call_cdf(it, "evcdf", &xs, &[pm.data[0], pm.data[1]])?)
            }
        }
    };
    let cvs: Vec<f64> = (0..27)
        .map(|i| match distr {
            "normal" => AD_NORM_A[i] + AD_NORM_B[i] / nf + AD_NORM_C[i] / (nf * nf),
            "exponential" => AD_EXP_A[i] + AD_EXP_B[i] / nf,
            _ => AD_EV_A[i] + AD_EV_B[i] / nf.sqrt(),
        })
        .collect();
    let la: Vec<f64> = AD_ALPHAS.iter().map(|v| v.ln()).collect();
    let mut cv = pchip_at(&la, &cvs, alpha.ln());
    if alpha < AD_ALPHAS[0] {
        cv = cvs[0];
        it.emit_warning("Significance level ALPHA is outside range of lookup table. You may get more accurate results using the MCTOL input argument. The critical value output may be inaccurate.".into(), "stats:adtest:BadAlpha3".into());
    } else if alpha > AD_ALPHAS[26] {
        cv = cvs[26];
        it.emit_warning("Significance level ALPHA is outside range of lookup table. You may get more accurate results using the MCTOL input argument. The critical value output may be inaccurate.".into(), "stats:adtest:BadAlpha3".into());
    }
    let pval;
    if adstat > cvs[0] {
        it.emit_warning("P is less than the smallest tabulated value, returning 0.0005.".into(), "stats:adtest:OutOfRangePLow".into());
        pval = AD_ALPHAS[0];
    } else if adstat < cvs[26] {
        it.emit_warning("P is greater than the largest tabulated value, returning 0.99.".into(), "stats:adtest:OutOfRangePHigh".into());
        pval = AD_ALPHAS[26];
    } else if adstat.is_nan() {
        pval = 0.0;
    } else {
        let i = cvs.iter().position(|&c| adstat > c).unwrap_or(1).max(1);
        let (mut lo, mut hi) = (la[i - 1], la[i]);
        // the pchip value decreases with the level: root of pchip(x) - adstat
        let f = |t: f64| pchip_at(&la, &cvs, t) - adstat;
        let (mut flo, _fhi) = (f(lo), f(hi));
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            let fm = f(mid);
            if fm == 0.0 {
                lo = mid;
                hi = mid;
                break;
            }
            if (fm > 0.0) == (flo > 0.0) {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
            if (hi - lo).abs() < 1e-16 * lo.abs().max(1e-300) {
                break;
            }
        }
        pval = (0.5 * (lo + hi)).exp();
    }
    let h = if adstat.is_nan() {
        true
    } else if alpha < AD_ALPHAS[0] || alpha > AD_ALPHAS[26] {
        pval < alpha
    } else {
        adstat > cv
    };
    let sx = |v: f64| Value::Num(finish(Matrix::scalar(v), single_x));
    let mut out = vec![Value::bool(h), sx(pval)];
    if nargout > 2 {
        out.push(sx(adstat));
    }
    if nargout > 3 {
        out.push(Value::Num(finish(Matrix::scalar(cv), single_alpha)));
    }
    Ok(Some(out))
}
