//! Univariate distributions: `<name>pdf`, `<name>cdf` (with `'upper'`), `<name>inv` and
//! `<name>stat` for the families in `FAMILIES`. Each family supplies its density, distribution
//! function, quantile and moments; this module handles the argument rules (defaults, scalar
//! expansion, size checks, single precision, MATLAB's messages).

use super::errtab::ERRS;
use super::util::*;
use crate::nd::Dims;
use crate::builtins::Ret;
use crate::builtins_specfun::{betainc_pair, betaincinv, erfc, gammainc_pq, gammaincinv, gammaln};
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;
use std::f64::consts::{LN_2, PI, SQRT_2};

type P<'a> = &'a [f64];

pub struct Family {
    pub name: &'static str,
    /// Number of parameters.
    pub np: usize,
    /// Value a missing trailing parameter takes (None: required).
    pub defaults: &'static [Option<f64>],
    pub valid: fn(P) -> bool,
    pub pdf: fn(f64, P) -> f64,
    /// (x, params, upper tail)
    pub cdf: fn(f64, P, bool) -> f64,
    pub inv: fn(f64, P) -> f64,
    pub stat: fn(P) -> (f64, f64),
}

pub const FAMILY_NAMES: &[&str] = &[
    "norm", "logn", "exp", "unif", "gam", "beta", "chi2", "t", "f", "wbl", "rayl", "ev", "gev", "gp", "bino", "poiss", "geo", "nbin", "hyge", "unid",
];

pub const NAMES: &[&str] = &[
    "normpdf", "normcdf", "norminv", "normstat", "lognpdf", "logncdf", "logninv", "lognstat",
    "exppdf", "expcdf", "expinv", "expstat", "unifpdf", "unifcdf", "unifinv", "unifstat",
    "gampdf", "gamcdf", "gaminv", "gamstat", "betapdf", "betacdf", "betainv", "betastat",
    "chi2pdf", "chi2cdf", "chi2inv", "chi2stat", "tpdf", "tcdf", "tinv", "tstat",
    "fpdf", "fcdf", "finv", "fstat", "wblpdf", "wblcdf", "wblinv", "wblstat",
    "raylpdf", "raylcdf", "raylinv", "raylstat", "evpdf", "evcdf", "evinv", "evstat",
    "gevpdf", "gevcdf", "gevinv", "gevstat", "gppdf", "gpcdf", "gpinv", "gpstat",
    "binopdf", "binocdf", "binoinv", "binostat", "poisspdf", "poisscdf", "poissinv", "geopdf",
    "geocdf", "geoinv", "geostat", "nbinpdf", "nbincdf", "nbininv", "nbinstat", "hygepdf",
    "hygecdf", "hygeinv", "hygestat", "unidpdf", "unidcdf", "unidinv", "unidstat",
];

// ---- helpers ----------------------------------------------------------------------------------

fn norm_cdf(z: f64) -> f64 {
    0.5 * erfc(-z / SQRT_2)
}

fn norm_inv(p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    -SQRT_2 * crate::builtins_specfun::erfcinv(2.0 * p)
}

fn gamma_fn(x: f64) -> f64 {
    if x > 0.0 { gammaln(x).exp() } else if x == 0.0 { f64::INFINITY } else { PI / ((PI * x).sin() * gammaln(1.0 - x).exp()) }
}

fn isint(x: f64) -> bool {
    x.fract() == 0.0 && x.is_finite()
}

fn lchoose(n: f64, k: f64) -> f64 {
    gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)
}

fn up<T>(upper: bool, lower: f64, upr: f64) -> f64 {
    let _ = std::marker::PhantomData::<T>;
    if upper { upr } else { lower }
}

fn tail(p: f64, upper: bool) -> f64 {
    if upper { 1.0 - p } else { p }
}

// ---- families ---------------------------------------------------------------------------------

fn norm_pdf(x: f64, p: P) -> f64 {
    let z = (x - p[0]) / p[1];
    (-0.5 * z * z).exp() / ((2.0 * PI).sqrt() * p[1])
}
fn norm_cdf_f(x: f64, p: P, upper: bool) -> f64 {
    let z = (x - p[0]) / p[1];
    norm_cdf(if upper { -z } else { z })
}

fn logn_pdf(x: f64, p: P) -> f64 {
    if x > 0.0 {
        let z = (x.ln() - p[0]) / p[1];
        (-0.5 * z * z).exp() / (x * p[1] * (2.0 * PI).sqrt())
    } else if x.is_nan() {
        f64::NAN
    } else {
        0.0
    }
}
fn logn_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x > 0.0 {
        let z = (x.ln() - p[0]) / p[1];
        norm_cdf(if upper { -z } else { z })
    } else if x.is_nan() {
        f64::NAN
    } else if upper {
        1.0
    } else {
        0.0
    }
}

fn exp_pdf(x: f64, p: P) -> f64 {
    if x >= 0.0 { (-x / p[0]).exp() / p[0] } else if x.is_nan() { f64::NAN } else { 0.0 }
}
fn exp_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if upper { (-x / p[0]).exp() } else { -(-x / p[0]).exp_m1() }
}
fn exp_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) { f64::NAN } else { -p[0] * (-q).ln_1p() }
}

fn unif_pdf(x: f64, p: P) -> f64 {
    if x.is_nan() {
        f64::NAN
    } else if x >= p[0] && x <= p[1] {
        1.0 / (p[1] - p[0])
    } else {
        0.0
    }
}
fn unif_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let v = ((x - p[0]) / (p[1] - p[0])).clamp(0.0, 1.0);
    if upper { 1.0 - v } else { v }
}

fn gam_pdf(x: f64, p: P) -> f64 {
    let (a, b) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return if a < 1.0 { f64::INFINITY } else if a == 1.0 { 1.0 / b } else { 0.0 };
    }
    ((a - 1.0) * x.ln() - x / b - gammaln(a) - a * b.ln()).exp()
}
fn gam_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if x == f64::INFINITY {
        return if upper { 0.0 } else { 1.0 };
    }
    let (l, u) = gammainc_pq(x / p[1], p[0]);
    if upper { u } else { l }
}
fn gam_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) { f64::NAN } else { p[1] * gammaincinv(q, p[0], false) }
}

fn beta_pdf(x: f64, p: P) -> f64 {
    let (a, b) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    if !(0.0..=1.0).contains(&x) {
        return 0.0;
    }
    if x == 0.0 {
        return if a < 1.0 { f64::INFINITY } else if a == 1.0 { b } else { 0.0 };
    }
    if x == 1.0 {
        return if b < 1.0 { f64::INFINITY } else if b == 1.0 { a } else { 0.0 };
    }
    ((a - 1.0) * x.ln() + (b - 1.0) * (-x).ln_1p() - crate::builtins_specfun::betaln(a, b)).exp()
}
fn beta_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if x >= 1.0 {
        return if upper { 0.0 } else { 1.0 };
    }
    let (l, u) = betainc_pair(x, p[0], p[1]);
    if upper { u } else { l }
}
fn beta_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) { f64::NAN } else { betaincinv(q, p[0], p[1], false) }
}

fn chi2_pdf(x: f64, p: P) -> f64 {
    gam_pdf(x, &[p[0] / 2.0, 2.0])
}
fn chi2_cdf(x: f64, p: P, upper: bool) -> f64 {
    gam_cdf(x, &[p[0] / 2.0, 2.0], upper)
}
fn chi2_inv(q: f64, p: P) -> f64 {
    gam_inv(q, &[p[0] / 2.0, 2.0])
}

fn t_pdf(x: f64, p: P) -> f64 {
    let v = p[0];
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return 0.0;
    }
    (gammaln((v + 1.0) / 2.0) - gammaln(v / 2.0) - 0.5 * (v * PI).ln() - (v + 1.0) / 2.0 * (x * x / v).ln_1p()).exp()
}
fn t_cdf(x: f64, p: P, upper: bool) -> f64 {
    let x = if upper { -x } else { x };
    super::corr::t_cdf(x, p[0])
}
fn t_inv(q: f64, p: P) -> f64 {
    let v = p[0];
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if q == 0.5 {
        return 0.0;
    }
    if q == 0.0 {
        return f64::NEG_INFINITY;
    }
    if q == 1.0 {
        return f64::INFINITY;
    }
    let (lo, neg) = if q < 0.5 { (q, true) } else { (1.0 - q, false) };
    // lo = 0.5 * I_{v/(v+t^2)}(v/2, 1/2)
    let z = betaincinv(2.0 * lo, v / 2.0, 0.5, false);
    let t = (v * (1.0 / z - 1.0)).sqrt();
    if neg { -t } else { t }
}

fn f_pdf(x: f64, p: P) -> f64 {
    let (a, b) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return if a < 2.0 { f64::INFINITY } else if a == 2.0 { 1.0 } else { 0.0 };
    }
    let lg = 0.5 * (a * (a * x).ln() + b * b.ln() - (a + b) * (a * x + b).ln()) - x.ln() - crate::builtins_specfun::betaln(a / 2.0, b / 2.0);
    let _ = lg;
    ((a / 2.0) * (a / b).ln() + (a / 2.0 - 1.0) * x.ln() - ((a + b) / 2.0) * (a * x / b).ln_1p() - crate::builtins_specfun::betaln(a / 2.0, b / 2.0)).exp()
}
fn f_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if x.is_infinite() {
        return if upper { 0.0 } else { 1.0 };
    }
    let z = p[0] * x / (p[0] * x + p[1]);
    if z.is_nan() {
        return f64::NAN;
    }
    let (l, u) = betainc_pair(z, p[0] / 2.0, p[1] / 2.0);
    if upper { u } else { l }
}
fn f_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if q == 1.0 {
        return f64::INFINITY;
    }
    let z = betaincinv(q, p[0] / 2.0, p[1] / 2.0, false);
    p[1] * z / (p[0] * (1.0 - z))
}

fn wbl_pdf(x: f64, p: P) -> f64 {
    let (a, b) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return if b < 1.0 { f64::INFINITY } else if b == 1.0 { 1.0 / a } else { 0.0 };
    }
    b / a * (x / a).powf(b - 1.0) * (-(x / a).powf(b)).exp()
}
fn wbl_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    let e = (x / p[0]).powf(p[1]);
    if upper { (-e).exp() } else { -(-e).exp_m1() }
}
fn wbl_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) { f64::NAN } else { p[0] * (-(-q).ln_1p()).powf(1.0 / p[1]) }
}

fn rayl_pdf(x: f64, p: P) -> f64 {
    let b = p[0];
    if x.is_nan() {
        f64::NAN
    } else if x < 0.0 {
        0.0
    } else {
        x / (b * b) * (-x * x / (2.0 * b * b)).exp()
    }
}
fn rayl_cdf(x: f64, p: P, upper: bool) -> f64 {
    let b = p[0];
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    let e = -x * x / (2.0 * b * b);
    if upper { e.exp() } else { -e.exp_m1() }
}
fn rayl_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) { f64::NAN } else { p[0] * (-2.0 * (-q).ln_1p()).sqrt() }
}

fn ev_pdf(x: f64, p: P) -> f64 {
    let z = (x - p[0]) / p[1];
    z.exp() * (-z.exp()).exp() / p[1]
}
fn ev_cdf(x: f64, p: P, upper: bool) -> f64 {
    let e = ((x - p[0]) / p[1]).exp();
    if upper { (-e).exp() } else { -(-e).exp_m1() }
}
fn ev_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) { f64::NAN } else { p[0] + p[1] * (-(-q).ln_1p()).ln() }
}

const EULER: f64 = 0.577_215_664_901_532_9;

fn gev_pdf(x: f64, p: P) -> f64 {
    let (k, s, m) = (p[0], p[1], p[2]);
    if x.is_nan() {
        return f64::NAN;
    }
    let z = (x - m) / s;
    if k == 0.0 {
        return (-z - (-z).exp()).exp() / s;
    }
    let t = 1.0 + k * z;
    if t > 0.0 {
        let tk = t.powf(-1.0 / k);
        t.powf(-1.0 - 1.0 / k) * (-tk).exp() / s
    } else if t == 0.0 && k < -1.0 {
        f64::INFINITY
    } else if t == 0.0 && k == -1.0 {
        1.0 / s
    } else {
        0.0
    }
}
fn gev_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (k, s, m) = (p[0], p[1], p[2]);
    if x.is_nan() {
        return f64::NAN;
    }
    let z = (x - m) / s;
    let lower = if k == 0.0 {
        (-(-z).exp()).exp()
    } else {
        let t = 1.0 + k * z;
        if t > 0.0 {
            (-t.powf(-1.0 / k)).exp()
        } else if k > 0.0 {
            0.0
        } else {
            1.0
        }
    };
    if upper {
        // 1 - exp(-w) computed without cancellation
        if k == 0.0 {
            return -(-(-z).exp()).exp_m1();
        }
        let t = 1.0 + k * z;
        if t > 0.0 {
            return -(-t.powf(-1.0 / k)).exp_m1();
        }
        return 1.0 - lower;
    }
    lower
}
fn gev_inv(q: f64, p: P) -> f64 {
    let (k, s, m) = (p[0], p[1], p[2]);
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if k == 0.0 {
        return m - s * (-q.ln()).ln();
    }
    m + s * ((-q.ln()).powf(-k) - 1.0) / k
}

fn gp_pdf(x: f64, p: P) -> f64 {
    let (k, s, th) = (p[0], p[1], p[2]);
    if x.is_nan() {
        return f64::NAN;
    }
    let z = (x - th) / s;
    if z < 0.0 {
        return 0.0;
    }
    if k == 0.0 {
        return (-z).exp() / s;
    }
    if k > 0.0 || z <= -1.0 / k {
        (1.0 + k * z).powf(-1.0 - 1.0 / k) / s
    } else {
        0.0
    }
}
fn gp_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (k, s, th) = (p[0], p[1], p[2]);
    if x.is_nan() {
        return f64::NAN;
    }
    let z = (x - th) / s;
    if z <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    let (lo, hi) = if k == 0.0 {
        (-(-z).exp_m1(), (-z).exp())
    } else if k > 0.0 || z < -1.0 / k {
        let sv = (1.0 + k * z).powf(-1.0 / k);
        (1.0 - sv, sv)
    } else {
        (1.0, 0.0)
    };
    if upper { hi } else { lo }
}
fn gp_inv(q: f64, p: P) -> f64 {
    let (k, s, th) = (p[0], p[1], p[2]);
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if k == 0.0 {
        return th - s * (-q).ln_1p();
    }
    th + s * ((1.0 - q).powf(-k) - 1.0) / k
}

// ---- discrete families ------------------------------------------------------------------------

fn bino_pdf(x: f64, p: P) -> f64 {
    let (n, pr) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    if !isint(x) || x < 0.0 || x > n {
        return 0.0;
    }
    if pr == 0.0 {
        return if x == 0.0 { 1.0 } else { 0.0 };
    }
    if pr == 1.0 {
        return if x == n { 1.0 } else { 0.0 };
    }
    (lchoose(n, x) + x * pr.ln() + (n - x) * (-pr).ln_1p()).exp()
}
fn bino_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (n, pr) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    let k = x.floor();
    if k < 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if k >= n {
        return if upper { 0.0 } else { 1.0 };
    }
    // P(X <= k) = I_{1-p}(n-k, k+1)
    let (l, u) = betainc_pair(1.0 - pr, n - k, k + 1.0);
    if upper { u } else { l }
}
fn bino_inv(q: f64, p: P) -> f64 {
    let (n, pr) = (p[0], p[1]);
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if q == 1.0 {
        return n;
    }
    let mut cum = bino_pdf(0.0, p);
    let mut x = 0.0;
    let _ = pr;
    while cum + 10.0 * f64::EPSILON * cum.max(f64::MIN_POSITIVE) < q && x < n {
        x += 1.0;
        let pdf = bino_pdf(x, p);
        cum += pdf;
        if cum + 10.0 * f64::EPSILON * pdf >= q {
            break;
        }
    }
    x
}

fn poiss_pdf(x: f64, p: P) -> f64 {
    let l = p[0];
    if x.is_nan() {
        return f64::NAN;
    }
    if !isint(x) || x < 0.0 {
        return 0.0;
    }
    if l == 0.0 {
        return if x == 0.0 { 1.0 } else { 0.0 };
    }
    (-l + x * l.ln() - gammaln(x + 1.0)).exp()
}
fn poiss_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let k = x.floor();
    if k < 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if x == f64::INFINITY {
        return if upper { 0.0 } else { 1.0 };
    }
    if p[0] == 0.0 {
        return if upper { 0.0 } else { 1.0 };
    }
    let (l, u) = gammainc_pq(p[0], k + 1.0);
    // gammainc(lambda, k+1, 'upper') is P(X <= k)
    if upper { l } else { u }
}
fn poiss_inv(q: f64, p: P) -> f64 {
    let l = p[0];
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if q == 1.0 {
        return if l > 0.0 { f64::INFINITY } else { 0.0 };
    }
    if q == 0.0 || l == 0.0 {
        return 0.0;
    }
    let mut x = if l > 10.0 { (l.sqrt() * norm_inv(q) + l).ceil().max(0.0) } else { 0.0 };
    while poiss_cdf(x, p, false) < q {
        x += 1.0;
    }
    while x > 0.0 && poiss_cdf(x - 1.0, p, false) >= q {
        x -= 1.0;
    }
    x
}

fn geo_pdf(x: f64, p: P) -> f64 {
    let pr = p[0];
    if x.is_nan() {
        return f64::NAN;
    }
    if !isint(x) || x < 0.0 {
        return 0.0;
    }
    pr * (x * (-pr).ln_1p()).exp()
}
fn geo_cdf(x: f64, p: P, upper: bool) -> f64 {
    let pr = p[0];
    if x.is_nan() {
        return f64::NAN;
    }
    let k = x.floor();
    if k < 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    let e = (k + 1.0) * (-pr).ln_1p();
    if upper { e.exp() } else { -e.exp_m1() }
}
fn geo_inv(q: f64, p: P) -> f64 {
    let pr = p[0];
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    let eps = f64::EPSILON;
    let numer = if q == 0.0 { 0.0 } else if q == 1.0 { f64::NEG_INFINITY } else { (-q).ln_1p() };
    let denom = if pr == 1.0 { f64::NEG_INFINITY } else if pr <= eps { -pr } else { (-pr).ln_1p() };
    let x = ((numer / denom).abs() - 1.0 - eps.sqrt()).ceil();
    if q == 1.0 && pr == 1.0 { 0.0 } else if x < 0.0 { 0.0 } else { x }
}

fn nbin_pdf(x: f64, p: P) -> f64 {
    let (r, pr) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    if !isint(x) || x < 0.0 {
        return 0.0;
    }
    (gammaln(r + x) - gammaln(r) - gammaln(x + 1.0) + r * pr.ln() + x * (-pr).ln_1p()).exp()
}
fn nbin_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (r, pr) = (p[0], p[1]);
    if x.is_nan() {
        return f64::NAN;
    }
    let k = x.floor();
    if k < 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if x == f64::INFINITY {
        return if upper { 0.0 } else { 1.0 };
    }
    let (l, u) = betainc_pair(pr, r, k + 1.0);
    if upper { u } else { l }
}
fn nbin_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    if q == 1.0 {
        return f64::INFINITY;
    }
    if q == 0.0 {
        return 0.0;
    }
    let mut x = 0.0;
    let mut cum = nbin_pdf(0.0, p);
    while cum + 10.0 * f64::EPSILON * cum < q {
        x += 1.0;
        cum = nbin_cdf(x, p, false);
    }
    x
}

fn hyge_pdf(x: f64, p: P) -> f64 {
    let (m, k, n) = (p[0], p[1], p[2]);
    if x.is_nan() {
        return f64::NAN;
    }
    if !isint(x) || x < 0.0 || x > n || x > k || n - x > m - k {
        return 0.0;
    }
    (lchoose(k, x) + lchoose(m - k, n - x) - lchoose(m, n)).exp()
}
fn hyge_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (m, k, n) = (p[0], p[1], p[2]);
    if x.is_nan() {
        return f64::NAN;
    }
    let top = x.floor().min(n).min(k);
    let lo = (n - (m - k)).max(0.0);
    let mut s = 0.0;
    let mut u = 0.0;
    let mut i = lo;
    while i <= n.min(k) {
        let v = hyge_pdf(i, p);
        if i <= top {
            s += v;
        } else {
            u += v;
        }
        i += 1.0;
    }
    if x < lo {
        return if upper { 1.0 } else { 0.0 };
    }
    if upper { u } else { s }
}
fn hyge_inv(q: f64, p: P) -> f64 {
    let (m, k, n) = (p[0], p[1], p[2]);
    if !(0.0..=1.0).contains(&q) {
        return f64::NAN;
    }
    let lo = (n - (m - k)).max(0.0);
    let hi = n.min(k);
    let mut cum = 0.0;
    let mut x = lo;
    while x < hi {
        cum += hyge_pdf(x, p);
        if cum + 10.0 * f64::EPSILON * cum >= q {
            return x;
        }
        x += 1.0;
    }
    hi
}

fn unid_pdf(x: f64, p: P) -> f64 {
    if isint(x) && x >= 1.0 && x <= p[0] {
        1.0 / p[0]
    } else {
        0.0
    }
}
fn unid_cdf(x: f64, p: P, upper: bool) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let v = (x.floor() / p[0]).clamp(0.0, 1.0);
    if upper { 1.0 - v } else { v }
}
fn unid_inv(q: f64, p: P) -> f64 {
    if q > 0.0 && q <= 1.0 { (p[0] * q).ceil() } else { f64::NAN }
}

// ---- validity and moments ---------------------------------------------------------------------

fn v_norm(p: P) -> bool {
    p[1] > 0.0 && !p[0].is_nan()
}
fn v_pos1(p: P) -> bool {
    p[0] > 0.0
}
fn v_unif(p: P) -> bool {
    p[0] < p[1]
}
fn v_pos2(p: P) -> bool {
    p[0] > 0.0 && p[1] > 0.0
}
fn v_logn(p: P) -> bool {
    p[1] > 0.0 && !p[0].is_nan()
}
fn v_ev(p: P) -> bool {
    p[1] > 0.0 && !p[0].is_nan()
}
fn v_gev(p: P) -> bool {
    p[1] > 0.0 && !p[0].is_nan() && !p[2].is_nan()
}
fn v_bino(p: P) -> bool {
    isint(p[0]) && p[0] >= 0.0 && (0.0..=1.0).contains(&p[1])
}
fn v_poiss(p: P) -> bool {
    p[0] >= 0.0
}
fn v_geo(p: P) -> bool {
    p[0] > 0.0 && p[0] <= 1.0
}
fn v_nbin(p: P) -> bool {
    p[0] > 0.0 && p[1] > 0.0 && p[1] <= 1.0
}
fn v_hyge(p: P) -> bool {
    isint(p[0]) && isint(p[1]) && isint(p[2]) && p[0] >= 0.0 && p[1] >= 0.0 && p[2] >= 0.0 && p[1] <= p[0] && p[2] <= p[0]
}
fn v_unid(p: P) -> bool {
    p[0] >= 1.0 && isint(p[0])
}

fn st_norm(p: P) -> (f64, f64) {
    (p[0], p[1] * p[1])
}
fn st_logn(p: P) -> (f64, f64) {
    ((p[0] + p[1] * p[1] / 2.0).exp(), (2.0 * p[0] + p[1] * p[1]).exp() * (p[1] * p[1]).exp_m1())
}
fn st_exp(p: P) -> (f64, f64) {
    (p[0], p[0] * p[0])
}
fn st_unif(p: P) -> (f64, f64) {
    ((p[0] + p[1]) / 2.0, (p[1] - p[0]).powi(2) / 12.0)
}
fn st_gam(p: P) -> (f64, f64) {
    (p[0] * p[1], p[0] * p[1] * p[1])
}
fn st_beta(p: P) -> (f64, f64) {
    let s = p[0] + p[1];
    (p[0] / s, p[0] * p[1] / (s * s * (s + 1.0)))
}
fn st_chi2(p: P) -> (f64, f64) {
    (p[0], 2.0 * p[0])
}
fn st_t(p: P) -> (f64, f64) {
    let v = p[0];
    (if v > 1.0 { 0.0 } else { f64::NAN }, if v > 2.0 { v / (v - 2.0) } else { f64::NAN })
}
fn st_f(p: P) -> (f64, f64) {
    let (a, b) = (p[0], p[1]);
    let m = if b > 2.0 { b / (b - 2.0) } else { f64::NAN };
    let v = if b > 4.0 { 2.0 * b * b * (a + b - 2.0) / (a * (b - 2.0).powi(2) * (b - 4.0)) } else { f64::NAN };
    (m, v)
}
fn st_wbl(p: P) -> (f64, f64) {
    let g1 = gamma_fn(1.0 + 1.0 / p[1]);
    let g2 = gamma_fn(1.0 + 2.0 / p[1]);
    (p[0] * g1, p[0] * p[0] * (g2 - g1 * g1))
}
fn st_rayl(p: P) -> (f64, f64) {
    (p[0] * (PI / 2.0).sqrt(), (4.0 - PI) / 2.0 * p[0] * p[0])
}
fn st_ev(p: P) -> (f64, f64) {
    (p[0] - EULER * p[1], PI * PI * p[1] * p[1] / 6.0)
}
fn st_gev(p: P) -> (f64, f64) {
    let (k, s, m) = (p[0], p[1], p[2]);
    let mean = if k == 0.0 {
        m + s * EULER
    } else if k < 1.0 {
        m + s * (gamma_fn(1.0 - k) - 1.0) / k
    } else {
        f64::INFINITY
    };
    let var = if k == 0.0 {
        s * s * PI * PI / 6.0
    } else if k < 0.5 {
        s * s * (gamma_fn(1.0 - 2.0 * k) - gamma_fn(1.0 - k).powi(2)) / (k * k)
    } else {
        f64::INFINITY
    };
    (mean, var)
}
fn st_gp(p: P) -> (f64, f64) {
    let (k, s, th) = (p[0], p[1], p[2]);
    let mean = if k < 1.0 { th + s / (1.0 - k) } else { f64::INFINITY };
    let var = if k < 0.5 { s * s / ((1.0 - k).powi(2) * (1.0 - 2.0 * k)) } else { f64::INFINITY };
    (mean, var)
}
fn st_bino(p: P) -> (f64, f64) {
    (p[0] * p[1], p[0] * p[1] * (1.0 - p[1]))
}
fn st_poiss(p: P) -> (f64, f64) {
    (p[0], p[0])
}
fn st_geo(p: P) -> (f64, f64) {
    ((1.0 - p[0]) / p[0], (1.0 - p[0]) / (p[0] * p[0]))
}
fn st_nbin(p: P) -> (f64, f64) {
    (p[0] * (1.0 - p[1]) / p[1], p[0] * (1.0 - p[1]) / (p[1] * p[1]))
}
fn st_hyge(p: P) -> (f64, f64) {
    let (m, k, n) = (p[0], p[1], p[2]);
    (n * k / m, n * (k / m) * (1.0 - k / m) * (m - n) / (m - 1.0))
}
fn st_unid(p: P) -> (f64, f64) {
    ((p[0] + 1.0) / 2.0, (p[0] * p[0] - 1.0) / 12.0)
}

pub const FAMILIES: &[Family] = &[
    Family { name: "norm", np: 2, defaults: &[Some(0.0), Some(1.0)], valid: v_norm, pdf: norm_pdf, cdf: norm_cdf_f, inv: |q, p| p[1] * norm_inv(q) + p[0], stat: st_norm },
    Family { name: "logn", np: 2, defaults: &[Some(0.0), Some(1.0)], valid: v_logn, pdf: logn_pdf, cdf: logn_cdf, inv: |q, p| norm_inv(q).mul_add(p[1], p[0]).exp(), stat: st_logn },
    Family { name: "exp", np: 1, defaults: &[Some(1.0)], valid: v_pos1, pdf: exp_pdf, cdf: exp_cdf, inv: exp_inv, stat: st_exp },
    Family { name: "unif", np: 2, defaults: &[Some(0.0), Some(1.0)], valid: v_unif, pdf: unif_pdf, cdf: unif_cdf, inv: |q, p| if (0.0..=1.0).contains(&q) { p[0] + q * (p[1] - p[0]) } else { f64::NAN }, stat: st_unif },
    Family { name: "gam", np: 2, defaults: &[None, Some(1.0)], valid: v_pos2, pdf: gam_pdf, cdf: gam_cdf, inv: gam_inv, stat: st_gam },
    Family { name: "beta", np: 2, defaults: &[None, None], valid: v_pos2, pdf: beta_pdf, cdf: beta_cdf, inv: beta_inv, stat: st_beta },
    Family { name: "chi2", np: 1, defaults: &[None], valid: v_pos1, pdf: chi2_pdf, cdf: chi2_cdf, inv: chi2_inv, stat: st_chi2 },
    Family { name: "t", np: 1, defaults: &[None], valid: v_pos1, pdf: t_pdf, cdf: t_cdf, inv: t_inv, stat: st_t },
    Family { name: "f", np: 2, defaults: &[None, None], valid: v_pos2, pdf: f_pdf, cdf: f_cdf, inv: f_inv, stat: st_f },
    Family { name: "wbl", np: 2, defaults: &[Some(1.0), Some(1.0)], valid: v_pos2, pdf: wbl_pdf, cdf: wbl_cdf, inv: wbl_inv, stat: st_wbl },
    Family { name: "rayl", np: 1, defaults: &[Some(1.0)], valid: v_pos1, pdf: rayl_pdf, cdf: rayl_cdf, inv: rayl_inv, stat: st_rayl },
    Family { name: "ev", np: 2, defaults: &[Some(0.0), Some(1.0)], valid: v_ev, pdf: ev_pdf, cdf: ev_cdf, inv: ev_inv, stat: st_ev },
    Family { name: "gev", np: 3, defaults: &[Some(0.0), Some(1.0), Some(0.0)], valid: v_gev, pdf: gev_pdf, cdf: gev_cdf, inv: gev_inv, stat: st_gev },
    Family { name: "gp", np: 3, defaults: &[Some(0.0), Some(1.0), Some(0.0)], valid: v_gev, pdf: gp_pdf, cdf: gp_cdf, inv: gp_inv, stat: st_gp },
    Family { name: "bino", np: 2, defaults: &[None, None], valid: v_bino, pdf: bino_pdf, cdf: bino_cdf, inv: bino_inv, stat: st_bino },
    Family { name: "poiss", np: 1, defaults: &[None], valid: v_poiss, pdf: poiss_pdf, cdf: poiss_cdf, inv: poiss_inv, stat: st_poiss },
    Family { name: "geo", np: 1, defaults: &[None], valid: v_geo, pdf: geo_pdf, cdf: geo_cdf, inv: geo_inv, stat: st_geo },
    Family { name: "nbin", np: 2, defaults: &[None, None], valid: v_nbin, pdf: nbin_pdf, cdf: nbin_cdf, inv: nbin_inv, stat: st_nbin },
    Family { name: "hyge", np: 3, defaults: &[None, None, None], valid: v_hyge, pdf: hyge_pdf, cdf: hyge_cdf, inv: hyge_inv, stat: st_hyge },
    Family { name: "unid", np: 1, defaults: &[None], valid: v_unid, pdf: unid_pdf, cdf: unid_cdf, inv: unid_inv, stat: st_unid },
];

// ---- entry point ------------------------------------------------------------------------------

fn table(fname: &str, situation: &str) -> Option<&'static str> {
    ERRS.iter().find(|(f, s, _)| *f == fname && *s == situation).map(|(_, _, m)| *m)
}

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    let (fam, kind) = FAMILIES
        .iter()
        .find_map(|f| name.strip_prefix(f.name).filter(|k| ["pdf", "cdf", "inv", "stat"].contains(k)).map(|k| (f, k)))
        .expect("family");
    let _ = up::<u8>;
    let is_stat = kind == "stat";
    let mut args: Vec<&Value> = a.iter().collect();
    let mut upper = false;
    if kind == "cdf" {
        if let Some(l) = args.last() {
            if l.is_text() && args.len() > 1 {
                if l.text()?.eq_ignore_ascii_case("upper") {
                    upper = true;
                    args.pop();
                } else {
                    return rt("Unrecognized value for last argument. Use 'upper' to compute the upper tail probability.");
                }
            }
        }
    }
    if args.is_empty() {
        let msg = table(name, "noargs").unwrap_or("Not enough input arguments.");
        return rt(msg);
    }
    let first = if is_stat { 0 } else { 1 };
    let given = args.len() - first;
    if given > fam.np {
        return rt("Too many input arguments.");
    }
    if given < fam.np && fam.defaults[given].is_none() {
        let sit = if given == 0 && !is_stat { "onlyx" } else { "noargs" };
        let msg = table(name, sit).or_else(|| table(name, "noargs")).unwrap_or("Not enough input arguments.");
        return rt(msg);
    }
    let mut mats: Vec<Matrix> = Vec::new();
    for v in &args {
        match v {
            Value::Num(m) => mats.push(m.clone()),
            other => return rt(format!("Undefined function '{name}' for input arguments of type '{}'.", other.class_name())),
        }
    }
    let single = mats.iter().any(is_single);
    if !is_stat && mats[0].kind == Kind::Char && kind == "pdf" {
        match fam.name {
            "beta" | "f" => return rt("Undefined function 'log' for input arguments of type 'char'."),
            "logn" => return rt("Non-scalar arguments must match in size."),
            _ => {}
        }
    }
    // Expand: all non-scalar arguments must share a size.
    let mut dims: Option<&Dims> = None;
    for m in &mats {
        if m.numel() != 1 {
            match dims {
                None => dims = Some(&m.dims),
                Some(d) if *d != m.dims => {
                    let sit = if upper { "mismatch_upper" } else { "mismatch" };
                    let msg = table(name, sit).or_else(|| table(name, "mismatch")).unwrap_or("Non-scalar arguments must match in size.");
                    return rt(msg);
                }
                _ => {}
            }
        }
    }
    let dims = dims.cloned().unwrap_or_else(|| vec![1, 1]);
    let n: usize = dims.iter().product();
    let at = |m: &Matrix, i: usize| if m.numel() == 1 { m.data[0] } else { m.data[i] };
    let pm = &mats[first..];
    let param = |i: usize, k: usize| -> f64 {
        if k < pm.len() { at(&pm[k], i) } else { fam.defaults[k].unwrap_or(f64::NAN) }
    };
    let mut out = Vec::with_capacity(n);
    let mut out2 = Vec::with_capacity(if is_stat { n } else { 0 });
    for i in 0..n {
        let p: Vec<f64> = (0..fam.np).map(|k| param(i, k)).collect();
        let ok = (fam.valid)(&p) && p.iter().all(|v| !v.is_nan());
        if is_stat {
            let (m, v) = if ok { (fam.stat)(&p) } else { (f64::NAN, f64::NAN) };
            out.push(m);
            out2.push(v);
            continue;
        }
        let x = at(&mats[0], i);
        out.push(if let Some(v) = override_invalid(fam.name, kind, x, &p, ok) {
            v
        } else if !ok || x.is_nan() {
            f64::NAN
        } else if single && fam.name == "unid" && kind == "inv" && x > 0.0 && x <= 1.0 {
            ((p[0] * x) as f32).ceil() as f64
        } else {
            match kind {
                "pdf" => (fam.pdf)(x, &p),
                "cdf" => (fam.cdf)(x, &p, upper),
                _ => (fam.inv)(x, &p),
            }
        });
    }
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(dims.clone(), out), single))];
    if is_stat && nargout > 1 {
        outs.push(Value::Num(finish(Matrix::from_dims(dims, out2), single)));
    }
    let _ = (LN_2, tail as fn(f64, bool) -> f64);
    Ok(Some(outs))
}

/// Results MATLAB gives even for parameters it calls invalid.
fn override_invalid(fam: &str, kind: &str, x: f64, p: P, ok: bool) -> Option<f64> {
    if ok {
        return if fam == "unid" && kind == "pdf" && x.is_nan() { Some(0.0) } else { None };
    }
    match (fam, kind) {
        ("gam", "pdf") if x == 0.0 && p[0] < 1.0 && p[1] > 0.0 => Some(f64::INFINITY),
        ("chi2", "pdf") if x == 0.0 && p[0] / 2.0 < 1.0 => Some(f64::INFINITY),
        ("rayl", "inv") if x == 1.0 => Some(f64::INFINITY),
        _ => None,
    }
}

// ---- scalar entry points for the other stats modules ------------------------------------------

pub fn tinv(q: f64, v: f64) -> f64 {
    if v.is_nan() || v <= 0.0 { f64::NAN } else { t_inv(q, &[v]) }
}
pub fn tcdf(x: f64, v: f64) -> f64 {
    super::corr::t_cdf(x, v)
}
pub fn norminv_std(q: f64) -> f64 {
    norm_inv(q)
}
pub fn normcdf_std(z: f64) -> f64 {
    norm_cdf(z)
}
/// Chi-square distribution function (`upper` for the upper tail).
pub fn chi2cdf(x: f64, v: f64, upper: bool) -> f64 {
    if v == 0.0 && !x.is_nan() {
        let lower = if x < 0.0 { 0.0 } else { 1.0 };
        return if upper { 1.0 - lower } else { lower };
    }
    if v.is_nan() || v <= 0.0 || x.is_nan() { f64::NAN } else { chi2_cdf(x, &[v], upper) }
}
pub fn chi2inv(q: f64, v: f64) -> f64 {
    if v == 0.0 && (0.0..=1.0).contains(&q) {
        return 0.0;
    }
    if v.is_nan() || v <= 0.0 { f64::NAN } else { chi2_inv(q, &[v]) }
}
pub fn fcdf(x: f64, d1: f64, d2: f64, upper: bool) -> f64 {
    if d1.is_nan() || d2.is_nan() || d1 <= 0.0 || d2 <= 0.0 || x.is_nan() { f64::NAN } else { f_cdf(x, &[d1, d2], upper) }
}
pub fn finv(q: f64, d1: f64, d2: f64) -> f64 {
    if d1.is_nan() || d2.is_nan() || d1 <= 0.0 || d2 <= 0.0 { f64::NAN } else { f_inv(q, &[d1, d2]) }
}
