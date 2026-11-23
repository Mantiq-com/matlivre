//! Noncentral chi-square, t and F distributions (series in Poisson weights of the central ones).

use crate::builtins_specfun::{betainc_pair, gammainc_pq, gammaln};

type P<'a> = &'a [f64];

/// Sum of `w_j * f(j)` over the Poisson(lambda) weights, from the mode outwards.
fn pois_sum(lambda: f64, f: &dyn Fn(usize) -> f64) -> f64 {
    if lambda == 0.0 {
        return f(0);
    }
    let lw = |j: usize| -lambda + j as f64 * lambda.ln() - gammaln(j as f64 + 1.0);
    let jm = lambda.floor() as usize;
    let mut sum = f(jm) * lw(jm).exp();
    let mut j = jm + 1;
    let mut wsum = lw(jm).exp();
    loop {
        let w = lw(j).exp();
        sum += w * f(j);
        wsum += w;
        if (w < 1e-18 * wsum.max(1e-300) && j > jm + 2) || j > jm + 200_000 {
            break;
        }
        j += 1;
    }
    let mut j = jm;
    while j > 0 {
        j -= 1;
        let w = lw(j).exp();
        sum += w * f(j);
        if w < 1e-18 * wsum.max(1e-300) {
            break;
        }
    }
    sum
}

fn chi2_pdf_k(x: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return if k == 2.0 {
            0.5
        } else if k < 2.0 {
            f64::INFINITY
        } else {
            0.0
        };
    }
    ((k / 2.0 - 1.0) * (x / 2.0).ln() - x / 2.0 - std::f64::consts::LN_2 - gammaln(k / 2.0)).exp()
}

fn gam_pq(x: f64, a: f64) -> (f64, f64) {
    if a == 0.0 {
        return if x > 0.0 { (1.0, 0.0) } else { (0.0, 1.0) };
    }
    if x <= 0.0 {
        return (0.0, 1.0);
    }
    if x == f64::INFINITY {
        return (1.0, 0.0);
    }
    gammainc_pq(x, a)
}

fn root(cdf: &dyn Fn(f64) -> f64, q: f64, mut lo: f64, mut hi: f64, two_sided: bool) -> f64 {
    let mut n = 0;
    while cdf(hi) < q && n < 2100 {
        hi = if hi == 0.0 { 1.0 } else { hi * 2.0 };
        n += 1;
    }
    if two_sided {
        n = 0;
        while cdf(lo) > q && n < 2100 {
            lo = if lo == 0.0 { -1.0 } else { lo * 2.0 };
            n += 1;
        }
    }
    for _ in 0..400 {
        let mid = 0.5 * (lo + hi);
        if mid == lo || mid == hi {
            break;
        }
        if cdf(mid) < q {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ---- noncentral chi-square ---------------------------------------------------------------------

pub fn v_ncx2(p: P) -> bool {
    p[0] >= 0.0 && p[1] >= 0.0
}

pub fn ncx2_pdf(x: f64, p: P) -> f64 {
    let (v, d) = (p[0], p[1]);
    if x < 0.0 || x == f64::INFINITY {
        return 0.0;
    }
    pois_sum(d / 2.0, &|j| chi2_pdf_k(x, v + 2.0 * j as f64))
}

pub fn ncx2_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (v, d) = (p[0], p[1]);
    if x == f64::INFINITY {
        return if upper { 0.0 } else { 1.0 };
    }
    if x == 0.0 && v == 0.0 {
        let lam = d / 2.0;
        return if upper { -(-lam).exp_m1() } else { (-lam).exp() };
    }
    if x < 0.0 || (x == 0.0 && v > 0.0) {
        return if upper { 1.0 } else { 0.0 };
    }
    pois_sum(d / 2.0, &|j| {
        let (pp, qq) = gam_pq(x / 2.0, v / 2.0 + j as f64);
        if upper {
            qq
        } else {
            pp
        }
    })
}

pub fn ncx2_inv(q: f64, p: P) -> f64 {
    let (v, d) = (p[0], p[1]);
    if !(0.0..=1.0).contains(&q) || q.is_nan() {
        return f64::NAN;
    }
    if q == 1.0 {
        return f64::INFINITY;
    }
    if q == 0.0 {
        return 0.0;
    }
    if v == 0.0 && ncx2_cdf(0.0, p, false) >= q {
        return 0.0;
    }
    root(&|x| ncx2_cdf(x, p, false), q, 0.0, (v + d).max(1.0), false)
}

pub fn st_ncx2(p: P) -> (f64, f64) {
    if p[0] == 0.0 {
        return (f64::NAN, f64::NAN);
    }
    (p[0] + p[1], 2.0 * (p[0] + 2.0 * p[1]))
}

// ---- noncentral t -------------------------------------------------------------------------------

pub fn v_nct(p: P) -> bool {
    p[0] > 0.0 && p[1].is_finite()
}

fn phi(z: f64) -> f64 {
    0.5 * crate::builtins_specfun::erfc(-z / std::f64::consts::SQRT_2)
}

fn nct_lower(t: f64, nu: f64, del: f64) -> f64 {
    if t == f64::INFINITY {
        return 1.0;
    }
    if t == f64::NEG_INFINITY {
        return 0.0;
    }
    if t < 0.0 {
        return 1.0 - nct_lower(-t, nu, -del);
    }
    if del == 0.0 {
        return super::corr::t_cdf(t, nu);
    }
    let x = t * t / (t * t + nu);
    let lam = del * del / 2.0;
    let sum = pois_sum(lam, &|j| {
        let jf = j as f64;
        let a = betainc_pair(x, jf + 0.5, nu / 2.0).0;
        let b = betainc_pair(x, jf + 1.0, nu / 2.0).0;
        a + del / std::f64::consts::SQRT_2 * (gammaln(jf + 1.0) - gammaln(jf + 1.5)).exp() * b
    });
    phi(-del) + 0.5 * sum
}

pub fn nct_cdf(x: f64, p: P, upper: bool) -> f64 {
    let lo = nct_lower(x, p[0], p[1]);
    if upper {
        1.0 - lo
    } else {
        lo
    }
}

pub fn nct_pdf(t: f64, p: P) -> f64 {
    let (nu, del) = (p[0], p[1]);
    if t.is_infinite() {
        return 0.0;
    }
    if t == 0.0 {
        return (gammaln((nu + 1.0) / 2.0) - gammaln(nu / 2.0) - 0.5 * (std::f64::consts::PI * nu).ln() - del * del / 2.0).exp();
    }
    nu / t * (nct_lower(t * (1.0 + 2.0 / nu).sqrt(), nu + 2.0, del) - nct_lower(t, nu, del))
}

pub fn nct_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) || q.is_nan() {
        return f64::NAN;
    }
    if q == 0.0 {
        return f64::NEG_INFINITY;
    }
    if q == 1.0 {
        return f64::INFINITY;
    }
    if p[1] == 0.0 {
        return super::dist::t_inv_pub(q, p[0]);
    }
    root(&|x| nct_lower(x, p[0], p[1]), q, -1.0, 1.0, true)
}

pub fn st_nct(p: P) -> (f64, f64) {
    let (nu, d) = (p[0], p[1]);
    let m = if nu > 1.0 { d * (nu / 2.0).sqrt() * (gammaln((nu - 1.0) / 2.0) - gammaln(nu / 2.0)).exp() } else { f64::NAN };
    let v = if nu > 2.0 { nu * (1.0 + d * d) / (nu - 2.0) - m * m } else { f64::NAN };
    (m, v)
}

// ---- noncentral F ---------------------------------------------------------------------------------

pub fn v_ncf(p: P) -> bool {
    p[0] > 0.0 && p[1] > 0.0 && p[2] >= 0.0
}

pub fn ncf_cdf(x: f64, p: P, upper: bool) -> f64 {
    let (v1, v2, d) = (p[0], p[1], p[2]);
    if x <= 0.0 {
        return if upper { 1.0 } else { 0.0 };
    }
    if x == f64::INFINITY {
        return if upper { 0.0 } else { 1.0 };
    }
    let y = v1 * x / (v2 + v1 * x);
    pois_sum(d / 2.0, &|j| {
        let (lo, up) = betainc_pair(y, v1 / 2.0 + j as f64, v2 / 2.0);
        if upper {
            up
        } else {
            lo
        }
    })
}

pub fn ncf_pdf(x: f64, p: P) -> f64 {
    let (v1, v2, d) = (p[0], p[1], p[2]);
    if x < 0.0 || x == f64::INFINITY {
        return 0.0;
    }
    if x == 0.0 {
        return if v1 > 2.0 {
            0.0
        } else if v1 == 2.0 {
            (-d / 2.0).exp()
        } else {
            f64::INFINITY
        };
    }
    let y = v1 * x / (v2 + v1 * x);
    let jac = v1 * v2 / ((v2 + v1 * x) * (v2 + v1 * x));
    let b = v2 / 2.0;
    pois_sum(d / 2.0, &|j| {
        let a = v1 / 2.0 + j as f64;
        let lb = (a - 1.0) * y.ln() + (b - 1.0) * (1.0 - y).ln() - (gammaln(a) + gammaln(b) - gammaln(a + b));
        lb.exp() * jac
    })
}

pub fn ncf_inv(q: f64, p: P) -> f64 {
    if !(0.0..=1.0).contains(&q) || q.is_nan() {
        return f64::NAN;
    }
    if q == 0.0 {
        return 0.0;
    }
    if q == 1.0 {
        return f64::INFINITY;
    }
    root(&|x| ncf_cdf(x, p, false), q, 0.0, 1.0, false)
}

pub fn st_ncf(p: P) -> (f64, f64) {
    let (v1, v2, d) = (p[0], p[1], p[2]);
    let m = if v2 > 2.0 { v2 * (v1 + d) / (v1 * (v2 - 2.0)) } else { f64::NAN };
    let v = if v2 > 4.0 { 2.0 * (v2 / v1).powi(2) * ((v1 + d).powi(2) + (v1 + 2.0 * d) * (v2 - 2.0)) / ((v2 - 2.0).powi(2) * (v2 - 4.0)) } else { f64::NAN };
    (m, v)
}
