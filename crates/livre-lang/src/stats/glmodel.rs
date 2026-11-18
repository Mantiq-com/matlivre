//! `fitglm` for matrix input and the `GeneralizedLinearModel` object (normal, binomial,
//! Poisson, gamma and inverse Gaussian responses; canonical and several other links).

use super::dist::{chi2cdf, fcdf, normcdf_std, norminv_std, tcdf, tinv};
use super::linmodel::{col, design, design_row_values, formula_text, g3, prop, strong_markup, table_of, term_name, terms_of};
use super::regress::{qr_pivot, solve_upper};
use super::util::*;
use crate::builtins::Ret;
use crate::builtins_specfun::{betaln, betainc_pair, gammaln};
use crate::class::{ClassDef, Object};
use crate::display::NumFormat;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub const NAMES: &[&str] = &["fitglm", "matlivre_glm_method"];

const PROPS: &[&str] = &[
    "Offset", "Residuals", "Fitted", "Diagnostics", "Dispersion", "DispersionEstimated", "Deviance", "Distribution", "Link", "SSE", "SST", "SSR", "NumEstimatedCoefficients", "Formula",
    "LogLikelihood", "DFE", "CoefficientCovariance", "CoefficientNames", "NumCoefficients", "Coefficients", "Rsquared", "ModelCriterion", "VariableInfo", "NumVariables",
    "VariableNames", "NumPredictors", "PredictorNames", "ResponseName", "NumObservations", "Steps", "ObservationInfo", "Variables", "ObservationNames",
];
const METHODS: &[&str] = &["predict", "feval", "coefCI"];

#[derive(Clone, Copy, PartialEq)]
enum Dist {
    Normal,
    Binomial,
    Poisson,
    Gamma,
    InvGaussian,
}

#[derive(Clone, Copy, PartialEq)]
enum Link {
    Identity,
    Log,
    Logit,
    Probit,
    Comploglog,
    Reciprocal,
    Loglog,
    Power2inv,
}

impl Dist {
    fn name(self) -> &'static str {
        match self {
            Dist::Normal => "Normal",
            Dist::Binomial => "Binomial",
            Dist::Poisson => "Poisson",
            Dist::Gamma => "Gamma",
            Dist::InvGaussian => "Inverse Gaussian",
        }
    }
    fn canonical(self) -> Link {
        match self {
            Dist::Normal => Link::Identity,
            Dist::Binomial => Link::Logit,
            Dist::Poisson => Link::Log,
            Dist::Gamma => Link::Reciprocal,
            Dist::InvGaussian => Link::Power2inv,
        }
    }
    fn sqrtvar(self, mu: f64, n: f64) -> f64 {
        match self {
            Dist::Normal => 1.0,
            Dist::Binomial => (mu * (1.0 - mu) / n).sqrt(),
            Dist::Poisson => mu.sqrt(),
            Dist::Gamma => mu,
            Dist::InvGaussian => mu.powf(1.5),
        }
    }
    fn dev(self, mu: f64, y: f64, n: f64) -> f64 {
        let ylogy = |a: f64, b: f64| if a == 0.0 { 0.0 } else { a * (a / b).ln() };
        match self {
            Dist::Normal => (y - mu).powi(2),
            Dist::Binomial => 2.0 * n * (ylogy(y, mu) + ylogy(1.0 - y, 1.0 - mu)),
            Dist::Poisson => 2.0 * (ylogy(y, mu) - (y - mu)),
            Dist::Gamma => 2.0 * (-(y / mu).ln() + (y - mu) / mu),
            Dist::InvGaussian => ((y - mu) / mu).powi(2) / y,
        }
    }
}

impl Link {
    fn name(self) -> &'static str {
        match self {
            Link::Identity => "identity",
            Link::Log => "log",
            Link::Logit => "logit",
            Link::Probit => "probit",
            Link::Comploglog => "comploglog",
            Link::Reciprocal => "reciprocal",
            Link::Loglog => "loglog",
            Link::Power2inv => "-2",
        }
    }
    fn link(self, mu: f64) -> f64 {
        match self {
            Link::Identity => mu,
            Link::Log => mu.ln(),
            Link::Logit => (mu / (1.0 - mu)).ln(),
            Link::Probit => norminv_std(mu),
            Link::Comploglog => (-(1.0 - mu).ln()).ln(),
            Link::Reciprocal => 1.0 / mu,
            Link::Loglog => (-mu.ln()).ln(),
            Link::Power2inv => mu.powi(-2),
        }
    }
    fn dlink(self, mu: f64) -> f64 {
        match self {
            Link::Identity => 1.0,
            Link::Log => 1.0 / mu,
            Link::Logit => 1.0 / (mu * (1.0 - mu)),
            Link::Probit => (2.0 * std::f64::consts::PI).sqrt() * (0.5 * norminv_std(mu).powi(2)).exp(),
            Link::Comploglog => 1.0 / ((mu - 1.0) * (1.0 - mu).ln()),
            Link::Reciprocal => -1.0 / (mu * mu),
            Link::Loglog => 1.0 / (mu * mu.ln()),
            Link::Power2inv => -2.0 * mu.powi(-3),
        }
    }
    fn ilink(self, eta: f64) -> f64 {
        match self {
            Link::Identity => eta,
            Link::Log => eta.exp(),
            Link::Logit => 1.0 / (1.0 + (-eta).exp()),
            Link::Probit => normcdf_std(eta),
            Link::Comploglog => 1.0 - (-eta.exp()).exp(),
            Link::Reciprocal => 1.0 / eta,
            Link::Loglog => (-(-eta).exp()).exp(),
            Link::Power2inv => eta.powf(-0.5),
        }
    }
    fn prefix(self) -> Option<&'static str> {
        match self {
            Link::Identity => None,
            Link::Power2inv => Some("power"),
            other => Some(other.name()),
        }
    }
}

fn class_def() -> Rc<ClassDef> {
    thread_local! {
        static DEF: RefCell<HashMap<u8, Rc<ClassDef>>> = RefCell::new(HashMap::new());
    }
    DEF.with(|d| {
        if let Some(c) = d.borrow().get(&0) {
            return c.clone();
        }
        let mut src = String::from("classdef GeneralizedLinearModel\n  properties\n");
        for p in PROPS {
            src += &format!("    {p}\n");
        }
        src += "  end\n  methods\n";
        for m in METHODS {
            src += &format!("    function varargout = {m}(obj, varargin)\n      varargout = matlivre_glm_method('{m}', nargout, obj, varargin{{:}});\n    end\n");
        }
        src += "  end\nend\n";
        let ast = crate::class::parse_classdef(&src).expect("built-in class source parses");
        let def = Rc::new(ClassDef::build(ast, &[]));
        d.borrow_mut().insert(0, def.clone());
        def
    })
}

struct Fit {
    beta: Vec<f64>,
    mu: Vec<f64>,
    eta: Vec<f64>,
    dev: f64,
    ssr: f64,
    covb_unscaled: Vec<Vec<f64>>,
}

/// Iteratively reweighted least squares, as `glmfit` runs it.
fn irls(xd: &[Vec<f64>], y: &[f64], nn: &[f64], dist: Dist, link: Link) -> Res<Fit> {
    let n = y.len();
    let p = xd.first().map_or(0, |r| r.len());
    let eps = f64::EPSILON;
    let mut mu: Vec<f64> = (0..n)
        .map(|i| match dist {
            Dist::Poisson => y[i] + 0.25,
            Dist::Binomial => (nn[i] * y[i] + 0.5) / (nn[i] + 1.0),
            Dist::Gamma | Dist::InvGaussian => y[i].max(eps),
            Dist::Normal => y[i],
        })
        .collect();
    let mut eta: Vec<f64> = mu.iter().map(|&m| link.link(m)).collect();
    let mut b = vec![0.0; p];
    let seps = eps.sqrt();
    let mut r_last: Vec<Vec<f64>> = vec![];
    let mut perm_last: Vec<usize> = vec![];
    for _iter in 0..=100 {
        let deta: Vec<f64> = mu.iter().map(|&m| link.dlink(m)).collect();
        let z: Vec<f64> = (0..n).map(|i| eta[i] + (y[i] - mu[i]) * deta[i]).collect();
        let sqrtirls: Vec<f64> = (0..n).map(|i| deta[i].abs() * dist.sqrtvar(mu[i], nn[i])).collect();
        let mut sqrtw: Vec<f64> = sqrtirls.iter().map(|s| 1.0 / s).collect();
        let wtol = sqrtw.iter().cloned().fold(0.0, f64::max) * eps.powf(2.0 / 3.0);
        for w in sqrtw.iter_mut() {
            if *w < wtol && *w != 0.0 {
                *w = wtol;
            }
        }
        let b_old = b.clone();
        let xw: Vec<Vec<f64>> = (0..n).map(|i| xd[i].iter().map(|v| v * sqrtw[i]).collect()).collect();
        let (q, r, perm) = qr_pivot(&xw, n, p);
        let rq: Vec<Vec<f64>> = r.iter().take(p).map(|row| row[..p].to_vec()).collect();
        let qty: Vec<f64> = (0..p).map(|c| (0..n).map(|i| q[i][c] * z[i] * sqrtw[i]).sum()).collect();
        let bp = solve_upper(&rq, &qty, p);
        b = vec![0.0; p];
        for (j, &pc) in perm.iter().enumerate() {
            b[pc] = bp[j];
        }
        r_last = rq;
        perm_last = perm;
        eta = (0..n).map(|i| xd[i].iter().zip(&b).map(|(a, c)| a * c).sum()).collect();
        mu = eta.iter().map(|&e| link.ilink(e)).collect();
        match dist {
            Dist::Binomial => mu.iter_mut().for_each(|m| *m = m.max(eps).min(1.0 - eps)),
            Dist::Poisson | Dist::Gamma | Dist::InvGaussian => mu.iter_mut().for_each(|m| *m = m.max(eps)),
            _ => {}
        }
        if !(0..p).any(|j| (b[j] - b_old[j]).abs() > 1e-6 * seps.max(b_old[j].abs())) {
            break;
        }
    }
    let dev: f64 = (0..n).map(|i| dist.dev(mu[i], y[i], nn[i])).sum();
    let ssr: f64 = (0..n)
        .map(|i| match dist {
            Dist::Normal => (y[i] - mu[i]).powi(2),
            Dist::Binomial => (y[i] - mu[i]).powi(2) / (mu[i] * (1.0 - mu[i]) / nn[i]),
            Dist::Poisson => (y[i] - mu[i]).powi(2) / mu[i],
            Dist::Gamma => ((y[i] - mu[i]) / mu[i]).powi(2),
            Dist::InvGaussian => ((y[i] - mu[i]) / mu[i].powf(1.5)).powi(2),
        })
        .sum();
    // (R'R)^-1 unpermuted
    let mut rinv = vec![vec![0.0; p]; p];
    for j in 0..p {
        let mut e = vec![0.0; p];
        e[j] = 1.0;
        let c = solve_upper(&r_last, &e, p);
        for i in 0..p {
            rinv[i][j] = c[i];
        }
    }
    let mut cp = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..p {
            cp[perm_last[i]][perm_last[j]] = (0..p).map(|m| rinv[i][m] * rinv[j][m]).sum();
        }
    }
    Ok(Fit { beta: b, mu, eta, dev, ssr, covb_unscaled: cp })
}

fn loglik(dist: Dist, mu: &[f64], y: &[f64], nn: &[f64], dev: f64, phi_hat: f64) -> f64 {
    let n = y.len() as f64;
    match dist {
        Dist::Normal => -0.5 * n * ((2.0 * std::f64::consts::PI * dev / n).ln() + 1.0),
        Dist::Binomial => (0..y.len())
            .map(|i| {
                let k = nn[i] * y[i];
                let lc = gammaln(nn[i] + 1.0) - gammaln(k + 1.0) - gammaln(nn[i] - k + 1.0);
                let t1 = if k == 0.0 { 0.0 } else { k * mu[i].ln() };
                let t2 = if nn[i] - k == 0.0 { 0.0 } else { (nn[i] - k) * (1.0 - mu[i]).ln() };
                lc + t1 + t2
            })
            .sum(),
        Dist::Poisson => (0..y.len()).map(|i| (if y[i] == 0.0 { 0.0 } else { y[i] * mu[i].ln() }) - mu[i] - gammaln(y[i] + 1.0)).sum(),
        Dist::Gamma => {
            let phi = phi_hat;
            let a = 1.0 / phi;
            (0..y.len()).map(|i| a * (a / mu[i]).ln() + (a - 1.0) * y[i].ln() - a * y[i] / mu[i] - gammaln(a)).sum()
        }
        Dist::InvGaussian => {
            let phi = phi_hat;
            (0..y.len()).map(|i| -0.5 * (2.0 * std::f64::consts::PI * phi * y[i].powi(3)).ln() - (y[i] - mu[i]).powi(2) / (2.0 * phi * mu[i] * mu[i] * y[i])).sum()
        }
    }
}

pub fn object_text(o: &Object) -> Option<String> {
    if o.class.name != "GeneralizedLinearModel" {
        return None;
    }
    let formula = prop(o, "Formula").text().unwrap_or_default();
    let Value::Table(t) = prop(o, "Coefficients") else { return None };
    let body = strong_markup(&t.body(NumFormat::Short), t.row_names.as_deref().unwrap_or(&[]));
    let n = prop(o, "NumObservations").scalar().ok()? as i64;
    let dfe = prop(o, "DFE").scalar().ok()? as i64;
    let disp = prop(o, "Dispersion").scalar().ok()?;
    let est = prop(o, "DispersionEstimated").scalar().ok()? != 0.0;
    let Value::Struct(d) = prop(o, "Distribution") else { return None };
    let dname = d.get("Name")?.text().ok()?;
    let nullstat = prop(o, "ModelFitVsNull");
    let _ = nullstat;
    let (stat, pv) = (null_stat(o)?, null_p(o)?);
    let mut s = format!("Generalized linear regression model:\n    {formula}\n    Distribution = {dname}\n\nEstimated Coefficients:\n{}\n\n\n{n} observations, {dfe} error degrees of freedom\n", body.trim_end_matches('\n'));
    let has_test = !stat.is_nan();
    if est {
        s += &format!("Estimated Dispersion: {}", g3(disp));
        if has_test {
            s += &format!("
F-statistic vs. constant model: {}, p-value = {}", g3(stat), g3(pv));
        }
    } else {
        s += &format!("Dispersion: {}", g3(disp));
        if has_test {
            s += &format!("
Chi^2-statistic vs. constant model: {}, p-value = {}", g3(stat), g3(pv));
        }
    }
    Some(s)
}

fn null_stat(o: &Object) -> Option<f64> {
    let Value::Struct(s) = prop(o, "Steps") else { return None };
    s.get("Stat")?.scalar().ok()
}
fn null_p(o: &Object) -> Option<f64> {
    let Value::Struct(s) = prop(o, "Steps") else { return None };
    s.get("Pvalue")?.scalar().ok()
}

fn fit(a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() < 2 {
        return rt("Y argument is required unless X is a dataset or table.");
    }
    let (Value::Num(x), Value::Num(ym)) = (&a[0], &a[1]) else { return rt("X and Y must be numeric matrices.") };
    let n = x.dims[0];
    let p = x.dims[1];
    let mut spec = "linear".to_string();
    let mut rest = &a[2..];
    if let Some(v) = rest.first() {
        if v.is_text() {
            let t = v.text()?;
            let low = t.to_lowercase();
            let named = ["distribution", "link", "intercept", "varnames", "dispersionflag", "offset", "weights"];
            if !named.iter().any(|nm| nm.starts_with(&low) && !low.is_empty()) {
                spec = t;
                rest = &rest[1..];
            }
        }
    }
    if rest.len() % 2 == 1 {
        return rt(format!("Formula character vector or string scalar is an invalid model formula: {}.", rest[0].text().unwrap_or_default()));
    }
    let mut dist = Dist::Normal;
    let mut link_opt: Option<Link> = None;
    let mut intercept = true;
    let mut var_names: Option<Vec<String>> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        match key.as_str() {
            "distribution" => {
                let t = pair[1].text().unwrap_or_default().to_lowercase();
                dist = match t.as_str() {
                    "normal" => Dist::Normal,
                    "binomial" => Dist::Binomial,
                    "poisson" => Dist::Poisson,
                    "gamma" => Dist::Gamma,
                    "inverse gaussian" => Dist::InvGaussian,
                    _ => return rt(format!("'{}' is not a valid value for the 'Distribution' argument.  Valid values are: 'Normal', 'Poisson', 'Binomial', 'Gamma', 'Inverse Gaussian'.", pair[1].text().unwrap_or_default())),
                };
            }
            "link" => {
                let t = pair[1].text().unwrap_or_default().to_lowercase();
                link_opt = Some(match t.as_str() {
                    "identity" => Link::Identity,
                    "log" => Link::Log,
                    "logit" => Link::Logit,
                    "probit" => Link::Probit,
                    "comploglog" => Link::Comploglog,
                    "reciprocal" => Link::Reciprocal,
                    "loglog" => Link::Loglog,
                    "canonical" => dist.canonical(),
                    _ => return rt("Unrecognized link function name."),
                });
            }
            "intercept" => intercept = pair[1].scalar().unwrap_or(1.0) != 0.0,
            "varnames" => {
                if let Value::Cell(c) = &pair[1] {
                    var_names = Some(c.data.iter().map(|v| v.text().unwrap_or_default()).collect());
                }
            }
            _ => return rt(format!("The '{}' option of fitglm is not implemented.", pair[0].text().unwrap_or_default())),
        }
    }
    let link = link_opt.unwrap_or_else(|| dist.canonical());
    // response: a vector, or [successes trials] for the binomial
    let (y, nn): (Vec<f64>, Vec<f64>) = if dist == Dist::Binomial && ym.dims.len() == 2 && ym.dims[1] == 2 && ym.dims[0] == n {
        let (c, t) = (&ym.data[..n], &ym.data[n..]);
        ((0..n).map(|i| c[i] / t[i]).collect(), t.to_vec())
    } else {
        if ym.numel() != n {
            return rt("Predictor and response variables must have the same length.");
        }
        (ym.data.clone(), vec![1.0; n])
    };
    match dist {
        Dist::Binomial if y.iter().any(|&v| !(0.0..=1.0).contains(&v)) => return rt("Y must contain values in the interval [0,N] for the binomial distribution."),
        Dist::Poisson if y.iter().any(|&v| v < 0.0) => return rt("Y must contain non-negative values for the Poisson distribution."),
        Dist::Gamma | Dist::InvGaussian if y.iter().any(|&v| v <= 0.0) => return rt(if dist == Dist::Gamma { "Y must contain positive values for the gamma distribution." } else { "Y must contain positive values for the inverse Gaussian distribution." }),
        _ => {}
    }
    let mut names: Vec<String> = (1..=p).map(|i| format!("x{i}")).collect();
    let mut resp = "y".to_string();
    if let Some(v) = var_names {
        if v.len() == p + 1 {
            names = v[..p].to_vec();
            resp = v[p].clone();
        }
    }
    let terms = terms_of(&spec, p, intercept).ok_or_else(|| crate::error::LangError::new(format!("Formula character vector or string scalar is an invalid model formula: {spec}."), 0))?;
    let has_int = terms.iter().any(|r| r.iter().all(|&e| e == 0));
    let cnames: Vec<String> = terms.iter().map(|r| term_name(r, &names)).collect();
    let xd = design(x, &terms);
    let k = terms.len();
    let f = irls(&xd, &y, &nn, dist, link)?;
    let dfe = (n as f64 - k as f64).max(0.0);
    let est = matches!(dist, Dist::Normal | Dist::Gamma | Dist::InvGaussian);
    let sfit2 = f.ssr / dfe;
    let dispersion = if est { sfit2 } else { 1.0 };
    let cov: Vec<Vec<f64>> = f.covb_unscaled.iter().map(|r| r.iter().map(|v| v * if est { sfit2 } else { 1.0 }).collect()).collect();
    let se: Vec<f64> = (0..k).map(|i| cov[i][i].sqrt()).collect();
    let tstat: Vec<f64> = f.beta.iter().zip(&se).map(|(b, s)| b / s).collect();
    let pval: Vec<f64> = tstat.iter().map(|t| if est { 2.0 * tcdf(-t.abs(), dfe) } else { 2.0 * normcdf_std(-t.abs()) }).collect();
    // null model: constant only
    let ones: Vec<Vec<f64>> = vec![vec![1.0]; n];
    let f0 = irls(&ones, &y, &nn, dist, link)?;
    let nf = n as f64;
    let ll = loglik(dist, &f.mu, &y, &nn, f.dev, sfit2);
    let ll0 = loglik(dist, &f0.mu, &y, &nn, f0.dev, f0.ssr / nf);
    let scale: Vec<f64> = if dist == Dist::Binomial { nn.clone() } else { vec![1.0; n] };
    let yv: Vec<f64> = (0..n).map(|i| y[i] * scale[i]).collect();
    let muv: Vec<f64> = (0..n).map(|i| f.mu[i] * scale[i]).collect();
    let ybar = yv.iter().sum::<f64>() / nf;
    let raw: Vec<f64> = (0..n).map(|i| yv[i] - muv[i]).collect();
    let sse: f64 = raw.iter().map(|v| v * v).sum();
    let sst: f64 = yv.iter().map(|v| (v - ybar).powi(2)).sum();
    let ssr: f64 = muv.iter().map(|m| (m - ybar).powi(2)).sum();
    let r2 = 1.0 - sse / sst;
    let r2adj = 1.0 - (1.0 - r2) * (nf - 1.0) / dfe;
    let llr = 1.0 - ll / ll0;
    let rdev = 1.0 - f.dev / f0.dev;
    let adjgen = (1.0 - (2.0 / nf * (ll0 - ll)).exp()) / (1.0 - (2.0 / nf * ll0).exp());
    let kk = k as f64;
    let aic = -2.0 * ll + 2.0 * kk;
    let aicc = aic + 2.0 * kk * (kk + 1.0) / (nf - kk - 1.0);
    let bic = -2.0 * ll + kk * nf.ln();
    let caic = -2.0 * ll + kk * (nf.ln() + 1.0);
    let df_model = kk - if has_int { 1.0 } else { 0.0 };
    let (stat, pv) = if !(has_int && k > 1) {
        (f64::NAN, f64::NAN)
    } else {
        let d = f0.dev - f.dev;
        if est {
            let fs = (d / df_model) / dispersion;
            (fs, fcdf(fs, df_model, dfe, true))
        } else {
            (d, chi2cdf(d, df_model, true))
        }
    };
    // residuals
    let pearson: Vec<f64> = (0..n).map(|i| (y[i] - f.mu[i]) / (dist.sqrtvar(f.mu[i], nn[i]) + if y[i] == f.mu[i] { 1.0 } else { 0.0 })).collect();
    let devres: Vec<f64> = (0..n).map(|i| (y[i] - f.mu[i]).signum() * dist.dev(f.mu[i], y[i], nn[i]).max(0.0).sqrt()).collect();
    let ansc: Vec<f64> = (0..n)
        .map(|i| match dist {
            Dist::Normal => y[i] - f.mu[i],
            Dist::Binomial => {
                let t = 2.0 / 3.0;
                let bt = betaln(t, t).exp();
                bt * (betainc_pair(y[i], t, t).0 - betainc_pair(f.mu[i], t, t).0) / ((f.mu[i] * (1.0 - f.mu[i])).powf(1.0 / 6.0) / nn[i].sqrt())
            }
            Dist::Poisson => 1.5 * ((y[i].powf(2.0 / 3.0) - f.mu[i].powf(2.0 / 3.0)) / f.mu[i].powf(1.0 / 6.0)),
            Dist::Gamma => 3.0 * (y[i].powf(1.0 / 3.0) - f.mu[i].powf(1.0 / 3.0)) / f.mu[i].powf(1.0 / 3.0),
            Dist::InvGaussian => (y[i].ln() - f.mu[i].ln()) / f.mu[i],
        })
        .collect();
    let rawr: Vec<f64> = if dist == Dist::Binomial { (0..n).map(|i| raw[i]).collect() } else { raw.clone() };
    let rows = Some(cnames.clone());
    let coef_t = table_of(&["Estimate", "SE", "tStat", "pValue"], vec![col(f.beta.clone()), col(se), col(tstat), col(pval)], rows)?;
    let resid_t = table_of(&["Raw", "LinearPredictor", "Pearson", "Anscombe", "Deviance"], vec![col(rawr), col((0..n).map(|i| (y[i] - f.mu[i]) * link.dlink(f.mu[i])).collect()), col(pearson), col(ansc), col(devres)], None)?;
    let fitted_t = if dist == Dist::Binomial {
        table_of(&["Response", "LinearPredictor", "Probability"], vec![col(muv.clone()), col(f.eta.clone()), col(f.mu.clone())], None)?
    } else {
        table_of(&["Response", "LinearPredictor"], vec![col(f.mu.clone()), col(f.eta.clone())], None)?
    };
    let mut vars = vec![];
    for j in 0..p {
        vars.push(col(x.data[j * n..(j + 1) * n].to_vec()));
    }
    vars.push(col(y.clone()));
    let mut vnames: Vec<&str> = names.iter().map(|s| s.as_str()).collect();
    vnames.push(resp.as_str());
    let variables_t = table_of(&vnames, vars, None)?;
    let mut rsq = StructArray::scalar();
    for (nm, v) in [("Ordinary", r2), ("Adjusted", r2adj), ("LLR", llr), ("Deviance", rdev), ("AdjGeneralized", adjgen)] {
        rsq.set(nm, Value::num(v));
    }
    let mut crit = StructArray::scalar();
    for (nm, v) in [("AIC", aic), ("AICc", aicc), ("BIC", bic), ("CAIC", caic)] {
        crit.set(nm, Value::num(v));
    }
    let mut dstruct = StructArray::scalar();
    dstruct.set("Name", Value::chars(dist.name()));
    let mut lstruct = StructArray::scalar();
    lstruct.set("Name", Value::chars(link.name()));
    let mut steps = StructArray::scalar();
    steps.set("Stat", Value::num(stat));
    steps.set("Pvalue", Value::num(pv));
    let covm = Value::Num(Matrix::from_dims(vec![k, k], (0..k).flat_map(|j| (0..k).map(|i| cov[i][j]).collect::<Vec<_>>()).collect()));
    let mut formula = formula_text(&resp, &terms, &names);
    if link == Link::Power2inv {
        // MATLAB prints the exponent of this link as the code 49 and drops the intercept's 1.
        formula = formula.replacen(&format!("{resp} ~ 1"), &format!("power({resp},49) ~ "), 1);
    } else if let Some(pref) = link.prefix() {
        formula = formula.replacen(&format!("{resp} ~"), &format!("{pref}({resp}) ~"), 1);
    }
    let cell_names = Value::Cell(CellArray::row(cnames.iter().map(|s| Value::chars(s)).collect()));
    let mut vals: Vec<Value> = Vec::new();
    for pn in PROPS {
        let v = match *pn {
            "Offset" => Value::num(0.0),
            "Residuals" => resid_t.clone(),
            "Fitted" => fitted_t.clone(),
            "Dispersion" => Value::num(dispersion),
            "DispersionEstimated" => Value::bool(est),
            "Deviance" => Value::num(f.dev),
            "Distribution" => Value::Struct(dstruct.clone()),
            "Link" => Value::Struct(lstruct.clone()),
            "SSE" => Value::num(sse),
            "SST" => Value::num(sst),
            "SSR" => Value::num(ssr),
            "NumEstimatedCoefficients" | "NumCoefficients" => Value::num(kk),
            "Formula" => Value::chars(&formula),
            "LogLikelihood" => Value::num(ll),
            "DFE" => Value::num(dfe),
            "CoefficientCovariance" => covm.clone(),
            "CoefficientNames" => cell_names.clone(),
            "Coefficients" => coef_t.clone(),
            "Rsquared" => Value::Struct(rsq.clone()),
            "ModelCriterion" => Value::Struct(crit.clone()),
            "NumVariables" => Value::num((p + 1) as f64),
            "VariableNames" => Value::Cell(CellArray::row(names.iter().chain(std::iter::once(&resp)).map(|s| Value::chars(s)).collect())),
            "NumPredictors" => Value::num((0..p).filter(|&j| terms.iter().any(|t| t[j] > 0)).count() as f64),
            "PredictorNames" => Value::Cell(CellArray::row(names.iter().map(|s| Value::chars(s)).collect())),
            "ResponseName" => Value::chars(&resp),
            "NumObservations" => Value::num(nf),
            "Steps" => Value::Struct(steps.clone()),
            "Variables" => variables_t.clone(),
            _ => Value::empty(),
        };
        vals.push(v);
    }
    let o = Object { class: class_def(), props: Rc::new(RefCell::new(vals)) };
    Ok(Some(vec![Value::Object(o)]))
}

fn method(_it: &mut Interpreter, a: &[Value]) -> Ret {
    let name = a[0].text().unwrap_or_default();
    let nargout = a[1].scalar().unwrap_or(1.0) as usize;
    if name == "feval" && nargout > 1 {
        return rt("Too many output arguments.");
    }
    let Value::Object(o) = &a[2] else { return rt("Invalid GeneralizedLinearModel.") };
    let rest = &a[3..];
    let Value::Table(coef) = prop(o, "Coefficients") else { return rt("Invalid GeneralizedLinearModel.") };
    let est: Vec<f64> = match coef.var("Estimate") {
        Some(Value::Num(m)) => m.data.clone(),
        _ => vec![],
    };
    let dfe = prop(o, "DFE").scalar().unwrap_or(f64::NAN);
    let dispest = prop(o, "DispersionEstimated").scalar().unwrap_or(0.0) != 0.0;
    let k = est.len();
    let covb: Vec<f64> = match prop(o, "CoefficientCovariance") {
        Value::Num(m) => m.data.clone(),
        _ => vec![],
    };
    let link_name = match prop(o, "Link") {
        Value::Struct(s) => s.get("Name").and_then(|v| v.text().ok()).unwrap_or_default(),
        _ => String::new(),
    };
    let link = match link_name.as_str() {
        "identity" => Link::Identity,
        "log" => Link::Log,
        "logit" => Link::Logit,
        "probit" => Link::Probit,
        "comploglog" => Link::Comploglog,
        "reciprocal" => Link::Reciprocal,
        "loglog" => Link::Loglog,
        _ => Link::Power2inv,
    };
    let crit = |alpha: f64| if dispest { tinv(1.0 - alpha / 2.0, dfe) } else { norminv_std(1.0 - alpha / 2.0) };
    match name.as_str() {
        "coefCI" => {
            let alpha = match rest.first() {
                Some(Value::Num(m)) if !m.is_empty() => m.data[0],
                _ => 0.05,
            };
            let t = tinv(1.0 - alpha / 2.0, dfe);
            let mut d = vec![0.0; 2 * k];
            for i in 0..k {
                let se = covb[i * k + i].sqrt();
                d[i] = est[i] - t * se;
                d[k + i] = est[i] + t * se;
            }
            Ok(Some(vec![Value::Cell(CellArray::row(vec![Value::Num(Matrix::from_dims(vec![k, 2], d))]))]))
        }
        _ => {
            let (xnew, opts): (Matrix, &[Value]) = if name == "feval" {
                let cols: Vec<Matrix> = rest.iter().take_while(|v| !v.is_text()).map(|v| v.as_matrix().cloned()).collect::<Result<_, _>>()?;
                let n = cols.first().map_or(0, |c| c.numel());
                let mut data = Vec::new();
                for c in &cols {
                    data.extend(c.data.iter().copied());
                }
                let used = cols.len();
                (Matrix::from_dims(vec![n, used], data), &rest[used..])
            } else {
                match rest.first() {
                    Some(Value::Num(m)) => (m.clone(), &rest[1..]),
                    _ => return rt("Not enough input arguments."),
                }
            };
            let mut alpha = 0.05;
            for pair in opts.chunks(2) {
                let key = pair[0].text().unwrap_or_default().to_lowercase();
                if key == "alpha" {
                    alpha = pair.get(1).and_then(|v| v.scalar().ok()).unwrap_or(0.05);
                } else {
                    return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default()));
                }
            }
            let (rows, _) = design_row_values(o, &xnew)?;
            let n = rows.len();
            let eta: Vec<f64> = rows.iter().map(|r| r.iter().zip(&est).map(|(a, b)| a * b).sum()).collect();
            let ypred: Vec<f64> = eta.iter().map(|&e| link.ilink(e)).collect();
            let t = crit(alpha);
            let mut ci = vec![0.0; 2 * n];
            for i in 0..n {
                let mut v = 0.0;
                for a1 in 0..k {
                    for b1 in 0..k {
                        v += rows[i][a1] * covb[b1 * k + a1] * rows[i][b1];
                    }
                }
                let (l, u) = (link.ilink(eta[i] - t * v.sqrt()), link.ilink(eta[i] + t * v.sqrt()));
                ci[i] = l.min(u);
                ci[n + i] = l.max(u);
            }
            Ok(Some(vec![Value::Cell(CellArray::row(vec![col(ypred), Value::Num(Matrix::from_dims(vec![n, 2], ci))]))]))
        }
    }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    if name == "matlivre_glm_method" {
        return Ok(Some(method(it, a)?.unwrap_or_default()));
    }
    fit(a)
}
