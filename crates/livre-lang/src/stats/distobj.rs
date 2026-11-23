//! `makedist` and `fitdist`: the univariate probability distribution objects (`NormalDistribution`,
//! `GammaDistribution`, ...) with their methods.

use super::util::*;
use crate::builtins::Ret;
use crate::class::{ClassDef, Object};
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub const NAMES: &[&str] = &["makedist", "fitdist", "matlivre_dist_method"];

#[derive(Clone, Copy, PartialEq)]
enum Tr {
    Id,
    Log,
    Logit,
}

struct Fam {
    name: &'static str,
    dname: &'static str,
    prefix: &'static str,
    params: &'static [&'static str],
    descr: &'static [&'static str],
    defaults: &'static [f64],
    tr: &'static [Tr],
    discrete: bool,
}

const FAMS: &[Fam] = &[
    Fam { name: "Normal", dname: "Normal", prefix: "norm", params: &["mu", "sigma"], descr: &["location", "scale"], defaults: &[0.0, 1.0], tr: &[Tr::Id, Tr::Log], discrete: false },
    Fam { name: "Lognormal", dname: "Lognormal", prefix: "logn", params: &["mu", "sigma"], descr: &["mean of logarithmic values", "standard deviation of logarithmic values"], defaults: &[0.0, 1.0], tr: &[Tr::Id, Tr::Log], discrete: false },
    Fam { name: "Exponential", dname: "Exponential", prefix: "exp", params: &["mu"], descr: &["mean"], defaults: &[1.0], tr: &[Tr::Log], discrete: false },
    Fam { name: "Gamma", dname: "Gamma", prefix: "gam", params: &["a", "b"], descr: &["shape", "scale"], defaults: &[1.0, 1.0], tr: &[Tr::Log, Tr::Log], discrete: false },
    Fam { name: "Weibull", dname: "Weibull", prefix: "wbl", params: &["a", "b"], descr: &["scale", "shape"], defaults: &[1.0, 1.0], tr: &[Tr::Log, Tr::Log], discrete: false },
    Fam { name: "Beta", dname: "Beta", prefix: "beta", params: &["a", "b"], descr: &["first shape", "second shape"], defaults: &[1.0, 1.0], tr: &[Tr::Log, Tr::Log], discrete: false },
    Fam { name: "Poisson", dname: "Poisson", prefix: "poiss", params: &["lambda"], descr: &["mean"], defaults: &[1.0], tr: &[Tr::Log], discrete: true },
    Fam { name: "Binomial", dname: "Binomial", prefix: "bino", params: &["N", "p"], descr: &["number of trials", "probability of success"], defaults: &[1.0, 0.5], tr: &[Tr::Id, Tr::Logit], discrete: true },
    Fam { name: "NegativeBinomial", dname: "Negative Binomial", prefix: "nbin", params: &["R", "p"], descr: &["number of successes", "probability of success in a single trial"], defaults: &[1.0, 0.5], tr: &[Tr::Log, Tr::Logit], discrete: true },
    Fam { name: "Rayleigh", dname: "Rayleigh", prefix: "rayl", params: &["B"], descr: &["scale"], defaults: &[1.0], tr: &[Tr::Log], discrete: false },
    Fam { name: "ExtremeValue", dname: "Extreme Value", prefix: "ev", params: &["mu", "sigma"], descr: &["location", "scale"], defaults: &[0.0, 1.0], tr: &[Tr::Id, Tr::Log], discrete: false },
    Fam { name: "GeneralizedExtremeValue", dname: "Generalized Extreme Value", prefix: "gev", params: &["k", "sigma", "mu"], descr: &["shape", "scale", "location"], defaults: &[0.0, 1.0, 0.0], tr: &[Tr::Id, Tr::Log, Tr::Id], discrete: false },
    Fam { name: "GeneralizedPareto", dname: "Generalized Pareto", prefix: "gp", params: &["k", "sigma", "theta"], descr: &["tail index (shape)", "scale", "threshold (location)"], defaults: &[1.0, 1.0, 1.0], tr: &[Tr::Id, Tr::Log, Tr::Id], discrete: false },
    Fam { name: "Uniform", dname: "Uniform", prefix: "unif", params: &["Lower", "Upper"], descr: &["lower endpoint (minimum)", "upper endpoint (maximum)"], defaults: &[0.0, 1.0], tr: &[Tr::Id, Tr::Id], discrete: false },
];

const COMMON: &[&str] = &["DistributionName", "NumParameters", "ParameterNames", "ParameterDescription", "ParameterValues", "Truncation", "IsTruncated", "ParameterCovariance", "ParameterIsFixed", "InputData"];
const METHODS: &[&str] = &["pdf", "cdf", "icdf", "mean", "var", "std", "median", "iqr", "negloglik", "paramci", "truncate"];

fn fam_of(name: &str) -> Option<&'static Fam> {
    FAMS.iter().find(|f| f.name.eq_ignore_ascii_case(name))
}

fn class_def(f: &'static Fam) -> Rc<ClassDef> {
    thread_local! {
        static DEF: RefCell<HashMap<&'static str, Rc<ClassDef>>> = RefCell::new(HashMap::new());
    }
    DEF.with(|d| {
        if let Some(c) = d.borrow().get(f.name) {
            return c.clone();
        }
        let mut src = format!("classdef {}Distribution\n  properties\n", f.name);
        for p in f.params {
            src += &format!("    {p}\n");
        }
        for p in COMMON {
            if f.name == "Uniform" && (*p == "ParameterCovariance" || *p == "ParameterIsFixed" || *p == "InputData") {
                continue;
            }
            src += &format!("    {p}\n");
        }
        src += "  end\n  methods\n";
        for m in METHODS {
            src += &format!("    function varargout = {m}(obj, varargin)\n      varargout = matlivre_dist_method('{m}', nargout, obj, varargin{{:}});\n    end\n");
        }
        src += "  end\nend\n";
        let ast = crate::class::parse_classdef(&src).expect("built-in class source parses");
        let def = Rc::new(ClassDef::build(ast, &[]));
        d.borrow_mut().insert(f.name, def.clone());
        def
    })
}

fn pnames(f: &Fam) -> Value {
    Value::Cell(CellArray::row(f.params.iter().map(|s| Value::chars(s)).collect()))
}

fn make_object(f: &'static Fam, vals: &[f64], cov: Option<Matrix>, fixed: Vec<bool>, input: Option<Value>, trunc: Option<(f64, f64)>) -> Value {
    let k = f.params.len();
    let class = class_def(f);
    let mut props: Vec<Value> = Vec::new();
    for pn in class.props.iter().map(|p| p.name.clone()).collect::<Vec<_>>() {
        let v = if let Some(i) = f.params.iter().position(|p| *p == pn) {
            Value::num(vals[i])
        } else {
            match pn.as_str() {
                "DistributionName" => Value::chars(f.dname),
                "NumParameters" => Value::num(k as f64),
                "ParameterNames" => pnames(f),
                "ParameterDescription" => Value::Cell(CellArray::row(f.descr.iter().map(|s| Value::chars(s)).collect())),
                "ParameterValues" => Value::Num(Matrix::from_dims(vec![1, k], vals.to_vec())),
                "Truncation" => match trunc {
                    Some((a, b)) => Value::Num(Matrix::from_dims(vec![1, 2], vec![a, b])),
                    None => Value::empty(),
                },
                "IsTruncated" => Value::bool(trunc.is_some()),
                "ParameterCovariance" => Value::Num(cov.clone().unwrap_or_else(|| Matrix::from_dims(vec![k, k], vec![0.0; k * k]))),
                "ParameterIsFixed" => Value::Num(Matrix::logical(1, k, fixed.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect())),
                "InputData" => input.clone().unwrap_or_else(Value::empty),
                _ => Value::empty(),
            }
        };
        props.push(v);
    }
    Value::Object(Object { class, props: Rc::new(RefCell::new(props)) })
}

fn check_param(f: &Fam, i: usize, v: &Value) -> Res<f64> {
    let name = f.params[i].to_uppercase();
    let x = match v {
        Value::Num(m) if m.numel() == 1 && m.kind != Kind::Char => m.data[0],
        _ => f64::NAN,
    };
    let bad = |what: &str| -> Res<f64> { rt(format!("The parameter {name} must be {what}.")) };
    match (f.name, f.params[i]) {
        ("Uniform", _) => {
            if !x.is_finite() {
                return bad("a real finite numeric scalar");
            }
        }
        ("Binomial", "N") => {
            if !(x >= 0.0 && x.fract() == 0.0 && x.is_finite()) {
                return bad("a nonnegative integer scalar");
            }
        }
        (_, "p") => {
            if !(0.0..=1.0).contains(&x) {
                return bad("a real numeric scalar between 0 and 1");
            }
        }
        (_, "sigma") if matches!(f.name, "Normal" | "Lognormal") => {
            if !(x >= 0.0 && x.is_finite()) {
                return bad("a nonnegative finite numeric scalar");
            }
        }
        (_, "mu") if matches!(f.name, "Normal" | "Lognormal" | "ExtremeValue" | "GeneralizedExtremeValue") => {
            if !x.is_finite() {
                return bad("a real finite numeric scalar");
            }
        }
        (_, "k") | (_, "theta") => {
            if !x.is_finite() {
                return bad("a real finite numeric scalar");
            }
        }
        _ => {
            if !(x > 0.0 && x.is_finite()) {
                return bad("a positive finite numeric scalar");
            }
        }
    }
    Ok(x)
}

fn makedist(a: &[Value]) -> Ret {
    if a.is_empty() {
        let names: Vec<Value> = FAMS.iter().map(|f| Value::chars(f.name)).collect();
        let n = names.len();
        return Ok(Some(vec![Value::Cell(CellArray::from_dims(vec![n, 1], names))]));
    }
    let name = a[0].text().unwrap_or_default();
    let Some(f) = fam_of(&name) else { return rt(format!("Distribution name '{name}' is unrecognized.")) };
    let rest = &a[1..];
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let mut vals = f.defaults.to_vec();
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        match f.params.iter().position(|p| p.eq_ignore_ascii_case(&key)) {
            Some(i) => vals[i] = check_param(f, i, &pair[1])?,
            None => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    if f.name == "Uniform" && vals[0] >= vals[1] {
        return rt("The Lower parameter must be less than the Upper parameter.");
    }
    let k = f.params.len();
    Ok(Some(vec![make_object(f, &vals, None, vec![true; k], None, None)]))
}

fn call_f(it: &mut Interpreter, name: &str, args: Vec<Value>) -> Res<Matrix> {
    let r = it.call_function(name, args, 1)?;
    match r.into_iter().next() {
        Some(Value::Num(m)) => Ok(m),
        _ => rt("distribution function failed"),
    }
}

fn pvec(f: &Fam, v: &[f64]) -> Vec<Value> {
    let _ = f;
    v.iter().map(|&x| Value::num(x)).collect()
}

fn fitdist(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.is_empty() {
        let names: Vec<Value> = FAMS.iter().filter(|f| f.name != "Uniform").map(|f| Value::chars(f.name)).collect();
        let n = names.len();
        return Ok(Some(vec![Value::Cell(CellArray::from_dims(vec![n, 1], names))]));
    }
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let name = a[1].text().unwrap_or_default();
    let Some(f) = fam_of(&name) else { return rt(format!("Distribution name '{name}' is unrecognized.")) };
    let Value::Num(xm) = &a[0] else { return rt("X must be a numeric column vector.") };
    if xm.dims.len() != 2 || xm.dims[1] != 1 || xm.dims[0] == 0 {
        return rt("X must be a numeric column vector.");
    }
    let rest = &a[2..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let (mut cens, mut freq, mut ntr): (Option<Vec<f64>>, Option<Vec<f64>>, Option<f64>) = (None, None, None);
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        match key.to_lowercase().as_str() {
            "censoring" => cens = Some(pair[1].as_matrix()?.data.clone()),
            "frequency" => freq = Some(pair[1].as_matrix()?.data.clone()),
            "ntrials" => ntr = Some(pair[1].scalar()?),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let n0 = xm.numel();
    if cens.as_ref().is_some_and(|c| c.len() != n0) || freq.as_ref().is_some_and(|c| c.len() != n0) {
        return rt("Sample data, frequency vector, and censoring vector must have the same size.");
    }
    if f.name == "Uniform" {
        return rt("The uniform distribution is not fittable.");
    }
    let keep: Vec<usize> = (0..n0).filter(|&i| !xm.data[i].is_nan()).collect();
    let x: Vec<f64> = keep.iter().map(|&i| xm.data[i]).collect();
    let cv: Option<Vec<f64>> = cens.as_ref().map(|c| keep.iter().map(|&i| c[i]).collect());
    let fv: Option<Vec<f64>> = freq.as_ref().map(|c| keep.iter().map(|&i| c[i]).collect());
    let n = x.len();
    let col = |v: &Vec<f64>| Value::Num(Matrix::from_dims(vec![v.len(), 1], v.clone()));
    let xcol = col(&x);
    let ccol = cv.as_ref().map(col).unwrap_or_else(Value::empty);
    let fcol = fv.as_ref().map(col).unwrap_or_else(Value::empty);
    let wts: Vec<f64> = fv.clone().unwrap_or_else(|| vec![1.0; n]);
    let wsum: f64 = wts.iter().sum();
    let fit_call = |it: &mut Interpreter, func: &str, with_cf: bool| -> Res<Vec<f64>> {
        let args = if with_cf { vec![xcol.clone(), Value::empty(), ccol.clone(), fcol.clone()] } else { vec![xcol.clone()] };
        Ok(call_f(it, func, args)?.data)
    };
    let mut ntrials = 1.0;
    let params: Vec<f64> = match f.name {
        "Normal" => {
            if cv.is_none() {
                let mu = x.iter().zip(&wts).map(|(v, w)| v * w).sum::<f64>() / wsum;
                let var = x.iter().zip(&wts).map(|(v, w)| w * (v - mu).powi(2)).sum::<f64>() / (wsum - 1.0);
                vec![mu, var.sqrt()]
            } else {
                fit_call(it, "normfit", true)?
            }
        }
        "Lognormal" => fit_call(it, "lognfit", true)?,
        "Exponential" => fit_call(it, "expfit", true)?,
        "Gamma" => fit_call(it, "gamfit", true)?,
        "Weibull" => fit_call(it, "wblfit", true)?,
        "Beta" => fit_call(it, "betafit", false)?,
        "Poisson" => vec![x.iter().zip(&wts).map(|(v, w)| v * w).sum::<f64>() / wsum],
        "Rayleigh" => vec![(x.iter().zip(&wts).map(|(v, w)| w * v * v).sum::<f64>() / (2.0 * wsum)).sqrt()],
        "ExtremeValue" => fit_call(it, "evfit", true)?,
        "NegativeBinomial" => fit_call(it, "nbinfit", false)?,
        "Binomial" => {
            ntrials = ntr.unwrap_or(1.0);
            if !(ntrials >= 0.0 && ntrials.fract() == 0.0) || x.iter().any(|&v| v > ntrials) {
                return rt("All N values must be non-negative integers at least as large as X.");
            }
            vec![ntrials, x.iter().zip(&wts).map(|(v, w)| v * w).sum::<f64>() / (ntrials * wsum)]
        }
        "GeneralizedExtremeValue" => fit_call(it, "gevfit", false)?,
        _ => fit_call(it, "gpfit", false)?,
    };
    let params: Vec<f64> = params.into_iter().take(f.params.len()).collect();
    let k = params.len();
    // covariance
    let pv = Matrix::from_dims(vec![1, k], params.clone());
    let like_cov = |it: &mut Interpreter, func: &str, p: Value| -> Res<Vec<f64>> {
        let args = vec![p, xcol.clone(), ccol.clone(), fcol.clone()];
        let r = it.call_function(func, args, 2)?;
        match r.get(1) {
            Some(Value::Num(m)) => Ok(m.data.clone()),
            _ => Ok(vec![f64::NAN; k * k]),
        }
    };
    let cov: Vec<f64> = match f.name {
        "Normal" => like_cov(it, "normlike", Value::Num(pv.clone()))?,
        "Lognormal" => like_cov(it, "lognlike", Value::Num(pv.clone()))?,
        "Exponential" => like_cov(it, "explike", Value::num(params[0]))?,
        "Gamma" => like_cov(it, "gamlike", Value::Num(pv.clone()))?,
        "Weibull" => like_cov(it, "wbllike", Value::Num(pv.clone()))?,
        "ExtremeValue" => like_cov(it, "evlike", Value::Num(pv.clone()))?,
        "Beta" => {
            let r = it.call_function("betalike", vec![Value::Num(pv.clone()), xcol.clone()], 2)?;
            match r.get(1) {
                Some(Value::Num(m)) => m.data.clone(),
                _ => vec![f64::NAN; 4],
            }
        }
        "GeneralizedExtremeValue" | "GeneralizedPareto" => {
            let func = if f.name == "GeneralizedExtremeValue" { "gevlike" } else { "gplike" };
            let r = it.call_function(func, vec![Value::Num(pv.clone()), xcol.clone()], 2)?;
            match r.get(1) {
                Some(Value::Num(m)) => m.data.clone(),
                _ => vec![f64::NAN; k * k],
            }
        }
        "Poisson" => vec![params[0] / wsum],
        "Rayleigh" => vec![params[0] * params[0] / (4.0 * wsum)],
        "Binomial" => vec![0.0, 0.0, 0.0, params[1] * (1.0 - params[1]) / (ntrials * wsum)],
        _ => {
            // negative binomial: observed information of (r, p)
            let (rr, pp) = (params[0], params[1]);
            let nn = x.len() as f64;
            let sumx: f64 = x.iter().sum();
            let mut tri = 0.0;
            for &v in &x {
                tri += it.call_function("psi", vec![Value::num(1.0), Value::num(rr + v)], 1)?[0].scalar()?;
            }
            tri -= nn * it.call_function("psi", vec![Value::num(1.0), Value::num(rr)], 1)?[0].scalar()?;
            let (h11, h12, h22) = (-tri, -nn / pp, nn * rr / (pp * pp) + sumx / (1.0 - pp).powi(2));
            let det = h11 * h22 - h12 * h12;
            vec![h22 / det, -h12 / det, -h12 / det, h11 / det]
        }
    };
    let mut fixed = vec![false; k];
    if f.name == "Binomial" {
        fixed[0] = true;
    }
    let mut input = StructArray::scalar();
    input.set("data", Value::Num(Matrix::from_dims(vec![n, 1], x.clone())));
    input.set("cens", cv.as_ref().map(col).unwrap_or_else(Value::empty));
    input.set("freq", fv.as_ref().map(col).unwrap_or_else(Value::empty));
    let covm = Matrix::from_dims(vec![k, k], if cov.len() == k * k { cov } else { vec![cov[0]; k * k] });
    Ok(Some(vec![make_object(f, &params, Some(covm), fixed, Some(Value::Struct(input)), None)]))
}

// ---- methods -----------------------------------------------------------------------------------

fn prop(o: &Object, name: &str) -> Value {
    let idx = o.class.props.iter().position(|p| p.name == name);
    idx.map(|i| o.props.borrow()[i].clone()).unwrap_or_else(Value::empty)
}

fn fam_of_obj(o: &Object) -> Option<&'static Fam> {
    FAMS.iter().find(|f| format!("{}Distribution", f.name) == o.class.name)
}

fn params_of(o: &Object) -> Vec<f64> {
    match prop(o, "ParameterValues") {
        Value::Num(m) => m.data.clone(),
        _ => vec![],
    }
}

fn trunc_of(o: &Object) -> Option<(f64, f64)> {
    match prop(o, "Truncation") {
        Value::Num(m) if m.numel() == 2 => Some((m.data[0], m.data[1])),
        _ => None,
    }
}

fn base_eval(it: &mut Interpreter, f: &Fam, kind: &str, x: &Matrix, p: &[f64], upper: bool) -> Res<Matrix> {
    let fname = format!("{}{}", f.prefix, kind);
    let mut args = vec![Value::Num(x.clone())];
    args.extend(pvec(f, p));
    if upper {
        args.push(Value::chars("upper"));
    }
    call_f(it, &fname, args)
}

fn cdf_at(it: &mut Interpreter, f: &Fam, p: &[f64], x: f64) -> Res<f64> {
    let r = base_eval(it, f, "cdf", &Matrix::scalar(x), p, false)?;
    Ok(r.data[0])
}

fn integrate2(g: &mut dyn FnMut(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
    // composite 16-point Gauss-Legendre
    const X: [f64; 8] = [0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499];
    const W: [f64; 8] = [0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767, 0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541];
    let panels = 64;
    let h = (b - a) / panels as f64;
    let (mut s1, mut s2) = (0.0, 0.0);
    for k in 0..panels {
        let (lo, hi) = (a + k as f64 * h, a + (k + 1) as f64 * h);
        let (c, r) = ((lo + hi) / 2.0, (hi - lo) / 2.0);
        for i in 0..8 {
            for xv in [c + r * X[i], c - r * X[i]] {
                let gv = g(xv);
                s1 += W[i] * r * xv * gv;
                s2 += W[i] * r * xv * xv * gv;
            }
        }
    }
    (s1, s2)
}

fn method(it: &mut Interpreter, a: &[Value]) -> Ret {
    let name = a[0].text().unwrap_or_default();
    let Value::Object(o) = &a[2] else { return rt("Invalid distribution object.") };
    let Some(f) = fam_of_obj(o) else { return rt("Invalid distribution object.") };
    let rest = &a[3..];
    let p = params_of(o);
    let tr = trunc_of(o);
    let single_out = |v: f64| Ok(Some(vec![Value::num(v)]));
    let mapm = |m: Matrix| Ok(Some(vec![Value::Num(m)]));
    let xarg = |v: Option<&Value>| -> Res<Matrix> {
        match v {
            Some(Value::Num(m)) => Ok(m.clone()),
            Some(_) => rt("X must be a numeric array."),
            None => rt("Not enough input arguments."),
        }
    };
    let (lo_c, hi_c) = match tr {
        Some((l, h)) => (cdf_at(it, f, &p, l)?, cdf_at(it, f, &p, h)?),
        None => (0.0, 1.0),
    };
    match name.as_str() {
        "pdf" => {
            let x = xarg(rest.first())?;
            let mut r = base_eval(it, f, "pdf", &x, &p, false)?;
            if let Some((l, h)) = tr {
                let z = hi_c - lo_c;
                for (v, xv) in r.data.iter_mut().zip(&x.data) {
                    *v = if *xv >= l && *xv <= h { *v / z } else { 0.0 };
                }
            }
            mapm(r)
        }
        "cdf" => {
            let x = xarg(rest.first())?;
            let upper = rest.get(1).map(|v| v.is_text() && v.text().map(|t| t.eq_ignore_ascii_case("upper")).unwrap_or(false)).unwrap_or(false);
            let mut r = base_eval(it, f, "cdf", &x, &p, false)?;
            if let Some((l, h)) = tr {
                let z = hi_c - lo_c;
                for (v, xv) in r.data.iter_mut().zip(&x.data) {
                    *v = if *xv < l { 0.0 } else if *xv > h { 1.0 } else { (*v - lo_c) / z };
                }
            }
            if upper {
                for v in r.data.iter_mut() {
                    *v = 1.0 - *v;
                }
            }
            mapm(r)
        }
        "icdf" => {
            let q = xarg(rest.first())?;
            let q2 = if tr.is_some() { Matrix::from_dims(q.dims.clone(), q.data.iter().map(|v| if (0.0..=1.0).contains(v) { lo_c + v * (hi_c - lo_c) } else { *v }).collect()) } else { q.clone() };
            mapm(base_eval(it, f, "inv", &q2, &p, false)?)
        }
        "mean" | "var" | "std" => {
            let (m, v) = if let Some((l, h)) = tr {
                let z = hi_c - lo_c;
                let pdf1 = |it: &mut Interpreter, x: f64| -> f64 { base_eval(it, f, "pdf", &Matrix::scalar(x), &p, false).map(|m| m.data[0]).unwrap_or(f64::NAN) };
                if f.discrete {
                    let (mut s1, mut s2) = (0.0, 0.0);
                    let mut x = l.max(0.0).ceil();
                    let top = if h.is_finite() { h } else { x + 100000.0 };
                    while x <= top {
                        let pm = pdf1(it, x) / z;
                        s1 += x * pm;
                        s2 += x * x * pm;
                        x += 1.0;
                    }
                    (s1, s2 - s1 * s1)
                } else {
                    let (a0, b0) = (if l.is_finite() { l } else { -40.0 * p.get(1).copied().unwrap_or(1.0).abs() + p[0] }, if h.is_finite() { h } else { 40.0 * p.get(1).copied().unwrap_or(1.0).abs() + p[0] });
                    let (m1, m2) = integrate2(&mut |x| pdf1(it, x) / z, a0, b0);
                    (m1, m2 - m1 * m1)
                }
            } else {
                let stat = if f.prefix == "poiss" { "poisstat".to_string() } else { format!("{}stat", f.prefix) };
                let r = it.call_function(&stat, pvec(f, &p), 2)?;
                (r[0].scalar()?, r[1].scalar()?)
            };
            single_out(match name.as_str() {
                "mean" => m,
                "var" => v,
                _ => v.sqrt(),
            })
        }
        "median" | "iqr" => {
            let q = if name == "median" { vec![0.5] } else { vec![0.25, 0.75] };
            let qm: Vec<f64> = q.iter().map(|v| if tr.is_some() { lo_c + v * (hi_c - lo_c) } else { *v }).collect();
            let r = base_eval(it, f, "inv", &Matrix::from_dims(vec![1, qm.len()], qm), &p, false)?;
            single_out(if name == "median" { r.data[0] } else { r.data[1] - r.data[0] })
        }
        "negloglik" => {
            let Value::Struct(inp) = prop(o, "InputData") else { return rt("This distribution has no input data.") };
            let data = match inp.get("data") {
                Some(Value::Num(m)) => m.clone(),
                _ => return rt("This distribution has no input data."),
            };
            let cens = match inp.get("cens") {
                Some(Value::Num(m)) if !m.is_empty() => m.data.clone(),
                _ => vec![0.0; data.numel()],
            };
            let freq = match inp.get("freq") {
                Some(Value::Num(m)) if !m.is_empty() => m.data.clone(),
                _ => vec![1.0; data.numel()],
            };
            let pdf = base_eval(it, f, "pdf", &data, &p, false)?;
            let sf = if cens.iter().any(|&c| c != 0.0) { Some(base_eval(it, f, "cdf", &data, &p, true)?) } else { None };
            let mut nll = 0.0;
            for i in 0..data.numel() {
                let l = if cens[i] != 0.0 { sf.as_ref().unwrap().data[i].ln() } else { pdf.data[i].ln() };
                nll -= freq[i] * l;
            }
            single_out(nll)
        }
        "paramci" => {
            let mut alpha = 0.05;
            for pair in rest.chunks(2) {
                let key = pair[0].text().unwrap_or_default();
                if key.eq_ignore_ascii_case("alpha") {
                    alpha = pair.get(1).and_then(|v| v.scalar().ok()).unwrap_or(0.05);
                    if !(alpha > 0.0 && alpha < 1.0) {
                        return rt("ALPHA must be a scalar between 0 and 1.");
                    }
                } else {
                    return rt(format!("Invalid parameter name: {key}."));
                }
            }
            let Value::Struct(inp) = prop(o, "InputData") else { return rt("This distribution has no input data.") };
            let get = |n: &str| match inp.get(n) {
                Some(Value::Num(m)) if !m.is_empty() => Value::Num(m.clone()),
                _ => Value::empty(),
            };
            let (data, cens, freq) = (get("data"), get("cens"), get("freq"));
            let k = p.len();
            let fitfn = match f.name {
                "Normal" => Some(("normfit", true)),
                "Lognormal" => Some(("lognfit", true)),
                "Exponential" => Some(("expfit", true)),
                "Gamma" => Some(("gamfit", true)),
                "Weibull" => Some(("wblfit", true)),
                "Beta" => Some(("betafit", false)),
                "Poisson" => Some(("poissfit", false)),
                "Rayleigh" => Some(("raylfit", false)),
                "ExtremeValue" => Some(("evfit", true)),
                "NegativeBinomial" => Some(("nbinfit", false)),
                "GeneralizedExtremeValue" => Some(("gevfit", false)),
                "GeneralizedPareto" => Some(("gpfit", false)),
                _ => None,
            };
            if let Some((func, cf)) = fitfn {
                let args = if cf { vec![data, Value::num(alpha), cens, freq] } else if func == "poissfit" || func == "raylfit" || func == "betafit" || func == "nbinfit" || func == "gevfit" || func == "gpfit" { vec![data, Value::num(alpha)] } else { vec![data] };
                if func == "normfit" || func == "lognfit" {
                    let r = it.call_function(func, args.clone(), 4)?;
                    if let (Some(Value::Num(a1)), Some(Value::Num(a2))) = (r.get(2), r.get(3)) {
                        let mut d = a1.data.clone();
                        d.extend(a2.data.iter().copied());
                        return mapm(Matrix::from_dims(vec![2, 2], d));
                    }
                }
                let r = it.call_function(func, args, 2)?;
                if let Some(Value::Num(ci)) = r.get(1) {
                    if ci.numel() == 2 * k {
                        return mapm(ci.clone());
                    }
                }
            }
            let Value::Num(cov) = prop(o, "ParameterCovariance") else { return rt("This distribution has no parameter covariance.") };
            let z = super::dist::norminv_std(1.0 - alpha / 2.0);
            let fixed = match prop(o, "ParameterIsFixed") {
                Value::Num(m) => m.data.clone(),
                _ => vec![0.0; k],
            };
            let mut out = vec![0.0; 2 * k];
            for i in 0..k {
                if fixed[i] != 0.0 {
                    out[i * 2] = p[i];
                    out[i * 2 + 1] = p[i];
                    continue;
                }
                let se = cov.data[i * k + i].sqrt();
                let (lo, hi) = match f.tr[i] {
                    Tr::Id => (p[i] - z * se, p[i] + z * se),
                    Tr::Log => {
                        let s = se / p[i];
                        ((p[i].ln() - z * s).exp(), (p[i].ln() + z * s).exp())
                    }
                    Tr::Logit => {
                        let s = se / (p[i] * (1.0 - p[i]));
                        let t = (p[i] / (1.0 - p[i])).ln();
                        (1.0 / (1.0 + (-(t - z * s)).exp()), 1.0 / (1.0 + (-(t + z * s)).exp()))
                    }
                };
                out[i * 2] = lo;
                out[i * 2 + 1] = hi;
            }
            mapm(Matrix::from_dims(vec![2, k], out))
        }
        _ => {
            // truncate
            if rest.len() < 2 {
                return rt("Not enough input arguments.");
            }
            let (l, h) = (rest[0].scalar()?, rest[1].scalar()?);
            if !(l < h) {
                return rt("The LOWER parameter must be less than the UPPER parameter.");
            }
            let (nl, nh) = match tr {
                Some((a0, b0)) => (l.max(a0), h.min(b0)),
                None => (l, h),
            };
            let cov = match prop(o, "ParameterCovariance") {
                Value::Num(m) => Some(m),
                _ => None,
            };
            let fixed = match prop(o, "ParameterIsFixed") {
                Value::Num(m) => m.data.iter().map(|&v| v != 0.0).collect(),
                _ => vec![true; p.len()],
            };
            let input = match prop(o, "InputData") {
                Value::Struct(s) => Some(Value::Struct(s)),
                _ => None,
            };
            Ok(Some(vec![make_object(f, &p, cov, fixed, input, Some((nl, nh)))]))
        }
    }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    match name {
        "makedist" => makedist(a),
        "fitdist" => fitdist(it, a),
        _ => {
            let outs = method(it, a)?.unwrap_or_default();
            Ok(Some(vec![Value::Cell(CellArray::row(outs))]))
        }
    }
}
