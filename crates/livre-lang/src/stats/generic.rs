//! The generic `pdf`, `cdf`, `icdf` and `random`: a distribution name routed to its family functions.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["pdf", "cdf", "icdf", "random"];

const DIST_NAMES: [&str; 25] = [
    "beta", "binomial", "chi-square", "extreme value", "exponential", "f", "gamma", "generalized extreme value", "generalized pareto", "geometric", "half normal", "hypergeometric", "lognormal", "negative binomial", "noncentral f",
    "noncentral t", "noncentral chi-square", "normal", "pearson", "poisson", "rayleigh", "t", "discrete uniform", "uniform", "weibull",
];

/// (canonical name, function prefix, number of distribution parameters).
const FAMILIES: [(&str, &str, usize); 23] = [
    ("beta", "beta", 2), ("binomial", "bino", 2), ("chi-square", "chi2", 1), ("exponential", "exp", 1), ("extreme value", "ev", 2), ("f", "f", 2), ("gamma", "gam", 2), ("generalized extreme value", "gev", 3), ("generalized pareto", "gp", 3), ("geometric", "geo", 1),
    ("hypergeometric", "hyge", 3), ("lognormal", "logn", 2), ("negative binomial", "nbin", 2), ("noncentral f", "ncf", 3), ("noncentral t", "nct", 2), ("noncentral chi-square", "ncx2", 2), ("normal", "norm", 2), ("poisson", "poiss", 1), ("rayleigh", "rayl", 1), ("t", "t", 1),
    ("discrete uniform", "unid", 1), ("uniform", "unif", 2), ("weibull", "wbl", 2),
];

fn canonical(name: &str) -> Res<String> {
    let n = name.to_lowercase();
    let hits: Vec<&&str> = DIST_NAMES.iter().filter(|d| !n.is_empty() && d.starts_with(&n)).collect();
    if hits.len() > 1 {
        return rt(format!("Ambiguous distribution name: '{name}'."));
    }
    if hits.len() == 1 {
        return Ok(hits[0].to_string());
    }
    Ok(match n.as_str() {
        "ev" => "extreme value",
        "gev" => "generalized extreme value",
        "gp" => "generalized pareto",
        "nbin" => "negative binomial",
        "wbl" => "weibull",
        "chi2" | "chisquare" => "chi-square",
        "hyge" => "hypergeometric",
        "ncf" => "noncentral f",
        "nct" => "noncentral t",
        "ncx2" => "noncentral chi-square",
        "unid" => "discrete uniform",
        _ => return Ok(n),
    }
    .to_string())
}

fn family(canon: &str) -> Option<(&'static str, usize)> {
    FAMILIES.iter().find(|f| f.0 == canon).map(|f| (f.1, f.2))
}

fn bad_name(who: &str, name: &str) -> LangErr {
    let _ = (who, name);
    unreachable!()
}
type LangErr = crate::error::LangError;

fn text_name(v: &Value, who: &str) -> Res<String> {
    if v.is_text() {
        return v.text();
    }
    rt(match who {
        "pdf" => "Requires the first input to be the name of a distribution.",
        "random" => "The NAME argument must be a distribution name.",
        _ => "First argument must be distribution name",
    })
}

pub fn call(it: &mut Interpreter, who: &str, a: &[Value], _nargout: usize) -> Ret {
    let _ = bad_name;
    if who == "random" {
        if a.is_empty() {
            return rt("Requires at least one input argument.");
        }
        let name = text_name(&a[0], who)?;
        let canon = canonical(&name)?;
        let Some((prefix, _)) = family(&canon) else { return rt(format!("Unrecognized distribution name: '{name}'.")) };
        if a[1..].iter().any(|v| !matches!(v, Value::Num(m) if m.kind != Kind::Char)) {
            return rt("Parameter values must be numeric.");
        }
        let r = it.call_function(&format!("{prefix}rnd"), a[1..].to_vec(), 1)?;
        return Ok(Some(r));
    }
    if a.len() < 2 {
        return rt(if who == "pdf" { "Requires at least two input arguments." } else { "Not enough input arguments" });
    }
    let name = text_name(&a[0], who)?;
    let canon = canonical(&name)?;
    let Some((prefix, np)) = family(&canon) else { return rt(format!("Distribution name '{name}' is unrecognized.")) };
    let mut rest: Vec<Value> = a[2..].to_vec();
    let mut tail: Option<Value> = None;
    if who == "cdf" {
        if let Some(l) = rest.last() {
            if l.is_text() && l.text()?.eq_ignore_ascii_case("upper") {
                tail = rest.pop();
            }
        }
    }
    let mut args = vec![a[1].clone()];
    for k in 0..np {
        args.push(rest.get(k).cloned().unwrap_or_else(|| Value::num(0.0)));
    }
    if let Some(t) = tail {
        args.push(t);
    }
    let fname = match who {
        "pdf" => format!("{prefix}pdf"),
        "cdf" => format!("{prefix}cdf"),
        _ => format!("{prefix}inv"),
    };
    Ok(Some(it.call_function(&fname, args, 1)?))
}
