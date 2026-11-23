//! `createns`, `ExhaustiveSearcher` and `KDTreeSearcher`: nearest-neighbour search objects. Their
//! `knnsearch` and `rangesearch` methods are the functions of the same names (see neighbors.rs).

use super::neighbors::{cols_stats, rows_of};
use super::util::*;
use crate::builtins::Ret;
use crate::class::{ClassDef, Object};
use crate::error::rt;
use crate::value::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub const NAMES: &[&str] = &["createns", "ExhaustiveSearcher", "KDTreeSearcher"];

const EXH_PROPS: &[&str] = &["X", "Distance", "DistParameter"];
const KD_PROPS: &[&str] = &["X", "Distance", "DistParameter", "BucketSize"];

fn class_def(class: &'static str, props: &'static [&'static str]) -> Rc<ClassDef> {
    thread_local! {
        static DEF: RefCell<HashMap<&'static str, Rc<ClassDef>>> = RefCell::new(HashMap::new());
    }
    DEF.with(|d| {
        if let Some(c) = d.borrow().get(class) {
            return c.clone();
        }
        let mut src = format!("classdef {class}\n  properties\n");
        for p in props {
            src += &format!("    {p}\n");
        }
        src += "  end\nend\n";
        let ast = crate::class::parse_classdef(&src).expect("built-in class source parses");
        let def = Rc::new(ClassDef::build(ast, &[]));
        d.borrow_mut().insert(class, def.clone());
        def
    })
}

const ALL: [&str; 11] = ["euclidean", "seuclidean", "cityblock", "chebychev", "mahalanobis", "minkowski", "cosine", "correlation", "spearman", "hamming", "jaccard"];
const KD: [&str; 4] = ["euclidean", "cityblock", "chebychev", "minkowski"];

fn build(kd: bool, a: &[Value], who: &str) -> Ret {
    let Some(Value::Num(x)) = a.first() else {
        return rt(if who == "createns" { "At least one input argument is required." } else { "At least one input is required" });
    };
    let rest = &a[1..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let mut kd = kd;
    let mut nsm: Option<String> = None;
    let (mut dist, mut pexp, mut cov, mut scale, mut bucket): (Option<String>, f64, Option<Matrix>, Option<Matrix>, Option<f64>) = (None, 2.0, None, None, None);
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        match key.to_lowercase().as_str() {
            "nsmethod" if who == "createns" => nsm = Some(pair[1].text().unwrap_or_default()),
            "distance" => dist = Some(pair[1].text().unwrap_or_default()),
            "p" => pexp = pair[1].scalar()?,
            "cov" => cov = Some(pair[1].as_matrix()?.clone()),
            "scale" => scale = Some(pair[1].as_matrix()?.clone()),
            "bucketsize" if kd || who == "createns" => bucket = Some(pair[1].scalar().unwrap_or(f64::NAN)),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    if who == "createns" {
        if let Some(m) = &nsm {
            match m.to_lowercase().as_str() {
                "kdtree" => kd = true,
                "exhaustive" => kd = false,
                _ => return rt(format!("Unknown 'NSMETHOD' parameter value: {m}.")),
            }
        } else {
            kd = match &dist {
                None => true,
                Some(d) => KD.iter().any(|n| n.starts_with(&d.to_lowercase()) && !d.is_empty()),
            } && x.dims[1] <= 10;
        }
    }
    let dname = match &dist {
        None => "euclidean".to_string(),
        Some(d) => {
            let l = d.to_lowercase();
            let pool: &[&str] = if kd { &KD } else { &ALL };
            match pool.iter().find(|n| **n == l) {
                Some(n) => n.to_string(),
                None => {
                    if ALL.iter().any(|n| **n == l) {
                        return rt(format!("'{d}' is not a valid value for the Distance argument.  Valid values are: 'euclidean', 'cityblock', 'chebychev', 'minkowski'."));
                    }
                    return rt(format!("'{d}' is not a valid value for the Distance argument."));
                }
            }
        }
    };
    if dname == "minkowski" && !(pexp > 0.0) {
        return rt("The exponent for Minkowski distance must be a positive scalar.");
    }
    if kd {
        if let Some(b) = bucket {
            if !(b >= 1.0 && b.fract() == 0.0) {
                return rt("The BucketSize argument must be a positive integer.");
            }
        }
    }
    let rows = rows_of(x);
    let p = x.dims[1];
    let dpar = match dname.as_str() {
        "minkowski" => Value::num(pexp),
        "mahalanobis" => match cov {
            Some(c) => Value::Num(c),
            None => {
                let ok: Vec<&Vec<f64>> = rows.iter().filter(|r| r.iter().all(|v| !v.is_nan())).collect();
                let n = ok.len() as f64;
                let mean: Vec<f64> = (0..p).map(|j| ok.iter().map(|r| r[j]).sum::<f64>() / n).collect();
                Value::Num(Matrix::from_dims(vec![p, p], (0..p * p).map(|idx| { let (i, j) = (idx % p, idx / p); ok.iter().map(|r| (r[i] - mean[i]) * (r[j] - mean[j])).sum::<f64>() / (n - 1.0) }).collect()))
            }
        },
        "seuclidean" => match scale {
            Some(s) => Value::Num(s),
            None => {
                let sd = cols_stats(&rows, p).1;
                Value::Num(Matrix::from_dims(vec![1, p], sd))
            }
        },
        _ => Value::empty(),
    };
    let (class, props): (&'static str, &'static [&'static str]) = if kd { ("KDTreeSearcher", KD_PROPS) } else { ("ExhaustiveSearcher", EXH_PROPS) };
    let vals: Vec<Value> = props
        .iter()
        .map(|pn| match *pn {
            "X" => Value::Num(x.clone()),
            "Distance" => Value::chars(&dname),
            "DistParameter" => dpar.clone(),
            _ => Value::num(bucket.unwrap_or(50.0)),
        })
        .collect();
    Ok(Some(vec![Value::Object(Object { class: class_def(class, props), props: Rc::new(RefCell::new(vals)) })]))
}

pub fn call(name: &str, a: &[Value]) -> Ret {
    match name {
        "createns" => build(false, a, "createns"),
        "ExhaustiveSearcher" => build(false, a, "ExhaustiveSearcher"),
        _ => build(true, a, "KDTreeSearcher"),
    }
}

/// The argument list of `knnsearch(obj,Y,...)` / `rangesearch(obj,Y,r,...)` rewritten as the
/// function form `knnsearch(X,Y,...)` with the object's distance options first.
pub fn expand_args(who: &str, a: &[Value]) -> Option<Vec<Value>> {
    let Some(Value::Object(o)) = a.first() else { return None };
    if o.class.name != "ExhaustiveSearcher" && o.class.name != "KDTreeSearcher" {
        return None;
    }
    let get = |n: &str| -> Value {
        let idx = o.class.props.iter().position(|p| p.name == n);
        idx.map(|i| o.props.borrow()[i].clone()).unwrap_or_else(Value::empty)
    };
    let dist = get("Distance");
    let par = get("DistParameter");
    let mut out = vec![get("X")];
    let fixed = if who == "knnsearch" { 2 } else { 3 };
    if a.len() < fixed {
        return Some(vec![get("X")]);
    }
    for i in 1..a.len().min(fixed) {
        out.push(a[i].clone());
    }
    out.push(Value::chars("Distance"));
    out.push(dist.clone());
    if let (Ok(d), Value::Num(pm)) = (dist.text(), &par) {
        if !pm.is_empty() {
            let key = match d.as_str() {
                "minkowski" => "P",
                "mahalanobis" => "Cov",
                "seuclidean" => "Scale",
                _ => "",
            };
            if !key.is_empty() {
                out.push(Value::chars(key));
                out.push(par.clone());
            }
        }
    }
    out.extend(a.iter().skip(fixed).cloned());
    Some(out)
}

#[allow(dead_code)]
fn _unused(_: &Matrix) -> bool {
    is_single(&Matrix::scalar(0.0))
}
