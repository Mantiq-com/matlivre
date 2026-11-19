//! Classification learners: `fitcknn` (ClassificationKNN), `fitcnb` (ClassificationNaiveBayes)
//! and `fitcdiscr` (ClassificationDiscriminant) with `predict`, `loss`, `margin`, `edge` and their
//! resubstitution forms.

use super::linmodel::prop;
use super::neighbors::{cols_stats, dist_to_all, make_metric, rows_of, Metric};
use super::tests2::{labels_of, Labels};
use super::util::*;
use crate::builtins::Ret;
use crate::class::{ClassDef, Object};
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub const NAMES: &[&str] = &["fitcknn", "fitcnb", "fitcdiscr", "matlivre_learner_method"];

const KNN_PROPS: &[&str] = &[
    "NumNeighbors", "Distance", "DistParameter", "IncludeTies", "DistanceWeight", "BreakTies", "NSMethod", "Mu", "Sigma", "Y", "X", "RowsUsed", "W", "ModelParameters", "NumObservations", "BinEdges", "HyperparameterOptimizationResults",
    "PredictorNames", "CategoricalPredictors", "ResponseName", "ExpandedPredictorNames", "ClassNames", "Prior", "Cost", "ScoreTransform",
];
const NB_PROPS: &[&str] = &[
    "Y", "X", "RowsUsed", "W", "ModelParameters", "NumObservations", "BinEdges", "HyperparameterOptimizationResults", "PredictorNames", "CategoricalPredictors", "ResponseName", "ExpandedPredictorNames", "ClassNames", "Prior", "Cost", "ScoreTransform",
    "DistributionNames", "Mu", "Sigma", "DistributionParameters", "CategoricalLevels", "Kernel", "Support", "Width",
];
const DA_PROPS: &[&str] = &[
    "XCentered", "Y", "X", "RowsUsed", "W", "ModelParameters", "NumObservations", "BinEdges", "HyperparameterOptimizationResults", "PredictorNames", "CategoricalPredictors", "ResponseName", "ExpandedPredictorNames", "ClassNames", "Prior", "Cost",
    "ScoreTransform", "DiscrimType", "Gamma", "Delta", "Coeffs", "Mu", "Sigma", "BetweenSigma", "LogDetSigma", "MinGamma", "DeltaPredictor",
];
const METHODS: &[&str] = &["predict", "resubPredict", "loss", "resubLoss", "margin", "resubMargin", "edge", "resubEdge"];

fn props_of(class: &str) -> &'static [&'static str] {
    match class {
        "ClassificationKNN" => KNN_PROPS,
        "ClassificationNaiveBayes" => NB_PROPS,
        _ => DA_PROPS,
    }
}

fn class_def(class: &'static str) -> Rc<ClassDef> {
    thread_local! {
        static DEF: RefCell<HashMap<&'static str, Rc<ClassDef>>> = RefCell::new(HashMap::new());
    }
    DEF.with(|d| {
        if let Some(c) = d.borrow().get(class) {
            return c.clone();
        }
        let mut src = format!("classdef {class}\n  properties\n");
        for p in props_of(class) {
            src += &format!("    {p}\n");
        }
        src += "  end\n  methods\n";
        for m in METHODS {
            src += &format!("    function varargout = {m}(obj, varargin)\n      varargout = matlivre_learner_method('{m}', nargout, obj, varargin{{:}});\n    end\n");
        }
        src += "  end\nend\n";
        let ast = crate::class::parse_classdef(&src).expect("built-in class source parses");
        let def = Rc::new(ClassDef::build(ast, &[]));
        d.borrow_mut().insert(class, def.clone());
        def
    })
}

fn col(v: Vec<f64>) -> Value {
    let n = v.len();
    Value::Num(Matrix::from_dims(vec![n, 1], v))
}

fn row(v: Vec<f64>) -> Value {
    let n = v.len();
    Value::Num(Matrix::from_dims(vec![1, n], v))
}

fn mat(rows: usize, cols: usize, d: Vec<f64>) -> Value {
    Value::Num(Matrix::from_dims(vec![rows, cols], d))
}

fn names_cell(v: &[String]) -> Value {
    Value::Cell(CellArray::row(v.iter().map(|s| Value::chars(s)).collect()))
}

struct Train {
    x: Vec<Vec<f64>>,
    cls: Vec<usize>,
    w: Vec<f64>,
    k: usize,
    p: usize,
    class_names: Value,
    y_val: Value,
    prior: Vec<f64>,
    cost: Vec<f64>,
    pred_names: Vec<String>,
    resp: String,
    rows_used: Value,
    x_mat: Value,
}

fn name_matches(key: &str, names: &[&str]) -> Option<usize> {
    names.iter().position(|n| n.eq_ignore_ascii_case(key))
}

/// The labels of `Y` as class keys, with the class list (sorted) in the type of `Y`.
fn class_keys(yv: &Value) -> Res<(Vec<Option<String>>, Vec<String>, Value, Value)> {
    let lab = labels_of(yv)?;
    match lab {
        Labels::Num(vals, kind) => {
            let key = |v: f64| if v.is_nan() { None } else { Some(format!("{v:?}")) };
            let mut u: Vec<f64> = vals.iter().copied().filter(|v| !v.is_nan()).collect();
            u.sort_by(|a, b| a.partial_cmp(b).unwrap());
            u.dedup();
            let mut cn = Matrix::from_dims(vec![u.len(), 1], u.clone());
            cn.kind = kind;
            let mut yc = Matrix::from_dims(vec![vals.len(), 1], vals.clone());
            yc.kind = kind;
            Ok((vals.iter().map(|&v| key(v)).collect(), u.iter().map(|&v| format!("{v:?}")).collect(), Value::Num(cn), Value::Num(yc)))
        }
        Labels::Text(items, _) => {
            let mut u: Vec<String> = items.iter().filter_map(|s| s.clone()).filter(|s| !s.is_empty()).collect();
            u.sort();
            u.dedup();
            let cn = Value::Cell(CellArray::from_dims(vec![u.len(), 1], u.iter().map(|s| Value::chars(s)).collect()));
            let yc = Value::Cell(CellArray::from_dims(vec![items.len(), 1], items.iter().map(|s| Value::chars(s.as_deref().unwrap_or(""))).collect()));
            Ok((items.iter().map(|s| s.clone().filter(|t| !t.is_empty())).collect(), u.clone(), cn, yc))
        }
    }
}

/// Class list given by the `ClassNames` option (keys in the same encoding as `class_keys`).
fn option_classes(v: &Value) -> Res<(Vec<String>, Value)> {
    match v {
        Value::Num(m) if m.kind != Kind::Char => {
            let mut cn = Matrix::from_dims(vec![m.numel(), 1], m.data.clone());
            cn.kind = m.kind;
            Ok((m.data.iter().map(|x| format!("{x:?}")).collect(), Value::Num(cn)))
        }
        Value::Cell(c) => {
            let s: Vec<String> = c.data.iter().filter_map(|e| e.text().ok()).collect();
            Ok((s.clone(), Value::Cell(CellArray::from_dims(vec![s.len(), 1], s.iter().map(|t| Value::chars(t)).collect()))))
        }
        other => rt(format!("ClassNames must be a vector or cell array, not {}.", other.class_name())),
    }
}

/// Shared argument handling of the three `fitc*` functions. `extra` lists the function's own
/// parameter names, whose values are returned untouched.
fn fit_common(a: &[Value], extra: &[&str], adjust: bool) -> Res<(Train, Vec<(usize, Value)>)> {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let Value::Num(xm) = &a[0] else { return rt("X must be a numeric matrix.") };
    if xm.dims.len() != 2 {
        return rt("X must be a numeric matrix.");
    }
    let (n0, p) = (xm.dims[0], xm.dims[1]);
    let rest = &a[2..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let common = ["ClassNames", "Cost", "Prior", "ResponseName", "PredictorNames", "Weights", "ScoreTransform", "CategoricalPredictors", "RowsUsed"];
    let mut opts: HashMap<&str, Value> = HashMap::new();
    let mut own: Vec<(usize, Value)> = Vec::new();
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        if let Some(i) = name_matches(&key, &common) {
            opts.insert(common[i], pair[1].clone());
        } else if let Some(i) = name_matches(&key, extra) {
            own.push((i, pair[1].clone()));
        } else {
            return rt(format!("{key} is not a valid parameter name."));
        }
    }
    let (keys, auto_classes, auto_cn, auto_y) = class_keys(&a[1])?;
    if keys.len() != n0 {
        return rt("X and Y do not have the same number of observations.");
    }
    let (class_keys_v, class_names, y_val) = match opts.get("ClassNames") {
        Some(v) => {
            let (ks, cn) = option_classes(v)?;
            (ks, cn, auto_y.clone())
        }
        None => (auto_classes.clone(), auto_cn.clone(), auto_y.clone()),
    };
    let k = class_keys_v.len();
    let weights: Vec<f64> = match opts.get("Weights") {
        Some(Value::Num(m)) => {
            if m.numel() != n0 {
                return rt("Weights must be a vector with one element per observation.");
            }
            m.data.clone()
        }
        _ => vec![1.0; n0],
    };
    let xrows = rows_of(xm);
    let mut keep: Vec<usize> = Vec::new();
    let mut cls: Vec<usize> = Vec::new();
    for i in 0..n0 {
        let Some(key) = &keys[i] else { continue };
        let Some(ci) = class_keys_v.iter().position(|c| c == key) else { continue };
        if xrows[i].iter().any(|v| v.is_nan()) || weights[i].is_nan() || weights[i] < 0.0 {
            continue;
        }
        keep.push(i);
        cls.push(ci);
    }
    let n = keep.len();
    let wsum: f64 = keep.iter().map(|&i| weights[i]).sum();
    let w0: Vec<f64> = keep.iter().map(|&i| weights[i] / wsum).collect();
    let x: Vec<Vec<f64>> = keep.iter().map(|&i| xrows[i].clone()).collect();
    // prior
    let mut ck = vec![0.0; k];
    for (j, &c) in cls.iter().enumerate() {
        ck[c] += w0[j];
    }
    let prior: Vec<f64> = match opts.get("Prior") {
        None => ck.clone(),
        Some(v) if v.is_text() => {
            let t = v.text()?.to_lowercase();
            if "empirical".starts_with(&t) && !t.is_empty() {
                ck.clone()
            } else if "uniform".starts_with(&t) && !t.is_empty() {
                vec![1.0 / k as f64; k]
            } else {
                return rt("Prior must be 'empirical', 'uniform', a vector or a structure.");
            }
        }
        Some(Value::Num(m)) => {
            if m.numel() != k {
                return rt(format!("Prior probabilities must be a vector of length {k}."));
            }
            let s: f64 = m.data.iter().sum();
            m.data.iter().map(|v| v / s).collect()
        }
        _ => return rt("Prior must be 'empirical', 'uniform', a vector or a structure."),
    };
    let w: Vec<f64> = {
        let mut tot = vec![0.0; k];
        for (j, &c) in cls.iter().enumerate() {
            tot[c] += w0[j];
        }
        w0.iter().enumerate().map(|(j, v)| if adjust && tot[cls[j]] > 0.0 { prior[cls[j]] * v / tot[cls[j]] } else { *v }).collect()
    };
    let cost: Vec<f64> = match opts.get("Cost") {
        None => (0..k * k).map(|i| if i % k == i / k { 0.0 } else { 1.0 }).collect(),
        Some(Value::Num(m)) => {
            if m.dims != [k, k] {
                return rt(format!("Misclassification cost matrix must be {k}-by-{k}."));
            }
            m.data.clone()
        }
        _ => return rt("Cost must be a numeric matrix."),
    };
    let pred_names: Vec<String> = match opts.get("PredictorNames") {
        Some(Value::Cell(c)) => c.data.iter().filter_map(|e| e.text().ok()).collect(),
        _ => (1..=p).map(|i| format!("x{i}")).collect(),
    };
    let resp = opts.get("ResponseName").and_then(|v| v.text().ok()).unwrap_or_else(|| "Y".into());
    let rows_used = if n == n0 {
        Value::empty()
    } else {
        let mut f = vec![0.0; n0];
        for &i in &keep {
            f[i] = 1.0;
        }
        Value::Num(Matrix::logical(n0, 1, f))
    };
    // Y of the rows used, in the type of the input
    let y_used = match &y_val {
        Value::Num(m) => {
            let mut o = Matrix::from_dims(vec![n, 1], keep.iter().map(|&i| m.data[i]).collect());
            o.kind = m.kind;
            Value::Num(o)
        }
        Value::Cell(c) => Value::Cell(CellArray::from_dims(vec![n, 1], keep.iter().map(|&i| c.data[i].clone()).collect())),
        v => v.clone(),
    };
    let mut xd = vec![0.0; n * p];
    for (r, row) in x.iter().enumerate() {
        for j in 0..p {
            xd[j * n + r] = row[j];
        }
    }
    let x_mat = Value::Num(Matrix::from_dims(vec![n, p], xd));
    Ok((Train { x, cls, w, k, p, class_names, y_val: y_used, prior, cost, pred_names, resp, rows_used, x_mat }, own))
}

fn build_object(class: &'static str, t: &Train, extra: Vec<(&str, Value)>) -> Value {
    let mut vals: Vec<Value> = Vec::new();
    let n = t.x.len();
    for pn in props_of(class) {
        let v = extra.iter().find(|(k, _)| k == pn).map(|(_, v)| v.clone()).unwrap_or_else(|| match *pn {
            "Y" => t.y_val.clone(),
            "X" => t.x_mat.clone(),
            "RowsUsed" => t.rows_used.clone(),
            "W" => col(t.w.clone()),
            "NumObservations" => Value::num(n as f64),
            "BinEdges" => Value::Cell(CellArray::from_dims(vec![0, 0], vec![])),
            "PredictorNames" | "ExpandedPredictorNames" => names_cell(&t.pred_names),
            "ResponseName" => Value::chars(&t.resp),
            "ClassNames" => t.class_names.clone(),
            "Prior" => row(t.prior.clone()),
            "Cost" => mat(t.k, t.k, t.cost.clone()),
            "ScoreTransform" => Value::chars("none"),
            "ModelParameters" => Value::Struct(StructArray::scalar()),
            _ => Value::empty(),
        });
        vals.push(v);
    }
    Value::Object(Object { class: class_def(class), props: Rc::new(RefCell::new(vals)) })
}

// ---- fitcknn ------------------------------------------------------------------------------------

fn fitcknn(it: &mut Interpreter, a: &[Value]) -> Ret {
    let extra = ["NumNeighbors", "Distance", "DistanceWeight", "Standardize", "IncludeTies", "BreakTies", "Exponent", "Cov", "Scale", "NSMethod", "BucketSize"];
    let (t, own) = fit_common(a, &extra, true)?;
    let mut o: HashMap<&str, Value> = HashMap::new();
    for (i, v) in own {
        o.insert(extra[i], v);
    }
    let k = match o.get("NumNeighbors") {
        None => 1.0,
        Some(v) => {
            let x = v.scalar().unwrap_or(f64::NAN);
            if !(x >= 1.0 && x.fract() == 0.0) {
                return rt("The NumNeighbors input must be a positive integer specifying the number of neighbors.");
            }
            x
        }
    };
    let dist_names = ["euclidean", "cityblock", "chebychev", "minkowski", "mahalanobis", "seuclidean", "cosine", "correlation", "spearman", "hamming", "jaccard"];
    let dist_v = o.get("Distance").cloned();
    let dname: Option<String> = match &dist_v {
        None => Some("euclidean".into()),
        Some(v) if v.is_text() => {
            let s = v.text()?.to_lowercase();
            match dist_names.iter().find(|d| d.starts_with(&s) && !s.is_empty()) {
                Some(d) => Some(d.to_string()),
                None => return rt(format!("Expected Distance to match one of these values:\n\n{}\n\nThe input, '{}', did not match any of the valid values.", dist_names.iter().map(|d| format!("'{d}'")).collect::<Vec<_>>().join(", "), v.text()?)),
            }
        }
        _ => None,
    };
    let weight = match o.get("DistanceWeight") {
        None => "equal".to_string(),
        Some(v) if v.is_text() => {
            let s = v.text()?.to_lowercase();
            match ["equal", "inverse", "squaredinverse"].iter().find(|d| d.starts_with(&s) && !s.is_empty()) {
                Some(d) => d.to_string(),
                None => return rt("The distance weight function must be a 'equal','inverse','squaredinverse', or a function handle."),
            }
        }
        Some(Value::Func(_)) => "custom".into(),
        _ => return rt("The distance weight function must be a 'equal','inverse','squaredinverse', or a function handle."),
    };
    let bt = match o.get("BreakTies") {
        None => "smallest".to_string(),
        Some(v) => {
            let s = v.text().unwrap_or_default();
            match ["smallest", "nearest", "random"].iter().find(|d| d.starts_with(&s.to_lowercase()) && !s.is_empty()) {
                Some(d) => d.to_string(),
                None => return rt(format!("'{s}' is not a valid value for the BreakTies argument.  Valid values are: 'smallest', 'nearest', 'random'.")),
            }
        }
    };
    let ties = o.get("IncludeTies").map(|v| v.scalar().unwrap_or(0.0) != 0.0).unwrap_or(false);
    let std = o.get("Standardize").map(|v| v.scalar().unwrap_or(0.0) != 0.0).unwrap_or(false);
    let (mu, sigma) = if std { cols_stats(&t.x, t.p) } else { (vec![], vec![]) };
    let dist_par = match dname.as_deref() {
        Some("minkowski") => Value::num(o.get("Exponent").and_then(|v| v.scalar().ok()).unwrap_or(2.0)),
        Some("mahalanobis") => match o.get("Cov") {
            Some(Value::Num(m)) => Value::Num(m.clone()),
            _ => {
                let n = t.x.len() as f64;
                let p = t.p;
                let mean: Vec<f64> = (0..p).map(|j| t.x.iter().map(|r| r[j]).sum::<f64>() / n).collect();
                mat(p, p, (0..p * p).map(|idx| { let (i, j) = (idx % p, idx / p); t.x.iter().map(|r| (r[i] - mean[i]) * (r[j] - mean[j])).sum::<f64>() / (n - 1.0) }).collect())
            }
        },
        Some("seuclidean") => match o.get("Scale") {
            Some(Value::Num(m)) => Value::Num(m.clone()),
            _ => row(cols_stats(&t.x, t.p).1),
        },
        _ => Value::empty(),
    };
    let ns = match o.get("NSMethod") {
        Some(v) if v.is_text() => v.text()?.to_lowercase(),
        _ => {
            if matches!(dname.as_deref(), Some("euclidean" | "cityblock" | "chebychev" | "minkowski")) && t.p <= 10 {
                "kdtree".into()
            } else {
                "exhaustive".into()
            }
        }
    };
    let _ = it;
    let distance_value = match (&dname, &dist_v) {
        (Some(d), _) => Value::chars(d),
        (None, Some(v)) => v.clone(),
        _ => Value::chars("euclidean"),
    };
    let obj = build_object(
        "ClassificationKNN",
        &t,
        vec![
            ("NumNeighbors", Value::num(k)),
            ("Distance", distance_value),
            ("DistParameter", dist_par),
            ("IncludeTies", Value::bool(ties)),
            ("DistanceWeight", Value::chars(&weight)),
            ("BreakTies", Value::chars(&bt)),
            ("NSMethod", Value::chars(&ns)),
            ("Mu", if std { row(mu) } else { Value::empty() }),
            ("Sigma", if std { row(sigma) } else { Value::empty() }),
        ],
    );
    Ok(Some(vec![obj]))
}

// ---- fitcnb ------------------------------------------------------------------------------------

fn fitcnb(a: &[Value]) -> Ret {
    let extra = ["DistributionNames", "Kernel", "Support", "Width"];
    let (t, own) = fit_common(a, &extra, true)?;
    for (i, v) in &own {
        if *i == 0 {
            let names: Vec<String> = match v {
                v if v.is_text() => vec![v.text()?; t.p],
                Value::Cell(c) => c.data.iter().filter_map(|e| e.text().ok()).collect(),
                _ => return rt("'DistributionNames' value must be a character vector or string scalar, or a string array or cell array whose length equals the number of predictors."),
            };
            if names.len() != t.p {
                return rt("'DistributionNames' value must be a character vector or string scalar, or a string array or cell array whose length equals the number of predictors.");
            }
            for nm in &names {
                if !["normal", "kernel", "mvmn", "mn"].iter().any(|d| d == nm) {
                    return rt(format!("Unknown distribution name '{nm}' passed for 'DistributionNames' argument."));
                }
                if nm != "normal" {
                    return rt("Only the normal distribution is supported.");
                }
            }
        }
    }
    let (k, p) = (t.k, t.p);
    let mut params = Vec::with_capacity(k * p);
    for j in 0..p {
        for c in 0..k {
            let idx: Vec<usize> = (0..t.x.len()).filter(|&i| t.cls[i] == c).collect();
            let w1: f64 = idx.iter().map(|&i| t.w[i]).sum();
            let w2: f64 = idx.iter().map(|&i| t.w[i] * t.w[i]).sum();
            let mu = idx.iter().map(|&i| t.w[i] * t.x[i][j]).sum::<f64>() / w1;
            let v = idx.iter().map(|&i| t.w[i] * (t.x[i][j] - mu).powi(2)).sum::<f64>() / (w1 - w2 / w1);
            if !(v > 0.0) || !v.is_finite() {
                return rt(format!("A normal distribution cannot be fit for the combination of class {} and predictor {}. The data has zero variance.", class_label(&t.class_names, c), t.pred_names[j]));
            }
            params.push(col(vec![mu, v.sqrt()]));
        }
    }
    let dp = Value::Cell(CellArray::from_dims(vec![k, p], params));
    let dnames = Value::Cell(CellArray::row((0..p).map(|_| Value::chars("normal")).collect()));
    let empties = || Value::Cell(CellArray::row((0..p).map(|_| Value::empty()).collect()));
    let obj = build_object(
        "ClassificationNaiveBayes",
        &t,
        vec![("DistributionNames", dnames), ("DistributionParameters", dp), ("CategoricalLevels", empties()), ("Kernel", empties()), ("Support", empties()), ("Width", mat(k, p, vec![f64::NAN; k * p]))],
    );
    Ok(Some(vec![obj]))
}

fn class_label(cn: &Value, i: usize) -> String {
    match cn {
        Value::Num(m) => {
            let v = m.data[i];
            if v == v.trunc() { format!("{}", v as i64) } else { format!("{v}") }
        }
        Value::Cell(c) => c.data[i].text().unwrap_or_default(),
        _ => String::new(),
    }
}

// ---- fitcdiscr ----------------------------------------------------------------------------------

fn fitcdiscr(it: &mut Interpreter, a: &[Value]) -> Ret {
    let extra = ["DiscrimType", "Gamma", "Delta", "SaveMemory", "FillCoeffs"];
    let (t, own) = fit_common(a, &extra, false)?;
    let types = ["linear", "quadratic", "diaglinear", "diagquadratic", "pseudolinear", "pseudoquadratic"];
    let mut ty = "linear".to_string();
    let mut gamma = 0.0;
    let mut delta = 0.0;
    for (i, v) in &own {
        match i {
            0 => {
                let s = v.text().unwrap_or_default().to_lowercase();
                match types.iter().find(|d| **d == s) {
                    Some(d) => ty = d.to_string(),
                    None => return rt("You must pass 'DiscrimType' as one of the following: linear quadratic diagLinear diagQuadratic pseudoLinear pseudoQuadratic."),
                }
            }
            1 => {
                gamma = v.scalar().unwrap_or(f64::NAN);
                if !(0.0..=1.0).contains(&gamma) {
                    return rt("You must pass 'gamma' as a numeric scalar between 0 and 1.");
                }
            }
            2 => {
                delta = v.scalar().unwrap_or(f64::NAN);
                if !(delta >= 0.0) {
                    return rt("You must pass 'delta' as a non-negative scalar.");
                }
            }
            _ => {}
        }
    }
    let (k, p, n) = (t.k, t.p, t.x.len());
    // class means (weighted), overall mean by prior
    let mut wk = vec![0.0; k];
    for (i, &c) in t.cls.iter().enumerate() {
        wk[c] += t.w[i];
    }
    let mut mu = vec![vec![0.0; p]; k];
    for (i, &c) in t.cls.iter().enumerate() {
        for j in 0..p {
            mu[c][j] += t.w[i] * t.x[i][j] / wk[c];
        }
    }
    let xcen: Vec<Vec<f64>> = (0..n).map(|i| (0..p).map(|j| t.x[i][j] - mu[t.cls[i]][j]).collect()).collect();
    let quad = ty.contains("quadratic");
    let diag = ty.starts_with("diag");
    if diag {
        gamma = 1.0;
    }
    let pseudo = ty.starts_with("pseudo");
    // within-class covariances
    let cov_of = |idx: &[usize]| -> Vec<f64> {
        let w1: f64 = idx.iter().map(|&i| t.w[i]).sum();
        let w2: f64 = idx.iter().map(|&i| t.w[i] * t.w[i]).sum();
        let mut c = vec![0.0; p * p];
        for &i in idx {
            for a1 in 0..p {
                for b1 in 0..p {
                    c[b1 * p + a1] += t.w[i] * xcen[i][a1] * xcen[i][b1];
                }
            }
        }
        c.iter().map(|v| v / (w1 - w2 / w1)).collect()
    };
    let mut sigmas: Vec<Vec<f64>> = Vec::new();
    if quad {
        for c in 0..k {
            let idx: Vec<usize> = (0..n).filter(|&i| t.cls[i] == c).collect();
            sigmas.push(cov_of(&idx));
        }
    } else {
        // pooled: sum of weighted scatter / (1 - sum prior-weights...)
        let mut c = vec![0.0; p * p];
        for i in 0..n {
            for a1 in 0..p {
                for b1 in 0..p {
                    c[b1 * p + a1] += t.w[i] * xcen[i][a1] * xcen[i][b1];
                }
            }
        }
        let w2: f64 = (0..k).map(|c| { let idx_w: Vec<f64> = (0..n).filter(|&i| t.cls[i] == c).map(|i| t.w[i]).collect(); idx_w.iter().map(|w| w * w).sum::<f64>() / wk[c] }).sum();
        let denom = 1.0 - w2;
        for j in 0..p {
            if c[j * p + j] <= 0.0 && !ty.starts_with("pseudo") {
                return rt(format!("Predictor {} has zero within-class variance. Either exclude this predictor or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.", t.pred_names[j]));
            }
        }
        sigmas.push(c.iter().map(|v| v / denom).collect());
    }
    if diag {
        for s in sigmas.iter_mut() {
            for a1 in 0..p {
                for b1 in 0..p {
                    if a1 != b1 {
                        s[b1 * p + a1] = 0.0;
                    }
                }
            }
        }
    }
    if gamma > 0.0 {
        for s in sigmas.iter_mut() {
            for a1 in 0..p {
                for b1 in 0..p {
                    if a1 != b1 {
                        s[b1 * p + a1] *= 1.0 - gamma;
                    }
                }
            }
        }
    }
    // zero-variance check
    for s in &sigmas {
        for j in 0..p {
            if s[j * p + j] <= 0.0 && !pseudo {
                return rt(format!("Predictor {} has zero within-class variance. Either exclude this predictor or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.", t.pred_names[j]));
            }
        }
    }
    let logdet = |it: &mut Interpreter, s: &Vec<f64>| -> Res<f64> {
        let r = it.call_function("det", vec![mat(p, p, s.clone())], 1)?;
        Ok(r[0].scalar()?.ln())
    };
    let mut logdets = Vec::new();
    for s in &sigmas {
        logdets.push(logdet(it, s)?);
    }
    let wj: Vec<f64> = { let s: f64 = wk.iter().sum(); wk.iter().map(|v| v / s).collect() };
    let mubar: Vec<f64> = (0..p).map(|j| (0..k).map(|c| wj[c] * mu[c][j]).sum()).collect();
    let w2s: f64 = wj.iter().map(|v| v * v).sum();
    let mut between = vec![0.0; p * p];
    for c in 0..k {
        for a1 in 0..p {
            for b1 in 0..p {
                between[b1 * p + a1] += wj[c] * (mu[c][a1] - mubar[a1]) * (mu[c][b1] - mubar[b1]) / (1.0 - w2s);
            }
        }
    }
    let mut xc = vec![0.0; n * p];
    for i in 0..n {
        for j in 0..p {
            xc[j * n + i] = xcen[i][j];
        }
    }
    let sigma_val = if diag && quad {
        let mut d = Vec::new();
        for s in &sigmas {
            d.extend((0..p).map(|j| s[j * p + j]));
        }
        Value::Num(Matrix::from_dims(vec![1, p, k], d))
    } else if diag {
        row((0..p).map(|j| sigmas[0][j * p + j]).collect())
    } else if quad {
        let mut d = Vec::new();
        for s in &sigmas {
            d.extend(s.iter().copied());
        }
        Value::Num(Matrix::from_dims(vec![p, p, k], d))
    } else {
        mat(p, p, sigmas[0].clone())
    };
    let ld_val = if quad { col(logdets.clone()) } else { Value::num(logdets[0]) };
    // linear coefficients of the centred class means, per predictor: max |Sigma^-1 (mu_c - mu)|
    let delta_pred: Vec<f64> = if !quad {
        let sd: Vec<f64> = (0..p).map(|j| sigmas[0][j * p + j].sqrt()).collect();
        let rmat: Vec<f64> = (0..p * p).map(|idx| sigmas[0][idx] / (sd[idx % p] * sd[idx / p])).collect();
        let inv = it.call_function("inv", vec![mat(p, p, rmat)], 1)?;
        let Value::Num(im) = &inv[0] else { return rt("inv failed") };
        (0..p).map(|j| (0..k).map(|c| (0..p).map(|q| im.data[j * p + q] * (mu[c][q] - mubar[q]) / sd[q]).sum::<f64>().abs()).fold(0.0, f64::max)).collect()
    } else {
        vec![0.0; p]
    };
    let mut mu_d = vec![0.0; k * p];
    for c in 0..k {
        for j in 0..p {
            mu_d[j * k + c] = mu[c][j];
        }
    }
    // coefficient structure
    let coeffs = discr_coeffs(it, &t, &ty, quad, &mu, &sigmas, &logdets)?;
    let tname = match ty.as_str() {
        "diaglinear" => "diagLinear",
        "diagquadratic" => "diagQuadratic",
        "pseudolinear" => "pseudoLinear",
        "pseudoquadratic" => "pseudoQuadratic",
        other => other,
    };
    let _ = tname;
    let obj = build_object(
        "ClassificationDiscriminant",
        &t,
        vec![
            ("XCentered", Value::Num(Matrix::from_dims(vec![n, p], xc))),
            ("DiscrimType", Value::chars(tname)),
            ("Gamma", Value::num(gamma)),
            ("Delta", Value::num(delta)),
            ("Coeffs", coeffs),
            ("Mu", mat(k, p, mu_d)),
            ("Sigma", sigma_val),
            ("BetweenSigma", mat(p, p, between)),
            ("LogDetSigma", ld_val),
            ("MinGamma", Value::num(0.0)),
            ("DeltaPredictor", row(delta_pred)),
        ],
    );
    Ok(Some(vec![obj]))
}

fn discr_coeffs(it: &mut Interpreter, t: &Train, ty: &str, quad: bool, mu: &[Vec<f64>], sigmas: &[Vec<f64>], logdets: &[f64]) -> Res<Value> {
    let (k, p) = (t.k, t.p);
    let mut fields: Vec<String> = vec!["DiscrimType".into(), "Const".into(), "Linear".into()];
    if quad {
        fields.push("Quadratic".into());
    }
    fields.push("Class1".into());
    fields.push("Class2".into());
    let mut st = StructArray { dims: vec![k, k], fields: fields.clone(), elems: vec![] };
    let cn = |i: usize| -> Value {
        match &t.class_names {
            Value::Num(m) => {
                let mut o = Matrix::scalar(m.data[i]);
                o.kind = m.kind;
                Value::Num(o)
            }
            Value::Cell(c) => c.data[i].clone(),
            v => v.clone(),
        }
    };
    for idx in 0..k * k {
        let (i, j) = (idx % k, idx / k);
        let mut e = vec![Value::Num(Matrix::from_dims(vec![0, 0], vec![])), Value::empty(), Value::empty()];
        if quad {
            e.push(Value::empty());
        }
        e.push(cn(i));
        e.push(cn(j));
        let _ = (i, j);
        st.elems.push(e);
    }
    let inv = |it: &mut Interpreter, c: &Vec<f64>| -> Res<Vec<f64>> {
        let r = it.call_function("inv", vec![mat(p, p, c.clone())], 1)?;
        let Value::Num(m) = &r[0] else { return rt("inv failed") };
        Ok(m.data.clone())
    };
    let tname = match ty {
        "diaglinear" => "diagLinear",
        "diagquadratic" => "diagQuadratic",
        "pseudolinear" => "pseudoLinear",
        "pseudoquadratic" => "pseudoQuadratic",
        o => o,
    };
    let invs: Vec<Vec<f64>> = sigmas.iter().map(|s| inv(it, s)).collect::<Res<_>>()?;
    for i in 0..k {
        for j in 0..k {
            if i == j {
                continue;
            }
            let (si, sj) = if quad { (&invs[i], &invs[j]) } else { (&invs[0], &invs[0]) };
            let (mi, mj) = (&mu[i], &mu[j]);
            let qf = |s: &Vec<f64>, m: &Vec<f64>| -> f64 { (0..p).map(|a1| (0..p).map(|b1| m[a1] * s[b1 * p + a1] * m[b1]).sum::<f64>()).sum() };
            let lin: Vec<f64> = (0..p).map(|a1| (0..p).map(|b1| si[b1 * p + a1] * mi[b1] - sj[b1 * p + a1] * mj[b1]).sum()).collect();
            let mut cst = -(0.5 * (qf(si, mi) - qf(sj, mj)));
            if quad {
                cst -= 0.5 * (logdets[i] - logdets[j]);
            }
            cst += t.prior[i].ln() - t.prior[j].ln();
            // MATLAB's Const is the negative of the boundary offset: Const + Linear'x + x'Qx = 0
            let mut e = vec![Value::chars(tname), Value::num(cst), col(lin)];
            if quad {
                e.push(mat(p, p, (0..p * p).map(|q| -0.5 * (si[q] - sj[q])).collect()));
            }
            e.push(cn(i));
            e.push(cn(j));
            st.elems[j * k + i] = e;
        }
    }
    Ok(Value::Struct(st))
}

// ---- scoring -------------------------------------------------------------------------------------

struct Model {
    class: String,
    k: usize,
    p: usize,
    prior: Vec<f64>,
    cost: Vec<f64>,
    class_names: Value,
}

fn model_of(o: &Object) -> Res<Model> {
    let cn = prop(o, "ClassNames");
    let k = match &cn {
        Value::Num(m) => m.numel(),
        Value::Cell(c) => c.data.len(),
        _ => 0,
    };
    let Value::Num(x) = prop(o, "X") else { return rt("Invalid classification object.") };
    let prior = match prop(o, "Prior") {
        Value::Num(m) => m.data.clone(),
        _ => vec![],
    };
    let cost = match prop(o, "Cost") {
        Value::Num(m) => m.data.clone(),
        _ => vec![],
    };
    Ok(Model { class: o.class.name.clone(), k, p: x.dims[1], prior, cost, class_names: cn })
}

fn class_index(cn: &Value, y: &Value) -> Res<Vec<usize>> {
    let (keys, _, _, _) = class_keys(y)?;
    let cks: Vec<String> = match cn {
        Value::Num(m) => m.data.iter().map(|v| format!("{v:?}")).collect(),
        Value::Cell(c) => c.data.iter().filter_map(|e| e.text().ok()).collect(),
        _ => vec![],
    };
    let mut out = Vec::with_capacity(keys.len());
    for (i, key) in keys.iter().enumerate() {
        match key.as_ref().and_then(|s| cks.iter().position(|c| c == s)) {
            Some(ix) => out.push(ix),
            None => {
                let shown = match y {
                    Value::Num(m) => {
                        let v = m.data[i];
                        if v == v.trunc() { format!("{}", v as i64) } else { format!("{v}") }
                    }
                    Value::Cell(c) => c.data[i].text().unwrap_or_default(),
                    _ => String::new(),
                };
                return rt(format!("You passed an unknown class '{shown}' of type {}.", if matches!(y, Value::Cell(_)) { "char" } else { "double" }));
            }
        }
    }
    Ok(out)
}

fn knn_scores(it: &mut Interpreter, o: &Object, m: &Model, q: &[Vec<f64>]) -> Res<(Vec<Vec<f64>>, Vec<Vec<usize>>)> {
    let Value::Num(xm) = prop(o, "X") else { return rt("Invalid classification object.") };
    let xs0 = rows_of(&xm);
    let yv = prop(o, "Y");
    let cls = class_index(&m.class_names, &yv)?;
    let w: Vec<f64> = match prop(o, "W") {
        Value::Num(wm) => wm.data.clone(),
        _ => vec![1.0; xs0.len()],
    };
    let k = prop(o, "NumNeighbors").scalar()? as usize;
    let ties = prop(o, "IncludeTies").scalar()? != 0.0;
    let dw = prop(o, "DistanceWeight").text()?;
    let (mu, sg): (Vec<f64>, Vec<f64>) = match (prop(o, "Mu"), prop(o, "Sigma")) {
        (Value::Num(a), Value::Num(b)) if !a.is_empty() => (a.data.clone(), b.data.clone()),
        _ => (vec![], vec![]),
    };
    let std_row = |r: &Vec<f64>| -> Vec<f64> { if mu.is_empty() { r.clone() } else { (0..r.len()).map(|j| (r[j] - mu[j]) / sg[j]).collect() } };
    let xs: Vec<Vec<f64>> = xs0.iter().map(std_row).collect();
    let xmat = Matrix::from_dims(vec![xs.len(), m.p], (0..m.p).flat_map(|j| xs.iter().map(move |r| r[j]).collect::<Vec<_>>()).collect());
    let dist_v = prop(o, "Distance");
    let dpar = prop(o, "DistParameter");
    let (pexp, cov, scale) = match (&dist_v, &dpar) {
        (d, Value::Num(pm)) if d.is_text() && d.text()? == "minkowski" => (pm.data.first().copied().unwrap_or(2.0), None, None),
        (d, Value::Num(pm)) if d.is_text() && d.text()? == "mahalanobis" => (2.0, Some(pm.clone()), None),
        (d, Value::Num(pm)) if d.is_text() && d.text()? == "seuclidean" => (2.0, None, Some(pm.clone())),
        _ => (2.0, None, None),
    };
    // the covariance and scale are in the units of the (standardized) training data
    let metric: Metric = make_metric(it, Some(dist_v.clone()), pexp, cov, scale, &xs, m.p)?;
    let mut scores = Vec::with_capacity(q.len());
    let mut nearest_cls: Vec<Vec<usize>> = Vec::with_capacity(q.len());
    let class_w: Vec<f64> = {
        let mut cw = vec![0.0; m.k];
        for (i, &c) in cls.iter().enumerate() {
            cw[c] += w[i];
        }
        cw
    };
    for row in q {
        if row.iter().any(|v| v.is_nan()) {
            scores.push(vec![f64::NAN; m.k]);
            nearest_cls.push(vec![]);
            continue;
        }
        let qs = std_row(row);
        let d = dist_to_all(it, &metric, &qs, &xs, &xmat)?;
        let mut cand: Vec<(f64, usize)> = d.iter().copied().enumerate().filter(|(_, v)| !v.is_nan()).map(|(i, v)| (v, i)).collect();
        cand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        let mut take = k.min(cand.len());
        if ties && take > 0 {
            let kth = cand[take - 1].0;
            while take < cand.len() && cand[take].0 <= kth {
                take += 1;
            }
        }
        let mut a = vec![0.0; m.k];
        let zero_hit = cand[..take].iter().any(|c| c.0 == 0.0) && (dw == "inverse" || dw == "squaredinverse");
        for &(dd, i) in &cand[..take] {
            let wd = if zero_hit {
                if dd == 0.0 { 1.0 } else { 0.0 }
            } else {
                match dw.as_str() {
                    "inverse" => 1.0 / dd,
                    "squaredinverse" => 1.0 / (dd * dd),
                    _ => 1.0,
                }
            };
            a[cls[i]] += wd * w[i];
        }
        let sc: Vec<f64> = (0..m.k).map(|c| m.prior[c] * a[c] / class_w[c]).collect();
        let tot: f64 = sc.iter().sum();
        scores.push(sc.iter().map(|v| v / tot).collect());
        nearest_cls.push(cand[..take].iter().map(|&(_, i)| cls[i]).collect());
    }
    Ok((scores, nearest_cls))
}

fn nb_scores(o: &Object, m: &Model, q: &[Vec<f64>]) -> Res<Vec<Vec<f64>>> {
    let Value::Cell(dp) = prop(o, "DistributionParameters") else { return rt("Invalid classification object.") };
    let at = |c: usize, j: usize| -> (f64, f64) {
        match &dp.data[j * m.k + c] {
            Value::Num(v) => (v.data[0], v.data[1]),
            _ => (f64::NAN, f64::NAN),
        }
    };
    Ok(q
        .iter()
        .map(|row| {
            let lp: Vec<f64> = (0..m.k)
                .map(|c| {
                    let mut s = m.prior[c].ln();
                    for j in 0..m.p {
                        if row[j].is_nan() {
                            continue;
                        }
                        let (mu, sd) = at(c, j);
                        s += -0.5 * ((row[j] - mu) / sd).powi(2) - sd.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln();
                    }
                    s
                })
                .collect();
            softmax(&lp)
        })
        .collect())
}

fn softmax(lp: &[f64]) -> Vec<f64> {
    let mx = lp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let e: Vec<f64> = lp.iter().map(|v| (v - mx).exp()).collect();
    let s: f64 = e.iter().sum();
    e.iter().map(|v| v / s).collect()
}

fn da_scores(it: &mut Interpreter, o: &Object, m: &Model, q: &[Vec<f64>]) -> Res<Vec<Vec<f64>>> {
    let ty = prop(o, "DiscrimType").text()?.to_lowercase();
    let quad = ty.contains("quadratic");
    let Value::Num(mu) = prop(o, "Mu") else { return rt("Invalid classification object.") };
    let Value::Num(sg) = prop(o, "Sigma") else { return rt("Invalid classification object.") };
    let p = m.p;
    let mut invs: Vec<Vec<f64>> = Vec::new();
    let mut lds: Vec<f64> = Vec::new();
    let nsig = if quad { m.k } else { 1 };
    let diag_form = sg.dims[0] == 1 && p > 1;
    for c in 0..nsig {
        let s = if diag_form {
            let mut d = vec![0.0; p * p];
            for j in 0..p {
                d[j * p + j] = sg.data[c * p + j];
            }
            Matrix::from_dims(vec![p, p], d)
        } else {
            Matrix::from_dims(vec![p, p], sg.data[c * p * p..(c + 1) * p * p].to_vec())
        };
        let r = it.call_function("inv", vec![Value::Num(s.clone())], 1)?;
        let Value::Num(im) = &r[0] else { return rt("inv failed") };
        invs.push(im.data.clone());
        let dt = it.call_function("det", vec![Value::Num(s)], 1)?[0].scalar()?;
        lds.push(dt.ln());
    }
    Ok(q
        .iter()
        .map(|row| {
            let lp: Vec<f64> = (0..m.k)
                .map(|c| {
                    let inv = if quad { &invs[c] } else { &invs[0] };
                    let ld = if quad { lds[c] } else { lds[0] };
                    let z: Vec<f64> = (0..p).map(|j| row[j] - mu.data[j * m.k + c]).collect();
                    let md: f64 = (0..p).map(|a| (0..p).map(|b| z[a] * inv[b * p + a] * z[b]).sum::<f64>()).sum();
                    m.prior[c].ln() - 0.5 * (md + ld)
                })
                .collect();
            softmax(&lp)
        })
        .collect())
}

fn predict_core(it: &mut Interpreter, o: &Object, q: &[Vec<f64>]) -> Res<(Vec<usize>, Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    let m = model_of(o)?;
    let (scores, near) = match m.class.as_str() {
        "ClassificationKNN" => knn_scores(it, o, &m, q)?,
        "ClassificationNaiveBayes" => (nb_scores(o, &m, q)?, vec![]),
        _ => (da_scores(it, o, &m, q)?, vec![]),
    };
    let breaks = if m.class == "ClassificationKNN" { prop(o, "BreakTies").text().unwrap_or_default() } else { "smallest".into() };
    let mut labels = Vec::new();
    let mut costs = Vec::new();
    for (i, s) in scores.iter().enumerate() {
        let c: Vec<f64> = (0..m.k).map(|j| (0..m.k).map(|kk| s[kk] * m.cost[j * m.k + kk]).sum()).collect();
        costs.push(c.clone());
        let mut best = None;
        for (j, &v) in c.iter().enumerate() {
            if v.is_nan() {
                continue;
            }
            match best {
                None => best = Some(j),
                Some(b) if v < c[b] => best = Some(j),
                _ => {}
            }
        }
        let mut lab = best.unwrap_or(0);
        if let (Some(b), "nearest", Some(nn)) = (best, breaks.as_str(), near.get(i)) {
            let tied: Vec<usize> = (0..m.k).filter(|&j| (c[j] - c[b]).abs() <= 1e-12 * c[b].abs().max(1e-300)).collect();
            if tied.len() > 1 {
                if let Some(&first) = nn.iter().find(|cl| tied.contains(cl)) {
                    lab = first;
                }
            }
        }
        labels.push(lab);
    }
    Ok((labels, scores, costs))
}

fn label_value(cn: &Value, idx: &[usize]) -> Value {
    match cn {
        Value::Num(m) => {
            let mut o = Matrix::from_dims(vec![idx.len(), 1], idx.iter().map(|&i| m.data[i]).collect());
            o.kind = m.kind;
            Value::Num(o)
        }
        Value::Cell(c) => Value::Cell(CellArray::from_dims(vec![idx.len(), 1], idx.iter().map(|&i| c.data[i].clone()).collect())),
        v => v.clone(),
    }
}

fn x_arg(o: &Object, m: &Model, v: Option<&Value>) -> Res<(Vec<Vec<f64>>, bool)> {
    let _ = o;
    let Some(Value::Num(x)) = v else { return rt("Not enough input arguments.") };
    if x.dims.len() != 2 || x.dims[1] != m.p {
        return match m.class.as_str() {
            "ClassificationKNN" => rt(format!("Y must be a matrix with {} columns.", m.p)),
            "ClassificationNaiveBayes" => rt(format!("You must pass X as a matrix with {} columns.", m.p)),
            _ => rt(format!("X must have {} columns.", m.p)),
        };
    }
    Ok((rows_of(x), is_single(x)))
}

fn loss_value(kind: &str, scores: &[Vec<f64>], costs: &[Vec<f64>], y: &[usize], pred: &[usize], w: &[f64], m: &Model) -> Res<f64> {
    let n = y.len();
    let margins: Vec<f64> = (0..n).map(|i| scores[i][y[i]]).collect();
    let per: Vec<f64> = match kind {
        "classiferror" => (0..n).map(|i| m.cost[pred[i] * m.k + y[i]]).collect(),
        "mincost" => (0..n).map(|i| m.cost[pred[i] * m.k + y[i]]).collect(),
        "binodeviance" => margins.iter().map(|v| (1.0 + (-2.0 * v).exp()).ln()).collect(),
        "exponential" => margins.iter().map(|v| (-v).exp()).collect(),
        "hinge" => margins.iter().map(|v| (1.0 - v).max(0.0)).collect(),
        "logit" => margins.iter().map(|v| (1.0 + (-v).exp()).ln()).collect(),
        "quadratic" => margins.iter().map(|v| (1.0 - v).powi(2)).collect(),
        _ => return rt("You passed an unknown 'LossFun'."),
    };
    let _ = costs;
    Ok((0..n).map(|i| w[i] * per[i]).sum())
}

fn prior_from_cost(it: &mut Interpreter, m: &Model) -> Res<Vec<f64>> {
    let k = m.k;
    let c = |i: usize, j: usize| m.cost[j * k + i];
    let pcost: Vec<f64> = if k == 1 {
        vec![1.0]
    } else if k == 2 {
        let t = c(0, 1) + c(1, 0);
        vec![c(0, 1) / t, c(1, 0) / t]
    } else {
        let nrows = k * (k - 1) / 2;
        let mut a = vec![0.0; nrows * k];
        let mut r = 0;
        for kk in 0..k - 1 {
            for l in kk + 1..k {
                a[kk * nrows + r] = c(l, kk);
                a[l * nrows + r] = -c(kk, l);
                r += 1;
            }
        }
        let am = Value::Num(Matrix::from_dims(vec![nrows, k], a));
        let rank = it.call_function("rank", vec![am.clone()], 1)?[0].scalar()?;
        let v: Vec<f64> = if (rank as usize) < k {
            let nl = it.call_function("null", vec![am], 1)?;
            let Value::Num(nm) = &nl[0] else { return rt("null failed") };
            nm.data[..k].to_vec()
        } else {
            (0..k).map(|i| (0..k).map(|j| c(i, j)).sum()).collect()
        };
        let s: f64 = v.iter().sum();
        v.iter().map(|x| x / s).collect()
    };
    let pr: Vec<f64> = (0..k).map(|i| m.prior[i] * pcost[i]).collect();
    let s: f64 = pr.iter().sum();
    Ok(pr.iter().map(|v| v / s).collect())
}

fn adj_weights(y: &[usize], w: &[f64], prior: &[f64], k: usize) -> Vec<f64> {
    let mut tot = vec![0.0; k];
    for (i, &c) in y.iter().enumerate() {
        tot[c] += w[i];
    }
    let adj: Vec<f64> = y.iter().enumerate().map(|(i, &c)| if tot[c] > 0.0 { prior[c] * w[i] / tot[c] } else { 0.0 }).collect();
    let s: f64 = adj.iter().sum();
    adj.iter().map(|v| v / s).collect()
}

fn method(it: &mut Interpreter, a: &[Value]) -> Ret {
    let name = a[0].text().unwrap_or_default();
    let Value::Object(o) = &a[2] else { return rt("Invalid classification object.") };
    let rest = &a[3..];
    let m = model_of(o)?;
    let resub = name.starts_with("resub");
    let base = name.trim_start_matches("resub").to_lowercase();
    let (q, single, ytrue, wobs, opts): (Vec<Vec<f64>>, bool, Option<Vec<usize>>, Vec<f64>, &[Value]) = if resub {
        let Value::Num(xm) = prop(o, "X") else { return rt("Invalid classification object.") };
        let y = class_index(&m.class_names, &prop(o, "Y"))?;
        let w = match prop(o, "W") {
            Value::Num(wm) => wm.data.clone(),
            _ => vec![],
        };
        (rows_of(&xm), false, Some(y), w, rest)
    } else if base == "predict" {
        let (q, s) = x_arg(o, &m, rest.first())?;
        (q, s, None, vec![], &rest[1.min(rest.len())..])
    } else {
        let (q, s) = x_arg(o, &m, rest.first())?;
        let Some(yv) = rest.get(1) else { return rt("Not enough input arguments.") };
        let (keys, ..) = class_keys(yv)?;
        if keys.len() != q.len() {
            return rt("X and Y must have the same number of observations.");
        }
        let y = class_index(&m.class_names, yv)?;
        (q, s, Some(y), vec![1.0; keys.len()], &rest[2..])
    };
    let (pred, scores, costs) = predict_core(it, o, &q)?;
    let sc = |v: Vec<Vec<f64>>, rows: usize, cols: usize| -> Value {
        let mut d = vec![0.0; rows * cols];
        for (i, r) in v.iter().enumerate() {
            for (j, x) in r.iter().enumerate() {
                d[j * rows + i] = *x;
            }
        }
        Value::Num(finish(Matrix::from_dims(vec![rows, cols], d), single))
    };
    if base == "predict" {
        let n = q.len();
        return Ok(Some(vec![label_value(&m.class_names, &pred), sc(scores, n, m.k), sc(costs, n, m.k)]));
    }
    let y = ytrue.unwrap();
    let mut lossfun = "classiferror".to_string();
    let mut wobs = wobs;
    for pair in opts.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        match (key.as_str(), pair.get(1)) {
            ("lossfun", Some(v)) if v.is_text() => lossfun = v.text()?.to_lowercase(),
            ("weights", Some(Value::Num(wm))) => wobs = wm.data.clone(),
            _ => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
        }
    }
    let pe = m.prior.clone();
    let _ = prior_from_cost;
    let w = adj_weights(&y, &wobs, &pe, m.k);
    let n = y.len();
    match base.as_str() {
        "loss" => {
            let v = loss_value(&lossfun, &scores, &costs, &y, &pred, &w, &m)?;
            Ok(Some(vec![Value::Num(finish(Matrix::scalar(v), single))]))
        }
        _ => {
            let marg: Vec<f64> = (0..n)
                .map(|i| {
                    let other = (0..m.k).filter(|&j| j != y[i]).map(|j| scores[i][j]).fold(f64::NEG_INFINITY, f64::max);
                    scores[i][y[i]] - other
                })
                .collect();
            if base == "margin" {
                Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![n, 1], marg), single))]))
            } else {
                let e: f64 = (0..n).map(|i| w[i] * marg[i]).sum();
                Ok(Some(vec![Value::Num(finish(Matrix::scalar(e), single))]))
            }
        }
    }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    match name {
        "fitcknn" => fitcknn(it, a),
        "fitcnb" => fitcnb(a),
        "fitcdiscr" => fitcdiscr(it, a),
        _ => {
            let outs = method(it, a)?.unwrap_or_default();
            Ok(Some(vec![Value::Cell(CellArray::row(outs))]))
        }
    }
}

/// Command-window text of the learner objects (`ClassificationKNN` ...).
pub fn object_text(o: &Object) -> Option<String> {
    let class = o.class.name.as_str();
    if !matches!(class, "ClassificationKNN" | "ClassificationNaiveBayes" | "ClassificationDiscriminant") {
        return None;
    }
    let mut lines: Vec<(String, String)> = Vec::new();
    let fmt_num_vec = |v: &Value| -> String {
        match v {
            Value::Num(m) if m.numel() == 0 => "[]".into(),
            Value::Num(m) => {
                let items: Vec<String> = m.data.iter().map(|x| if *x == x.trunc() { format!("{}", *x as i64) } else { format!("{x}") }).collect();
                format!("[{}]", items.join(" "))
            }
            Value::Cell(c) => format!("{{{}}}", c.data.iter().map(|e| format!("'{}'", e.text().unwrap_or_default())).collect::<Vec<_>>().join(" ")),
            _ => String::new(),
        }
    };
    lines.push(("ResponseName".into(), format!("'{}'", prop(o, "ResponseName").text().unwrap_or_default())));
    lines.push(("CategoricalPredictors".into(), fmt_num_vec(&prop(o, "CategoricalPredictors"))));
    lines.push(("ClassNames".into(), fmt_num_vec(&prop(o, "ClassNames"))));
    lines.push(("ScoreTransform".into(), format!("'{}'", prop(o, "ScoreTransform").text().unwrap_or_default())));
    lines.push(("NumObservations".into(), format!("{}", prop(o, "NumObservations").scalar().unwrap_or(0.0) as i64)));
    match class {
        "ClassificationKNN" => {
            let d = prop(o, "Distance");
            lines.push(("Distance".into(), if d.is_text() { format!("'{}'", d.text().unwrap_or_default()) } else { "@function".into() }));
            lines.push(("NumNeighbors".into(), format!("{}", prop(o, "NumNeighbors").scalar().unwrap_or(1.0) as i64)));
        }
        "ClassificationNaiveBayes" => {
            let dn = prop(o, "DistributionNames");
            lines.push(("DistributionNames".into(), match &dn {
                Value::Cell(c) => format!("{{{}}}", c.data.iter().map(|e| format!("'{}'", e.text().unwrap_or_default())).collect::<Vec<_>>().join("  ")),
                _ => String::new(),
            }));
            if let Value::Cell(c) = prop(o, "DistributionParameters") {
                lines.push(("DistributionParameters".into(), format!("{{{}x{} cell}}", c.dims[0], c.dims[1])));
            }
        }
        _ => {
            lines.push(("DiscrimType".into(), format!("'{}'", prop(o, "DiscrimType").text().unwrap_or_default())));
            if let Value::Num(m) = prop(o, "Mu") {
                lines.push(("Mu".into(), format!("[{}x{} double]", m.dims[0], m.dims[1])));
            }
            if let Value::Struct(s) = prop(o, "Coeffs") {
                lines.push(("Coeffs".into(), format!("[{}x{} struct]", s.dims[0], s.dims[1])));
            }
        }
    }
    let w = lines.iter().map(|(k, _)| k.len()).max().unwrap_or(0);
    let mut s = format!("  {class}\n");
    for (k, v) in &lines {
        s += &format!("    {:>w$}: {}\n", k, v, w = w);
    }
    Some(s)
}
