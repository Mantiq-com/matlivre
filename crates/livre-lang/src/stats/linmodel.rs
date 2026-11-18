//! `fitlm` for matrix input and the `LinearModel` object: coefficient table, fit statistics,
//! residuals, diagnostics, `predict`, `feval`, `coefCI`, `anova` and the display.

use super::dist::{fcdf, tcdf, tinv};
use super::regress::{qr_pivot, solve_upper};
use super::util::*;
use crate::builtins::Ret;
use crate::class::{ClassDef, Object};
use crate::display::NumFormat;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::table::Table;
use crate::value::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub const NAMES: &[&str] = &["fitlm", "matlivre_lm_method"];

const PROPS: &[&str] = &[
    "ModelFitVsNullModel", "Residuals", "Fitted", "Diagnostics", "MSE", "Robust", "RMSE", "SSE", "SST", "SSR", "NumEstimatedCoefficients", "Formula", "LogLikelihood", "DFE",
    "CoefficientCovariance", "CoefficientNames", "NumCoefficients", "Coefficients", "Rsquared", "ModelCriterion", "VariableInfo", "NumVariables", "VariableNames", "NumPredictors",
    "PredictorNames", "ResponseName", "NumObservations", "Steps", "ObservationInfo", "Variables", "ObservationNames",
];
const METHODS: &[&str] = &["predict", "feval", "coefCI", "anova", "disp"];

fn class_def() -> Rc<ClassDef> {
    thread_local! {
        static DEF: RefCell<HashMap<u8, Rc<ClassDef>>> = RefCell::new(HashMap::new());
    }
    DEF.with(|d| {
        if let Some(c) = d.borrow().get(&0) {
            return c.clone();
        }
        let mut src = String::from("classdef LinearModel\n  properties\n");
        for p in PROPS {
            src += &format!("    {p}\n");
        }
        src += "  end\n  methods\n";
        for m in METHODS {
            if *m == "disp" {
                continue;
            }
            src += &format!("    function varargout = {m}(obj, varargin)\n      varargout = matlivre_lm_method('{m}', nargout, obj, varargin{{:}});\n    end\n");
        }
        src += "  end\nend\n";
        let ast = crate::class::parse_classdef(&src).expect("built-in class source parses");
        let def = Rc::new(ClassDef::build(ast, &[]));
        d.borrow_mut().insert(0, def.clone());
        def
    })
}

fn prop(o: &Object, name: &str) -> Value {
    o.class.prop_index(name).map(|i| o.props.borrow()[i].clone()).unwrap_or_else(Value::empty)
}

fn col(v: Vec<f64>) -> Value {
    let n = v.len();
    Value::Num(Matrix::from_dims(vec![n, 1], v))
}

fn cellstr_col(v: &[String]) -> Value {
    Value::Cell(CellArray::col(v.iter().map(|s| Value::chars(s)).collect()))
}

fn table_of(names: &[&str], vars: Vec<Value>, rows: Option<Vec<String>>) -> Res<Value> {
    let mut t = Table::new(names.iter().map(|s| s.to_string()).collect(), vars)?;
    t.row_names = rows;
    Ok(Value::Table(Box::new(t)))
}

/// Terms of a model: exponent rows over the predictors.
fn terms_of(spec: &str, p: usize, intercept: bool) -> Option<Vec<Vec<usize>>> {
    let zero = vec![0usize; p];
    let mut t: Vec<Vec<usize>> = Vec::new();
    let low = spec.to_lowercase();
    let linear = |t: &mut Vec<Vec<usize>>| {
        for i in 0..p {
            let mut r = zero.clone();
            r[i] = 1;
            t.push(r);
        }
    };
    let inter = |t: &mut Vec<Vec<usize>>| {
        for c in 0..p {
            for r in c + 1..p {
                let mut row = zero.clone();
                row[r] = 1;
                row[c] = 1;
                t.push(row);
            }
        }
    };
    let squares = |t: &mut Vec<Vec<usize>>| {
        for i in 0..p {
            let mut r = zero.clone();
            r[i] = 2;
            t.push(r);
        }
    };
    if intercept {
        t.push(zero.clone());
    }
    match low.as_str() {
        "constant" => {
            t.clear();
            t.push(zero.clone());
        }
        "linear" => linear(&mut t),
        "interactions" => {
            linear(&mut t);
            inter(&mut t);
        }
        "purequadratic" => {
            linear(&mut t);
            squares(&mut t);
        }
        "quadratic" => {
            linear(&mut t);
            inter(&mut t);
            squares(&mut t);
        }
        s if s.starts_with("poly") && s.len() > 4 && s[4..].chars().all(|c| c.is_ascii_digit()) => {
            let degs: Vec<usize> = s[4..].chars().map(|c| c as usize - '0' as usize).collect();
            if degs.len() != p {
                return None;
            }
            let maxd = *degs.iter().max().unwrap();
            // all exponent vectors with each exponent within its degree and total degree <= max
            let mut all: Vec<Vec<usize>> = vec![zero.clone()];
            for i in 0..p {
                let mut next = Vec::new();
                for r in &all {
                    for e in 0..=degs[i] {
                        let mut rr = r.clone();
                        rr[i] = e;
                        next.push(rr);
                    }
                }
                all = next;
            }
            all.retain(|r| r.iter().sum::<usize>() <= maxd && r.iter().sum::<usize>() > 0);
            all.sort_by_key(|r| (r.iter().sum::<usize>(), r.iter().filter(|&&e| e > 0).count() as i64 * 0, std::cmp::Reverse(r.clone())));
            for r in all {
                t.push(r);
            }
            // order: linear terms first, then by total degree, interactions last within a degree
            t.sort_by(|a, b| {
                let (sa, sb) = (a.iter().sum::<usize>(), b.iter().sum::<usize>());
                sa.cmp(&sb)
            });
        }
        _ => return None,
    }
    Some(t)
}

fn term_name(r: &[usize], names: &[String]) -> String {
    if r.iter().all(|&e| e == 0) {
        return "(Intercept)".into();
    }
    let parts: Vec<String> = r
        .iter()
        .enumerate()
        .filter(|(_, &e)| e > 0)
        .map(|(i, &e)| if e == 1 { names[i].clone() } else { format!("{}^{}", names[i], e) })
        .collect();
    parts.join(":")
}

fn formula_text(resp: &str, terms: &[Vec<usize>], names: &[String]) -> String {
    let p = names.len();
    let has_main = |i: usize| terms.iter().any(|r| r.iter().enumerate().all(|(j, &e)| if j == i { e == 1 } else { e == 0 }));
    let mut parts: Vec<String> = Vec::new();
    let mut absorbed = vec![false; p];
    for r in terms {
        let idx: Vec<usize> = (0..p).filter(|&i| r[i] > 0).collect();
        if idx.len() == 2 && r[idx[0]] == 1 && r[idx[1]] == 1 && has_main(idx[0]) && has_main(idx[1]) {
            absorbed[idx[0]] = true;
            absorbed[idx[1]] = true;
        }
    }
    for r in terms {
        let idx: Vec<usize> = (0..p).filter(|&i| r[i] > 0).collect();
        if idx.is_empty() {
            parts.push("1".into());
        } else if idx.len() == 1 && r[idx[0]] == 1 {
            if !absorbed[idx[0]] {
                parts.push(names[idx[0]].clone());
            }
        } else if idx.len() == 2 && r[idx[0]] == 1 && r[idx[1]] == 1 && absorbed[idx[0]] && has_main(idx[0]) && has_main(idx[1]) {
            parts.push(format!("{}*{}", names[idx[0]], names[idx[1]]));
        } else {
            parts.push(idx.iter().map(|&i| if r[i] == 1 { names[i].clone() } else { format!("{}^{}", names[i], r[i]) }).collect::<Vec<_>>().join("*"));
        }
    }
    format!("{resp} ~ {}", parts.join(" + "))
}

fn design(x: &Matrix, terms: &[Vec<usize>]) -> Vec<Vec<f64>> {
    let n = x.dims[0];
    (0..n).map(|i| terms.iter().map(|t| t.iter().enumerate().map(|(j, &e)| if e == 0 { 1.0 } else { x.data[j * n + i].powi(e as i32) }).product()).collect()).collect()
}

fn g3(v: f64) -> String {
    // %.3g
    if v == 0.0 {
        return "0".into();
    }
    if v.is_nan() {
        return "NaN".into();
    }
    let exp = v.abs().log10().floor() as i32;
    let s = if !(-5..3).contains(&exp) {
        let m = v / 10f64.powi(exp);
        let mut ms = format!("{:.2}", m);
        if ms.contains('.') {
            ms = ms.trim_end_matches('0').trim_end_matches('.').to_string();
        }
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{ms}e{sign}{:02}", exp.abs())
    } else {
        let decimals = (2 - exp).max(0) as usize;
        let mut ms = format!("{:.*}", decimals, v);
        if ms.contains('.') {
            ms = ms.trim_end_matches('0').trim_end_matches('.').to_string();
        }
        ms
    };
    s
}

/// MATLAB's `<strong>` markup around the headers and row names of a displayed table.
fn strong_markup(body: &str, rows: &[String]) -> String {
    let width = rows.iter().map(|r| r.chars().count()).max().unwrap_or(0);
    let mut out = Vec::new();
    for (i, line) in body.lines().enumerate() {
        if i < 2 {
            let mut s = String::new();
            let mut tok = String::new();
            for ch in line.chars() {
                if ch == ' ' {
                    if !tok.is_empty() {
                        s += &format!("<strong>{tok}</strong>");
                        tok.clear();
                    }
                    s.push(' ');
                } else {
                    tok.push(ch);
                }
            }
            if !tok.is_empty() {
                s += &format!("<strong>{tok}</strong>");
            }
            out.push(s);
        } else if i >= 3 && !line.trim().is_empty() {
            let chars: Vec<char> = line.chars().collect();
            let lead = 4.min(chars.len());
            let end = (lead + width).min(chars.len());
            let field: String = chars[lead..end].iter().collect();
            let rest: String = chars[end..].iter().collect();
            out.push(format!("{}<strong>{field}</strong>{rest}", " ".repeat(lead)));
        } else {
            out.push(line.to_string());
        }
    }
    out.join("\n")
}

pub fn object_text(o: &Object) -> Option<String> {
    if o.class.name != "LinearModel" {
        return None;
    }
    let formula = prop(o, "Formula").text().unwrap_or_default();
    let coef = prop(o, "Coefficients");
    let Value::Table(t) = &coef else { return None };
    let body = strong_markup(&t.body(NumFormat::Short), t.row_names.as_deref().unwrap_or(&[]));
    let n = prop(o, "NumObservations").scalar().ok()? as i64;
    let dfe = prop(o, "DFE").scalar().ok()? as i64;
    let rmse = prop(o, "RMSE").scalar().ok()?;
    let Value::Struct(r2) = prop(o, "Rsquared") else { return None };
    let ord = r2.get("Ordinary")?.scalar().ok()?;
    let adj = r2.get("Adjusted")?.scalar().ok()?;
    let Value::Struct(f) = prop(o, "ModelFitVsNullModel") else { return None };
    let (fs, fp) = (f.get("Fstat")?.scalar().ok()?, f.get("Pvalue")?.scalar().ok()?);
    let mut s = format!("Linear regression model:\n    {formula}\n\nEstimated Coefficients:\n{}\n\n\nNumber of observations: {n}, Error degrees of freedom: {dfe}\nRoot Mean Squared Error: {}", body.trim_end_matches('\n'), g3(rmse));
    if !fs.is_nan() {
        s += &format!("\nR-squared: {},  Adjusted R-Squared: {}", g3(ord), g3(adj));
    }
    if !fs.is_nan() {
        s += &format!("\nF-statistic vs. constant model: {}, p-value = {}", g3(fs), g3(fp));
    }
    Some(s)
}

fn fit(it: &mut Interpreter, a: &[Value]) -> Ret {
    let _ = it;
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() < 2 {
        return rt("Y argument is required unless X is a dataset or table.");
    }
    let (Value::Num(x), Value::Num(y)) = (&a[0], &a[1]) else { return rt("X and Y must be numeric matrices.") };
    let n = x.dims[0];
    let p = x.dims[1];
    if y.numel() != n || !(y.dims.len() == 2 && (y.dims[0] == 1 || y.dims[1] == 1)) {
        return rt("Predictor and response variables must have the same length.");
    }
    let mut spec = "linear".to_string();
    let mut intercept = true;
    let mut rest = &a[2..];
    let mut var_names: Option<Vec<String>> = None;
    if let Some(v) = rest.first() {
        if v.is_text() {
            let t = v.text()?;
            let low = t.to_lowercase();
            let named = ["intercept", "varnames", "responsevar", "predictorvars", "weights", "exclude", "robustopts", "categoricalvars"];
            if !named.iter().any(|nm| nm.starts_with(&low) && !low.is_empty()) {
                spec = t;
                rest = &rest[1..];
            }
        }
    }
    if rest.len() % 2 == 1 {
        return rt(format!("Formula character vector or string scalar is an invalid model formula: {}.", rest[0].text().unwrap_or_default()));
    }
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        match key.as_str() {
            "intercept" => intercept = pair[1].scalar().unwrap_or(1.0) != 0.0,
            "varnames" => {
                if let Value::Cell(c) = &pair[1] {
                    var_names = Some(c.data.iter().map(|v| v.text().unwrap_or_default()).collect());
                }
            }
            _ => return rt(format!("The '{}' option of fitlm is not implemented.", pair[0].text().unwrap_or_default())),
        }
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
    let intercept = terms.iter().any(|r| r.iter().all(|&e| e == 0));
    let cnames: Vec<String> = terms.iter().map(|r| term_name(r, &names)).collect();
    let xd = design(x, &terms);
    let k = terms.len();
    let single = is_single(x) || is_single(y);
    let yv = y.data.clone();
    let (q, r, perm) = qr_pivot(&xd, n, k);
    let kd = n.min(k);
    let rank = (0..kd).filter(|&i| r[i][i].abs() > n.max(k) as f64 * eps_of(r[0][0])).count();
    if rank < k {
        return rt("Rank-deficient design matrices are not implemented by fitlm.");
    }
    let rq: Vec<Vec<f64>> = r.iter().take(k).map(|row| row[..k].to_vec()).collect();
    let qty: Vec<f64> = (0..k).map(|c| (0..n).map(|i| q[i][c] * yv[i]).sum()).collect();
    let bp = solve_upper(&rq, &qty, k);
    let mut beta = vec![0.0; k];
    for (j, &pc) in perm.iter().enumerate() {
        beta[pc] = bp[j];
    }
    // (X'X)^-1 = R^-1 R^-T, unpermuted
    let mut rinv = vec![vec![0.0; k]; k];
    for j in 0..k {
        let mut e = vec![0.0; k];
        e[j] = 1.0;
        let c = solve_upper(&rq, &e, k);
        for i in 0..k {
            rinv[i][j] = c[i];
        }
    }
    let mut xtxi_p = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in 0..k {
            xtxi_p[i][j] = (0..k).map(|m| rinv[i][m] * rinv[j][m]).sum();
        }
    }
    let mut xtxi = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in 0..k {
            xtxi[perm[i]][perm[j]] = xtxi_p[i][j];
        }
    }
    let fitted: Vec<f64> = xd.iter().map(|row| row.iter().zip(&beta).map(|(a, b)| a * b).sum()).collect();
    let res: Vec<f64> = yv.iter().zip(&fitted).map(|(a, b)| a - b).collect();
    let sse: f64 = res.iter().map(|v| v * v).sum();
    let ybar = yv.iter().sum::<f64>() / n as f64;
    let sst: f64 = yv.iter().map(|v| (v - ybar).powi(2)).sum();
    let dfe = n as f64 - k as f64;
    let mse = sse / dfe;
    let rmse = mse.sqrt();
    let ssr: f64 = fitted.iter().map(|v| (v - ybar).powi(2)).sum();
    let covb: Vec<Vec<f64>> = xtxi.iter().map(|r| r.iter().map(|v| v * mse).collect()).collect();
    let se: Vec<f64> = (0..k).map(|i| covb[i][i].sqrt()).collect();
    let tstat: Vec<f64> = beta.iter().zip(&se).map(|(b, s)| b / s).collect();
    let pval: Vec<f64> = tstat.iter().map(|t| 2.0 * tcdf(-t.abs(), dfe)).collect();
    let ord = 1.0 - sse / sst;
    let adj = 1.0 - (n as f64 - 1.0) / dfe * (sse / sst);
    let (fstat, fp) = if intercept && k > 1 {
        let f = (ssr / (k as f64 - 1.0)) / mse;
        (f, fcdf(f, k as f64 - 1.0, dfe, true))
    } else {
        (f64::NAN, f64::NAN)
    };
    let nf = n as f64;
    let ll = -0.5 * nf * (2.0 * std::f64::consts::PI * sse / nf).ln() - 0.5 * nf;
    let kk = k as f64;
    let aic = -2.0 * ll + 2.0 * kk;
    let aicc = aic + 2.0 * kk * (kk + 1.0) / (nf - kk - 1.0);
    let bic = -2.0 * ll + kk * nf.ln();
    let caic = -2.0 * ll + kk * (nf.ln() + 1.0);
    // hat matrix and diagnostics
    let mut hat = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for a2 in 0..k {
                for b2 in 0..k {
                    s += xd[i][a2] * xtxi[a2][b2] * xd[j][b2];
                }
            }
            hat[j * n + i] = s;
        }
    }
    let h: Vec<f64> = (0..n).map(|i| hat[i * n + i]).collect();
    let s2i: Vec<f64> = (0..n).map(|i| (sse - res[i] * res[i] / (1.0 - h[i])) / (dfe - 1.0)).collect();
    let pearson: Vec<f64> = res.iter().map(|r| r / rmse).collect();
    let stud: Vec<f64> = (0..n).map(|i| res[i] / (s2i[i] * (1.0 - h[i])).sqrt()).collect();
    let stdz: Vec<f64> = (0..n).map(|i| res[i] / (rmse * (1.0 - h[i]).sqrt())).collect();
    let cook: Vec<f64> = (0..n).map(|i| res[i] * res[i] / (kk * mse) * h[i] / (1.0 - h[i]).powi(2)).collect();
    let dffits: Vec<f64> = (0..n).map(|i| stud[i] * (h[i] / (1.0 - h[i])).sqrt()).collect();
    let covratio: Vec<f64> = (0..n).map(|i| 1.0 / (((dfe - 1.0 + stud[i] * stud[i]) / dfe).powf(kk) * (1.0 - h[i]))).collect();
    let mut dfbetas = vec![0.0; n * k];
    for i in 0..n {
        for j in 0..k {
            let dfb: f64 = (0..k).map(|m| xtxi[j][m] * xd[i][m]).sum::<f64>() * res[i] / (1.0 - h[i]);
            dfbetas[j * n + i] = dfb / (s2i[i].sqrt() * xtxi[j][j].sqrt());
        }
    }
    let fin = |v: Vec<f64>| -> Value { col(v) };
    let _ = fin;
    let rows = Some(cnames.clone());
    let coef_t = table_of(&["Estimate", "SE", "tStat", "pValue"], vec![col(beta.clone()), col(se.clone()), col(tstat), col(pval)], rows)?;
    let resid_t = table_of(&["Raw", "Pearson", "Studentized", "Standardized"], vec![col(res.clone()), col(pearson), col(stud), col(stdz)], None)?;
    let diag_t = table_of(
        &["Leverage", "CooksDistance", "Dffits", "S2_i", "CovRatio", "Dfbetas", "HatMatrix"],
        vec![col(h), col(cook), col(dffits), col(s2i), col(covratio), Value::Num(Matrix::from_dims(vec![n, k], dfbetas)), Value::Num(Matrix::from_dims(vec![n, n], hat))],
        None,
    )?;
    let mut vars = vec![];
    let mut vnames: Vec<&str> = names.iter().map(|s| s.as_str()).collect();
    for j in 0..p {
        vars.push(col(x.data[j * n..(j + 1) * n].to_vec()));
    }
    vars.push(col(yv.clone()));
    vnames.push(resp.as_str());
    let variables_t = table_of(&vnames, vars, None)?;
    let mut rsq = StructArray::scalar();
    rsq.set("Ordinary", Value::num(ord));
    rsq.set("Adjusted", Value::num(adj));
    let mut crit = StructArray::scalar();
    crit.set("AIC", Value::num(aic));
    crit.set("AICc", Value::num(aicc));
    crit.set("BIC", Value::num(bic));
    crit.set("CAIC", Value::num(caic));
    let mut fvn = StructArray::scalar();
    fvn.set("Fstat", Value::num(fstat));
    fvn.set("Pvalue", Value::num(fp));
    fvn.set("NullModel", Value::chars("constant"));
    let covm = Value::Num(Matrix::from_dims(vec![k, k], (0..k).flat_map(|j| (0..k).map(|i| covb[i][j]).collect::<Vec<_>>()).collect()));
    let formula = formula_text(&resp, &terms, &names);
    let cell_names = Value::Cell(CellArray::row(cnames.iter().map(|s| Value::chars(s)).collect()));
    let _ = cellstr_col;
    let mut vals: Vec<Value> = Vec::new();
    for pn in PROPS {
        let v = match *pn {
            "ModelFitVsNullModel" => Value::Struct(fvn.clone()),
            "Residuals" => resid_t.clone(),
            "Fitted" => col(fitted.clone()),
            "Diagnostics" => diag_t.clone(),
            "MSE" => Value::num(mse),
            "RMSE" => Value::num(rmse),
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
            "Variables" => variables_t.clone(),
            _ => Value::empty(),
        };
        let v = if single { crate::signal::to_single(v) } else { v };
        vals.push(v);
    }
    let o = Object { class: class_def(), props: Rc::new(RefCell::new(vals)) };
    Ok(Some(vec![Value::Object(o)]))
}

fn design_row_values(o: &Object, xnew: &Matrix) -> Res<(Vec<Vec<f64>>, Vec<usize>)> {
    let p = prop(o, "NumPredictors").scalar().unwrap_or(0.0) as usize;
    if xnew.dims[1] != p {
        return rt(format!("X must have {p} columns."));
    }
    let names: Vec<String> = match prop(o, "CoefficientNames") {
        Value::Cell(c) => c.data.iter().map(|v| v.text().unwrap_or_default()).collect(),
        _ => vec![],
    };
    let pnames: Vec<String> = match prop(o, "PredictorNames") {
        Value::Cell(c) => c.data.iter().map(|v| v.text().unwrap_or_default()).collect(),
        _ => vec![],
    };
    let n = xnew.dims[0];
    let rows: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            names
                .iter()
                .map(|nm| {
                    if nm == "(Intercept)" {
                        return 1.0;
                    }
                    nm.split(':')
                        .map(|part| {
                            let (base, e) = match part.split_once('^') {
                                Some((b, e)) => (b, e.parse::<i32>().unwrap_or(1)),
                                None => (part, 1),
                            };
                            let j = pnames.iter().position(|pn| pn == base).unwrap_or(0);
                            xnew.data[j * n + i].powi(e)
                        })
                        .product()
                })
                .collect()
        })
        .collect();
    Ok((rows, vec![]))
}

fn method(it: &mut Interpreter, a: &[Value]) -> Ret {
    let name = a[0].text().unwrap_or_default();
    let nargout = a[1].scalar().unwrap_or(1.0) as usize;
    if name == "feval" && nargout > 1 {
        return rt("Too many output arguments.");
    }
    let Value::Object(o) = &a[2] else { return rt("Invalid LinearModel.") };
    let rest = &a[3..];
    let coef = match prop(o, "Coefficients") {
        Value::Table(t) => t,
        _ => return rt("Invalid LinearModel."),
    };
    let est: Vec<f64> = match coef.var("Estimate") {
        Some(Value::Num(m)) => m.data.clone(),
        _ => vec![],
    };
    let dfe = prop(o, "DFE").scalar().unwrap_or(f64::NAN);
    let mse = prop(o, "MSE").scalar().unwrap_or(f64::NAN);
    let k = est.len();
    let covb: Vec<f64> = match prop(o, "CoefficientCovariance") {
        Value::Num(m) => m.data.clone(),
        _ => vec![],
    };
    let _ = it;
    match name.as_str() {
        "coefCI" => {
            let alpha = match rest.first() {
                Some(Value::Num(m)) if !m.is_empty() => m.data[0],
                _ => 0.05,
            };
            let t = tinv(1.0 - alpha / 2.0, dfe);
            let se: Vec<f64> = (0..k).map(|i| covb[i * k + i].sqrt()).collect();
            let mut d = vec![0.0; 2 * k];
            for i in 0..k {
                d[i] = est[i] - t * se[i];
                d[k + i] = est[i] + t * se[i];
            }
            Ok(Some(vec![Value::Cell(CellArray::row(vec![Value::Num(Matrix::from_dims(vec![k, 2], d))]))]))
        }
        "predict" | "feval" => {
            let (xnew, opts): (Matrix, &[Value]) = if name == "feval" {
                // one argument per predictor
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
            let mut observation = false;
            for pair in opts.chunks(2) {
                let key = pair[0].text().unwrap_or_default().to_lowercase();
                match key.as_str() {
                    "alpha" => alpha = pair.get(1).and_then(|v| v.scalar().ok()).unwrap_or(0.05),
                    "prediction" => observation = pair.get(1).and_then(|v| v.text().ok()).map_or(false, |t| t.to_lowercase().starts_with("obs")),
                    _ => return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default())),
                }
            }
            let (rows, _) = design_row_values(o, &xnew)?;
            let n = rows.len();
            let ypred: Vec<f64> = rows.iter().map(|r| r.iter().zip(&est).map(|(a, b)| a * b).sum()).collect();
            let t = tinv(1.0 - alpha / 2.0, dfe);
            let mut ci = vec![0.0; 2 * n];
            for i in 0..n {
                let mut v = 0.0;
                for a1 in 0..k {
                    for b1 in 0..k {
                        v += rows[i][a1] * covb[b1 * k + a1] * rows[i][b1];
                    }
                }
                if observation {
                    v += mse;
                }
                ci[i] = ypred[i] - t * v.sqrt();
                ci[n + i] = ypred[i] + t * v.sqrt();
            }
            Ok(Some(vec![Value::Cell(CellArray::row(vec![col(ypred), Value::Num(Matrix::from_dims(vec![n, 2], ci))]))]))
        }
        "anova" => rt("anova of a LinearModel is not implemented."),
        _ => rt("Unknown LinearModel method."),
    }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    if name == "matlivre_lm_method" {
        let r = method(it, a)?.unwrap_or_default();
        // The method hands back a cell holding its outputs.
        return Ok(Some(r));
    }
    fit(it, a)
}
