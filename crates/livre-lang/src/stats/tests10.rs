//! `bootstrp` and `candgen`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["bootstrp", "candgen"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "bootstrp" => bootstrp(it, a, nargout),
        "candgen" => candgen(it, a, nargout),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn eval_fun(it: &mut Interpreter, f: &Value, args: Vec<Value>) -> Result<Vec<Value>, crate::error::LangError> {
    if f.is_text() {
        let name = f.text().unwrap_or_default();
        it.call_function(&name, args, 1)
    } else if let Value::Func(fv) = f {
        it.call_value(fv, args, 1)
    } else {
        rt("BOOTFUN must be a function handle or a function name.")
    }
}

/// The rows `idx` (1-based) of a matrix.
fn take_rows(m: &Matrix, idx: &[usize]) -> Matrix {
    let (r, c) = (m.dims[0], m.dims[1..].iter().product::<usize>());
    let mut data = Vec::with_capacity(idx.len() * c);
    for j in 0..c {
        for &i in idx {
            data.push(m.data[j * r + i - 1]);
        }
    }
    let mut out = Matrix::from_dims(vec![idx.len(), c], data);
    out.kind = m.kind;
    out
}

fn bootstrp(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("BOOTSTRP requires at least two arguments.");
    }
    let nboot = match &a[0] {
        Value::Num(m) if m.data.len() == 1 && m.kind != Kind::Char => m.data[0],
        _ => f64::NAN,
    };
    if !(nboot > 0.0) || nboot != nboot.round() {
        return rt("NBOOT must be a positive integer.");
    }
    let nboot = nboot as usize;
    let bootfun = a[1].clone();
    let mut rest: Vec<Value> = a[2..].to_vec();
    // 'weights' and 'options' name/value pairs after the data
    let mut weights_given = false;
    if rest.len() > 1 {
        let is_special: Vec<bool> = rest.iter().map(|v| matches!(v, Value::Num(m) if m.kind == Kind::Char && m.dims[0] == 1)).collect();
        if let Some(first) = is_special.iter().position(|&b| b) {
            let special = rest.split_off(first);
            if special.len() % 2 != 0 {
                return rt("Wrong number of arguments.");
            }
            for pair in special.chunks(2) {
                let key = pair[0].text().unwrap_or_default().to_lowercase();
                if !key.is_empty() && "weights".starts_with(&key) {
                    if !is_empty_val(&pair[1]) {
                        weights_given = true;
                    }
                } else if !key.is_empty() && "options".starts_with(&key) {
                } else {
                    return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default()));
                }
            }
        }
    }
    if weights_given {
        return rt("The 'Weights' option is not implemented.");
    }
    // the data arguments
    if rest.is_empty() {
        return rt("BOOTFUN must have at least one argument.");
    }
    let la = rest.len();
    let mut scalard = vec![false; la];
    let mut n = 1usize;
    for k in 0..la {
        let Value::Num(m) = &rest[k] else { return rt("Invalid data type. First argument must be numeric or logical.") };
        let (mut row, col) = (m.dims[0], m.dims[1..].iter().product::<usize>());
        if row.max(col) == 1 {
            scalard[k] = true;
        }
        if row == 1 && col != 1 {
            row = col;
            let mut cm = m.clone();
            cm.dims = vec![col, 1];
            rest[k] = Value::Num(cm);
        }
        if n > 1 && row > 1 && row != n {
            return rt("Nonscalar arguments to BOOTFUN must have the same number of rows.");
        }
        n = n.max(row);
    }
    if n < 2 {
        return rt("BOOTFUN must have at least one non-scalar argument.");
    }
    let sample_args = |rest: &[Value], idx: &[usize]| -> Vec<Value> {
        (0..la)
            .map(|k| {
                if scalard[k] {
                    rest[k].clone()
                } else if let Value::Num(m) = &rest[k] {
                    Value::Num(take_rows(m, idx))
                } else {
                    rest[k].clone()
                }
            })
            .collect()
    };
    let have_all = nargout >= 2;
    if is_empty_val(&bootfun) {
        let mut sam = Vec::with_capacity(n * nboot);
        for _ in 0..nboot {
            let r = it.call_function("randi", vec![Value::num(n as f64), Value::num(n as f64), Value::num(1.0)], 1)?;
            if let Some(Value::Num(m)) = r.first() {
                sam.extend(m.data.iter().copied());
            }
        }
        let mut out = vec![Value::Num(Matrix::from_dims(vec![nboot, 0], vec![]))];
        if have_all {
            out.push(Value::Num(Matrix::from_dims(vec![n, nboot], sam)));
        }
        return Ok(Some(out));
    }
    // the statistic on the original data: it fixes the width
    let first = match eval_fun(it, &bootfun, rest.clone()) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };
    let width = match first.first() {
        Some(Value::Num(m)) => m.data.len(),
        _ => return rt("BOOTFUN must return numeric values."),
    };
    let mut stat = vec![0.0; nboot * width];
    let mut sam = Vec::with_capacity(n * nboot);
    for b in 0..nboot {
        let r = it.call_function("randi", vec![Value::num(n as f64), Value::num(n as f64), Value::num(1.0)], 1)?;
        let idx: Vec<usize> = match r.first() {
            Some(Value::Num(m)) => m.data.iter().map(|&d| d as usize).collect(),
            _ => return rt("randi returned no value"),
        };
        sam.extend(idx.iter().map(|&i| i as f64));
        let res = eval_fun(it, &bootfun, sample_args(&rest, &idx))?;
        let Some(Value::Num(m)) = res.first() else { return rt("BOOTFUN must return numeric values.") };
        if m.data.len() != width {
            return rt("BOOTFUN must return the same number of values for every sample.");
        }
        for (j, v) in m.data.iter().enumerate() {
            stat[j * nboot + b] = *v;
        }
    }
    let mut out = vec![Value::Num(Matrix::from_dims(vec![nboot, width], stat))];
    if have_all {
        out.push(Value::Num(Matrix::from_dims(vec![n, nboot], sam)));
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// candgen

fn candgen(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    let nfactors = match &a[0] {
        Value::Num(m) if m.data.len() == 1 && m.data[0] >= 1.0 && m.data[0] == m.data[0].round() => m.data[0] as usize,
        Value::Num(m) if m.data.len() == 1 && m.data[0] == m.data[0].round() && m.data[0] < 1.0 => return rt("Requires a vector input."),
        Value::Num(m) if m.data.len() == 1 && m.data[0] != m.data[0].round() => return rt("Size inputs must be integers."),
        _ => return rt("Requires a vector input."),
    };
    let model = a.get(1).cloned().filter(|v| !is_empty_val(v)).unwrap_or_else(|| Value::chars("linear"));
    let rest = if a.len() > 2 { &a[2..] } else { &[][..] };
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let names = ["bounds", "levels", "categorical"];
    let mut bounds_v: Option<Value> = None;
    let mut levels_v: Option<Value> = None;
    let mut categ_v: Option<Value> = None;
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        let v = Some(pair[1].clone());
        match pick {
            Some(0) => bounds_v = v,
            Some(1) => levels_v = v,
            Some(2) => categ_v = v,
            None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let nums = |v: &Value| -> Vec<f64> {
        match v {
            Value::Num(m) => m.data.clone(),
            _ => vec![],
        }
    };
    let categ: Vec<usize> = match &categ_v {
        Some(v) if !is_empty_val(v) => {
            let c = nums(v);
            if !c.iter().all(|&d| d >= 1.0 && d <= nfactors as f64 && d == d.round()) {
                return rt("The 'categorical' value must be a vector of factor numbers.");
            }
            c.iter().map(|&d| d as usize).collect()
        }
        _ => vec![],
    };
    let mut levels: Option<Vec<f64>> = levels_v.as_ref().filter(|v| !is_empty_val(v)).map(nums);
    // bounds: a 2-by-nfactors matrix, or a cell array of vectors of values
    let mut bnd_cells: Option<Vec<Vec<f64>>> = None;
    let mut bnd_mat: Vec<f64> = Vec::new();
    match &bounds_v {
        None => {}
        Some(v) if is_empty_val(v) => {}
        Some(Value::Cell(c)) => {
            if c.data.len() != nfactors {
                return rt("The 'bounds' cell array must have one entry per factor.");
            }
            let cells: Vec<Vec<f64>> = c.data.iter().map(nums).collect();
            let obs: Vec<f64> = cells.iter().map(|v| v.len() as f64).collect();
            if obs.iter().any(|&l| l < 2.0) {
                return rt("The 'bounds' array must have at least 2 levels for each factor.");
            }
            if let Some(l) = &levels {
                if l.len() != obs.len() && l.len() != 1 || l.iter().zip(&obs).any(|(x, y)| x != y) {
                    it.emit_warning("The number of levels requested does not match the number of levels in BOUNDS; using the values in BOUNDS.".into(), "stats:doptargchk:BadLevels".into());
                }
            }
            levels = Some(obs);
            bnd_cells = Some(cells);
        }
        Some(Value::Num(m)) => {
            if m.dims != vec![2, nfactors] {
                return rt("The 'bounds' matrix must have 2 rows, and 1 column per factor.");
            }
            bnd_mat = m.data.clone();
        }
        Some(_) => return rt("The 'bounds' matrix must have 2 rows, and 1 column per factor."),
    }
    if bnd_cells.is_none() && bnd_mat.is_empty() {
        bnd_mat = (0..nfactors).flat_map(|_| [-1.0, 1.0]).collect();
        if !categ.is_empty() {
            let lev: Vec<f64> = match &levels {
                None => vec![2.0],
                Some(l) if l.len() == 1 => l.clone(),
                Some(l) => categ.iter().map(|&c| l.get(c - 1).copied().unwrap_or(2.0)).collect(),
            };
            for (k, &c) in categ.iter().enumerate() {
                bnd_mat[(c - 1) * 2] = 1.0;
                bnd_mat[(c - 1) * 2 + 1] = if lev.len() == 1 { lev[0] } else { lev[k] };
            }
        }
    }
    let levels_given = levels.is_some();
    let mut lev: Vec<f64> = match levels {
        Some(l) => {
            if l.len() == 1 {
                vec![l[0]; nfactors]
            } else {
                l
            }
        }
        None => Vec::new(),
    };
    if levels_given {
        if lev.len() != nfactors {
            return rt("The 'levels' vector must have one element for each factor.");
        }
        if lev.iter().any(|&d| d != d.floor() || d < 2.0) {
            return rt("The 'levels' values must be integers greater than 1.");
        }
    } else {
        match &model {
            Value::Num(m) if m.kind != Kind::Char => {
                let (r, c) = (m.dims[0], m.dims[1]);
                if c != nfactors {
                    return rt("MODEL matrix must have one column per factor.");
                }
                lev = (0..c).map(|j| (1.0 + (0..r).map(|i| m.data[j * r + i]).fold(f64::NEG_INFINITY, f64::max)).max(2.0)).collect();
            }
            v => {
                let t = v.text().unwrap_or_default().to_lowercase();
                let pre = |name: &str| !t.is_empty() && name.starts_with(&t);
                let isquad = pre("quadratic") || pre("purequadratic");
                let islin = pre("linear") || pre("interaction") || pre("additive");
                let l = if isquad { 3.0 } else if islin { 2.0 } else { 5.0 };
                lev = vec![l; nfactors];
            }
        }
        for &c in &categ {
            lev[c - 1] = 2.0;
        }
    }
    let fv = it.call_function("fullfact", vec![Value::Num(Matrix::from_dims(vec![1, nfactors], lev.clone()))], 1)?;
    let Some(Value::Num(ff)) = fv.first() else { return rt("fullfact returned no value") };
    let nr = ff.dims[0];
    let mut x = vec![0.0; nr * nfactors];
    for j in 0..nfactors {
        let u: Vec<f64> = match &bnd_cells {
            Some(c) => c[j].clone(),
            None => {
                let (lo, hi) = (bnd_mat[j * 2], bnd_mat[j * 2 + 1]);
                let l = lev[j] as usize;
                (0..l).map(|k| if k + 1 == l { hi } else { lo + (hi - lo) * k as f64 / (l as f64 - 1.0) }).collect()
            }
        };
        for i in 0..nr {
            x[j * nr + i] = u[ff.data[j * nr + i] as usize - 1];
        }
    }
    let xm = Value::Num(Matrix::from_dims(vec![nr, nfactors], x));
    let mut out = vec![xm.clone()];
    if nargout > 1 {
        let mut args = vec![xm, model];
        if let Some(c) = &categ_v {
            if !is_empty_val(c) {
                args.push(c.clone());
            }
        }
        let fx = it.call_function("x2fx", args, 1)?;
        out.push(fx[0].clone());
    }
    Ok(Some(out))
}
