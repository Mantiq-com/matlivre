//! `anovan`: N-way analysis of variance with dummy variables under sum-to-zero constraints,
//! continuous covariates and sequential, hierarchical or constrained sums of squares.

use super::anova::{grp2idx, group_labels};
use super::dist::fcdf;
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["anovan"];

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

/// A matrix kept as rows plus an explicit column count (a constraint matrix can have no rows).
#[derive(Clone)]
struct Mat {
    rows: Vec<Vec<f64>>,
    cols: usize,
}

fn eye(n: usize) -> Mat {
    Mat { rows: (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect(), cols: n }
}

fn kron(a: &Mat, b: &Mat) -> Mat {
    let mut rows = Vec::with_capacity(a.rows.len() * b.rows.len());
    for ra in &a.rows {
        for rb in &b.rows {
            let mut r = Vec::with_capacity(a.cols * b.cols);
            for &x in ra {
                for &y in rb {
                    r.push(x * y);
                }
            }
            rows.push(r);
        }
    }
    Mat { rows, cols: a.cols * b.cols }
}

fn vstack(a: Mat, b: Mat) -> Mat {
    let cols = a.cols;
    let mut rows = a.rows;
    rows.extend(b.rows);
    Mat { rows, cols }
}

/// `length` of a matrix: the larger dimension, or 0 when it has no element.
fn mat_len(m: &Mat) -> usize {
    if m.rows.is_empty() || m.cols == 0 { 0 } else { m.rows.len().max(m.cols) }
}

/// One grouping variable after the missing observations are dropped.
struct Var {
    cont: bool,
    /// Dummy columns: one indicator column per level, or the values.
    dum: Vec<Vec<f64>>,
    /// Sum-to-zero constraint of the variable.
    constr: Mat,
    nlev: usize,
    names: Vec<String>,
    df: f64,
    mean: f64,
}

struct Term {
    dum: Vec<Vec<f64>>,
    constr: Mat,
    name: String,
    df0: f64,
}

/// `termcross`: the products of every column of `a` (varying fastest) with every column of `b`.
fn termcross(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return b.to_vec();
    }
    if b.is_empty() {
        return a.to_vec();
    }
    let mut out = Vec::with_capacity(a.len() * b.len());
    for cb in b {
        for ca in a {
            out.push(ca.iter().zip(cb).map(|(x, y)| x * y).collect());
        }
    }
    out
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Rank-revealing least squares by column-pivoted Gram-Schmidt: (rank, explained sum of squares,
/// coefficients with the dependent columns set to zero).
fn ls(cols: &[Vec<f64>], y: &[f64]) -> (usize, f64, Vec<f64>) {
    let n = y.len();
    let p = cols.len();
    let mut a: Vec<Vec<f64>> = cols.to_vec();
    let maxnorm = a.iter().map(|c| dot(c, c).sqrt()).fold(0.0, f64::max);
    let tol = n.max(p) as f64 * f64::EPSILON * maxnorm;
    let mut perm: Vec<usize> = (0..p).collect();
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut r = vec![vec![0.0; p]; p];
    let mut rank = 0;
    for k in 0..p {
        let mut best = k;
        let mut bn = -1.0;
        for j in k..p {
            let nn = dot(&a[perm[j]], &a[perm[j]]).sqrt();
            if nn > bn {
                bn = nn;
                best = j;
            }
        }
        if bn <= tol {
            break;
        }
        perm.swap(k, best);
        for row in r.iter_mut().take(k) {
            row.swap(k, best);
        }
        let qk: Vec<f64> = a[perm[k]].iter().map(|v| v / bn).collect();
        r[k][k] = bn;
        for j in k + 1..p {
            let pj = perm[j];
            let c = dot(&qk, &a[pj]);
            for i in 0..n {
                a[pj][i] -= c * qk[i];
            }
            let c2 = dot(&qk, &a[pj]);
            for i in 0..n {
                a[pj][i] -= c2 * qk[i];
            }
            r[k][j] = c + c2;
        }
        q.push(qk);
        rank += 1;
    }
    let z: Vec<f64> = q.iter().map(|qk| dot(qk, y)).collect();
    let ss = dot(&z, &z);
    let mut sol = vec![0.0; rank];
    for k in (0..rank).rev() {
        let mut s = z[k];
        for j in k + 1..rank {
            s -= r[k][j] * sol[j];
        }
        sol[k] = s / r[k][k];
    }
    let mut beta = vec![0.0; p];
    for k in 0..rank {
        beta[perm[k]] = sol[k];
    }
    (rank, ss, beta)
}

/// An orthonormal basis (as columns) of the null space of the constraints.
fn null_basis(rows: &[Vec<f64>], p: usize) -> Vec<Vec<f64>> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for r in rows {
        let rn = dot(r, r).sqrt();
        if rn == 0.0 {
            continue;
        }
        let mut v = r.clone();
        for _ in 0..2 {
            for b in &basis {
                let c = dot(b, &v);
                for i in 0..p {
                    v[i] -= c * b[i];
                }
            }
        }
        let nv = dot(&v, &v).sqrt();
        if nv > 1e-10 * rn {
            basis.push(v.iter().map(|x| x / nv).collect());
        }
    }
    let want = p - basis.len();
    let mut out: Vec<Vec<f64>> = Vec::new();
    let mut all = basis.clone();
    for _ in 0..want {
        let mut best = 0;
        let mut bv = -1.0;
        for i in 0..p {
            let proj: f64 = all.iter().map(|b| b[i] * b[i]).sum();
            let rem = 1.0 - proj;
            if rem > bv {
                bv = rem;
                best = i;
            }
        }
        let mut v = vec![0.0; p];
        v[best] = 1.0;
        for _ in 0..2 {
            for b in &all {
                let c = dot(b, &v);
                for i in 0..p {
                    v[i] -= c * b[i];
                }
            }
        }
        let nv = dot(&v, &v).sqrt();
        let u: Vec<f64> = v.iter().map(|x| x / nv).collect();
        all.push(u.clone());
        out.push(u);
    }
    out
}

/// The terms of an integer model: every term of order up to `k`, by order and then
/// lexicographically.
fn terms_up_to(ng: usize, k: usize) -> Vec<Vec<f64>> {
    fn rec(start: usize, ng: usize, left: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if left == 0 {
            out.push(cur.clone());
            return;
        }
        for i in start..ng {
            cur.push(i);
            rec(i + 1, ng, left - 1, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    for ord in 1..=k.min(ng) {
        let mut combos = Vec::new();
        rec(0, ng, ord, &mut Vec::new(), &mut combos);
        for c in combos {
            let mut row = vec![0.0; ng];
            for i in c {
                row[i] = 1.0;
            }
            out.push(row);
        }
    }
    out
}

fn text_list(v: &Value) -> Option<Vec<String>> {
    match v {
        Value::Num(m) if m.kind == Kind::Char => {
            let (r, c) = (m.dims[0], m.dims[1]);
            Some((0..r).map(|i| (0..c).map(|j| char::from_u32(m.data[j * r + i] as u32).unwrap_or(' ')).collect::<String>().trim_end().to_string()).collect())
        }
        Value::Cell(c) => c.data.iter().map(|e| if e.is_text() { e.text().ok() } else { None }).collect(),
        Value::Str(s) => Some(s.data.iter().map(|d| d.clone().unwrap_or_default()).collect()),
        _ => None,
    }
}

fn nrows_of(v: &Value) -> usize {
    match v {
        Value::Num(m) => m.dims[0],
        Value::Cell(c) => c.dims[0],
        Value::Str(s) => s.dims[0],
        _ => 1,
    }
}

pub fn call(_it: &mut Interpreter, _name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    // Y
    let (y_all, single_y): (Vec<f64>, bool) = match &a[0] {
        Value::Num(m) if m.kind == Kind::Double || m.kind == Kind::Single => {
            if m.im.is_some() {
                return rt("Expected Y to be real.");
            }
            if !(m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1)) {
                return rt("Expected Y to be a vector.");
            }
            (m.data.clone(), m.kind == Kind::Single)
        }
        Value::Num(m) if m.kind == Kind::Char => return rt("Expected Y to be one of these types:\n\ndouble, single"),
        Value::Num(m) => return rt(format!("Expected Y to be one of these types:\n\ndouble, single\n\nInstead its type was {}.", m.kind.name())),
        v => return rt(format!("Expected Y to be one of these types:\n\ndouble, single\n\nInstead its type was {}.", v.class_name())),
    };
    let n_all = y_all.len();
    // GROUP
    let groups: Vec<Value> = match &a[1] {
        Value::Cell(c) => c.data.clone(),
        Value::Num(m) if m.kind != Kind::Char && m.dims.len() == 2 && m.dims[0] == n_all => {
            (0..m.dims[1]).map(|j| Value::Num(Matrix::from_dims(vec![n_all, 1], m.data[j * n_all..(j + 1) * n_all].to_vec()))).collect()
        }
        _ => return rt("GROUP must be a cell array or matrix of grouping variables."),
    };
    let ng = groups.len();
    // options: the name/value form, or the old positional one
    let rest = &a[2..];
    let mut model_v: Option<Value> = None;
    let mut sstype_v: Option<Value> = None;
    let mut varnames_v: Option<Value> = None;
    let mut alpha_v: Option<Value> = None;
    let mut random_v: Option<Value> = None;
    let mut continuous_v: Option<Value> = None;
    let mut nested_v: Option<Value> = None;
    let mut newstyle = true;
    if !rest.is_empty() && rest.len() <= 4 {
        let va1 = &rest[0];
        if !va1.is_text() || is_empty_val(va1) {
            newstyle = false;
        } else if let Ok(t) = va1.text() {
            if ["l", "i", "f", "linear", "interaction", "full"].contains(&t.as_str()) {
                newstyle = false;
            }
        }
    }
    if newstyle {
        if rest.len() % 2 != 0 {
            return rt("Wrong number of arguments.");
        }
        let names = ["model", "sstype", "varnames", "display", "alpha", "random", "continuous", "nested"];
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
                Some(0) => model_v = v,
                Some(1) => sstype_v = v,
                Some(2) => varnames_v = v,
                Some(3) => {}
                Some(4) => alpha_v = v,
                Some(5) => random_v = v,
                Some(6) => continuous_v = v,
                Some(7) => nested_v = v,
                _ if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
                _ => return rt(format!("Invalid parameter name: {key}.")),
            }
        }
    } else {
        let n = rest.len();
        if n >= 1 {
            model_v = Some(rest[0].clone());
        }
        if n >= 2 {
            sstype_v = Some(rest[1].clone());
        }
        if n >= 3 {
            varnames_v = Some(rest[2].clone());
        }
    }
    // argcheck
    let varnames: Vec<String> = match &varnames_v {
        None => (1..=ng).map(|i| format!("X{i}")).collect(),
        Some(v) if is_empty_val(v) => (1..=ng).map(|i| format!("X{i}")).collect(),
        Some(v) => {
            let nr = match v {
                Value::Cell(c) => c.data.len(),
                other => nrows_of(other),
            };
            if nr != ng {
                return rt(format!("VARNAMES must have {ng} rows"));
            }
            match text_list(v) {
                Some(t) => t,
                None => return rt("'VarNames' value must be a character array, string array, or cell array of character vectors."),
            }
        }
    };
    let mut hier = false;
    let sstype: usize = match &sstype_v {
        None => 3,
        Some(v) if is_empty_val(v) => 3,
        Some(v) if v.is_text() => {
            let t = v.text().unwrap_or_default().to_lowercase();
            match t.as_str() {
                "one" => 1,
                "two" => 2,
                "three" => 3,
                "h" | "hierarchical" => {
                    hier = true;
                    2
                }
                _ => return rt("SSTYPE must be 1, 2, or 3."),
            }
        }
        Some(Value::Num(m)) if m.data.len() == 1 && (m.data[0] == 1.0 || m.data[0] == 2.0 || m.data[0] == 3.0) => m.data[0] as usize,
        _ => return rt("SSTYPE must be 1, 2, or 3."),
    };
    let mut continuous_col = false;
    let mut continuous = vec![false; ng];
    let bad_cont = format!("CONTINUOUS must be a vector of indices between 1 and {ng}.");
    if let Some(v) = &continuous_v {
        if !is_empty_val(v) {
            match v {
                Value::Num(m) if m.kind == Kind::Logical && m.data.len() == ng && m.numel() == ng => {
                    continuous_col = true;
                    for i in 0..ng {
                        continuous[i] = m.data[i] != 0.0;
                    }
                }
                Value::Num(m) if m.kind != Kind::Char && m.data.iter().all(|&d| d >= 1.0 && d <= ng as f64 && d == d.floor()) => {
                    for &d in &m.data {
                        continuous[d as usize - 1] = true;
                    }
                }
                _ => return rt(bad_cont),
            }
        }
    }
    if let Some(v) = &nested_v {
        if !is_empty_val(v) {
            let ok = matches!(v, Value::Num(m) if m.dims == vec![ng, ng] && m.data.iter().all(|&d| d == 0.0 || d == 1.0));
            if !ok {
                return rt(format!("The NESTED parameter must be a {ng}-by-{ng} matrix of 0's and 1's."));
            }
            if let Value::Num(m) = v {
                if (0..ng).any(|i| (0..ng).any(|j| m.data[j * ng + i] != 0.0 && m.data[i * ng + j] != 0.0)) {
                    return rt("MODEL is invalid.");
                }
            }
        }
    }
    if let Some(v) = &alpha_v {
        match v {
            Value::Num(m) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] > 0.0 && m.data[0] < 1.0 => {}
            _ => return rt("ALPHA must be a scalar between 0 and 1."),
        }
    }
    let mut randomvar = vec![false; ng];
    if let Some(v) = &random_v {
        match v {
            Value::Num(m) if m.kind == Kind::Logical && m.data.len() == ng => {
                for i in 0..ng {
                    randomvar[i] = m.data[i] != 0.0;
                }
            }
            Value::Num(m) if m.kind != Kind::Char && m.data.iter().all(|&d| d >= 1.0 && d <= ng as f64 && d == d.floor()) => {
                for &d in &m.data {
                    randomvar[d as usize - 1] = true;
                }
            }
            _ => return rt(format!("RANDOMVAR must be a vector of indices between 1 and {ng}.")),
        }
    }
    if (0..ng).any(|j| randomvar[j] && continuous[j]) {
        return rt("You cannot specify a continuous predictor as having a random effect");
    }
    let termlist: Vec<Vec<f64>> = {
        let mv = model_v.clone().filter(|v| !is_empty_val(v));
        match mv {
            None => terms_up_to(ng, 1),
            Some(v) if v.is_text() => {
                let t = v.text().unwrap_or_default();
                match t.as_str() {
                    "linear" | "l" => terms_up_to(ng, 1),
                    "interaction" | "i" => terms_up_to(ng, 2),
                    "full" | "f" => terms_up_to(ng, ng),
                    _ => return rt("MODEL is invalid."),
                }
            }
            Some(Value::Num(m)) if m.dims.len() == 2 && (m.dims[0] > 1 && m.dims[1] > 1 || (m.dims[0] == 1 && m.dims[1] == ng && ng > 1)) => {
                if m.dims[1] != ng {
                    return rt(format!("MODEL matrix must have {ng} columns."));
                }
                let r = m.dims[0];
                (0..r).map(|i| (0..ng).map(|j| m.data[j * r + i]).collect::<Vec<f64>>()).collect()
            }
            Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] >= 1.0 && m.data[0] == m.data[0].floor() => terms_up_to(ng, m.data[0] as usize),
            Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() > 1 && m.data.iter().all(|&d| d >= 1.0 && d == d.floor()) => {
                // a vector of interaction orders: the terms of exactly those orders
                let mut orders: Vec<usize> = m.data.iter().map(|&d| d as usize).collect();
                orders.sort();
                orders.dedup();
                let all = terms_up_to(ng, *orders.last().unwrap());
                all.into_iter().filter(|t| orders.contains(&(t.iter().filter(|&&v| v != 0.0).count()))).collect()
            }
            Some(_) => return rt("MODEL is invalid."),
        }
    };
    let termlist: Vec<Vec<f64>> = termlist.into_iter().filter(|t| t.iter().any(|&d| d != 0.0)).collect();
    if termlist.is_empty() {
        return rt("MODEL is invalid.");
    }
    for t in &termlist {
        for j in 0..ng {
            if !continuous[j] && t[j] > 1.0 {
                return rt("MODEL cannot have squared or higher terms for categorical factors.");
            }
        }
    }
    if randomvar.iter().any(|&r| r) || nested_v.as_ref().map(|v| !is_empty_val(v)).unwrap_or(false) {
        return rt("The 'random' and 'nested' options are not implemented.");
    }
    // removenans
    let mut keep = vec![true; n_all];
    for (i, v) in y_all.iter().enumerate() {
        if v.is_nan() {
            keep[i] = false;
        }
    }
    let mut labels_all: Vec<Option<Vec<super::anova::Label>>> = Vec::with_capacity(ng);
    let mut cont_vals: Vec<Option<Vec<f64>>> = Vec::with_capacity(ng);
    for (j, g) in groups.iter().enumerate() {
        let nr = match g {
            Value::Num(m) if m.dims.len() == 2 && m.dims[0] == 1 => m.dims[1],
            Value::Cell(c) if c.dims.len() == 2 && c.dims[0] == 1 => c.data.len(),
            Value::Str(s) if s.dims.len() == 2 && s.dims[0] == 1 => s.data.len(),
            other => nrows_of(other),
        };
        if nr != n_all {
            return rt(format!("Group variable {} must have {} elements.", j + 1, n_all));
        }
        let is_vec = match g {
            Value::Num(m) if m.kind == Kind::Char => true,
            Value::Num(m) => m.dims.len() == 2 && (m.dims[0] == 1 || m.dims[1] == 1),
            Value::Cell(c) => c.dims.len() == 2 && (c.dims[0] == 1 || c.dims[1] == 1),
            Value::Str(s) => s.dims.len() == 2 && (s.dims[0] == 1 || s.dims[1] == 1),
            _ => true,
        };
        if !is_vec {
            return rt("Grouping variable must be a vector or a character array.");
        }
        if continuous[j] {
            let vals = match g {
                Value::Num(m) if m.kind == Kind::Double || m.kind == Kind::Single || m.kind == Kind::Logical => m.data.clone(),
                _ => return rt(format!("Factor '{}' is specified as continuous but is not a double or single vector.", varnames[j])),
            };
            for (i, v) in vals.iter().enumerate() {
                if v.is_nan() {
                    keep[i] = false;
                }
            }
            labels_all.push(None);
            cont_vals.push(Some(vals));
        } else {
            let labs = group_labels(g)?;
            let (idx, _) = grp2idx(&labs);
            for (i, v) in idx.iter().enumerate() {
                if v.is_nan() {
                    keep[i] = false;
                }
            }
            labels_all.push(Some(labs));
            cont_vals.push(None);
        }
    }
    let y: Vec<f64> = y_all.iter().zip(&keep).filter(|(_, &k)| k).map(|(v, _)| *v).collect();
    let n = y.len();
    let mut vars: Vec<Var> = Vec::with_capacity(ng);
    for j in 0..ng {
        if continuous[j] {
            let vals: Vec<f64> = cont_vals[j].as_ref().unwrap().iter().zip(&keep).filter(|(_, &k)| k).map(|(v, _)| *v).collect();
            let mean = vals.iter().sum::<f64>() / n.max(1) as f64;
            vars.push(Var { cont: true, dum: vec![vals], constr: Mat { rows: vec![], cols: 1 }, nlev: 1, names: vec![String::new()], df: 1.0, mean });
        } else {
            let labs: Vec<super::anova::Label> = labels_all[j].as_ref().unwrap().iter().zip(&keep).filter(|(_, &k)| k).map(|(l, _)| l.clone()).collect();
            let (idx, names) = grp2idx(&labs);
            let l = names.len();
            let dum: Vec<Vec<f64>> = (0..l).map(|k| idx.iter().map(|&i| if i as usize == k + 1 { 1.0 } else { 0.0 }).collect()).collect();
            vars.push(Var { cont: false, dum, constr: Mat { rows: vec![vec![1.0; l]], cols: l }, nlev: l, names, df: l as f64 - 1.0, mean: 0.0 });
        }
    }
    // centered response, in two passes
    let mut mu = y.iter().sum::<f64>() / n.max(1) as f64;
    let yc0: Vec<f64> = y.iter().map(|v| v - mu).collect();
    let mu2 = yc0.iter().sum::<f64>() / n.max(1) as f64;
    mu += mu2;
    let yc: Vec<f64> = yc0.iter().map(|v| v - mu2).collect();
    let sst = dot(&yc, &yc);
    // the terms: dummy variables and constraints, built in sorted term order
    let nterms = termlist.len();
    let mut sindex: Vec<usize> = (0..nterms).collect();
    sindex.sort_by(|&x, &y| termlist[x].partial_cmp(&termlist[y]).unwrap_or(std::cmp::Ordering::Equal));
    let mut terms: Vec<Option<Term>> = (0..nterms).map(|_| None).collect();
    for jj in 0..nterms {
        let sj = sindex[jj];
        let tm0 = &termlist[sj];
        let tvars: Vec<usize> = (0..ng).filter(|&v| tm0[v] != 0.0).collect();
        let mut pwrs: Vec<f64> = tvars.iter().map(|&v| tm0[v]).collect();
        let mut tm: Vec<bool> = tm0.iter().map(|&v| v > 0.0).collect();
        let mut df0 = 1.0;
        let mut tdum: Vec<Vec<f64>> = vec![];
        let mut tconstr = Mat { rows: vec![vec![1.0]], cols: 1 };
        let mut tn = String::new();
        for varidx in 0..tvars.len() {
            let varnum = tvars[varidx];
            let thispwr = pwrs.remove(0);
            tm[varnum] = false;
            df0 *= vars[varnum].df;
            let mut g = vars[varnum].dum.clone();
            let mut thisname = varnames[varnum].clone();
            if thispwr > 1.0 {
                for col in g.iter_mut() {
                    for v in col.iter_mut() {
                        *v = v.powf(thispwr);
                    }
                }
                thisname = format!("{}^{}", thisname, thispwr);
            }
            tdum = termcross(&g, &tdum);
            let vc = &vars[varnum].constr;
            if tn.is_empty() {
                tn = thisname;
                tconstr = vc.clone();
            } else {
                tn = format!("{tn}:{thisname}");
                let ncol = tconstr.cols;
                let first = kron(vc, &eye(ncol));
                let second = kron(&eye(mat_len(vc)), &tconstr);
                tconstr = vstack(first, second);
            }
            if varidx + 1 < tvars.len() && jj > 0 {
                // an earlier term made of the remaining variables
                let true_cols: Vec<usize> = (0..ng).filter(|&c| tm[c]).collect();
                let mut found: Option<usize> = None;
                for q in 0..jj {
                    let pt = &termlist[sindex[q]];
                    let same = true_cols.iter().zip(&pwrs).all(|(&c, &pw)| pt[c] == pw);
                    let outside_zero = (0..ng).filter(|&c| !tm[c]).all(|c| pt[c] == 0.0);
                    if same && outside_zero {
                        found = Some(sindex[q]);
                        break;
                    }
                }
                if let Some(k) = found {
                    let old = terms[k].as_ref().unwrap();
                    tdum = termcross(&old.dum, &tdum);
                    let oc = &old.constr;
                    let first = kron(&tconstr, &eye(oc.cols));
                    let second = kron(&eye(tconstr.cols), oc);
                    tconstr = vstack(first, second);
                    tn = format!("{tn}:{}", old.name);
                    df0 *= old.df0;
                    break;
                }
            }
        }
        terms[sj] = Some(Term { dum: tdum, constr: tconstr, name: tn, df0 });
    }
    let terms: Vec<Term> = terms.into_iter().map(|t| t.unwrap()).collect();
    let ones = vec![1.0; n];
    // fit of a set of terms: (rank, explained sum of squares, full-length coefficients)
    let fit_terms = |sel: &[usize]| -> (usize, f64, Vec<f64>, Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut cols: Vec<Vec<f64>> = vec![ones.clone()];
        let mut offsets = vec![];
        for &t in sel {
            offsets.push(cols.len());
            cols.extend(terms[t].dum.iter().cloned());
        }
        let p = cols.len();
        let mut crows: Vec<Vec<f64>> = vec![];
        for (si, &t) in sel.iter().enumerate() {
            for r in &terms[t].constr.rows {
                let mut row = vec![0.0; p];
                for (c, &v) in r.iter().enumerate() {
                    row[offsets[si] + c] = v;
                }
                crows.push(row);
            }
        }
        let q0 = null_basis(&crows, p);
        let proj: Vec<Vec<f64>> = q0
            .iter()
            .map(|qc| {
                let mut v = vec![0.0; n];
                for (c, col) in cols.iter().enumerate() {
                    if qc[c] != 0.0 {
                        for i in 0..n {
                            v[i] += col[i] * qc[c];
                        }
                    }
                }
                v
            })
            .collect();
        let (rank, ss, gamma) = ls(&proj, &yc);
        let mut beta = vec![0.0; p];
        for (k, qc) in q0.iter().enumerate() {
            for c in 0..p {
                beta[c] += gamma[k] * qc[c];
            }
        }
        (rank, ss, beta, cols, q0, proj)
    };
    let all: Vec<usize> = (0..nterms).collect();
    let (dfx, ssx, beta_full, cols_full, q0_full, proj_full) = fit_terms(&all);
    let sse = sst - ssx;
    let dfe = n as f64 - dfx as f64;
    // which terms are not contained in each other (for the type 2 sums of squares)
    let notcontained = |i: usize, j: usize| -> bool {
        let mut m = false;
        for v in 0..ng {
            if !continuous[v] && termlist[i][v] != 0.0 && termlist[j][v] == 0.0 {
                m = true;
            }
            if continuous[v] {
                let (x, y) = (termlist[i][v], termlist[j][v]);
                if (hier && x > y) || (!hier && x != y) {
                    m = true;
                }
            }
        }
        m || i == j
    };
    let mut ssterm = vec![0.0; nterms];
    let mut dfterm = vec![0.0; nterms];
    for k in 0..nterms {
        let with: Vec<usize> = match sstype {
            1 => (0..=k).collect(),
            2 => (0..nterms).filter(|&j| notcontained(k, j)).collect(),
            _ => all.clone(),
        };
        let without: Vec<usize> = with.iter().copied().filter(|&j| j != k).collect();
        let (dw, sw, _, _, _, _) = fit_terms(&with);
        let (dwo, swo, _, _, _, _) = fit_terms(&without);
        ssterm[k] = (sw - swo).max(0.0);
        dfterm[k] = dw as f64 - dwo as f64;
        if dfterm[k] == 0.0 {
            ssterm[k] = 0.0;
        }
    }
    let msterm: Vec<f64> = (0..nterms).map(|k| ssterm[k] / dfterm[k].max(1.0)).collect();
    let mse = if dfe > 0.0 { sse / dfe.max(1.0) } else { 0.0 };
    let mut fstat = vec![0.0; nterms];
    let mut pval = vec![f64::NAN; nterms];
    let mut sing = vec![0.0; nterms];
    for k in 0..nterms {
        fstat[k] = if mse > 0.0 { msterm[k] / mse } else if msterm[k] > 0.0 { f64::INFINITY } else { f64::NAN };
        if dfe > 0.0 && dfterm[k] > 0.0 {
            pval[k] = fcdf(fstat[k], dfterm[k], dfe, true);
        }
        if dfterm[k] < terms[k].df0 {
            sing[k] = 1.0;
        }
    }
    let cast = |v: f64| if single_y { v as f32 as f64 } else { v };
    let mut outs: Vec<Value> = vec![Value::Num(Matrix::col(pval.clone()))];
    if nargout > 1 {
        let num = |v: f64| Value::Num(finish(Matrix::scalar(v), single_y));
        let heads = ["Source", "Sum Sq.", "d.f.", "Singular?", "Mean Sq.", "F", "Prob>F"];
        let nr = nterms + 3;
        let mut data: Vec<Value> = vec![Value::empty(); nr * 7];
        let mut set = |r: usize, c: usize, v: Value| data[c * nr + r] = v;
        for (c, h) in heads.iter().enumerate() {
            set(0, c, Value::chars(h));
        }
        for k in 0..nterms {
            set(k + 1, 0, Value::chars(&terms[k].name));
            set(k + 1, 1, num(ssterm[k]));
            set(k + 1, 2, num(dfterm[k]));
            set(k + 1, 3, num(sing[k]));
            set(k + 1, 4, num(msterm[k]));
            set(k + 1, 5, num(fstat[k]));
            set(k + 1, 6, num(pval[k]));
        }
        set(nterms + 1, 0, Value::chars("Error"));
        set(nterms + 1, 1, num(sse));
        set(nterms + 1, 2, num(dfe));
        set(nterms + 1, 3, num(0.0));
        set(nterms + 1, 4, num(mse));
        set(nterms + 2, 0, Value::chars("Total"));
        set(nterms + 2, 1, num(sst));
        set(nterms + 2, 2, num(n as f64 - 1.0));
        set(nterms + 2, 3, num(0.0));
        outs.push(Value::Cell(CellArray::from_dims(vec![nr, 7], data)));
    }
    if nargout > 2 {
        let sv = |v: Matrix| Value::Num(finish(v, single_y));
        // coefficients: the constant (with the mean put back), then every cell of every term
        let mut coeffs = vec![mu + beta_full[0]];
        let mut codes: Vec<Vec<f64>> = vec![vec![0.0; ng]];
        let mut names = vec!["Constant".to_string()];
        let mut termcols = vec![1.0];
        let mut pos = 1;
        for (k, t) in termlist.iter().enumerate() {
            let tv: Vec<usize> = (0..ng).filter(|&j| t[j] != 0.0).collect();
            let ncell = terms[k].dum.len();
            let dims: Vec<usize> = tv.iter().map(|&j| vars[j].nlev).collect();
            for ci in 0..ncell {
                coeffs.push(beta_full[pos + ci]);
                // var 1 slowest
                let mut rem = ci;
                let mut stride = ncell;
                let mut code = vec![0.0; ng];
                let mut parts = Vec::new();
                for (d, &j) in tv.iter().enumerate() {
                    stride /= dims[d];
                    let ix = rem / stride;
                    rem %= stride;
                    code[j] = (ix + 1) as f64;
                    if vars[j].cont {
                        parts.push(varnames[j].clone());
                    } else {
                        parts.push(format!("{}={}", varnames[j], vars[j].names[ix]));
                    }
                }
                names.push(parts.join(" : "));
                codes.push(code);
            }
            termcols.push(ncell as f64);
            pos += ncell;
        }
        let ncoef = coeffs.len();
        let mut fitted = vec![0.0; n];
        for (c, col) in cols_full.iter().enumerate() {
            if beta_full[c] != 0.0 {
                for i in 0..n {
                    fitted[i] += beta_full[c] * col[i];
                }
            }
        }
        let resid_valid: Vec<f64> = (0..n).map(|i| yc[i] - fitted[i]).collect();
        let resid_out = if n < n_all {
            let mut r = vec![f64::NAN; n_all];
            let mut k = 0;
            for i in 0..n_all {
                if keep[i] {
                    r[i] = resid_valid[k];
                    k += 1;
                }
            }
            Matrix::from_dims(vec![1, n_all], r)
        } else {
            Matrix::col(resid_valid)
        };
        let mut st = StructArray::scalar();
        st.set("source", Value::chars("anovan"));
        st.set("resid", sv(resid_out));
        st.set("coeffs", sv(Matrix::col(coeffs)));
        let p_all = cols_full.len();
        let k_red = q0_full.len();
        let full_rank = dfx == k_red && k_red > 0;
        if full_rank {
            // R of the reduced design (rows of Rtr are its columns), a basis of the row space of the design
            let mut qs: Vec<Vec<f64>> = Vec::new();
            let mut rmat = vec![vec![0.0; k_red]; k_red];
            for j in 0..k_red {
                let mut v = proj_full[j].clone();
                for _ in 0..2 {
                    for (i, qi) in qs.iter().enumerate() {
                        let c = dot(qi, &v);
                        rmat[i][j] += c;
                        for t in 0..n {
                            v[t] -= c * qi[t];
                        }
                    }
                }
                let nv = dot(&v, &v).sqrt();
                rmat[j][j] = nv;
                qs.push(v.iter().map(|x| x / nv).collect());
            }
            let mut rtr = vec![0.0; k_red * k_red];
            for i in 0..k_red {
                for j in 0..k_red {
                    rtr[j * k_red + i] = rmat[j][i];
                }
            }
            st.set("Rtr", sv(Matrix::from_dims(vec![k_red, k_red], rtr)));
            let mut rowb: Vec<Vec<f64>> = Vec::new();
            for i in 0..n {
                let mut v: Vec<f64> = cols_full.iter().map(|c| c[i]).collect();
                let rn = dot(&v, &v).sqrt();
                for _ in 0..2 {
                    for b in &rowb {
                        let c = dot(b, &v);
                        for t in 0..p_all {
                            v[t] -= c * b[t];
                        }
                    }
                }
                let nv = dot(&v, &v).sqrt();
                if nv > 1e-10 * rn.max(1e-300) && nv > 0.0 {
                    rowb.push(v.iter().map(|x| x / nv).collect());
                }
            }
            let rr = rowb.len();
            st.set("rowbasis", sv(Matrix::from_dims(vec![rr, p_all], (0..p_all).flat_map(|c| rowb.iter().map(move |r| r[c])).collect())));
        }
        st.set("dfe", Value::num(dfe));
        st.set("mse", Value::Num(finish(Matrix::scalar(mse), single_y)));
        if full_rank {
            st.set("nullproject", sv(Matrix::from_dims(vec![p_all, k_red], q0_full.iter().flat_map(|c| c.iter().copied()).collect())));
        }
        let tmat = Matrix::from_dims(vec![nterms, ng], (0..ng).flat_map(|j| termlist.iter().map(move |t| t[j])).collect());
        st.set("terms", Value::Num(tmat.clone()));
        st.set("nlevels", sv(Matrix::col(vars.iter().map(|v| v.nlev as f64).collect())));
        let mut cm = if continuous_col { Matrix::from_dims(vec![ng, 1], vec![0.0; ng]) } else { Matrix::from_dims(vec![1, ng], vec![0.0; ng]) };
        for j in 0..ng {
            cm.data[j] = if continuous[j] { 1.0 } else { 0.0 };
        }
        cm.kind = Kind::Logical;
        st.set("continuous", Value::Num(cm));
        st.set("vmeans", sv(Matrix::col(vars.iter().map(|v| v.mean).collect())));
        st.set("termcols", Value::Num(Matrix::col(termcols)));
        st.set("coeffnames", Value::Cell(CellArray::col(names.iter().map(|s| Value::chars(s)).collect())));
        st.set("vars", Value::Num(Matrix::from_dims(vec![ncoef, ng], (0..ng).flat_map(|j| codes.iter().map(move |c| c[j])).collect::<Vec<f64>>())));
        st.set("varnames", Value::Cell(CellArray::col(varnames.iter().map(|s| Value::chars(s)).collect())));
        let gn: Vec<Value> = vars.iter().map(|v| if v.cont { Value::Cell(CellArray::col(vec![Value::chars("")])) } else { Value::Cell(CellArray::col(v.names.iter().map(|s| Value::chars(s)).collect())) }).collect();
        st.set("grpnames", Value::Cell(CellArray::col(gn)));
        st.set("vnested", Value::empty());
        for f in ["denom", "dfdenom", "msdenom", "varest", "varci", "txtdenom", "txtems", "rtnames"] {
            st.set(f, Value::empty());
        }
        outs.push(Value::Struct(st));
    }
    if nargout > 3 {
        outs.push(Value::Num(Matrix::from_dims(vec![nterms, ng], (0..ng).flat_map(|j| termlist.iter().map(move |t| t[j])).collect())));
    }
    let _ = cast;
    Ok(Some(outs))
}
