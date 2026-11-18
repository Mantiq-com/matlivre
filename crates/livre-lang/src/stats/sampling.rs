//! `randsample` and `datasample`, drawing from the global stream the way MATLAB's do.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::random::global_uniform;
use crate::value::*;

pub const NAMES: &[&str] = &["randsample", "datasample"];

/// `randi(n,1,k)`
pub fn randi_vec(it: &mut Interpreter, n: usize, k: usize) -> Vec<usize> {
    global_uniform(it, k).iter().map(|u| ((u * n as f64).floor() as usize + 1).min(n)).collect()
}

/// `histc(rand(1,k), edges)` for the cumulative probabilities of `w`.
pub fn weighted_with_replacement(it: &mut Interpreter, w: &[f64], k: usize) -> Vec<usize> {
    let sum: f64 = w.iter().sum();
    let mut edges = vec![0.0];
    let mut acc = 0.0;
    for &v in w {
        acc += v / sum;
        edges.push(acc.min(1.0));
    }
    *edges.last_mut().unwrap() = 1.0;
    global_uniform(it, k)
        .iter()
        .map(|&u| {
            let mut idx = edges.len() - 1;
            for j in 0..edges.len() - 1 {
                if u >= edges[j] && u < edges[j + 1] {
                    idx = j + 1;
                    break;
                }
            }
            idx
        })
        .collect()
}

fn randperm(it: &mut Interpreter, n: usize, k: Option<usize>) -> Res<Vec<usize>> {
    let mut args = vec![Value::num(n as f64)];
    if let Some(k) = k {
        args.push(Value::num(k as f64));
    }
    let out = it.call_function("randperm", args, 1)?;
    match out.first() {
        Some(Value::Num(m)) => Ok(m.data.iter().map(|&v| v as usize).collect()),
        _ => rt("randperm failed"),
    }
}

/// Indices of `k` draws without replacement, the way `randsample(n,k)` makes them.
pub fn sample_without_replacement(it: &mut Interpreter, n: usize, k: usize) -> Res<Vec<usize>> {
    if 4 * k > n {
        let rp = randperm(it, n, None)?;
        Ok(rp[..k].to_vec())
    } else {
        let mut flags = vec![false; n];
        let mut sumx = 0;
        while sumx < k {
            for i in randi_vec(it, n, k - sumx) {
                flags[i - 1] = true;
            }
            sumx = flags.iter().filter(|&&f| f).count();
        }
        let y: Vec<usize> = (0..n).filter(|&i| flags[i]).map(|i| i + 1).collect();
        let p = randperm(it, k, None)?;
        Ok(p.iter().map(|&q| y[q - 1]).collect())
    }
}

fn randsample(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let Value::Num(first) = &a[0] else { return rt("N must be a positive integer.") };
    let (n, population): (usize, Option<&Matrix>) = if first.numel() == 1 {
        let v = first.data[0];
        if !(v >= 0.0 && v.fract() == 0.0) {
            return rt("N must be a nonnegative integer, or a vector of population values.");
        }
        (v as usize, None)
    } else {
        (first.numel(), Some(first))
    };
    let k = match &a[1] {
        Value::Num(m) if m.numel() == 1 && m.data[0] >= 0.0 && m.data[0].fract() == 0.0 => m.data[0] as usize,
        _ => return rt("K must be a nonnegative integer."),
    };
    let replace = match a.get(2) {
        None => false,
        Some(Value::Num(m)) if m.numel() == 1 => m.data[0] != 0.0,
        Some(v) if v.is_text() => v.text()?.eq_ignore_ascii_case("true"),
        _ => false,
    };
    let w: Option<Vec<f64>> = match a.get(3) {
        Some(Value::Num(m)) if !m.is_empty() => {
            if m.numel() != n {
                return rt(format!("W must have length {n}."));
            }
            let s: f64 = m.data.iter().sum();
            if !(s > 0.0) || m.data.iter().any(|&v| !(v >= 0.0)) {
                return rt("W must be a vector of non-negative values with a positive sum.");
            }
            Some(m.data.clone())
        }
        _ => None,
    };
    let idx: Vec<usize> = if replace {
        if n == 0 {
            if k == 0 {
                vec![]
            } else {
                return rt("The population is empty, so K must be 0.");
            }
        } else if let Some(w) = &w {
            weighted_with_replacement(it, w, k)
        } else {
            randi_vec(it, n, k)
        }
    } else {
        if k > n {
            return rt(format!("K must be less than or equal to N, the population size ({n}), when sampling without replacement."));
        }
        if w.is_some() {
            return rt("Weights are not supported when sampling without replacement.");
        }
        sample_without_replacement(it, n, k)?
    };
    let vals: Vec<f64> = match population {
        Some(m) => idx.iter().map(|&i| m.data[i - 1]).collect(),
        None => idx.iter().map(|&i| i as f64).collect(),
    };
    let len = vals.len();
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![len, 1], vals))]))
}

fn datasample(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let Value::Num(x) = &a[0] else { return rt("DATA must be numeric.") };
    let k = match &a[1] {
        Value::Num(m) if m.numel() == 1 && m.data[0] >= 0.0 && m.data[0].fract() == 0.0 => m.data[0] as usize,
        _ => return rt("K must be a nonnegative integer."),
    };
    let mut rest = &a[2..];
    let mut dim = x.dims.iter().position(|&d| d != 1).map_or(1, |p| p + 1);
    if let Some(Value::Num(d)) = rest.first() {
        if d.kind != Kind::Char {
            dim = d.data[0] as usize;
            rest = &rest[1..];
        }
    }
    let n = x.dims.get(dim - 1).copied().unwrap_or(1);
    let mut replace = true;
    let mut wgts: Option<Vec<f64>> = None;
    if rest.len() % 2 == 1 {
        return rt("Wrong number of arguments.");
    }
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default().to_lowercase();
        if !key.is_empty() && "replace".starts_with(&key) {
            replace = match &pair[1] {
                Value::Num(m) if m.numel() == 1 && m.kind == Kind::Logical => m.data[0] != 0.0,
                _ => return rt("REPLACE must be a logical scalar."),
            };
        } else if !key.is_empty() && "weights".starts_with(&key) {
            if let Value::Num(m) = &pair[1] {
                if !m.is_empty() {
                    if m.numel() != n {
                        return rt("Weights must be a vector with one element for each element of the data along the sampling dimension.");
                    }
                    wgts = Some(m.data.clone());
                }
            }
        } else {
            return rt(format!("Invalid parameter name: {}.", pair[0].text().unwrap_or_default()));
        }
    }
    let idx: Vec<usize> = if replace {
        match &wgts {
            None => randi_vec(it, n, k),
            Some(w) => weighted_with_replacement(it, w, k),
        }
    } else {
        if k > n {
            return rt("K must be less than or equal to the number of elements of DATA along the sampling dimension when sampling without replacement.");
        }
        match &wgts {
            None => randperm(it, n, Some(k))?,
            Some(w) => {
                // successive draws proportional to the remaining weights
                let mut rem = w.clone();
                let mut out = Vec::new();
                for _ in 0..k {
                    let s: f64 = rem.iter().sum();
                    let u = global_uniform(it, 1)[0] * s;
                    let mut acc = 0.0;
                    let mut pick = rem.len() - 1;
                    for (j, &v) in rem.iter().enumerate() {
                        acc += v;
                        if u < acc && v > 0.0 {
                            pick = j;
                            break;
                        }
                    }
                    out.push(pick + 1);
                    rem[pick] = 0.0;
                }
                out
            }
        }
    };
    // y = x(i,:) for dim 1, x(:,i) for dim 2
    let (r, c) = (x.dims[0], x.dims[1]);
    let y = if dim == 1 {
        let mut d = vec![0.0; k * c];
        for j in 0..c {
            for (t, &i) in idx.iter().enumerate() {
                d[j * k + t] = x.data[j * r + i - 1];
            }
        }
        Matrix::from_dims(vec![k, c], d)
    } else if dim == 2 {
        let mut d = vec![0.0; r * k];
        for (t, &i) in idx.iter().enumerate() {
            for j in 0..r {
                d[t * r + j] = x.data[(i - 1) * r + j];
            }
        }
        Matrix::from_dims(vec![r, k], d)
    } else {
        return rt("datasample supports only dimensions 1 and 2.");
    };
    let mut y = y;
    y.kind = x.kind;
    let mut outs = vec![Value::Num(y)];
    if nargout > 1 {
        outs.push(Value::Num(Matrix::from_dims(vec![1, k], idx.iter().map(|&i| i as f64).collect())));
    }
    Ok(Some(outs))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if name == "randsample" { randsample(it, a) } else { datasample(it, a, nargout) }
}
