//! `dbscan`, `procrustes`, the hidden-Markov-model functions `hmmdecode`, `hmmviterbi` and
//! `hmmestimate`, `barttest` and `bbdesign`.

use super::neighbors::{dist_to_all, make_metric, rows_of};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["dbscan", "procrustes", "hmmdecode", "hmmviterbi", "hmmestimate", "barttest", "bbdesign"];

fn num<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(_) => rt("Input must be numeric."),
        None => rt("Not enough input arguments."),
    }
}

fn matches_name(given: &str, full: &str) -> bool {
    !given.is_empty() && full.to_lowercase().starts_with(&given.to_lowercase())
}

// ---- dbscan --------------------------------------------------------------------------------------

fn dbscan(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    let x = num(a, 0)?;
    if x.is_empty() {
        return rt("Expected X to be nonempty.");
    }
    let eps = a[1].scalar()?;
    if !(eps >= 0.0) {
        return rt("Expected epsilon to be nonnegative.");
    }
    let minpts = a[2].scalar()?;
    if !(minpts >= 1.0 && minpts.fract() == 0.0) {
        return rt("Expected minpts to be positive.");
    }
    let rest = &a[3..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let (mut dist, mut pexp, mut cov, mut scale): (Option<Value>, f64, Option<Matrix>, Option<Matrix>) = (None, 2.0, None, None);
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        match key.to_lowercase().as_str() {
            "distance" => dist = Some(pair[1].clone()),
            "p" => pexp = pair[1].scalar()?,
            "cov" => cov = Some(pair[1].as_matrix()?.clone()),
            "scale" => scale = Some(pair[1].as_matrix()?.clone()),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let rows = rows_of(x);
    let n = rows.len();
    let p = x.dims[1];
    let metric = make_metric(it, dist, pexp, cov, scale, &rows, p)?;
    let mut nbrs: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        if rows[i].iter().any(|v| v.is_nan()) {
            nbrs.push(vec![]);
            continue;
        }
        let d = dist_to_all(it, &metric, &rows[i], &rows, x)?;
        nbrs.push((0..n).filter(|&j| d[j] <= eps).collect());
    }
    let core: Vec<bool> = nbrs.iter().map(|v| v.len() as f64 >= minpts).collect();
    let mut label = vec![0i64; n]; // 0 = unvisited, -1 = noise, >0 cluster
    let mut cluster = 0i64;
    for i in 0..n {
        if label[i] != 0 {
            continue;
        }
        if !core[i] {
            label[i] = -1;
            continue;
        }
        cluster += 1;
        label[i] = cluster;
        let mut queue: std::collections::VecDeque<usize> = nbrs[i].iter().copied().filter(|&j| j != i).collect();
        while let Some(j) = queue.pop_front() {
            if label[j] == -1 {
                label[j] = cluster;
                continue;
            }
            if label[j] != 0 {
                continue;
            }
            label[j] = cluster;
            if core[j] {
                for &k in &nbrs[j] {
                    if label[k] == 0 || label[k] == -1 {
                        queue.push_back(k);
                    }
                }
            }
        }
    }
    let idx: Vec<f64> = label.iter().map(|&v| v as f64).collect();
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![n, 1], idx)), Value::Num(Matrix::logical(n, 1, core.iter().map(|&c| if c { 1.0 } else { 0.0 }).collect()))]))
}

// ---- procrustes ----------------------------------------------------------------------------------

fn procrustes(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let (x, y) = (num(a, 0)?, num(a, 1)?);
    let single = is_single(x) || is_single(y);
    let rest = &a[2..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let mut scaling = true;
    let mut reflection: Option<bool> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        if matches_name(&key, "scaling") {
            scaling = pair[1].scalar()? != 0.0;
        } else if matches_name(&key, "reflection") {
            if pair[1].is_text() {
                if !pair[1].text()?.eq_ignore_ascii_case("best") {
                    return rt("REFLECTION must be a scalar, or 'best'.");
                }
            } else {
                reflection = Some(pair[1].scalar()? != 0.0);
            }
        } else {
            return rt(format!("Invalid parameter name: {key}."));
        }
    }
    let (n, m) = (x.dims[0], x.dims[1]);
    let (ny, my) = (y.dims[0], y.dims[1]);
    if ny != n {
        return rt("X and Y must have the same number of rows (points).");
    }
    if my > m {
        return rt("Y cannot have more columns (variables) than X.");
    }
    let mu = |m: &Matrix, c: usize| (0..m.dims[0]).map(|i| m.data[c * m.dims[0] + i]).sum::<f64>() / m.dims[0] as f64;
    let mux: Vec<f64> = (0..m).map(|c| mu(x, c)).collect();
    let muy: Vec<f64> = (0..my).map(|c| mu(y, c)).collect();
    let x0: Vec<f64> = (0..n * m).map(|i| x.data[i] - mux[i / n]).collect();
    let y0: Vec<f64> = (0..n * my).map(|i| y.data[i] - muy[i / n]).collect();
    let ssq_cols_x: Vec<f64> = (0..m).map(|c| (0..n).map(|i| x0[c * n + i].powi(2)).sum()).collect();
    let ssq_cols_y: Vec<f64> = (0..my).map(|c| (0..n).map(|i| y0[c * n + i].powi(2)).sum()).collect();
    let epsx = if single { f32::EPSILON as f64 } else { f64::EPSILON };
    let constx = (0..m).all(|c| ssq_cols_x[c] <= (epsx * n as f64 * mux[c]).abs().powi(2));
    let consty = (0..my).all(|c| ssq_cols_y[c] <= (epsx * n as f64 * muy[c]).abs().powi(2));
    let (ssqx, ssqy): (f64, f64) = (ssq_cols_x.iter().sum(), ssq_cols_y.iter().sum());
    let mkst = |t: Matrix, b: f64, c: Matrix| -> Value {
        let mut st = StructArray::scalar();
        st.set("T", Value::Num(finish(t, single)));
        st.set("b", Value::Num(finish(Matrix::scalar(b), single)));
        st.set("c", Value::Num(finish(c, single)));
        Value::Struct(st)
    };
    let rep = |v: &[f64], rows: usize| -> Matrix {
        let k = v.len();
        Matrix::from_dims(vec![rows, k], (0..rows * k).map(|i| v[i / rows]).collect())
    };
    if constx || consty {
        let z = rep(&mux, n);
        let t = Matrix::from_dims(vec![my, m], (0..my * m).map(|i| if i % my == i / my { 1.0 } else { 0.0 }).collect());
        return Ok(Some(vec![Value::Num(finish(Matrix::scalar(if constx { 0.0 } else { 1.0 }), single)), Value::Num(finish(z.clone(), single)), mkst(t, 0.0, z)]));
    }
    let (nx, ny_) = (ssqx.sqrt(), ssqy.sqrt());
    let x0n: Vec<f64> = x0.iter().map(|v| v / nx).collect();
    let mut y0n: Vec<f64> = y0.iter().map(|v| v / ny_).collect();
    if my < m {
        y0n.extend(vec![0.0; n * (m - my)]);
    }
    // A = X0' * Y0 (m x m)
    let amat: Vec<f64> = (0..m * m).map(|idx| { let (i, j) = (idx % m, idx / m); (0..n).map(|r| x0n[i * n + r] * y0n[j * n + r]).sum() }).collect();
    let sv = it.call_function("svd", vec![Value::Num(Matrix::from_dims(vec![m, m], amat))], 3)?;
    let (Value::Num(l), Value::Num(d), Value::Num(mm)) = (&sv[0], &sv[1], &sv[2]) else { return rt("svd failed") };
    let mut mvec = mm.data.clone();
    let mut dvec: Vec<f64> = (0..m).map(|i| d.data[i * m + i]).collect();
    let tmat = |mv: &Vec<f64>| -> Vec<f64> { (0..m * m).map(|idx| { let (i, j) = (idx % m, idx / m); (0..m).map(|k| mv[k * m + i] * l.data[k * m + j]).sum() }).collect() };
    let mut t = tmat(&mvec);
    if let Some(want) = reflection {
        let det = it.call_function("det", vec![Value::Num(Matrix::from_dims(vec![m, m], t.clone()))], 1)?[0].scalar()?;
        if want != (det < 0.0) {
            for i in 0..m {
                mvec[(m - 1) * m + i] = -mvec[(m - 1) * m + i];
            }
            dvec[m - 1] = -dvec[m - 1];
            t = tmat(&mvec);
        }
    }
    let trace_ta: f64 = dvec.iter().sum();
    let (b, dd) = if scaling { (trace_ta * nx / ny_, 1.0 - trace_ta * trace_ta) } else { (1.0, 1.0 + ssqy / ssqx - 2.0 * trace_ta * ny_ / nx) };
    let scale = if scaling { nx * trace_ta } else { ny_ };
    // Z = scale * Y0 * T + muX
    let mut z = vec![0.0; n * m];
    for c in 0..m {
        for r in 0..n {
            z[c * n + r] = scale * (0..m).map(|k| y0n[k * n + r] * t[c * m + k]).sum::<f64>() + mux[c];
        }
    }
    let t_out: Matrix = if my < m { Matrix::from_dims(vec![my, m], (0..my * m).map(|idx| { let (i, j) = (idx % my, idx / my); t[j * m + i] }).collect()) } else { Matrix::from_dims(vec![m, m], t.clone()) };
    let cvec: Vec<f64> = (0..m).map(|j| mux[j] - b * (0..my).map(|k| muy[k] * t_out.data[j * my + k]).sum::<f64>()).collect();
    Ok(Some(vec![Value::Num(finish(Matrix::scalar(dd), single)), Value::Num(finish(Matrix::from_dims(vec![n, m], z), single)), mkst(t_out, b, rep(&cvec, n))]))
}

// ---- hidden Markov models -------------------------------------------------------------------------

struct Hmm {
    a: Vec<Vec<f64>>,
    e: Vec<Vec<f64>>,
    seq: Vec<usize>,
    names: Option<Value>,
}

fn hmm_inputs(a: &[Value], who: &str) -> Res<Hmm> {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    let (tr, em) = (num(a, 1)?, num(a, 2)?);
    let l = tr.dims[0];
    if tr.dims.len() != 2 || tr.dims[1] != l {
        return rt("TRANSITION matrix must be square.");
    }
    if em.dims[0] != l {
        return rt("EMISSIONS matrix must have the same number of rows as TRANSITIONS.");
    }
    let nsym = em.dims[1];
    let mut symbols: Option<Value> = None;
    let mut names = None;
    let rest = &a[3..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        if matches_name(&key, "symbols") {
            symbols = Some(pair[1].clone());
        } else if matches_name(&key, "statenames") && who == "hmmviterbi" {
            names = Some(pair[1].clone());
        } else {
            return rt(format!("Invalid parameter name: {key}."));
        }
    }
    let seq_idx: Vec<usize> = match (&a[0], &symbols) {
        (Value::Num(s), None) => {
            if s.is_empty() && who == "hmmviterbi" {
                return rt("Output argument \"logP\" (and possibly others) not assigned a value in the execution with \"hmmviterbi\" function.");
            }
            let mut out = Vec::new();
            for &v in &s.data {
                if v.fract() != 0.0 || v < 1.0 || v > nsym as f64 {
                    return rt(if who == "hmmdecode" { format!("SEQ must be a vector of integers between 1 and {nsym}.") } else { format!("SEQ must consist of integers between 1 and {nsym}.") });
                }
                out.push(v as usize - 1);
            }
            out
        }
        (Value::Num(s), Some(Value::Num(sym))) => {
            let mut out = Vec::new();
            for &v in &s.data {
                match sym.data.iter().position(|&x| x == v) {
                    Some(p) => out.push(p),
                    None => return rt(format!("SEQ must consist of integers between 1 and {nsym}.")),
                }
            }
            out
        }
        (Value::Cell(c), Some(Value::Cell(sym))) => {
            let syms: Vec<String> = sym.data.iter().filter_map(|e| e.text().ok()).collect();
            let mut out = Vec::new();
            for e in &c.data {
                let t = e.text()?;
                match syms.iter().position(|x| *x == t) {
                    Some(p) => out.push(p),
                    None => return rt("SEQ contains a symbol that is not in SYMBOLS."),
                }
            }
            out
        }
        _ => return rt("SEQ must be a numeric vector or a cell array of symbols."),
    };
    let at = |m: &Matrix, i: usize, j: usize| m.data[j * m.dims[0] + i];
    Ok(Hmm { a: (0..l).map(|i| (0..l).map(|j| at(tr, i, j)).collect()).collect(), e: (0..l).map(|i| (0..nsym).map(|j| at(em, i, j)).collect()).collect(), seq: seq_idx, names })
}

fn hmmdecode(a: &[Value]) -> Ret {
    let h = hmm_inputs(a, "hmmdecode")?;
    let (l, len) = (h.a.len(), h.seq.len());
    let mut fs = vec![vec![0.0; l]; len + 1];
    let mut s = vec![1.0; len + 1];
    fs[0][0] = 1.0;
    for t in 0..len {
        let mut f: Vec<f64> = (0..l).map(|j| h.e[j][h.seq[t]] * (0..l).map(|i| h.a[i][j] * fs[t][i]).sum::<f64>()).collect();
        let sc: f64 = f.iter().sum();
        s[t + 1] = sc;
        f.iter_mut().for_each(|v| *v /= sc);
        fs[t + 1] = f;
    }
    let mut bs = vec![vec![1.0; l]; len + 1];
    for t in (0..len).rev() {
        bs[t] = (0..l).map(|i| (0..l).map(|j| h.a[i][j] * h.e[j][h.seq[t]] * bs[t + 1][j]).sum::<f64>() / s[t + 1]).collect();
    }
    let logp: f64 = s[1..].iter().map(|v| v.ln()).sum();
    let mut ps = vec![0.0; l * len];
    for t in 0..len {
        for i in 0..l {
            ps[t * l + i] = fs[t + 1][i] * bs[t + 1][i];
        }
    }
    let to_m = |v: &Vec<Vec<f64>>| Value::Num(Matrix::from_dims(vec![l, len + 1], v.iter().flat_map(|c| c.iter().copied()).collect()));
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![l, len], ps)), Value::num(logp), to_m(&fs), to_m(&bs), Value::Num(Matrix::from_dims(vec![1, len + 1], s))]))
}

fn hmmviterbi(a: &[Value]) -> Ret {
    let h = hmm_inputs(a, "hmmviterbi")?;
    let (l, len) = (h.a.len(), h.seq.len());
    let lg = |v: f64| if v > 0.0 { v.ln() } else { f64::NEG_INFINITY };
    let mut cur: Vec<f64> = (0..l).map(|j| if j == 0 { 0.0 } else { f64::NEG_INFINITY }).collect();
    let mut ptr = vec![vec![0usize; l]; len];
    for t in 0..len {
        let mut nxt = vec![f64::NEG_INFINITY; l];
        for j in 0..l {
            let mut best = (f64::NEG_INFINITY, 0usize);
            for i in 0..l {
                let v = cur[i] + lg(h.a[i][j]);
                if v > best.0 {
                    best = (v, i);
                }
            }
            nxt[j] = best.0 + lg(h.e[j][h.seq[t]]);
            ptr[t][j] = best.1;
        }
        cur = nxt;
    }
    let mut best = (f64::NEG_INFINITY, 0usize);
    for (j, &v) in cur.iter().enumerate() {
        if v > best.0 {
            best = (v, j);
        }
    }
    let mut states = vec![0usize; len];
    let mut st = best.1;
    for t in (0..len).rev() {
        states[t] = st;
        st = ptr[t][st];
    }
    let out = match &h.names {
        Some(Value::Cell(c)) => Value::Cell(CellArray::from_dims(vec![1, len], states.iter().map(|&s| c.data[s].clone()).collect())),
        _ => Value::Num(Matrix::from_dims(vec![1, len], states.iter().map(|&s| s as f64 + 1.0).collect())),
    };
    Ok(Some(vec![out, Value::num(best.0)]))
}

fn hmmestimate(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let rest = &a[2..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let mut symbols: Option<Value> = None;
    let mut pe: Option<Matrix> = None;
    let mut pt: Option<Matrix> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        if matches_name(&key, "symbols") {
            symbols = Some(pair[1].clone());
        } else if matches_name(&key, "pseudoemissions") {
            pe = Some(pair[1].as_matrix()?.clone());
        } else if matches_name(&key, "pseudotransitions") {
            pt = Some(pair[1].as_matrix()?.clone());
        } else if matches_name(&key, "statenames") {
            return rt("Input A of class double and input B of class cell must be cell arrays of character vectors, unless one is a character vector.");
        } else {
            return rt(format!("Invalid parameter name: {key}."));
        }
    }
    let Value::Num(st) = &a[1] else { return rt("STATES must be numeric.") };
    let states: Vec<usize> = st.data.iter().map(|&v| v as usize - 1).collect();
    let seq: Vec<usize> = match (&a[0], &symbols) {
        (Value::Num(s), None) => s.data.iter().map(|&v| v as usize - 1).collect(),
        (Value::Num(s), Some(Value::Num(sym))) => s.data.iter().map(|v| sym.data.iter().position(|x| x == v).unwrap_or(0)).collect(),
        (Value::Cell(_), _) => return rt("Unable to use a value of type cell as an index."),
        _ => return rt("SEQ must be numeric."),
    };
    if seq.len() != states.len() {
        return rt("Input sequence and states must be the same length.");
    }
    let nstates = states.iter().copied().max().map_or(0, |m| m + 1).max(pt.as_ref().map_or(0, |m| m.dims[0])).max(pe.as_ref().map_or(0, |m| m.dims[0]));
    let nsym = match (&symbols, &pe) {
        (Some(Value::Num(s)), _) => s.numel(),
        (_, Some(p)) => p.dims[1],
        _ => seq.iter().copied().max().map_or(0, |m| m + 1),
    };
    if seq.iter().any(|&s| s >= nsym) {
        return rt("There are more symbols in SEQ than in PSEUDO.");
    }
    let mut tr = vec![0.0; nstates * nstates];
    let mut em = vec![0.0; nstates * nsym];
    if let Some(p) = &pt {
        for i in 0..p.dims[0].min(nstates) {
            for j in 0..p.dims[1].min(nstates) {
                tr[j * nstates + i] += p.data[j * p.dims[0] + i];
            }
        }
    }
    if let Some(p) = &pe {
        for i in 0..p.dims[0].min(nstates) {
            for j in 0..p.dims[1].min(nsym) {
                em[j * nstates + i] += p.data[j * p.dims[0] + i];
            }
        }
    }
    for w in states.windows(2) {
        tr[w[1] * nstates + w[0]] += 1.0;
    }
    for (s, y) in states.iter().zip(&seq) {
        em[y * nstates + s] += 1.0;
    }
    let norm = |m: &mut Vec<f64>, rows: usize, cols: usize| {
        for i in 0..rows {
            let t: f64 = (0..cols).map(|j| m[j * rows + i]).sum();
            for j in 0..cols {
                m[j * rows + i] = if t > 0.0 { m[j * rows + i] / t } else { 0.0 };
            }
        }
    };
    norm(&mut tr, nstates, nstates);
    norm(&mut em, nstates, nsym);
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![nstates, nstates], tr)), Value::Num(Matrix::from_dims(vec![nstates, nsym], em))]))
}

// ---- barttest / bbdesign -------------------------------------------------------------------------

fn barttest(it: &mut Interpreter, a: &[Value]) -> Ret {
    let x = num(a, 0)?;
    let single = is_single(x);
    let alpha = match a.get(1) {
        Some(v) => {
            let al = v.scalar().unwrap_or(f64::NAN);
            if !(al > 0.0 && al < 1.0) || v.as_matrix().map(|m| m.numel() != 1).unwrap_or(true) {
                return rt("ALPHA must be a scalar between 0 and 1.");
            }
            al
        }
        None => 0.05,
    };
    let (nu, n) = (x.dims[0], x.dims[1]);
    if n <= 1 || nu <= 1 {
        return rt("X must have more than one row and column.");
    }
    let xc: Vec<f64> = (0..nu * n).map(|i| { let c = i / nu; x.data[i] - (0..nu).map(|r| x.data[c * nu + r]).sum::<f64>() / nu as f64 }).collect();
    let sv = it.call_function("svd", vec![Value::Num(Matrix::from_dims(vec![nu, n], xc))], 1)?;
    let Value::Num(s) = &sv[0] else { return rt("svd failed") };
    let mut latent: Vec<f64> = s.data.iter().map(|v| v * v / (nu as f64 - 1.0)).collect();
    latent.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let nu1 = (nu - 1) as f64;
    let lnl: Vec<f64> = latent.iter().map(|v| v.ln()).collect();
    // flipud(cumsum(log(latent)))
    let mut cum = vec![0.0; n];
    let mut acc = 0.0;
    for i in 0..n {
        acc += lnl[i];
        cum[i] = acc;
    }
    let lnlatent: Vec<f64> = cum.iter().rev().copied().collect();
    let k: Vec<f64> = (0..n - 1).map(|v| v as f64).collect();
    let pk: Vec<f64> = k.iter().map(|v| n as f64 - v).collect();
    // logsum = log(flipud((latent(1)+cumsum(latent(2:n)))./flipud(pk)))
    let mut cs = vec![0.0; n - 1];
    let mut acc2 = 0.0;
    for i in 1..n {
        acc2 += latent[i];
        cs[i - 1] = latent[0] + acc2;
    }
    let pk_rev: Vec<f64> = pk.iter().rev().copied().collect();
    let ratio: Vec<f64> = (0..n - 1).map(|i| cs[i] / pk_rev[i]).collect();
    let logsum: Vec<f64> = ratio.iter().rev().map(|v| v.ln()).collect();
    let chisq: Vec<f64> = (0..n - 1).map(|i| (pk[i] * logsum[i] - lnlatent[i]) * nu1).collect();
    let df: Vec<f64> = pk.iter().map(|p| (p - 1.0) * (p + 2.0) / 2.0).collect();
    let mut prob = Vec::new();
    for i in 0..n - 1 {
        prob.push(it.call_function("chi2cdf", vec![Value::num(chisq[i]), Value::num(df[i]), Value::chars("upper")], 1)?[0].scalar()?);
    }
    let col = |v: Vec<f64>| Value::Num(finish(Matrix::from_dims(vec![v.len(), 1], v), single));
    let ndim = match prob.iter().position(|&p| p > alpha) {
        None => n as f64,
        Some(0) => f64::NAN,
        Some(d) => d as f64,
    };
    Ok(Some(vec![Value::num(ndim), col(prob), col(chisq)]))
}

fn bbdesign(a: &[Value]) -> Ret {
    let nf = a.first().map(|v| v.scalar()).transpose()?;
    let Some(nf) = nf else { return rt("Not enough input arguments.") };
    if nf.fract() != 0.0 || nf < 2.0 {
        return rt("Size arguments must be real integers.");
    }
    let nf = nf as usize;
    if nf < 3 {
        return rt(format!("Cannot generate Box-Behnken design for {nf} factors."));
    }
    let rest = &a[1..];
    if rest.len() % 2 != 0 {
        return rt("Invalid parameter name/value pair.");
    }
    let mut center: Option<usize> = None;
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        if matches_name(&key, "center") {
            center = Some(pair[1].scalar()? as usize);
        } else if matches_name(&key, "blocksize") {
            return rt("The 'blocksize' option is not implemented.");
        } else {
            return rt(format!("Invalid parameter name: {key}."));
        }
    }
    if nf > 5 {
        return rt("Box-Behnken designs for more than 5 factors are not implemented.");
    }
    let mut pairs: Vec<(usize, usize, usize)> = Vec::new(); // (i, j, block)
    let nblocks;
    if nf == 4 {
        nblocks = 3;
        for (k, &(i, j)) in [(0, 1), (2, 3), (0, 3), (1, 2), (0, 2), (1, 3)].iter().enumerate() {
            pairs.push((i, j, k / 2 + 1));
        }
    } else if nf == 5 {
        nblocks = 2;
        for (k, &(i, j)) in [(0, 1), (2, 3), (1, 4), (0, 2), (3, 4), (1, 2), (0, 3), (2, 4), (0, 4), (1, 3)].iter().enumerate() {
            pairs.push((i, j, k / 5 + 1));
        }
    } else {
        nblocks = 1;
        for i in 0..nf {
            for j in i + 1..nf {
                pairs.push((i, j, 1));
            }
        }
    }
    let ncen = center.unwrap_or(match nf {
        3 => 3,
        4 => 3,
        _ => 6,
    });
    let nruns = pairs.len() * 4 + ncen;
    let mut d = vec![0.0; nruns * nf];
    let mut blocks = vec![1.0; nruns];
    let signs = [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)];
    for (pi, &(i, j, blk)) in pairs.iter().enumerate() {
        for (k, &(s1, s2)) in signs.iter().enumerate() {
            let r = pi * 4 + k;
            d[i * nruns + r] = s1;
            d[j * nruns + r] = s2;
            blocks[r] = blk as f64;
        }
    }
    for c in 0..ncen {
        blocks[pairs.len() * 4 + c] = if nblocks == 1 { 1.0 } else if nf == 5 { (c * 2 / ncen + 1) as f64 } else { (c % nblocks + 1) as f64 };
    }
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![nruns, nf], d)), Value::Num(Matrix::from_dims(vec![nruns, 1], blocks))]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    match name {
        "dbscan" => dbscan(it, a),
        "procrustes" => procrustes(it, a),
        "hmmdecode" => hmmdecode(a),
        "hmmviterbi" => hmmviterbi(a),
        "hmmestimate" => hmmestimate(a),
        "barttest" => barttest(it, a),
        _ => bbdesign(a),
    }
}
