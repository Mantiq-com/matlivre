//! `kmeans` for the squared Euclidean, cosine and correlation distances: batch reassignment
//! followed (optionally) by the online phase, with `Start`, `Replicates`, `EmptyAction`,
//! `OnlinePhase` and `MaxIter`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::random::global_uniform;
use crate::value::*;

pub const NAMES: &[&str] = &["kmeans"];

#[derive(Clone, Copy, PartialEq)]
enum Dist {
    Sqe,
    Cos,
}

fn prefix(given: &str, list: &[&str]) -> Option<usize> {
    let g = given.to_lowercase();
    if g.is_empty() {
        return None;
    }
    let hits: Vec<usize> = list.iter().enumerate().filter(|(_, c)| c.starts_with(&g)).map(|(i, _)| i).collect();
    if hits.len() == 1 { Some(hits[0]) } else { None }
}

fn dist_to(x: &[f64], c: &[f64], d: Dist, cnorm: f64) -> f64 {
    match d {
        Dist::Sqe => x.iter().zip(c).map(|(a, b)| (a - b) * (a - b)).sum(),
        Dist::Cos => 1.0 - x.iter().zip(c).map(|(a, b)| a * b / cnorm).sum::<f64>(),
    }
}

struct State<'a> {
    x: &'a [Vec<f64>],
    n: usize,
    k: usize,
    dist: Dist,
}

impl State<'_> {
    fn column_dist(&self, c: &[f64]) -> Result<Vec<f64>, crate::error::LangError> {
        let nrm = c.iter().map(|v| v * v).sum::<f64>().sqrt();
        if self.dist == Dist::Cos && nrm < f64::EPSILON {
            return rt("Zero cluster centroid created at iteration.");
        }
        Ok(self.x.iter().map(|r| dist_to(r, c, self.dist, if self.dist == Dist::Cos { nrm } else { 1.0 })).collect())
    }
}

fn centroid(x: &[Vec<f64>], idx: &[usize], cl: usize, p: usize) -> (Vec<f64>, f64) {
    let mut s = vec![0.0; p];
    let mut cnt = 0.0;
    for (i, &g) in idx.iter().enumerate() {
        if g == cl {
            cnt += 1.0;
            for j in 0..p {
                s[j] += x[i][j];
            }
        }
    }
    if cnt == 0.0 {
        return (vec![f64::NAN; p], 0.0);
    }
    (s.iter().map(|v| v / cnt).collect(), cnt)
}

struct Rep {
    total: f64,
    sumd: Vec<f64>,
    idx: Vec<usize>,
    c: Vec<Vec<f64>>,
    d: Vec<Vec<f64>>,
}

fn run_rep(x: &[Vec<f64>], k: usize, p: usize, dist: Dist, c0: Vec<Vec<f64>>, empty: usize, online: bool, maxit: usize, rep: usize, reps: usize) -> Result<Rep, crate::error::LangError> {
    let n = x.len();
    let st = State { x, n, k, dist };
    let mut c = c0;
    let mut d: Vec<Vec<f64>> = vec![vec![0.0; k]; n];
    for j in 0..k {
        let col = st.column_dist(&c[j])?;
        for i in 0..n {
            d[i][j] = col[i];
        }
    }
    let mut idx: Vec<usize> = (0..n).map(|i| (0..k).fold(0, |b, j| if d[i][j] < d[i][b] { j } else { b })).collect();
    let mut m: Vec<f64> = (0..k).map(|j| idx.iter().filter(|&&g| g == j).count() as f64).collect();
    let rowd = |d: &Vec<Vec<f64>>, idx: &[usize]| -> Vec<f64> { (0..n).map(|i| d[i][idx[i]]).collect() };
    // ---- batch phase
    let mut changed: Vec<usize> = (0..k).collect();
    let mut previdx = vec![usize::MAX; n];
    let mut prevtot = f64::INFINITY;
    let mut iter = 0usize;
    let mut totsum;
    loop {
        iter += 1;
        for &cl in &changed {
            let (cc, cnt) = centroid(x, &idx, cl, p);
            c[cl] = cc;
            m[cl] = cnt;
        }
        for &cl in &changed {
            if m[cl] > 0.0 || dist == Dist::Sqe {
                if m[cl] == 0.0 {
                    for i in 0..n {
                        d[i][cl] = f64::NAN;
                    }
                } else {
                    let col = st.column_dist(&c[cl])?;
                    for i in 0..n {
                        d[i][cl] = col[i];
                    }
                }
            } else {
                for i in 0..n {
                    d[i][cl] = f64::NAN;
                }
            }
        }
        let empties: Vec<usize> = changed.iter().copied().filter(|&cl| m[cl] == 0.0).collect();
        if !empties.is_empty() {
            match empty {
                0 => {
                    return if reps == 1 { rt(format!("Empty cluster created at iteration {iter}.")) } else { rt(format!("Empty cluster created at iteration {iter} during replicate {rep}.")) };
                }
                1 => {
                    changed.retain(|&cl| m[cl] > 0.0);
                }
                _ => {
                    for &i in &empties {
                        let dd = rowd(&d, &idx);
                        let mut lonely = 0;
                        for t in 0..n {
                            if dd[t] > dd[lonely] || dd[lonely].is_nan() && !dd[t].is_nan() {
                                lonely = t;
                            }
                        }
                        let mut from = idx[lonely];
                        if m[from] < 2.0 {
                            from = (0..k).find(|&j| m[j] > 1.0).unwrap_or(from);
                            lonely = idx.iter().position(|&g| g == from).unwrap_or(lonely);
                        }
                        c[i] = x[lonely].clone();
                        m[i] = 1.0;
                        idx[lonely] = i;
                        let col = st.column_dist(&c[i])?;
                        for t in 0..n {
                            d[t][i] = col[t];
                        }
                        let (cc, cnt) = centroid(x, &idx, from, p);
                        c[from] = cc;
                        m[from] = cnt;
                        let col = st.column_dist(&c[from])?;
                        for t in 0..n {
                            d[t][from] = col[t];
                        }
                        if !changed.contains(&from) {
                            changed.push(from);
                        }
                        changed.sort();
                        changed.dedup();
                    }
                }
            }
        }
        totsum = rowd(&d, &idx).iter().sum::<f64>();
        if prevtot <= totsum {
            idx = previdx.clone();
            for &cl in &changed {
                let (cc, cnt) = centroid(x, &idx, cl, p);
                c[cl] = cc;
                m[cl] = cnt;
            }
            iter -= 1;
            break;
        }
        if iter >= maxit {
            break;
        }
        previdx = idx.clone();
        prevtot = totsum;
        let nidx: Vec<usize> = (0..n).map(|i| (0..k).fold(0, |b, j| if d[i][j] < d[i][b] || d[i][b].is_nan() && !d[i][j].is_nan() { j } else { b })).collect();
        let moved: Vec<usize> = (0..n).filter(|&i| nidx[i] != previdx[i] && d[i][previdx[i]] > d[i][nidx[i]]).collect();
        if moved.is_empty() {
            break;
        }
        for &i in &moved {
            idx[i] = nidx[i];
        }
        let mut ch: Vec<usize> = moved.iter().flat_map(|&i| [idx[i], previdx[i]]).collect();
        ch.sort();
        ch.dedup();
        changed = ch;
    }
    let _ = totsum;
    // ---- online phase
    if online {
        let mut changed: Vec<usize> = (0..k).filter(|&j| m[j] > 0.0).collect();
        let mut lastmoved = 0usize;
        let mut nummoved = 0usize;
        let iter1 = iter;
        let mut del = vec![vec![f64::NAN; k]; n];
        let mut totsumd = rowd(&d, &idx).iter().sum::<f64>();
        while iter < maxit {
            match dist {
                Dist::Sqe => {
                    for &i in &changed {
                        for t in 0..n {
                            let mbr = idx[t] == i;
                            let sgn = if m[i] == 1.0 && mbr { 0.0 } else if mbr { -1.0 } else { 1.0 };
                            let dd: f64 = x[t].iter().zip(&c[i]).map(|(a, b)| (a - b) * (a - b)).sum();
                            del[t][i] = (m[i] / (m[i] + sgn)) * dd;
                        }
                    }
                }
                Dist::Cos => {
                    let norms: Vec<f64> = c.iter().map(|cc| cc.iter().map(|v| v * v).sum::<f64>().sqrt()).collect();
                    if norms.iter().any(|&v| v < f64::EPSILON) {
                        return rt(format!("Zero cluster centroid created at iteration {iter}."));
                    }
                    for &i in &changed {
                        for t in 0..n {
                            let xc: f64 = c[i].iter().zip(&x[t]).map(|(a, b)| a * b).sum();
                            let sgn = if idx[t] == i { -1.0 } else { 1.0 };
                            let mn = m[i] * norms[i];
                            del[t][i] = 1.0 + sgn * (mn - (mn * mn + 2.0 * sgn * m[i] * xc + 1.0).sqrt());
                        }
                    }
                }
            }
            let previdx2 = idx.clone();
            let nidx: Vec<usize> = (0..n).map(|i| (0..k).fold(0, |b, j| if del[i][j] < del[i][b] || del[i][b].is_nan() && !del[i][j].is_nan() { j } else { b })).collect();
            let mut moved: Vec<usize> = (0..n).filter(|&i| previdx2[i] != nidx[i] && m[previdx2[i]] != 1.0).collect();
            moved.retain(|&i| del[i][previdx2[i]] > del[i][nidx[i]]);
            if moved.is_empty() {
                if iter == iter1 || nummoved > 0 {
                    iter += 1;
                }
                break;
            }
            // first moved index at or after lastmoved (cyclically), 1-based arithmetic
            let nn = n as i64;
            let lm = lastmoved as i64;
            let best = moved.iter().map(|&i| (((i as i64 + 1) - lm - 1).rem_euclid(nn)) + lm).min().unwrap();
            let mv1 = best.rem_euclid(nn) + 1;
            let mv = (mv1 - 1) as usize;
            if mv1 <= lm {
                iter += 1;
                if iter >= maxit {
                    break;
                }
                nummoved = 0;
            }
            nummoved += 1;
            lastmoved = mv1 as usize;
            let (o, nw) = (idx[mv], nidx[mv]);
            totsumd += del[mv][nw] - del[mv][o];
            idx[mv] = nw;
            m[nw] += 1.0;
            m[o] -= 1.0;
            for j in 0..p {
                c[nw][j] += (x[mv][j] - c[nw][j]) / m[nw];
                c[o][j] -= (x[mv][j] - c[o][j]) / m[o];
            }
            changed = if o < nw { vec![o, nw] } else { vec![nw, o] };
        }
        let _ = totsumd;
    }
    let nonempty: Vec<usize> = (0..k).filter(|&j| m[j] > 0.0).collect();
    for &cl in &nonempty {
        let col = st.column_dist(&c[cl])?;
        for i in 0..n {
            d[i][cl] = col[i];
        }
    }
    let dd = rowd(&d, &idx);
    let mut sumd = vec![0.0; k];
    for i in 0..n {
        sumd[idx[i]] += dd[i];
    }
    let total: f64 = nonempty.iter().map(|&j| sumd[j]).sum();
    Ok(Rep { total, sumd, idx, c, d })
}

pub fn call(it: &mut Interpreter, _name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("At least two input arguments required.");
    }
    let Some(Value::Num(x0)) = a.first() else { return rt("X must be numeric.") };
    if x0.is_complex() {
        return rt("X must be real.");
    }
    let single = is_single(x0);
    let (n0, p) = (x0.dims[0], x0.dims[1]);
    let all_rows: Vec<Vec<f64>> = (0..n0).map(|i| (0..p).map(|j| x0.data[j * n0 + i]).collect()).collect();
    let wasnan: Vec<bool> = all_rows.iter().map(|r| r.iter().any(|v| v.is_nan())).collect();
    let mut x: Vec<Vec<f64>> = all_rows.iter().zip(&wasnan).filter(|(_, &w)| !w).map(|(r, _)| r.clone()).collect();
    let n = x.len();
    let k_arg = match &a[1] {
        Value::Num(m) if m.numel() == 1 => Some(m.data[0]),
        Value::Num(m) if m.is_empty() => None,
        _ => Some(f64::NAN),
    };
    let rest = &a[2..];
    if rest.len() % 2 == 1 {
        return rt("Wrong number of arguments.");
    }
    let (mut distance, mut start, mut reps, mut empty, mut online, mut maxit): (usize, Option<&Value>, Option<f64>, usize, bool, f64) = (0, None, None, 2, false, 100.0);
    for pair in rest.chunks(2) {
        let key = pair[0].text().unwrap_or_default();
        match prefix(&key, &["distance", "start", "replicates", "emptyaction", "onlinephase", "options", "maxiter", "display"]) {
            Some(0) => {
                let t = pair[1].text().unwrap_or_default();
                match prefix(&t, &["sqeuclidean", "cityblock", "cosine", "correlation", "hamming"]) {
                    Some(i) => distance = i,
                    None => return rt(format!("'{t}' is not a valid value for the 'Distance' argument.  Valid values are: 'sqeuclidean', 'cityblock', 'cosine', 'correlation', 'hamming'.")),
                }
            }
            Some(1) => start = Some(&pair[1]),
            Some(2) => reps = pair[1].scalar().ok(),
            Some(3) => {
                let t = pair[1].text().unwrap_or_default();
                match prefix(&t, &["error", "drop", "singleton"]) {
                    Some(i) => empty = i,
                    None => return rt(format!("'{t}' is not a valid value for the 'EmptyAction' argument.  Valid values are: 'error', 'drop', 'singleton'.")),
                }
            }
            Some(4) => {
                let t = pair[1].text().unwrap_or_default();
                match prefix(&t, &["on", "off"]) {
                    Some(i) => online = i == 0,
                    None => return rt(format!("'{t}' is not a valid value for the 'OnlinePhase' argument.  Valid values are: 'on', 'off'.")),
                }
            }
            Some(5) => {
                if let Value::Struct(s) = &pair[1] {
                    if let Some(Value::Num(m)) = s.get("MaxIter") {
                        if !m.is_empty() {
                            maxit = m.data[0];
                        }
                    }
                }
            }
            Some(6) => maxit = pair[1].scalar().unwrap_or(100.0),
            Some(7) => {}
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    if distance == 4 && x.iter().any(|r| r.iter().any(|&v| v != 0.0 && v != 1.0)) {
        return rt("Non-binary data cannot be clustered using Hamming distance.");
    }
    if distance == 1 || distance == 4 {
        return rt("The cityblock and hamming distances of kmeans are not implemented.");
    }
    let dist = if distance == 0 { Dist::Sqe } else { Dist::Cos };
    match distance {
        2 => {
            let norms: Vec<f64> = x.iter().map(|r| r.iter().map(|v| v * v).sum::<f64>().sqrt()).collect();
            let mx = norms.iter().cloned().fold(0.0, f64::max);
            if norms.iter().cloned().fold(f64::INFINITY, f64::min) <= crate::stats::util::eps_of(mx) {
                return rt("Some points have small relative magnitudes, making them effectively zero.  Either remove those points, or choose a distance other than 'cosine'.");
            }
            for (r, nm) in x.iter_mut().zip(&norms) {
                r.iter_mut().for_each(|v| *v /= nm);
            }
        }
        3 => {
            for r in x.iter_mut() {
                let mean = r.iter().sum::<f64>() / p as f64;
                r.iter_mut().for_each(|v| *v -= mean);
            }
            let norms: Vec<f64> = x.iter().map(|r| r.iter().map(|v| v * v).sum::<f64>().sqrt()).collect();
            let mx = norms.iter().cloned().fold(0.0, f64::max);
            if norms.iter().cloned().fold(f64::INFINITY, f64::min) <= crate::stats::util::eps_of(mx) {
                return rt("Some points are essentially constant, making them effectively zero after centering.  Either remove those points, or choose a distance other than 'correlation'.");
            }
            for (r, nm) in x.iter_mut().zip(&norms) {
                r.iter_mut().for_each(|v| *v /= nm);
            }
        }
        _ => {}
    }
    let dist = if distance == 3 { Dist::Cos } else { dist };
    let start_v = start.cloned().unwrap_or_else(|| Value::chars("plus"));
    let mut k = k_arg;
    let mut cc: Option<Matrix> = None;
    let mut start_kind = "plus".to_string();
    if start_v.is_text() {
        let t = start_v.text()?;
        match prefix(&t, &["uniform", "sample", "cluster", "plus", "kmeans++"]) {
            Some(i) => start_kind = ["uniform", "sample", "cluster", "plus", "plus"][i].to_string(),
            None => return rt(format!("Unknown 'Start' parameter value:  {t}.")),
        }
    } else if let Value::Num(m) = &start_v {
        start_kind = "numeric".into();
        let (kr, kc) = (m.dims[0], m.dims[1]);
        match k {
            None => k = Some(kr as f64),
            Some(kv) if kv != kr as f64 => return rt("The 'Start' matrix must have K rows."),
            _ => {}
        }
        if kc != p {
            return rt("The 'Start' matrix must have the same number of columns as X.");
        }
        let r3 = m.dims.get(2).copied().unwrap_or(1);
        match reps {
            None => reps = Some(r3 as f64),
            Some(rv) if rv != r3 as f64 => return rt("The third dimension of the 'Start' array must match the 'replicates' parameter value."),
            _ => {}
        }
        cc = Some(m.clone());
    } else {
        return rt("The 'Start' parameter value must be a string or a numeric matrix or array.");
    }
    let kf = match k {
        Some(v) if v.is_finite() && v > 0.0 && v.fract() == 0.0 => v,
        _ => return rt("K must be a positive integer value."),
    };
    let kk = kf as usize;
    if n < kk {
        return rt("X must have more rows than the number of clusters.");
    }
    let reps_n = match reps {
        None => 1,
        Some(r) if r >= 1.0 && r.fract() == 0.0 => r as usize,
        _ => return rt("The 'replicates' parameter value must be a positive integer."),
    };
    let maxit = maxit as usize;
    let mut best: Option<Rep> = None;
    let mut empty_err = 0;
    for rep in 0..reps_n {
        let c0: Vec<Vec<f64>> = match start_kind.as_str() {
            "numeric" => {
                let m = cc.as_ref().unwrap();
                let page = kk * p;
                let mut c: Vec<Vec<f64>> = (0..kk).map(|i| (0..p).map(|j| m.data[rep * page + j * kk + i]).collect()).collect();
                if distance == 3 {
                    for r in c.iter_mut() {
                        let mean = r.iter().sum::<f64>() / p as f64;
                        r.iter_mut().for_each(|v| *v -= mean);
                    }
                }
                c
            }
            "uniform" => {
                let mins: Vec<f64> = (0..p).map(|j| x.iter().map(|r| r[j]).fold(f64::INFINITY, f64::min)).collect();
                let maxs: Vec<f64> = (0..p).map(|j| x.iter().map(|r| r[j]).fold(f64::NEG_INFINITY, f64::max)).collect();
                let u = global_uniform(it, kk * p);
                let mut c: Vec<Vec<f64>> = (0..kk).map(|i| (0..p).map(|j| mins[j] + u[j * kk + i] * (maxs[j] - mins[j])).collect()).collect();
                if distance == 3 {
                    // C - mean(C,1): per-column centering of the p-by-k centroid matrix
                    for j in 0..p {
                        let mean = c.iter().map(|r| r[j]).sum::<f64>() / kk as f64;
                        c.iter_mut().for_each(|r| r[j] -= mean);
                    }
                }
                c
            }
            "sample" => {
                let ix = super::sampling::sample_without_replacement(it, n, kk)?;
                ix.iter().map(|&i| x[i - 1].clone()).collect()
            }
            "cluster" => {
                let k0 = (kk as f64).max(0.1 * n as f64).floor() as usize;
                let ix = super::sampling::sample_without_replacement(it, n, k0)?;
                let sub: Vec<Vec<f64>> = ix.iter().map(|&i| x[i - 1].clone()).collect();
                if kk < k0 {
                    let s0 = super::sampling::sample_without_replacement(it, k0, kk)?;
                    let c0: Vec<Vec<f64>> = s0.iter().map(|&i| sub[i - 1].clone()).collect();
                    run_rep(&sub, kk, p, dist, c0, empty, online, maxit, 1, 1)?.c
                } else {
                    sub
                }
            }
            "plus" => {
                let mut c: Vec<Vec<f64>> = Vec::with_capacity(kk);
                let first = super::sampling::randi_vec(it, n, 1)[0];
                c.push(x[first - 1].clone());
                let mut min_d = vec![f64::INFINITY; n];
                for ii in 1..kk {
                    let last = c[ii - 1].clone();
                    let nrm = last.iter().map(|v| v * v).sum::<f64>().sqrt();
                    for t in 0..n {
                        let dd = dist_to(&x[t], &last, dist, if dist == Dist::Cos { nrm } else { 1.0 });
                        if dd < min_d[t] {
                            min_d[t] = dd;
                        }
                    }
                    let denom: f64 = min_d.iter().sum();
                    if denom == 0.0 || denom.is_infinite() || denom.is_nan() {
                        let ix = {
                            let args = vec![Value::num(n as f64), Value::num((kk - ii) as f64)];
                            let out = it.call_function("randperm", args, 1)?;
                            match out.first() {
                                Some(Value::Num(m)) => m.data.iter().map(|&v| v as usize).collect::<Vec<_>>(),
                                _ => return rt("randperm failed"),
                            }
                        };
                        for &i in &ix {
                            c.push(x[i - 1].clone());
                        }
                        break;
                    }
                    let pick = super::sampling::weighted_with_replacement(it, &min_d, 1)[0];
                    c.push(x[pick - 1].clone());
                }
                c
            }
            _ => return rt("This 'Start' option of kmeans is not implemented."),
        };
        match run_rep(&x, kk, p, dist, c0, empty, online, maxit, rep + 1, reps_n) {
            Ok(r) => {
                if best.as_ref().map_or(true, |b| r.total < b.total) {
                    best = Some(r);
                }
            }
            Err(e) => {
                if reps_n == 1 || !e.message.starts_with("Empty cluster") {
                    return Err(e);
                }
                empty_err += 1;
                if empty_err == reps_n {
                    return rt("Empty clusters were created in all replicates.");
                }
            }
        }
    }
    let b = best.unwrap();
    let mut idx_full = vec![f64::NAN; n0];
    let mut d_full = vec![vec![f64::NAN; kk]; n0];
    let mut t = 0;
    for i in 0..n0 {
        if !wasnan[i] {
            idx_full[i] = (b.idx[t] + 1) as f64;
            d_full[i] = b.d[t].clone();
            t += 1;
        }
    }
    let cmat = {
        let mut data = vec![0.0; kk * p];
        for i in 0..kk {
            for j in 0..p {
                data[j * kk + i] = b.c[i][j];
            }
        }
        Matrix::from_dims(vec![kk, p], data)
    };
    let dmat = {
        let mut data = vec![0.0; n0 * kk];
        for i in 0..n0 {
            for j in 0..kk {
                data[j * n0 + i] = d_full[i][j];
            }
        }
        Matrix::from_dims(vec![n0, kk], data)
    };
    let outs = vec![
        Value::Num(Matrix::from_dims(vec![n0, 1], idx_full)),
        Value::Num(finish(cmat, single)),
        Value::Num(finish(Matrix::from_dims(vec![kk, 1], b.sumd), single)),
        Value::Num(finish(dmat, single)),
    ];
    Ok(Some(outs.into_iter().take(nargout.max(1)).collect()))
}
