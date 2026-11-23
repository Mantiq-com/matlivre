//! Hierarchical clustering: `linkage`, `cluster`, `clusterdata`, `cophenet`, `inconsistent`.

use super::dist_matrix::pdist_core;
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["linkage", "cluster", "clusterdata", "cophenet", "inconsistent"];

const LINKAGE_METHODS: &[&str] = &["single", "nearest", "complete", "farthest", "average", "upgma", "weighted", "wpgma", "centroid", "upgmc", "median", "wpgmc", "ward's", "incremental"];

fn prefix_index(given: &str, list: &[&str]) -> Option<usize> {
    let g = given.to_lowercase();
    if g.is_empty() {
        return None;
    }
    if let Some(i) = list.iter().position(|c| *c == g) {
        return Some(i);
    }
    let hits: Vec<usize> = list.iter().enumerate().filter(|(_, c)| c.starts_with(&g)).map(|(i, _)| i).collect();
    if hits.len() == 1 { Some(hits[0]) } else { None }
}

/// Agglomerates the condensed distances `y` of `m` observations.
fn agglomerate(y: &[f64], m: usize, method: &str) -> Vec<[f64; 3]> {
    let squared = matches!(method, "ce" | "me" | "wa");
    let mut d = vec![vec![0.0; m]; m];
    let mut k = 0;
    for i in 0..m {
        for j in i + 1..m {
            let v = if squared { y[k] * y[k] } else { y[k] };
            d[i][j] = v;
            d[j][i] = v;
            k += 1;
        }
    }
    // Positions are kept in order; each holds a cluster id and its size.
    let mut ids: Vec<usize> = (1..=m).collect();
    let mut size: Vec<f64> = vec![1.0; m];
    let mut rows: Vec<[f64; 3]> = Vec::with_capacity(m - 1);
    for s in 0..m - 1 {
        let cur = ids.len();
        // Minimum over pairs (i<j) in row-major order; the first minimum wins.
        let (mut bi, mut bj, mut bv) = (0usize, 1usize, f64::NAN);
        let mut found = false;
        for i in 0..cur {
            for j in i + 1..cur {
                let v = if method == "av" { d[i][j] / (size[i] * size[j]) } else { d[i][j] };
                if v.is_nan() {
                    continue;
                }
                if !found || v <= bv {
                    bv = v;
                    bi = i;
                    bj = j;
                    found = true;
                }
            }
        }
        let (ni, nj) = (size[bi], size[bj]);
        rows.push([ids[bi] as f64, ids[bj] as f64, bv]);
        for u in 0..cur {
            if u == bi || u == bj {
                continue;
            }
            let (dui, duj) = (d[u][bi], d[u][bj]);
            let nu = size[u];
            let nv = match method {
                "si" => dui.min(duj),
                "co" => dui.max(duj),
                "av" => dui + duj,
                "we" => (dui + duj) / 2.0,
                "ce" => {
                    let kk = ni + nj;
                    (ni * dui + nj * duj - ni * nj * bv / kk) / kk
                }
                "me" => (dui + duj) / 2.0 - bv / 4.0,
                _ => ((nu + ni) * dui + (nu + nj) * duj - nu * bv) / (ni + nj + nu),
            };
            d[u][bi] = nv;
            d[bi][u] = nv;
        }
        size[bi] = ni + nj;
        ids[bi] = m + s + 1;
        // remove position bj
        ids.remove(bj);
        size.remove(bj);
        d.remove(bj);
        for r in d.iter_mut() {
            r.remove(bj);
        }
    }
    if squared {
        for r in rows.iter_mut() {
            r[2] = r[2].sqrt();
        }
    }
    for r in rows.iter_mut() {
        if r[0] > r[1] {
            r.swap(0, 1);
        }
    }
    rows
}

/// MATLAB's re-ordering of the rows by height for average, complete, weighted and Ward linkage.
fn rearrange(z: &mut Vec<[f64; 3]>) {
    let n = z.len();
    let mut idx: Vec<usize> = (0..n - 1).collect();
    idx.sort_by(|&a, &b| z[a][2].partial_cmp(&z[b][2]).unwrap_or(std::cmp::Ordering::Equal));
    idx.push(n - 1);
    let mut revidx = vec![0usize; n];
    for (k, &i) in idx.iter().enumerate() {
        revidx[i] = k + 1;
    }
    let map = |v: f64| -> f64 {
        let v = v as usize;
        if v <= n + 1 { v as f64 } else { (n + 1 + revidx[v - n - 2]) as f64 }
    };
    let mut out: Vec<[f64; 3]> = idx.iter().map(|&i| z[i]).collect();
    for r in out.iter_mut() {
        r[0] = map(r[0]);
        r[1] = map(r[1]);
    }
    // Non-chronological trees: a row that refers to a node created later is swapped down.
    let nl = n + 1;
    if out.iter().enumerate().any(|(i, r)| r[0].max(r[1]) > (nl + i) as f64) {
        let mut last = nl as f64;
        for i in 0..n - 1 {
            let tn = (nl + i + 1) as f64;
            if out[i][0] > last || out[i][1] > last {
                let h = (i + 1..n).find(|&h| out[h][0] <= last && out[h][1] <= last).unwrap_or(i + 1);
                out.swap(i, h);
                let nn = (nl + h + 1) as f64;
                // positions of the references to tn and nn in the (column-major) pair columns
                let mut to_nn = None;
                let mut to_tn = None;
                'search: for c in 0..2 {
                    for r in 0..n {
                        if to_nn.is_none() && out[r][c] == tn {
                            to_nn = Some((r, c));
                        }
                        if to_tn.is_none() && out[r][c] == nn {
                            to_tn = Some((r, c));
                        }
                        if to_nn.is_some() && to_tn.is_some() {
                            break 'search;
                        }
                    }
                }
                if let Some((r, c)) = to_nn {
                    out[r][c] = nn;
                }
                if let Some((r, c)) = to_tn {
                    out[r][c] = tn;
                }
            }
            last = tn;
        }
    }
    for r in out.iter_mut() {
        if r[0] > r[1] {
            r.swap(0, 1);
        }
    }
    *z = out;
}

fn z_matrix(z: &[[f64; 3]], single: bool) -> Value {
    let m = z.len();
    let mut data = vec![0.0; m * 3];
    for (i, r) in z.iter().enumerate() {
        for c in 0..3 {
            data[c * m + i] = r[c];
        }
    }
    Value::Num(finish(Matrix::from_dims(vec![m, 3], data), single))
}

fn linkage(it: &mut Interpreter, a: &[Value]) -> Ret {
    let Some(Value::Num(y)) = a.first() else { return rt("Not enough input arguments.") };
    let (k, n) = (y.dims[0], y.dims[1]);
    let m = ((2.0 * n as f64).sqrt()).ceil() as usize;
    let nargin = a.len();
    let data_input = k > 1;
    let nargs = if data_input { nargin.max(3) } else { nargin };
    let single = is_single(y);
    if nargs >= 3 && !data_input {
        if k < 2 {
            return rt("You have to have at least one distance (two observations) to do a linkage.");
        }
    }
    let passx = nargs >= 3;
    if !passx {
        if n < 1 {
            return rt("You have to have at least one distance (two observations) to do a linkage.");
        }
        if k != 1 || m * (m - 1) / 2 != n {
            return rt("The first input does not appear to be a distance matrix because its size is not compatible with the output of the PDIST function. A data matrix input must have more than one row.");
        }
    }
    let method_str: &str = if nargs == 1 {
        "single"
    } else {
        let Some(v) = a.get(1) else { return rt("Not enough input arguments.") };
        let t = if v.is_text() { v.text()? } else { String::new() };
        match prefix_index(&t, LINKAGE_METHODS) {
            Some(i) => LINKAGE_METHODS[i & !1],
            None => return rt(format!("'{t}' is not a valid value for the METHOD argument.")),
        }
    };
    let method = &method_str[..2];
    if passx && nargin > 3 && (nargin - 3) % 2 == 1 {
        return rt("Wrong number of arguments.");
    }
    let dist: Vec<f64> = if passx {
        let spec: Vec<Value> = match a.get(2) {
            None => vec![Value::chars("euclidean")],
            Some(Value::Func(f)) => vec![Value::Func(f.clone())],
            Some(v) if v.is_text() => vec![v.clone()],
            Some(Value::Cell(c)) => c.data.clone(),
            _ => return rt("The pdist arguments must be a character vector, function handle or cell array."),
        };
        let (d, _) = pdist_core(it, y, &spec)?;
        d
    } else {
        y.data.clone()
    };
    let mm = if passx { y.dims[0] } else { m };
    let mut z = agglomerate(&dist, mm, method);
    if matches!(method, "av" | "wa" | "co" | "we") && z.len() > 1 {
        rearrange(&mut z);
    }
    Ok(Some(vec![z_matrix(&z, single)]))
}

fn inconsistent_rows(z: &[[f64; 3]], depth: f64) -> Vec<[f64; 4]> {
    let m = z.len();
    let leaves = m + 1;
    z.iter()
        .enumerate()
        .map(|(k, row)| {
            let (mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0);
            let mut klist = vec![k];
            let mut dlist = vec![depth];
            let mut cur = 0;
            while cur < klist.len() {
                let kk = klist[cur];
                let dd = dlist[cur];
                s1 += z[kk][2];
                s2 += z[kk][2] * z[kk][2];
                s3 += 1.0;
                if dd > 1.0 {
                    for c in 0..2 {
                        let i = z[kk][c] as usize;
                        if i > leaves {
                            klist.push(i - leaves - 1 + 0);
                            dlist.push(dd - 1.0);
                        }
                    }
                }
                cur += 1;
            }
            let mean = s1 / s3;
            let v = (s2 - s1 * s1 / s3) / (s3 - if s3 != 1.0 { 1.0 } else { 0.0 });
            let sd = v.max(0.0).sqrt();
            let crit = if sd > 0.0 { (row[2] - mean) / sd } else { 0.0 };
            [mean, sd, s3, crit]
        })
        .collect()
}

fn read_z(v: Option<&Value>) -> Res<Vec<[f64; 3]>> {
    let Some(Value::Num(m)) = v else { return rt("Not enough input arguments.") };
    if m.dims.len() != 2 || m.dims[1] != 3 {
        return Err(crate::error::LangError::new("The first input must be a valid linkage matrix (a cluster tree) with three columns.", 0));
    }
    let n = m.dims[0];
    Ok((0..n).map(|i| [m.data[i], m.data[n + i], m.data[2 * n + i]]).collect())
}

fn inconsistent(a: &[Value]) -> Ret {
    if let Some(Value::Num(m)) = a.first() {
        if m.dims.len() == 2 && m.dims[1] < 3 && m.dims[0] > 0 {
            return rt(format!("Index in position 2 exceeds array bounds. Index must not exceed {}.", m.dims[1]));
        }
    }
    let z = read_z(a.first())?;
    let depth = match a.get(1) {
        Some(Value::Num(d)) if d.numel() == 1 => d.data[0],
        None => 2.0,
        _ => return rt("The depth must be a scalar."),
    };
    let single = matches!(a.first(), Some(Value::Num(m)) if is_single(m));
    let rows = inconsistent_rows(&z, depth);
    let m = rows.len();
    let mut data = vec![0.0; m * 4];
    for (i, r) in rows.iter().enumerate() {
        for c in 0..4 {
            data[c * m + i] = r[c];
        }
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![m, 4], data), single))]))
}

fn checkcut(z: &[[f64; 3]], cutoff: f64, crit: &[f64]) -> Vec<bool> {
    let n = z.len() + 1;
    let mut conn: Vec<bool> = crit.iter().map(|&c| c <= cutoff).collect();
    let mut todo: Vec<bool> = (0..z.len()).map(|i| conn[i] && (z[i][0] as usize > n || z[i][1] as usize > n)).collect();
    while todo.iter().any(|&t| t) {
        let rows: Vec<usize> = (0..z.len()).filter(|&i| todo[i]).collect();
        let mut cdone = vec![[true, true]; rows.len()];
        let mut newconn = conn.clone();
        for j in 0..2 {
            for (ri, &r) in rows.iter().enumerate() {
                let c = z[r][j] as usize;
                if c > n {
                    let child = c - n - 1;
                    cdone[ri][j] = !todo[child];
                    newconn[r] = newconn[r] && conn[child];
                }
            }
        }
        conn = newconn;
        for (ri, &r) in rows.iter().enumerate() {
            if cdone[ri][0] && cdone[ri][1] {
                todo[r] = false;
            }
        }
    }
    conn
}

fn labeltree(z: &[[f64; 3]], conn_in: &[bool]) -> Vec<f64> {
    let n = z.len();
    let leaves = n + 1;
    let mut conn = conn_in.to_vec();
    let mut t = vec![1.0; n + 1];
    let mut todo = vec![true; n];
    // clustlist(r,1) = r, clustlist(r,2) = n + r (1-based values)
    let mut cl = vec![[0.0; 2]; n];
    for r in 0..n {
        cl[r] = [(r + 1) as f64, (n + r + 1) as f64];
    }
    while todo.iter().any(|&x| x) {
        let rows: Vec<usize> = (0..n).filter(|&i| todo[i] && !conn[i]).collect();
        if rows.is_empty() {
            break;
        }
        for j in 0..2 {
            let mut set_children: Vec<(usize, f64)> = Vec::new();
            for &r in &rows {
                let c = z[r][j] as usize;
                if c <= leaves {
                    t[c - 1] = cl[r][j];
                } else if conn[c - leaves - 1] {
                    set_children.push((c - leaves - 1, cl[r][j]));
                }
            }
            for (child, num) in set_children {
                cl[child] = [num, num];
                conn[child] = false;
            }
        }
        for &r in &rows {
            todo[r] = false;
        }
    }
    // [~,~,T] = unique(T): rank of each value among the sorted distinct values
    let mut vals: Vec<f64> = t.clone();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vals.dedup();
    t.iter().map(|v| (vals.iter().position(|u| u == v).unwrap() + 1) as f64).collect()
}

fn clusternum(z: &[[f64; 3]], t: &mut [f64], k: usize, c: f64) {
    let m = z.len() + 1;
    let mut list = vec![k];
    while !list.is_empty() {
        let mut next = Vec::new();
        for kk in list {
            for col in 0..2 {
                let child = z[kk - 1][col] as usize;
                if child <= m {
                    t[child - 1] = c;
                } else {
                    next.push(child - m);
                }
            }
        }
        list = next;
    }
}

fn cluster(z: &[[f64; 3]], rest: &[Value]) -> Ret {
    if rest.is_empty() {
        return rt("Not enough input arguments.");
    }
    let mut depth = 2.0;
    let mut criterion = "inconsistent".to_string();
    let mut maxclust: Vec<f64> = vec![f64::INFINITY];
    let mut cutoff: Vec<f64> = vec![f64::NAN];
    let mut usecutoff = true;
    let nums = |v: &Value| -> Vec<f64> {
        match v {
            Value::Num(m) => m.data.clone(),
            _ => vec![],
        }
    };
    if !rest[0].is_text() {
        cutoff = nums(&rest[0]);
        if rest.len() >= 2 {
            depth = nums(&rest[1]).first().copied().unwrap_or(2.0);
        }
        if rest.len() < 3 {
            if cutoff.iter().all(|&c| c >= 2.0 && c == c.trunc()) {
                maxclust = cutoff.clone();
                cutoff = vec![f64::NAN];
                usecutoff = false;
            }
        } else {
            let flag = rest[2].text().unwrap_or_default();
            match prefix_index(&flag, &["clusters", "inconsistent", "distance"]) {
                None => return rt("The FLAG must be 'clusters', 'inconsistent', or 'distance'."),
                Some(0) => {
                    maxclust = cutoff.clone();
                    cutoff = vec![f64::NAN];
                    usecutoff = false;
                }
                Some(2) => criterion = "distance".into(),
                _ => {}
            }
        }
    } else {
        if rest.len() % 2 == 1 {
            return rt("Invalid parameter name or missing value.");
        }
        for pair in rest.chunks(2) {
            let key = pair[0].text().unwrap_or_default();
            match prefix_index(&key, &["cutoff", "maxclust", "criterion", "depth"]) {
                Some(0) => cutoff = nums(&pair[1]),
                Some(1) => maxclust = nums(&pair[1]),
                Some(2) => {
                    let t = pair[1].text().unwrap_or_default();
                    match prefix_index(&t, &["inconsistent", "distance"]) {
                        Some(i) => criterion = ["inconsistent", "distance"][i].into(),
                        None => return rt(format!("'{t}' is not a valid value for the Criterion argument.")),
                    }
                }
                Some(3) => depth = nums(&pair[1]).first().copied().unwrap_or(2.0),
                _ => return rt(format!("Invalid parameter name: {key}.")),
            }
        }
    }
    if maxclust.iter().any(|&m| m < 1.0 || (m.is_finite() && m.fract() != 0.0)) {
        return rt("The MAXCLUST value must be an integer number 1 or larger.");
    }
    if !depth.is_finite() || depth < 1.0 {
        return rt("The DEPTH value must be a finite number 1 or larger.");
    }
    if cutoff.iter().any(|c| c.is_nan()) {
        if maxclust.iter().any(|m| !m.is_finite()) {
            return rt("You must specify either the CUTOFF or the MAXCLUST value.");
        }
        usecutoff = false;
    } else if maxclust.iter().any(|m| !m.is_infinite()) {
        return rt("You cannot specify both the CUTOFF and the MAXCLUST values.");
    } else if cutoff.iter().any(|&c| c <= 0.0) {
        return rt("CUTOFF must be positive.");
    }
    let m = z.len() + 1;
    let mut cols: Vec<Vec<f64>> = Vec::new();
    if usecutoff {
        let crit: Vec<f64> = if criterion == "inconsistent" { inconsistent_rows(z, depth).iter().map(|r| r[3]).collect() } else { z.iter().map(|r| r[2]).collect() };
        for &c in &cutoff {
            let conn = checkcut(z, c, &crit);
            cols.push(labeltree(z, &conn));
        }
    } else {
        for &mc in &maxclust {
            let mcn = mc as usize;
            let mut t = vec![0.0; m];
            if m <= mcn {
                t = (1..=m).map(|v| v as f64).collect();
            } else if mcn == 1 {
                t = vec![1.0; m];
            } else {
                let mut clsnum = 1.0;
                for k in (m - mcn + 1)..m {
                    for col in 0..2 {
                        let i = z[k - 1][col] as usize;
                        if i <= m {
                            t[i - 1] = clsnum;
                            clsnum += 1.0;
                        } else if i < 2 * m - mcn + 1 {
                            clusternum(z, &mut t, i - m, clsnum);
                            clsnum += 1.0;
                        }
                    }
                }
            }
            cols.push(t);
        }
    }
    let nc = cols.len();
    let data: Vec<f64> = cols.into_iter().flatten().collect();
    Ok(Some(vec![Value::Num(Matrix::from_dims(vec![m, nc], data))]))
}

fn cophenet(a: &[Value], nargout: usize) -> Ret {
    let z = read_z(a.first()).map_err(|_| crate::error::LangError::new("Z must be a valid linkage matrix.", 0))?;
    let Some(Value::Num(yv)) = a.get(1) else { return rt("Not enough input arguments.") };
    let n = z.len() + 1;
    let ny = yv.numel();
    let mut link = vec![0usize; n + 1];
    let mut listhead: Vec<usize> = (0..=2 * n).map(|i| if i >= 1 && i <= n { i } else { 0 }).collect();
    let (mut sum1, mut sum2, mut s11, mut s22, mut s12) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut dout = vec![0.0; ny];
    let nf = n as f64;
    for (k, row) in z.iter().enumerate() {
        let (i, j, t) = (row[0] as usize, row[1] as usize, row[2]);
        let mut m1 = listhead[i];
        let mut msav = 0;
        while m1 > 0 {
            let mut m = listhead[j];
            while m > 0 {
                let (m2, m3) = (m1.min(m), m1.max(m));
                let idx = ((m2 as f64 - 1.0) * (nf - m2 as f64 / 2.0) + (m3 - m2) as f64) as usize;
                if idx < 1 || idx > ny {
                    return rt(format!("Index exceeds the number of array elements. Index must not exceed {ny}."));
                }
                dout[idx - 1] = t;
                let u = yv.data[idx - 1];
                sum1 += t;
                sum2 += u;
                s11 += t * t;
                s22 += u * u;
                s12 += t * u;
                msav = m;
                m = link[m];
            }
            m1 = link[m1];
        }
        link[msav] = listhead[i];
        listhead[n + k + 1] = listhead[j];
    }
    let t = 2.0 / (nf * (nf - 1.0));
    s11 -= sum1 * sum1 * t;
    s22 -= sum2 * sum2 * t;
    s12 -= sum1 * sum2 * t;
    let c = s12 / (s11 * s22).sqrt();
    let single = is_single(yv) || matches!(a.first(), Some(Value::Num(m)) if is_single(m));
    let mut outs = vec![Value::Num(finish(Matrix::scalar(c), single))];
    if nargout > 1 {
        outs.push(Value::Num(finish(Matrix::from_dims(yv.dims.clone(), dout), single)));
    }
    Ok(Some(outs))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "linkage" => linkage(it, a),
        "inconsistent" => inconsistent(a),
        "cophenet" => cophenet(a, nargout),
        "cluster" => {
            let z = read_z(a.first()).map_err(|_| crate::error::LangError::new("The first input must be a valid linkage matrix with three columns of cluster tree data.", 0));
            if a.len() < 2 {
                return rt("Not enough input arguments.");
            }
            cluster(&z?, &a[1..])
        }
        _ => {
            // clusterdata(X, ...)
            let Some(Value::Num(x)) = a.first() else { return rt("Not enough input arguments.") };
            if a.len() < 2 {
                return rt("Not enough input arguments.");
            }
            let mut pdist_args: Vec<Value> = Vec::new();
            let mut link_args: Vec<Value> = Vec::new();
            let mut rest: Vec<Value> = Vec::new();
            if a.len() == 2 {
                pdist_args.push(Value::chars("euclid"));
                link_args.push(Value::chars("single"));
                let c = match &a[1] {
                    Value::Num(m) if m.data.iter().all(|&v| v > 0.0) && !m.is_empty() => m.clone(),
                    _ => return rt("The cutoff must be a positive number."),
                };
                if c.data.iter().all(|&v| v < 2.0) {
                    rest = vec![Value::chars("cutoff"), Value::Num(c)];
                } else if c.data.iter().all(|&v| v.fract() == 0.0) {
                    rest = vec![Value::chars("maxclust"), Value::Num(c)];
                } else {
                    return rt("The cutoff must be a positive number.");
                }
            } else {
                let args = &a[1..];
                let mut done = vec![false; args.len()];
                let mut j = 0;
                while j + 1 < args.len() {
                    let pname = if args[j].is_text() { args[j].text()?.to_lowercase() } else { String::new() };
                    if !pname.is_empty() && "distance".starts_with(&pname) {
                        pdist_args.push(args[j + 1].clone());
                        done[j] = true;
                        done[j + 1] = true;
                        if args[j + 1].text().map(|t| t.eq_ignore_ascii_case("minkowski")).unwrap_or(false) {
                            if let Some(Value::Num(e)) = args.get(j + 2) {
                                if e.kind != Kind::Char {
                                    pdist_args.push(args[j + 2].clone());
                                    done[j + 2] = true;
                                    j += 1;
                                }
                            }
                        }
                    } else if !pname.is_empty() && "linkage".starts_with(&pname) {
                        link_args.push(args[j + 1].clone());
                        done[j] = true;
                        done[j + 1] = true;
                    }
                    j += 2;
                }
                for (i, v) in args.iter().enumerate() {
                    if !done[i] {
                        rest.push(v.clone());
                    }
                }
                if rest.is_empty() {
                    return rt("Missing parameter: specify the CUTOFF or the MAXCLUST value.");
                }
                if pdist_args.is_empty() {
                    pdist_args.push(Value::chars("euclid"));
                }
                if link_args.is_empty() {
                    link_args.push(Value::chars("single"));
                }
            }
            let mut largs = vec![Value::Num(x.clone()), link_args[0].clone(), Value::Cell(CellArray::row(pdist_args))];
            let z = linkage(it, &largs)?.unwrap();
            largs.clear();
            let Value::Num(zm) = &z[0] else { return rt("internal error") };
            let n = zm.dims[0];
            let zr: Vec<[f64; 3]> = (0..n).map(|i| [zm.data[i], zm.data[n + i], zm.data[2 * n + i]]).collect();
            cluster(&zr, &rest)
        }
    }
}
