//! Linear regression helpers: `regress`, `dummyvar`, `x2fx`.

use super::dist::{fcdf, tinv};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["regress", "dummyvar", "x2fx"];

const MTIMES_MSG: &str = "Incorrect dimensions for matrix multiplication. Check that the number of columns in the first matrix matches the number of rows in the second matrix. To operate on each element of the matrix individually, use TIMES (.*) for elementwise multiplication.";

type Mat = Vec<Vec<f64>>;

/// Householder QR with column pivoting (economy size): returns Q (n x k), R (k x k) with
/// k = min(n, p), and the column permutation.
fn qr_pivot(x: &Mat, n: usize, p: usize) -> (Mat, Mat, Vec<usize>) {
    let k = n.min(p);
    let mut a: Mat = x.clone();
    let mut perm: Vec<usize> = (0..p).collect();
    let mut vs: Vec<Vec<f64>> = Vec::new();
    for j in 0..k {
        // pivot: the remaining column with the largest norm (first on ties)
        let mut best = j;
        let mut bn = -1.0;
        for c in j..p {
            let nn: f64 = (j..n).map(|i| a[i][c] * a[i][c]).sum();
            if nn > bn * (1.0 + 1e-14) {
                bn = nn;
                best = c;
            }
        }
        if best != j {
            for row in a.iter_mut() {
                row.swap(j, best);
            }
            perm.swap(j, best);
        }
        let mut v: Vec<f64> = (j..n).map(|i| a[i][j]).collect();
        let norm = v.iter().map(|t| t * t).sum::<f64>().sqrt();
        if norm == 0.0 {
            vs.push(vec![0.0; n - j]);
            continue;
        }
        let alpha = if v[0] > 0.0 { -norm } else { norm };
        v[0] -= alpha;
        let vn: f64 = v.iter().map(|t| t * t).sum::<f64>();
        if vn > 0.0 {
            for c in j..p {
                let dot: f64 = (j..n).map(|i| v[i - j] * a[i][c]).sum();
                let f = 2.0 * dot / vn;
                for i in j..n {
                    a[i][c] -= f * v[i - j];
                }
            }
        }
        vs.push(v);
    }
    let r: Mat = (0..k).map(|i| (0..p).map(|c| if c >= i { a[i][c] } else { 0.0 }).collect()).collect();
    // Q = H1 H2 ... Hk applied to the first k columns of the identity
    let mut q: Mat = (0..n).map(|i| (0..k).map(|c| if i == c { 1.0 } else { 0.0 }).collect()).collect();
    for j in (0..k).rev() {
        let v = &vs[j];
        let vn: f64 = v.iter().map(|t| t * t).sum();
        if vn == 0.0 {
            continue;
        }
        for c in 0..k {
            let dot: f64 = (j..n).map(|i| v[i - j] * q[i][c]).sum();
            let f = 2.0 * dot / vn;
            for i in j..n {
                q[i][c] -= f * v[i - j];
            }
        }
    }
    (q, r, perm)
}

fn solve_upper(r: &Mat, b: &[f64], k: usize) -> Vec<f64> {
    let mut x = vec![0.0; k];
    for i in (0..k).rev() {
        let mut s = b[i];
        for j in i + 1..k {
            s -= r[i][j] * x[j];
        }
        x[i] = s / r[i][i];
    }
    x
}

fn regress(a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("REGRESS requires at least two input arguments.");
    }
    let (Some(Value::Num(y)), Some(Value::Num(x))) = (a.first(), a.get(1)) else { return rt("Y and X must be numeric.") };
    let alpha = match a.get(2) {
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        _ => 0.05,
    };
    let (n, ncol) = (x.dims[0], x.dims[1]);
    let is_vec = y.dims.len() == 2 && (y.dims[0] == 1 || y.dims[1] == 1);
    if !is_vec || y.numel() != n {
        return rt("Y must be a vector and must have the same number of rows as X.");
    }
    if y.dims[0] == 1 && y.dims[1] > 1 {
        return rt(MTIMES_MSG);
    }
    let single = is_single(x) || is_single(y);
    let wasnan: Vec<bool> = (0..n).map(|i| y.data[i].is_nan() || (0..ncol).any(|c| x.data[c * n + i].is_nan())).collect();
    let keep: Vec<usize> = (0..n).filter(|&i| !wasnan[i]).collect();
    let nn = keep.len();
    let yv: Vec<f64> = keep.iter().map(|&i| y.data[i]).collect();
    let xr: Mat = keep.iter().map(|&i| (0..ncol).map(|c| x.data[c * n + i]).collect()).collect();
    let (q, r, perm) = qr_pivot(&xr, nn, ncol);
    let kdim = nn.min(ncol);
    let p = if kdim == 0 {
        0
    } else if ncol == 1 || kdim == 1 {
        if r[0][0].abs() > 0.0 { 1 } else { 0 }
    } else {
        let tol = nn.max(ncol) as f64 * eps_of(r[0][0]);
        (0..kdim).filter(|&i| r[i][i].abs() > tol).count()
    };
    let rq: Mat = r.iter().take(p).map(|row| row[..p].to_vec()).collect();
    let qk: Mat = q.iter().map(|row| row[..p].to_vec()).collect();
    let qty: Vec<f64> = (0..p).map(|c| (0..nn).map(|i| qk[i][c] * yv[i]).sum()).collect();
    let bp = solve_upper(&rq, &qty, p);
    let mut b = vec![0.0; ncol];
    for (k, &pc) in perm.iter().take(p).enumerate() {
        b[pc] = bp[k];
    }
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(vec![ncol, 1], b.clone()), single))];
    if nargout < 2 {
        return Ok(Some(outs));
    }
    // inverse of R
    let mut ri = vec![vec![0.0; p]; p];
    for j in 0..p {
        let mut e = vec![0.0; p];
        e[j] = 1.0;
        let col = solve_upper(&rq, &e, p);
        for i in 0..p {
            ri[i][j] = col[i];
        }
    }
    let nu = (nn as f64 - p as f64).max(0.0);
    let yhat: Vec<f64> = xr.iter().map(|row| row.iter().zip(&b).map(|(u, v)| u * v).sum()).collect();
    let res: Vec<f64> = yv.iter().zip(&yhat).map(|(u, v)| u - v).collect();
    let normr = res.iter().map(|t| t * t).sum::<f64>().sqrt();
    let (rmse, mut tval) = if nu != 0.0 { (normr / nu.sqrt(), tinv(1.0 - alpha / 2.0, nu)) } else { (f64::NAN, 0.0) };
    let s2 = rmse * rmse;
    let mut se = vec![0.0; ncol];
    for k in 0..p {
        se[perm[k]] = rmse * (0..p).map(|j| ri[k][j] * ri[k][j]).sum::<f64>().sqrt();
    }
    let mut bint = vec![0.0; 2 * ncol];
    for i in 0..ncol {
        bint[i] = b[i] - tval * se[i];
        bint[ncol + i] = b[i] + tval * se[i];
    }
    outs.push(Value::Num(finish(Matrix::from_dims(vec![ncol, 2], bint), single)));
    let full = |v: &[f64]| -> Vec<f64> {
        let mut t = vec![f64::NAN; n];
        for (k, &i) in keep.iter().enumerate() {
            t[i] = v[k];
        }
        t
    };
    if nargout >= 3 {
        outs.push(Value::Num(finish(Matrix::from_dims(vec![n, 1], full(&res)), single)));
    }
    if nargout >= 4 {
        let mut hat: Vec<f64> = qk.iter().map(|row| row.iter().map(|t| t * t).sum()).collect();
        let thr = if single { (f32::EPSILON as f64).sqrt() } else { f64::EPSILON.sqrt() };
        let ok: Vec<bool> = hat.iter().map(|&h| 1.0 - h > thr).collect();
        for (h, &o) in hat.iter_mut().zip(&ok) {
            if !o {
                *h = 1.0;
            }
        }
        let ser: Vec<f64>;
        if nu > 1.0 {
            let sigmai: Vec<f64> = (0..nn).map(|i| if ok[i] { (0f64).max(nu * s2 / (nu - 1.0) - res[i] * res[i] / ((nu - 1.0) * (1.0 - hat[i]))).sqrt() } else { 0.0 }).collect();
            ser = (0..nn).map(|i| if ok[i] { (1.0 - hat[i]).sqrt() * sigmai[i] } else { f64::INFINITY }).collect();
            tval = tinv(1.0 - alpha / 2.0, nu - 1.0);
        } else if nu == 1.0 {
            ser = (0..nn).map(|i| if ok[i] { (1.0 - hat[i]).sqrt() * rmse } else { f64::INFINITY }).collect();
        } else {
            ser = vec![rmse; nn];
        }
        let lo: Vec<f64> = (0..nn).map(|i| res[i] - tval * ser[i]).collect();
        let hi: Vec<f64> = (0..nn).map(|i| res[i] + tval * ser[i]).collect();
        let mut rint = full(&lo);
        rint.extend(full(&hi));
        outs.push(Value::Num(finish(Matrix::from_dims(vec![n, 2], rint), single)));
    }
    if nargout >= 5 {
        let ybar = yv.iter().sum::<f64>() / nn as f64;
        let sse = normr * normr;
        let rss: f64 = yhat.iter().map(|t| (t - ybar).powi(2)).sum();
        let tss: f64 = yv.iter().map(|t| (t - ybar).powi(2)).sum();
        let r2 = 1.0 - sse / tss;
        let f = if p > 1 { (rss / (p as f64 - 1.0)) / s2 } else { f64::NAN };
        let prob = fcdf(f, p as f64 - 1.0, nu, true);
        outs.push(Value::Num(finish(Matrix::from_dims(vec![1, 4], vec![r2, f, prob, s2]), single)));
    }
    Ok(Some(outs))
}

fn dummyvar(a: &[Value]) -> Ret {
    if a.len() > 1 {
        return rt("Too many input arguments.");
    }
    let (m, n, group): (usize, usize, Vec<f64>) = match a.first() {
        None => return rt("Not enough input arguments."),
        Some(Value::Cell(c)) => {
            let nc = c.data.len();
            let mut cols: Vec<Vec<f64>> = Vec::new();
            let mut m = 0;
            for (j, v) in c.data.iter().enumerate() {
                let Value::Num(g) = v else { return rt("Number of elements must not change. Use [] as one of the size inputs to automatically calculate the appropriate size for that dimension.") };
                if g.kind == Kind::Char {
                    return rt("Number of elements must not change. Use [] as one of the size inputs to automatically calculate the appropriate size for that dimension.");
                }
                // grp2idx of a numeric vector: rank among the sorted distinct values (NaN stays NaN)
                let mut vals: Vec<f64> = g.data.iter().copied().filter(|x| !x.is_nan()).collect();
                vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
                vals.dedup();
                let idx: Vec<f64> = g.data.iter().map(|x| if x.is_nan() { f64::NAN } else { (vals.iter().position(|u| u == x).unwrap() + 1) as f64 }).collect();
                if j == 0 {
                    m = idx.len();
                } else if idx.len() != m {
                    return rt("All grouping variables must have the same number of rows.");
                }
                cols.push(idx);
            }
            (m, nc, cols.into_iter().flatten().collect())
        }
        Some(Value::Num(g)) => {
            if g.dims.len() > 2 {
                return rt("GROUP must be a vector or matrix.");
            }
            if matches!(g.kind, Kind::Int(_)) {
                return rt("Integers can only be combined with integers of the same class, or scalar doubles.");
            }
            let (mut m, mut n) = (g.dims[0], g.dims[1]);
            if m == 1 {
                m = n;
                n = 1;
            }
            (m, n, g.data.clone())
        }
        Some(_) => return rt("GROUP must be a vector or matrix."),
    };
    if group.iter().any(|&v| v.fract() != 0.0 && !v.is_nan()) || group.iter().any(|&v| v <= 0.0) {
        return rt("Each element of GROUP must be a positive integer.");
    }
    let single = matches!(a.first(), Some(Value::Num(g)) if is_single(g));
    let mut maxg = vec![0.0f64; n];
    for j in 0..n {
        maxg[j] = group[j * m..(j + 1) * m].iter().cloned().filter(|v| !v.is_nan()).fold(0.0, f64::max);
    }
    let total: usize = maxg.iter().map(|&v| v as usize).sum();
    let mut d = vec![0.0; m * total];
    let mut start = 0;
    for j in 0..n {
        let w = maxg[j] as usize;
        for i in 0..m {
            let g = group[j * m + i];
            if g.is_nan() {
                for c in 0..w {
                    d[(start + c) * m + i] = f64::NAN;
                }
            } else {
                d[(start + g as usize - 1) * m + i] = 1.0;
            }
        }
        start += w;
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![m, total], d), single))]))
}

fn x2fx(a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(x0)) = a.first() else { return rt("Not enough input arguments.") };
    let (m, n) = (x0.dims[0], x0.dims[1]);
    let single = is_single(x0);
    let mut x: Mat = (0..m).map(|i| (0..n).map(|j| x0.data[j * m + i]).collect()).collect();
    let categ: Vec<usize> = match a.get(2) {
        Some(Value::Num(c)) => c.data.iter().map(|&v| v as usize).collect(),
        _ => vec![],
    };
    let catbool: Vec<bool> = (1..=n).map(|j| categ.contains(&j)).collect();
    // The model: a name or a matrix of exponents.
    let model: Mat = match a.get(1) {
        None => model_of_name("linear", n, &categ)?,
        Some(Value::Num(mm)) if mm.is_empty() => model_of_name("linear", n, &categ)?,
        Some(v) if v.is_text() => {
            let t = v.text()?;
            if t.is_empty() {
                model_of_name("linear", n, &categ)?
            } else {
                match model_of_name(&t, n, &categ) {
                    Ok(mm) => mm,
                    Err(_) => return rt(format!("Undefined function '{t}' for input arguments of type 'double'.")),
                }
            }
        }
        Some(Value::Num(mm)) => {
            if mm.dims[1] != n {
                return rt("A numeric second argument must have the same number of columns as the first argument.");
            }
            (0..mm.dims[0]).map(|i| (0..n).map(|j| mm.data[j * mm.dims[0] + i]).collect()).collect()
        }
        Some(_) => return rt("The model must be a name or a matrix."),
    };
    let mut vardf = vec![1.0f64; n];
    match a.get(3) {
        Some(Value::Num(cl)) if !cl.is_empty() => {
            if cl.numel() != categ.len() {
                return rt("Unable to perform assignment because the left and right sides have a different number of elements.");
            }
            for (k, &cj) in categ.iter().enumerate() {
                let lv = *cl.data.get(k).unwrap_or(&1.0);
                vardf[cj - 1] = lv - 1.0;
            }
            for (k, &cj) in categ.iter().enumerate() {
                let lv = *cl.data.get(k).unwrap_or(&1.0);
                if x.iter().any(|r| !(r[cj - 1] >= 1.0 && r[cj - 1] <= lv && r[cj - 1].fract() == 0.0)) {
                    return rt(format!("Values of categorical variable {cj} must be integers between 1 and {lv}."));
                }
            }
        }
        _ => {
            for &cj in &categ {
                let mut u: Vec<f64> = x.iter().map(|r| r[cj - 1]).collect();
                u.sort_by(|p, q| p.partial_cmp(q).unwrap());
                u.dedup();
                vardf[cj - 1] = u.len() as f64 - 1.0;
                for r in x.iter_mut() {
                    r[cj - 1] = (u.iter().position(|t| *t == r[cj - 1]).unwrap() + 1) as f64;
                }
            }
        }
    }
    let rows = model.len();
    if model.iter().any(|r| r.iter().enumerate().any(|(j, &v)| catbool[j] && v > 1.0)) {
        return rt("The model matrix cannot request a power above 1 for a categorical variable.");
    }
    let termdf: Vec<usize> = model
        .iter()
        .map(|r| r.iter().enumerate().map(|(j, &v)| if v > 0.0 { vardf[j].max(1.0) } else { 1.0 }).product::<f64>() as usize)
        .collect();
    let mut termend = Vec::new();
    let mut acc = 0;
    for &t in &termdf {
        acc += t;
        termend.push(acc);
    }
    let termstart: Vec<usize> = termend.iter().zip(&termdf).map(|(e, d)| e - d + 1).collect();
    let total = acc;
    let mut d = vec![0.0; m * total];
    for idx in 0..rows {
        let pw = &model[idx];
        let t: Vec<bool> = pw.iter().map(|&v| v > 0.0).collect();
        let mut c: Vec<f64> = vec![1.0; m];
        if t.iter().any(|&b| b) {
            if pw.iter().enumerate().any(|(j, &v)| !catbool[j] && v != 0.0) {
                c = vec![1.0; m];
                for (j, &e) in pw.iter().enumerate() {
                    if !catbool[j] && e > 0.0 {
                        for i in 0..m {
                            c[i] *= if e == 1.0 { x[i][j] } else { x[i][j].powf(e) };
                        }
                    }
                }
            }
        }
        let cats: Vec<usize> = (0..n).filter(|&j| pw[j] > 0.0 && catbool[j]).collect();
        let cols_n = termdf[idx];
        let mut z: Option<Vec<f64>> = None;
        if !cats.is_empty() {
            let mut zz = vec![0.0; m * cols_n];
            let mut keep = vec![true; m];
            let mut colnum: Vec<f64> = (0..m).map(|i| x[i][cats[0]]).collect();
            let mut cumdf = 1.0;
            for i in 0..m {
                keep[i] = x[i][cats[0]] <= vardf[cats[0]];
            }
            for k in 1..cats.len() {
                cumdf *= vardf[cats[k - 1]];
                for i in 0..m {
                    keep[i] = keep[i] && x[i][cats[k]] <= vardf[cats[k]];
                    colnum[i] += cumdf * (x[i][cats[k]] - 1.0);
                }
            }
            for i in 0..m {
                if keep[i] {
                    zz[(colnum[i] as usize - 1) * m + i] = c[i];
                }
            }
            z = Some(zz);
        }
        let start = termstart[idx] - 1;
        match z {
            Some(zz) => {
                for k in 0..cols_n {
                    for i in 0..m {
                        d[(start + k) * m + i] = zz[k * m + i];
                    }
                }
            }
            None => {
                for i in 0..m {
                    d[start * m + i] = c[i];
                }
            }
        }
    }
    let mdl: Vec<f64> = (0..n).flat_map(|j| model.iter().map(move |r| r[j])).collect();
    let mut outs = vec![Value::Num(finish(Matrix::from_dims(vec![m, total], d), single))];
    if nargout > 1 {
        outs.push(Value::Num(Matrix::from_dims(vec![rows, n], mdl)));
    }
    if nargout > 2 {
        outs.push(Value::Num(Matrix::from_dims(vec![rows, 1], termstart.iter().map(|&v| v as f64).collect())));
    }
    if nargout > 3 {
        outs.push(Value::Num(Matrix::from_dims(vec![rows, 1], termend.iter().map(|&v| v as f64).collect())));
    }
    Ok(Some(outs))
}

fn model_of_name(name: &str, n: usize, categ: &[usize]) -> Result<Mat, crate::error::LangError> {
    let low = name.to_lowercase();
    let is = |full: &str| !low.is_empty() && full.starts_with(&low);
    let (hasint, hasquad) = if is("linear") || is("additive") {
        (false, false)
    } else if is("interactions") {
        (true, false)
    } else if is("quadratic") {
        (true, true)
    } else if is("purequadratic") {
        (false, true)
    } else {
        return rt("unknown model");
    };
    let mut m: Mat = vec![vec![0.0; n]];
    for i in 0..n {
        let mut r = vec![0.0; n];
        r[i] = 1.0;
        m.push(r);
    }
    if hasint && n > 1 {
        for c in 0..n {
            for r in c + 1..n {
                let mut row = vec![0.0; n];
                row[r] = 1.0;
                row[c] = 1.0;
                m.push(row);
            }
        }
    }
    if hasquad {
        for i in 0..n {
            if !categ.contains(&(i + 1)) {
                let mut r = vec![0.0; n];
                r[i] = 2.0;
                m.push(r);
            }
        }
    }
    Ok(m)
}

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "regress" => regress(a, nargout),
        "dummyvar" => dummyvar(a),
        _ => x2fx(a, nargout),
    }
}
