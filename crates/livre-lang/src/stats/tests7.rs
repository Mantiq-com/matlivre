//! `regstats` and `stepwisefit`.

use super::dist::{fcdf, tcdf};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::{Interpreter, OutputKind};
use crate::value::*;

pub const NAMES: &[&str] = &["regstats", "stepwisefit"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "regstats" => regstats(it, a, nargout),
        "stepwisefit" => stepwisefit(it, a, nargout),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn col_major(cols: &[Vec<f64>]) -> Vec<f64> {
    cols.iter().flat_map(|c| c.iter().copied()).collect()
}

fn mat_cols(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1..].iter().product::<usize>());
    (0..c).map(|j| m.data[j * r..(j + 1) * r].to_vec()).collect()
}

fn tcdf_f(x: f64, df: f64) -> f64 {
    tcdf(x, df)
}

fn tcell(mat: Matrix) -> Value {
    Value::Num(mat)
}

// ---------------------------------------------------------------------------------------------
// regstats

const REG_NAMES: [&str; 23] = [
    "Q", "R", "beta", "covb", "yhat", "r", "mse", "rsquare", "adjrsquare", "leverage", "hatmat", "s2_i", "beta_i", "standres", "studres", "dfbetas", "dffit", "dffits", "covratio", "cookd", "tstat", "fstat", "dwstat",
];

fn regstats(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("At least two arguments are required.");
    }
    let model = a.get(2).cloned().filter(|v| !is_empty_val(v)).unwrap_or_else(|| Value::chars("linear"));
    let which_v = a.get(3).cloned().filter(|v| !is_empty_val(v) && !matches!(v, Value::Cell(c) if c.data.is_empty())).unwrap_or_else(|| Value::chars("all"));
    let (Value::Num(resp), Value::Num(data)) = (&a[0], &a[1]) else { return rt("RESPONSES and DATA must be numeric.") };
    if nargout == 0 && a.len() < 4 {
        return rt("The interactive form of REGSTATS is not available.");
    }
    let mut y: Vec<f64> = resp.data.clone();
    let (yr, yc) = (resp.dims[0], resp.dims[1]);
    let yrows = if yr == 1 { yc } else { yr };
    let ycols = if yr == 1 { 1 } else { yc };
    let (xr, xc) = (data.dims[0], data.dims[1]);
    let (xrows, xcols, mut dcols) = if xr == 1 { (xc, 1, vec![data.data.clone()]) } else { (xr, xc, mat_cols(data)) };
    if yrows != xrows {
        return rt("RESPONSES and DATA must have the same number of rows.");
    }
    if ycols > 1 {
        return rt("RESPONSES must have a single column.");
    }
    let _ = xcols;
    let n_all = yrows;
    let wasnan: Vec<bool> = (0..n_all).map(|i| y[i].is_nan() || dcols.iter().any(|c| c[i].is_nan())).collect();
    let havenans = wasnan.iter().any(|&b| b);
    if havenans {
        y = y.iter().zip(&wasnan).filter(|(_, &w)| !w).map(|(v, _)| *v).collect();
        for c in dcols.iter_mut() {
            *c = c.iter().zip(&wasnan).filter(|(_, &w)| !w).map(|(v, _)| *v).collect();
        }
    }
    let nobs = y.len();
    let dm = Matrix::from_dims(vec![nobs, dcols.len()], col_major(&dcols));
    let xv = it.call_function("x2fx", vec![Value::Num(dm), model], 1)?;
    let Some(Value::Num(xm)) = xv.first() else { return rt("x2fx returned no value") };
    let p = xm.dims[1];
    if xm.dims[0] < p {
        return rt("The design matrix has more predictor variables than observations.");
    }
    // which statistics
    let names: Vec<String> = match &which_v {
        Value::Cell(c) => {
            if !c.data.iter().all(|e| e.is_text()) {
                return rt("WHICHSTATS argument must be one or more statistic names.");
            }
            c.data.iter().map(|e| e.text().unwrap_or_default()).collect()
        }
        v if v.is_text() => vec![v.text().unwrap_or_default()],
        _ => return rt("WHICHSTATS argument must be one or more statistic names."),
    };
    let mut idx = [false; 23];
    for nm in &names {
        if nm == "all" {
            idx = [true; 23];
            break;
        }
        match REG_NAMES.iter().position(|v| v == nm) {
            Some(k) => idx[k] = true,
            None => return rt(format!("Invalid statistic name '{nm}'.")),
        }
    }
    if !idx.iter().any(|&b| b) {
        return Ok(Some(vec![]));
    }
    let qr = it.call_function("qr", vec![Value::Num(xm.clone()), Value::num(0.0)], 2)?;
    let (Some(Value::Num(qm)), Some(Value::Num(rm))) = (qr.first(), qr.get(1)) else { return rt("qr returned no value") };
    let qcols = mat_cols(qm);
    let n = nobs;
    // R is upper triangular p x p, stored column-major
    let rr = |i: usize, j: usize| rm.data[j * p + i];
    let qty: Vec<f64> = qcols.iter().map(|q| dot(q, &y)).collect();
    let mut beta = vec![0.0; p];
    for i in (0..p).rev() {
        let mut s = qty[i];
        for j in i + 1..p {
            s -= rr(i, j) * beta[j];
        }
        beta[i] = s / rr(i, i);
    }
    let xcols = mat_cols(xm);
    let yhat: Vec<f64> = (0..n).map(|i| (0..p).map(|j| xcols[j][i] * beta[j]).sum()).collect();
    let resid: Vec<f64> = (0..n).map(|i| y[i] - yhat[i]).collect();
    let dfe = (n - p) as f64;
    let dft = n as f64 - 1.0;
    let ybar = y.iter().sum::<f64>() / n as f64;
    let sse = dot(&resid, &resid);
    let ssr: f64 = yhat.iter().map(|v| (v - ybar).powi(2)).sum();
    let sst: f64 = y.iter().map(|v| (v - ybar).powi(2)).sum();
    let mse = sse / dfe;
    let h: Vec<f64> = (0..n).map(|i| qcols.iter().map(|q| q[i] * q[i]).sum()).collect();
    let s2i: Vec<f64> = (0..n).map(|i| (dfe * mse - resid[i] * resid[i] / (1.0 - h[i])) / (dfe - 1.0)).collect();
    let ei: Vec<f64> = (0..n).map(|i| resid[i] / (s2i[i] * (1.0 - h[i])).sqrt()).collect();
    // inverse of R (upper triangular)
    let mut ri = vec![vec![0.0; p]; p];
    for c in 0..p {
        for i in (0..p).rev() {
            let mut s = if i == c { 1.0 } else { 0.0 };
            for j in i + 1..p {
                s -= rr(i, j) * ri[j][c];
            }
            ri[i][c] = s / rr(i, i);
        }
    }
    let xtxi: Vec<Vec<f64>> = (0..p).map(|i| (0..p).map(|j| (0..p).map(|k| ri[i][k] * ri[j][k]).sum()).collect()).collect();
    let covb: Vec<Vec<f64>> = xtxi.iter().map(|r| r.iter().map(|v| v * mse).collect()).collect();
    let mut b_i: Vec<Vec<f64>> = vec![];
    if idx[12] || idx[15] {
        for i in 0..n {
            let stde = resid[i] / (1.0 - h[i]);
            let col: Vec<f64> = (0..p).map(|r| beta[r] - (0..p).map(|k| ri[r][k] * qcols[k][i] * stde).sum::<f64>()).collect();
            b_i.push(col);
        }
    }
    let fix = |v: Vec<f64>| -> Vec<f64> {
        if !havenans {
            return v;
        }
        let mut out = vec![f64::NAN; n_all];
        let mut k = 0;
        for i in 0..n_all {
            if !wasnan[i] {
                out[i] = v[k];
                k += 1;
            }
        }
        out
    };
    let colv = |v: Vec<f64>| Value::Num(Matrix::col(v));
    let mut st = StructArray::scalar();
    st.set("source", Value::chars("regstats"));
    if idx[0] {
        st.set("Q", Value::Num(qm.clone()));
    }
    if idx[1] {
        st.set("R", Value::Num(rm.clone()));
    }
    if idx[2] {
        st.set("beta", colv(beta.clone()));
    }
    if idx[3] {
        st.set("covb", Value::Num(Matrix::from_dims(vec![p, p], (0..p).flat_map(|j| (0..p).map(move |i| (i, j))).map(|(i, j)| covb[i][j]).collect())));
    }
    if idx[4] {
        st.set("yhat", colv(fix(yhat.clone())));
    }
    if idx[5] {
        st.set("r", colv(fix(resid.clone())));
    }
    if idx[6] {
        st.set("mse", Value::num(mse));
    }
    if idx[7] {
        st.set("rsquare", Value::num(1.0 - sse / sst));
    }
    if idx[8] {
        st.set("adjrsquare", Value::num(1.0 - (sse / sst) * (dft / dfe)));
    }
    if idx[9] {
        st.set("leverage", colv(fix(h.clone())));
    }
    if idx[10] {
        let mut hm = vec![0.0; n_all * n_all];
        let kept: Vec<usize> = (0..n_all).filter(|&i| !wasnan[i]).collect();
        if havenans {
            for &i in &kept {
                for &j in &kept {
                    let (ii, jj) = (kept.iter().position(|&k| k == i).unwrap(), kept.iter().position(|&k| k == j).unwrap());
                    hm[j * n_all + i] = qcols.iter().map(|q| q[ii] * q[jj]).sum();
                }
            }
            for (i, &w) in wasnan.iter().enumerate() {
                if w {
                    hm[i * n_all + i] = f64::NAN;
                }
            }
        } else {
            for i in 0..n {
                for j in 0..n {
                    hm[j * n + i] = qcols.iter().map(|q| q[i] * q[j]).sum();
                }
            }
        }
        st.set("hatmat", Value::Num(Matrix::from_dims(vec![n_all, n_all], hm)));
    }
    if idx[11] {
        st.set("s2_i", colv(fix(s2i.clone())));
    }
    if idx[12] {
        let mut out = vec![0.0; p * n_all];
        let mut k = 0;
        for i in 0..n_all {
            for r in 0..p {
                out[i * p + r] = if havenans && wasnan[i] { beta[r] } else { b_i[k][r] };
            }
            if !(havenans && wasnan[i]) {
                k += 1;
            }
        }
        st.set("beta_i", Value::Num(Matrix::from_dims(vec![p, n_all], out)));
    }
    if idx[13] {
        let v: Vec<f64> = (0..n).map(|i| resid[i] / (mse * (1.0 - h[i])).sqrt()).collect();
        st.set("standres", colv(fix(v)));
    }
    if idx[14] {
        st.set("studres", colv(fix(ei.clone())));
    }
    if idx[15] {
        let rtri: Vec<f64> = (0..p).map(|i| xtxi[i][i].sqrt()).collect();
        let mut out = vec![0.0; p * n_all];
        let mut k = 0;
        for i in 0..n_all {
            if havenans && wasnan[i] {
                continue;
            }
            for r in 0..p {
                out[i * p + r] = (beta[r] - b_i[k][r]) / (s2i[k].sqrt() * rtri[r]);
            }
            k += 1;
        }
        st.set("dfbetas", Value::Num(Matrix::from_dims(vec![p, n_all], out)));
    }
    if idx[16] {
        let v: Vec<f64> = (0..n).map(|i| h[i] * resid[i] / (1.0 - h[i])).collect();
        st.set("dffit", colv(fix(v)));
    }
    if idx[17] {
        let v: Vec<f64> = (0..n).map(|i| (h[i] / (1.0 - h[i])).sqrt() * ei[i]).collect();
        st.set("dffits", colv(fix(v)));
    }
    if idx[18] {
        let pf = p as f64;
        let v: Vec<f64> = (0..n).map(|i| 1.0 / ((((n as f64 - pf - 1.0 + ei[i] * ei[i]) / (n as f64 - pf)).powf(pf)) * (1.0 - h[i]))).collect();
        st.set("covratio", colv(fix(v)));
    }
    if idx[19] {
        let v: Vec<f64> = (0..n).map(|i| resid[i] * resid[i] * (h[i] / ((1.0 - h[i]) * (1.0 - h[i]))) / (p as f64 * mse)).collect();
        st.set("cookd", colv(fix(v)));
    }
    if idx[20] {
        let mut d = StructArray::scalar();
        let se: Vec<f64> = (0..p).map(|i| covb[i][i].sqrt()).collect();
        let t: Vec<f64> = (0..p).map(|i| beta[i] / se[i]).collect();
        d.set("beta", colv(beta.clone()));
        d.set("se", colv(se));
        d.set("t", colv(t.clone()));
        d.set("pval", colv(t.iter().map(|v| 2.0 * tcdf_f(-v.abs(), dfe)).collect()));
        d.set("dfe", Value::num(dfe));
        st.set("tstat", Value::Struct(d));
    }
    if idx[21] {
        let mut d = StructArray::scalar();
        let dfr = p as f64 - 1.0;
        let f = (ssr / dfr) / (sse / dfe);
        d.set("sse", Value::num(sse));
        d.set("dfe", Value::num(dfe));
        d.set("dfr", Value::num(dfr));
        d.set("ssr", Value::num(ssr));
        d.set("f", Value::num(f));
        d.set("pval", Value::num(fcdf(f, dfr, dfe, true)));
        st.set("fstat", Value::Struct(d));
    }
    if idx[22] {
        let mut d = StructArray::scalar();
        let r = it.call_function("dwtest", vec![colv(resid.clone()), Value::Num(xm.clone())], 2)?;
        d.set("dw", r[1].clone());
        d.set("pval", r[0].clone());
        st.set("dwstat", Value::Struct(d));
    }
    let _ = tcell;
    Ok(Some(vec![Value::Struct(st)]))
}

// ---------------------------------------------------------------------------------------------
// stepwisefit

/// Column-pivoted Gram-Schmidt: the orthonormal columns, the pivoted R (rank x rank, row-major
/// rows), the pivot order and the rank under the tolerance `tol * |R(1,1)|`.
fn pivoted_qr(cols: &[Vec<f64>], tol: f64) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<usize>) {
    let n = cols.first().map(|c| c.len()).unwrap_or(0);
    let p = cols.len();
    let mut a: Vec<Vec<f64>> = cols.to_vec();
    let mut perm: Vec<usize> = (0..p).collect();
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut r = vec![vec![0.0; p]; p];
    let mut r11 = 0.0;
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
        if k == 0 {
            r11 = bn;
        }
        if !(bn > tol * r11) || bn == 0.0 {
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
    }
    let rank = q.len();
    let rr: Vec<Vec<f64>> = (0..rank).map(|i| r[i][..rank].to_vec()).collect();
    (q, rr, perm[..rank].to_vec())
}

fn safedivide(num: f64, den: f64) -> f64 {
    if den == 0.0 {
        let s = if num > 0.0 { 1.0 } else if num < 0.0 { -1.0 } else { 0.0 };
        f64::INFINITY * s
    } else {
        num / den
    }
}

struct StepCalc {
    b: Vec<f64>,
    se: Vec<f64>,
    pval: Vec<f64>,
    tstat: Vec<f64>,
    dfe: f64,
    df0: f64,
    sstotal: f64,
    ssresid: f64,
    fstat: f64,
    pval_model: f64,
    rmse: f64,
    xr: Vec<Vec<f64>>,
    yr: Vec<f64>,
    covb: Vec<Vec<f64>>,
    intercept: f64,
}

fn stepcalc(allx: &[Vec<f64>], y: &[f64], inn: &[bool]) -> StepCalc {
    let n = y.len();
    let pp = inn.len();
    let mut cols: Vec<Vec<f64>> = vec![vec![1.0; n]];
    for (j, c) in allx.iter().enumerate() {
        if inn[j] {
            cols.push(c.clone());
        }
    }
    let nin = cols.len();
    let tol = n.max(pp + 1) as f64 * f64::EPSILON;
    let outc: Vec<usize> = (0..pp).filter(|&j| !inn[j]).collect();
    let sumxsq: Vec<f64> = outc.iter().map(|&j| dot(&allx[j], &allx[j])).collect();
    let (q, r, perm) = pivoted_qr(&cols, tol);
    let rank = q.len();
    let mut qb: Vec<f64> = q.iter().map(|qk| dot(qk, y)).collect();
    let maxqb = qb.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    for v in qb.iter_mut() {
        if v.abs() < tol * maxqb {
            *v = 0.0;
        }
    }
    let mut sol = vec![0.0; rank];
    for i in (0..rank).rev() {
        let mut s = qb[i];
        for j in i + 1..rank {
            s -= r[i][j] * sol[j];
        }
        sol[i] = s / r[i][i];
    }
    let mut b = vec![0.0; nin];
    for (k, &pk) in perm.iter().enumerate() {
        b[pk] = sol[k];
    }
    let mut resid: Vec<f64> = (0..n).map(|i| y[i] - (0..nin).map(|j| cols[j][i] * b[j]).sum::<f64>()).collect();
    let dfe = n as f64 - rank as f64;
    let df0 = rank as f64 - 1.0;
    let ym = y.iter().sum::<f64>() / n as f64;
    let sstotal: f64 = y.iter().map(|v| (v - ym).powi(2)).sum();
    let mut ssresid = dot(&resid, &resid);
    let perfect = dfe == 0.0 || ssresid < tol * sstotal;
    if perfect {
        ssresid = 0.0;
        resid.iter_mut().for_each(|v| *v = 0.0);
    }
    let rmse = safedivide(ssresid, dfe).sqrt();
    // inverse of R
    let mut rinv = vec![vec![0.0; rank]; rank];
    for c in 0..rank {
        for i in (0..rank).rev() {
            let mut s = if i == c { 1.0 } else { 0.0 };
            for j in i + 1..rank {
                s -= r[i][j] * rinv[j][c];
            }
            rinv[i][c] = s / r[i][i];
        }
    }
    let mut covb_full = vec![vec![0.0; nin]; nin];
    for i in 0..rank {
        for j in 0..rank {
            covb_full[perm[i]][perm[j]] = rmse * rmse * (0..rank).map(|k| rinv[i][k] * rinv[j][k]).sum::<f64>();
        }
    }
    let se_full: Vec<f64> = (0..nin).map(|i| covb_full[i][i].sqrt()).collect();
    // remove the effect of the included predictors from the others
    let mut xr: Vec<Vec<f64>> = outc.iter().map(|&j| allx[j].clone()).collect();
    for c in xr.iter_mut() {
        let proj: Vec<f64> = q.iter().map(|qk| dot(qk, c)).collect();
        for (k, qk) in q.iter().enumerate() {
            for i in 0..n {
                c[i] -= proj[k] * qk[i];
            }
        }
    }
    let yr = resid.clone();
    let mut xx: Vec<f64> = xr.iter().map(|c| dot(c, c)).collect();
    for (k, c) in xr.iter_mut().enumerate() {
        if xx[k] <= tol * sumxsq[k] {
            c.iter_mut().for_each(|v| *v = 0.0);
            xx[k] = 1.0;
        }
    }
    let b2: Vec<f64> = xr.iter().zip(&xx).map(|(c, &x)| safedivide(dot(&yr, c), x)).collect();
    let df2 = (dfe - 1.0).max(0.0);
    let s2: Vec<f64> = xr
        .iter()
        .enumerate()
        .map(|(k, c)| {
            let r2: f64 = (0..n).map(|i| (yr[i] - c[i] * b2[k]).powi(2)).sum();
            safedivide(safedivide(r2, df2).sqrt(), xx[k].sqrt())
        })
        .collect();
    let mut bb = vec![0.0; pp];
    let mut se = vec![0.0; pp];
    let mut covb = vec![vec![0.0; pp]; pp];
    let mut k_in = 1;
    let inidx: Vec<usize> = (0..pp).filter(|&j| inn[j]).collect();
    for &j in &inidx {
        bb[j] = b[k_in];
        se[j] = se_full[k_in];
        k_in += 1;
    }
    for (k, &j) in outc.iter().enumerate() {
        bb[j] = b2[k];
        se[j] = s2[k];
    }
    for (ia, &ja) in inidx.iter().enumerate() {
        for (ib, &jb) in inidx.iter().enumerate() {
            covb[ja][jb] = covb_full[ia + 1][ib + 1];
        }
    }
    let mut pval = vec![0.0; pp];
    let mut tstat = vec![0.0; pp];
    for &j in &inidx {
        let t = safedivide(bb[j], se[j]);
        tstat[j] = t;
        pval[j] = 2.0 * tcdf_f(-t.abs(), dfe);
    }
    for &j in &outc {
        if dfe > 1.0 {
            let t = safedivide(bb[j], se[j]);
            tstat[j] = t;
            pval[j] = 2.0 * tcdf_f(-t.abs(), dfe - 1.0);
        } else {
            tstat[j] = f64::NAN;
            pval[j] = f64::NAN;
        }
    }
    let msexpl = safedivide(sstotal - ssresid, df0);
    let fstat = safedivide(msexpl, rmse * rmse);
    let pval_model = fcdf(fstat, df0, dfe, true);
    StepCalc { b: bb, se, pval, tstat, dfe, df0, sstotal, ssresid, fstat, pval_model, rmse, xr, yr, covb, intercept: b[0] }
}

fn vec_mat(v: &[f64]) -> Value {
    Value::Num(Matrix::col(v.to_vec()))
}

fn stepwisefit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    let rest = &a[2..];
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let names = ["inmodel", "penter", "premove", "display", "maxiter", "keep", "scale"];
    let mut vals: Vec<Option<&Value>> = vec![None; names.len()];
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match pick {
            Some(i) => vals[i] = Some(&pair[1]),
            None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            None => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let (Value::Num(xm), Value::Num(ym)) = (&a[0], &a[1]) else { return rt("Both X and Y must be real-valued.") };
    let pcount = xm.dims[1];
    let bad_in = "The 'InModel' parameter value must be a logical vector with one value for each column of X, or a list of column numbers of X.";
    let bad_keep = "The 'Keep' parameter value must be a logical vector with one value for each column of X, or a list of column numbers of X.";
    let to_mask = |v: Option<&Value>, msg: &str| -> Result<Vec<bool>, crate::error::LangError> {
        match v {
            None => Ok(vec![false; pcount]),
            Some(v) if is_empty_val(v) => Ok(vec![false; pcount]),
            Some(Value::Num(m)) if m.kind == Kind::Logical => {
                if m.data.len() != pcount {
                    return rt(msg);
                }
                Ok(m.data.iter().map(|&d| d != 0.0).collect())
            }
            Some(Value::Num(m)) => {
                if m.data.iter().any(|&d| !(d >= 1.0 && d <= pcount as f64 && d == d.floor())) {
                    return rt(msg);
                }
                let mut mask = vec![false; pcount];
                for &d in &m.data {
                    mask[d as usize - 1] = true;
                }
                Ok(mask)
            }
            Some(_) => rt(msg),
        }
    };
    let mut inn = to_mask(vals[0], bad_in)?;
    let keep = to_mask(vals[5], bad_keep)?;
    let get_p = |v: Option<&Value>| -> f64 {
        match v {
            Some(Value::Num(m)) if !m.is_empty() && m.kind != Kind::Char && m.data.len() == 1 => m.data[0],
            _ => f64::NAN,
        }
    };
    let pe_given = vals[1].map(|v| !is_empty_val(v)).unwrap_or(false);
    let pr_given = vals[2].map(|v| !is_empty_val(v)).unwrap_or(false);
    let (mut pe, mut pr) = (get_p(vals[1]), get_p(vals[2]));
    if !pe_given && !pr_given {
        pe = 0.05;
        pr = 0.10;
    } else if !pe_given {
        pe = pr.min(0.05);
    } else if !pr_given {
        pr = pe.max(0.10);
    }
    if !(pe > 0.0 && pe < 1.0) {
        return rt("The 'PEnter' parameter value must be a scalar strictly between zero and one.");
    }
    if !(pr > 0.0 && pr < 1.0) {
        return rt("The 'PRemove' parameter value must be a scalar strictly between zero and one.");
    }
    if pe > pr {
        return rt("'PEnter' must be less than or equal to 'PRemove'.");
    }
    if ym.dims[1] != 1 {
        return rt("Y must be a column vector.");
    }
    if ym.dims[0] != xm.dims[0] {
        return rt("X and Y must have the same number of rows.");
    }
    if xm.im.is_some() || ym.im.is_some() {
        return rt("Both X and Y must be real-valued.");
    }
    let n_all = ym.dims[0];
    let mut xcols = mat_cols(xm);
    let mut y = ym.data.clone();
    let wasnan: Vec<bool> = (0..n_all).map(|i| y[i].is_nan() || xcols.iter().any(|c| c[i].is_nan())).collect();
    if wasnan.iter().any(|&w| w) {
        y = y.iter().zip(&wasnan).filter(|(_, &w)| !w).map(|(v, _)| *v).collect();
        for c in xcols.iter_mut() {
            *c = c.iter().zip(&wasnan).filter(|(_, &w)| !w).map(|(v, _)| *v).collect();
        }
    }
    let doscale = matches!(vals[6], Some(v) if v.is_text() && v.text().unwrap_or_default() == "on");
    let dodisp = match vals[3] {
        None => true,
        Some(v) => v.is_text() && v.text().unwrap_or_default() == "on",
    };
    let maxiter = match vals[4] {
        None => f64::INFINITY,
        Some(v) if is_empty_val(v) => f64::INFINITY,
        Some(Value::Num(m)) if m.data.len() == 1 => m.data[0],
        _ => f64::INFINITY,
    };
    let nobs = y.len();
    let sx: Vec<f64> = xcols
        .iter()
        .map(|c| {
            let m = c.iter().sum::<f64>() / nobs as f64;
            let s = (c.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (nobs as f64 - 1.0)).sqrt();
            if s == 0.0 { 1.0 } else { s }
        })
        .collect();
    let scaled: Vec<Vec<f64>> = xcols.iter().zip(&sx).map(|(c, &s)| c.iter().map(|v| v / s).collect()).collect();
    let cols_text = |inn: &[bool]| -> String {
        if !inn.iter().any(|&b| b) {
            "none".to_string()
        } else {
            inn.iter().enumerate().filter(|(_, &b)| b).map(|(j, _)| format!("{} ", j + 1)).collect()
        }
    };
    if dodisp {
        it.emit(OutputKind::Text(format!("Initial columns included: {}\n", cols_text(&inn))));
    }
    let mut hist_b: Vec<Vec<f64>> = vec![];
    let mut hist_rmse: Vec<f64> = vec![];
    let mut hist_df0: Vec<f64> = vec![];
    let mut hist_in: Vec<Vec<bool>> = vec![];
    let mut jstep = 0.0;
    let mut sc;
    let mut nextstep;
    loop {
        sc = stepcalc(&scaled, &y, &inn);
        let mut b = sc.b.clone();
        let mut se = sc.se.clone();
        if !doscale {
            for j in 0..pcount {
                b[j] /= sx[j];
                se[j] /= sx[j];
            }
            sc.b = b.clone();
            sc.se = se.clone();
            let outc: Vec<usize> = (0..pcount).filter(|&j| !inn[j]).collect();
            for (k, c) in sc.xr.iter_mut().enumerate() {
                c.iter_mut().for_each(|v| *v *= sx[outc[k]]);
            }
            for i in 0..pcount {
                for j in 0..pcount {
                    sc.covb[i][j] /= sx[i] * sx[j];
                }
            }
        }
        // the next move
        let mut swap = 0usize;
        let mut pinout = f64::NAN;
        let termsout: Vec<usize> = (0..pcount).filter(|&j| !inn[j] && !keep[j]).collect();
        if !termsout.is_empty() {
            let mut kmin = termsout[0];
            let mut pmin = sc.pval[kmin];
            for &t in &termsout {
                if sc.pval[t] < pmin || pmin.is_nan() && !sc.pval[t].is_nan() {
                    pmin = sc.pval[t];
                    kmin = t;
                }
            }
            if pmin < pe {
                swap = kmin + 1;
                pinout = pmin;
            }
        }
        if swap == 0 {
            let termsin: Vec<usize> = (0..pcount).filter(|&j| inn[j] && !keep[j]).collect();
            if !termsin.is_empty() {
                let bad: Vec<usize> = termsin.iter().copied().filter(|&j| sc.pval[j].is_nan()).collect();
                if !bad.is_empty() {
                    if let Some(&j) = bad.iter().find(|&&j| b[j].is_nan()) {
                        swap = j + 1;
                    } else {
                        let mut bj = bad[0];
                        for &j in &bad {
                            if b[j].abs() < b[bj].abs() {
                                bj = j;
                            }
                        }
                        swap = bj + 1;
                    }
                    pinout = f64::NAN;
                } else {
                    let mut kmax = termsin[0];
                    let mut pmax = sc.pval[kmax];
                    for &t in &termsin {
                        if sc.pval[t] > pmax {
                            pmax = sc.pval[t];
                            kmax = t;
                        }
                    }
                    if pmax > pr {
                        swap = kmax + 1;
                        pinout = pmax;
                    }
                }
            }
        }
        nextstep = swap;
        if nargout >= 7 && jstep > 0.0 {
            let mut col = vec![0.0; pcount];
            for j in 0..pcount {
                if inn[j] {
                    col[j] = b[j];
                }
            }
            hist_b.push(col);
            hist_rmse.push(sc.rmse);
            hist_df0.push(sc.df0);
            hist_in.push(inn.clone());
        }
        if jstep >= maxiter {
            break;
        }
        jstep += 1.0;
        if nextstep == 0 {
            break;
        } else if dodisp {
            let word = if inn[nextstep - 1] { "removed" } else { "added" };
            it.emit(OutputKind::Text(format!("Step {}, {} column {}, p={}\n", jstep as i64, word, nextstep, fmt_g(pinout))));
        }
        inn[nextstep - 1] = !inn[nextstep - 1];
    }
    if dodisp {
        it.emit(OutputKind::Text(format!("Final columns included:  {}\n", cols_text(&inn))));
        let mut data: Vec<Value> = Vec::new();
        let rows = pcount + 1;
        data.resize(rows * 4, Value::empty());
        for (c, h) in ["Coeff", "Std.Err.", "Status", "P"].iter().enumerate() {
            data[c * rows] = Value::chars(h);
        }
        for j in 0..pcount {
            data[j + 1] = Value::num(sc.b[j]);
            data[rows + j + 1] = Value::num(sc.se[j]);
            data[2 * rows + j + 1] = Value::chars(if inn[j] { "In" } else { "Out" });
            data[3 * rows + j + 1] = Value::num(sc.pval[j]);
        }
        let cell = Value::Cell(CellArray::from_dims(vec![rows, 4], data));
        it.call_function("disp", vec![cell], 0)?;
    }
    let single_in = matches!(&a[0], Value::Num(m) if is_single(m)) || matches!(&a[1], Value::Num(m) if is_single(m));
    let vec_mat = |v: &[f64]| -> Value { Value::Num(finish(Matrix::col(v.to_vec()), single_in)) };
    let sn = |v: f64| -> Value { Value::Num(finish(Matrix::scalar(v), single_in)) };
    let mut out = vec![vec_mat(&sc.b), vec_mat(&sc.se), Value::Num(Matrix::col(sc.pval.clone()))];
    let mut inm = Matrix::from_dims(vec![1, pcount], inn.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect());
    inm.kind = Kind::Logical;
    out.push(Value::Num(inm));
    if nargout >= 5 {
        let mut st = StructArray::scalar();
        st.set("source", Value::chars("stepwisefit"));
        st.set("dfe", Value::num(sc.dfe));
        st.set("df0", Value::num(sc.df0));
        st.set("SStotal", sn(sc.sstotal));
        st.set("SSresid", sn(sc.ssresid));
        st.set("fstat", sn(sc.fstat));
        st.set("pval", sn(sc.pval_model));
        st.set("rmse", sn(sc.rmse));
        let nout = sc.xr.len();
        st.set("xr", Value::Num(finish(Matrix::from_dims(vec![nobs, nout], col_major(&sc.xr)), single_in)));
        st.set("yr", vec_mat(&sc.yr));
        st.set("B", vec_mat(&sc.b));
        st.set("SE", vec_mat(&sc.se));
        st.set("TSTAT", vec_mat(&sc.tstat));
        st.set("PVAL", vec_mat(&sc.pval));
        st.set("covb", Value::Num(finish(Matrix::from_dims(vec![pcount, pcount], (0..pcount).flat_map(|j| (0..pcount).map(move |i| (i, j))).map(|(i, j)| sc.covb[i][j]).collect()), single_in)));
        st.set("intercept", sn(sc.intercept));
        let mut wn = Matrix::from_dims(vec![n_all, 1], wasnan.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect());
        wn.kind = Kind::Logical;
        st.set("wasnan", Value::Num(wn));
        out.push(Value::Struct(st));
    }
    if nargout >= 6 {
        out.push(Value::num(nextstep as f64));
    }
    if nargout >= 7 {
        let mut h = StructArray::scalar();
        let nst = hist_b.len();
        h.set("B", Value::Num(Matrix::from_dims(vec![pcount, nst], col_major(&hist_b))));
        if nst == 0 {
            h.set("rmse", Value::empty());
            h.set("df0", Value::empty());
        } else {
            h.set("rmse", Value::Num(Matrix::from_dims(vec![1, nst], hist_rmse)));
            h.set("df0", Value::Num(Matrix::from_dims(vec![1, nst], hist_df0)));
        }
        let mut hin = Matrix::from_dims(vec![nst, pcount], (0..pcount).flat_map(|j| hist_in.iter().map(move |r| if r[j] { 1.0 } else { 0.0 })).collect());
        hin.kind = Kind::Logical;
        h.set("in", Value::Num(hin));
        out.push(Value::Struct(h));
    }
    Ok(Some(out))
}

/// `%g` formatting.
fn fmt_g(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    if v == 0.0 {
        return "0".into();
    }
    let exp = v.abs().log10().floor() as i32;
    if exp < -5 || exp >= 6 {
        let s = format!("{:.5e}", v);
        let (m, e) = s.split_once('e').unwrap();
        let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.').to_string() } else { m.to_string() };
        let e: i32 = e.parse().unwrap();
        format!("{}e{}{:02}", m, if e < 0 { "-" } else { "+" }, e.abs())
    } else {
        let decimals = (5 - exp).max(0) as usize;
        let s = format!("{:.*}", decimals, v);
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s }
    }
}
