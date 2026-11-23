//! `poisstat`, `canoncorr`, `nbinfit`, `gevfit` and `gpfit`.

use super::regress::{qr_pivot, solve_upper};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;


const HELPERS: &str = "function v = matlivre_gev_nll(p, x)
  k = p(1); lnsigma = p(2); sigma = exp(lnsigma); mu = p(3);
  n = numel(x); z = (x - mu) ./ sigma;
  if abs(k) > eps
    u = 1 + k.*z;
    if min(u) > 0
      lnu = log1p(k.*z); t = exp(-(1/k)*lnu);
      v = n*lnsigma + sum(t) + (1+1/k)*sum(lnu);
    else
      v = Inf;
    end
  else
    v = n*lnsigma + sum(exp(-z) + z);
  end
end
function v = matlivre_gp_nll(p, x)
  k = p(1); lnsigma = p(2); sigma = exp(lnsigma);
  n = numel(x); z = x ./ sigma;
  if abs(k) > eps
    if k > 0 || max(z) < -1/k
      v = n*lnsigma + (1+1/k) * sum(log1p(k.*z));
    else
      v = Inf;
    end
  else
    v = n*lnsigma + sum(z);
  end
end
function v = matlivre_nb_nll(r, x, n, sumx, tolbnd)
  if r < tolbnd
    v = Inf;
  else
    xbar = sumx / n;
    v = -sum(gammaln(r+x)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r));
  end
end
";

fn ensure_helpers(it: &mut Interpreter) -> Res<()> {
    if it.has_function("matlivre_gev_nll") {
        return Ok(());
    }
    let prog = crate::parser::parse(HELPERS).map_err(|e| crate::error::LangError::new(e.message.clone(), 0))?;
    it.define_functions(&prog);
    Ok(())
}

pub const NAMES: &[&str] = &["poisstat", "canoncorr", "nbinfit", "gevfit", "gpfit"];

fn num<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(_) => rt("Input must be numeric."),
        None => rt("Not enough input arguments."),
    }
}

fn lit(v: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&format!("{:e}", x));
    }
    s.push(']');
    s
}

fn anon(it: &mut Interpreter, src: &str) -> Res<Value> {
    Ok(it.call_function("str2func", vec![Value::chars(src)], 1)?.into_iter().next().unwrap())
}

fn optset(it: &mut Interpreter, tolx: f64, tolfun: f64, maxfun: f64, maxiter: f64) -> Res<Value> {
    let args = vec![Value::chars("Display"), Value::chars("off"), Value::chars("TolX"), Value::num(tolx), Value::chars("TolFun"), Value::num(tolfun), Value::chars("MaxFunEvals"), Value::num(maxfun), Value::chars("MaxIter"), Value::num(maxiter)];
    Ok(it.call_function("optimset", args, 1)?.into_iter().next().unwrap())
}

fn options_of(a: &[Value], idx: usize) -> (f64, f64) {
    let (mut tolx, mut tolfun) = (1e-6, 1e-6);
    if let Some(Value::Struct(s)) = a.get(idx) {
        if let Some(Value::Num(t)) = s.get("TolX") {
            if !t.is_empty() {
                tolx = t.data[0];
            }
        }
        if let Some(Value::Num(t)) = s.get("TolFun") {
            if !t.is_empty() {
                tolfun = t.data[0];
            }
        }
    }
    (tolx, tolfun)
}

fn poisstat(a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Requires one input argument.");
    }
    let m = num(a, 0)?;
    if m.kind == Kind::Char {
        return rt("Argument following 'like' must be a variable with a data type that supports ZEROS, such as a double or single variable.");
    }
    let single = is_single(m);
    let v: Vec<f64> = m.data.iter().map(|&l| if l <= 0.0 { f64::NAN } else if l > 0.0 { l } else { 0.0 }).collect();
    let o = Value::Num(finish(Matrix::from_dims(m.dims.clone(), v), single));
    Ok(Some(vec![o.clone(), o]))
}

fn pf_upper(it: &mut Interpreter, name: &str, x: f64, p: &[f64]) -> Res<f64> {
    let mut args = vec![Value::num(x)];
    args.extend(p.iter().map(|&v| Value::num(v)));
    args.push(Value::chars("upper"));
    it.call_function(name, args, 1)?[0].scalar()
}

fn canoncorr(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Requires two arguments.");
    }
    let (x, y) = (num(a, 0)?, num(a, 1)?);
    let single = is_single(x) || is_single(y);
    let (n, p1) = (x.dims[0], x.dims[1]);
    if y.dims[0] != n {
        return rt("X and Y must have the same number of rows.");
    }
    if n == 1 {
        return rt("X and Y must have at least two rows.");
    }
    let p2 = y.dims[1];
    let center = |m: &Matrix| -> Vec<Vec<f64>> {
        let (r, c) = (m.dims[0], m.dims[1]);
        let mu: Vec<f64> = (0..c).map(|j| (0..r).map(|i| m.data[j * r + i]).sum::<f64>() / r as f64).collect();
        (0..r).map(|i| (0..c).map(|j| m.data[j * r + i] - mu[j]).collect()).collect()
    };
    let (xc, yc) = (center(x), center(y));
    type Red = (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<usize>, usize);
    let reduce = |rows: &Vec<Vec<f64>>, p: usize, name: &str| -> Res<Red> {
        let (q, r, perm) = qr_pivot(rows, n, p);
        let tol = crate::stats::util::eps_of(r[0][0].abs()) * n.max(p) as f64;
        let rank = (0..p).filter(|&j| r[j][j].abs() > tol).count();
        if rank == 0 {
            return rt(format!("{name} must contain at least one non-constant column."));
        }
        let qk: Vec<Vec<f64>> = (0..n).map(|i| q[i][..rank].to_vec()).collect();
        let tk: Vec<Vec<f64>> = (0..rank).map(|i| r[i][..rank].to_vec()).collect();
        Ok((qk, tk, perm, rank))
    };
    let (q1, t11, perm1, rank_x) = reduce(&xc, p1, "X")?;
    let (q2, t22, perm2, rank_y) = reduce(&yc, p2, "Y")?;
    let d = rank_x.min(rank_y);
    let q1tq2: Vec<f64> = (0..rank_x * rank_y)
        .map(|idx| {
            let (i, j) = (idx % rank_x, idx / rank_x);
            (0..n).map(|r| q1[r][i] * q2[r][j]).sum()
        })
        .collect();
    let sv = it.call_function("svd", vec![Value::Num(Matrix::from_dims(vec![rank_x, rank_y], q1tq2)), Value::num(0.0)], 3)?;
    let (Value::Num(l), Value::Num(dm), Value::Num(m)) = (&sv[0], &sv[1], &sv[2]) else { return rt("svd failed") };
    let scale = ((n - 1) as f64).sqrt();
    let solve = |t: &Vec<Vec<f64>>, k: usize, cols: &Matrix| -> Vec<Vec<f64>> {
        (0..d)
            .map(|c| {
                let rhs: Vec<f64> = (0..k).map(|i| cols.data[c * cols.dims[0] + i] * scale).collect();
                solve_upper(t, &rhs, k)
            })
            .collect()
    };
    let a_cols = solve(&t11, rank_x, l);
    let b_cols = solve(&t22, rank_y, m);
    let unperm = |cols: &Vec<Vec<f64>>, perm: &Vec<usize>, rank: usize, p: usize| -> Vec<f64> {
        let mut out = vec![0.0; p * d];
        for c in 0..d {
            for (j, &pc) in perm.iter().take(rank).enumerate() {
                out[c * p + pc] = cols[c][j];
            }
        }
        out
    };
    let a_mat = unperm(&a_cols, &perm1, rank_x, p1);
    let b_mat = unperm(&b_cols, &perm2, rank_y, p2);
    let r: Vec<f64> = (0..d).map(|i| dm.data[i * dm.dims[0] + i].max(0.0).min(1.0)).collect();
    let mut outs = vec![
        Value::Num(finish(Matrix::from_dims(vec![p1, d], a_mat.clone()), single)),
        Value::Num(finish(Matrix::from_dims(vec![p2, d], b_mat.clone()), single)),
        Value::Num(finish(Matrix::from_dims(vec![1, d], r.clone()), single)),
    ];
    if nargout > 3 {
        let uv = |rows: &Vec<Vec<f64>>, coef: &Vec<f64>, p: usize| -> Vec<f64> {
            let mut o = vec![0.0; n * d];
            for c in 0..d {
                for i in 0..n {
                    o[c * n + i] = (0..p).map(|j| rows[i][j] * coef[c * p + j]).sum();
                }
            }
            o
        };
        outs.push(Value::Num(finish(Matrix::from_dims(vec![n, d], uv(&xc, &a_mat, p1)), single)));
        outs.push(Value::Num(finish(Matrix::from_dims(vec![n, d], uv(&yc, &b_mat, p2)), single)));
    }
    if nargout > 5 {
        let k: Vec<f64> = (0..d).map(|v| v as f64).collect();
        let d1k: Vec<f64> = k.iter().map(|v| rank_x as f64 - v).collect();
        let d2k: Vec<f64> = k.iter().map(|v| rank_y as f64 - v).collect();
        let mut log_lambda = vec![f64::NEG_INFINITY; d];
        let nondegen: Vec<usize> = (0..d).filter(|&i| r[i] < 1.0).collect();
        let mut acc = 0.0;
        for &i in nondegen.iter().rev() {
            acc += (1.0 - r[i] * r[i]).ln();
            log_lambda[i] = acc;
        }
        let wilks: Vec<f64> = log_lambda.iter().map(|v| v.exp()).collect();
        let mut s = vec![1.0; d];
        for i in 0..d {
            if d1k[i] * d2k[i] > 2.0 {
                s[i] = ((d1k[i].powi(2) * d2k[i].powi(2) - 4.0) / (d1k[i].powi(2) + d2k[i].powi(2) - 5.0)).sqrt();
            }
        }
        let df1: Vec<f64> = (0..d).map(|i| d1k[i] * d2k[i]).collect();
        let nn = n as f64;
        let df2: Vec<f64> = (0..d).map(|i| (nn - 0.5 * (rank_x + rank_y + 3) as f64) * s[i] - 0.5 * d1k[i] * d2k[i] + 1.0).collect();
        let mut ratio = vec![f64::INFINITY; d];
        for &i in &nondegen {
            let pl = wilks[i].powf(1.0 / s[i]);
            ratio[i] = (1.0 - pl) / pl;
        }
        let fstat: Vec<f64> = (0..d).map(|i| ratio[i] * df2[i] / df1[i]).collect();
        let mut pf = Vec::new();
        for i in 0..d {
            pf.push(pf_upper(it, "fcdf", fstat[i], &[df1[i], df2[i]])?);
        }
        let mut cum = vec![0.0; d];
        let mut c = 0.0;
        for i in 0..d {
            if i > 0 {
                c += (1.0 / r[i - 1]).powi(2);
            }
            cum[i] = c;
        }
        let chisq: Vec<f64> = (0..d).map(|i| -(nn - k[i] - 0.5 * (rank_x + rank_y + 3) as f64 + cum[i]) * log_lambda[i]).collect();
        let mut pc = Vec::new();
        for i in 0..d {
            pc.push(pf_upper(it, "chi2cdf", chisq[i], &[df1[i]])?);
        }
        let row = |v: Vec<f64>| Value::Num(finish(Matrix::from_dims(vec![1, v.len()], v), single));
        let mut st = StructArray::scalar();
        st.set("Wilks", row(wilks));
        st.set("df1", row(df1.clone()));
        st.set("df2", row(df2));
        st.set("F", row(fstat));
        st.set("pF", row(pf));
        st.set("chisq", row(chisq));
        st.set("pChisq", row(pc.clone()));
        st.set("dfe", row(df1));
        st.set("p", row(pc));
        outs.push(Value::Struct(st));
    }
    Ok(Some(outs))
}

fn nbinfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = num(a, 0)?;
    if x.dims.len() > 2 || (x.dims[0] > 1 && x.dims[1] > 1) {
        return rt("The first argument must be a vector.");
    }
    if x.dims == [0, 0] {
        return rt("Non-scalar arguments must match in size.");
    }
    let alpha = match a.get(1) {
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        _ => 0.05,
    };
    if x.data.iter().any(|&v| v < 0.0 || !v.is_finite() || v != v.round()) {
        return rt("Negative or non-integer data.");
    }
    if is_single(x) {
        return rt("FMINSEARCH accepts inputs only of data type double.");
    }
    let n = x.numel() as f64;
    let xbar = x.data.iter().sum::<f64>() / n;
    let s2 = x.data.iter().map(|v| (v - xbar).powi(2)).sum::<f64>() / (n - 1.0);
    if s2 <= xbar {
        return Ok(Some(vec![Value::Num(Matrix::from_dims(vec![1, 2], vec![f64::INFINITY, 1.0])), Value::Num(Matrix::from_dims(vec![2, 2], vec![f64::INFINITY, f64::INFINITY, 1.0, 1.0]))]));
    }
    let r0 = xbar * xbar / (s2 - xbar);
    let sumx: f64 = x.data.iter().sum();
    let (tolx, tolfun) = options_of(a, 2);
    ensure_helpers(it)?;
    let body = format!("@(r) matlivre_nb_nll(r, {xs}, {n:e}, {sumx:e}, 1e-6)", xs = lit(&x.data));
    let f2 = anon(it, &body)?;
    let opts = optset(it, tolx, tolfun, 400.0, 200.0)?;
    let res = it.call_function("fminsearch", vec![f2, Value::num(r0), opts], 1)?;
    let rhat = res[0].scalar()?;
    let phat = rhat / (xbar + rhat);
    let ph = Value::Num(Matrix::from_dims(vec![1, 2], vec![rhat, phat]));
    if nargout < 2 {
        return Ok(Some(vec![ph]));
    }
    let tri_sum: f64 = {
        let mut s = 0.0;
        for &v in &x.data {
            s += it.call_function("psi", vec![Value::num(1.0), Value::num(rhat + v)], 1)?[0].scalar()?;
        }
        s - n * it.call_function("psi", vec![Value::num(1.0), Value::num(rhat)], 1)?[0].scalar()?
    };
    let (l11, l12, l22) = (tri_sum, n / phat, -n * rhat / (phat * phat) - sumx / (1.0 - phat).powi(2));
    let (h11, h12, h22) = (-l11, -l12, -l22);
    let det = h11 * h22 - h12 * h12;
    let (v11, v22) = (h22 / det, h11 / det);
    let z = super::dist::norminv_std(1.0 - alpha / 2.0);
    let (s1, s2v) = (v11.sqrt(), v22.sqrt());
    let ci = vec![rhat - z * s1, rhat + z * s1, phat - z * s2v, phat + z * s2v];
    Ok(Some(vec![ph, Value::Num(Matrix::from_dims(vec![2, 2], ci))]))
}

fn prepare(a: &[Value], who: &str) -> Res<(Vec<f64>, f64, bool)> {
    let x = num(a, 0)?;
    if x.dims.len() > 2 || (x.dims[0] > 1 && x.dims[1] > 1) || x.is_empty() && x.dims[0] != 0 || (x.dims[0] == 0 && x.dims[1] == 0) {
        return rt("X must be a vector.");
    }
    if who == "gpfit" && x.data.iter().any(|&v| v <= 0.0) {
        return rt("The data in X must be positive.");
    }
    let alpha = match a.get(1) {
        Some(Value::Num(m)) if !m.is_empty() => m.data[0],
        _ => 0.05,
    };
    let mut v = x.data.clone();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Ok((v, alpha, is_single(x)))
}

fn gevfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let (x, alpha, single) = prepare(a, "gevfit")?;
    let n = x.len();
    if n == 0 {
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 3], vec![f64::NAN; 3]), single)), Value::Num(finish(Matrix::from_dims(vec![2, 3], vec![f64::NAN; 6]), single))]));
    }
    let (xmin, xmax) = (x[0], x[n - 1]);
    let rangex = xmax - xmin;
    if !rangex.is_finite() {
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 3], vec![f64::NAN; 3]), single)), Value::Num(finish(Matrix::from_dims(vec![2, 3], vec![f64::NAN; 6]), single))]));
    }
    if rangex < f64::MIN_POSITIVE {
        let ph = vec![0.0, 0.0, x[0]];
        let ci = if n == 1 { vec![f64::NEG_INFINITY, f64::INFINITY, 0.0, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY] } else { vec![0.0, 0.0, 0.0, 0.0, x[0], x[0]] };
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 3], ph), single)), Value::Num(finish(Matrix::from_dims(vec![2, 3], ci), single))]));
    }
    let (tolx, tolfun) = options_of(a, 2);
    let nf = n as f64;
    let fv: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / nf).collect();
    let xs = lit(&x);
    let fl = lit(&fv);
    let fun0 = anon(it, &format!("@(k) 1-corr({xs}', gevinv({fl}',k,1,0))"))?;
    let opts0 = it.call_function("optimset", vec![], 1)?.into_iter().next().unwrap();
    let k0 = it.call_function("fminsearch", vec![fun0, Value::num(0.0), opts0], 1)?[0].scalar()?;
    let gi = it.call_function("gevinv", vec![Value::Num(Matrix::from_dims(vec![n, 1], fv.clone())), Value::num(k0), Value::num(1.0), Value::num(0.0)], 1)?;
    let pf = it.call_function("polyfit", vec![gi[0].clone(), Value::Num(Matrix::from_dims(vec![n, 1], x.clone())), Value::num(1.0)], 1)?;
    let Value::Num(pm) = &pf[0] else { return rt("polyfit failed") };
    let (mut sigma0, mut mu0) = (pm.data[0], pm.data[1]);
    let mut k0 = k0;
    if (k0 < 0.0 && xmax > -sigma0 / k0 + mu0) || (k0 > 0.0 && xmin < -sigma0 / k0 + mu0) {
        k0 = 0.0;
        let ev = it.call_function("evfit", vec![Value::Num(Matrix::from_dims(vec![n, 1], x.clone()))], 1)?;
        let Value::Num(em) = &ev[0] else { return rt("evfit failed") };
        sigma0 = em.data[1];
        mu0 = em.data[0];
    }
    let start = Value::Num(Matrix::from_dims(vec![1, 3], vec![k0, sigma0.ln(), mu0]));
    ensure_helpers(it)?;
    let body = format!("@(p) matlivre_gev_nll(p, {xs}')", xs = xs);
    let fun = anon(it, &body)?;
    let opts = optset(it, tolx, tolfun, 400.0, 200.0)?;
    let res = it.call_function("fminsearch", vec![fun, start, opts], 1)?;
    let Value::Num(rp) = &res[0] else { return rt("fminsearch failed") };
    let ph = [rp.data[0], rp.data[1].exp(), rp.data[2]];
    let tolbnd = 1e-6;
    let at_boundary = (ph[0] < 0.0 && xmax > -ph[1] / ph[0] + ph[2] - tolbnd) || (ph[0] > 0.0 && xmin < -ph[1] / ph[0] + ph[2] + tolbnd);
    let mut ci = vec![f64::NAN; 6];
    if nargout > 1 && !at_boundary {
        let lk = it.call_function("gevlike", vec![Value::Num(Matrix::from_dims(vec![1, 3], ph.to_vec())), Value::Num(Matrix::from_dims(vec![n, 1], x.clone()))], 2)?;
        let Value::Num(acov) = &lk[1] else { return rt("gevlike failed") };
        let se: Vec<f64> = (0..3).map(|i| acov.data[i * 3 + i].sqrt()).collect();
        let z = super::dist::norminv_std(1.0 - alpha / 2.0);
        let k_ci = [ph[0] - z * se[0], ph[0] + z * se[0]];
        let s_ci = [(ph[1].ln() - z * se[1] / ph[1]).exp(), (ph[1].ln() + z * se[1] / ph[1]).exp()];
        let m_ci = [ph[2] - z * se[2], ph[2] + z * se[2]];
        ci = vec![k_ci[0], k_ci[1], s_ci[0], s_ci[1], m_ci[0], m_ci[1]];
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 3], ph.to_vec()), single)), Value::Num(finish(Matrix::from_dims(vec![2, 3], ci), single))]))
}

fn gpfit(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let (x, alpha, single) = prepare(a, "gpfit")?;
    let n = x.len();
    if n == 0 {
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 2], vec![f64::NAN; 2]), single)), Value::Num(finish(Matrix::from_dims(vec![2, 2], vec![f64::NAN; 4]), single))]));
    }
    let xmax = x[n - 1];
    let rangex = xmax - x[0];
    if !rangex.is_finite() {
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 2], vec![f64::NAN; 2]), single)), Value::Num(finish(Matrix::from_dims(vec![2, 2], vec![f64::NAN; 4]), single))]));
    }
    if rangex < f64::MIN_POSITIVE {
        let ph = if xmax <= f64::MAX.sqrt() { vec![f64::NAN, 0.0] } else { vec![f64::NEG_INFINITY, f64::INFINITY] };
        let ci = vec![ph[0], ph[0], ph[1], ph[1]];
        return Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 2], ph), single)), Value::Num(finish(Matrix::from_dims(vec![2, 2], ci), single))]));
    }
    let (tolx, tolfun) = options_of(a, 2);
    let nf = n as f64;
    let xbar = x.iter().sum::<f64>() / nf;
    let s2 = x.iter().map(|v| (v - xbar).powi(2)).sum::<f64>() / (nf - 1.0);
    let mut k0 = -0.5 * (xbar * xbar / s2 - 1.0);
    let mut sigma0 = 0.5 * xbar * (xbar * xbar / s2 + 1.0);
    if k0 < 0.0 && xmax >= -sigma0 / k0 {
        k0 = 0.0;
        sigma0 = xbar;
    }
    let xs = lit(&x);
    ensure_helpers(it)?;
    let body = format!("@(p) matlivre_gp_nll(p, {xs}')", xs = xs);
    let fun = anon(it, &body)?;
    let opts = optset(it, tolx, tolfun, 400.0, 200.0)?;
    let start = Value::Num(Matrix::from_dims(vec![1, 2], vec![k0, sigma0.ln()]));
    let res = it.call_function("fminsearch", vec![fun, start, opts], 1)?;
    let Value::Num(rp) = &res[0] else { return rt("fminsearch failed") };
    let ph = [rp.data[0], rp.data[1].exp()];
    let at_boundary = (ph[0] < 0.0 && xmax > -ph[1] / ph[0] - 1e-6) || ph[0] <= -0.5;
    let mut ci = vec![f64::NAN; 4];
    if nargout > 1 && !at_boundary {
        let lk = it.call_function("gplike", vec![Value::Num(Matrix::from_dims(vec![1, 2], ph.to_vec())), Value::Num(Matrix::from_dims(vec![n, 1], x.clone()))], 2)?;
        let Value::Num(acov) = &lk[1] else { return rt("gplike failed") };
        let se: Vec<f64> = (0..2).map(|i| acov.data[i * 2 + i].sqrt()).collect();
        let z = super::dist::norminv_std(1.0 - alpha / 2.0);
        ci = vec![ph[0] - z * se[0], ph[0] + z * se[0], (ph[1].ln() - z * se[1] / ph[1]).exp(), (ph[1].ln() + z * se[1] / ph[1]).exp()];
    }
    Ok(Some(vec![Value::Num(finish(Matrix::from_dims(vec![1, 2], ph.to_vec()), single)), Value::Num(finish(Matrix::from_dims(vec![2, 2], ci), single))]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "poisstat" => poisstat(a),
        "canoncorr" => canoncorr(it, a, nargout),
        "nbinfit" => nbinfit(it, a, nargout),
        "gevfit" => gevfit(it, a, nargout),
        _ => gpfit(it, a, nargout),
    }
}
