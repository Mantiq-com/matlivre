//! `plsregress`: partial least squares regression (SIMPLS).

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["plsregress"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "plsregress" => plsregress(it, a, nargout),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn cols_of(m: &Matrix) -> Vec<Vec<f64>> {
    let (r, c) = (m.dims[0], m.dims[1..].iter().product::<usize>());
    (0..c).map(|j| m.data[j * r..(j + 1) * r].to_vec()).collect()
}

fn mat(cols: &[Vec<f64>], rows: usize, single: bool) -> Value {
    Value::Num(finish(Matrix::from_dims(vec![rows, cols.len()], cols.iter().flat_map(|c| c.iter().copied()).collect()), single))
}

fn on_off(v: &Value, what: &str) -> Result<bool, crate::error::LangError> {
    match v {
        Value::Num(m) if m.data.len() == 1 && m.kind != Kind::Char => Ok(m.data[0] != 0.0),
        v if v.is_text() => match v.text().unwrap_or_default().to_lowercase().as_str() {
            "on" | "true" => Ok(true),
            "off" | "false" => Ok(false),
            _ => rt(format!("'{}' is not a valid value for the {what} argument.  Valid values are: 'on', 'off'.", v.text().unwrap_or_default())),
        },
        _ => rt(format!("Invalid {what} data type. Argument must be a character vector or string scalar, or a logical or numeric scalar.")),
    }
}

fn plsregress(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Requires at least two inputs.");
    }
    let (Value::Num(xm), Value::Num(ym)) = (&a[0], &a[1]) else { return rt("X and Y must be numeric.") };
    let (n, dx) = (xm.dims[0], xm.dims[1]);
    let (ny, dy) = (ym.dims[0], ym.dims[1]);
    if ny != n {
        return rt("X and Y must have the same number of rows.");
    }
    let maxn = (n as f64 - 1.0).min(dx as f64) as usize;
    let ncomp: usize = if a.len() < 3 {
        maxn
    } else {
        match &a[2] {
            Value::Num(m) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] == m.data[0].round() && m.data[0] > 0.0 => {
                if m.data[0] as usize > maxn {
                    return rt(format!("NCOMP must be less than or equal to {maxn} for these data."));
                }
                m.data[0] as usize
            }
            _ => return rt("NCOMP must be a positive integer value."),
        }
    };
    // options
    let rest = if a.len() > 3 { &a[3..] } else { &[][..] };
    if rest.len() % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    let names = ["cv", "mcreps", "options", "intercept"];
    let mut cv: Option<&Value> = None;
    let mut mcreps: Option<&Value> = None;
    let mut intercept_v: Option<&Value> = None;
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match pick {
            Some(0) => cv = Some(&pair[1]),
            Some(1) => mcreps = Some(&pair[1]),
            Some(2) => {}
            Some(3) => intercept_v = Some(&pair[1]),
            None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            _ => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    let resub = match cv {
        None => true,
        Some(v) if v.is_text() && v.text().unwrap_or_default() == "resubstitution" => true,
        Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] == m.data[0].round() && m.data[0] > 0.0 && m.data[0] <= n as f64 => false,
        Some(_) => return rt("CV must be 'resubstitution', a positive integer no larger than the number of rows in X, or a partition created with CVPARTITION."),
    };
    let mc = match mcreps {
        None => 1.0,
        Some(Value::Num(m)) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] == m.data[0].round() && m.data[0] > 0.0 => m.data[0],
        Some(_) => return rt("MCReps must be a positive integer."),
    };
    if mc > 1.0 && resub {
        return rt("MCReps must be 1 for resubstitution.");
    }
    let with_intercept = match intercept_v {
        None => true,
        Some(v) => on_off(v, "Intercept")?,
    };
    if !resub && nargout > 6 {
        return rt("Cross-validated mean squared errors are not implemented.");
    }
    let single = is_single(xm) || is_single(ym);
    let mut x0 = cols_of(xm);
    let mut y0 = cols_of(ym);
    let mut mean_x = vec![0.0; dx];
    let mut mean_y = vec![0.0; dy];
    if with_intercept {
        for (j, c) in x0.iter_mut().enumerate() {
            let m = c.iter().sum::<f64>() / n as f64;
            mean_x[j] = m;
            c.iter_mut().for_each(|v| *v -= m);
        }
        for (j, c) in y0.iter_mut().enumerate() {
            let m = c.iter().sum::<f64>() / n as f64;
            mean_y[j] = m;
            c.iter_mut().for_each(|v| *v -= m);
        }
    }
    // SIMPLS
    let mut xload = vec![vec![0.0; dx]; ncomp];
    let mut yload = vec![vec![0.0; dy]; ncomp];
    let mut xscore = vec![vec![0.0; n]; ncomp];
    let mut yscore = vec![vec![0.0; n]; ncomp];
    let mut weights = vec![vec![0.0; dx]; ncomp];
    let mut vmat: Vec<Vec<f64>> = Vec::new();
    // Cov = X0' * Y0  (dx x dy), kept as columns
    let mut cov: Vec<Vec<f64>> = (0..dy).map(|j| (0..dx).map(|i| dot(&x0[i], &y0[j])).collect()).collect();
    for i in 0..ncomp {
        let cm = Value::Num(Matrix::from_dims(vec![dx, dy], cov.iter().flat_map(|c| c.iter().copied()).collect()));
        let r = it.call_function("svd", vec![cm, Value::chars("econ")], 3)?;
        let (Value::Num(um), Value::Num(sm), Value::Num(vm)) = (&r[0], &r[1], &r[2]) else { return rt("svd returned no value") };
        let ri: Vec<f64> = um.data[..dx].to_vec();
        let ci: Vec<f64> = vm.data[..dy].to_vec();
        let si = sm.data[0];
        let mut ti: Vec<f64> = (0..n).map(|k| (0..dx).map(|j| x0[j][k] * ri[j]).sum()).collect();
        let normti = dot(&ti, &ti).sqrt();
        ti.iter_mut().for_each(|v| *v /= normti);
        for j in 0..dx {
            xload[i][j] = dot(&x0[j], &ti);
        }
        let qi: Vec<f64> = ci.iter().map(|c| si * c / normti).collect();
        yload[i] = qi.clone();
        xscore[i] = ti;
        yscore[i] = (0..n).map(|k| (0..dy).map(|j| y0[j][k] * qi[j]).sum()).collect();
        weights[i] = ri.iter().map(|v| v / normti).collect();
        let mut vi = xload[i].clone();
        for _ in 0..2 {
            for vj in &vmat {
                let c = dot(vj, &vi);
                for k in 0..dx {
                    vi[k] -= c * vj[k];
                }
            }
        }
        let nv = dot(&vi, &vi).sqrt();
        vi.iter_mut().for_each(|v| *v /= nv);
        vmat.push(vi.clone());
        // Cov = Cov - vi*(vi'*Cov)
        for c in cov.iter_mut() {
            let d = dot(&vi, c);
            for k in 0..dx {
                c[k] -= vi[k] * d;
            }
        }
        // Cov = Cov - Vi*(Vi'*Cov)
        for c in cov.iter_mut() {
            let ds: Vec<f64> = vmat.iter().map(|v| dot(v, c)).collect();
            for (v, d) in vmat.iter().zip(&ds) {
                for k in 0..dx {
                    c[k] -= v[k] * d;
                }
            }
        }
    }
    for i in 0..ncomp {
        let mut ui = yscore[i].clone();
        for _ in 0..2 {
            for j in 0..i {
                let c = dot(&xscore[j], &ui);
                for k in 0..n {
                    ui[k] -= c * xscore[j][k];
                }
            }
        }
        yscore[i] = ui;
    }
    let mut out = vec![mat(&xload, dx, single), mat(&yload, dy, single)];
    if nargout > 2 {
        out.push(mat(&xscore, n, single));
        out.push(mat(&yscore, n, single));
    }
    if nargout > 4 {
        // beta = Weights * Yloadings'  (dx x dy)
        let beta_cols: Vec<Vec<f64>> = (0..dy).map(|j| (0..dx).map(|i| (0..ncomp).map(|k| weights[k][i] * yload[k][j]).sum()).collect()).collect();
        let (beta, rows) = if with_intercept {
            let cols: Vec<Vec<f64>> = (0..dy)
                .map(|j| {
                    let mut c = vec![mean_y[j] - (0..dx).map(|i| mean_x[i] * beta_cols[j][i]).sum::<f64>()];
                    c.extend(beta_cols[j].iter().copied());
                    c
                })
                .collect();
            (cols, dx + 1)
        } else {
            (beta_cols, dx)
        };
        out.push(mat(&beta, rows, single));
        let sumsq = |cols: &[Vec<f64>]| -> f64 { cols.iter().map(|c| dot(c, c)).sum() };
        if nargout > 5 {
            let tx = sumsq(&x0);
            let ty = sumsq(&y0);
            let mut pv = vec![0.0; 2 * ncomp];
            for k in 0..ncomp {
                pv[2 * k] = dot(&xload[k], &xload[k]) / tx;
                pv[2 * k + 1] = dot(&yload[k], &yload[k]) / ty;
            }
            out.push(Value::Num(finish(Matrix::from_dims(vec![2, ncomp], pv), single)));
        }
        // reconstructions
        let recon = |scores: &[Vec<f64>], load: &[Vec<f64>], upto: usize, dcols: usize| -> Vec<Vec<f64>> {
            (0..dcols).map(|j| (0..n).map(|k| (0..upto).map(|i| scores[i][k] * load[i][j]).sum()).collect()).collect()
        };
        if nargout > 6 {
            let mut mse = vec![0.0; 2 * (ncomp + 1)];
            mse[0] = sumsq(&x0);
            mse[1] = sumsq(&y0);
            for i in 1..=ncomp {
                let xr = recon(&xscore, &xload, i, dx);
                let yr = recon(&xscore, &yload, i, dy);
                let ex: f64 = (0..dx).map(|j| (0..n).map(|k| (x0[j][k] - xr[j][k]).powi(2)).sum::<f64>()).sum();
                let ey: f64 = (0..dy).map(|j| (0..n).map(|k| (y0[j][k] - yr[j][k]).powi(2)).sum::<f64>()).sum();
                mse[2 * i] = ex;
                mse[2 * i + 1] = ey;
            }
            mse.iter_mut().for_each(|v| *v /= n as f64);
            out.push(Value::Num(finish(Matrix::from_dims(vec![2, ncomp + 1], mse), single)));
        }
        if nargout > 7 {
            let mut st = StructArray::scalar();
            st.set("W", mat(&weights, dx, single));
            let vars: Vec<f64> = (0..ncomp)
                .map(|i| {
                    let m = xscore[i].iter().sum::<f64>() / n as f64;
                    xscore[i].iter().map(|v| (v - m).powi(2)).sum::<f64>() / (n as f64 - 1.0)
                })
                .collect();
            let t2: Vec<f64> = (0..n).map(|k| (0..ncomp).map(|i| xscore[i][k] * xscore[i][k] / vars[i]).sum()).collect();
            st.set("T2", Value::Num(finish(Matrix::col(t2), single)));
            let xr = recon(&xscore, &xload, ncomp, dx);
            let yr = recon(&xscore, &yload, ncomp, dy);
            let xres: Vec<Vec<f64>> = (0..dx).map(|j| (0..n).map(|k| x0[j][k] - xr[j][k]).collect()).collect();
            let yres: Vec<Vec<f64>> = (0..dy).map(|j| (0..n).map(|k| y0[j][k] - yr[j][k]).collect()).collect();
            st.set("Xresiduals", mat(&xres, n, single));
            st.set("Yresiduals", mat(&yres, n, single));
            out.push(Value::Struct(st));
        }
    }
    Ok(Some(out))
}
