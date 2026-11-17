//! `tachorpm`: rotational speed from a tachometer pulse train.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["tachorpm"];

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    if nargout > 3 {
        return rt("Too many output arguments.");
    }
    let x = match &a[0] {
        Value::Num(m) if m.kind != Kind::Char && m.kind != Kind::Logical => m.clone(),
        Value::Num(_) => return rt("Expected X to be one of these types:  double, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64"),
        v => return rt(format!("Expected X to be one of these types:  double, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64  Instead its type was {}.", v.class_name())),
    };
    let fsm = match &a[1] {
        Value::Num(m) if m.kind != Kind::Char => m.clone(),
        _ => return rt("Expected Fs to be positive."),
    };
    if fsm.numel() != 1 {
        return rt("Expected Fs to be a scalar.");
    }
    let fs = fsm.data[0];
    if !(fs > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    if (a.len() - 2) % 2 != 0 {
        return rt("Optional name-value arguments must be in pairs.");
    }
    let mut ppr = 1.0;
    let mut ofs = fs;
    let mut levels: Option<Matrix> = None;
    let mut fit = "smooth".to_string();
    let mut fpts = 10.0;
    let mut i = 2;
    while i < a.len() {
        let s = a[i].text().unwrap_or_default();
        let v = &a[i + 1];
        let m = v.to_matrix().unwrap_or_else(|_| Matrix::from_dims(vec![0, 0], vec![]));
        match s.to_lowercase().as_str() {
            "pulsesperrev" => {
                if m.numel() != 1 {
                    return rt("Expected PPR to be a scalar.");
                }
                if !(m.data[0] > 0.0) {
                    return rt("Expected PPR to be positive.");
                }
                ppr = m.data[0];
            }
            "outputfs" => {
                if m.numel() != 1 {
                    return rt("Expected OFs to be a scalar.");
                }
                if !(m.data[0] > 0.0) {
                    return rt("Expected OFs to be positive.");
                }
                ofs = m.data[0];
            }
            "statelevels" => {
                if m.numel() != 2 {
                    return rt("Expected SL to be an array with number of elements equal to 2.");
                }
                if !(m.data[1] > m.data[0]) {
                    return rt("The state levels must be in increasing order.");
                }
                levels = Some(m);
            }
            "fittype" => {
                let t = v.text().unwrap_or_default();
                match t.to_lowercase().as_str() {
                    "linear" | "smooth" => fit = t.to_lowercase(),
                    _ => return rt(format!("Expected FT to match one of these values:  'linear', 'smooth'  The input, '{t}', did not match any of the valid values.")),
                }
            }
            "fitpoints" => {
                if m.numel() != 1 {
                    return rt("Expected FP to be a scalar.");
                }
                if !(m.data[0] > 0.0) {
                    return rt("Expected FP to be positive.");
                }
                fpts = m.data[0];
            }
            _ => return rt(format!("'{s}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    if !(x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1)) {
        return rt("Expected X to be a vector.");
    }
    if x.im.is_some() {
        return rt("Expected X to be real.");
    }
    let n = x.numel();
    let xv = real_col(x.data.clone());
    let sl = match &levels {
        Some(m) => m.data.clone(),
        None => it.call_function("statelevels", vec![xv.clone()], 1)?.remove(0).to_matrix()?.data,
    };
    if sl.len() < 2 || !(sl[1] > sl[0]) {
        return rt("The state levels must be in increasing order.");
    }
    let lv = real_row(vec![sl[0], sl[1]]);
    let edge = |it: &mut Interpreter, name: &str| -> R<Vec<f64>> {
        let r = it.call_function(name, vec![xv.clone(), scalar_value(fs), Value::string("StateLevels"), lv.clone()], 3);
        match r {
            Ok(v) => {
                let lt = v[1].to_matrix()?.data;
                let ut = v[2].to_matrix()?.data;
                Ok(lt.iter().zip(ut.iter()).map(|(a, b)| (a + b) / 2.0).collect())
            }
            Err(_) => Ok(vec![]),
        }
    };
    let rises = edge(it, "risetime")?;
    let falls = edge(it, "falltime")?;
    let mut edges: Vec<f64> = rises.iter().chain(falls.iter()).cloned().collect();
    edges.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
    let tp: Vec<f64> = edges.chunks(2).filter(|c| c.len() == 2).map(|c| (c[0] + c[1]) / 2.0).collect();
    if tp.is_empty() {
        return rt("No tachometer pulses were detected. Try adjusting the state levels.");
    }
    if tp.len() < 3 {
        return rt("Tachometer pulse vector X must contain at least 3 pulses.");
    }
    let tm: Vec<f64> = tp.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();
    let rpm_pts: Vec<f64> = tp.windows(2).map(|w| 60.0 / ((w[1] - w[0]) * ppr)).collect();
    let nout = (n as f64 * ofs / fs).floor() as usize;
    let grid: Vec<f64> = (0..nout).map(|k| k as f64 / ofs).collect();
    let _ = fpts;
    if fit == "smooth" {
        return rt("tachorpm: the 'smooth' fit is not implemented (its least-squares spline was not recovered); use FitType 'linear'.");
    }
    let r = it
        .call_function("interp1", vec![real_col(tm), real_col(rpm_pts), real_col(grid.clone()), Value::string("linear"), Value::string("extrap")], 1)?
        .remove(0)
        .to_matrix()?;
    let mut rm = Matrix::from_dims(vec![nout, 1], r.data.clone());
    if x.kind == Kind::Single {
        rm.kind = Kind::Single;
    }
    out(vec![Value::Num(rm), real_col(grid), real_col(tp)].into_iter().take(nargout.max(1)).collect())
}
