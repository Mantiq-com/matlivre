//! `tfridge`: time-frequency ridges by penalized dynamic programming.

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;

pub const NAMES: &[&str] = &["tfridge"];

pub fn call(_it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    if nargout > 3 {
        return rt("Too many output arguments.");
    }
    let tfm = match &a[0] {
        Value::Num(m) if matches!(m.kind, Kind::Double | Kind::Single) => m.clone(),
        _ => return rt("Expected input number 1, TFM, to be one of these types:  single, double"),
    };
    if tfm.dims.len() != 2 || tfm.dims[0] == 1 || tfm.dims[1] == 1 {
        return rt("Expected input number 1, TFM, to be a matrix.");
    }
    let fv = a[1].to_matrix()?;
    let (nr, nc) = (tfm.dims[0], tfm.dims[1]);
    let mut penalty = 0.0;
    let mut nridges = 1usize;
    let mut nbins = 4usize;
    let mut i = 2;
    if let Some(Value::Num(m)) = a.get(2) {
        if m.kind != Kind::Char {
            if m.numel() != 1 {
                return rt("Expected input number 3, PENALTY, to be a scalar.");
            }
            if m.data[0] < 0.0 {
                return rt("Expected input number 3, PENALTY, to be nonnegative.");
            }
            penalty = m.data[0];
            i = 3;
        }
    }
    while i < a.len() {
        let key = a[i].text().unwrap_or_default();
        let Some(v) = a.get(i + 1) else { return rt(format!("No value was given for '{key}'. Name-value pair arguments require a name followed by a value.")) };
        match key.to_lowercase().as_str() {
            "numridges" => {
                let x = v.scalar().unwrap_or(f64::NAN);
                if !(x > 0.0) {
                    return rt("The value of 'NumRidges' is invalid. Expected NR to be positive.");
                }
                if x.fract() != 0.0 {
                    return rt("The value of 'NumRidges' is invalid. Expected NR to be integer-valued.");
                }
                nridges = x as usize;
            }
            "numfrequencybins" => {
                let x = v.scalar().unwrap_or(f64::NAN);
                if !(x > 0.0) {
                    return rt("The value of 'NumFrequencyBins' is invalid. Expected NFB to be positive.");
                }
                if x.fract() != 0.0 {
                    return rt("The value of 'NumFrequencyBins' is invalid. Expected NFB to be integer-valued.");
                }
                nbins = x as usize;
            }
            _ => return rt(format!("'{key}' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.")),
        }
        i += 2;
    }
    let im = tfm.im.clone();
    let mut mag: Vec<f64> = (0..nr * nc).map(|k| match &im {
        Some(v) => tfm.data[k].hypot(v[k]),
        None => tfm.data[k].abs(),
    }).collect();
    let mut iridge = vec![0usize; nc * nridges];
    for rr in 0..nridges {
        let path: Vec<usize> = if penalty == 0.0 {
            (0..nc)
                .map(|c| {
                    let col = &mag[c * nr..(c + 1) * nr];
                    let mut best = 0;
                    for r in 1..nr {
                        if col[r] > col[best] {
                            best = r;
                        }
                    }
                    best
                })
                .collect()
        } else {
            let e = |r: usize, c: usize| -> f64 { -mag[c * nr + r].ln() };
            let mut cost = vec![0.0; nr * nc];
            let mut from = vec![0usize; nr * nc];
            for r in 0..nr {
                cost[r] = e(r, 0);
            }
            for c in 1..nc {
                for r in 0..nr {
                    let mut bk = 0;
                    let mut bv = f64::INFINITY;
                    for k in 0..nr {
                        let d = r as f64 - k as f64;
                        let v = cost[(c - 1) * nr + k] + penalty * d * d;
                        if v < bv {
                            bv = v;
                            bk = k;
                        }
                    }
                    cost[c * nr + r] = e(r, c) + bv;
                    from[c * nr + r] = bk;
                }
            }
            let mut last = 0;
            for r in 1..nr {
                if cost[(nc - 1) * nr + r] < cost[(nc - 1) * nr + last] {
                    last = r;
                }
            }
            let mut path = vec![0usize; nc];
            path[nc - 1] = last;
            for c in (1..nc).rev() {
                path[c - 1] = from[c * nr + path[c]];
            }
            path
        };
        for c in 0..nc {
            iridge[rr * nc + c] = path[c];
            let lo = path[c].saturating_sub(nbins);
            let hi = (path[c] + nbins).min(nr - 1);
            for r in lo..=hi {
                mag[c * nr + r] = 0.0;
            }
        }
    }
    if iridge.iter().any(|&r| r >= fv.numel()) {
        let n = fv.numel();
        return rt(format!("Index exceeds the number of array elements. Index must not exceed {n}."));
    }
    let single = tfm.kind == Kind::Single;
    let mut fr = Matrix::from_dims(vec![nc, nridges], iridge.iter().map(|&r| fv.data[r]).collect());
    let mut ir = Matrix::from_dims(vec![nc, nridges], iridge.iter().map(|&r| r as f64 + 1.0).collect());
    let mut lr = Matrix::from_dims(vec![nc, nridges], (0..nc * nridges).map(|k| (k % nc * nr + iridge[k] + 1) as f64).collect());
    if single {
        ir.kind = Kind::Single;
        lr.kind = Kind::Single;
    }
    out(vec![Value::Num(fr), Value::Num(ir), Value::Num(lr)].into_iter().take(nargout.max(1)).collect())
}
