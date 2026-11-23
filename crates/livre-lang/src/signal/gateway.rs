//! `window` (the window function gateway), `db` (convert to decibels across a load) and
//! `scaleFilterSections` (scale the sections of a cascaded transfer function).

use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::scalar_value;
use crate::value::*;

pub const NAMES: &[&str] = &["window", "db", "scaleFilterSections", "hanning"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "window" => window(it, a, nargout),
        "db" => db(it, a),
        "hanning" => hanning(it, a),
        _ => scale_sections(a),
    }
}

fn window(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 2 {
        return rt("Not enough input arguments.");
    }
    if nargout >= 2 {
        return rt("Output argument \"varargout{2}\" (and possibly others) not assigned a value in the execution with \"window\" function.");
    }
    let rest: Vec<Value> = a[1..].to_vec();
    let r = match &a[0] {
        Value::Func(f) => it.call_value(f, rest, 1)?,
        other => {
            let Ok(n) = other.text() else {
                return rt("Function to evaluate must be represented as a string scalar, character vector, or function_handle object.");
            };
            let class = match &a[1] {
                Value::Num(m) => m.kind.name().to_string(),
                _ => a[1].class_name().to_string(),
            };
            match it.call_function(&n, rest, 1) {
                Ok(v) => v,
                Err(e) if e.message.starts_with("Unrecognized function or variable") || e.message.starts_with("Undefined function") => {
                    return rt(format!("Undefined function '{n}' for input arguments of type '{class}'."));
                }
                Err(e) => return Err(e),
            }
        }
    };
    Ok(Some(r.into_iter().take(1).collect()))
}

fn db(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() > 3 {
        return rt("Too many input arguments.");
    }
    let unit_err = |given: Option<String>| -> Ret {
        let tail = match given {
            Some(s) => format!("The input, '{s}', did not match any of the valid values."),
            None => "The input did not match any of the valid values.".to_string(),
        };
        rt(format!("Expected UNIT to match one of these values:  'voltage', 'power'  {tail}"))
    };
    let is_char = |v: &Value| matches!(v, Value::Num(m) if m.kind == Kind::Char) || matches!(v, Value::Str(_));
    let (unit, load): (Option<&Value>, Option<&Value>) = match a.len() {
        1 => (None, None),
        2 if is_char(&a[1]) => (Some(&a[1]), None),
        2 => (None, Some(&a[1])),
        _ => (Some(&a[1]), Some(&a[2])),
    };
    let mut power = false;
    if let Some(u) = unit {
        if !is_char(u) {
            return unit_err(None);
        }
        let s = u.text().unwrap_or_default().to_lowercase();
        if s.is_empty() || !("voltage".starts_with(&s) || "power".starts_with(&s)) {
            return unit_err(Some(u.text().unwrap_or_default()));
        }
        power = "power".starts_with(&s);
    }
    let x = a[0].clone();
    let cls = |v: &Value| match v {
        Value::Num(m) => m.kind.name().to_string(),
        other => other.class_name().to_string(),
    };
    let undefined = |f: &str, v: &Value| -> Ret { rt(format!("Undefined function '{f}' for input arguments of type '{}'.", cls(v))) };
    let non_numeric = |v: &Value| !matches!(v, Value::Num(_));
    let v = if power {
        if non_numeric(&x) {
            return undefined("log10", &x);
        }
        if matches!(&x, Value::Num(m) if m.kind.is_int() || m.kind == Kind::Logical) {
            return undefined("log10", &x);
        }
        if let Value::Num(m) = &x {
            if m.im.is_none() && m.data.iter().any(|v| *v < 0.0) {
                return rt("Power cannot be negative.");
            }
        }
        x
    } else {
        if non_numeric(&x) {
            return undefined("abs", &x);
        }
        let ax = it.call_function("abs", vec![x.clone()], 1)?.remove(0);
        let sq = it.call_function("power", vec![ax, scalar_value(2.0)], 1)?.remove(0);
        let v = match load {
            Some(r) => {
                if non_numeric(r) {
                    return undefined("rdivide", r);
                }
                it.call_function("rdivide", vec![sq, r.clone()], 1)?.remove(0)
            }
            None => sq,
        };
        if matches!(&x, Value::Num(m) if m.kind.is_int()) {
            return undefined("log10", &x);
        }
        v
    };
    let l = it.call_function("log10", vec![v], 1)?.remove(0);
    Ok(Some(it.call_function("times", vec![scalar_value(10.0), l], 1)?))
}

/// `hanning(N)`: the N points of a Hann window without its two zero end points (`hann(N+2)`
/// inside); `'periodic'` gives the periodic Hann window.
fn hanning(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if a.len() > 2 {
        return rt("Too many input arguments.");
    }
    const TYPES: &str = "Expected N to be one of these types:

double, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64";
    let m = match &a[0] {
        Value::Num(m) if m.kind == Kind::Char => return rt(TYPES),
        Value::Num(m) if m.kind == Kind::Logical => return rt(format!("{TYPES}

Instead its type was logical.")),
        Value::Num(m) => m.clone(),
        _ => return rt(TYPES),
    };
    if m.im.is_some() {
        return rt("Expected N to be real.");
    }
    if m.numel() != 1 {
        return rt("Expected N to be a scalar.");
    }
    let v = m.data[0];
    if !v.is_finite() {
        return rt("Expected N to be finite.");
    }
    if v < 0.0 {
        return rt("Expected N to be nonnegative.");
    }
    if v.fract() != 0.0 {
        it.emit(crate::interp::OutputKind::Warning("Rounding order to nearest integer.".into()));
    }
    let n = v.round();
    if a.len() == 2 {
        // The periodic flag is validated (and applied) by hann itself.
        let probe = it.call_function("hann", vec![scalar_value(3.0), a[1].clone()], 1)?.remove(0);
        let _ = probe;
        let periodic = a[1].text().map(|t| !t.is_empty() && "periodic".starts_with(&t.to_lowercase())).unwrap_or(false);
        if periodic {
            return Ok(Some(it.call_function("hann", vec![scalar_value(n), a[1].clone()], 1)?));
        }
    }
    let w = it.call_function("hann", vec![scalar_value(n + 2.0)], 1)?.remove(0).to_matrix()?;
    let inner: Vec<f64> = w.data[1..1 + n as usize].to_vec();
    Ok(Some(vec![Value::Num(Matrix::col(inner))]))
}

fn scale_sections(a: &[Value]) -> Ret {
    if a.len() < 2 {
        return rt("Invalid argument list. Function requires 1 more input(s).");
    }
    if a.len() > 2 {
        return rt("Too many input arguments.");
    }
    let numeric = |v: &Value| matches!(v, Value::Num(m) if m.kind != Kind::Char && m.kind != Kind::Logical);
    if !numeric(&a[0]) {
        return rt("Invalid argument at position 1. Value must be numeric.");
    }
    if !numeric(&a[1]) {
        return rt("Invalid argument at position 2. Value must be numeric.");
    }
    let (n, g) = (a[0].to_matrix()?, a[1].to_matrix()?);
    if g.dims.len() != 2 || (g.dims[0] != 1 && g.dims[1] != 1) {
        return rt("Invalid argument at position 2. Value must be a 1-by-n vector or an n-by-1 vector.");
    }
    if g.im.as_ref().is_some_and(|im| im.iter().any(|v| *v != 0.0)) {
        return rt("Invalid argument at position 2. Value must be real.");
    }
    // An N-D array is a matrix with the trailing dimensions folded into the columns.
    let rows = n.dims[0];
    let cols = if rows == 0 { n.dims[1..].iter().product::<usize>() } else { n.numel() / rows };
    if n.numel() == 0 && cols == 0 {
        return rt("Input argument must be an L-by-N+1 matrix, where L is the number of sections and N is the order of the section.");
    }
    let l = rows;
    if g.numel() != 1 && g.numel() != l + 1 {
        return rt(format!("Invalid number of scale values. Specify either a scalar or a vector of length equal to the number of sections + 1, {}.", l + 1));
    }
    if n.kind.is_int() && g.numel() != 1 {
        return rt("Integers can only be combined with integers of the same class, or scalar doubles.");
    }
    // The last (or only) scale value is spread evenly over the sections: g^(1/L), principal root.
    let last = g.data[g.numel() - 1];
    let (ur, ui) = if l == 0 {
        (1.0, 0.0)
    } else if last >= 0.0 || last.is_nan() {
        (last.powf(1.0 / l as f64), 0.0)
    } else {
        let m = (-last).powf(1.0 / l as f64);
        let th = std::f64::consts::PI / l as f64;
        (m * th.cos(), m * th.sin())
    };
    let complex = ui != 0.0 || n.im.is_some();
    let total = rows * cols;
    let mut re = vec![0.0; total];
    let mut im = vec![0.0; total];
    for j in 0..cols {
        for i in 0..rows {
            let gi = if g.numel() == 1 { 1.0 } else { g.data[i] };
            let (sr, si) = (gi * ur, gi * ui);
            let (xr, xi) = (n.data[j * rows + i], n.im.as_ref().map(|v| v[j * rows + i]).unwrap_or(0.0));
            re[j * rows + i] = xr * sr - xi * si;
            im[j * rows + i] = xr * si + xi * sr;
            if gi == 0.0 {
                // A zero scale value zeroes the section even where the root is not a number.
                re[j * rows + i] = 0.0;
                im[j * rows + i] = 0.0;
            }
        }
    }
    let dims = vec![rows, cols];
    let mut m = if complex { Matrix::complex(dims, re, im) } else { Matrix::from_dims(dims, re) };
    if n.kind == Kind::Single || g.kind == Kind::Single {
        m.kind = Kind::Single;
        for v in m.data.iter_mut() {
            *v = *v as f32 as f64;
        }
        if let Some(im) = m.im.as_mut() {
            for v in im.iter_mut() {
                *v = *v as f32 as f64;
            }
        }
    } else if n.kind.is_int() {
        m.kind = n.kind;
        for v in m.data.iter_mut() {
            *v = n.kind.cast(*v);
        }
        m.im = None;
    }
    Ok(Some(vec![Value::Num(m)]))
}
