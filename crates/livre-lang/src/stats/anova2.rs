//! `anova2` (two-way ANOVA with balanced replication) and `friedman`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["anova2", "friedman"];

fn upper(it: &mut Interpreter, f: &str, x: f64, p: &[f64]) -> Res<f64> {
    let mut args = vec![Value::num(x)];
    args.extend(p.iter().map(|&v| Value::num(v)));
    args.push(Value::chars("upper"));
    it.call_function(f, args, 1)?[0].scalar()
}

struct Core {
    css: f64,
    rss: f64,
    iss: f64,
    sse: f64,
    tss: f64,
    df1: f64,
    df2: f64,
    idf: f64,
    edf: f64,
    colf: f64,
    rowf: f64,
    intf: f64,
    p: Vec<f64>,
    mse: f64,
    colmean: Vec<f64>,
    rowmean: Vec<f64>,
    m: usize,
    c: usize,
    r: usize,
    reps: usize,
}

fn core(it: &mut Interpreter, x: &Matrix, reps: usize) -> Res<Core> {
    let (r, c) = (x.dims[0], x.dims[1]);
    let at = |i: usize, j: usize| x.data[j * r + i];
    let m = r / reps;
    let y: Vec<Vec<f64>> = (0..m).map(|i| (0..c).map(|j| (0..reps).map(|k| at(i * reps + k, j)).sum::<f64>() / reps as f64).collect()).collect();
    let colmean: Vec<f64> = (0..c).map(|j| y.iter().map(|row| row[j]).sum::<f64>() / m as f64).collect();
    let rowmean: Vec<f64> = (0..m).map(|i| y[i].iter().sum::<f64>() / c as f64).collect();
    let gm = colmean.iter().sum::<f64>() / c as f64;
    let (df1, df2) = (c as f64 - 1.0, m as f64 - 1.0);
    let (edf, idf) = if reps == 1 { ((c as f64 - 1.0) * (r as f64 - 1.0), 0.0) } else { ((c * m * (reps - 1)) as f64, (c as f64 - 1.0) * (m as f64 - 1.0)) };
    let css = (m * reps) as f64 * colmean.iter().map(|v| (v - gm).powi(2)).sum::<f64>();
    let rss = (c * reps) as f64 * rowmean.iter().map(|v| (v - gm).powi(2)).sum::<f64>();
    let correction = (c * m * reps) as f64 * gm * gm;
    let tss = x.data.iter().map(|v| v * v).sum::<f64>() - correction;
    let iss = reps as f64 * y.iter().map(|row| row.iter().map(|v| v * v).sum::<f64>()).sum::<f64>() - correction - css - rss;
    let sse = if reps == 1 { iss } else { tss - css - rss - iss };
    let (mut colf, mut rowf, mut intf) = (f64::NAN, f64::NAN, f64::NAN);
    let (mut p, mse);
    if sse != 0.0 {
        mse = sse / edf;
        colf = (css / df1) / mse;
        rowf = (rss / df2) / mse;
        p = vec![upper(it, "fcdf", colf, &[df1, edf])?, upper(it, "fcdf", rowf, &[df2, edf])?];
        if reps > 1 {
            intf = (iss / idf) / mse;
            p.push(upper(it, "fcdf", intf, &[idf, edf])?);
        }
    } else {
        mse = if edf > 0.0 { 0.0 } else { f64::NAN };
        let (cp, rp);
        if css == 0.0 {
            colf = f64::NAN;
            cp = f64::NAN;
        } else {
            colf = f64::INFINITY;
            cp = 0.0;
        }
        if rss == 0.0 {
            rowf = f64::NAN;
            rp = f64::NAN;
        } else {
            rowf = f64::INFINITY;
            rp = 0.0;
        }
        p = vec![cp, rp];
        if reps > 1 && iss == 0.0 {
            intf = f64::NAN;
            p.push(f64::NAN);
        } else if reps > 1 {
            intf = f64::INFINITY;
            p.push(0.0);
        }
    }
    Ok(Core { css, rss, iss, sse, tss, df1, df2, idf, edf, colf, rowf, intf, p, mse, colmean, rowmean, m, c, r, reps })
}

fn cell(rows: Vec<Vec<Value>>) -> Value {
    let (nr, nc) = (rows.len(), rows[0].len());
    let mut data = vec![Value::empty(); nr * nc];
    for (i, r) in rows.into_iter().enumerate() {
        for (j, v) in r.into_iter().enumerate() {
            data[j * nr + i] = v;
        }
    }
    Value::Cell(CellArray::from_dims(vec![nr, nc], data))
}

fn s(t: &str) -> Value {
    Value::chars(t)
}

fn n(v: f64) -> Value {
    Value::num(v)
}

fn anova2_table(c: &Core) -> Vec<Vec<Value>> {
    let e = Value::empty;
    let mut rows = vec![vec![s("Source"), s("SS"), s("df"), s("MS"), s("F"), s("Prob>F")]];
    rows.push(vec![s("Columns"), n(c.css), n(c.df1), n(c.css / c.df1), n(c.colf), n(c.p[0])]);
    rows.push(vec![s("Rows"), n(c.rss), n(c.df2), n(c.rss / c.df2), n(c.rowf), n(c.p[1])]);
    if c.reps > 1 {
        rows.push(vec![s("Interaction"), n(c.iss), n(c.idf), n(c.iss / c.idf), n(c.intf), n(c.p[2])]);
    }
    rows.push(vec![s("Error"), n(c.sse), n(c.edf), n(c.sse / c.edf), e(), e()]);
    rows.push(vec![s("Total"), n(c.tss), n((c.r * c.c) as f64 - 1.0), e(), e(), e()]);
    rows
}

fn display_opt(a: &[Value], idx: usize) -> Res<()> {
    if let Some(v) = a.get(idx) {
        let t = v.text().unwrap_or_default();
        if t != "on" && t != "off" {
            return rt("DISPLAYOPT must be 'on' or 'off'.");
        }
    }
    Ok(())
}

fn anova2(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Requires at least one input argument.");
    }
    let Value::Num(x) = &a[0] else { return rt("X must be numeric.") };
    if x.data.iter().any(|v| v.is_nan()) {
        return rt("NaN values in input not allowed.  Use anovan instead.");
    }
    let (r, _c) = (x.dims[0], x.dims[1]);
    let reps = match a.get(1) {
        None => 1.0,
        Some(v) => v.scalar()?,
    };
    if !reps.is_finite() || reps < 1.0 {
        return rt("NaN and Inf not allowed.");
    }
    if reps.fract() != 0.0 || (r as f64 / reps).fract() != 0.0 {
        return rt("The number of rows must be a multiple of reps.");
    }
    let reps = reps as usize;
    let core = core(it, x, reps)?;
    let table = cell(anova2_table(&core));
    let mut st = StructArray::scalar();
    st.set("source", s("anova2"));
    st.set("sigmasq", n(core.mse));
    st.set("colmeans", Value::Num(Matrix::from_dims(vec![1, core.c], core.colmean.clone())));
    st.set("coln", n((core.m * reps) as f64));
    st.set("rowmeans", Value::Num(Matrix::from_dims(vec![1, core.m], core.rowmean.clone())));
    st.set("rown", n((core.c * reps) as f64));
    st.set("inter", Value::bool(reps > 1));
    st.set("pval", n(if reps > 1 { core.p.get(2).copied().unwrap_or(f64::NAN) } else { f64::NAN }));
    st.set("df", n(core.edf));
    let p = Value::Num(Matrix::from_dims(vec![1, core.p.len()], core.p.clone()));
    Ok(Some(vec![p, table, Value::Struct(st)]))
}

fn friedman(it: &mut Interpreter, a: &[Value]) -> Ret {
    if a.is_empty() {
        return rt("Requires at least one input argument.");
    }
    let Value::Num(x) = &a[0] else { return rt("X must be numeric.") };
    let reps_v = match a.get(1) {
        None => 1.0,
        Some(v) => v.scalar()?,
    };
    if x.data.iter().any(|v| v.is_nan()) {
        return rt("NaN values in input not allowed.");
    }
    let (rr, c) = (x.dims[0], x.dims[1]);
    let mut r = rr as f64;
    if reps_v > 1.0 {
        r /= reps_v;
        if r.floor() != r {
            return rt("The number of rows must be a multiple of reps.");
        }
    }
    if r <= 1.0 || c <= 1 {
        return rt("Must have at least two rows and columns.");
    }
    display_opt(a, 2)?;
    if !reps_v.is_finite() || reps_v < 1.0 {
        return rt("NaN and Inf not allowed.");
    }
    if reps_v.fract() != 0.0 {
        return rt("Size arguments must be real integers.");
    }
    let reps = reps_v as usize;
    let r = r as usize;
    let at = |i: usize, j: usize| x.data[j * rr + i];
    let mut ranks = vec![0.0; rr * c];
    let mut sumta = 0.0;
    for j in 0..r {
        let mut vals: Vec<(f64, usize, usize)> = Vec::new();
        for col in 0..c {
            for k in 0..reps {
                vals.push((at(reps * j + k, col), k, col));
            }
        }
        let mut idx: Vec<usize> = (0..vals.len()).collect();
        idx.sort_by(|&p, &q| vals[p].0.partial_cmp(&vals[q].0).unwrap());
        let mut i = 0;
        while i < idx.len() {
            let mut e = i;
            while e + 1 < idx.len() && vals[idx[e + 1]].0 == vals[idx[i]].0 {
                e += 1;
            }
            let avg = (i + e) as f64 / 2.0 + 1.0;
            for t in i..=e {
                let (_, k, col) = vals[idx[t]];
                ranks[col * rr + reps * j + k] = avg;
            }
            let t = (e - i + 1) as f64;
            sumta += t * t * t - t;
            i = e + 1;
        }
    }
    let rm = Matrix::from_dims(vec![rr, c], ranks);
    let core = core(it, &rm, reps)?;
    let cf = c as f64;
    let rf = reps as f64;
    let mut sigmasq = cf * rf * (rf * cf + 1.0) / 12.0;
    if sumta > 0.0 {
        sigmasq -= sumta / (12.0 * r as f64 * (rf * cf - 1.0));
    }
    let mut chistat = core.css;
    if chistat > 0.0 {
        chistat /= sigmasq;
    }
    let p = upper(it, "chi2cdf", chistat, &[cf - 1.0])?;
    let e = Value::empty;
    let mut rows = vec![vec![s("Source"), s("SS"), s("df"), s("MS"), s("Chi-sq"), s("Prob>Chi-sq")]];
    rows.push(vec![s("Columns"), n(core.css), n(core.df1), n(core.css / core.df1), n(chistat), n(p)]);
    if reps > 1 {
        rows.push(vec![s("Interaction"), n(core.iss), n(core.idf), n(core.iss / core.idf), e(), e()]);
    }
    rows.push(vec![s("Error"), n(core.sse), n(core.edf), n(core.sse / core.edf), e(), e()]);
    rows.push(vec![s("Total"), n(core.tss), n((core.r * core.c) as f64 - 1.0), e(), e(), e()]);
    let mut st = StructArray::scalar();
    st.set("source", s("friedman"));
    st.set("n", n(r as f64));
    let means: Vec<f64> = (0..c).map(|j| (0..rr).map(|i| rm.data[j * rr + i]).sum::<f64>() / rr as f64).collect();
    st.set("meanranks", Value::Num(Matrix::from_dims(vec![1, c], means)));
    st.set("sigma", n(sigmasq.sqrt()));
    Ok(Some(vec![n(p), cell(rows), Value::Struct(st)]))
}

fn sg(v: Value, single: bool) -> Value {
    match v {
        Value::Num(m) if single && m.kind == Kind::Double && !m.is_empty() => Value::Num(m.cast_to(Kind::Single)),
        other => other,
    }
}

fn single_cols(v: Value, cols: &[usize], single: bool) -> Value {
    // converts the numeric entries of the given columns of a table cell (and nothing else)
    match v {
        Value::Cell(mut c) if single => {
            let nr = c.dims[0];
            for &j in cols {
                for i in 0..nr {
                    let idx = j * nr + i;
                    c.data[idx] = sg(std::mem::replace(&mut c.data[idx], Value::empty()), true);
                }
            }
            Value::Cell(c)
        }
        other => other,
    }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], _nargout: usize) -> Ret {
    let single = matches!(a.first(), Some(Value::Num(m)) if m.kind == Kind::Single);
    let r = match name {
        "anova2" => anova2(it, a)?,
        _ => friedman(it, a)?,
    };
    Ok(r.map(|mut v| {
        if single {
            v[0] = sg(v[0].clone(), true);
            v[1] = single_cols(v[1].clone(), &[1, 3, 4, 5], true);
            if let Value::Struct(st) = &mut v[2] {
                for fld in ["sigmasq", "colmeans", "rowmeans", "meanranks"] {
                    if let Some(i) = st.field_index(fld) {
                        st.elems[0][i] = sg(st.elems[0][i].clone(), true);
                    }
                }
            }
        }
        v
    }))
}
