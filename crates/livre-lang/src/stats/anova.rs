//! One-way analysis of variance: `anova1` and `kruskalwallis` (the table and statistics
//! outputs; the figures are not drawn).

use super::dist::{chi2cdf, fcdf};
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["anova1", "kruskalwallis"];

#[derive(Clone, PartialEq)]
enum Label {
    Num(f64),
    Str(String),
}

fn group_labels(v: &Value) -> Res<Vec<Label>> {
    match v {
        Value::Num(m) if m.kind == Kind::Char => {
            let (r, c) = (m.dims[0], m.dims[1]);
            Ok((0..r).map(|i| Label::Str((0..c).map(|j| char::from_u32(m.data[j * r + i] as u32).unwrap_or(' ')).collect::<String>().trim_end().to_string())).collect())
        }
        Value::Num(m) => Ok(m.data.iter().map(|&d| Label::Num(d)).collect()),
        Value::Cell(c) => c
            .data
            .iter()
            .map(|e| match e {
                Value::Num(m) if m.kind == Kind::Char => Ok(Label::Str(m.to_string_lossy())),
                _ => rt("GROUP must be a vector, character array or cell array of strings."),
            })
            .collect(),
        Value::Str(s) => Ok(s.data.iter().map(|d| Label::Str(d.clone().unwrap_or_default())).collect()),
        _ => rt("GROUP must be a vector, character array or cell array of strings."),
    }
}

/// `grp2idx`: group numbers (NaN for missing) and group names.
fn grp2idx(labels: &[Label]) -> (Vec<f64>, Vec<String>) {
    let numeric = labels.iter().all(|l| matches!(l, Label::Num(_)));
    if numeric {
        let mut vals: Vec<f64> = labels.iter().filter_map(|l| if let Label::Num(v) = l { Some(*v) } else { None }).filter(|v| !v.is_nan()).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        vals.dedup();
        let idx = labels
            .iter()
            .map(|l| match l {
                Label::Num(v) if !v.is_nan() => (vals.iter().position(|u| u == v).unwrap() + 1) as f64,
                _ => f64::NAN,
            })
            .collect();
        let names = vals.iter().map(|v| num_name(*v)).collect();
        (idx, names)
    } else {
        let mut names: Vec<String> = Vec::new();
        let mut idx = Vec::new();
        for l in labels {
            let s = match l {
                Label::Str(s) => s.clone(),
                Label::Num(v) => num_name(*v),
            };
            if s.is_empty() {
                idx.push(f64::NAN);
                continue;
            }
            let pos = match names.iter().position(|n| *n == s) {
                Some(p) => p,
                None => {
                    names.push(s);
                    names.len() - 1
                }
            };
            idx.push((pos + 1) as f64);
        }
        (idx, names)
    }
}

fn num_name(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e15 { format!("{}", v as i64) } else { format!("{v}") }
}

fn char_rows(rows: &[String]) -> Value {
    let w = rows.iter().map(|s| s.chars().count()).max().unwrap_or(0);
    let r = rows.len();
    let mut data = vec![32.0; r * w];
    for (i, s) in rows.iter().enumerate() {
        for (j, ch) in s.chars().enumerate() {
            data[j * r + i] = ch as u32 as f64;
        }
    }
    let mut m = Matrix::from_dims(vec![r, w], data);
    m.kind = Kind::Char;
    Value::Num(m)
}

pub fn call(_it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    let mut args: Vec<Value> = a.to_vec();
    let classical = name == "anova1";
    if args.is_empty() {
        return rt("Not enough input arguments.");
    }
    if args.len() > 4 {
        return rt("Too many input arguments.");
    }
    let Value::Num(x0) = args[0].clone() else { return rt("X must be numeric.") };
    let single = is_single(&x0);
    let group_v = args.get(1).cloned();
    let _ = &mut args;
    let mut group: Vec<Label> = match &group_v {
        Some(Value::Num(m)) if m.is_empty() && m.kind != Kind::Char => vec![],
        Some(v) => group_labels(v)?,
        None => vec![],
    };
    let is_vec = x0.dims.len() == 2 && (x0.dims[0] == 1 || x0.dims[1] == 1);
    let mismatch = || -> Ret { rt("X and GROUP must have the same length.") };
    let mut x: Vec<f64> = x0.data.clone();
    let mut grouped = !group.is_empty();
    let mut matrix_dims: Option<(usize, usize)> = None;
    let mut named_cols: Option<Vec<String>> = None;
    if !is_vec {
        let (r, c) = (x0.dims[0], x0.dims[1..].iter().product::<usize>());
        let has_nan = x0.data.iter().any(|v| v.is_nan());
        let unique_groups = {
            let mut seen: Vec<&Label> = Vec::new();
            group.iter().all(|g| {
                if seen.contains(&g) {
                    false
                } else {
                    seen.push(g);
                    true
                }
            })
        };
        let need_vector = has_nan || (!group.is_empty() && !unique_groups);
        if need_vector {
            let gi: Vec<usize> = (0..r * c).map(|i| i / r).collect();
            group = if group.is_empty() {
                gi.iter().map(|&g| Label::Num((g + 1) as f64)).collect()
            } else if group.len() == c {
                gi.iter().map(|&g| group[g].clone()).collect()
            } else {
                return mismatch();
            };
            grouped = true;
        } else if !group.is_empty() {
            if c != group.len() {
                return mismatch();
            }
            let (_, gn) = grp2idx(&group);
            let _ = gn;
            named_cols = Some(group.iter().map(|g| match g { Label::Str(s) => s.clone(), Label::Num(v) => num_name(*v) }).collect());
            grouped = false;
            matrix_dims = Some((r, c));
        } else {
            matrix_dims = Some((r, c));
        }
    } else if !group.is_empty() && group.len() != x.len() {
        return mismatch();
    }
    let (xm, countx, mu, df1, df2, rss, tss, lx, tieadj, gnames): (Vec<f64>, Vec<f64>, f64, f64, f64, f64, f64, usize, f64, Option<Vec<String>>);
    if grouped {
        let (idx0, names) = grp2idx(&group);
        // drop NaN data, then NaN groups
        let mut xs = Vec::new();
        let mut gs = Vec::new();
        for (v, g) in x.iter().zip(&idx0) {
            if !v.is_nan() && !g.is_nan() {
                xs.push(*v);
                gs.push(*g as usize);
            }
        }
        // after dropping, regroup only the groups that remain, as grp2idx would on the kept labels
        let kept_labels: Vec<Label> = group
            .iter()
            .zip(x.iter())
            .filter(|(_, v)| !v.is_nan())
            .map(|(l, _)| l.clone())
            .collect();
        let (idx1, names1) = grp2idx(&kept_labels);
        let mut gs2 = Vec::new();
        let mut xs2 = Vec::new();
        for (v, g) in x.iter().filter(|v| !v.is_nan()).zip(&idx1) {
            if !g.is_nan() {
                xs2.push(*v);
                gs2.push(*g as usize);
            }
        }
        let _ = (gs, xs, names);
        let maxi = names1.len();
        let n_obs = xs2.len();
        let (xr, tadj): (Vec<f64>, f64) = if classical {
            let m = mean_of(&xs2);
            (xs2.iter().map(|v| v - m).collect(), 0.0)
        } else {
            let (r, adj) = super::descr::tiedrank_vec(&xs2, false, false);
            (r, adj[0])
        };
        let mut sums = vec![0.0; maxi];
        let mut cnt = vec![0.0; maxi];
        for (v, &g) in xr.iter().zip(&gs2) {
            sums[g - 1] += v;
            cnt[g - 1] += 1.0;
        }
        let means: Vec<f64> = sums.iter().zip(&cnt).map(|(s, c)| s / c).collect();
        let gm = mean_of(&xr);
        let d1 = cnt.iter().filter(|&&c| c > 0.0).count() as f64 - 1.0;
        let d2 = n_obs as f64 - d1 - 1.0;
        let r_ss: f64 = cnt.iter().zip(&means).map(|(c, m)| if *c > 0.0 { c * (m - gm).powi(2) } else { 0.0 }).sum();
        let t_ss: f64 = xr.iter().map(|v| (v - gm).powi(2)).sum();
        xm = means;
        countx = cnt;
        mu = if classical { mean_of(&xs2) } else { 0.0 };
        df1 = d1;
        df2 = d2;
        rss = r_ss;
        tss = t_ss;
        lx = n_obs;
        tieadj = tadj;
        gnames = Some(names1);
    } else {
        let (r, c) = matrix_dims.unwrap_or((x.len(), 1));
        let n_all = r * c;
        let (xr, tadj): (Vec<f64>, f64) = if classical {
            let m = mean_of(&x);
            (x.iter().map(|v| v - m).collect(), 0.0)
        } else {
            let (rk, adj) = super::descr::tiedrank_vec(&x, false, false);
            (rk, adj[0])
        };
        let cols: Vec<f64> = (0..c).map(|j| mean_of(&xr[j * r..(j + 1) * r])).collect();
        let gm = mean_of(&cols);
        let r_ss = r as f64 * cols.iter().map(|m| (m - gm).powi(2)).sum::<f64>();
        let t_ss: f64 = xr.iter().map(|v| (v - gm).powi(2)).sum();
        xm = cols;
        countx = vec![r as f64; c];
        mu = if classical { mean_of(&x) } else { 0.0 };
        df1 = c as f64 - 1.0;
        df2 = c as f64 * (r as f64 - 1.0);
        rss = r_ss;
        tss = t_ss;
        lx = n_all;
        tieadj = tadj;
        gnames = named_cols.clone();
        let _ = &mut x;
    }
    let sse = tss - rss;
    let mse = if df2 > 0.0 { sse / df2 } else { f64::NAN };
    let (fstat, p);
    if classical {
        if sse != 0.0 {
            fstat = (rss / df1) / mse;
            p = fcdf(fstat, df1, df2, true);
        } else if rss == 0.0 {
            fstat = f64::NAN;
            p = f64::NAN;
        } else {
            fstat = f64::INFINITY;
            p = 0.0;
        }
    } else {
        let lxf = lx as f64;
        let mut f = (12.0 * rss) / (lxf * (lxf + 1.0));
        if tieadj > 0.0 {
            f /= 1.0 - 2.0 * tieadj / (lxf.powi(3) - lxf);
        }
        fstat = f;
        p = chi2cdf(f, df1, true);
    }
    let fin = |v: f64| Value::Num(finish(Matrix::scalar(v), single));
    let mut outs = vec![fin(p)];
    if nargout > 1 {
        let heads: Vec<&str> = if classical { vec!["Source", "SS", "df", "MS", "F", "Prob>F"] } else { vec!["Source", "SS", "df", "MS", "Chi-sq", "Prob>Chi-sq"] };
        let row1 = if grouped { "Groups" } else { "Columns" };
        let rows = [row1, "Error", "Total"];
        let num = |v: f64| if v.is_infinite() { Value::empty() } else { Value::num(if single { v as f32 as f64 } else { v }) };
        let none = Value::empty();
        let table: Vec<Vec<Value>> = vec![
            heads.iter().map(|h| Value::chars(h)).collect(),
            vec![Value::chars(rows[0]), num(rss), num(df1), num(rss / df1), num(fstat), num(p)],
            vec![Value::chars(rows[1]), num(sse), num(df2), num(mse), none.clone(), none.clone()],
            vec![Value::chars(rows[2]), num(tss), num(df1 + df2), none.clone(), none.clone(), none],
        ];
        let mut data = Vec::new();
        for c in 0..6 {
            for r in 0..4 {
                data.push(table[r][c].clone());
            }
        }
        outs.push(Value::Cell(CellArray::new(4, 6, data)));
    }
    if nargout > 2 {
        let mut st = StructArray::scalar();
        let m = xm.len();
        let gn_val = match &gnames {
            Some(names) => Value::Cell(CellArray::col(names.iter().map(|s| Value::chars(s)).collect())),
            None => char_rows(&(1..=m).map(|i| i.to_string()).collect::<Vec<_>>().iter().map(|s| format!("{s:>w$}", w = m.to_string().len())).map(|s| s.trim_start().to_string()).collect::<Vec<_>>()),
        };
        st.set("gnames", gn_val);
        st.set("n", Value::Num(Matrix::from_dims(vec![1, m], countx.clone())));
        if classical {
            st.set("source", Value::chars("anova1"));
            st.set("means", Value::Num(finish(Matrix::from_dims(vec![1, m], xm.iter().map(|v| v + mu).collect()), single)));
            st.set("df", Value::num(df2));
            st.set("s", fin(mse.sqrt()));
        } else {
            st.set("source", Value::chars("kruskalwallis"));
            st.set("meanranks", Value::Num(finish(Matrix::from_dims(vec![1, m], xm.clone()), single)));
            st.set("sumt", fin(2.0 * tieadj));
        }
        outs.push(Value::Struct(st));
    }
    Ok(Some(outs))
}
