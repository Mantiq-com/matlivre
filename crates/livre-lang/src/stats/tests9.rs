//! `hmmtrain`: training of a hidden Markov model by the Baum-Welch or the Viterbi algorithm.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::{Interpreter, OutputKind};
use crate::value::*;

pub const NAMES: &[&str] = &["hmmtrain"];

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "hmmtrain" => hmmtrain(it, a, nargout),
        _ => rt(format!("Undefined function '{name}'.")),
    }
}

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn scalar_of(v: &Value) -> Option<f64> {
    match v {
        Value::Num(m) if m.data.len() == 1 => Some(m.data[0]),
        _ => None,
    }
}

/// The symbol names as strings (numbers print as `%g`-like text).
fn symbol_names(v: &Value) -> Option<Vec<String>> {
    match v {
        Value::Num(m) if m.kind == Kind::Char => {
            if m.dims.len() == 2 && (m.dims[0] == 1) {
                Some(m.data.iter().map(|&c| char::from_u32(c as u32).unwrap_or(' ').to_string()).collect())
            } else {
                let (r, c) = (m.dims[0], m.dims[1]);
                Some((0..r).map(|i| (0..c).map(|j| char::from_u32(m.data[j * r + i] as u32).unwrap_or(' ')).collect::<String>().trim_end().to_string()).collect())
            }
        }
        Value::Num(m) => Some(m.data.iter().map(|d| format!("{d}")).collect()),
        Value::Cell(c) => c.data.iter().map(|e| e.text().ok()).collect(),
        _ => None,
    }
}

fn seq_items(v: &Value) -> Option<Vec<String>> {
    symbol_names(v)
}

fn hmmtrain(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    if a.len() < 3 {
        return rt("Not enough input arguments.");
    }
    let (Value::Num(trm), Value::Num(em)) = (&a[1], &a[2]) else { return rt("Transition and emission guesses must be numeric.") };
    let (mut num_states, chk_tr) = (trm.dims[0], trm.dims[1]);
    if chk_tr != num_states {
        return rt("TRANSITION matrix must be square.");
    }
    let (chk_e, mut num_emis) = (em.dims[0], em.dims[1]);
    if chk_e != num_states {
        return rt("EMISSIONS matrix must have the same number of rows as TRANSITIONS.");
    }
    if num_states == 0 || num_emis == 0 {
        if nargout > 2 {
            return rt("Output argument \"logliks\" (and possibly others) not assigned a value in the execution with \"hmmtrain\" function.");
        }
        return Ok(Some(vec![Value::empty(), Value::empty()]));
    }
    let rest = &a[3..];
    if !rest.is_empty() && a.len() % 2 == 0 {
        return rt("Incorrect number of arguments to hmmtrain.");
    }
    let names = ["symbols", "tolerance", "pseudoemissions", "pseudotransitions", "maxiterations", "verbose", "algorithm", "trtol", "etol"];
    let mut opt: Vec<Option<&Value>> = vec![None; names.len()];
    for pair in rest.chunks(2) {
        let key = match pair[0].text() {
            Ok(k) if pair[0].is_text() => k,
            _ => return rt("Parameter name must be text."),
        };
        let kl = key.to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, nm)| !kl.is_empty() && nm.starts_with(&kl)).map(|(i, _)| i).collect();
        let pick = names.iter().position(|nm| *nm == kl).or(if hits.len() == 1 { Some(hits[0]) } else { None });
        match pick {
            Some(i) => opt[i] = Some(&pair[1]),
            None if hits.len() > 1 => return rt(format!("Ambiguous parameter name: {key}.")),
            None => return rt(format!("Invalid parameter name: {key}.")),
        }
    }
    // the sequences
    let mut seqs: Vec<Vec<f64>> = Vec::new();
    let symbols = opt[0].filter(|v| !is_empty_val(v));
    let sym_names: Option<Vec<String>> = match symbols {
        Some(v) => {
            let nm = symbol_names(v).unwrap_or_default();
            if nm.len() != num_emis {
                return rt(format!("SYMBOLS must be a vector with length equal to the number of emissions ({num_emis})."));
            }
            Some(nm)
        }
        None => None,
    };
    let map_items = |items: &[String]| -> Result<Vec<f64>, crate::error::LangError> {
        let names = sym_names.as_ref().unwrap();
        items
            .iter()
            .map(|s| names.iter().position(|n| n == s).map(|p| (p + 1) as f64).ok_or_else(|| crate::error::LangError::new("At least one entry in SEQ does not appear in SYMBOLS.", 0)))
            .collect()
    };
    match &a[0] {
        Value::Cell(c) => {
            for e in &c.data {
                if sym_names.is_some() {
                    let items = seq_items(e).unwrap_or_default();
                    seqs.push(map_items(&items)?);
                } else if let Value::Num(m) = e {
                    seqs.push(m.data.clone());
                } else {
                    return rt("SEQS must be a numeric matrix or a cell array of sequences.");
                }
            }
        }
        Value::Num(m) if m.kind != Kind::Char || sym_names.is_some() => {
            if matches!(symbols, Some(Value::Cell(_))) && m.kind == Kind::Char {
                return rt("At least one entry in SEQ does not appear in SYMBOLS.");
            }
            if sym_names.is_some() {
                let items: Vec<String> = if m.kind == Kind::Char {
                    let (r, c) = (m.dims[0], m.dims[1]);
                    (0..r * c).map(|k| char::from_u32(m.data[k] as u32).unwrap_or(' ').to_string()).collect()
                } else {
                    m.data.iter().map(|d| format!("{d}")).collect()
                };
                let mapped = map_items(&items)?;
                let (r, c) = (m.dims[0], m.dims[1]);
                for i in 0..r {
                    seqs.push((0..c).map(|j| mapped[j * r + i]).collect());
                }
            } else {
                let (r, c) = (m.dims[0], m.dims[1]);
                for i in 0..r {
                    seqs.push((0..c).map(|j| m.data[j * r + i]).collect());
                }
            }
        }
        _ => return rt("Seqs must be cell array or numerical array."),
    }
    let get_mat = |v: Option<&Value>| -> Option<Matrix> {
        match v {
            Some(Value::Num(m)) if !m.is_empty() => Some(m.clone()),
            _ => None,
        }
    };
    let pseudo_e = get_mat(opt[2]);
    let pseudo_tr = get_mat(opt[3]);
    if let Some(pe) = &pseudo_e {
        if pe.dims[0] < num_states {
            return rt("There are more states in GUESSTR than in PSEUDOE.");
        }
        if pe.dims[1] < num_emis {
            return rt("There are more symbols in SEQ than in PSEUDOE.");
        }
        num_states = pe.dims[0];
        num_emis = pe.dims[1];
    }
    if let Some(pt) = &pseudo_tr {
        if pt.dims[0] != pt.dims[1] {
            return rt("PSEUDOTRANSITIONS matrix must be square.");
        }
        if pt.dims[0] < num_states {
            return rt("There are more states in GUESSTR than in PSEUDOTR.");
        }
        num_states = pt.dims[0];
    }
    let (ns0, ne0) = (trm.dims[0], em.dims[1]);
    if num_states != ns0 || num_emis != ne0 {
        return rt("The pseudocount matrices must have the same size as the guesses.");
    }
    let verbose = match opt[5] {
        Some(v) if v.is_text() => ["on", "true", "yes"].contains(&v.text().unwrap_or_default().to_lowercase().as_str()),
        Some(v) => scalar_of(v).map(|d| d != 0.0).unwrap_or(false),
        None => false,
    };
    let mut baumwelch = true;
    if let Some(v) = opt[6] {
        if !is_empty_val(v) {
            let t = if v.is_text() { v.text().unwrap_or_default().to_lowercase() } else { String::new() };
            let choices = ["baumwelch", "viterbi"];
            let hits: Vec<usize> = choices.iter().enumerate().filter(|(_, c)| !t.is_empty() && c.starts_with(&t)).map(|(i, _)| i).collect();
            if hits.len() != 1 {
                return rt(format!("'{}' is not a valid value for the Algorithm argument.  Valid values are: 'baumwelch', 'viterbi'.", if v.is_text() { v.text().unwrap_or_default() } else { String::new() }));
            }
            baumwelch = hits[0] == 0;
        }
    }
    let tol = opt[1].filter(|v| !is_empty_val(v)).and_then(scalar_of).unwrap_or(1e-6);
    let trtol = opt[7].filter(|v| !is_empty_val(v)).and_then(scalar_of).unwrap_or(tol);
    let etol = opt[8].filter(|v| !is_empty_val(v)).and_then(scalar_of).unwrap_or(tol);
    let maxiter = opt[4].filter(|v| !is_empty_val(v)).and_then(scalar_of).unwrap_or(500.0) as usize;
    let ptr: Vec<f64> = pseudo_tr.as_ref().map(|m| m.data.clone()).unwrap_or_else(|| vec![0.0; num_states * num_states]);
    let pe: Vec<f64> = pseudo_e.as_ref().map(|m| m.data.clone()).unwrap_or_else(|| vec![0.0; num_states * num_emis]);
    let mut guess_tr = trm.data.clone();
    let mut guess_e = em.data.clone();
    let mut tr = vec![0.0; num_states * num_states];
    let mut e = vec![0.0; num_states * num_emis];
    let mut converged = false;
    let mut loglik = 1.0f64;
    let mut logliks: Vec<f64> = Vec::new();
    let tr_at = |m: &Vec<f64>, k: usize, l: usize| m[l * num_states + k];
    let e_at = |m: &Vec<f64>, k: usize, s: usize| m[s * num_states + k];
    let mut iteration = 0;
    while iteration < maxiter {
        iteration += 1;
        let old_ll = loglik;
        loglik = 0.0;
        let old_e = guess_e.clone();
        let old_tr = guess_tr.clone();
        for seq in &seqs {
            let l = seq.len();
            let trv = Value::Num(Matrix::from_dims(vec![num_states, num_states], guess_tr.clone()));
            let ev = Value::Num(Matrix::from_dims(vec![num_states, num_emis], guess_e.clone()));
            let sv = Value::Num(Matrix::from_dims(vec![1, l], seq.clone()));
            if baumwelch {
                let r = it.call_function("hmmdecode", vec![sv, trv, ev], 5)?;
                let log_pseq = r[1].scalar()?;
                loglik += log_pseq;
                let (Value::Num(fs), Value::Num(bs), Value::Num(sc)) = (&r[2], &r[3], &r[4]) else { return rt("hmmdecode returned no value") };
                let ncolf = fs.dims[1];
                let fsv = |k: usize, i: usize| fs.data[i * num_states + k];
                let bsv = |k: usize, i: usize| bs.data[i * num_states + k];
                let _ = ncolf;
                for k in 0..num_states {
                    for ll in 0..num_states {
                        for i in 0..l {
                            // fs(k,i) is the forward value after i symbols (column 0 is the start)
                            let s = seq[i] as usize;
                            let lf = fsv(k, i + 1 - 1).ln();
                            let _ = lf;
                            let term = (fsv(k, i).ln() + tr_at(&guess_tr, k, ll).ln() + e_at(&guess_e, ll, s - 1).ln() + bsv(ll, i + 1).ln()).exp() / sc.data[i + 1];
                            tr[ll * num_states + k] += term;
                        }
                    }
                }
                for k in 0..num_states {
                    for s in 1..=num_emis {
                        let mut acc = 0.0;
                        for pos in 0..l {
                            if seq[pos] as usize == s {
                                acc += (fsv(k, pos + 1).ln() + bsv(k, pos + 1).ln()).exp();
                            }
                        }
                        e[(s - 1) * num_states + k] += acc;
                    }
                }
            } else {
                let r = it.call_function("hmmviterbi", vec![sv.clone(), trv, ev], 2)?;
                loglik += r[1].scalar()?;
                let mut args = vec![sv, r[0].clone(), Value::chars("pseudoe"), Value::Num(Matrix::from_dims(vec![num_states, num_emis], pe.clone())), Value::chars("pseudoTR"), Value::Num(Matrix::from_dims(vec![num_states, num_states], ptr.clone()))];
                let _ = &mut args;
                let r2 = it.call_function("hmmestimate", args, 2)?;
                let (Value::Num(itr), Value::Num(ie)) = (&r2[0], &r2[1]) else { return rt("hmmestimate returned no value") };
                for (d, s) in tr.iter_mut().zip(&itr.data) {
                    *d += if s.is_nan() { 0.0 } else { *s };
                }
                for (d, s) in e.iter_mut().zip(&ie.data) {
                    *d += if s.is_nan() { 0.0 } else { *s };
                }
            }
        }
        let tot_e: Vec<f64> = (0..num_states).map(|k| (0..num_emis).map(|s| e_at(&e, k, s)).sum()).collect();
        let tot_t: Vec<f64> = (0..num_states).map(|k| (0..num_states).map(|ll| tr_at(&tr, k, ll)).sum()).collect();
        for k in 0..num_states {
            for s in 0..num_emis {
                guess_e[s * num_states + k] = e[s * num_states + k] / tot_e[k];
            }
            for ll in 0..num_states {
                guess_tr[ll * num_states + k] = tr[ll * num_states + k] / tot_t[k];
            }
        }
        for k in 0..num_states {
            if tot_t[k] == 0.0 {
                for ll in 0..num_states {
                    guess_tr[ll * num_states + k] = 0.0;
                }
                guess_tr[k * num_states + k] = 1.0;
            }
        }
        guess_tr.iter_mut().for_each(|v| {
            if v.is_nan() {
                *v = 0.0;
            }
        });
        guess_e.iter_mut().for_each(|v| {
            if v.is_nan() {
                *v = 0.0;
            }
        });
        let norm_inf = |a: &Vec<f64>, b: &Vec<f64>| -> f64 {
            (0..num_states).map(|k| (0..a.len() / num_states).map(|c| (a[c * num_states + k] - b[c * num_states + k]).abs()).sum::<f64>()).fold(0.0, f64::max)
        };
        if verbose {
            if iteration == 1 {
                it.emit(OutputKind::Text("Relative Changes in Log Likelihood, Transition Matrix and Emission Matrix\n   Iteration       Log Lik    Transition     Emmission\n".to_string()));
            } else {
                it.emit(OutputKind::Text(format!(
                    "  {:6}      {:12}  {:12}  {:12}\n",
                    iteration,
                    (loglik - old_ll).abs() / (1.0 + old_ll.abs()),
                    norm_inf(&guess_tr, &old_tr) / num_states as f64,
                    norm_inf(&guess_e, &old_e) / num_emis as f64
                )));
            }
        }
        logliks.push(loglik);
        if (loglik - old_ll).abs() / (1.0 + old_ll.abs()) < tol && norm_inf(&guess_tr, &old_tr) / (num_states as f64) < trtol && norm_inf(&guess_e, &old_e) / (num_emis as f64) < etol {
            converged = true;
            break;
        }
        e = if pseudo_e.is_some() { pe.clone() } else { vec![0.0; num_states * num_emis] };
        tr = if pseudo_tr.is_some() { ptr.clone() } else { vec![0.0; num_states * num_states] };
    }
    if !converged {
        it.emit_warning(format!("Algorithm did not converge with tolerance {} in {} iterations.", fmt_short(tol), maxiter), "stats:hmmtrain:NoConvergence".into());
    }
    logliks.retain(|&v| v != 0.0);
    let mut out = vec![Value::Num(Matrix::from_dims(vec![num_states, num_states], guess_tr)), Value::Num(Matrix::from_dims(vec![num_states, num_emis], guess_e))];
    if nargout > 2 {
        let n = logliks.len();
        out.push(Value::Num(Matrix::from_dims(vec![1, n], logliks)));
    }
    Ok(Some(out))
}

/// `num2str`-like short text of a tolerance.
fn fmt_short(v: f64) -> String {
    let s = format!("{:e}", v);
    let (m, e) = s.split_once('e').unwrap();
    let ex: i32 = e.parse().unwrap();
    if (-4..5).contains(&ex) {
        let t = format!("{}", v);
        return t;
    }
    let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.').to_string() } else { m.to_string() };
    format!("{}e{}{:02}", m, if ex < 0 { "-" } else { "+" }, ex.abs())
}
