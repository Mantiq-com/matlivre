//! `designfilt` and the `digitalFilter` object (FIR responses: window, equiripple and least squares).

use crate::builtins::Ret;
use crate::class::{ClassDef, Object};
use crate::display::NumFormat;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::signal::poly::*;
use crate::value::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub const NAMES: &[&str] = &["designfilt", "matlivre_df_method"];

const KNOWN: &[&str] = &[
    "FilterOrder", "CutoffFrequency", "CutoffFrequency1", "CutoffFrequency2", "PassbandFrequency", "PassbandFrequency1", "PassbandFrequency2", "StopbandFrequency",
    "StopbandFrequency1", "StopbandFrequency2", "PassbandRipple", "PassbandRipple1", "PassbandRipple2", "StopbandAttenuation", "StopbandAttenuation1",
    "StopbandAttenuation2", "PassbandWeight", "StopbandWeight", "StopbandWeight1", "StopbandWeight2", "SampleRate", "DesignMethod", "Window", "ScalePassband", "MinOrder", "MinOrder",
    "HalfPowerFrequency", "HalfPowerFrequency1", "HalfPowerFrequency2",
];

/// Properties that are used by the design but not listed among the specifications.
const HIDDEN: &[&str] = &["PassbandWeight", "StopbandWeight", "StopbandWeight1", "StopbandWeight2", "Window", "ScalePassband", "MinOrder"];

const RESPONSES: &[&str] = &["lowpassfir", "highpassfir", "bandpassfir", "bandstopfir", "lowpassiir", "highpassiir", "bandpassiir", "bandstopiir"];

fn class_for(props: &[String]) -> Rc<ClassDef> {
    thread_local! {
        static DEFS: RefCell<HashMap<String, Rc<ClassDef>>> = RefCell::new(HashMap::new());
    }
    let key = props.join(",");
    DEFS.with(|d| {
        if let Some(c) = d.borrow().get(&key) {
            return c.clone();
        }
        let mut src = String::from("classdef digitalFilter\n  properties\n");
        for p in props {
            src += &format!("    {p}\n");
        }
        src += "  end\n  methods\n";
        for m in METHODS {
            src += &format!("    function varargout = {m}(obj, varargin)\n      varargout = matlivre_df_method('{m}', nargout, obj, varargin{{:}});\n    end\n");
        }
        src += "  end\nend\n";
        let ast = crate::class::parse_classdef(&src).expect("built-in class source parses");
        let def = Rc::new(ClassDef::build(ast, &[]));
        d.borrow_mut().insert(key, def.clone());
        def
    })
}

fn short(v: f64) -> String {
    crate::display::disp_value(&Value::num(v), NumFormat::Short).trim().to_string()
}

/// The display of a `digitalFilter` (the part after `d = `).
pub fn object_text(o: &Object) -> Option<String> {
    if o.class.name != "digitalFilter" {
        return None;
    }
    let vals = o.props.borrow().clone();
    let names: Vec<String> = o.class.props.iter().map(|p| p.name.clone()).collect();
    let width = names.iter().map(|n| n.len()).max().unwrap_or(0);
    let line = |name: &str, text: &str| format!("{}{}: {}", " ".repeat(width + 4 - name.len()), name, text);
    let coef = vals[0].to_matrix().ok()?;
    let (r, c) = (coef.dims[0], coef.numel() / coef.dims[0].max(1));
    let prefix = format!("{}{}: ", " ".repeat(width + 4 - "Coefficients".len()), "Coefficients");
    let ctext = if r > 1 {
        format!("[{r}x{c} double]")
    } else {
        let els: Vec<String> = coef.data.iter().map(|v| short(*v)).collect();
        let full = format!("[{}]", els.join(" "));
        if prefix.len() + full.len() <= 179 {
            full
        } else {
            let suffix = format!(" ... ] ({r}x{c} double)");
            let mut shown = String::from("[");
            for (i, e) in els.iter().enumerate() {
                let cand = if i == 0 { format!("{shown}{e}") } else { format!("{shown} {e}") };
                if prefix.len() + cand.len() + suffix.len() > 179 {
                    break;
                }
                shown = cand;
            }
            format!("{shown}{suffix}")
        }
    };
    let mut s = String::from(" digitalFilter with properties:\n\n");
    s += &format!("{prefix}{ctext}\n\n   Specifications:\n");
    for (i, n) in names.iter().enumerate().skip(1) {
        let text = match &vals[i] {
            Value::Num(m) if m.kind == Kind::Char => format!("'{}'", vals[i].text().unwrap_or_default()),
            Value::Num(m) if m.numel() == 1 => short(m.data[0]),
            other => other.text().unwrap_or_default(),
        };
        s += &line(n, &text);
        s.push('\n');
    }
    s += "\n Use fvtool to visualize filter\n Use designfilt to edit filter\n Use filter to filter data";
    Some(s)
}

fn canon(name: &str) -> Option<&'static str> {
    KNOWN.iter().find(|k| k.eq_ignore_ascii_case(name)).copied()
}

fn sets_for(resp: &str) -> Vec<Vec<&'static str>> {
    match resp {
        "lowpassfir" | "highpassfir" => vec![
            vec!["FilterOrder", "HalfPowerFrequency"],
            vec!["FilterOrder", "CutoffFrequency"],
            vec!["FilterOrder", "CutoffFrequency", "PassbandRipple", "StopbandAttenuation"],
            vec!["FilterOrder", "PassbandFrequency", "StopbandFrequency"],
            vec!["PassbandFrequency", "StopbandFrequency", "PassbandRipple", "StopbandAttenuation"],
        ],
        "bandpassfir" => vec![
            vec!["FilterOrder", "CutoffFrequency1", "CutoffFrequency2"],
            vec!["FilterOrder", "StopbandFrequency1", "PassbandFrequency1", "PassbandFrequency2", "StopbandFrequency2"],
            vec!["StopbandFrequency1", "PassbandFrequency1", "PassbandFrequency2", "StopbandFrequency2", "StopbandAttenuation1", "PassbandRipple", "StopbandAttenuation2"],
        ],
        _ => vec![
            vec!["FilterOrder", "CutoffFrequency1", "CutoffFrequency2"],
            vec!["FilterOrder", "PassbandFrequency1", "StopbandFrequency1", "StopbandFrequency2", "PassbandFrequency2"],
            vec!["PassbandFrequency1", "StopbandFrequency1", "StopbandFrequency2", "PassbandFrequency2", "PassbandRipple1", "StopbandAttenuation", "PassbandRipple2"],
        ],
    }
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if name == "matlivre_df_method" {
        // The class methods take their outputs from a cell (varargout).
        let r = method(it, a)?;
        return out(vec![Value::Cell(CellArray::row(r.unwrap_or_default()))]);
    }
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    if nargout > 1 {
        return rt("Too many output arguments.");
    }
    let resp = match &a[0] {
        Value::Num(m) if m.kind == Kind::Char => a[0].text().unwrap_or_default(),
        Value::Object(_) => return rt("designfilt: editing a digitalFilter is not implemented."),
        _ => return rt("Filter response is not valid."),
    };
    let resp_l = resp.to_lowercase();
    if a.len() == 1 {
        if RESPONSES.contains(&resp_l.as_str()) {
            return rt("Filter Design Assistant requires a MATLAB desktop. Call designfilt with a complete set of name-value pairs to design your filter.");
        }
        return rt("Filter response is not valid.");
    }
    if (a.len() - 1) % 2 != 0 {
        return rt("All inputs except the filter response must be in the form of name-value pairs.");
    }
    if !RESPONSES.contains(&resp_l.as_str()) {
        let other = ["lowpassiir", "highpassiir", "bandpassiir", "bandstopiir", "differentiatorfir", "hilbertfir", "arbmagfir", "notchiir", "peakiir"];
        if other.contains(&resp_l.as_str()) {
            return rt(format!("designfilt: the response '{resp}' is not implemented."));
        }
        return rt("Filter response is not valid.");
    }
    let mut given: Vec<(&'static str, Value)> = vec![];
    let mut i = 1;
    while i < a.len() {
        let key = a[i].text().unwrap_or_default();
        let Some(c) = canon(&key) else {
            let have: Vec<&str> = (1..a.len()).step_by(2).filter_map(|k| canon(&a[k].text().unwrap_or_default())).collect();
            return rt(format!("Parameter {key} is not recognized, and you have specified too few parameters for '{resp_l}'.\nThe following are valid parameter sets that are close to your inputs:\n{}\n", sets_for(&resp_l).iter().filter(|s| have.iter().all(|n| s.contains(n))).map(|s| format!("  - {}", s.join(", "))).collect::<Vec<_>>().join("\n")));
        };
        given.push((c, a[i + 1].clone()));
        i += 2;
    }
    if resp_l.ends_with("iir") {
        return design_iir(it, &resp_l, &given);
    }
    let has = |n: &str| given.iter().any(|(k, _)| *k == n);
    let get = |n: &str| -> Option<Value> { given.iter().find(|(k, _)| *k == n).map(|(_, v)| v.clone()) };
    let spec_names: Vec<&str> = given.iter().map(|(k, _)| *k).filter(|k| !HIDDEN.contains(k) && *k != "SampleRate" && *k != "DesignMethod").collect();
    // Which parameter set is it?
    let sets = sets_for(&resp_l);
    let found = sets.iter().position(|s| s.len() == spec_names.len() && s.iter().all(|n| spec_names.contains(n)));
    let Some(si) = found else {
        let close: Vec<&Vec<&str>> = sets.iter().filter(|s| spec_names.iter().all(|n| s.contains(n))).collect();
        let list = close.iter().map(|s| format!("  - {}", s.join(", "))).collect::<Vec<_>>().join("\n");
        return rt(format!("You have specified too few parameters for '{resp_l}'.\nThe following are valid parameter sets that are close to your inputs:\n{list}\n"));
    };
    let set = &sets[si];
    if set.contains(&"HalfPowerFrequency") || set.contains(&"CutoffFrequency") && set.contains(&"PassbandRipple") {
        return rt("designfilt: the maxflat and kaiserwin designs are not implemented.");
    }
    let num = |n: &str| -> R<f64> {
        let v = get(n).unwrap().to_matrix()?;
        if v.numel() != 1 {
            return rt(format!("Expected input to be a scalar."));
        }
        Ok(v.data[0])
    };
    let fs = if has("SampleRate") { num("SampleRate")? } else { 2.0 };
    if !(fs > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    let method_set = if set.contains(&"CutoffFrequency") || set.contains(&"CutoffFrequency1") { "window" } else { "equiripple" };
    let minorder = !set.contains(&"FilterOrder");
    let method = match get("DesignMethod") {
        Some(v) => {
            let m = v.text().unwrap_or_default().to_lowercase();
            let ok = match method_set {
                "window" => m == "window",
                _ => m == "equiripple" || (m == "ls" && !minorder) || (m == "kaiserwin" && minorder),
            };
            if !ok {
                return rt(format!("'{}' is an invalid design method for {resp_l} designs with {} specifications.", v.text().unwrap_or_default(), set.join(", ")));
            }
            m
        }
        None => method_set.to_string(),
    };
    let mut order = 0.0;
    if !minorder {
        order = num("FilterOrder")?;
        if order.fract() != 0.0 {
            return rt("Expected input to be integer-valued.");
        }
        if order < 1.0 {
            return rt("Expected input to be positive.");
        }
    }
    let nyq = fs / 2.0;
    let mut freqs: HashMap<&str, f64> = HashMap::new();
    for k in set.iter() {
        if k.contains("Frequency") {
            let f = num(k)?;
            if f < 0.0 || f > nyq {
                return rt(if has("SampleRate") { format!("Frequency specifications must be between 0 and {}.", short(nyq)) } else { "Frequency specifications must be between 0 and 1.".to_string() });
            }
            freqs.insert(k, f / nyq);
        }
    }
    let fq = |k: &str| freqs.get(k).copied().unwrap_or(0.0);
    // The edges must increase in the order of the response.
    let order_names: Vec<&str> = match resp_l.as_str() {
        "lowpassfir" => if method == "window" { vec!["CutoffFrequency"] } else { vec!["PassbandFrequency", "StopbandFrequency"] },
        "highpassfir" => if method == "window" { vec!["CutoffFrequency"] } else { vec!["StopbandFrequency", "PassbandFrequency"] },
        "bandpassfir" => if method == "window" { vec!["CutoffFrequency1", "CutoffFrequency2"] } else { vec!["StopbandFrequency1", "PassbandFrequency1", "PassbandFrequency2", "StopbandFrequency2"] },
        _ => if method == "window" { vec!["CutoffFrequency1", "CutoffFrequency2"] } else { vec!["PassbandFrequency1", "StopbandFrequency1", "StopbandFrequency2", "PassbandFrequency2"] },
    };
    if order_names.windows(2).any(|w| !(fq(w[0]) < fq(w[1]))) {
        let shown: Vec<&str> = if order_names.len() == 2 { order_names.clone() } else { order_names.clone() };
        return rt(format!("The frequency specifications {} must have increasing values.", shown.join(", ")));
    }
    // Edges, desired amplitudes and weights per band.
    let (edges, amps): (Vec<f64>, Vec<f64>) = match resp_l.as_str() {
        "lowpassfir" => (vec![0.0, fq("PassbandFrequency"), fq("StopbandFrequency"), 1.0], vec![1.0, 1.0, 0.0, 0.0]),
        "highpassfir" => (vec![0.0, fq("StopbandFrequency"), fq("PassbandFrequency"), 1.0], vec![0.0, 0.0, 1.0, 1.0]),
        "bandpassfir" => (vec![0.0, fq("StopbandFrequency1"), fq("PassbandFrequency1"), fq("PassbandFrequency2"), fq("StopbandFrequency2"), 1.0], vec![0.0, 0.0, 1.0, 1.0, 0.0, 0.0]),
        _ => (vec![0.0, fq("PassbandFrequency1"), fq("StopbandFrequency1"), fq("StopbandFrequency2"), fq("PassbandFrequency2"), 1.0], vec![1.0, 1.0, 0.0, 0.0, 1.0, 1.0]),
    };
    let wt = |n: &str| -> f64 { get(n).and_then(|v| v.scalar().ok()).unwrap_or(1.0) };
    let weights: Vec<f64> = match resp_l.as_str() {
        "lowpassfir" => vec![wt("PassbandWeight"), wt("StopbandWeight")],
        "highpassfir" => vec![wt("StopbandWeight"), wt("PassbandWeight")],
        "bandpassfir" => vec![wt("StopbandWeight1"), wt("PassbandWeight"), wt("StopbandWeight2")],
        _ => vec![wt("PassbandWeight"), wt("StopbandWeight"), wt("PassbandWeight")],
    };
    let coef: Matrix = if method == "window" {
        let wn: Vec<f64> = match resp_l.as_str() {
            "lowpassfir" | "highpassfir" => vec![num("CutoffFrequency")? / nyq],
            _ => vec![num("CutoffFrequency1")? / nyq, num("CutoffFrequency2")? / nyq],
        };
        let kind = match resp_l.as_str() {
            "highpassfir" => "high",
            "bandstopfir" => "stop",
            "bandpassfir" => "bandpass",
            _ => "low",
        };
        let n = order as usize;
        if n % 2 == 1 && resp_l == "bandstopfir" {
            return rt("designfilt: odd-order bandstop window designs are not implemented.");
        }
        let win = match get("Window") {
            Some(v) if matches!(&v, Value::Num(m) if m.kind == Kind::Char) => {
                let name = v.text().unwrap_or_default();
                it.call_function(&name, vec![scalar_value(n as f64 + 1.0)], 1)?.remove(0)
            }
            Some(v) => v,
            None => it.call_function("hamming", vec![scalar_value(n as f64 + 1.0)], 1)?.remove(0),
        };
        let scale = get("ScalePassband").and_then(|v| v.scalar().ok()).map(|v| v != 0.0).unwrap_or(true);
        if n % 2 == 1 && resp_l == "highpassfir" {
            // An odd order highpass is the antisymmetric type IV filter: a lowpass of cutoff 1 - wc, modulated by -(-1)^k.
            let r = it.call_function("fir1", vec![scalar_value(n as f64), scalar_value(1.0 - wn[0]), Value::string("low"), win, Value::string(if scale { "scale" } else { "noscale" })], 1)?.remove(0).to_matrix()?;
            let data: Vec<f64> = r.data.iter().enumerate().map(|(k, v)| if k % 2 == 0 { -v } else { *v }).collect();
            Matrix::from_dims(vec![1, data.len()], data)
        } else {
            let wnv = if wn.len() == 1 { scalar_value(wn[0]) } else { real_row(wn) };
            let r = it.call_function("fir1", vec![scalar_value(n as f64), wnv, Value::string(kind), win, Value::string(if scale { "scale" } else { "noscale" })], 1)?.remove(0);
            r.to_matrix()?
        }
    } else if !minorder {
        let n = order as usize;
        let name = if method == "ls" { "firls" } else { "firpm" };
        let wv = {
            // Bands come in pairs: one weight per band.
            let mut w = vec![];
            for k in 0..edges.len() / 2 {
                w.push(weights[(k).min(weights.len() - 1)]);
            }
            w
        };
        it.call_function(name, vec![scalar_value(n as f64), real_row(edges.clone()), real_row(amps.clone()), real_row(wv)], 1)?.remove(0).to_matrix()?
    } else {
        // Minimum order equiripple: firpmord from the ripples and attenuations.
        let dev = |rp: f64| (10f64.powf(rp / 20.0) - 1.0) / (10f64.powf(rp / 20.0) + 1.0);
        let att = |a: f64| 10f64.powf(-a / 20.0);
        let (fe, av, dv): (Vec<f64>, Vec<f64>, Vec<f64>) = match resp_l.as_str() {
            "lowpassfir" => (vec![fq("PassbandFrequency"), fq("StopbandFrequency")], vec![1.0, 0.0], vec![dev(num("PassbandRipple")?), att(num("StopbandAttenuation")?)]),
            "highpassfir" => (vec![fq("StopbandFrequency"), fq("PassbandFrequency")], vec![0.0, 1.0], vec![att(num("StopbandAttenuation")?), dev(num("PassbandRipple")?)]),
            "bandpassfir" => (
                vec![fq("StopbandFrequency1"), fq("PassbandFrequency1"), fq("PassbandFrequency2"), fq("StopbandFrequency2")],
                vec![0.0, 1.0, 0.0],
                vec![att(num("StopbandAttenuation1")?), dev(num("PassbandRipple")?), att(num("StopbandAttenuation2")?)],
            ),
            _ => (
                vec![fq("PassbandFrequency1"), fq("StopbandFrequency1"), fq("StopbandFrequency2"), fq("PassbandFrequency2")],
                vec![1.0, 0.0, 1.0],
                vec![dev(num("PassbandRipple1")?), att(num("StopbandAttenuation")?), dev(num("PassbandRipple2")?)],
            ),
        };
        let fe_hz: Vec<f64> = fe.iter().map(|v| v * nyq).collect();
        let r = it.call_function("firpmord", vec![real_row(fe_hz), real_row(av.clone()), real_row(dv.clone()), scalar_value(fs)], 4)?;
        let n0 = r[0].scalar().unwrap_or(1.0) as i64;
        let fo = r[1].to_matrix()?.data;
        let wgt: Vec<f64> = dv.iter().map(|d| 1.0 / d).collect();
        // The smallest order whose design meets every deviation on a dense grid.
        let grid: Vec<f64> = (0..2048).map(|k| k as f64 / 2047.0).collect();
        let meets = |h: &[f64]| -> bool {
            for (bi, d) in dv.iter().enumerate() {
                let (lo, hi) = (fo[2 * bi], fo[2 * bi + 1]);
                let (mut mx, mut mn) = (0.0f64, f64::INFINITY);
                for &g in grid.iter().filter(|g| **g >= lo && **g <= hi) {
                    let (mut re, mut im) = (0.0, 0.0);
                    for (k, c) in h.iter().enumerate() {
                        let ph = -std::f64::consts::PI * g * k as f64;
                        re += c * ph.cos();
                        im += c * ph.sin();
                    }
                    let m = re.hypot(im);
                    mx = mx.max(m);
                    mn = mn.min(m);
                }
                if av[bi] > 0.5 {
                    // The passband ripple is the peak-to-peak variation, wherever it sits.
                    if mn.is_finite() && mx / mn > (1.0 + d) / (1.0 - d) {
                        return false;
                    }
                } else if mx > *d {
                    return false;
                }
            }
            true
        };
        let mut found: Option<Matrix> = None;
        let mut n = (n0 - 8).max(1);
        if method == "kaiserwin" {
            let stop_dev = dv.iter().cloned().fold(f64::INFINITY, f64::min);
            let a_att = -20.0 * stop_dev.log10();
            let beta = if a_att > 50.0 { 0.1102 * (a_att - 8.7) } else if a_att >= 21.0 { 0.5842 * (a_att - 21.0).powf(0.4) + 0.07886 * (a_att - 21.0) } else { 0.0 };
            let (wn, ftype): (Value, &str) = match resp_l.as_str() {
                "lowpassfir" => (scalar_value((fe[0] + fe[1]) / 2.0), "low"),
                "highpassfir" => (scalar_value((fe[0] + fe[1]) / 2.0), "high"),
                "bandpassfir" => (real_row(vec![(fe[0] + fe[1]) / 2.0, (fe[2] + fe[3]) / 2.0]), "bandpass"),
                _ => (real_row(vec![(fe[0] + fe[1]) / 2.0, (fe[2] + fe[3]) / 2.0]), "stop"),
            };
            let even_only = resp_l == "highpassfir" || resp_l == "bandstopfir" || get("MinOrder").and_then(|v| v.text().ok()).map(|t| t.to_lowercase() == "even").unwrap_or(false);
            let kr = it.call_function("kaiserord", vec![real_row(fe.clone()), real_row(av.clone()), real_row(dv.clone()), scalar_value(2.0)], 1)?;
            let mut n = (kr[0].scalar().unwrap_or(2.0) as i64).max(2);
            if even_only && n % 2 == 1 {
                n += 1;
            }
            while n <= n0 + 400 {
                let win = it.call_function("kaiser", vec![scalar_value(n as f64 + 1.0), scalar_value(beta)], 1)?.remove(0);
                let h = it.call_function("fir1", vec![scalar_value(n as f64), wn.clone(), Value::string(ftype), win, Value::string("scale")], 1);
                if let Ok(mut hv) = h {
                    let hm = hv.remove(0).to_matrix()?;
                    if meets(&hm.data) {
                        found = Some(hm);
                        break;
                    }
                }
                n += if even_only { 2 } else { 1 };
            }
        }
        while found.is_none() && method != "kaiserwin" && n <= n0 + 40 {
            let h = it.call_function("firpm", vec![scalar_value(n as f64), r[1].clone(), r[2].clone(), real_row(wgt.clone())], 1);
            if let Ok(mut hv) = h {
                let hm = hv.remove(0).to_matrix()?;
                if meets(&hm.data) {
                    found = Some(hm);
                    break;
                }
            }
            n += 1;
        }
        match found {
            Some(m) => m,
            None => it.call_function("firpm", vec![r[0].clone(), r[1].clone(), r[2].clone(), r[3].clone()], 1)?.remove(0).to_matrix()?,
        }
    };
    let n_coef = coef.numel();
    let coef = Matrix::from_dims(vec![1, n_coef], coef.data.clone());
    // The object: Coefficients, FrequencyResponse, ImpulseResponse, SampleRate, the specifications, DesignMethod.
    let mut names: Vec<String> = vec!["Coefficients".into(), "FrequencyResponse".into(), "ImpulseResponse".into(), "SampleRate".into()];
    let mut vals: Vec<Value> = vec![
        Value::Num(coef),
        Value::chars(resp_l.trim_end_matches("fir").trim_end_matches("iir")),
        Value::chars("fir"),
        scalar_value(fs),
    ];
    for (k, v) in given.iter() {
        if HIDDEN.contains(k) || *k == "SampleRate" || *k == "DesignMethod" {
            continue;
        }
        names.push(k.to_string());
        vals.push(v.clone());
    }
    names.push("DesignMethod".into());
    vals.push(Value::chars(&method));
    let o = Object { class: class_for(&names), props: Rc::new(RefCell::new(vals)) };
    out(vec![Value::Object(o)])
}


/// The IIR responses: butter, cheby1, cheby2 and ellip, for a given order or the minimum order.
fn design_iir(it: &mut Interpreter, resp: &str, given: &[(&'static str, Value)]) -> Ret {
    let has = |n: &str| given.iter().any(|(k, _)| *k == n);
    let get = |n: &str| -> Option<Value> { given.iter().find(|(k, _)| *k == n).map(|(_, v)| v.clone()) };
    let spec: Vec<&str> = given.iter().map(|(k, _)| *k).filter(|k| !HIDDEN.contains(k) && *k != "SampleRate" && *k != "DesignMethod").collect();
    let kind = &resp[..resp.len() - 3];
    let band = kind == "bandpass" || kind == "bandstop";
    let n = |name: &str| -> R<f64> {
        let v = get(name).unwrap().to_matrix()?;
        if v.numel() != 1 {
            return rt("Expected input to be a scalar.");
        }
        Ok(v.data[0])
    };
    let fs = if has("SampleRate") { n("SampleRate")? } else { 2.0 };
    if !(fs > 0.0) {
        return rt("Expected Fs to be positive.");
    }
    let nyq = fs / 2.0;
    let sfx = |base: &str| -> Vec<String> {
        if band {
            vec![format!("{base}1"), format!("{base}2")]
        } else {
            vec![base.to_string()]
        }
    };
    let all_in = |names: &[String]| names.iter().all(|nm| has(nm));
    let fixed = has("FilterOrder");
    let hp = sfx("HalfPowerFrequency");
    let pf = sfx("PassbandFrequency");
    let sf = sfx("StopbandFrequency");
    let mut want: Vec<String> = vec![];
    let method_default: &str;
    if fixed {
        want.push("FilterOrder".into());
        if all_in(&hp) {
            want.extend(hp.clone());
            method_default = "butter";
        } else if all_in(&pf) && has("PassbandRipple") && has("StopbandAttenuation") {
            want.extend(pf.clone());
            want.push("PassbandRipple".into());
            want.push("StopbandAttenuation".into());
            method_default = "ellip";
        } else if all_in(&pf) && has("PassbandRipple") {
            want.extend(pf.clone());
            want.push("PassbandRipple".into());
            method_default = "cheby1";
        } else if all_in(&sf) && has("StopbandAttenuation") {
            want.extend(sf.clone());
            want.push("StopbandAttenuation".into());
            method_default = "cheby2";
        } else {
            let sets: Vec<Vec<String>> = {
                let mut v = vec![];
                let mut a = vec!["FilterOrder".to_string()];
                a.extend(hp.clone());
                v.push(a);
                let mut b2 = vec!["FilterOrder".to_string()];
                b2.extend(pf.clone());
                b2.push("PassbandRipple".into());
                v.push(b2.clone());
                let mut c = b2.clone();
                c.push("StopbandAttenuation".into());
                v.push(c);
                let mut d = vec!["FilterOrder".to_string()];
                d.extend(sf.clone());
                d.push("StopbandAttenuation".into());
                v.push(d);
                v
            };
            let list = sets.iter().filter(|s| spec.iter().all(|n| s.iter().any(|x| x == *n))).map(|s| format!("  - {}", s.join(", "))).collect::<Vec<_>>().join("\n");
            return rt(format!("You have specified too few parameters for '{resp}'.\nThe following are valid parameter sets that are close to your inputs:\n{list}\n"));
        }
        if band && has("PassbandRipple") && has("StopbandAttenuation") {
            return rt(format!("You have specified too many parameters for '{resp}'. A solution is to remove the following parameter from your inputs:\n\n  - StopbandAttenuation\n"));
        }
    } else {
        method_default = "butter";
        let names: Vec<&str> = match kind {
            "lowpass" | "highpass" => vec!["PassbandFrequency", "StopbandFrequency", "PassbandRipple", "StopbandAttenuation"],
            "bandpass" => vec!["StopbandFrequency1", "PassbandFrequency1", "PassbandFrequency2", "StopbandFrequency2", "StopbandAttenuation1", "PassbandRipple", "StopbandAttenuation2"],
            _ => vec!["PassbandFrequency1", "StopbandFrequency1", "StopbandFrequency2", "PassbandFrequency2", "PassbandRipple1", "StopbandAttenuation", "PassbandRipple2"],
        };
        if !names.iter().all(|nm| has(nm)) {
            return rt(format!("You have specified too few parameters for '{resp}'."));
        }
        want.extend(names.iter().map(|s| s.to_string()));
    }
    if let Some(extra) = spec.iter().find(|k| !want.iter().any(|w| w == **k)) {
        return rt(format!("You have specified too many parameters for '{resp}'. A solution is to remove the following parameter from your inputs:\n\n  - {extra}\n"));
    }
    let method = match get("DesignMethod") {
        Some(v) => {
            let m = v.text().unwrap_or_default().to_lowercase();
            let ok = if fixed { m == method_default } else { ["butter", "cheby1", "cheby2", "ellip"].contains(&m.as_str()) };
            if !ok {
                return rt(format!("'{}' is an invalid design method for {resp} designs with {} specifications.", v.text().unwrap_or_default(), want.join(", ")));
            }
            m
        }
        None => method_default.to_string(),
    };
    let mut fr: HashMap<String, f64> = HashMap::new();
    for nm in &want {
        if nm.contains("Frequency") {
            let f = n(nm)?;
            if f < 0.0 || f > nyq {
                return rt(if has("SampleRate") { format!("Frequency specifications must be between 0 and {}.", short(nyq)) } else { "Frequency specifications must be between 0 and 1.".to_string() });
            }
            fr.insert(nm.clone(), f / nyq);
        }
    }
    let f = |nm: &str| fr[nm];
    if band && fixed && method == "butter" && !(f("HalfPowerFrequency1") < f("HalfPowerFrequency2")) {
        return rt("Expected Wcutoff to be positive.");
    }
    let vecv = |names: &[String]| -> Value {
        if names.len() == 1 {
            scalar_value(f(&names[0]))
        } else {
            real_row(names.iter().map(|x| f(x)).collect())
        }
    };
    let ty = match kind {
        "highpass" => Some("high"),
        "bandstop" => Some("stop"),
        _ => None,
    };
    let order: f64;
    let wn_val: Value;
    let rp: f64;
    let rs: f64;
    if fixed {
        let ord = n("FilterOrder")?;
        if ord.fract() != 0.0 {
            return rt("Expected input to be integer-valued.");
        }
        if ord < 1.0 {
            return rt("Expected input to be positive.");
        }
        let edges = if method == "butter" { hp.clone() } else if method == "cheby2" { sf.clone() } else { pf.clone() };
        wn_val = vecv(&edges);
        order = if band { ord / 2.0 } else { ord };
        rp = if has("PassbandRipple") { n("PassbandRipple")? } else { 0.0 };
        rs = if has("StopbandAttenuation") { n("StopbandAttenuation")? } else { 0.0 };
    } else {
        let (wp, ws, rpv, rsv): (Value, Value, f64, f64) = match kind {
            "lowpass" | "highpass" => (scalar_value(f("PassbandFrequency")), scalar_value(f("StopbandFrequency")), n("PassbandRipple")?, n("StopbandAttenuation")?),
            "bandpass" => (
                real_row(vec![f("PassbandFrequency1"), f("PassbandFrequency2")]),
                real_row(vec![f("StopbandFrequency1"), f("StopbandFrequency2")]),
                n("PassbandRipple")?,
                n("StopbandAttenuation1")?.max(n("StopbandAttenuation2")?),
            ),
            _ => (
                real_row(vec![f("PassbandFrequency1"), f("PassbandFrequency2")]),
                real_row(vec![f("StopbandFrequency1"), f("StopbandFrequency2")]),
                n("PassbandRipple1")?.min(n("PassbandRipple2")?),
                n("StopbandAttenuation")?,
            ),
        };
        rp = rpv;
        rs = rsv;
        let ordname = match method.as_str() {
            "butter" => "buttord",
            "cheby1" => "cheb1ord",
            "cheby2" => "cheb2ord",
            _ => "ellipord",
        };
        let r = it.call_function(ordname, vec![wp, ws, scalar_value(rp), scalar_value(rs)], 2)?;
        order = r[0].scalar().unwrap_or(1.0);
        wn_val = r[1].clone();
    }
    let mut args: Vec<Value> = vec![scalar_value(order)];
    match method.as_str() {
        "butter" => {}
        "cheby1" => args.push(scalar_value(rp)),
        "cheby2" => args.push(scalar_value(rs)),
        _ => {
            args.push(scalar_value(rp));
            args.push(scalar_value(rs));
        }
    }
    args.push(wn_val);
    if let Some(t) = ty {
        args.push(Value::string(t));
    }
    let zpk = it.call_function(&method, args, 3)?;
    let sos_out = it.call_function("zp2sos", vec![zpk[0].clone(), zpk[1].clone(), zpk[2].clone(), Value::string("down"), Value::string("inf")], 2)?;
    let mut sos = sos_out[0].to_matrix()?;
    let g = sos_out[1].scalar().unwrap_or(1.0);
    let l = sos.dims[0];
    // The real zeros of a Butterworth or Chebyshev I bandpass sit at +1 and -1: every section takes one of each.
    if kind == "bandpass" && (method == "butter" || method == "cheby1") {
        for i in 0..l {
            let lead = sos.data[i];
            sos.data[l + i] = 0.0;
            sos.data[2 * l + i] = -lead;
        }
    }
    // An even-order elliptic design puts the sections of the poles nearest the unit circle last.
    if method == "ellip" && order % 2.0 == 0.0 {
        let m = if band { 2 } else { 1 }.min(l - 1);
        for c in 0..6 {
            sos.data[c * l..(c + 1) * l].rotate_left(m);
        }
    }
    // The leftover gain goes on the last section.
    for c in 0..3 {
        sos.data[c * l + (l - 1)] *= g;
    }
    let mut names: Vec<String> = vec!["Coefficients".into(), "FrequencyResponse".into(), "ImpulseResponse".into(), "SampleRate".into()];
    let mut vals: Vec<Value> = vec![Value::Num(sos), Value::chars(kind), Value::chars("iir"), scalar_value(fs)];
    for (k, v) in given.iter() {
        if HIDDEN.contains(k) || *k == "SampleRate" || *k == "DesignMethod" {
            continue;
        }
        names.push(k.to_string());
        vals.push(v.clone());
    }
    names.push("DesignMethod".into());
    vals.push(Value::chars(&method));
    let o = Object { class: class_for(&names), props: Rc::new(RefCell::new(vals)) };
    out(vec![Value::Object(o)])
}

const METHODS: &[&str] = &[
    "filter", "filtfilt", "fftfilt", "freqz", "impz", "impzlength", "stepz", "zplane", "grpdelay", "phasez", "phasedelay", "zerophase", "filtord", "firtype", "info",
    "isfir", "islinphase", "isminphase", "ismaxphase", "isallpass", "isstable", "isdouble", "issingle", "double", "single", "tf", "zpk", "ss", "eq", "ne",
];

fn is_fir(o: &Object) -> bool {
    o.class.prop_index("ImpulseResponse").map(|i| o.props.borrow()[i].text().unwrap_or_default() == "fir").unwrap_or(true)
}

/// The numerator and denominator (row vectors) of a digital filter.
fn tf_of(it: &mut Interpreter, o: &Object) -> R<(Value, Value)> {
    let c = o.props.borrow()[0].to_matrix()?;
    if is_fir(o) {
        Ok((Value::Num(c), scalar_value(1.0)))
    } else {
        let r = it.call_function("sos2tf", vec![Value::Num(c)], 2)?;
        Ok((r[0].clone(), r[1].clone()))
    }
}

fn prop_num(o: &Object, name: &str) -> Option<f64> {
    o.class.prop_index(name).and_then(|i| o.props.borrow()[i].scalar().ok())
}

/// The object functions (the methods the class forwards here).
fn method(it: &mut Interpreter, a: &[Value]) -> Ret {
    let name = a[0].text().unwrap_or_default();
    let nargout = a[1].scalar().unwrap_or(1.0) as usize;
    let o = match &a[2] {
        Value::Object(o) => o.clone(),
        _ => return rt("Invalid digitalFilter."),
    };
    let rest: Vec<Value> = a[3..].to_vec();
    let fir = is_fir(&o);
    let coef = o.props.borrow()[0].clone();
    let boolv = |b: bool| out(vec![Value::bool(b)]);
    match name.as_str() {
        "isfir" => boolv(fir),
        "isdouble" | "issingle" => {
            let single = coef.to_matrix()?.kind == Kind::Single;
            boolv(if name == "issingle" { single } else { !single })
        }
        "double" | "single" => {
            let mut m = coef.to_matrix()?;
            m.kind = if name == "single" { Kind::Single } else { Kind::Double };
            let mut vals = o.props.borrow().clone();
            vals[0] = Value::Num(m);
            out(vec![Value::Object(Object { class: o.class.clone(), props: Rc::new(RefCell::new(vals)) })])
        }
        "eq" | "ne" => {
            let same = match rest.first() {
                Some(Value::Object(p)) => p.class.name == o.class.name && *p.props.borrow() == *o.props.borrow(),
                _ => false,
            };
            boolv(if name == "eq" { same } else { !same })
        }
        "filter" => {
            if nargout > 1 {
                return rt("Too many output arguments.");
            }
            let x = rest.first().cloned().ok_or_else(|| crate::error::LangError::new("Not enough input arguments.", 0))?;
            let r = if fir {
                it.call_function("filter", vec![coef.clone(), scalar_value(1.0), x], 1)?
            } else {
                it.call_function("sosfilt", vec![coef.clone(), x], 1)?
            };
            out(r)
        }
        "fftfilt" | "filtfilt" => {
            if !fir && name == "filtfilt" {
                let x = rest.first().cloned().ok_or_else(|| crate::error::LangError::new("Not enough input arguments.", 0))?;
                let m = coef.to_matrix()?;
                // One section is a plain transfer function; several are filtered section by section.
                let args = if m.dims[0] >= 2 {
                    vec![coef.clone(), scalar_value(1.0), x]
                } else {
                    vec![real_row((0..3).map(|c| m.data[c]).collect()), real_row((3..6).map(|c| m.data[c]).collect()), x]
                };
                return out(it.call_function("filtfilt", args, 1)?);
            }
            if !fir {
                return rt("fftfilt: the filter must be an FIR digitalFilter.");
            }
            let x = rest.first().cloned().ok_or_else(|| crate::error::LangError::new("Not enough input arguments.", 0))?;
            let args = if name == "fftfilt" { vec![coef.clone(), x] } else { vec![coef.clone(), scalar_value(1.0), x] };
            out(it.call_function(&name, args, 1)?)
        }
        "firtype" => {
            if !fir {
                return rt("The firtype function is only valid for FIR filters.");
            }
            out(it.call_function("firtype", vec![coef.clone()], 1)?)
        }
        "filtord" => {
            if fir {
                out(vec![scalar_value(coef.to_matrix()?.numel() as f64 - 1.0)])
            } else {
                let (b, den) = tf_of(it, &o)?;
                out(it.call_function("filtord", vec![b, den], 1)?)
            }
        }
        "tf" => {
            let (b, den) = tf_of(it, &o)?;
            out(vec![b, den].into_iter().take(nargout.max(1)).collect())
        }
        "zpk" => {
            if fir {
                let (b, den) = tf_of(it, &o)?;
                out(it.call_function("tf2zpk", vec![b, den], nargout.max(1))?)
            } else {
                // Zeros and poles section by section, in the order of the sections.
                let m = coef.to_matrix()?;
                let l = m.dims[0];
                let (mut zr, mut zi, mut pr, mut pi_) = (vec![], vec![], vec![], vec![]);
                let mut k = 1.0;
                for i in 0..l {
                    k *= m.data[i];
                    for (part, re, im) in [(0usize, &mut zr, &mut zi), (3usize, &mut pr, &mut pi_)] {
                        let row = real_row((0..3).map(|c| m.data[(part + c) * l + i]).collect());
                        let r = it.call_function("roots", vec![row], 1)?.remove(0).to_matrix()?;
                        let imv = r.im.clone().unwrap_or_else(|| vec![0.0; r.numel()]);
                        let mut idx: Vec<usize> = (0..r.numel()).collect();
                        if r.numel() == 2 && imv.iter().all(|v| *v == 0.0) {
                            idx.sort_by(|&x, &y| r.data[x].partial_cmp(&r.data[y]).unwrap_or(std::cmp::Ordering::Equal));
                        }
                        for q in idx {
                            re.push(r.data[q]);
                            im.push(imv[q]);
                        }
                    }
                }
                let zm = Matrix::complex(vec![zr.len(), 1], zr, zi);
                let pm = Matrix::complex(vec![pr.len(), 1], pr, pi_);
                let outs = vec![Value::Num(zm), Value::Num(pm), scalar_value(k)];
                out(outs.into_iter().take(nargout.max(1)).collect())
            }
        }
        "ss" => {
            if fir {
                let m = coef.to_matrix()?;
                let mut den = vec![0.0; m.numel()];
                den[0] = 1.0;
                out(it.call_function("tf2ss", vec![coef.clone(), real_row(den)], nargout.max(1))?)
            } else {
                out(it.call_function("sos2ss", vec![coef.clone()], nargout.max(1))?)
            }
        }
        "isminphase" if !fir => {
            // Zeros and poles of every section on or inside the unit circle.
            let m = coef.to_matrix()?;
            let l = m.dims[0];
            let mut ok = true;
            for i in 0..l {
                for part in [0usize, 3] {
                    let row = real_row((0..3).map(|c| m.data[(part + c) * l + i]).collect());
                    let r = it.call_function("roots", vec![row], 1)?.remove(0).to_matrix()?;
                    let im = r.im.clone().unwrap_or_else(|| vec![0.0; r.numel()]);
                    if (0..r.numel()).any(|k| r.data[k].hypot(im[k]) > 1.0 + 1e-8) {
                        ok = false;
                    }
                }
            }
            out(vec![Value::bool(ok)])
        }
        "ismaxphase" if !fir => {
            let m = coef.to_matrix()?;
            let l = m.dims[0];
            let mut ok = true;
            for i in 0..l {
                let row = real_row((0..3).map(|c| m.data[c * l + i]).collect());
                let r = it.call_function("roots", vec![row], 1)?.remove(0).to_matrix()?;
                let im = r.im.clone().unwrap_or_else(|| vec![0.0; r.numel()]);
                if (0..r.numel()).any(|k| r.data[k].hypot(im[k]) < 1.0 - 1e-8) {
                    ok = false;
                }
            }
            out(vec![Value::bool(ok)])
        }
        "grpdelay" | "phasez" | "phasedelay" | "zerophase" if !fir => {
            // Sums over the sections: group delay and (unwrapped) phase add.
            let m = coef.to_matrix()?;
            let l = m.dims[0];
            let mut acc: Option<Vec<f64>> = None;
            let mut wv: Option<Value> = None;
            for i in 0..l {
                let bsec = real_row((0..3).map(|c| m.data[c * l + i]).collect());
                let asec = real_row((0..3).map(|c| m.data[(3 + c) * l + i]).collect());
                let call_name = if name == "grpdelay" { "grpdelay" } else { "phasez" };
                let mut args = vec![bsec, asec];
                args.extend(rest.iter().cloned());
                let r = it.call_function(call_name, args, 2)?;
                let v = r[0].to_matrix()?.data;
                acc = Some(match acc {
                    None => v,
                    Some(prev) => prev.iter().zip(v.iter()).map(|(x, y)| x + y).collect(),
                });
                wv = Some(r[1].clone());
            }
            let w = wv.unwrap();
            let wm = w.to_matrix()?;
            let ph = acc.unwrap();
            let n = ph.len();
            let res: Vec<f64> = match name.as_str() {
                "phasedelay" => (0..n).map(|k| -ph[k] / wm.data[k]).collect(),
                _ => ph.clone(),
            };
            if name == "zerophase" {
                let hr = it.call_function("freqz", vec![coef.clone(), Value::Num(Matrix::from_dims(vec![0, 0], vec![]))], 1);
                let _ = hr;
                // Hr = real(H exp(-j phi)) with H from the transfer function.
                let (bb, aa) = tf_of(it, &o)?;
                let mut args = vec![bb, aa];
                args.extend(rest.iter().cloned());
                let rr = it.call_function("freqz", args, 2)?;
                let h = rr[0].to_matrix()?;
                let him = h.im.clone().unwrap_or_else(|| vec![0.0; n]);
                let hrv: Vec<f64> = (0..n)
                    .map(|k| {
                        let p = if ph[k].is_nan() { 0.0 } else { ph[k] };
                        h.data[k] * p.cos() + him[k] * p.sin()
                    })
                    .collect();
                let outs = vec![Value::Num(Matrix::from_dims(vec![n, 1], hrv)), w, Value::Num(Matrix::from_dims(vec![n, 1], ph))];
                return out(outs.into_iter().take(nargout.max(1)).collect());
            }
            let outs = vec![Value::Num(Matrix::from_dims(vec![n, 1], res)), w];
            out(outs.into_iter().take(nargout.max(1)).collect())
        }
        "info" => rt("info: not implemented."),
        _ => {
            // Responses and tests of the transfer function: the same function on (b, a, ...).
            let (b, den) = tf_of(it, &o)?;
            let mut args = vec![b, den];
            args.extend(rest);
            out(it.call_function(&name, args, nargout)?)
        }
    }
}
