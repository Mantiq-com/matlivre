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

pub const NAMES: &[&str] = &["designfilt"];

const KNOWN: &[&str] = &[
    "FilterOrder", "CutoffFrequency", "CutoffFrequency1", "CutoffFrequency2", "PassbandFrequency", "PassbandFrequency1", "PassbandFrequency2", "StopbandFrequency",
    "StopbandFrequency1", "StopbandFrequency2", "PassbandRipple", "PassbandRipple1", "PassbandRipple2", "StopbandAttenuation", "StopbandAttenuation1",
    "StopbandAttenuation2", "PassbandWeight", "StopbandWeight", "StopbandWeight1", "StopbandWeight2", "SampleRate", "DesignMethod", "Window", "ScalePassband",
    "HalfPowerFrequency",
];

/// Properties that are used by the design but not listed among the specifications.
const HIDDEN: &[&str] = &["PassbandWeight", "StopbandWeight", "StopbandWeight1", "StopbandWeight2", "Window", "ScalePassband"];

const RESPONSES: &[&str] = &["lowpassfir", "highpassfir", "bandpassfir", "bandstopfir"];

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

pub fn call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
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
                _ => m == "equiripple" || (m == "ls" && !minorder),
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
                for &g in grid.iter().filter(|g| **g >= lo && **g <= hi) {
                    let (mut re, mut im) = (0.0, 0.0);
                    for (k, c) in h.iter().enumerate() {
                        let ph = -std::f64::consts::PI * g * k as f64;
                        re += c * ph.cos();
                        im += c * ph.sin();
                    }
                    let m = re.hypot(im);
                    let e = if av[bi] > 0.5 { (m - 1.0).abs() } else { m };
                    if e > *d {
                        return false;
                    }
                }
            }
            true
        };
        let mut found: Option<Matrix> = None;
        let mut n = (n0 - 8).max(1);
        while n <= n0 + 40 {
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
