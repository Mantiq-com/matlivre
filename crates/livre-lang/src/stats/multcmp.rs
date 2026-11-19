//! `multcompare` (pairwise comparisons of group means after `anova1`, `anova2`,
//! `kruskalwallis` or `friedman`) and the studentized range distribution it needs.

use super::dist::{fcdf, finv, norminv_std, tcdf, tinv};
use super::tests3::param_val;
use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::{Interpreter, OutputKind};
use crate::value::*;

pub const NAMES: &[&str] = &["multcompare"];

fn normcdf(x: f64) -> f64 {
    0.5 * crate::builtins_specfun::erfc(-x / std::f64::consts::SQRT_2)
}

// ---------------------------------------------------------------------------------------------
// Studentized range distribution: the grid sum of Lund and Lund (Applied Statistics 190)

/// Distribution function (or the upper tail when `upper`) of the studentized range of `r` samples
/// with `v` degrees of freedom at `q`.
pub fn stdrcdf(q: f64, v: f64, r: f64, upper: bool) -> f64 {
    const JMAX: i32 = 15;
    const KMAX: i32 = 15;
    const STEP: f64 = 0.45;
    const VMAX: f64 = 120.0;
    if !(v >= 1.0 && r > 1.0 && r == r.round()) {
        return f64::NAN;
    }
    if !(q > 0.0) {
        return if upper { 1.0 } else { 0.0 };
    }
    let g = STEP / r.powf(0.2);
    let big_v = v > VMAX;
    let (c, h) = if big_v {
        ((r * g / (2.0 * std::f64::consts::PI).sqrt()).ln(), 0.0)
    } else {
        let h = STEP / v.sqrt();
        let v2 = v * 0.5;
        (((2.0 / std::f64::consts::PI).sqrt() * r * g * h).ln() - v2 + v2 * v2.ln() - crate::builtins_specfun::gammaln(v2), h)
    };
    let nk = (2 * KMAX + 1) as usize;
    let gk: Vec<f64> = (0..nk).map(|i| 0.5 * r.ln() + g * ((i as i32 - KMAX) as f64)).collect();
    let w0: Vec<f64> = gk.iter().map(|x| c - 0.5 * x * x).collect();
    let pz: Vec<f64> = gk.iter().map(|&x| normcdf(-x)).collect();
    let rm1 = (r - 1.0) as i32;
    let mut xx = 0.0;
    if big_v {
        for k in 0..nk {
            let rect = (normcdf(q - gk[k]) - pz[k]).powi(rm1);
            let term = if upper { (1.0 - pz[k]).powi(rm1) - rect } else { rect };
            xx += w0[k].exp() * term;
        }
    } else {
        let nj = (2 * JMAX + 1) as usize;
        let hj: Vec<f64> = (0..nj).map(|i| h * ((i as i32 - JMAX) as f64)).collect();
        let ehj: Vec<f64> = hj.iter().map(|x| x.exp()).collect();
        let qw: Vec<f64> = ehj.iter().map(|e| q * e).collect();
        let vw: Vec<f64> = hj.iter().zip(&ehj).map(|(x, e)| v * (x + 0.5 * (1.0 - e * e))).collect();
        for k in 0..nk {
            let xhq = (1.0 - pz[k]).powi(rm1);
            let mut col = 0.0;
            for j in 0..nj {
                let rect = (normcdf(qw[j] - gk[k]) - pz[k]).powi(rm1);
                let term = if upper { xhq - rect } else { rect };
                col += (w0[k] + vw[j]).exp() * term;
            }
            xx += col;
        }
    }
    xx
}

/// Initial quantile of AS 190.2.
fn qtrng0(p: f64, v: f64, r: f64) -> f64 {
    let mut t = norminv_std(0.5 + 0.5 * p);
    if v < 120.0 {
        t += 0.25 * (t * t * t + t) / v;
    }
    let mut q = 0.8843 - 0.2368 * t;
    if v < 120.0 {
        q = q - (1.214 / v) + (1.208 * t / v);
    }
    t * (q * (r - 1.0).ln() + 1.4142)
}

/// Inverse of the studentized range distribution (the secant iteration of AS 190.1).
pub fn stdrinv(p: f64, v: f64, r: f64) -> f64 {
    const JMAX: usize = 20;
    const PCUT: f64 = 0.00001;
    const TINY: f64 = 0.000001;
    if !(v >= 1.0 && r > 1.0 && r == r.round() && p < 1.0) {
        return f64::NAN;
    }
    if !(p > 0.0) {
        return 0.0;
    }
    let upper = p > 0.99;
    let p0 = if upper { 1.0 - p } else { p };
    let mut q1 = qtrng0(p, v, r);
    let mut p1 = stdrcdf(q1, v, r, upper);
    let mut xx = q1;
    if (p1 - p0).abs() >= PCUT * p0 {
        let mut p2 = f64::NAN;
        if p1 > p0 {
            p2 = (0.75 * p0).max(p0 - 0.75 * (p1 - p0));
        }
        if p1 < p0 {
            p2 = p0 + (p0 - p1) * (1.0 - p0) / (1.0 - p1) * 0.75;
        }
        let mut q2 = if upper { qtrng0(1.0 - p2, v, r) } else { qtrng0(p2, v, r) };
        for _ in 2..=JMAX {
            p2 = stdrcdf(q2, v, r, upper);
            let e1 = p1 - p0;
            let e2 = p2 - p0;
            let d = e2 - e1;
            xx = (q1 + q2) / 2.0;
            if d.abs() > TINY * p0 {
                xx = (e2 * q1 - e1 * q2) / d;
            }
            if e1.abs() >= e2.abs() {
                q1 = q2;
                p1 = p2;
            }
            if (p1 - p0).abs() < PCUT * p0 {
                break;
            }
            q2 = xx;
        }
    }
    xx
}

// ---------------------------------------------------------------------------------------------
// multcompare

const BAD_STATS: &str = "STATS must be a stats output structure from ANOVA1, ANOVA2, ANOVAN, AOCTOOL, KRUSKALWALLIS, or FRIEDMAN.";
const BAD_CTYPE: &str = "Unrecognized value for 'CriticalValueType' argument.";

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn field_vec(st: &StructArray, name: &str) -> Res<Vec<f64>> {
    match st.get(name) {
        Some(Value::Num(m)) => Ok(m.data.clone()),
        Some(_) => rt(BAD_STATS),
        None => rt(format!("Unrecognized field name \"{name}\".")),
    }
}

fn field_scalar(st: &StructArray, name: &str) -> Res<f64> {
    let v = field_vec(st, name)?;
    match v.first() {
        Some(&x) => Ok(x),
        None => rt(BAD_STATS),
    }
}

/// The critical value and p-values of one critical-value type.
fn crit_one(onetype: &str, alpha: f64, df: f64, ng: usize, t: &[f64], alternative: &str) -> Res<(f64, Vec<f64>)> {
    let kstar = (ng * (ng - 1) / 2) as f64;
    let mut onetype = onetype.to_string();
    if onetype.chars().count() == 1 {
        onetype = match onetype.as_str() {
            "t" | "h" => "tukey-kramer",
            "d" => "dunn-sidak",
            "b" => "bonferroni",
            "s" => "scheffe",
            "l" => "lsd",
            o => o,
        }
        .to_string();
    }
    if onetype == "hsd" {
        onetype = "tukey-kramer".into();
    }
    let ngf = ng as f64;
    let (crit1, mut pval1): (f64, Vec<f64>) = match onetype.as_str() {
        "tukey-kramer" => {
            let (ub, _) = crit_one("dunn-sidak", alpha, df, ng, &[], alternative)?;
            let crit1 = ub.min(stdrinv(1.0 - alpha, df, ngf) / std::f64::consts::SQRT_2);
            (crit1, t.iter().map(|x| stdrcdf(std::f64::consts::SQRT_2 * x.abs(), df, ngf, true)).collect())
        }
        "dunn-sidak" => {
            let alf = 1.0 - (1.0 - alpha).powf(1.0 / kstar);
            (tinv(1.0 - alf / 2.0, df), t.iter().map(|x| 1.0 - (1.0 - 2.0 * tcdf(-x.abs(), df)).powf(kstar)).collect())
        }
        "bonferroni" => (tinv(1.0 - alpha / (2.0 * kstar), df), t.iter().map(|x| 2.0 * kstar * tcdf(-x.abs(), df)).collect()),
        "lsd" => (tinv(1.0 - alpha / 2.0, df), t.iter().map(|x| 2.0 * tcdf(-x.abs(), df)).collect()),
        "scheffe" => {
            let tmp = finv(1.0 - alpha, ngf - 1.0, df);
            ((((ngf - 1.0) * tmp).sqrt()), t.iter().map(|x| fcdf((x * x) / (ngf - 1.0), ngf - 1.0, df, true)).collect())
        }
        "dunnett" => {
            let df = df.min(1e4);
            if t.len() != 1 {
                return rt("The 'dunnett' critical value type is only available here for two groups.");
            }
            match alternative {
                "notequal" => (-tinv(alpha / 2.0, df), vec![2.0 * tcdf(-t[0].abs(), df)]),
                "less" => (tinv(alpha, df), vec![tcdf(t[0], df)]),
                "greater" => (-tinv(alpha, df), vec![1.0 - tcdf(t[0], df)]),
                _ => return rt("Output argument \"crit\" (and possibly others) not assigned a value in the execution with \"multcompare>getcritDunnett\" function."),
            }
        }
        _ => return rt(BAD_CTYPE),
    };
    for p in pval1.iter_mut() {
        if *p > 1.0 {
            *p = 1.0;
        }
    }
    Ok((crit1, pval1))
}

fn getcrit(ctype: &str, alpha: f64, df: f64, ng: usize, t: &[f64], alternative: &str) -> Res<(f64, Vec<f64>)> {
    let ctype = if ctype.is_empty() { "tukey-kramer" } else { ctype };
    let types: Vec<&str> = ctype.split_whitespace().collect();
    if types.len() > 1 && types.contains(&"dunnett") {
        return rt(BAD_CTYPE);
    }
    let mut crit = f64::INFINITY;
    let mut pval = vec![f64::NAN; t.len()];
    if types.is_empty() {
        return rt(BAD_CTYPE);
    }
    for ty in types {
        let (crit1, pval1) = crit_one(ty, alpha, df, ng, t, alternative)?;
        if !crit.is_finite() || crit1 < crit {
            crit = crit.min(crit1);
            pval = pval1;
        }
    }
    Ok((crit, pval))
}

struct Model {
    gmeans: Vec<f64>,
    /// Variances of the group means (the covariance matrix is diagonal for these sources).
    gvar: Vec<f64>,
    /// The full covariance of the means, when it is not diagonal.
    gcov: Option<Vec<Vec<f64>>>,
    df: f64,
    ng_crit: usize,
    gnames: Value,
    single: bool,
    note: Option<String>,
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

fn number_names(ng: usize) -> Value {
    char_rows(&(1..=ng).map(|i| i.to_string()).collect::<Vec<_>>())
}

fn single_of(st: &StructArray, name: &str) -> bool {
    matches!(st.get(name), Some(Value::Num(m)) if is_single(m))
}

fn model(st: &StructArray, source: &str, estimate: &Option<String>, dim_v: &Option<Value>, dodisp: bool, alpha: f64) -> Res<Model> {
    match source {
        "anovan" => anovan_model(st, dim_v),
        "anova1" => {
            let gmeans = field_vec(st, "means")?;
            let n = field_vec(st, "n")?;
            let df = field_scalar(st, "df")?;
            let s = field_scalar(st, "s")?;
            let ng = n.iter().filter(|&&c| c > 0.0).count();
            if df < 1.0 {
                return rt("There is not enough data to compute multiple comparisons: at least one group needs more than one observation.");
            }
            let gvar = n.iter().map(|&c| s * s / c).collect();
            let gnames = st.get("gnames").cloned().unwrap_or_else(|| number_names(gmeans.len()));
            Ok(Model { gmeans, gvar, gcov: None, df, ng_crit: ng, gnames, single: single_of(st, "means"), note: None })
        }
        "anova2" => {
            let docols = estimate.as_deref() != Some("row");
            let (gmeans, n) = if docols { (field_vec(st, "colmeans")?, field_scalar(st, "coln")?) } else { (field_vec(st, "rowmeans")?, field_scalar(st, "rown")?) };
            let sigmasq = field_scalar(st, "sigmasq")?;
            let df = field_scalar(st, "df")?;
            if df < 1.0 {
                return rt("There is not enough data to compute multiple comparisons: at least one group needs more than one observation.");
            }
            let ng = gmeans.len();
            let sigma = sigmasq.sqrt();
            let inter = matches!(st.get("inter"), Some(v) if v.scalar().map(|x| x != 0.0).unwrap_or(false));
            let note = if inter && dodisp {
                let pv = field_scalar(st, "pval")?;
                Some(if pv < alpha {
                    "Note: Your model includes an interaction term that is significant at the level you specified. Testing main effects under these conditions is questionable.".to_string()
                } else {
                    "Note: Your model includes an interaction term.  A test of main effects can be difficult to interpret when the model includes interactions.".to_string()
                })
            } else {
                None
            };
            Ok(Model { gmeans: gmeans.clone(), gcov: None, gvar: vec![sigma * sigma / n; ng], df, ng_crit: ng, gnames: number_names(ng), single: single_of(st, if docols { "colmeans" } else { "rowmeans" }), note })
        }
        "kruskalwallis" => {
            let gmeans = field_vec(st, "meanranks")?;
            let n = field_vec(st, "n")?;
            let sumt = field_scalar(st, "sumt")?;
            let ng = n.len();
            let big_n: f64 = n.iter().sum();
            let gvar = n.iter().map(|&c| ((big_n * (big_n + 1.0) / 12.0) - (sumt / (12.0 * (big_n - 1.0)))) / c).collect();
            let gnames = st.get("gnames").cloned().unwrap_or_else(|| number_names(gmeans.len()));
            Ok(Model { gmeans, gvar, gcov: None, df: f64::INFINITY, ng_crit: ng, gnames, single: single_of(st, "meanranks"), note: Some("Note: Intervals can be used for testing but are not simultaneous confidence intervals.".into()) })
        }
        "friedman" => {
            let gmeans = field_vec(st, "meanranks")?;
            let n = field_scalar(st, "n")?;
            let sigma = field_scalar(st, "sigma")?;
            let ng = gmeans.len();
            Ok(Model { gmeans: gmeans.clone(), gcov: None, gvar: vec![sigma * sigma / n; ng], df: f64::INFINITY, ng_crit: ng, gnames: number_names(ng), single: single_of(st, "meanranks"), note: Some("Note: Intervals can be used for testing but are not simultaneous confidence intervals.".into()) })
        }
        _ => rt(BAD_STATS),
    }
}

/// Canonical parameter names: (alias list, canonical).
const PARAMS: [(&[&str], &str); 8] = [
    (&["alpha"], "alpha"),
    (&["displayopt"], "displayopt"),
    (&["ctype", "criticalvaluetype", "c"], "ctype"),
    (&["estimate"], "estimate"),
    (&["dimension"], "dimension"),
    (&["controlgroup"], "controlgroup"),
    (&["approximate"], "approximate"),
    (&["alternative"], "alternative"),
];

fn match_param(key: &str) -> Res<&'static str> {
    let k = key.to_lowercase();
    for (aliases, canon) in PARAMS.iter() {
        if aliases.contains(&k.as_str()) {
            return Ok(canon);
        }
    }
    let hits: Vec<&'static str> = PARAMS.iter().filter(|(al, _)| !k.is_empty() && al.iter().any(|a| a.starts_with(&k))).map(|(_, c)| *c).collect();
    match hits.len() {
        1 => Ok(hits[0]),
        0 => rt(format!("Invalid parameter name: {key}.")),
        _ => rt(format!("Ambiguous parameter name: {key}.")),
    }
}

pub fn call(it: &mut Interpreter, _name: &str, a: &[Value], nargout: usize) -> Ret {
    if a.is_empty() {
        return rt("Not enough input arguments.");
    }
    let rest = &a[1..];
    let mut alpha_v = Value::num(0.05);
    let mut disp_v = Value::chars("on");
    let mut ctype_v = Value::chars("");
    let mut est_v = Value::chars("");
    let mut dim_v: Option<Value> = None;
    let mut control_set = false;
    let mut control_v = Value::num(1.0);
    let mut alt_v = Value::chars("notequal");
    let mut approx_v: Option<Value> = None;
    let named = rest.is_empty() || (!is_empty_val(&rest[0]) && rest[0].is_text());
    if named {
        if rest.len() % 2 != 0 {
            return rt("Wrong number of arguments.");
        }
        for pair in rest.chunks(2) {
            let key = match pair[0].text() {
                Ok(k) => k,
                Err(_) => return rt("Parameter name must be text."),
            };
            let v = pair[1].clone();
            match match_param(&key)? {
                "alpha" => alpha_v = v,
                "displayopt" => disp_v = v,
                "ctype" => ctype_v = v,
                "estimate" => est_v = v,
                "dimension" => dim_v = Some(v),
                "controlgroup" => {
                    control_set = true;
                    control_v = v;
                }
                "approximate" => approx_v = Some(v),
                _ => alt_v = v,
            }
        }
    } else {
        let n = rest.len();
        if n > 0 {
            alpha_v = rest[0].clone();
        }
        if n > 1 {
            disp_v = rest[1].clone();
        }
        if n > 2 {
            ctype_v = rest[2].clone();
        }
        if n > 3 {
            est_v = rest[3].clone();
        }
    }
    let Value::Struct(st) = &a[0] else { return rt(BAD_STATS) };
    let source = match st.get("source") {
        Some(v) if v.is_text() => v.text().unwrap_or_default(),
        Some(_) => String::new(),
        None => return rt(BAD_STATS),
    };
    let alpha = match &alpha_v {
        Value::Num(m) if m.kind != Kind::Char && m.kind != Kind::Logical && m.data.len() == 1 && m.data[0].is_finite() && m.data[0] > 0.0 && m.data[0] < 1.0 && (m.kind == Kind::Double || is_single(m)) => m.data[0],
        _ => return rt("ALPHA must be a scalar between 0 and 1, exclusive."),
    };
    let dodisp = match disp_v.text() {
        Ok(t) if disp_v.is_text() && t == "on" => true,
        Ok(t) if disp_v.is_text() && t == "off" => false,
        _ => return rt("DISPLAYOPT must be 'on' or 'off'."),
    };
    let ctype_text = if ctype_v.is_text() { Some(ctype_v.text().unwrap_or_default()) } else { None };
    let is_dunnett = ctype_text.as_deref().map(|t| t.eq_ignore_ascii_case("dunnett")).unwrap_or(false);
    if !is_dunnett && control_set {
        return rt("'Control' argument is valid only when 'CriticalValueType' is set to 'dunnett'.");
    }
    let mut alternative = "notequal".to_string();
    if !is_empty_val(&alt_v) {
        let i = param_val(&alt_v, &["notequal", "less", "greater"], "'Alternative'", true)?;
        let _ = i;
        alternative = alt_v.text().unwrap_or_default();
        if !alternative.eq_ignore_ascii_case("notequal") && !is_dunnett {
            return rt("'Alternative' argument is valid only when 'CriticalValueType' is set to 'dunnett'.");
        }
    }
    if let Some(v) = &approx_v {
        if !is_empty_val(v) {
            let ok = matches!(v, Value::Num(m) if m.kind != Kind::Char);
            if !ok {
                return rt("Expected Approximate to be one of these types:\n\ndouble, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64");
            }
            if let Value::Num(m) = v {
                if m.kind != Kind::Logical {
                    if m.data.len() != 1 {
                        return rt("Expected Approximate to be a scalar.");
                    }
                    let x = m.data[0];
                    if x != x.round() || x.is_infinite() {
                        return rt("Expected Approximate to be integer-valued.");
                    }
                    if x < 0.0 {
                        return rt("Expected Approximate to be a scalar with value >= 0.");
                    }
                    if x > 1.0 {
                        return rt("Expected Approximate to be a scalar with value <= 1.");
                    }
                }
            }
        }
    }
    let estimate = if is_empty_val(&est_v) || (est_v.is_text() && est_v.text().unwrap_or_default().is_empty()) {
        None
    } else if source == "anova2" {
        if !est_v.is_text() {
            return rt("Invalid ESTIMATE data type. Argument must be a character vector or string scalar.");
        }
        let i =param_val(&est_v, &["row", "column"], "ESTIMATE", true)?;
        Some(["row", "column"][i].to_string())
    } else {
        None
    };
    let md = model(st, &source, &estimate, &dim_v, dodisp, alpha)?;
    if !ctype_v.is_text() && !is_empty_val(&ctype_v) {
        return rt(BAD_CTYPE);
    }
    let ctype = ctype_text.unwrap_or_default();
    let ng = md.gmeans.len();
    // gvar of groups with a NaN mean does not enter the comparisons
    let nan_g: Vec<bool> = md.gmeans.iter().map(|v| v.is_nan()).collect();
    let cov0: Vec<Vec<f64>> = match &md.gcov {
        Some(c) => c.clone(),
        None => (0..ng).map(|i| (0..ng).map(|j| if i == j { md.gvar[i] } else { 0.0 }).collect()).collect(),
    };
    let mut cov = cov0.clone();
    for i in 0..ng {
        if nan_g[i] {
            for j in 0..ng {
                cov[i][j] = 0.0;
                cov[j][i] = 0.0;
            }
        }
    }
    let (g1, g2): (Vec<usize>, Vec<usize>) = if is_dunnett {
        let control = match &control_v {
            Value::Num(m) if m.kind != Kind::Char && m.data.len() == 1 && m.data[0] >= 1.0 && m.data[0] <= ng as f64 && m.data[0] == m.data[0].round() => m.data[0] as usize,
            _ => return rt(format!("'Control' value must be an integer between 1 and {ng}.")),
        };
        let g1: Vec<usize> = (1..=ng).filter(|&g| g != control).collect();
        let g2 = vec![control; g1.len()];
        (g1, g2)
    } else {
        let mut g1 = Vec::new();
        let mut g2 = Vec::new();
        for i in 1..=ng {
            for j in i + 1..=ng {
                g1.push(i);
                g2.push(j);
            }
        }
        (g1, g2)
    };
    let mn: Vec<f64> = g1.iter().zip(&g2).map(|(&i, &j)| md.gmeans[i - 1] - md.gmeans[j - 1]).collect();
    let se: Vec<f64> = g1.iter().zip(&g2).map(|(&i, &j)| (cov[i - 1][i - 1] + cov[j - 1][j - 1] - 2.0 * cov[i - 1][j - 1]).sqrt()).collect();
    let t: Vec<f64> = mn.iter().zip(&se).map(|(m, s)| m / s).collect();
    let (crit, pval) = getcrit(&ctype, alpha, md.df, md.ng_crit, &t, &alternative)?;
    let np = mn.len();
    let mut m = vec![0.0; np * 6];
    for i in 0..np {
        let delta = crit * se[i];
        m[i] = g1[i] as f64;
        m[np + i] = g2[i] as f64;
        m[2 * np + i] = mn[i] - delta;
        m[3 * np + i] = mn[i];
        m[4 * np + i] = mn[i] + delta;
        m[5 * np + i] = pval[i];
    }
    match alternative.as_str() {
        "less" => (0..np).for_each(|i| m[2 * np + i] = f64::NEG_INFINITY),
        "greater" => (0..np).for_each(|i| m[4 * np + i] = f64::INFINITY),
        _ => {}
    }
    if dodisp {
        if let Some(note) = &md.note {
            if source == "anova2" {
                it.emit(OutputKind::Text(format!("{note}\n")));
            } else {
                it.emit(OutputKind::Text(format!("{note}\n")));
            }
        }
    }
    let mut mm = vec![0.0; ng * 2];
    for i in 0..ng {
        mm[i] = md.gmeans[i];
        mm[ng + i] = if nan_g[i] { f64::NAN } else { cov0[i][i].sqrt() };
    }
    let mut out = vec![Value::Num(finish(Matrix::from_dims(vec![np, 6], m), md.single))];
    if nargout > 1 {
        out.push(Value::Num(finish(Matrix::from_dims(vec![ng, 2], mm), md.single)));
    }
    if nargout > 2 {
        out.push(Value::empty());
    }
    if nargout > 3 {
        out.push(md.gnames);
    }
    Ok(Some(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1.0)
    }

    #[test]
    fn studentized_range_matches_matlab_values() {
        // Values printed by R2023b for the same grid-sum algorithm.
        assert!(close(stdrcdf(3.5, 10.0, 3.0, false), 0.92289668896292587));
        assert!(close(stdrcdf(3.5, 10.0, 3.0, true), 0.077103304075219037));
        assert!(close(stdrcdf(2.0, 2.0, 2.0, false), 0.70710786078337895));
        assert!(close(stdrcdf(0.5, 1000.0, 2.0, false), 0.27632639016823701));
        assert!(close(stdrcdf(10.0, f64::INFINITY, 20.0, true), 1.8358978886418918e-10));
        assert!(close(stdrinv(0.95, 10.0, 4.0), 4.3265821216854299));
        assert!(close(stdrinv(0.5, 2.0, 2.0), 1.1546939740310438));
        assert!(close(stdrinv(0.999, 100.0, 8.0), 6.0837567032283468));
        assert!(stdrinv(1.0, 10.0, 3.0).is_nan());
        assert_eq!(stdrinv(0.0, 10.0, 3.0), 0.0);
        assert!(stdrcdf(1.0, 0.5, 3.0, false).is_nan());
    }
}

/// The group means and their covariance after `anovan`, for the combinations of the levels of the
/// factors in `Dimension`.
fn anovan_model(st: &StructArray, dim_v: &Option<Value>) -> Res<Model> {
    let nlevels = field_vec(st, "nlevels")?;
    let nvars = nlevels.len();
    let terms_m = match st.get("terms") {
        Some(Value::Num(m)) => m.clone(),
        _ => return rt(BAD_STATS),
    };
    let nt = terms_m.dims[0];
    let termcols = field_vec(st, "termcols")?;
    let coeffs = field_vec(st, "coeffs")?;
    let continuous: Vec<bool> = match st.get("continuous") {
        Some(Value::Num(m)) => m.data.iter().map(|&d| d != 0.0).collect(),
        _ => vec![false; nvars],
    };
    let vmeans = field_vec(st, "vmeans")?;
    let dfe = field_scalar(st, "dfe")?;
    let mse0 = field_scalar(st, "mse")?;
    let single = single_of(st, "coeffs");
    if let Some(Value::Num(m)) = st.get("vnested") {
        if m.data.iter().any(|&d| d != 0.0) {
            return rt("MULTCOMPARE does not support nested models.");
        }
    }
    let missing = "The STATS structure is missing the information needed for the comparisons.";
    let nullp = match st.get("nullproject") {
        Some(Value::Num(m)) => m.clone(),
        _ => return rt(missing),
    };
    let rtr = match st.get("Rtr") {
        Some(Value::Num(m)) => m.clone(),
        _ => return rt(missing),
    };
    let rowb = match st.get("rowbasis") {
        Some(Value::Num(m)) => m.clone(),
        _ => return rt(missing),
    };
    let bad_dim = format!("Values of DIM must be integers between 1 and {nvars}.");
    let mut dim: Vec<usize> = match dim_v {
        Some(Value::Num(m)) if !m.is_empty() => {
            if m.data.iter().any(|&d| d < 1.0 || d > nvars as f64 || d != d.round()) {
                return rt(bad_dim);
            }
            m.data.iter().map(|&d| d as usize).collect()
        }
        _ => match nlevels.iter().position(|&l| l > 1.0) {
            Some(p) => vec![p + 1],
            None => return rt(bad_dim),
        },
    };
    dim.sort();
    dim.dedup();
    if dim.iter().any(|&d| nlevels[d - 1] < 2.0) {
        return rt("DIM must specify only categorical factors with 2 or more degrees of freedom.");
    }
    // every combination of the levels (the first factor varies fastest)
    let dl: Vec<usize> = dim.iter().map(|&d| nlevels[d - 1] as usize).collect();
    let nrows: usize = dl.iter().product();
    let mut ffd = vec![vec![0usize; dim.len()]; nrows];
    for r in 0..nrows {
        let mut rem = r;
        for (c, &l) in dl.iter().enumerate() {
            ffd[r][c] = rem % l;
            rem /= l;
        }
    }
    // dummy variables of every factor: the chosen ones, or the average over the levels
    let dums: Vec<Vec<Vec<f64>>> = (0..nvars)
        .map(|j| {
            let l = nlevels[j] as usize;
            if let Some(c) = dim.iter().position(|&d| d == j + 1) {
                (0..l).map(|k| (0..nrows).map(|r| if ffd[r][c] == k { 1.0 } else { 0.0 }).collect()).collect()
            } else if continuous[j] {
                vec![vec![vmeans[j]; nrows]]
            } else {
                (0..l).map(|_| vec![1.0 / l as f64; nrows]).collect()
            }
        })
        .collect();
    let ncols: usize = termcols.iter().sum::<f64>() as usize;
    let mut x = vec![vec![0.0; ncols]; nrows];
    let mut t0 = 0usize;
    for j in 0..=nt {
        let width = termcols[j] as usize;
        if j == 0 {
            for r in 0..nrows {
                x[r][0] = 1.0;
            }
        } else {
            let mut x0: Vec<Vec<f64>> = vec![];
            for k in (0..nvars).rev() {
                if terms_m.data[k * nt + (j - 1)] != 0.0 {
                    let d = &dums[k];
                    if x0.is_empty() {
                        x0 = d.clone();
                    } else {
                        let mut next = Vec::with_capacity(x0.len() * d.len());
                        for cb in d {
                            for ca in &x0 {
                                next.push(ca.iter().zip(cb).map(|(a, b)| a * b).collect());
                            }
                        }
                        x0 = next;
                    }
                }
            }
            for (c, col) in x0.iter().enumerate().take(width) {
                for r in 0..nrows {
                    x[r][t0 + c] = col[r];
                }
            }
        }
        t0 += width;
    }
    let mut gmeans: Vec<f64> = x.iter().map(|row| row.iter().zip(&coeffs).map(|(a, b)| a * b).sum()).collect();
    // covariance: mse * (L^-1 (x P0)')' (L^-1 (x P0)')
    let kk = nullp.dims[1];
    let pp = nullp.dims[0];
    let mse = if dfe == 0.0 { f64::NAN } else { mse0.max(0.0) };
    let mut tmp: Vec<Vec<f64>> = Vec::with_capacity(nrows);
    for r in 0..nrows {
        let v: Vec<f64> = (0..kk).map(|c| (0..pp).map(|i| x[r][i] * nullp.data[c * pp + i]).sum()).collect();
        // forward substitution with the lower-triangular Rtr
        let mut w = vec![0.0; kk];
        for i in 0..kk {
            let mut sacc = v[i];
            for j in 0..i {
                sacc -= rtr.data[j * kk + i] * w[j];
            }
            w[i] = sacc / rtr.data[i * kk + i];
        }
        tmp.push(w);
    }
    let gcov: Vec<Vec<f64>> = (0..nrows).map(|i| (0..nrows).map(|j| mse * dot_v(&tmp[i], &tmp[j])).collect()).collect();
    // means outside the row space of the design are not estimable
    let nb = rowb.dims[0];
    for r in 0..nrows {
        let mut resid: Vec<f64> = x[r].clone();
        for b in 0..nb {
            let c: f64 = (0..pp).map(|i| rowb.data[i * nb + b] * x[r][i]).sum();
            for i in 0..pp {
                resid[i] -= c * rowb.data[i * nb + b];
            }
        }
        let sum_abs: f64 = resid.iter().map(|v| v.abs()).sum();
        if sum_abs > f64::EPSILON.sqrt() * nrows as f64 {
            gmeans[r] = f64::NAN;
        }
    }
    // names of the groups
    let varnames: Vec<String> = match st.get("varnames") {
        Some(Value::Cell(c)) => c.data.iter().map(|v| v.text().unwrap_or_default()).collect(),
        _ => (1..=nvars).map(|i| format!("X{i}")).collect(),
    };
    let grpnames: Vec<Vec<String>> = match st.get("grpnames") {
        Some(Value::Cell(c)) => c.data.iter().map(|v| if let Value::Cell(cc) = v { cc.data.iter().map(|e| e.text().unwrap_or_default()).collect() } else { vec![] }).collect(),
        _ => vec![vec![]; nvars],
    };
    let gnames: Vec<Value> = (0..nrows)
        .map(|r| {
            let parts: Vec<String> = dim.iter().enumerate().map(|(c, &d)| format!("{}={}", varnames[d - 1], grpnames[d - 1][ffd[r][c]])).collect();
            Value::chars(&parts.join(","))
        })
        .collect();
    let gvar: Vec<f64> = (0..nrows).map(|i| gcov[i][i]).collect();
    Ok(Model { gmeans, gvar, gcov: Some(gcov), df: dfe, ng_crit: nrows, gnames: Value::Cell(CellArray::col(gnames)), single, note: None })
}

fn dot_v(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
