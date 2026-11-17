//! Warning states: `warning('on'|'off'|'error'|'query', id)`, restoring a saved state with
//! `warning(s)`, and the identifiers of built-in warnings.

use crate::error::LangError;
use crate::interp::{Interpreter, OutputKind};
use crate::value::*;

/// The state of every warning: a default (`all`) and the identifiers set apart from it.
#[derive(Clone)]
pub(crate) struct WarningStates {
    all: String,
    over: Vec<(String, String)>,
    /// `backtrace` and `verbose` are modes: `warning('off', 'backtrace')`.
    modes: Vec<(String, String)>,
}

impl Default for WarningStates {
    fn default() -> Self {
        Self { all: "on".into(), over: Vec::new(), modes: vec![("backtrace".into(), "on".into()), ("verbose".into(), "off".into())] }
    }
}

impl WarningStates {
    pub(crate) fn state_of(&self, id: &str) -> &str {
        if let Some((_, s)) = self.modes.iter().find(|(i, _)| i == id) {
            return s;
        }
        self.over.iter().find(|(i, _)| i == id).map_or(self.all.as_str(), |(_, s)| s.as_str())
    }

    fn set(&mut self, id: &str, state: &str) {
        if let Some(e) = self.modes.iter_mut().find(|(i, _)| i == id) {
            e.1 = state.to_string();
        } else if id == "all" {
            self.all = state.to_string();
            self.over.clear();
        } else if state == self.all {
            // Back at the default: nothing is set apart.
            self.over.retain(|(i, _)| i != id);
        } else if let Some(e) = self.over.iter_mut().find(|(i, _)| i == id) {
            e.1 = state.to_string();
        } else {
            self.over.push((id.to_string(), state.to_string()));
        }
    }

    /// `all` first, then the identifiers set apart from it, in the order they were set.
    fn all_states(&self) -> Vec<(String, String)> {
        let mut v = vec![("all".to_string(), self.all.clone())];
        v.extend(self.over.iter().cloned());
        v
    }
}

/// The identifier MATLAB gives a built-in warning, from its text.
pub(crate) fn builtin_id(msg: &str) -> &'static str {
    if msg.starts_with("Matrix is singular to working precision") {
        "MATLAB:singularMatrix"
    } else if msg.starts_with("Matrix is close to singular or badly scaled") {
        "MATLAB:nearlySingularMatrix"
    } else {
        ""
    }
}

fn states_struct(v: &[(String, String)]) -> Value {
    let mut s = StructArray { dims: vec![v.len(), 1], fields: vec!["identifier".into(), "state".into()], elems: Vec::new() };
    for (id, st) in v {
        s.elems.push(vec![Value::chars(id), Value::chars(st)]);
    }
    Value::Struct(s)
}

impl Interpreter {
    /// Raises a warning: it is remembered for `lastwarn` and shown unless it is turned off.
    pub(crate) fn emit_warning(&mut self, msg: String, id: String) {
        self.last_warning = (msg.clone(), id.clone());
        if self.warning_states.state_of(&id) == "off" {
            return;
        }
        self.push_output(OutputKind::Warning(msg));
    }

    /// A warning raised by `warning(...)`: an error when its identifier is set to `error`.
    pub(crate) fn raise_warning(&mut self, msg: String, id: String) -> Result<(), LangError> {
        if self.warning_states.state_of(&id) == "error" {
            self.last_warning = (msg.clone(), id.clone());
            return Err(LangError::new(msg, 0).with_id(id).raised());
        }
        self.emit_warning(msg, id);
        Ok(())
    }
}

/// `warning` called with a state ('on', 'off', 'error', 'query'), or a saved state structure:
/// the values it returns. None for a call that raises a warning instead.
pub(crate) fn state_call(it: &mut Interpreter, a: &[Value], nargout: usize) -> Result<Option<Vec<Value>>, LangError> {
    let out = |it: &Interpreter, v: Vec<(String, String)>| -> Option<Vec<Value>> {
        let _ = it;
        Some(if nargout > 0 { vec![states_struct(&v)] } else { Vec::new() })
    };
    let Some(first) = a.first() else {
        return Ok(out(it, it.warning_states.all_states()));
    };
    if let Value::Struct(s) = first {
        let (Some(ik), Some(sk)) = (s.field_index("identifier"), s.field_index("state")) else {
            return Err(LangError::new("Invalid warning state structure.", 0));
        };
        let mut changes = Vec::new();
        for e in &s.elems {
            let (id, st) = (e[ik].text()?, e[sk].text()?);
            if !matches!(st.as_str(), "on" | "off" | "error") {
                return Err(LangError::new(format!("Invalid warning state '{st}'."), 0));
            }
            changes.push((id, st));
        }
        let before = it.warning_states.all_states();
        for (id, st) in changes {
            it.warning_states.set(&id, &st);
        }
        return Ok(out(it, before));
    }
    let Ok(cmd) = first.text() else { return Ok(None) };
    if !matches!(cmd.as_str(), "on" | "off" | "error" | "query") {
        return Ok(None);
    }
    if a.len() > 2 {
        return Err(LangError::new("Too many input arguments.", 0));
    }
    let id = match a.get(1) {
        Some(v) => v.text()?,
        None => "all".to_string(),
    };
    if cmd == "query" {
        if id == "all" {
            return Ok(out(it, it.warning_states.all_states()));
        }
        let st = it.warning_states.state_of(&id).to_string();
        if nargout == 0 {
            it.emit(OutputKind::Text(format!("The state of warning '{id}' is '{st}'.
")));
        }
        return Ok(out(it, vec![(id, st)]));
    }
    // The state before the call, so that `warning(s)` puts it back.
    let before = if id == "all" { it.warning_states.all_states() } else { vec![(id.clone(), it.warning_states.state_of(&id).to_string())] };
    it.warning_states.set(&id, &cmd);
    Ok(out(it, before))
}
