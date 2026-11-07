//! The bridge to a page that embeds the editor in a frame (`index.html?host=parent`): the
//! page sends the document and requests, the editor reports its state and hands back the
//! bytes to store. Messages travel with `postMessage` between the frame and its parent,
//! both on the same origin, under the prefix `matlivre:`.
//!
//! From the frame: `matlivre:ready { version }`, `matlivre:failed { error }`,
//! `matlivre:dirty { dirty }`, `matlivre:save { name, bytes, saveAs }` and
//! `matlivre:response { id, reply }`, where `reply` is `{ ok: true, result }` or
//! `{ ok: false, error }`. From the page: `matlivre:open { name, bytes }` and
//! `matlivre:request { id, method, params }`.

use std::cell::RefCell;
use std::collections::VecDeque;
use std::rc::Rc;

use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

const PREFIX: &str = "matlivre:";

/// What the page asked for, read by the app once a frame.
pub enum HostRequest {
    Open { name: String, bytes: Vec<u8> },
    Request { id: JsValue, method: String, params: serde_json::Value },
}

pub struct HostLink {
    parent: web_sys::Window,
    origin: String,
    inbox: Rc<RefCell<VecDeque<HostRequest>>>,
    _listener: Closure<dyn FnMut(web_sys::MessageEvent)>,
}

impl HostLink {
    /// Listens to the parent page when the editor was opened by one, which `?host=parent`
    /// says; a page opened on its own has no host.
    pub fn connect(ctx: eframe::egui::Context) -> Option<Self> {
        let window = web_sys::window()?;
        let search = window.location().search().ok()?;
        let hosted = web_sys::UrlSearchParams::new_with_str(&search).ok()?.get("host").as_deref() == Some("parent");
        let parent = window.parent().ok()??;
        if !hosted || parent == window {
            return None;
        }
        let origin = window.location().origin().ok()?;
        let inbox: Rc<RefCell<VecDeque<HostRequest>>> = Rc::new(RefCell::new(VecDeque::new()));
        let listener = {
            let inbox = inbox.clone();
            let origin = origin.clone();
            let parent = parent.clone();
            Closure::<dyn FnMut(web_sys::MessageEvent)>::new(move |event: web_sys::MessageEvent| {
                // Only the page that opened the frame, on this origin, is listened to. The
                // parent's window belongs to another realm, so it is compared by identity
                // rather than cast: a cross-realm window fails `instanceof Window`.
                if event.origin() != origin {
                    return;
                }
                let from_parent = event.source().is_some_and(|source| {
                    let source: &JsValue = source.as_ref();
                    let parent: &JsValue = parent.as_ref();
                    source == parent
                });
                if !from_parent {
                    return;
                }
                if let Some(request) = read_request(&event.data()) {
                    inbox.borrow_mut().push_back(request);
                    ctx.request_repaint();
                }
            })
        };
        window.add_event_listener_with_callback("message", listener.as_ref().unchecked_ref()).ok()?;
        Some(Self { parent, origin, inbox, _listener: listener })
    }

    /// The requests received since the last call, oldest first.
    pub fn take(&self) -> Vec<HostRequest> {
        self.inbox.borrow_mut().drain(..).collect()
    }

    pub fn ready(&self, version: &str) {
        let message = object("ready");
        set(&message, "version", &version.into());
        self.post(&message, None);
    }

    pub fn failed(&self, error: &str) {
        let message = object("failed");
        set(&message, "error", &error.into());
        self.post(&message, None);
    }

    pub fn dirty(&self, dirty: bool) {
        let message = object("dirty");
        set(&message, "dirty", &dirty.into());
        self.post(&message, None);
    }

    /// Hands the page the bytes of the document under `name`, which it stores; `save_as`
    /// tells it the document took this name just now.
    pub fn save(&self, name: &str, bytes: &[u8], save_as: bool) {
        let message = object("save");
        set(&message, "name", &name.into());
        let buffer = js_sys::Uint8Array::from(bytes).buffer();
        set(&message, "bytes", &buffer);
        set(&message, "saveAs", &save_as.into());
        self.post(&message, Some(&buffer));
    }

    pub fn respond(&self, id: JsValue, reply: Result<serde_json::Value, String>) {
        let message = object("response");
        set(&message, "id", &id);
        let answer = js_sys::Object::new();
        match reply {
            Ok(result) => {
                set(&answer, "ok", &true.into());
                let json = js_sys::JSON::parse(&result.to_string()).unwrap_or(JsValue::NULL);
                set(&answer, "result", &json);
            }
            Err(error) => {
                set(&answer, "ok", &false.into());
                set(&answer, "error", &error.into());
            }
        }
        set(&message, "reply", &answer);
        self.post(&message, None);
    }

    fn post(&self, message: &js_sys::Object, transfer: Option<&js_sys::ArrayBuffer>) {
        let result = match transfer {
            Some(buffer) => self.parent.post_message_with_transfer(message, &self.origin, &js_sys::Array::of1(buffer)),
            None => self.parent.post_message(message, &self.origin),
        };
        if let Err(e) = result {
            web_sys::console::error_2(&"MatLivre could not reach its host page".into(), &e);
        }
    }
}

fn object(kind: &str) -> js_sys::Object {
    let message = js_sys::Object::new();
    set(&message, "type", &format!("{PREFIX}{kind}").into());
    message
}

fn set(target: &js_sys::Object, key: &str, value: &JsValue) {
    let _ = js_sys::Reflect::set(target, &key.into(), value);
}

fn get(data: &JsValue, key: &str) -> Option<JsValue> {
    js_sys::Reflect::get(data, &key.into()).ok().filter(|v| !v.is_undefined() && !v.is_null())
}

fn read_request(data: &JsValue) -> Option<HostRequest> {
    let kind = get(data, "type")?.as_string()?;
    match kind.strip_prefix(PREFIX)? {
        "open" => {
            let name = get(data, "name")?.as_string()?;
            let bytes = get(data, "bytes")?;
            let bytes = if let Some(buffer) = bytes.dyn_ref::<js_sys::ArrayBuffer>() {
                js_sys::Uint8Array::new(buffer).to_vec()
            } else {
                js_sys::Uint8Array::new(&bytes).to_vec()
            };
            Some(HostRequest::Open { name, bytes })
        }
        "request" => {
            let id = get(data, "id")?;
            let method = get(data, "method")?.as_string()?;
            let params = match get(data, "params") {
                Some(p) => js_sys::JSON::stringify(&p).ok().and_then(|s| s.as_string()).and_then(|s| serde_json::from_str(&s).ok()).unwrap_or(serde_json::Value::Null),
                None => serde_json::Value::Null,
            };
            Some(HostRequest::Request { id, method, params })
        }
        _ => None,
    }
}
