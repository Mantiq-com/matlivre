//! The browser build: the editor drawn on the page's canvas, and the entry the engine's
//! worker calls. Built with trunk from `index.html`; see `README.md` in this directory.

#[cfg(target_arch = "wasm32")]
pub mod web {
    use wasm_bindgen::prelude::*;
    use wasm_bindgen::JsCast;

    /// Called by `engine-worker.js` once the worker has loaded this module: from then on
    /// the worker runs the interpreter.
    #[wasm_bindgen]
    pub fn engine_worker_entry() {
        livre_app::engine::worker_main();
    }

    pub fn start() {
        wasm_bindgen_futures::spawn_local(async {
            let document = web_sys::window().and_then(|w| w.document()).expect("the page has a document");
            let canvas = document
                .get_element_by_id("matlivre_canvas")
                .and_then(|c| c.dyn_into::<web_sys::HtmlCanvasElement>().ok())
                .expect("the page has the matlivre_canvas canvas");
            let result = eframe::WebRunner::new()
                .start(
                    canvas,
                    eframe::WebOptions::default(),
                    Box::new(|cc| Ok(Box::new(livre_app::App::new(cc, livre_app::Launch::Welcome, false)))),
                )
                .await;
            let loading = document.get_element_by_id("matlivre_loading");
            match result {
                Ok(()) => {
                    if let Some(l) = loading {
                        l.remove();
                    }
                }
                Err(e) => {
                    let message = e.as_string().unwrap_or_else(|| format!("{e:?}"));
                    if let Some(l) = loading {
                        l.set_text_content(Some(&format!("MatLivre could not start: {message}")));
                    }
                    web_sys::console::error_1(&e);
                    failed(&message);
                }
            }
        });
    }

    /// Tells a hosting page that the editor did not start, in the host protocol's words.
    fn failed(message: &str) {
        let Some(window) = web_sys::window() else { return };
        let Ok(Some(parent)) = window.parent() else { return };
        if parent == window {
            return;
        }
        let Ok(origin) = window.location().origin() else { return };
        let data = js_sys::Object::new();
        let _ = js_sys::Reflect::set(&data, &"type".into(), &"matlivre:failed".into());
        let _ = js_sys::Reflect::set(&data, &"error".into(), &message.into());
        let _ = parent.post_message(&data, &origin);
    }
}

#[cfg(target_arch = "wasm32")]
fn main() {
    // The engine worker loads this module too, but has no page to draw.
    if web_sys::window().is_some() {
        web::start();
    }
}

/// Only the WebAssembly target has a browser to run in.
#[cfg(not(target_arch = "wasm32"))]
fn main() {}
