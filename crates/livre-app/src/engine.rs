//! The engine owns the interpreter (and so the base workspace) and runs apart from the
//! window: on its own thread, or in a Web Worker in the browser. The UI sends run requests
//! and receives outputs as they are produced (spec section 3).

use std::cell::Cell;
use std::rc::Rc;
use std::sync::mpsc::Receiver;
use std::time::Duration;

use livre_lang::clock::Instant;
use livre_lang::{Interpreter, Output, OutputKind, VarInfo};

#[derive(serde::Serialize, serde::Deserialize)]
pub enum Command {
    /// Run `blocks` in order. `all_code` holds every code block of the document so local
    /// functions defined anywhere are callable.
    Run { blocks: Vec<(u64, String)>, all_code: Vec<String>, folder: Option<std::path::PathBuf> },
    /// A command-window line, run in the base workspace without resetting figures or
    /// functions. Its outputs come back under [`COMMAND_ID`].
    Eval(String),
    ClearWorkspace,
}

/// The block id the command window's outputs are reported under.
pub const COMMAND_ID: u64 = u64::MAX;

#[derive(serde::Serialize, serde::Deserialize)]
pub enum Event {
    BlockStarted(u64),
    Output(u64, Output),
    Finished { workspace: Vec<VarInfo>, ok: bool, elapsed: Duration },
}

/// The interpreter and the order it handles commands in, the same on a thread and in a
/// worker; `emit` carries each event to the window.
struct Runner {
    it: Interpreter,
    current: Rc<Cell<u64>>,
    emit: Rc<dyn Fn(Event)>,
}

impl Runner {
    fn new(emit: Rc<dyn Fn(Event)>) -> Self {
        let mut it = Interpreter::new();
        let current = Rc::new(Cell::new(0u64));
        {
            let emit = emit.clone();
            let current = current.clone();
            it.set_output_sink(move |outs| {
                for o in outs {
                    emit(Event::Output(current.get(), o));
                }
            });
        }
        Self { it, current, emit }
    }

    fn handle(&mut self, cmd: Command) {
        match cmd {
            Command::ClearWorkspace => {
                self.it.clear_workspace();
                self.it.close_all_figures();
                (self.emit)(Event::Finished { workspace: self.it.workspace(), ok: true, elapsed: Duration::ZERO });
            }
            Command::Eval(src) => {
                self.current.set(COMMAND_ID);
                let started = Instant::now();
                let ok = match livre_lang::parser::parse(&src) {
                    Ok(prog) => {
                        self.it.define_functions(&prog);
                        self.it.run(&prog.body).is_ok()
                    }
                    Err(e) => {
                        (self.emit)(Event::Output(COMMAND_ID, Output { line: 1, kind: OutputKind::Error(e.message) }));
                        false
                    }
                };
                (self.emit)(Event::Finished { workspace: self.it.workspace(), ok, elapsed: started.elapsed() });
            }
            Command::Run { blocks, all_code, folder } => {
                // The script's folder is on the search path (function and classdef files);
                // resetting it each run picks up edited files.
                self.it.set_search_path(folder.into_iter().collect());
                let started = Instant::now();
                // Each run draws into fresh figures, like re-running a live script.
                self.it.close_all_figures();
                self.it.clear_functions();
                for src in &all_code {
                    if let Ok(p) = livre_lang::parser::parse(src) {
                        self.it.define_functions(&p);
                    }
                }
                let mut ok = true;
                for (id, src) in blocks {
                    self.current.set(id);
                    (self.emit)(Event::BlockStarted(id));
                    match livre_lang::parser::parse(&src) {
                        Ok(prog) => {
                            if self.it.run(&prog.body).is_err() {
                                ok = false;
                                break;
                            }
                        }
                        Err(e) => {
                            let out = Output { line: e.line.max(1), kind: OutputKind::Error(e.message) };
                            (self.emit)(Event::Output(id, out));
                            ok = false;
                            break;
                        }
                    }
                }
                (self.emit)(Event::Finished { workspace: self.it.workspace(), ok, elapsed: started.elapsed() });
            }
        }
    }
}

pub struct Engine {
    pub rx: Receiver<Event>,
    link: Link,
}

// ---- the engine thread -----------------------------------------------------------------

#[cfg(not(target_arch = "wasm32"))]
struct Link {
    tx: std::sync::mpsc::Sender<Command>,
    cancel: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

#[cfg(not(target_arch = "wasm32"))]
impl Engine {
    pub fn start(ctx: eframe::egui::Context) -> Self {
        use std::sync::atomic::AtomicBool;
        use std::sync::mpsc::channel;
        use std::sync::Arc;

        let (tx, cmd_rx) = channel::<Command>();
        let (ev_tx, rx) = channel::<Event>();
        let cancel = Arc::new(AtomicBool::new(false));
        let flag = cancel.clone();
        std::thread::Builder::new()
            .name("livre-engine".into())
            .stack_size(64 * 1024 * 1024)
            .spawn(move || {
                let emit_ctx = ctx.clone();
                let mut runner = Runner::new(Rc::new(move |ev| {
                    let _ = ev_tx.send(ev);
                    emit_ctx.request_repaint();
                }));
                runner.it.set_cancel_flag(flag);
                runner.it.enable_persistent_settings();
                install_exporter(&mut runner.it, ctx);
                while let Ok(cmd) = cmd_rx.recv() {
                    runner.handle(cmd);
                }
            })
            .expect("failed to start engine thread");
        Self { rx, link: Link { tx, cancel } }
    }

    pub fn send(&self, cmd: Command) {
        if matches!(cmd, Command::Run { .. }) {
            self.link.cancel.store(false, std::sync::atomic::Ordering::Relaxed);
        }
        let _ = self.link.tx.send(cmd);
    }

    /// Asks the interpreter to stop between statements; the workspace is kept.
    pub fn stop(&self) {
        self.link.cancel.store(true, std::sync::atomic::Ordering::Relaxed);
    }
}

/// saveas/print/exportgraphics render with the light skin at MATLAB's figure size.
#[cfg(not(target_arch = "wasm32"))]
pub fn install_exporter(it: &mut Interpreter, ctx: eframe::egui::Context) {
    it.set_figure_exporter(move |fig, format, scale| {
        let skin = crate::theme::Skin::light();
        let size = eframe::egui::vec2(560.0, 420.0);
        match format {
            livre_lang::ExportFormat::Svg => Ok(crate::plot::to_svg(&ctx, &skin, fig, size).into_bytes()),
            livre_lang::ExportFormat::Png => crate::plot::to_png(&ctx, &skin, fig, size, scale),
        }
    });
}

// ---- the engine worker -----------------------------------------------------------------

/// The window's end of the Web Worker running the engine. Commands sent before the worker
/// has loaded the module wait in `queued`; a stop ends the worker and starts another, since
/// a running script cannot be interrupted without shared memory.
#[cfg(target_arch = "wasm32")]
struct Link {
    worker: std::cell::RefCell<Option<WorkerLink>>,
    tx: std::sync::mpsc::Sender<Event>,
    ctx: eframe::egui::Context,
}

#[cfg(target_arch = "wasm32")]
struct WorkerLink {
    worker: web_sys::Worker,
    ready: Rc<Cell<bool>>,
    queued: Rc<std::cell::RefCell<Vec<Command>>>,
    _onmessage: wasm_bindgen::closure::Closure<dyn FnMut(web_sys::MessageEvent)>,
}

#[cfg(target_arch = "wasm32")]
impl Engine {
    pub fn start(ctx: eframe::egui::Context) -> Self {
        let (tx, rx) = std::sync::mpsc::channel::<Event>();
        let link = Link { worker: std::cell::RefCell::new(None), tx, ctx };
        link.spawn();
        Self { rx, link }
    }

    pub fn send(&self, cmd: Command) {
        let worker = self.link.worker.borrow();
        let Some(w) = worker.as_ref() else {
            return;
        };
        if w.ready.get() {
            w.post(&cmd);
        } else {
            w.queued.borrow_mut().push(cmd);
        }
    }

    /// Ends the worker mid-run; the workspace it held is gone, and a fresh engine takes over.
    pub fn stop(&self) {
        if let Some(w) = self.link.worker.borrow_mut().take() {
            w.worker.terminate();
        }
        let _ = self.link.tx.send(Event::Finished { workspace: Vec::new(), ok: false, elapsed: Duration::ZERO });
        self.link.spawn();
        self.link.ctx.request_repaint();
    }
}

#[cfg(target_arch = "wasm32")]
impl Link {
    fn spawn(&self) {
        use wasm_bindgen::prelude::*;
        use wasm_bindgen::JsCast;

        let options = web_sys::WorkerOptions::new();
        options.set_type(web_sys::WorkerType::Module);
        let worker = match web_sys::Worker::new_with_options("./engine-worker.js", &options) {
            Ok(w) => w,
            Err(e) => {
                web_sys::console::error_2(&"MatLivre could not start its engine worker".into(), &e);
                return;
            }
        };
        let ready = Rc::new(Cell::new(false));
        let queued: Rc<std::cell::RefCell<Vec<Command>>> = Rc::new(std::cell::RefCell::new(Vec::new()));
        let onmessage = {
            let tx = self.tx.clone();
            let ctx = self.ctx.clone();
            let ready = ready.clone();
            let queued = queued.clone();
            let worker = worker.clone();
            Closure::<dyn FnMut(web_sys::MessageEvent)>::new(move |event: web_sys::MessageEvent| {
                let data = event.data();
                if data.as_string().as_deref() == Some("ready") {
                    ready.set(true);
                    for cmd in queued.borrow_mut().drain(..) {
                        post_command(&worker, &cmd);
                    }
                    return;
                }
                let bytes = js_sys::Uint8Array::new(&data).to_vec();
                match postcard::from_bytes::<Event>(&bytes) {
                    Ok(ev) => {
                        let _ = tx.send(ev);
                        ctx.request_repaint();
                    }
                    Err(e) => web_sys::console::error_1(&format!("MatLivre could not read an engine event: {e}").into()),
                }
            })
        };
        worker.set_onmessage(Some(onmessage.as_ref().unchecked_ref()));
        let init = js_sys::Object::new();
        let (module, wasm) = module_urls();
        let _ = js_sys::Reflect::set(&init, &"type".into(), &"init".into());
        let _ = js_sys::Reflect::set(&init, &"module".into(), &module.into());
        let _ = js_sys::Reflect::set(&init, &"wasm".into(), &wasm.into());
        let _ = worker.post_message(&init);
        *self.worker.borrow_mut() = Some(WorkerLink { worker, ready, queued, _onmessage: onmessage });
    }
}

#[cfg(target_arch = "wasm32")]
impl WorkerLink {
    fn post(&self, cmd: &Command) {
        post_command(&self.worker, cmd);
    }
}

#[cfg(target_arch = "wasm32")]
fn post_command(worker: &web_sys::Worker, cmd: &Command) {
    match postcard::to_allocvec(cmd) {
        Ok(bytes) => {
            let array = js_sys::Uint8Array::from(bytes.as_slice());
            let buffer = array.buffer();
            let transfer = js_sys::Array::of1(&buffer);
            let _ = worker.post_message_with_transfer(&buffer, &transfer);
        }
        Err(e) => web_sys::console::error_1(&format!("MatLivre could not send a command to its engine: {e}").into()),
    }
}

/// The URLs of the window's own module and WebAssembly file, which the worker loads too.
/// The bundler names them with a hash, so they are read off the page: the preload links it
/// writes, or else the import in its loading script.
#[cfg(target_arch = "wasm32")]
fn module_urls() -> (String, String) {
    let href = |selector: &str| -> Option<String> {
        let document = web_sys::window()?.document()?;
        let link = document.query_selector(selector).ok()??;
        link.get_attribute("href")
    };
    let module = href("link[rel='modulepreload']");
    let wasm = href("link[rel='preload'][as='fetch']");
    if let (Some(m), Some(w)) = (module, wasm) {
        return (m, w);
    }
    let script = web_sys::window()
        .and_then(|w| w.document())
        .and_then(|d| d.query_selector("script[type='module']").ok().flatten())
        .and_then(|s| s.text_content())
        .unwrap_or_default();
    let quoted = |after: &str| -> Option<String> {
        let rest = &script[script.find(after)? + after.len()..];
        let rest = rest.trim_start_matches(['"', '\'']);
        let end = rest.find(['"', '\''])?;
        Some(rest[..end].to_string())
    };
    (quoted("from ").unwrap_or_default(), quoted("module_or_path: ").unwrap_or_default())
}

/// The worker's end: it runs the interpreter in the worker's own copy of the module. The
/// worker script calls this once the module is loaded; from then on every message is a
/// command, and every event goes back to the window.
#[cfg(target_arch = "wasm32")]
pub fn worker_main() {
    use wasm_bindgen::prelude::*;
    use wasm_bindgen::JsCast;

    let scope: web_sys::DedicatedWorkerGlobalScope = js_sys::global().unchecked_into();
    let post = {
        let scope = scope.clone();
        Rc::new(move |ev: Event| match postcard::to_allocvec(&ev) {
            Ok(bytes) => {
                let array = js_sys::Uint8Array::from(bytes.as_slice());
                let buffer = array.buffer();
                let transfer = js_sys::Array::of1(&buffer);
                let _ = scope.post_message_with_transfer(&buffer, &transfer);
            }
            Err(e) => web_sys::console::error_1(&format!("MatLivre could not send an engine event: {e}").into()),
        })
    };
    let runner = Rc::new(std::cell::RefCell::new(Runner::new(post)));
    let onmessage = Closure::<dyn FnMut(web_sys::MessageEvent)>::new(move |event: web_sys::MessageEvent| {
        let bytes = js_sys::Uint8Array::new(&event.data()).to_vec();
        match postcard::from_bytes::<Command>(&bytes) {
            Ok(cmd) => runner.borrow_mut().handle(cmd),
            Err(e) => web_sys::console::error_1(&format!("MatLivre could not read a command: {e}").into()),
        }
    });
    scope.set_onmessage(Some(onmessage.as_ref().unchecked_ref()));
    onmessage.forget();
    let _ = scope.post_message(&JsValue::from_str("ready"));
}
