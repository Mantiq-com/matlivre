#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result<()> {
    use eframe::egui;
    use livre_app::{app, icons, mcp, updater, App, Launch};

    // Velopack's install, update and uninstall hooks: answered and quit before anything else.
    velopack::VelopackApp::build().run();
    // An MCP server over stdio for Claude Code (`claude mcp add matlivre -- matlivre --mcp`).
    if std::env::args().nth(1).as_deref() == Some("--mcp") {
        mcp::stdio();
        return Ok(());
    }
    updater::check_in_background();
    let mut launch = Launch::Welcome;
    let mut at = None;
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--tab" => launch = args.next().map(|p| Launch::Tab(p.into())).unwrap_or(Launch::Welcome),
            "--at" => {
                at = args.next().and_then(|s| {
                    let (x, y) = s.split_once(',')?;
                    Some(egui::pos2(x.parse().ok()?, y.parse().ok()?))
                })
            }
            _ => launch = Launch::Open(a.into()),
        }
    }
    let mut viewport = egui::ViewportBuilder::default()
        .with_title("MatLivre")
        .with_inner_size(app::DEFAULT_SIZE)
        .with_min_inner_size(app::MIN_SIZE)
        .with_icon(std::sync::Arc::new(icons::window_icon()))
        .with_decorations(false)
        .with_resizable(true);
    if let Some(p) = at {
        viewport = viewport.with_position(p);
    }
    // A tab torn off to a new window opens where it was dropped; otherwise the window is
    // centred on the screen the pointer is on.
    let place = at.is_none();
    let options = eframe::NativeOptions { viewport, ..Default::default() };
    eframe::run_native("matlivre", options, Box::new(move |cc| Ok(Box::new(App::new(cc, launch, place)))))
}

/// The desktop program has no WebAssembly build; `apps/matlivre-web` is the browser one.
#[cfg(target_arch = "wasm32")]
fn main() {}
