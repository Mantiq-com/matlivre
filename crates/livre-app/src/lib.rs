//! The MatLivre editor: the document view, the engine and the window chrome, shared by
//! the desktop program (`src/main.rs`) and the browser build (`apps/matlivre-web`).

pub mod app;
pub mod engine;
pub mod highlight;
pub mod icons;
pub mod plot;
pub mod theme;
pub mod widgets;

#[cfg(not(target_arch = "wasm32"))]
mod chat;
#[cfg(not(target_arch = "wasm32"))]
mod claude;
#[cfg(not(target_arch = "wasm32"))]
mod complete;
#[cfg(not(target_arch = "wasm32"))]
mod files;
#[cfg(target_os = "macos")]
mod macos;
#[cfg(not(target_arch = "wasm32"))]
pub mod mcp;
#[cfg(not(target_arch = "wasm32"))]
mod os;
/// The browser build has no native window; only the shortcut labels remain, shown as written.
#[cfg(target_arch = "wasm32")]
mod os {
    pub fn keys(s: &str) -> String {
        s.to_string()
    }
}
#[cfg(not(target_arch = "wasm32"))]
mod transfer;
#[cfg(not(target_arch = "wasm32"))]
pub mod updater;

#[cfg(target_arch = "wasm32")]
pub mod host;

pub use app::{App, DEFAULT_SIZE, MIN_SIZE};

/// What the window starts with.
pub enum Launch {
    Welcome,
    #[cfg(not(target_arch = "wasm32"))]
    Open(std::path::PathBuf),
    /// A tab dragged out of another window (a manifest from `transfer::write`).
    #[cfg(not(target_arch = "wasm32"))]
    Tab(std::path::PathBuf),
}
