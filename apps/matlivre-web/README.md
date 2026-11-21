# MatLivre in the browser

A build of the editor for a web page, made with [trunk](https://trunkrs.dev) from
`index.html`. The page draws the editor on a canvas; the engine runs in a Web Worker
(`engine-worker.js`) that loads the same module, so a long script never freezes the editor.
The module starts the page UI only in a window; the worker starts the interpreter
through its exported entry point.
The interpreter runs without the native compiler, threads, the clipboard, file functions or
figure export, which a browser does not give a program.

## Source and retained browser changes

The upstream target is [`master` at `e099eee04472959ec5f40d0d441c75ab25cc209b`](https://github.com/campuslivre/matlivre/commit/e099eee04472959ec5f40d0d441c75ab25cc209b),
retrieved on 2026-11-21. The latest stable release at retrieval is
[v0.2.0](https://github.com/campuslivre/matlivre/releases/tag/v0.2.0), commit
`29ed3a19f7f14a1e34aac8d7dbf1ff7e2cc9daa3`. This fork continues tracking
the upstream default branch.
The previous published browser release is `host-bridge-2` at
`65799972e1abeeeac73162fa4901ee94ca3a424a`. That published tag contains changes beyond
the older `web-host` branch and supplies the browser changes carried here.

| Fork change | Required behavior | Selected upstream evidence | Decision |
| --- | --- | --- | --- |
| `livre-lang` feature gates and `clock` module | Run the interpreter in WebAssembly without native compilation, native threads or process clocks | The upstream manifest requires native dependencies, and its interpreter uses `std::time::Instant` | Keep the browser feature selection and clocks |
| `livre-doc/src/mlx.rs` and its browser clock dependency | Save an MLX package in the browser | Upstream timestamps packages with `std::time::SystemTime` | Keep the browser clock; preserve the package format |
| `livre-app/src/lib.rs`, platform gates and `engine.rs` | Draw the same editor on desktop and in a browser, with a worker that can be stopped | Upstream provides a desktop binary and a thread-based engine | Keep the shared library and worker transport, including variable inspection |
| `livre-app/src/host.rs` and `apps/matlivre-web` | Exchange MLX/M files, commands, dirty state and saved bytes with a same-origin parent | The upstream checkout has no browser application or parent protocol | Keep the host bridge and browser package |
| `.github/workflows/web.yml` | Publish an immutable browser archive without replacing an existing release's bytes | Upstream release workflows package desktop applications | Keep the browser packaging workflow |
| `signal_known_deviations.tsv`, case 6948 | Require the complex FIR callback result to match or track platform differences cleanly | Upstream now marks case 6948 `platform-dependent:` in `signal_known_deviations.tsv` and handles it in `signal_matlab.rs` | Removed the fork override; track upstream's file directly |

An upstream replacement must provide these browser behaviors and pass the same hosted
open, edit, save, reopen, export, worker-stop and host-validation checks. A desktop
build or native interpreter test alone does not qualify a replacement.

This upstream checkout still declares workspace version `0.1.6` even after the `v0.2.0`
tag. The ready message therefore cannot identify this source revision on its own; identify
an installed bundle by its source commit and archive checksum as well.

### Browser validation on the selected source

The release package builds from source commit `735772695c5de2d150d30e8372a5fb941525a7d3`
with Rust 1.99.1, Trunk 0.22.2 and wasm-bindgen 0.2.129. The local archive
`matlivre-web-upstream-e099eee.tar.gz` has SHA-256
`1479c77f5552a631e6c37252d30c737468d0985b8003c12a5b2a958d191cb677` and serves
`matlivre-web-5b32c1c0e12b5093_bg.wasm`. The build reports 24 warnings in unchanged
upstream `livre-lang` code.

Hosted MLX and M creation, editing, Save As, ordinary save, export, `ui.hostTheme`
and reopening pass. Saved MLX packages (2,475 bytes) preserve their title, explanatory
text and code blocks. Export keeps the source filename and its unsaved changes. The new
bundle and the published `host-bridge-2` bundle read each other's saved MLX files with
all blocks preserved.

Three fresh workers start cleanly with no startup panic. In the hosted editor UI,
running `disp(6 * 7);` returns `42`, running an infinite loop blocks concurrent host
edits, pressing Stop terminates the worker and starts a replacement worker, and the
replacement worker returns `56` and completes a 65,536-point FFT of ones with zero
error. All observed runtime resources load from the host's own origin under a policy
that blocks external connections, with zero console or page errors. The host rejects
messages from another origin or source window; disposal unloads the iframe and rejects
further requests. These synthetic checks do not qualify authenticated application
storage, autosave acknowledgements or permissions.

The native workspace test suite (`cargo test --locked --workspace -- --test-threads=2`)
passes all unit and sample tests (`livre-app`: 39 passed, 4 ignored; `matlivre`: 2
passed; `livre-doc`: 12 passed + 1 sample test; `livre-lang` lib: 17 passed) and finishes
`livre-lang`'s integration suite without stack overflow: 412 passed, 2 failed, 2 ignored.
The two failures are:

- `signal_matlab::signal_functions_match_matlab` (`crates/livre-lang/tests/all/signal_matlab.rs:203`):
  2 of 10,311 cases differ from MATLAB (9,765 match): case 2138 (`ellip(7,0.1,20,0.45,'high')`)
  and case 8269 (`yulewalk(11,[0 .6 .65 1],[0.001 0.001 1 1])`).
- `stats_matlab::stats_functions_match_matlab` (`crates/livre-lang/tests/all/stats_matlab.rs:204`):
  1 of 4,672 cases differs from MATLAB (4,612 match): case 3862 (`gevfit(xg*5)`).

An independent run on the unmodified selected upstream commit
`e099eee04472959ec5f40d0d441c75ab25cc209b` reproduces both failures with the exact same
cases and counts. All fixtures and assertions remain unchanged.

`cargo build --locked --workspace --bins` and `python3 tools/linux/smoke.py` pass,
including CLI linear-system solving, failing script exit status and headless MCP PNG/SVG
exports.

### Previous validation on the 2026-11-17 master revision

The following results apply to upstream `87764fc0fbdae0999eabe163cfd10d022d5c4e8f`
and do not qualify the newer source above.

Upstream rewrote its history between 2026-11-15 and 2026-11-17: the old base
`e66a9f6d67c0ffb015c19d967847835d0e94206b` and its replacement `8fe0390` have the
identical tree `919a448a36da57e4922ff4b9430b45d05a156d54`.

The release package builds with Rust 1.99.1 and wasm-bindgen 0.2.129. The
`jit_log` field follows the same `jit` feature gate as its initializer, so
the browser build excludes the native compiler's diagnostics state. The
build reports eleven warnings in unchanged upstream signal-processing code.

Hosted MLX and M creation, editing, Save As, ordinary save, export and reopening
pass. Saved MLX packages preserve their title, explanatory text and code blocks.
Export keeps the source filename and its unsaved changes. The new bundle and
the published `host-bridge-2` bundle read each other's saved MLX files with the
same blocks.

Three fresh workers each return 42 for `disp(6 * 7);` with no startup failure.
The Stop shortcut terminates an infinite loop, and a subsequent run returns 42
in the editor. All observed runtime resources load from the host's own origin
under a policy that blocks external connections. The host rejects messages from
another origin or source window; disposal unloads the iframe and rejects further
requests. These synthetic checks do not qualify authenticated application storage,
autosave acknowledgements or permissions.

The native workspace test command,
`cargo test --locked --workspace -- --test-threads=2`, initially stops with
403 passes, one failure and one ignored test in the language integration suite.
Case 6948 matches the recorded MATLAB result but remains listed as a known
deviation; the unmodified upstream target produces the same failure and counts.
Removing only that exception makes the test require a match for that case.

Three full reruns after the removal remain unsuccessful. Two complete with
403 passes, one failure and one ignored test: case 2138, an elliptic high-pass
filter, returns two conjugate zero pairs in a different order from the MATLAB
fixture. The pristine upstream CLI reproduces that ordering. One run
aborts with a stack overflow before the fixture assertions; its cause remains
unresolved. A debugger run of the unmodified upstream suite does not reproduce
the overflow and stops at the resolved exception. The case 2138 input, recorded
output and assertions remain unchanged.

`cargo build --locked --workspace --bins` and `python tools/linux/smoke.py`
pass, including CLI success and failure exits and headless MCP PNG/SVG exports.

### Previous validation on the 2026-11-15 master revision

The following results apply to upstream `e66a9f6d67c0ffb015c19d967847835d0e94206b`
and do not qualify the newer source above.

The browser build passes hosted MLX and M creation, editing, Save As, export and
reopening checks. MLX files pass in both directions with `host-bridge-2`. An export
keeps the source name and its unsaved changes. Messages from another origin or
source window are ignored, and disposing the host bridge unloads the frame.

The shared module's entry point now starts the UI only in a window. Loading it in
a worker previously also attempted to access a page document and panicked; the
published `host-bridge-2` bundle reproduces that failure. Three fresh worker starts
with the fix return a calculation without a startup error. The editor's Stop
control terminates an infinite loop, and its replacement worker runs another
calculation. A 65,536-point FFT of a vector of ones returns the expected result
with zero error. These checks use synthetic files and the actual host bridge;
they do not establish authenticated storage or permissions in a consuming app.

The native workspace gate is not green. Case 2138 of
`signal_matlab::signal_functions_match_matlab`, an elliptic high-pass filter,
differs from the MATLAB fixture in both this fork and the unmodified upstream
commit. Three repeated language-suite runs on each revision all abort with a
stack overflow after reporting that mismatch. An earlier upstream run and a
fork run under a debugger completed with only the filter failure, so the
overflow is intermittent and also affects upstream. Its cause remains
unresolved. No fixture or assertion was changed.

### Validation of the release-based snapshot

The following results apply to the browser adaptation based on upstream v0.1.8,
commit `53089b98ec340a7c488c0eab00c9a816df4364e3`. They do not qualify the current
`upstream/master` base; the checks above apply to that base.

That browser build passed hosted checks for MLX and M opening, edits, creation,
Save As, export and reopening. A file saved by this build also reopens in
`host-bridge-2`. The worker runs calculations, returns variable inspection data,
and restarts after stopping an infinite loop. The host rejects messages with a
different origin or source window. These checks use synthetic documents and
local assets under a policy that blocks external connections; they do not
establish authenticated storage, autosave or permission behavior in a consuming app.

The native workspace suite failed in both the release-based fork and unmodified
upstream v0.1.8. Case 6948 in `signal_known_deviations.tsv` agreed with the
MATLAB fixture, while `signal_matlab.rs` requires every listed deviation to
remain a mismatch. The fixture and assertion remain unchanged. The new
`variable_inspection_survives_worker_transport` test passes.
`cargo build --locked --workspace --bins` and `tools/linux/smoke.py` also pass,
including CLI success and failure exits and headless MCP PNG/SVG output.

```bash
rustup target add wasm32-unknown-unknown
cargo install trunk
cd apps/matlivre-web
trunk serve          # http://127.0.0.1:8080, the editor on its own
./package.sh host-bridge-1   # dist/ and pkg/matlivre-web-host-bridge-1.tar.gz
```

Opened on its own, the page shows the welcome script and keeps nothing. Opened by another
page in a frame as `index.html?host=parent`, it follows that page's instructions: the page
sends the document and asks for its bytes back, through `postMessage` on the same origin.

## The host protocol

Every message is an object with a `type` under the prefix `matlivre:`.

| From the editor | Meaning |
| --- | --- |
| `matlivre:ready { version }` | The editor runs and takes messages. |
| `matlivre:failed { error }` | The editor could not start. |
| `matlivre:dirty { dirty }` | The document gained or lost unsaved changes. |
| `matlivre:save { name, bytes, saveAs }` | The document's bytes to store under `name` (`bytes` is an `ArrayBuffer`); `saveAs` is true when the document took this name just now. |
| `matlivre:response { id, reply }` | The answer to a request: `{ ok: true, result }` or `{ ok: false, error }`. |

| From the page | Meaning |
| --- | --- |
| `matlivre:open { name, bytes }` | Show this document: a `.mlx` package, or the text of a `.m` script, by the name's extension. |
| `matlivre:request { id, method, params }` | Call `document.inspect`, `engine.commands` or `engine.execute { command, params }`. |

`document.inspect` answers `{ name, kind, blocks: [{ index, kind, text }] }`, `kind` being
`live` or `script` and a block's kind `title`, `heading`, `heading2`, `text`, `code` or
`sectionBreak`. `engine.commands` lists the document commands with their parameters.
`engine.execute` runs one command:

| Command | Effect |
| --- | --- |
| `file.new { template }` | A new, empty live script. |
| `file.save {}` | Sends `matlivre:save` under the current name. |
| `file.saveAs { path, keepSource }` | Writes `.m` or `.mlx` by the extension of `path` and sends `matlivre:save` under that name. The document takes the name, with `saveAs: true`, unless `keepSource` is true: then the copy is an export and the document keeps its name and its unsaved changes. |
| `file.info {}` | `{ name, kind, dirty }`. |
| `ui.hostTheme { appearance }` | The light or dark skin. |
| `block.setText { index, text }` | Replaces one block's text. |
| `block.insert { at, kind, text }` | Inserts a block. |
| `block.delete { index }` | Removes a block. |
| `text.set { text }` | Replaces the whole document from its `.m` form. |
| `text.replaceAll { text, with }` | Replaces every occurrence in the `.m` form; answers `{ replaced }`. |

The page's own file menu offers Save only while a page hosts the editor: the page names
copies and exports, creates and opens documents.
