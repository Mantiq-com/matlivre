# MatLivre in the browser

A build of the editor for a web page, made with [trunk](https://trunkrs.dev) from
`index.html`. The page draws the editor on a canvas; the engine runs in a Web Worker
(`engine-worker.js`) that loads the same module, so a long script never freezes the editor.
The interpreter runs without the native compiler, threads, the clipboard, file functions or
figure export, which a browser does not give a program.

## Source and retained browser changes

The upstream target is [`master` at `e66a9f6d67c0ffb015c19d967847835d0e94206b`](https://github.com/campuslivre/matlivre/commit/e66a9f6d67c0ffb015c19d967847835d0e94206b).
The previous browser release is `host-bridge-2` at
`65799972e1abeeeac73162fa4901ee94ca3a424a`. That published tag contains changes beyond
the older `web-host` branch and supplies the browser changes carried here.

| Fork change | Required behavior | Selected upstream evidence | Decision |
| --- | --- | --- | --- |
| `livre-lang` feature gates and `clock` module | Run the interpreter in WebAssembly without native compilation, native threads or process clocks | The upstream manifest requires native dependencies, and its interpreter uses `std::time::Instant` | Keep the browser feature selection and clocks |
| `livre-doc/src/mlx.rs` and its browser clock dependency | Save an MLX package in the browser | Upstream timestamps packages with `std::time::SystemTime` | Keep the browser clock; preserve the package format |
| `livre-app/src/lib.rs`, platform gates and `engine.rs` | Draw the same editor on desktop and in a browser, with a worker that can be stopped | Upstream provides a desktop binary and a thread-based engine | Keep the shared library and worker transport, including variable inspection |
| `livre-app/src/host.rs` and `apps/matlivre-web` | Exchange MLX/M files, commands, dirty state and saved bytes with a same-origin parent | The upstream checkout has no browser application or parent protocol | Keep the host bridge and browser package |
| `.github/workflows/web.yml` | Publish an immutable browser archive without replacing an existing release's bytes | Upstream release workflows package desktop applications | Keep the browser packaging workflow |

An upstream replacement must provide these browser behaviors and pass the same hosted
open, edit, save, reopen, export, worker-stop and host-validation checks. A desktop
build or native interpreter test alone does not qualify a replacement.

This upstream checkout declares workspace version `0.1.6`. The ready message therefore
cannot identify this source revision on its own; identify an installed bundle by its
source commit and archive checksum as well.

### Validation of the release-based snapshot

The following results apply to the browser adaptation based on upstream v0.1.8,
commit `53089b98ec340a7c488c0eab00c9a816df4364e3`. They do not qualify the current
`upstream/master` base, which still needs build and browser checks.

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
