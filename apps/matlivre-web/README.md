# MatLivre in the browser

A build of the editor for a web page, made with [trunk](https://trunkrs.dev) from
`index.html`. The page draws the editor on a canvas; the engine runs in a Web Worker
(`engine-worker.js`) that loads the same module, so a long script never freezes the editor.
The interpreter runs without the native compiler, threads, the clipboard, file functions or
figure export, which a browser does not give a program.

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
