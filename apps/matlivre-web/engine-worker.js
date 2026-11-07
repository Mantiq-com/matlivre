// The engine's Web Worker. It loads the same WebAssembly module as the window, so the
// interpreter runs apart from the editor and a long script never freezes it. The window
// sends the module's URLs first; every message after the module has answered "ready" is a
// command, handled by the module itself.
self.onmessage = async (event) => {
    const data = event.data;
    if (!data || data.type !== 'init') return;
    const module = await import(data.module);
    await module.default({ module_or_path: data.wasm });
    module.engine_worker_entry();
};
