//! Local functions are visible only inside their own file: two files may each define a
//! local function of the same name.

use livre_lang::{Interpreter, OutputKind};

fn run(files: &[(&str, &str)], src: &str) -> Vec<String> {
    let dir = std::env::temp_dir().join(format!("matlivre-local-fns-{}-{}", std::process::id(), src.len()));
    std::fs::create_dir_all(&dir).unwrap();
    for (name, text) in files {
        std::fs::write(dir.join(format!("{name}.m")), text).unwrap();
    }
    let mut it = Interpreter::new();
    it.set_search_path(vec![dir.clone()]);
    let _ = it.run_source(src);
    let _ = std::fs::remove_dir_all(&dir);
    it.take_outputs()
        .into_iter()
        .filter_map(|o| match o.kind {
            OutputKind::Value { text, .. } => Some(text),
            OutputKind::Text(t) => Some(t),
            OutputKind::Error(e) => Some(format!("Error: {e}")),
            _ => None,
        })
        .collect()
}

const FA: &str = "function r = fa()\nr = helper();\nend\nfunction r = helper()\nr = 'helper of fa';\nend\n";
const FB: &str = "function r = fb()\nr = helper();\nend\nfunction r = helper()\nr = 'helper of fb';\nend\n";

#[test]
fn same_named_local_functions_in_two_files() {
    let out = run(&[("fa", FA), ("fb", FB)], "disp(fa()); disp(fb()); disp(fa());");
    assert_eq!(out, ["helper of fa\n", "helper of fb\n", "helper of fa\n"]);
}

#[test]
fn numeric_local_functions_with_loops_and_handles() {
    let ga = "function r = ga(n)\nr = 0;\nfor k = 1:n\n  r = r + h(k);\nend\nend\nfunction y = h(x)\ny = x;\nend\n";
    let gb = "function r = gb(n)\nr = 0;\nfor k = 1:n\n  r = r + h(k);\nend\nr = r + sum(arrayfun(@h, 1:n));\nend\nfunction y = h(x)\ny = 10*x;\nend\n";
    let gc = "function r = gc(x)\nr = ga(x) + gb(x);\nend\n";
    let out = run(&[("ga", ga), ("gb", gb), ("gc", gc)], "a = ga(3)\nb = gb(3)\nc = gc(3)\nd = ga(3)");
    assert_eq!(out, ["a = 6", "b = 120", "c = 126", "d = 6"]);
}
