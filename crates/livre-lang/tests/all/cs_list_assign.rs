//! Multiple assignment into a comma-separated list, `[varargout{1:nargout}] = f(...)`: output k
//! goes to element k (it used to give every element the first output). Expected values from
//! MATLAB R2023b; the script forms are in `tests/matlab/cases/basic/cs_list_assign.m`.

use livre_lang::{Interpreter, OutputKind};

fn run(files: &[(&str, &str)], src: &str) -> Vec<String> {
    let dir = std::env::temp_dir().join(format!("matlivre-cs-list-{}-{}", std::process::id(), src.len()));
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

const VOUT: &str = "function varargout = vout()\n[varargout{1:nargout}] = size(zeros(2,3));\nend\n";
const FIELDS: &str = "function s = fields3()\ns = struct('f', {0, 0, 0});\n[s.f] = deal('a', 'b', 'c');\nend\n";

#[test]
fn varargout_takes_one_output_per_element() {
    let out = run(&[("vout", VOUT)], "[a, b] = vout();\ndisp(a), disp(b)\n[a, b, c] = vout();\ndisp([a b c])\na = vout();\ndisp(a)");
    // One output: varargout{1} takes size's single output, [2 3].
    assert_eq!(out, ["     2\n", "     3\n", "     2     3     1\n", "     2     3\n"]);
}

#[test]
fn a_list_target_is_shown_once() {
    let out = run(&[], "[c{1:2}] = deal(1, 2)");
    assert_eq!(out, ["c = 1×2 cell\n\n    {[1]}    {[2]}"]);
}

#[test]
fn struct_array_fields_inside_a_function() {
    let out = run(&[("fields3", FIELDS)], "s = fields3();\ndisp([s.f])");
    assert_eq!(out, ["abc\n"]);
}

#[test]
fn script_forms() {
    let out = run(
        &[],
        concat!(
            "[c{1:2}] = deal(1, 2); disp(c)\n",
            "[e{1:2}] = size(zeros(2, 3)); disp(e)\n",
            "d = cell(1, 3); [d{:}] = deal(1, 2, 3); disp(d)\n",
            "[x, g{2:3}, y] = deal(1, 2, 3, 4); disp(g)\n",
            "s = struct('f', {1, 2, 3}); [s(2:3).f] = deal(20, 30); disp([s.f])\n",
        ),
    );
    assert_eq!(
        out,
        [
            "    {[1]}    {[2]}\n",
            "    {[2]}    {[3]}\n",
            "    {[1]}    {[2]}    {[3]}\n",
            "    {0×0 double}    {[2]}    {[3]}\n",
            "     1    20    30\n",
        ]
    );
}
