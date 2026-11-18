"""Builds tests/all/data/stats_cases.m from tools/stats/cases.txt.

The script defines a few test data sets, runs every case inside try/catch and prints each result
with %.17g (class, size, real parts, imaginary parts) or the error message. The same script runs
in MATLAB, to record stats_matlab.txt, and in MatLivre, from the test stats_matlab.rs:

    python tools/stats/gen_cases.py
    matlab -batch "run('crates/livre-lang/tests/all/data/stats_cases.m')" > out.txt
    # keep only the result lines (those starting with <case>[.<n>]|) as stats_matlab.txt

Cases may assign any variable except the ones the loop and the data sets use (ci_, ki_, it_, w_,
w2_, ME_, cases, r, x, y, A, P, N, v).
"""
import os

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = os.path.join(root, "tools", "stats", "cases.txt")
dst = os.path.join(root, "crates", "livre-lang", "tests", "all", "data", "stats_cases.m")

cases = [line.rstrip("\n") for line in open(src, encoding="utf-8")]
cases = [c for c in cases if c.strip() and not c.startswith("#")]

LOOP = r"""for ci_ = 1:numel(cases)
  try
    r = 'unset';
    eval(cases{ci_});
    if iscell(r), it_ = r; else, it_ = {r}; end
    for ki_ = 1:numel(it_)
      w_ = it_{ki_};
      if isnumeric(w_) || islogical(w_) || ischar(w_)
        w2_ = double(w_(:)).';
        fprintf('%d.%d|%s|%s|%s|%s\n', ci_, ki_, class(w_), sprintf('%d ', size(w_)), sprintf('%.17g ', real(w2_)), sprintf('%.17g ', imag(w2_)));
      else
        fprintf('%d.%d|%s\n', ci_, ki_, class(w_));
      end
    end
  catch ME_
    fprintf('%d|ERR|%s\n', ci_, strrep(ME_.message, char(10), ' '));
  end
end"""

out = [
    "x = sin(0.3*(1:100)) + 0.5*cos(1.1*(1:100));",
    "y = cos(0.9*(1:100)) + 0.3*sin(2.4*(1:100));",
    "A = reshape(round(1e4*(sin((1:30).^1.5)*3+2))/1e4, 6, 5);",
    "P = abs(A) + 0.5;",
    "N = A; N(2,1) = NaN; N(5,3) = NaN; N(1,5) = NaN;",
    "v = [2 4 4 4 5 5 7 9]';",
    "cases = {",
]
out += ["'" + c.replace("'", "''") + "'" for c in cases]
out.append("};")
out.append(LOOP)
with open(dst, "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(out) + "\n")
print(len(cases), "cases ->", dst)
