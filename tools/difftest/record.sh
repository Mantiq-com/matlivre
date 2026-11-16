#!/usr/bin/env bash
# Record MATLAB's snapshots of the difftest cases into
# tests/matlab/expected/<release>/<case>.json (the release is read from the MATLAB that ran them).
#
#   tools/difftest/record.sh                 record every case under tests/matlab/cases
#   tools/difftest/record.sh tests/matlab/cases/basic/workspace.m ...
#   tools/difftest/record.sh --release R2023b ...   refuse to record unless MATLAB is that release
#
# One MATLAB session runs all the cases; difftest_run.m resets the state between them.
#
# MATLAB runs on this machine unless MATLAB_SSH_HOST is set:
#   local:  MATLAB_BIN, the MATLAB executable (default: `matlab` on PATH, else the R2023b
#           install under C:\Program Files on Windows)
#   remote: the same settings as tools/matlab-remote.sh (MATLAB_SSH_HOST, MATLAB_REMOTE_DIR,
#           MATLAB_REMOTE_BIN, MATLAB_REMOTE_OS = windows|posix, default posix); see
#           AGENTS.local.md.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"

want_release=""
if [ "${1:-}" = "--release" ]; then want_release="$2"; shift 2; fi

if [ $# -gt 0 ]; then cases=("$@"); else
  mapfile -t cases < <(cd "$root" && find tests/matlab/cases -name '*.m' | sort)
fi
[ ${#cases[@]} -gt 0 ] || { echo "no cases"; exit 1; }

local_tmp="$(mktemp -d)"
host="${MATLAB_SSH_HOST:-}"
work="dt_$$"

if [ -z "$host" ]; then
  # MATLAB on this machine: the work folder is a local temporary one.
  where="this machine"
  bin="${MATLAB_BIN:-}"
  if [ -z "$bin" ]; then
    if command -v matlab >/dev/null 2>&1; then bin=matlab
    elif [ -x "/c/Program Files/MATLAB/R2023b/bin/matlab.exe" ]; then bin="/c/Program Files/MATLAB/R2023b/bin/matlab"
    else echo "MATLAB not found: set MATLAB_BIN, or MATLAB_SSH_HOST for a remote one" >&2; exit 1
    fi
  fi
  wdir="$local_tmp/$work"
  mkdir -p "$wdir"
  # MATLAB on Windows wants C:/... paths, not Git Bash's /c/...
  mpath() { if command -v cygpath >/dev/null 2>&1; then cygpath -m "$1"; else printf '%s' "$1"; fi; }
  put() { cp "$1" "$wdir/$2"; }
  get() { cp "$wdir/$1" "$2" 2>/dev/null; }
  run_driver() { "$bin" -batch "run('$(mpath "$wdir")/dt_driver.m')"; }
  cleanup() { :; }
else
  where="$host"
  rdir="${MATLAB_REMOTE_DIR:-/tmp}"
  bin="${MATLAB_REMOTE_BIN:-matlab}"
  os="${MATLAB_REMOTE_OS:-posix}"
  quiet() { grep -v -e 'post-quantum' -e 'store now' -e 'openssh.com/pq' || true; }
  rsh() { ssh "$host" "$@" 2> >(quiet >&2); }
  case "$rdir" in /*) scp_dir="$rdir" ;; *) scp_dir="/$rdir" ;; esac
  if [ "$os" = windows ]; then
    wdir_win="${rdir//\//\\}\\$work"
    rsh "mkdir \"$wdir_win\"" >/dev/null
    cleanup() { rsh "rmdir /s /q \"$wdir_win\"" >/dev/null 2>&1 || true; }
  else
    rsh "mkdir -p \"$rdir/$work\"" >/dev/null
    cleanup() { rsh "rm -rf \"$rdir/$work\"" >/dev/null 2>&1 || true; }
  fi
  put() { scp -q "$1" "$host:$scp_dir/$work/$2"; }
  get() { scp -q "$host:$scp_dir/$work/$1" "$2" 2>/dev/null; }
  # The script's own folder is the current folder while it runs, so difftest_run.m is found.
  run_driver() { rsh "\"$bin\" -batch \"run('$rdir/$work/dt_driver.m')\""; }
fi
trap 'cleanup; rm -rf "$local_tmp"' EXIT

# Send the harness and the cases (c001.m, c002.m, ...), plus a driver that runs them in order.
put "$here/difftest_run.m" difftest_run.m
driver="$local_tmp/dt_driver.m"
{
  echo "here = fileparts(mfilename('fullpath')); addpath(here);"
  echo "fid = fopen(fullfile(here,'release.txt'),'w'); fprintf(fid,'%s',version('-release')); fclose(fid);"
} > "$driver"
i=0
for c in "${cases[@]}"; do
  i=$((i+1)); id=$(printf 'c%03d' "$i")
  put "$root/$c" "$id.m"
  echo "try, cd(tempdir); difftest_run(fullfile(here,'$id.m'), fullfile(here,'$id.json')); catch ME, disp(['FAILED $id: ' ME.message]); end" >> "$driver"
done
put "$driver" dt_driver.m

run_driver

get release.txt "$local_tmp/release.txt"
release="R$(tr -d '\r\n' < "$local_tmp/release.txt")"
if [ -n "$want_release" ] && [ "$want_release" != "$release" ]; then
  echo "MATLAB on $where is $release, not $want_release: nothing recorded" >&2; exit 1
fi
echo "Recording against MATLAB $release"

i=0
for c in "${cases[@]}"; do
  i=$((i+1)); id=$(printf 'c%03d' "$i")
  rel="${c#tests/matlab/cases/}"; rel="${rel%.m}"
  out="$root/tests/matlab/expected/$release/$rel.json"
  mkdir -p "$(dirname "$out")"
  if get "$id.json" "$out"; then
    tr -d '\r' < "$out" > "$out.tmp" && mv "$out.tmp" "$out"
    echo "  $rel"
  else
    echo "  $rel: no snapshot produced" >&2
  fi
done
