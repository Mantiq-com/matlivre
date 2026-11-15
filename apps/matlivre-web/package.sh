#!/usr/bin/env bash
# Builds the browser bundle and packs it as the artifact an embedding application installs:
# matlivre-web-<tag>/{index.html, engine-worker.js, …} inside matlivre-web-<tag>.tar.gz, in
# this directory's pkg/. Prints the archive's SHA-256, which the application pins.
set -euo pipefail

tag="${1:?usage: package.sh <release-tag>}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="matlivre-web-${tag}"

(cd "${here}" && trunk build --release)
rm -rf "${here}/pkg/${name}"
mkdir -p "${here}/pkg/${name}"
cp -r "${here}/dist/." "${here}/pkg/${name}/"
rm -f "${here}/pkg/${name}/.trunk-lock"
tar --create --gzip --file "${here}/pkg/${name}.tar.gz" --directory "${here}/pkg" "${name}"
sha256sum "${here}/pkg/${name}.tar.gz"
