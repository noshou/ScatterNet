#!/usr/bin/env bash
# One-shot environment setup for a fresh clone: patched owl, this
# project's OCaml deps, and the Python side (numpy/xraydb).
set -euo pipefail
cd "$(dirname "$0")"

# 1. owl needs a local patch to build at all (see _patches/README.md) -
#    pin it before the deps-only install below, so that install doesn't
#    waste time building upstream's broken release first. Skipped if
#    already pinned.
if ! opam list owl 2>/dev/null | grep -q "pinned"; then
    echo "==> pinning patched owl"
    OWL_SRC="$(mktemp -d)"
    git clone --depth 1 --branch 1.2 https://github.com/owlbarn/owl.git "$OWL_SRC"
    patch -p1 -d "$OWL_SRC" < _patches/owl-1.2-exponpow-args.patch
    opam pin add owl "$OWL_SRC" -y
else
    echo "==> owl already pinned, skipping"
fi

# 2. everything else this project depends on (OCaml side)
echo "==> installing OCaml dependencies"
opam install . --deps-only -y

# 3. Python side - numpy/xraydb, used by scattering/form_factor_xraydb
echo "==> installing Python dependencies"
pip install -r scattering/form_factor_xraydb/requirements.txt

echo "==> done"
