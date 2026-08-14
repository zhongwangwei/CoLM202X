#!/usr/bin/env bash
# Numerical regression for the methane host-water partition (conservation).
# Source-string tests cannot see a conservation error; this runs the arithmetic.
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"
fc=${FC:-$(command -v gfortran || true)}
[[ -n "$fc" ]] || { echo "ERROR: gfortran required" >&2; exit 1; }
[[ -f .bld/MOD_Precision.o ]] || { echo "ERROR: build the model first" >&2; exit 1; }
wd=$(mktemp -d /tmp/colm-methane-host-water.XXXXXX); trap 'rm -rf "$wd"' EXIT
"$fc" -ffree-form -fdefault-real-8 -ffree-line-length-0 -I.bld -J"$wd" \
  -o "$wd/mhw" tests/methane_host_water_harness.F90 .bld/MOD_Precision.o
"$wd/mhw"
