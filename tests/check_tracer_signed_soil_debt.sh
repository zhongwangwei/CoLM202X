#!/usr/bin/env bash
set -euo pipefail
# Uses an existing verified full build; never rebuilds or edits that build.
BUILD=$(cd "${1:?usage: $0 BUILD OUTPUT [SOILWATER_OBJECT]}" && pwd)
OUT=${2:?output directory required}
SRC=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
objects=()
includes=(-I"$BUILD/.bld" -I"$BUILD/include")
if [[ -n ${3:-} ]]; then
  objects+=("$3")
  includes=(-I"$(dirname "$3")" "${includes[@]}")
fi
${FC:-mpif90} -cpp -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fcheck=all \
  "${includes[@]}" "$SRC/tracer_signed_soil_debt.f90" "${objects[@]}" "$BUILD/lib/libcolm.a" \
  ${NETCDF_LIBS:-$(nf-config --flibs)} ${EXTRA_LIBS:--llapack -lblas} -o "$OUT/tracer_signed_soil_debt"
failed=0
for mode in uniform nondebt atlimit partial reverse throughdown throughup gradientdown gradientup different overlimit enriched nan nan-coupling overflow phase diffusion reaction provider solubility lifecycle; do
  if "$OUT/tracer_signed_soil_debt" "$mode" >"$OUT/$mode.log" 2>&1; then
    status=0
  else
    status=$?
  fi
  marker=''
  case "$mode" in
    different) marker='incompatible inventory signs' ;;
    overlimit) marker='oversized carrier debt' ;;
    enriched) marker='concentration overshoot' ;;
    nan) marker='nonfinite inventory' ;;
    nan-coupling) marker='nonfinite coupling' ;;
    lifecycle) marker='full timestep; direct-kernel experiment only' ;;
    overflow) marker='ratio overflow' ;;
    phase) marker='nonclosed or nonliquid step' ;;
    diffusion) marker='snow or diffusion' ;;
    reaction|provider|solubility) marker='non-passive species' ;;
  esac
  if [[ -z $marker ]]; then
    if [[ $status == 0 ]] && grep -q "PASS $mode" "$OUT/$mode.log"; then
      echo "PASS $mode"
    else
      echo "FAIL $mode (exit $status)"; failed=1
    fi
  else
    # A generic crash or harness assertion is not an accepted safety rejection.
    if [[ $status != 0 ]] && grep -Fq "signed soil debt: rejected $marker" "$OUT/$mode.log"; then
      echo "PASS $mode (explicit safety rejection)"
    else
      echo "FAIL $mode (expected signed soil debt safety rejection, exit $status)"; failed=1
    fi
  fi
done
exit "$failed"
