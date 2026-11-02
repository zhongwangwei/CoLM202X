#!/usr/bin/env bash
set -euo pipefail
export PATH=/home/zhwei/miniconda3/bin:$PATH
mode=${1:?MTON or MTOFF}
out=${OBREF_JOB_RESULTS:?runner results directory}
if [[ -n "${OBREF_INPUT_BLD:-}" ]]; then
  cp -a "$OBREF_INPUT_BLD" .bld
  cp -p "$OBREF_INPUT_HEADER" include/define.h
  printf '%s\n' 'Reuse prior private build; compile changed sources and relink.' > "$out/config.log"
else
  python3 tests/create_material_defineh.py GRID LULC_IGBP_PC URBANOFF vanGenu CaMaOFF BGCON CROPOFF "$mode" > "$out/config.log"
fi
cp include/define.h "$out/define.h"
make colm.x mkinidata.x mksrfdata.x FF=mpif90 NETCDF_INC="$(nf-config --includedir)" \
  NETCDF_LIB="$(nc-config --libdir)" \
  FOPTS='-fdefault-real-8 -fdefault-double-8 -ffree-form -O0 -g -cpp -ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium -fcheck=bounds -ffpe-trap=invalid,zero,overflow' \
  > "$out/build.log" 2>&1
test -x run/colm.x
test -x run/mkinidata.x
test -x run/mksrfdata.x
sha256sum run/colm.x run/mkinidata.x run/mksrfdata.x > "$out/binaries.sha256"
printf '%s\n' "PASS complete GRID+BGC $mode compile and link" | tee "$out/result.txt"
