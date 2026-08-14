#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Step 4 end-to-end acceptance: shards -> aggregator -> single file.
#
# Builds a real shard set with the shard harness (which always includes a
# zero-length worker and, at >= 2 groups, an IO group that owns nothing), runs
# river_hist_concatenate.x on it, and checks that every unit-catchment value
# landed in the grid cell its global id maps to and that pathway columns follow
# pth_global_id rather than rank order.
#
# Requires a prior 'make all' and 'make postprocess.x'.
# ---------------------------------------------------------------------------
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
[[ -n "$compiler" ]] || { echo "ERROR: mpif90 required" >&2; exit 1; }
[[ -n "$launcher" ]] || { echo "ERROR: mpiexec/mpirun required" >&2; exit 1; }

blddir=${COLM_BLD_DIR:-$repo_root/.bld}
libcolm=${COLM_LIB:-$repo_root/lib/libcolm.a}
aggx=$repo_root/run/river_hist_concatenate.x
[[ -d "$blddir" && -f "$libcolm" ]] || { echo "ERROR: build the model first" >&2; exit 1; }
[[ -x "$aggx" ]] || { echo "ERROR: build 'make postprocess.x' first" >&2; exit 1; }

if [[ -n "${NETCDF_PREFIX:-}" ]]; then
  nc_inc=(-I"$NETCDF_PREFIX/include"); nc_lib=(-L"$NETCDF_PREFIX/lib" -lnetcdff -lnetcdf)
else
  nc_inc=(); nc_lib=()
  command -v nf-config >/dev/null 2>&1 && {
    read -r -a _f <<< "$(nf-config --fflags)"; nc_inc+=("${_f[@]}")
    read -r -a _l <<< "$(nf-config --flibs)";  nc_lib+=("${_l[@]}"); }
  command -v nc-config >/dev/null 2>&1 && {
    read -r -a _c <<< "$(nc-config --libs)";   nc_lib+=("${_c[@]}"); }
fi
nc_rpath=()
for f in "${nc_lib[@]}"; do [[ "$f" == -L* ]] && nc_rpath+=("-Wl,-rpath,${f#-L}"); done

wd=$(mktemp -d /tmp/colm-river-concat-e2e.XXXXXX)
trap 'rm -rf "$wd"' EXIT

# A namelist the aggregator can actually read on this machine: the shipped ones
# point their forcing namelist at an HPC path.
nml=${RH_E2E_NAMELIST:-}
if [[ -z "$nml" ]]; then
  for cand in run/Amazon.nml run/*.nml; do [[ -f "$cand" ]] && { nml=$cand; break; }; done
  sed "s#DEF_forcing_namelist = .*#DEF_forcing_namelist = '$repo_root/run/forcing/JRA3Q.nml'#" \
    "$nml" > "$wd/test.nml"
  nml=$wd/test.nml
fi

flags=(-fopenmp -fdefault-real-8 -ffree-form -cpp -ffree-line-length-0
       -fallow-argument-mismatch -fbacktrace -w)

echo "== building shard harness"
"$compiler" "${flags[@]}" -DUSEMPI -I"$blddir" "${nc_inc[@]}" -J"$wd" \
  -c tests/river_hist_shard_harness.F90 -o "$wd/h.o"
"$compiler" "${flags[@]}" "$wd/h.o" "$libcolm" "${nc_lib[@]}" "${nc_rpath[@]}" \
  -llapack -lblas -o "$wd/h"

NUCAT=${RH_TOTALNUMUCAT:-1200}; NLON=${RH_NLON:-40}; NLAT=${RH_NLAT:-30}
NPTH=${RH_TOTALNPTHOUT:-131}; NLEV=${RH_NPTHLEV:-3}

echo "== writing shards (7 ranks, 3 IO groups; one group owns nothing)"
( cd "$wd" && RH_NGROUP=3 RH_TOTALNUMUCAT=$NUCAT RH_NLON=$NLON RH_NLAT=$NLAT \
    RH_TOTALNPTHOUT=$NPTH RH_NPTHLEV=$NLEV RH_OUT="$wd/e2e.nc" \
    env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        OMPI_MCA_rmaps_base_oversubscribe=1 \
    "$launcher" -n 7 "$wd/h" ) > "$wd/harness.log" 2>&1
grep -q "RHSHARD PASS" "$wd/harness.log" || { echo "harness FAILED"; cat "$wd/harness.log"; exit 1; }
echo "   $(ls "$wd"/e2e_shard*.nc | wc -l | tr -d ' ') shards"

echo "== aggregating"
"$launcher" -n 1 "$aggx" "$nml" "$wd/e2e.nc" > "$wd/agg.log" 2>&1 || {
  echo "aggregator FAILED"; tail -20 "$wd/agg.log"; exit 1; }
grep -E "shards, identity consistent|rebuilt" "$wd/agg.log" | sed 's/^/   /'
[[ -f "$wd/e2e.nc" && -f "$wd/e2e.nc.complete" ]] || {
  echo "missing output or .complete marker"; exit 1; }

echo "== verifying every id landed in the right place"
python3 - "$wd/e2e.nc" "$NUCAT" "$NLON" "$NLAT" "$NPTH" "$NLEV" <<'PY'
import sys
import numpy as np
from netCDF4 import Dataset

path, nuc, nlon, nlat, npth, nlev = sys.argv[1], *map(int, sys.argv[2:])
fail = 0
with Dataset(path) as ds:
    g = np.asarray(ds.variables['f_ucat_shard'][0])
    if ds.variables['f_ucat_shard'].dimensions != ('time', 'lat_ucat', 'lon_ucat'):
        print("   FAIL: unexpected dimension order", ds.variables['f_ucat_shard'].dimensions)
        fail += 1
    b = np.asarray(ds.variables['f_bifflw_lev'][0])

exp = np.full((nlat, nlon), -1e36)
for gid in range(1, nuc + 1):
    exp[((gid - 1) // nlon) % nlat, (gid - 1) % nlon] = gid + 0.5
bad = int((~np.isclose(g, exp, rtol=0, atol=1e-9)).sum())
print(f"   unit catchments: {bad} of {g.size} cells wrong")
fail += bad != 0

expb = np.empty((npth, nlev))
for gid in range(1, npth + 1):
    for l in range(1, nlev + 1):
        expb[gid - 1, l - 1] = gid + 0.001 * l
badb = int((~np.isclose(b, expb, rtol=0, atol=1e-9)).sum())
print(f"   pathways       : {badb} of {b.size} cells wrong")
fail += badb != 0

sys.exit(1 if fail else 0)
PY

echo "river history concatenate e2e: PASS"
