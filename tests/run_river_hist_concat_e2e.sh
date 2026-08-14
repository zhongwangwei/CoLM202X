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

# The directory name contains a space ON PURPOSE. The aggregator moves its
# .tmp onto the target with a shell command; unquoted, that split the path into
# separate arguments, the move failed, and the .complete success marker was
# written anyway. A tidy path would never have shown it.
wd=$(mktemp -d "/tmp/colm river concat e2e.XXXXXX")
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
NRESV=${RH_TOTALNUMRESV:-37}

write_segment () {   # $1 = ranks, $2 = groups, $3 = day, $4 = log
  ( cd "$wd" && RH_NGROUP=$2 RH_DAY=$3 RH_TOTALNUMUCAT=$NUCAT RH_NLON=$NLON RH_NLAT=$NLAT \
      RH_TOTALNPTHOUT=$NPTH RH_NPTHLEV=$NLEV RH_TOTALNUMRESV=$NRESV RH_OUT="$wd/e2e.nc" \
      env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
          OMPI_MCA_rmaps_base_oversubscribe=1 \
      "$launcher" -n "$1" "$wd/h" ) > "$4" 2>&1
  grep -q "RHSHARD PASS" "$4" || { echo "harness FAILED"; cat "$4"; exit 1; }
}

echo "== segment 1 (7 ranks, 3 IO groups; one group owns nothing)"
write_segment 7 3 1 "$wd/seg1.log"
echo "   $(ls "$wd"/e2e_seg*_shard*.nc | wc -l | tr -d ' ') shards"

# A restart inside the same history period, with a DIFFERENT IO-group count.
# The previous segment's shards stay valid; this writes its own complete set.
echo "== segment 2: restart, 5 ranks / 2 IO groups, next day"
write_segment 5 2 2 "$wd/seg2.log"
echo "   $(ls "$wd"/e2e_seg*_shard*.nc | wc -l | tr -d ' ') shards total"

echo "== aggregating"
"$launcher" -n 1 "$aggx" "$nml" "$wd/e2e.nc" > "$wd/agg.log" 2>&1 || {
  echo "aggregator FAILED"; tail -20 "$wd/agg.log"; exit 1; }
grep -E "shards, identity consistent|rebuilt" "$wd/agg.log" | sed 's/^/   /'
[[ -f "$wd/e2e.nc" && -f "$wd/e2e.nc.complete" ]] || {
  echo "missing output or .complete marker"; exit 1; }
# The marker must never outlive a failed promotion: if the target is absent the
# marker has to be absent too, or an operator is told to trust a file that is
# not there.
if [[ ! -f "$wd/e2e.nc" && -f "$wd/e2e.nc.complete" ]]; then
  echo "FAIL: .complete written without a target file"; exit 1
fi

echo "== verifying every id landed in the right place"
python3 - "$wd/e2e.nc" "$NUCAT" "$NLON" "$NLAT" "$NPTH" "$NLEV" "${RH_TOTALNUMRESV:-37}" <<'PY'
import sys
import numpy as np
from netCDF4 import Dataset

path, nuc, nlon, nlat, npth, nlev, nresv = sys.argv[1], *map(int, sys.argv[2:])
fail = 0
with Dataset(path) as ds:
    nt = len(ds.dimensions['time'])
    if nt != 2:
        print(f"   FAIL: expected 2 merged time records, got {nt}")
        fail += 1
    if ds.variables['f_ucat_shard'].dimensions != ('time', 'lat_ucat', 'lon_ucat'):
        print("   FAIL: unexpected dimension order", ds.variables['f_ucat_shard'].dimensions)
        fail += 1
    grids = [np.asarray(ds.variables['f_ucat_shard'][t]) for t in range(nt)]
    bifs  = [np.asarray(ds.variables['f_bifflw_lev'][t]) for t in range(nt)]

# Each segment encodes its day into the values, so a record taken from the
# wrong segment is visible rather than plausible.
for t in range(len(grids)):
    day = t + 1
    exp = np.full((nlat, nlon), -1e36)
    for gid in range(1, nuc + 1):
        exp[((gid - 1) // nlon) % nlat, (gid - 1) % nlon] = gid + 0.5 + 1000.0 * day
    bad = int((~np.isclose(grids[t], exp, rtol=0, atol=1e-9)).sum())
    print(f"   unit catchments t{t}: {bad} of {grids[t].size} cells wrong")
    fail += bad != 0

    expb = np.empty((npth, nlev))
    for gid in range(1, npth + 1):
        for l in range(1, nlev + 1):
            expb[gid - 1, l - 1] = gid + 0.001 * l + 1000.0 * day
    badb = int((~np.isclose(bifs[t], expb, rtol=0, atol=1e-9)).sum())
    print(f"   pathways       t{t}: {badb} of {bifs[t].size} cells wrong")
    fail += badb != 0

# Reservoir fields must SURVIVE aggregation. Excluding them from the
# unit-catchment sweep stopped the wrong reconstruction but silently dropped
# them, which is why this check exists at all.
with Dataset(path) as ds:
    for name in ("volresv", "qresv_in", "qresv_out"):
        if name not in ds.variables:
            print(f"   FAIL: {name} missing from the aggregate (silently dropped)")
            fail += 1
            continue
        for t in range(len(grids)):
            expr = np.arange(1, nresv + 1, dtype=np.float64) + 0.25 + 1000.0 * (t + 1)
            r = np.asarray(ds.variables[name][t], dtype=np.float64)
            if r.shape != expr.shape:
                print(f"   FAIL: {name} shape {r.shape} != {expr.shape}")
                fail += 1
                continue
            badr = int((~np.isclose(r, expr, rtol=0, atol=1e-9)).sum())
            print(f"   {name:<12} t{t}: {badr} of {r.size} values wrong")
            fail += badr != 0

# Static per-reservoir metadata. Time-varying reservoir fields surviving does
# not prove this one did: it takes a separate rebuild path (no time axis), and
# the harness writes 70000+gid so a mis-scattered value is visible, not just a
# missing variable.
with Dataset(path) as ds:
    if "resv_GRAND_ID" not in ds.variables:
        print("   FAIL: resv_GRAND_ID missing from the aggregate")
        fail += 1
    else:
        g = np.asarray(ds.variables["resv_GRAND_ID"][:], dtype=np.int64)
        expg = 70000 + np.arange(1, nresv + 1, dtype=np.int64)
        if g.shape != expg.shape:
            print(f"   FAIL: resv_GRAND_ID shape {g.shape} != {expg.shape}")
            fail += 1
        else:
            badg = int((g != expg).sum())
            print(f"   resv_GRAND_ID    : {badg} of {g.size} values wrong")
            fail += badg != 0

sys.exit(1 if fail else 0)
PY

echo "river history concatenate e2e: PASS"
