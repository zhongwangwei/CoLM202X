#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Step 5 negative acceptance: the aggregator must REFUSE bad shard sets.
#
# Plan step 5 items 7 and 14: missing / duplicated / out-of-range global ids,
# and schema-version / period-key / grid-fingerprint / shard-count conflicts,
# must fail rather than produce a plausible-looking file. A guard that is only
# present in the source is not a guard; each one is provoked here.
#
# Requires a prior 'make all' and 'make postprocess.x'.
# ---------------------------------------------------------------------------
set -uo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
aggx=$repo_root/run/river_hist_concatenate.x
blddir=${COLM_BLD_DIR:-$repo_root/.bld}
libcolm=${COLM_LIB:-$repo_root/lib/libcolm.a}
[[ -n "$compiler" && -n "$launcher" ]] || { echo "ERROR: MPI toolchain required" >&2; exit 1; }
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
nc_rpath=(); for f in "${nc_lib[@]}"; do [[ "$f" == -L* ]] && nc_rpath+=("-Wl,-rpath,${f#-L}"); done

wd=$(mktemp -d /tmp/colm-river-concat-neg.XXXXXX)
trap 'rm -rf "$wd"' EXIT

nml=$wd/test.nml
sed "s#DEF_forcing_namelist = .*#DEF_forcing_namelist = '$repo_root/run/forcing/JRA3Q.nml'#" \
  run/Amazon.nml > "$nml"

flags=(-fopenmp -fdefault-real-8 -ffree-form -cpp -ffree-line-length-0
       -fallow-argument-mismatch -fbacktrace -w)
"$compiler" "${flags[@]}" -DUSEMPI -I"$blddir" "${nc_inc[@]}" -J"$wd" \
  -c tests/river_hist_shard_harness.F90 -o "$wd/h.o" >/dev/null
"$compiler" "${flags[@]}" "$wd/h.o" "$libcolm" "${nc_lib[@]}" "${nc_rpath[@]}" \
  -llapack -lblas -o "$wd/h" >/dev/null

# 8 ranks / 2 groups: master + 2 IO + 5 workers, worker 0 deliberately empty,
# which still leaves both groups holding data -- needed so the duplicate-id
# case has two non-empty shards to work with.
make_shards () {   # $1 = destination stem
  ( cd "$wd" && RH_NGROUP=2 RH_TOTALNUMUCAT=600 RH_NLON=30 RH_NLAT=20 \
      RH_TOTALNPTHOUT=61 RH_OUT="$1" \
      env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
          OMPI_MCA_rmaps_base_oversubscribe=1 \
      "$launcher" -n 8 "$wd/h" ) >/dev/null 2>&1
}

# Runs the aggregator and requires it to FAIL.
expect_fail () {   # $1 = case name, $2 = target
  local out rc
  out=$("$launcher" -n 1 "$aggx" "$nml" "$2" 2>&1); rc=$?
  if [[ $rc -eq 0 ]]; then
    echo "   FAIL: $1 -- aggregator accepted a bad shard set"
    return 1
  fi
  echo "   ok: $1"
  return 0
}

status=0

echo "== sanity: an untouched shard set must aggregate"
make_shards "$wd/good.nc"
if "$launcher" -n 1 "$aggx" "$nml" "$wd/good.nc" >/dev/null 2>&1; then
  echo "   ok: clean set accepted"
else
  echo "   FAIL: clean set was rejected"; status=1
fi

echo "== schema version mismatch"
make_shards "$wd/ver.nc"
python3 - $(ls "$wd"/ver_seg*_shard00001.nc | head -1) <<'PY'
import sys
from netCDF4 import Dataset
with Dataset(sys.argv[1], 'a') as ds:
    ds.river_hist_shard_schema_version = 99
PY
expect_fail "schema version" "$wd/ver.nc" || status=1

echo "== history period mismatch"
make_shards "$wd/per.nc"
python3 - $(ls "$wd"/per_seg*_shard00001.nc | head -1) <<'PY'
import sys
from netCDF4 import Dataset
with Dataset(sys.argv[1], 'a') as ds:
    ds.history_period_key = '1999-12'
PY
expect_fail "period key" "$wd/per.nc" || status=1

echo "== grid fingerprint mismatch"
make_shards "$wd/grid.nc"
python3 - $(ls "$wd"/grid_seg*_shard00001.nc | head -1) <<'PY'
import sys
from netCDF4 import Dataset
with Dataset(sys.argv[1], 'a') as ds:
    ds.grid_fingerprint = 'a-different-grid'
PY
expect_fail "grid fingerprint" "$wd/grid.nc" || status=1

echo "== incomplete shard set (one shard removed)"
make_shards "$wd/miss.nc"
rm -f "$wd"/miss_seg*_shard00001.nc
expect_fail "missing shard" "$wd/miss.nc" || status=1

echo "== duplicated global id across shards"
make_shards "$wd/dup.nc"
python3 - "$wd" dup <<'PY'
import glob, sys
import numpy as np
from netCDF4 import Dataset
# An IO group can legitimately own nothing, so pick two shards that actually
# hold ids -- mutating an empty one would test nothing.
files = sorted(glob.glob(f"{sys.argv[1]}/{sys.argv[2]}_seg*_shard*.nc"))
nonempty = []
for f in files:
    with Dataset(f) as ds:
        if np.asarray(ds.variables['ucat_ucid'][:]).size:
            nonempty.append(f)
assert len(nonempty) >= 2, f"need two non-empty shards, got {len(nonempty)}"
with Dataset(nonempty[1]) as b:
    other = int(np.asarray(b.variables['ucat_ucid'][:])[0])
with Dataset(nonempty[0], 'a') as a:
    ids = np.asarray(a.variables['ucat_ucid'][:])
    ids[0] = other
    a.variables['ucat_ucid'][:] = ids
print(f"  (duplicated id {other} into {nonempty[0].split('/')[-1]})")
PY
expect_fail "duplicate id" "$wd/dup.nc" || status=1

echo "== out-of-range global id"
make_shards "$wd/range.nc"
python3 - "$wd" range <<'PY'
import glob, sys
import numpy as np
from netCDF4 import Dataset
for f in sorted(glob.glob(f"{sys.argv[1]}/{sys.argv[2]}_seg*_shard*.nc")):
    with Dataset(f, 'a') as ds:
        ids = np.asarray(ds.variables['ucat_ucid'][:])
        if ids.size:
            ids[0] = 999999
            ds.variables['ucat_ucid'][:] = ids
            break
PY
expect_fail "id out of range" "$wd/range.nc" || status=1

echo
if [[ $status -eq 0 ]]; then
  echo "river history concatenate negative tests: PASS"
else
  echo "river history concatenate negative tests: FAILED" >&2
fi
exit $status
