#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Step 2 acceptance driver for .omx/plans/river-history-sharded-output.md
#
# Runs the IO-group shard writers across several rank / IO-group counts,
# including layouts where a whole IO group owns no data, and fails if any
# global id is missing, duplicated, out of range, or paired with the wrong
# value.
#
# Requires a prior model build (lib/libcolm.a + .bld).
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
[[ -d "$blddir" && -f "$libcolm" ]] || {
  echo "ERROR: model build not found ($blddir, $libcolm); build first." >&2; exit 1; }

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
for _flag in "${nc_lib[@]}"; do
  [[ "$_flag" == -L* ]] && nc_rpath+=("-Wl,-rpath,${_flag#-L}")
done

workdir=$(mktemp -d /tmp/colm-river-hist-shard.XXXXXX)
trap 'rm -rf "$workdir"' EXIT

flags=(-fopenmp -fdefault-real-8 -ffree-form -cpp -ffree-line-length-0
       -fallow-argument-mismatch -fbacktrace -w)

echo "== building shard harness"
"$compiler" "${flags[@]}" -DUSEMPI -I"$blddir" "${nc_inc[@]}" -J"$workdir" \
  -c tests/river_hist_shard_harness.F90 -o "$workdir/shard.o"
"$compiler" "${flags[@]}" "$workdir/shard.o" "$libcolm" \
  "${nc_lib[@]}" "${nc_rpath[@]}" -llapack -lblas -o "$workdir/shard"

# ranks:groups -- the last case leaves an entire IO group without data
cases=${RH_SHARD_CASES:-"3:1 4:1 6:2 8:2 9:3"}

status=0
for case in $cases; do
  ranks=${case%%:*}
  groups=${case##*:}
  out="$workdir/case_${ranks}_${groups}.nc"
  echo "-- ranks=$ranks groups=$groups"
  if ! RH_NGROUP=$groups RH_OUT="$out" \
      env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
          OMPI_MCA_rmaps_base_oversubscribe=1 \
        "$launcher" -n "$ranks" "$workdir/shard" > "$workdir/log_${ranks}_${groups}.txt" 2>&1
  then
    echo "   FAILED (nonzero exit)"; cat "$workdir/log_${ranks}_${groups}.txt"; status=1; continue
  fi
  if grep -q "RHSHARD PASS" "$workdir/log_${ranks}_${groups}.txt"; then
    echo "   PASS  ($(ls "$workdir"/case_${ranks}_${groups}_seg*_shard*.nc 2>/dev/null | wc -l | tr -d ' ') shards)"
  else
    echo "   FAILED"; cat "$workdir/log_${ranks}_${groups}.txt"; status=1
  fi
done

[[ $status -eq 0 ]] && echo "river history shard evaluator: PASS"
exit $status
