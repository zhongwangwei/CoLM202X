#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
make_cmd=${GMAKE:-$(command -v gmake || command -v make || true)}
for dependency in "$compiler" "$launcher" "$make_cmd" "$(command -v nf-config || true)" \
  "$(command -v ncks || true)" "$(command -v ncap2 || true)"; do
  if [[ -z "$dependency" ]]; then
    echo "ERROR: mpif90, mpiexec, GNU make, nf-config, ncks, and ncap2 are required" >&2
    exit 1
  fi
done
if [[ ! -d .bld ]]; then
  echo "ERROR: .bld dependency modules are required" >&2
  exit 1
fi

build_dir=$(mktemp -d /tmp/colm-river-bif-restart.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT

flags=(
  -O0
  -g
  -fopenmp
  -fdefault-real-8
  -ffree-form
  -fcheck=all
  -ffpe-trap=invalid,zero,overflow
  -fbacktrace
  -cpp
  -ffree-line-length-0
  -fallow-argument-mismatch
  -w
)
read -r -a netcdf_fflags <<< "$(nf-config --fflags)"
read -r -a netcdf_libs <<< "$(nf-config --flibs)"
if command -v pkg-config >/dev/null 2>&1; then
  read -r -a netcdf_c_libs <<< "$(pkg-config --libs netcdf)"
elif command -v nc-config >/dev/null 2>&1; then
  read -r -a netcdf_c_libs <<< "$(nc-config --libs)"
else
  netcdf_c_libs=()
fi
includes=(-I"$build_dir" -Iinclude -I.bld -J"$build_dir" "${netcdf_fflags[@]}")

# The evaluator exercises the river/BIF restart unit only. Keep optional
# tracer code out so the link surface cannot silently depend on stale tracer
# objects from .bld.
printf '%s\n' '#define USEMPI' '#define GridRiverLakeFlow' > "$build_dir/define.h"

compile() {
  local source=$1
  local object=$2
  "$compiler" "${flags[@]}" "${includes[@]}" -c "$source" -o "$build_dir/$object"
}

# These are the production modules under test. Compile every one freshly so a
# stale .bld object or .mod cannot hide an interface or restart-format change.
compile share/MOD_NetCDFSerial.F90 MOD_NetCDFSerial.o
compile share/MOD_WorkerPushData.F90 MOD_WorkerPushData.o
compile main/HYDRO/MOD_Vector_ReadWrite.F90 MOD_Vector_ReadWrite.o
compile main/HYDRO/MOD_Grid_RiverLakeNetwork.F90 MOD_Grid_RiverLakeNetwork.o
compile main/HYDRO/MOD_Grid_Reservoir.F90 MOD_Grid_Reservoir.o
compile main/HYDRO/MOD_Grid_RiverLakeLevee.F90 MOD_Grid_RiverLakeLevee.o
compile main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90 MOD_Grid_RiverLakeBifurcation.o
compile main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90 MOD_Grid_RiverLakeTimeVars.o
compile tests/river_bif_restart_mpi_harness.F90 river_bif_restart_mpi_harness.o
compile tests/river_restart_netcdf_mutator.F90 river_restart_netcdf_mutator.o

shared_line=$("$make_cmd" --no-print-directory \
  --eval='print-shared:;@echo $(OBJS_SHARED_T)' print-shared)
read -r -a shared_objects <<< "$shared_line"
link_objects=()
for object in "${shared_objects[@]}"; do
  case "$object" in
    .bld/MOD_NetCDFSerial.o|.bld/MOD_WorkerPushData.o) ;;
    *) link_objects+=("$object") ;;
  esac
done

fresh_objects=(
  "$build_dir/MOD_NetCDFSerial.o"
  "$build_dir/MOD_WorkerPushData.o"
  "$build_dir/MOD_Vector_ReadWrite.o"
  "$build_dir/MOD_Grid_RiverLakeNetwork.o"
  "$build_dir/MOD_Grid_Reservoir.o"
  "$build_dir/MOD_Grid_RiverLakeLevee.o"
  "$build_dir/MOD_Grid_RiverLakeBifurcation.o"
  "$build_dir/MOD_Grid_RiverLakeTimeVars.o"
)

"$compiler" -fopenmp -o "$build_dir/river_bif_restart_mpi_harness" \
  "${link_objects[@]}" \
  "${fresh_objects[@]}" "$build_dir/river_bif_restart_mpi_harness.o" \
  "${netcdf_libs[@]}" "${netcdf_c_libs[@]}" -llapack -lblas

"$compiler" -fopenmp -o "$build_dir/river_restart_netcdf_mutator" \
  "${link_objects[@]}" "$build_dir/MOD_NetCDFSerial.o" \
  "$build_dir/river_restart_netcdf_mutator.o" \
  "${netcdf_libs[@]}" "${netcdf_c_libs[@]}" -llapack -lblas

run_mpi() {
  local ranks=$1
  shift
  local command=(env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 "$launcher" -n "$ranks" "$@")
  if command -v timeout >/dev/null 2>&1; then
    timeout 90s "${command[@]}"
  else
    "${command[@]}"
  fi
}

restart="$build_dir/bif-write-3r.nc"
repartitioned="$build_dir/bif-readwrite-5r.nc"
run_mpi 3 "$build_dir/river_bif_restart_mpi_harness" write "$restart"
"$build_dir/river_restart_netcdf_mutator" check-bif-nonzero "$restart"
run_mpi 5 "$build_dir/river_bif_restart_mpi_harness" read-valid "$restart" "$repartitioned"
"$build_dir/river_restart_netcdf_mutator" compare-bif-state "$restart" "$repartitioned"

check_cold_start() {
  local label=$1
  local damaged=$2
  local rewritten="$build_dir/${label}-cold.nc"
  run_mpi 5 "$build_dir/river_bif_restart_mpi_harness" read-cold "$damaged" "$rewritten"
  "$build_dir/river_restart_netcdf_mutator" check-bif-cold "$rewritten"
}

check_declared_bif_corruption() {
  local label=$1
  local damaged=$2
  local expected=${3:-'GridRiverLake restart declares bifurcation enabled but pathway state is incomplete'}
  local log="$build_dir/${label}.log"
  set +e
  run_mpi 5 "$build_dir/river_bif_restart_mpi_harness" read-valid \
    "$damaged" "$build_dir/${label}-unused.nc" >"$log" 2>&1
  local status=$?
  set -e
  if [[ $status -eq 0 ]] || ! grep -q "$expected" "$log" || \
     grep -q 'Floating point exception' "$log"; then
    echo "ERROR: declared BIF corruption $label was not rejected" >&2
    sed -n '1,200p' "$log" >&2
    exit 1
  fi
}

missing_previous="$build_dir/missing-wdsrf-ucat-prev.nc"
ncks -O -x -v wdsrf_ucat_prev "$restart" "$missing_previous"
check_declared_bif_corruption missing-wdsrf-ucat-prev "$missing_previous" \
  'GridRiverLake restart transaction is missing wdsrf_ucat_prev'

nan_previous="$build_dir/nan-wdsrf-ucat-prev.nc"
cp "$restart" "$nan_previous"
"$build_dir/river_restart_netcdf_mutator" previous-depth-nan "$nan_previous"
check_declared_bif_corruption nan-wdsrf-ucat-prev "$nan_previous" \
  'GridRiverLake restart has invalid wdsrf_ucat_prev'

negative_previous="$build_dir/negative-wdsrf-ucat-prev.nc"
cp "$restart" "$negative_previous"
"$build_dir/river_restart_netcdf_mutator" previous-depth-negative "$negative_previous"
check_declared_bif_corruption negative-wdsrf-ucat-prev "$negative_previous" \
  'GridRiverLake restart has invalid wdsrf_ucat_prev'

legacy_missing_previous="$build_dir/legacy-missing-wdsrf-ucat-prev.nc"
ncks -O -x -v gridriver_restart_schema,gridriver_restart_complete,gridriver_ucatch_identity,\
gridriver_restart_feature_bifurcation,gridriver_restart_feature_levee,wdsrf_ucat_prev \
  "$restart" "$legacy_missing_previous"
check_cold_start legacy-missing-wdsrf-ucat-prev "$legacy_missing_previous"

missing_momentum="$build_dir/missing-pth-momen.nc"
ncks -O -x -v pth_momen "$restart" "$missing_momentum"
check_declared_bif_corruption missing-pth-momen "$missing_momentum"

# Schema-v1 predates feature declarations.  Missing optional BIF state must
# remain backward-compatible by atomically cold-starting momentum and depth.
schema1_momentum="$build_dir/schema1-missing-pth-momen-base.nc"
missing_momentum_v1="$build_dir/schema1-missing-pth-momen.nc"
ncap2 -O -s 'gridriver_restart_schema=1' "$restart" "$schema1_momentum"
ncks -O -x -v gridriver_restart_feature_bifurcation,gridriver_restart_feature_levee,pth_momen \
  "$schema1_momentum" "$missing_momentum_v1"
check_cold_start schema1-missing-pth-momen "$missing_momentum_v1"

missing_signature="$build_dir/missing-bif-signature.nc"
ncks -O -x -v bif_path_signature "$restart" "$missing_signature"
check_declared_bif_corruption missing-bif-signature "$missing_signature"

nan_momentum="$build_dir/nan-pth-momen.nc"
cp "$restart" "$nan_momentum"
"$build_dir/river_restart_netcdf_mutator" bif-momentum-nan "$nan_momentum"
check_declared_bif_corruption nan-pth-momen "$nan_momentum" \
  'GridRiverLake restart declares bifurcation enabled but pathway state is invalid'

schema1_nan="$build_dir/schema1-nan-pth-momen.nc"
ncap2 -O -s 'gridriver_restart_schema=1' "$restart" "$schema1_nan"
ncks -O -x -v gridriver_restart_feature_bifurcation,gridriver_restart_feature_levee \
  "$schema1_nan" "$build_dir/schema1-nan-pth-momen-no-manifest.nc"
schema1_nan="$build_dir/schema1-nan-pth-momen-no-manifest.nc"
"$build_dir/river_restart_netcdf_mutator" bif-momentum-nan "$schema1_nan"
check_cold_start schema1-nan-pth-momen "$schema1_nan"

nan_signature="$build_dir/nan-bif-signature.nc"
cp "$restart" "$nan_signature"
"$build_dir/river_restart_netcdf_mutator" bif-signature-nan "$nan_signature"
set +e
run_mpi 5 "$build_dir/river_bif_restart_mpi_harness" read-valid \
  "$nan_signature" "$build_dir/unused.nc" >"$build_dir/nan-signature.log" 2>&1
nan_signature_status=$?
set -e
if [[ $nan_signature_status -eq 0 ]] || \
   ! grep -q 'Refusing to load bifurcation momentum for a different pathway network' \
     "$build_dir/nan-signature.log" || \
   grep -q 'Floating point exception' "$build_dir/nan-signature.log"; then
  echo "ERROR: non-finite BIF signature was not rejected trap-safely" >&2
  sed -n '1,200p' "$build_dir/nan-signature.log" >&2
  exit 1
fi

echo "river BIF restart MPI evaluator: PASS"
