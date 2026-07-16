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

build_dir=$(mktemp -d /tmp/colm-river-levee-restart.XXXXXX)
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

printf '%s\n' '#define USEMPI' '#define GridRiverLakeFlow' > "$build_dir/define.h"

compile() {
  local source=$1
  local object=$2
  "$compiler" "${flags[@]}" "${includes[@]}" -c "$source" -o "$build_dir/$object"
}

# Compile every river/levee restart module freshly.  This keeps the harness
# independent of stale .bld objects while retaining unrelated shared support.
compile share/MOD_NetCDFSerial.F90 MOD_NetCDFSerial.o
compile share/MOD_WorkerPushData.F90 MOD_WorkerPushData.o
compile main/HYDRO/MOD_Vector_ReadWrite.F90 MOD_Vector_ReadWrite.o
compile main/HYDRO/MOD_Grid_RiverLakeNetwork.F90 MOD_Grid_RiverLakeNetwork.o
compile main/HYDRO/MOD_Grid_Reservoir.F90 MOD_Grid_Reservoir.o
compile main/HYDRO/MOD_Grid_RiverLakeLevee.F90 MOD_Grid_RiverLakeLevee.o
compile main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90 MOD_Grid_RiverLakeBifurcation.o
compile main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90 MOD_Grid_RiverLakeTimeVars.o
compile tests/river_levee_restart_mpi_harness.F90 river_levee_restart_mpi_harness.o
compile tests/river_levee_restart_netcdf_mutator.F90 river_levee_restart_netcdf_mutator.o

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

"$compiler" -fopenmp -o "$build_dir/river_levee_restart_mpi_harness" \
  "${link_objects[@]}" \
  "${fresh_objects[@]}" "$build_dir/river_levee_restart_mpi_harness.o" \
  "${netcdf_libs[@]}" "${netcdf_c_libs[@]}" -llapack -lblas

"$compiler" -fopenmp -o "$build_dir/river_levee_restart_netcdf_mutator" \
  "${link_objects[@]}" "$build_dir/MOD_NetCDFSerial.o" \
  "$build_dir/river_levee_restart_netcdf_mutator.o" \
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

check_rejected() {
  local label=$1
  local damaged=$2
  local expected=$3
  local log="$build_dir/${label}.log"
  set +e
  run_mpi 5 "$build_dir/river_levee_restart_mpi_harness" read-enabled \
    "$damaged" "$build_dir/${label}-unused.nc" >"$log" 2>&1
  local status=$?
  set -e
  if [[ $status -eq 0 ]] || ! grep -q "$expected" "$log" || \
     grep -q 'Floating point exception' "$log"; then
    echo "ERROR: levee restart corruption $label was not rejected trap-safely" >&2
    sed -n '1,200p' "$log" >&2
    exit 1
  fi
}

restart="$build_dir/levee-enabled-3r.nc"
repartitioned="$build_dir/levee-enabled-5r.nc"
run_mpi 3 "$build_dir/river_levee_restart_mpi_harness" write-enabled "$restart"
run_mpi 5 "$build_dir/river_levee_restart_mpi_harness" read-enabled "$restart" "$repartitioned"

missing="$build_dir/missing-levsto.nc"
ncks -O -x -v levsto "$restart" "$missing"
check_rejected missing-levsto "$missing" \
  'GridRiverLake restart declares levee enabled but levsto is missing'

nan_levsto="$build_dir/nan-levsto.nc"
cp "$restart" "$nan_levsto"
"$build_dir/river_levee_restart_netcdf_mutator" nan "$nan_levsto"
check_rejected nan-levsto "$nan_levsto" \
  'levsto contains a negative or non-finite value'

negative_levsto="$build_dir/negative-levsto.nc"
cp "$restart" "$negative_levsto"
"$build_dir/river_levee_restart_netcdf_mutator" negative "$negative_levsto"
check_rejected negative-levsto "$negative_levsto" \
  'levsto contains a negative or non-finite value'

# A disabled manifest owns no levsto state.  Even if an old/stale variable is
# present, the deferred reader must leave the current levee pool cold at zero.
manifest_disabled="$build_dir/manifest-levee-disabled-with-stale-levsto.nc"
ncap2 -O -s 'gridriver_restart_feature_levee=0' "$restart" "$manifest_disabled"
run_mpi 5 "$build_dir/river_levee_restart_mpi_harness" \
  read-manifest-disabled "$manifest_disabled"

# An enabled restart read under a currently levee-disabled configuration must
# move protected storage into the persisted visible pool without losing water.
run_mpi 5 "$build_dir/river_levee_restart_mpi_harness" read-fold-disabled "$restart"

echo "river levee restart MPI evaluator: PASS"
