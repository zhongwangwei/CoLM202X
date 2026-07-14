#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
for dependency in "$compiler" "$launcher" "$(command -v nf-config || true)" "$(command -v ncks || true)"; do
  if [[ -z "$dependency" ]]; then
    echo "ERROR: mpif90, mpiexec, nf-config, and ncks are required" >&2
    exit 1
  fi
done

build_dir=$(mktemp -d /tmp/colm-river-restart-mpi.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT

flags=(
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
includes=(-I"$build_dir" -J"$build_dir" "${netcdf_fflags[@]}")

printf '%s\n' '#define USEMPI' '#define GridRiverLakeFlow' > "$build_dir/define.h"

compile() {
  local source=$1
  local object=$2
  "$compiler" "${flags[@]}" "${includes[@]}" -c "$source" -o "$build_dir/$object"
}

compile share/MOD_Precision.F90 MOD_Precision.o
compile share/MOD_SPMD_Task.F90 MOD_SPMD_Task.o
compile tests/river_restart_mpi_test_support.F90 river_restart_mpi_test_support.o
compile share/MOD_NetCDFSerial.F90 MOD_NetCDFSerial.o
compile main/HYDRO/MOD_Vector_ReadWrite.F90 MOD_Vector_ReadWrite.o
compile main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90 MOD_Grid_RiverLakeTimeVars.o
compile tests/river_restart_mpi_harness.F90 river_restart_mpi_harness.o
compile tests/river_restart_netcdf_mutator.F90 river_restart_netcdf_mutator.o

common_objects=(
  "$build_dir/MOD_Precision.o"
  "$build_dir/MOD_SPMD_Task.o"
  "$build_dir/river_restart_mpi_test_support.o"
  "$build_dir/MOD_NetCDFSerial.o"
)
"$compiler" "${flags[@]}" "${common_objects[@]}" \
  "$build_dir/MOD_Vector_ReadWrite.o" "$build_dir/MOD_Grid_RiverLakeTimeVars.o" \
  "$build_dir/river_restart_mpi_harness.o" "${netcdf_libs[@]}" "${netcdf_c_libs[@]}" \
  -o "$build_dir/river_restart_mpi_harness"
"$compiler" "${flags[@]}" "${common_objects[@]}" \
  "$build_dir/river_restart_netcdf_mutator.o" "${netcdf_libs[@]}" "${netcdf_c_libs[@]}" \
  -o "$build_dir/river_restart_netcdf_mutator"

run_mpi() {
  local ranks=$1
  shift
  local command=(env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 "$launcher" -n "$ranks" "$@")
  if command -v timeout >/dev/null 2>&1; then
    timeout 60s "${command[@]}"
  else
    "${command[@]}"
  fi
}

expect_reject() {
  local label=$1
  local damaged=$2
  local expected=$3
  local log="$build_dir/${label}.log"
  set +e
  run_mpi 5 "$build_dir/river_restart_mpi_harness" read "$damaged" >"$log" 2>&1
  local status=$?
  set -e
  if [[ $status -eq 0 ]] || ! grep -q "$expected" "$log" || grep -q 'Floating point exception' "$log"; then
    echo "ERROR: $label was not rejected trap-safely" >&2
    sed -n '1,180p' "$log" >&2
    exit 1
  fi
}

restart="$build_dir/roundtrip.nc"
run_mpi 3 "$build_dir/river_restart_mpi_harness" write "$restart"
run_mpi 5 "$build_dir/river_restart_mpi_harness" read "$restart"

ncks -m "$restart" >"$build_dir/roundtrip-metadata.txt"
if ! grep -q 'gridriver_reservoir_identity' "$build_dir/roundtrip-metadata.txt"; then
  echo "ERROR: active-reservoir schema-v2 restart lacks reservoir identity" >&2
  exit 1
fi
if grep -q 'wdsrf_ucat_prev' "$build_dir/roundtrip-metadata.txt"; then
  echo "ERROR: BIF-disabled schema-v2 restart unexpectedly contains wdsrf_ucat_prev" >&2
  exit 1
fi

for field in wdsrf_ucat veloc_riv acctime_rnof acc_rnof_uc volwater_ucat volresv; do
  missing="$build_dir/missing-${field}.nc"
  ncks -O -x -v "$field" "$restart" "$missing"
  expect_reject "missing-${field}" "$missing" 'schema-v2 restart is missing required base state'
done

missing_reservoir_identity="$build_dir/missing-reservoir-identity.nc"
ncks -O -x -v gridriver_reservoir_identity "$restart" "$missing_reservoir_identity"
expect_reject "missing-reservoir-identity" "$missing_reservoir_identity" \
  'schema-v2 restart is missing required reservoir identity'

for mode in reservoir-identity-swap reservoir-identity-nan; do
  damaged="$build_dir/${mode}.nc"
  cp "$restart" "$damaged"
  "$build_dir/river_restart_netcdf_mutator" "$mode" "$damaged"
  expect_reject "$mode" "$damaged" 'different reservoir configuration'
done

for mode in wdsrf-nan wdsrf-negative velocity-nan velocity-out-of-range \
            acctime-nan acctime-negative acc-runoff-nan volwater-nan volwater-negative \
            volresv-nan volresv-negative; do
  damaged="$build_dir/${mode}.nc"
  cp "$restart" "$damaged"
  "$build_dir/river_restart_netcdf_mutator" "$mode" "$damaged"
  expect_reject "$mode" "$damaged" 'invalid GridRiverLake restart base state'
done

for mode in manifest-bif-invalid manifest-levee-invalid; do
  damaged="$build_dir/${mode}.nc"
  cp "$restart" "$damaged"
  "$build_dir/river_restart_netcdf_mutator" "$mode" "$damaged"
  expect_reject "$mode" "$damaged" 'invalid GridRiverLake restart feature manifest'
done

missing_manifest="$build_dir/transaction-missing-feature-manifest.nc"
ncks -O -x -v gridriver_restart_feature_bifurcation "$restart" "$missing_manifest"
set +e
run_mpi 5 "$build_dir/river_restart_mpi_harness" read "$missing_manifest" >"$build_dir/missing-manifest.log" 2>&1
missing_manifest_status=$?
set -e
if [[ $missing_manifest_status -eq 0 ]] || \
   ! grep -q 'GridRiverLake schema-v2 restart is missing feature manifest' \
     "$build_dir/missing-manifest.log"; then
  echo "ERROR: schema-v2 restart missing a feature manifest variable was not rejected" >&2
  sed -n '1,160p' "$build_dir/missing-manifest.log" >&2
  exit 1
fi

legacy="$build_dir/legacy-no-metadata-or-prev.nc"
ncks -O -x -v gridriver_restart_schema,gridriver_restart_complete,gridriver_ucatch_identity,\
gridriver_reservoir_identity,gridriver_restart_feature_bifurcation,gridriver_restart_feature_levee \
  "$restart" "$legacy"
run_mpi 5 "$build_dir/river_restart_mpi_harness" read-legacy-prev "$legacy"

tampered="$build_dir/tampered.nc"
cp "$restart" "$tampered"
"$build_dir/river_restart_netcdf_mutator" identity "$tampered"
set +e
run_mpi 5 "$build_dir/river_restart_mpi_harness" read "$tampered" >"$build_dir/tampered.log" 2>&1
tampered_status=$?
set -e
if [[ $tampered_status -eq 0 ]] || ! grep -q 'Refusing to load GridRiverLake state for a different ucatch network' "$build_dir/tampered.log"; then
  echo "ERROR: tampered identity was not rejected as expected" >&2
  sed -n '1,160p' "$build_dir/tampered.log" >&2
  exit 1
fi

nan_identity="$build_dir/nan-identity.nc"
cp "$restart" "$nan_identity"
"$build_dir/river_restart_netcdf_mutator" identity-nan "$nan_identity"
set +e
run_mpi 5 "$build_dir/river_restart_mpi_harness" read "$nan_identity" >"$build_dir/nan-identity.log" 2>&1
nan_identity_status=$?
set -e
if [[ $nan_identity_status -eq 0 ]] || ! grep -q 'Refusing to load GridRiverLake state for a different ucatch network' "$build_dir/nan-identity.log"; then
  echo "ERROR: non-finite ucatch identity was not rejected trap-safely" >&2
  sed -n '1,160p' "$build_dir/nan-identity.log" >&2
  exit 1
fi

incomplete="$build_dir/incomplete.nc"
run_mpi 3 "$build_dir/river_restart_mpi_harness" write-incomplete "$incomplete"
set +e
run_mpi 5 "$build_dir/river_restart_mpi_harness" read "$incomplete" >"$build_dir/incomplete.log" 2>&1
incomplete_status=$?
set -e
if [[ $incomplete_status -eq 0 ]] || ! grep -q 'GridRiverLake restart transaction is not complete' "$build_dir/incomplete.log"; then
  echo "ERROR: incomplete restart was not rejected as expected" >&2
  sed -n '1,160p' "$build_dir/incomplete.log" >&2
  exit 1
fi

empty="$build_dir/empty.nc"
run_mpi 3 "$build_dir/river_restart_mpi_harness" write-empty "$empty"
run_mpi 5 "$build_dir/river_restart_mpi_harness" read-empty "$empty"
ncks -m "$empty" >"$build_dir/empty-metadata.txt"
if grep -q 'volresv' "$build_dir/empty-metadata.txt" || \
   grep -q 'gridriver_reservoir_identity' "$build_dir/empty-metadata.txt" || \
   grep -q 'reservoir =' "$build_dir/empty-metadata.txt"; then
  echo "ERROR: zero-ucatch restart unexpectedly contains reservoir state" >&2
  exit 1
fi

echo "river restart MPI evaluator: PASS"
