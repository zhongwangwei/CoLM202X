#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
[[ -n "$compiler" && -n "$launcher" ]] || {
  echo "ERROR: mpif90 and mpiexec/mpirun are required" >&2; exit 1; }

build_dir=$(mktemp -d /tmp/colm-river-reservoir-mpi.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT

flags=(-fopenmp -fdefault-real-8 -ffree-form -fcheck=all
       -ffpe-trap=invalid,zero,overflow -fbacktrace -cpp
       -ffree-line-length-0 -fallow-argument-mismatch -w)
read -r -a netcdf_fflags <<< "$(nf-config --fflags)"
read -r -a netcdf_libs <<< "$(nf-config --flibs)"
read -r -a netcdf_c_libs <<< "$(nc-config --libs)"
printf '%s\n' '#define USEMPI' '#define GridRiverLakeFlow' > "$build_dir/define.h"
includes=(-I"$build_dir" -J"$build_dir" "${netcdf_fflags[@]}")

compile() {
  "$compiler" "${flags[@]}" "${includes[@]}" -c "$1" -o "$build_dir/$2"
}
compile share/MOD_Precision.F90 MOD_Precision.o
compile share/MOD_SPMD_Task.F90 MOD_SPMD_Task.o
compile tests/river_reservoir_mpi_test_support.F90 support.o
compile share/MOD_NetCDFSerial.F90 MOD_NetCDFSerial.o
compile share/MOD_Utils.F90 MOD_Utils.o
compile main/HYDRO/MOD_Grid_Reservoir.F90 MOD_Grid_Reservoir.o
compile tests/river_reservoir_mpi_harness.F90 harness.o

"$compiler" "${flags[@]}" "$build_dir"/*.o "${netcdf_libs[@]}" \
  "${netcdf_c_libs[@]}" -o "$build_dir/harness"
env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 \
  "$launcher" -n 5 "$build_dir/harness" "$build_dir/reservoir.nc"
