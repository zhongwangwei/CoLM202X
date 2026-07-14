#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
make_cmd=${GMAKE:-$(command -v gmake || command -v make || true)}
if [[ -z "$compiler" || -z "$launcher" || -z "$make_cmd" ]]; then
  echo "ERROR: mpif90, mpiexec/mpirun, and GNU make are required" >&2
  exit 1
fi
if [[ ! -d .bld ]]; then
  echo "ERROR: .bld is required for the production-module evaluator" >&2
  exit 1
fi

moddir=$(mktemp -d /tmp/colm-river-bif-physics.XXXXXX)
trap 'rm -rf "$moddir"' EXIT

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
includes=(-Iinclude -I"$moddir" -I.bld -J"$moddir")

# Compile the modules under test from the current source tree.  Keeping the
# freshly generated module directory before .bld prevents a stale .mod file
# from hiding an interface change in WorkerPushData, Network, Levee, or BIF.
"$compiler" "${flags[@]}" "${includes[@]}" -c \
  share/MOD_WorkerPushData.F90 -o "$moddir/MOD_WorkerPushData.o"
"$compiler" "${flags[@]}" "${includes[@]}" -c \
  main/HYDRO/MOD_Grid_RiverLakeNetwork.F90 -o "$moddir/MOD_Grid_RiverLakeNetwork.o"
"$compiler" "${flags[@]}" "${includes[@]}" -c \
  main/HYDRO/MOD_Grid_RiverLakeLevee.F90 -o "$moddir/MOD_Grid_RiverLakeLevee.o"
"$compiler" "${flags[@]}" "${includes[@]}" -c \
  main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90 -o "$moddir/MOD_Grid_RiverLakeBifurcation.o"
"$compiler" "${flags[@]}" "${includes[@]}" -c \
  tests/river_bif_physics_harness.F90 -o "$moddir/river_bif_physics_harness.o"

shared_line=$("$make_cmd" --no-print-directory \
  --eval='print-shared:;@echo $(OBJS_SHARED_T)' print-shared)
read -r -a shared_objects <<< "$shared_line"
link_objects=()
for object in "${shared_objects[@]}"; do
  if [[ "$object" != .bld/MOD_WorkerPushData.o ]]; then
    link_objects+=("$object")
  fi
done

netcdf_flags=()
if command -v nf-config >/dev/null 2>&1; then
  read -r -a netcdf_flags <<< "$(nf-config --flibs)"
else
  netcdf_flags=(-lnetcdff -lnetcdf)
fi
# Homebrew's nf-config may omit the separate netcdf-c search directory.
if command -v brew >/dev/null 2>&1; then
  netcdf_flags=(-L"$(brew --prefix netcdf)/lib" "${netcdf_flags[@]}")
fi

"$compiler" -fopenmp -o "$moddir/river_bif_physics_harness" \
  "${link_objects[@]}" \
  .bld/MOD_Vector_ReadWrite.o .bld/MOD_Grid_Reservoir.o \
  "$moddir/MOD_WorkerPushData.o" \
  "$moddir/MOD_Grid_RiverLakeNetwork.o" \
  "$moddir/MOD_Grid_RiverLakeLevee.o" \
  "$moddir/MOD_Grid_RiverLakeBifurcation.o" \
  "$moddir/river_bif_physics_harness.o" \
  "${netcdf_flags[@]}" -llapack -lblas

for ranks in 1 2 4; do
  run=(env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 "$launcher" -n "$ranks" \
    "$moddir/river_bif_physics_harness")
  if command -v timeout >/dev/null 2>&1; then
    timeout 60s "${run[@]}"
  else
    "${run[@]}"
  fi
done

echo "river BIF physics evaluator: PASS"
