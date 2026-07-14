#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

pytest -q

bash tests/run_river_mpi_evaluator.sh
bash tests/run_river_restart_mpi_evaluator.sh
bash tests/run_river_levee_restart_mpi_evaluator.sh
bash tests/run_river_bif_restart_mpi_evaluator.sh

compiler=${MPIFC:-$(command -v mpif90 || true)}
if [[ -z "$compiler" ]]; then
  echo "ERROR: mpif90 is required for the fresh-module compile" >&2
  exit 1
fi
if [[ ! -d .bld ]]; then
  echo "ERROR: .bld dependency modules are required for the fresh-module compile" >&2
  exit 1
fi

moddir=$(mktemp -d /tmp/colm-river-mpi-phase2.XXXXXX)
trap 'rm -rf "$moddir"' EXIT

flags=(
  -fopenmp
  -fdefault-real-8
  -ffree-form
  -cpp
  -ffree-line-length-0
  -fallow-argument-mismatch
  -w
)
includes=(-I"$moddir" -Iinclude -I.bld -J"$moddir")
if command -v nf-config >/dev/null 2>&1; then
  includes+=(-I"$(nf-config --includedir)")
fi

sources=(
  share/MOD_WorkerPushData.F90
  main/HYDRO/MOD_Grid_RiverLakeNetwork.F90
  main/HYDRO/MOD_Grid_Reservoir.F90
  main/HYDRO/MOD_Grid_RiverLakeLevee.F90
  main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90
  main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90
  main/HYDRO/MOD_Grid_RiverLakeHistState.F90
  main/TRACER/MOD_Tracer_RiverLake.F90
  main/HYDRO/MOD_Grid_RiverLakeFlow.F90
  main/MOD_Vars_TimeVariables.F90
)

for source in "${sources[@]}"; do
  object="$moddir/$(basename "${source%.*}").o"
  "$compiler" "${flags[@]}" "${includes[@]}" -c "$source" -o "$object"
done

echo "river fresh-module compile: PASS"
echo "river MPI phase-2 evaluator: PASS"
