#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

pytest -q

compiler=${MPIFC:-$(command -v mpif90 || true)}
if [[ -z "$compiler" ]]; then
  echo "ERROR: mpif90 is required for the MPI evaluator" >&2
  exit 1
fi

moddir=$(mktemp -d /tmp/colm-river-mpi-evaluator.XXXXXX)
trap 'rm -rf "$moddir"' EXIT

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

"$compiler" "${flags[@]}" -fsyntax-only -Iinclude -I.bld -J "$moddir" \
  share/MOD_WorkerPushData.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/HYDRO/MOD_Grid_RiverLakeNetwork.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/HYDRO/MOD_Grid_RiverLakeLevee.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/HYDRO/MOD_Grid_RiverLakeHistState.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/TRACER/MOD_Tracer_RiverLake.F90
"$compiler" "${flags[@]}" -fsyntax-only -I"$moddir" -Iinclude -I.bld -J "$moddir" \
  main/HYDRO/MOD_Grid_RiverLakeFlow.F90

echo "river MPI evaluator: PASS"
