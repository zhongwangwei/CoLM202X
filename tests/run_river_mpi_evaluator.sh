#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
if [[ -z "$compiler" ]]; then
  echo "ERROR: mpif90 is required for the MPI evaluator" >&2
  exit 1
fi

launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
if [[ -z "$launcher" ]]; then
  echo "ERROR: mpiexec or mpirun is required for the MPI evaluator" >&2
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

cat > "$moddir/define.h" <<'EOF'
#define USEMPI
EOF

"$compiler" "${flags[@]}" -DUSEMPI -I"$moddir" -J"$moddir" -c \
  tests/river_mpi_test_support.F90 -o "$moddir/river_mpi_test_support.o"
"$compiler" "${flags[@]}" -DUSEMPI -I"$moddir" -J"$moddir" -c \
  share/MOD_WorkerPushData.F90 -o "$moddir/MOD_WorkerPushData.o"
"$compiler" "${flags[@]}" -DUSEMPI -I"$moddir" -J"$moddir" -c \
  tests/river_mpi_harness.F90 -o "$moddir/river_mpi_harness.o"
"$compiler" "${flags[@]}" \
  "$moddir/river_mpi_test_support.o" "$moddir/MOD_WorkerPushData.o" \
  "$moddir/river_mpi_harness.o" -o "$moddir/river_mpi_harness"

for ranks in 2 4; do
  run=(env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 "$launcher" -n "$ranks" \
    "$moddir/river_mpi_harness")
  if command -v timeout >/dev/null 2>&1; then
    timeout 60s "${run[@]}"
  else
    "${run[@]}"
  fi
done

echo "river MPI evaluator: PASS"
