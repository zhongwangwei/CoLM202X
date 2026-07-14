#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

compiler=${MPIFC:-$(command -v mpif90 || true)}
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
if [[ -z "$compiler" || -z "$launcher" ]]; then
  echo "ERROR: mpif90 and mpiexec/mpirun are required" >&2
  exit 1
fi

moddir=$(mktemp -d /tmp/colm-worker-push-benchmark.XXXXXX)
trap 'rm -rf "$moddir"' EXIT

flags=(
  -O2
  -fdefault-real-8
  -ffree-form
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
  tests/worker_push_mapping_benchmark.F90 -o "$moddir/worker_push_mapping_benchmark.o"
"$compiler" "${flags[@]}" \
  "$moddir/river_mpi_test_support.o" "$moddir/MOD_WorkerPushData.o" \
  "$moddir/worker_push_mapping_benchmark.o" -o "$moddir/worker_push_mapping_benchmark"

ranks_to_run=${WORKER_PUSH_BENCH_RANKS:-"2 4 8 16 32"}
if [[ ${WORKER_PUSH_BENCH_INCLUDE_64:-0} == 1 ]]; then
  ranks_to_run="$ranks_to_run 64"
fi

launcher_args=()
if "$launcher" --version 2>&1 | grep -qi 'Open MPI'; then
  launcher_args+=(--oversubscribe)
fi

for ranks in $ranks_to_run; do
  timeout_seconds=${WORKER_PUSH_BENCH_TIMEOUT_SECONDS:-300}
  run=(env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 "$launcher" "${launcher_args[@]}" -n "$ranks" \
    "$moddir/worker_push_mapping_benchmark")
  if command -v timeout >/dev/null 2>&1; then
    timeout "${timeout_seconds}s" "${run[@]}"
  else
    "${run[@]}"
  fi
done
