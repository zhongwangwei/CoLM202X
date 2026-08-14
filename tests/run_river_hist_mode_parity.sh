#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Step 5 items 8-12: one vs block+postprocess parity on a REAL case.
#
#   8.  DEF_HIST_mode='one' compatibility regression
#   9.  block shards + postprocess compared with the 'one' file, variable by
#       variable
#   10. restart inside a history window, same MPI layout
#   11. restart inside a history window, DIFFERENT worker/IO-group layout
#   12. the variable sets, dimensions and attributes are enumerated from the
#       files themselves, never from a hand-maintained expected list
#
# These need land data and forcing, so they cannot run on a laptop -- this is
# the script to run on the target machine. Everything it needs is a working
# case: the same namelist, run twice (plus twice more for the restart cases).
#
# Usage:
#   tests/run_river_hist_mode_parity.sh <namelist> <workdir> [ranks] [ranks_alt]
#
# The namelist is copied and edited per phase; the original is not modified.
# ---------------------------------------------------------------------------
set -uo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

if [[ $# -lt 2 ]]; then
  sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  exit 2
fi

namelist=$1
workdir=$2
ranks=${3:-8}
ranks_alt=${4:-6}

colm=$repo_root/run/colm.x
wrapper=$repo_root/run/scripts/concatenate_history
launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}

for f in "$colm" "$wrapper" "$repo_root/run/river_hist_concatenate.x"; do
  [[ -x "$f" ]] || { echo "ERROR: $f missing; build 'make all' and 'make postprocess.x'" >&2; exit 1; }
done
[[ -f "$namelist" ]] || { echo "ERROR: no such namelist: $namelist" >&2; exit 1; }
mkdir -p "$workdir"

# --- helpers ---------------------------------------------------------------
mk_namelist () {   # $1 = out, $2 = mode, $3 = history dir, $4 = restart dir
  sed -e "s#DEF_HIST_mode *=.*#DEF_HIST_mode = '$2'#" \
      -e "s#DEF_dir_history *=.*#DEF_dir_history = '$3'#" \
      -e "s#DEF_dir_restart *=.*#DEF_dir_restart = '$4'#" \
      "$namelist" > "$1"
  grep -q "DEF_HIST_mode" "$1" || echo "   WARNING: DEF_HIST_mode not present in namelist" >&2
}

run_case () {      # $1 = namelist, $2 = ranks, $3 = label
  echo "-- running $3 on $2 ranks"
  if ! "$launcher" -n "$2" "$colm" "$1" > "$workdir/$3.log" 2>&1; then
    echo "   FAILED (see $workdir/$3.log)" >&2
    return 1
  fi
}

status=0

# --- phase 1: reference 'one' run -----------------------------------------
mk_namelist "$workdir/one.nml"   one   "$workdir/hist_one"   "$workdir/rest_one"
run_case "$workdir/one.nml" "$ranks" one || exit 1

# --- phase 2: block run + aggregation --------------------------------------
mk_namelist "$workdir/block.nml" block "$workdir/hist_block" "$workdir/rest_block"
run_case "$workdir/block.nml" "$ranks" block || exit 1

echo "-- checking the pending manifest exists before aggregation"
if ! ls "$workdir"/hist_block/*.pending >/dev/null 2>&1; then
  echo "   FAIL: block run left no .pending note" >&2; status=1
else
  echo "   ok"
fi

echo "-- aggregating"
"$wrapper" "$workdir/block.nml" "$workdir/hist_block" || { echo "   aggregation FAILED" >&2; exit 1; }

echo "-- checking the pending manifest was cleared"
if ls "$workdir"/hist_block/*.pending >/dev/null 2>&1; then
  echo "   FAIL: .pending survived a successful aggregation" >&2; status=1
else
  echo "   ok"
fi

# --- phase 3: restart mid-window, same layout ------------------------------
# (Driven by the namelist's own restart settings; the caller is expected to
#  point DEF_simulation_time at a window that straddles a history boundary.)
if [[ "${RH_PARITY_RESTART:-1}" == "1" ]]; then
  mk_namelist "$workdir/block_rs.nml" block "$workdir/hist_block_rs" "$workdir/rest_block"
  run_case "$workdir/block_rs.nml" "$ranks" block_restart_same || status=1
  "$wrapper" "$workdir/block_rs.nml" "$workdir/hist_block_rs" || status=1

  # --- phase 4: restart mid-window, different layout ----------------------
  mk_namelist "$workdir/block_rl.nml" block "$workdir/hist_block_rl" "$workdir/rest_block"
  run_case "$workdir/block_rl.nml" "$ranks_alt" block_restart_layout || status=1
  "$wrapper" "$workdir/block_rl.nml" "$workdir/hist_block_rl" || status=1
fi

# --- comparison ------------------------------------------------------------
echo
echo "== comparing against the 'one' reference"
compare () {   # $1 = candidate dir, $2 = label
  python3 "$repo_root/tests/river_hist_compare.py" \
    "$workdir/hist_one" "$1" --label "$2" || return 1
}
compare "$workdir/hist_block"    "block+postprocess"     || status=1
if [[ "${RH_PARITY_RESTART:-1}" == "1" ]]; then
  compare "$workdir/hist_block_rs" "restart, same layout"  || status=1
  compare "$workdir/hist_block_rl" "restart, other layout" || status=1
fi

echo
if [[ $status -eq 0 ]]; then
  echo "river history mode parity: PASS"
else
  echo "river history mode parity: FAILED" >&2
fi
exit $status
