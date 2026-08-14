#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Step 1 driver for .omx/plans/river-history-sharded-output.md
#
# Builds tests/river_hist_baseline_harness.F90 against the REAL model objects
# (lib/libcolm.a + .bld module directory) so the reference file and the timings
# come from the actual `one`-mode writer, then runs it at two worker/IO-group
# scales with repeats and reports medians (plan 1.3, 1.4).
#
# Requires a prior model build, e.g.
#   ln -sf Makeoptions.Mac-arm include/Makeoptions && gmake all TRACER_ENABLED=YES
#
# Usage:
#   tests/run_river_hist_baseline.sh [small|large]
# Environment overrides are forwarded to the harness (RH_* variables).
# ---------------------------------------------------------------------------
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

profile=${1:-small}

compiler=${MPIFC:-$(command -v mpif90 || true)}
if [[ -z "$compiler" ]]; then
  echo "ERROR: mpif90 is required for the river history baseline" >&2
  exit 1
fi

launcher=${MPIEXEC:-$(command -v mpiexec || command -v mpirun || true)}
if [[ -z "$launcher" ]]; then
  echo "ERROR: mpiexec or mpirun is required for the river history baseline" >&2
  exit 1
fi

blddir=${COLM_BLD_DIR:-$repo_root/.bld}
libcolm=${COLM_LIB:-$repo_root/lib/libcolm.a}

if [[ ! -d "$blddir" || ! -f "$libcolm" ]]; then
  echo "ERROR: model build not found ($blddir, $libcolm)." >&2
  echo "       Build the model first, then re-run this script." >&2
  exit 1
fi

# netcdf-fortran and netcdf-c frequently live under different prefixes
# (Homebrew Cellar), so prefer the *-config helpers over a single prefix.
if [[ -n "${NETCDF_PREFIX:-}" ]]; then
  nc_inc=(-I"$NETCDF_PREFIX/include")
  nc_lib=(-L"$NETCDF_PREFIX/lib" -lnetcdff -lnetcdf)
else
  nc_inc=()
  nc_lib=()
  if command -v nf-config >/dev/null 2>&1; then
    read -r -a _fflags <<< "$(nf-config --fflags)"; nc_inc+=("${_fflags[@]}")
    read -r -a _flibs  <<< "$(nf-config --flibs)";  nc_lib+=("${_flibs[@]}")
  fi
  if command -v nc-config >/dev/null 2>&1; then
    read -r -a _clibs <<< "$(nc-config --libs)";    nc_lib+=("${_clibs[@]}")
  fi
  if [[ ${#nc_lib[@]} -eq 0 ]]; then
    nc_inc=(-I/usr/local/include); nc_lib=(-L/usr/local/lib -lnetcdff -lnetcdf)
  fi
fi

# The *-config helpers emit -L but no -rpath, so the shared netCDF libraries
# are not found at run time on macOS. Mirror every -L as an -rpath.
nc_rpath=()
for _flag in "${nc_lib[@]}"; do
  if [[ "$_flag" == -L* ]]; then
    nc_rpath+=("-Wl,-rpath,${_flag#-L}")
  fi
done

workdir=$(mktemp -d /tmp/colm-river-hist-baseline.XXXXXX)
trap 'rm -rf "$workdir"' EXIT

flags=(
  -fopenmp
  -fdefault-real-8
  -ffree-form
  -cpp
  -ffree-line-length-0
  -fallow-argument-mismatch
  -w
)

echo "== building baseline harness against $libcolm"
"$compiler" "${flags[@]}" -DUSEMPI \
  -I"$blddir" "${nc_inc[@]}" -J"$workdir" \
  -c tests/river_hist_baseline_harness.F90 -o "$workdir/harness.o"
"$compiler" "${flags[@]}" \
  "$workdir/harness.o" "$libcolm" \
  "${nc_lib[@]}" "${nc_rpath[@]}" -llapack -lblas \
  -o "$workdir/river_hist_baseline"

case "$profile" in
  small) : "${RH_TOTALNUMUCAT:=4000}";  : "${RH_NLON:=200}"; : "${RH_NLAT:=100}"
         : "${RH_TOTALNPTHOUT:=500}";   : "${RH_NVAR:=8}" ;;
  large) : "${RH_TOTALNUMUCAT:=200000}"; : "${RH_NLON:=1440}"; : "${RH_NLAT:=720}"
         : "${RH_TOTALNPTHOUT:=20000}";  : "${RH_NVAR:=16}" ;;
  *) echo "ERROR: unknown profile '$profile' (expected small|large)" >&2; exit 1 ;;
esac
: "${RH_NPTHLEV:=3}"
: "${RH_NTIME:=2}"
: "${RH_REPEAT:=1}"
: "${RH_RUNS:=3}"
export RH_TOTALNUMUCAT RH_NLON RH_NLAT RH_TOTALNPTHOUT RH_NVAR RH_NPTHLEV RH_NTIME RH_REPEAT

# Two scales, so rank scaling is measured rather than inferred (plan 1.4).
scales=${RH_SCALES:-"4 8"}

summary="$workdir/summary.txt"
: > "$summary"

for ranks in $scales; do
  # one master + one IO rank per group; keep >=2 workers at every scale.
  groups=$(( ranks / 4 )); [[ $groups -lt 1 ]] && groups=1
  for run in $(seq 1 "${RH_RUNS}"); do
    out="$workdir/ref_r${ranks}_run${run}.nc"
    RH_NGROUP=$groups RH_OUT="$out" \
    env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
        OMPI_MCA_rmaps_base_oversubscribe=1 \
      "$launcher" -n "$ranks" "$workdir/river_hist_baseline" \
      | tee "$workdir/log_r${ranks}_run${run}.txt" \
      | sed -n 's/^RHBASE /RHBASE /p' >> "$summary" || {
        echo "ERROR: harness failed at ranks=$ranks run=$run" >&2
        cat "$workdir/log_r${ranks}_run${run}.txt" >&2
        exit 1
      }
    # keep the last reference file of each scale for schema locking
    cp -f "$out" "$repo_root/tests/artifacts/river_hist_ref_r${ranks}.nc" 2>/dev/null || {
      mkdir -p "$repo_root/tests/artifacts"
      cp -f "$out" "$repo_root/tests/artifacts/river_hist_ref_r${ranks}.nc"
    }
  done
done

echo
echo "== medians over ${RH_RUNS} runs per scale (profile=$profile)"
python3 - "$workdir" $scales <<'PY'
import re, statistics, sys, pathlib
workdir = pathlib.Path(sys.argv[1])
scales = sys.argv[2:]
keys = ["gather_only_s", "full_write_s", "netcdf_est_s", "bif_matrix_s"]
print(f"{'ranks':>6} " + " ".join(f"{k:>16}" for k in keys))
prev = {}
for ranks in scales:
    vals = {k: [] for k in keys}
    mem = {}
    for log in sorted(workdir.glob(f"log_r{ranks}_run*.txt")):
        for line in log.read_text().splitlines():
            m = re.match(r"RHBASE TIME (\S+)\s+(\S+)", line)
            if m and m.group(1) in vals:
                vals[m.group(1)].append(float(m.group(2)))
            m = re.match(r"RHBASE MEM (\S+)\s+(\S+)", line)
            if m:
                mem[m.group(1)] = int(m.group(2))
    med = {k: statistics.median(v) if v else float('nan') for k, v in vals.items()}
    print(f"{ranks:>6} " + " ".join(f"{med[k]:16.6e}" for k in keys))
    prev[ranks] = med
    if mem:
        print(f"       root lower-bound: ucat={mem.get('ucat_bytes',0)/1e6:.1f} MB "
              f"bif={mem.get('bif_bytes',0)/1e6:.1f} MB")
if len(scales) == 2:
    lo, hi = scales
    for k in keys:
        a, b = prev[lo][k], prev[hi][k]
        if a > 0:
            print(f"  scaling {k}: ranks {lo}->{hi} factor {b/a:.2f}x "
                  f"({'grows' if b > a else 'shrinks'})")
PY

cp -f "$summary" "$repo_root/tests/artifacts/river_hist_baseline_summary.txt" 2>/dev/null || true

# Verify the artifacts we just produced, rather than trusting that they are
# fine because the runs exited zero: the schema lock is the thing that would
# notice an id landing in the wrong cell.
echo
echo "== validating the generated reference files"
python3 - "$repo_root" $scales <<PYEOF
import sys, pathlib
root = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(root / "tests"))
from river_hist_schema_lock import (describe, compare_schema,
                                    check_ucat_values, check_bif_values)
scales = sys.argv[2:]
files = [root / f"tests/artifacts/river_hist_ref_r{r}.nc" for r in scales]
problems = []
if len(files) == 2:
    problems += [f"r{scales[0]} vs r{scales[1]}: {p}" for p in
                 compare_schema(describe(files[0]), describe(files[1]))]
for f, r in zip(files, scales):
    for ivar in range(1, ${RH_NVAR} + 1):
        for itime in range(1, ${RH_NTIME} + 1):
            problems += [f"r{r}: {p}" for p in check_ucat_values(
                str(f), f"f_ucat_synth{ivar:02d}", ${RH_TOTALNUMUCAT},
                ${RH_NLON}, ${RH_NLAT}, ivar, itime)]
    for itime in range(1, ${RH_NTIME} + 1):
        problems += [f"r{r}: {p}" for p in check_bif_values(
            str(f), "f_bifflw_lev", ${RH_NPTHLEV}, ${RH_TOTALNPTHOUT}, itime)]
for p in problems[:20]:
    print("   FAIL", p)
if problems:
    sys.exit(1)
print("   layout-invariant schema; every global id in the right cell")
PYEOF
rc=$?
if [[ $rc -ne 0 ]]; then
  echo "river history baseline: FAILED validation" >&2
  exit 1
fi

echo
echo "river history baseline: PASS (reference files in tests/artifacts/)"
