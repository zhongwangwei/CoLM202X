#!/bin/bash
# methane10 two-year W1 runbook script (one annual job per submit; not launched
# in this task -- Lead runs it after choosing the methane mode).
#
# W1 = 1995-01-01 03:00 -> 1996-12-31 18:00, 35,070 steps (zero external
# CRUJRA boundary files; production time-index probe).  Two jobs:
#   Y1: mkinidata at 1995-01-01 03:00 + colm to 1996-01-01 00:00 (17,514 steps)
#   Y2: colm from the Y1 YEARLY restart to 1996-12-31 18:00 (17,556 steps)
#
# Required inputs (declare with --input-remote):
#   LANDDATA FORCING RT [RAWROOT] BUILD [MKINI_RESTART] [MKINI_XML]
# Env:
#   CASE_TAG  y1|y2 (default y1)
#   CASE_MODE wetwat (runnable, default) | hybrid (blocked by host-water guard
#             in the delivered code, see EVIDENCE.md sec.5)
#   CASE_TOL  host_water_tolerance (default 0.20 for wetwat)
#   CASE_RICE true|false (default false)
#   OBREF_INPUT_BUILD          frozen build directory (contains run/colm.x)
#   OBREF_INPUT_MKINI_RESTART  Y1: unused; Y2: Y1 restart directory (read-only)
#   OBREF_INPUT_MKINI_XML      optional: existing mkinidata case.nml for Y1
set -o pipefail
export PATH=/home/zhwei/miniconda3/bin:$PATH
export OMP_NUM_THREADS=1

TAG=${CASE_TAG:-y1}
CMODE=${CASE_MODE:-wetwat}
CTOL=${CASE_TOL:-0.20}
CRICE=${CASE_RICE:-false}
COLM_NP=${COLM_NP:-64}
MKINI_NP=${MKINI_NP:-64}

if [ "$TAG" = y1 ]; then
  START_Y=1995 START_M=1 START_D=1 START_S=10800
  END_Y=1996   END_M=1 END_D=1 END_S=0
  EXPECT_STEPS=17514 EXPECT_LAST="1995-12-31-84600"
else
  START_Y=1996 START_M=1 START_D=1 START_S=0
  END_Y=1996   END_M=12 END_D=31 END_S=64800
  EXPECT_STEPS=17556 EXPECT_LAST="1996-12-31-63000"
fi

ROOT=$(pwd)
W=$PWD
OUT=${OBREF_JOB_RESULTS:-$W/results}
mkdir -p "$OUT/logs" "$OUT/phases"
LANDDATA=${OBREF_INPUT_LANDDATA:?missing input}
FORCING=${OBREF_INPUT_FORCING:?missing input}
RT=${OBREF_INPUT_RT:?missing input}
RAWROOT=${OBREF_INPUT_RAWROOT:-$RT}
BINROOT=${OBREF_INPUT_BUILD:?missing build input}
export ROOT W OUT LANDDATA FORCING RT RAWROOT BINROOT COLM_NP MKINI_NP

source "$ROOT/lib.sh"
export -f mkparam case_nml seed_restart run_segment phase_json

{
  echo "date: $(date -Is)"
  echo "tag: $TAG mode: $CMODE tol: $CTOL rice: $CRICE"
  echo "window: $START_Y-$START_M-$START_D sec=$START_S -> $END_Y-$END_M-$END_D sec=$END_S"
  echo "OBREF_JOB_ID: ${OBREF_JOB_ID:-unknown}"
  env | grep '^OBREF_' | sort
  sha256sum "$BINROOT/run/colm.x" "$BINROOT/run/mkinidata.x"
  du -sh "$LANDDATA" "$FORCING"
  df -h "$W" "$OUT" | head -4
  uptime
} > "$OUT/logs/environment.txt" 2>&1
cat "$OUT/logs/environment.txt"

HISTVARS='f_methane_surf_flux_tot,f_methane_surf_flux_tot_active,f_methane_surf_flux_active_total_without_lake,f_methane_surf_flux_global_total_with_lake,f_methane_surf_flux_global_phys_with_lake,f_methane_balance_residual_global_with_lake,f_methane_ch4_clip_credit_global_with_lake,f_methane_surf_flux_tot_phys,f_methane_balance_residual,f_methane_ch4_clip_credit,f_methane_surf_flux_wetland,f_methane_surf_flux_soil,f_methane_surf_flux_lake,f_methane_surf_flux_rice,f_methane_surf_flux_rice_intensive,f_methane_prod_tot,f_methane_prod_tot_soil,f_methane_prod_tot_rice,f_methane_oxid_tot,f_methane_oxid_tot_soil,f_methane_oxid_tot_rice,f_co2_hr_methane,f_totcol_methane,f_methane_finundated,f_methane_soil_finundated,f_methane_soil_zwt,f_o2_cap_loss,f_o2_cap_gain'
STD="$ROOT/new/run/standard_methane_parameter.nml"
YEAR="DEF_WRST_FREQ = 'YEARLY'
DEF_HIST_FREQ = 'MONTHLY'
DEF_HIST_groupby = 'YEAR'"

mkparam "$STD" "$OUT/param_$TAG.nml" "$CRICE" false "$HISTVARS" "$CTOL" "$CMODE" || exit 1

# mkinidata only in Y1
if [ "$TAG" = y1 ]; then
  case_nml "mk_$TAG" "$START_Y" "$START_M" "$START_D" "$START_S" "$END_Y" "$END_M" "$END_D" "$END_S" \
      "$OUT/param_$TAG.nml" "$YEAR" || { phase_json "${TAG}_mkini" FAIL "case_nml"; exit 1; }
  if run_segment "mk_$TAG" mkinidata "$MKINI_NP" && \
     grep -q 'CoLM Initialization Execution Completed' "$W/pilot/mk_$TAG/case/log.mkinidata"; then
    phase_json "${TAG}_mkini" PASS
    MKINI_RESTART=$W/pilot/mk_$TAG/case/restart
  else
    phase_json "${TAG}_mkini" FAIL "log=$W/pilot/mk_$TAG/case/log.mkinidata"
    exit 1
  fi
else
  MKINI_RESTART=${OBREF_INPUT_MKINI_RESTART:?Y2 needs the Y1 restart tree as --input-remote MKINI_RESTART=...}
fi
export MKINI_RESTART

case_nml "$TAG" "$START_Y" "$START_M" "$START_D" "$START_S" "$END_Y" "$END_M" "$END_D" "$END_S" \
    "$OUT/param_$TAG.nml" "$YEAR" || { phase_json "$TAG" FAIL "case_nml"; exit 1; }
seed_restart "$TAG" || { phase_json "$TAG" FAIL "seed_restart"; exit 1; }

if run_segment "$TAG" colm "$COLM_NP"; then
  python3 "$ROOT/validate.py" --single "$W/pilot/$TAG/case" --tool colm --steps "$EXPECT_STEPS" \
      --last "$EXPECT_LAST" --out "$OUT/phases/$TAG.json" || true
  python3 - "$OUT/phases/$TAG.json" "$W/pilot/$TAG/case/log.colm" <<'PYX'
import json, re, sys
rep, log = sys.argv[1], sys.argv[2]
d = json.load(open(rep))
txt = open(log, 'r', errors='replace').read()
worst = re.findall(r'worst disagreement over all ranks\s+([0-9.Ee+-]+) of pore volume', txt)
d['host_water_worst'] = max((float(x) for x in worst), default=None)
d['n_inundation_abort'] = txt.count('ERROR: methane inundation exceeds host soil water')
d['n_o2_roundoff_warnings'] = txt.count('O2 roundoff adjustment')
json.dump(d, open(rep, 'w'), indent=1)
PYX
  if python3 "$ROOT/validate.py" --gate "$OUT/phases/$TAG.json"; then
    phase_json "$TAG" PASS
  else
    phase_json "$TAG" FAIL "gate"
  fi
else
  phase_json "$TAG" FAIL "log=$W/pilot/$TAG/case/log.colm"
fi

cat > "$OUT/expectations.json" <<JSON
{ "$TAG": {"tool": "colm", "expect_steps": $EXPECT_STEPS, "expect_last": "$EXPECT_LAST"} }
JSON
python3 "$ROOT/validate.py" --all "$W" --out "$OUT" 2>&1 | tee "$OUT/logs/validation.log"
echo "=== validation rc=${PIPESTATUS[0]}"
du -sh "$W/pilot/$TAG/case/history" "$W/pilot/$TAG/case/restart" 2>/dev/null
exit 0
