#!/bin/bash
# p9s long-window transparency probe: re-run only the target month from an
# existing monthly restart with private trace instrumentation, verify that the
# overlapping history is bit-identical to the baseline, and keep the trace text
# for the Decimal50 replay on the Mac.
set -o pipefail
export PATH=/home/zhwei/miniconda3/bin:$PATH OMP_NUM_THREADS=1
ROOT=$(pwd); OUT=${OBREF_JOB_RESULTS:-/tmp/out}
mkdir -p "$OUT/build" "$OUT/bld" "$OUT/runs" "$OUT/trace"
P9S=$OBREF_INPUT_P9S ISOGSM=$OBREF_INPUT_ISOGSM
RAW=$OBREF_INPUT_RAW RT=$OBREF_INPUT_RT JRA=$OBREF_INPUT_JRA UCAT=$OBREF_INPUT_UCAT
cp "$ROOT/targets.json" "$ROOT/probes.json" "$ROOT/compare_subset.py" "$OUT/"
[ -f "$ROOT/selection.json" ] && cp "$ROOT/selection.json" "$OUT/"

d=$OUT/bld/new; cp -r "$ROOT/new" "$d" || exit 1; cd "$d" || exit 1
bash .github/workflows/create_defineh.bash GRID LULC_IGBP URBANOFF vanGenu CaMaOFF BGCOFF CROPOFF MTON > "$OUT/build/defines.log"
grep -q '^#define GridRiverLakeFlow' include/define.h || sed -i 's/^#undef CaMa_Flood$/#undef CaMa_Flood\n#define GridRiverLakeFlow/' include/define.h
grep -q '^#define GridRiverLakeFlow' include/define.h || exit 1
grep -q '^#define MaterialTransport' include/define.h || exit 1
! grep -q '^#define LULCC' include/define.h || exit 1
make all FF=mpif90 NETCDF_INC="$(nf-config --includedir)" NETCDF_LIB="$(nc-config --libdir)" \
  FOPTS='-fdefault-real-8 -fdefault-double-8 -ffree-form -O2 -cpp -ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium' \
  > "$OUT/build/new.log" 2>&1
test -x run/colm.x || { echo "BUILD FAIL"; tail -40 "$OUT/build/new.log"; exit 1; }
echo "BUILD PASS $(date -Is)"
grep -q 'trace_write' main/MATERIAL/MOD_MT_Defs.F90 || { echo "instrumentation missing"; exit 1; }
echo "INSTRUMENTATION PRESENT"
cd "$ROOT" || exit 1

# Reproduce the reported extremes from the baseline history files themselves.
if [ -f "$ROOT/selection.json" ] && [ -f "$ROOT/report_p9s.json" ]; then
  mkdir -p "$OUT/logs"
  python3 "$ROOT/verify_extremes.py" "$P9S" "$ROOT/selection.json" > "$OUT/logs/extreme-verify.json" 2>&1 || {
     echo "EXTREME VERIFICATION FAILED"; cat "$OUT/logs/extreme-verify.json"; exit 1; }
  echo "EXTREMES VERIFIED FROM BASELINE HISTORY"
  cat "$OUT/logs/extreme-verify.json"
else
  echo "no selection.json/report_p9s.json: extreme reproduction check skipped"
fi

sed 's/DEF_ISO_USE_FRACTIONATION *=.*/DEF_ISO_USE_FRACTIONATION = .true./' "$OUT/bld/new/run/standard_isotope_parameter.nml" > "$OUT/iso_frac.nml"
sed "s#DEF_dir_forcing *=.*#DEF_dir_forcing = '$ISOGSM/'#; s#DEF_forcing%startyr *=.*#DEF_forcing%startyr = 2000#; s#DEF_forcing%endyr *=.*#DEF_forcing%endyr = 2004#" \
  "$OUT/bld/new/run/forcing/IsoGSM.nml" | sed -n '/^&nl_colm_forcing/,/^\/$/p' > "$OUT/forcing_isogsm.nml"

python3 - > "$OUT/segments.txt" <<'EOF'
import json
for p in json.load(open('probes.json')):
    print('|'.join([p['name'], p['baseline_case'], str(p['start_year']), str(p['start_month']),
                    str(p['start_day']), str(p['end_year']), str(p['end_month']), str(p['end_day'])]))
EOF

ptask() {
  IFS='|' read -r name case sy sm sd ey em ed <<< "$1"
  local c=$OUT/runs/$name src=$P9S/runs/$case
  mkdir -p "$c/case" || return 1
  cp -r "$src/case/landdata" "$c/case/" || return 1
  cp -r "$src/case/restart" "$c/case/" || return 1
  cp "$src/forcing.nml" "$c/forcing.nml" || return 1
  sed -e "s#^DEF_dir_output  *=.*#DEF_dir_output  = '$c/'#" \
      -e "s#^DEF_rawdata_namelist *=.*#DEF_rawdata_namelist = '$OUT/bld/new/run/rawdata/colm2024.nml'#" \
      -e "s#^DEF_forcing_namelist *=.*#DEF_forcing_namelist = '$c/forcing.nml'#" \
      -e "s#^DEF_ISOTOPE_PARAM_FILE *=.*#DEF_ISOTOPE_PARAM_FILE = '$OUT/iso_frac.nml'#" \
      -e "s#^DEF_simulation_time%start_year *=.*#DEF_simulation_time%start_year    = $sy#" \
      -e "s#^DEF_simulation_time%start_month *=.*#DEF_simulation_time%start_month   = $sm#" \
      -e "s#^DEF_simulation_time%start_day *=.*#DEF_simulation_time%start_day     = $sd#" \
      -e "s#^DEF_simulation_time%start_sec *=.*#DEF_simulation_time%start_sec     = 0#" \
      -e "s#^DEF_simulation_time%end_year *=.*#DEF_simulation_time%end_year      = $ey#" \
      -e "s#^DEF_simulation_time%end_month *=.*#DEF_simulation_time%end_month     = $em#" \
      -e "s#^DEF_simulation_time%end_day *=.*#DEF_simulation_time%end_day       = $ed#" \
      -e "s#^DEF_simulation_time%end_sec *=.*#DEF_simulation_time%end_sec       = 86400#" \
      "$src/case.nml" > "$c/case.nml" || return 1
  ( cd "$c" && mpiexec -n 8 "$OUT/bld/new/run/colm.x" case.nml > "log.colm" 2>&1 ) || {
     echo "PROBE FAIL $name run"; tail -20 "$c/log.colm"; return 1; }
  grep -q 'CoLM Execution Completed' "$c/log.colm" || { echo "PROBE FAIL $name no END"; return 1; }
  echo "PROBE END $name $(date -Is)"
  python3 "$ROOT/compare_subset.py" "$c/case/history" "$P9S/runs/$case/case/history" > "$c/compare.json" || {
     echo "PROBE FAIL $name history not transparent"; return 1; }
  echo "TRANSPARENCY PASS $name"
  mkdir -p "$OUT/trace/$name"
  ls "$c"/trace_rank_*.txt >/dev/null 2>&1 || { echo "PROBE FAIL $name no trace files"; return 1; }
  cp "$c"/trace_rank_*.txt "$OUT/trace/$name/" || return 1
  cp "$c/log.colm" "$c/case.nml" "$c/forcing.nml" "$c/compare.json" "$OUT/trace/$name/"
  echo "TRACE COPIED $name $(ls "$c"/trace_rank_*.txt | wc -l) rank files"
}
export -f ptask
export ROOT OUT P9S ISOGSM RAW RT JRA UCAT
rc=0
if [ -n "$PROBE_ONLY" ]; then
  grep -F "$PROBE_ONLY" "$OUT/segments.txt" > "$OUT/segments.filtered" || { echo "no probe matches PROBE_ONLY=$PROBE_ONLY"; exit 1; }
  echo "PROBE_ONLY=$PROBE_ONLY selected:"; cat "$OUT/segments.filtered"
  mv "$OUT/segments.filtered" "$OUT/segments.txt"
fi
xargs -P 6 -I{} bash -c 'ptask "$1"' _ '{}' < "$OUT/segments.txt" || rc=1
python3 - "$OUT" <<'EOF'
import json, re, sys
from pathlib import Path
out = Path(sys.argv[1]); rep = {'probes': {}, 'failures': []}
for d in sorted((out / 'runs').glob('*')):
    c = d / 'log.colm'
    if not c.exists():
        rep['failures'].append(f'{d.name}: no log'); continue
    s = c.read_text(errors='replace')
    steps = re.findall(r'TIMESTEP\s*=\s*(\d+)\s*\|\s*DATE\s*=\s*([\d-]+)', s)
    ends = 'CoLM Execution Completed' in s
    cmpf = d / 'compare.json'
    tracedir = out / 'trace' / d.name
    rep['probes'][d.name] = dict(
        end=ends, steps=len(steps), last=steps[-1] if steps else None,
        trace_files=len(list(tracedir.glob('trace_rank_*.txt'))) if tracedir.exists() else 0,
        transparency=json.load(open(cmpf)) if cmpf.exists() else None,
        bad=[x for x in ('TRC_SIG step report', 'TRC_BAL step report', 'residual hard failure',
                         'ERROR STOP', 'MPI_Abort') if x in s])
    if not ends:
        rep['failures'].append(f'{d.name}: no END')
    if rep['probes'][d.name]['transparency'] is None:
        rep['failures'].append(f'{d.name}: no transparency comparison')
(out / 'report.json').write_text(json.dumps(rep, indent=2) + '\n')
print(json.dumps(rep, indent=2))
sys.exit(bool(rep['failures']))
EOF
r2=$?
[ $r2 = 0 ] || rc=1
[ $rc = 0 ] && echo "PROBE JOB PASS" || echo "PROBE JOB FAIL"
exit $rc
