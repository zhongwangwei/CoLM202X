#!/bin/bash
# Validate the carrierless-isotope quarantine fix:
#  1) build LULC_IGBP (no LULCC) and LULC_IGBP+LULCC MT+GridRiverLakeFlow binaries from the modified tree
#  2) run the 12-step/1-day flood-feedback + IsoGSM isotope trigger case on both
#  3) check model END, zero conservation diagnostics, and that the isotope surface
#     residue is now a finite conservative inventory instead of a numerical residual
#  4) run the repository material-transport self checks
set -o pipefail
export PATH=/home/zhwei/miniconda3/bin:$PATH OMP_NUM_THREADS=1
ROOT=$(pwd); OUT=${OBREF_JOB_RESULTS:-/tmp/out}; SRC=$ROOT/src
mkdir -p "$OUT/build" "$OUT/runs" "$OUT/logs" "$OUT/self"
RAW=$OBREF_INPUT_RAW; RT=$OBREF_INPUT_RT; UCAT=$OBREF_INPUT_UCAT
P6B=$OBREF_INPUT_P6B; ISOGSM=$OBREF_INPUT_ISOGSM2
{
  echo "date: $(date -Is)"; echo "host: $(hostname)"; echo "OBREF_JOB_ID: ${OBREF_JOB_ID:-unknown}"
  env | grep '^OBREF_' | sort; mpif90 --version 2>&1 | head -1; gfortran --version 2>&1 | head -1
} > "$OUT/logs/environment.txt" 2>&1
echo "=== source diff vs branch head (information only) ===" > "$OUT/logs/source-note.txt"
grep -n "tracer_is_isotope" "$SRC/main/MATERIAL/GENERIC_WATER/MOD_MT_SoilWater.F90" | wc -l >> "$OUT/logs/source-note.txt"

build() {  # name lulc [lulcc]
  local name=$1 lulc=$2 lulcc=${3:-}
  local d=$OUT/build/$name; rm -rf "$d"; cp -r "$SRC" "$d" || return 1; cd "$d" || return 1
  if [ "$lulcc" = lulcc ]; then
    bash .github/workflows/create_defineh.bash GRID "$lulc" URBANOFF vanGenu CaMaOFF BGCOFF CROPOFF MTON LULCCON > "$OUT/logs/define-$name.log"
  else
    bash .github/workflows/create_defineh.bash GRID "$lulc" URBANOFF vanGenu CaMaOFF BGCOFF CROPOFF MTON > "$OUT/logs/define-$name.log"
  fi
  grep -q '^#define GridRiverLakeFlow' include/define.h || sed -i 's/^#undef CaMa_Flood$/#undef CaMa_Flood\n#define GridRiverLakeFlow/' include/define.h
  grep -q '^#define GridRiverLakeFlow' include/define.h || { echo "BUILD $name: GridRiverLakeFlow missing"; return 1; }
  grep -q '^#define MaterialTransport' include/define.h || { echo "BUILD $name: MaterialTransport missing"; return 1; }
  if [ "$lulcc" = lulcc ]; then grep -q '^#define LULCC' include/define.h || { echo "BUILD $name: LULCC missing"; return 1; }
  else ! grep -q '^#define LULCC' include/define.h || { echo "BUILD $name: unexpected LULCC"; return 1; }; fi
  make all FF=mpif90 NETCDF_INC="$(nf-config --includedir)" NETCDF_LIB="$(nc-config --libdir)" \
    FOPTS='-fdefault-real-8 -fdefault-double-8 -ffree-form -O2 -cpp -ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium' \
    > "$OUT/logs/build-$name.log" 2>&1
  if [ -x run/colm.x ] && [ -x run/mkinidata.x ]; then
    cp include/define.h "$OUT/logs/define-$name.h"
    sha256sum run/colm.x run/mkinidata.x > "$OUT/logs/binary-$name.sha256"
    echo "BUILD PASS $name $(date -Is)"
    cd "$ROOT"; return 0
  fi
  echo "BUILD FAIL $name"; tail -30 "$OUT/logs/build-$name.log"; cd "$ROOT"; return 1
}

build plain LULC_IGBP
RC_PLAIN=$?
RC_LULCC=1
echo "build rc: plain=$RC_PLAIN (lulcc skipped for diagnosis)"

sed 's/DEF_ISO_USE_FRACTIONATION *=.*/DEF_ISO_USE_FRACTIONATION = .true./' \
    "$SRC/run/standard_isotope_parameter.nml" > "$OUT/iso_frac.nml"
sed "s#DEF_dir_forcing *=.*#DEF_dir_forcing = '$ISOGSM/'#; \
     s#DEF_forcing%startyr *=.*#DEF_forcing%startyr = 2005#; \
     s#DEF_forcing%endyr *=.*#DEF_forcing%endyr = 2006#" \
    "$SRC/run/forcing/IsoGSM.nml" | sed -n '/^&nl_colm_forcing/,/^\//p' > "$OUT/forcing_isogsm.nml"

run_case() {  # name build [lulcc]
  local name=$1 bld=$2 lulcc=${3:-}
  local c=$OUT/runs/$name; rm -rf "$c"; mkdir -p "$c/case"
  cp -r "$P6B/srf/LULC_IGBP/NEC/case/landdata" "$c/case/" || { echo "FAIL $name landdata"; return 1; }
  cp "$OUT/forcing_isogsm.nml" "$c/forcing.nml"
  local extra=""
  [ "$lulcc" = lulcc ] && extra="DEF_LULCC_SCHEME = 2"
  cat > "$c/case.nml" <<NML
&nl_colm
DEF_CASE_NAME = 'case'
DEF_domain%edgew = 125
DEF_domain%edgee = 127.5
DEF_domain%edges = 47.5
DEF_domain%edgen = 50
DEF_simulation_time%greenwich = .true.
DEF_GRIDBASED_lon_res = 0.25
DEF_GRIDBASED_lat_res = 0.25
DEF_nx_blocks = 72
DEF_ny_blocks = 36
DEF_PIO_groupsize = 4
DEF_HIST_FREQ    = 'DAILY'
DEF_HIST_groupby = 'MONTH'
DEF_HIST_lon_res = 0.25
DEF_HIST_lat_res = 0.25
DEF_simulation_time%start_year    = 2005
DEF_simulation_time%start_month   = 12
DEF_simulation_time%start_day     = 25
DEF_simulation_time%start_sec     = 0
DEF_simulation_time%end_year      = 2005
DEF_simulation_time%end_month     = 12
DEF_simulation_time%end_day       = 25
DEF_simulation_time%end_sec       = 21600
DEF_simulation_time%spinup_year   = 2005
DEF_simulation_time%spinup_month  = 12
DEF_simulation_time%spinup_day    = 25
DEF_simulation_time%spinup_sec    = 0
DEF_simulation_time%spinup_repeat = 0
DEF_simulation_time%timestep      = 1800.
DEF_dir_rawdata = '$RAW/'
DEF_dir_runtime = '$RT/'
DEF_dir_output  = '$c/'
DEF_rawdata_namelist = '$OUT/build/$bld/run/rawdata/colm2024.nml'
DEF_LAI_MONTHLY       = .true.
DEF_LAI_CHANGE_YEARLY = .false.
DEF_forcing_namelist = '$c/forcing.nml'
DEF_WRST_FREQ = 'DAILY'
DEF_UnitCatchment_file = '$UCAT'
DEF_GridRiverLake_FloodFeedback = .true.
$extra
/

&nl_colm_material
DEF_USE_ISOTOPE = .true.
DEF_ISOTOPE_PARAM_FILE = '$OUT/iso_frac.nml'
DEF_MT_AQUIFER_MIXING_WATER_MM = 1000.
DEF_MT_BALANCE_ABORT_NBAD = 0
/
NML
  cd "$c" || return 1
  mpiexec -n 5 "$OUT/build/$bld/run/mkinidata.x" case.nml > log.mkinidata 2>&1 || { echo "FAIL $name mkinidata"; cd "$ROOT"; return 1; }
  ( cd "$c" && mpiexec -n 5 "$OUT/build/$bld/run/colm.x" case.nml > log.colm 2>&1 )
  echo "RUN DONE $name rc=$? $(date -Is)"
  cd "$ROOT"
}
if [ "$RC_PLAIN" = 0 ]; then run_case plain plain; fi
if [ "$RC_LULCC" = 0 ]; then run_case lulcc lulcc lulcc; fi

echo "=== validation ==="
python3 - "$OUT" "$SRC" <<'PY'
import json, re, sys, glob
from pathlib import Path
import numpy as np
import netCDF4 as nc
out = Path(sys.argv[1]); src = Path(sys.argv[2])
rep = {'runs': {}, 'self_checks': {}, 'failures': []}
for run in sorted((out / 'runs').glob('*')):
    name = run.name
    logf = run / 'log.colm'
    if not logf.exists():
        rep['failures'].append(f'{name}: no model log'); continue
    log = logf.read_text(errors='replace')
    steps = re.findall(r'TIMESTEP\s*=\s*(\d+)\s*\|\s*DATE\s*=\s*([\d-]+)', log)
    bad = [x for x in ('TRC_SIG step report', 'TRC_BAL step report', 'residual hard failure',
                       'ERROR STOP', 'MPI_Abort', 'Water balance violation', 'Energy balance violation')
           if x in log]
    entry = dict(end='CoLM Execution Completed' in log, steps=len(steps),
                 last=steps[-1] if steps else None, bad=bad)
    # restart: isotope surface residue inventory
    rdir = run / 'case/restart'
    dates = sorted([p for p in rdir.iterdir() if p.is_dir() and p.name[:1].isdigit()]) if rdir.exists() else []
    mres = {}
    for dd in dates[-1:]:
        for f in dd.rglob('*.nc'):
            with nc.Dataset(f) as ds:
                for key in ds.variables:
                    if 'surface_residue' not in key and 'numerical_residue' not in key:
                        continue
                    a = np.asarray(ds[key][:])
                    if a.size == 0: continue
                    mres[key] = dict(shape=list(a.shape), max_abs=float(np.max(np.abs(a))),
                                     nonzero=int(np.count_nonzero(a)))
    entry['restart_residue_fields'] = mres
    # history: residue variables and gross sanity
    hdir = run / 'case/history'
    hvars = {}
    for f in sorted(hdir.glob('*material*.nc')):
        with nc.Dataset(f) as ds:
            for key in ds.variables:
                if 'surface_residue' not in key: continue
                a = np.ma.filled(ds[key][:], np.nan)
                hvars[key] = dict(vmin=float(np.nanmin(a)), vmax=float(np.nanmax(a)),
                                  finite_masked=int(np.count_nonzero(~np.isfinite(a))))
    entry['history_residue_fields'] = hvars
    rep['runs'][name] = entry
    if not entry['end']: rep['failures'].append(f'{name}: no model END')
    if entry['bad']: rep['failures'].append(f'{name}: diagnostics {entry["bad"]}')
    if entry['steps'] != 12: rep['failures'].append(f'{name}: steps {entry["steps"]} != 12')

checks = {
    'compare': ['tests/material_transport/compare.py', '--self-test'],
    'iso_check': ['tests/material_transport/iso_check.py', '--self-test'],
    'forcing_cache': ['tests/material_transport/check_forcing_cache.py', '--self-test'],
    'forcing_edges': ['tests/material_transport/check_forcing_edges.py', str(src)],
    'singlebucket': ['tests/material_transport/check_singlebucket_precip.py'],
    'orphan_leaf': ['tests/material_transport/check_orphan_leaf_inventory.py'],
    'frcsat': ['tests/material_transport/check_frcsat_history.py', str(src)],
    'regionclip': ['tests/material_transport/check_regionclip_paths.py'],
    'defineh': ['tests/material_transport/check_defineh.py'],
}
import subprocess
for key, arg in checks.items():
    p = subprocess.run([sys.executable, str(src / arg[0])] + arg[1:], capture_output=True, text=True)
    rep['self_checks'][key] = dict(rc=p.returncode, tail=(p.stdout + p.stderr).strip().splitlines()[-3:])
    if p.returncode != 0:
        rep['failures'].append(f'self-check {key}: rc={p.returncode}')
(out / 'validation.json').write_text(json.dumps(rep, indent=2) + '\n')
print(json.dumps(rep, indent=2))
sys.exit(bool(rep['failures']))
PY
VRC=$?
echo "validation rc=$VRC"
exit "$VRC"
