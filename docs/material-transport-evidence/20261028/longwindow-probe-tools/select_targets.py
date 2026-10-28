#!/usr/bin/env python
"""Select real p9s long-window extremes as probe targets.

Reads the p9s `report.json` (audit output), keeps only the evaporation-ratio
diagnostics that the private trace instrumentation covers
(`f_mt_delta_soilevap_*`, `f_mt_delta_wetland_evap_*`), and turns the reported
min/max records into:

* targets.json  - lat/lon/doy/name/host for instrument.py and the replay;
* probes.json   - one segment per (case, month): start at the month first day,
                  ending on the latest selected target day in that month;
* selection.json- the full ranked evidence used to choose them.

The selection is by the reported delta magnitude, never by convenience: every
target keeps its source file, record index, delta, lat and lon.
"""
import datetime as dt
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent
MONTH_START_DOY = {1: 1, 2: 32, 3: 61, 4: 92, 5: 122, 6: 153,
                   7: 184, 8: 215, 9: 246, 10: 276, 11: 307, 12: 337}


def record_date(fname, time_index):
    stem = Path(fname).stem                      # case_hist_material_2004-02
    year, month = (int(x) for x in stem.split('_')[-1].split('-'))
    return dt.date(year, month, 1) + dt.timedelta(days=int(time_index))


def doy(date):
    return int(date.strftime('%j'))


def main(report_path, per_case, out_dir):
    report = json.loads(Path(report_path).read_text())
    candidates = []
    for case, data in report['cases'].items():
        for var, item in data['variables'].items():
            host = 'soil' if 'soilevap' in var else 'wet' if 'wetland_evap' in var else None
            if host is None:
                continue
            species = 'H2_18O' if var.endswith('H2_18O') else 'HDO' if var.endswith('HDO') else None
            if species is None:
                continue
            for ext in item['extremes']:
                if len(ext['index']) != 3:
                    continue  # need (time, lat, lon) to fix the date/cell
                date = record_date(ext['file'], ext['index'][0])
                doy_ = doy(date)
                candidates.append(dict(
                    case=case, var=var, host=host, species=species, mode=ext['mode'],
                    delta=float(ext['delta']), lat=float(ext['lat']), lon=float(ext['lon']),
                    file=ext['file'], index=[int(i) for i in ext['index']],
                    date=date.isoformat(), doy=doy_,
                    month=f'{date.year:04d}-{date.month:02d}',
                    restart=f'{date.year:04d}-{MONTH_START_DOY[date.month]:03d}-00000',
                    start_day=1, end_day=date.day))
    # rank by reported magnitude, but keep every record for the evidence table
    candidates.sort(key=lambda c: abs(c['delta']), reverse=True)
    chosen = []
    for case in sorted(report['cases']):
        picked = []
        for cand in candidates:
            if cand['case'] != case or cand in picked:
                continue
            if len(picked) >= per_case:
                break
            # prefer coverage of both isotopes, both hosts and both signs
            picked.append(cand)
        chosen.extend(picked)
    targets, probes, seen_names = [], [], set()
    for cand in chosen:
        name = f"{cand['case']}-{cand['host']}-{cand['species']}-{cand['mode']}-{cand['date']}"
        if name in seen_names:
            continue
        seen_names.add(name)
        targets.append(dict(name=name, lat=cand['lat'], lon=cand['lon'], doy=cand['doy'],
                            host=cand['host'], case=cand['case'], species=cand['species'],
                            expected_delta=cand['delta'], source=cand))
    by_segment = {}
    for t in targets:
        c = t['source']
        key = (c['case'], c['month'])
        seg = by_segment.setdefault(key, dict(
            name=f"{c['case']}-{c['month']}", baseline_case=c['case'],
            start_year=int(c['month'][:4]), start_month=int(c['month'][5:7]), start_day=1,
            end_year=int(c['date'][:4]), end_month=int(c['date'][5:7]), end_day=c['end_day'],
            targets=[]))
        seg['end_day'] = max(seg['end_day'], c['end_day'])
        seg['targets'].append(t['name'])
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'targets.json').write_text(json.dumps(targets, indent=2) + '\n')
    (out / 'probes.json').write_text(json.dumps(list(by_segment.values()), indent=2) + '\n')
    (out / 'selection.json').write_text(json.dumps(
        dict(report=str(report_path), per_case=per_case, chosen=targets,
             ranked_candidates=candidates[:40]), indent=2) + '\n')
    print(json.dumps(dict(n_targets=len(targets), n_probes=len(by_segment),
                          segments={k[0] + '/' + k[1]: (v['start_day'], v['end_day'])
                                    for k, v in by_segment.items()},
                          targets=[(t['name'], t['lat'], t['lon'], t['doy']) for t in targets]),
                     indent=1))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2, sys.argv[3] if len(sys.argv) > 3 else '.')
