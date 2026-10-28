#!/usr/bin/env python
"""Run the Decimal50 replay for every fetched probe and aggregate the results.

Usage: run_replay_all.py FETCHED_DIR  (contains one subdirectory per probe name,
with trace_rank_*.txt copied from the job's results/trace/<probe>/)

Writes replay/<probe>/ results plus an aggregate summary.json / summary.md and
targets-<probe>.json subset files.
"""
import json
import subprocess
import sys
from pathlib import Path

root = Path(__file__).parent
fetched = Path(sys.argv[1])
outroot = root / 'replay-results'
outroot.mkdir(exist_ok=True)
probes = json.loads((root / 'probes.json').read_text())
targets = json.loads((root / 'targets.json').read_text())
by_name = {t['name']: t for t in targets}

def dedupe(subset):
    """Merge records that share one traced (cell, doy, host): both isotope grids are
    traced for the cell, while the audit may report separate extremes per isotope."""
    order, merged = [], {}
    for t in subset:
        key = (t['lat'], t['lon'], t['doy'], t['host'])
        if key not in merged:
            merged[key] = dict(name=t['name'], lat=t['lat'], lon=t['lon'], doy=t['doy'],
                               host=t['host'], expected={}, records=[])
            order.append(key)
        merged[key]['expected'][t.get('species')] = t['expected_delta']
        merged[key]['records'].append(dict(name=t['name'], species=t.get('species'),
                                           expected=t['expected_delta'],
                                           var=t['source']['var'], date=t['name'][-10:]))
    return [merged[k] for k in order]


agg = {'probes': {}, 'failures': [], 'targets': {}, 'expected': {}}
for probe in probes:
    name = probe['name']
    trace = fetched / name
    if not trace.exists() or not list(trace.glob('trace_rank_*.txt')):
        agg['failures'].append(f'{name}: no trace files at {trace}')
        continue
    subset = [by_name[n] for n in probe['targets'] if n in by_name]
    if not subset:
        agg['failures'].append(f'{name}: no matching targets')
        continue
    unique = dedupe(subset)
    for u in unique:
        agg['expected'][u['name']] = u
    tpath = outroot / f'targets-{name}.json'
    tpath.write_text(json.dumps(unique, indent=2) + '\n')
    odir = outroot / name
    cmd = [sys.executable, str(root / 'replay_longwindow.py'), str(trace), str(odir), str(tpath)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    (outroot / f'{name}.stderr.txt').write_text(p.stderr)
    if p.returncode != 0:
        agg['failures'].append(f"{name}: replay rc={p.returncode}: {p.stderr.strip()[:500]}")
        continue
    summary = json.loads((odir / 'summary.json').read_text())
    daily = json.loads((odir / 'finite_daily_replay.json').read_text())
    agg['probes'][name] = dict(summary=summary, daily=daily)
    for g in summary['grids_detail']:
        agg['targets'].setdefault(g['target'], {})[g['var']] = dict(
            model_delta=g['model_delta'], decimal_delta=g['decimal_delta'],
            mask_reason=g['mask_reason'], signed_mass=g['signed_mass'],
            gross_water=g['gross_water'], patch_reconstruction_delta_error=g['patch_reconstruction_delta_error'],
            patches=g['patches'], eligible_patches=g['eligible_patches'])

    for t in unique:
        got = agg['targets'].get(t['name'], {})
        if not got:
            agg['failures'].append(f"{name}: target {t['name']} produced no grid record")
        for sp, exp in t['expected'].items():
            var = [k for k in got if k.endswith(sp)]
            if not var:
                agg['failures'].append(f"{name}: target {t['name']} missing {sp} grid")
            elif abs(float(got[var[0]]['model_delta']) - float(exp)) > 1e-9:
                agg['failures'].append(f"{name}: target {t['name']} {sp} model {got[var[0]]['model_delta']} "
                                       f"!= reported extreme {exp}")

(outroot / 'summary.json').write_text(json.dumps(agg, indent=2) + '\n')

lines = ['# p9s long-window probe replay aggregate', '']
if agg['failures']:
    lines += ['## Failures', ''] + [f'- {f}' for f in agg['failures']] + ['']
lines += ['| target | expected | variable | model delta | decimal delta | mask | gross water | signed mass | patches (eligible) |',
          '| --- | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: |']
for tname, exp in agg['expected'].items():
    for sp in ('H2_18O', 'HDO'):
        got = {k: v for k, v in agg['targets'].get(tname, {}).items() if k.endswith(sp)}
        if not got:
            continue
        var, g = sorted(got.items())[0]
        lines.append(f"| {tname} | {exp['expected'].get(sp, '')} | {var.split('f_mt_delta_')[-1]} | "
                     f"{float(g['model_delta']):.6f} | {float(g['decimal_delta']):.6f} | {g['mask_reason']} | "
                     f"{float(g['gross_water']):.6e} | {float(g['signed_mass']):.6e} | "
                     f"{g['patches']} ({g['eligible_patches']}) |")
lines += ['', '## Per-probe replay summary', '']
for name, d in agg['probes'].items():
    s = d['summary']
    lines.append(f"- **{name}**: {s['finite_pool_events']} finite-pool events "
                 f"({s['negative_finite_pool_events']} signed negative), "
                 f"max mass replay error {s['max_finite_pool_mass_abs_error']}, "
                 f"max affine slope/offset error {s['max_affine_slope_abs_error']}/{s['max_affine_offset_abs_error']}, "
                 f"daily water/signed book residual {s['max_daily_water_book_residual']}/{s['max_daily_signed_mass_book_residual']}, "
                 f"mask reasons {s['mask_reasons']}")
(outroot / 'summary.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
sys.exit(1 if agg['failures'] else 0)
