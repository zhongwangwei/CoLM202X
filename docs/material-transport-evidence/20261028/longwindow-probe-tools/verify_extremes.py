#!/usr/bin/env python
"""Reproduce the reported p9s extremes directly from the baseline history.

For every selected target this reads the exact baseline file/variable/index that
the p9s audit used and asserts that the stored value equals the reported delta.
This makes the probe target selection traceable to the actual baseline output,
not to a copied number.

usage: verify_extremes.py P9S_RESULTS_ROOT selection.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def main(p9s_root, selection_path):
    p9s_root = Path(p9s_root)
    sel = json.loads(Path(selection_path).read_text())
    failures = []
    checked = []
    for target in sel['chosen']:
        src = target['source']
        case = src['case']
        path = p9s_root / 'runs' / case / 'case' / 'history' / src['file']
        if not path.exists():
            failures.append(f"{target['name']}: missing {path}")
            continue
        with Dataset(path) as ds:
            if src['var'] not in ds.variables:
                failures.append(f"{target['name']}: variable {src['var']} absent")
                continue
            var = ds[src['var']]
            if len(src['index']) != len(var.shape):
                failures.append(f"{target['name']}: index rank {len(src['index'])} != {len(var.shape)}")
                continue
            a = np.asarray(np.ma.filled(np.ma.asarray(var[:]).astype('f8'), np.nan))
            value = float(a[tuple(src['index'])])
        ok = bool(np.isclose(value, float(src['delta']), rtol=0.0, atol=1e-12))
        checked.append(dict(target=target['name'], case=case, file=src['file'],
                            var=src['var'], index=src['index'],
                            reported=float(src['delta']), stored=value, equal=ok))
        if not ok:
            failures.append(f"{target['name']}: reported {src['delta']} != stored {value}")
    out = dict(checked=len(checked), failures=failures, detail=checked)
    print(json.dumps(out, indent=1))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
