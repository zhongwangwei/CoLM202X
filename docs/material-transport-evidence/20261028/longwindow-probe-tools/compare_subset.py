#!/usr/bin/env python
"""Compare a probe segment's history against the overlapping part of the baseline.

The probe deliberately runs fewer days than the baseline month, so the full-file
comparator cannot be used.  This comparator still requires:

* every probe history file has a baseline counterpart with the same name;
* identical dimension sets (time may be longer in the baseline) and, except for
  the time dimension, identical shapes;
* identical variable sets and dtypes;
* identical time coordinate values for the overlapping records;
* raw variable bytes for every overlapping time record, and every non-time
  variable in full.

usage: compare_subset.py PROBE_HIST_DIR BASELINE_HIST_DIR
Prints one JSON object; exits nonzero on any mismatch.
"""
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def raw_equal(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    if a.dtype.hasobject or b.dtype.hasobject:
        return np.array_equal(a, b)
    return a.tobytes() == b.tobytes()


def compare(probe, base):
    pf = {p.name: p for p in Path(probe).glob('*.nc')}
    if not pf:
        raise AssertionError(f'probe history is empty: {probe}')
    bf = {p.name: p for p in Path(base).glob('*.nc')}
    missing = sorted(set(pf) - set(bf))
    assert not missing, f'missing baseline history files: {missing}'
    result = {}
    for name in sorted(pf):
        with Dataset(pf[name]) as x, Dataset(bf[name]) as y:
            x.set_auto_maskandscale(False)
            y.set_auto_maskandscale(False)
            assert set(x.dimensions) == set(y.dimensions), f'{name}: DIMSET'
            for k in x.dimensions:
                if k == 'time':
                    continue
                assert len(x.dimensions[k]) == len(y.dimensions[k]), f'{name}: DIM {k}'
            n = len(x.dimensions['time'])
            m = len(y.dimensions['time'])
            assert n <= m, f'{name}: probe has more records ({n}) than baseline ({m})'
            assert set(x.variables) == set(y.variables), f'{name}: VARSET'
            nvar = 0
            for k, u in x.variables.items():
                v = y[k]
                assert u.dimensions == v.dimensions, f'{name}:{k}: DIMENSIONS'
                if 'time' in u.dimensions:
                    assert u.shape[1:] == v.shape[1:], f'{name}:{k}: SHAPE'
                    assert u.shape[0] == n, f'{name}:{k}: probe time length'
                    assert u.dtype == v.dtype, f'{name}:{k}: DTYPE'
                    if k == 'time':
                        assert raw_equal(u[...], v[:n]), f'{name}:{k}: TIME VALUES'
                    else:
                        assert raw_equal(u[...], v[:n]), f'{name}:{k}: BITS'
                else:
                    assert u.shape == v.shape, f'{name}:{k}: SHAPE'
                    assert u.dtype == v.dtype, f'{name}:{k}: DTYPE'
                    assert raw_equal(u[...], v[...]), f'{name}:{k}: BITS'
                nvar += 1
            result[name] = dict(records=n, baseline_records=m, variables=nvar)
    out = dict(probe=str(probe), baseline=str(base),
               files=len(result), variables=sum(v['variables'] for v in result.values()),
               records={k: v['records'] for k, v in result.items()},
               detail=result)
    return out


if __name__ == '__main__':
    probe, base = sys.argv[1], sys.argv[2]
    print(json.dumps(compare(probe, base), indent=1))
