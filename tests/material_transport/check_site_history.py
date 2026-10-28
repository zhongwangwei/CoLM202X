#!/usr/bin/env python3
"""Check daily SinglePoint host/material history; optionally compare a reference.

Requires the same netCDF4/numpy dependencies as the integration report tools.
Example: check_site_history.py CASE/history --days 33 --reference DIRECT/history
"""
import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def read_streams(directory):
    streams = {'host': [], 'material': []}
    for path in sorted(directory.glob('*.nc')):
        stream = 'material' if '_hist_material_' in path.name else 'host'
        with Dataset(path) as data:
            assert 'time' in data.variables, f'{path}: missing time'
            times = np.asarray(data['time'][:])
            assert len(times) and np.all(np.diff(times) == 1440), path
            variables = {}
            for name, variable in data.variables.items():
                values = np.ma.asarray(variable[:])
                if values.dtype.kind in 'fi':
                    assert np.all(np.isfinite(values.compressed())), (path, name)
                variables[name] = values
            streams[stream].append((path.name.split('_hist_', 1)[-1], times, variables))
    return streams


def check(directory, days, reference=None):
    streams = read_streams(directory)
    times = {}
    for stream, files in streams.items():
        assert files, f'No {stream} history in {directory}'
        times[stream] = np.concatenate([item[1] for item in files])
        assert len(times[stream]) == days, (stream, len(times[stream]), days)
        assert np.all(np.diff(times[stream]) == 1440), f'{stream}: missing/duplicate daily time'
    np.testing.assert_array_equal(times['host'], times['material'])
    if reference:
        other = read_streams(reference)
        for stream, files in streams.items():
            assert len(files) == len(other[stream]), stream
            for (name, _, variables), (refname, _, refvars) in zip(files, other[stream]):
                assert name == refname and variables.keys() == refvars.keys(), name
                for varname, values in variables.items():
                    ref = refvars[varname]
                    np.testing.assert_array_equal(np.ma.getmaskarray(values), np.ma.getmaskarray(ref), err_msg=varname)
                    np.testing.assert_array_equal(values.compressed(), ref.compressed(), err_msg=f'{name}:{varname}')
    print(f'PASS: {days} daily records per stream, matching times, finite values'
          + ('; all fields equal reference' if reference else ''))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--days', type=int, required=True)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    check(args.directory, args.days, args.reference)
