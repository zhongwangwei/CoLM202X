#!/usr/bin/env python3
"""Check every forced-isotope checkpoint block, with an explicit cold-start date.

The exception is only for host-only or committed empty-material files at the declared initial checkpoint.
It never permits a missing cache in a later checkpoint or in an initialized
material restart. Requires NumPy and netCDF4; --self-test uses tiny fixtures.
"""
import argparse
import json
import re
import tempfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

DATE = re.compile(r'^\d{4}-\d{3}-\d{5}$')
FIELDS = ('mt_forcing_cache_schema', 'mt_forcing_cache_count',
          'mt_forcing_cache_identity', 'mt_forcing_precip_last',
          'mt_forcing_vapor_last', 'mt_forcing_precip_reliable')


def values(var):
    a = np.ma.asarray(var[:])
    assert not np.ma.getmaskarray(a).any(), f'{var.name}: missing values'
    a = np.asarray(a)
    assert np.isfinite(a).all(), f'{var.name}: non-finite values'
    return a


def check(root, initial_date=None):
    root = Path(root)
    dirs = sorted(p for p in root.iterdir() if p.is_dir())
    assert all(DATE.fullmatch(p.name) or p.name in ('const', 'ParaOpt') for p in dirs), 'unknown restart directory'
    dates = [p for p in dirs if DATE.fullmatch(p.name)]
    assert dates, 'no dated checkpoints'
    if initial_date is not None:
        assert dates[0].name == initial_date, 'declared initial checkpoint is not the first checkpoint'
    report = []
    for date in dates:
        files = sorted(date.glob(f'*_restart_{date.name}_lc*.nc'))
        assert files, f'{date}: no land restart blocks'
        cold = 0
        blocks = []
        for f in files:
            with Dataset(f) as ds:
                present = {n for n in ds.variables if n.startswith('mt_forcing_')}
                if not present:
                    assert date.name == initial_date, f'{f}: missing cache after cold start'
                    material = {n for n in ds.variables if n.startswith('mt_')}
                    if material:
                        expected = {'mt_land_restart_complete': 1, 'mt_land_restart_schema': 5,
                                    'mt_land_transport_count': 0}
                        assert material == set(expected) | {'mt_aquifer_mixing_water_mm'}, f'{f}: initialized material payload without cache'
                        for name, value in expected.items():
                            var = ds[name]
                            assert var.shape == () and var.dtype.kind in 'iu', f'{f}: invalid empty transaction metadata'
                            assert int(values(var)) == value, f'{f}: invalid empty transaction {name}'
                        assert ds['mt_aquifer_mixing_water_mm'].shape == (), f'{f}: non-scalar aquifer metadata'
                        values(ds['mt_aquifer_mixing_water_mm'])
                    cold += 1
                    continue
                assert all(n in ds.variables for n in FIELDS), f'{f}: incomplete forcing cache'
                schema, count = (ds[n] for n in FIELDS[:2])
                assert schema.shape == count.shape == (), f'{f}: schema/count must be scalar'
                assert schema.dtype.kind in 'iu' and count.dtype.kind in 'iu', f'{f}: schema/count type'
                assert int(values(schema)) == 2, f'{f}: expected schema 2'
                nvars = int(values(count))
                assert nvars > 0, f'{f}: forced case has no forcing variables'
                identity = ds[FIELDS[2]]
                assert identity.dimensions == ('mt_forcing_id_var', 'mt_forcing_id_field'), f'{f}: identity dimensions'
                assert identity.shape[0] == nvars + 1 and identity.shape[1] > 0, f'{f}: identity shape'
                values(identity)
                arrays = []
                for name in FIELDS[3:]:
                    var = ds[name]
                    assert var.dimensions == ('patch', 'mt_forcing_species'), f'{f}: {name} dimensions'
                    arrays.append(values(var))
                rain, vapor, reliable = arrays
                assert rain.shape == vapor.shape == reliable.shape and all(rain.shape), f'{f}: cache shape'
                assert (rain > 0).all() and (vapor > 0).all(), f'{f}: invalid isotope ratio'
                assert np.isin(reliable, [0., 1.]).all(), f'{f}: invalid reliability flag'
                blocks.append(dict(file=f.name, patches=rain.shape[0], species=rain.shape[1],
                                   forcing_variables=nvars, reliable=int(reliable.sum())))
        assert not (cold and blocks), f'{date}: mixed cold/initialized land blocks'
        report.append(dict(date=date.name, cold_host_blocks=cold, cache_blocks=blocks))
    assert any(r['cache_blocks'] for r in report), 'no initialized forced checkpoints'
    return report


def self_test():
    tested = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        initial, later = '2005-359-00000', '2005-360-00000'
        for date in (initial, later):
            (root/date).mkdir()
        def fixture(date, cache, suffix='e125_n45'):
            f = root/date/f'case_restart_{date}_lc2005_{suffix}.nc'
            with Dataset(f, 'w') as ds:
                ds.createDimension('patch', 2)
                ds.createVariable('wdsrf', 'f8', ('patch',))[:] = 1.
                if cache:
                    for n, v in zip(FIELDS[:2], (2, 6)):
                        ds.createVariable(n, 'i4')[:] = v
                    ds.createDimension('mt_forcing_species', 2)
                    ds.createDimension('mt_forcing_id_var', 7)
                    ds.createDimension('mt_forcing_id_field', 3)
                    ds.createVariable(FIELDS[2], 'i4', ('mt_forcing_id_var', 'mt_forcing_id_field'))[:] = 1
                    for n in FIELDS[3:]:
                        ds.createVariable(n, 'f8', ('patch', 'mt_forcing_species'))[:] = 1. if n == FIELDS[-1] else .002
            return f
        fixture(initial, False)
        warm = fixture(later, True)
        assert check(root, initial)[0]['cold_host_blocks'] == 1
        tested.append('explicit host-only initial accepted')
        def reject(name, initial_arg=initial):
            try:
                check(root, initial_arg)
            except AssertionError:
                tested.append(name)
            else:
                raise AssertionError(f'negative fixture accepted: {name}')
        reject('initial requires explicit exception', None)
        reject('wrong initial date', later)
        fixture(later, False)
        reject('later missing cache')
        fixture(later, True)
        bad = fixture(later, False, 'e130_n45')
        reject('one later block missing cache')
        bad.unlink()
        for field, bad_value, name in ((FIELDS[0], 1, 'schema 1'), (FIELDS[1], 0, 'empty forced cache'),
                                      (FIELDS[-1], .5, 'invalid flag'), (FIELDS[3], np.nan, 'NaN ratio'),
                                      (FIELDS[4], 0., 'zero vapor ratio')):
            fixture(later, True)
            with Dataset(warm, 'a') as ds:
                ds[field][:] = bad_value
            reject(name)
        fixture(later, True)
        with Dataset(warm, 'a') as ds:
            ds.renameVariable(FIELDS[-1], 'removed_flag')
        reject('incomplete cache')
        fixture(later, True)
        with Dataset(root/initial/f'case_restart_{initial}_lc2005_e125_n45.nc', 'a') as ds:
            ds.createVariable('mt_land_restart_complete', 'i4')[:] = 1
        reject('incomplete initial material transaction')
        initial_file = root/initial/f'case_restart_{initial}_lc2005_e125_n45.nc'
        with Dataset(initial_file, 'a') as ds:
            ds.createVariable('mt_land_restart_schema', 'i4')[:] = 5
            ds.createVariable('mt_land_transport_count', 'i4')[:] = 0
            ds.createVariable('mt_aquifer_mixing_water_mm', 'f8')[:] = -1.
        assert check(root, initial)[0]['cold_host_blocks'] == 1
        tested.append('committed empty initial material transaction accepted')
        for field, bad_value in [('mt_land_restart_complete', 0), ('mt_land_transport_count', 2),
                                 ('mt_land_restart_schema', 4)]:
            with Dataset(initial_file, 'a') as ds:
                good = ds[field][:].copy()
                ds[field][:] = bad_value
            reject('initial invalid ' + field)
            with Dataset(initial_file, 'a') as ds:
                ds[field][:] = good
        with Dataset(initial_file, 'a') as ds:
            ds.createVariable('mt_wliq_soisno', 'f8', ('patch',))[:] = 1.
        reject('initial actual material pool without cache')
    print('PASS forcing cache self-test:', len(tested), 'positive/negative behaviors')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('restart', nargs='?')
    parser.add_argument('--initial-date')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        parser.error('restart directory required') if args.restart is None else None
        print(json.dumps(check(args.restart, args.initial_date), indent=2))
