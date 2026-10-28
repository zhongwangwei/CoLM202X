#!/usr/bin/env python
"""Isotope sanity check of material history files.

Every isotope variable is turned into delta (permil vs VSMOW): f_mt_delta_*
are already delta, f_mt_conc_* / f_mt_levdelta_* carry R (or delta, read
from units).  For each variable: valid count, NaN/Inf, percentiles and the
points outside a generous plausible range, with the host water amount at
that point (so a near-empty pool can be told from a real error).

usage: iso_check.py CASEDIR [CASEDIR ...]
"""
import sys, glob, os
import numpy as np
import netCDF4 as nc

RSTD = {'H2_18O': 2.0052e-3, 'HDO': 1.5576e-4}
RANGE = {'H2_18O': (-60., 40.), 'HDO': (-450., 250.)}
# host water that goes with each material pool / flux
PAIR = [('soisno', 'f_wliq_soisno'), ('ldew', 'f_ldew'), ('snowpack', 'f_scv'),
        ('scv', 'f_scv'), ('wdsrf', 'f_wdsrf'), ('wetwat', 'f_wetwat'),
        ('wa', 'f_wa'), ('precip', 'f_xy_prc'), ('runoff', 'f_rnof'),
        ('soilevap', 'f_fevpg'), ('canopyevap', 'f_fevpl'), ('transp', 'f_etr'),
        ('leaf', 'f_etr'), ('evap', 'f_fevpa'), ('subl', 'f_fevpa')]


# material masses and fluxes, not ratios: no delta
NOT_RATIO = ('f_mt_flux_', 'f_mt_evap_mass_', 'f_mt_vapor_exchange_')


def species(name):
    if name.startswith(NOT_RATIO):
        return None
    for s in RSTD:
        if name.endswith('_' + s):
            return s
    return None


def to_delta(v, var, sp):
    u = getattr(var, 'units', '')
    if 'permil' in u.lower() or 'delta' in var.name or 'levdelta' in var.name:
        return v
    return (v / RSTD[sp] - 1.) * 1000.


def host_of(name):
    for key, h in PAIR:
        if key in name:
            return h
    return None


def check(case):
    mfiles = sorted(glob.glob(os.path.join(case, 'history', '*hist_material*.nc')))
    rfiles = sorted(glob.glob(os.path.join(case, 'history', '*unitcat*.nc')))
    hfiles = {os.path.basename(f).replace('_material', ''): f for f in mfiles}
    print('#### %s: %d material files, %d river files' % (case, len(mfiles), len(rfiles)))
    stats = {}
    for f in mfiles + rfiles:
        ds = nc.Dataset(f)
        hostf = os.path.join(case, 'history', os.path.basename(f).replace('_material', ''))
        hs = nc.Dataset(hostf) if (f in mfiles and os.path.exists(hostf) and hostf != f) else None
        for name, var in ds.variables.items():
            sp = species(name)
            if sp is None or var.ndim < 2:
                continue
            a = var[:]
            data = np.ma.filled(a.astype('f8'), np.nan)
            if hasattr(a, 'mask'):
                data[np.ma.getmaskarray(a)] = np.nan
            d = to_delta(data, var, sp)
            st = stats.setdefault(name, {'n': 0, 'nonfinite': 0, 'vals': {}, 'bad': [], 'sp': sp,
                                         'units': getattr(var, 'units', '')})
            fin = np.isfinite(d)
            st['n'] += int(fin.sum())
            st['nonfinite'] += int((~np.ma.getmaskarray(a) & ~fin).sum())
            st['vals'][f] = (np.flatnonzero(fin), d[fin], d.shape)
            lo, hi = RANGE[sp]
            badmask = fin & ((d < lo) | (d > hi))
            if badmask.any():
                idx = np.argwhere(badmask)
                h = host_of(name)
                hv = None
                if hs is not None and h in hs.variables:
                    hv = np.ma.filled(hs.variables[h][:].astype('f8'), np.nan)
                order = np.argsort(-np.abs(d[badmask]))[:8]
                for k in order:
                    ix = tuple(idx[k])
                    w = None
                    if hv is not None and hv.shape == d.shape:
                        w = hv[ix]
                    elif hv is not None and hv.ndim == d.ndim - 0 and hv.shape[-2:] == d.shape[-2:]:
                        try:
                            w = hv[ix]
                        except Exception:
                            w = None
                    st['bad'].append((abs(d[ix]), d[ix], os.path.basename(f), ix, w, h))
        ds.close()
        if hs is not None:
            hs.close()
    print('%-38s %-6s %8s %5s %9s %9s %9s %9s %9s %9s %9s %6s' %
          ('variable', 'sp', 'n', 'nonf', 'min', 'p0.1', 'p1', 'p50', 'p99', 'p99.9', 'max', 'nout'))
    for name in sorted(stats):
        st = stats[name]
        v = np.concatenate([a[1] for a in st['vals'].values()]) if st['vals'] else np.array([])
        v = v[np.isfinite(v)]
        if v.size == 0:
            print('%-38s %-6s %8d %5d   (no valid value)' % (name, st['sp'], 0, st['nonfinite']))
            continue
        lo, hi = RANGE[st['sp']]
        nout = int(((v < lo) | (v > hi)).sum())
        q = np.percentile(v, [0.1, 1, 50, 99, 99.9])
        print('%-38s %-6s %8d %5d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %6d' %
              (name, st['sp'], v.size, st['nonfinite'], v.min(), q[0], q[1], q[2], q[3], q[4], v.max(), nout))
    print('-- worst points outside the plausible range (delta, file, index, host water)')
    for name in sorted(stats):
        bad = sorted(stats[name]['bad'], key=lambda x: -x[0])[:5]
        for b in bad:
            print('   %-36s %10.2f  %s %s  %s=%s' % (name, b[1], b[2], b[3], b[5], b[4]))
    # d-excess where both species exist
    print('-- d-excess (dD - 8 d18O) of the same variable')
    paired = {}
    for name in sorted(stats):
        if not name.endswith('_H2_18O'):
            continue
        other = name[:-len('H2_18O')] + 'HDO'
        if other not in stats:
            continue
        chunks = []
        for f in sorted(stats[name]['vals'].keys() & stats[other]['vals'].keys()):
            ia, a, shape_a = stats[name]['vals'][f]
            ib, b, shape_b = stats[other]['vals'][f]
            assert shape_a == shape_b, (f, name, other)
            _, pa, pb = np.intersect1d(ia, ib, assume_unique=True, return_indices=True)
            chunks.append(b[pb] - 8. * a[pa])
        dx = np.concatenate(chunks) if chunks else np.array([])
        if dx.size == 0:
            continue
        paired[name] = dx
        q = np.percentile(dx, [0.1, 1, 50, 99, 99.9])
        print('   %-36s p0.1 %8.2f p1 %8.2f p50 %8.2f p99 %8.2f p99.9 %8.2f' % (name[:-7], *q))

    assert not any(st['nonfinite'] for st in stats.values()), 'unmasked NaN/Inf in isotope output'
    assert stats, 'no isotope history variables found'
    return stats, paired


def self_test():
    import tempfile, contextlib, io
    with tempfile.TemporaryDirectory() as root:
        os.mkdir(os.path.join(root, 'history'))
        f = os.path.join(root, 'history', 'case_hist_material.nc')
        with nc.Dataset(f, 'w') as ds:
            ds.createDimension('time', 1); ds.createDimension('patch', 3)
            a = ds.createVariable('f_mt_delta_precip_H2_18O', 'f8', ('time', 'patch'), fill_value=-999.)
            b = ds.createVariable('f_mt_delta_precip_HDO', 'f8', ('time', 'patch'), fill_value=-999.)
            a[:] = [[-999., 2., 3.]]; b[:] = [[10., -999., 30.]]
        with contextlib.redirect_stdout(io.StringIO()):
            stats, pairs = check(root)
        assert list(pairs.values())[0].tolist() == [6.], pairs
        with nc.Dataset(f, 'a') as ds:
            ds['f_mt_delta_precip_H2_18O'][0, 0] = np.nan
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                check(root)
        except AssertionError as e:
            assert 'NaN/Inf' in str(e), e
        else:
            raise AssertionError('unmasked NaN was not detected')
    print('iso_check self-test PASS (joint masks, unmasked NaN)')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        for c in sys.argv[1:]:
            check(c)
