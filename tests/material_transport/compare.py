"""Compare every NetCDF variable bit for bit; attributes are not compared."""
import sys
from pathlib import Path
import tempfile
import numpy as np
from netCDF4 import Dataset


def compare(a, b):
    a, b = Path(a), Path(b)
    fa = {f.relative_to(a) for f in a.rglob('*.nc')}
    fb = {f.relative_to(b) for f in b.rglob('*.nc')}
    if not fa or fa != fb:
        return False, f'FILESET a={len(fa)} b={len(fb)} unmatched={len(fa ^ fb)}'
    nvar = 0
    for f in sorted(fa):
        with Dataset(a / f) as x, Dataset(b / f) as y:
            x.set_auto_maskandscale(False)
            y.set_auto_maskandscale(False)
            if set(x.variables) != set(y.variables):
                return False, f'VARIABLESET {f}'
            for name in x.variables:
                u, v = x[name], y[name]
                if u.dimensions != v.dimensions or u.shape != v.shape or u.dtype != v.dtype:
                    return False, f'SHAPE/TYPE {f}:{name}'
                p, q = np.asarray(u[:]), np.asarray(v[:])
                same = np.array_equal(p, q) if p.dtype.hasobject else p.tobytes() == q.tobytes()
                if not same:
                    return False, f'BITS {f}:{name}'
                nvar += 1
    return True, f'IDENTICAL files={len(fa)} vars={nvar} (variable bits; not attributes)'


def self_test():
    with tempfile.TemporaryDirectory() as d:
        a, b = Path(d) / 'a', Path(d) / 'b'
        a.mkdir(); b.mkdir()
        assert not compare(a, b)[0]
        for p in (a, b):
            with Dataset(p / 'x.nc', 'w') as ds:
                ds.createDimension('n', 1)
                ds.createVariable('v', 'f8', ('n',))[:] = 0.
        assert compare(a, b)[0]
        with Dataset(b / 'x.nc', 'a') as ds:
            ds['v'][:] = -0.
        assert not compare(a, b)[0]
        with Dataset(b / 'x.nc', 'a') as ds:
            ds['v'][:] = 0.
            ds.createVariable('extra', 'i4', ('n',))[:] = 1
        assert not compare(a, b)[0]
        (b / 'x.nc').unlink()
        assert not compare(a, b)[0]
    print('compare self-test PASS (empty, equal, signed zero, missing variable/file)')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        ok, message = compare(*sys.argv[1:])
        print(message)
        sys.exit(0 if ok else 1)
