"""Check dependency preparation ordering without compilers, downloads or libraries."""
from pathlib import Path
import os
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    (root / 'tools').mkdir()
    (root / 'tools/netcdf.mk').write_bytes((ROOT / 'tools/netcdf.mk').read_bytes())
    (root / 'tools/prepare_netcdf.py').write_text('''
import sys
from pathlib import Path
assert sys.argv[sys.argv.index('--fc') + 1] == 'fake-fortran'
if '--check-fc' in sys.argv:
    assert sys.argv[sys.argv.index('--check-fc') + 1] == 'cama-fortran'
if '--key' in sys.argv:
    print('test-key')
else:
    import os
    if os.environ.get('FAIL_PROBE'):
        sys.exit('cached probe failed')
    p = Path('.netcdf/test-key/config.mk')
    p.parent.mkdir(parents=True, exist_ok=True)
    if not p.exists():
        p.write_text('NETCDF_F_INCLUDE := /matching/include\\nNETCDF_F_LDFLAGS := -L/matching/lib\\n')
    with Path('checks').open('a') as stream:
        stream.write('probe\\n')
''')
    (root / 'Makefile').write_text('''
NETCDF_AUTO = 1
FC = other-fortran
FF = fake-fortran
NETCDF_INC = /system/include
NETCDF_LIB = /system/lib
LDFLAGS = -L/system/lib -lnetcdff
all:
	@test -f .netcdf/test-key/config.mk
	@echo $(LDFLAGS)
include tools/netcdf.mk
clean:
	@echo clean
''')
    def make(*args):
        return subprocess.check_output(['make', '--no-print-directory', *args], cwd=root, text=True, timeout=10)
    make('-n', 'all')
    make('clean')
    assert not (root / 'checks').exists(), 'dry-run/clean prepared dependencies'
    result = make('all')
    assert '-L/matching/lib -L/system/lib -lnetcdff' in result, result
    count = (root / 'checks').read_text().count('probe')
    assert count == 2, ('fresh make must restart only once', count)
    make('all')
    assert (root / 'checks').read_text().count('probe') == count + 1
    make('all', 'CaMa=YES', 'FCMP=cama-fortran')
    result = make('all', 'LDFLAGS=-L/custom/lib -lnetcdff')
    assert '-L/matching/lib -L/custom/lib -lnetcdff' in result, result
    result = make('all', 'LDFLAGS=-static -L/custom/lib -lnetcdff -lnetcdf')
    assert '-L/matching/lib -static -L/custom/lib -lnetcdff -lnetcdf' in result, result
    assert result.count('-lnetcdff') == 1, result
    failed = subprocess.run(['make', '--no-print-directory', 'all'], cwd=root,
                            env=dict(os.environ, FAIL_PROBE='1'), capture_output=True, text=True)
    assert failed.returncode != 0, failed.stdout
    assert '-L/matching/lib' not in failed.stdout, 'model target ran after failed cache probe'
print('PASS prepare ordering, cache revalidation, library precedence and no dry-run/clean side effects')
