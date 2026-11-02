"""Keep OpenMP on consumers without imposing its storage rules on the host."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    shutil.copy2(ROOT / 'Makefile', root / 'Makefile')
    (root / 'include').mkdir()
    shutil.copy2(ROOT / 'include/Makeoptions.github', root / 'include/Makeoptions')
    for name in ('main', 'share', 'mkinidata', 'mksrfdata', 'postprocess', 'extends', 'tools'):
        (root / name).symlink_to(ROOT / name, target_is_directory=True)
    custom = root / 'custom-mpif90'
    custom.write_text('#!/bin/sh\necho 14\n')
    custom.chmod(0o755)
    for mode, cama in [('GRID', 'CaMaOFF'), ('GRID', 'CaMaON'), ('UNSTRUCTURED', 'CaMaOFF')]:
        subprocess.run(['bash', str(ROOT / '.github/workflows/create_defineh.bash'),
                        mode, 'LULC_IGBP_PFT', 'URBANOFF', 'vanGenu', cama, 'BGCON', 'CROPOFF'],
                       cwd=root, check=True, stdout=subprocess.DEVNULL)
        def commands(target, *options):
            return subprocess.check_output(['make', '-n', target, 'GCC_VERSION=14', *options],
                                           cwd=root, text=True).splitlines()
        lines = commands('colm.x')
        host = next(line for line in lines if '-o .bld/MOD_Precision.o' in line)
        assert ('-fopenmp' in host) == (cama == 'CaMaON'), host
        link = next(line for line in lines if '-o run/colm.x' in line)
        assert ('-fopenmp' in link) == (cama == 'CaMaON'), link
        if cama == 'CaMaON':
            coupled = next(line for line in lines if '-o .bld/cmf_ctrl_maps_mod.o' in line)
            assert '-fopenmp' in coupled, coupled
        lines = commands('post_vector2grid.x')
        for target in ('.bld/MOD_Vector2Grid.o', '.bld/POST_Vector2Grid.o', 'run/post_vector2grid.x'):
            line = next(line for line in lines if '-o ' + target in line)
            assert '-fopenmp' in line, line
        lines = commands('colm.x', f'FF={custom} -fopenmp')
        host = next(line for line in lines if '-o .bld/MOD_Precision.o' in line)
        assert host.startswith(f'{custom} -fopenmp '), host
print('PASS host serial flags, CaMa/Vector2Grid OpenMP and explicit compiler overrides')
