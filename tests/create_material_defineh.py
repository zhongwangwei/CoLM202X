#!/usr/bin/env python3
"""Generate material test configurations on top of the unchanged upstream CI helper."""
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / '.github/workflows/create_defineh.bash'


def configure(args, root):
    if not 7 <= len(args) <= 9:
        raise ValueError('expected seven upstream arguments, optional MTON/MTOFF and LULCCON/LULCCOFF')
    args = list(args) + ['MTOFF', 'LULCCOFF'][len(args) - 7:]
    choices = [('GRID', 'CATCHMENT', 'UNSTRUCTURED', 'SinglePoint'),
               ('LULC_USGS', 'LULC_IGBP', 'LULC_IGBP_PFT', 'LULC_IGBP_PC'),
               ('URBANON', 'URBANOFF'), ('Campbell', 'vanGenu'), ('CaMaON', 'CaMaOFF'),
               ('BGCON', 'BGCOFF'), ('CROPON', 'CROPOFF'), ('MTON', 'MTOFF'),
               ('LULCCON', 'LULCCOFF')]
    for value, allowed in zip(args, choices):
        if value not in allowed:
            raise ValueError(f'{value!r}: expected one of {allowed}')
    if args[7] == 'MTON' and (args[3] != 'vanGenu' or args[0] == 'CATCHMENT' or args[2] == 'URBANON'):
        raise ValueError('material transport requires vanGenu and excludes CATCHMENT/URBANON')
    subprocess.run(['bash', str(GENERATOR), *args[:7]], cwd=root, check=True)
    header = root / 'include/define.h'
    text = header.read_text()
    guard = '#ifndef LULC_IGBP_PFT\n#undef BGC\n#endif'
    if text.count(guard) != 1:
        raise ValueError('upstream BGC guard changed; review material test configuration')
    text = text.replace(guard, '#if !defined(LULC_IGBP_PFT) && !defined(LULC_IGBP_PC)\n#undef BGC\n#endif')
    text += '\n! Material regression test options.\n#undef MaterialTransport\n#undef LULCC\n'
    if args[7] == 'MTON':
        text += '#define MaterialTransport\n'
    if args[8] == 'LULCCON':
        text += '#define LULCC\n'
    header.write_text(text)
    macros = subprocess.check_output(['cc', '-E', '-x', 'c', '-dM', '-'], input=text, text=True)
    defined = {line.split()[1] for line in macros.splitlines() if line.startswith('#define ')}
    bgc = args[5] == 'BGCON' and args[1] in ('LULC_IGBP_PFT', 'LULC_IGBP_PC')
    expected = {'MaterialTransport': args[7] == 'MTON', 'LULCC': args[8] == 'LULCCON',
                'BGC': bgc, 'CROP': bgc and args[6] == 'CROPON',
                'USEMPI': args[0] != 'SinglePoint', args[1]: True,
                {'GRID': 'GRIDBASED'}.get(args[0], args[0]): True}
    for macro, enabled in expected.items():
        if (macro in defined) != enabled:
            raise ValueError(f'{macro}: effective macro differs from requested test configuration')


def self_test():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / 'include').mkdir()
        for grid, land, bgc, crop, mt, lulcc in [
            ('GRID', 'LULC_IGBP_PC', 'BGCON', 'CROPOFF', 'MTON', 'LULCCOFF'),
            ('GRID', 'LULC_IGBP_PFT', 'BGCON', 'CROPOFF', 'MTON', 'LULCCOFF'),
            ('GRID', 'LULC_IGBP', 'BGCOFF', 'CROPOFF', 'MTOFF', 'LULCCOFF'),
            ('GRID', 'LULC_IGBP', 'BGCOFF', 'CROPOFF', 'MTON', 'LULCCON'),
            ('SinglePoint', 'LULC_IGBP', 'BGCOFF', 'CROPOFF', 'MTON', 'LULCCOFF'),
            ('UNSTRUCTURED', 'LULC_IGBP_PC', 'BGCON', 'CROPOFF', 'MTON', 'LULCCOFF'),
            ('GRID', 'LULC_IGBP_PFT', 'BGCON', 'CROPON', 'MTON', 'LULCCOFF'),
        ]:
            configure([grid, land, 'URBANOFF', 'vanGenu', 'CaMaOFF', bgc, crop, mt, lulcc], root)
        base = ['GRID', 'LULC_IGBP_PC', 'URBANOFF', 'vanGenu', 'CaMaOFF', 'BGCON', 'CROPOFF']
        configure(base, root)
        configure(base + ['MTON'], root)
        for args in [base + ['INVALID'], base[:3] + ['Campbell'] + base[4:] + ['MTON']]:
            try:
                configure(args, root)
            except ValueError:
                continue
            raise AssertionError('invalid configuration accepted')
    print('PASS 9 effective macro configurations and invalid configuration rejection')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        configure(sys.argv[1:], Path.cwd())
