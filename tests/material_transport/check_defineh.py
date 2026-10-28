#!/usr/bin/env python3
"""Check effective CI configuration macros without changing the working header."""
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
GENERATOR = ROOT / '.github/workflows/create_defineh.bash'


def check(args):
    with tempfile.TemporaryDirectory() as directory:
        header = pathlib.Path(directory) / 'include/define.h'
        header.parent.mkdir()
        subprocess.run(['bash', str(GENERATOR), *args], cwd=directory,
                       check=True, capture_output=True, text=True)
        output = subprocess.run(['cc', '-E', '-dM', '-x', 'c', str(header)], check=True,
                                capture_output=True, text=True).stdout
    macros = {line.split()[1] for line in output.splitlines()
              if line.startswith('#define ')}
    bgc = args[5] == 'BGCON' and args[1] in ('LULC_IGBP_PFT', 'LULC_IGBP_PC')
    expected = {
        args[1]: True,
        'BGC': bgc,
        'CROP': bgc and args[6] == 'CROPON',
        'MaterialTransport': len(args) >= 8 and args[7] == 'MTON',
        'LULCC': len(args) >= 9 and args[8] == 'LULCCON',
    }
    for macro, enabled in expected.items():
        assert (macro in macros) == enabled, (args, macro, enabled)


def main():
    base = ['GRID', 'LULC_IGBP_PFT', 'URBANOFF', 'vanGenu',
            'CaMaOFF', 'BGCOFF', 'CROPOFF']
    for suffix in ([], ['MTON'], ['MTOFF'], ['MTON', 'LULCCOFF'],
                   ['MTON', 'LULCCON'], ['MTOFF', 'LULCCON']):
        check(base + suffix)
    with tempfile.TemporaryDirectory() as directory:
        pathlib.Path(directory, 'include').mkdir()
        result = subprocess.run(['bash', str(GENERATOR), *base, 'MTON', 'INVALID'],
                                cwd=directory, capture_output=True, text=True)
        assert result.returncode != 0, 'Invalid LULCC argument was accepted'
        assert 'Error in argument 9' in result.stdout
    for land in ('LULC_IGBP_PC', 'LULC_IGBP_PFT', 'LULC_IGBP'):
        args = base.copy()
        args[1], args[5] = land, 'BGCON'
        check(args + ['MTON'])

    cases = [line.split() for line in
             (ROOT / '.github/workflows/TestCaseLists').read_text().splitlines()
             if line.strip()]
    names = [case[0] for case in cases]
    assert len(names) == len(set(names)), 'Duplicate test case names'
    coverage = set()
    for name, *args in cases:
        assert 7 <= len(args) <= 9, (name, args)
        assert not (args[5] == 'BGCON' and args[-1] == 'LULCCON'), name
        check(args)
        if len(args) >= 8 and args[7] == 'MTON':
            coverage.add((args[1], args[5], args[8] if len(args) == 9 else 'LULCCOFF'))
    assert ('LULC_IGBP_PC', 'BGCON', 'LULCCOFF') in coverage
    for land in ('LULC_IGBP', 'LULC_IGBP_PFT'):
        assert (land, 'BGCOFF', 'LULCCON') in coverage
    print(f'PASS: optional arguments and effective macros for {len(cases)} unique CI cases')


if __name__ == '__main__':
    main()
