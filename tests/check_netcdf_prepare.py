#!/usr/bin/env python3
"""Check that dependency fallback and cache validation fail closed."""
import argparse
import importlib.util
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prepare_netcdf', ROOT / 'tools/prepare_netcdf.py')
netcdf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(netcdf)


class PreparationTests(unittest.TestCase):
    def test_checksum_failure_removes_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / 'source.tar.gz'
            with patch.object(netcdf.urllib.request, 'urlopen', return_value=io.BytesIO(b'wrong source')):
                with self.assertRaisesRegex(RuntimeError, 'checksum mismatch'):
                    netcdf.download('4.6.4', archive)
            self.assertFalse(archive.exists())

    def test_failed_upstream_check_prevents_install_and_probe(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            commands = []

            def run(argv, cwd, log, env=None):
                commands.append(argv)
                if argv == ['make', '-j2', 'check']:
                    raise subprocess.CalledProcessError(2, argv)

            def output(argv):
                return {'--version': 'netCDF 4.9.3', '--cflags': '-I/system/include',
                        '--libdir': '/system/lib', '--libs': '-lnetcdf'}[argv[-1]]

            with patch.object(netcdf, 'download'), patch.object(netcdf, 'run', side_effect=run), patch.object(netcdf, 'output', side_effect=output), patch.object(netcdf, 'probe') as probe:
                with self.assertRaises(subprocess.CalledProcessError):
                    netcdf.build(['gfortran'], directory)
                self.assertIn(['make', '-j2', 'check'], commands)
                self.assertNotIn(['make', 'install'], commands)
                probe.assert_not_called()

    def test_system_success_never_builds_and_cached_selection_is_rechecked(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(cache_dir=tmp, include='/system/include', lib='/system/lib')
            with patch.object(netcdf, 'output', return_value='-lnetcdf'), patch.object(netcdf, 'probe') as probe, patch.object(netcdf, 'build') as build:
                config = netcdf.prepare(args, 'key', ['gfortran'])
                stamp = config.stat().st_mtime_ns
                netcdf.prepare(args, 'key', ['gfortran'])
                self.assertEqual(probe.call_count, 2)
                build.assert_not_called()
                self.assertEqual(config.stat().st_mtime_ns, stamp)
                self.assertIn('NETCDF_F_LDFLAGS := \n', config.read_text())

    def test_cama_probe_failure_blocks_new_and_cached_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(cache_dir=tmp, include='/system/include', lib='/system/lib',
                                      check_fc='cama-fortran -fopenmp')
            error = subprocess.CalledProcessError(1, ['cama-fortran'])
            directory = Path(tmp) / 'key'
            with patch.object(netcdf, 'output', return_value='-lnetcdf'), patch.object(netcdf, 'resolve_compiler', return_value=['cama-fortran', '-fopenmp']), patch.object(netcdf, 'probe', side_effect=[None, error]), patch.object(netcdf, 'build') as build:
                with self.assertRaisesRegex(RuntimeError, 'CaMa compiler cannot use selected NetCDF'):
                    netcdf.prepare(args, 'key', ['gfortran'])
                build.assert_not_called()
            self.assertFalse((directory / 'config.mk').exists())
            self.assertFalse((directory / 'selection.json').exists())
            with patch.object(netcdf, 'output', return_value='-lnetcdf'), patch.object(netcdf, 'resolve_compiler', return_value=['cama-fortran']), patch.object(netcdf, 'probe'):
                netcdf.prepare(args, 'key', ['gfortran'])
            with patch.object(netcdf, 'resolve_compiler', return_value=['cama-fortran']), patch.object(netcdf, 'probe', side_effect=[None, error]), patch.object(netcdf, 'build') as build:
                with self.assertRaisesRegex(RuntimeError, 'align FF/FCMP/FC'):
                    netcdf.prepare(args, 'key', ['gfortran'])
                build.assert_not_called()
            self.assertFalse((directory / 'config.mk').exists())

    def test_cama_on_off_reuses_primary_compiler_cache_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(fc='gfortran', cache_dir=tmp, include=tmp, lib=tmp,
                                      check_fc='cama-fortran -fopenmp')
            with patch.object(netcdf, 'output', return_value='version'), patch.object(netcdf.shutil, 'which', side_effect=lambda command: '/compiler/' + command):
                first, _ = netcdf.identity(args)
                args.check_fc = None
                second, _ = netcdf.identity(args)
                self.assertEqual(first, second)

    def test_failed_system_probe_builds_and_publishes_only_after_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(cache_dir=tmp, include='/system/include', lib='/system/lib')
            error = subprocess.CalledProcessError(1, ['compiler'])
            with patch.object(netcdf, 'output', return_value='-lnetcdf'), patch.object(netcdf, 'probe', side_effect=error), patch.object(netcdf, 'build', side_effect=error):
                with self.assertRaises(subprocess.CalledProcessError):
                    netcdf.prepare(args, 'key', ['gfortran'])
            self.assertFalse((Path(tmp) / 'key/config.mk').exists())
            self.assertFalse((Path(tmp) / 'key/selection.json').exists())
            with patch.object(netcdf, 'output', return_value='-lnetcdf'), patch.object(netcdf, 'probe', side_effect=error), patch.object(netcdf, 'build', return_value=(Path('/private/include'), ['-L/private/lib', '-lnetcdff'])) as build:
                config = netcdf.prepare(args, 'key', ['gfortran'])
                build.assert_called_once()
                self.assertIn('NETCDF_F_INCLUDE := /private/include', config.read_text())
                self.assertIn('NETCDF_F_LDFLAGS := -L/private/lib\n', config.read_text())
                self.assertNotIn('-lnetcdff', config.read_text())
            with patch.object(netcdf, 'probe', side_effect=error), patch.object(netcdf, 'build') as build:
                with self.assertRaisesRegex(RuntimeError, 'cached NetCDF failed'):
                    netcdf.prepare(args, 'key', ['gfortran'])
                build.assert_not_called()
            self.assertFalse(config.exists())


if __name__ == '__main__':
    unittest.main()
