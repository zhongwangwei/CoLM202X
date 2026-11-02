#!/usr/bin/env python3
"""Verify system NetCDF-Fortran, rebuilding it for FC only when incompatible."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import urllib.request

RELEASES = {
    '4.6.4': '98159c1e0f63b3b59bb5eda12f2d80126f5b1aad93032d1490989a5752e0df99',
    '4.5.4': '0a19b26a2b6e29fab5d29d7d7e08c24e87712d09a5cafeea90e16e0a2ab86b81',
}
PROBE = '''program probe
  use netcdf
  implicit none
  integer :: ncid, dimid, varid, scalarid
  real(kind=8) :: values(2), scalar
  call checked(nf90_create('probe.nc', ior(nf90_clobber, nf90_netcdf4), ncid))
  call checked(nf90_def_dim(ncid, 'sample', 2, dimid))
  call checked(nf90_def_var(ncid, 'values', nf90_double, [dimid], varid))
  call checked(nf90_def_var(ncid, 'scalar', nf90_double, scalarid))
  call checked(nf90_enddef(ncid))
  call checked(nf90_put_var(ncid, varid, [42.5d0, -1.25d0]))
  call checked(nf90_put_var(ncid, scalarid, 3.5d0))
  call checked(nf90_close(ncid))
  call checked(nf90_open('probe.nc', nf90_nowrite, ncid))
  call checked(nf90_inq_varid(ncid, 'values', varid))
  call checked(nf90_inq_varid(ncid, 'scalar', scalarid))
  call checked(nf90_get_var(ncid, varid, values))
  call checked(nf90_get_var(ncid, scalarid, scalar))
  call checked(nf90_close(ncid))
  if (any(values /= [42.5d0, -1.25d0]) .or. scalar /= 3.5d0) stop 2
contains
  subroutine checked(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, nf90_strerror(status)
      stop 1
    endif
  end subroutine
end program
'''


def output(argv):
    return subprocess.check_output(argv, text=True, stderr=subprocess.STDOUT).strip()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def resolve_compiler(command):
    compiler = shlex.split(command)
    if not compiler or not shutil.which(compiler[0]):
        raise RuntimeError(f'Fortran compiler not found: {command}')
    compiler[0] = str(Path(shutil.which(compiler[0])).absolute())
    return compiler


def identity(args):
    compiler = resolve_compiler(args.fc)
    include, lib = Path(args.include).resolve(), Path(args.lib).resolve()
    libs = {str(p.resolve()): digest(p) for p in sorted(lib.glob('libnetcdff.*')) if p.is_file()}
    nc = output(['nc-config', '--all'])
    info = {'compiler': compiler, 'version': output(compiler + ['--version']),
            'binary': str(Path(compiler[0]).resolve()),
            'include': str(include), 'module': digest(include / 'netcdf.mod'),
            'lib': str(lib), 'libraries': libs, 'nc_config': nc,
            'script': digest(Path(__file__)), 'cache_dir': str(Path(args.cache_dir).resolve()),
            'backends': {key: os.environ.get(key) for key in
                         ['OMPI_FC', 'MPICH_FC', 'I_MPI_F90', 'I_MPI_FC']}}
    return hashlib.sha256(json.dumps(info, sort_keys=True).encode()).hexdigest(), compiler


def run(argv, cwd, log, env=None):
    with log.open('a') as stream:
        stream.write('$ ' + shlex.join(map(str, argv)) + '\n')
        stream.flush()
        result = subprocess.run(argv, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        print(f'Failed command; log: {log}\n' + log.read_text(errors='replace')[-16000:], file=sys.stderr)
        if 'check' in argv:
            for report in sorted(Path(cwd).rglob('test-suite.log')):
                print(f'\n{report}:\n' + report.read_text(errors='replace')[-16000:], file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, argv)


def probe(compiler, include, libraries, directory, log):
    with tempfile.TemporaryDirectory(prefix='probe-', dir=directory) as tmp:
        work = Path(tmp)
        (work / 'probe.f90').write_text(PROBE)
        run(compiler + ['-I' + str(include), 'probe.f90', '-o', 'probe'] + libraries, work, log)
        run([str(work / 'probe')], work, log)


def download(version, destination):
    url = f'https://downloads.unidata.ucar.edu/netcdf-fortran/{version}/netcdf-fortran-{version}.tar.gz'
    with urllib.request.urlopen(url, timeout=180) as response, destination.open('wb') as stream:
        shutil.copyfileobj(response, stream)
    if digest(destination) != RELEASES[version]:
        destination.unlink()
        raise RuntimeError(f'NetCDF-Fortran {version}: source checksum mismatch')


def build(compiler, directory):
    version_text = output(['nc-config', '--version'])
    match = re.search(r'(\d+)\.(\d+)\.(\d+)', version_text)
    if not match or tuple(map(int, match.groups())) < (4, 7, 4):
        raise RuntimeError(f'unsupported NetCDF-C version: {version_text}')
    version = '4.6.4' if tuple(map(int, match.groups())) >= (4, 9, 2) else '4.5.4'
    archive = directory / 'source.tar.gz'
    download(version, archive)
    work = directory / 'build'
    work.mkdir()
    log = directory / 'build.log'
    run(['tar', '-xzf', str(archive), '-C', str(work)], directory, log)
    source = work / f'netcdf-fortran-{version}'
    prefix = directory / 'install'
    env = os.environ.copy()
    # Model precision flags and an installed netcdf.mod must not enter this build.
    for key in ['FFLAGS', 'FCFLAGS', 'F77', 'CPPFLAGS', 'CFLAGS', 'LDFLAGS', 'LIBS',
                'CPATH', 'C_INCLUDE_PATH', 'CPLUS_INCLUDE_PATH', 'MAKEFLAGS', 'MFLAGS', 'MAKELEVEL']:
        env.pop(key, None)
    # NetCDF's serial tests use large arrays; -fopenmp moves them onto the stack.
    library_compiler = [arg for arg in compiler if arg != '-fopenmp']
    env.update(FC=shlex.join(library_compiler), F77=shlex.join(library_compiler),
               CPPFLAGS='-I' + str(source / 'fortran') + ' ' + output(['nc-config', '--cflags']),
               LDFLAGS='-L' + output(['nc-config', '--libdir']),
               LIBS=output(['nc-config', '--libs']))
    run(['./configure', '--prefix=' + str(prefix), '--libdir=' + str(prefix / 'lib'),
         '--disable-shared', '--enable-static'], source, log, env)
    run(['make', '-j2'], source, log, env)
    run(['make', '-j2', 'check'], source, log, env)
    run(['make', 'install'], source, log, env)
    flags = shlex.split(output([str(prefix / 'bin/nf-config'), '--flibs']))
    probe(compiler, prefix / 'include', flags, directory, directory / 'probe.log')
    return prefix / 'include', flags


def make_value(value):
    if any(char in value for char in '\n\r'):
        raise RuntimeError('newline in generated Make configuration')
    return value.replace('$', '$$').replace('#', '\\#')


def prepare(args, key, compiler):
    directory = Path(args.cache_dir).resolve() / key
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        config = directory / 'config.mk'
        manifest = directory / 'selection.json'
        if manifest.exists():
            try:
                saved = json.loads(manifest.read_text())
                probe(compiler, saved['include'], saved['flags'], directory, directory / 'probe.log')
            except (subprocess.CalledProcessError, OSError, ValueError, KeyError) as error:
                config.unlink(missing_ok=True)
                raise RuntimeError(f'cached NetCDF failed validation; inspect {directory / "probe.log"}; remove {directory} before retrying') from error
            include, flags = Path(saved['include']), saved['flags']
        else:
            include = Path(args.include).resolve()
            lib = Path(args.lib).resolve()
            flags = ['-L' + str(lib), '-Wl,-rpath,' + str(lib), '-lnetcdff'] + shlex.split(output(['nc-config', '--libs']))
            try:
                probe(compiler, include, flags, directory, directory / 'system-probe.log')
            except subprocess.CalledProcessError:
                print(f'System NetCDF-Fortran is incompatible; rebuilding for {shlex.join(compiler)}. Logs: {directory}', file=sys.stderr)
                if (directory / 'build').exists():
                    raise RuntimeError(f'previous dependency build failed; inspect {directory / "build.log"}; remove {directory} before retrying')
                include, flags = build(compiler, directory)
                saved = {'include': str(include), 'flags': flags, 'fallback': True}
            else:
                saved = {'include': str(include), 'flags': flags, 'fallback': False}
        if getattr(args, 'check_fc', None):
            try:
                probe(resolve_compiler(args.check_fc), include, flags, directory, directory / 'cama-probe.log')
            except (subprocess.CalledProcessError, OSError, RuntimeError) as error:
                config.unlink(missing_ok=True)
                raise RuntimeError(f'CaMa compiler cannot use selected NetCDF; align FF/FCMP/FC; inspect {directory / "cama-probe.log"}') from error
        if not manifest.exists():
            manifest.write_text(json.dumps(saved, indent=2) + '\n')
        link = shlex.join(['-L' + str(include.parent / 'lib')]) if saved['fallback'] else ''
        content = f'NETCDF_F_INCLUDE := {make_value(str(include))}\nNETCDF_F_LDFLAGS := {make_value(link)}\n'
        if not config.exists() or config.read_text() != content:
            temporary = directory / 'config.mk.tmp'
            temporary.write_text(content)
            temporary.replace(config)
        return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', required=True)
    parser.add_argument('--check-fc', help='also require the active CaMa compiler to use the selected library')
    parser.add_argument('--include', required=True)
    parser.add_argument('--lib', required=True)
    parser.add_argument('--cache-dir', required=True)
    parser.add_argument('--key', action='store_true')
    args = parser.parse_args()
    key, compiler = identity(args)
    print(key if args.key else prepare(args, key, compiler))


if __name__ == '__main__':
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError, ValueError) as error:
        sys.exit(f'NetCDF preparation failed: {error}')
