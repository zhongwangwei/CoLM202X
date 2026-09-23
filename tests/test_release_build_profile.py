"""Release flags stay explicit and cannot reuse checked-build objects."""
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def test_release_flags_and_profile_switch_guard(tmp_path):
    shutil.copy(ROOT / 'Makefile', tmp_path)
    (tmp_path / 'include').mkdir()
    for name in ('Makeoptions', 'define.h'):
        shutil.copy(ROOT / 'include' / name, tmp_path / 'include' / name)

    def make(*args, input=None):
        return subprocess.run(['make', '-s', *args], cwd=tmp_path,
                              input=input, text=True, capture_output=True)

    probe = "include Makefile\nflags:\n\t@echo $(FOPTS)\n"
    debug = make('-f', '-', 'flags', input=probe)
    release = make('-f', '-', 'BUILD_PROFILE=release', 'flags', input=probe)
    assert debug.returncode == release.returncode == 0
    assert '-fcheck=all' in debug.stdout and '-ffpe-trap=' in debug.stdout
    assert '-O2' in release.stdout
    assert '-fcheck=' not in release.stdout and '-ffpe-trap=' not in release.stdout
    assert make('mkdir_build').returncode == 0
    assert make('BUILD_PROFILE=release', 'mkdir_build').returncode != 0
    assert make('clean').returncode == 0
    assert make('BUILD_PROFILE=release', 'mkdir_build').returncode == 0
    assert make('mkdir_build').returncode != 0
    assert make('clean').returncode == 0
    (tmp_path / '.bld').mkdir()
    (tmp_path / '.bld' / 'legacy.o').touch()
    assert make('BUILD_PROFILE=release', 'mkdir_build').returncode != 0
