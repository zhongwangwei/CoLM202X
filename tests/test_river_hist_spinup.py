"""The spinup hist_out early return drops diagnostics, not pending routing inputs."""
from pathlib import Path
import os
import shutil
import subprocess

import pytest

from fortran_test_support import netcdf_fortran_flags

ROOT = Path(__file__).resolve().parents[1]
FLAGS = ["-fopenmp", "-fdefault-real-8", "-ffree-form", "-cpp",
         "-ffree-line-length-0", "-fallow-argument-mismatch", "-w", "-g", "-fbounds-check", "-DUSEMPI"]


def _run(tmp_path, hist_source):
    compiler, launcher = shutil.which("mpif90"), shutil.which("mpiexec") or shutil.which("mpirun")
    bld = Path(os.environ.get("COLM_BLD_DIR", "")).resolve()
    lib = Path(os.environ.get("COLM_LIB", "")).resolve()
    if not compiler or not launcher or not bld.is_dir() or not lib.is_file():
        pytest.skip("MPI compiler and COLM_BLD_DIR/COLM_LIB build required")
    includes, libs = netcdf_fortran_flags()
    objects = []
    for name, source in (("namelist", ROOT / "share/MOD_Namelist.F90"),
                         ("hist", hist_source),
                         ("harness", ROOT / "tests/river_hist_spinup_harness.F90")):
        obj = tmp_path / f"{name}.o"
        built = subprocess.run([compiler, *FLAGS, f"-I{ROOT / 'include'}", f"-I{tmp_path}",
                                f"-I{bld}", *includes, f"-J{tmp_path}", "-c", str(source),
                                "-o", str(obj)], capture_output=True, text=True, timeout=180)
        assert built.returncode == 0, built.stdout + built.stderr
        objects.append(str(obj))
    exe = tmp_path / "spinup"
    linked = subprocess.run([compiler, *FLAGS, *reversed(objects), str(lib), *libs,
                             *(["-framework", "Accelerate"] if os.uname().sysname == "Darwin"
                               else ["-llapack", "-lblas"]), "-o", str(exe)],
                            capture_output=True, text=True, timeout=180)
    assert linked.returncode == 0, linked.stdout + linked.stderr
    return subprocess.run([launcher, "-n", "1", str(exe)], capture_output=True,
                          text=True, timeout=30)


def test_spinup_flushes_only_route_diagnostics(tmp_path):
    result = _run(tmp_path, ROOT / "main/MOD_Hist.F90")
    assert result.returncode == 0 and "SPINUP_ROUTE_HISTORY_OK" in result.stdout, result.stdout + result.stderr


def test_spinup_harness_fails_without_route_flush(tmp_path):
    source = (ROOT / "main/MOD_Hist.F90").read_text(encoding="utf-8")
    assert source.count("CALL flush_acc_fluxes_riverlake ()") == 1
    variant = tmp_path / "MOD_Hist_no_route_flush.F90"
    variant.write_text(source.replace("CALL flush_acc_fluxes_riverlake ()", "CONTINUE"), encoding="utf-8")
    result = _run(tmp_path, variant)
    assert result.returncode != 0 and "SPINUP_ROUTE_HISTORY_OK" not in result.stdout
