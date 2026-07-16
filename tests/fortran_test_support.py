from pathlib import Path
import shutil
import subprocess

import pytest


SMOKE_TIMEOUT = 5


def require_runnable_fortran_compiler(workdir: Path) -> str:
    compiler = shutil.which("gfortran") or shutil.which("mpif90")
    if compiler is None:
        pytest.skip("Fortran compiler is not available")

    source = workdir / "compiler_smoke.f90"
    executable = workdir / "compiler_smoke"
    source.write_text(
        "program compiler_smoke\nprint '(A)', 'FORTRAN_OK'\nend program compiler_smoke\n",
        encoding="utf-8",
    )
    try:
        compiled = subprocess.run(
            [compiler, str(source), "-o", str(executable)],
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=SMOKE_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        pytest.fail("Fortran compiler timed out while building the smoke test")
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)

    try:
        ran = subprocess.run(
            [str(executable)],
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=SMOKE_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        pytest.skip("Fortran executables cannot start in this test environment")
    if ran.returncode != 0:
        pytest.skip("Fortran executables cannot start in this test environment")
    if ran.stdout.strip() != "FORTRAN_OK":
        pytest.fail("Fortran compiler smoke test returned unexpected output")
    return compiler
