from pathlib import Path
import shutil
import subprocess

import pytest


# Generous on purpose.  These helpers gate every compile-and-run numerical
# test in the suite, and a skip reads as a pass in pytest's summary -- so a
# timeout that is merely tight silently deletes coverage instead of reporting
# it.  Observed on a developer machine at load average ~5: a 5 s budget made
# 65-127 tests vanish while the run still said "passed", and the only clue was
# the wall time going from 24 s to 177 s.  A busy CI runner would do the same.
SMOKE_TIMEOUT = 30


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
        # Distinguish this from a sandbox that forbids exec: a trivial print
        # program that needs more than SMOKE_TIMEOUT seconds means the machine
        # is loaded, not that Fortran is unusable.  Say so, because the result
        # is a skip and skips are easy to read as passes.
        pytest.skip(
            f"Fortran smoke executable did not finish within {SMOKE_TIMEOUT}s; "
            "the host is likely overloaded, so compile-and-run coverage is "
            "being skipped rather than verified"
        )
    if ran.returncode != 0:
        pytest.skip(
            "Fortran executables cannot start in this test environment "
            f"(exit code {ran.returncode})"
        )
    if ran.stdout.strip() != "FORTRAN_OK":
        pytest.fail("Fortran compiler smoke test returned unexpected output")
    return compiler
