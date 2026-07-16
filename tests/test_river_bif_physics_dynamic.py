from pathlib import Path
import shutil
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = REPO_ROOT / "tests" / "run_river_bif_physics_evaluator.sh"


def test_river_bif_physics_dynamic_harness() -> None:
    if shutil.which("mpif90") is None or shutil.which("mpiexec") is None:
        pytest.skip("MPI Fortran compiler/launcher unavailable")
    if not (REPO_ROOT / ".bld").is_dir():
        pytest.skip("production .bld dependency objects unavailable")

    try:
        launcher = subprocess.run(
            ["mpiexec", "-n", "1", "/usr/bin/true"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except subprocess.TimeoutExpired:
        pytest.skip("MPI launcher cannot start processes in this environment")
    if launcher.returncode != 0:
        pytest.skip("MPI launcher cannot start processes in this environment")

    subprocess.run(["bash", str(RUNNER)], cwd=REPO_ROOT, check=True, timeout=180)
