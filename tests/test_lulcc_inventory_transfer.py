"""Execute production LULCC transfer and Tracer pool remappers on synthetic areas."""

from pathlib import Path
import os
import shutil
import subprocess
import sys

import pytest

from fortran_test_support import netcdf_fortran_flags

ROOT = Path(__file__).resolve().parents[1]


def test_lulcc_inventory_transfer(tmp_path):
    driver = (ROOT / "main/LULCC/MOD_Lulcc_Driver.F90").read_text()
    initialize = (ROOT / "main/LULCC/MOD_Lulcc_Initialize.F90").read_text()
    snapshot = (ROOT / "main/LULCC/MOD_Lulcc_Vars_TimeInvariants.F90").read_text().split(
        "SUBROUTINE SAVE_LulccTimeInvariants", 1)[1]
    assert driver.index("CALL lulcc_check_inventory_transfer") < driver.index("CALL LulccMassEnergyConserve()")
    assert "inventory_trace, new_patch_area, old_patch_area" in driver
    assert "landpatch%pctshared, landpatch_%pctshared" not in driver
    assert initialize.index("old_has_patches = numpatch > 0") < initialize.index("CALL mesh_free_mem")
    assert initialize.index("IF (old_has_patches .neqv. (numpatch > 0))") < initialize.index(
        "CALL deallocate_TimeInvariants")
    assert snapshot.index("numpatch_ = numpatch") < snapshot.index("IF (numpatch > 0) THEN")
    assert snapshot.index("numelm_ = numelm") < snapshot.index("IF (numpatch > 0) THEN")
    build = Path(os.environ.get("COLM_BLD_DIR", ROOT / ".bld")).resolve()
    library = Path(os.environ.get("COLM_LIB", ROOT / "libcolm.a")).resolve()
    compiler = shutil.which("mpif90")
    if not compiler or not library.exists() or not (build / "mod_tracer_defs.mod").exists():
        pytest.skip("MPI CoLM Tracer build required")

    include_dir = tmp_path / "include"
    module_dir = tmp_path / "modules"
    include_dir.mkdir()
    module_dir.mkdir()
    definition = (ROOT / "include/define.h").read_text()
    (include_dir / "define.h").write_text(definition.replace("#undef LULCC", "#define LULCC", 1))
    includes, libs = netcdf_fortran_flags()
    flags = ["-cpp", "-fopenmp", "-fdefault-real-8", "-ffree-form", "-ffree-line-length-0",
             "-fallow-argument-mismatch", "-fcheck=all", "-ffunction-sections", "-fdata-sections",
             "-I" + str(include_dir), "-I" + str(ROOT / "include"), "-I" + str(module_dir),
             "-I" + str(build), "-J" + str(module_dir), *includes]
    objects = []
    for path in (ROOT / "main/LULCC/MOD_Lulcc_Vars_TimeInvariants.F90",
                 ROOT / "main/LULCC/MOD_Lulcc_Initialize.F90",
                 ROOT / "main/LULCC/MOD_Lulcc_MassEnergyConserve.F90",
                 ROOT / "main/LULCC/MOD_Lulcc_Driver.F90"):
        obj = tmp_path / (path.stem + ".o")
        result = subprocess.run([compiler, *flags, "-c", str(path), "-o", str(obj)],
                                capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        objects.append(obj)
    probe = ROOT / "tests/lulcc_inventory_transfer_probe.F90"
    binary = tmp_path / "lulcc_inventory_probe"
    tail = (["-framework", "Accelerate", "-Wl,-dead_strip"] if sys.platform == "darwin"
            else ["-llapack", "-lblas", "-Wl,--gc-sections"])
    result = subprocess.run([compiler, *flags, str(probe), str(objects[-1]), str(library),
                             *libs, *tail, "-o", str(binary)],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr

    good = subprocess.run([str(binary)], capture_output=True, text=True, timeout=30)
    assert good.returncode == 0, good.stdout + good.stderr
    for scenario, diagnostic in (
        ("mismatch", "source-class physical area mismatch"),
        ("missing_donor", "source-class physical area mismatch"),
        ("missing_target", "source-class physical area mismatch"),
        ("footprint", "element physical footprint changed"),
        ("bad_row", "incomplete source trace"),
        ("nan_area", "invalid inventory area or transfer trace"),
    ):
        bad = subprocess.run([str(binary), scenario], capture_output=True, text=True, timeout=30)
        assert bad.returncode != 0, (scenario, bad.stdout, bad.stderr)
        assert diagnostic in bad.stdout + bad.stderr, (scenario, bad.stdout, bad.stderr)
