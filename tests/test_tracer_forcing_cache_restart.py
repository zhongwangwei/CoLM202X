"""The remembered precipitation/vapor composition must survive a hot restart."""

from pathlib import Path
from datetime import date
import os
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def source(path: str) -> str:
    return (ROOT / path).read_text()


def test_hot_restart_restores_cache_after_forcing_allocation():
    driver = source("main/CoLM.F90")
    init = source("main/TRACER/MOD_Tracer_LandPhase.F90")
    assert driver.index("CALL land_tracer_init") < driver.index("CALL tracer_forcing_init")
    assert driver.index("CALL tracer_forcing_init") < driver.index(
        "CALL tracer_forcing_read_restart"
    )
    assert "loaded_restart=tracer_loaded_restart" in driver
    assert "IF (present(loaded_restart)) loaded_restart = found_restart" in init


def test_cache_is_part_of_physical_transaction_not_history_window():
    rest = source("main/TRACER/MOD_Tracer_Rest.F90")
    writer = rest.split("SUBROUTINE write_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE write_land_tracer_restart", 1
    )[0]
    assert writer.index("write_land_tracer_transaction_marker(file_restart, 0)") < writer.index(
        "tracer_forcing_write_restart(file_restart)"
    ) < writer.rindex("write_land_tracer_transaction_marker(file_restart, 1)")


def test_cache_schema_covers_mixed_provider_and_exact_forcing_configuration():
    forcing = source("main/TRACER/MOD_Tracer_Forcing.F90")
    writer = forcing.split("SUBROUTINE tracer_forcing_write_restart", 1)[1].split(
        "END SUBROUTINE tracer_forcing_write_restart", 1
    )[0]
    identity = forcing.split("SUBROUTINE tracer_forcing_identity", 1)[1].split(
        "END SUBROUTINE tracer_forcing_identity", 1
    )[0]
    reader = forcing.split("SUBROUTINE tracer_forcing_read_restart", 1)[1].split(
        "END SUBROUTINE tracer_forcing_read_restart", 1
    )[0]

    # Unlike the compact generic land-water rows, this dimension includes
    # provider-owned CH4/particle rows and cannot truncate mixed registries.
    assert "'trc_forcing_species', ntracers" in writer
    assert "'trc_forcing_precip_last', 'trc_forcing_species', ntracers" in writer
    assert "'trc_forcing_vapor_last', 'trc_forcing_species', ntracers" in writer
    for name in (
        "DEF_forcing%dataset", "DEF_forcing%groupby", "DEF_forcing%startyr",
        "DEF_forcing%startmo", "DEF_forcing%leapyear", "DEF_dir_forcing",
        "trc_var_stream", "trc_var_itrc", "trc_var_mode", "trc_var_total",
        "trc_var_dtime", "trc_var_offset", "trc_var_fprefix",
        "trc_var_vname", "trc_var_tintalgo", "trc_var_timelog",
    ):
        assert name in identity
    assert "tracer_forcing_dim_names_match" in reader
    assert reader.index("tracer_forcing_dim_names_match(fileblock") < reader.index(
        "CALL ncio_read_vector(file_restart, 'trc_forcing_precip_last'"
    )
    assert "any(.not. ieee_is_finite(precip))" in reader
    assert "old tracer restart lacks last-valid forcing cache" in reader
    assert "IF (.not. loaded_restart .or. ntracers <= 0) RETURN" in reader


def test_lulcc_remaps_before_restart_and_preserves_cache_after_init():
    driver = source("main/CoLM.F90")
    transition = driver.split("! DO land use and land cover change simulation", 1)[1].split(
        "! Get leaf area index", 1
    )[0]
    assert transition.index("CALL tracer_forcing_lulcc_save") < transition.index(
        "CALL tracer_forcing_final") < transition.index("CALL LulccDriver")
    assert transition.index("CALL LulccDriver") < transition.index(
        "CALL tracer_forcing_init (gforc, numpatch)") < transition.index(
            "CALL tracer_forcing_lulcc_restore")
    lulcc = source("main/LULCC/MOD_Lulcc_Driver.F90")
    assert lulcc.index("CALL tracer_forcing_lulcc_remap") < lulcc.index(
        "CALL deallocate_LulccTransferTrace")
    assert "CALL tracer_forcing_read_restart" not in transition


def test_lulcc_cache_split_merge_and_post_reinit_roundtrip(tmp_path):
    """Compile the actual Fortran LULCC routines against tiny lifecycle stubs."""
    compiler = shutil.which("gfortran")
    if not compiler:
        pytest.skip("gfortran unavailable")
    forcing = source("main/TRACER/MOD_Tracer_Forcing.F90")
    start = forcing.index("   SUBROUTINE tracer_forcing_lulcc_save")
    end = forcing.index("   END SUBROUTINE tracer_forcing_lulcc_map", start) + len(
        "   END SUBROUTINE tracer_forcing_lulcc_map")
    routines = forcing[start:end]
    probe = r"""
MODULE cache_probe
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  INTEGER, PARAMETER :: r8 = KIND(1.0D0)
  INTEGER, PARAMETER :: N_land_classification = 17, CROPLAND = 12
  INTEGER :: ntracers = 3
  LOGICAL :: p_is_worker = .true., trc_runtime_forcing_enabled = .true.
  LOGICAL :: DEF_USE_PFT = .false., DEF_SOLO_PFT = .false., DEF_FAST_PC = .false.
  INTEGER :: patchtypes(N_land_classification) = 0
  TYPE element_patch_type
    REAL(r8), ALLOCATABLE :: subfrc(:)
  END TYPE
  TYPE(element_patch_type) :: elm_patch
  REAL(r8), ALLOCATABLE :: trc_forc_precip_value(:,:), trc_forc_vapor_value(:,:)
  REAL(r8), ALLOCATABLE :: lulcc_precip_old(:,:), lulcc_vapor_old(:,:)
  REAL(r8), ALLOCATABLE :: lulcc_precip_new(:,:), lulcc_vapor_new(:,:)
  REAL(r8), ALLOCATABLE :: lulcc_old_patch_area(:)
CONTAINS
  SUBROUTINE CoLM_stop(message)
    CHARACTER(*), INTENT(IN) :: message
    PRINT *, message
    ERROR STOP 1
  END SUBROUTINE
  REAL(r8) FUNCTION tracer_precip_default_ratio(itrc)
    INTEGER, INTENT(IN) :: itrc
    tracer_precip_default_ratio = -REAL(itrc,r8)
  END FUNCTION
  REAL(r8) FUNCTION tracer_vapor_default_ratio(itrc)
    INTEGER, INTENT(IN) :: itrc
    tracer_vapor_default_ratio = -10._r8*REAL(itrc,r8)
  END FUNCTION
""" + routines + r"""
END MODULE
PROGRAM test_cache
  USE cache_probe
  IMPLICIT NONE
  INTEGER :: bad, unit, i
  INTEGER :: cold(4), cnew(4)
  INTEGER(KIND=8) :: eold(4), enew(4)
  REAL(r8) :: trace(4,0:17), recovered(3,4), vapor_recovered(3,4), old_precip(3,4)
  cold = [1,1,2,1]
  cnew = [3,3,1,1]
  eold = [10_8,10_8,10_8,20_8]
  enew = [10_8,10_8,20_8,30_8]
  ALLOCATE(elm_patch%subfrc(4),trc_forc_precip_value(3,4),trc_forc_vapor_value(3,4))
  elm_patch%subfrc = [.1_r8,.3_r8,.6_r8,1._r8]
  trc_forc_precip_value(1,:) = [10._r8,20._r8,30._r8,40._r8]
  trc_forc_precip_value(2,:) = 10._r8*trc_forc_precip_value(1,:)
  trc_forc_precip_value(3,:) = 100._r8*trc_forc_precip_value(1,:)
  old_precip = trc_forc_precip_value
  trc_forc_vapor_value = 2._r8*trc_forc_precip_value
  CALL tracer_forcing_lulcc_save()
  IF (ALLOCATED(trc_forc_precip_value)) ERROR STOP 2
  trace = 0._r8
  trace(1,1:2) = .5_r8
  trace(2,1) = 1._r8
  trace(3,1) = 1._r8
  CALL tracer_forcing_lulcc_remap(cnew,enew,cold,eold,trace)
  ALLOCATE(trc_forc_precip_value(3,4),trc_forc_vapor_value(3,4))
  trc_forc_precip_value = -99._r8
  trc_forc_vapor_value = -99._r8
  CALL tracer_forcing_lulcc_restore()
  IF (ABS(trc_forc_precip_value(1,1)-23.75_r8) > 1.e-12_r8) ERROR STOP 3
  IF (ABS(trc_forc_precip_value(1,2)-17.5_r8) > 1.e-12_r8) ERROR STOP 4
  IF (ABS(trc_forc_precip_value(1,3)-40._r8) > 1.e-12_r8) ERROR STOP 5
  IF (trc_forc_precip_value(1,4) /= -1._r8) ERROR STOP 6
  IF (ABS(trc_forc_precip_value(3,1)-2375._r8) > 1.e-10_r8) ERROR STOP 7
  IF (ABS(trc_forc_vapor_value(1,1)-47.5_r8) > 1.e-12_r8) ERROR STOP 8
  IF (trc_forc_vapor_value(1,4) /= -10._r8) ERROR STOP 9
  ! Verify that the restored cache is serializable without value changes.
  ! The full-model NetCDF writer is separately covered by lifecycle checks.
  OPEN(NEWUNIT=unit,STATUS='SCRATCH',FORM='UNFORMATTED')
  WRITE(unit) trc_forc_precip_value,trc_forc_vapor_value
  REWIND(unit)
  READ(unit) recovered,vapor_recovered
  IF (ANY(recovered /= trc_forc_precip_value)) ERROR STOP 10
  IF (ANY(vapor_recovered /= trc_forc_vapor_value)) ERROR STOP 11
  CLOSE(unit)
  ! SAT has no trace: preserve same class first, else same-element area mean.
  recovered = -9._r8
  CALL tracer_forcing_lulcc_map(old_precip,recovered, &
       [2,3,1,1],enew,cold,eold,bad,old_patch_area=elm_patch%subfrc)
  IF (ABS(recovered(1,1)-30._r8) > 1.e-12_r8) ERROR STOP 12
  IF (ABS(recovered(1,2)-25._r8) > 1.e-12_r8) ERROR STOP 13
  IF (recovered(1,3) /= 40._r8 .or. recovered(1,4) /= -9._r8) ERROR STOP 14
  IF (bad /= 1) ERROR STOP 15
  ! Zero-area class is not a donor despite a positive transfer fraction.
  recovered = -9._r8
  CALL tracer_forcing_lulcc_map(old_precip,recovered,cnew,enew,cold,eold,bad, &
       trace,[.1_r8,.3_r8,0._r8,1._r8])
  IF (ABS(recovered(1,1)-17.5_r8) > 1.e-12_r8) ERROR STOP 16
  recovered = -9._r8
  CALL tracer_forcing_lulcc_map(old_precip,recovered, &
       [2,3,1,1],enew,cold,eold,bad,old_patch_area=[.1_r8,.3_r8,0._r8,1._r8])
  IF (ABS(recovered(1,1)-17.5_r8) > 1.e-12_r8) ERROR STOP 19
  ! PFT folds raw vegetation class 2 into old grouped soil class 1.
  DEF_USE_PFT = .true.
  patchtypes(11) = 2
  trace = 0._r8
  trace(1,2) = 1._r8
  ! Recreate a new LULCC cycle with the original old values.
  DEALLOCATE(trc_forc_precip_value,trc_forc_vapor_value)
  ALLOCATE(trc_forc_precip_value(3,4),trc_forc_vapor_value(3,4))
  trc_forc_precip_value = old_precip
  trc_forc_vapor_value = 2._r8*old_precip
  CALL tracer_forcing_lulcc_save()
  CALL tracer_forcing_lulcc_remap(cnew,enew,cold,eold,trace)
  ALLOCATE(trc_forc_precip_value(3,4),trc_forc_vapor_value(3,4))
  CALL tracer_forcing_lulcc_restore()
  IF (ABS(trc_forc_precip_value(1,1)-17.5_r8) > 1.e-12_r8) ERROR STOP 17
  ! FAST_PC combines raw 12 and 14 into the crop donor; special stays raw.
  DEF_USE_PFT = .false.
  DEF_FAST_PC = .true.
  patchtypes = 0
  patchtypes(11) = 2
  cold = [CROPLAND,CROPLAND,11,1]
  trace = 0._r8
  trace(1,12) = .25_r8
  trace(1,14) = .25_r8
  trace(1,11) = .5_r8
  DEALLOCATE(trc_forc_precip_value,trc_forc_vapor_value)
  ALLOCATE(trc_forc_precip_value(3,4),trc_forc_vapor_value(3,4))
  trc_forc_precip_value = old_precip
  trc_forc_vapor_value = 2._r8*old_precip
  CALL tracer_forcing_lulcc_save()
  CALL tracer_forcing_lulcc_remap(cnew,enew,cold,eold,trace)
  ALLOCATE(trc_forc_precip_value(3,4),trc_forc_vapor_value(3,4))
  CALL tracer_forcing_lulcc_restore()
  IF (ABS(trc_forc_precip_value(1,1)-23.75_r8) > 1.e-12_r8) ERROR STOP 18
  ! Workers without patches neither snapshot nor remap an unallocated map.
  DEALLOCATE(trc_forc_precip_value,trc_forc_vapor_value,elm_patch%subfrc)
  ALLOCATE(trc_forc_precip_value(3,0),trc_forc_vapor_value(3,0))
  CALL tracer_forcing_lulcc_save()
  IF (ALLOCATED(lulcc_precip_old)) ERROR STOP 20
  CALL tracer_forcing_lulcc_restore()
END PROGRAM
"""
    file = tmp_path / "cache_probe.f90"
    binary = tmp_path / "cache_probe"
    file.write_text(probe)
    build = subprocess.run(
        [compiler, "-O0", "-fcheck=all", "-ffree-line-length-none", str(file), "-o", str(binary)],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert build.returncode == 0, build.stdout + build.stderr
    run = subprocess.run([str(binary)], capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr


def _set_nml(text: str, key: str, value: str) -> str:
    pattern = rf"(?im)^(\s*{re.escape(key)}\s*=\s*)[^\n]+"
    assert re.search(pattern, text), key
    return re.sub(pattern, lambda m: m[1] + value, text)


def _get_nml(text: str, key: str) -> str:
    match = re.search(rf"(?im)^\s*{re.escape(key)}\s*=\s*([^\n]+)", text)
    assert match, key
    value = match[1].strip()
    if value.startswith(("'", '"')):
        return value.split(value[0], 2)[1]
    return value.split("!", 1)[0].strip()


def _runtime_case(tmp_path: Path):
    fixture = Path(os.environ["COLM_FORCING_CACHE_FIXTURE_ROOT"]) / "split"
    binary = Path(os.environ["COLM_FORCING_CACHE_TEST_BINARY"])
    nml = (fixture / "resume.nml").read_text()
    case = _get_nml(nml, "DEF_CASE_NAME")
    for part in ("landdata", "restart"):
        shutil.copytree(fixture / case / part, tmp_path / case / part)
    nml = _set_nml(nml, "DEF_dir_output", repr(str(tmp_path) + "/"))
    (tmp_path / "resume.nml").write_text(nml)
    return binary, case, nml


def _run(binary: Path, nml: Path):
    run = subprocess.run(
        [os.environ.get("MPIEXEC", "mpiexec"), "-n", "3", str(binary), str(nml)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=600,
        env={**os.environ, "OMP_NUM_THREADS": "1", "OMPI_MCA_rmaps_base_oversubscribe": "1"},
    )
    return run.returncode, run.stdout


@pytest.mark.skipif(
    not (os.environ.get("COLM_FORCING_CACHE_TEST_BINARY") and
         os.environ.get("COLM_FORCING_CACHE_FIXTURE_ROOT")),
    reason="set COLM_FORCING_CACHE_TEST_BINARY and COLM_FORCING_CACHE_FIXTURE_ROOT",
)
@pytest.mark.parametrize(
    ("damage", "expected"),
    [
        ("missing_vapor", "incomplete or malformed tracer forcing cache restart"),
        ("nan_precip", "non-finite or invalid tracer forcing cache restart"),
        ("changed_config", "tracer forcing cache configuration differs from restart"),
    ],
)
def test_real_mpi_restart_rejects_corrupt_or_changed_cache(tmp_path, damage, expected):
    nc = pytest.importorskip("netCDF4")
    binary, case, nml = _runtime_case(tmp_path)
    assert binary.is_file()
    if damage == "changed_config":
        forcing = Path(_get_nml(nml, "DEF_forcing_namelist"))
        forcing_text = forcing.read_text()
        current = _get_nml(forcing_text, "DEF_forcing%groupby")
        changed = "month" if current == "year" else "year"
        copy = tmp_path / "forcing_changed.nml"
        copy.write_text(_set_nml(forcing_text, "DEF_forcing%groupby", repr(changed)))
        (tmp_path / "resume.nml").write_text(_set_nml(nml, "DEF_forcing_namelist", repr(str(copy))))
    else:
        checkpoints = sorted((tmp_path / case / "restart").glob("*/" + case + "_restart_*_lc*_w*_s*.nc"))
        blocks = []
        for file in checkpoints:
            with nc.Dataset(file) as data:
                if "trc_forcing_cache_schema" in data.variables:
                    blocks.append(file)
        assert blocks, "fixture has no forcing-cache physical checkpoint"
        # The fixture's resume start is midnight, not the later final save.
        checkpoint = next(p for p in blocks if "-00000/" in str(p))
        with nc.Dataset(checkpoint, "r+") as data:
            if damage == "missing_vapor":
                data.renameVariable("trc_forcing_vapor_last", "trc_forcing_vapor_missing")
            else:
                data.variables["trc_forcing_precip_last"][0, 0] = float("nan")
    code, log = _run(binary, tmp_path / "resume.nml")
    assert code != 0 and "CoLM Execution Completed." not in log
    assert expected in log


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "COLM_FORCING_CACHE_TEST_BINARY", "COLM_FORCING_CACHE_FIXTURE_ROOT",
        "COLM_FORCING_CACHE_CH4_PARAM",
    )),
    reason="set forcing-cache binary, fixture root, and CH4 parameter file",
)
def test_real_mpi_mixed_isotopes_and_ch4_cache_roundtrip(tmp_path):
    nc = pytest.importorskip("netCDF4")
    binary, case, _ = _runtime_case(tmp_path)
    nml = (Path(os.environ["COLM_FORCING_CACHE_FIXTURE_ROOT"]) / "split" / "run.nml").read_text()
    nml = _set_nml(nml, "DEF_dir_output", repr(str(tmp_path) + "/"))
    for key, value in (
        ("DEF_TRACER_NUM", "3"),
        ("DEF_TRACER_NAMES", repr("H2_18O,HDO,CH4")),
        ("DEF_TRACER_TYPES", repr("isotope,isotope,gas")),
        ("DEF_USE_PLANTHYDRAULICS", ".false."),
        ("DEF_simulation_time%end_day", _get_nml(nml, "DEF_simulation_time%start_day")),
        ("DEF_simulation_time%end_sec", "3600"),
    ):
        nml = _set_nml(nml, key, value)
    files = _get_nml(nml, "DEF_TRACER_PARAM_FILES")
    nml = _set_nml(nml, "DEF_TRACER_PARAM_FILES", repr(files + ",CH4:" +
        os.environ["COLM_FORCING_CACHE_CH4_PARAM"]))
    nml = nml.replace("&nl_colm", "&nl_colm\n   DEF_USE_Dynamic_Wetland = .true.", 1)
    (tmp_path / "run.nml").write_text(nml)
    code, log = _run(binary, tmp_path / "run.nml")
    assert code == 0 and "CoLM Execution Completed." in log
    year = int(_get_nml(nml, "DEF_simulation_time%start_year"))
    month = int(_get_nml(nml, "DEF_simulation_time%start_month"))
    day = int(_get_nml(nml, "DEF_simulation_time%start_day"))
    stamp = f"{year}-{date(year, month, day).timetuple().tm_yday:03d}-03600"
    blocks = []
    for file in (tmp_path / case / "restart" / stamp).glob(case + "_restart_*_lc*_w*_s*.nc"):
        with nc.Dataset(file) as data:
            if "trc_forcing_precip_last" in data.variables:
                blocks.append(file)
    assert blocks
    for file in blocks:
        with nc.Dataset(file) as data:
            assert len(data.dimensions["trc_forcing_species"]) == 3
            assert len(data.dimensions["trc_land_transport"]) == 2
            for name in ("trc_forcing_precip_last", "trc_forcing_vapor_last"):
                assert data.variables[name].dimensions[-1] == "trc_forcing_species"
    nml = _set_nml(nml, "DEF_simulation_time%start_sec", "3600")
    nml = _set_nml(nml, "DEF_simulation_time%end_sec", "7200")
    (tmp_path / "resume.nml").write_text(nml)
    code, log = _run(binary, tmp_path / "resume.nml")
    assert code == 0 and "CoLM Execution Completed." in log
