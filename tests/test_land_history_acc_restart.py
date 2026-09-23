"""Guard the sidecar's field inventory and the write/read ordering."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_generic_history_accumulator_inventory():
    module = (ROOT / "main/MOD_Vars_1DAccFluxes.F90").read_text()
    sidecar = (ROOT / "include/land_history_restart.inc").read_text()
    declarations = set(re.findall(
        r"real\(r8\),\s*allocatable\s*::\s*(a_\w+)\s*\(",
        module.split("CONTAINS", 1)[0], flags=re.I,
    )) | {"nac_ln", "nac_dt"}
    manifest = set(re.findall(r"names\(i\)='(\w+)'", sidecar))
    writes = set(re.findall(r"CALL history_acc_write\(file, '(\w+)'", sidecar))
    reads = set(re.findall(r"CALL ncio_read_vector\(file, '(a_\w+|nac_ln|nac_dt)'", sidecar))
    assert declarations == manifest == writes == reads
    assert "filter_dt" not in manifest
    assert "decomp_vr_tmp" not in manifest


def test_generic_land_tracer_history_inventory_and_transaction():
    module = (ROOT / "main/TRACER/MOD_Tracer_Vars.F90").read_text()
    tracer = (ROOT / "include/tracer_land_history_restart.inc").read_text()
    sidecar = (ROOT / "include/land_history_restart.inc").read_text()
    declarations = set(re.findall(
        r"real\(r8\),\s*allocatable\s*::\s*(a_(?:trc|water)\w+)\s*\(",
        module.split("CONTAINS", 1)[0], flags=re.I,
    ))
    manifest = set(re.findall(r"^\s*'(a_(?:trc|water)\w+)'(?:, &|\])", tracer, flags=re.M))
    writes = set(re.findall(r"CALL ncio_write_vector\(file, '(a_(?:trc|water)\w+)'", tracer))
    reads = set(re.findall(r"CALL ncio_read_vector\(file, '(a_(?:trc|water)\w+)'", tracer))
    assert len(declarations) == 44
    assert declarations == manifest == writes == reads
    assert "DO i=1,TRC_HISTORY_FIELDS" in tracer
    assert "IF (ranks(i)>1) CALL ncio_define_dimension_vector" in tracer
    assert "IF (ranks(i)>2) CALL ncio_define_dimension_vector" in tracer
    assert tracer.index("CALL ncio_define_dimension_vector(file, landpatch") < tracer.index(
        "CALL ncio_write_vector(file, 'a_trc_precip'"
    )
    assert "CALL tracer_history_write(file)" in sidecar
    assert sidecar.index("CALL tracer_history_write(file)") < sidecar.index(
        "CALL ncio_write_vector(file, 'history_complete'"
    )
    assert "CALL tracer_history_restore(file)" in sidecar
    assert "generic land-tracer history missing from active restart window" in sidecar
    assert "trc_hist_descriptor" in tracer


def test_real_land_tracer_history_continuous_matches_midday_resume():
    """Opt-in real MPI fixture: compare every final tracer-history field and mask."""
    import os

    fixture = os.environ.get("COLM_LAND_TRACER_HISTORY_ROOT")
    if not fixture:
        pytest.skip("set COLM_LAND_TRACER_HISTORY_ROOT to the real day/split fixture")
    np = pytest.importorskip("numpy")
    netcdf4 = pytest.importorskip("netCDF4")
    root = Path(fixture)
    case = "AB_Amazon_2003"
    continuous = root / "continuous" / case / "history"
    resumed = root / "split" / case / "history"
    files = sorted(continuous.glob("*_hist_tracer_*.nc"))
    assert files, f"no real tracer-history files in {continuous}"
    compared = set()
    for source in files:
        target = resumed / source.name
        assert target.exists(), target
        with netcdf4.Dataset(source) as before, netcdf4.Dataset(target) as after:
            fields = [name for name in before.variables
                      if name.startswith(("f_", "history_window"))]
            assert fields, source
            for name in fields:
                assert name in after.variables, (target, name)
                lhs, rhs = before[name][-1], after[name][-1]
                np.testing.assert_array_equal(np.ma.getmaskarray(lhs), np.ma.getmaskarray(rhs), err_msg=name)
                np.testing.assert_array_equal(np.ma.filled(lhs, 0), np.ma.filled(rhs, 0), err_msg=name)
            compared.update(fields)
    assert any(name.startswith("f_trc_conc_wa") for name in compared)
    assert any(name.startswith("f_trc_delta_precip") for name in compared)


def test_checkpoint_restore_boundaries():
    driver = (ROOT / "main/CoLM.F90").read_text()
    history = (ROOT / "main/MOD_Hist.F90").read_text()
    sidecar = (ROOT / "include/land_history_restart.inc").read_text()
    assert driver.index("CALL hist_init (dir_hist)") < driver.index("CALL read_history_acc_restart")
    assert history.index("CALL accumulate_fluxes ()") < history.index(
        "CALL write_history_acc_restart(restart_date, casename, dir_restart)",
        history.index("CALL accumulate_fluxes ()"),
    ) < history.index("IF (lwrite) THEN")
    assert sidecar.index("CALL history_acc_transfer(file, .true.") < sidecar.index(
        "CALL ncio_write_vector(file, 'history_complete'"
    )
    assert "marker=0._r8\n      CALL ncio_write_vector(file, 'history_complete'" in sidecar
    assert sidecar.index("any_sidecar=history_acc_any_file(file)") < sidecar.index(
        "IF (.not.any(markers)) THEN"
    )
    assert "Warning: legacy land restart lacks history accumulators" in sidecar
    assert "ERROR: incomplete land-history restart sidecar" in sidecar
    assert "ERROR: incomplete land-history sidecar commit marker" in sidecar
    assert "ERROR: invalid land-history sample count" in sidecar
    assert "ieee_is_finite(marker)" in sidecar
    assert "IF (.not.history_saved_raw)" in driver
    assert "IF (.not. (itstamp < etstamp) .and. .not.natural_boundary)" in history
    assert "history_sidecar_required" in sidecar
    assert "ERROR: land restart requires a missing history sidecar" in sidecar
    assert driver.index("CALL write_history_acc_restart (jdate") < driver.index(
        "CALL WRITE_TimeVariables (jdate, lc_year"
    ) < driver.index("CALL mark_history_acc_restart (jdate") < driver.index(
        "CALL complete_history_acc_restart (jdate"
    )


def test_dry_lake_frcsat_is_not_an_undefined_history_sample():
    source = (ROOT / "main/CoLMMAIN.F90").read_text()
    assert "IF (is_dry_lake) frcsat = spval" in source


def test_river_history_partial_restart_restores_complete_window():
    driver = (ROOT / "main/CoLM.F90").read_text()
    sidecar = (ROOT / "include/land_history_restart.inc").read_text()
    route = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90").read_text()
    assert driver.index("CALL grid_riverlake_flow_init") < driver.index(
        "CALL restore_river_history_acc_restart"
    )
    assert "CALL write_gridriverlake_hist_restart(river_file)" in sidecar
    assert "CALL write_tracer_history_acc_restart(river_file)" in sidecar
    assert "CALL sediment_history_acc_sidecar(river_file, .true.)" in sidecar
    assert "primary_identity_validated=.true., strict=.true." in sidecar
    assert "history_river_required" in sidecar
    assert "history_window_seconds" in route
    assert "history_window_end_minutes" in route


def test_real_river_history_midday_resume_has_full_day_window():
    """Opt-in three-rank real-input 24 h versus 12+12 h regression."""
    import os
    np = pytest.importorskip("numpy")
    netcdf4 = pytest.importorskip("netCDF4")
    continuous_dir = os.environ.get("COLM_RIVER_HISTORY_CONTINUOUS_DIR")
    split_dir = os.environ.get("COLM_RIVER_HISTORY_SPLIT_DIR")
    if not continuous_dir or not split_dir:
        pytest.skip("set both COLM_RIVER_HISTORY_*_DIR fixture paths")
    name = "AB_Amazon_2003_hist_unitcat_2003-01-01.nc"
    with netcdf4.Dataset(Path(continuous_dir) / name) as continuous, \
         netcdf4.Dataset(Path(split_dir) / name) as resumed:
        np.testing.assert_allclose(continuous["history_window_seconds"][:], [86400.], rtol=0, atol=1e-9)
        np.testing.assert_allclose(resumed["history_window_seconds"][:], [43200., 86400.], rtol=0, atol=1e-9)
        continuous_start = (continuous["history_window_end_minutes"][:] -
                            continuous["history_window_seconds"][:] / 60)
        resumed_start = (resumed["history_window_end_minutes"][:] -
                         resumed["history_window_seconds"][:] / 60)
        np.testing.assert_allclose(resumed_start, continuous_start[0], rtol=0, atol=1e-9)
        for field in ("f_discharge", "f_rivsto", "f_bifflw_lev",
                      "f_sedout_1", "f_trc_flux_H2_18O", "f_trc_flux_HDO"):
            np.testing.assert_array_equal(continuous[field][-1], resumed[field][-1], err_msg=field)
    for suffix in ("hist", "hist_tracer"):
        name = f"AB_Amazon_2003_{suffix}_2003-01-01.nc"
        with netcdf4.Dataset(Path(continuous_dir) / name) as continuous, \
             netcdf4.Dataset(Path(split_dir) / name) as resumed:
            np.testing.assert_allclose(continuous["history_window_seconds"][:], [86400.], rtol=0, atol=1e-9)
            # Generic gridded files are recreated on resume; the final full-day
            # record replaces the provisional terminal partial record.
            np.testing.assert_allclose(resumed["history_window_seconds"][:], [86400.], rtol=0, atol=1e-9)
            np.testing.assert_array_equal(continuous["history_window_end_minutes"][-1],
                                          resumed["history_window_end_minutes"][-1])


def test_river_history_sidecar_missing_numerator_is_rejected(tmp_path):
    """An apparently valid denominator must not hide a missing discharge sum."""
    import os
    import shutil
    import subprocess

    source = os.environ.get("COLM_RIVER_HISTORY_SPLIT_CASE_DIR")
    binary = os.environ.get("COLM_RIVER_HISTORY_COLM_X")
    if not source or not binary:
        pytest.skip("set COLM_RIVER_HISTORY_SPLIT_CASE_DIR and COLM_RIVER_HISTORY_COLM_X")
    netcdf4 = pytest.importorskip("netCDF4")
    source = Path(source)
    copied = tmp_path / "split"
    shutil.copytree(source, copied)
    namelist = copied / "resume.nml"
    namelist.write_text(re.sub(
        r"(?im)^(\s*DEF_dir_output\s*=\s*)[^\n]+",
        lambda match: match[1] + repr(str(copied) + "/"), namelist.read_text(),
    ))
    checkpoints = list((copied / "AB_Amazon_2003/restart").glob("*-43200"))
    assert len(checkpoints) == 1, "fixture must contain one noon restart"
    checkpoint = checkpoints[0]
    sidecars = list(checkpoint.glob("*.nc.river"))
    assert len(sidecars) == 1
    with netcdf4.Dataset(sidecars[0], "r+") as dataset:
        dataset.renameVariable("hist_discharge", "hist_discharge_removed")
    env = dict(os.environ, OMP_NUM_THREADS="1", OMPI_MCA_rmaps_base_oversubscribe="1")
    result = subprocess.run(["mpiexec", "-n", "3", binary, str(namelist)],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=180, env=env, check=False)
    assert result.returncode != 0
    assert "ERROR: missing river-history sidecar field hist_discharge" in result.stdout


@pytest.mark.parametrize("env, day", [
    ("COLM_LAND_HISTORY_RESTART_EDGE", "2003-04-02"),
    ("COLM_LAND_HISTORY_MIDNIGHT", "2003-04-03"),
])
def test_real_resume_matches_all_land_history_fields(env, day):
    """Compare continuous/resumed land history at forced and natural boundaries."""
    import os

    fixture = os.environ.get(env)
    if not fixture:
        pytest.skip(f"set {env} for the real restart fixture")
    np = pytest.importorskip("numpy")
    netcdf4 = pytest.importorskip("netCDF4")
    case = "AB_Amazon_2003"
    name = f"{case}_hist_{day}.nc"
    paths = {kind: Path(fixture) / kind / case / "history" / name
             for kind in ("continuous", "split", "resumed")}
    with netcdf4.Dataset(paths["continuous"]) as continuous, \
         netcdf4.Dataset(paths["resumed"]) as resumed:
        np.testing.assert_array_equal(continuous["time"][:], resumed["time"][:])
        checked = 0
        for field, baseline in continuous.variables.items():
            if field == "time" or "time" not in baseline.dimensions:
                continue
            assert field in resumed.variables, field
            checked += 1
            expected, actual = baseline[:], resumed[field][:]
            np.testing.assert_array_equal(
                np.ma.getmaskarray(expected), np.ma.getmaskarray(actual),
                err_msg=f"{field} mask",
            )
            np.testing.assert_allclose(
                np.ma.filled(expected, 0), np.ma.filled(actual, 0),
                rtol=2e-6, atol=2e-5,
                err_msg=field,
            )
        assert checked > 300
    split_day = "2003-04-02" if env == "COLM_LAND_HISTORY_MIDNIGHT" else day
    split_path = paths["split"].with_name(f"{case}_hist_{split_day}.nc")
    with netcdf4.Dataset(split_path) as split:
        assert len(split.dimensions["time"]) > 0  # terminal record remains


@pytest.mark.parametrize("damage, expected, succeeds", [
    ("missing", "ERROR: land restart requires a missing history sidecar", False),
    ("dimension", "ERROR: incompatible land-history sidecar shape", False),
    ("incomplete", "ERROR: incomplete land-history sidecar commit marker", False),
    ("markerless", "ERROR: incomplete land-history restart sidecar", False),
    ("legacy", "Warning: legacy land restart lacks history accumulators", True),
])
def test_sidecar_missing_or_corrupt_restart(tmp_path, damage, expected, succeeds):
    """Use a copy of a real split checkpoint; never mutate the fixture."""
    import os
    import shutil
    import subprocess

    fixture = os.environ.get("COLM_LAND_HISTORY_RESTART_EDGE")
    binary = os.environ.get("COLM_LAND_HISTORY_BINARY")
    if not fixture or not binary:
        pytest.skip("set COLM_LAND_HISTORY_RESTART_EDGE and COLM_LAND_HISTORY_BINARY")
    netcdf4 = pytest.importorskip("netCDF4")
    source = Path(fixture) / "resumed"
    copied = tmp_path / "resumed"
    shutil.copytree(source, copied)
    namelist = copied / "run.nml"
    namelist.write_text(re.sub(
        r"(?im)^(\s*DEF_dir_output\s*=\s*)[^\n]+",
        lambda match: match[1] + repr(str(copied) + "/"), namelist.read_text(),
    ))
    checkpoint = copied / "AB_Amazon_2003/restart/2003-092-43200"
    sidecars = sorted(checkpoint.glob("*_restart_hist_*.nc"))
    assert sidecars
    if damage in {"incomplete", "markerless", "legacy"}:
        for primary in checkpoint.glob("*_restart_2003-092-43200_lc*.nc"):
            with netcdf4.Dataset(primary, "r+") as dataset:
                dataset.renameVariable("history_sidecar_required", "legacy_marker")
    if damage == "dimension":
        with netcdf4.Dataset(sidecars[0], "r+") as dataset:
            dataset.renameDimension("patch", "wrong_patch")
    elif damage == "incomplete":
        for sidecar in sidecars:
            with netcdf4.Dataset(sidecar, "r+") as dataset:
                dataset["history_complete"][:] = 0.0
    elif damage == "markerless":
        for sidecar in sidecars:
            sidecar.unlink()
            with netcdf4.Dataset(sidecar, "w"):
                pass
    else:
        for sidecar in sidecars:
            sidecar.unlink()
    result = subprocess.run(
        ["mpirun", "-np", "3", binary, str(namelist)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=120, check=False,
    )
    assert (result.returncode == 0) is succeeds, result.stdout[-3000:]
    assert expected in result.stdout
