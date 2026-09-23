"""Regression checks for the optional real CaMa split/resume fixture.

Set COLM_CAMA_RESTART_EDGE to the output directory after running the local
continuous, split, and resumed integration jobs.
"""

import os
from pathlib import Path

import pytest


def test_restart_carries_raw_daily_history_and_checks_metadata() -> None:
    cama = Path(__file__).resolve().parents[1] / "extends/CaMa/src"
    source = (cama / "cmf_ctrl_restart_mod.F90").read_text().lower()
    diagnostics = (cama / "yos_cmf_diag.F90").read_text().lower()
    for field in (
        "history_elapsed", "history_rivout", "history_fldout", "history_outflw",
        "history_rivvel", "history_pthout", "history_gdwrtn", "history_runoff",
        "history_rofsub", "history_outflw_max", "history_rivdph_max",
        "history_storge_max", "history_daminf", "history_wevapex",
        "history_winfiltex", "history_pthflw",
    ):
        assert field in source + diagnostics
    read = source.split("subroutine read_rest_cdf", 1)[1].split("end subroutine read_rest_cdf", 1)[0]
    assert read.index("checkpoint time does not match") < read.index("nf90_inq_varid(ncid,'rivsto'")
    assert read.index("lon/lat grid mismatch") < read.index("nf90_inq_varid(ncid,'rivsto'")
    assert "legacy checkpoint lacks a safe daily history boundary" in read
    assert ".not.llegacy_daily_history .or. mod(kminstart,1440_jpim)/=0" in read
    coupler = (cama / "MOD_CaMa_colmCaMa.F90").read_text().lower()
    assert "llegacy_daily_history = trim(adjustl(def_hist_freq)) == 'daily'" in coupler
    assert "start_date(3)<=def_simulation_time%spinup_sec" in coupler
    history = (cama / "MOD_CaMa_Vars.F90").read_text().lower()
    for cadence in ("hourly", "daily", "monthly", "yearly"):
        assert f"case ('{cadence}')" in history
    assert "preserve_partial_history" in history


def test_real_midday_resume_matches_continuous_history() -> None:
    root = os.environ.get("COLM_CAMA_RESTART_EDGE")
    if not root:
        pytest.skip("set COLM_CAMA_RESTART_EDGE after running split/resume fixture")
    np = pytest.importorskip("numpy")
    netcdf4 = pytest.importorskip("netCDF4")
    root = Path(root)
    case = "AB_Amazon_2003"
    for kind, pattern in (
        ("restart", "restart/CaMa/2003-093-00000/*.nc"),
        ("history", "history/*hist_cama*.nc"),
    ):
        expected = next((root / "continuous" / case).glob(pattern))
        actual = next((root / "resumed" / case).glob(pattern))
        with netcdf4.Dataset(expected) as a, netcdf4.Dataset(actual) as b:
            for name in a.variables:
                if name in {"time", "lon", "lat", "lat_cama", "lon_cama"}:
                    continue
                assert name in b.variables, (kind, name)
                np.testing.assert_allclose(a[name][:], b[name][:], rtol=2e-6, atol=2e-5, err_msg=f"{kind}/{name}")
    wrong_log = (root / "wrong_date" / "run.log").read_text()
    assert "checkpoint time does not match simulation start" in wrong_log
    assert not list((root / "wrong_date" / case / "history").glob("*hist_cama*.nc"))
    grid_log = (root / "wrong_grid" / "run.log").read_text()
    assert "lon/lat grid mismatch" in grid_log
    assert not list((root / "wrong_grid" / case / "history").glob("*hist_cama*.nc"))
    assert "legacy day-boundary checkpoint; history starts empty" in (root / "legacy_daily" / "run.log").read_text()
    assert list((root / "legacy_daily" / case / "history").glob("*hist_cama*.nc"))
    assert "legacy checkpoint lacks a safe daily history boundary" in (root / "legacy_monthly" / "run.log").read_text()
    assert not list((root / "legacy_monthly" / case / "history").glob("*hist_cama*.nc"))
    assert "legacy checkpoint lacks a safe daily history boundary" in (root / "legacy_spinup" / "run.log").read_text()
    for case_name, diagnostic in (
        ("wrong_grid_nan", "non-finite lon/lat grid coordinate"),
        ("wrong_units", "unsupported time units"),
        ("wrong_elapsed_negative", "invalid history elapsed time"),
        ("wrong_elapsed_nan", "non-finite history elapsed time"),
    ):
        log = (root / case_name / "run.log").read_text()
        assert diagnostic in log
        assert "STOP 9" in log


def test_real_monthly_midday_resume_matches_continuous_history() -> None:
    root = os.environ.get("COLM_CAMA_RESTART_MONTHLY")
    if not root:
        pytest.skip("set COLM_CAMA_RESTART_MONTHLY after running monthly split/resume fixture")
    np = pytest.importorskip("numpy")
    netcdf4 = pytest.importorskip("netCDF4")
    root = Path(root)
    case = "AB_Amazon_2003"
    checkpoint = next((root / "split" / case).glob("restart/CaMa/*/restart2003040212.nc"))
    with netcdf4.Dataset(checkpoint) as dataset:
        assert float(dataset["history_elapsed"][...]) == 43200.0
    for kind, pattern in (
        ("restart", "restart/CaMa/2003-093-00000/*.nc"),
        ("history", "history/*hist_cama*.nc"),
    ):
        expected = next((root / "continuous" / case).glob(pattern))
        actual = next((root / "resumed" / case).glob(pattern))
        with netcdf4.Dataset(expected) as a, netcdf4.Dataset(actual) as b:
            for name in a.variables:
                if name in {"time", "lon", "lat", "lat_cama", "lon_cama"}:
                    continue
                assert name in b.variables, (kind, name)
                np.testing.assert_allclose(a[name][:], b[name][:], rtol=2e-6, atol=2e-5, err_msg=f"{kind}/{name}")
