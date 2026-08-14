"""Prove the Step-1 schema/value lock catches the failures the plan names.

Plan 1.1 acceptance: 测试能识别变量缺失、维度转置和 ID 错位。

The reference file itself needs a working MPI runtime to produce, so these
tests build equivalent synthetic files with netCDF4 and then corrupt them one
way at a time. That keeps the detector itself under test even where the MPI
harness cannot run.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from river_hist_schema_lock import (  # noqa: E402
    SPVAL,
    check_bif_values,
    check_ucat_values,
    compare_schema,
    describe,
    expected_ucat_grid,
)

TOTAL_UCAT = 60
NLON, NLAT = 10, 8
NPTHLEV, NPTHOUT = 3, 7
NTIME = 2
VARNAME = "f_ucat_synth01"
BIFNAME = "f_bifflw_lev"


def _write(path: pathlib.Path, *, drop_var=False, transpose=False, shuffle_ids=False,
           shuffle_pathways=False) -> pathlib.Path:
    with Dataset(path, "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("lat_ucat", NLAT)
        ds.createDimension("lon_ucat", NLON)
        ds.createDimension("bifurcation_level", NPTHLEV)
        ds.createDimension("bifurcation_pathway", NPTHOUT)

        ds.createVariable("lat_ucat", "f8", ("lat_ucat",))[:] = np.linspace(90, -90, NLAT)
        ds.createVariable("lon_ucat", "f8", ("lon_ucat",))[:] = np.linspace(-180, 180, NLON)
        tvar = ds.createVariable("time", "f8", ("time",))
        tvar[:] = np.arange(1, NTIME + 1, dtype=np.float64)

        if not drop_var:
            dims = ("time", "lon_ucat", "lat_ucat") if transpose else ("time", "lat_ucat", "lon_ucat")
            var = ds.createVariable(VARNAME, "f8", dims)
            var.long_name = "baseline synthetic unitcat field"
            var.units = "m"
            var.missing_value = SPVAL
            for itime in range(1, NTIME + 1):
                grid = expected_ucat_grid(TOTAL_UCAT, NLON, NLAT, 1, itime)
                if shuffle_ids:
                    # Same values, wrong cells: exactly what rank-order
                    # reassembly of shards would produce.
                    flat = grid.ravel().copy()
                    valid = np.nonzero(flat != SPVAL)[0]
                    flat[valid] = flat[valid][::-1]
                    grid = flat.reshape(grid.shape)
                var[itime - 1, :, :] = grid.T if transpose else grid

        # Same axis order the real writer produces: CoLM passes dimension names
        # fastest-first, so a Fortran (level, pathway) array lands on disk as
        # (time, pathway, level).
        bif = ds.createVariable(BIFNAME, "f8",
                                ("time", "bifurcation_pathway", "bifurcation_level"))
        bif.long_name = "baseline synthetic bifurcation pathway flow"
        bif.units = "m^3/s"
        bif.missing_value = SPVAL
        for itime in range(1, NTIME + 1):
            mat = np.empty((NPTHLEV, NPTHOUT), dtype=np.float64)
            for gid in range(1, NPTHOUT + 1):
                for lev in range(1, NPTHLEV + 1):
                    mat[lev - 1, gid - 1] = gid + 1e-3 * lev + 1e-6 * itime
            if shuffle_pathways:
                mat = mat[:, ::-1]
            bif[itime - 1, :, :] = mat.T
    return path


@pytest.fixture
def reference(tmp_path: pathlib.Path) -> pathlib.Path:
    return _write(tmp_path / "ref.nc")


def test_reference_matches_itself(reference: pathlib.Path) -> None:
    assert compare_schema(describe(reference), describe(reference)) == []
    for itime in (1, 2):
        assert check_ucat_values(reference, VARNAME, TOTAL_UCAT, NLON, NLAT, 1, itime) == []
        assert check_bif_values(reference, BIFNAME, NPTHLEV, NPTHOUT, itime) == []


def test_missing_variable_is_detected(reference: pathlib.Path, tmp_path: pathlib.Path) -> None:
    bad = _write(tmp_path / "drop.nc", drop_var=True)
    problems = compare_schema(describe(reference), describe(bad))
    assert any("missing variable" in p and VARNAME in p for p in problems)
    assert check_ucat_values(bad, VARNAME, TOTAL_UCAT, NLON, NLAT, 1, 1)


def test_transposed_dimensions_are_detected(reference: pathlib.Path, tmp_path: pathlib.Path) -> None:
    bad = _write(tmp_path / "transpose.nc", transpose=True)
    problems = compare_schema(describe(reference), describe(bad))
    assert any("dimension order differs" in p and VARNAME in p for p in problems)


def test_misplaced_unitcat_ids_are_detected(reference: pathlib.Path, tmp_path: pathlib.Path) -> None:
    bad = _write(tmp_path / "shuffle.nc", shuffle_ids=True)
    # The schema is untouched -- only the id->cell mapping is wrong, which is
    # the failure a schema-only comparison would miss.
    assert compare_schema(describe(reference), describe(bad)) == []
    problems = check_ucat_values(bad, VARNAME, TOTAL_UCAT, NLON, NLAT, 1, 1)
    assert problems and "wrong unit-catchment value" in problems[0]


def test_rank_ordered_pathways_are_detected(reference: pathlib.Path, tmp_path: pathlib.Path) -> None:
    bad = _write(tmp_path / "pathways.nc", shuffle_pathways=True)
    assert compare_schema(describe(reference), describe(bad)) == []
    problems = check_bif_values(bad, BIFNAME, NPTHLEV, NPTHOUT, 1)
    assert problems and "global-id order" in problems[0]


def test_time_axis_difference_is_detected(reference: pathlib.Path, tmp_path: pathlib.Path) -> None:
    bad = _write(tmp_path / "time.nc")
    with Dataset(bad, "a") as ds:
        ds.variables["time"][1] = 999.0
    problems = compare_schema(describe(reference), describe(bad))
    assert any(p.startswith("time[1] differs") for p in problems)


def test_missing_value_attribute_difference_is_detected(
    reference: pathlib.Path, tmp_path: pathlib.Path
) -> None:
    bad = _write(tmp_path / "attr.nc")
    with Dataset(bad, "a") as ds:
        ds.variables[VARNAME].missing_value = -9.0e33
    problems = compare_schema(describe(reference), describe(bad))
    assert any("attr missing_value differs" in p for p in problems)
