"""Prove the one-vs-block comparator catches what it must.

A comparator that only ever prints "identical" is worse than none: it makes the
parity run look like evidence. Each failure mode the plan names is provoked
here against synthetic pairs, so the tool is under test before it is trusted to
sign off a real case.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from river_hist_compare import main as compare_main  # noqa: E402

NLON, NLAT, NT = 6, 4, 2


def _write(path: pathlib.Path, *, drop_var=False, extra_var=False, transpose=False,
           perturb=0.0, bad_attr=False, dim_size=NLAT) -> None:
    with Dataset(path, "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("lat_ucat", dim_size)
        ds.createDimension("lon_ucat", NLON)
        ds.createVariable("time", "f8", ("time",))[:] = np.arange(NT, dtype=np.float64)
        if not drop_var:
            dims = ("time", "lon_ucat", "lat_ucat") if transpose else ("time", "lat_ucat", "lon_ucat")
            v = ds.createVariable("f_rivsto", "f8", dims)
            v.long_name = "river channel storage"
            v.units = "m^3" if not bad_attr else "m3"
            v.missing_value = -1.0e36
            base = np.arange(dim_size * NLON, dtype=np.float64).reshape(dim_size, NLON)
            for t in range(NT):
                data = base + t + perturb
                v[t, :, :] = data.T if transpose else data
        if extra_var:
            e = ds.createVariable("f_only_here", "f8", ("time", "lat_ucat", "lon_ucat"))
            e[:] = 0.0


@pytest.fixture
def pair(tmp_path: pathlib.Path):
    ref, cand = tmp_path / "ref", tmp_path / "cand"
    ref.mkdir(); cand.mkdir()
    return ref, cand


def _run(ref, cand) -> int:
    return compare_main([str(ref), str(cand), "--label", "test"])


def test_identical_directories_pass(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc")
    assert _run(ref, cand) == 0


def test_missing_variable_is_caught(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", drop_var=True)
    assert _run(ref, cand) == 1


def test_extra_variable_is_caught(pair) -> None:
    """A field present only in block mode is as wrong as a missing one."""
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", extra_var=True)
    assert _run(ref, cand) == 1


def test_transposed_dimensions_are_caught(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", transpose=True)
    assert _run(ref, cand) == 1


def test_value_difference_is_caught(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", perturb=1.0e-6)
    assert _run(ref, cand) == 1


def test_tiny_difference_is_still_reported_not_waved_through(pair) -> None:
    """Rearranged data must reproduce bitwise; 'close enough' is a finding."""
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", perturb=1.0e-15)
    assert _run(ref, cand) == 1


def test_attribute_difference_is_caught(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", bad_attr=True)
    assert _run(ref, cand) == 1


def test_dimension_size_difference_is_caught(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc", dim_size=NLAT + 1)
    assert _run(ref, cand) == 1


def test_missing_file_is_caught(pair) -> None:
    ref, cand = pair
    _write(ref / "h.nc")
    assert _run(ref, cand) == 1


def test_shard_and_marker_files_are_ignored(pair) -> None:
    """The sharded path's own artefacts are not part of the 'one' schema."""
    ref, cand = pair
    _write(ref / "h.nc")
    _write(cand / "h.nc")
    _write(cand / "h_shard00000.nc")
    (cand / "h.nc.complete").write_text("done")
    (cand / "h.nc.pending").write_text("todo")
    assert _run(ref, cand) == 0
