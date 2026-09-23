"""mksrfdata/MOD_UnitCatchmentSubset against the reference implementation.

The Fortran module that mksrfdata uses and tools/subset_unitcatchment.py are two
implementations of one operation.  Each test cuts the same synthetic network with
both and requires identical results, variable by variable, so a fault in either
shows up as a disagreement.
"""
from pathlib import Path
import subprocess

import numpy as np
import pytest

from fortran_test_support import netcdf_fortran_flags, require_runnable_fortran_compiler
from test_subset_unitcatchment import make_mask, make_network, tool

netCDF4 = pytest.importorskip("netCDF4")

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "mksrfdata/MOD_UnitCatchmentSubset.F90"
HARNESS = ROOT / "tests/unitcatchment_subset_harness.F90"

COUNT_ATTRS = ("nseqall", "nseqmax", "nseqriv", "npthout", "dam_ndams", "source_nseqmax", "subset_bif_mode")


@pytest.fixture(scope="module")
def harness(tmp_path_factory):
    work = tmp_path_factory.mktemp("fsubset")
    compiler = require_runnable_fortran_compiler(work)
    includes, libs = netcdf_fortran_flags()
    exe = work / "subset_harness"
    built = subprocess.run(
        [compiler, "-cpp", "-Wall", "-fcheck=all", "-ffree-line-length-none", *includes, f"-J{work}",
         str(MODULE), str(HARNESS), "-o", str(exe), *libs],
        capture_output=True, text=True, timeout=180)
    assert built.returncode == 0, built.stdout + built.stderr
    assert "Warning" not in built.stderr, built.stderr
    return exe


def cut_both(harness, tmp_path, cells, closure):
    net = tmp_path / "net.nc"
    if not net.exists():
        make_network(net)
    cells_file = tmp_path / "cells.txt"
    cells_file.write_text("".join(f"{x} {y}\n" for x, y in cells))
    fortran_out = tmp_path / "fortran.nc"
    ran = subprocess.run([str(harness), str(net), str(fortran_out), "1" if closure else "0", str(cells_file)],
                         capture_output=True, text=True, timeout=60)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "SUBSET_OK" in ran.stdout

    make_mask(tmp_path / "mask.nc", cells)
    python_out = tmp_path / "python.nc"
    tool.main(["network", str(net), str(python_out), "--mask", str(tmp_path / "mask.nc"),
               "--bif", "closure" if closure else "drop"])
    return fortran_out, python_out


def assert_same_file(fortran_out, python_out):
    with netCDF4.Dataset(fortran_out) as f, netCDF4.Dataset(python_out) as p:
        f.set_auto_mask(False)
        p.set_auto_mask(False)
        assert set(f.variables) == set(p.variables)
        assert list(f.dimensions) == list(p.dimensions)
        for name in f.dimensions:
            assert len(f.dimensions[name]) == len(p.dimensions[name]), name
        for name in f.variables:
            a, b = f[name][:], p[name][:]
            assert a.dtype == b.dtype and a.shape == b.shape, name
            assert np.array_equal(a, b), name
            assert f[name].dimensions == p[name].dimensions, name
        for attr in COUNT_ATTRS:
            assert f.getncattr(attr) == p.getncattr(attr), attr


@pytest.mark.parametrize("cells,closure", [
    ([(1, 1)], True),
    ([(1, 1)], False),
    ([(5, 3), (7, 4)], True),
    ([(1, 1), (3, 2), (5, 3), (7, 4)], True),
])
def test_fortran_and_reference_implementations_agree(harness, tmp_path, cells, closure):
    fortran_out, python_out = cut_both(harness, tmp_path, cells, closure)
    assert_same_file(fortran_out, python_out)


def test_fortran_cut_has_the_expected_content(harness, tmp_path):
    fortran_out, _ = cut_both(harness, tmp_path, [(1, 1)], True)
    with netCDF4.Dataset(fortran_out) as f:
        f.set_auto_mask(False)
        assert list(f["seq_src_index"][:]) == [1, 2, 3, 8]
        assert list(f["seq_next"][:]) == [3, 3, -9, -10]
        assert list(f["bifurcation_upst"][:]) == [3, 4] and list(f["bifurcation_down"][:]) == [2, 2]
        assert list(f["dam_seq"][:]) == [2, 4]
        assert (f.nseqmax, f.npthout, f.dam_ndams, f.nseqriv) == (4, 2, 2, 3)


def test_domain_without_runoff_stops_the_program(harness, tmp_path):
    net = tmp_path / "net.nc"
    make_network(net)
    cells = tmp_path / "cells.txt"
    cells.write_text("8 4\n")
    ran = subprocess.run([str(harness), str(net), str(tmp_path / "o.nc"), "1", str(cells)],
                         capture_output=True, text=True, timeout=60)
    assert ran.returncode != 0
    assert "no unit catchment receives runoff" in ran.stdout + ran.stderr
