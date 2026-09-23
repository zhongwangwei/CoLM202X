"""Real-MPI run of the river-history dispatcher with the master included.

tests/river_hist_route_master_harness.F90 calls route_hist_begin /
route_hist_write_ucat / _resv / _bif_matrix / route_hist_end on EVERY rank, as
the model does, in DEF_HIST_mode='block', and writes two history files in one
process.  It therefore covers what the other shard tests cannot:

* the master must stay out of the group collectives (route_hist_write_ucat and
  route_hist_write_resv are called by every rank);
* the bifurcation dimensions and pth_global_id must exist in EVERY history
  file, not just the first one of a run.

It links against a built model, like tests/run_river_hist_baseline.sh, so it is
skipped unless the build is named explicitly:

    gmake all TRACER_ENABLED=YES ...
    COLM_BLD_DIR=.bld COLM_LIB=libcolm.a python3 -m pytest tests/test_river_hist_route_mpi.py

Two negative controls rebuild the route module with the old behaviour and
require the harness to fail, so a harness that cannot fail is caught.
"""
from pathlib import Path
import glob
import os
import shutil
import subprocess

import pytest

from fortran_test_support import netcdf_fortran_flags

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tests/river_hist_route_master_harness.F90"
ROUTE = ROOT / "main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90"

TOTAL_UCAT, TOTAL_PTH, TOTAL_RESV = 37, 11, 5


def _build_dirs():
    bld, lib = os.environ.get("COLM_BLD_DIR"), os.environ.get("COLM_LIB")
    if not bld or not lib:
        pytest.skip("set COLM_BLD_DIR and COLM_LIB to a model build to run the MPI harness")
    bld, lib = Path(bld).resolve(), Path(lib).resolve()
    if not bld.is_dir() or not lib.is_file():
        pytest.skip(f"model build not found: {bld} / {lib}")
    return bld, lib


def _toolchain():
    compiler, launcher = shutil.which("mpif90"), shutil.which("mpiexec") or shutil.which("mpirun")
    if not compiler or not launcher:
        pytest.skip("MPI compiler/runtime unavailable")
    return compiler, launcher


FLAGS = ["-fopenmp", "-fdefault-real-8", "-ffree-form", "-cpp", "-ffree-line-length-0",
         "-fallow-argument-mismatch", "-w", "-g", "-fbounds-check"]


def _compile_harness(compiler, bld, lib, work, route_source=None):
    includes, libs = netcdf_fortran_flags()
    objects = []
    for source, name in ((ROOT / "share/MOD_Namelist.F90", "namelist"),
                         (route_source or ROUTE, "route")):
        obj = work / f"{name}.o"
        built = subprocess.run([compiler, *FLAGS, "-DUSEMPI", f"-I{ROOT / 'include'}",
                                f"-I{work}", f"-I{bld}", *includes, f"-J{work}",
                                "-c", str(source), "-o", str(obj)],
                               capture_output=True, text=True, timeout=180)
        assert built.returncode == 0, built.stdout + built.stderr
        objects.append(str(obj))
    harness_o = work / "harness.o"
    subprocess.run([compiler, *FLAGS, "-DUSEMPI", f"-I{work}", f"-I{bld}", *includes, f"-J{work}",
                    "-c", str(HARNESS), "-o", str(harness_o)],
                   check=True, capture_output=True, text=True, timeout=180)
    objects.insert(0, str(harness_o))
    exe = work / "hr"
    linked = subprocess.run([compiler, *FLAGS, *objects, str(lib), *libs,
                             *(["-framework", "Accelerate"] if os.uname().sysname == "Darwin"
                               else ["-llapack", "-lblas"]), "-o", str(exe)],
                            capture_output=True, text=True, timeout=180)
    assert linked.returncode == 0, linked.stdout + linked.stderr
    return exe


def _run(launcher, exe, work, ranks, groups, empty_last_group=False):
    out = work / f"out_{ranks}"
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir()
    env = dict(os.environ, RH_OUTDIR=str(out), RH_NGROUP=str(groups),
               RH_EMPTY_LAST_GROUP=str(int(empty_last_group)),
               OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
               OMPI_MCA_rmaps_base_oversubscribe="1")
    try:
        result = subprocess.run([launcher, "-n", str(ranks), str(exe)], env=env,
                                capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return out, None, "HUNG"
    return out, result.returncode, result.stdout + result.stderr


def _check_shards(out):
    netCDF4 = pytest.importorskip("netCDF4")
    for day in (1, 2):
        shards = sorted(glob.glob(f"{out}/hrtest_hist_unitcat_{day:03d}_seg*_shard*.nc"))
        assert shards, f"no shards for day {day}"
        ids, pids, rids = [], [], []
        for path in shards:
            with netCDF4.Dataset(path) as ds:
                # per-file state: every history file needs its own bifurcation dims
                assert "bifurcation_pathway_local" in ds.dimensions, (path, "no bif dim")
                assert "bifurcation_level" in ds.dimensions, (path, "no level dim")
                assert "pth_global_id" in ds.variables, (path, "no pth_global_id")
                these = list(ds["ucat_ucid"][:])
                assert list(ds["x_ucat"][:]) == [(gid - 1) % 12 + 1 for gid in these]
                assert list(ds["y_ucat"][:]) == [(gid - 1) // 12 + 1 for gid in these]
                values = list(ds["f_test"][0, :])
                for gid, value in zip(these, values):
                    assert abs(value - (gid + 0.5 + 1000 * day)) < 1e-9, (path, gid, value)
                ids += these
                pids += list(ds["pth_global_id"][:])
                if "resv_global_index" in ds.variables:
                    rids += list(ds["resv_global_index"][:])
        assert sorted(ids) == list(range(1, TOTAL_UCAT + 1)), (day, "unit-catchment ids")
        assert sorted(pids) == list(range(1, TOTAL_PTH + 1)), (day, "pathway ids")
        assert sorted(rids) == list(range(1, TOTAL_RESV + 1)), (day, "reservoir ids")


@pytest.mark.parametrize("ranks,groups,empty_last_group", [(4, 1, False), (6, 2, False),
                                                            (6, 2, True)])
def test_master_included_block_mode_writes_two_complete_files(tmp_path, ranks, groups,
                                                               empty_last_group):
    compiler, launcher = _toolchain()
    bld, lib = _build_dirs()
    exe = _compile_harness(compiler, bld, lib, tmp_path)
    out, code, output = _run(launcher, exe, tmp_path, ranks, groups, empty_last_group)
    assert code == 0 and "HISTROUTE_OK" in output, output[-1500:]
    _check_shards(out)


def _variant(tmp_path, name, old, new, count=1):
    source = ROUTE.read_text(encoding="utf-8")
    assert source.count(old) >= count, f"pre-fix pattern for {name} not found"
    path = tmp_path / f"variant_{name}.F90"
    path.write_text(source.replace(old, new), encoding="utf-8")
    return path


def test_harness_fails_without_the_master_guard(tmp_path):
    compiler, launcher = _toolchain()
    bld, lib = _build_dirs()
    variant = _variant(
        tmp_path, "noguard",
        "         IF (p_is_io .or. p_is_worker) THEN\n            IF (with_time) THEN",
        "         IF (.true.) THEN\n            IF (with_time) THEN")
    exe = _compile_harness(compiler, bld, lib, tmp_path, route_source=variant)
    _, code, output = _run(launcher, exe, tmp_path, 4, 1)
    assert code != 0 and "HISTROUTE_OK" not in output


def test_harness_fails_with_bif_dimensions_defined_once_per_run(tmp_path):
    compiler, launcher = _toolchain()
    bld, lib = _build_dirs()
    source = ROUTE.read_text(encoding="utf-8")
    start = source.index("            IF (.not. rh_bif_layout_built) THEN\n               CALL route_shard_layout_build (rh_bif_layout")
    end = source.index("            CALL route_shard_write_matrix (rh_bif_layout", start)
    old_block = source[start:end]
    old_behaviour = (
        "            IF (.not. rh_bif_layout_built) THEN\n"
        "               CALL route_shard_layout_build (rh_bif_layout, ncol_local, global_id)\n"
        "               rh_bif_layout_built = .true.\n"
        "               IF (p_is_io) CALL define_bif_shard_dims (nrow)\n"
        "            ENDIF\n")
    variant = _variant(tmp_path, "oldbif", old_block, old_behaviour)
    exe = _compile_harness(compiler, bld, lib, tmp_path, route_source=variant)
    out, code, output = _run(launcher, exe, tmp_path, 4, 1)
    assert code != 0 and "HISTROUTE_OK" not in output
