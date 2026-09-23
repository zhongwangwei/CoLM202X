"""tools/subset_unitcatchment.py on a small synthetic network.

Layout (grid 8 x 4, 1-degree cells, west=0 north=4; cell k has inpmat cells below):

    system A : 1 -> 3 (mouth, -9), 2 -> 3           input cells (1,1) (2,1)
    system B : 4 -> 5 (mouth)                       input cell  (5,3)
    system C : 6 -> 7 (mouth)                       input cell  (7,4)
    system D : 8 (inland depression, -10)           input cell  (3,2)

    pathways : P1 5->7 (B-C)   P2 3->2 (inside A)   P3 8->2 (D-A)
    dams     : at cells 2, 4, 6, 8
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("subset_unitcatchment", ROOT / "tools/subset_unitcatchment.py")
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)

N = 8
NEXT = [3, 3, -9, 5, -9, 7, -9, -10]
INP = {1: [(1, 1), (2, 1)], 2: [(1, 1)], 3: [(2, 1)], 4: [(5, 3)], 5: [(5, 3)],
       6: [(7, 4)], 7: [(7, 4)], 8: [(3, 2)]}
FILL = -9999


def make_network(path):
    with netCDF4.Dataset(path, "w") as f:
        for name, size in dict(nx=8, ny=4, nseqmax=N, upnmax=3, nlfp=2, npthout=3, npthlev=1,
                               inpn=2, dam_ndams=4, dam_namelen=6).items():
            f.createDimension(name, size)
        f.title = "synthetic"
        f.west, f.east, f.north, f.south = 0.0, 8.0, 4.0, 0.0
        f.nseqriv, f.nseqall, f.nseqmax = 7, N, N
        f.npthout, f.dam_ndams = 3, 4

        def put(name, dtype, dims, data):
            var = f.createVariable(name, dtype, dims)
            var[:] = data
            return var

        put("lon", "f8", ("nx",), np.arange(8) + 0.5)
        put("lat", "f8", ("ny",), 3.5 - np.arange(4))
        put("seq", "i4", ("nseqmax",), np.arange(1, N + 1))
        put("seq_x", "i4", ("nseqmax",), np.arange(1, N + 1))
        put("seq_y", "i4", ("nseqmax",), [1, 1, 1, 3, 3, 4, 4, 2])
        put("seq_next", "i4", ("nseqmax",), NEXT)
        upst = np.full((3, N), FILL)
        upst[:2, 2] = [1, 2]
        upst[0, 4] = 4
        upst[0, 6] = 6
        put("seq_upst", "i4", ("upnmax", "nseqmax"), upst)
        put("seq_upn", "i4", ("nseqmax",), (upst > 0).sum(axis=0))
        put("topo_area", "f8", ("nseqmax",), 100.0 * np.arange(1, N + 1))
        put("topo_fldhgt", "f8", ("nseqmax", "nlfp"), np.outer(np.arange(1, N + 1), [1.0, 2.0]))
        inpx = np.zeros((N, 2), dtype=int)
        inpy = np.zeros((N, 2), dtype=int)
        for k, cells in INP.items():
            for j, (x, y) in enumerate(cells):
                inpx[k - 1, j], inpy[k - 1, j] = x, y
        put("inpmat_x", "i4", ("nseqmax", "inpn"), inpx)
        put("inpmat_y", "i4", ("nseqmax", "inpn"), inpy)
        put("bifurcation_upst", "i4", ("npthout",), [5, 3, 8])
        put("bifurcation_down", "i4", ("npthout",), [7, 2, 2])
        put("bifurcation_distance", "f8", ("npthout",), [10.0, 20.0, 30.0])
        put("bifurcation_manning", "f8", ("npthlev",), [0.03])
        put("dam_seq", "i4", ("dam_ndams",), [2, 4, 6, 8])
        names = netCDF4.stringtochar(np.array(["d2", "d4", "d6", "d8"], dtype="S6"))
        put("dam_DamName", "S1", ("dam_ndams", "dam_namelen"), names)


def make_mask(path, cell_lonlat):
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("nlon", 8)
        f.createDimension("nlat", 4)
        lon = f.createVariable("longitude", "f8", ("nlon",))
        lat = f.createVariable("latitude", "f8", ("nlat",))
        lon[:] = np.arange(8) + 0.5
        lat[:] = 3.5 - np.arange(4)
        field = f.createVariable("elmindex", "i4", ("nlat", "nlon"))
        data = np.zeros((4, 8), dtype=int)
        for x, y in cell_lonlat:
            data[y - 1, x - 1] = 1
        field[:] = data


def run(tmp_path, *args):
    net = tmp_path / "net.nc"
    if not net.exists():
        make_network(net)
    out = tmp_path / "sub.nc"
    tool.main(["network", str(net), str(out), *args])
    return out


def test_bif_closure_keeps_systems_joined_by_pathways(tmp_path):
    make_mask(tmp_path / "mask.nc", [(1, 1)])
    out = run(tmp_path, "--mask", str(tmp_path / "mask.nc"))
    with netCDF4.Dataset(out) as f:
        f.set_auto_mask(False)
        # domain touches A; P3 pulls in D; B and C stay out
        assert list(f["seq_src_index"][:]) == [1, 2, 3, 8]
        assert list(f["seq"][:]) == [1, 2, 3, 4]
        assert list(f["seq_next"][:]) == [3, 3, -9, -10]
        assert list(f["seq_upst"][0, :3]) == [FILL, FILL, 1] and list(f["seq_upst"][1, :3]) == [FILL, FILL, 2]
        assert list(f["seq_upn"][:]) == [0, 0, 2, 0]
        assert list(f["topo_area"][:]) == [100.0, 200.0, 300.0, 800.0]
        assert f["topo_fldhgt"][:].tolist() == [[1, 2], [2, 4], [3, 6], [8, 16]]
        assert list(f["inpmat_x"][:, 0]) == [1, 1, 2, 3]
        # pathways P2 (3->2) and P3 (8->2) survive, renumbered; P1 (B-C) is dropped
        assert list(f["bifurcation_upst"][:]) == [3, 4] and list(f["bifurcation_down"][:]) == [2, 2]
        assert list(f["bifurcation_distance"][:]) == [20.0, 30.0]
        assert list(f["bifurcation_manning"][:]) == [0.03]
        # dams at 2 and 8 survive
        assert list(f["dam_seq"][:]) == [2, 4]
        assert [b"".join(row).decode().strip("\x00") for row in f["dam_DamName"][:]] == ["d2", "d8"]
        assert (f.nseqmax, f.nseqall, f.nseqriv, f.npthout, f.dam_ndams) == (4, 4, 3, 2, 2)
        assert len(f.dimensions["nseqmax"]) == 4 and len(f.dimensions["npthout"]) == 2
        assert f.source_nseqmax == N
        assert f.subset_bif_mode == "closure"


def test_bif_drop_discards_straddling_pathways(tmp_path):
    make_mask(tmp_path / "mask.nc", [(1, 1)])
    out = run(tmp_path, "--mask", str(tmp_path / "mask.nc"), "--bif", "drop")
    with netCDF4.Dataset(out) as f:
        f.set_auto_mask(False)
        assert list(f["seq_src_index"][:]) == [1, 2, 3]
        assert list(f["bifurcation_upst"][:]) == [3] and list(f["bifurcation_down"][:]) == [2]
        assert list(f["dam_seq"][:]) == [2]


def test_bounding_box_selects_the_same_systems(tmp_path):
    out = run(tmp_path, "--bbox", "0.1", "1.9", "3.1", "3.9")
    with netCDF4.Dataset(out) as f:
        assert list(f["seq_src_index"][:]) == [1, 2, 3, 8]


def test_domain_over_two_systems(tmp_path):
    make_mask(tmp_path / "mask.nc", [(5, 3), (7, 4)])
    out = run(tmp_path, "--mask", str(tmp_path / "mask.nc"))
    with netCDF4.Dataset(out) as f:
        assert list(f["seq_src_index"][:]) == [4, 5, 6, 7]
        assert list(f["seq_next"][:]) == [2, -9, 4, -9]
        assert list(f["bifurcation_upst"][:]) == [2] and list(f["bifurcation_down"][:]) == [4]


def test_domain_without_runoff_is_refused(tmp_path):
    make_mask(tmp_path / "mask.nc", [(8, 4)])
    with pytest.raises(SystemExit, match="no unit catchment receives runoff"):
        run(tmp_path, "--mask", str(tmp_path / "mask.nc"))


def make_restart(path, extra_dimension=False, size=N):
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("ucatch", size)
        f.createDimension("gridriver_ucatch_identity_field", 4)
        f.createVariable("gridriver_restart_schema", "i4")[...] = 2
        f.createVariable("acctime_rnof", "f8")[...] = 5.0
        ident = f.createVariable("gridriver_ucatch_identity", "f8",
                                 ("ucatch", "gridriver_ucatch_identity_field"))
        ident[:] = np.column_stack([np.full(size, 7.0), np.arange(size) + 1, np.arange(size) + 1,
                                    np.array((NEXT + [0] * size)[:size], dtype=float)])
        f.createVariable("wdsrf_ucat", "f8", ("ucatch",))[:] = 10.0 * (np.arange(size) + 1)
        if extra_dimension:
            f.createDimension("npthout", 3)
            f.createVariable("bifflw", "f8", ("npthout",))[:] = 1.0


def test_restart_is_cut_to_the_subset_and_identity_follows_the_new_network(tmp_path):
    make_mask(tmp_path / "mask.nc", [(1, 1)])
    sub = run(tmp_path, "--mask", str(tmp_path / "mask.nc"))
    make_restart(tmp_path / "rst.nc")
    tool.main(["restart", str(sub), str(tmp_path / "rst.nc"), str(tmp_path / "rst_sub.nc")])
    with netCDF4.Dataset(tmp_path / "rst_sub.nc") as f, netCDF4.Dataset(sub) as s:
        f.set_auto_mask(False)
        assert len(f.dimensions["ucatch"]) == 4
        assert list(f["wdsrf_ucat"][:]) == [10.0, 20.0, 30.0, 80.0]
        assert f["acctime_rnof"][...] == 5.0 and f["gridriver_restart_schema"][...] == 2
        ident = f["gridriver_ucatch_identity"][:]
        assert list(ident[:, 0]) == [7.0] * 4                 # version carried over
        assert list(ident[:, 1]) == list(s["seq_x"][:])
        assert list(ident[:, 2]) == list(s["seq_y"][:])
        assert list(ident[:, 3]) == [3.0, 3.0, -9.0, -10.0]   # renumbered seq_next


def test_restart_of_another_network_or_with_extra_state_is_refused(tmp_path):
    make_mask(tmp_path / "mask.nc", [(1, 1)])
    sub = run(tmp_path, "--mask", str(tmp_path / "mask.nc"))
    make_restart(tmp_path / "wrong_size.nc", size=5)
    with pytest.raises(SystemExit, match="restart has 5 unit catchments"):
        tool.main(["restart", str(sub), str(tmp_path / "wrong_size.nc"), str(tmp_path / "o1.nc")])
    make_restart(tmp_path / "with_bif.nc", extra_dimension=True)
    with pytest.raises(SystemExit, match="cannot be converted"):
        tool.main(["restart", str(sub), str(tmp_path / "with_bif.nc"), str(tmp_path / "o2.nc")])
