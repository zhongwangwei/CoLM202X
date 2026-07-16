from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NETWORK = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeNetwork.F90").read_text(
    encoding="utf-8"
)


def routine(name: str) -> str:
    return NETWORK.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_bif_width_and_layer_geometry_contracts_are_explicit() -> None:
    reader = routine("read_bifurcation_global_arrays")
    assert "any(bif_wdth_all < 0._r8)" in reader
    assert "bifurcation width must be non-negative" in reader
    assert "active bifurcation layers must be contiguous from layer 1" in reader
    assert (
        "active bifurcation layer elevation must be non-decreasing from layer 1"
        in reader
    )
    assert "bif_wdth_all(ilev-1, ip) <= 0._r8" in reader
    assert "bif_elev_all(ilev, ip) < bif_elev_all(ilev-1, ip)" in reader


def test_equal_active_layer_elevations_remain_valid() -> None:
    elevations = [12.0, 12.0, 13.5]
    widths = [30.0, 20.0, 5.0]
    assert all(width >= 0.0 for width in widths)
    assert all(width > 0.0 for width in widths)
    assert all(
        elevations[level] >= elevations[level - 1]
        for level in range(1, len(elevations))
    )


def test_network_final_releases_owned_arrays_and_river_communicator() -> None:
    finalizer = routine("riverlake_network_final")
    assert "IF (allocated(wts_ups" in finalizer
    assert "IF (allocated(lake_type" in finalizer
    assert "IF (p_comm_rivsys /= MPI_COMM_NULL) THEN" in finalizer
    assert "CALL mpi_comm_free (p_comm_rivsys, p_err)" in finalizer
    assert "p_comm_rivsys = MPI_COMM_NULL" in finalizer

    declaration = NETWORK.split("CONTAINS", 1)[0]
    assert "integer :: p_comm_rivsys = MPI_COMM_NULL" in declaration


def test_network_final_releases_owned_mapping_and_grid_objects() -> None:
    finalizer = routine("riverlake_network_final")

    for push in (
        "push_inpm2ucat",
        "push_ucat2inpm",
        "push_ucat2grid",
        "allreduce_inpm",
        "push_next2ucat",
        "push_ups2ucat",
        "push_bif_dn2pth",
        "push_bif_influx",
    ):
        assert finalizer.count(f"CALL worker_pushdata_free_mem ({push})") == 1

    assert (
        finalizer.count(
            "CALL worker_remapdata_free_mem (remap_patch2inpm)"
        )
        == 1
    )
    assert finalizer.count("CALL grid_free_mem (griducat)") == 1
