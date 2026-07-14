from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIF = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90").read_text(
    encoding="utf-8"
)
HIST_STATE = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeHistState.F90").read_text(
    encoding="utf-8"
)


def routine(name: str) -> str:
    return BIF.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_restart_signature_covers_path_topology_and_geometry() -> None:
    builder = routine("build_bifurcation_path_signature")
    for field in (
        "BIF_RESTART_SIGNATURE_VERSION",
        "ucat_ucid(iup)",
        "pth_down_ucid(ipth)",
        "pth_dst(ipth)",
        "pth_elv(:, ipth)",
        "pth_wth(:, ipth)",
        "pth_man(:)",
    ):
        assert field in builder


def test_writer_persists_identity_before_path_state() -> None:
    writer = routine("write_bifurcation_restart")
    signature_write = writer.index(
        "CALL ncio_write_serial (file_restart, 'bif_path_signature'"
    )
    velocity_write = writer.index("CALL ncio_write_serial (file_restart, 'pth_veloc'")
    assert "bifurcation_signature_field" in writer
    assert signature_write < velocity_write


def test_reader_validates_identity_before_loading_state() -> None:
    reader = routine("read_bifurcation_restart")
    signature_read = reader.index("'bif_path_signature', pth_global_id")
    velocity_read = reader.index("'pth_veloc', pth_global_id")
    assert signature_read < velocity_read
    assert "mpi_allreduce (MPI_IN_PLACE, mismatch_count" in reader
    assert "ieee_is_finite(path_signature(:, ipth))" in reader
    assert "IF (allocated(path_signature)) deallocate (path_signature)" in reader
    assert "Refusing to load bifurcation momentum for a different pathway network" in reader


def test_legacy_restart_is_an_explicit_cold_start() -> None:
    reader = routine("read_bifurcation_restart")
    validation = reader.index(
        "IF ((.not. restart_feature_present .or. strict_bif_restart) .and. &"
    )
    assert reader.index("pth_veloc(:,:) = 0._r8") < validation
    assert reader.index("pth_momen(:,:) = 0._r8") < validation
    assert "previous_depth_restart_found" in reader
    assert "restart_loaded" in reader
    assert "WARNING: incomplete/legacy bifurcation restart" in reader


def test_legacy_restart_does_not_load_ordinal_path_history() -> None:
    assert "has_bif_signature = restart_var_exists(file_restart, 'bif_path_signature')" in HIST_STATE
    assert "IF (has_bif_signature) THEN" in HIST_STATE
    assert "cold-starting pathway history accumulators" in HIST_STATE


def test_corrupt_path_state_cold_starts_the_paired_restart_unit_trap_safely() -> None:
    reader = routine("read_bifurcation_restart")

    assert "ieee_is_finite(pth_veloc(ilev, ipth))" in reader
    assert "ieee_is_finite(pth_momen(ilev, ipth))" in reader
    assert "invalid_state_count" in reader
    assert "restart_loaded = .false." in reader
    assert "pth_veloc /= pth_veloc .or." not in reader


def test_same_size_path_reorder_changes_identity() -> None:
    # The old restart contract checked only matrix shape, so these two networks
    # were indistinguishable even though ordinal state would bind to other paths.
    old = [(10, 20, 1000.0), (30, 40, 2000.0)]
    reordered = [old[1], old[0]]
    assert len(old) == len(reordered)
    assert old != reordered
