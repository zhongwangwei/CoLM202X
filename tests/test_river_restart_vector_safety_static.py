from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text(
    encoding="utf-8"
)
VECTOR_IO = (ROOT / "main/HYDRO/MOD_Vector_ReadWrite.F90").read_text(
    encoding="utf-8"
)


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_restart_vectors_are_allocated_on_every_rank() -> None:
    allocator = routine(TIMEVARS, "allocate_GridRiverLakeTimeVars")

    for allocation in (
        "allocate (wdsrf_ucat (ncell_state))",
        "allocate (veloc_riv  (ncell_state))",
        "allocate (momen_riv  (ncell_state))",
        "allocate (acc_rnof_uc (ncell_state))",
        "allocate (volresv (nresv_state))",
    ):
        assert allocation in allocator


def test_zero_global_vectors_return_before_gather_or_write() -> None:
    writer = routine(VECTOR_IO, "vector_gather_and_write")
    mapped_writer = routine(VECTOR_IO, "vector_gather_map2grid_and_write")

    for body in (writer, mapped_writer):
        guard = body.index("IF (totalvlen <= 0) RETURN")
        assert guard < body.index("CALL vector_gather_to_master")


def test_gather_initializes_unassigned_global_slots_to_missing_value() -> None:
    gather = routine(VECTOR_IO, "vector_gather_to_master")
    assert "USE MOD_Vars_Global, only: spval" in gather
    assert "wdata = spval" in gather


def test_vector_reader_rejects_shape_and_address_mismatch() -> None:
    reader = routine(VECTOR_IO, "vector_read_and_scatter")
    read = reader.index("CALL ncio_read_serial (filein, varname, rdata)")
    length_check = reader.index("size(rdata) /= expected_length")
    address_check = reader.index("data_address(iwork)%val < 1 .or. &")
    scatter = reader.index("rcache = rdata(data_address(iwork)%val)")

    assert read < length_check < scatter
    assert read < address_check < scatter
    assert "restart vector length mismatch" in reader
    assert "restart vector address out of range" in reader


def test_each_restart_read_starts_from_explicit_defaults() -> None:
    reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")
    first_probe = reader.index("ncio_var_exist(file_restart, 'acctime_rnof'")

    for reset in (
        "acctime_rnof = 0._r8",
        "IF (allocated(acc_rnof_uc)) acc_rnof_uc = 0._r8",
        "IF (allocated(volwater_ucat)) volwater_ucat = 0._r8",
        "volwater_ucat_valid = .false.",
    ):
        assert reader.index(reset) < first_probe


def test_timevars_lifecycle_clears_deferred_restart_path() -> None:
    allocator = routine(TIMEVARS, "allocate_GridRiverLakeTimeVars")
    finalizer = routine(TIMEVARS, "deallocate_GridRiverLakeTimeVars")

    assert "gridriver_restart_file = ''" in allocator
    assert "gridriver_restart_file = ''" in finalizer
