from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text(
    encoding="utf-8"
)
TOPLEVEL = (ROOT / "main/MOD_Vars_TimeVariables.F90").read_text(
    encoding="utf-8"
)
RESERVOIR = (ROOT / "main/HYDRO/MOD_Grid_Reservoir.F90").read_text(
    encoding="utf-8"
)


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_identity_signature_covers_version_location_and_downstream() -> None:
    builder = routine(TIMEVARS, "build_gridriver_ucatch_identity")

    for field in (
        "GRIDRIVER_UCATCH_IDENTITY_VERSION",
        "x_ucat",
        "y_ucat",
        "ucat_next",
    ):
        assert field in builder
    assert "allocate (identity(4, numucat))" in builder


def test_writer_opens_transaction_and_persists_identity_before_state() -> None:
    writer = routine(TIMEVARS, "WRITE_GridRiverLakeTimeVars")

    created = writer.index("CALL ncio_create_file")
    schema = writer.index("'gridriver_restart_schema', GRIDRIVER_RESTART_SCHEMA_VERSION")
    incomplete = writer.index("'gridriver_restart_complete', 0")
    identity = writer.index("CALL write_gridriver_ucatch_identity")
    first_state = writer.index("wdsrf_ucat, numucat")

    assert created < schema < incomplete < identity < first_state


def test_reader_validates_transaction_and_identity_before_any_state() -> None:
    reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")

    transaction = reader.index("CALL validate_gridriver_restart_transaction")
    identity = reader.index("CALL validate_gridriver_ucatch_identity")
    first_state = reader.index("'wdsrf_ucat', ucat_data_address")
    assert transaction < identity < first_state


def test_zero_global_ucatch_skips_identity_matrix_and_required_vectors() -> None:
    identity_writer = routine(TIMEVARS, "write_gridriver_ucatch_identity")
    identity_reader = routine(TIMEVARS, "validate_gridriver_ucatch_identity")
    state_reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")

    assert "IF (totalnumucat <= 0) RETURN" in identity_writer
    assert "IF (totalnumucat <= 0) RETURN" in identity_reader
    zero_guard = state_reader.index("IF (totalnumucat > 0) THEN")
    first_state = state_reader.index("'wdsrf_ucat', ucat_data_address")
    second_state = state_reader.index("'veloc_riv',  ucat_data_address")
    assert zero_guard < first_state < second_state


def test_transaction_accepts_only_complete_new_or_all_absent_legacy() -> None:
    validator = routine(TIMEVARS, "validate_gridriver_restart_transaction")

    for name in (
        "gridriver_restart_schema",
        "gridriver_restart_complete",
        "gridriver_ucatch_identity",
    ):
        assert f"ncio_var_exist(file_restart, '{name}'" in validator

    assert "legacy GridRiverLake restart has no transaction or ucatch identity metadata" in validator
    assert "incomplete GridRiverLake restart transaction metadata" in validator
    assert "unsupported GridRiverLake restart schema" in validator
    assert "GridRiverLake restart transaction is not complete" in validator


def test_identity_is_read_collectively_and_compared_globally() -> None:
    validator = routine(TIMEVARS, "validate_gridriver_ucatch_identity")

    assert "CALL vector_read_matrix_and_scatter" in validator
    assert "'gridriver_ucatch_identity'" in validator
    assert "CALL build_gridriver_ucatch_identity" in validator
    assert "ieee_is_finite(identity_restart(:, i))" in validator
    assert "mpi_allreduce (MPI_IN_PLACE, mismatch_count" in validator
    assert "GridRiverLake restart ucatch identity mismatch" in validator


def test_reservoir_identity_binds_ordinal_state_to_dam_ucatch() -> None:
    builder = routine(TIMEVARS, "build_gridriver_reservoir_identity")
    writer = routine(TIMEVARS, "write_gridriver_reservoir_identity")
    validator = routine(TIMEVARS, "validate_gridriver_reservoir_identity")
    reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")

    assert "GRIDRIVER_RESERVOIR_IDENTITY_VERSION" in builder
    assert "resv_ucid" in builder
    assert "resv_global_id" in writer
    assert "gridriver_reservoir_identity" in writer
    assert "schema-v2 restart is missing required reservoir identity" in validator
    assert "ieee_is_finite(identity_restart(:, i))" in validator
    assert "different reservoir configuration" in validator
    assert reader.index("CALL validate_gridriver_reservoir_identity") < reader.index(
        "'volresv', resv_data_address"
    )


def test_reservoir_parameter_ordinals_are_total_and_one_to_one() -> None:
    initializer = routine(RESERVOIR, "reservoir_init")

    assert "ucat2resv = 0" in initializer
    assert "ordinal_count(resv_global_id(irsv))" in initializer
    assert "CALL mpi_allreduce" in initializer
    assert "any(ordinal_count /= 1)" in initializer
    assert "reservoir parameter rows must map one-to-one to active ucatch IDs" in initializer


def test_top_level_commits_only_after_all_gridriver_writers_return() -> None:
    writer = routine(TOPLEVEL, "WRITE_TimeVariables")

    state = writer.index("CALL WRITE_GridRiverLakeTimeVars")
    history = writer.index("CALL write_gridriverlake_hist_restart")
    tracer = writer.index("CALL write_tracer_restart(file_restart)")
    commit = writer.index("CALL commit_GridRiverLakeRestart")
    assert state < history < tracer < commit

    committer = routine(TIMEVARS, "commit_GridRiverLakeRestart")
    assert "CALL mpi_barrier" in committer
    assert "'gridriver_restart_complete', 1" in committer
