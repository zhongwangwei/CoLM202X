from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _source(relative: str) -> str:
    return (ROOT / relative).read_text()


def test_land_restart_dimension_probe_is_global_before_vector_reads():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    probe = source.split("logical FUNCTION tracer_dim_matches", 1)[1].split(
        "END FUNCTION tracer_dim_matches", 1
    )[0]
    reader = source.split("SUBROUTINE read_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_land_tracer_restart", 1
    )[0]

    assert "MPI_SUM, p_comm_io" in probe
    assert "p_root, p_comm_group" in probe
    assert probe.index("MPI_SUM, p_comm_io") < probe.index("p_root, p_comm_group")
    assert "size(varsize) == expected_rank" in probe
    assert "varsize(expected_rank) == landpatch%vecgs%vlen(iblk,jblk)" in probe
    assert "expected_rank = merge(3, 2, present(expect_soilsnow))" in probe
    assert reader.index("tracer_dim_matches(file_restart, 'trc_ldew_rain'") < reader.index(
        "CALL ncio_read_vector(file_restart, 'trc_ldew_rain'"
    )
    # A collective-bearing function must not be hidden in a short-circuitable
    # OR chain: every IO/worker rank must execute the same probe sequence.
    required_probe = reader.split("! Reject the restart", 1)[1].split(
        "CALL ncio_read_vector(file_restart, 'trc_ldew_rain'", 1
    )[0]
    assert ".or." not in required_probe.lower()


def test_complete_vector_reader_rejects_partial_block_variables():
    source = _source("share/MOD_NetCDFVector.F90")

    assert "INTERFACE ncio_read_vector_complete" in source
    assert "ncio_require_complete_vector_var" in source
    check = source.split("SUBROUTINE ncio_require_complete_vector_var", 1)[1].split(
        "END SUBROUTINE ncio_require_complete_vector_var", 1
    )[0]
    assert "ncio_var_exist" in check
    assert "MPI_SUM, p_comm_io" in check
    assert "p_root, p_comm_group" in check
    assert "counts(2) > 0" in check
    assert "counts(2) < counts(1)" in check
    assert "CALL CoLM_stop" in check


def test_complete_vector_preflight_skips_the_control_master_singleton_group():
    source = _source("share/MOD_NetCDFVector.F90")
    check = source.split("SUBROUTINE ncio_require_complete_vector_var", 1)[1].split(
        "END SUBROUTINE ncio_require_complete_vector_var", 1
    )[0]

    assert "p_is_worker" in check
    assert "IF (.not. (p_is_io .or. p_is_worker)) RETURN" in check


def test_complete_vector_reader_rejects_missing_required_and_wrong_shape_before_read():
    source = _source("share/MOD_NetCDFVector.F90")
    check = source.split("SUBROUTINE ncio_require_complete_vector_var", 1)[1].split(
        "END SUBROUTINE ncio_require_complete_vector_var", 1
    )[0]
    wrapper_1d = source.split("SUBROUTINE ncio_read_vector_complete_real8_1d", 1)[1].split(
        "END SUBROUTINE ncio_read_vector_complete_real8_1d", 1
    )[0]
    wrapper_2d = source.split("SUBROUTINE ncio_read_vector_complete_real8_2d", 1)[1].split(
        "END SUBROUTINE ncio_read_vector_complete_real8_2d", 1
    )[0]

    assert "allow_missing" in check
    assert "counts(2) == 0 .and. .not. allow_missing" in check
    assert "ncio_inquire_varsize" in check
    assert "size(varsize) == expected_rank" in check
    assert "varsize(expected_rank) == pixelset%vecgs%vlen(iblk,jblk)" in check
    assert "varsize(1) == expected_dim1" in check
    assert "counts(3) > 0" in check
    assert wrapper_1d.index("ncio_require_complete_vector_var") < wrapper_1d.index(
        "ncio_read_vector_real8_1d"
    )
    assert "expected_rank=1" in wrapper_1d
    assert "allow_missing=allow_field_missing" in wrapper_1d
    assert wrapper_2d.index("ncio_require_complete_vector_var") < wrapper_2d.index(
        "ncio_read_vector_real8_2d"
    )
    assert "expected_rank=2" in wrapper_2d
    assert "expected_dim1=ndim1" in wrapper_2d
    assert "allow_missing=allow_field_missing" in wrapper_2d


def test_all_ch4_restart_readers_opt_into_complete_block_policy():
    for relative, routine in (
        ("main/TRACER/MOD_Tracer_Reactive_Methane_State.F90", "read_methane_restart"),
        ("main/TRACER/MOD_Tracer_Reactive_Methane_Microbes.F90", "read_methane_microbes_restart"),
        ("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90", "read_methane_accflux_restart"),
    ):
        source = _source(relative)
        body = source.split(f"SUBROUTINE {routine}", 1)[1].split(
            f"END SUBROUTINE {routine}", 1
        )[0]
        if routine == "read_methane_accflux_restart":
            assert "INTERFACE ncio_read_vector" in source
            assert source.count(
                "USE MOD_NetCDFVector, only: ncio_read_vector_complete"
            ) >= 2
        else:
            assert (
                "ncio_read_vector => ncio_read_vector_complete" in body
            ), f"{routine} must reject variables present in only some restart blocks"


def test_ch4_restart_uses_schema_in_progress_and_commit_markers():
    source = _source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    reader = source.split("SUBROUTINE ch4_reactive_read_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_read_restart", 1
    )[0]
    writer = source.split("SUBROUTINE ch4_reactive_write_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_write_restart", 1
    )[0]
    validator = source.split("SUBROUTINE validate_methane_restart_transaction", 1)[1].split(
        "END SUBROUTINE validate_methane_restart_transaction", 1
    )[0]

    assert "ch4_restart_schema" in validator
    assert "ch4_restart_complete" in validator
    assert "#ifdef USEMPI" in validator
    assert "#else\n      USE MOD_SPMD_Task, only: p_is_worker, p_is_master" in validator
    assert "IF (.not. has_schema .and. .not. has_commit) RETURN" in validator
    assert "IF (.not. has_schema .or. .not. has_commit) THEN" in validator
    assert "logical, intent(out) :: strict_restart" in validator
    assert "strict_restart = .true." in validator
    assert "ncio_set_complete_require_present (strict_restart)" in reader
    assert "ncio_set_complete_require_present (.false.)" in reader
    assert "CALL ncio_vector_group_presence" in reader
    assert "n_microbe_state_fields = count(microbe_field_present(1:4))" in reader
    assert "partial/inconsistent microbial feature group" in reader
    assert "IF (file_has_pools) THEN" in reader
    assert "IF (DEF_METHANE%use_microbial_pools) THEN" in reader
    assert "validate_methane_microbes_restart_values" in reader
    assert reader.index("validate_methane_restart_transaction") < reader.index(
        "read_methane_restart"
    )
    assert writer.count("write_methane_restart_marker") == 4
    assert writer.index("write_methane_restart_marker(file_restart, 'ch4_restart_complete', 0._r8") < writer.index(
        "write_methane_restart(file_restart"
    )
    assert writer.index("write_methane_restart_marker(file_restart, 'ch4_history_accumulation_mode'") < writer.index(
        "write_methane_restart(file_restart"
    )
    assert writer.rindex("write_methane_restart_marker(file_restart, 'ch4_restart_complete', 1._r8") > writer.index(
        "write_methane_microbes_restart(file_restart"
    )


def test_ch4_presence_probes_use_block_aware_vector_metadata():
    methane = _source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    accflux = _source("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90")
    reader = methane.split("SUBROUTINE ch4_reactive_read_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_read_restart", 1
    )[0]
    acc_reader = accflux.split("SUBROUTINE read_methane_accflux_restart", 1)[1].split(
        "END SUBROUTINE read_methane_accflux_restart", 1
    )[0]

    assert "ncio_vector_group_presence" in reader
    assert "ncio_var_exist(file_restart" not in reader
    assert "ncio_vector_var_present" in acc_reader
    assert "ncio_var_exist(file_restart" not in acc_reader
