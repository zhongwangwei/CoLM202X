from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "share/MOD_NetCDFVector.F90").read_text()
METHANE = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane.F90").read_text()
ACCFLUX = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90").read_text()
MICROBES = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Microbes.F90").read_text()
STATE = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_State.F90").read_text()


def _routine(name: str) -> str:
    return SOURCE.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_complete_restart_preflight_result_drives_the_read_path():
    preflight = _routine("ncio_require_complete_vector_var")

    assert "logical, intent(out) :: all_present" in preflight
    assert "all_present = counts(1) > 0 .and. counts(2) == counts(1)" in preflight

    for rank in ("1d", "2d"):
        wrapper = _routine(f"ncio_read_vector_complete_real8_{rank}")
        reader = _routine(f"ncio_read_vector_real8_{rank}")

        assert "logical :: field_present" in wrapper
        assert "all_present=field_present" in wrapper
        assert "known_present=field_present" in wrapper
        assert "logical, intent(in), optional :: known_present" in reader
        assert "IF (present(known_present)) THEN" in reader
        assert "IF (known_present) THEN" in reader


def test_optional_feature_group_presence_is_batched_into_one_collective_scan():
    probe = _routine("ncio_vector_group_presence")

    assert "DO ivar = 1, size(datanames)" in probe
    assert "ncio_var_exist" in probe
    assert probe.count("mpi_allreduce") == 2
    assert "present_flags(ivar)" in probe
    assert "missing from some vector blocks" in probe


def test_strict_read_skips_duplicate_presence_probe_and_reduction():
    for rank in ("1d", "2d"):
        reader = _routine(f"ncio_read_vector_real8_{rank}")

        # The ordinary reader retains one probe for non-strict callers, while
        # strict callers enter the known-present branch before that fallback.
        assert reader.count("ncio_var_exist") == 1
        assert reader.index("IF (present(known_present)) THEN") < reader.index(
            "ncio_var_exist"
        )
        assert (
            "IF (.not. present(known_present)) THEN CALL mpi_allreduce"
            in " ".join(reader.split())
        )


def test_complete_restart_semantics_remain_strict_and_backward_compatible():
    preflight = _routine("ncio_require_complete_vector_var")

    assert "counts(2) > 0 .and. counts(2) < counts(1)" in preflight
    assert "counts(2) == 0 .and. .not. allow_missing" in preflight
    assert "counts(3) > 0" in preflight
    assert "size(varsize) == expected_rank" in preflight
    assert "varsize(expected_rank) == pixelset%vecgs%vlen(iblk,jblk)" in preflight
    assert "block_shape_ok = varsize(1) == expected_dim1" in preflight

    for rank in ("1d", "2d"):
        wrapper = _routine(f"ncio_read_vector_complete_real8_{rank}")
        reader = _routine(f"ncio_read_vector_real8_{rank}")

        assert "allow_missing=allow_field_missing" in wrapper
        assert (
            "allow_field_missing = present(defval) .and. "
            "(.not. complete_require_present)"
            in " ".join(wrapper.split())
        )
        assert (
            "rdata, & defval, known_present=field_present"
            in " ".join(wrapper.split())
        )
        assert "ELSEIF (present(defval)) THEN" in reader
        assert "sbuff" in reader and "= defval" in reader


def test_committed_methane_schema_requires_every_requested_restart_field():
    callback = _routine_from(METHANE, "ch4_reactive_read_restart")
    validator = _routine_from(METHANE, "validate_methane_restart_transaction")
    acc_reader = _routine_from(ACCFLUX, "read_methane_accflux_restart")

    assert "logical, private, save :: complete_require_present = .false." in SOURCE
    assert "SUBROUTINE ncio_set_complete_require_present" in SOURCE
    assert "logical, intent(out) :: strict_restart" in validator
    assert "strict_restart = .false." in validator
    assert "strict_restart = .true." in validator
    assert "ncio_set_complete_require_present (strict_restart)" in callback
    assert "read_methane_restart (file_restart, strict_restart, restart_schema)" in callback
    assert callback.index("CALL ncio_vector_group_presence") < callback.index(
        "ncio_set_complete_require_present (strict_restart)"
    )
    assert callback.index("ncio_set_complete_require_present (strict_restart)") < callback.index(
        "read_methane_restart"
    )
    assert (
        "read_methane_accflux_restart (file_restart, file_has_pools, &"
        in callback
    )
    assert "IF (file_has_pools) THEN" in callback
    assert "IF (DEF_METHANE%use_microbial_pools) THEN" in callback
    assert "validate_methane_microbes_restart_values" in callback
    assert callback.index("read_methane_microbes_restart") < callback.index(
        "ncio_set_complete_require_present (.false.)"
    )
    assert "logical, intent(in), optional :: checkpoint_has_microbe_pools" in acc_reader
    assert "logical, intent(in), optional :: checkpoint_has_microbe_accumulators" in acc_reader
    assert "read_microbe_fields = checkpoint_has_microbe_pools" in acc_reader
    assert "has_microbe_accumulator .neqv. read_microbe_fields" in acc_reader
    assert "IF (read_microbe_fields) THEN" in acc_reader
    assert "IF (allocated(a_methane_acc_num_microbe)) THEN" in acc_reader
    assert "discarded_microbe_accumulator_2d" in acc_reader
    assert "discarded_microbe_accumulator_1d" in acc_reader


def test_ignored_microbe_restart_validator_uses_zero_extent_on_nonworkers():
    validator = _routine_from(
        MICROBES, "validate_methane_microbes_restart_values"
    )

    assert "IF (p_is_worker) THEN" in validator
    assert "allocate(values(nl_soil, landpatch%nset))" in validator
    assert "allocate(values(nl_soil, 0))" in validator


def test_restart_writer_invalidates_commit_marker_before_payload():
    writer = _routine_from(METHANE, "ch4_reactive_write_restart")

    invalidate = "'ch4_restart_complete', 0._r8"
    schema = "'ch4_restart_schema'"
    commit = "'ch4_restart_complete', 1._r8"
    assert invalidate in writer
    assert commit in writer
    assert writer.index(invalidate) < writer.index(schema) < writer.index(commit)


def test_legacy_restart_without_lake_soilc_preserves_surface_initialization():
    reader = _routine_from(STATE, "read_methane_restart")
    normalized = " ".join(reader.split())

    assert "lake_soilc_present = strict_restart_active" in reader
    assert (
        "IF (.not. strict_restart_active) lake_soilc_present = & "
        "ncio_vector_var_present(file_restart, 'ch4_lake_soilc', landpatch)"
        in normalized
    )
    assert (
        "IF (lake_soilc_present) CALL ncio_read_vector (file_restart, "
        "'ch4_lake_soilc'"
        in normalized
    )
    assert "'ch4_lake_soilc',       nl_soil, landpatch, lake_soilc,       defval = 0._r8" not in reader


def _routine_from(text: str, name: str) -> str:
    return text.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]
