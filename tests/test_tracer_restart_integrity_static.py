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

    assert "MPI_LAND, p_comm_io" in probe
    assert "p_root, p_comm_group" in probe
    assert probe.index("MPI_LAND, p_comm_io") < probe.index("p_root, p_comm_group")
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
    assert "present_count > 0" in check
    assert "present_count < block_count" in check
    assert "CALL CoLM_stop" in check


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
    assert "present_count == 0 .and. .not. allow_missing" in check
    assert "ncio_inquire_varsize" in check
    assert "size(varsize) == expected_rank" in check
    assert "varsize(expected_rank) == pixelset%vecgs%vlen(iblk,jblk)" in check
    assert "varsize(1) == expected_dim1" in check
    assert "shape_error_count > 0" in check
    assert wrapper_1d.index("ncio_require_complete_vector_var") < wrapper_1d.index(
        "ncio_read_vector_real8_1d"
    )
    assert "expected_rank=1" in wrapper_1d
    assert "allow_missing=present(defval)" in wrapper_1d
    assert wrapper_2d.index("ncio_require_complete_vector_var") < wrapper_2d.index(
        "ncio_read_vector_real8_2d"
    )
    assert "expected_rank=2" in wrapper_2d
    assert "expected_dim1=ndim1" in wrapper_2d
    assert "allow_missing=present(defval)" in wrapper_2d


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
        assert (
            "ncio_read_vector => ncio_read_vector_complete" in body
        ), f"{routine} must reject variables present in only some restart blocks"
