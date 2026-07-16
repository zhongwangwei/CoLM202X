from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90"
).read_text(encoding="utf-8")


def routine(name: str) -> str:
    return SOURCE.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def normalized(text: str) -> str:
    return " ".join(text.replace("&", " ").split())


def test_every_persisted_accumulator_uses_the_validating_read_generic() -> None:
    declarations = SOURCE.split("CONTAINS", 1)[0]
    writer = normalized(routine("write_methane_accflux_restart"))
    reader = normalized(routine("read_methane_accflux_restart"))

    declared = set(
        re.findall(
            r"real\(r8\),\s*allocatable\s*::\s*(a_[A-Za-z0-9_]+)",
            declarations,
            flags=re.IGNORECASE,
        )
    )
    written = set(
        re.findall(
            r"ncio_write_vector\s*\([^)]*?\b(a_[A-Za-z0-9_]+)\s*,\s*compress\s*\)",
            writer,
            flags=re.IGNORECASE,
        )
    )
    read = set(
        re.findall(
            r"ncio_read_vector\s*\([^)]*?\b(a_[A-Za-z0-9_]+)\s*,\s*defval\s*=",
            reader,
            flags=re.IGNORECASE,
        )
    )

    # Do not freeze a field count: the safety property is complete symmetry
    # between every declared, written, and validating-read accumulator.
    assert declared
    assert written == declared
    assert read == declared
    assert "INTERFACE ncio_read_vector" in SOURCE
    assert "ncio_read_vector => ncio_read_vector_complete" not in routine(
        "read_methane_accflux_restart"
    )

    # Optional microbial diagnostics are still read into one discard buffer
    # when the current configuration disables their runtime allocations.
    for field in (
        "ch4_a_B_methanogen",
        "ch4_a_B_methanotroph",
        "ch4_a_B_methanogen_dormant",
        "ch4_a_B_methanotroph_dormant",
        "ch4_a_f_T_methanogen",
        "ch4_a_f_S_methanogen",
        "ch4_a_f_O2_methanogen",
        "ch4_a_f_T_methanotroph",
        "ch4_a_methanogen_growth_rate",
        "ch4_a_methanotroph_growth_rate",
        "ch4_a_microbial_prod_potential",
        "ch4_a_microbial_oxid_potential",
        "ch4_a_methane_acc_num_microbe",
    ):
        expected = 3 if field == "ch4_a_methane_acc_num_microbe" else 2
        assert reader.count(f"'{field}'") == expected
    assert "discarded_microbe_accumulator_2d" in reader
    assert "discarded_microbe_accumulator_1d" in reader


def test_both_restart_vector_ranks_reject_nonfinite_and_sentinel_values() -> None:
    predicate = SOURCE.split(
        "ELEMENTAL LOGICAL FUNCTION invalid_accflux_restart_value", 1
    )[1].split("END FUNCTION invalid_accflux_restart_value", 1)[0]

    assert "ieee_is_finite(x)" in predicate
    assert "abs(x) >= 0.5_r8 * abs(spval)" in normalized(predicate)
    for name in (
        "read_accflux_restart_vector_1d",
        "read_accflux_restart_vector_2d",
    ):
        wrapper = routine(name)
        assert wrapper.index("CALL ncio_read_vector_complete") < wrapper.index(
            "any(invalid_accflux_restart_value(rdata))"
        )
        assert "p_is_worker" in wrapper


def test_committed_schema_validation_is_collective_and_precedes_any_reset() -> None:
    reader = routine("read_methane_accflux_restart")
    compact = normalized(reader)

    for counter in (
        "a_methane_acc_num",
        "a_methane_acc_num_unsat",
        "a_methane_acc_num_sat",
        "a_methane_acc_num_lake",
        "a_methane_acc_num_extra",
        "a_annavg_finrw_acc_num",
        "a_methane_acc_num_microbe",
    ):
        assert re.search(rf"any\({counter}\s*<\s*0\._r8\)", compact)

    assert "MPI_LOGICAL" in reader
    assert "MPI_LOR" in reader
    assert "CALL mpi_allreduce" in reader
    assert "any(discarded_microbe_accumulator_1d < 0._r8)" in compact
    strict_guard = (
        "strict_restart_active .and. restart_schema_active >= 1 .and. "
        "corrupt_accumulator_values"
    )
    assert strict_guard in compact
    assert reader.index("CALL mpi_allreduce") < reader.index(
        "IF (strict_restart_active .and. restart_schema_active >= 1"
    )
    assert reader.index("CALL CoLM_stop()") < reader.index(
        "CALL flush_methane_acc_fluxes ()"
    )


def test_legacy_and_non_strict_corruption_start_a_clean_history_window() -> None:
    reader = normalized(routine("read_methane_accflux_restart"))

    assert "IF (restart_schema_active < 3) THEN" in reader
    assert "legacy methane restart resets the in-progress CH4 history accumulation window" in reader
    assert "ELSEIF (corrupt_accumulator_values) THEN" in reader
    assert "non-strict methane restart resets corrupt CH4 history accumulators" in reader
    assert "history selection changed across restart" in reader
