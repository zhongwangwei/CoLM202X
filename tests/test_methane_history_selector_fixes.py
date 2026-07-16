from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Hist.F90"
).read_text(encoding="utf-8")
CONST_SOURCE = (
    ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90"
).read_text(encoding="utf-8")
ACCFLUX_SOURCE = (
    ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90"
).read_text(encoding="utf-8")
METHANE_SOURCE = (
    ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane.F90"
).read_text(encoding="utf-8")


def test_custom_history_selector_fails_when_no_real_query_matches():
    assert "methane_history_match_count = 0" in SOURCE
    assert "IF (mhist_on) methane_history_match_count = methane_history_match_count + 1" in SOURCE
    assert "ch4_history_vars selects no known methane history variable" in SOURCE
    assert "USE MOD_Tracer_Reactive_Methane_Const, only: methane_history_enabled" in SOURCE


def test_history_match_tracking_wraps_the_existing_selector():
    start = SOURCE.index("logical FUNCTION mhist_on (varname)")
    end = SOURCE.index("END FUNCTION mhist_on", start)
    wrapper = SOURCE[start:end]
    assert "mhist_on = methane_history_enabled(varname)" in wrapper

    # The wrapper updates module state, so keep at most one invocation in
    # each Fortran statement (evaluation order between function references
    # is otherwise unspecified).
    flattened = re.sub(r"&\s*\n\s*", " ", SOURCE)
    assert all(line.count("mhist_on(") <= 1 for line in flattened.splitlines())


def test_every_writable_history_selector_reaches_an_accumulation_gate():
    pattern = re.compile(r"mhist_on\('([^']+)'\)")
    writable = set(pattern.findall(SOURCE))
    accumulated = set(pattern.findall(ACCFLUX_SOURCE))

    assert writable == accumulated


def test_restart_resets_partial_window_when_accumulation_mode_changes():
    mode = CONST_SOURCE.split(
        "integer FUNCTION methane_history_accumulation_mode", 1
    )[1].split("END FUNCTION methane_history_accumulation_mode", 1)[0]
    assert "IF (.not. DEF_METHANE%write_ch4_history) RETURN" in mode
    assert "CASE ('none', 'off', 'false', '.false.')" in mode
    assert "CASE ('core', 'default', 'minimal', 'fast')" in mode
    assert "methane_history_accumulation_mode = 2" in mode
    assert "history_mode = methane_history_accumulation_mode ()" in ACCFLUX_SOURCE

    writer = METHANE_SOURCE.split("SUBROUTINE ch4_reactive_write_restart", 1)[
        1
    ].split("END SUBROUTINE ch4_reactive_write_restart", 1)[0]
    reader = METHANE_SOURCE.split("SUBROUTINE ch4_reactive_read_restart", 1)[
        1
    ].split("END SUBROUTINE ch4_reactive_read_restart", 1)[0]
    validator = METHANE_SOURCE.split(
        "SUBROUTINE validate_methane_restart_transaction", 1
    )[1].split("END SUBROUTINE validate_methane_restart_transaction", 1)[0]

    assert "'ch4_history_accumulation_mode'" in writer
    assert "ncio_vector_group_presence" in validator
    assert "restart_metadata_fields(4)" in validator
    assert "ncio_read_vector_complete" in validator
    assert "history_mode_changed = .not. has_history_mode" in validator
    assert "MPI_LOR" in validator
    assert "restart_schema >= 2 .and. history_mode_changed" in reader
    assert "CALL flush_methane_acc_fluxes ()" in reader
