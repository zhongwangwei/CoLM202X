from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
STATE = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_State.F90").read_text(
    encoding="utf-8"
)
CONST = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90").read_text(
    encoding="utf-8"
)


def routine(text: str, name: str) -> str:
    return text.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


CONDUCTANCES = (
    "grnd_methane_cond",
    "grnd_methane_cond_unsat",
    "grnd_methane_cond_sat",
    "grnd_methane_cond_lake",
)


def test_inactive_reset_is_closed_under_committed_restart_domain() -> None:
    reset = routine(STATE, "reset_methane_inactive_lake_diagnostics")
    writer = routine(STATE, "write_methane_restart")
    reader = routine(STATE, "read_methane_restart")

    for field in CONDUCTANCES:
        assert re.search(
            rf"{field}\(ipatch\)\s*=\s*DEF_METHANE%grnd_methane_cond_default",
            reset,
            flags=re.IGNORECASE,
        )
        assert f"'ch4_{field}'" in writer
        assert re.search(
            rf"invalid_restart_value\({field}\)\s*\.or\.\s*{field}\s*<=\s*0\._r8",
            reader,
            flags=re.IGNORECASE,
        )


def test_reset_value_satisfies_the_strict_reader_semantically() -> None:
    match = re.search(
        r"grnd_methane_cond_default\s*=\s*([0-9.eE+-]+)_r8", CONST
    )
    assert match is not None
    default_conductance = float(match.group(1))

    # This is the exact finite/positive domain enforced for committed state.
    assert default_conductance > 0.0
    assert default_conductance < float("inf")
    reset_values = {field: default_conductance for field in CONDUCTANCES}
    assert all(value > 0.0 and value < float("inf") for value in reset_values.values())


def test_namelist_validation_guards_the_reset_source_value() -> None:
    validator = routine(CONST, "validate_methane_namelist")
    assert re.search(
        r"DEF_METHANE%grnd_methane_cond_default\s*<=\s*0\._r8",
        validator,
        flags=re.IGNORECASE,
    )
    assert "grnd_methane_cond_default must be > 0" in validator
