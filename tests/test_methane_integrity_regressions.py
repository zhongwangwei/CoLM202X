from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def source(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def routine_body(text: str, name: str, kind: str = "SUBROUTINE") -> str:
    match = re.search(
        rf"{kind}\s+{name}\b(?P<body>.*?)END\s+{kind}\s+{name}",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    assert match is not None, f"{kind.lower()} {name} not found"
    return match.group("body")


def test_committed_restart_rejects_negative_inventory_and_legacy_warns() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    reader = routine_body(state, "read_methane_restart")

    assert "corrupt_prognostic_values" in reader
    assert "IF (strict_restart_active) THEN" in reader
    assert "invalid or negative prognostic state values" in reader
    assert "legacy methane restart contains" in reader
    assert "CALL CoLM_stop()" in reader
    for field in (
        "conc_o2",
        "conc_o2_unsat",
        "conc_o2_sat",
        "conc_o2_lake",
        "conc_methane",
        "conc_methane_unsat",
        "conc_methane_sat",
        "conc_methane_lake",
        "totcol_methane",
        "totcol_methane_unsat",
        "totcol_methane_sat",
        "totcol_methane_lake",
        "lake_water_ch4_stock",
        "lake_water_o2_stock",
        "lake_soilc",
    ):
        assert f"invalid_restart_value({field})" in reader
    assert "lake_water_ch4_present = ncio_vector_var_present" in reader
    assert "lake_water_o2_present = ncio_vector_var_present" in reader
    for field in (
        "grnd_methane_cond",
        "annavg_agnpp",
        "annavg_bgnpp",
        "annavg_somhr",
        "tempavg_agnpp",
        "tempavg_bgnpp",
        "tempavg_somhr",
        "tempavg_finrw",
        "annsum_counter",
        "methane_dfsat_tot",
        "f_h2osfc",
        "f_inund_levee_patch",
        "f_inund_flood_patch",
        "f_inund_flood_depth_patch",
    ):
        assert f"invalid_restart_value({field})" in reader
    for field in ("annavg_finrw", "fsat_bef", "finundated_lag", "layer_sat_lag"):
        assert f"invalid_restart_fraction_or_sentinel({field})" in reader
    assert "restart_ch4_clip_credit_mass" not in state
    assert "credit_methane_restart_clip" not in state


def test_restart_lag_states_reject_nan_before_flux_arithmetic() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    physics = source("MOD_Tracer_Reactive_Methane_Physics.F90")
    reader = routine_body(state, "read_methane_restart")
    methane = routine_body(physics, "methane")

    for field in ("fsat_bef", "finundated_lag", "layer_sat_lag"):
        assert re.search(
            rf"WHERE\s*\(\s*invalid_restart_fraction_or_sentinel\s*\(\s*{field}\s*\)\s*\)\s*"
            rf"{field}\s*=\s*spval",
            reader,
            flags=re.IGNORECASE,
        )
        assert re.search(
            rf"ieee_is_nan\s*\(\s*{field}(?:\(j\))?\s*\)",
            methane,
            flags=re.IGNORECASE,
        )

    assert methane.index("ieee_is_nan(fsat_bef)") < methane.index("dfsat =")
    assert methane.index("ieee_is_nan(finundated_lag)") < methane.index(
        "finundated_lag = finundated_lag *"
    )


def test_lake_restart_total_and_components_are_consistency_checked() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    reader = routine_body(state, "read_methane_restart")

    assert "logical, intent(in), optional :: strict_restart" in reader
    assert "lake_component_total = totcol_methane_lake(ipatch) +" in reader
    assert "lake_water_ch4_stock(ipatch)" in reader
    assert re.search(
        r"lake_inventory_tolerance\s*=\s*1\.e-12_r8\s*\+\s*1\.e-10_r8",
        reader,
        flags=re.IGNORECASE,
    )
    assert "invalid_lake_inventory = 1" in reader
    assert "totcol_methane(ipatch) = lake_component_total" in reader
    assert "fsat_bef(ipatch) = spval" in reader
    assert "MPI_INTEGER, MPI_SUM" in reader
    assert "CALL CoLM_stop()" in reader

    def inconsistent(total: float, sediment: float, water: float) -> bool:
        components = sediment + water
        tolerance = 1.0e-12 + 1.0e-10 * max(abs(total), abs(components))
        return abs(total - components) > tolerance

    assert not inconsistent(0.100000000005, 0.08, 0.02)
    assert inconsistent(0.12, 0.08, 0.02)


def test_patch_ph_rejects_nan_and_infinity() -> None:
    ph = source("MOD_Tracer_Reactive_Methane_pH.F90")
    reader = routine_body(ph, "read_methane_ph_patch")

    assert "ieee_is_finite" in ph
    assert re.search(
        r"\.not\.\s*ieee_is_finite\s*\(\s*methane_ph_patch\s*\(\s*ipatch\s*\)\s*\)",
        reader,
        flags=re.IGNORECASE,
    )


def test_rice_calendar_uses_current_and_previous_leap_year_lengths() -> None:
    bgc = source("MOD_Tracer_Reactive_Methane_BgcLink.F90")

    for name in ("rice_days_past_planting", "rice_days_since_harvest"):
        body = routine_body(bgc, name, kind="FUNCTION")
        assert re.search(r"integer\s*,\s*intent\(in\)\s*,\s*optional\s*::\s*year", body, re.I)
        assert "isleapyear(year)" in body
        assert "isleapyear(year - 1)" in body
        assert "event_dayspyr + jday" in body
        assert "integer, parameter :: dayspyr = 365" not in body


def test_dormancy_transfer_preserves_active_seed_and_total_pool_is_bounded() -> None:
    microbes = source("MOD_Tracer_Reactive_Methane_Microbes.F90")
    step = routine_body(microbes, "methane_microbes_step")
    cap = routine_body(microbes, "cap_microbe_biomass_layer")

    assert re.search(
        r"B_methanogen_comp\(j,component,ipatch\)\s*-\s*&?\s*DEF_METHANE%B_min_methanogen",
        step,
        flags=re.IGNORECASE,
    )
    assert re.search(
        r"B_methanotroph_comp\(j,component,ipatch\)\s*-\s*&?\s*DEF_METHANE%B_min_methanotroph",
        step,
        flags=re.IGNORECASE,
    )
    assert step.count("CALL cap_microbe_biomass_layer") == 2
    assert "CALL cap_microbe_biomass_layer" in step.split("CYCLE", 1)[0]
    assert "seed_total" in cap
    assert "allowed_surplus" in cap
    assert re.search(r"scale\s*=\s*allowed_surplus\s*/\s*surplus", cap, re.I)


def test_uncoupled_rice_substrate_multiplier_is_hard_disabled() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    physics = source("MOD_Tracer_Reactive_Methane_Physics.F90")

    assert re.search(
        r"abs\s*\(\s*DEF_METHANE%rice_substrate_boost\s*-\s*1\._r8\s*\)",
        const,
        flags=re.IGNORECASE,
    )
    assert "rice_substrate_boost_eff" not in physics
    assert "methane_prod_depth_unsat(:) * rice_substrate_boost" not in physics
    assert "methane_prod_depth_sat(:) * rice_substrate_boost" not in physics
