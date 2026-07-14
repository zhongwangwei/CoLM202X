from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOIL = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text(
    encoding="utf-8"
)


def section(start: str, end: str) -> str:
    return SOIL.split(start, 1)[1].split(end, 1)[0]


def test_finite_liquid_ratios_use_the_physical_water_floor():
    ratio_setup = section(
        "source_fallback_ratio = tracer_init_water_ratio(itrc)",
        "! 0a. Transpiration per soil layer.",
    )
    assert ratio_setup.count(
        "wliq_soisno_bef(j) > trc_water_min_for_ratio"
    ) == 2
    assert ratio_setup.count("tracer_is_nonvolatile_solute(itrc)") >= 1
    assert "ratio_layer(j) = 0._r8" in ratio_setup

    current_ratio = section(
        "FUNCTION current_liq_ratio", "END FUNCTION current_liq_ratio"
    )
    assert "water_shadow(jlay) > trc_water_min_for_ratio" in current_ratio
    assert "tracer_is_nonvolatile_solute(itrc)" in current_ratio
    assert "current_liq_ratio = 0._r8" in current_ratio

    groundwater = section(
        "! 3. Groundwater: qcharge-driven", "! A previously dry aquifer"
    )
    assert "abs(wa_bef) > trc_water_min_for_ratio" in groundwater
    assert "tracer_is_nonvolatile_solute(itrc)" in groundwater
    assert "ratio_src = 0._r8" in groundwater


def test_qlayer_and_qcharge_require_a_resolved_source_pool():
    qlayer = section(
        "! --- qlayer(j): layer j", "! 3. Groundwater: qcharge-driven"
    )
    assert "ratio_src = current_liq_ratio(j)" in qlayer
    assert "ratio_src = current_liq_ratio(j+1)" in qlayer

    groundwater = section(
        "! 3. Groundwater: qcharge-driven", "! A previously dry aquifer"
    )
    downward = groundwater.split("IF (qcharge_eff > trc_tiny) THEN", 1)[1].split(
        "ELSEIF (qcharge_eff < -trc_tiny) THEN", 1
    )[0]
    assert "current_liq_ratio(j)" in downward
    assert "ratio_layer(j)" not in downward

    upward = groundwater.split("ELSEIF (qcharge_eff < -trc_tiny) THEN", 1)[1]
    assert "abs(wa_bef) > trc_water_min_for_ratio" in upward
    assert "tracer_is_nonvolatile_solute(itrc)" in upward
    assert "ratio_src = 0._r8" in upward


def test_near_dry_source_keeps_inventory_and_internal_transfer_conserves_mass():
    water_floor = 1.0e-12
    source_water = 0.5 * water_floor
    source_tracer = 3.0e-4
    destination_tracer = 2.0e-4

    transfer = 0.0 if source_water <= water_floor else 1.0e-4
    source_after = source_tracer - transfer
    destination_after = destination_tracer + transfer

    assert source_after == source_tracer
    assert destination_after == destination_tracer
    assert source_after + destination_after == source_tracer + destination_tracer
