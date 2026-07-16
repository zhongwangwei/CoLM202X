from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PHYSICS = ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Physics.F90"


def _source() -> str:
    return PHYSICS.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"&\s*", " ", text.lower()))


def _is_leap(year: int) -> bool:
    return year % 400 == 0 or (year % 4 == 0 and year % 100 != 0)


def _split_step_at_year_end(
    end_year: int, end_day: int, end_second: int, dt: float
) -> list[tuple[int, float]]:
    """Reference split for CoLM's end-of-step timestamp convention."""
    end_sec_of_year = (end_day - 1) * 86400 + end_second
    if end_sec_of_year >= dt:
        return [(end_year, dt)]
    return [(end_year - 1, dt - end_sec_of_year), (end_year, end_sec_of_year)]


def test_surface_water_depth_uses_host_or_routing_water_only() -> None:
    source = _flat(_source())
    depth_block = source.split("surface water depth must come from", 1)[1].split(
        "jwt_sat = 0", 1
    )[0]

    assert "area_forcing_min_water_depth_mm" not in depth_block
    assert "rice_ponded_depth_mm" not in depth_block
    assert "max(finundated, 0.01_r8) * routing_depth_mm" in depth_block

    inundated = 0.05
    local_depth_mm = 100.0
    patch_mean_depth_mm = max(inundated, 0.01) * local_depth_mm
    assert patch_mean_depth_mm / max(inundated, 0.01) == local_depth_mm


def test_conditional_columns_partition_host_water_without_creation() -> None:
    source = _flat(_source())
    partition = source.split("partition host soil water between conditional columns", 1)[
        1
    ].split("loop", 1)[0]
    assert "vtot < finundated * pore_volume" in partition
    assert "methane inundation exceeds host soil water" in partition
    assert "wliq_soisno_unsat(j) = max(0._r8, (vliq - finundated * vliq_sat_alloc)" in partition
    assert "physical column is dry" not in partition

    flooded = 0.25
    pore, host_liq, host_ice = 0.4, 0.18, 0.02
    host_total = host_liq + host_ice
    sat_liq = pore * host_liq / host_total
    sat_ice = pore * host_ice / host_total
    dry_liq = (host_liq - flooded * sat_liq) / (1.0 - flooded)
    dry_ice = (host_ice - flooded * sat_ice) / (1.0 - flooded)
    assert flooded * sat_liq + (1.0 - flooded) * dry_liq == host_liq
    assert flooded * sat_ice + (1.0 - flooded) * dry_ice == host_ice


def test_plant_fluxes_are_limited_only_in_joint_competition() -> None:
    source = _flat(_source())
    aere = source.split("subroutine methane_aere", 1)[1].split(
        "end subroutine methane_aere", 1
    )[0]
    tran = source.split("subroutine methane_tran", 1)[1].split(
        "end subroutine methane_tran", 1
    )[0]

    assert "aeretran_requested" not in aere
    assert "methane_aere_depth = aere" in aere
    assert "methane_tran_depth = tranloss" in aere
    assert "methane_supply = conc_methane(j) / deltim + methane_prod_depth(j)" in tran
    assert "- min(methane_aere_depth(j), 0._r8)" in tran
    assert "- min(methane_tran_depth(j), 0._r8)" in tran

    # Raw oxidation is O2-limited from 8 to 2.  A one-pass allocation then
    # leaves the full plant request of 8, instead of pre-clipping it to 2.
    available, oxidation, plant, o2_stress = 10.0, 8.0, 8.0, 0.25
    plant_after = min(plant, available - o2_stress * oxidation)
    assert plant_after == 8.0

    # Atmospheric uptake through a negative plant flux is a source available
    # to the positive sinks in the same layer and timestep.
    concentration_supply, production, aere_flux = 1.0, 0.0, -3.0
    methane_supply = concentration_supply + production - min(aere_flux, 0.0)
    assert methane_supply == 4.0


def test_annual_accumulator_uses_end_timestamp_and_skips_lakes() -> None:
    source = _flat(_source())
    call_site = source.split("do ch4 annual averages", 1)[1].split("call henry_law", 1)[0]
    annual = source.split("subroutine methane_annualupdate", 1)[1].split(
        "end subroutine methane_annualupdate", 1
    )[0]

    assert "if (patchtype /= 4) then" in call_site
    assert "end_sec_of_year" in annual
    assert "previous_year" in annual
    assert "remaining_this_year" not in annual
    assert "end_sec_of_year >= secsperyear .or. annsum_counter >= secsperyear" in annual

    assert _split_step_at_year_end(2023, 365, 86400, 1800.0) == [(2023, 1800.0)]
    assert _split_step_at_year_end(2024, 1, 600, 1800.0) == [
        (2023, 1200.0),
        (2024, 600),
    ]
    assert _split_step_at_year_end(2025, 1, 600, 1800.0) == [
        (2024, 1200.0),
        (2025, 600),
    ]
    # Exact year-end must publish even for a partial first calendar year;
    # otherwise a mid-year cold start would mix the next year's samples into it.
    assert (365 - 1) * 86400 + 86400 == 365 * 86400
    assert (366 - 1) * 86400 + 86400 == 366 * 86400
    assert _is_leap(2024) and not _is_leap(2023)
