from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PHYSICS = ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Physics.F90"


def _source() -> str:
    return PHYSICS.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text.lower())


def _distribute_grid_flood(grid_fraction: float, wetland_fraction: float) -> tuple[float, float]:
    """Wetlands fill first; return wetland- and soil-relative inundation."""
    flood = min(max(grid_fraction, 0.0), 1.0)
    wetland = min(max(wetland_fraction, 0.0), 1.0)
    wet = min(flood / wetland, 1.0) if wetland > 0.0 else 0.0
    soil = max(flood - wetland, 0.0) / (1.0 - wetland) if wetland < 1.0 else 0.0
    return wet, soil


def test_grid_flood_partition_closes_area_and_handles_endpoints():
    for flood, wetland in ((0.10, 0.20), (0.70, 0.20), (0.25, 0.0), (1.0, 1.0)):
        wet, soil = _distribute_grid_flood(flood, wetland)
        assert 0.0 <= wet <= 1.0
        assert 0.0 <= soil <= 1.0
        assert abs(wetland * wet + (1.0 - wetland) * soil - flood) < 1.0e-14


def test_absolute_grid_forcing_and_hybrid_soil_use_wetland_first_partition():
    source = _flat(_source())
    assert source.count("methane_distribute_grid_finundation(") >= 4

    scheme5 = source.split("def_wetland_finundation_scheme == 5", 1)[1].split(
        "def_wetland_finundation_scheme == 6", 1
    )[0]
    scheme6 = source.split("def_wetland_finundation_scheme == 6", 1)[1].split(
        "def_wetland_finundation_scheme == 7", 1
    )[0]
    scheme7 = source.split("def_wetland_finundation_scheme == 7", 1)[1].split(
        "unsupported def_wetland_finundation_scheme", 1
    )[0]
    for branch in (scheme5, scheme7):
        assert "methane_distribute_grid_finundation(" in branch

    # Hybrid deliberately keeps wetland sigmoid seasonality.  Its soil branch
    # must nevertheless consume only the residual grid flood, avoiding overlap.
    assert "methane_distribute_grid_finundation(" in scheme6

    helper = source.split("function methane_distribute_grid_finundation", 1)[1].split(
        "end function methane_distribute_grid_finundation", 1
    )[0]
    assert "flood / wetland_fraction" in helper
    assert "max(flood - wetland_fraction, 0._r8) / & (1._r8 - wetland_fraction)" in helper


def test_subsurface_ebullition_is_reinjected_above_water_table():
    tran = _flat(_source().split("subroutine methane_tran", 1)[1].split("end subroutine methane_tran", 1)[0])
    assert re.search(
        r"if \(jwt /= 0\) then source\(jwt,1\) = source\(jwt,1\) \+ "
        r"methane_ebul_tot/dz_soisno\(jwt\)",
        tran,
    )
    assert re.search(
        r"if \(jwt /= 0\) then.*?methane_surf_ebul = 0\._r8.*?else.*?"
        r"methane_surf_ebul = methane_ebul_tot",
        tran,
    )


def test_wetland_dry_column_is_prognostic_and_area_changes_close_mass():
    source = _flat(_source())
    assert "wetland_inactive_dry_area" not in source
    assert "zwt_unsat = zi_soisno(nl_soil) + 1.e-3_r8" in source
    assert "jwt_unsat = nl_soil" in source

    # Growing inundation mixes the newly flooded dry stock into the saturated pool.
    f0, f1, sat, dry = 0.2, 0.5, 4.0, 1.0
    before = f0 * sat + (1.0 - f0) * dry
    sat_after = (f0 * sat + (f1 - f0) * dry) / f1
    after = f1 * sat_after + (1.0 - f1) * dry
    assert abs(after - before) < 1.0e-14

    # Shrinking inundation transfers the retired saturated stock into the
    # enlarged unsaturated column; area bookkeeping must not create a flux.
    f0, f1 = 0.5, 0.2
    before = f0 * sat + (1.0 - f0) * dry
    dry_after = ((1.0 - f0) * dry + (f0 - f1) * sat) / (1.0 - f1)
    after = f1 * sat + (1.0 - f1) * dry_after
    assert abs(after - before) < 1.0e-14
    assert "conc_methane_unsat(j) = ((1._r8 - fsat_bef)" in source
    assert "methane_dfsat_tot = methane_dfsat_tot -" not in source


def test_decomposition_o2_follows_methanogenesis_stoichiometry_and_root_co2():
    prod = _flat(
        _source().split("subroutine methane_prod", 1)[1].split(
            "end subroutine methane_prod", 1
        )[0]
    )
    assert "o2_decomp_depth(j) = max(0._r8, carbon_decomp_depth - & 2._r8 * methane_prod_depth(j))" in prod
    assert "if (patchtype /= 4) then o2_decomp_depth(j) = o2_decomp_depth(j) + rr_vr(j)/catomw/dz_soisno(j)" in prod
    assert "j <= jwt" not in prod.split("root respiration contributes", 1)[1].split(
        "add oxygen demand for nitrification", 1
    )[0]

    decomp, methane, root = 8.0, 3.0, 0.5
    o2_demand = decomp - 2.0 * methane + root
    co2 = decomp - methane + root
    assert o2_demand == 2.5
    assert co2 == 5.5


def test_cold_start_uses_atmospheric_henry_equilibrium_and_ice_is_not_storage():
    source = _flat(_source())
    assert source.count("if (methane_cold_start") >= 3
    assert "c_atm(1) * & (vol_gas_init + k_h_cc(j,1) * vol_liq_init)" in source
    assert "c_atm(2) * & (vol_gas_init + k_h_cc(j,2) * vol_liq_init)" in source
    assert "if (lake_liquid_depth_init <= 1.e-12_r8) lake_liquid_depth_init = lake_total_depth_init" not in source
    assert "lake_water_ch4_stock = max(lake_liquid_depth_init, 0._r8)" in source

    dz_lake = (0.5, 1.5)
    ice_fraction = (1.0, 1.0)
    liquid_depth = sum(
        dz * (1.0 - min(max(ice, 0.0), 1.0))
        for dz, ice in zip(dz_lake, ice_fraction)
    )
    assert liquid_depth == 0.0
    assert liquid_depth * 1.4 * 1.9e-6 == 0.0

    split = source.split("subroutine split_ch4_o2_phases", 1)[1].split(
        "end subroutine split_ch4_o2_phases", 1
    )[0]
    assert "vol_ch4_storage(j) = vol_aqu(j)" in split
    assert "vol_ch4_storage(j) = vol_aqu(j) + vol_ice" not in split
    assert "vol_ch4_storage(j) = vol_ice" not in split


def test_frozen_lake_retains_inventory_without_using_ice_as_dissolved_storage():
    source = _flat(_source())
    tran = _flat(
        _source().split("subroutine methane_tran", 1)[1].split(
            "end subroutine methane_tran", 1
        )[0]
    )
    setup = tran.split("lake_ch4_conc =", 1)[0]
    assert "lake_water_storage_depth = max(lake_liquid_depth, smallnumber)" in setup
    assert "lake_water_storage_depth = max(lake_total_depth, smallnumber)" not in setup
    assert "lake_frozen_ch4_stock" in source
    assert "lake_frozen_o2_stock" in source
    assert "lake_liquid_fraction_prev" in source

    def phase_transfer(mobile, frozen, previous, current):
        if current < previous:
            moved = mobile * (previous - current) / previous
            return mobile - moved, frozen + moved
        if current > previous:
            moved = frozen * (current - previous) / (1.0 - previous)
            return mobile + moved, frozen - moved
        return mobile, frozen

    # Freeze and thaw without reactions: total inventory and liquid
    # concentration stay continuous; fully frozen inventory is all immobile.
    mobile, frozen, previous, total = 2.5, 0.0, 1.0, 2.5
    for current in (0.75, 0.25, 0.0, 0.25, 0.75, 1.0):
        mobile, frozen = phase_transfer(mobile, frozen, previous, current)
        assert abs(mobile + frozen - total) < 1.0e-14
        assert abs(frozen - total * (1.0 - current)) < 1.0e-14
        if current > 0.0:
            assert abs(mobile / current - total) < 1.0e-14
        else:
            assert abs(mobile) < 1.0e-14
        previous = current

    # At fixed ice fraction the frozen pool is not consumed by liquid reactions.
    mobile, frozen = phase_transfer(total, 0.0, 1.0, 0.5)
    mobile -= 0.1
    mobile_after, frozen_after = phase_transfer(mobile, frozen, 0.5, 0.5)
    assert mobile_after == mobile
    assert frozen_after == frozen == 1.25

    # Tiny residual liquid must still transfer completely at full freeze.
    mobile, frozen = phase_transfer(2.5e-15, total - 2.5e-15, 1.0e-15, 0.0)
    assert mobile == 0.0
    assert frozen == total


def test_soil_o2_uses_positive_backward_euler_and_allows_supersaturation():
    source = _flat(_source())
    o2_postsolve = source.split("elseif (s == 2)", 1)[1].split("endif ! species", 1)[0]
    assert "logical, parameter :: backward_euler_transport = .true." in source
    assert "if (backward_euler_transport) then" in o2_postsolve
    assert "negative o2 after backward-euler transport" in o2_postsolve
    assert "conc_o2_rel(j) = max (conc_o2_rel(j), 0._r8)" in o2_postsolve
    assert "conc_o2_rel(j) = min" not in o2_postsolve
    assert "o2_rel_cap = c_atm(2)" not in o2_postsolve
    assert "o2_rel_cap = k_h_cc(j,2) * c_atm(2)" not in o2_postsolve


def test_driver_uses_the_same_soil_flood_residual():
    driver = _flat(
        (ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Driver.F90").read_text(
            encoding="utf-8"
        )
    )
    assert "f_inund_flood_patch(i) > wetland_frac_per_patch(i)" in driver
