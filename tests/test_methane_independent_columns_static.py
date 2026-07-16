from pathlib import Path
import re

import pytest


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def text(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_soil_and_rice_keep_independent_prognostic_and_microbe_state() -> None:
    state = text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = text("MOD_Tracer_Reactive_Methane_Microbes.F90")

    for field in (
        "conc_methane_unsat_component",
        "conc_methane_sat_component",
        "conc_o2_unsat_component",
        "conc_o2_sat_component",
        "layer_sat_lag_component",
        "fsat_bef_component",
        "finundated_lag_component",
    ):
        assert field in state
    for field in (
        "B_methanogen_comp",
        "B_methanotroph_comp",
        "B_methanogen_dormant_comp",
        "B_methanotroph_dormant_comp",
    ):
        assert field in microbes
    step = routine(microbes, "methane_microbes_step")
    assert "B_methanogen(j,ipatch)" not in step
    assert "B_methanogen_comp(j,component,ipatch)" in step


def test_microbe_component_aggregation_handles_endpoints_and_exact_repartition() -> None:
    microbes = text("MOD_Tracer_Reactive_Methane_Microbes.F90")
    aggregate = routine(microbes, "aggregate_field")
    repartition = routine(microbes, "repartition_methane_microbes")

    assert "IF (rice_weight <= 0._r8) THEN" in aggregate
    assert "aggregate_values = component_values(:,METHANE_COMP_SOIL,ipatch)" in aggregate
    assert "ELSEIF (rice_weight >= 1._r8) THEN" in aggregate
    assert "aggregate_values = component_values(:,METHANE_COMP_RICE,ipatch)" in aggregate
    assert "IF (.not. (new_rice > old_rice .or. new_rice < old_rice)) RETURN" in repartition
    assert "epsilon(1._r8)" not in repartition


def test_driver_solves_pure_columns_then_weights_once_and_finalizes_once() -> None:
    driver = text("MOD_Tracer_Reactive_Methane_Driver.F90")
    physics = text("MOD_Tracer_Reactive_Methane_Physics.F90")

    soil = "CALL run_methane_component(METHANE_COMP_SOIL, .false., soil_column)"
    rice = "CALL run_methane_component(METHANE_COMP_RICE, .true., rice_column)"
    aggregate = "CALL aggregate_methane_columns(soil_column, rice_column, rice_weight)"
    assert soil in driver and rice in driver and aggregate in driver
    assert "column_fraction = merge(1._r8, 0._r8, rice_column_active)" in driver
    assert "methane_surf_flux_tot(i) = ws*soil%surf_flux + wr*rice%surf_flux" in driver
    assert "methane_surf_flux_soil(i) = ws*soil%surf_flux" in driver
    assert "methane_surf_flux_rice(i) = wr*rice%surf_flux" in driver
    assert driver.index(aggregate) < driver.index(
        "CALL tracer_ch4_bgc_finalize_step", driver.index(aggregate)
    )
    assert "tracer_ch4_bgc_finalize_step" not in physics


def test_pure_rice_soil_inundation_reports_the_pre_override_default() -> None:
    physics = text("MOD_Tracer_Reactive_Methane_Physics.F90")
    driver = text("MOD_Tracer_Reactive_Methane_Driver.F90")

    assert "finundated_default_out = finundated_default" in physics
    assert "finundated_default_out=finundated_default_used" in driver
    assert "result%finundated_default = finundated_default_used" in driver
    # In a pure-rice patch soil_column is copied from the solved rice column;
    # the retained field is nevertheless the default before the paddy override.
    assert "soil_column = rice_column" in driver
    assert "methane_soil_finundated(i) = soil%finundated_default" in driver


def test_component_pft_forcing_is_intensive_and_uses_one_local_scan() -> None:
    bgc = text("MOD_Tracer_Reactive_Methane_BgcLink.F90")
    driver = text("MOD_Tracer_Reactive_Methane_Driver.F90")
    helper = routine(bgc, "tracer_ch4_bgc_component_veg_inputs")

    assert helper.count("DO m = ps, pe") == 1
    assert "component_fraction(METHANE_COMP_SOIL) = 1._r8 - rice_area" in helper
    for field in (
        "lai_component(component)",
        "agnpp_component(component)",
        "bgnpp_component(component)",
        "rr_component(component)",
        "annsum_npp_component(component)",
    ):
        assert f"{field} = {field} / component_fraction(component)" in helper
    assert "crootfr_component(:,component) = crootfr_component(:,component) / profile_sum" in helper
    assert "CALL get_rice_veg_proxy(column_lai, i, 1._r8)" in driver
    assert "column_agnpp,column_bgnpp" in driver
    assert "column_crootfr" in driver
    assert "mpi_" not in helper.lower()


def test_component_intensive_weighting_reconstructs_patch_forcing_with_bare_ground() -> None:
    # Rice covers 1/4 of the patch, non-rice vegetation 1/2, and bare ground
    # the remaining 1/4.  Bare ground belongs to the soil column denominator.
    rice_area = 0.25
    nonrice_area = 0.50
    soil_area = 1.0 - rice_area
    rice_pft_value = 11.0
    nonrice_pft_value = 4.0

    rice_intensive = rice_area * rice_pft_value / rice_area
    soil_intensive = nonrice_area * nonrice_pft_value / soil_area
    reconstructed = rice_area * rice_intensive + soil_area * soil_intensive
    direct_patch_sum = rice_area * rice_pft_value + nonrice_area * nonrice_pft_value

    assert rice_intensive == pytest.approx(rice_pft_value)
    assert soil_intensive == pytest.approx(8.0 / 3.0)
    assert reconstructed == pytest.approx(direct_patch_sum)


@pytest.mark.parametrize("rice_fraction", [0.0, 0.35, 1.0])
def test_conditional_sat_unsat_aggregation_reconstructs_total_inventory(
    rice_fraction: float,
) -> None:
    soil_h, rice_h = 0.15, 0.90
    soil_unsat, soil_sat = 2.0, 11.0
    rice_unsat, rice_sat = 7.0, 19.0
    soil_weight = 1.0 - rice_fraction
    patch_h = soil_weight * soil_h + rice_fraction * rice_h

    if patch_h > 0.0:
        patch_sat = (
            soil_weight * soil_h * soil_sat
            + rice_fraction * rice_h * rice_sat
        ) / patch_h
    else:
        patch_sat = soil_weight * soil_sat + rice_fraction * rice_sat
    if patch_h < 1.0:
        patch_unsat = (
            soil_weight * (1.0 - soil_h) * soil_unsat
            + rice_fraction * (1.0 - rice_h) * rice_unsat
        ) / (1.0 - patch_h)
    else:
        patch_unsat = soil_weight * soil_unsat + rice_fraction * rice_unsat

    direct_total = soil_weight * (
        (1.0 - soil_h) * soil_unsat + soil_h * soil_sat
    ) + rice_fraction * ((1.0 - rice_h) * rice_unsat + rice_h * rice_sat)
    assert (1.0 - patch_h) * patch_unsat + patch_h * patch_sat == pytest.approx(
        direct_total
    )


def test_floodplain_gate_requires_residual_after_static_wetland() -> None:
    driver = text("MOD_Tracer_Reactive_Methane_Driver.F90")
    component = routine(driver, "run_methane_component")

    assert "allocated(wetland_frac_per_patch)" in component
    assert "f_inund_flood_patch(i) > DEF_METHANE%hybrid_soil_threshold" in component
    assert "f_inund_flood_patch(i) > wetland_frac_per_patch(i)" in component
    assert "get_biome_redoxlag" in component
    assert "is_floodplain_active" in component.split("get_biome_redoxlag", 1)[1]


def test_schema3_component_groups_are_atomic_and_batch_probed() -> None:
    facade = text("MOD_Tracer_Reactive_Methane.F90")
    state = text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = text("MOD_Tracer_Reactive_Methane_Microbes.F90")

    assert "METHANE_RESTART_SCHEMA_VERSION = 4" in facade
    assert "nmicrobe_restart_fields = 25" in facade
    assert "n_microbe_component_fields = count(microbe_field_present(5:12))" in facade
    assert "restart_schema >= 3 .and. (file_has_pools .neqv. file_has_component_pools)" in facade
    assert "ncio_vector_group_presence(file_restart, component_restart_fields" in state
    microbe_reader = routine(microbes, "read_methane_microbes_restart")
    assert "ncio_vector_group_presence" in microbe_reader
    assert "ncio_vector_var_present" not in microbe_reader
    assert microbe_reader.index("invalid component microbial biomass values") < microbe_reader.index(
        "WHERE (ieee_is_nan(B_methanogen_comp)"
    )


def test_rice_history_uses_same_window_raw_sums() -> None:
    history = text("MOD_Tracer_Reactive_Methane_Hist.F90")
    accflux = text("MOD_Tracer_Reactive_Methane_AccFlux.F90")

    assert "hist_ch4_rice_flux_mean(ipatch) = a_methane_surf_flux_rice(ipatch)" in history
    assert "hist_ch4_rice_area_frac(ipatch) = a_methane_rice_fraction(ipatch)" in history
    assert "paddy_rice_fraction(ipatch)" not in history
    assert "'ch4_a_methane_rice_fraction'" in routine(
        accflux, "write_methane_accflux_restart"
    )
    assert "'ch4_a_methane_rice_fraction'" in routine(
        accflux, "read_methane_accflux_restart"
    )


@pytest.mark.parametrize(
    "old_fraction,new_fraction", [(0.0, 0.3), (0.7, 0.2), (0.4, 1.0), (1.0, 0.0)]
)
def test_rice_area_transfer_conserves_phase_inventory(
    old_fraction: float, new_fraction: float
) -> None:
    soil_u, soil_s, soil_h = 2.0, 8.0, 0.25
    rice_u, rice_s, rice_h = 5.0, 13.0, 0.8
    old = (1.0 - old_fraction) * ((1.0 - soil_h) * soil_u + soil_h * soil_s)
    old += old_fraction * ((1.0 - rice_h) * rice_u + rice_h * rice_s)

    if new_fraction > old_fraction:
        delta = new_fraction - old_fraction
        mixed = (
            old_fraction * ((1.0 - rice_h) * rice_u + rice_h * rice_s)
            + delta * ((1.0 - soil_h) * soil_u + soil_h * soil_s)
        ) / new_fraction
        rice_u = rice_s = mixed
    elif new_fraction < old_fraction:
        delta = old_fraction - new_fraction
        mixed = (
            (1.0 - old_fraction) * ((1.0 - soil_h) * soil_u + soil_h * soil_s)
            + delta * ((1.0 - rice_h) * rice_u + rice_h * rice_s)
        ) / (1.0 - new_fraction)
        soil_u = soil_s = mixed

    new = (1.0 - new_fraction) * ((1.0 - soil_h) * soil_u + soil_h * soil_s)
    new += new_fraction * ((1.0 - rice_h) * rice_u + rice_h * rice_s)
    assert new == pytest.approx(old)


def test_lulcc_component_and_target_area_weights_conserve_inventory() -> None:
    old_area = [40.0, 60.0]
    old_rice = [0.1, 0.8]
    soil_inventory = [2.0, 4.0]
    rice_inventory = [10.0, 20.0]
    new_area = 200.0

    old_mass = sum(
        area * ((1.0 - rice) * soil + rice * rice_column)
        for area, rice, soil, rice_column in zip(
            old_area, old_rice, soil_inventory, rice_inventory
        )
    )
    transferred_area = sum(old_area)
    rice_area = sum(area * rice for area, rice in zip(old_area, old_rice))
    soil_area = transferred_area - rice_area
    new_rice = rice_area / transferred_area
    new_soil_column = sum(
        area * (1.0 - rice) * soil
        for area, rice, soil in zip(old_area, old_rice, soil_inventory)
    ) / soil_area
    new_rice_column = sum(
        area * rice * rice_column
        for area, rice, rice_column in zip(old_area, old_rice, rice_inventory)
    ) / rice_area
    unscaled_column = (
        (1.0 - new_rice) * new_soil_column + new_rice * new_rice_column
    )
    target_column = old_mass / new_area
    scale = target_column / unscaled_column
    assert new_area * scale * unscaled_column == pytest.approx(old_mass)

    state = text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = text("MOD_Tracer_Reactive_Methane_Microbes.F90")
    assert "CALL remap_component_phase" in state
    assert "w = component_transfer_weight(np, link, component)" in state
    assert "scale = target/column" in state
    for source in (state, microbes):
        base_weight = source.split("FUNCTION component_base_weight", 1)[1].split(
            "END FUNCTION component_base_weight", 1
        )[0]
        assert "w = map_mass_weight(link)" in base_weight
        assert "w = map_source_weight(link)" in base_weight


def test_nonsoil_lulcc_targets_mirror_hidden_component_state() -> None:
    rice_fraction = 0.4
    soil_pool = 2.0
    rice_pool = 10.0
    aggregate_pool = (1.0 - rice_fraction) * soil_pool + rice_fraction * rice_pool

    normalized_soil = aggregate_pool
    normalized_rice = aggregate_pool
    normalized_fraction = 0.0
    returned_inventory = (
        (1.0 - normalized_fraction) * normalized_soil
        + normalized_fraction * normalized_rice
    )
    assert returned_inventory == pytest.approx(aggregate_pool)

    state = text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = text("MOD_Tracer_Reactive_Methane_Microbes.F90")
    assert "patchtype(np) /= 0" in state
    assert "spread(conc_methane_unsat(:,np), 2, N_METHANE_COMP)" in state
    assert "CALL mirror_microbe_components_from_aggregate(np)" in microbes


def test_rice_process_diagnostics_flow_directly_to_accumulators_and_history() -> None:
    driver = text("MOD_Tracer_Reactive_Methane_Driver.F90")
    accflux = text("MOD_Tracer_Reactive_Methane_AccFlux.F90")
    history = text("MOD_Tracer_Reactive_Methane_Hist.F90")

    paths = {
        "surf_aere": "surf_aere",
        "surf_ebul": "surf_ebul",
        "surf_diff": "surf_diff",
        "prod_tot": "prod_tot",
        "oxid_tot": "oxid_tot",
    }
    for field, column_member in paths.items():
        assert f"methane_{field}_rice(i) = wr*rice%{column_member}" in driver
        assert re.search(
            rf"CALL\s+acc1d\s*\(methane_{field}_rice\s*,\s*"
            rf"a_methane_{field}_rice\s*\)",
            accflux,
        )
        assert f"'f_methane_{field}_rice'" in history
