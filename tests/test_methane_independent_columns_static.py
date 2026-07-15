from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def _text(name: str) -> str:
    return (TRACER / name).read_text()


def test_soil_and_rice_have_separate_prognostic_state():
    state = _text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = _text("MOD_Tracer_Reactive_Methane_Microbes.F90")

    assert "conc_methane_unsat_component" in state
    assert "conc_methane_sat_component" in state
    assert "fsat_bef_component" in state
    assert "finundated_lag_component" in state
    assert "layer_sat_lag_component" in state
    assert "B_methanogen_comp" in microbes
    assert "B_methanotroph_comp" in microbes
    assert "methane_microbes_step(ipatch, component" in microbes


def test_driver_solves_pure_columns_before_area_aggregation():
    driver = _text("MOD_Tracer_Reactive_Methane_Driver.F90")

    assert "CALL run_methane_component(METHANE_COMP_SOIL, .false., soil_column)" in driver
    assert "CALL run_methane_component(METHANE_COMP_RICE, .true., rice_column)" in driver
    assert "column_fraction = merge(1._r8, 0._r8, rice_column_active)" in driver
    assert "methane_surf_flux_tot(i) = ws*soil%surf_flux + wr*rice%surf_flux" in driver
    assert "methane_surf_flux_soil(i) = ws*soil%surf_flux" in driver
    assert "methane_surf_flux_rice(i) = wr*rice%surf_flux" in driver


def test_pure_rice_soil_inundation_uses_pre_override_default():
    physics = _text("MOD_Tracer_Reactive_Methane_Physics.F90")
    driver = _text("MOD_Tracer_Reactive_Methane_Driver.F90")

    assert "finundated_default_out = finundated_default" in physics
    assert "finundated_default_out=finundated_default_used" in driver
    assert "result%finundated_default = finundated_default_used" in driver
    assert "methane_soil_finundated(i) = soil%finundated_default" in driver


def test_component_restart_is_backward_compatible():
    state = _text("MOD_Tracer_Reactive_Methane_State.F90")
    facade = _text("MOD_Tracer_Reactive_Methane.F90")
    microbes = _text("MOD_Tracer_Reactive_Methane_Microbes.F90")

    assert "ch4_conc_ch4_unsat_soil" in state
    assert "ch4_conc_ch4_unsat_rice" in state
    assert "ch4_fsat_bef_soil" in state
    assert "ch4_fsat_bef_rice" in state
    assert "conc_methane_unsat_component(:,component,:) = conc_methane_unsat" in state
    assert "component_fields_present(33)" in state
    assert "component_fields_present(33) = ncio_vector_var_present" in state
    assert "'ch4_rice_fraction_prev'" in state
    assert "require_components .and. .not. all(component_fields_present)" in state
    assert "'ch4_restart_schema', 2._r8" in facade
    assert "CALL read_methane_restart (file_restart, component_state_required)" in facade
    assert (
        "CALL read_methane_microbes_restart (file_restart, component_state_required .and. file_has_pools)"
        in facade
    )
    assert "require_components .and. .not. all(component_fields_present)" in microbes
    assert "abs(schema-1._r8) > epsilon(1._r8)" in facade
    assert "abs(schema-2._r8) > epsilon(1._r8)" in facade


def test_rice_area_change_uses_conservative_donor_transfer():
    driver = _text("MOD_Tracer_Reactive_Methane_Driver.F90")
    microbes = _text("MOD_Tracer_Reactive_Methane_Microbes.F90")

    assert "CALL repartition_methane_column_state" in driver
    assert "CALL repartition_phase_state" in driver
    assert "rice_unsat(j) = mixed" in driver
    assert "rice_sat(j) = mixed" in driver
    assert "SUBROUTINE repartition_methane_microbes" in microbes


@pytest.mark.parametrize("old_fraction,new_fraction", [(0.0, 0.3), (0.7, 0.2), (0.4, 1.0), (1.0, 0.0)])
def test_rice_area_change_conserves_phase_inventory_numerically(old_fraction, new_fraction):
    soil_unsat, soil_sat, soil_h = 2.0, 8.0, 0.25
    rice_unsat, rice_sat, rice_h = 5.0, 13.0, 0.8
    old_inventory = (1.0 - old_fraction) * (
        (1.0 - soil_h) * soil_unsat + soil_h * soil_sat
    ) + old_fraction * ((1.0 - rice_h) * rice_unsat + rice_h * rice_sat)

    if new_fraction > old_fraction:
        delta = new_fraction - old_fraction
        mixed = (
            old_fraction * ((1.0 - rice_h) * rice_unsat + rice_h * rice_sat)
            + delta * ((1.0 - soil_h) * soil_unsat + soil_h * soil_sat)
        ) / new_fraction
        rice_unsat = rice_sat = mixed
    elif new_fraction < old_fraction:
        delta = old_fraction - new_fraction
        mixed = (
            (1.0 - old_fraction) * ((1.0 - soil_h) * soil_unsat + soil_h * soil_sat)
            + delta * ((1.0 - rice_h) * rice_unsat + rice_h * rice_sat)
        ) / (1.0 - new_fraction)
        soil_unsat = soil_sat = mixed

    new_inventory = (1.0 - new_fraction) * (
        (1.0 - soil_h) * soil_unsat + soil_h * soil_sat
    ) + new_fraction * ((1.0 - rice_h) * rice_unsat + rice_h * rice_sat)
    assert new_inventory == pytest.approx(old_inventory)


@pytest.mark.parametrize("rice_fraction", [0.0, 0.35, 1.0])
def test_branch_conditioning_closes_for_different_component_inundation(rice_fraction):
    soil_h, rice_h = 0.15, 0.9
    soil_unsat, soil_sat = 2.0, 11.0
    rice_unsat, rice_sat = 7.0, 19.0
    soil_weight = 1.0 - rice_fraction
    patch_h = soil_weight * soil_h + rice_fraction * rice_h

    if patch_h > 0.0:
        patch_sat = (
            soil_weight * soil_h * soil_sat + rice_fraction * rice_h * rice_sat
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
    assert (1.0 - patch_h) * patch_unsat + patch_h * patch_sat == pytest.approx(direct_total)


def test_lulcc_component_weighting_and_area_scaling_conserve_inventory():
    old_area = [40.0, 60.0]
    old_rice = [0.1, 0.8]
    soil_inventory = [2.0, 4.0]
    rice_inventory = [10.0, 20.0]
    new_area = 200.0

    old_mass = sum(
        area * ((1.0 - rice) * soil + rice * rice_col)
        for area, rice, soil, rice_col in zip(
            old_area, old_rice, soil_inventory, rice_inventory
        )
    )
    transferred_area = sum(old_area)
    rice_area = sum(area * rice for area, rice in zip(old_area, old_rice))
    soil_area = transferred_area - rice_area
    new_rice = rice_area / transferred_area
    new_soil_col = sum(
        area * (1.0 - rice) * soil
        for area, rice, soil in zip(old_area, old_rice, soil_inventory)
    ) / soil_area
    new_rice_col = sum(
        area * rice * rice_col
        for area, rice, rice_col in zip(old_area, old_rice, rice_inventory)
    ) / rice_area
    unscaled_column = (1.0 - new_rice) * new_soil_col + new_rice * new_rice_col
    target_column = old_mass / new_area
    scale = target_column / unscaled_column

    reconstructed_mass = new_area * scale * unscaled_column
    assert reconstructed_mass == pytest.approx(old_mass)

    state = _text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = _text("MOD_Tracer_Reactive_Methane_Microbes.F90")
    assert "CALL remap_component_phase" in state
    assert "scale = target/column" in state
    assert microbes.count("w = component_base_weight(ipnew, link)") == 2
    assert "w = map_mass_weight(link)" in microbes


def test_soil_to_nonsoil_to_soil_normalization_does_not_lose_rice_component():
    rice_fraction = 0.4
    soil_pool = 2.0
    rice_pool = 10.0
    aggregate_pool = (1.0 - rice_fraction) * soil_pool + rice_fraction * rice_pool

    # Entering a non-soil patch collapses the former subarea columns into the
    # live aggregate, mirrors that value to both hidden components, and sets r=0.
    normalized_soil = aggregate_pool
    normalized_rice = aggregate_pool
    normalized_fraction = 0.0
    returned_inventory = (
        (1.0 - normalized_fraction) * normalized_soil
        + normalized_fraction * normalized_rice
    )
    assert returned_inventory == pytest.approx(aggregate_pool)

    state = _text("MOD_Tracer_Reactive_Methane_State.F90")
    microbes = _text("MOD_Tracer_Reactive_Methane_Microbes.F90")
    assert "patchtype(np) /= 0" in state
    assert "CALL mirror_microbe_components_from_aggregate(np)" in microbes


def test_rice_process_diagnostics_are_direct_column_contributions():
    driver = _text("MOD_Tracer_Reactive_Methane_Driver.F90")
    accflux = _text("MOD_Tracer_Reactive_Methane_AccFlux.F90")
    history = _text("MOD_Tracer_Reactive_Methane_Hist.F90")

    assert "methane_surf_aere_rice(i) = wr*rice%surf_aere" in driver
    assert "methane_surf_ebul_rice(i) = wr*rice%surf_ebul" in driver
    assert "methane_surf_diff_rice(i) = wr*rice%surf_diff" in driver
    assert "methane_prod_tot_rice(i) = wr*rice%prod_tot" in driver
    assert "methane_oxid_tot_rice(i) = wr*rice%oxid_tot" in driver
    assert "CALL acc1d (methane_surf_aere_rice, a_methane_surf_aere_rice)" in accflux
    assert "CALL acc1d (methane_surf_ebul_rice, a_methane_surf_ebul_rice)" in accflux
    assert "'f_methane_surf_aere_rice'" in history
    assert "'f_methane_surf_ebul_rice'" in history
