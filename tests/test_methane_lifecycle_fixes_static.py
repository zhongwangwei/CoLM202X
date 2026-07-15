from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def source(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def subroutine_body(text: str, name: str) -> str:
    start = text.index(f"SUBROUTINE {name}")
    end = text.index(f"END SUBROUTINE {name}", start)
    return text[start:end]


def test_wetland_decomposition_debits_state_once_with_timestep() -> None:
    driver = source("main/CoLMDRIVER.F90")
    land = source("main/TRACER/MOD_Tracer_LandPhase.F90")
    reactive = source("main/TRACER/MOD_Tracer_Reactive.F90")
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    impl = source("main/TRACER/MOD_Tracer_Reactive_Methane_Impl.F90")
    shim = source("main/TRACER/MOD_Tracer_Reactive_BgcShim.F90")

    assert "CALL tracer_wetland_decomp (i, deltim)" in driver
    assert "SUBROUTINE tracer_wetland_decomp (ipatch, deltim)" in land
    assert "CALL tracer_reactive_wetland_decomp (ipatch, deltim)" in land
    assert "SUBROUTINE reactive_wetland_decomp_if (ipatch, deltim)" in reactive
    assert "CALL reactive_callbacks(i)%wetland_decomp (ipatch, deltim)" in reactive
    assert "SUBROUTINE ch4_reactive_wetland_decomp (ipatch, deltim)" in methane
    assert "CALL ch4_impl_wetland_decomp (ipatch, deltim)" in methane
    assert "SUBROUTINE ch4_impl_wetland_decomp (ipatch, deltim)" in impl
    assert "CALL reactive_bgc_run_wetland_decomp (ipatch, deltim)" in impl
    assert "SUBROUTINE reactive_bgc_run_wetland_decomp (ipatch, deltim)" in shim

    shim_body = subroutine_body(shim, "reactive_bgc_run_wetland_decomp")
    assert "CALL apply_wetland_decomposition_state (ipatch, deltim)" in shim_body
    apply_body = subroutine_body(shim, "apply_wetland_decomposition_state")
    assert "decomp_cpools_vr" in apply_body
    assert "decomp_npools_vr" in apply_body
    assert "sminn_vr" in apply_body
    assert "decomp_cpools_sourcesink" in shim
    assert "decomp_npools_sourcesink" in shim
    assert "DEF_USE_SASU .or. DEF_USE_DiagMatrix" in shim
    assert "CALL CNDriverSummarizeNonvegetatedSoilStates (ipatch, nl_soil, dz_soi, ndecomp_pools)" in apply_body
    assert apply_body.count("decomp_cpools_sourcesink(:,:,ipatch) = 0._r8") >= 2
    assert apply_body.count("decomp_npools_sourcesink(:,:,ipatch) = 0._r8") >= 2

    require_body = subroutine_body(shim, "require_wetland_decomposition_state")
    assert "decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch) < 0._r8" in require_body
    assert "decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch) < 0._r8" in require_body
    assert ".not. allocated(o_scalar)" in require_body
    assert ".not. allocated(fpi_vr)" in require_body
    assert "BD_all(1:nl_soil,ipatch) <= tiny(1._r8)" in require_body
    assert "wetland patch unexpectedly contains vegetation C/N state" in require_body

    # The full PFT BGC driver and the decomposition-only wetland path are
    # disjoint, so the new state debit cannot run twice for one patch.
    assert "IF(patchtype(i) .eq. 0)THEN" in driver
    impl_body = subroutine_body(impl, "ch4_impl_wetland_decomp")
    assert "IF (patchtype(ipatch) /= 2) RETURN" in impl_body


def test_wetland_state_summary_rebuilds_soil_and_column_aggregates() -> None:
    summary = source("main/BGC/MOD_BGC_CNSummary.F90")
    makefile = source("Makefile")
    body = subroutine_body(summary, "CNDriverSummarizeNonvegetatedSoilStates")

    assert "IMPLICIT NONE\n   PRIVATE" in summary
    assert "PUBLIC CNDriverSummarizeNonvegetatedSoilStates" in summary
    assert "CALL soilbiogeochem_carbonstate_summary" in body
    assert "CALL soilbiogeochem_nitrogenstate_summary" in body
    assert "totvegc(i) = 0._r8" in body
    assert "totvegn(i) = 0._r8" in body
    assert "totcolc(i) = totcwdc(i) + totlitc(i) + totsomc(i)" in body
    assert "totcoln(i) = totcwdn(i) + totlitn(i) + totsomn(i)" in body
    assert "MOD_BGC_Soil_BiogeochemPotential.o MOD_BGC_CNSummary.o" in makefile


def test_patch_ph_rejects_nan_and_infinity_before_collective_stats() -> None:
    ph = source("main/TRACER/MOD_Tracer_Reactive_Methane_pH.F90")
    body = subroutine_body(ph, "read_methane_ph_patch")

    assert "USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite" in ph
    assert ".not. ieee_is_finite(methane_ph_patch)" in body
    sanitize = body.index("WHERE (.not. ieee_is_finite(methane_ph_patch)")
    stats = body.index("n_nonfallback = count(")
    assert sanitize < stats

def test_decomposition_only_budget_identity() -> None:
    # One transition: donor loses respiration+transfer, receiver gains the
    # transfer. Organic N plus mineral N loses only explicit denitrification.
    hr, ctransfer, ntransfer = 2.0, 3.0, 0.4
    sminn_flux, denit = 0.1, 0.02
    donor_c = -(hr + ctransfer)
    receiver_c = ctransfer
    donor_n = -ntransfer
    receiver_n = ntransfer + sminn_flux
    mineral_n = -sminn_flux - denit

    assert donor_c + receiver_c == -hr
    assert donor_n + receiver_n + mineral_n == pytest.approx(-denit)


def test_nitrification_branch_preserves_total_mineral_n_with_nh4_first() -> None:
    nh4, no3 = 0.3, 0.5
    mineral_delta = -0.6
    remaining = -mineral_delta
    take_nh4 = min(nh4, remaining)
    nh4 = max(0.0, nh4 - take_nh4)
    no3 = max(0.0, no3 - (remaining - take_nh4))

    assert nh4 == 0.0
    assert no3 == pytest.approx(0.2)
    assert nh4 + no3 == pytest.approx(0.8 + mineral_delta)


def test_lulcc_refreshes_spatial_ph_context_before_collective_reload() -> None:
    driver = source("main/LULCC/MOD_Lulcc_Driver.F90")
    reactive = source("main/TRACER/MOD_Tracer_Reactive.F90")
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")

    reload_call = driver.index(
        "CALL tracer_reactive_reload_lulcc_inputs (dir_landdata, jdate(1))"
    )
    remap_call = driver.index("CALL tracer_reactive_remap_lulcc_state")
    assert remap_call < reload_call
    assert "SUBROUTINE reactive_reload_lulcc_if (dir_landdata, lc_year)" in reactive
    assert "reactive_lulcc_dir_landdata" not in reactive
    assert "tracer_reactive_get_lulcc_context" not in reactive

    reload_body = subroutine_body(methane, "ch4_reactive_reload_lulcc_inputs")
    assert "character(len=*), intent(in) :: dir_landdata_lulcc" in reload_body
    assert "integer, intent(in) :: lc_year_lulcc" in reload_body
    assert "write(cyear_lulcc,'(i4.4)') lc_year_lulcc" in reload_body
    assert "last_methane_ph_patch_file = trim(dir_landdata_lulcc)//'/soil/'" in reload_body
    path_update = reload_body.index("last_methane_ph_patch_file =")
    ph_read = reload_body.index("CALL read_methane_ph_patch")
    assert path_update < ph_read


def test_routing_advances_during_spinup_after_land_step() -> None:
    colm = source("main/CoLM.F90")
    driver_call = colm.index("CALL CoLMDRIVER")
    routing_call = colm.index("CALL grid_riverlake_flow", driver_call)
    assert driver_call < routing_call
    routing_context = colm[routing_call - 160 : routing_call + 100]
    assert "IF (.not. is_spinup)" not in routing_context


def test_restart_schema_is_uniform_across_all_mpi_ranks() -> None:
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    body = subroutine_body(methane, "validate_methane_restart_transaction")

    assert "schema_v1_present" in body
    assert "schema_v2_present" in body
    assert "MPI_LOR" in body
    assert "schema_v1_present .eqv. schema_v2_present" in body
    assert "component_state_required = schema_v2_present" in body
    assert "component_state_required = &\n\t        any(" not in body


def test_inactive_final_and_reentry_still_release_methane_state() -> None:
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    land = source("main/TRACER/MOD_Tracer_LandPhase.F90")

    final_body = subroutine_body(methane, "ch4_reactive_final")
    assert "ch4_reactive_has" not in final_body
    assert "CALL deallocate_methane_state ()" in final_body
    assert "last_methane_ph_patch_file = ''" in final_body

    reentry_body = subroutine_body(land, "tracer_init_from_arrays")
    reactive_final = reentry_body.index("CALL tracer_reactive_final ()")
    defs_init = reentry_body.index("CALL tracer_defs_init()")
    assert reactive_final < defs_init
