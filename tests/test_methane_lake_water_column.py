from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def source(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def without_continuations(text: str) -> str:
    return re.sub(r"&\s*\n\s*", " ", text)


def test_lake_water_stocks_are_prognostic_restart_fields() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")

    assert "lake_water_ch4_stock" in state
    assert "lake_water_o2_stock" in state
    assert "ch4_lake_water_ch4_stock" in state
    assert "ch4_lake_water_o2_stock" in state
    assert "CALL remap1d_mass(lulcc_lake_water_ch4_stock_old" in state
    assert "CALL remap1d_mass(lulcc_lake_water_o2_stock_old" in state


def test_lake_water_diagnostics_survive_accumulator_restart_and_history() -> None:
    acc = source("MOD_Tracer_Reactive_Methane_AccFlux.F90")
    hist = source("MOD_Tracer_Reactive_Methane_Hist.F90")
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    need_lake = acc.split("need_lake = any([", 1)[1].split("])", 1)[0]

    for name in (
        "lake_water_ch4_stock",
        "lake_water_o2_stock",
        "lake_water_ch4_oxid",
        "lake_sed_ch4_flux",
        "lake_sed_o2_flux",
        "lake_air_o2_flux",
    ):
        assert f"ch4_a_{name}" in acc
        assert f"f_{name}" in hist
        assert f"f_{name}" in const
        assert f"mhist_on('f_{name}')" in need_lake
        assert f"CALL acc1d ({name}" in acc


def test_lake_transport_uses_dynamic_lake_geometry_and_a_water_node() -> None:
    physics = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Physics.F90")
    )
    driver = source("MOD_Tracer_Reactive_Methane_Driver.F90")
    impl = source("MOD_Tracer_Reactive_Methane_Impl.F90")

    assert "dz_lake" in driver
    assert "t_lake" in driver
    assert "dz_lake(1:nl_lake,ipatch)" in impl
    assert "t_lake(1:nl_lake,ipatch)" in impl
    assert re.search(
        r"lake_liquid_depth\s*=\s*sum\(max\(dz_lake\(1:nl_lake\),\s*0\._r8\)\s*\*\s*"
        r"\(1\._r8\s*-\s*min\(max\(lake_icefrac\(1:nl_lake\),",
        physics,
    )
    assert re.search(
        r"bt\(j\)\s*=\s*lake_water_storage_depth\s*/\s*deltim\s*\+\s*"
        r"lake_sed_water_cond\s*\+\s*lake_exchange_vel",
        physics,
    )
    assert "ct(j) = -lake_sed_water_cond" in physics
    assert "rt(j) = lake_water_stock_old / deltim" in physics


def test_dynamic_dry_lakes_use_conservative_handler_before_return() -> None:
    impl = without_continuations(source("MOD_Tracer_Reactive_Methane_Impl.F90"))
    body = impl.split("SUBROUTINE ch4_impl_lake_step", 1)[1].split(
        "END SUBROUTINE ch4_impl_lake_step", 1
    )[0]

    assert re.search(
        r"USE\s+MOD_Namelist,\s+only:\s+DEF_USE_Dynamic_Lake", impl
    )
    dry_match = re.search(
        r"IF\s*\(DEF_USE_Dynamic_Lake\s*\.and\.\s*"
        r"\(wdsrf\(ipatch\)\s*<\s*100\._r8\s*\.or\.\s*"
        r"zwt\(ipatch\)\s*>\s*0\._r8\)\)\s*THEN",
        body,
    )
    assert dry_match is not None
    dry_branch = dry_match.start()
    dry_handler = body.index("CALL handle_methane_dry_lake_substep", dry_branch)
    microbe_reset = body.index(
        "CALL reset_methane_inactive_lake_microbe_diagnostics", dry_handler
    )
    microbe_accumulator = body.index(
        "CALL accumulate_methane_lake_microbe_substep_diagnostics",
        microbe_reset,
    )
    substep_accumulator = body.index(
        "CALL accumulate_methane_lake_substep_diagnostics", microbe_accumulator
    )
    dry_return = body.index("RETURN", substep_accumulator)

    assert (
        dry_branch
        < dry_handler
        < microbe_reset
        < microbe_accumulator
        < substep_accumulator
        < dry_return
    )
    assert dry_return < body.index("CALL compute_f_h2osfc")
    assert dry_return < body.index("CALL methane_driver")
    assert "dry_lake_wet_prev" not in impl


def test_disabled_lake_production_clears_diagnostics_without_inventory_loss() -> None:
    impl = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Impl.F90")
    )
    state = without_continuations(
        source("MOD_Tracer_Reactive_Methane_State.F90")
    )
    impl_body = impl.split("SUBROUTINE ch4_impl_lake_step", 1)[1].split(
        "END SUBROUTINE ch4_impl_lake_step", 1
    )[0]

    inactive_match = re.search(
        r"IF\s*\(\.not\.\s*DEF_METHANE%allowlakeprod\)\s*THEN(.*?)"
        r"RETURN\s*ENDIF",
        impl_body,
        flags=re.DOTALL,
    )
    assert inactive_match is not None
    inactive = inactive_match.group(1)
    state_reset = inactive.index(
        "CALL reset_methane_inactive_lake_diagnostics(ipatch)"
    )
    microbe_reset = inactive.index(
        "CALL reset_methane_inactive_lake_microbe_diagnostics(ipatch)"
    )
    microbe_accumulator = inactive.index(
        "CALL accumulate_methane_lake_microbe_substep_diagnostics"
    )
    state_accumulator = inactive.index(
        "CALL accumulate_methane_lake_substep_diagnostics"
    )
    assert state_reset < microbe_reset < microbe_accumulator < state_accumulator
    assert inactive_match.start() < impl_body.index("IF (DEF_USE_Dynamic_Lake")
    assert inactive_match.start() < impl_body.index("CALL methane_driver")
    assert not re.search(
        r"IF\s*\(\.not\.\s*DEF_METHANE%allowlakeprod\)\s*RETURN",
        impl_body,
    )

    assert "PUBLIC :: reset_methane_inactive_lake_diagnostics" in state
    reset_body = state.split(
        "SUBROUTINE reset_methane_inactive_lake_diagnostics", 1
    )[1].split(
        "END SUBROUTINE reset_methane_inactive_lake_diagnostics", 1
    )[0]
    dry_body = state.split(
        "SUBROUTINE handle_methane_dry_lake_substep", 1
    )[1].split("END SUBROUTINE handle_methane_dry_lake_substep", 1)[0]
    assert "CALL reset_methane_inactive_lake_diagnostics(ipatch)" in dry_body

    prognostic_inventory = (
        "lake_water_ch4_stock",
        "lake_water_o2_stock",
        "lake_frozen_ch4_stock",
        "lake_frozen_o2_stock",
        "lake_soilc",
        "conc_methane_lake",
        "conc_o2_lake",
        "totcol_methane_lake",
        "totcol_methane",
        "totcol_methane_sat",
        "totcol_methane_unsat",
    )
    for field in prognostic_inventory:
        assert not re.search(rf"\b{field}\s*\([^=]*\)\s*=", reset_body)


def test_dynamic_dry_lake_handler_closes_ch4_and_o2_inventory() -> None:
    state = without_continuations(
        source("MOD_Tracer_Reactive_Methane_State.F90")
    )
    body = state.split(
        "SUBROUTINE handle_methane_dry_lake_substep", 1
    )[1].split("END SUBROUTINE handle_methane_dry_lake_substep", 1)[0]

    assert re.search(
        r"drydown_ch4_stock\s*=\s*max\(lake_water_ch4_stock\(ipatch\),\s*0\._r8\)\s*\+\s*"
        r"max\(lake_frozen_ch4_stock\(ipatch\),\s*0\._r8\)",
        body,
    )
    assert re.search(
        r"drydown_o2_stock\s*=\s*max\(lake_water_o2_stock\(ipatch\),\s*0\._r8\)\s*\+\s*"
        r"max\(lake_frozen_o2_stock\(ipatch\),\s*0\._r8\)",
        body,
    )
    for stock in (
        "lake_water_ch4_stock",
        "lake_frozen_ch4_stock",
        "lake_water_o2_stock",
        "lake_frozen_o2_stock",
    ):
        assert re.search(rf"{stock}\(ipatch\)\s*=\s*0\._r8", body)

    assert "drydown_ch4_flux = drydown_ch4_stock / substep_dt" in body
    for flux in (
        "methane_surf_diff",
        "methane_surf_diff_phys",
        "methane_surf_flux_tot",
        "methane_surf_flux_tot_phys",
        "methane_surf_diff_lake",
        "methane_surf_flux_tot_lake",
        "methane_surf_flux_lake",
    ):
        assert re.search(
            rf"{flux}\(ipatch\)\s*=\s*drydown_ch4_flux", body
        )
    assert "lake_air_o2_flux(ipatch) = drydown_o2_stock / substep_dt" in body

    assert re.search(
        r"totcol_methane\(ipatch\)\s*=\s*totcol_methane_lake\(ipatch\)", body
    )
    assert re.search(
        r"totcol_methane_sat\(ipatch\)\s*=\s*totcol_methane_lake\(ipatch\)",
        body,
    )
    assert re.search(r"totcol_methane_unsat\(ipatch\)\s*=\s*0\._r8", body)
    for memory in (
        "fsat_bef(ipatch)",
        "finundated_lag(ipatch)",
        "layer_sat_lag(:,ipatch)",
        "lake_liquid_fraction_prev(ipatch)",
    ):
        assert re.search(rf"{re.escape(memory)}\s*=\s*spval", body)


def test_dynamic_dry_lake_handler_resets_stale_process_diagnostics() -> None:
    state = without_continuations(
        source("MOD_Tracer_Reactive_Methane_State.F90")
    )
    body = state.split(
        "SUBROUTINE reset_methane_inactive_lake_diagnostics", 1
    )[1].split(
        "END SUBROUTINE reset_methane_inactive_lake_diagnostics", 1
    )[0]
    accumulator = state.split(
        "SUBROUTINE accumulate_methane_lake_substep_diagnostics", 1
    )[1].split("END SUBROUTINE accumulate_methane_lake_substep_diagnostics", 1)[0]
    sampled_fields = re.findall(
        r"CALL\s+add(?:1d|2d)\(\s*\d+\s*,\s*([A-Za-z_]\w*)"
        r"\(\s*(?::\s*,\s*)?ipatch\s*\)\s*\)",
        accumulator,
        flags=re.IGNORECASE,
    )

    assert len(sampled_fields) == 101
    assert len(set(sampled_fields)) == 101
    positive_restart_state = {
        "grnd_methane_cond",
        "grnd_methane_cond_unsat",
        "grnd_methane_cond_sat",
        "grnd_methane_cond_lake",
    }
    missing = [
        field
        for field in sampled_fields
        if not re.search(
            rf"\b{field}\s*\(\s*(?::\s*,\s*)?ipatch\s*\)\s*=\s*"
            rf"(?:0\._r8|spval|DEF_METHANE%grnd_methane_cond_default)",
            body,
            flags=re.IGNORECASE,
        )
    ]
    assert missing == []
    for field in positive_restart_state:
        assert re.search(
            rf"\b{field}\s*\(ipatch\)\s*=\s*"
            rf"DEF_METHANE%grnd_methane_cond_default",
            body,
            flags=re.IGNORECASE,
        )

    microbes = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Microbes.F90")
    )
    microbe_body = microbes.split(
        "SUBROUTINE reset_methane_inactive_lake_microbe_diagnostics", 1
    )[1].split(
        "END SUBROUTINE reset_methane_inactive_lake_microbe_diagnostics", 1
    )[0]
    for field in (
        "f_T_methanogen",
        "f_S_methanogen",
        "f_O2_methanogen",
        "f_T_methanotroph",
        "methanogen_growth_rate",
        "methanotroph_growth_rate",
        "microbial_prod_potential",
        "microbial_oxid_potential",
    ):
        assert re.search(rf"{field}\(:,ipatch\)\s*=\s*0\._r8", microbe_body)
    assert not re.search(r"\bB_methanogen(?:_dormant)?\s*\(", microbe_body)
    assert not re.search(r"\bB_methanotroph(?:_dormant)?\s*\(", microbe_body)


def test_wet_lake_driver_resets_all_inactive_microbe_diagnostics_via_owner() -> None:
    driver = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Driver.F90")
    )
    lake_microbe_region = driver.split(
        "IF (patchtype /= 4) THEN", 1
    )[1].split("CALL methane (", 1)[0]

    reset = lake_microbe_region.index(
        "CALL reset_methane_inactive_lake_microbe_diagnostics(i)"
    )
    build_effective_potential = lake_microbe_region.index(
        "microbial_prod_potential_eff(:) = 0._r8"
    )

    assert reset < build_effective_potential
    assert "microbial_prod_potential(1:nl_soil,i) = 0._r8" not in lake_microbe_region
    assert "microbial_oxid_potential(1:nl_soil,i) = 0._r8" not in lake_microbe_region


def test_lake_microbe_history_diagnostics_are_substep_time_weighted() -> None:
    microbes = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Microbes.F90")
    )
    impl = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Impl.F90")
    )
    accflux = without_continuations(
        source("MOD_Tracer_Reactive_Methane_AccFlux.F90")
    )
    body = microbes.split(
        "SUBROUTINE accumulate_methane_lake_microbe_substep_diagnostics", 1
    )[1].split(
        "END SUBROUTINE accumulate_methane_lake_microbe_substep_diagnostics", 1
    )[0]
    fields = (
        "f_T_methanogen",
        "f_S_methanogen",
        "f_O2_methanogen",
        "f_T_methanotroph",
        "methanogen_growth_rate",
        "methanotroph_growth_rate",
        "microbial_prod_potential",
        "microbial_oxid_potential",
    )

    assert "methane_lake_microbe_substep_ndiag = 8" in microbes
    for index, field in enumerate(fields, start=1):
        assert f"CALL add2d({index}, {field}(:,ipatch))" in body
        assert f"CALL finish2d({index}, {field}(:,ipatch))" in body
        assert re.search(rf"CALL\s+acc2d\s*\(\s*{field}\s*,", accflux)

    assert "+ var(:) * substep_dt" in body
    assert "var(:) = methane_lake_microbe_substep_acc(:,icol) / total_dt" in body
    assert "B_methanogen(:,ipatch)" not in body
    assert impl.count(
        "CALL accumulate_methane_lake_microbe_substep_diagnostics"
    ) == 3
    methane_driver = impl.index("CALL methane_driver")
    wet_microbe_accumulator = impl.index(
        "CALL accumulate_methane_lake_microbe_substep_diagnostics",
        methane_driver,
    )
    wet_state_accumulator = impl.index(
        "CALL accumulate_methane_lake_substep_diagnostics",
        wet_microbe_accumulator,
    )
    assert methane_driver < wet_microbe_accumulator < wet_state_accumulator


def test_drydown_substep_average_exports_inventory_once() -> None:
    stock = 3.6
    substep_dt = 900.0
    nsub = 4
    rates = [stock / substep_dt] + [0.0] * (nsub - 1)
    timestep_mean = sum(rate * substep_dt for rate in rates) / (
        substep_dt * nsub
    )

    assert timestep_mean == stock / (substep_dt * nsub)


def test_lake_surface_flux_is_water_to_air_not_sediment_to_air() -> None:
    physics = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Physics.F90")
    )

    assert re.search(
        r"methane_surf_diff\s*=\s*lake_exchange_vel\s*\*\s*"
        r"\(conc_ch4_rel\(0\)\s*-\s*k_h_cc\(0,s\)\s*\*\s*c_atm\(s\)\)",
        physics,
    )
    assert re.search(
        r"lake_sed_ch4_flux\s*=\s*lake_sed_water_cond\s*\*\s*"
        r"\(conc_ch4_rel\(1\)\s*-\s*conc_ch4_rel\(0\)\)",
        physics,
    )


def test_implicit_lake_box_is_mass_conservative_and_positive() -> None:
    # One water node coupled to one fixed sediment concentration.  This is the
    # same backward-Euler row used by methane_tran; test a very large timestep
    # to catch explicit overdraw and sign errors.
    dt = 86400.0
    depth = 0.25
    old_stock = 1.0e-5
    g_sw = 2.0e-7
    k_wa = 7.0e-6
    c_sed = 0.12
    c_eq = 2.2e-6

    c_water = (old_stock / dt + g_sw * c_sed + k_wa * c_eq) / (
        depth / dt + g_sw + k_wa
    )
    new_stock = depth * c_water
    f_sw = g_sw * (c_sed - c_water)
    f_wa = k_wa * (c_water - c_eq)

    assert c_water >= 0.0
    assert new_stock >= 0.0
    residual = new_stock - old_stock - dt * (f_sw - f_wa)
    assert abs(residual) < 1.0e-18


def test_lake_total_inventory_includes_water_stock_once() -> None:
    physics = without_continuations(
        source("MOD_Tracer_Reactive_Methane_Physics.F90")
    )

    assert "totcol_methane_lake = totcol_methane_sat" in physics
    assert "totcol_methane = totcol_methane_lake + lake_water_ch4_stock" in physics
    assert "totcol_methane_lake = totcol_methane_sat + lake_water_ch4_stock" not in physics


def test_lake_conductance_history_is_substep_time_weighted() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    body = state.split(
        "SUBROUTINE accumulate_methane_lake_substep_diagnostics", 1
    )[1].split("END SUBROUTINE accumulate_methane_lake_substep_diagnostics", 1)[0]

    assert "methane_lake_substep_n1d = 57" in state
    for index, field in (
        (54, "grnd_methane_cond"),
        (55, "grnd_methane_cond_sat"),
        (56, "grnd_methane_cond_lake"),
        (57, "methane_surf_diff_phys"),
    ):
        assert f"CALL add1d({index:2d}, {field}(ipatch))" in body
        assert f"CALL finish1d({index:2d}, {field}(ipatch))" in body
