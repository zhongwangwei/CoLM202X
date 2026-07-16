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

    assert "methane_lake_substep_n1d = 56" in state
    for index, field in (
        (54, "grnd_methane_cond"),
        (55, "grnd_methane_cond_sat"),
        (56, "grnd_methane_cond_lake"),
    ):
        assert f"CALL add1d({index:2d}, {field}(ipatch))" in body
        assert f"CALL finish1d({index:2d}, {field}(ipatch))" in body
