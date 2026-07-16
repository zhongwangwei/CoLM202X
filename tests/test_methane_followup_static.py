from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def source(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def test_lulcc_preserves_routing_inundation_inputs() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    for name in (
        "f_inund_levee_patch",
        "f_inund_flood_patch",
        "f_inund_flood_depth_patch",
    ):
        old = f"lulcc_{name}_old"
        assert f"allocate({old}" in state
        assert f"CALL remap1d({old}, {name})" in state
        assert f"allocated({old})" in state


def test_unimplemented_online_methane_is_rejected() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    validation = const.split("SUBROUTINE validate_methane_namelist", 1)[1].split(
        "END SUBROUTINE validate_methane_namelist", 1
    )[0]
    assert "IF (.not. DEF_METHANE%methane_offline) THEN" in validation
    assert "online atmosphere/NEE coupling that is not implemented" in validation


def test_core_history_exposes_numerical_flux_corrections() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    core = const.split("logical FUNCTION methane_history_is_core", 1)[1].split(
        "END FUNCTION methane_history_is_core", 1
    )[0]
    for name in (
        "f_methane_surf_flux_tot_active",
        "f_methane_surf_flux_tot_phys",
        "f_methane_balance_residual",
        "f_methane_ch4_clip_credit",
        "f_o2_cap_loss",
        "f_o2_cap_gain",
    ):
        assert f"'{name}'" in core


def test_default_flux_identity_uses_one_area_basis() -> None:
    hist = source("MOD_Tracer_Reactive_Methane_Hist.F90")
    for name in (
        "f_methane_surf_flux_tot_active",
        "f_methane_surf_flux_tot_phys",
        "f_methane_balance_residual",
        "f_methane_ch4_clip_credit",
    ):
        call = hist.split(f"'{name}'", 1)[1].split("acc_num=", 1)[1]
        assert call.lstrip().startswith("a_methane_acc_num")


def test_physical_flux_and_effective_conductance_include_all_physics() -> None:
    physics = source("MOD_Tracer_Reactive_Methane_Physics.F90")
    assert re.search(
        r"methane_surf_flux_tot_phys\s*=\s*methane_surf_diff_phys\s*\+\s*"
        r"methane_dfsat_tot",
        physics,
        re.DOTALL,
    )
    diff_phys = physics.split("methane_surf_diff_phys =", 1)[1].split(
        "methane_balance_residual =", 1
    )[0]
    assert "methane_dfsat_tot" not in diff_phys
    assert re.search(
        r"grnd_methane_cond\s*=\s*grnd_methane_cond_sat\s*\*\s*finundated\s*\+.*?"
        r"grnd_methane_cond_unsat",
        physics,
        re.DOTALL,
    )
    assert "grnd_methane_cond = grnd_methane_cond_lake" in physics


def test_diagnostic_flux_selector_keeps_corrected_active_identity() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    diagnostic = const.split("logical FUNCTION methane_history_is_diagnostic", 1)[1].split(
        "END FUNCTION methane_history_is_diagnostic", 1
    )[0]
    for name in (
        "f_methane_surf_flux_tot_active",
        "f_methane_surf_flux_tot_phys",
        "f_methane_balance_residual",
        "f_methane_ch4_clip_credit",
    ):
        assert f"'{name}'" in diagnostic


def test_lulcc_flood_sentinels_are_cleared_before_clamp() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    remap = state.split("SUBROUTINE remap_methane_lulcc_state", 1)[1]
    for name in (
        "f_inund_levee_patch",
        "f_inund_flood_patch",
        "f_inund_flood_depth_patch",
    ):
        assert f"ieee_is_nan({name})" in remap
