from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PHYSICS = ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Physics.F90"
CONST = ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Const.F90"
DRIVER = ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Driver.F90"
IMPL = ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Impl.F90"
STANDARD_PARAM = ROOT / "run" / "standard_ch4_parameter.nml"


def source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def without_continuations(text: str) -> str:
    return re.sub(r"&\s*\n\s*", " ", text)


def test_transpiration_uses_true_aqueous_concentration() -> None:
    physics = without_continuations(source(PHYSICS))

    assert "tranloss(j) = conc_ch4_aqu_porsl(j) * rootr(j)*etr" in physics
    assert "tranloss(j) = conc_ch4_aqu(j) * rootr(j)*etr" not in physics


def test_fsat_decrease_conservatively_regrids_oxygen() -> None:
    physics = without_continuations(source(PHYSICS))

    assert re.search(
        r"conc_o2_unsat\(j\)\s*=\s*\(\(1\._r8\s*-\s*fsat_bef\)\s*\*\s*"
        r"conc_o2_unsat\(j\)\s*-\s*dfsat\s*\*\s*conc_o2_sat\(j\)\)\s*/\s*"
        r"\(1\._r8\s*-\s*finundated\)",
        physics,
        re.DOTALL,
    )
    assert "wetland_inactive_dry_area" not in physics
    assert "zwt_unsat   = zi_soisno(nl_soil) + 1.e-3_r8" in physics
    assert "jwt_unsat   = nl_soil" in physics

    fsat_old, fsat_new = 0.7, 0.25
    o2_sat, o2_unsat = 0.2, 1.1
    dfsat = fsat_new - fsat_old
    o2_unsat_new = (
        (1.0 - fsat_old) * o2_unsat - dfsat * o2_sat
    ) / (1.0 - fsat_new)
    inventory_old = fsat_old * o2_sat + (1.0 - fsat_old) * o2_unsat
    inventory_new = fsat_new * o2_sat + (1.0 - fsat_new) * o2_unsat_new
    assert abs(inventory_new - inventory_old) < 1.0e-14


def test_rice_flood_fraction_does_not_create_surface_water() -> None:
    physics = without_continuations(source(PHYSICS))

    assert "rice_ponded_depth_mm" not in physics
    assert "rice_paddy_min_finundated" not in physics
    assert "finundated_rice" not in physics
    assert "wdsrf_sat = methane_wetland_water_depth(patchtype, wdsrf, wetwat)" in physics


def test_surface_conductance_uses_current_monin_obukhov_state_with_fallback() -> None:
    physics = source(PHYSICS)
    driver = source(DRIVER)
    impl = source(IMPL)

    assert "ustar_in" in physics
    assert "fq_in" in physics
    assert re.search(
        r"grnd_methane_cond_base\s*=\s*vonkar\s*\*\s*ustar_in\s*/\s*fq_in",
        physics,
    )
    assert "patchtype == 4 .or. lai + sai <= 1.e-6_r8" in without_continuations(physics)
    assert "fq-fqt plus canopy/ground resistances" in physics
    assert re.search(
        r"grnd_methane_cond_base\s*=\s*DEF_METHANE%grnd_methane_cond_default",
        physics,
    )
    assert "ustar_in=ustar, fq_in=fq" in without_continuations(driver)
    assert impl.count("sai(ipatch), rootr(1:nl_soil,ipatch)") == 2
    assert impl.count("ustar(ipatch), fq(ipatch)") == 2


def test_walter_renormalization_claim_matches_implemented_weights() -> None:
    physics = source(PHYSICS)
    walter = physics.split("Walter & Heimann 2001 depth-attenuation pre-pass", 1)[1]
    walter = walter.split("column loop to partition decomposition_rate", 1)[0]

    assert "preserves the unmodified HR partition-weight sum" in walter
    assert "preserves total column CH4 production" not in walter


def test_rice_calendar_is_used_for_physiology_not_hydrology() -> None:
    physics = source(PHYSICS)
    driver = source(DRIVER)

    assert "rice_days_past_planting(ipatch" not in physics
    assert "rice_days_since_harvest(ipatch" not in physics
    assert "rice_days_since_harvest(i, idate(2), idate(1))" in driver


def test_material_numerical_corrections_fail_fast_in_production() -> None:
    physics = without_continuations(source(PHYSICS))
    const = source(CONST)
    standard = source(STANDARD_PARAM)

    assert "real(r8) :: numerical_correction_fatal_threshold = -1._r8" in const
    assert "DEF_METHANE%numerical_correction_fatal_threshold" in standard
    assert "= 1.0e-3" in standard

    threshold = "DEF_METHANE%numerical_correction_fatal_threshold"
    assert physics.count(f"{threshold} > 0._r8") == 4
    assert f"abs(err_methane) > {threshold}" in physics
    assert "ch4_clip_deficit_col = 0._r8" in physics
    assert "ch4_clip_deficit_col = ch4_clip_deficit_col + deficit" in physics
    assert f"ch4_clip_deficit_col > {threshold}" in physics
    assert f"o2_cap_gain_col > {threshold}" in physics
    assert "CH4 column-budget correction exceeds fatal threshold" in physics
    assert "CH4 nonnegative column clip exceeds fatal threshold" in physics
    assert "O2 nonnegative floor exceeds fatal threshold" in physics
    assert "CH4 transport closure residual exceeds fatal threshold" in physics
