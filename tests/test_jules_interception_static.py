from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "extends" / "interception" / "MOD_LeafInterception_Extended.F90"
).read_text(encoding="utf-8")
JULES = SOURCE.split("SUBROUTINE LEAF_interception_JULES", 1)[1].split(
    "END SUBROUTINE LEAF_interception_JULES", 1
)[0]
LEAF_TEMPERATURE = (
    ROOT / "extends" / "interception" / "MOD_LeafTemperature_Extended.F90"
).read_text(encoding="utf-8")
LEAF_TEMPERATURE_PC = (
    ROOT / "extends" / "interception" / "MOD_LeafTemperaturePC_Extended.F90"
).read_text(encoding="utf-8")


def compact(source: str) -> str:
    return " ".join(source.lower().replace("&", "").split())


def test_jules_uses_full_precipitation_without_generic_prefilter():
    jules = compact(JULES)

    assert "default_interception_alpha" not in jules
    assert "fpi_pre" not in jules
    assert "rain_free_tf" not in jules
    assert "snow_free_tf" not in jules


def test_jules_rebuilds_total_rain_and_uses_upstream_convective_area():
    jules = compact(JULES)

    assert (
        "r_rain_con = 0.1_r8 * max(0.0_r8, prc_rain + prl_rain)" in jules
    )
    assert (
        "r_rain_ls = 0.9_r8 * max(0.0_r8, prc_rain + prl_rain) "
        "+ max(0.0_r8, qflx_irrig_sprinkler)" in jules
    )
    assert "area = 0.3_r8" in jules
    assert "area = 0.1_r8" not in jules


def test_jules_uses_default_zero_background_snow_unloading():
    jules = compact(JULES)

    assert "unload_backgrnd = 0.0_r8" in jules
    assert "2.31e-6" not in jules
    assert "5.56e-7" not in jules


def test_jules_warm_rain_evaporation_uses_step_start_storage():
    single = compact(LEAF_TEMPERATURE)
    pc = compact(LEAF_TEMPERATURE_PC)

    assert (
        "ldew_jules_rain_evap = min(ldew_rain, "
        "max(0._r8, ldew_rain - qintr_rain * deltim))" in single
    )
    assert "ldew_vic_evap = ldew_jules_rain_evap" in single
    assert "ldew_jules_rain_evap / sigf_safe_jules" in single

    assert (
        "ldew_jules_rain_evap(i) = min(ldew_rain(i), "
        "max(0._r8, ldew_rain(i) - qintr_rain(i) * deltim))" in pc
    )
    assert "ldew_vic_evap(i) = ldew_jules_rain_evap(i)" in pc
    assert "ldew_jules_rain_evap(i) / sigf_safe_jules" in pc
