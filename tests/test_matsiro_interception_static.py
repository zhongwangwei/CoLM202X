from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "extends" / "interception" / "MOD_LeafInterception_Extended.F90"
).read_text(encoding="utf-8")
LEAF_TEMPERATURE = (
    ROOT / "extends" / "interception" / "MOD_LeafTemperature_Extended.F90"
).read_text(encoding="utf-8")
LEAF_TEMPERATURE_PC = (
    ROOT / "extends" / "interception" / "MOD_LeafTemperaturePC_Extended.F90"
).read_text(encoding="utf-8")
MATSIRO = SOURCE.split("SUBROUTINE LEAF_interception_MATSIRO", 1)[1].split(
    "END SUBROUTINE LEAF_interception_MATSIRO", 1
)[0]


def compact(source: str) -> str:
    return " ".join(source.lower().replace("&", "").split())


def test_matsiro_rebuilds_total_precipitation_as_10_90_split():
    matsiro = compact(MATSIRO)

    assert (
        "rain_conv_clamp = 0.1_r8 * max(0.0_r8, prc_rain + prl_rain)"
        in matsiro
    )
    assert (
        "rain_strat_clamp = 0.9_r8 * max(0.0_r8, prc_rain + prl_rain) "
        "+ max(0.0_r8, qflx_irrig_sprinkler)"
        in matsiro
    )
    assert (
        "snow_conv_clamp = 0.1_r8 * max(0.0_r8, prc_snow + prl_snow)"
        in matsiro
    )
    assert (
        "snow_strat_clamp = 0.9_r8 * max(0.0_r8, prc_snow + prl_snow)"
        in matsiro
    )
    assert "rain_conv_clamp = max(0.0_r8, prc_rain)" not in matsiro
    assert "rain_strat_clamp = max(0.0_r8, prl_rain" not in matsiro


def test_matsiro_uses_full_precipitation_without_generic_prefilter():
    matsiro = compact(MATSIRO)

    assert "if (lai > 1e-6_r8) then" in matsiro
    assert "default_interception_alpha" not in matsiro
    assert "fpi_pre" not in matsiro
    assert "rain_free_tf" not in matsiro
    assert "snow_free_tf" not in matsiro


def test_matsiro_evaporation_uses_step_start_storage_and_fwet_weight():
    leaf = compact(LEAF_TEMPERATURE)
    leaf_pc = compact(LEAF_TEMPERATURE_PC)
    marker = "elseif (def_interception_scheme == 5) then"
    leaf_scheme5 = leaf.split(marker, 1)[1].split(" else", 1)[0]
    leaf_pc_scheme5 = leaf_pc.rsplit(marker, 1)[1].split(" else", 1)[0]

    assert "ldew - (qintr_rain + qintr_snow) * deltim" in leaf
    assert "ldew(i) - (qintr_rain(i) + qintr_snow(i)) * deltim" in leaf_pc
    assert "evp_weight = fwet" in leaf_scheme5
    assert "evp_weight = fwet(i)" in leaf_pc_scheme5
    assert "1._r8 - delta" not in leaf_scheme5
    assert "1._r8 - delta" not in leaf_pc_scheme5
    assert "cfw = fwet*wet_cond_cfw" in leaf
    assert "cfw(i) = fwet(i)*wet_cond_cfw" in leaf_pc
    assert "ldew_vic_evap = min(ldew, ldew_vic_evap)" in leaf
    assert "ldew_vic_evap(i) = min(ldew(i), ldew_vic_evap(i))" in leaf_pc
