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
NOAHMP = SOURCE.split("SUBROUTINE LEAF_interception_NOAHMP", 1)[1].split(
    "END SUBROUTINE LEAF_interception_NOAHMP", 1
)[0]


def compact(source: str) -> str:
    return " ".join(source.lower().replace("&", "").split())


def test_noahmp_receives_full_precipitation_before_vegetation_gap_partition():
    noahmp = compact(NOAHMP)

    assert "default_interception_alpha" not in noahmp
    assert "fpi_pre" not in noahmp
    assert "rain_free_tf" not in noahmp
    assert "snow_free_tf" not in noahmp
    assert (
        "rain_clamp = max(0.0_r8, prc_rain + prl_rain + qflx_irrig_sprinkler)"
        in noahmp
    )
    assert "snow_clamp = max(0.0_r8, prc_snow + prl_snow)" in noahmp
    assert "p0 = (rain_clamp + snow_clamp) * deltim" in noahmp
    assert (
        "ppc = 0.1_r8 * max(0.0_r8, prc_rain + prc_snow + prl_rain + prl_snow) * deltim"
        in noahmp
    )
    assert "ppl = max(0.0_r8, p0 - ppc)" in noahmp


def test_noahmp_precipitation_area_fraction_matches_upstream_formula():
    noahmp = compact(NOAHMP)

    assert "precipareafrac = p0 / (10.0_r8*ppc + ppl)" in noahmp
    sprinkler_branch = noahmp.index("if (qflx_irrig_sprinkler > 0._r8) then")
    uniform_area = noahmp.index("precipareafrac = 1.0_r8", sprinkler_branch)
    formula = noahmp.index("precipareafrac = p0 / (10.0_r8*ppc + ppl)")
    assert sprinkler_branch < uniform_area < formula

    def precip_area_fraction(convective: float, large_scale: float) -> float:
        total = convective + large_scale
        return total / (10.0 * convective + large_scale)

    assert precip_area_fraction(1.0, 0.0) == 0.1
    assert precip_area_fraction(0.0, 1.0) == 1.0
    assert precip_area_fraction(1.0, 1.0) == 2.0 / 11.0


def test_noahmp_drains_existing_storage_above_current_capacity():
    noahmp = compact(NOAHMP)

    assert "xsc_rain = max(0._r8, ldew_rain - satcap_rain)" in noahmp
    assert "xsc_snow = max(0._r8, ldew_snow - satcap_snow)" in noahmp
    assert (
        "xsc_rain = xsc_rain + max(0._r8, ldew_rain - satcap_rain)"
        in noahmp
    )
    assert (
        "xsc_snow = xsc_snow + max(0._r8, ldew_snow - satcap_snow)"
        in noahmp
    )


def test_noahmp_snow_unloading_source_uses_existing_storage_only():
    noahmp = compact(NOAHMP)

    assert "icedrip = max(0._r8, ldew_snow) * (fv+ft)" in noahmp
    assert "max(0._r8, ldew_snow) + int_snow*deltim" not in noahmp


def test_noahmp_processes_every_positive_precipitation_amount():
    noahmp = compact(NOAHMP)

    assert "if (p0 > 0._r8) then" in noahmp
    assert "if (p0 > 1.e-8) then" not in noahmp


def test_noahmp_wet_canopy_flux_keeps_outer_vegetation_fraction():
    single = compact(LEAF_TEMPERATURE)
    pc = compact(LEAF_TEMPERATURE_PC)

    assert (
        "wet_area = max(0.05_r8, 1._r8 - exp(-0.52_r8 * max(lai + sai, 0._r8))) * wet_area_cfw"
        in single
    )
    assert (
        "wet_area = max(0.05_r8, 1._r8 - exp(-0.52_r8 * max(lsai(i), 0._r8))) * wet_area_cfw"
        in pc
    )


def test_noahmp_wet_fraction_uses_ice_when_any_canopy_ice_is_present():
    for source in (LEAF_TEMPERATURE, LEAF_TEMPERATURE_PC):
        text = compact(source)

        assert "if (ldew_snow > 0._r8) then" in text
        assert "(ldew_snow / satcap_rain_eff)**.666666666666_r8" in text
        assert "elseif (ldew_rain > 0._r8) then" in text
        assert "(ldew_rain / satcap_rain_eff)**.666666666666_r8" in text
        assert "ldew_snow >= ldew_rain" not in text
