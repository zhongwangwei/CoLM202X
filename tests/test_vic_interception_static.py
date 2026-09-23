from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INTERCEPTION = (
    ROOT / "extends" / "interception" / "MOD_LeafInterception_Extended.F90"
).read_text(encoding="utf-8")
LEAF_TEMPERATURE = (
    ROOT / "extends" / "interception" / "MOD_LeafTemperature_Extended.F90"
).read_text(encoding="utf-8")
LEAF_TEMPERATURE_PC = (
    ROOT / "extends" / "interception" / "MOD_LeafTemperaturePC_Extended.F90"
).read_text(encoding="utf-8")
VIC = INTERCEPTION.split("SUBROUTINE LEAF_interception_VIC", 1)[1].split(
    "END SUBROUTINE LEAF_interception_VIC", 1
)[0]


def compact(source: str) -> str:
    return " ".join(source.lower().replace("&", "").split())


def test_vic_uses_full_precipitation_and_current_mixed_phase_capacity():
    vic = compact(VIC)

    assert "fpi_pre" not in vic
    assert "rain_free_tf" not in vic
    assert "snow_free_tf" not in vic
    assert "p0 = (rain_clamp + snow_clamp) * deltim" in vic
    assert "maxwaterint = 0.035_r8 * ldew_snow + maxint" in vic
    assert "tex_rain = max(0.0,ldew_rain-maxwaterint)" in vic


def test_vic_wet_fraction_keeps_phase_aware_dewfraction_result():
    for source in (LEAF_TEMPERATURE, LEAF_TEMPERATURE_PC):
        text = compact(source)

        assert "wetfrac_vic = wetfrac_vic_base" in text
        assert "sigf_safe_vic" not in text
        assert "maxint_vic" not in text
        assert (
            "(def_veg_snow .or. def_interception_scheme == 6) .and. "
            "def_interception_scheme /= 2"
        ) in text


def test_vic_snow_capacity_and_unloading_use_canopy_conditions():
    vic = compact(VIC)
    source = compact(INTERCEPTION)

    assert "if (tleaf > 272.15_r8)" in vic
    assert "lr = 1.5_r8*(tleaf - 273.15_r8) + 5.5_r8" in vic
    assert "if (tair > 272.15_r8)" not in vic
    assert "wind = vic_canopy_wind_speed" in vic
    assert "if (tleaf-tfrz < -3.0_r8" in vic
    assert "wind_attenuation = 0.5_r8" in source
    assert "z0m_height_ratio = 0.1_r8" in source


def test_vic_dynamic_snow_capacity_only_limits_new_interception():
    vic = compact(VIC)

    # VIC's MaxSnowInt limits DeltaSnowInt; existing canopy snow is retained
    # unless the separate structural limit Imax1 is exceeded.
    assert "xsc_snow = max(0., ldew_snow-satcap_snow)" not in vic
    assert "xsc_snow = xsc_snow + max(0., ldew_snow-satcap_snow)" not in vic
    assert "tex_snow = max(0.0,ldew_snow-satcap_snow)" not in vic

    # Structural unloading must also run when no new precipitation arrives.
    precip_end = vic.index("endif ! end precipitation-dependent interception")
    overload = vic.index("if (ldew_rain + ldew_snow > imax1) then")
    assert overload > precip_end


def test_vic_processes_every_positive_precipitation_amount():
    vic = compact(VIC)

    assert "if (p0 > 0._r8) then" in vic
    assert "if (p0 > 1.e-8) then" not in vic


def test_vic_liquid_thin_storage_cutoff_requires_snow_context():
    vic = compact(VIC)

    assert "vic_snow_intercept_context =" in vic
    assert "fsno > 1.e-10_r8" in vic
    assert (
        "if (vic_snow_intercept_context .and. rain_clamp < 1.e-8_r8 "
        ".and. ldew_rain < 1.0_r8"
    ) in vic


def test_vic_runs_inside_the_vegetation_tile_with_f_equal_one():
    vic = compact(VIC)

    assert "ldew_rain = ldew_rain / sigf" not in vic
    assert "ldew_snow = ldew_snow / sigf" not in vic
    assert "sigf_safe" not in vic
    assert "tti_snow = snow - deltasnowint" in vic
    assert "tti_rain = 0._r8" in vic
    assert "pg_rain = (xsc_rain + thru_rain) / deltim" in vic
    assert "gross_intr_rain = actual_rain_int / deltim" in vic

    for source in (LEAF_TEMPERATURE, LEAF_TEMPERATURE_PC):
        text = compact(source)
        assert "satcap_rain_eff = 0.1_r8 * max(lai_perveg, 0._r8)" in text
        assert "satcap_snow_eff = 0.5_r8 * lr * max(lai_perveg, 0._r8)" in text


def test_vic_all_timestep_rain_evaporation_uses_step_start_storage():
    for source in (LEAF_TEMPERATURE, LEAF_TEMPERATURE_PC):
        text = compact(source)

        assert "def_interception_scheme == 6 .and. deltim < 86400._r8" not in text
        assert "fsno <= 1.e-10_r8" in text
        assert "abs(canopy_phase_heat" in text
        assert "ldew_vic_evap" in text
        assert "ldew_rain" in text and "- qintr_rain" in text
        assert "if (deltim >= 86400._r8) ldew_vic_evap" in text
        assert "canopy_rain_capacity_for_fwet" in text
        assert "elwmax = ldew_vic_evap" in text
