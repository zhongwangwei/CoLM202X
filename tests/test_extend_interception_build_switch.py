from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAKEFILE = (ROOT / "Makefile").read_text(encoding="utf-8")
COLMMAIN = (ROOT / "main" / "CoLMMAIN.F90").read_text(encoding="utf-8")
THERMAL_EXTENDED = (
    ROOT / "extends" / "interception" / "MOD_Thermal_CanopyPhase_Extended.F90"
).read_text(encoding="utf-8")
INTERCEPTION_EXTENDED = (
    ROOT / "extends" / "interception" / "MOD_LeafInterception_Extended.F90"
).read_text(encoding="utf-8")


def test_makefile_selects_all_four_extended_interception_modules():
    assert "EXTENDED_INTERCEPTION_ENABLED" in MAKEFILE
    assert "MOD_LeafInterception_Extended.F90" in MAKEFILE
    assert "MOD_LeafTemperature_Extended.F90" in MAKEFILE
    assert "MOD_LeafTemperaturePC_Extended.F90" in MAKEFILE
    assert "MOD_Thermal_CanopyPhase_Extended.F90" in MAKEFILE
    assert "$(filter-out $(INTERCEPTION_CORE_OBJS),$(OBJS_MAIN))" in MAKEFILE


def test_colmmain_passes_canopy_phase_heat_only_to_extended_thermal():
    start = COLMMAIN.index("CALL THERMAL")
    call = COLMMAIN[start : COLMMAIN.index("#ifdef TRACER", start)]

    assert "#ifdef extend_interception" in call
    assert ",canopy_phase_heat, canopy_phase_heat_p" in call
    assert ",canopy_phase_heat" in call


def test_colmmain_passes_ground_snow_fraction_only_to_extended_interception():
    start = COLMMAIN.index("CALL LEAF_interception_wrap")
    call = COLMMAIN[start : COLMMAIN.index("canopy_phase_heat)", start)]

    assert "#ifdef extend_interception" in call
    assert "fsno," in call


def test_extended_thermal_initializes_empty_pft_end_index():
    loop = THERMAL_EXTENDED.index("DO i = ps, pe", THERMAL_EXTENDED.index("pn = ps - 1"))
    assert THERMAL_EXTENDED.index("pn = ps - 1") < loop


def test_extended_thermal_guards_ground_phase_outputs_with_tracer():
    start = THERMAL_EXTENDED.index("CALL GroundTemperature")
    call = THERMAL_EXTENDED[start : THERMAL_EXTENDED.index(")", start) + 1]
    assert "#ifdef TRACER" in call
    assert "qphs_thaw_lay = qphs_thaw_lay_th" in call


def test_extended_colm2014_keeps_legacy_urban_call_compatible():
    declaration = INTERCEPTION_EXTENDED[
        INTERCEPTION_EXTENDED.index("SUBROUTINE LEAF_interception_CoLM2014") :
        INTERCEPTION_EXTENDED.index("IF (lai+sai", INTERCEPTION_EXTENDED.index("SUBROUTINE LEAF_interception_CoLM2014"))
    ]
    assert declaration.count("intent(out), optional") == 7
