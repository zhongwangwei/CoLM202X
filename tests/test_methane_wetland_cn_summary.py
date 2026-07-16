from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def src(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_nonvegetated_wetland_cn_summary_is_called_after_state_updates() -> None:
    summary = src("main/BGC/MOD_BGC_CNSummary.F90")
    bridge = src("main/TRACER/MOD_Tracer_Reactive_Methane_BgcLink.F90")

    assert "PUBLIC CNDriverSummarizeNonvegetatedSoilStates" in summary
    wrapper = summary.split(
        "SUBROUTINE CNDriverSummarizeNonvegetatedSoilStates", 1
    )[1].split("END SUBROUTINE CNDriverSummarizeNonvegetatedSoilStates", 1)[0]
    assert "CALL soilbiogeochem_carbonstate_summary" in wrapper
    assert "CALL soilbiogeochem_nitrogenstate_summary" in wrapper
    assert "totvegc(i) = 0._r8" in wrapper
    assert "totcolc(i) = totcwdc(i) + totlitc(i) + totsomc(i) + ctrunc_soil(i)" in wrapper
    assert "totcoln(i) = totcwdn(i) + totlitn(i) + totsomn(i) + sminn(i) + ntrunc_soil(i)" in wrapper

    finalize = bridge.split("SUBROUTINE tracer_ch4_bgc_finalize_step", 1)[1].split(
        "END SUBROUTINE tracer_ch4_bgc_finalize_step", 1
    )[0]
    assert finalize.index("CDecompStateUpdate") < finalize.index(
        "CNDriverSummarizeNonvegetatedSoilStates"
    )
    assert finalize.index("SoilBiogeochemNDecompStateUpdate") < finalize.index(
        "CNDriverSummarizeNonvegetatedSoilStates"
    )


def test_wetland_cn_summary_is_initialized_and_make_dependencies_are_explicit() -> None:
    init = src("mkinidata/MOD_IniTimeVariable.F90")
    makefile = src("Makefile")

    assert "CNDriverSummarizeNonvegetatedSoilStates" in init
    assert "IF (patchtype == 2)" in init
    target = "MOD_Tracer_Reactive_Methane_BgcLink.o: MOD_Tracer_Reactive_Methane_Const.o"
    assert makefile.count(target) == 1
    dependency = makefile.split(target, 1)[1].split("\nMOD_", 1)[0]
    assert "MOD_BGC_CNSummary.o" in dependency


def test_wetland_competition_uses_callers_actual_positive_timestep() -> None:
    shim = src("main/TRACER/MOD_Tracer_Reactive_BgcShim.F90")

    assert "ieee_is_finite(deltim)" in shim
    assert "deltim <= 0._r8" in shim
    assert "SoilBiogeochemCompetitionNoPlant (ipatch, deltim" in shim
    assert "DEF_simulation_time%timestep" not in shim
