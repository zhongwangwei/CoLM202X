from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def test_campbell_tracer_conflict_fails_fast_without_silently_undefining_tracer():
    define = text("include/define.h")
    conflict = define.split("#if (defined TRACER) && (defined Campbell_SOIL_MODEL)", 1)[1].split("#endif", 1)[0]
    assert "#error" in conflict
    assert "#undef TRACER" not in conflict
    assert "This repository template currently enables TRACER" in define
    assert "Default is OFF" not in define


def test_isogsm_forcing_template_uses_current_per_species_tracer_forcing_schema():
    isogsm = text("run/forcing/IsoGSM.nml")
    for retired in (
        "precipitation_O18",
        "precipitation_H2",
        "water_vapor_O18",
        "water_vapor_H2",
    ):
        assert retired not in isogsm
    for current in ("IsoGSM_temperature", "IsoGSM_Q", "IsoGSM_prate"):
        assert current in isogsm
    for parameter_file in (
        "run/standard_O18_parameter.nml",
        "run/standard_HDO_parameter.nml",
    ):
        parameters = text(parameter_file)
        assert "forcing_role       = 'precip', 'vapor'" in parameters
        assert "forcing_fprefix    = 'IsoGSM_prate', 'IsoGSM_Q'" in parameters


def test_readme_states_ch4_is_offline_only_for_atmosphere_not_land_bgc():
    readme = text("main/TRACER/README.md")
    assert "offline only with respect to atmospheric CH4 feedback" in readme
    flat = " ".join(readme.split())
    assert "coupled to land BGC carbon inputs" in flat
    assert "updates methane-specific land state" in flat


def test_amazon_template_has_no_nonexistent_runtime_tracer_switches():
    amazon = text("run/Amazon.nml")
    assert "DEF_USE_TRACER" not in amazon
    assert "DEF_USE_SEDIMENT" not in amazon


def test_ci_define_generator_and_matrix_cover_current_tracer_contract():
    generator = text(".github/workflows/create_defineh.bash")
    cases = text(".github/workflows/TestCaseLists")
    assert "#define extend_interception" in generator
    assert "#if (defined TRACER) && (defined Campbell_SOIL_MODEL)" in generator
    assert '#error "TRACER requires vanGenuchten_Mualem_SOIL_MODEL' in generator
    tracer_only = next(line for line in cases.splitlines() if "GRID_PFT_TRACER_ONLY" in line)
    assert "BGCOFF" in tracer_only
    assert tracer_only.endswith("TRACERON")


def test_tracer_forcing_final_releases_parsed_parameter_state():
    forcing = text("main/TRACER/MOD_Tracer_Forcing.F90")
    final = forcing.split("SUBROUTINE tracer_forcing_final", 1)[1].split(
        "END SUBROUTINE tracer_forcing_final", 1
    )[0]
    assert "CALL tracer_forcing_input_final()" in final
