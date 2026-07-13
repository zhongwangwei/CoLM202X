from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFS = (ROOT / "main" / "TRACER" / "MOD_Tracer_Defs.F90").read_text(
    encoding="utf-8"
)
FORCING = (ROOT / "main" / "TRACER" / "MOD_Tracer_Forcing.F90").read_text(
    encoding="utf-8"
)
FORCING_INPUT = (
    ROOT / "main" / "TRACER" / "MOD_Tracer_ForcingInput.F90"
).read_text(encoding="utf-8")
METHANE_REGISTRY = (
    ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Registry.F90"
).read_text(encoding="utf-8")
CONSERVATION = (
    ROOT / "main" / "TRACER" / "MOD_Tracer_Conservation.F90"
).read_text(encoding="utf-8")
CHLORIDE = (ROOT / "run" / "standard_chloride_parameter.nml").read_text(
    encoding="utf-8"
)


def function_body(source, name):
    return source.split(f"FUNCTION {name}", 1)[1].split(f"END FUNCTION {name}", 1)[0]


def subroutine_body(source, name):
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_solute_initial_and_forcing_defaults_are_separate():
    for field in ("init_conc", "precip_default_conc", "vapor_default_conc"):
        assert f"real(r8)          :: {field}" in DEFS
        assert f"DEF_TRACER%{field}" in CHLORIDE

    prepare = subroutine_body(FORCING, "tracer_forcing_prepare_step")
    assert "tracer_precip_default_ratio(itrc)" in prepare
    assert "tracer_vapor_default_ratio(itrc)" in prepare
    assert "tracer_init_water_ratio(itrc)" not in prepare


def test_legacy_init_delta_is_the_per_field_fallback():
    reader = subroutine_body(DEFS, "read_tracer_parameter_file")
    compact = "".join(reader.split())
    for field in ("init_conc", "precip_default_conc", "vapor_default_conc"):
        assert f"tracers(itrc)%{field}=DEF_TRACER%init_delta" in compact
        assert f"tracers(itrc)%{field}=DEF_TRACER%{field}" in compact


def test_nonvolatile_capability_is_tied_to_generic_land_water_state():
    helper = function_body(DEFS, "tracer_is_nonvolatile_solute")
    assert "tracer_uses_land_water_transport(itrc)" in helper
    assert "tracers(itrc)%is_nonvolatile" in helper
    assert "tracer_is_conservative(itrc)" not in helper

    validator = subroutine_body(DEFS, "validate_tracer_descriptor")
    assert "is only implemented for category=conservative" not in validator
    assert "tracer_is_conservative(itrc) .or. tracer_is_reactive(itrc)" in validator
    assert "requires generic land-water transport" in validator


def test_reactive_decay_includes_dry_residue_in_source_sink_ledger():
    body = subroutine_body(CONSERVATION, "tracer_apply_reactive_processes")
    assert "CALL decay_pool(trc_surface_residue(itrc, ipatch)" in body
    assert "CALL decay_pool(trc_subsurface_residue(itrc, ipatch)" in body
    assert body.index("CALL decay_pool(trc_subsurface_residue(itrc, ipatch)") < body.index(
        "trc_reactive_source_step(itrc, ipatch)"
    )


def test_public_category_helpers_guard_allocation_and_bounds_explicitly():
    for name in (
        "tracer_is_isotope",
        "tracer_is_conservative",
        "tracer_is_reactive",
        "tracer_is_particle",
        "tracer_is_nonvolatile_solute",
    ):
        helper = function_body(DEFS, name)
        assert "IF (.not. allocated(tracers)) RETURN" in helper
        assert "IF (itrc < 1 .or. itrc > ntracers) RETURN" in helper
        assert "allocated(tracers) .and." not in helper


def test_forcing_namelist_absence_and_parse_failure_are_distinguished():
    loader = subroutine_body(FORCING_INPUT, "tracer_forcing_input_load")
    assert "tracer_forcing_group_present" in loader
    assert "iomsg=iomsg" in loader.replace(" ", "")
    assert "CALL CoLM_stop()" in loader
    assert "IF (ierr /= 0) CYCLE" not in loader


def test_forcing_dtime_is_validated_before_use():
    loader = subroutine_body(FORCING_INPUT, "tracer_forcing_input_load")
    add_var = subroutine_body(FORCING, "tracer_forcing_add_var")
    stamp_lb = subroutine_body(FORCING, "tracer_forcing_setstamp_LB")
    assert "forcing_dtime(k) <= 0" in loader
    assert "dtime <= 0" in add_var
    assert "tracer_forcing_validate_var_dtime" in stamp_lb


def test_methane_registry_has_no_unwired_o2_or_co2_placeholders():
    assert "igas_" + "o2" not in METHANE_REGISTRY
    assert "igas_" + "co2" not in METHANE_REGISTRY
