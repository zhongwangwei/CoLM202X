"""Non-isotope species must not index the isotope registry at zero."""

from pathlib import Path


SOURCE = (
    Path(__file__).resolve().parents[1]
    / "main/TRACER/MOD_Tracer_Isotope_Registry.F90"
).read_text()


def test_unregistered_species_returns_before_dereferencing_registry():
    for name in (
        "isotope_alpha_liq_vap",
        "isotope_alpha_ice_vap",
        "isotope_diffusivity_ratio_air",
        "isotope_leaf_liquid_diffusivity",
    ):
        body = SOURCE.split(f"FUNCTION {name}", 1)[1].split(f"END FUNCTION {name}", 1)[0]
        assert body.index("IF (idx <= 0) RETURN") < body.index(
            "associated(isotope_physics(idx)%"
        )


def test_zero_ratio_hint_is_skipped_before_division():
    body = SOURCE.split("FUNCTION find_isotope_physics", 1)[1].split(
        "END FUNCTION find_isotope_physics", 1
    )[0]
    assert body.index("IF (isotope_physics(i)%ref_ratio_hint <= trc_tiny) CYCLE") < body.index(
        "/ &\n                isotope_physics(i)%ref_ratio_hint"
    )
