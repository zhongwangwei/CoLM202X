"""Alignment of the land isotope fractionation physics with the IsoGSM
reference model (codes/IsoGSM: gsml/ISOTOPE/{freq,frkin,eqm_deg}.F,
gsml/CLD1/lrgscl.F, gsml/moninp.F).

The equilibrium coefficients already match IsoGSM exactly.  These tests pin
down the kinetic/diffusive parts, where the two models had genuinely
different physics.
"""

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


def text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


NAMELIST = text("share/MOD_Namelist.F90")
AMAZON = text("run/Amazon_iso_0p5_2020.nml")
FRAC = text("main/TRACER/MOD_Tracer_Frac.F90")
REGISTRY = text("main/TRACER/MOD_Tracer_Isotope_Registry.F90")
O18 = text("main/TRACER/MOD_Tracer_Isotope_O18.F90")
HDO = text("main/TRACER/MOD_Tracer_Isotope_HDO.F90")
SPECIAL = text("main/TRACER/MOD_Tracer_SpecialPatches.F90")
FORCING = text("main/TRACER/MOD_Tracer_Forcing.F90")
MAIN = text("main/CoLMMAIN.F90")


# ---------------------------------------------------------------------------
# Diffusivity-ratio scheme default (IsoGSM lrgscl.F:103-105 uses the
# Merlivat 1978 pair dif18o=1.02849 / difhdo=1.02512).
# ---------------------------------------------------------------------------


def test_default_kinetic_scheme_reproduces_the_isogsm_reference_experiment():
    """CAPPA2003 gives eps_k(HDO)/eps_k(18O)=0.51 vs Merlivat's 0.88.

    The shipped default retains the Merlivat pair for reproducibility with
    IsoGSM.  Cappa remains a valid land-surface sensitivity choice; the
    atmospheric forcing producer does not physically constrain this knob.
    """
    decl = re.search(
        r"character\(len=16\)\s*::\s*DEF_TRACER_KINETIC_SCHEME\s*=\s*'([A-Za-z0-9]+)'",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_KINETIC_SCHEME declaration not found"
    assert decl.group(1) == "MERLIVAT1978"


def test_merlivat_default_is_justified_in_place():
    """The choice is a science decision, not a typo: keep the reason next
    to the declaration so it survives future edits."""
    decl_idx = NAMELIST.index("DEF_TRACER_KINETIC_SCHEME =")
    context = NAMELIST[max(0, decl_idx - 600) : decl_idx]
    assert "IsoGSM" in context
    assert "d-excess" in context or "deuterium excess" in context


# ---------------------------------------------------------------------------
# Vapour -> ice deposition: Jouzel & Merlivat (1984) supersaturation kinetics
# (IsoGSM lrgscl.F:196-224).  The numerical properties of the closure live in
# test_tracer_isotope_frac_runtime.py; these tests pin the wiring, which a
# standalone-compiled function cannot see.
# ---------------------------------------------------------------------------


def test_ice_deposition_uses_supersaturation_alpha_not_bare_equilibrium():
    """The from_ice branch of the deposition ratio must go through the
    JM84 effective alpha; using alpha_ice_vap directly over-enriches every
    frost/rime event below 0 C."""
    branch = FRAC.split("FUNCTION tracer_equilibrium_deposition_ratio", 1)[1].split(
        "END FUNCTION tracer_equilibrium_deposition_ratio", 1
    )[0]
    from_ice = branch.split("IF (from_ice) THEN", 1)[1].split("ELSE", 1)[0]
    assert "tracer_alpha_ice_vap_deposition(itrc, temp_k)" in from_ice
    assert "tracer_alpha_ice_vap(itrc, temp_k) * vapor_ratio" not in from_ice


def test_supersaturation_slope_is_namelist_controlled_with_isogsm_default():
    decl = re.search(
        r"real\(r8\)\s*::\s*DEF_TRACER_ICE_SUPERSAT_SLOPE\s*=\s*([0-9._r8]+)",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_ICE_SUPERSAT_SLOPE declaration not found"
    assert float(decl.group(1).replace("_r8", "")) == 0.003
    assert "DEF_TRACER_ICE_SUPERSAT_SLOPE" in FRAC, "closure must read the namelist knob"


def test_supersaturation_slope_is_validated_and_broadcast():
    """A negative slope would invert the kinetic term; and a knob that is
    not broadcast leaves MPI workers silently on the default value."""
    assert "DEF_TRACER_ICE_SUPERSAT_SLOPE < 0._r8" in NAMELIST
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_ICE_SUPERSAT_SLOPE\s*,\s*1\s*,\s*mpi_double_precision",
        NAMELIST,
    ), "DEF_TRACER_ICE_SUPERSAT_SLOPE is not broadcast to MPI workers"
    assert re.search(
        r"DEF_TRACER_ICE_SUPERSAT_SLOPE,\s*&", NAMELIST
    ), "DEF_TRACER_ICE_SUPERSAT_SLOPE is not in the namelist group"


# ---------------------------------------------------------------------------
# Craig-Gordon: the humidity cap and the sign of the ratio.
# Numerical properties live in test_tracer_isotope_frac_runtime.py.
# ---------------------------------------------------------------------------


def test_one_sided_craig_gordon_floor_is_gone():
    """The old 0.75*R_eq floor and the trailing max(...,0) were physical
    restrictions disguised as guards: together they made net heavy-isotope
    uptake unrepresentable.  Neither may come back."""
    assert "craig_gordon_min_net_ratio_frac" not in FRAC
    body = FRAC.split("FUNCTION tracer_craig_gordon_evap_ratio", 1)[1].split(
        "END FUNCTION tracer_craig_gordon_evap_ratio", 1
    )[0]
    assert "max(tracer_craig_gordon_evap_ratio, 0._r8)" not in body, (
        "a non-negative clamp on the Craig-Gordon ratio is back; it silently "
        "discards net heavy-isotope uptake at high humidity"
    )


def test_craig_gordon_bound_is_symmetric():
    """The remaining magnitude bound must clamp both signs, otherwise it is
    the old one-sided floor under a new name."""
    core = FRAC.split("FUNCTION tracer_craig_gordon_ratio_core", 1)[1].split(
        "END FUNCTION tracer_craig_gordon_ratio_core", 1
    )[0]
    assert "craig_gordon_max_ratio_amplification" in core
    assert re.search(
        r"min\(max\(tracer_craig_gordon_ratio_core,\s*-bound\),\s*bound\)", core
    ), "the Craig-Gordon bound is not symmetric about zero"


def test_humidity_cap_is_namelist_controlled_near_saturation():
    decl = re.search(
        r"real\(r8\)\s*::\s*DEF_TRACER_CG_RELHUM_MAX\s*=\s*([0-9._r8]+)",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_CG_RELHUM_MAX declaration not found"
    value = float(decl.group(1).replace("_r8", ""))
    assert value == 0.99, f"expected a near-saturation cap, got {value}"
    assert "DEF_TRACER_CG_RELHUM_MAX" in FRAC, "the closure must read the namelist cap"
    assert "0.95_r8" not in FRAC.split("FUNCTION tracer_craig_gordon_ratio_core", 1)[1].split(
        "END FUNCTION tracer_craig_gordon_ratio_core", 1
    )[0], "a hard-coded 0.95 cap survives inside the Craig-Gordon core"


def test_humidity_cap_is_validated_and_broadcast():
    assert "DEF_TRACER_CG_RELHUM_MAX <= 0._r8" in NAMELIST
    assert "DEF_TRACER_CG_RELHUM_MAX >= 1._r8" in NAMELIST
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_CG_RELHUM_MAX\s*,\s*1\s*,\s*mpi_double_precision",
        NAMELIST,
    ), "DEF_TRACER_CG_RELHUM_MAX is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_CG_RELHUM_MAX,\s*&", NAMELIST)


# ---------------------------------------------------------------------------
# Open water bodies get their own kinetic law (Merlivat & Jouzel 1979).
# Numerical properties live in test_tracer_isotope_frac_runtime.py.
# ---------------------------------------------------------------------------


def test_species_mj79_scaling_is_registered_not_hard_coded_in_the_closure():
    """The 0.88 HDO factor is a published species property, so it belongs in
    the per-species registration next to the diffusivity ratio -- not as a
    branch on tracer name inside the shared closure."""
    assert "mj79_relative_factor" in REGISTRY
    assert "isotope_mj79_relative_factor" in REGISTRY
    assert re.search(r"mj79_relative_factor\s*=\s*1\.0_r8", O18), (
        "18O must register the reference MJ79 factor 1.0"
    )
    assert re.search(r"mj79_relative_factor\s*=\s*0\.88_r8", HDO), (
        "HDO must register the published MJ79 factor 0.88"
    )
    assert "0.88_r8" not in FRAC, (
        "the HDO scaling must come from the registry, not be a literal in "
        "MOD_Tracer_Frac (mentioning it in a comment is fine)"
    )


def test_water_body_evaporation_uses_mj79_and_not_the_soil_exponent():
    """tracer_waterbody_patch is the lake / water-body path.  Its LIQUID
    kinetic factor must come from the open-water law; the ICE factor stays
    on the stagnant-diffusion exponent because sublimation is not an
    open-water process."""
    body = SPECIAL.split("SUBROUTINE tracer_waterbody_patch", 1)[1].split(
        "END SUBROUTINE tracer_waterbody_patch", 1
    )[0]
    assert "tracer_alpha_kinetic_open_water(itrc" in body, (
        "water-body liquid evaporation still uses the soil/generic exponent"
    )
    assert "alpha_k_liq = tracer_alpha_kinetic_craig_gordon(itrc, .false.)" not in body
    assert "alpha_k_ice = tracer_alpha_kinetic_craig_gordon(itrc, .true.)" in body, (
        "sublimation from a frozen water body must keep the ice exponent"
    )


def test_glacier_patch_keeps_the_exponent_law():
    """A glacier surface is not open water; only the water-body path changes."""
    body = SPECIAL.split("SUBROUTINE tracer_glacier_patch", 1)[1].split(
        "END SUBROUTINE tracer_glacier_patch", 1
    )[0]
    assert "alpha_k_liq = tracer_alpha_kinetic_craig_gordon(itrc, .false.)" in body
    assert "tracer_alpha_kinetic_open_water" not in body


def test_wind_forcing_is_actually_plumbed_to_the_water_body_patch():
    """MJ79 is wind-dependent, so a call site that does not pass wind would
    silently pin every lake to the calm-regime value."""
    call = MAIN.split("CALL tracer_waterbody_patch", 1)[1].split(")", 1)[0] + ")"
    assert "forc_us" in call and "forc_vs" in call, (
        f"forc_us/forc_vs are not passed to tracer_waterbody_patch: {call!r}"
    )
    signature = SPECIAL.split("SUBROUTINE tracer_waterbody_patch", 1)[1][:1200]
    assert "forc_us" in signature and "forc_vs" in signature


def test_open_water_scheme_is_selectable_with_mj79_default():
    decl = re.search(
        r"character\(len=16\)\s*::\s*DEF_TRACER_OPEN_WATER_KINETIC\s*=\s*'([A-Za-z0-9]+)'",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_OPEN_WATER_KINETIC declaration not found"
    assert decl.group(1) == "MJ79"
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_OPEN_WATER_KINETIC\s*,\s*16\s*,\s*mpi_character",
        NAMELIST,
    ), "DEF_TRACER_OPEN_WATER_KINETIC is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_OPEN_WATER_KINETIC,\s*&", NAMELIST)
    # An unknown value must fail loudly rather than fall through to a default.
    assert "is invalid; use MJ79 or EXPONENT" in NAMELIST


# ---------------------------------------------------------------------------
# Sublimation exchanges only through a surface skin.
# Numerical properties live in test_tracer_isotope_frac_runtime.py.
# ---------------------------------------------------------------------------

EVAPO = text("main/TRACER/MOD_Tracer_Evapo.F90")
SOIL = text("main/TRACER/MOD_Tracer_SoilWater.F90")
EVAPLIMIT = text("main/TRACER/MOD_Tracer_EvapLimit.F90")


def test_sublimation_skin_is_namelist_controlled_validated_and_broadcast():
    decl = re.search(
        r"real\(r8\)\s*::\s*DEF_TRACER_SUBL_SKIN_MM\s*=\s*([0-9._r8]+)", NAMELIST
    )
    assert decl is not None, "DEF_TRACER_SUBL_SKIN_MM declaration not found"
    assert float(decl.group(1).replace("_r8", "")) == 5.0
    assert "DEF_TRACER_SUBL_SKIN_MM < 0._r8" in NAMELIST
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_SUBL_SKIN_MM\s*,\s*1\s*,\s*mpi_double_precision",
        NAMELIST,
    ), "DEF_TRACER_SUBL_SKIN_MM is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_SUBL_SKIN_MM,\s*&", NAMELIST)


def test_only_the_ice_branch_takes_the_skin_limiter():
    """Soil evaporation is a genuinely mixed-pool process; the skin argument
    is specifically about a snow LAYER not being well mixed.  Applying it to
    the liquid branch as well would silently damp soil evaporative
    enrichment for the wrong reason."""
    wrapper = EVAPO.split("FUNCTION evaporative_tracer_loss", 1)[1].split(
        "END FUNCTION evaporative_tracer_loss", 1
    )[0]
    # The wrapper's outer IF handles nonvolatile solutes; the phase split is
    # the inner "IF (from_ice) THEN".
    phase_split = wrapper.split("IF (from_ice) THEN", 1)
    assert len(phase_split) == 2, "no from_ice branch in evaporative_tracer_loss"
    ice_branch, liquid_branch = phase_split[1].split("ELSE", 1)
    assert "tracer_skin_limited_tracer_loss" in ice_branch
    assert "DEF_TRACER_SUBL_SKIN_MM" in ice_branch
    assert "tracer_skin_limited_tracer_loss" not in liquid_branch
    assert "tracer_evaporative_tracer_loss" in liquid_branch


def test_snow_layer_sublimation_routes_through_the_skin_limiter():
    """Both tracer_soil_water and tracer_wetland reach sublimation through
    atmospheric_loss_tracer, so both call sites must pass the skin mass."""
    assert SOIL.count("skin_mass = DEF_TRACER_SUBL_SKIN_MM") == 2, (
        "expected both atmospheric_loss_tracer wrappers to pass the skin mass"
    )
    guard = EVAPLIMIT.split("FUNCTION tracer_atmospheric_tracer_loss", 1)[1].split(
        "END FUNCTION tracer_atmospheric_tracer_loss", 1
    )[0]
    assert "from_ice .and. present(skin_mass)" in guard, (
        "the skin limiter must apply only to the ice branch, and only when a "
        "caller actually supplies a skin mass"
    )


# ---------------------------------------------------------------------------
# Soil evaporation: resistance-weighted kinetic law.
# Numerical properties live in test_tracer_isotope_frac_runtime.py.
# ---------------------------------------------------------------------------


def test_resistance_weighted_soil_law_is_actually_reachable():
    """tracer_alpha_kinetic_soil used to be dead code while soil evaporation
    ran on the fixed exponent.  It must now be called, and it must be fed
    CoLM's own resistances."""
    assert "tracer_alpha_kinetic_soil(itrc, ra_frac, rss_frac)" in SOIL
    assert "rss_frac = rss" in MAIN, "rss is not plumbed from CoLMMAIN"
    assert MAIN.count("rss_frac = rss") == 2, (
        "both tracer_soil_water call sites (snl<0 and snl==0) must pass rss"
    )
    assert re.search(
        r"real\(r8\),\s*intent\(in\),\s*optional\s*::\s*rss_frac", SOIL
    ), "rss_frac must be an optional dummy so third-party callers keep working"


def test_soil_law_applies_only_to_bare_soil_liquid_evaporation():
    """rss is a soil pore resistance: applying it to snow-layer or
    surface-water evaporation would be physically meaningless."""
    guard = SOIL.split("IF (kinetic_on_soil_surface", 1)[1].split("ENDIF", 1)[0]
    assert ".not. from_ice" in guard
    assert "DEF_TRACER_SOIL_KINETIC" in guard
    assert "present(ra_frac)" in guard and "present(rss_frac)" in guard
    # The flag must be armed and disarmed around each soil-surface call.
    assert SOIL.count("kinetic_on_soil_surface = .true.") == 2
    assert SOIL.count("kinetic_on_soil_surface = .false.") == 3, (
        "expected one reset per armed call plus the per-call initialisation"
    )


def test_soil_surface_flag_is_not_implicitly_saved():
    """CoLM builds with -fopenmp.  An initialised local variable acquires an
    implicit SAVE and would be SHARED between threads, so the flag must be
    declared bare and assigned at runtime."""
    assert re.search(
        r"logical\s+::\s+kinetic_on_soil_surface\s*(?:!.*)?$",
        SOIL,
        re.MULTILINE,
    ), "kinetic_on_soil_surface must be declared without an initialiser"
    assert "logical  :: kinetic_on_soil_surface = " not in SOIL


def test_soil_kinetic_scheme_is_selectable_with_resistance_default():
    decl = re.search(
        r"character\(len=16\)\s*::\s*DEF_TRACER_SOIL_KINETIC\s*=\s*'([A-Za-z]+)'",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_SOIL_KINETIC declaration not found"
    assert decl.group(1) == "RESISTANCE"
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_SOIL_KINETIC\s*,\s*16\s*,\s*mpi_character",
        NAMELIST,
    ), "DEF_TRACER_SOIL_KINETIC is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_SOIL_KINETIC,\s*&", NAMELIST)
    assert "is invalid; use RESISTANCE or EXPONENT" in NAMELIST


def test_soil_segments_use_their_own_transport_exponents():
    """The aerodynamic segment is turbulent (n=2/3) and the pore segment is
    stagnant diffusion (n=1).  Using one exponent for both collapses the
    wet/dry contrast the change exists to represent."""
    core = FRAC.split("FUNCTION tracer_soil_kinetic_alpha_core", 1)[1].split(
        "END FUNCTION tracer_soil_kinetic_alpha_core", 1
    )[0]
    assert "craig_gordon_kinetic_exponent_liquid" in core
    assert "craig_gordon_kinetic_exponent_ice" in core


# ---------------------------------------------------------------------------
# Liquid-phase soil diffusion.
# Numerical properties live in test_tracer_isotope_frac_runtime.py.
# ---------------------------------------------------------------------------


def test_soil_diffusion_is_switchable_and_broadcast():
    assert re.search(
        r"logical\s+::\s+DEF_TRACER_SOIL_DIFFUSION\s*=\s*\.true\.", NAMELIST
    ), "soil isotope diffusion should ship enabled"
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_SOIL_DIFFUSION\s*,\s*1\s*,\s*mpi_logical",
        NAMELIST,
    ), "DEF_TRACER_SOIL_DIFFUSION is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_SOIL_DIFFUSION,\s*&", NAMELIST)


def test_diffusion_needs_geometry_and_is_gated_on_it():
    """dz and porosity are optional dummies, so the sweep must check they are
    present -- otherwise a third-party caller crashes on unassociated args."""
    assert re.search(
        r"real\(r8\),\s*intent\(in\),\s*optional\s*::\s*dz_soi_frac\(1:nl_soil\)", SOIL
    )
    assert re.search(
        r"real\(r8\),\s*intent\(in\),\s*optional\s*::\s*porsl_frac\(1:nl_soil\)", SOIL
    )
    assert "DEF_TRACER_SOIL_DIFFUSION .or. DEF_TRACER_SOIL_VAPOR_DIFFUSION" in SOIL
    assert "present(dz_soi_frac)" in SOIL
    assert "present(porsl_frac)" in SOIL
    assert MAIN.count("dz_soi_frac = dz_soisno(1:nl_soil)") == 2
    assert MAIN.count("porsl_frac = porsl(1:nl_soil)") == 2


def test_diffusion_is_an_internal_exchange_with_no_flux_accumulator():
    """Diffusion moves tracer between layers only.  If it ever booked into an
    a_trc_* accumulator it would be double-counted as a boundary flux and the
    balance check would break."""
    sweep = SOIL.split("liquid-phase molecular diffusion between soil layers", 1)[1]
    sweep = sweep.split("firn vapour diffusion between snow layers", 1)[0]
    assert "a_trc_" not in sweep, "diffusion must not touch flux accumulators"
    assert "trc_numerical_residual_step" not in sweep, (
        "diffusion is exactly conserving; it is not a numerical residual"
    )
    # Equal and opposite: the same signed amount leaves j and enters j+1.
    assert "trc_wliq_soisno(itrc, j,   ipatch) - diff_transfer" in sweep
    assert "trc_wliq_soisno(itrc, j+1, ipatch) + diff_transfer" in sweep


def test_diffusion_runs_on_the_post_advection_water_state():
    """Using the pre-advection water amounts would make the ratios
    inconsistent with the tracer inventory the sweep is reading."""
    sweep = SOIL.split("liquid-phase molecular diffusion between soil layers", 1)[1]
    sweep = sweep.split("firn vapour diffusion between snow layers", 1)[0]
    assert "water_shadow(j)" in sweep and "water_shadow(j+1)" in sweep
    assert "wliq_soisno_bef" not in sweep


def _snapshot_diffuse(inv, faces):
    inv = list(inv)
    out = [0.0] * len(inv)
    for i, f in enumerate(faces):
        if f > 0.0:
            out[i] += f
        else:
            out[i + 1] -= f
    scale = [min(1.0, stock / out_i) if out_i > 0.0 else 1.0 for stock, out_i in zip(inv, out)]
    for i, f in enumerate(faces):
        f = f * (scale[i] if f > 0.0 else scale[i + 1])
        inv[i] -= f
        inv[i + 1] += f
    return inv


def test_diffusion_snapshot_limiter_is_conservative_nonnegative_and_symmetric():
    # Middle layer tries to donate 0.6 to each neighbour from a 0.2 inventory.
    # Snapshot limiting scales the two outgoing faces together, not one at a time.
    result = _snapshot_diffuse([1.0, 0.2, 1.0], [-0.6, 0.6])
    assert sum(result) == pytest.approx(2.2)
    assert min(result) >= -1e-15
    assert result == pytest.approx(list(reversed(result)))
    assert result[1] == pytest.approx(0.0)


def test_soil_and_firn_diffusion_compute_faces_before_applying():
    for start, stop in (
        ("liquid-phase molecular diffusion between soil layers", "firn vapour diffusion between snow layers"),
        ("firn vapour diffusion between snow layers", "3. Groundwater"),
    ):
        sweep = SOIL.split(start, 1)[1].split(stop, 1)[0]
        assert "diff_face" in sweep and "diff_out" in sweep and "diff_scale" in sweep
        assert sweep.index("diff_face(j) = tracer_soil_diffusive_transfer") < sweep.rindex(
            "trc_w"
        )


# ---------------------------------------------------------------------------
# Wet-leaf two-way equilibrium exchange with ambient vapour.
# Numerical properties live in test_tracer_isotope_frac_runtime.py.
# ---------------------------------------------------------------------------

VARS = text("main/TRACER/MOD_Tracer_Vars.F90")
CONSERVATION = text("main/TRACER/MOD_Tracer_Conservation.F90")


def test_canopy_equilibration_ships_disabled_and_is_validated():
    """Unlike the other terms added here, the equilibration degree has no
    reference land-surface implementation to calibrate against, so enabling
    it must be an explicit choice."""
    decl = re.search(
        r"real\(r8\)\s*::\s*DEF_TRACER_CANOPY_EQUILIBRATION\s*=\s*([0-9._r8]+)",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_CANOPY_EQUILIBRATION declaration not found"
    assert float(decl.group(1).replace("_r8", "")) == 0.0
    assert "DEF_TRACER_CANOPY_EQUILIBRATION < 0._r8" in NAMELIST
    assert "DEF_TRACER_CANOPY_EQUILIBRATION > 1._r8" in NAMELIST
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_CANOPY_EQUILIBRATION\s*,\s*1\s*,\s*mpi_double_precision",
        NAMELIST,
    ), "DEF_TRACER_CANOPY_EQUILIBRATION is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_CANOPY_EQUILIBRATION,\s*&", NAMELIST)


def test_vapour_exchange_has_its_own_accumulator_with_a_full_lifecycle():
    """A boundary flux that is allocated but never zeroed or deallocated
    leaks across runs and LULCC remaps."""
    assert "a_trc_vapor_exchange" in VARS
    assert "allocate(a_trc_vapor_exchange" in VARS
    assert "deallocate(a_trc_vapor_exchange)" in VARS
    assert "a_trc_vapor_exchange = 0._r8" in VARS
    assert "a_trc_vapor_exchange(itrc, :) = 0._r8" in VARS, (
        "per-tracer reset (LULCC / tracer re-init path) is missing"
    )
    assert "PUBLIC :: a_trc_vapor_exchange" in VARS


def test_vapour_exchange_enters_the_balance_as_a_signed_input():
    """Zero water flux plus non-zero tracer flux means the balance check must
    see it explicitly, or every exchanging step is reported as a leak."""
    assert "snap_vapor_exchange" in CONSERVATION
    assert "snap_vapor_exchange(itrc, ipatch) = a_trc_vapor_exchange(itrc, ipatch)" in CONSERVATION
    assert "step_input_check = step_input_check + step_vapor_exchange" in CONSERVATION


def test_fixed_signature_fluxes_are_checked_without_hiding_actual_accounting():
    """Mass balance must use booked fluxes; the theoretical fixed signature is
    a separate contract. Otherwise replacing input/evaporation by water*R_init
    can make an accounting bug disappear from the conservation residual."""
    accounting = CONSERVATION.split("step_input_check = step_input\n", 1)[1].split(
        "! Conservation:", 1
    )[0]
    assert "step_output_check = step_output" in accounting
    assert "water_input_in * R_init" not in accounting
    assert "water_evap_in  * R_init" not in accounting
    assert "step_input_check = step_input_check + step_vapor_exchange" in accounting

    signature = CONSERVATION.split("signature_error = 0._r8", 1)[1].split(
        "IF (abs(check_err) > balance_tol)", 1
    )[0]
    for term in (
        "abs(in_minus_water_R)",
        "abs(evap_minus_water_R)",
        "abs(rnof_minus_water_R)",
    ):
        assert term in signature
    assert "fixed_signature_step .and. signature_error > signature_tol" in signature
    assert "TRC_SIG step report" in CONSERVATION


def test_canopy_equilibration_applies_to_liquid_only_and_books_the_flux():
    block = EVAPO.split("Wet-leaf two-way equilibrium exchange", 1)[1].split(
        "Soil+snow layers", 1
    )[0]
    assert "tracer_equilibration_exchange(" in block
    assert "tracer_alpha_liq_vap(itrc, canopy_temp())" in block, (
        "the canopy pool must relax toward the LIQUID-vapour equilibrium"
    )
    assert "ldew_rain" in block and "ldew_snow" not in block, (
        "exchange with a frozen canopy pool is far slower and is out of scope"
    )
    assert "a_trc_vapor_exchange(itrc, ipatch) + equil_gain" in block, (
        "the exchange must be booked, or the balance check reports a leak"
    )
    assert "a_trc_precip(itrc, ipatch) =" not in block, (
        "this flux carries no water and must not be written into the "
        "precipitation accumulator, which is cross-checked against "
        "water_input * R_init (mentioning it in a comment is fine)"
    )
    assert "tracer_is_nonvolatile_solute(itrc)" in block, (
        "nonvolatile solutes do not exchange through the vapour phase"
    )


def test_fortran_smoke_budget_is_not_tight_enough_to_delete_coverage():
    """Every compile-and-run numerical test is gated by this timeout, and a
    skip reads as a pass.  At 5 s and load average ~5 this silently removed
    65-127 tests from a run that still reported success -- the only symptom
    was the wall time.  Keep the budget generous."""
    support = text("tests/fortran_test_support.py")
    budget = re.search(r"^SMOKE_TIMEOUT\s*=\s*(\d+)", support, re.MULTILINE)
    assert budget is not None, "SMOKE_TIMEOUT not found"
    assert int(budget.group(1)) >= 30, (
        f"SMOKE_TIMEOUT is {budget.group(1)}s; too tight a budget turns load "
        "spikes into silent coverage loss"
    )
    # A real timeout must say it was a timeout, not blame the environment.
    assert "the host is likely overloaded" in support


def test_snowmelt_exchange_targets_ice_water_equilibrium_with_the_reciprocal():
    """The percolating water must relax toward R_ice / alpha_ice_liq, so the
    closure has to be handed 1/alpha_ice_liq.  Passing alpha_ice_liq itself
    would enrich the meltwater instead of depleting it -- the opposite of the
    observed early-melt signal."""
    block = SOIL.split("meltwater <-> layer ice isotopic exchange", 1)[1]
    block = block.split("IF (qout_snow > trc_tiny", 1)[0]
    assert "tracer_equilibration_exchange(" in block
    assert re.search(
        r"1\._r8\s*/\s*max\(tracer_alpha_ice_liq\(itrc,\s*layer_temp\(j\)\)", block
    ), "the equilibrium target is not the reciprocal of alpha_ice_liq"
    # Exchange is with the layer's own ice, and it is internal.
    assert "trc_wice_soisno(itrc, j, ipatch) - melt_exchange" in block
    assert "trc_before_flow + melt_exchange" in block
    assert "a_trc_" not in block, "internal exchange must not book a boundary flux"
    assert "trc_numerical_residual_step" not in block
    assert "tracer_is_nonvolatile_solute(itrc)" in block


def test_snowmelt_equilibration_is_switchable_validated_and_broadcast():
    decl = re.search(
        r"real\(r8\)\s*::\s*DEF_TRACER_SNOWMELT_EQUILIBRATION\s*=\s*([0-9._r8]+)",
        NAMELIST,
    )
    assert decl is not None, "DEF_TRACER_SNOWMELT_EQUILIBRATION not declared"
    assert float(decl.group(1).replace("_r8", "")) == 0.0
    assert "DEF_TRACER_SNOWMELT_EQUILIBRATION < 0._r8" in NAMELIST
    assert "DEF_TRACER_SNOWMELT_EQUILIBRATION > 1._r8" in NAMELIST
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_SNOWMELT_EQUILIBRATION\s*,\s*1\s*,\s*mpi_double_precision",
        NAMELIST,
    ), "DEF_TRACER_SNOWMELT_EQUILIBRATION is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_SNOWMELT_EQUILIBRATION,\s*&", NAMELIST)


def test_firn_diffusion_acts_on_the_ice_carrier_not_the_liquid_film():
    """In a snowpack the isotopes sit in the ice; the liquid film is a minor,
    transient carrier.  Diffusing trc_wliq there would move almost nothing and
    leave the actual signal untouched."""
    sweep = SOIL.split("firn vapour diffusion between snow layers", 1)[1]
    sweep = sweep.split("3. Groundwater", 1)[0]
    assert "trc_wice_soisno(itrc, j,   ipatch) - diff_transfer" in sweep
    assert "trc_wice_soisno(itrc, j+1, ipatch) + diff_transfer" in sweep
    assert "trc_wliq_soisno" not in sweep, (
        "the firn sweep must move the ice inventory, not the liquid film"
    )
    assert "tracer_alpha_ice_vap(itrc, layer_temp(j))" in sweep, (
        "pore vapour over a snowpack equilibrates with ICE"
    )
    # Snow layers only, and only when they exist.
    assert "snl < 0" in sweep and "DO j = lb, -1" in sweep
    assert "present(dz_sno_frac)" in sweep
    # Same internal-exchange discipline as the soil sweep.
    assert "a_trc_" not in sweep
    assert "trc_numerical_residual_step" not in sweep


def test_snow_layer_thickness_is_plumbed_only_where_snow_exists():
    """dz_sno_frac is dimensioned snl+1:0, so passing it from the snl==0 call
    site would hand over a zero-size array for a sweep that cannot run."""
    assert MAIN.count("dz_sno_frac = dz_soisno(snl+1:0)") == 2, (
        "expected the ordinary and wetland snl<0 call sites to pass snow thickness"
    )
    assert re.search(
        r"real\(r8\),\s*intent\(in\),\s*optional\s*::\s*dz_sno_frac\(snl\+1:0\)", SOIL
    )


def test_soil_vapour_diffusion_is_added_to_the_liquid_term_not_replacing_it():
    """Liquid and vapour pathways have opposite moisture dependences, so the
    column only stays connected at all wetnesses if BOTH are present."""
    sweep = SOIL.split("liquid-phase molecular diffusion between soil layers", 1)[1]
    sweep = sweep.split("firn vapour diffusion between snow layers", 1)[0]
    assert "d_eff_up = d_eff_up + tracer_soil_vapor_equivalent_diffusivity(" in sweep
    assert "d_eff_dn = d_eff_dn + tracer_soil_vapor_equivalent_diffusivity(" in sweep
    assert "DEF_TRACER_SOIL_VAPOR_DIFFUSION .and. present(forc_psrf_frac)" in sweep, (
        "the vapour term needs surface pressure and must be gated on it"
    )
    # Folding it into the diffusivity is what lets the verified flux form and
    # overshoot limiter stay untouched.
    assert sweep.count("tracer_soil_diffusive_transfer(") == 1


def test_soil_vapour_diffusion_is_switchable_and_broadcast():
    assert re.search(
        r"logical\s+::\s+DEF_TRACER_SOIL_VAPOR_DIFFUSION\s*=\s*\.true\.", NAMELIST
    )
    assert re.search(
        r"mpi_bcast\s*\(\s*DEF_TRACER_SOIL_VAPOR_DIFFUSION\s*,\s*1\s*,\s*mpi_logical",
        NAMELIST,
    ), "DEF_TRACER_SOIL_VAPOR_DIFFUSION is not broadcast to MPI workers"
    assert re.search(r"DEF_TRACER_SOIL_VAPOR_DIFFUSION,\s*&", NAMELIST)


def test_differing_ra_exponents_for_leaf_and_soil_are_documented():
    """The leaf form weights ra by 1 and the soil form by D**(2/3).  That is a
    physics choice (CoLM's soil ra absorbs the quasi-laminar sublayer, so
    weighting it by 1 would zero out wet-soil kinetic fractionation), and it
    must stay explained where someone would otherwise 'fix' it."""
    note = FRAC.split("NOTE on the `ra` exponent", 1)
    assert len(note) == 2, "the leaf/soil ra exponent difference is undocumented"
    note = note[1][:1400]
    assert "Farquhar" in note
    assert "quasi-laminar sublayer" in note
    assert "zero kinetic fractionation" in note or "alpha_k = 1" in note


def test_vapour_exchange_flux_is_diagnosable():
    """It is already a term in the balance equation, so leaving it out of
    history makes it a black box the moment anyone enables it.

    It must NOT go through the delta writer: that divides tracer mass by a
    paired water flux, and this flux has none (zero net water by
    construction), so a delta is undefined for it.
    """
    hist = text("main/TRACER/MOD_Tracer_Hist.F90")
    assert "a_trc_vapor_exchange" in hist, (
        "the vapour-exchange flux has no history output"
    )
    assert "f_trc_vapor_exchange_" in hist
    block = hist.split("f_trc_vapor_exchange_", 1)[1][:900]
    assert "write_history_variable_2d" in block, (
        "a flux with no paired water flux must be written as a mass"
    )
    assert "write_history_tracer_delta_2d" not in block, (
        "the delta writer needs a paired water flux; this flux has none"
    )


def test_leaf_nss_shares_the_craig_gordon_humidity_cap():
    """The NSS end-member is delta_es = ... + h*(delta_v - eps_k - delta_x),
    i.e. LINEAR in h with no 1/(1-h) singularity.  Its former hard-coded 0.95
    was therefore not a numerical guard either -- it truncated exactly the
    same real depletion the Craig-Gordon cap did, and left the leaf and
    surface closures disagreeing about how humid the air is allowed to get.
    One namelist knob must own the cap.
    """
    nss = FRAC.split("SUBROUTINE tracer_transpiration_nss_ratio", 1)[1].split(
        "END SUBROUTINE tracer_transpiration_nss_ratio", 1
    )[0]
    assert "h = min(max(relhum, 0._r8), 0.95_r8)" not in nss, (
        "the leaf NSS humidity cap is still hard-coded at 0.95"
    )
    assert "DEF_TRACER_CG_RELHUM_MAX" in nss
    # The unrelated storage-tendency limiter also uses 0.95 and must survive.
    assert "max_storage_tendency_fraction = 0.95_r8" in FRAC


def test_isogsm_driven_amazon_case_uses_the_merlivat_pair():
    """run/Amazon_iso_0p5_2020.nml is driven by IsoGSM forcing, so it must
    not override the default back to Cappa."""
    assert "IsoGSM" in AMAZON, "expected the IsoGSM-forced Amazon configuration"
    setting = re.search(
        r"DEF_TRACER_KINETIC_SCHEME\s*=\s*'([A-Za-z0-9]+)'",
        AMAZON,
    )
    assert setting is not None
    assert setting.group(1) == "MERLIVAT1978"


# ---------------------------------------------------------------------------
# Cross-patch consistency: special surfaces and wetlands must use the same
# finite-pool and phase-specific process implementations as ordinary land.
# ---------------------------------------------------------------------------


def _subroutine(source: str, name: str) -> str:
    match = re.search(
        rf"SUBROUTINE\s+{name}\b(?P<body>.*?)END\s+SUBROUTINE\s+{name}\b",
        source,
        flags=re.IGNORECASE | re.DOTALL,
    )
    assert match is not None, f"{name} not found"
    return match.group("body")


def _compact(source: str) -> str:
    return " ".join(source.lower().replace("&", " ").split())


def test_special_patches_share_the_finite_pool_atmospheric_loss():
    assert "USE MOD_Tracer_EvapLimit, only: tracer_atmospheric_tracer_loss" in SPECIAL
    for name in ("tracer_glacier_patch", "tracer_waterbody_patch"):
        body = _subroutine(SPECIAL, name)
        assert body.count("tracer_atmospheric_tracer_loss") == 2
        assert not re.search(r"R_evap_(?:liq|ice)\s*=\s*min\s*\(", body, re.IGNORECASE)
        assert re.search(
            r"trc_evap_ice\s*=\s*tracer_atmospheric_tracer_loss\(.*?"
            r"skin_mass\s*=\s*DEF_TRACER_SUBL_SKIN_MM",
            body,
            flags=re.DOTALL,
        )


def test_special_patch_ice_loss_sees_the_pool_after_liquid_loss():
    for name in ("tracer_glacier_patch", "tracer_waterbody_patch"):
        body = _subroutine(SPECIAL, name)
        assert "trc_after_liq = trc_available - trc_evap_liq" in body
        assert "water_after_liq = water_before_output - evap_liq_mass" in body
        assert re.search(
            r"trc_evap_ice\s*=\s*tracer_atmospheric_tracer_loss\(\s*"
            r"trc_after_liq,\s*water_after_liq",
            body,
            flags=re.DOTALL,
        )


def test_wetland_liquid_loss_uses_wind_dependent_open_water_kinetics():
    body = _compact(_subroutine(SOIL, "tracer_wetland"))
    assert "forc_us_frac" in body and "forc_vs_frac" in body
    assert (
        "alpha_k = tracer_alpha_kinetic_open_water(itrc, sqrt(max("
        "forc_us_frac*forc_us_frac + forc_vs_frac*forc_vs_frac, 0._r8)))"
    ) in body
    assert "alpha_k = tracer_alpha_kinetic_craig_gordon(itrc, .true.)" in body
    assert "alpha_k = tracer_alpha_kinetic_craig_gordon(itrc, from_ice)" not in body


def test_both_wetland_call_sites_pass_wind_components():
    calls = MAIN.split("CALL tracer_wetland")[1:]
    assert len(calls) == 2
    for call in calls:
        normalized = _compact(call.split("forc_psrf_frac = forc_psrf", 1)[0])
        assert "forc_us, forc_vs, waterstorage_trc_ground" in normalized


def test_wetland_snow_percolation_applies_snowmelt_equilibration():
    body = _subroutine(SOIL, "tracer_wetland")
    percolation = body.split("! Snow percolation → trc_gwat_snow_local", 1)[1].split(
        "trc_gwat_snow_local = trc_qin_snow", 1
    )[0]
    compact = _compact(percolation)
    assert "def_tracer_snowmelt_equilibration > 0._r8" in compact
    assert "tracer_equilibration_exchange(trc_before_flow" in compact
    assert (
        "1._r8 / max(tracer_alpha_ice_liq(itrc, layer_temp(j)), trc_tiny)"
        in compact
    )
    assert (
        "trc_wice_soisno(itrc, j, ipatch) = max("
        "trc_wice_soisno(itrc, j, ipatch) - melt_exchange, 0._r8)"
        in compact
    )


def test_wetland_transpiration_uses_leaf_nss_storage_when_plumbed():
    body = _subroutine(SOIL, "tracer_wetland")
    loss = body.split("transp_source_tracer = q_etr_out * pool_ratio_loss", 1)[1].split(
        "trc_loss = trc_evap_loss + trc_subl_loss + trc_etr_loss", 1
    )[0]
    assert "tracer_transpiration_nss_ratio" in loss
    assert "trc_leaf_iso_storage(itrc, ipatch)" in loss
    assert "transp_source_tracer - transp_output_tracer" in loss
    assert "trc_etr_loss = transp_source_tracer" in body
    booking = body.split("TRC_EVAP_KIND_TRANSP", 1)[0]
    assert "transp_output_tracer" in booking


def test_wetland_signed_fluxes_are_booked_before_net_cancellation():
    body = _subroutine(SOIL, "tracer_wetland")
    booking = body.split("Book signed components independently", 1)[1].split(
        "pool_tracer = pool_tracer - trc_loss", 1
    )[0]
    assert "IF (abs(trc_loss)" not in booking
    assert "trc_evap_loss + trc_subl_loss" in booking
    assert "transp_output_tracer" in booking
    assert "a_trc_transp_src(itrc, ipatch) = a_trc_transp_src(itrc, ipatch) + trc_etr_loss" in booking

    # A negative atmospheric exchange can exactly cancel source-water
    # transpiration while NSS still changes the external signature.
    atmospheric, source, output = -1.0, 1.0, 0.8
    net_pool_loss = atmospheric + source
    leaf_storage_change = source - output
    booked_external_loss = atmospheric + output
    assert net_pool_loss == 0.0
    assert leaf_storage_change == pytest.approx(-booked_external_loss)


def test_wetland_leaf_off_releases_leaf_nss_storage_like_ordinary_path():
    body = _subroutine(SOIL, "tracer_wetland")
    setup = body.split("IF (tracer_has_dissolved_limit(itrc))", 1)[0]
    assert "release_leaf_iso_storage(itrc, ipatch, nl_soil" in setup
    assert "wliq_soisno_bef(1:nl_soil), wa_bef" in setup
    assert "lai_frac <= trc_tiny" in setup
    assert "trc_leaf_water_moles(itrc, ipatch) <= trc_tiny" in setup


def test_wetland_snow_path_runs_firn_diffusion_and_gets_snow_thickness():
    body = _subroutine(SOIL, "tracer_wetland")
    assert "dz_sno_frac(snl+1:0)" in body
    firn = body.split("Wetland snowpacks keep the same firn vapour diffusion", 1)[1].split(
        "2) Strip wresi", 1
    )[0]
    assert "tracer_snow_vapor_equivalent_diffusivity" in firn
    assert "trc_wice_soisno(itrc, j,   ipatch) - diff_transfer" in firn
    snow_call = MAIN.split("CALL tracer_wetland", 1)[1].split("ELSE", 1)[0]
    assert "dz_sno_frac = dz_soisno(snl+1:0)" in snow_call


def test_fractionating_vapor_fallback_is_reported():
    log = FORCING.split("SUBROUTINE tracer_forcing_log_ranges", 1)[1].split(
        "END SUBROUTINE tracer_forcing_log_ranges", 1
    )[0]
    assert "tracer_fractionation_active(itrc)" in log
    assert "vmiss(itrc) = vmiss(itrc) + 1" in log
    assert "WARNING vapor " in log
    assert "fallback patches=" in log
