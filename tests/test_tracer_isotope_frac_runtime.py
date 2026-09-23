"""Runtime (behavioural) tests for the isotope fractionation closures.

Why this file exists
--------------------
The fractionation coefficients are pure functions of temperature, so the
mathematical properties they must satisfy -- continuity across the
piecewise branches, reduction to the equilibrium limit when a switch is
off, the correct *sign* of the kinetic effect -- can be checked against
the really shipped formulas rather than against their source text.

A mistyped sign or a swapped equilibrium/effective alpha in the
supersaturation branch leaves the code compiling and every static
assertion green, while silently over-enriching deposited frost at every
sub-freezing timestep.  Static tests cannot see that; these can.

Reference: Jouzel & Merlivat (1984), as implemented in the IsoGSM
reference model at codes/IsoGSM/gsm/src/gsml/CLD1/lrgscl.F:196-224.

Scope and limits
----------------
The functions under test are extracted verbatim from the shipped source
and compiled standalone, so the *formulas and constants are the real
ones*.  What is not covered here is whether the callers pass the right
temperature or select the ice branch correctly -- that is the job of the
static tests in test_tracer_isotope_isogsm_alignment.py.
"""

from pathlib import Path
import math
import re
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
FRAC_SOURCE = ROOT / "main/TRACER/MOD_Tracer_Frac.F90"
O18_SOURCE = ROOT / "main/TRACER/MOD_Tracer_Isotope_O18.F90"
HDO_SOURCE = ROOT / "main/TRACER/MOD_Tracer_Isotope_HDO.F90"
SUBPROCESS_TIMEOUT = 60

# Merlivat (1978) 18O air diffusivity ratio -- the shipped default pair.
DIFF_RATIO_O18 = 1.0285
# Jouzel & Merlivat (1984) supersaturation slope, S = 1 - slope * T[C].
JM84_SLOPE = 0.003


def _extract_function(source: Path, name: str) -> str:
    """Return the verbatim text of a module function, header through END."""
    text = source.read_text(encoding="utf-8")
    start = re.search(
        rf"^\s*(?:real\(r8\)\s+)?FUNCTION\s+{name}\s*\(",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if start is None:
        raise AssertionError(
            f"{name} not found in {source}. If it was renamed, update this test "
            "rather than deleting it -- the property still needs an owner."
        )
    end = re.search(
        rf"^\s*END\s+FUNCTION\s+{name}\b",
        text[start.start():],
        re.MULTILINE | re.IGNORECASE,
    )
    if end is None:
        raise AssertionError(f"unterminated FUNCTION {name} in {source}")
    return text[start.start(): start.start() + end.end()]


def _module_parameter(source: Path, name: str) -> str:
    """Read a module-level real parameter value from the shipped source.

    Extracted rather than hard-coded so that editing the source moves the
    test with it instead of silently decoupling.
    """
    text = source.read_text(encoding="utf-8")
    # Capture to end of line (minus any trailing comment) so that parameters
    # written as expressions -- e.g. `2._r8 / 3._r8` -- survive intact.  A
    # token-only match would silently inject `2._r8` and turn an exponent of
    # 2/3 into 2.
    found = re.search(
        rf"real\(r8\),\s*parameter\s*::\s*{name}\s*=\s*([^!\n]+)",
        text,
    )
    if found is None:
        raise AssertionError(f"module parameter {name} not found in {source}")
    return found.group(1).strip()


@pytest.fixture(scope="module")
def ice_alpha_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_ice_alpha")
    compiler = require_runnable_fortran_compiler(workdir)

    tfrz = _module_parameter(FRAC_SOURCE, "tfrz")
    tcold = _module_parameter(FRAC_SOURCE, "jm84_full_kinetic_temp")

    bodies = "\n".join(
        (
            _extract_function(O18_SOURCE, "o18_alpha_ice_vap"),
            _extract_function(FRAC_SOURCE, "tracer_jm84_effective_alpha"),
            _extract_function(FRAC_SOURCE, "tracer_ice_deposition_alpha"),
        )
    )

    (workdir / "closure.f90").write_text(
        "module tracer_ice_alpha\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        f"  real(r8), parameter :: tfrz = {tfrz}\n"
        f"  real(r8), parameter :: jm84_full_kinetic_temp = {tcold}\n"
        "CONTAINS\n"
        + bodies
        + "\nend module tracer_ice_alpha\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program ice_alpha_driver
  use tracer_ice_alpha
  implicit none
  real(r8) :: temp_k, slope, diff_ratio
  real(r8) :: a_t, a_frz, a_cold
  integer  :: stat

  do
    read(*, *, iostat=stat) temp_k, slope, diff_ratio
    if (stat /= 0) exit
    a_t    = o18_alpha_ice_vap(temp_k)
    a_frz  = o18_alpha_ice_vap(tfrz)
    a_cold = o18_alpha_ice_vap(jm84_full_kinetic_temp)
    write(*, '(2ES24.16)') &
      tracer_ice_deposition_alpha(temp_k, slope, diff_ratio, a_t, a_frz, a_cold), a_t
  end do
end program ice_alpha_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "ice_alpha"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _evaluate(executable: Path, cases) -> list:
    """cases: iterable of (temp_k, slope, diff_ratio) -> [(alpha_eff, alpha_eq)]"""
    cases = list(cases)
    payload = "\n".join(f"{t!r} {s!r} {d!r}" for t, s, d in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    numbers = [float(tok) for tok in ran.stdout.split()]
    assert len(numbers) == 2 * len(cases), ran.stdout + ran.stderr
    return list(zip(numbers[0::2], numbers[1::2]))


def test_zero_slope_reduces_to_pure_equilibrium(ice_alpha_driver):
    """slope=0 is the documented off switch: no supersaturation kinetics."""
    cases = [(t, 0.0, DIFF_RATIO_O18) for t in (233.15, 253.15, 263.15, 273.15)]
    for (alpha_eff, alpha_eq), (temp_k, _, _) in zip(
        _evaluate(ice_alpha_driver, cases), cases
    ):
        assert alpha_eff == pytest.approx(alpha_eq, rel=1e-12), (
            f"slope=0 must disable the kinetic term, but at T={temp_k} K "
            f"alpha_eff={alpha_eff:.9g} != alpha_eq={alpha_eq:.9g}"
        )


def test_at_and_above_freezing_there_is_no_supersaturation_effect(ice_alpha_driver):
    """S = 1 - slope*T[C] is unity at 0 C, so the branch must be inert there."""
    cases = [(t, JM84_SLOPE, DIFF_RATIO_O18) for t in (273.15, 275.0, 280.0)]
    for (alpha_eff, alpha_eq), (temp_k, _, _) in zip(
        _evaluate(ice_alpha_driver, cases), cases
    ):
        assert alpha_eff == pytest.approx(alpha_eq, rel=1e-12), (
            f"at T={temp_k} K (>= 0 C) deposition must stay at equilibrium"
        )


def test_kinetic_term_weakens_enrichment_below_freezing(ice_alpha_driver):
    """Supersaturation reduces alpha below the equilibrium value.

    Getting this sign backwards would *increase* frost enrichment, which is
    the failure mode the whole branch exists to prevent.
    """
    cases = [(t, JM84_SLOPE, DIFF_RATIO_O18) for t in (233.15, 243.15, 253.15, 263.15)]
    for (alpha_eff, alpha_eq), (temp_k, _, _) in zip(
        _evaluate(ice_alpha_driver, cases), cases
    ):
        assert alpha_eff < alpha_eq, (
            f"at T={temp_k} K the effective alpha ({alpha_eff:.9g}) must be BELOW "
            f"the equilibrium alpha ({alpha_eq:.9g}); the kinetic term has the "
            "wrong sign and is over-enriching deposited ice"
        )
        assert alpha_eff > 1.0, (
            f"at T={temp_k} K alpha collapsed to {alpha_eff:.9g}; ice-vapour "
            "fractionation must remain > 1"
        )


@pytest.mark.parametrize("boundary", ("jm84_full_kinetic_temp", "tfrz"))
def test_piecewise_branches_join_continuously(ice_alpha_driver, boundary):
    """The -20 C / 0 C blend must not tear.

    Between -20 C and 0 C the scheme linearly blends the equilibrium alpha
    at 0 C with the *effective* (supersaturated) alpha at -20 C, exactly as
    IsoGSM lrgscl.F:210-224 does.  Blending against the equilibrium value at
    -20 C instead would leave a visible step at the cold end.
    """
    temp = float(
        _module_parameter(FRAC_SOURCE, boundary).replace("_r8", "").replace("d", "e")
    )
    below, above = temp - 1.0e-6, temp + 1.0e-6
    (eff_below, _), (eff_above, _) = _evaluate(
        ice_alpha_driver,
        [(below, JM84_SLOPE, DIFF_RATIO_O18), (above, JM84_SLOPE, DIFF_RATIO_O18)],
    )
    assert eff_below == pytest.approx(eff_above, rel=1e-9), (
        f"discontinuity at {boundary}={temp} K: {eff_below:.12g} -> {eff_above:.12g}. "
        "The blend endpoints no longer match the branch values."
    )


@pytest.fixture(scope="module")
def hdo_alpha_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_hdo_alpha")
    compiler = require_runnable_fortran_compiler(workdir)
    bodies = "\n".join(
        (
            _extract_function(HDO_SOURCE, "hdo_alpha_liq_vap"),
            _extract_function(HDO_SOURCE, "hdo_alpha_ice_vap"),
        )
    )
    (workdir / "closure.f90").write_text(
        "module hdo_alpha\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        "contains\n"
        + bodies
        + "\nend module hdo_alpha\n",
        encoding="utf-8",
    )
    (workdir / "driver.f90").write_text(
        """
program driver
  use hdo_alpha
  implicit none
  real(r8) :: temp_k
  integer :: stat
  do
    read(*, *, iostat=stat) temp_k
    if (stat /= 0) exit
    write(*, '(2ES24.16)') hdo_alpha_liq_vap(temp_k), hdo_alpha_ice_vap(temp_k)
  end do
end program driver
""",
        encoding="utf-8",
    )
    executable = workdir / "hdo_alpha"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


@pytest.mark.parametrize("temp_k", (140.0, 233.15, 253.15, 273.15, 298.15))
def test_hdo_equilibrium_alphas_match_reference_formulas(hdo_alpha_driver, temp_k):
    ran = subprocess.run(
        [str(hdo_alpha_driver)],
        input=f"{temp_k}\n",
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    liquid, ice = (float(value) for value in ran.stdout.split())
    tk = max(temp_k, 150.0)
    assert liquid == pytest.approx(
        math.exp(24844.0 / tk**2 - 76.248 / tk + 0.052612), rel=1e-12
    )
    assert ice == pytest.approx(math.exp(16289.0 / tk**2 - 0.0945), rel=1e-12)


# ---------------------------------------------------------------------------
# Craig-Gordon evaporation ratio.
#
# R_E = (R_s/alpha_eq - h*R_a) / (alpha_k*(1-h))
#
# The 1/(1-h) factor is not a singularity in the *flux*: the water flux
# itself carries a (1-h) factor, so the product E*R_E stays finite.  IsoGSM
# gets this for free by computing the isotope flux from the same bulk
# gradient as the water flux (moninp.F:369-373).  CoLM computes the ratio,
# so the cancellation has to survive the humidity cap instead of being
# truncated by it.
# ---------------------------------------------------------------------------

# Merlivat n=2/3 kinetic factor for 18O over a liquid surface.
ALPHA_K_O18 = 1.0285 ** (2.0 / 3.0)
R_SMOW_RELATIVE = 1.0  # work in R/R_SMOW units; delta = (R-1)*1000


def _delta_to_r(delta: float) -> float:
    return 1.0 + delta / 1000.0


def _majoube_alpha_liq_vap(temp_k: float) -> float:
    """Majoube (1971) liquid-vapour 18O fractionation, independent restatement."""
    return math.exp(1137.0 / temp_k**2 - 0.4156 / temp_k - 0.0020667)


def _craig_gordon_closed_form(rs: float, ra: float, temp_k: float,
                             alpha_k: float, h: float) -> float:
    """R_E = (R_s/alpha_eq - h*R_a) / (alpha_k*(1-h)), computed independently."""
    r_eq = rs / _majoube_alpha_liq_vap(temp_k)
    return (r_eq - h * ra) / (alpha_k * (1.0 - h))


@pytest.fixture(scope="module")
def craig_gordon_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_craig_gordon")
    compiler = require_runnable_fortran_compiler(workdir)

    amplification = _module_parameter(
        FRAC_SOURCE, "craig_gordon_max_ratio_amplification"
    )

    bodies = "\n".join(
        (
            _extract_function(O18_SOURCE, "o18_alpha_liq_vap"),
            _extract_function(FRAC_SOURCE, "tracer_craig_gordon_ratio_core"),
        )
    )

    (workdir / "closure.f90").write_text(
        "module tracer_craig_gordon\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        f"  real(r8), parameter :: craig_gordon_max_ratio_amplification = {amplification}\n"
        "CONTAINS\n"
        + bodies
        + "\nend module tracer_craig_gordon\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program craig_gordon_driver
  use tracer_craig_gordon
  implicit none
  real(r8) :: rs, ra, temp_k, alpha_k, relhum, relhum_max
  real(r8) :: alpha_eq
  integer  :: stat

  do
    read(*, *, iostat=stat) rs, ra, temp_k, alpha_k, relhum, relhum_max
    if (stat /= 0) exit
    alpha_eq = o18_alpha_liq_vap(temp_k)
    write(*, '(2ES24.16)') &
      tracer_craig_gordon_ratio_core(rs, ra, alpha_eq, alpha_k, relhum, relhum_max), &
      rs / alpha_eq
  end do
end program craig_gordon_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "craig_gordon"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _craig_gordon(executable: Path, cases) -> list:
    """cases: (R_s, R_a, T, alpha_k, h, h_max) -> [(R_E, R_equilibrium)]"""
    cases = list(cases)
    payload = "\n".join(" ".join(repr(float(v)) for v in c) for c in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    numbers = [float(tok) for tok in ran.stdout.split()]
    assert len(numbers) == 2 * len(cases), ran.stdout + ran.stderr
    return list(zip(numbers[0::2], numbers[1::2]))


def test_dry_limit_matches_the_closed_form(craig_gordon_driver):
    """At h=0 Craig-Gordon collapses to R_s/(alpha_eq*alpha_k)."""
    (value, r_eq), = _craig_gordon(
        craig_gordon_driver,
        [(_delta_to_r(-5.0), _delta_to_r(-12.0), 298.15, ALPHA_K_O18, 0.0, 0.99)],
    )
    assert value == pytest.approx(r_eq / ALPHA_K_O18, rel=1e-12)


def test_evaporating_vapour_gets_lighter_as_air_gets_wetter(craig_gordon_driver):
    """Monotone in h -- the physical content of the Craig-Gordon closure."""
    humidities = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.98, 0.99]
    results = _craig_gordon(
        craig_gordon_driver,
        [
            (_delta_to_r(-5.0), _delta_to_r(-12.0), 298.15, ALPHA_K_O18, h, 0.99)
            for h in humidities
        ],
    )
    values = [v for v, _ in results]
    for (h_dry, v_dry), (h_wet, v_wet) in zip(
        zip(humidities, values), zip(humidities[1:], values[1:])
    ):
        assert v_wet < v_dry, (
            f"R_E must decrease with humidity, but h={h_dry} gives {v_dry:.9g} "
            f"and h={h_wet} gives {v_wet:.9g}"
        )


def test_humid_regime_is_not_truncated_at_the_old_0p95_cap(craig_gordon_driver):
    """The h=0.95 cap cost ~180 permil of depletion in the humid regime.

    With soil water at -5 permil and vapour at -12 permil at 25 C the true
    h=0.99 evaporate is -250 permil, while capping h at 0.95 reports only
    -74 permil -- a 176 permil error that under-enriches residual soil
    water at every humid timestep.  The cap must now sit at 0.99, so
    evaluating at h=0.99 must NOT return the h=0.95 answer, and must match
    the closed form.
    """
    rs, ra, temp_k = _delta_to_r(-5.0), _delta_to_r(-12.0), 298.15
    (at_95, _), (at_99, _) = _craig_gordon(
        craig_gordon_driver,
        [
            (rs, ra, temp_k, ALPHA_K_O18, 0.95, 0.99),
            (rs, ra, temp_k, ALPHA_K_O18, 0.99, 0.99),
        ],
    )
    delta_95 = (at_95 - 1.0) * 1000.0
    delta_99 = (at_99 - 1.0) * 1000.0
    assert delta_99 < delta_95 - 100.0, (
        f"h=0.99 must be much more depleted than h=0.95, got {delta_99:.1f} vs "
        f"{delta_95:.1f} permil -- the humidity cap is still truncating at 0.95"
    )
    for h, got in ((0.95, at_95), (0.99, at_99)):
        expected = _craig_gordon_closed_form(rs, ra, temp_k, ALPHA_K_O18, h)
        assert got == pytest.approx(expected, rel=1e-12), (
            f"at h={h} the shipped closure gives {got!r} but the Craig-Gordon "
            f"closed form gives {expected!r}"
        )


def test_cap_is_applied_at_the_supplied_maximum_not_a_hidden_constant(craig_gordon_driver):
    """Above relhum_max the result must equal the value AT relhum_max."""
    for cap in (0.95, 0.99):
        (capped, _), (at_cap, _) = _craig_gordon(
            craig_gordon_driver,
            [
                (_delta_to_r(-5.0), _delta_to_r(-12.0), 298.15, ALPHA_K_O18, 0.9999, cap),
                (_delta_to_r(-5.0), _delta_to_r(-12.0), 298.15, ALPHA_K_O18, cap, cap),
            ],
        )
        assert capped == pytest.approx(at_cap, rel=1e-12), (
            f"h above the cap must clamp to relhum_max={cap}"
        )


def test_net_heavy_isotope_uptake_is_representable(craig_gordon_driver):
    """When h*R_a exceeds the equilibrium vapour ratio the net isotope flux
    reverses sign, even though the net water flux is still evaporation.

    The old 0.75*R_eq floor made this impossible to express: it reported a
    positive, near-equilibrium loss where the physics wants a gain.
    """
    (value, r_eq), = _craig_gordon(
        craig_gordon_driver,
        [(_delta_to_r(-5.0), _delta_to_r(0.0), 298.15, ALPHA_K_O18, 0.99, 0.99)],
    )
    assert value < 0.0, (
        f"expected a negative (net uptake) ratio, got {value:.9g}; a hard "
        "non-negative floor is still clamping the Craig-Gordon numerator"
    )
    assert value > -0.75 * r_eq, "sanity: this case should not hit the amplification bound"


def test_extreme_ratio_contrast_is_bounded_not_unbounded(craig_gordon_driver):
    """1/(1-h) amplifies any inconsistency between the host water flux and
    the humidity used here, so a magnitude bound must remain -- but as a
    symmetric numerical guard, not as a one-sided physical floor."""
    amplification = float(
        _module_parameter(FRAC_SOURCE, "craig_gordon_max_ratio_amplification")
        .replace("_r8", "")
    )
    (value, r_eq), = _craig_gordon(
        craig_gordon_driver,
        [(_delta_to_r(-5.0), _delta_to_r(200.0), 298.15, ALPHA_K_O18, 0.99, 0.99)],
    )
    assert value == pytest.approx(-amplification * r_eq, rel=1e-9), (
        f"expected the symmetric bound -{amplification}*R_eq, got {value:.9g}"
    )


# ---------------------------------------------------------------------------
# Finite-pool evaporation limiter.
#
# Craig-Gordon may now return a NEGATIVE ratio (net heavy-isotope uptake
# during net water loss).  The limiter sits between that ratio and the
# prognostic pools, so it has to carry the sign through -- including on its
# sub-stepping path, which exists because enrichment is non-linear as a pool
# dries down.
# ---------------------------------------------------------------------------

EVAPLIMIT_SOURCE = ROOT / "main/TRACER/MOD_Tracer_EvapLimit.F90"


def _module_int_parameter(source: Path, name: str) -> str:
    text = source.read_text(encoding="utf-8")
    found = re.search(rf"integer,\s*parameter\s*::\s*{name}\s*=\s*([^\s!]+)", text)
    if found is None:
        raise AssertionError(f"integer parameter {name} not found in {source}")
    return found.group(1)


@pytest.fixture(scope="module")
def evaplimit_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_evaplimit")
    compiler = require_runnable_fortran_compiler(workdir)

    max_frac = _module_parameter(EVAPLIMIT_SOURCE, "evaplimit_default_max_loss_fraction")
    max_sub = _module_int_parameter(EVAPLIMIT_SOURCE, "evaplimit_default_max_substeps")
    body = "\n".join(
        (
            _extract_function(EVAPLIMIT_SOURCE, "tracer_evaporative_tracer_loss"),
            _extract_function(EVAPLIMIT_SOURCE, "tracer_skin_limited_tracer_loss"),
        )
    )

    (workdir / "closure.f90").write_text(
        "module tracer_evaplimit\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        f"  real(r8), parameter :: evaplimit_default_max_loss_fraction = {max_frac}\n"
        f"  integer,  parameter :: evaplimit_default_max_substeps = {max_sub}\n"
        "  real(r8) :: test_flux_ratio = 0._r8\n"
        "  real(r8) :: test_flux_scale = -1._r8\n"
        "  abstract interface\n"
        "     real(r8) FUNCTION tracer_evap_ratio_callback (source_ratio, temp_k, from_ice)\n"
        "        IMPORT :: r8\n"
        "        real(r8), intent(in) :: source_ratio\n"
        "        real(r8), intent(in) :: temp_k\n"
        "        logical,  intent(in) :: from_ice\n"
        "     END FUNCTION tracer_evap_ratio_callback\n"
        "  end interface\n"
        "CONTAINS\n"
        + body
        + "\n"
        "  real(r8) FUNCTION fixed_ratio_cb (source_ratio, temp_k, from_ice)\n"
        "    real(r8), intent(in) :: source_ratio\n"
        "    real(r8), intent(in) :: temp_k\n"
        "    logical,  intent(in) :: from_ice\n"
        "    IF (test_flux_scale >= 0._r8) THEN\n"
        "       fixed_ratio_cb = test_flux_scale * source_ratio\n"
        "    ELSE\n"
        "       fixed_ratio_cb = test_flux_ratio\n"
        "    ENDIF\n"
        "  END FUNCTION fixed_ratio_cb\n"
        "end module tracer_evaplimit\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program evaplimit_driver
  use tracer_evaplimit
  implicit none
  real(r8) :: pool_trc, pool_water, water_loss, ratio, r_max, skin, scale
  integer  :: mode, stat

  do
    read(*, *, iostat=stat) mode, pool_trc, pool_water, water_loss, ratio, r_max, skin, scale
    if (stat /= 0) exit
    test_flux_ratio = ratio
    test_flux_scale = scale
    if (mode == 0) then
      write(*, '(ES24.16)') tracer_evaporative_tracer_loss(pool_trc, pool_water, &
        water_loss, 293.15_r8, .false., fixed_ratio_cb, 1.e-30_r8, r_max)
    else
      write(*, '(ES24.16)') tracer_skin_limited_tracer_loss(pool_trc, pool_water, &
        water_loss, skin, 293.15_r8, .true., fixed_ratio_cb, 1.e-30_r8, r_max)
    end if
  end do
end program evaplimit_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "evaplimit"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _evaplimit(executable: Path, cases, mode: int = 0, scale: float = -1.0) -> list:
    """cases: (pool_trc, pool_water, water_loss, ratio, r_max[, skin_mass]).

    mode=0 exercises tracer_evaporative_tracer_loss, mode=1 the skin-limited
    wrapper.  scale >= 0 makes the stub callback return scale*source_ratio
    (a fractionating evaporate) instead of the fixed `ratio`.
    """
    cases = list(cases)
    rows = []
    for case in cases:
        fields = list(case)
        skin = fields.pop(5) if len(fields) > 5 else 0.0
        rows.append(
            f"{mode} "
            + " ".join(repr(float(v)) for v in fields)
            + f" {float(skin)!r} {float(scale)!r}"
        )
    payload = "\n".join(rows) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    values = [float(tok) for tok in ran.stdout.split()]
    assert len(values) == len(cases), ran.stdout + ran.stderr
    return values


def test_limiter_passes_a_negative_ratio_through_as_a_pool_gain(evaplimit_driver):
    """A negative Craig-Gordon ratio must produce a negative 'loss', i.e. the
    caller's `pool - loss` adds tracer.  A non-negative clamp anywhere in the
    limiter would silently discard net heavy-isotope uptake."""
    pool_trc, pool_water, water_loss, ratio = 2.0, 100.0, 1.0, -0.4
    (value,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, ratio, 0.0)]
    )
    assert value == pytest.approx(water_loss * ratio, rel=1e-12), (
        f"expected {water_loss * ratio} (a gain), got {value}"
    )


def test_negative_ratio_can_uptake_into_zero_tracer_inventory(evaplimit_driver):
    (value,) = _evaplimit(evaplimit_driver, [(0.0, 100.0, 1.0, -0.4, 0.0)])
    assert value == pytest.approx(-0.4, rel=1e-12)

    (skinned,) = _evaplimit(
        evaplimit_driver, [(0.0, 100.0, 1.0, -0.4, 0.0, 2.0)], mode=1
    )
    assert skinned == pytest.approx(-0.4, rel=1e-12)


def test_substepping_and_fast_path_agree_for_a_constant_ratio(evaplimit_driver):
    """With a ratio independent of the pool state, the sub-stepped integral
    must equal the single-shot answer -- for both signs.

    This pins the sub-step loop's bookkeeping: a wrong EXIT condition on the
    negative branch would stop integrating early and silently lose part of
    the flux.
    """
    pool_trc, pool_water = 2.0, 100.0
    for ratio in (-0.4, 0.02):
        # water_loss below 10% of the pool takes the fast path; above it the
        # limiter sub-steps.
        (fast,) = _evaplimit(evaplimit_driver, [(pool_trc, pool_water, 1.0, ratio, 0.0)])
        (stepped,) = _evaplimit(
            evaplimit_driver, [(pool_trc, pool_water, 50.0, ratio, 0.0)]
        )
        assert fast == pytest.approx(1.0 * ratio, rel=1e-12)
        assert stepped == pytest.approx(50.0 * ratio, rel=1e-9), (
            f"sub-stepped integral for ratio={ratio} gave {stepped}, expected "
            f"{50.0 * ratio}; the loop is not conserving the constant ratio"
        )


def test_full_pool_evaporation_still_removes_exactly_the_inventory(evaplimit_driver):
    """The all-pool branch must stay exact: leaving tracer in a zero-water
    pool breaks conservation regardless of the ratio's sign."""
    for ratio in (-0.4, 0.02):
        (value,) = _evaplimit(evaplimit_driver, [(2.0, 100.0, 100.0, ratio, 0.0)])
        assert value == pytest.approx(2.0, rel=1e-12)

def test_evaplimit_clamps_single_step_at_rmax(evaplimit_driver):
    pool_water = 100.0
    r_max = 0.020
    pool_trc = 1.95
    water_loss = 5.0
    (loss,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, 0.0, r_max)], scale=0.5
    )
    residual_ratio = (pool_trc - loss) / (pool_water - water_loss)
    assert residual_ratio <= r_max + 1e-12
    assert residual_ratio == pytest.approx(r_max, rel=1e-12)


def test_evaplimit_rmax_boundary_is_non_fractionating(evaplimit_driver):
    pool_water = 100.0
    r_max = 0.020
    pool_trc = pool_water * r_max
    water_loss = 5.0
    (loss,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, 0.0, r_max)], scale=0.5
    )
    assert loss == pytest.approx(water_loss * r_max, rel=1e-12)
    assert (pool_trc - loss) / (pool_water - water_loss) == pytest.approx(r_max, rel=1e-12)


# ---------------------------------------------------------------------------
# Open-water kinetic fractionation: Merlivat & Jouzel (1979).
#
# IsoGSM gsml/ISOTOPE/frkin.F:12-23:
#   smooth regime (u <  7 m/s): k = 0.006
#   rough  regime (u >= 7 m/s): k = 0.000285*u + 0.00082
#   HDO takes 0.88 of the 18O value.
# The former fixed n=2/3 exponent gave eps_k(18O) ~ 21 permil, 3-6x the
# MJ79 value, over-enriching every lake/water-body evaporation event.
# ---------------------------------------------------------------------------

HDO_MJ79_RELATIVE = 0.88


@pytest.fixture(scope="module")
def mj79_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_mj79")
    compiler = require_runnable_fortran_compiler(workdir)

    params = "\n".join(
        f"  real(r8), parameter :: {name} = {_module_parameter(FRAC_SOURCE, name)}"
        for name in (
            "mj79_smooth_k",
            "mj79_rough_slope",
            "mj79_rough_offset",
            "mj79_wind_threshold",
        )
    )
    body = _extract_function(FRAC_SOURCE, "tracer_mj79_kinetic_alpha")

    (workdir / "closure.f90").write_text(
        "module tracer_mj79\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        + params
        + "\nCONTAINS\n"
        + body
        + "\nend module tracer_mj79\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program mj79_driver
  use tracer_mj79
  implicit none
  real(r8) :: wind, relative_factor
  integer  :: stat

  do
    read(*, *, iostat=stat) wind, relative_factor
    if (stat /= 0) exit
    write(*, '(ES24.16)') tracer_mj79_kinetic_alpha(wind, relative_factor)
  end do
end program mj79_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "mj79"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _mj79(executable: Path, cases) -> list:
    cases = list(cases)
    payload = "\n".join(f"{u!r} {f!r}" for u, f in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    values = [float(tok) for tok in ran.stdout.split()]
    assert len(values) == len(cases), ran.stdout + ran.stderr
    return values


def _mj79_closed_form(wind: float, relative_factor: float) -> float:
    k = 0.000285 * wind + 0.00082 if wind >= 7.0 else 0.006
    return 1.0 / (1.0 - k * relative_factor)


@pytest.mark.parametrize("wind", (0.0, 3.0, 6.999, 7.0, 10.0, 15.0, 25.0))
def test_mj79_matches_the_published_piecewise_form(mj79_driver, wind):
    (value,) = _mj79(mj79_driver, [(wind, 1.0)])
    assert value == pytest.approx(_mj79_closed_form(wind, 1.0), rel=1e-12)


def test_rough_regime_fractionates_less_than_smooth_regime(mj79_driver):
    """The defining feature of MJ79: crossing into the rough (windy) regime
    at 7 m/s DROPS the kinetic factor, because added turbulence shrinks the
    diffusive sublayer.  Implementing this as a smooth/monotone increase
    with wind -- the intuitive but wrong reading -- would invert it.
    """
    smooth, rough = _mj79(mj79_driver, [(6.999, 1.0), (7.0, 1.0)])
    assert rough < smooth, (
        f"at the 7 m/s threshold the kinetic factor must drop, got "
        f"{smooth:.9g} -> {rough:.9g}"
    )
    eps_smooth = (smooth - 1.0) * 1000.0
    eps_rough = (rough - 1.0) * 1000.0
    assert eps_smooth == pytest.approx(6.04, abs=0.1), f"got {eps_smooth:.2f} permil"
    assert eps_rough == pytest.approx(2.82, abs=0.1), f"got {eps_rough:.2f} permil"


def test_open_water_kinetic_effect_is_far_below_the_soil_exponent_value(mj79_driver):
    """Sanity against the value being replaced: n=2/3 on the Merlivat
    diffusivity ratio gives ~18.9 permil, MJ79 gives ~6 permil at low wind."""
    (calm,) = _mj79(mj79_driver, [(2.0, 1.0)])
    eps_mj79 = (calm - 1.0) * 1000.0
    eps_exponent = (1.0285 ** (2.0 / 3.0) - 1.0) * 1000.0
    assert eps_mj79 < 0.4 * eps_exponent, (
        f"MJ79 should be several times smaller than the n=2/3 exponent value: "
        f"{eps_mj79:.2f} vs {eps_exponent:.2f} permil"
    )


def test_hdo_keeps_the_published_0p88_ratio_to_o18(mj79_driver):
    """MJ79 fixes eps_k(HDO)/eps_k(18O) = 0.88 by construction -- the same
    0.88 that makes the Merlivat diffusivity pair self-consistent."""
    for wind in (0.0, 10.0, 20.0):
        o18, hdo = _mj79(mj79_driver, [(wind, 1.0), (wind, HDO_MJ79_RELATIVE)])
        ratio = (hdo - 1.0) / (o18 - 1.0)
        assert ratio == pytest.approx(HDO_MJ79_RELATIVE, abs=0.002), (
            f"at u={wind} m/s the HDO/18O kinetic ratio is {ratio:.4f}, "
            f"expected {HDO_MJ79_RELATIVE}"
        )


def test_negative_wind_is_treated_as_calm(mj79_driver):
    calm, negative = _mj79(mj79_driver, [(0.0, 1.0), (-5.0, 1.0)])
    assert negative == pytest.approx(calm, rel=1e-12)


# ---------------------------------------------------------------------------
# Sublimation: only a thin surface skin can exchange isotopically.
#
# Applying Craig-Gordon to a whole snow LAYER treats tens of mm of ice as one
# well-mixed reservoir, which over-enriches the pack: in reality only the top
# few mm exchange with the atmosphere and the rest is removed bodily as the
# snow surface retreats.  IsoGSM sidesteps this by not fractionating land/ice
# evaporation at all (moninp.F:384-389).  The skin mass interpolates between
# those two limits: skin -> infinity recovers the layer-mixed behaviour, skin
# -> 0 recovers IsoGSM's non-fractionating sublimation.
# ---------------------------------------------------------------------------

# A fractionating evaporate: the vapour leaves lighter than the source, so the
# residual pool enriches.  0.9 exaggerates a real alpha to keep the signal
# well clear of round-off.
FRACTIONATING_SCALE = 0.9


def test_skin_larger_than_the_loss_reproduces_layer_mixed_behaviour(evaplimit_driver):
    """When the whole sublimated mass fits inside the skin, nothing changes --
    the wrapper must be a no-op relative to the existing limiter."""
    pool_trc, pool_water, water_loss = 2.0, 100.0, 5.0
    (reference,) = _evaplimit(
        evaplimit_driver,
        [(pool_trc, pool_water, water_loss, 0.0, 0.0)],
        mode=0,
        scale=FRACTIONATING_SCALE,
    )
    (skinned,) = _evaplimit(
        evaplimit_driver,
        [(pool_trc, pool_water, water_loss, 0.0, 0.0, water_loss * 2.0)],
        mode=1,
        scale=FRACTIONATING_SCALE,
    )
    assert skinned == pytest.approx(reference, rel=1e-12)


def test_zero_skin_makes_sublimation_non_fractionating(evaplimit_driver):
    """skin=0 is the documented IsoGSM limit: the vapour carries the pool's
    own ratio, so the residual pack does not enrich at all."""
    pool_trc, pool_water, water_loss = 2.0, 100.0, 5.0
    (value,) = _evaplimit(
        evaplimit_driver,
        [(pool_trc, pool_water, water_loss, 0.0, 0.0, 0.0)],
        mode=1,
        scale=FRACTIONATING_SCALE,
    )
    assert value == pytest.approx(water_loss * pool_trc / pool_water, rel=1e-12), (
        "with no exchanging skin the loss must be exactly water-matched"
    )


def test_enrichment_grows_monotonically_with_skin_mass(evaplimit_driver):
    """More exchanging mass => more fractionation => less tracer leaves =>
    more residual enrichment.  A sign slip here would make thicker skins
    enrich less, which is the opposite of the intended physics."""
    pool_trc, pool_water, water_loss = 2.0, 100.0, 20.0
    skins = [0.0, 2.0, 5.0, 10.0, 20.0]
    losses = _evaplimit(
        evaplimit_driver,
        [(pool_trc, pool_water, water_loss, 0.0, 0.0, s) for s in skins],
        mode=1,
        scale=FRACTIONATING_SCALE,
    )
    for (s_small, l_small), (s_big, l_big) in zip(
        zip(skins, losses), zip(skins[1:], losses[1:])
    ):
        assert l_big < l_small, (
            f"skin={s_big} mm must retain MORE tracer (smaller loss) than "
            f"skin={s_small} mm, got {l_big:.9g} vs {l_small:.9g}"
        )


def test_skin_limited_loss_never_exceeds_the_inventory(evaplimit_driver):
    for skin in (0.0, 1.0, 1.0e6):
        (value,) = _evaplimit(
            evaplimit_driver,
            [(2.0, 100.0, 99.9, 0.0, 0.0, skin)],
            mode=1,
            scale=FRACTIONATING_SCALE,
        )
        assert value <= 2.0 + 1e-12, f"skin={skin}: loss {value} exceeds inventory"
        assert value > 0.0


def test_skin_limiter_carries_a_negative_ratio_through(evaplimit_driver):
    """The two features interact: sublimation now goes through the skin
    limiter, and Craig-Gordon may return a negative ratio (net heavy-isotope
    uptake).  Only the skin fraction can take up; the remainder leaves the
    pool bodily at the pool's own ratio, so the two parts must combine with
    their own signs rather than one clobbering the other.
    """
    pool_trc, pool_water, water_loss = 2.0, 100.0, 20.0
    ratio = -0.4

    # skin >= loss: the whole loss fractionates, so this must equal the plain
    # limiter's answer (a pure gain).
    (all_skin,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, ratio, 0.0, 50.0)], mode=1
    )
    (reference,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, ratio, 0.0)], mode=0
    )
    assert all_skin == pytest.approx(reference, rel=1e-12)
    assert all_skin < 0.0

    # skin = 0: nothing fractionates, so the uptake disappears entirely and
    # the pool simply loses tracer in proportion to the water removed.
    (no_skin,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, ratio, 0.0, 0.0)], mode=1
    )
    assert no_skin == pytest.approx(water_loss * pool_trc / pool_water, rel=1e-12)
    assert no_skin > 0.0

    # A partial skin must land strictly between the two, not outside them.
    (partial,) = _evaplimit(
        evaplimit_driver, [(pool_trc, pool_water, water_loss, ratio, 0.0, 5.0)], mode=1
    )
    assert all_skin < partial < no_skin, (
        f"partial skin gave {partial:.6g}, outside the bracket "
        f"[{all_skin:.6g}, {no_skin:.6g}]"
    )


def test_skin_limited_full_pool_sublimation_is_exact(evaplimit_driver):
    """Draining the pool must move the whole inventory regardless of skin,
    otherwise tracer is stranded in a zero-water layer."""
    for skin in (0.0, 3.0, 1.0e6):
        (value,) = _evaplimit(
            evaplimit_driver,
            [(2.0, 100.0, 100.0, 0.0, 0.0, skin)],
            mode=1,
            scale=FRACTIONATING_SCALE,
        )
        assert value == pytest.approx(2.0, rel=1e-12)


# ---------------------------------------------------------------------------
# Soil evaporation kinetic factor: resistance-weighted rather than a fixed
# exponent.
#
# A fixed n=2/3 pins soil evaporation at ~19 permil regardless of how dry the
# surface is.  Physically the evaporation front retreats into the pores as the
# soil dries, lengthening the purely diffusive path, so the kinetic effect
# should migrate from the turbulent limit (n=2/3, ~19 permil) toward the
# stagnant-diffusion limit (n=1, ~28.5 permil).  Weighting by CoLM's own
# aerodynamic and soil-surface resistances does exactly that, and matches the
# form already used for leaves (tracer_alpha_kinetic_leaf).
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def soil_kinetic_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_soil_kinetic")
    compiler = require_runnable_fortran_compiler(workdir)

    params = "\n".join(
        f"  real(r8), parameter :: {name} = {_module_parameter(FRAC_SOURCE, name)}"
        for name in (
            "craig_gordon_kinetic_exponent_liquid",
            "craig_gordon_kinetic_exponent_ice",
        )
    )
    body = _extract_function(FRAC_SOURCE, "tracer_soil_kinetic_alpha_core")

    (workdir / "closure.f90").write_text(
        "module tracer_soil_kinetic\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        + params
        + "\nCONTAINS\n"
        + body
        + "\nend module tracer_soil_kinetic\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program soil_kinetic_driver
  use tracer_soil_kinetic
  implicit none
  real(r8) :: ra, rs, diff_ratio
  integer  :: stat

  do
    read(*, *, iostat=stat) ra, rs, diff_ratio
    if (stat /= 0) exit
    write(*, '(ES24.16)') tracer_soil_kinetic_alpha_core(ra, rs, diff_ratio)
  end do
end program soil_kinetic_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "soil_kinetic"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _soil_kinetic(executable: Path, cases) -> list:
    cases = list(cases)
    payload = "\n".join(" ".join(repr(float(v)) for v in c) for c in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    values = [float(tok) for tok in ran.stdout.split()]
    assert len(values) == len(cases), ran.stdout + ran.stderr
    return values


DIFF_O18 = 1.0285


def test_wet_soil_recovers_the_turbulent_two_thirds_limit(soil_kinetic_driver):
    """No soil-surface resistance means a purely turbulent path: the factor
    must equal D**(2/3), i.e. exactly the value the fixed exponent gave."""
    (value,) = _soil_kinetic(soil_kinetic_driver, [(50.0, 0.0, DIFF_O18)])
    assert value == pytest.approx(DIFF_O18 ** (2.0 / 3.0), rel=1e-12)


def test_dry_soil_approaches_the_stagnant_diffusion_limit(soil_kinetic_driver):
    """A large soil-surface resistance means the path is dominated by pore
    diffusion: the factor must approach the bare diffusivity ratio D."""
    (value,) = _soil_kinetic(soil_kinetic_driver, [(50.0, 1.0e6, DIFF_O18)])
    assert value == pytest.approx(DIFF_O18, rel=1e-4)
    eps = (value - 1.0) * 1000.0
    assert eps == pytest.approx(28.5, abs=0.2), f"got {eps:.2f} permil"


def test_kinetic_factor_increases_monotonically_as_soil_dries(soil_kinetic_driver):
    """The whole point of the change: drying must strengthen the kinetic
    effect instead of leaving it pinned at the wet-soil value."""
    resistances = [0.0, 10.0, 50.0, 200.0, 1000.0, 10000.0]
    values = _soil_kinetic(
        soil_kinetic_driver, [(50.0, rs, DIFF_O18) for rs in resistances]
    )
    for (rs_wet, v_wet), (rs_dry, v_dry) in zip(
        zip(resistances, values), zip(resistances[1:], values[1:])
    ):
        assert v_dry > v_wet, (
            f"rss={rs_dry} s/m must fractionate more than rss={rs_wet} s/m, "
            f"got {v_dry:.9g} vs {v_wet:.9g}"
        )
    assert values[0] == pytest.approx(DIFF_O18 ** (2.0 / 3.0), rel=1e-12)
    assert values[-1] < DIFF_O18


def test_degenerate_resistances_give_no_fractionation_rather_than_nan(
    soil_kinetic_driver,
):
    (value,) = _soil_kinetic(soil_kinetic_driver, [(0.0, 0.0, DIFF_O18)])
    assert value == pytest.approx(1.0, rel=1e-12)
    (negative,) = _soil_kinetic(soil_kinetic_driver, [(-10.0, -10.0, DIFF_O18)])
    assert negative == pytest.approx(1.0, rel=1e-12)


# ---------------------------------------------------------------------------
# Liquid-phase isotope diffusion in the soil column.
#
# Both CoLM and IsoGSM previously moved soil-water isotopes by advection only.
# Without molecular diffusion the model cannot produce the peak-shaped delta
# profile around an evaporation front (Barnes & Allison 1983), which is the
# single most-compared feature of soil water isotope observations.
# ---------------------------------------------------------------------------

# Wang et al. (1953) liquid self-diffusivity of H2-18O near 25 C.
D_LIQ_O18_298 = 2.28e-9


@pytest.fixture(scope="module")
def soil_diffusion_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_soil_diffusion")
    compiler = require_runnable_fortran_compiler(workdir)

    bodies = "\n".join(
        (
            _extract_function(FRAC_SOURCE, "tracer_soil_effective_diffusivity"),
            _extract_function(FRAC_SOURCE, "tracer_soil_diffusive_transfer"),
        )
    )

    (workdir / "closure.f90").write_text(
        "module tracer_soil_diffusion\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        "CONTAINS\n"
        + bodies
        + "\nend module tracer_soil_diffusion\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program soil_diffusion_driver
  use tracer_soil_diffusion
  implicit none
  real(r8) :: ru, rl, wu, wl, dzu, dzl, por, dliq, dt
  real(r8) :: du, dl
  integer  :: stat

  do
    read(*, *, iostat=stat) ru, rl, wu, wl, dzu, dzl, por, dliq, dt
    if (stat /= 0) exit
    du = tracer_soil_effective_diffusivity(wu, dzu, por, dliq)
    dl = tracer_soil_effective_diffusivity(wl, dzl, por, dliq)
    write(*, '(3ES24.16)') &
      tracer_soil_diffusive_transfer(ru, rl, wu, wl, du, dl, dzu, dzl, dt), du, dl
  end do
end program soil_diffusion_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "soil_diffusion"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _diffuse(executable: Path, cases) -> list:
    """cases: (R_up, R_low, w_up, w_low, dz_up, dz_low, porsl, D_liq, dt)
    -> [(transfer_up_to_low, D_eff_up, D_eff_low)]"""
    cases = list(cases)
    payload = "\n".join(" ".join(repr(float(v)) for v in c) for c in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    nums = [float(tok) for tok in ran.stdout.split()]
    assert len(nums) == 3 * len(cases), ran.stdout + ran.stderr
    return list(zip(nums[0::3], nums[1::3], nums[2::3]))


# A representative pair of 10 cm layers holding 30 mm each (theta = 0.30).
BASE = dict(w=30.0, dz=0.1, por=0.45, dliq=D_LIQ_O18_298, dt=1800.0)


def _case(r_up, r_low, **kw):
    p = {**BASE, **kw}
    return (r_up, r_low, p["w"], p["w"], p["dz"], p["dz"], p["por"], p["dliq"], p["dt"])


def test_no_gradient_means_no_diffusion(soil_diffusion_driver):
    (transfer, _, _), = _diffuse(soil_diffusion_driver, [_case(2.0e-3, 2.0e-3)])
    assert transfer == pytest.approx(0.0, abs=1e-30)


def test_diffusion_runs_from_heavy_to_light(soil_diffusion_driver):
    """Sign convention: a positive return moves tracer from upper to lower."""
    heavy_above, light_above = _diffuse(
        soil_diffusion_driver,
        [_case(2.05e-3, 2.00e-3), _case(2.00e-3, 2.05e-3)],
    )
    assert heavy_above[0] > 0.0, "an enriched upper layer must export downward"
    assert light_above[0] < 0.0, "an enriched lower layer must export upward"
    assert heavy_above[0] == pytest.approx(-light_above[0], rel=1e-12), (
        "the closure must be antisymmetric in the gradient"
    )


def test_effective_diffusivity_uses_millington_quirk_tortuosity(soil_diffusion_driver):
    """D_eff = D_liq * theta**(7/3) / porsl**2, and it must fall steeply as
    the soil dries -- that steepness is what keeps dry layers from mixing."""
    (_, d_wet, _), = _diffuse(soil_diffusion_driver, [_case(2.0e-3, 2.0e-3, w=45.0)])
    (_, d_mid, _), = _diffuse(soil_diffusion_driver, [_case(2.0e-3, 2.0e-3, w=30.0)])
    (_, d_dry, _), = _diffuse(soil_diffusion_driver, [_case(2.0e-3, 2.0e-3, w=3.0)])

    def expected(theta):
        return D_LIQ_O18_298 * theta ** (7.0 / 3.0) / 0.45**2

    assert d_wet == pytest.approx(expected(0.45), rel=1e-10)
    assert d_mid == pytest.approx(expected(0.30), rel=1e-10)
    assert d_dry == pytest.approx(expected(0.03), rel=1e-10)
    assert d_dry < d_mid < d_wet


def test_bone_dry_layer_does_not_diffuse(soil_diffusion_driver):
    (transfer, d_up, _), = _diffuse(
        soil_diffusion_driver,
        [(2.05e-3, 2.0e-3, 0.0, 30.0, 0.1, 0.1, 0.45, D_LIQ_O18_298, 1800.0)],
    )
    assert d_up == pytest.approx(0.0, abs=1e-30)
    assert transfer == pytest.approx(0.0, abs=1e-30)


def test_transfer_never_overshoots_equilibrium(soil_diffusion_driver):
    """Explicit diffusion must not create a new extremum, even for a thin
    layer plus a long timestep where the plain CFL condition is violated."""
    r_up, r_low, w = 2.10e-3, 2.00e-3, 5.0
    # 1 mm layers and a 30-minute step: strongly CFL-violating.
    (transfer, _, _), = _diffuse(
        soil_diffusion_driver,
        [(r_up, r_low, w, w, 0.001, 0.001, 0.45, D_LIQ_O18_298, 1800.0)],
    )
    equalise = (r_up - r_low) / (1.0 / w + 1.0 / w)
    assert 0.0 < transfer <= equalise * (1.0 + 1e-12), (
        f"transfer {transfer:.6g} must not exceed the equalising amount "
        f"{equalise:.6g}; the limiter is missing and the scheme will ring"
    )
    new_r_up = (r_up * w - transfer) / w
    new_r_low = (r_low * w + transfer) / w
    assert new_r_up >= new_r_low - 1e-18, "layers crossed over: overshoot"


def test_diffusive_timescale_matches_the_analytic_estimate(soil_diffusion_driver):
    """Order-of-magnitude guard on units.  For theta=0.3 in 10 cm layers the
    e-folding time L**2/D_eff is ~170 days; a unit slip (m vs mm, or a missing
    1000) would move this by three orders of magnitude."""
    r_up, r_low = 2.05e-3, 2.00e-3
    (transfer, d_eff, _), = _diffuse(soil_diffusion_driver, [_case(r_up, r_low)])
    # d(R_up)/dt from this step, converted to an e-folding time of the contrast.
    d_ratio_per_step = transfer / BASE["w"]
    tau_seconds = (r_up - r_low) / (d_ratio_per_step / BASE["dt"])
    tau_days = tau_seconds / 86400.0
    assert 20.0 < tau_days < 2000.0, (
        f"diffusive e-folding time came out as {tau_days:.3g} days "
        f"(D_eff={d_eff:.3g} m2/s); expected O(100) days"
    )


# ---------------------------------------------------------------------------
# Two-way equilibrium exchange between a wet surface and ambient vapour.
#
# A wet leaf keeps exchanging molecules with the surrounding vapour even when
# the NET water flux is zero, which relaxes its ratio toward the equilibrium
# value alpha_eq*R_vapor.  IsoGSM carries the analogous idea for falling
# raindrops (a constant equilibration degree eqf = 0.95, plus the full
# Marshall-Palmer model in eqm_deg.F); CoLM had no zero-water-flux exchange
# anywhere.  Because water does not move but tracer does, this is a genuine
# boundary flux and needs its own accumulator.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def equilibration_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_equilibration")
    compiler = require_runnable_fortran_compiler(workdir)

    body = _extract_function(FRAC_SOURCE, "tracer_equilibration_exchange")

    (workdir / "closure.f90").write_text(
        "module tracer_equilibration\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        "CONTAINS\n"
        + body
        + "\nend module tracer_equilibration\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program equilibration_driver
  use tracer_equilibration
  implicit none
  real(r8) :: pool_trc, pool_water, vapor_ratio, alpha_eq, frac
  integer  :: stat

  do
    read(*, *, iostat=stat) pool_trc, pool_water, vapor_ratio, alpha_eq, frac
    if (stat /= 0) exit
    write(*, '(ES24.16)') tracer_equilibration_exchange(pool_trc, pool_water, &
      vapor_ratio, alpha_eq, frac)
  end do
end program equilibration_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "equilibration"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _equilibrate(executable: Path, cases) -> list:
    """cases: (pool_trc, pool_water, R_vapor, alpha_eq, fraction) -> [gain]"""
    cases = list(cases)
    payload = "\n".join(" ".join(repr(float(v)) for v in c) for c in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    values = [float(tok) for tok in ran.stdout.split()]
    assert len(values) == len(cases), ran.stdout + ran.stderr
    return values


ALPHA_EQ_298 = 1.0093736282405357  # Majoube at 298.15 K
R_VAP = 0.988                       # ambient vapour, delta = -12 permil
W_CANOPY = 0.4                      # mm of intercepted water


def test_zero_fraction_disables_the_exchange(equilibration_driver):
    """The shipped default is 0, so the term must be an exact no-op then."""
    (gain,) = _equilibrate(
        equilibration_driver,
        [(W_CANOPY * 0.995, W_CANOPY, R_VAP, ALPHA_EQ_298, 0.0)],
    )
    assert gain == 0.0


def test_full_equilibration_lands_exactly_on_the_equilibrium_ratio(
    equilibration_driver,
):
    (gain,) = _equilibrate(
        equilibration_driver,
        [(W_CANOPY * 0.995, W_CANOPY, R_VAP, ALPHA_EQ_298, 1.0)],
    )
    new_ratio = (W_CANOPY * 0.995 + gain) / W_CANOPY
    assert new_ratio == pytest.approx(ALPHA_EQ_298 * R_VAP, rel=1e-12)


def test_a_pool_already_in_equilibrium_does_not_move(equilibration_driver):
    equilibrium_ratio = ALPHA_EQ_298 * R_VAP
    inventory = W_CANOPY * equilibrium_ratio
    (gain,) = _equilibrate(
        equilibration_driver,
        [(inventory, W_CANOPY, R_VAP, ALPHA_EQ_298, 0.95)],
    )
    # Tolerance scaled to the inventory: forming ratio_old = trc/water and
    # comparing it back against alpha_eq*R_vapor costs a couple of ulp.
    assert abs(gain) <= 1.0e-14 * inventory, (
        f"a pool at equilibrium must not move, got gain={gain:.6g} against an "
        f"inventory of {inventory:.6g}"
    )


def test_exchange_pulls_toward_equilibrium_from_both_sides(equilibration_driver):
    """Sign check: a pool lighter than equilibrium gains heavy isotopes and a
    heavier pool loses them.  Getting this backwards would drive the canopy
    away from equilibrium instead of toward it."""
    equilibrium_ratio = ALPHA_EQ_298 * R_VAP
    light = equilibrium_ratio - 0.01
    heavy = equilibrium_ratio + 0.01
    gain_light, gain_heavy = _equilibrate(
        equilibration_driver,
        [
            (W_CANOPY * light, W_CANOPY, R_VAP, ALPHA_EQ_298, 0.5),
            (W_CANOPY * heavy, W_CANOPY, R_VAP, ALPHA_EQ_298, 0.5),
        ],
    )
    assert gain_light > 0.0, "a pool lighter than equilibrium must gain"
    assert gain_heavy < 0.0, "a pool heavier than equilibrium must lose"
    assert gain_light == pytest.approx(-gain_heavy, rel=1e-12)


def test_exchange_scales_with_pool_size_and_fraction(equilibration_driver):
    args = (R_VAP, ALPHA_EQ_298)
    ratio = 0.995
    base, double_water = _equilibrate(
        equilibration_driver,
        [(1.0 * ratio, 1.0, *args, 0.5), (2.0 * ratio, 2.0, *args, 0.5)],
    )
    assert double_water == pytest.approx(2.0 * base, rel=1e-12)

    half, full = _equilibrate(
        equilibration_driver,
        [(1.0 * ratio, 1.0, *args, 0.5), (1.0 * ratio, 1.0, *args, 1.0)],
    )
    assert full == pytest.approx(2.0 * half, rel=1e-12)


def test_fraction_outside_zero_one_is_clamped(equilibration_driver):
    args = (W_CANOPY * 0.995, W_CANOPY, R_VAP, ALPHA_EQ_298)
    at_one, above_one = _equilibrate(equilibration_driver, [(*args, 1.0), (*args, 4.0)])
    assert above_one == pytest.approx(at_one, rel=1e-12), (
        "an over-unity equilibration degree must clamp, not overshoot past "
        "equilibrium and oscillate"
    )
    (negative,) = _equilibrate(equilibration_driver, [(*args, -1.0)])
    assert negative == 0.0


def test_meltwater_ice_equilibration_depletes_the_percolating_water(
    equilibration_driver,
):
    """Snowmelt percolation reuses the same closure, but the target is the
    ice-water equilibrium R_ice / alpha_ice_liq, supplied as
    alpha_eq = 1/alpha_ice_liq with the ICE ratio in the vapour slot.

    alpha_ice_liq > 1, so ice is enriched relative to water in equilibrium.
    The percolating water must therefore come out LIGHTER than the ice it
    passes through -- that is the observed early-melt depletion (Taylor et
    al. 2001).  Passing alpha_ice_liq instead of its reciprocal would flip
    this and enrich the meltwater.
    """
    temp_k = 271.15
    alpha_ice_liq = _alpha_ice_vap(temp_k) / _majoube_alpha_liq_vap(temp_k)
    assert alpha_ice_liq > 1.0, "sanity: ice is enriched relative to liquid"

    r_ice = 1.002
    water_mm = 3.0
    # Water starts at the ice ratio, then exchanges toward equilibrium.
    (gain,) = _equilibrate(
        equilibration_driver,
        [(water_mm * r_ice, water_mm, r_ice, 1.0 / alpha_ice_liq, 1.0)],
    )
    r_new = (water_mm * r_ice + gain) / water_mm
    assert gain < 0.0, "percolating meltwater must lose heavy isotopes to the ice"
    assert r_new == pytest.approx(r_ice / alpha_ice_liq, rel=1e-12)
    assert r_new < r_ice
    depletion_permil = (r_new / r_ice - 1.0) * 1000.0
    assert -6.0 < depletion_permil < -1.0, (
        f"expected a few permil of depletion, got {depletion_permil:.2f}"
    )


def test_dry_pool_exchanges_nothing(equilibration_driver):
    (gain,) = _equilibrate(
        equilibration_driver, [(0.0, 0.0, R_VAP, ALPHA_EQ_298, 0.95)]
    )
    assert gain == 0.0


# ---------------------------------------------------------------------------
# Vapour-phase isotope diffusion in soil pores.
#
# Liquid diffusion alone still cannot make the Barnes & Allison (1983)
# profile: once the surface dries, the evaporation front retreats into the
# soil and transport above it is through the PORE AIR, not the water film.
# Vapour diffusivity is ~1e4 times the liquid value while pore vapour density
# is ~2e-5 of liquid water, so the two are comparable overall -- and their
# theta dependences are OPPOSITE, so vapour takes over exactly where liquid
# shuts down.
#
# Implemented as an equivalent liquid diffusivity so the existing flux and
# overshoot limiter are reused unchanged:
#   D_vap_equiv = D_vap_air * tau_gas / diff_ratio * rho_v_sat / (alpha * rho_w)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def vapor_diffusion_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_vapor_diffusion")
    compiler = require_runnable_fortran_compiler(workdir)

    tfrz = _module_parameter(FRAC_SOURCE, "tfrz")
    bodies = "\n".join(
        (
            _extract_function(FRAC_SOURCE, "tracer_saturation_vapor_pressure"),
            _extract_function(FRAC_SOURCE, "tracer_soil_effective_diffusivity"),
            _extract_function(FRAC_SOURCE, "tracer_soil_vapor_equivalent_diffusivity"),
        )
    )

    (workdir / "closure.f90").write_text(
        "module tracer_vapor_diffusion\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        f"  real(r8), parameter :: tfrz = {tfrz}\n"
        "CONTAINS\n"
        + bodies
        + "\nend module tracer_vapor_diffusion\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program vapor_diffusion_driver
  use tracer_vapor_diffusion
  implicit none
  real(r8) :: w, ice, dz, por, temp_k, psrf, dr, alpha, dliq
  integer  :: stat

  do
    read(*, *, iostat=stat) w, ice, dz, por, temp_k, psrf, dr, alpha, dliq
    if (stat /= 0) exit
    write(*, '(2ES24.16)') &
      tracer_soil_vapor_equivalent_diffusivity(w, ice, dz, por, temp_k, psrf, dr, alpha), &
      tracer_soil_effective_diffusivity(w, dz, por, dliq)
  end do
end program vapor_diffusion_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "vapor_diffusion"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


VD_POR = 0.45
VD_DZ = 0.1
VD_T = 298.15
VD_P = 101325.0
VD_ALPHA = 1.0093736282405357


def _vapor_diff(executable: Path, cases) -> list:
    """cases: (water_mm, ice_mm, dz, porsl, T, psrf, diff_ratio, alpha, D_liq)
    -> [(D_vap_equiv, D_liq_eff)]"""
    cases = list(cases)
    payload = "\n".join(" ".join(repr(float(v)) for v in c) for c in cases) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    nums = [float(tok) for tok in ran.stdout.split()]
    assert len(nums) == 2 * len(cases), ran.stdout + ran.stderr
    return list(zip(nums[0::2], nums[1::2]))


def _vd_case(theta, theta_ice=0.0, **kw):
    p = dict(por=VD_POR, dz=VD_DZ, T=VD_T, psrf=VD_P, dr=1.0285, alpha=VD_ALPHA)
    p.update(kw)
    water = theta * 1000.0 * p["dz"]
    ice = theta_ice * 1000.0 * p["dz"]
    return (water, ice, p["dz"], p["por"], p["T"], p["psrf"], p["dr"], p["alpha"],
            D_LIQ_O18_298)


def test_saturated_soil_closes_the_vapour_pathway(vapor_diffusion_driver):
    """No air-filled porosity means no vapour transport."""
    (vap, liq), = _vapor_diff(vapor_diffusion_driver, [_vd_case(VD_POR)])
    assert vap == pytest.approx(0.0, abs=1e-30)
    assert liq > 0.0, "the liquid pathway must be wide open when saturated"


def test_vapour_and_liquid_have_opposite_moisture_dependence(vapor_diffusion_driver):
    """This is the property that makes vapour worth adding: it takes over
    exactly where the liquid film shuts down.  If both increased with theta
    the term would be redundant."""
    thetas = [0.05, 0.10, 0.20, 0.30, 0.40]
    results = _vapor_diff(vapor_diffusion_driver, [_vd_case(t) for t in thetas])
    vap = [v for v, _ in results]
    liq = [l for _, l in results]
    for (t_dry, v_dry), (t_wet, v_wet) in zip(zip(thetas, vap), zip(thetas[1:], vap[1:])):
        assert v_dry > v_wet, (
            f"vapour diffusivity must DECREASE with wetness: theta={t_dry} gave "
            f"{v_dry:.3e}, theta={t_wet} gave {v_wet:.3e}"
        )
    for (t_dry, l_dry), (t_wet, l_wet) in zip(zip(thetas, liq), zip(thetas[1:], liq[1:])):
        assert l_dry < l_wet, "liquid diffusivity must INCREASE with wetness"


def test_vapour_dominates_in_dry_soil_and_is_negligible_when_wet(
    vapor_diffusion_driver,
):
    """Order-of-magnitude guard on the unit chain.  A missing rho_w or a
    swapped diffusivity would move this ratio by orders of magnitude and the
    dry-soil crossover -- the whole point of the term -- would vanish."""
    (v_dry, l_dry), = _vapor_diff(vapor_diffusion_driver, [_vd_case(0.05)])
    (v_wet, l_wet), = _vapor_diff(vapor_diffusion_driver, [_vd_case(0.40)])
    assert v_dry / l_dry > 5.0, (
        f"vapour must dominate in dry soil, got ratio {v_dry / l_dry:.3g}"
    )
    assert v_wet / l_wet < 0.05, (
        f"vapour must be negligible in wet soil, got ratio {v_wet / l_wet:.3g}"
    )


def test_vapour_diffusivity_matches_the_closed_form(vapor_diffusion_driver):
    (vap, _), = _vapor_diff(vapor_diffusion_driver, [_vd_case(0.20)])

    theta = 0.20
    d_vap_air = 2.12e-5 * (VD_T / 273.15) ** 2 * (101325.0 / VD_P)
    tau_gas = (VD_POR - theta) ** (7.0 / 3.0) / VD_POR**2
    esat = 611.2 * math.exp(17.67 * (VD_T - 273.15) / (243.5 + VD_T - 273.15))
    rho_v = esat / (461.5 * VD_T)
    expected = d_vap_air * tau_gas / 1.0285 * rho_v / (VD_ALPHA * 1000.0)
    assert vap == pytest.approx(expected, rel=1e-9)


def test_heavier_isotope_diffuses_more_slowly_in_the_vapour_phase(
    vapor_diffusion_driver,
):
    (light, _), = _vapor_diff(vapor_diffusion_driver, [_vd_case(0.20, dr=1.0)])
    (heavy, _), = _vapor_diff(vapor_diffusion_driver, [_vd_case(0.20, dr=1.0285)])
    assert heavy < light
    assert light / heavy == pytest.approx(1.0285, rel=1e-12)


def test_vapour_transport_strengthens_with_temperature_and_weakens_with_pressure(
    vapor_diffusion_driver,
):
    """rho_v_sat follows Clausius-Clapeyron, and gas diffusivity scales as
    1/P -- both are easy to drop when assembling the unit chain."""
    cold, warm = _vapor_diff(
        vapor_diffusion_driver, [_vd_case(0.20, T=283.15), _vd_case(0.20, T=303.15)]
    )
    assert warm[0] > 2.0 * cold[0], (
        f"a 20 K warming should more than double pore vapour transport, got "
        f"{cold[0]:.3e} -> {warm[0]:.3e}"
    )
    sea, alt = _vapor_diff(
        vapor_diffusion_driver, [_vd_case(0.20), _vd_case(0.20, psrf=50000.0)]
    )
    assert alt[0] == pytest.approx(sea[0] * 101325.0 / 50000.0, rel=1e-9)


def test_pore_ice_blocks_the_vapour_pathway(vapor_diffusion_driver):
    """Ice occupies pore space just as water does.  Charging air-filled
    porosity with only the LIQUID content overstates vapour transport in
    frozen soil by more than an order of magnitude (porsl=0.45,
    theta_liq=0.10, theta_ice=0.25 -> 18.6x), exactly where a seasonally
    frozen column is trying to hold an isotope gradient.
    """
    (no_ice, _), = _vapor_diff(vapor_diffusion_driver, [_vd_case(0.10)])
    (some_ice, _), = _vapor_diff(
        vapor_diffusion_driver, [_vd_case(0.10, theta_ice=0.25)]
    )
    assert some_ice < no_ice, "pore ice must reduce vapour diffusion"

    # Air-filled porosity is porsl - theta_liq - theta_ice, so the ratio is
    # set purely by the Millington-Quirk exponent on those two values.
    expected_ratio = ((0.45 - 0.10) / (0.45 - 0.10 - 0.25)) ** (7.0 / 3.0)
    assert no_ice / some_ice == pytest.approx(expected_ratio, rel=1e-9), (
        f"ice is not being subtracted from air-filled porosity: ratio "
        f"{no_ice / some_ice:.3g}, expected {expected_ratio:.3g}"
    )

    # Ice filling the remaining pore space shuts the pathway completely.
    (frozen, _), = _vapor_diff(
        vapor_diffusion_driver, [_vd_case(0.10, theta_ice=0.35)]
    )
    assert frozen == pytest.approx(0.0, abs=1e-30)


# ---------------------------------------------------------------------------
# Vapour-phase isotope diffusion within snow / firn.
#
# Isotopes barely move through the ice lattice; a snowpack smooths its own
# vertical signal through the PORE VAPOUR (Whillans & Grootes 1985; Johnsen et
# al. 2000).  Two things differ from the soil case: pore air is saturated with
# respect to ICE, and the carrier the gradient acts on is the ice, so the
# equilibrium coefficient is alpha_ice_vap.  Porosity comes from snow density
# against the density of ice, not from a soil porosity parameter.
# ---------------------------------------------------------------------------

RHO_ICE = 917.0


@pytest.fixture(scope="module")
def firn_diffusion_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("tracer_firn_diffusion")
    compiler = require_runnable_fortran_compiler(workdir)

    tfrz = _module_parameter(FRAC_SOURCE, "tfrz")
    bodies = "\n".join(
        (
            _extract_function(FRAC_SOURCE, "tracer_saturation_vapor_pressure"),
            _extract_function(FRAC_SOURCE, "tracer_snow_vapor_equivalent_diffusivity"),
        )
    )

    (workdir / "closure.f90").write_text(
        "module tracer_firn_diffusion\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        f"  real(r8), parameter :: tfrz = {tfrz}\n"
        "CONTAINS\n"
        + bodies
        + "\nend module tracer_firn_diffusion\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program firn_diffusion_driver
  use tracer_firn_diffusion
  implicit none
  real(r8) :: wliq, wice, dz, temp_k, psrf, dr, alpha
  integer  :: stat

  do
    read(*, *, iostat=stat) wliq, wice, dz, temp_k, psrf, dr, alpha
    if (stat /= 0) exit
    write(*, '(ES24.16)') tracer_snow_vapor_equivalent_diffusivity( &
      wliq, wice, dz, temp_k, psrf, dr, alpha)
  end do
end program firn_diffusion_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "firn_diffusion"
    compiled = subprocess.run(
        [compiler, "closure.f90", "driver.f90", "-o", str(executable)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    if compiled.returncode != 0:
        pytest.fail(compiled.stdout + compiled.stderr)
    return executable


def _firn(executable: Path, cases) -> list:
    """cases: (wliq_mm, wice_mm, dz, T, psrf, diff_ratio, alpha_ice_vap)"""
    cases = list(cases)
    payload = "\n".join(" ".join(repr(float(v)) for v in c) for c in cases) + "\n"
    ran = subprocess.run(
        [str(executable)], input=payload, capture_output=True, text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    values = [float(tok) for tok in ran.stdout.split()]
    assert len(values) == len(cases), ran.stdout + ran.stderr
    return values


FIRN_T = 263.15
FIRN_DZ = 0.05


def _alpha_ice_vap(temp_k):
    return math.exp(11.839 / temp_k - 0.028224)


def _firn_case(rho_snow, liquid_fraction=0.0, dz=FIRN_DZ, temp_k=FIRN_T,
               psrf=101325.0, dr=1.0285):
    """Snow of bulk density rho_snow, optionally with part of it as liquid."""
    total_mm = rho_snow * dz          # kg/m2 == mm water equivalent
    wliq = total_mm * liquid_fraction
    wice = total_mm - wliq
    return (wliq, wice, dz, temp_k, psrf, dr, _alpha_ice_vap(temp_k))


def test_firn_diffusivity_is_in_the_published_range(firn_diffusion_driver):
    """Johnsen et al. (2000) put firn vapour diffusivity around
    1e-11 to 5e-11 m2/s for densities of a few hundred kg/m3.  This is the
    end-to-end unit check on the whole chain."""
    (value,) = _firn(firn_diffusion_driver, [_firn_case(300.0)])
    assert 5.0e-12 < value < 1.0e-10, f"got {value:.3e} m2/s, outside firn range"


def test_denser_snow_diffuses_more_slowly(firn_diffusion_driver):
    densities = [150.0, 300.0, 500.0, 700.0]
    values = _firn(firn_diffusion_driver, [_firn_case(r) for r in densities])
    for (r_light, v_light), (r_dense, v_dense) in zip(
        zip(densities, values), zip(densities[1:], values[1:])
    ):
        assert v_dense < v_light, (
            f"rho={r_dense} must diffuse slower than rho={r_light}, got "
            f"{v_dense:.3e} vs {v_light:.3e}"
        )


def test_snow_at_ice_density_closes_the_pathway(firn_diffusion_driver):
    """Solid ice has no pore space, so there is no vapour pathway left."""
    (value,) = _firn(firn_diffusion_driver, [_firn_case(RHO_ICE)])
    assert value == pytest.approx(0.0, abs=1e-30)


def test_porosity_uses_ice_density_for_the_ice_fraction(firn_diffusion_driver):
    """Ice is 917 kg/m3, not 1000.  Charging its volume at 1000 would
    understate the ice fraction by ~8% and overstate porosity."""
    rho = 500.0
    (value,) = _firn(firn_diffusion_driver, [_firn_case(rho)])
    phi = 1.0 - rho / RHO_ICE
    d_vap_air = 2.12e-5 * (FIRN_T / 273.15) ** 2
    tc = FIRN_T - 273.15
    esat_ice = 611.2 * math.exp(22.46 * tc / (272.62 + tc))
    rho_v = esat_ice / (461.5 * FIRN_T)
    expected = (d_vap_air * phi ** (7.0 / 3.0) / 1.0285
                * rho_v / (_alpha_ice_vap(FIRN_T) * 1000.0))
    assert value == pytest.approx(expected, rel=1e-9)


def test_liquid_water_in_the_pack_also_blocks_pores(firn_diffusion_driver):
    """Wet snow at the same bulk density has less pore air than dry snow,
    because liquid at 1000 kg/m3 occupies less volume per mm than ice does --
    so the two fractions must be converted separately, not lumped."""
    dry, wet = _firn(
        firn_diffusion_driver,
        [_firn_case(400.0), _firn_case(400.0, liquid_fraction=0.5)],
    )
    assert wet > dry, (
        "at fixed bulk mass, replacing ice with liquid frees pore volume "
        f"(liquid is denser), so porosity rises: got dry={dry:.3e} wet={wet:.3e}"
    )


def test_firn_diffusion_uses_saturation_over_ice_not_water(firn_diffusion_driver):
    """Over ice, esat is lower than over supercooled water at the same
    sub-freezing temperature.  Using the water curve would overstate pore
    vapour density, here by ~10% at -10 C."""
    (value,) = _firn(firn_diffusion_driver, [_firn_case(300.0)])
    tc = FIRN_T - 273.15
    esat_ice = 611.2 * math.exp(22.46 * tc / (272.62 + tc))
    esat_water = 611.2 * math.exp(17.67 * tc / (243.5 + tc))
    assert esat_ice < esat_water, "sanity: ice curve must sit below water curve"
    over_water = value * esat_water / esat_ice
    assert value < over_water * 0.95, (
        "the closure appears to be using saturation over water"
    )


def test_firn_diffusion_strengthens_with_temperature(firn_diffusion_driver):
    cold, warm = _firn(
        firn_diffusion_driver,
        [_firn_case(300.0, temp_k=253.15), _firn_case(300.0, temp_k=270.15)],
    )
    assert warm > 2.0 * cold, (
        f"pore vapour transport must rise steeply with temperature, got "
        f"{cold:.3e} -> {warm:.3e}"
    )


def test_kinetic_weakening_deepens_as_temperature_drops(ice_alpha_driver):
    """Colder air holds more supersaturation, so the kinetic damping grows."""
    temps = [213.15, 223.15, 233.15, 243.15, 253.15]
    results = _evaluate(
        ice_alpha_driver, [(t, JM84_SLOPE, DIFF_RATIO_O18) for t in temps]
    )
    damping = [eff / eq for eff, eq in results]
    for (t_cold, d_cold), (t_warm, d_warm) in zip(
        zip(temps, damping), zip(temps[1:], damping[1:])
    ):
        assert d_cold < d_warm, (
            f"damping alpha_eff/alpha_eq must strengthen with cooling, but "
            f"T={t_cold} K gives {d_cold:.9g} and T={t_warm} K gives {d_warm:.9g}"
        )
