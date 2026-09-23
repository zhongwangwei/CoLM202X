"""Runtime (behavioural) tests for the sediment critical-shear closure.

Why this file exists
--------------------
``calc_critical_shear_vel_sq`` implements Iwagaki's piecewise fit for the
critical shear velocity.  The fit is *designed* to be continuous at its
breakpoints -- that is the property that pins down the exponents.  A single
mistyped exponent (``31/32`` instead of ``31/22``) leaves the code compiling,
leaves every static/text assertion green, and silently raises the entrainment
threshold of the 1.18-3.03 mm classes -- the bed-load workhorse -- by a factor
of 1.7-2.6.

Static tests cannot see that.  These tests compile the *real* shipped function
and check the mathematical properties the fit must satisfy.

Scope and limits
----------------
``MOD_Tracer_Particle_Sediment`` is ``PRIVATE`` with 46 ``USE`` dependencies, so
stubbing the whole module to reach one pure leaf function is not worth the
coupling.  Instead the function body is extracted verbatim from the shipped
source and compiled standalone.  That means:

* the *constants and formula under test are the real ones* -- editing the source
  changes what these tests see, which is the point;
* what is **not** covered here is whether callers use the result correctly
  (unit handling, ``sqrt``/``* 0.01`` conversion at the call sites).

If the function is renamed or removed, extraction fails and the tests error out
rather than silently passing.
"""

from pathlib import Path
import re
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
SEDIMENT_SOURCE = ROOT / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"
FUNCTION_NAME = "calc_critical_shear_vel_sq"
SUBPROCESS_TIMEOUT = 60

# Iwagaki branch boundaries, in metres, as written in the shipped source.
BREAKPOINTS_M = (0.00303, 0.00118, 0.000565, 0.000065)

# The piecewise fit is continuous by construction.  With the correct exponents
# every breakpoint matches to ~2%; the mistyped 31/32 exponent produces 73% and
# 161% jumps, so this threshold separates the two cases by a wide margin.
CONTINUITY_TOL = 0.05


def _extract_function(source_text: str, name: str) -> str:
    """Return the verbatim text of a module function, header through END."""
    start = re.search(
        rf"^\s*(?:real\(r8\)\s+)?FUNCTION\s+{name}\s*\(", source_text, re.MULTILINE | re.IGNORECASE
    )
    if start is None:
        raise AssertionError(
            f"{name} not found in {SEDIMENT_SOURCE}. If it was renamed or removed, "
            "update this test rather than deleting it -- the continuity property "
            "still needs an owner."
        )
    end = re.search(
        rf"^\s*END\s+FUNCTION\s+{name}\b", source_text[start.start():], re.MULTILINE | re.IGNORECASE
    )
    if end is None:
        raise AssertionError(f"unterminated FUNCTION {name} in {SEDIMENT_SOURCE}")
    return source_text[start.start(): start.start() + end.end()]


@pytest.fixture(scope="module")
def critical_shear_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("sediment_critical_shear")
    compiler = require_runnable_fortran_compiler(workdir)

    body = _extract_function(SEDIMENT_SOURCE.read_text(encoding="utf-8"), FUNCTION_NAME)

    (workdir / "closure.f90").write_text(
        "module sediment_critical_shear\n"
        "  implicit none\n"
        "  integer, parameter :: r8 = selected_real_kind(12)\n"
        "  public :: " + FUNCTION_NAME + "\n"
        "CONTAINS\n"
        + body
        + "\nend module sediment_critical_shear\n",
        encoding="utf-8",
    )

    (workdir / "driver.f90").write_text(
        """
program critical_shear_driver
  use sediment_critical_shear, only: calc_critical_shear_vel_sq
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
  real(r8) :: diam
  integer  :: stat

  do
    read(*, *, iostat=stat) diam
    if (stat /= 0) exit
    write(*, '(ES24.16)') calc_critical_shear_vel_sq(diam)
  end do
end program critical_shear_driver
""",
        encoding="utf-8",
    )

    executable = workdir / "critical_shear"
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


def _evaluate(executable: Path, diameters) -> list:
    payload = "\n".join(repr(float(d)) for d in diameters) + "\n"
    ran = subprocess.run(
        [str(executable)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    values = [float(line) for line in ran.stdout.split()]
    assert len(values) == len(list(diameters)), ran.stdout + ran.stderr
    return values


@pytest.mark.parametrize("breakpoint_m", BREAKPOINTS_M)
def test_piecewise_fit_is_continuous_at_each_breakpoint(critical_shear_driver, breakpoint_m):
    """Iwagaki's branches must join up; a wrong exponent tears them apart."""
    below = breakpoint_m * (1.0 - 1.0e-9)
    above = breakpoint_m * (1.0 + 1.0e-9)
    value_below, value_above = _evaluate(critical_shear_driver, (below, above))

    scale = max(abs(value_below), abs(value_above))
    assert scale > 0.0, f"degenerate critical shear at d={breakpoint_m} m"
    relative_jump = abs(value_above - value_below) / scale
    assert relative_jump < CONTINUITY_TOL, (
        f"discontinuity at d={breakpoint_m} m: {value_below:.6g} -> {value_above:.6g} "
        f"({relative_jump:.1%} jump). The Iwagaki branch exponents no longer join; "
        "check the exponent on the 0.00118-0.00303 m branch (must be 31/22, not 31/32)."
    )


def test_critical_shear_has_no_large_inversion_in_grain_size(critical_shear_driver):
    """Coarser grains are harder to entrain -- no branch may badly invert that.

    Iwagaki's branches are an empirical fit, so they join to ~2% rather than
    exactly; crossing a breakpoint can dip slightly.  That small seam is the fit
    itself and is tolerated here.  A mistyped exponent is a different animal: the
    31/32 typo drops the value ~42% across the 0.00303 m breakpoint, far outside
    the seam.
    """
    lo, hi, count = 1.0e-5, 1.0e-2, 400
    step = (hi / lo) ** (1.0 / (count - 1))
    diameters = [lo * step**i for i in range(count)]
    values = _evaluate(critical_shear_driver, diameters)

    for (d_prev, v_prev), (d_next, v_next) in zip(
        zip(diameters, values), zip(diameters[1:], values[1:])
    ):
        assert v_next >= v_prev * (1.0 - CONTINUITY_TOL), (
            f"critical shear fell {1.0 - v_next / v_prev:.1%} with increasing grain size: "
            f"d={d_prev:.6g} m -> {v_prev:.6g}, d={d_next:.6g} m -> {v_next:.6g}. "
            "That exceeds the Iwagaki branch-join seam; check the branch exponents."
        )


def test_bed_load_classes_are_finite_and_positive(critical_shear_driver):
    """The default sediment classes must all yield a usable threshold."""
    # Spans the shipped default grain sizes plus the coarse bed-load classes
    # that sit on the previously-broken branch.
    diameters = (1.0e-5, 6.5e-5, 5.65e-4, 1.18e-3, 2.0e-3, 3.03e-3, 1.0e-2)
    values = _evaluate(critical_shear_driver, diameters)
    for diameter, value in zip(diameters, values):
        assert value == value, f"NaN critical shear at d={diameter} m"
        assert 0.0 < value < 1.0e6, f"implausible critical shear {value} at d={diameter} m"


def test_two_millimetre_class_sits_on_the_repaired_branch(critical_shear_driver):
    """Regression guard for the exact class the 31/32 typo mis-scaled.

    2 mm falls in the 0.00118-0.00303 m branch.  Anchoring it against the
    neighbouring linear branch (80.9*d, which has no disputed exponent) pins the
    exponent without restating it.
    """
    anchor_d = 0.00303
    (anchor_value,) = _evaluate(critical_shear_driver, (anchor_d,))
    (value_2mm,) = _evaluate(critical_shear_driver, (0.002,))

    # 80.9 * (0.303 cm) is the unambiguous linear-branch value at the boundary.
    assert anchor_value == pytest.approx(80.9 * 0.303, rel=1.0e-12)
    # Continuity + monotonicity force the 2 mm value into this window; the
    # 31/32 typo lands it near 30, far outside.
    assert 12.0 < value_2mm < anchor_value, (
        f"2 mm critical shear {value_2mm:.6g} (cm/s)^2 is outside the range implied "
        f"by the adjacent linear branch ({anchor_value:.6g})"
    )
