"""Guards the TRACER conservation tolerance against being tightened past the
floating-point noise floor.

Background
----------
``tracer_balance_check`` is a *per-step* check::

    err = storage_end - storage_beg - step_input + step_output
    balance_tol = trc_balance_abs_tol + trc_balance_rel_tol * balance_scale

with ``balance_scale`` the largest magnitude among the storage and flux terms.
Because the check does not carry state across steps, its error does not grow
with run length -- but it *does* carry the round-off of summing the soil/snow
layers and the flux terms within one step.

Paired with ``DEF_TRACER_BALANCE_ABORT_NBAD = 0`` (a single violation aborts the
run), a tolerance set below that noise floor would abort a physically correct
multi-year run.  Such a failure needs a long integration to surface, so it is
exactly the class of regression the audit found the test suite blind to.

These tests pin the headroom as a numeric property rather than pinning the
constants themselves, so the values stay free to change for physical reasons
while staying above the noise floor.
"""

from pathlib import Path
import re

import pytest


ROOT = Path(__file__).resolve().parents[1]
CONSERVATION_SOURCE = ROOT / "main/TRACER/MOD_Tracer_Conservation.F90"

# IEEE-754 binary64 unit round-off.
DOUBLE_EPS = 2.220446049250313e-16

# Terms summed into one balance residual: soil layers + snow layers + the
# canopy/aquifer/surface pools + the per-process flux terms.  Deliberately
# generous -- the point is a bound, not an exact count.
BALANCE_TERM_COUNT = 64

# Required margin between the tolerance and the round-off bound.
REQUIRED_HEADROOM = 100.0

# Magnitudes spanning the shipped tracer families: trace isotopic signatures
# through bulk solute and sediment loads.
MAGNITUDE_SWEEP = (1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4, 1.0e6)


def _read_parameter(name: str) -> float:
    text = CONSERVATION_SOURCE.read_text(encoding="utf-8")
    match = re.search(
        rf"^\s*real\(r8\),\s*parameter\s*::\s*{name}\s*=\s*([0-9.eEdD+_-]+?)_r8\s*$",
        text,
        re.MULTILINE,
    )
    if match is None:
        raise AssertionError(
            f"{name} not found in {CONSERVATION_SOURCE}. If the tolerance moved or "
            "was renamed, update this test -- the headroom property still needs an owner."
        )
    return float(match.group(1).replace("d", "e").replace("D", "e"))


@pytest.fixture(scope="module")
def tolerances() -> tuple:
    return _read_parameter("trc_balance_abs_tol"), _read_parameter("trc_balance_rel_tol")


@pytest.mark.parametrize("scale", MAGNITUDE_SWEEP)
def test_balance_tolerance_clears_the_roundoff_floor(tolerances, scale):
    """A correct run must not trip the balance check on arithmetic noise alone."""
    abs_tol, rel_tol = tolerances
    tolerance = abs_tol + rel_tol * scale
    roundoff = BALANCE_TERM_COUNT * DOUBLE_EPS * scale

    assert tolerance >= REQUIRED_HEADROOM * roundoff, (
        f"at magnitude {scale:g} the balance tolerance {tolerance:g} leaves only "
        f"{tolerance / roundoff:.1f}x headroom over the {roundoff:g} round-off floor "
        f"(need {REQUIRED_HEADROOM:g}x). With DEF_TRACER_BALANCE_ABORT_NBAD = 0 this "
        "aborts physically correct long runs."
    )


def test_relative_term_governs_above_unit_scale(tolerances):
    """The tolerance must track magnitude, not sit at a fixed absolute value.

    A purely absolute tolerance is simultaneously too loose for trace species and
    too tight for bulk loads -- the defect the 5e-6 absolute tolerance had.
    """
    abs_tol, rel_tol = tolerances
    assert rel_tol > 0.0, "a relative term is required for scale-aware tolerance"

    crossover = abs_tol / rel_tol
    assert crossover <= 1.0e4, (
        f"the relative term only takes over above magnitude {crossover:g}; below that "
        "the check degenerates to a fixed absolute tolerance across every tracer family"
    )


def test_tolerance_is_not_so_loose_that_real_leaks_hide(tolerances):
    """The other side of the trade: keep the check able to see a real leak.

    A per-step relative tolerance at or above 1e-6 would let a systematic
    0.0001%-per-step drift pass unnoticed, which compounds to a large error over a
    multi-year integration.
    """
    _, rel_tol = tolerances
    assert rel_tol <= 1.0e-6, (
        f"relative tolerance {rel_tol:g} is too permissive for a per-step check; "
        "a systematic per-step drift below it would accumulate undetected"
    )
