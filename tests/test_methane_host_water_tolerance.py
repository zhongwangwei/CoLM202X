"""Pin the host-water/inundation tolerance behaviour.

Replaces the machine-epsilon guard that aborted global runs on their first
timestep, without adopting the silent clamp that PR #11 marked TEMP: the abort
survives, the threshold becomes a configurable fraction of pore volume, and the
tolerated band is reported rather than hidden.
"""

from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONST = ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90"
PHYSICS = ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Physics.F90"


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


def _partition_block() -> str:
    src = PHYSICS.read_text(encoding="utf-8")
    start = src.index("Partition host soil water between conditional columns")
    # Anchor on the actual statement, not the prose: the explanatory comment
    # inside this block quotes the same expression without spaces.
    end = src.index("vliq_sat_alloc = pore_volume * vliq / max", start)
    return _flat(src[start:end])


def test_tolerance_is_a_namelist_parameter_with_a_sane_default() -> None:
    const = CONST.read_text(encoding="utf-8")
    assert "real(r8) :: host_water_tolerance = 0.05_r8" in const
    # bounds-checked: 0 restores the aborting behaviour, >=1 disables the guard
    assert "host_water_tolerance < 0._r8" in const
    assert "host_water_tolerance >= 1._r8" in const
    # carried through the NaN/inf sweep like every other real parameter
    assert "DEF_METHANE%host_water_tolerance, &" in const


def test_threshold_scales_with_pore_volume_not_machine_epsilon() -> None:
    block = _partition_block()
    assert "host_water_tol = DEF_METHANE%host_water_tolerance * pore_volume" in block
    # the old flat 1e-10 absolute tolerance is gone
    assert "1.e-10_r8 * max(pore_volume, 1._r8)" not in block


def test_the_guard_still_aborts_beyond_tolerance() -> None:
    block = _partition_block()
    assert block.count("CALL CoLM_stop") == 2
    assert "if (vtot > pore_volume + host_water_tol) then" in block
    assert "if (vtot < finundated * pore_volume - host_water_tol) then" in block
    # and says what to set, not just that it failed
    assert "tolerance =" in block


def test_tolerated_band_is_reported_not_silent() -> None:
    block = _partition_block()
    assert "host_water_max_resid" in block
    assert "host_water_warned" in block
    assert "WARNING: methane host-water disagreement" in block


def test_no_clamping_so_host_mass_balance_is_preserved() -> None:
    """The module must neither create nor discard host water.

    PR #11's TEMP clamp rescaled vtot without rescaling vliq/vice, which both
    broke `vliq + vice == vtot` and silently changed the saturated allocation.
    Rescaling in the other direction would invent water outright.
    """
    block = _partition_block()
    for forbidden in (
        "vtot = pore_volume",
        "vtot = finundated * pore_volume",
        "vliq = vliq * pore_volume",
        "vice = vice * pore_volume",
    ):
        assert forbidden not in block, forbidden


def test_allocation_cannot_exceed_pore_volume_without_a_clamp() -> None:
    """vliq_sat_alloc = pore_volume*vliq/vtot with vliq <= vtot = vliq+vice.

    This is why no clamp is needed; assert the algebra the argument rests on
    is still what the source computes.
    """
    src = _flat(PHYSICS.read_text(encoding="utf-8"))
    assert "vtot = vliq + vice" in src
    assert "vliq_sat_alloc = pore_volume * vliq / max(vtot, 1.e-12_r8)" in src
    assert "vice_sat_alloc = pore_volume * vice / max(vtot, 1.e-12_r8)" in src
    # and the unsaturated branch stays dry rather than inventing water
    assert "wliq_soisno_unsat(j) = max(0._r8," in src


def test_temp_marker_from_pr11_is_not_present() -> None:
    assert "TEMP-METHANE-TOL" not in PHYSICS.read_text(encoding="utf-8")
