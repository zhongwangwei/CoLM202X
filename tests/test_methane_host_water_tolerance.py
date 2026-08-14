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
    end = src.index("wliq_soisno_sat(j) =", start)
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
    assert "if (vtot > pore_volume + host_water_excess_tol) then" in block
    assert "if (vtot < finundated * pore_volume - host_water_tol) then" in block
    # and says what to set, not just that it failed
    assert "tolerance =" in block


def test_tolerated_band_is_reported_not_silent() -> None:
    block = _partition_block()
    assert "host_water_max_resid" in block
    assert "host_water_warned" in block
    assert "WARNING: methane host-water disagreement" in block


def test_host_water_itself_is_never_modified() -> None:
    """PR #11's TEMP clamp rescaled vtot without rescaling vliq/vice, breaking
    `vliq + vice == vtot`. The host state must be read, not rewritten."""
    block = _partition_block()
    for forbidden in (
        "vtot = pore_volume",
        "vtot = finundated * pore_volume",
        "vliq = vliq * pore_volume",
        "vice = vice * pore_volume",
    ):
        assert forbidden not in block, forbidden


def test_unsaturated_floor_is_not_load_bearing() -> None:
    """max(0,...) must be a safety net, not the conservation mechanism.

    It was previously relied on to 'leave the column dry rather than invent
    water' -- but zeroing the unsaturated side does not undo an over-allocated
    saturated side, which is exactly how water was being created. The cap on
    the allocation is what conserves; this asserts both are present.
    """
    src = _flat(PHYSICS.read_text(encoding="utf-8"))
    assert "vtot = vliq + vice" in src
    assert "wliq_soisno_unsat(j) = max(0._r8," in src
    block = _partition_block()
    assert "host_water_scale" in block


def test_temp_marker_from_pr11_is_not_present() -> None:
    assert "TEMP-METHANE-TOL" not in PHYSICS.read_text(encoding="utf-8")


def test_saturated_allocation_is_capped_so_no_water_is_created() -> None:
    """The conservation bug this file previously missed.

    With sat = pore*v/vtot and no cap, a column holding less than
    finundated*pore_volume gets an allocation the unsaturated side cannot pay
    for; its max(0,...) floor stops it going negative but does not undo the
    over-allocation, so the area-weighted total lands on finundated*pore_volume
    and (finundated*pore_volume - vtot) of water is created per layer.

    Numerical proof lives in tests/methane_host_water_harness.F90; this pins
    the source so the two cannot drift apart.
    """
    block = _partition_block()
    assert "host_water_scale = min(pore_volume / max(vtot, 1.e-12_r8)" in block
    assert "1._r8 / max(finundated, 1.e-12_r8))" in block
    assert "vliq_sat_alloc = vliq * host_water_scale" in block
    assert "vice_sat_alloc = vice * host_water_scale" in block
    # the un-capped form must not come back
    assert "vliq_sat_alloc = pore_volume * vliq" not in block


def test_numerical_harness_mirrors_the_source_formula() -> None:
    """A mirrored formula is only useful while it still mirrors."""
    harness = (ROOT / "tests/methane_host_water_harness.F90").read_text(encoding="utf-8")
    assert "min(pore_volume / max(vtot, 1.e-12_r8)" in harness
    assert "1._r8 / max(finundated, 1.e-12_r8))" in harness
    assert "MHW PASS" in harness


def test_report_collective_matches_the_ranks_that_reach_it() -> None:
    """CoLM.F90 calls CoLMDRIVER inside IF (p_is_worker), so tracer_report --
    and this reduction -- is worker-only. A p_comm_glb collective there
    deadlocks every MPI run with methane enabled."""
    src = PHYSICS.read_text(encoding="utf-8")
    body = src[src.index("subroutine methane_host_water_report"):]
    body = body[: body.index("end subroutine methane_host_water_report")]
    code = _flat(re.sub(r"!.*", "", body))
    assert "p_comm_worker" in code
    assert "p_comm_glb" not in code, "global collective on a worker-only path"
    # and it must not report from a rank that never arrives
    assert "p_is_master" not in code
    assert "p_iam_worker == 0" in body


def test_excess_side_keeps_a_roundoff_only_tolerance() -> None:
    """vtot > pore_volume cannot be partitioned conservatively, so it must not
    share the wide deficit tolerance."""
    block = _partition_block()
    assert "host_water_excess_tol = 1.e-9_r8 * pore_volume" in block
    assert "if (vtot > pore_volume + host_water_excess_tol) then" in block
    # the wide tolerance stays on the deficit side only
    assert "if (vtot < finundated * pore_volume - host_water_tol) then" in block
