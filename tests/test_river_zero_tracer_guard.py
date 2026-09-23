"""Regression for a TRACER build with DEF_TRACER_NUM=0.

The river tracer initializer allocates no state in that configuration, but
water routing and levee handling still call these entry points.
"""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_zero_tracer_entry_points_return_before_unallocated_state():
    source = (ROOT / "main/TRACER/MOD_Tracer_RiverLake.F90").read_text()
    for name, first_access in (
        ("tracer_substep", "CALL ensure_tracer_substep_workspace"),
        ("tracer_diag_accumulate_substep", "CALL get_cell_volume"),
    ):
        body = routine(source, name)
        assert body.index("IF (ntracers <= 0) RETURN") < body.index(first_access)


def test_levee_call_never_slices_unallocated_pending_pool():
    source = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text()
    assert re.search(
        r"IF \(ntracers > 0\) CALL levee_tracer_repartition\(i,[\s\S]*?"
        r"pending_trc_pool = trc_inp_buf\(:, i\)\)",
        source,
    )
