"""Routing retains the original 60-second upper bound before local limiters."""
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90"


def test_hydrodynamic_substep_has_60_second_upper_bound():
    src = SOURCE.read_text()
    step = src.split("dt_all(:) = min(dt_res(:), 60._r8)", 1)
    assert len(step) == 2
    assert "dt_all(:) = dt_res(:)" not in src
    assert step[1].index("DO i = 1, numucat") < step[1].index("dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))")
    assert step[1].index("dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))") < step[1].index("dt_res = dt_res - dt_all")
