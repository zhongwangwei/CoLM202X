"""Calendar-year inputs must not inherit the previous-day restart timestamp."""
from pathlib import Path
import re


def test_year_based_initializers_use_unshifted_calendar_year():
    source = (Path(__file__).resolve().parents[1] / "main/CoLM.F90").read_text()
    assert re.search(r"s_year\s*=\s*DEF_simulation_time%start_year", source)
    assert "CALL adj2end(sdate)" in source
    for routine in (
        "init_ndep_data_annually", "init_ndep_data_monthly",
        "init_fire_data", "grid_riverlake_flow_init",
    ):
        calls = re.findall(r"CALL\s+" + routine + r"\s*\(([^\n]*)", source)
        assert calls and all(re.match(r"s_year\s*[,)]", args) for args in calls), routine
