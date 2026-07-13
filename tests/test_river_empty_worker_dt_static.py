from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NETWORK = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeNetwork.F90").read_text()
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text()
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text()


def test_empty_worker_river_system_metadata_is_allocated() -> None:
    assert "numrivsys = 0\n         allocate (irivsys (numucat))" in NETWORK
    assert "IF (numucat > 0) allocate (irivsys (numucat))" not in NETWORK


def test_empty_worker_primary_routing_state_is_allocated() -> None:
    for declaration in (
        "allocate (wdsrf_ucat (numucat))",
        "allocate (veloc_riv  (numucat))",
        "allocate (momen_riv  (numucat))",
    ):
        assert declaration in TIMEVARS

    assert "IF (numucat > 0)  allocate (wdsrf_ucat" not in TIMEVARS
    assert "IF (numucat > 0)  allocate (veloc_riv" not in TIMEVARS
    assert "IF (numucat > 0)  allocate (momen_riv" not in TIMEVARS


def test_positive_sub_ten_second_dt_is_not_raised() -> None:
    assert "ieee_is_finite" in FLOW
    assert "IF (dt_global < ROUTING_MIN_ADAPTIVE_DT" not in FLOW
    assert "dt_all = min(ROUTING_PATHOLOGICAL_DT_FALLBACK, dt_res)" in FLOW
