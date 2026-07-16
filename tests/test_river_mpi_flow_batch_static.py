from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text(
    encoding="utf-8"
)
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text(
    encoding="utf-8"
)


def test_downstream_state_uses_one_single_mapping_batch() -> None:
    assert "type(worker_push_real8_field_type) :: downstream_state_fields(2)" in FLOW
    assert (
        "CALL worker_push_data (push_next2ucat, downstream_state_fields)" in FLOW
    )
    assert "push_next2ucat, wdsrf_ucat, wdsrf_next" not in FLOW
    assert "push_next2ucat, veloc_riv,  veloc_next" not in FLOW


def test_reservoir_flux_uses_one_multi_mapping_batch() -> None:
    assert "type(worker_push_real8_field_type) :: reservoir_flux_fields(2)" in FLOW
    assert (
        "CALL worker_push_data (push_ups2ucat, reservoir_flux_fields, mode = 'sum')"
        in FLOW
    )
    assert "push_ups2ucat, hflux_resv, hflux_sumups" not in FLOW
    assert "push_ups2ucat, mflux_resv, mflux_sumups" not in FLOW


def test_batched_flow_state_has_stable_pointer_targets() -> None:
    assert "allocatable, target :: wdsrf_ucat" in TIMEVARS
    assert "allocatable, target :: veloc_riv" in TIMEVARS
    assert "allocatable, target :: wdsrf_next" in FLOW
    assert "allocatable, target :: veloc_next" in FLOW
    assert "allocatable, target :: hflux_resv" in FLOW
    assert "allocatable, target :: mflux_resv" in FLOW
