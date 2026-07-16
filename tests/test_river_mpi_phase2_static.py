import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text(
    encoding="utf-8"
)
PUSH = (ROOT / "share/MOD_WorkerPushData.F90").read_text(encoding="utf-8")


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_bif_substep_uses_one_collective_for_dt_and_continuation() -> None:
    loop_entry = FLOW.split("DO WHILE (loop_active)", 1)[0].rsplit(
        "dt_res(:) = acctime_rnof", 1
    )[1]
    substep_loop = FLOW.split("DO WHILE (loop_active)", 1)[1].split(
        "! Keep restart-visible state", 1
    )[0]
    sync = routine(FLOW, "sync_global_routing_dt")

    assert "mpi_allreduce (MPI_IN_PLACE, loop_active" not in loop_entry
    assert "mpi_allreduce (MPI_IN_PLACE, loop_active" not in substep_loop
    assert "next_loop_active" in substep_loop
    assert "dt_reduce(2)" in sync
    assert "dt_reduce, 2, MPI_REAL8, MPI_MIN" in sync
    assert sync.count("CALL mpi_allreduce") == 1


def test_bif_loop_entry_uses_the_globally_identical_accumulator() -> None:
    loop_entry = FLOW.split("DO WHILE (loop_active)", 1)[0].rsplit(
        "dt_res(:) = acctime_rnof", 1
    )[1]

    assert "IF (.not. ieee_is_finite(acctime_rnof)) THEN" in loop_entry
    assert "loop_active = .true." in loop_entry
    assert "loop_active = acctime_rnof > 0._r8" in loop_entry
    assert "acctime_rnof > 0._r8 .or." not in loop_entry
    assert "mpi_allreduce (MPI_IN_PLACE, loop_active" not in loop_entry


def test_nonfinite_residual_cannot_skip_the_collective_fail_fast() -> None:
    loop_prefix = FLOW.split("DO WHILE (loop_active)", 1)[0]
    bif_entry = loop_prefix.rsplit("dt_res(:) = acctime_rnof", 1)[1]
    assert "IF (.not. ieee_is_finite(acctime_rnof)) THEN" in bif_entry
    assert "loop_active = .true." in bif_entry
    assert "loop_active = acctime_rnof > 0._r8" in bif_entry
    assert "ELSE\n            loop_active = any(dt_res > 0._r8)" in bif_entry
    sync = routine(FLOW, "sync_global_routing_dt")
    assert "dt_reduce(1) = -huge(1._r8)" in sync
    assert "non-finite routing residual" in sync


def test_dt_guards_do_not_compare_nan_under_invalid_traps() -> None:
    sync = routine(FLOW, "sync_global_routing_dt")

    assert "pathological_mask" in sync
    assert "IF (.not. ieee_is_finite(dt_all(i))) THEN" in sync
    assert "ELSEIF (dt_all(i) <= 0._r8) THEN" in sync
    assert "dt_all <= 0._r8 .or." not in sync
    assert ".not. ieee_is_finite(dt_reduce(1))" in sync
    assert "ieee_is_finite(dt_global) .or." not in sync


def test_combined_dt_reduction_matches_explicit_residual_update() -> None:
    cases = (
        ([90.0, 25.0], [10.0, 5.0]),
        ([25.0, 0.0], [25.0, 0.0]),
        ([math.nextafter(1.0, math.inf), 1.0], [1.0, 1.0]),
        ([1.0e-12, 2.0e-12, 0.0], [5.0e-13, 1.0e-12, 0.0]),
    )

    for residuals, candidates in cases:
        active = [i for i, residual in enumerate(residuals) if residual > 0.0]
        dt_global = min(candidates[i] for i in active)
        max_residual = max(residuals[i] for i in active)
        combined_next = (max_residual - dt_global) > 0.0
        explicit_next = any(
            (residuals[i] - dt_global) > 0.0 for i in active
        )
        assert combined_next == explicit_next


def test_three_coupled_flux_fields_use_one_batch_push() -> None:
    assert "worker_push_real8_field_type" in PUSH
    assert "worker_push_data_multi_real8_batch" in PUSH
    assert "type(worker_push_real8_field_type) :: upstream_flux_fields(3)" in FLOW
    assert "CALL worker_push_data (push_ups2ucat, upstream_flux_fields, mode = 'sum')" in FLOW

    for old_call in (
        "CALL worker_push_data (push_ups2ucat, hflux_fc, hflux_sumups",
        "CALL worker_push_data (push_ups2ucat, mflux_fc, mflux_sumups",
        "CALL worker_push_data (push_ups2ucat, zgrad_dn, zgrad_sumups",
    ):
        assert FLOW.count(old_call) == 0


def test_batch_push_preserves_field_order_and_single_message_phase() -> None:
    batch = routine(PUSH, "worker_push_data_multi_real8_batch")
    assert batch.count("CALL mpi_isend") == 1
    assert batch.count("CALL mpi_irecv") == 1
    assert "fields(ifield)%send" in batch
    assert "fields(ifield)%recv" in batch
    assert "DO ifield = 1, size(fields)" in batch
    assert "batch_capacity" in PUSH


def test_zero_request_mapping_uses_allocated_zero_size_actuals() -> None:
    single = routine(PUSH, "build_worker_pushdata_single")
    multi = routine(PUSH, "build_worker_pushdata_multi")
    assert single.index("allocate (ids_req_uniq (num_req))") < single.index(
        "IF (num_req > 0) THEN"
    )
    assert multi.index("allocate (ids_req_uniq (ndim1*num_req))") < multi.index(
        "IF (num_req > 0) THEN"
    )
