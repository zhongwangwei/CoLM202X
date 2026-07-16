from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PUSH = (ROOT / "share/MOD_WorkerPushData.F90").read_text(encoding="utf-8")


def routine(name: str) -> str:
    return PUSH.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_batch_dispatch_uses_explicit_mapping_kind() -> None:
    assert "mapping_kind = WORKER_PUSH_MAPPING_UNSET" in PUSH
    assert "pushdata%mapping_kind = WORKER_PUSH_MAPPING_SINGLE" in routine(
        "build_worker_pushdata_single"
    )
    assert "pushdata%mapping_kind = WORKER_PUSH_MAPPING_MULTI" in routine(
        "build_worker_pushdata_multi"
    )

    batch = routine("worker_push_data_multi_real8_batch")
    assert "character(len=*), intent(in), optional :: mode" in batch
    assert "SELECT CASE (pushdata%mapping_kind)" in batch
    assert "fields(ifield)%recv = pushdata%recv_uniq_real8(pushdata%addr_single)" in batch
    assert "CASE (WORKER_PUSH_MAPPING_MULTI)" in batch


def test_batch_contract_checks_mapping_and_field_shapes() -> None:
    validate = routine("validate_worker_push_real8_batch")
    assert "batch used with an untyped push mapping" in validate
    assert "null field in real8 batch" in validate
    assert "short send field in real8 batch" in validate
    assert "wrong receive field size in real8 batch" in validate
    assert "mode is invalid for a batched single mapping" in validate
    assert "mode is required for a batched multi mapping" in validate
    assert "pushdata%required_send_size" in validate
    assert "DO iworker = 0, p_np_worker-1" not in validate
    assert "maxval(pushdata%to_other" not in validate


def test_required_send_size_is_cached_with_mapping_lifecycle() -> None:
    assert "integer :: required_send_size = 0" in PUSH
    cache = routine("update_worker_push_required_send_size")
    assert "pushdata%required_send_size = 0" in cache
    assert "pushdata%to_other(iworker)%val" in cache
    assert "CALL update_worker_push_required_send_size (pushdata)" in routine(
        "build_worker_pushdata_uniq"
    )
    assert "CALL update_worker_push_required_send_size (pushdata_out)" in routine(
        "build_worker_pushdata_subset"
    )
    finalizer = routine("worker_pushdata_free_mem")
    assert "this%required_send_size = 0" in finalizer


def test_batch_still_uses_one_message_per_peer() -> None:
    batch = routine("worker_push_data_multi_real8_batch")
    assert batch.count("CALL mpi_isend") == 1
    assert batch.count("CALL mpi_irecv") == 1
    assert batch.count("CALL mpi_waitall") == 2


def test_multi_outputs_are_filled_before_the_zero_unique_request_guard() -> None:
    scalar = routine("worker_push_data_multi_real8")
    scalar_fill = scalar.index("vec_recv(:) = fillvalue")
    scalar_guard = scalar.index("IF (pushdata%num_req_uniq > 0) THEN")
    assert scalar_fill < scalar_guard

    batch = routine("worker_push_data_multi_real8_batch")
    batch_fill = batch.index("fields(ifield)%recv(:) = fields(ifield)%fillvalue")
    batch_guard = batch.index("IF (pushdata%num_req_uniq > 0) THEN")
    assert batch_fill < batch_guard


def test_multi_average_skips_nonpositive_area_before_accumulation() -> None:
    scalar = routine("worker_push_data_multi_real8")
    batch = routine("worker_push_data_multi_real8_batch")

    for implementation in (scalar, batch):
        assert "IF (pushdata%area_multi(i,j) <= 0._r8) CYCLE" in implementation
        assert "sumarea > 0._r8" in implementation


def test_grid_to_pset_average_skips_nonpositive_area() -> None:
    remap = routine("worker_remap_data_grid2pset_real8")
    assert "IF (area <= 0._r8) CYCLE" in remap
    assert "sumarea > 0._r8" in remap
