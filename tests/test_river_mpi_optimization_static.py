from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text(
    encoding="utf-8"
)
PUSH = (ROOT / "share/MOD_WorkerPushData.F90").read_text(encoding="utf-8")
EVALUATOR = ROOT / "tests/run_river_mpi_evaluator.sh"
GNU_WORKFLOW = ROOT / ".github/workflows/build_CoLM_gnu.yml"
BIF = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90").read_text(
    encoding="utf-8"
)
NETWORK = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeNetwork.F90").read_text(
    encoding="utf-8"
)
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text(
    encoding="utf-8"
)


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_bif_substep_has_one_global_dt_reduction() -> None:
    substep_loop = FLOW.split("DO WHILE (loop_active)", 1)[1].split(
        "! Keep restart-visible state", 1
    )[0]
    assert substep_loop.count("CALL sync_global_routing_dt(dt_res, dt_all)") == 1


def test_real8_worker_push_reuses_communication_scratch() -> None:
    push_type = PUSH.split("type :: worker_pushdata_type", 1)[1].split(
        "END type worker_pushdata_type", 1
    )[0]
    for component in (
        "recv_uniq_real8",
        "sendcache_real8",
        "recvcache_real8",
        "req_send_real8",
        "req_recv_real8",
        "send_peer_real8",
        "recv_peer_real8",
    ):
        assert component in push_type

    real8_push = routine(PUSH, "worker_push_data_uniq_real8")
    assert "allocate (sendcache" not in real8_push
    assert "allocate (recvcache" not in real8_push
    assert "allocate (req_send" not in real8_push
    assert "allocate (req_recv" not in real8_push
    assert "pushdata%sendcache_real8" in real8_push
    assert "pushdata%recvcache_real8" in real8_push
    assert "pushdata%send_peer_real8" in real8_push
    assert "pushdata%recv_peer_real8" in real8_push
    assert "DO iworker = 0, p_np_worker-1" not in real8_push


def test_zero_data_worker_does_not_search_unallocated_ids() -> None:
    builder = routine(PUSH, "build_worker_pushdata_uniq")
    self_lookup = builder.split("IF (n_req_uniq > 0) THEN", 1)[1].split(
        "pushdata%nself = count", 1
    )[0]
    assert "IF (num_me > 0) THEN" in self_lookup
    assert "find_in_sorted_list1" in self_lookup


def test_bif_does_not_send_globally_identical_destination_dt() -> None:
    calc = routine(BIF, "bifurcation_calc")
    assert "dt_ucat_buf" not in BIF
    assert "dt_dn_pth_buf" not in BIF
    assert "worker_push_data (push_bif_dn2pth, dt_ucat" not in calc
    assert "dt_dn = dt" in calc


def test_no_levee_path_skips_split_storage_pushes() -> None:
    calc = routine(BIF, "bifurcation_calc")
    split_push = calc.index(
        "CALL worker_push_data (push_bif_dn2pth, visible_storage_ucat"
    )
    guard = calc.rfind("IF (DEF_USE_LEVEE) THEN", 0, split_push)
    fallback = calc.index("visible_storage_dn_pth(:) = storage_dn_pth(:)", split_push)
    assert guard >= 0
    assert split_push < fallback


def test_levee_restart_reuses_preloaded_visible_volume() -> None:
    levee = (
        ROOT / "main/HYDRO/MOD_Grid_RiverLakeLevee.F90"
    ).read_text(encoding="utf-8")
    reader = routine(levee, "read_levee_restart")
    assert "visible_volume_preloaded" in reader
    assert "IF (.not. visible_volume_preloaded) THEN" in reader


def test_multi_rank_river_system_partition_invariant_is_checked() -> None:
    assert "invalid_rivsys_partition = minval(rivermouth) /= maxval(rivermouth)" in NETWORK
    assert "invalid river-system MPI partition" in NETWORK


def test_restart_read_resets_visible_volume_validity() -> None:
    reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")
    reset = reader.index("volwater_ucat_valid = .false.")
    probe = reader.index("ncio_var_exist(file_restart, 'volwater_ucat'")
    assert reset < probe


def test_mpi_evaluator_executes_multi_rank_harness() -> None:
    script = EVALUATOR.read_text(encoding="utf-8")
    harness = (ROOT / "tests/river_mpi_harness.F90").read_text(encoding="utf-8")
    assert EVALUATOR.stat().st_mode & 0o111
    assert "mpiexec" in script or "mpirun" in script
    assert "for ranks in 2 4" in script
    assert "river_mpi_harness" in script
    assert "MPI_Comm_split" in harness
    assert "is_master = rank == nranks - 1" in harness
    assert "is_io = rank == 0" in harness
    assert "p_is_worker = .not. is_master .and. .not. is_io" in harness
    assert "MPI_Barrier(MPI_COMM_WORLD" in harness


def test_ci_runs_python_and_mpi_tests_for_shell_changes() -> None:
    workflow = GNU_WORKFLOW.read_text(encoding="utf-8")
    assert "python -m pytest -q" in workflow
    assert "./tests/run_river_mpi_evaluator.sh" in workflow
    assert "'**/**.sh'" not in workflow
