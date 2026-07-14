from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HARNESS = (ROOT / "tests/worker_push_mapping_benchmark.F90").read_text(encoding="utf-8")
RUNNER = (ROOT / "tests/run_worker_push_mapping_benchmark.sh").read_text(encoding="utf-8")


def test_mapping_benchmark_uses_three_trial_medians() -> None:
    assert "integer, parameter :: timing_trials = 3" in HARNESS
    assert "build_time = median_of_three(build_samples)" in HARNESS
    assert "push_time = median_of_three(push_samples)" in HARNESS


def test_mapping_benchmark_covers_empty_and_duplicate_requests() -> None:
    for case in (
        "single-zero",
        "single-empty-rank",
        "single-duplicate",
        "multi-zero",
        "multi-empty-rank",
        "multi-duplicate",
    ):
        assert case in HARNESS


def test_mapping_benchmark_defaults_through_32_ranks_and_makes_64_optional() -> None:
    assert 'WORKER_PUSH_BENCH_RANKS:-"2 4 8 16 32"' in RUNNER
    assert "WORKER_PUSH_BENCH_INCLUDE_64" in RUNNER
    assert "launcher_args+=(--oversubscribe)" in RUNNER


def test_mapping_benchmark_has_a_test_only_directory_ab_baseline() -> None:
    assert "CALL run_suite(.false.)" in HARNESS
    assert "CALL run_suite(.true.)" in HARNESS
    assert "MPI_Allgatherv" in HARNESS
    assert "MPI_Alltoallv" in HARNESS
    assert "directory_temp_bytes_max" in HARNESS
