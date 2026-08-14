"""Step-1 wiring and formula tests for the river-history sharding plan.

Covers the parts of ``.omx/plans/river-history-sharded-output.md`` step 1 that
do not need an MPI launcher, so they run in CI:

* the memory formulas the go/no-go quotes are pinned and agree with the
  numbers the Fortran harness prints;
* the harness and its driver exist, are wired together, and exercise the
  properties the plan requires (zero-length worker, two rank scales, gather
  and NetCDF timed separately, real writer rather than a reimplementation).
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from river_hist_memory_baseline import (  # noqa: E402
    bif_root_bytes,
    even_split,
    unitcat_root_bytes,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
HARNESS = ROOT / "tests/river_hist_baseline_harness.F90"
DRIVER = ROOT / "tests/run_river_hist_baseline.sh"
# The plan itself lives outside the repository (.omx/ is locally excluded), so
# the committed report is what the gate criteria are asserted against.
REPORT = ROOT / "tests/river_hist_step1_baseline_report.md"


def test_memory_formulas_match_the_plan() -> None:
    # plan 1.5: unitcat 8*(totalnumucat + nlon*nlat + max_local_ucat)
    assert unitcat_root_bytes(4000, 200, 100, 4000) == 8 * (4000 + 200 * 100 + 4000)
    # plan 1.5: BIF 8*npthlev*(totalnpthout + max_local_path)
    assert bif_root_bytes(3, 500, 500) == 8 * 3 * (500 + 500)


def test_memory_formulas_agree_with_the_fortran_harness() -> None:
    """The harness prints its own lower bound; the tool must reproduce it.

    The harness leaves worker 0 empty on purpose, so its effective worker count
    is one less than p_np_worker -- hence workers=1 for the 2-worker run.
    """
    assert unitcat_root_bytes(4000, 200, 100, even_split(4000, 1)) == 224000
    assert bif_root_bytes(3, 500, even_split(500, 1)) == 24000


def test_even_split_is_a_ceiling_and_respects_imbalance() -> None:
    assert even_split(4000, 3) == 1334
    assert even_split(4000, 0) == 4000
    assert even_split(1000, 4, imbalance=2.0) == 500
    # imbalance can never claim more than the whole domain
    assert even_split(1000, 2, imbalance=100.0) == 1000


def test_harness_uses_the_real_writer_not_a_reimplementation() -> None:
    src = HARNESS.read_text(encoding="utf-8")
    assert "USE MOD_Vector_ReadWrite" in src
    for entry in (
        "vector_gather_to_master",
        "vector_gather_map2grid_and_write",
        "vector_gather_matrix_to_master",
    ):
        assert entry in src, entry
    # real SPMD layout, not a hand-rolled rank split
    assert "CALL spmd_init ()" in src
    assert "CALL divide_processes_into_groups (ngroup)" in src


def test_harness_always_exercises_a_zero_length_worker() -> None:
    src = HARNESS.read_text(encoding="utf-8")
    assert "allocate (ucat_data_address(iw)%val (0))" in src
    assert "IF (iw == 0) THEN" in src


def test_harness_separates_mpi_and_netcdf_timing() -> None:
    src = HARNESS.read_text(encoding="utf-8")
    assert "RHBASE TIME gather_only_s" in src
    assert "RHBASE TIME full_write_s" in src
    assert "RHBASE TIME netcdf_est_s" in src
    assert "RHBASE TIME bif_matrix_s" in src


def test_harness_avoids_the_promoted_double_precision_mpi_wtime_trap() -> None:
    """-fdefault-real-8 without -fdefault-double-8 promotes mpif.h's
    DOUBLE PRECISION MPI_WTIME to real(16), which silently reads as 0."""
    src = HARNESS.read_text(encoding="utf-8")
    assert "system_clock" in src
    assert "t0 = mpi_wtime" not in src


def test_driver_measures_two_scales_with_repeats() -> None:
    script = DRIVER.read_text(encoding="utf-8")
    assert DRIVER.stat().st_mode & 0o111
    assert "mpiexec" in script or "mpirun" in script
    assert 'scales=${RH_SCALES:-"4 8"}' in script
    assert 'RH_RUNS:=3' in script
    assert "statistics.median" in script
    # links against the built model rather than recompiling a subset
    assert "libcolm.a" in script


def test_report_records_the_go_no_go_verdict_and_criteria() -> None:
    report = REPORT.read_text(encoding="utf-8")
    # all three gate criteria are stated, with their thresholds
    assert "go/no-go" in report
    assert "5%" in report
    assert "25%" in report
    # and the verdict is explicit, not left to the reader
    assert "暂不进入第 2 步" in report


def test_report_records_the_measurements_the_gate_still_needs() -> None:
    """A stopped gate is only actionable if it names what would open it."""
    report = REPORT.read_text(encoding="utf-8")
    for needed in ("端到端墙钟", "gather", "routing compute"):
        assert needed in report, needed
