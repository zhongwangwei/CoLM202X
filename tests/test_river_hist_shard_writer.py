"""Step-2 static checks for the IO-group route-history shard writer.

Plan step 2 requires that the block path never funnels route history through
the global master, and that a shard's identity lives in its attributes rather
than in its filename. Those are structural properties, so they are asserted
against the source; the runtime behaviour (zero-data groups, complete and
mutually exclusive id coverage) is covered by
tests/run_river_hist_shard_evaluator.sh.
"""

from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
VRW = ROOT / "main/HYDRO/MOD_Vector_ReadWrite.F90"
HARNESS = ROOT / "tests/river_hist_shard_harness.F90"
EVALUATOR = ROOT / "tests/run_river_hist_shard_evaluator.sh"

MASTER_GATHERERS = (
    "vector_gather_to_master",
    "vector_gather_map2grid_and_write",
    "vector_gather_matrix_to_master",
)


def _routine(name: str) -> str:
    src = VRW.read_text(encoding="utf-8")
    start = src.index(f"SUBROUTINE {name} ")
    end = src.index(f"END SUBROUTINE {name}", start)
    return src[start:end]


def _shard_routines() -> dict[str, str]:
    return {
        name: _routine(name)
        for name in (
            "route_shard_layout_build",
            "route_shard_write_vector",
            "route_shard_write_matrix",
        )
    }


def test_shard_writers_are_public() -> None:
    src = VRW.read_text(encoding="utf-8")
    for name in (
        "route_shard_layout_build",
        "route_shard_layout_free",
        "route_shard_write_vector",
        "route_shard_write_matrix",
        "route_shard_write_identity",
        "route_shard_filename",
        "route_shard_grid_fingerprint",
    ):
        assert f"PUBLIC :: {name}" in src, name


def test_shard_path_never_gathers_to_the_global_master() -> None:
    for name, body in _shard_routines().items():
        for gatherer in MASTER_GATHERERS:
            assert gatherer not in body, f"{name} calls {gatherer}"
        # the group communicator, never the global one
        assert "p_comm_glb" not in body, f"{name} uses p_comm_glb"


def test_shard_collectives_run_on_the_group_communicator() -> None:
    for name, body in _shard_routines().items():
        if "mpi_gather" in body:
            assert "p_comm_group" in body, name
            assert "p_root" in body, name


def test_master_takes_no_part() -> None:
    """The master forms a singleton group, so it must not enter a group
    collective; every branch is gated on p_is_io or p_is_worker."""
    for name, body in _shard_routines().items():
        assert "p_is_master" not in body, f"{name} branches on p_is_master"


def test_no_global_sized_buffers_on_the_shard_path() -> None:
    """IO ranks hold only their own group's share."""
    for name, body in _shard_routines().items():
        for forbidden in ("totalnumucat", "totalvlen", "totalnpthout"):
            assert forbidden not in body, f"{name} sizes a buffer with {forbidden}"


def test_zero_length_ranks_still_enter_the_collective() -> None:
    """A rank that owns nothing must present a valid send buffer and call
    gatherv anyway, or the collective deadlocks."""
    build = _routine("route_shard_layout_build")
    assert "allocate (sbuf (1))" in build
    vec = _routine("route_shard_write_vector")
    assert "allocate (sbuff (1));" in vec


def test_shard_identity_is_written_and_versioned() -> None:
    src = VRW.read_text(encoding="utf-8")
    assert "ROUTE_SHARD_SCHEMA_VERSION = 1" in src
    identity = _routine("route_shard_write_identity")
    for attr in (
        "river_hist_shard_schema_version",
        "target_history_basename",
        "case_name",
        "history_period_key",
        "segment_id",
        "shard_index",
        "shard_count",
        "grid_fingerprint",
        "time_first",
        "time_last",
        "time_count",
    ):
        assert attr in identity, attr
    # one open/close cycle, not one per attribute
    assert identity.count("nf90_open") == 1
    assert identity.count("nf90_close") == 1


def test_filename_is_not_the_source_of_truth() -> None:
    src = VRW.read_text(encoding="utf-8")
    assert "never has to infer correctness from a filename" in src
    fname = _routine("route_shard_filename")
    assert "_shard" in fname


def test_shard_count_follows_io_groups_not_workers() -> None:
    harness = HARNESS.read_text(encoding="utf-8")
    assert "shard_count = p_np_io" in harness
    assert "p_np_worker" not in re.sub(r"!.*", "", harness).split("shard_count =")[1][:200]


def test_evaluator_covers_zero_data_groups_and_several_layouts() -> None:
    script = EVALUATOR.read_text(encoding="utf-8")
    assert EVALUATOR.stat().st_mode & 0o111
    assert "mpiexec" in script or "mpirun" in script
    # several rank:group layouts, including more groups than data-owning workers
    assert '"3:1 4:1 6:2 8:2 9:3"' in script
    assert "RHSHARD PASS" in script


def test_harness_checks_coverage_is_complete_and_exclusive() -> None:
    harness = HARNESS.read_text(encoding="utf-8")
    assert "duplicate ucat id" in harness
    assert "duplicate pathway id" in harness
    assert "out of range" in harness
    assert "coverage broken" in harness
    # values must travel with their id, not with their rank
    assert "value misplaced" in harness


ROUTE_MOD = ROOT / "main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90"


def test_production_actually_writes_shard_identity() -> None:
    """The aggregator refuses a shard without identity attributes, so a writer
    that exists but is never called from production makes block mode unusable
    while every test still passes -- the harness wrote them itself.

    Assert the CALL SITE, not the routine's existence.
    """
    route = re.sub(r"!.*", "", ROUTE_MOD.read_text(encoding="utf-8"))
    assert "CALL route_shard_write_identity" in route, (
        "production never stamps shard identity; the aggregator will reject "
        "every real block-mode shard"
    )
    # and it must happen where the shard is created
    begin = route[route.index("SUBROUTINE route_hist_begin"):]
    begin = begin[: begin.index("END SUBROUTINE route_hist_begin")]
    assert "write_shard_identity" in begin


def test_shard_identity_carries_a_run_segment() -> None:
    """A restart inside one period writes a second segment of the same period;
    without a distinguishing id the aggregator cannot tell them apart."""
    route = ROUTE_MOD.read_text(encoding="utf-8")
    assert "rh_seg_date" in route
    assert "segment_id" in route


def test_restart_writes_a_new_segment_rather_than_reusing_shards() -> None:
    """A restart inside one history period must not reuse the previous run's
    shards: if the IO-group count changed they partition the domain
    differently and stale ones would linger. Each continuous run segment gets
    its own complete set, and the aggregator merges segments by time record.
    """
    route = ROUTE_MOD.read_text(encoding="utf-8")
    assert "rh_segment_id" in route
    # the segment must reach the filename, not only the attributes
    assert "route_shard_filename" in route
    # the call is continued across lines, so flatten before matching
    call = re.sub(r"&\s*\n\s*", " ", route[route.index("CALL route_shard_filename"):])
    call = call[: call.index("\n")]
    assert "rh_segment_id" in call, "segment must be part of the shard filename"

    vrw = (ROOT / "main/HYDRO/MOD_Vector_ReadWrite.F90").read_text(encoding="utf-8")
    sig = vrw[vrw.index("SUBROUTINE route_shard_filename"):]
    sig = sig[: sig.index("END SUBROUTINE route_shard_filename")]
    assert "segment" in sig
    assert "'_shard'" in sig



def test_harness_stamps_identity_through_the_production_routine() -> None:
    """The harness must not hand-roll identity: that is exactly what hid the
    fact that production never wrote it."""
    h = (ROOT / "tests/river_hist_shard_harness.F90").read_text(encoding="utf-8")
    assert "CALL route_shard_write_identity" in h
    assert "route_shard_grid_fingerprint" in h, "use the real fingerprint helper"
    assert "max(p_np_io,1)" in h, "shard_count must come from the real IO group count"
