"""Step-4 checks: the shard aggregator and the workflow that must invoke it.

The aggregator's numerical correctness is covered end-to-end by
tests/run_river_hist_concat_e2e.sh, which builds a real shard set and checks
every id landed in the right cell. These are the structural guarantees:

* shards are merged on their identity attributes, never on their filename;
* a mismatched schema version is refused rather than misread;
* a half-written aggregate can never take the target name;
* a block-mode run leaves a discoverable note that aggregation is still owed.
"""

from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
AGG = ROOT / "postprocess/RiverHistConcatenate.F90"
MAIN = ROOT / "postprocess/HistConcatenate.F90"
ROUTE = ROOT / "main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90"
WRAPPER = ROOT / "run/scripts/concatenate_history"
MAKEFILE = ROOT / "Makefile"


def _src(p: pathlib.Path) -> str:
    return p.read_text(encoding="utf-8")


def test_main_history_concatenate_covers_the_river_remap_variables() -> None:
    main = _src(MAIN)
    for var in (
        "f_wdpth_ucat_regrid",
        "f_veloc_riv_regrid",
        "f_discharge",
        "f_discharge_rivermouth_regrid",
        "f_floodfrc",
        "f_floodarea",
    ):
        assert f"'{var}'" in main, var


def test_shards_are_merged_on_identity_not_filename() -> None:
    agg = _src(AGG)
    for attr in (
        "river_hist_shard_schema_version",
        "case_name",
        "history_period_key",
        "grid_fingerprint",
        "target_history_basename",
        "shard_count",
    ):
        assert attr in agg, attr
    assert "identity_error" in agg
    assert "never by the\n!    filename" in agg or "never by the" in agg


def test_unsupported_schema_version_is_refused() -> None:
    agg = _src(AGG)
    assert "SUPPORTED_SCHEMA_VERSION = 1" in agg
    assert "unsupported shard schema" in agg


def test_incomplete_or_overlapping_coverage_stops_the_job() -> None:
    agg = _src(AGG)
    for failure in (
        "id out of range",
        "duplicate id",
        "incomplete coverage",
        "duplicate pathway id",
        "incomplete pathway coverage",
        "incomplete shard set",
    ):
        assert failure in agg, failure


def test_pathways_are_keyed_on_global_id() -> None:
    agg = _src(AGG)
    assert "pth_global_id" in agg
    body = agg[agg.index("SUBROUTINE rebuild_bif_matrix"):]
    body = body[: body.index("END SUBROUTINE rebuild_bif_matrix")]
    assert "mat(l, ids(k))" in body, "pathway column must be indexed by global id"


def test_a_half_written_aggregate_cannot_take_the_target_name() -> None:
    agg = _src(AGG)
    assert "'.tmp'" in agg
    assert "'.complete'" in agg
    promote = agg[agg.index("SUBROUTINE verify_and_promote"):]
    promote = promote[: promote.index("END SUBROUTINE verify_and_promote")]
    # verify the temporary opens, only then rename onto the target
    assert promote.index("nf90_open") < promote.index("CALL rename_file")
    assert "shards were NOT deleted" in promote


def test_block_mode_leaves_a_pending_note_removed_only_on_success() -> None:
    route = _src(ROUTE)
    assert "route_hist_write_pending_manifest" in route
    assert "concatenate_history" in route
    assert "not yet aggregated" in route
    agg = _src(AGG)
    assert "'.pending'" in agg
    assert "status='old'" in agg and "status='delete'" in agg


def test_wrapper_handles_both_halves_and_reports_failure() -> None:
    assert WRAPPER.stat().st_mode & 0o111
    w = _src(WRAPPER)
    assert "hist_concatenate.x" in w
    assert "river_hist_concatenate.x" in w
    assert "exit $status" in w
    assert "MUST RUN THIS BEFORE ANALYSING" in w
    # shards are never deleted automatically
    assert "NOT deleted" in w or "left in place" in w


def test_aggregator_is_wired_into_the_build() -> None:
    mk = _src(MAKEFILE)
    assert "RiverHistConcatenate.o" in mk
    assert "river_hist_concatenate.x" in mk
    assert re.search(r"postprocess\.x :.*river_hist_concatenate\.x", mk)
    # and removed by the clean target, like its siblings
    assert "run/river_hist_concatenate.x" in mk


def test_reservoir_fields_are_not_treated_as_unit_catchment() -> None:
    """Reservoir shards are (reservoir_local, time) -- the same shape as a
    unit-catchment field. Matching on 'two dims, second is time' swept them up
    and rebuilt them against ucat_ucid/x_ucat/y_ucat, whose length is
    unitcat_local. The first dimension must be matched by name.
    """
    agg = _src(AGG)
    body = agg[agg.index("SUBROUTINE collect_varnames"):]
    body = body[: body.index("END SUBROUTINE collect_varnames")]
    assert "nf90_inq_dimid (nc, 'unitcat_local', udim)" in body
    assert "IF (dimids(1) /= udim) CYCLE" in body


def test_reservoir_fields_are_rebuilt_not_merely_excluded() -> None:
    """Excluding reservoir from the unit-catchment sweep stopped the wrong
    reconstruction but silently dropped the variables from the aggregate --
    worse than being visibly wrong, because nothing reports it."""
    agg = _src(AGG)
    assert "CALL rebuild_resv_variables ()" in agg
    body = agg[agg.index("SUBROUTINE rebuild_resv_variables"):]
    body = body[: body.index("END SUBROUTINE rebuild_resv_variables")]
    assert "resv_global_index" in body, "reservoirs must be keyed on their global index"
    assert "'reservoir'" in body
    for failure in ("reservoir id out of range", "duplicate reservoir id",
                    "incomplete reservoir coverage"):
        assert failure in body, failure


def test_e2e_covers_reservoir_survival() -> None:
    e2e = (ROOT / "tests/run_river_hist_concat_e2e.sh").read_text(encoding="utf-8")
    for name in ("volresv", "qresv_in", "qresv_out"):
        assert name in e2e, name
    assert "silently dropped" in e2e
