"""Step-3 checks: route history obeys DEF_HIST_mode through a single dispatcher.

The plan's requirement is not merely that block mode works, but that adding a
new route variable cannot wire it into one mode and silently miss the other:

    不维护两份手工变量清单。每个河网变量只保留一个写调用。

So the property under test is structural -- every route-history write in
MOD_Grid_RiverLakeHist goes through route_hist_write_*, and the mode branch
exists in exactly one place.
"""

from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
HIST = ROOT / "main/HYDRO/MOD_Grid_RiverLakeHist.F90"
ROUTE = ROOT / "main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90"
MAKEFILE = ROOT / "Makefile"

LOW_LEVEL_WRITERS = (
    "vector_gather_map2grid_and_write",
    "vector_gather_and_write",
    "vector_gather_matrix_to_master",
)


def _hist_out() -> str:
    src = HIST.read_text(encoding="utf-8")
    start = src.index("SUBROUTINE hist_grid_riverlake_out")
    end = src.index("END SUBROUTINE hist_grid_riverlake_out", start)
    return src[start:end]


def _strip_comments(text: str) -> str:
    return "\n".join(re.sub(r"!.*$", "", line) for line in text.split("\n"))


def test_no_route_variable_calls_a_low_level_writer_directly() -> None:
    body = _strip_comments(_hist_out())
    for writer in LOW_LEVEL_WRITERS:
        assert writer not in body, (
            f"hist_grid_riverlake_out still calls {writer} directly; every route "
            "variable must go through the route_hist_write_* dispatcher"
        )


def test_every_route_variable_goes_through_the_dispatcher() -> None:
    body = _strip_comments(_hist_out())
    n = sum(body.count(f"CALL route_hist_write_{k}") for k in ("ucat", "resv", "bif_matrix"))
    # 15 unit-catchment fields + 3 reservoir fields + the BIF matrix
    assert n == 19, f"expected 19 dispatched route writes, found {n}"


def test_mode_branch_lives_in_exactly_one_module() -> None:
    """DEF_HIST_mode must not be tested in the history routine itself."""
    assert "DEF_HIST_mode" not in _strip_comments(_hist_out())
    # the branch itself, not the name appearing in an error string
    route = _strip_comments(ROUTE.read_text(encoding="utf-8"))
    assert route.count("trim(DEF_HIST_mode) ==") == 1


def test_file_creation_and_time_axis_are_shared_by_both_modes() -> None:
    """Both skeletons live in route_hist_begin, so they cannot drift."""
    body = _strip_comments(_hist_out())
    for leaked in ("ncio_create_file", "ncio_write_time"):
        assert leaked not in body, f"{leaked} escaped route_hist_begin"
    route = ROUTE.read_text(encoding="utf-8")
    assert route.count("ncio_write_time") == 2      # one per mode
    assert route.count("ncio_create_file") == 2


def test_unmigrated_writers_fail_loudly_under_block_mode() -> None:
    """A shard set missing tracer variables would still look complete."""
    body = _hist_out()
    assert body.count("CALL route_hist_require_migrated") == 2
    for what in ("tracer_lifecycle_route_write_history", "write_tracer_history"):
        assert f"route_hist_require_migrated ('{what}')" in body
    guard = ROUTE.read_text(encoding="utf-8")
    guard_body = guard[guard.index("SUBROUTINE route_hist_require_migrated"):]
    guard_body = guard_body[: guard_body.index("END SUBROUTINE route_hist_require_migrated")]
    assert "CoLM_stop" in guard_body
    assert "IF (rh_block)" in guard_body


def test_dispatcher_layering_keeps_the_network_out_of_the_generic_utility() -> None:
    vrw = (ROOT / "main/HYDRO/MOD_Vector_ReadWrite.F90").read_text(encoding="utf-8")
    assert "MOD_Grid_RiverLakeNetwork" not in vrw
    assert "MOD_Grid_Reservoir" not in vrw
    route = ROUTE.read_text(encoding="utf-8")
    assert "USE MOD_Grid_RiverLakeNetwork" in route
    assert "USE MOD_Vector_ReadWrite" in route
    # and the dispatcher must not depend back on the history module
    assert "USE MOD_Grid_RiverLakeHist" not in route.replace(
        "MOD_Grid_RiverLakeHistRoute", "")


def test_dispatcher_is_wired_into_the_build() -> None:
    mk = MAKEFILE.read_text(encoding="utf-8")
    assert "MOD_Grid_RiverLakeHistRoute.o" in mk
    assert re.search(
        r"MOD_Grid_RiverLakeHistRoute\.o:.*MOD_Vector_ReadWrite\.o", mk)
    assert re.search(
        r"MOD_Grid_RiverLakeHist\.o:.*MOD_Grid_RiverLakeHistRoute\.o", mk)


def test_filename_is_derived_on_every_rank() -> None:
    """IO ranks need it to build their shard name, not just the master."""
    body = _hist_out()
    derive = body.index("file_hist_ucat = file_hist(1:i)")
    begin = body.index("CALL route_hist_begin")
    assert derive < begin
    # the derivation must not sit inside a master-only branch
    preceding = body[:derive]
    assert preceding.rstrip().rsplit("\n", 1)[-1].strip() != "IF (p_is_master) THEN"


def test_flush_ordering_is_unchanged() -> None:
    """Each history window must still be zeroed exactly once, after the write."""
    body = _hist_out()
    assert body.count("CALL flush_acc_fluxes_riverlake ()") == 1
    assert body.index("CALL route_hist_end ()") < body.index(
        "CALL flush_acc_fluxes_riverlake ()")
