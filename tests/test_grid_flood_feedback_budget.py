"""Small conservative geometry checks for the grid flood-credit ledger."""

from pathlib import Path


FLOW = Path(__file__).resolve().parents[1] / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90"


def exchange(overlap, patch_overlap, visible, protected, patch_ratios):
    donor_area = [sum(row) for row in overlap]
    route_area = [sum(row[g] for row in overlap) for g in range(len(overlap[0]))]
    covered_area = [sum(row[g] for row in patch_overlap) for g in range(len(route_area))]
    grid_area = [max(a, b) for a, b in zip(route_area, covered_area)]
    grid_credit = [
        sum((visible[d] + protected[d]) * overlap[d][g] / donor_area[d]
            for d in range(len(overlap)) if donor_area[d])
        for g in range(len(route_area))
    ]
    patch_credit = [
        sum(grid_credit[g] / grid_area[g] * row[g] for g in range(len(grid_area)) if grid_area[g])
        / sum(row)
        for row in patch_overlap
    ]
    grid_ratio = [
        sum(patch_ratios[p] * patch_overlap[p][g] for p in range(len(patch_overlap))) / grid_area[g]
        if grid_area[g] else 0
        for g in range(len(grid_area))
    ]
    donor_ratio = [
        sum(grid_ratio[g] * overlap[d][g] for g in range(len(grid_area))) / donor_area[d]
        if donor_area[d] else 0
        for d in range(len(overlap))
    ]
    return patch_credit, grid_credit, donor_ratio


def test_two_donors_two_grids_two_crossing_patches_close_by_source():
    overlap = [[10, 5], [10, 15]]
    patches = [[12, 8], [8, 12]]
    visible, protected = [20, 15], [10, 10]
    ratios = [0.2, 0.6]
    credit, _, donor_ratio = exchange(overlap, patches, visible, protected, ratios)
    used_by_patch = sum(credit[p] * sum(patches[p]) * ratios[p] for p in range(2))
    used_by_source = sum((visible[d] + protected[d]) * donor_ratio[d] for d in range(2))
    assert abs(used_by_patch - used_by_source) < 1e-12
    assert all(0 <= f <= 1 for f in donor_ratio)
    assert all(visible[d] * donor_ratio[d] <= visible[d] for d in range(2))
    assert all(protected[d] * donor_ratio[d] <= protected[d] for d in range(2))


def test_partial_route_coverage_and_dry_grid_cannot_overpublish():
    # 10 m2 route support, 20 m2 patch support: credit is 10 m3, not 20.
    credit, grid_credit, donor_ratio = exchange([[10, 0]], [[20, 0]], [10], [0], [1])
    assert credit == [0.5]
    assert credit[0] * 20 == 10
    assert donor_ratio == [1]
    # A patch crossing wet/dry grids cannot debit the dry grid.
    credit, grid_credit, donor_ratio = exchange([[10, 0]], [[10, 10]], [10], [0], [1])
    assert grid_credit == [10, 0]
    assert credit[0] * 20 == 10
    assert donor_ratio == [1]


def test_debit_is_before_routing_early_return_and_uses_adjoint_weights():
    step = FLOW.read_text().split("SUBROUTINE grid_riverlake_flow (", 1)[1]
    debit = step.index("CALL debit_flood_feedback(totalflood_evap, totalflood_infil)")
    early_return = step.index("IF (acctime_rnof+0.01 < acctime_rnof_max) THEN")
    assert debit < early_return
    helper = FLOW.read_text().split("SUBROUTINE debit_flood_feedback(", 1)[1]
    assert "filter=flood_credit_patch>0._r8" in helper
    assert "ratio_grid = ratio_grid / max(flood_grid_area" in helper
    assert "worker_push_data(push_inpm2ucat, ratio_grid, ratio_uc" in helper
