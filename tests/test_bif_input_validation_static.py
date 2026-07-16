from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NETWORK = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeNetwork.F90").read_text()


def test_bif_parameter_shapes_and_finite_values_are_validated() -> None:
    assert "bifurcation parameter dimensions are inconsistent" in NETWORK
    assert "ieee_is_finite(bif_dist_all)" in NETWORK
    assert "ieee_is_finite(bif_elev_all)" in NETWORK
    assert "ieee_is_finite(bif_wdth_all)" in NETWORK
    assert "ieee_is_finite(bif_mann_all)" in NETWORK


def test_bif_active_geometry_requires_positive_distance_width_and_manning() -> None:
    assert "any(bif_dist_all <= 0._r8)" in NETWORK
    finite_distance = NETWORK.index("any(.not. ieee_is_finite(bif_dist_all))")
    positive_distance = NETWORK.index("any(bif_dist_all <= 0._r8)")
    assert finite_distance < positive_distance
    assert "ieee_is_finite(bif_dist_all)) .or." not in NETWORK
    assert "active bifurcation layer requires a positive Manning coefficient" in NETWORK
    assert "bifurcation pathway has no active positive-width layer" in NETWORK
