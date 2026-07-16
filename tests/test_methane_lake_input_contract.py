from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATE = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_State.F90").read_text()
INITIALIZE = (ROOT / "mkinidata/MOD_Initialize.F90").read_text()
AGGREGATION = (ROOT / "mksrfdata/Aggregation_LakeSoilC.F90").read_text()


def _subroutine(source: str, name: str) -> str:
    start = source.index(f"SUBROUTINE {name}")
    end = source.index(f"END SUBROUTINE {name}", start)
    return source[start:end]


def test_lulcc_remap_initializes_new_lakes_from_current_surface_data() -> None:
    remap_lulcc = _subroutine(STATE, "remap_methane_lulcc_state")
    no_snapshot = remap_lulcc.index("IF (.not. methane_lulcc_snapshot_valid) THEN")
    build_map = remap_lulcc.index("CALL build_lulcc_remap_map")
    assert "initialize_methane_lake_soilc_from_surface" in remap_lulcc
    assert "lake_soilc_srf" in remap_lulcc
    assert "initialize_lake_from_surface" in remap_lulcc
    assert no_snapshot < build_map


def test_enabled_lake_production_rejects_any_missing_lake_patch() -> None:
    initialize = _subroutine(STATE, "initialize_methane_lake_soilc_from_surface")
    assert "IF (lake_soilc_nlake == 0) RETURN" in initialize
    assert "lake_soilc_missing > 0" in initialize
    assert "CALL CoLM_stop" in initialize
    assert "WARNING: lake CH4 production is enabled" not in initialize


def test_single_point_lake_uses_existing_site_organic_matter_input() -> None:
    assert "lake_soilc_srf(:,ipatch)" in INITIALIZE
    assert "OM_density(:,ipatch)" in INITIALIZE


def test_gridded_lake_carbon_uses_area_weights_and_no_depth_proxy() -> None:
    assert AGGREGATION.index("IF (lake_global > 0) THEN") < AGGREGATION.index(
        "lake CH4 production requires raw lake_soilc.nc"
    )
    assert "no lake patches in domain; writing zero lake sediment carbon" in AGGREGATION
    assert "area = area_one" in AGGREGATION
    assert "median(lake_soilc_one" not in AGGREGATION
    assert "src_layer" not in AGGREGATION
    assert "vf_om_s_l" not in AGGREGATION
    assert "CALL CoLM_stop" in AGGREGATION
