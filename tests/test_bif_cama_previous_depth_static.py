from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIF = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90").read_text(
    encoding="utf-8"
)
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text(
    encoding="utf-8"
)
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text(
    encoding="utf-8"
)


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_previous_channel_depth_is_restart_visible_and_legacy_safe() -> None:
    allocator = routine(TIMEVARS, "allocate_GridRiverLakeTimeVars")
    reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")
    writer = routine(TIMEVARS, "WRITE_GridRiverLakeTimeVars")

    assert "wdsrf_ucat_prev" in TIMEVARS
    assert "wdsrf_ucat_prev_valid" in TIMEVARS
    assert "allocate (wdsrf_ucat_prev (ncell_state))" in allocator
    assert "ncio_var_exist(file_restart, 'wdsrf_ucat_prev'" in reader
    assert "wdsrf_ucat_prev = wdsrf_ucat" in reader
    assert "wdsrf_ucat_prev_restart_found" in reader
    assert ".not. legacy_restart" in reader
    assert "transaction is missing wdsrf_ucat_prev" in reader
    assert "ieee_is_finite(wdsrf_ucat_prev(i))" in reader
    assert "invalid wdsrf_ucat_prev" in reader
    assert "'wdsrf_ucat_prev', 'ucatch'" in writer


def test_bif_depth_rule_selects_cama_levee_and_nonlevee_variants() -> None:
    calc = routine(BIF, "bifurcation_calc")

    assert "wdsrf_ucat_prev" in calc
    assert "wdsrf_prev_dn_pth" in calc
    assert "h_face_current = max(height_up, height_dn)" in calc
    assert "h_face_previous = max(height_up_prev, height_dn_prev)" in calc
    assert "h_face = sqrt(h_face_current * h_face_previous)" in calc
    assert "IF (h_face <= 0._r8) h_face = h_face_current" in calc
    assert "IF (DEF_USE_LEVEE) THEN" in calc
    assert "sqrt(h_face_current * 0.01_r8)" in calc
    assert "IF (h_face <= BIFMIN) THEN" in calc
    assert "IF (ilev == 1) THEN" in calc
    assert "h_face = h_face_current" in calc


def test_previous_depth_shares_dynamic_push_and_advances_after_flux_calculation() -> None:
    calc = routine(BIF, "bifurcation_calc")
    flow = FLOW.split("SUBROUTINE grid_riverlake_flow (year, deltime)", 1)[1].split(
        "END SUBROUTINE grid_riverlake_flow", 1
    )[0]

    assert "dynamic_state_fields(2)%send => wdsrf_ucat_prev" in calc
    assert "dynamic_state_fields(2)%recv => wdsrf_prev_dn_pth" in calc

    fallback = flow.index(".not. wdsrf_ucat_prev_valid")
    call = flow.index("CALL bifurcation_calc")
    advance = flow.index("wdsrf_ucat_prev = wdsrf_ucat", call)
    assert fallback < call < advance


def test_previous_depth_and_path_momentum_restore_as_one_state_unit() -> None:
    reader = routine(BIF, "read_bifurcation_restart")
    init = routine(FLOW, "grid_riverlake_flow_init")

    assert "previous_depth_restart_found" in reader
    assert "restart_loaded" in reader
    condition = (
        "has_pth_veloc .and. has_pth_momen .and. has_path_signature .and. &"
    )
    assert condition in reader
    assert "previous_depth_restart_found" in reader.split(condition, 1)[1]
    assert "CALL read_bifurcation_restart(gridriver_restart_file, &" in init
    assert "wdsrf_ucat_prev_restart_found" in init
    assert "IF (.not. bif_restart_loaded)" in init
    assert "wdsrf_ucat_prev = wdsrf_ucat" in init
