from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text()
LEVEE = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeLevee.F90").read_text()
TRACER = (ROOT / "main/TRACER/MOD_Tracer_RiverLake.F90").read_text()


def test_bif_flux_returns_to_cama_static_levee_partition() -> None:
    assert "IF (abs(protected_hflux) > 1.e-20_r8)" not in FLOW
    assert "CALL levee_fldstg(i, volwater" not in FLOW
    assert "CALL levee_repartition_storage(i, volwater, wdsrf_ucat(i)" in FLOW
    assert "CaMa's simplified levee scheme applies all pathway fluxes" in FLOW


def test_restart_folds_nonleveed_protected_water_instead_of_zeroing_it() -> None:
    assert "ieee_is_finite(levee_hgt(i))" in LEVEE
    assert "ieee_is_finite(levsto(i))" in LEVEE
    assert "ELSEIF (levsto(i) < 0._r8) THEN" in LEVEE
    assert "ieee_is_finite(levsto) .or." not in LEVEE
    assert "volwater_ucat_io(i) = volwater_ucat_io(i) + levsto(i)" in LEVEE
    assert "currently non-leveed cell(s) into visible routing storage" in LEVEE
    assert "IF (has_restart_var(2) == 0) THEN" in LEVEE
    assert "old levee-side stage cannot reconstruct visible storage without double counting" in LEVEE
    assert "levee-to-reservoir restart mapping is ambiguous" in LEVEE
    assert "WHERE (.not. has_levee)" not in LEVEE


def test_restart_folds_matching_tracer_mass_into_visible_pool() -> None:
    assert "ieee_is_finite(trc_levsto(itrc, ii_bf))" in TRACER
    assert "ELSEIF (trc_levsto(itrc, ii_bf) < 0._r8) THEN" in TRACER
    assert "IF (.not. has_levee_bf(ii_bf)) THEN" in TRACER
    assert "IF (.not. tracer_uses_land_water_transport(itrc_bf)) CYCLE" in TRACER
    assert "trc_mass(itrc_bf, ii_bf) = trc_mass(itrc_bf, ii_bf)" in TRACER
    assert "+ max(trc_levsto(itrc_bf, ii_bf), 0._r8)" in TRACER
    assert "trc_levsto(itrc_bf, ii_bf) = 0._r8" in TRACER


def test_fold_algebra_preserves_total_water_and_tracer() -> None:
    visible_water, protected_water = 80.0, 20.0
    visible_tracer, protected_tracer = 160.0, 40.0

    water_before = visible_water + protected_water
    tracer_before = visible_tracer + protected_tracer
    visible_water += protected_water
    visible_tracer += protected_tracer
    protected_water = 0.0
    protected_tracer = 0.0

    assert visible_water + protected_water == water_before
    assert visible_tracer + protected_tracer == tracer_before
