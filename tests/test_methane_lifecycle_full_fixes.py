from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def source(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def routine(text: str, name: str) -> str:
    return text.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_methane_restart_schema4_preserves_legacy_migration_boundaries() -> None:
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    state = source("main/TRACER/MOD_Tracer_Reactive_Methane_State.F90")
    accflux = source("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90")
    callback = routine(methane, "ch4_reactive_read_restart")
    validator = routine(methane, "validate_methane_restart_transaction")
    state_reader = routine(state, "read_methane_restart")
    acc_reader = routine(accflux, "read_methane_accflux_restart")

    assert "METHANE_RESTART_SCHEMA_VERSION = 4" in methane
    assert "real(METHANE_RESTART_SCHEMA_VERSION, r8)" in routine(
        methane, "ch4_reactive_write_restart"
    )
    assert "integer, intent(out) :: restart_schema" in validator
    normalized_validator = " ".join(validator.replace("&", " ").split())
    assert "restart_schema /= 1 .and. restart_schema /= 2 .and. restart_schema /= 3" in normalized_validator
    assert "any(schema == 3._r8)" in validator
    assert "restart_schema /= METHANE_RESTART_SCHEMA_VERSION" in validator
    assert "read_methane_restart (file_restart, strict_restart, restart_schema)" in callback
    normalized_callback = " ".join(callback.replace("&", "").split())
    assert (
        "read_methane_accflux_restart (file_restart, file_has_pools, "
        "file_has_microbe_accumulators, strict_restart, restart_schema)"
        in normalized_callback
    )

    state_migration = state_reader.split(
        "IF (strict_restart_active .and. restart_schema_active == 1) THEN", 1
    )[1]
    state_disable = state_migration.index("ncio_set_complete_require_present (.false.)")
    state_restore = state_migration.index("ncio_set_complete_require_present (.true.)")
    for field in ("ch4_lake_water_ch4_stock", "ch4_lake_water_o2_stock"):
        assert state_disable < state_migration.index(field) < state_restore

    phase_migration = state_reader.split(
        "IF (strict_restart_active .and. restart_schema_active < 4)", 1
    )[1]
    phase_disable = phase_migration.index("ncio_set_complete_require_present (.false.)")
    phase_restore = phase_migration.index("ncio_set_complete_require_present (.true.)")
    for field in (
        "ch4_lake_frozen_ch4_stock",
        "ch4_lake_frozen_o2_stock",
        "ch4_lake_liquid_fraction_prev",
    ):
        assert phase_disable < phase_migration.index(field) < phase_restore
    assert "any(lake_phase_fields_present) .and. .not. all(lake_phase_fields_present)" in state_reader

    # Schema 3 added the current accumulator transaction.  Schemas 1/2 may
    # default missing fields only while loading, then the whole partial window
    # is flushed so old sums cannot be mixed with new zero denominators.
    assert "strict_restart_active .and. restart_schema_active < 3" in acc_reader
    assert "ncio_set_complete_require_present (.false.)" in acc_reader
    assert "ncio_set_complete_require_present (.true.)" in acc_reader
    assert "IF (restart_schema_active < 3) THEN" in acc_reader


def test_restart_write_read_names_remain_symmetric() -> None:
    for relative, writer_name, reader_name in (
        (
            "main/TRACER/MOD_Tracer_Reactive_Methane_State.F90",
            "write_methane_restart",
            "read_methane_restart",
        ),
        (
            "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90",
            "write_methane_accflux_restart",
            "read_methane_accflux_restart",
        ),
    ):
        text = source(relative)
        writer = routine(text, writer_name)
        reader = routine(text, reader_name)
        written = set(re.findall(r"'((?:ch4_|CH4_)[^']+)'", writer))
        read = set(re.findall(r"'((?:ch4_|CH4_)[^']+)'", reader))
        assert written == read


def test_lake_frozen_phase_state_covers_restart_lulcc_and_history() -> None:
    state = source("main/TRACER/MOD_Tracer_Reactive_Methane_State.F90")
    physics = source("main/TRACER/MOD_Tracer_Reactive_Methane_Physics.F90")
    accflux = source("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90")

    for field in ("lake_frozen_ch4_stock", "lake_frozen_o2_stock"):
        assert f"allocate ({field}" in state
        assert f"ch4_{field}" in routine(state, "write_methane_restart")
        assert f"ch4_{field}" in routine(state, "read_methane_restart")
        assert f"CALL remap1d_mass(lulcc_{field}_old, {field})" in state
        assert f"deallocate ({field})" in state
    assert "ch4_lake_liquid_fraction_prev" in routine(state, "write_methane_restart")
    assert "CALL remap1d(lulcc_lake_liquid_fraction_prev_old, lake_liquid_fraction_prev)" in state
    assert "totcol_methane_lake + lake_water_ch4_stock + lake_frozen_ch4_stock" in physics
    assert "lake_water_ch4_stock + lake_frozen_ch4_stock" in accflux
    assert "lake_water_o2_stock + lake_frozen_o2_stock" in accflux


def test_schema1_migration_cold_starts_only_lakes_and_resets_mixed_window() -> None:
    state = source("main/TRACER/MOD_Tracer_Reactive_Methane_State.F90")
    accflux = source("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90")
    state_reader = routine(state, "read_methane_restart")
    acc_reader = routine(accflux, "read_methane_accflux_restart")

    cold_start = state_reader.split(
        "Missing water-column fields identify an old restart", 1
    )[1].split("WHERE (invalid_restart_value(totcol_methane))", 1)[0]
    assert "patchtype == PATCHTYPE_LAKE" in cold_start
    assert "PATCHTYPE_LAKE = 4" in state
    assert "patchtype == WATERBODY" not in state
    # WATERBODY remains the separate land-cover class used by LULCC remapping.
    assert "patchclass_new(np) == WATERBODY" in state
    assert "fsat_bef = spval" in cold_start

    migration_reset = acc_reader.split(
        "Schemas 1 and 2 predate part of the current accumulator transaction", 1
    )[1]
    assert "restart_schema_active < 3" in migration_reset
    assert "CALL flush_methane_acc_fluxes ()" in migration_reset


def test_initial_history_init_preserves_reactive_restart_accumulators() -> None:
    hist = source("main/MOD_Hist.F90")
    accflux = source("main/MOD_Vars_1DAccFluxes.F90")
    hist_init = routine(hist, "hist_init")
    flush = routine(accflux, "FLUSH_acc_fluxes")

    assert "flush_reactive = .false." in hist_init
    assert "IF (present(lulcc_call)) flush_reactive = lulcc_call" in hist_init
    assert "CALL FLUSH_acc_fluxes (flush_reactive=flush_reactive)" in hist_init
    assert "logical, intent(in), optional :: flush_reactive" in flush
    assert "IF (flush_reactive_active) CALL tracer_flush_acc_fluxes ()" in flush


def test_lulcc_reload_uses_current_land_year_path() -> None:
    reactive = source("main/TRACER/MOD_Tracer_Reactive.F90")
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    driver = source("main/LULCC/MOD_Lulcc_Driver.F90")
    dispatch = routine(reactive, "tracer_reactive_reload_lulcc_inputs")
    reload = routine(methane, "ch4_reactive_reload_lulcc_inputs")

    assert "SUBROUTINE reactive_reload_lulcc_if (lc_year, dir_landdata)" in reactive
    assert "CALL tracer_reactive_reload_lulcc_inputs (jdate(1), dir_landdata)" in driver
    assert "integer, intent(in) :: lc_year" in dispatch
    assert "character(len=*), intent(in) :: dir_landdata" in dispatch
    assert "%reload_lulcc (lc_year, dir_landdata)" in dispatch
    assert "write(cyear_lulcc,'(i4.4)') lc_year" in reload
    assert re.search(
        r"last_methane_ph_patch_file\s*=\s*trim\(dir_landdata\).*?"
        r"trim\(cyear_lulcc\).*?methane_ph_patches\.nc",
        reload,
        re.DOTALL,
    )
    assert reload.index("last_methane_ph_patch_file =") < reload.index(
        "CALL read_methane_ph_patch"
    )


def test_accflux_skips_disabled_history_and_prunes_core_selector() -> None:
    accflux = source("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90")
    body = routine(accflux, "accumulate_methane_fluxes")
    first_acc = body.index("CALL acc1d")

    assert body.index("history_mode = methane_history_accumulation_mode ()") < first_acc
    assert body.index("IF (history_mode == 0) RETURN") < first_acc
    assert body.index("core_history_only = history_mode == 1") < first_acc

    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    mode = const.split("FUNCTION methane_history_accumulation_mode", 1)[1].split(
        "END FUNCTION methane_history_accumulation_mode", 1
    )[0]
    assert ".not. DEF_METHANE%write_ch4_history" in mode
    assert "CASE ('none', 'off', 'false', '.false.')" in mode
    assert "CASE ('core', 'default', 'minimal', 'fast')" in mode

    core = body.split("IF (core_history_only) THEN", 1)[1].split(
        "ENDIF", 1
    )[0]
    required = (
        "methane_surf_flux_tot",
        "methane_surf_flux_tot_phys",
        "methane_balance_residual",
        "methane_ch4_clip_credit",
        "o2_cap_loss",
        "o2_cap_gain",
        "methane_prod_tot",
        "methane_oxid_tot",
        "totcol_methane",
        "methane_surf_flux_wetland",
        "methane_surf_flux_soil",
        "methane_surf_flux_lake",
        "methane_surf_flux_rice",
        "methane_surf_flux_tot_lake",
        "a_methane_acc_num",
        "a_methane_acc_num_lake",
    )
    for name in required:
        assert name in core
    assert "RETURN" in core
    assert "methane_prod_depth" not in core
