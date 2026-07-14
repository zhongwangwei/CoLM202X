from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TIMEVARS = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90").read_text()
BIF = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90").read_text()
LEVEE = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeLevee.F90").read_text()


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_schema_v2_writer_persists_explicit_feature_bits() -> None:
    writer = routine(TIMEVARS, "WRITE_GridRiverLakeTimeVars")

    assert "GRIDRIVER_RESTART_SCHEMA_VERSION = 2" in TIMEVARS
    assert "'gridriver_restart_feature_bifurcation'" in writer
    assert "merge(1, 0, DEF_USE_BIFURCATION)" in writer
    assert "'gridriver_restart_feature_levee'" in writer
    assert "merge(1, 0, DEF_USE_LEVEE)" in writer


def test_reader_requires_manifest_only_for_schema_v2() -> None:
    validator = routine(TIMEVARS, "validate_gridriver_restart_transaction")

    assert "restart_schema == 1" in validator
    assert "restart_schema == GRIDRIVER_RESTART_SCHEMA_VERSION" in validator
    assert "GridRiverLake schema-v2 restart is missing feature manifest" in validator
    assert "invalid GridRiverLake restart feature manifest" in validator


def test_bif_manifest_controls_strict_restore_vs_legacy_cold_start() -> None:
    reader = routine(BIF, "read_bifurcation_restart")
    time_reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")
    time_writer = routine(TIMEVARS, "WRITE_GridRiverLakeTimeVars")

    assert "restart_transaction_validated_in" in reader
    assert "restart_feature_manifest_present_in" in reader
    assert "restart_bifurcation_enabled_in" in reader
    assert "requires validated transaction context" in reader
    assert reader.index("requires validated transaction context") < reader.index(
        "totalnpthout <= 0"
    )
    assert "gridriver_restart_schema" not in reader
    assert "gridriver_restart_feature_bifurcation" not in reader
    assert "strict_bif_restart" in reader
    assert "declares bifurcation enabled but pathway state is incomplete" in reader
    assert "declares bifurcation enabled but pathway state is invalid" in reader
    assert "restart declares bifurcation disabled; cold-starting pathway state" in reader
    assert "incomplete/legacy bifurcation restart" in reader
    assert "totalnpthout > 0 .and. npthlev_bif > 0" in time_reader
    assert "IF (restart_feature_manifest_present .and. .not. restore_previous_depth) has_var = .false." in time_reader
    assert "IF (DEF_USE_BIFURCATION .and. totalnpthout > 0 .and. npthlev_bif > 0) THEN" in time_writer


def test_levee_manifest_controls_required_levsto_payload() -> None:
    reader = routine(LEVEE, "read_levee_restart")

    assert "restart_transaction_validated_in" in reader
    assert "restart_feature_manifest_present_in" in reader
    assert "restart_levee_enabled_in" in reader
    assert "requires validated transaction context" in reader
    assert reader.index("requires validated transaction context") < reader.index(
        "has_restart_var = 0"
    )
    assert "gridriver_restart_schema" not in reader
    assert "gridriver_restart_feature_levee" not in reader
    assert "strict_levee_restart" in reader
    assert "declares levee enabled but levsto is missing" in reader
    assert "strict_levee_restart .and. totalnumucat > 0" in reader
    assert "IF (.not. restart_levee_enabled_in) has_restart_var(1) = 0" in reader


def test_timevars_owns_validated_restart_context_for_deferred_readers() -> None:
    validator = routine(TIMEVARS, "validate_gridriver_restart_transaction")

    assert "restart_transaction_validated = .false." in validator
    assert validator.count("restart_transaction_validated = .true.") == 3
    assert "legacy_restart" in validator
    assert "restart_schema == 1" in validator


def test_flow_passes_the_validated_context_to_both_deferred_readers() -> None:
    flow = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text()
    init = routine(flow, "grid_riverlake_flow_init")

    for name in (
        "restart_transaction_validated",
        "restart_feature_manifest_present",
        "restart_levee_enabled",
        "restart_bifurcation_enabled",
    ):
        assert name in init


def test_schema_v2_base_payload_is_complete_and_finite_first() -> None:
    reader = routine(TIMEVARS, "READ_GridRiverLakeTimeVars")

    for field in (
        "wdsrf_ucat",
        "veloc_riv",
        "acctime_rnof",
        "acc_rnof_uc",
        "volwater_ucat",
        "volresv",
    ):
        assert f"'{field}'" in reader

    assert "schema-v2 restart is missing required base state" in reader
    assert "ieee_is_finite(wdsrf_ucat(i))" in reader
    assert "ELSEIF (wdsrf_ucat(i) < 0._r8)" in reader
    assert "ieee_is_finite(veloc_riv(i))" in reader
    assert "ELSEIF (abs(veloc_riv(i)) > 50._r8)" in reader
    assert "ieee_is_finite(acctime_rnof)" in reader
    assert "ieee_is_finite(acc_rnof_uc(i))" in reader
    assert "ieee_is_finite(volwater_ucat(i))" in reader
    assert "ELSEIF (volwater_ucat(i) < 0._r8)" in reader
    assert "ieee_is_finite(volresv(i))" in reader
    assert "ELSEIF (volresv(i) < 0._r8)" in reader
    assert "invalid GridRiverLake restart base state" in reader
