from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEDIMENT = (ROOT / "main/TRACER/MOD_Tracer_Particle_Sediment.F90").read_text()
LIFECYCLE = (ROOT / "main/TRACER/MOD_Tracer_Lifecycle.F90").read_text()


def test_hydro_handoff_keeps_gross_directions_until_morphology() -> None:
    bif = SEDIMENT.split("SUBROUTINE sediment_bif_accumulate(", 1)[1].split(
        "END SUBROUTINE sediment_bif_accumulate", 1
    )[0]
    levee = SEDIMENT.split("SUBROUTINE sediment_levee_repartition(", 1)[1].split(
        "END SUBROUTINE sediment_levee_repartition", 1
    )[0]
    calc = SEDIMENT.split("SUBROUTINE grid_sediment_calc(", 1)[1].split(
        "END SUBROUTINE grid_sediment_calc", 1
    )[0]

    assert "max(bif_hflux_lev(:,ipth), 0._r8) * dt" in bif
    assert "max(-bif_hflux_lev(:,ipth), 0._r8) * dt" in bif
    assert "max(transfer, 0._r8)" in levee
    assert "max(-transfer, 0._r8)" in levee
    assert "sed_acc_bif_forward = 0._r8" in calc
    assert "sed_acc_to_protected = 0._r8" in calc
    assert "calc_sediment_advection" not in bif + levee
    assert "route_sediment_bif    => sediment_bif_accumulate" in SEDIMENT
    assert "route_sediment_levee  => sediment_levee_repartition" in SEDIMENT
    assert "tracer_lifecycle_route_sediment_bif_accumulate" in LIFECYCLE
    assert "tracer_lifecycle_route_sediment_levee_repartition" in LIFECYCLE


def test_combined_modes_are_not_rejected_before_parameter_validation() -> None:
    validator = SEDIMENT.split("SUBROUTINE validate_sediment_parameters()", 1)[1].split(
        "END SUBROUTINE validate_sediment_parameters", 1
    )[0]
    assert "sediment bifurcation transport is not yet implemented" not in validator
    assert "sediment levee transport is not yet implemented" not in validator
    assert "IF (nsed <= 0) THEN" in validator


def test_levee_restart_requires_separate_class_stocks_and_carrier_history() -> None:
    reader = SEDIMENT.split("SUBROUTINE read_sediment_restart(", 1)[1].split(
        "END SUBROUTINE read_sediment_restart", 1
    )[0]
    writer = SEDIMENT.split("SUBROUTINE write_sediment_restart(", 1)[1].split(
        "END SUBROUTINE write_sediment_restart", 1
    )[0]
    assert "(DEF_USE_LEVEE .or. DEF_USE_BIFURCATION) .and." in reader
    assert "schema_version /= SED_RESTART_SCHEMA_VERSION" in reader
    for name in (
        "sedsto_protected_", "sedbed_protected_", "sed_acc_protected_start",
        "sed_acc_pre_repartition_start", "sed_acc_protected_end",
        "sed_acc_to_protected", "sed_acc_from_protected",
    ):
        assert name in reader
        assert name in writer
    assert "merge(SED_RESTART_SCHEMA_VERSION, 2, DEF_USE_LEVEE .or. DEF_USE_BIFURCATION)" in writer


def test_combined_mode_cold_start_requires_absent_river_tracer_commit() -> None:
    reader = SEDIMENT.split("SUBROUTINE read_sediment_restart(", 1)[1].split(
        "END SUBROUTINE read_sediment_restart", 1
    )[0]
    assert "trc_river_restart_complete" in reader
    assert "(DEF_USE_LEVEE .or. DEF_USE_BIFURCATION) .and. has_river_tracer_commit" in reader
    assert "schema_version /= SED_RESTART_SCHEMA_VERSION" in reader


def test_bif_path_restart_and_shared_donor_are_morphology_only() -> None:
    prepare = SEDIMENT.split("SUBROUTINE prepare_bif_sediment(", 1)[1].split(
        "END SUBROUTINE prepare_bif_sediment", 1
    )[0]
    calc = SEDIMENT.split("SUBROUTINE grid_sediment_calc(", 1)[1].split(
        "END SUBROUTINE grid_sediment_calc", 1
    )[0]
    assert calc.index("CALL prepare_bif_sediment") < calc.index("CALL calc_sediment_advection")
    assert calc.index("CALL apply_bif_sediment_credits") > calc.index("CALL calc_sediment_advection")
    assert "push_bif_influx" in prepare and "push_bif_dn2pth" in prepare
    assert "main_sed_demand(ised,:) + levee_to_demand(ised,:)" in prepare
    assert "joint_sediment_scale(available_visible(i), recv_visible(i))" in prepare
    assert "joint_sediment_scale(available_protected(i), recv_protected(i))" in prepare
    for name in (
        "sed_acc_bif_forward", "sed_acc_bif_reverse",
        "sed_acc_bif_forward_time", "sed_acc_bif_reverse_time",
    ):
        assert f"{name} = 0._r8" in SEDIMENT.split("SUBROUTINE read_sediment_restart(", 1)[1]
        assert f"any({name} /= 0._r8)" in SEDIMENT.split("SUBROUTINE write_sediment_restart(", 1)[1]
