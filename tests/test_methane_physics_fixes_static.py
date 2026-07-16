from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def source(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def test_soil_ice_is_excluded_from_mobile_ch4_storage() -> None:
    physics = source("MOD_Tracer_Reactive_Methane_Physics.F90")
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    split = physics.split("subroutine split_ch4_o2_phases", 1)[1].split(
        "end subroutine split_ch4_o2_phases", 1
    )[0]

    assert "vol_ch4_storage" in split
    assert "vol_ch4_storage(j) = vol_aqu(j)" in split
    assert "methane_frzout is retired" in const


def test_fully_frozen_ch4_storage_is_not_zeroed_by_transport() -> None:
    physics = source("MOD_Tracer_Reactive_Methane_Physics.F90")
    transport = physics.split("subroutine methane_tran", 1)[1].split(
        "end subroutine methane_tran", 1
    )[0]

    assert "vol_ch4_storage" in transport
    assert "conc_ch4_rel(j) = conc_methane(j)/epsilon_t(j,1)" in transport
    assert "max(vol_ch4_storage(j), smallnumber)" in transport


def test_slpratio_has_one_ratio_contract_without_value_based_unit_guessing() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    compute = state.split("SUBROUTINE compute_f_h2osfc", 1)[1].split(
        "END SUBROUTINE compute_f_h2osfc", 1
    )[0]

    assert "slope_arg = max(slpratio_in, 0._r8)" in compute
    assert "slope_angle = atan(slope_arg)" in compute
    assert "ELSEIF (slope_arg" not in compute
    assert "slope_arg / 100._r8" not in compute
    assert "slope_arg * PI / 180._r8" not in compute


def test_retired_hydrology_parameters_are_inert_and_fail_closed() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    standard = (ROOT / "run" / "standard_ch4_parameter.nml").read_text(
        encoding="utf-8"
    )

    assert "vdcf and pc remain in the namelist only for backward compatibility" in const
    assert "retired methane hydrology knobs vdcf/pc must retain defaults" in const
    assert "rice_paddy_min_wdsrf_mm" not in const
    assert "rice_paddy_min_wdsrf_mm" not in standard
    assert "DEF_METHANE_hydrology%pc" not in state


def test_uncoupled_microbial_pools_and_dependent_options_fail_fast() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    validation = const.split("SUBROUTINE validate_methane_namelist", 1)[1].split(
        "END SUBROUTINE validate_methane_namelist", 1
    )[0]
    standard = (ROOT / "run" / "standard_ch4_parameter.nml").read_text(
        encoding="utf-8"
    )

    assert "IF (DEF_METHANE%use_microbial_pools) THEN" in validation
    assert "microbial pools are disabled until biomass growth/loss has donor/sink carbon coupling" in validation
    assert "use_microbial_flux_override .and." in validation
    assert "use_microbial_dormancy .and." in validation
    assert "ieee_is_finite(DEF_METHANE%B_max_fraction_methanogen)" in validation
    assert "ieee_is_finite(DEF_METHANE%B_max_fraction_methanotroph)" in validation
    assert "DEF_METHANE%use_microbial_pools         = .false." in standard
    assert "DEF_METHANE%use_microbial_flux_override = .false." in standard
    assert "DEF_METHANE%use_microbial_dormancy      = .false." in standard
    assert "B_max_fraction_methanogen  = -1.0" in standard
    assert "B_max_fraction_methanotroph = -1.0" in standard


def test_lulcc_remap_builds_and_reuses_one_patch_mapping() -> None:
    state = source("MOD_Tracer_Reactive_Methane_State.F90")
    remap = state.split("SUBROUTINE remap_methane_lulcc_state", 1)[1].split(
        "END SUBROUTINE remap_methane_lulcc_state", 1
    )[0]

    assert remap.index("CALL build_lulcc_remap_map") < remap.index("CALL remap2d")
    assert "map_start" in remap
    assert "map_old" in remap
    assert "map_source_weight" in remap
    assert "map_mass_weight" in remap
    # The new/old Cartesian scan is paid once while building the mapping,
    # rather than once for every remapped state field.
    assert remap.count("IF (eindex_old(op) /= eindex_new(np)) CYCLE") == 1
    # Link weights must use precomputed class totals.  Calling any of these
    # historical helpers per link hides another old/new scan inside the map
    # construction and restores cubic-like work for dense patch groups.
    assert "old_group_class_area" in remap
    assert "target_group_class_area" in remap
    assert "lulcc_source_weight" not in remap
    assert "lulcc_mass_transfer_area" not in remap
    assert "lulcc_target_class_area" not in remap
    builder = remap.split("SUBROUTINE build_lulcc_remap_map", 1)[1].split(
        "END SUBROUTINE build_lulcc_remap_map", 1
    )[0]
    assert builder.count("DO op = 1, nold") == 2
    assert "DO oq" not in builder
    assert "DO nq" not in builder


def test_microbe_lulcc_fields_reuse_one_precomputed_patch_mapping() -> None:
    microbes = source("MOD_Tracer_Reactive_Methane_Microbes.F90")
    remap = microbes.split(
        "SUBROUTINE remap_methane_microbes_lulcc_state", 1
    )[1].split("END SUBROUTINE remap_methane_microbes_lulcc_state", 1)[0]

    assert remap.index("CALL build_lulcc_remap_map") < remap.index("CALL remap2d")
    for name in (
        "map_start",
        "map_old",
        "map_source_weight",
        "map_mass_weight",
        "fallback_map",
        "old_group_class_area",
        "target_group_class_area",
    ):
        assert name in remap
    assert "lulcc_source_weight" not in remap
    assert "lulcc_mass_transfer_area" not in remap
    assert "lulcc_target_class_area" not in remap
    builder = remap.split("SUBROUTINE build_lulcc_remap_map", 1)[1].split(
        "END SUBROUTINE build_lulcc_remap_map", 1
    )[0]
    assert builder.count("DO op = 1, nold") == 2
    assert "DO oq" not in builder
    assert "DO nq" not in builder
    remap2d = remap.split("SUBROUTINE remap2d", 1)[1].split(
        "END SUBROUTINE remap2d", 1
    )[0]
    remap3d = remap.split("SUBROUTINE remap_component3d", 1)[1].split(
        "END SUBROUTINE remap_component3d", 1
    )[0]
    assert "DO link = map_start" in remap2d
    assert "DO link = map_start" in remap3d
    assert "DO op = 1, min(size(old" not in remap2d
    assert "DO op = 1, min(size(old" not in remap3d


def test_precomputed_lulcc_weight_formula_matches_legacy_and_conserves() -> None:
    old_eindex = [10, 10, 10, 20]
    old_class = [0, 0, 1, 0]
    old_area = [20.0, 30.0, 50.0, 40.0]
    old_value = [2.0, 4.0, 8.0, 3.0]
    new_eindex = [10, 10, 20, 30]
    new_area = [100.0, 50.0, 80.0, 60.0]
    lcc = [[0.6, 0.4], [0.2, 0.8], [1.0, 0.0], [0.5, 0.5]]

    old_group_class_area = {}
    for eindex, cls, area in zip(old_eindex, old_class, old_area):
        old_group_class_area[eindex, cls] = (
            old_group_class_area.get((eindex, cls), 0.0) + area
        )

    target = {}
    target_group_class_area = {}
    for np, (eindex, area, fractions) in enumerate(
        zip(new_eindex, new_area, lcc)
    ):
        class_sum = sum(max(value, 0.0) for value in fractions)
        for cls, fraction in enumerate(fractions):
            value = area * max(fraction, 0.0) / class_sum
            target[np, cls] = value
            target_group_class_area[eindex, cls] = (
                target_group_class_area.get((eindex, cls), 0.0) + value
            )

    new_mass = [0.0] * len(new_eindex)
    for np, eindex in enumerate(new_eindex):
        for op, (old_index, cls, area, value) in enumerate(
            zip(old_eindex, old_class, old_area, old_value)
        ):
            if old_index != eindex:
                continue
            legacy_class_area = sum(
                source_area
                for source_index, source_class, source_area in zip(
                    old_eindex, old_class, old_area
                )
                if source_index == eindex and source_class == cls
            )
            legacy_source = lcc[np][cls] * area / legacy_class_area
            cached_source = lcc[np][cls] * area / old_group_class_area[eindex, cls]
            assert cached_source == pytest.approx(legacy_source)

            legacy_target = new_area[np] * lcc[np][cls] / sum(lcc[np])
            legacy_target_group = sum(
                candidate_area * lcc[nq][cls] / sum(lcc[nq])
                for nq, (candidate_index, candidate_area) in enumerate(
                    zip(new_eindex, new_area)
                )
                if candidate_index == eindex
            )
            legacy_mass = area * legacy_target / legacy_target_group
            cached_mass = (
                area * target[np, cls] / target_group_class_area[eindex, cls]
            )
            assert cached_mass == pytest.approx(legacy_mass)
            new_mass[np] += cached_mass * value

    for eindex in (10, 20):
        old_total = sum(
            area * value
            for index, area, value in zip(old_eindex, old_area, old_value)
            if index == eindex
        )
        remapped_total = sum(
            mass for index, mass in zip(new_eindex, new_mass) if index == eindex
        )
        assert remapped_total == pytest.approx(old_total)
    # A new eindex with no old links remains a zero-link/fallback case.
    assert new_mass[3] == 0.0


def test_reentry_resets_transient_atmospheric_methane_cache() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    reader = const.split("SUBROUTINE read_methane_namelist", 1)[1].split(
        "END SUBROUTINE read_methane_namelist", 1
    )[0]

    assert "DEF_METHANE = Methane_type()" in reader
    assert "DEF_METHANE_hydrology = Methane_hydrology_type()" in reader
    assert "atm_ch4_file_loaded = .false." in reader
    assert "atm_ch4_file_molmol(:,:) = -1._r8" in reader
    assert "atm_ch4_file_warned" not in const


def test_offline_only_contract_still_fails_fast() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    validation = const.split("SUBROUTINE validate_methane_namelist", 1)[1].split(
        "END SUBROUTINE validate_methane_namelist", 1
    )[0]

    assert "IF (.not. DEF_METHANE%methane_offline) THEN" in validation
    assert "online atmosphere/NEE coupling that is not implemented" in validation


def test_default_history_exposes_global_physical_and_residual_fluxes() -> None:
    const = source("MOD_Tracer_Reactive_Methane_Const.F90")
    core = const.split("logical FUNCTION methane_history_is_core", 1)[1].split(
        "END FUNCTION methane_history_is_core", 1
    )[0]

    assert "'f_methane_surf_flux_global_phys_with_lake'" in core
    assert "'f_methane_balance_residual_global_with_lake'" in core


def test_standard_comments_match_enabled_soil_and_prognostic_lake_scope() -> None:
    standard = (ROOT / "run" / "standard_ch4_parameter.nml").read_text(
        encoding="utf-8"
    )

    assert "routing-connected soil columns" in standard
    assert "DEF_METHANE%only_wetland      = .false." in standard
    assert "explicit well-mixed water-column stock" in standard
    assert "DEF_METHANE%allowlakeprod         = .true." in standard
    assert "DEF_METHANE%ch4_history_vars  = 'core'" in standard
