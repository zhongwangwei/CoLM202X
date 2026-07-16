from pathlib import Path
import math
import re


ROOT = Path(__file__).resolve().parents[1]


def source(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def routine(text: str, name: str) -> str:
    return text.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def function(text: str, name: str) -> str:
    return text.split(f"FUNCTION {name}", 1)[1].split(
        f"END FUNCTION {name}", 1
    )[0]


def test_routed_floodplain_redoxlag_uses_floodplain_parameter_and_keeps_rice_blend() -> None:
    driver = source("main/TRACER/MOD_Tracer_Reactive_Methane_Driver.F90")
    bgc = source("main/TRACER/MOD_Tracer_Reactive_Methane_BgcLink.F90")
    redox = function(bgc, "get_biome_redoxlag")

    call = driver.split(
        "biome_redoxlag_patch(i) = get_biome_redoxlag", 1
    )[1].split("ENDIF", 1)[0]
    floodplain_assignment = driver.index("is_floodplain_active = .false.")
    biome_f_guard = driver.index("IF (allocated(biome_f_methane_patch)")
    assert floodplain_assignment < biome_f_guard
    assert "is_floodplain_active" in call
    assert "is_floodplain" in redox
    assert re.search(
        r"IF \(floodplain_active\) THEN\s+"
        r"nonrice_redoxlag\s*=\s*DEF_METHANE%redoxlag_tropical_floodplain",
        redox,
    )
    assert "(1._r8 - rice_frac) * nonrice_redoxlag" in redox
    assert "rice_frac * DEF_METHANE%redoxlag_rice_paddy" in redox


def test_spatial_ph_requires_complete_vector_and_rejects_invalid_values() -> None:
    ph = source("main/TRACER/MOD_Tracer_Reactive_Methane_pH.F90")
    read = routine(ph, "read_methane_ph_patch")

    assert "ncio_read_vector => ncio_read_vector_complete" in read
    call = read.split("CALL ncio_read_vector", 1)[1].split(")", 1)[0]
    assert "defval" not in call
    assert "ieee_is_finite" in read
    assert "n_invalid" in read
    assert "CALL CoLM_Stop" in read
    assert "methane_ph_patch = methane_ph_fallback" not in read
    assert "abs(methane_ph_patch - methane_ph_fallback)" not in read
    assert "ph_active = .true." in read


def test_mksrfdata_ph_is_sparse_area_weighted_and_fail_closed() -> None:
    aggregate = source("mksrfdata/Aggregation_MethanePH.F90")

    assert aggregate.index("IF (relevant_global > 0) THEN") < aggregate.index(
        "required PHH2O raw file not found"
    )
    assert "no soil or wetland patches in domain; writing fallback pH only" in aggregate

    missing_file = aggregate.split("IF (.not. raw_exists) THEN", 1)[1].split(
        "ENDIF", 1
    )[0]
    assert "CALL CoLM_Stop" in missing_file
    assert "writing fallback" not in missing_file

    for message in (
        "PHH2O open failed",
        "PHH2O missing or invalid lat/lon/depth dimensions",
        "PHH2O latitude read failed",
        "PHH2O longitude read failed",
        "PHH2O variable read failed",
    ):
        error_branch = aggregate.split(message, 1)[1][:300]
        assert "CALL CoLM_Stop" in error_branch

    assert "logical, allocatable :: methane_ph_valid(:)" in aggregate
    assert "patch_valid(iset) = .true." in aggregate
    assert "requires_spatial_ph .and. .not. methane_ph_valid(ipatch)" in aggregate
    assert "incomplete methane spatial-pH aggregation" in aggregate

    assert "type(methane_ph_mapping_type) :: map_ph" in aggregate
    assert (
        "CALL build_methane_ph_areal_mapping(map_ph, grid_ph, landpatch)"
        in aggregate
    )
    ph_mapping = source(
        "main/TRACER/MOD_Tracer_Reactive_Methane_PHMapping.F90"
    )
    spatial_mapping = source("share/MOD_SpatialMapping.F90")
    assert "build_grid_sumarea" not in spatial_mapping
    assert "SUBROUTINE spatial_mapping_build_arealweighted (this, fgrid, pixelset)" in spatial_mapping
    assert "type, PUBLIC :: methane_ph_mapping_type" in ph_mapping
    assert "SUBROUTINE build_methane_ph_areal_mapping" in ph_mapping
    assert "area = areaquad(lat_s, lat_n, lon_w, lon_e)" in ph_mapping
    assert "areagrid" not in ph_mapping
    assert "allocate_block_data" not in ph_mapping
    assert "get_sumarea" not in ph_mapping
    assert "patch_center_deg" not in aggregate
    assert "native_patch_ph" not in aggregate
    assert aggregate.count("nf90_open") == 2

    sparse = routine(aggregate, "aggregate_sparse_ph")
    assert "mapping%areapart(iset)%val(ipart)" in sparse
    assert "10._r8 ** (-value) * area * depth_weight" in sparse
    assert "valid_area_depth = valid_area_depth + area * depth_weight" in sparse
    assert "patch_ph(iset) = -log10(sum_h_activity_area / valid_area_depth)" in sparse
    assert "CALL mpi_send (grid_values, ng, MPI_REAL8" in sparse
    assert "CALL mpi_send (grid_weights, ng, MPI_REAL8" in sparse
    assert "CALL mpi_recv (recv_values(iproc)%val, ng, MPI_REAL8" in sparse
    assert "CALL mpi_recv (recv_weights(iproc)%val, ng, MPI_REAL8" in sparse

    runs = routine(aggregate, "read_ph_runs")
    assert "ilon(iend+1) /= ilon(iend) + 1" in runs
    assert "count=[nrun, 1, ndepth]" in runs
    assert "sum_h_activity_weight" in runs
    assert "valid_weight = valid_weight + layer_weight(d)" in runs
    assert "valid_weights(i) = valid_weight" in runs
    assert "expected_depth_bottom_cm(4) = [4.5_r8, 9.1_r8, 16.6_r8, 28.9_r8]" in aggregate
    assert "nf90_inq_varid(ncid, 'depth', depth_vid)" in aggregate
    assert "nf90_get_att(ncid, depth_vid, 'units', depth_units)" in aggregate
    assert "depth_g(2:4) <= depth_g(1:3)" in aggregate
    assert "ph_layer_thickness_cm(2:4) = depth_g(2:4) - depth_g(1:3)" in aggregate
    assert "PHH2O top-four depth coordinate is incompatible" in aggregate

    metadata = routine(aggregate, "validate_ph_variable_metadata")
    assert "ndims_loc /= 3" in metadata
    assert "xtype_loc /= NF90_BYTE" in metadata
    assert "PHH2O Fortran dimension order must be lon,lat,depth" in metadata
    assert "dimlen(1) /= nlon" in metadata
    assert "dimlen(2) /= nlat" in metadata
    assert "dimlen(3) /= ndepth" in metadata
    assert "nf90_get_att(ncid, vid, 'units', ph_units)" in metadata
    assert "scale_factor - 0.1_r8" in metadata
    assert "attr_value /= int(missing_byte)" in metadata
    assert "PHH2O missing marker is required" in metadata

    values_x10 = (40.0, 50.0, 60.0, 70.0)
    weights = (4.5, 4.6, 7.5, 12.3)
    expected = -math.log10(
        sum(10 ** (-0.1 * v) * w for v, w in zip(values_x10, weights))
        / sum(weights)
    )
    assert abs(expected - 4.757838723782229) < 1.0e-12

    horizontal = -math.log10((10**-4 + 10**-8) / 2.0)
    assert abs(horizontal - 4.300986568387119) < 1.0e-12

    # A cell with one valid depth layer must not have the same horizontal
    # influence as an equal-area cell with four valid layers.
    missing_aware_horizontal = -math.log10(
        (10**-4 * 1.0 + 10**-8 * 4.0) / (1.0 + 4.0)
    )
    assert abs(missing_aware_horizontal - 4.698796321277554) < 1.0e-12
    assert missing_aware_horizontal != horizontal


def test_inundation_mode_resets_only_hydrology_owned_state_before_dispatch() -> None:
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    configure = routine(const, "configure_methane_inundation_mode")
    reset = configure.split("SELECT CASE", 1)[0]

    for assignment in (
        "DEF_METHANE%enable_wetwat_finundated_override = .false.",
        "DEF_METHANE%wetland_dry_unsat_branch = .true.",
        "DEF_METHANE%use_routing_for_soil = .false.",
    ):
        assert assignment in reset

    for independent_control in (
        "DEF_METHANE%use_biome_f_methane",
        "DEF_METHANE%use_biome_redoxlag",
        "DEF_METHANE%z0_methane_prod",
        "DEF_METHANE%hybrid_soil_threshold",
    ):
        assert independent_control not in configure

    standard = source("run/standard_ch4_parameter.nml")
    assert "DEF_METHANE%use_biome_f_methane            = .true." in standard
    assert "DEF_METHANE%use_biome_redoxlag           = .true." in standard
    assert "DEF_METHANE%z0_methane_prod = 0.30" in standard
    assert "DEF_METHANE%hybrid_soil_threshold = 0.05" in standard


def test_independent_hydrology_and_depth_controls_are_range_checked() -> None:
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    validation = routine(const, "validate_methane_namelist")

    assert "ieee_is_finite" in validation
    for control in (
        "wtd_inflection_soil",
        "wtd_steepness_soil",
        "hybrid_soil_threshold",
        "z0_methane_prod",
    ):
        assert f"DEF_METHANE%{control}" in validation
    assert "hybrid_soil_threshold must be finite and in [0,1]" in validation
    assert "z0_methane_prod must be finite and >= 0 m" in validation


def test_retired_methane_freeze_out_flag_fails_fast() -> None:
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    validation = routine(const, "validate_methane_namelist")

    assert "IF (DEF_METHANE%methane_frzout) THEN" in validation
    assert "methane_frzout is retired" in validation


def test_hydrology_uses_declared_slope_ratio_and_rejects_retired_knob_overrides() -> None:
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    state = source("main/TRACER/MOD_Tracer_Reactive_Methane_State.F90")
    validation = routine(const, "validate_methane_namelist")
    surface = routine(state, "compute_f_h2osfc")

    assert "DEF_METHANE_hydrology%slopebeta >= 0._r8" in validation
    assert "DEF_METHANE_hydrology%vdcf - 2._r8" in validation
    assert "DEF_METHANE_hydrology%pc - 0.4_r8" in validation
    assert "retired methane hydrology knobs" in validation

    assert "slope_angle = atan(slope_arg)" in surface
    assert "slope_arg <= 1._r8" not in surface
    assert "slope_arg * PI / 180._r8" not in surface
    assert "atan(slope_arg / 100._r8)" not in surface
    assert "ieee_is_finite(wdsrf_in)" in surface
    assert "ieee_is_finite(slpratio_in)" in surface


def test_explicit_transient_atmospheric_ch4_cannot_silently_fallback() -> None:
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    validation = routine(const, "validate_methane_namelist")
    lookup = function(const, "methane_atm_mixing_ratio")
    loader = routine(const, "load_methane_atm_file")

    assert "use_transient_atm_methane" in validation
    assert "transient atmospheric CH4 requires atm_methane_file" in validation
    assert "incomplete transient atmospheric CH4 table" in lookup
    assert "cannot open transient atmospheric CH4 file" in loader
    assert "malformed transient atmospheric CH4 table" in loader
    assert "using DEF_METHANE%atm_methane" not in loader


def test_every_real_methane_namelist_parameter_is_checked_for_finiteness() -> None:
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    validation = routine(const, "validate_methane_namelist")
    finite_guard = validation.split("all(ieee_is_finite([", 1)[1].split("])))", 1)[0]

    for type_name, prefix in (
        ("Methane_type", "DEF_METHANE"),
        ("Methane_hydrology_type", "DEF_METHANE_hydrology"),
    ):
        body = re.search(
            rf"type\s+{type_name}\b(.*?)END\s+type\s+{type_name}",
            const,
            flags=re.IGNORECASE | re.DOTALL,
        )
        assert body is not None
        real_fields = re.findall(
            r"real\s*\(r8\)\s*::\s*([A-Za-z0-9_]+)",
            body.group(1),
            flags=re.IGNORECASE,
        )
        for field in real_fields:
            assert f"{prefix}%{field}" in finite_guard

    assert "biome redox lags must be >= 0 days" in validation
    assert "oxinhib must be >= 0" in validation
    assert "microbial biomass fractions must be <= 1" in validation


def test_standard_config_conserves_lake_carbon_and_example_declares_dependencies() -> None:
    standard = source("run/standard_ch4_parameter.nml")
    example = source("run/examples/Global_Grid_2x2_PFT_VG_BGC_CH4.nml")

    assert "DEF_METHANE%replenishlakec        = .false." in standard
    const = source("main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90")
    assert "logical :: replenishlakec = .false." in const
    assert "falls back to pH=6.2" not in standard
    assert "DEF_USE_Dynamic_Wetland = .true." in example
    methane_comment = example.split("! ----- Methane reactive tracer -----", 1)[1]
    assert "CROP" in methane_comment.split("DEF_TRACER_NUM", 1)[0]
