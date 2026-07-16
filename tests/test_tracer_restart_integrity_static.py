from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _source(relative: str) -> str:
    return (ROOT / relative).read_text()


def test_land_restart_dimension_probe_is_global_before_vector_reads():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    probe = source.split("logical FUNCTION tracer_dim_matches", 1)[1].split(
        "END FUNCTION tracer_dim_matches", 1
    )[0]
    reader = source.split("SUBROUTINE read_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_land_tracer_restart", 1
    )[0]

    assert "MPI_SUM, p_comm_glb" in probe
    assert "size(varsize) == expected_rank" in probe
    assert "varsize(expected_rank) == landpatch%vecgs%vlen(iblk,jblk)" in probe
    assert "expected_rank = merge(3, 2, present(expect_soilsnow))" in probe
    assert "varsize(1) == ntransport" in probe
    assert reader.index("tracer_dim_matches(file_restart, 'trc_ldew_rain'") < reader.index(
        "CALL read_transport_patch_field(file_restart, 'trc_ldew_rain'"
    )
    # A collective-bearing function must not be hidden in a short-circuitable
    # OR chain: every IO/worker rank must execute the same probe sequence.
    required_probe = reader.split("! A committed current transaction", 1)[1].split(
        "CALL read_transport_patch_field(file_restart, 'trc_ldew_rain'", 1
    )[0]
    assert ".or." not in required_probe.lower()


def test_land_restart_rejects_descriptor_mismatch_before_state_reads():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    defs = _source("main/TRACER/MOD_Tracer_Defs.F90")
    reader = source.split("SUBROUTINE read_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_land_tracer_restart", 1
    )[0]
    writer = source.split("SUBROUTINE write_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE write_land_tracer_restart", 1
    )[0]
    identity = defs.split("SUBROUTINE tracer_build_descriptor_identity", 1)[1].split(
        "END SUBROUTINE tracer_build_descriptor_identity", 1
    )[0]
    metadata_writer = source.split(
        "SUBROUTINE write_land_tracer_descriptor_metadata", 1
    )[1].split("END SUBROUTINE write_land_tracer_descriptor_metadata", 1)[0]

    assert "LAND_TRACER_RESTART_SCHEMA_VERSION" in source
    assert "TRACER_DESCRIPTOR_IDENTITY_WIDTH" in source
    assert "SUBROUTINE build_land_tracer_descriptor_identity" not in source
    assert "CALL tracer_build_descriptor_identity(expected, transport_only=.true.)" in source
    assert "write_land_tracer_descriptor_metadata(file_restart)" in writer
    assert "trc_land_restart_schema" in metadata_writer
    assert "trc_land_transport_count" in metadata_writer
    assert "trc_land_descriptor_identity" in metadata_writer
    assert "transport_only=.true." in metadata_writer
    assert "IF (size(identity, 2) > 0) THEN" in metadata_writer
    assert "trim(tracers(itrc)%name)" in identity
    assert "trim(tracers(itrc)%category)" in identity
    assert "trim(tracers(itrc)%unit_kind)" in identity
    assert "iachar(descriptor(k:k))" in identity
    assert "tracers(itrc)%family_id" in identity
    assert "tracers(itrc)%state_owner" in identity
    assert "tracers(itrc)%reaction_mode" in identity
    assert "tracers(itrc)%charge" in identity
    for field in (
        "mol_weight",
        "ref_ratio",
        "init_delta",
        "init_conc",
        "precip_default_conc",
        "vapor_default_conc",
        "max_dissolved_conc",
        "reactive_decay_rate",
    ):
        assert f"tracers(itrc)%{field}" in identity
    assert reader.index("land_tracer_descriptor_matches(file_restart)") < reader.index(
        "CALL read_transport_patch_field(file_restart, 'trc_ldew_rain'"
    )
    identity_check = source.split(
        "logical FUNCTION land_tracer_descriptor_matches", 1
    )[1].split("END FUNCTION land_tracer_descriptor_matches", 1)[0]
    assert identity_check.index("ncio_inquire_varsize") < identity_check.index(
        "ncio_read_serial(fileblock, 'trc_land_descriptor_identity'"
    )
    assert "size(varsize) /= 2" in identity_check


def test_land_restart_is_a_compact_atomic_transport_transaction():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    reader = source.split("SUBROUTINE read_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_land_tracer_restart", 1
    )[0]
    writer = source.split("SUBROUTINE write_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE write_land_tracer_restart", 1
    )[0]
    packer = source.split("SUBROUTINE pack_transport_patch", 1)[1].split(
        "END SUBROUTINE pack_transport_patch", 1
    )[0]
    unpacker = source.split("SUBROUTINE read_transport_patch_field", 1)[1].split(
        "END SUBROUTINE read_transport_patch_field", 1
    )[0]

    assert "land_transport_tracer_count()" in reader
    assert "land_transport_tracer_count()" in writer
    assert "trc_land_transport" in writer
    assert "'tracer', ntracers" not in writer
    assert "IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE" in packer
    assert "IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE" in unpacker
    assert "target = 0._r8" in unpacker
    assert reader.index("tracer_dim_matches(file_restart, 'trc_leaf_iso_storage'") < reader.index(
        "read_transport_patch_field(file_restart, 'trc_ldew_rain'"
    )
    assert writer.index("write_land_tracer_transaction_marker(file_restart, 0)") < writer.index(
        "ncio_write_vector(file_restart, 'trc_ldew_rain'"
    )
    assert writer.rindex("write_land_tracer_transaction_marker(file_restart, 1)") > writer.rindex(
        "write_land_tracer_descriptor_metadata(file_restart)"
    )


def test_land_restart_zero_transport_avoids_zero_length_vector_io():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    reader = source.split("SUBROUTINE read_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_land_tracer_restart", 1
    )[0]
    writer = source.split("SUBROUTINE write_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE write_land_tracer_restart", 1
    )[0]

    read_zero = reader.split("IF (ntransport <= 0) THEN", 1)[1].split("ENDIF", 1)[0]
    write_zero = writer.split("IF (ntransport <= 0) THEN", 1)[1].split("ENDIF", 1)[0]
    assert reader.index("descriptor_matches = land_tracer_descriptor_matches(file_restart)") < reader.index(
        "IF (ntransport <= 0) THEN"
    )
    assert "found_restart = .true." in read_zero
    assert "RETURN" in read_zero
    assert "write_land_tracer_transaction_marker(file_restart, 0)" in write_zero
    assert "write_land_tracer_descriptor_metadata(file_restart)" in write_zero
    assert "write_land_tracer_transaction_marker(file_restart, 1)" in write_zero
    assert "ncio_write_vector" not in write_zero
    assert "RETURN" in write_zero
    assert "IF (ntracers <= 0) RETURN" not in reader
    assert "IF (ntracers <= 0) RETURN" not in writer


def test_land_restart_empty_transaction_is_explicit_and_cannot_revive_stale_state():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    descriptor = source.split(
        "logical FUNCTION land_tracer_descriptor_matches", 1
    )[1].split("END FUNCTION land_tracer_descriptor_matches", 1)[0]
    metadata = source.split(
        "SUBROUTINE write_land_tracer_descriptor_metadata", 1
    )[1].split("END SUBROUTINE write_land_tracer_descriptor_metadata", 1)[0]

    assert "LAND_TRACER_RESTART_SCHEMA_VERSION = 4" in source
    assert "'trc_land_transport_count', size(identity, 2)" in metadata
    assert metadata.index("IF (size(identity, 2) > 0) THEN") < metadata.index(
        "'trc_land_descriptor_identity'"
    )
    assert "transport_count >= 0 .and. transport_count <= 1000" in descriptor
    assert "transport_count > 0" in descriptor
    assert "ELSEIF (transport_count == 0) THEN" in descriptor
    # An explicit zero count dominates any stale descriptor variable left by a
    # reused NetCDF target; only positive committed counts read identity bytes.
    positive_identity_gate = descriptor.split("transport_count > 0", 1)[1].split(
        "ELSEIF (block_ok .and. schema /=", 1
    )[0]
    assert "ncio_read_serial(fileblock, 'trc_land_descriptor_identity'" in positive_identity_gate
    assert "complete /= 1" in descriptor
    assert "(counts(3) > 0 .and. counts(4) > 0)" in descriptor


def test_land_restart_requires_one_exact_generation_across_all_blocks():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    descriptor = source.split(
        "logical FUNCTION land_tracer_descriptor_matches", 1
    )[1].split("END FUNCTION land_tracer_descriptor_matches", 1)[0]

    # Every IO rank compares all local blocks exactly, then all ranks compare
    # their exact representatives against one broadcast reference.  A hash is
    # deliberately insufficient because collisions could admit mixed state.
    assert "schema /= local_reference_schema" in descriptor
    assert "transport_count /= local_reference_count" in descriptor
    assert "any(ondisk /= local_reference_identity)" in descriptor
    assert "CALL mpi_bcast(reference_metadata" in descriptor
    assert "CALL mpi_bcast(reference_identity" in descriptor
    assert "any(local_reference_identity /= reference_identity)" in descriptor
    assert "generation_mismatch .or. counts(5) > 0" in descriptor

    # The executable truth table documents the formerly missed cases.  Count
    # zero has no identity key, so stale identity bytes cannot affect equality.
    identity_a = ((65, 66), (67, 68))
    identity_b = ((65, 66), (67, 69))

    def uniform(blocks):
        return all(block == blocks[0] for block in blocks[1:])

    assert not uniform(((4, 0, None), (4, 1, identity_a)))
    assert not uniform(((3, 1, identity_a), (2, 1, identity_a)))
    assert not uniform(((4, 1, identity_a), (4, 1, identity_b)))

    # A coherent descriptor for a different configuration remains safe to
    # cold-start as one unit rather than being reclassified as corruption.
    current = (4, 1, identity_a)
    same_noncurrent = ((4, 1, identity_b), (4, 1, identity_b))
    assert uniform(same_noncurrent)
    assert same_noncurrent[0] != current
    assert "land_tracer_descriptor_matches = counts(1) > 0 .and. counts(3) == counts(1)" in descriptor


def test_land_restart_validation_enforces_physical_domains_collectively():
    source = _source("main/TRACER/MOD_Tracer_Rest.F90")
    validator = source.split("SUBROUTINE validate_land_tracer_restart_state", 1)[1].split(
        "END SUBROUTINE validate_land_tracer_restart_state", 1
    )[0]
    signed = source.split("SUBROUTINE count_signed_patch", 1)[1].split(
        "END SUBROUTINE count_signed_patch", 1
    )[0]
    nonnegative = source.split("SUBROUTINE count_nonnegative_patch", 1)[1].split(
        "END SUBROUTINE count_nonnegative_patch", 1
    )[0]
    peclet = source.split("SUBROUTINE count_peclet_patch", 1)[1].split(
        "END SUBROUTINE count_peclet_patch", 1
    )[0]

    assert "MPI_INTEGER, MPI_SUM, p_comm_glb" in validator
    assert "count_signed_patch(trc_wa" in validator
    assert "count_signed_patch(trc_leaf_iso_storage" in validator
    assert "count_signed_patch(trc_leaf_delta_e" in validator
    assert "count_signed_patch(trc_leaf_delta_b" in validator
    assert "count_nonnegative_patch(trc_leaf_water_moles" in validator
    for field in ("trc_ldew_rain", "trc_wdsrf", "trc_surface_residue", "trc_canopy_solid"):
        assert f"count_nonnegative_patch({field}" in validator
    assert "ieee_is_finite" in signed
    assert "< -LAND_TRACER_RESTART_NEGATIVE_DUST" in nonnegative
    assert "field(itrc, ip) < 0._r8 .or. field(itrc, ip) > 1._r8" in peclet
    assert "invalid generic land tracer restart state" in validator


def test_provider_restart_is_independent_of_generic_land_descriptor_gate():
    land_phase = _source("main/TRACER/MOD_Tracer_LandPhase.F90")
    methane = _source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    init = land_phase.split("SUBROUTINE tracer_init_from_arrays", 1)[1].split(
        "END SUBROUTINE tracer_init_from_arrays", 1
    )[0]
    reader = methane.split("SUBROUTINE ch4_reactive_read_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_read_restart", 1
    )[0]

    assert init.count("tracer_lifecycle_land_read_restart(file_restart)") == 1
    assert "ENDIF\n      ! Provider-owned state" in init
    assert "IF (present(file_restart)) CALL tracer_lifecycle_land_read_restart" in init
    assert "provider_restart_probe_fields" in reader
    assert "IF (.not. any(provider_restart_probe_present)) RETURN" in reader
    assert "metadata exists without mandatory ch4_conc_methane" in reader
    assert reader.index("ncio_vector_group_presence") < reader.index(
        "validate_methane_restart_transaction"
    )


def test_tracer_forcing_window_advances_until_it_brackets_model_time():
    source = _source("main/TRACER/MOD_Tracer_Forcing.F90")
    reader = source.split("SUBROUTINE tracer_forcing_read_LBUB", 1)[1].split(
        "END SUBROUTINE tracer_forcing_read_LBUB", 1
    )[0]

    assert "DO WHILE" in reader
    loop = reader.split("DO WHILE", 1)[1].split("ENDDO", 1)[0]
    assert "trc_tstamp_UB(iv) <= mtstamp" in loop
    assert "block_data_copy(trc_forcn_UB(iv), trc_forcn_LB(iv))" in loop
    assert "tracer_forcing_setstamp_UB" in loop


def test_complete_vector_reader_rejects_partial_block_variables():
    source = _source("share/MOD_NetCDFVector.F90")

    assert "INTERFACE ncio_read_vector_complete" in source
    assert "ncio_require_complete_vector_var" in source
    check = source.split("SUBROUTINE ncio_require_complete_vector_var", 1)[1].split(
        "END SUBROUTINE ncio_require_complete_vector_var", 1
    )[0]
    assert "ncio_var_exist" in check
    assert "MPI_SUM, p_comm_io" in check
    assert "p_root, p_comm_group" in check
    assert "counts(2) > 0" in check
    assert "counts(2) < counts(1)" in check
    assert "CALL CoLM_stop" in check


def test_complete_vector_preflight_skips_the_control_master_singleton_group():
    source = _source("share/MOD_NetCDFVector.F90")
    check = source.split("SUBROUTINE ncio_require_complete_vector_var", 1)[1].split(
        "END SUBROUTINE ncio_require_complete_vector_var", 1
    )[0]

    assert "p_is_worker" in check
    assert "IF (.not. (p_is_io .or. p_is_worker)) RETURN" in check


def test_complete_vector_reader_rejects_missing_required_and_wrong_shape_before_read():
    source = _source("share/MOD_NetCDFVector.F90")
    check = source.split("SUBROUTINE ncio_require_complete_vector_var", 1)[1].split(
        "END SUBROUTINE ncio_require_complete_vector_var", 1
    )[0]
    wrapper_1d = source.split("SUBROUTINE ncio_read_vector_complete_real8_1d", 1)[1].split(
        "END SUBROUTINE ncio_read_vector_complete_real8_1d", 1
    )[0]
    wrapper_2d = source.split("SUBROUTINE ncio_read_vector_complete_real8_2d", 1)[1].split(
        "END SUBROUTINE ncio_read_vector_complete_real8_2d", 1
    )[0]

    assert "allow_missing" in check
    assert "counts(2) == 0 .and. .not. allow_missing" in check
    assert "ncio_inquire_varsize" in check
    assert "size(varsize) == expected_rank" in check
    assert "varsize(expected_rank) == pixelset%vecgs%vlen(iblk,jblk)" in check
    assert "varsize(1) == expected_dim1" in check
    assert "counts(3) > 0" in check
    assert wrapper_1d.index("ncio_require_complete_vector_var") < wrapper_1d.index(
        "ncio_read_vector_real8_1d"
    )
    assert "expected_rank=1" in wrapper_1d
    assert "allow_missing=allow_field_missing" in wrapper_1d
    assert wrapper_2d.index("ncio_require_complete_vector_var") < wrapper_2d.index(
        "ncio_read_vector_real8_2d"
    )
    assert "expected_rank=2" in wrapper_2d
    assert "expected_dim1=ndim1" in wrapper_2d
    assert "allow_missing=allow_field_missing" in wrapper_2d


def test_all_ch4_restart_readers_opt_into_complete_block_policy():
    for relative, routine in (
        ("main/TRACER/MOD_Tracer_Reactive_Methane_State.F90", "read_methane_restart"),
        ("main/TRACER/MOD_Tracer_Reactive_Methane_Microbes.F90", "read_methane_microbes_restart"),
        ("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90", "read_methane_accflux_restart"),
    ):
        source = _source(relative)
        body = source.split(f"SUBROUTINE {routine}", 1)[1].split(
            f"END SUBROUTINE {routine}", 1
        )[0]
        if routine == "read_methane_accflux_restart":
            assert "INTERFACE ncio_read_vector" in source
            assert source.count(
                "USE MOD_NetCDFVector, only: ncio_read_vector_complete"
            ) >= 2
        else:
            assert (
                "ncio_read_vector => ncio_read_vector_complete" in body
            ), f"{routine} must reject variables present in only some restart blocks"


def test_ch4_restart_uses_schema_in_progress_and_commit_markers():
    source = _source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    reader = source.split("SUBROUTINE ch4_reactive_read_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_read_restart", 1
    )[0]
    writer = source.split("SUBROUTINE ch4_reactive_write_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_write_restart", 1
    )[0]
    validator = source.split("SUBROUTINE validate_methane_restart_transaction", 1)[1].split(
        "END SUBROUTINE validate_methane_restart_transaction", 1
    )[0]

    assert "ch4_restart_schema" in validator
    assert "ch4_restart_complete" in validator
    assert "#ifdef USEMPI" in validator
    assert "#else\n      USE MOD_SPMD_Task, only: p_is_worker, p_is_master" in validator
    assert "IF (.not. has_schema .and. .not. has_commit) THEN" in validator
    assert "IF (has_history_mode .or. has_microbe_feature) THEN" in validator
    assert "orphaned methane restart transaction metadata" in validator
    assert "IF (.not. has_schema .or. .not. has_commit) THEN" in validator
    assert "logical, intent(out) :: strict_restart" in validator
    assert "strict_restart = .true." in validator
    assert "ncio_set_complete_require_present (strict_restart)" in reader
    assert "ncio_set_complete_require_present (.false.)" in reader
    assert "CALL ncio_vector_group_presence" in reader
    assert "n_microbe_state_fields = count(microbe_field_present(1:4))" in reader
    assert "partial/inconsistent microbial feature group" in reader
    assert "IF (file_has_pools) THEN" in reader
    assert "IF (DEF_METHANE%use_microbial_pools) THEN" in reader
    assert "validate_methane_microbes_restart_values" in reader
    assert reader.index("validate_methane_restart_transaction") < reader.index(
        "read_methane_restart"
    )
    assert writer.count("write_methane_restart_marker") == 5
    assert "'ch4_feature_microbe_pools'" in writer
    assert writer.index("write_methane_restart_marker(file_restart, 'ch4_restart_complete', 0._r8") < writer.index(
        "write_methane_restart(file_restart"
    )
    assert writer.index("write_methane_restart_marker(file_restart, 'ch4_history_accumulation_mode'") < writer.index(
        "write_methane_restart(file_restart"
    )
    assert writer.index("write_methane_restart_marker(file_restart, 'ch4_feature_microbe_pools'") < writer.index(
        "write_methane_restart(file_restart"
    )
    assert writer.rindex("write_methane_restart_marker(file_restart, 'ch4_restart_complete', 1._r8") > writer.index(
        "write_methane_microbes_restart(file_restart"
    )


def test_ch4_presence_probes_use_block_aware_vector_metadata():
    methane = _source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    accflux = _source("main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90")
    reader = methane.split("SUBROUTINE ch4_reactive_read_restart", 1)[1].split(
        "END SUBROUTINE ch4_reactive_read_restart", 1
    )[0]
    acc_reader = accflux.split("SUBROUTINE read_methane_accflux_restart", 1)[1].split(
        "END SUBROUTINE read_methane_accflux_restart", 1
    )[0]

    assert "ncio_vector_group_presence" in reader
    assert "ncio_var_exist(file_restart" not in reader
    assert "ncio_vector_var_present" in acc_reader
    assert "ncio_var_exist(file_restart" not in acc_reader
