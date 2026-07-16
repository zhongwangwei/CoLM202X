from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
RIVER = ROOT / "main" / "TRACER" / "MOD_Tracer_RiverLake.F90"
FLOW = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeFlow.F90"


def test_all_generic_river_loops_use_the_transport_contract():
    source = RIVER.read_text(encoding="utf-8")

    assert "tracer_uses_land_water_transport" in source.split("IMPLICIT NONE", 1)[0]
    assert "legacy_tracer_count_meta" not in source
    assert "legacy_tracer_namehash_meta" not in source

    tracer_loops = re.findall(
        r"DO\s+itrc\w*\s*=\s*1\s*,\s*ntracers",
        source,
        flags=re.IGNORECASE,
    )
    guarded_loops = re.findall(
        r"DO\s+(?P<index>itrc\w*)\s*=\s*1\s*,\s*ntracers\s*"
        r"IF\s*\(\.not\.\s*tracer_uses_land_water_transport\("
        r"(?P=index)\)\)\s*CYCLE",
        source,
        flags=re.IGNORECASE,
    )
    assert tracer_loops
    assert len(guarded_loops) == len(tracer_loops)


def test_restart_writer_persists_the_full_shared_transport_descriptor():
    source = RIVER.read_text(encoding="utf-8")
    writer = source.split("SUBROUTINE write_tracer_restart", 1)[1]
    writer = writer.split("END SUBROUTINE write_tracer_restart", 1)[0]

    assert (
        "tracer_build_descriptor_identity(descriptor_identity, transport_only=.true.)"
        in writer
    )
    assert "TRACER_DESCRIPTOR_IDENTITY_WIDTH" in writer
    assert "trc_river_restart_schema" in writer
    assert "trc_river_restart_complete" in writer
    assert "trc_river_descriptor_count" in writer
    assert "trc_river_descriptor_identity" in writer
    assert "trc_n_meta" not in source
    assert "trc_namehash_meta" not in source
    assert "tracer_count_meta" not in source
    assert "tracer_namehash_meta" not in source
    assert "trc_ucat_gdid_meta" in writer
    assert "trc_ucat_next_meta" in writer
    assert "ucat_gdid" in writer
    assert "ucat_next" in writer
    assert writer.count("trc_river_restart_complete") == 4
    assert writer.rindex("'trc_river_restart_complete', 0") < writer.index(
        "CALL vector_gather_and_write"
    )
    assert writer.rindex("'trc_river_restart_complete', 1") > writer.rindex(
        "vector_gather_and_write"
    )
    assert writer.rindex("'trc_river_restart_schema'") < writer.rindex(
        "'trc_river_restart_complete', 1"
    )
    assert writer.index("CALL validate_riverlake_restart_state()") < writer.rindex(
        "'trc_river_restart_complete', 0"
    )
    zero_transport = writer.split("IF (size(descriptor_identity, 2) <= 0) THEN", 1)[1].split(
        "RETURN", 1
    )[0]
    assert "'trc_river_restart_complete', 0" in zero_transport
    assert "'trc_river_descriptor_count', 0" in zero_transport
    assert "RIVER_TRACER_RESTART_SCHEMA_VERSION" in zero_transport
    assert "'trc_river_restart_complete', 1" in zero_transport
    assert "trc_river_descriptor_identity" not in zero_transport
    assert "ncio_define_dimension" not in zero_transport
    assert "deallocate(descriptor_identity)" in zero_transport


def test_restart_reader_uses_the_full_descriptor_as_the_only_tracer_identity():
    source = RIVER.read_text(encoding="utf-8")
    descriptor = source.split(
        "logical FUNCTION riverlake_restart_descriptor_compatible", 1
    )[1].split("END FUNCTION riverlake_restart_descriptor_compatible", 1)[0]
    reader = source.split("SUBROUTINE read_tracer_restart", 1)[1]
    reader = reader.split("END SUBROUTINE read_tracer_restart", 1)[0]

    assert "trc_n_meta" not in source
    assert "trc_namehash_meta" not in source
    assert "tracer_count_meta" not in source
    assert "tracer_namehash_meta" not in source
    assert "tracer_is_particle" not in source.split("IMPLICIT NONE", 1)[0]
    assert "flags(4) == 0" in descriptor
    assert "flags(1) /= 1 .or. flags(2) /= 1" in descriptor
    assert "ELSEIF (flags(3) /= 1)" in descriptor
    assert "trc_river_restart_complete" in descriptor
    assert "IF (complete /= 1) status = -1" in descriptor
    marker_rank_probe = descriptor.split(
        "ncio_inquire_varsize(file_restart, 'trc_river_restart_complete'", 1
    )[1].split("IF (status == 1)", 1)[0]
    assert "size(varsize) /= 0" in marker_rank_probe
    assert "identity_count /= descriptor_count" in descriptor
    assert "descriptor_count == 0" in descriptor
    assert "size(expected, 2) == 0" in descriptor
    empty_branch = descriptor.split("ELSEIF (descriptor_count == 0) THEN", 1)[1].split(
        "ELSEIF (flags(3) /= 1)", 1
    )[0]
    assert "trc_river_descriptor_identity" not in empty_branch
    assert "schema /= RIVER_TRACER_RESTART_SCHEMA_VERSION" in descriptor
    assert "descriptor_count /= size(expected, 2)" in descriptor
    assert "status = 0" in descriptor
    assert "status = -1" in descriptor
    assert reader.index("riverlake_restart_descriptor_compatible(file_restart)") < reader.index(
        "vector_read_and_scatter"
    )
    cold_start = reader.split("IF (.not. descriptor_compatible) THEN", 1)[1].split(
        "! A committed current transaction", 1
    )[0]
    assert "missing_mask(itrc) = .true." in cold_start
    assert "found_restart = .false." in cold_start

    assert "single tracer-identity source of truth" in reader
    assert "trc_numucat_meta" in reader


def test_provider_only_reader_validates_transaction_before_ignoring_generic_state():
    source = RIVER.read_text(encoding="utf-8")
    descriptor = source.split(
        "logical FUNCTION riverlake_restart_descriptor_compatible", 1
    )[1].split("END FUNCTION riverlake_restart_descriptor_compatible", 1)[0]
    reader = source.split("SUBROUTINE read_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_tracer_restart", 1
    )[0]

    validation = (
        "descriptor_compatible = "
        "riverlake_restart_descriptor_compatible(file_restart)"
    )
    provider_gate = "IF (.not. has_transport_tracer) THEN"
    first_vector_probe = "CALL probe_riverlake_restart_vector"

    assert reader.count("riverlake_restart_descriptor_compatible(file_restart)") == 1
    assert reader.index(validation) < reader.index(provider_gate)
    assert reader.index(provider_gate) < reader.index(first_vector_probe)
    provider_exit = reader.split(provider_gate, 1)[1].split("ENDIF", 1)[0]
    assert "found_restart = .true." in provider_exit
    assert "RETURN" in provider_exit
    assert "probe_riverlake_restart_vector" not in provider_exit
    assert "vector_read_and_scatter" not in provider_exit
    assert "IF (.not. descriptor_compatible) THEN" in reader

    # Legacy/incompatible metadata returns false and is harmless when there are
    # no generic rows. A present but uncommitted/partial transaction is damage
    # and the shared validator must stop before the provider-only return.
    assert "flags(4) == 0" in descriptor
    assert "IF (complete /= 1) status = -1" in descriptor
    assert "flags(1) /= 1 .or. flags(2) /= 1" in descriptor
    assert "malformed/incomplete river tracer restart descriptor metadata" in descriptor



def test_restart_state_is_collectively_validated_before_transport():
    source = RIVER.read_text(encoding="utf-8")
    validator = source.split("SUBROUTINE validate_riverlake_restart_state", 1)[1].split(
        "END SUBROUTINE validate_riverlake_restart_state", 1
    )[0]
    reader = source.split("SUBROUTINE read_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_tracer_restart", 1
    )[0]

    for state in (
        "trc_mass(itrc, i)",
        "trc_inp_buf(itrc, i)",
        "acc_trc_inp(itrc, i)",
        "acc_rnof_ref(i)",
        "trc_levsto(itrc, i)",
    ):
        assert f"ieee_is_finite({state})" in validator
    assert "trc_mass(itrc, i) < -TRC_RESTART_NEGATIVE_DUST" in validator
    assert "trc_levsto(itrc, i) < -TRC_RESTART_NEGATIVE_DUST" in validator
    assert "MPI_INTEGER, MPI_SUM" in validator
    assert "invalid river/lake tracer restart state" in validator
    assert "trc_inp_buf < 0" not in validator
    assert "acc_trc_inp < 0" not in validator
    assert reader.count("CALL validate_riverlake_restart_state()") == 2
    assert reader.rindex("CALL validate_riverlake_restart_state()") < reader.index(
        "deallocate (tmpvec)"
    )


def test_restart_rejects_same_size_different_network_identity():
    source = RIVER.read_text(encoding="utf-8")
    reader = source.split("SUBROUTINE read_tracer_restart", 1)[1]
    reader = reader.split("END SUBROUTINE read_tracer_restart", 1)[0]

    for field, current in (
        ("trc_ucat_gdid_meta", "ucat_gdid"),
        ("trc_ucat_next_meta", "ucat_next"),
    ):
        assert field in reader
        assert current in reader
    assert "network_meta_matches" in reader
    preflight = reader.split("! A committed current transaction", 1)[1].split(
        "! vector_read_and_scatter", 1
    )[0]
    assert "probe_riverlake_restart_vector(file_restart, 'trc_ucat_gdid_meta', .true." in preflight
    assert "probe_riverlake_restart_vector(file_restart, 'trc_ucat_next_meta', .true." in preflight
    assert reader.index("probe_riverlake_restart_vector(file_restart, 'trc_ucat_gdid_meta'") < reader.index(
        "vector_read_and_scatter(file_restart, tmpvec"
    )
    assert "belongs to a different catchment network" in reader


def test_runoff_tracer_is_added_once_without_substep_smoothing():
    source = RIVER.read_text(encoding="utf-8")
    flow = FLOW.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]

    assert "trc_inp_buf(itrc, :) = trc_inp_buf(itrc, :) + acc_trc_inp(itrc, :)" in flow
    assert "acc_trc_inp(itrc, :) = 0._r8" in flow
    assert "release = max(trc_inp_buf(itrc, i), 0._r8)" in substep
    release_block = substep.split("release = max(trc_inp_buf(itrc, i), 0._r8)", 1)[1].split(
        "CALL update_tracer_concentration", 1
    )[0]
    assert "IF (volwater > trc_v_dry_off) THEN" not in release_block
    assert "inj_frac" not in substep
    assert "m_tau" not in substep

    def add_input(buffer: float, runoff_mass: float) -> tuple[float, float]:
        buffer += runoff_mass
        runoff_mass = 0.0
        release = max(buffer, 0.0)
        return buffer - release, release + runoff_mass

    for substeps in (1, 2, 8, 64):
        buffer, delivered = add_input(0.0, 7.5)
        for _ in range(substeps - 1):
            buffer, extra = add_input(buffer, 0.0)
            delivered += extra
        assert delivered == 7.5
        assert buffer == 0.0
