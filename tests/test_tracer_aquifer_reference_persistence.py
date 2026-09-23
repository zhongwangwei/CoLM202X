from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def source(name):
    return (ROOT / name).read_text()


def test_reference_is_explicit_and_broadcast():
    nml = source("share/MOD_Namelist.F90")
    assert "DEF_TRACER_AQUIFER_MIXING_WATER_MM = -1._r8" in nml
    assert nml.count("DEF_TRACER_AQUIFER_MIXING_WATER_MM") >= 3
    assert "mpi_bcast (DEF_TRACER_AQUIFER_MIXING_WATER_MM" in nml


def test_restart_persists_reference_and_rejects_incompatible_isotope_state():
    rest = source("main/TRACER/MOD_Tracer_Rest.F90")
    assert "LAND_TRACER_RESTART_SCHEMA_VERSION = 5" in rest
    for field in ("trc_aquifer_ref_water", "trc_aquifer_ref_mass"):
        assert f"tracer_dim_matches(file_restart, '{field}'" in rest
        assert f"'trc_aquifer_ref_{field.rsplit('_', 1)[-1]}'" in rest
        assert f"'{field}'," in rest
    assert "ondisk_mixing == DEF_TRACER_AQUIFER_MIXING_WATER_MM" in rest
    assert "IF (active_isotope .and. state_present)" in rest
    assert "old/incompatible isotope tracer restart" in rest
    assert "(wa(ip) + trc_aquifer_ref_water(ip)) * R_init" in rest


def test_provider_only_schema4_restart_is_not_old_generic_isotope_inventory():
    rest = source("main/TRACER/MOD_Tracer_Rest.F90")
    descriptor = rest.split("logical FUNCTION land_tracer_descriptor_matches", 1)[1].split(
        "END FUNCTION land_tracer_descriptor_matches", 1
    )[0]
    assert "schema == 4" in descriptor
    assert "IF (block_ok .and. transport_count > 0) counts(6)" in descriptor
    assert "IF (ncio_var_exist(fileblock, 'trc_wa'" in descriptor
    assert "IF (has_commit .or. has_schema" not in descriptor
    assert "require_isotope_reference .and. transport_count > 0" in descriptor


def test_schema4_nonisotope_generic_hot_state_remains_readable():
    rest = source("main/TRACER/MOD_Tracer_Rest.F90")
    descriptor = rest.split("logical FUNCTION land_tracer_descriptor_matches", 1)[1].split(
        "END FUNCTION land_tracer_descriptor_matches", 1
    )[0]
    reader = rest.split("SUBROUTINE read_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE read_land_tracer_restart", 1
    )[0]
    assert "(schema /= 4 .or. require_isotope_reference)" in descriptor
    assert "IF (restart_schema == LAND_TRACER_RESTART_SCHEMA_VERSION) THEN" in reader
    assert "trc_aquifer_ref_water = 0._r8" in reader
    assert "trc_aquifer_ref_mass = 0._r8" in reader
    assert "CALL read_transport_patch_field(file_restart, 'trc_wa', trc_wa)" in reader


def test_history_keeps_isotope_carrier_separate_from_solute_denominator():
    hist = source("main/TRACER/MOD_Tracer_Hist.F90")
    assert "a_water_aquifer_actual(ipatch)" in hist
    assert "a_trc_aquifer_actual_mass(itrc, ipatch)" in hist
    assert "IF (tracer_is_isotope(itrc_loc)) THEN" in hist
    assert "a_trc_wa_mass(itrc_loc, :), a_water_wa" in hist


def test_lulcc_remaps_reference_and_rejects_transfer_to_special_patch():
    vars_ = source("main/TRACER/MOD_Tracer_Vars.F90")
    assert "CALL remap2d_mass(lulcc_trc_aquifer_ref_mass_old, trc_aquifer_ref_mass)" in vars_
    assert "CALL remap2d_mass(lulcc_trc_aquifer_ref_water_old, remapped_ref_water)" in vars_
    assert "LULCC cannot transfer isotope aquifer reference to special patch" in vars_
    assert "LULCC cannot create soil or wetland isotope aquifer without reference" in vars_
