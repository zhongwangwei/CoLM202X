from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def source(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def routine(text: str, name: str) -> str:
    return text.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_residue_pools_have_complete_state_lifecycle():
    variables = source("MOD_Tracer_Vars.F90")
    for name in ("trc_surface_residue", "trc_subsurface_residue"):
        assert f"allocate({name}" in variables
        assert f"allocated({name})" in variables
        assert f"lulcc_{name}_old" in variables
        assert f"CALL remap2d_mass(lulcc_{name}_old" in variables
        assert f"CALL accumulate_lulcc_mass_2d({name}" in variables


def test_residue_restart_is_required_in_committed_transaction_and_round_trips():
    restart = source("MOD_Tracer_Rest.F90")
    reader = routine(restart, "read_land_tracer_restart")
    writer = routine(restart, "write_land_tracer_restart")
    preflight = reader.split("! A committed current transaction", 1)[1].split(
        "CALL read_transport_patch_field", 1
    )[0]
    for name in ("trc_surface_residue", "trc_subsurface_residue"):
        assert f"tracer_dim_matches(file_restart, '{name}')" in preflight
        assert f"read_transport_patch_field(file_restart, '{name}'" in reader
        assert f"ncio_write_vector(file_restart, '{name}'" in writer
        assert f"pack_transport_patch({name}, restart_patch)" in writer
    assert "incomplete or malformed committed generic land tracer restart" in reader


def test_residue_is_a_conserved_storage_component():
    conservation = source("MOD_Tracer_Conservation.F90")
    assert "n_storage_diag = 12" in conservation
    assert "storage_comp(10) = trc_surface_residue" in conservation
    assert "storage_comp(11) = trc_subsurface_residue" in conservation
    assert "storage_comp_end(10) = trc_surface_residue" in conservation
    assert "storage_comp_end(11) = trc_subsurface_residue" in conservation


def test_all_exhausted_surface_phases_use_shared_nonvolatile_guard():
    soil = source("MOD_Tracer_SoilWater.F90")
    helper = routine(soil, "exhaust_surface_phase")
    assert soil.count("CALL exhaust_surface_phase") == 5
    assert "tracer_is_nonvolatile_solute(itrc)" in helper
    nonvolatile_branch = helper.split(
        "IF (tracer_is_nonvolatile_solute(itrc)) THEN", 1
    )[1].split("ENDIF", 1)[0]
    assert "RETURN" in nonvolatile_branch
    assert "phase_tracer = 0._r8" not in nonvolatile_branch
    assert "trc_surface_residue" not in helper
    assert "CALL tracer_book_evap_loss" in helper
    assert not re.search(r"trc_flux\s*=\s*max\(trc_w(?:ice|liq)_soisno", soil)


def test_surface_and_aquifer_collapse_transfer_without_numerical_sink():
    soil = source("MOD_Tracer_SoilWater.F90")
    assert re.search(
        r"tracer_is_nonvolatile_solute\(itrc\).*?"
        r"trc_surface_residue\(itrc, ipatch\).*?trc_pool_total",
        soil,
        re.DOTALL,
    )
    assert re.search(
        r"wa_bef > trc_water_min_for_ratio.*?"
        r"trc_wa\(itrc, ipatch\).*?trc_subsurface_residue",
        soil,
        re.DOTALL,
    )
    assert re.search(
        r"abs\(wa\) <= trc_water_min_for_ratio.*?tracer_is_nonvolatile_solute\(itrc\).*?"
        r"trc_subsurface_residue\(itrc, ipatch\).*?trc_wa\(itrc, ipatch\)",
        soil,
        re.DOTALL,
    )
    assert soil.count("wa > trc_water_min_for_ratio") >= 2


def test_special_patch_drydown_preserves_surface_residue():
    special = source("MOD_Tracer_SpecialPatches.F90")
    assert special.count("surface_residue_beg + trc_final") == 2
    assert special.count(
        "trc_subsurface_residue(itrc, ipatch) + surface_residue_beg"
    ) == 2
    assert special.count("trc_surface_residue(itrc, ipatch) = trc_surface_residue(itrc, ipatch) +") >= 2
    assert special.count("trc_ldew_rain(itrc, ipatch) = 0._r8") == 2
    assert special.count("trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) +") == 2
    assert "trc_wliq_soisno(itrc, snl+1, ipatch)" in special
    assert special.count("surface_residue_export = min(surface_residue_beg") == 2
    assert special.count("trc_rnof = trc_rnof + surface_residue_export") == 2
    assert "wliq_soisno(snl+1) > trc_water_min_for_ratio" in special
    assert "wliq_soisno(snl+1) + wice_soisno(snl+1)" not in special
    assert "trc_held_storage + trc_input" in special

    residue = 4.0
    runoff = 3.0
    surface_liquid_end = 1.0
    exported = residue * runoff / (runoff + surface_liquid_end)
    assert exported == 3.0
    assert residue - exported == 1.0


def test_special_patch_subsurface_residue_has_explicit_surface_destination():
    special = source("MOD_Tracer_SpecialPatches.F90")
    helper = special.split("SUBROUTINE move_special_subsurface_residue_to_surface", 1)[1].split(
        "END SUBROUTINE move_special_subsurface_residue_to_surface", 1
    )[0]
    assert "trc_subsurface_residue(itrc, ipatch) <= trc_tiny" in helper
    assert "trc_surface_solid(itrc, ipatch) = trc_surface_solid(itrc, ipatch) +" in helper
    assert "trc_surface_residue(itrc, ipatch) = trc_surface_residue(itrc, ipatch) +" in helper
    assert helper.index("trc_subsurface_residue(itrc, ipatch)") < helper.rindex(
        "trc_subsurface_residue(itrc, ipatch) = 0._r8"
    )

    for name in ("tracer_glacier_patch", "tracer_waterbody_patch"):
        body = routine(special, name)
        assert body.index("CALL move_special_subsurface_residue_to_surface") < body.index(
            "CALL tracer_save_storage"
        )

    waterbody = routine(special, "tracer_waterbody_patch")
    assert "trc_held_storage = trc_held_storage - trc_subsurface_residue" not in waterbody


def test_nonvolatile_wblc_ice_sink_stays_in_layer():
    soil = source("MOD_Tracer_SoilWater.F90")
    wblc = soil.split("IF (wblc_ice_sink(j) > trc_tiny) THEN", 1)[1].split(
        "d_wice = d_wice + wblc_ice_sink(j)", 1
    )[0]
    guard = wblc.split("IF (.not. tracer_is_nonvolatile_solute(itrc)) THEN", 1)[1]
    assert "trc_wice_soisno(itrc, j, ipatch) =" in guard
    assert "CALL tracer_book_evap_loss" in guard


def test_dry_layer_inventory_is_not_counted_as_dissolved_concentration():
    history = source("MOD_Tracer_Hist.F90")
    variables = source("MOD_Tracer_Vars.F90")
    assert "layer_water <= trc_water_min_for_ratio" in history
    assert history.count("a_trc_layer_dry_mass(itrc, ipatch)") >= 4
    assert "f_trc_layer_dry_inventory_" in history
    assert "allocate(a_trc_layer_dry_mass" in variables
    assert "allocated(a_trc_layer_dry_mass)" in variables


def test_residue_history_is_inventory_not_concentration():
    # Check each variable's own unit rather than counting how often the unit
    # string appears in the file -- that count moves whenever an unrelated
    # inventory-style output is added, which says nothing about these three.
    history = source("MOD_Tracer_Hist.F90")
    for var in (
        "f_trc_surface_residue_",
        "f_trc_subsurface_residue_",
        "f_trc_layer_dry_inventory_",
    ):
        assert var in history
        block = history.split(var, 1)[1][:800]
        assert "'tracer amount/m2'" in block, f"{var} must be written as an inventory"
        assert "'permil'" not in block, f"{var} must not be written as a delta"


def test_internal_transfer_identity():
    mobile_before = 3.0e-4
    residue_before = 2.0e-4
    mobile_after_dry = 0.0
    residue_after_dry = residue_before + mobile_before
    assert mobile_before + residue_before == mobile_after_dry + residue_after_dry

    mobile_after_rewet = mobile_after_dry + residue_after_dry
    residue_after_rewet = 0.0
    assert mobile_after_rewet + residue_after_rewet == residue_after_dry


def test_vapor_exchange_history_is_not_solute_only():
    history = source("MOD_Tracer_Hist.F90")
    block = history.split("f_trc_vapor_exchange_", 1)[0][-900:]
    assert "tracer_is_nonvolatile_solute" not in block, (
        "isotope vapor exchange must be reachable without pretending the tracer is a solute"
    )
    assert "tracer_uses_delta_diagnostics(itrc_loc)" in block
    out = history.split("f_trc_vapor_exchange_", 1)[1][:700]
    assert "write_history_variable_2d" in out
    assert "'tracer amount/m2'" in out
    assert "write_history_tracer_delta_2d" not in out


def test_signed_evaporation_mass_history_keeps_negative_exchange_observable():
    history = source("MOD_Tracer_Hist.F90")
    delta = history.split("f_trc_delta_evap_", 1)[1].split("f_trc_delta_soilevap_", 1)[0]
    assert "write_history_tracer_delta_2d" in delta
    assert "f_trc_evap_mass_" in delta
    mass = delta.split("f_trc_evap_mass_", 1)[1]
    assert "write_history_variable_2d" in mass
    assert "a_trc_evap(itrc_loc, :)" in mass
    assert "'tracer amount/m2'" in mass
    assert "write_history_tracer_delta_2d" not in mass


def test_vapor_forcing_retains_last_valid_value_instead_of_resetting_to_default_each_step():
    forcing = source("MOD_Tracer_Forcing.F90")
    prep = forcing.split("SUBROUTINE tracer_forcing_prepare_step", 1)[1].split(
        "END SUBROUTINE tracer_forcing_prepare_step", 1
    )[0]
    assert "trc_forc_precip_value(itrc, :) = tracer_precip_default_ratio(itrc)" not in prep
    assert "trc_forc_vapor_value(itrc, :) = tracer_vapor_default_ratio(itrc)" not in prep
    alloc = forcing.split("SUBROUTINE tracer_forcing_allocate_state", 1)[1].split(
        "END SUBROUTINE tracer_forcing_allocate_state", 1
    )[0]
    assert "trc_forc_vapor_value(itrc, :) = tracer_vapor_default_ratio(itrc)" in alloc
    accessor = forcing.split("FUNCTION tracer_forcing_vapor_value", 1)[1].split(
        "END FUNCTION tracer_forcing_vapor_value", 1
    )[0]
    assert "tracer_forcing_vapor_value = trc_forc_vapor_value(itrc, ipatch)" in accessor
    assert "IF (trc_forc_has_vapor" not in accessor


def test_vapor_forcing_diagnostics_split_near_dry_from_wet_invalid_missing():
    forcing = source("MOD_Tracer_Forcing.F90")
    assert "DECODE_NEAR_DRY" in forcing
    decode = forcing.split("SUBROUTINE tracer_forcing_decode_value", 1)[1].split(
        "END SUBROUTINE tracer_forcing_decode_value", 1
    )[0]
    assert re.search(r"total <= min_total.*?DECODE_NEAR_DRY", decode, re.DOTALL)
    log = forcing.split("SUBROUTINE tracer_forcing_log_ranges", 1)[1].split(
        "END SUBROUTINE tracer_forcing_log_ranges", 1
    )[0]
    assert "near-dry fallback patches=" in log
    assert "invalid/missing wet fallback patches=" in log
    assert "invalid/near-dry isotope forcing" not in log
    assert "previous valid/default vapor ratio" in log
    assert "trc_forc_vapor_warned_dry" in log
    assert "trc_forc_vapor_warned_invalid" in log
    has_vapor = forcing.split("FUNCTION tracer_forcing_has_vapor", 1)[1].split(
        "END FUNCTION tracer_forcing_has_vapor", 1
    )[0]
    assert "tracer_forcing_vapor_configured(itrc)" in has_vapor


def test_raw_over_total_forcing_requires_matching_temporal_config():
    forcing = source("MOD_Tracer_Forcing.F90")
    add_var = forcing.split("SUBROUTINE tracer_forcing_add_var", 1)[1].split(
        "END SUBROUTINE tracer_forcing_add_var", 1
    )[0]
    assert "raw/total forcing temporal config mismatch" in add_var
    assert "dtime /= trc_var_dtime(total_idx)" in add_var
    assert "offset /= trc_var_offset(total_idx)" in add_var
    assert "tracer_lower(tintalgo)" in add_var
    assert "trc_var_timelog(n_trc_forc_vars) = trc_var_timelog(total_idx)" in add_var


def test_precip_forcing_retains_last_valid_value_instead_of_resetting_to_default_each_step():
    forcing = source("MOD_Tracer_Forcing.F90")
    prep = forcing.split("SUBROUTINE tracer_forcing_prepare_step", 1)[1].split(
        "END SUBROUTINE tracer_forcing_prepare_step", 1
    )[0]
    assert "trc_forc_precip_value(itrc, :) = tracer_precip_default_ratio(itrc)" not in prep
    alloc = forcing.split("SUBROUTINE tracer_forcing_allocate_state", 1)[1].split(
        "END SUBROUTINE tracer_forcing_allocate_state", 1
    )[0]
    assert "trc_forc_precip_value(itrc, :) = tracer_precip_default_ratio(itrc)" in alloc
    accessor = forcing.split("FUNCTION tracer_forcing_precip_value", 1)[1].split(
        "END FUNCTION tracer_forcing_precip_value", 1
    )[0]
    assert "tracer_forcing_precip_value = trc_forc_precip_value(itrc, ipatch)" in accessor
    assert "IF (trc_forc_has_precip" not in accessor


def test_precip_forcing_diagnostics_split_near_dry_from_wet_invalid_missing():
    forcing = source("MOD_Tracer_Forcing.F90")
    assert "trc_forc_precip_status" in forcing
    decode = forcing.split("SUBROUTINE tracer_forcing_decode_value", 1)[1].split(
        "END SUBROUTINE tracer_forcing_decode_value", 1
    )[0]
    near_dry = decode.index("total <= min_total")
    raw_valid = decode.index("tracer_forcing_valid_value(raw)")
    assert near_dry < raw_valid, "near-dry total must win over invalid raw isotope numerator"
    log = forcing.split("SUBROUTINE tracer_forcing_log_ranges", 1)[1].split(
        "END SUBROUTINE tracer_forcing_log_ranges", 1
    )[0]
    assert "WARNING precip " in log
    assert "near-dry fallback patches=" in log
    assert "invalid/missing wet fallback patches=" in log
    assert "previous valid/default precip ratio" in log
    assert "trc_forc_precip_warned_dry" in log
    assert "trc_forc_precip_warned_invalid" in log
    has_precip = forcing.split("FUNCTION tracer_forcing_has_precip", 1)[1].split(
        "END FUNCTION tracer_forcing_has_precip", 1
    )[0]
    assert "tracer_forcing_precip_configured(itrc)" in has_precip
