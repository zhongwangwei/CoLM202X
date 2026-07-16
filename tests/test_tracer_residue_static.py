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


def test_residue_restart_is_optional_and_round_trips():
    restart = source("MOD_Tracer_Rest.F90")
    required = restart.split("Tracer restart ntracers/soilsnow mismatch", 1)[0]
    for name in ("trc_surface_residue", "trc_subsurface_residue"):
        assert f"'{name}'" not in required
        assert f"IF (tracer_dim_matches(file_restart, '{name}'))" in restart
        assert f"CALL ncio_read_vector(file_restart, '{name}'" in restart
        assert f"CALL ncio_write_vector(file_restart, '{name}'" in restart


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
    history = source("MOD_Tracer_Hist.F90")
    assert "f_trc_surface_residue_" in history
    assert "f_trc_subsurface_residue_" in history
    assert "f_trc_layer_dry_inventory_" in history
    assert history.count("'tracer amount/m2'") == 4


def test_internal_transfer_identity():
    mobile_before = 3.0e-4
    residue_before = 2.0e-4
    mobile_after_dry = 0.0
    residue_after_dry = residue_before + mobile_before
    assert mobile_before + residue_before == mobile_after_dry + residue_after_dry

    mobile_after_rewet = mobile_after_dry + residue_after_dry
    residue_after_rewet = 0.0
    assert mobile_after_rewet + residue_after_rewet == residue_after_dry
