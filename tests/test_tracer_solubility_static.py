from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"
SNOW_TOPOLOGY = ROOT / "main" / "MOD_SnowLayersCombineDivide.F90"
COLM_MAIN = ROOT / "main" / "CoLMMAIN.F90"


def source(name: str) -> str:
    return (TRACER / name).read_text(encoding="utf-8")


def test_solid_pools_have_location_preserving_state_lifecycle():
    variables = source("MOD_Tracer_Vars.F90")
    for name in (
        "trc_solid_soisno",
        "trc_canopy_solid",
        "trc_surface_solid",
        "trc_subsurface_solid",
        "trc_waterstorage_solid",
    ):
        assert f"allocate({name}" in variables
        assert f"allocated({name})" in variables
        assert f"lulcc_{name}_old" in variables

    assert "CALL remap3d_mass(lulcc_trc_solid_soisno_old" in variables
    for name in (
        "trc_canopy_solid",
        "trc_surface_solid",
        "trc_subsurface_solid",
        "trc_waterstorage_solid",
    ):
        assert f"CALL remap2d_mass(lulcc_{name}_old" in variables


def test_solid_pools_are_restartable_conserved_inventories():
    restart = source("MOD_Tracer_Rest.F90")
    conservation = source("MOD_Tracer_Conservation.F90")
    history = source("MOD_Tracer_Hist.F90")

    for name in (
        "trc_solid_soisno",
        "trc_canopy_solid",
        "trc_surface_solid",
        "trc_subsurface_solid",
        "trc_waterstorage_solid",
    ):
        assert f"'{name}'" in restart
        assert name in conservation
        assert name in history

    assert ".not. tracer_has_dissolved_limit(itrc)" in conservation


def test_all_land_transport_paths_use_the_shared_solubility_equilibrium():
    for filename in (
        "MOD_Tracer_Precip.F90",
        "MOD_Tracer_Evapo.F90",
        "MOD_Tracer_Snow.F90",
        "MOD_Tracer_SoilWater.F90",
        "MOD_Tracer_SpecialPatches.F90",
    ):
        assert "tracer_equilibrate_dissolved" in source(filename)


def test_special_patch_ratios_use_the_physical_water_floor():
    special = source("MOD_Tracer_SpecialPatches.F90")
    for expression in (
        "water_before_output > trc_tiny",
        "water_after_evap > trc_tiny",
        "water_beg > trc_tiny",
    ):
        assert expression not in special


def test_snow_topology_moves_solid_inventory_with_liquid_and_ice_tracer():
    topology = SNOW_TOPOLOGY.read_text(encoding="utf-8")
    colm_main = COLM_MAIN.read_text(encoding="utf-8")

    for routine in (
        "SUBROUTINE snowlayerscombine",
        "SUBROUTINE snowlayersdivide",
        "SUBROUTINE SnowLayersCombine_snicar",
        "SUBROUTINE SnowLayersDivide_snicar",
    ):
        block = topology.split(routine, 1)[1].split("END SUBROUTINE", 1)[0]
        assert "trc_solid" in block

    assert colm_main.count("trc_solid = trc_solid_soisno") == 4


def test_near_dry_transport_denominators_use_the_physical_water_floor():
    restart = source("MOD_Tracer_Rest.F90")
    soilwater = source("MOD_Tracer_SoilWater.F90")

    assert "water_ref > trc_tiny" not in restart
    assert "surface_water > trc_tiny" not in restart
    assert "water_before_flow > trc_tiny" not in soilwater
    assert "wliq_soisno_bef(j) > trc_tiny" not in soilwater.split(
        "SUBROUTINE tracer_wetland", 1
    )[1]
