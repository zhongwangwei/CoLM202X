"""DEF_UnitCatchment_regional: wiring of the regional unit-catchment network.

mksrfdata cuts the routing network to the river systems that receive runoff from
the land domain; mkinidata and the model then read the cut file.  The cut file
renumbers unit catchments, so everything that is keyed by unit-catchment number
must either come from the cut file or be translated.  These tests pin that.
"""
from pathlib import Path
import os
import re

import pytest

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


# ------------------------------------------------------------------- namelist

def test_switch_is_off_by_default_and_reaches_every_rank():
    src = read("share/MOD_Namelist.F90")
    assert re.search(r"logical\s+::\s+DEF_UnitCatchment_regional\s*=\s*\.false\.", src)
    assert re.search(r"^\s*DEF_UnitCatchment_regional,\s*&", src, re.M), "missing from the namelist group"
    assert re.search(r"mpi_bcast \(DEF_UnitCatchment_regional\s+,1\s+,mpi_logical", src)


def test_switch_is_refused_where_the_regional_network_cannot_be_used():
    src = read("share/MOD_Namelist.F90")
    block = src.split("IF (DEF_UnitCatchment_regional) THEN", 1)[1].split("IF (.not. ieee_is_finite(DEF_simulation_time%timestep))", 1)[0]
    # CaMa-Flood, catchment-based flow and single point keep the full map: define.h
    # undefines GridRiverLakeFlow for all of them.
    assert re.search(r"#ifndef GridRiverLakeFlow.*?CoLM_Stop", block, re.S)
    assert re.search(r"#ifdef LULCC.*?CoLM_Stop", block, re.S)
    assert "trim(DEF_UnitCatchment_file) == 'null'" in block
    assert "REGIONAL_UNITCATCHMENT_SUFFIX" in block, "path length check"
    assert "GridRiverLakeFlow" in read("include/define.h")


def test_file_in_use_is_the_regional_file_only_when_requested():
    src = read("share/MOD_Namelist.F90")
    assert "REGIONAL_UNITCATCHMENT_SUFFIX = '/riverlake/unitcatchment_regional.nc'" in src
    body = src.split("FUNCTION get_unitcatchment_file", 1)[1].split("END FUNCTION get_unitcatchment_file", 1)[0]
    assert re.search(r"IF \(DEF_UnitCatchment_regional\) THEN\s+fname = regional_unitcatchment_file \(\)\s+ELSE\s+fname = DEF_UnitCatchment_file", body)


# ------------------------------------------------------------------ consumers

def code_lines(src):
    """Source lines without comments and string literals."""
    for line in src.splitlines():
        if not line.strip().startswith("!"):
            yield re.sub(r"'[^']*'", "''", line)


def test_routing_code_reads_the_network_file_in_use():
    for rel in ("main/HYDRO/MOD_Grid_RiverLakeNetwork.F90",
                "main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90",
                "main/TRACER/MOD_Tracer_Particle_Sediment.F90"):
        src = read(rel)
        assert "get_unitcatchment_file" in src, rel
        if rel.endswith("RiverLakeNetwork.F90"):
            # the one direct use: the check against the network the file was cut from
            outside = src.replace(src.split("SUBROUTINE verify_regional_network", 1)[1]
                                  .split("END SUBROUTINE verify_regional_network", 1)[0], "")
            src = outside
        assert not any("DEF_UnitCatchment_file" in line for line in code_lines(src)), rel


def test_cama_flood_keeps_reading_the_full_network():
    for rel in ("extends/CaMa/src/MOD_CaMa_colmCaMa.F90", "extends/CaMa/src/cmf_ctrl_maps_mod.F90"):
        src = read(rel)
        assert "DEF_UnitCatchment_file" in src
        assert "get_unitcatchment_file" not in src and "DEF_UnitCatchment_regional" not in src, rel


def test_network_build_refuses_a_regional_file_from_another_source():
    src = read("main/HYDRO/MOD_Grid_RiverLakeNetwork.F90")
    assert "IF (DEF_UnitCatchment_regional) CALL verify_regional_network (parafile, x_ucat, y_ucat)" in src
    body = src.split("SUBROUTINE verify_regional_network", 1)[1].split("END SUBROUTINE verify_regional_network", 1)[0]
    for needle in ("seq_src_index", "x_source(src_index) == x_regional", "y_source(src_index) == y_regional", "CoLM_stop"):
        assert needle in body, needle


def test_reservoir_catalogue_is_translated_before_it_is_matched():
    """dam_seq numbers the unit catchments of the full network; matching it against the
    renumbered ucat_ucid would silently attach dams to the wrong cells."""
    src = read("main/HYDRO/MOD_Grid_Reservoir.F90")
    body = src.split("SUBROUTINE reservoir_init", 1)[1]
    translate = body.index("IF (DEF_UnitCatchment_regional) THEN")
    assert "'seq_src_index'" in body[translate:translate + 600]
    assert translate < body.index("CALL quicksort (nresv_catalogue, dam_seq, order)")
    assert translate < body.index("find_in_sorted_list1 (ucat_ucid(i), nresv_catalogue, dam_seq)")
    # dams outside the region must never match and must stay distinct (duplicate check follows)
    assert "dam_seq(i) = -i" in body


# ------------------------------------------------------------------- mksrfdata

def test_mksrfdata_cuts_the_network_after_the_land_patches_exist():
    src = read("mksrfdata/MKSRFDATA.F90")
    assert "USE MOD_UnitCatchmentRegional, only: unitcatchment_regional_build" in src
    call = src.index("CALL unitcatchment_regional_build ()")
    assert re.search(r"IF \(DEF_UnitCatchment_regional\) THEN\s+CALL unitcatchment_regional_build", src)
    assert src.index("CALL landpatch_build(lc_year)") < call
    assert src.index("CALL landpft_build  (lc_year)") < call


def test_regional_step_finds_input_cells_the_way_the_model_does():
    reg = read("mksrfdata/MOD_UnitCatchmentRegional.F90")
    net = read("main/HYDRO/MOD_Grid_RiverLakeNetwork.F90")
    # same grid and same patch-to-grid mapping as build_riverlake_network
    assert "CALL gridro%define_by_ndims (nlon, nlat)" in reg
    assert "CALL build_worker_remapdata (landpatch, gridro, remap)" in reg
    assert "CALL build_worker_remapdata (landpatch, griducat, remap_patch2inpm)" in net
    # same cell numbering: id = (y-1)*nlon + x
    assert "idmap_gd2uc = (idmap_y-1)*nlon_ucat + idmap_x" in net
    assert "ix = mod(ids(i)-1, nlon) + 1" in reg and "iy = (ids(i)-1) / nlon + 1" in reg
    # bifurcation-joined river systems are always kept together
    assert "touched, .true., nkeep, nsystem" in reg


def test_new_modules_are_built_in_dependency_order():
    make = read("Makefile")
    objs = make.split("OBJS_MKSRFDATA = ", 1)[1].split("MKSRFDATA.o\n", 1)[0]
    assert objs.index("MOD_UnitCatchmentSubset.o") < objs.index("MOD_UnitCatchmentRegional.o")
    assert "MOD_UnitCatchmentRegional.o: MOD_UnitCatchmentSubset.o" in make
    assert "MKSRFDATA.o: MOD_UnitCatchmentRegional.o" in make


# ---------------------------------------------------------- what gets renumbered

FORTRAN_RENUMBERED = {"seq", "seq_next", "seq_upst", "bifurcation_upst", "bifurcation_down", "dam_seq"}


def test_fortran_and_reference_tool_renumber_the_same_variables():
    fortran = read("mksrfdata/MOD_UnitCatchmentSubset.F90")
    listed = re.search(r"CASE \(([^)]*'dam_seq'[^)]*)\)", fortran).group(1)
    assert set(re.findall(r"'([^']+)'", listed)) == FORTRAN_RENUMBERED
    tool = read("tools/subset_unitcatchment.py")
    block = tool.split("index_vars = {", 1)[1].split("}", 1)[0]
    assert set(re.findall(r'"([^"]+)":', block)) == FORTRAN_RENUMBERED


GLOBAL_NETWORK = Path(os.environ.get(
    "COLM_UNITCATCHMENT_GLOBAL", "/Volumes/Data01/Data/CoLMruntime/unitcatchment/grid_routing_data_15min.nc"))


def test_every_variable_of_the_network_file_that_holds_cell_numbers_is_renumbered():
    """Audit of the real file: an integer variable whose long_name speaks of a
    'sequence' (unit-catchment number) but is not a grid index must be in the
    renumbered set.  Fails when a future network file adds such a variable."""
    netCDF4 = pytest.importorskip("netCDF4")
    if not GLOBAL_NETWORK.is_file():
        pytest.skip(f"global network file not available: {GLOBAL_NETWORK}")
    with netCDF4.Dataset(GLOBAL_NETWORK) as f:
        numbered = set()
        for name, var in f.variables.items():
            text = getattr(var, "long_name", "").lower()
            if var.dtype.kind in "iu" and "sequence" in text and "x-index" not in text and "y-index" not in text:
                numbered.add(name)
        # seq_upn is a count, not a cell number
        numbered.discard("seq_upn")
    assert numbered == FORTRAN_RENUMBERED
