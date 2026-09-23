"""Regression guards for the defects found in the two-day review.

These are source-structure checks, so each one states the failure it prevents.
Where a behaviour can be compiled and run cheaply it is (see the last tests).
"""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


def routine(src, name):
    match = re.search(
        rf"^\s*SUBROUTINE {name}\b.*?^\s*END SUBROUTINE {name}\b",
        src, re.M | re.S | re.I,
    )
    assert match, f"routine {name} not found"
    return match.group()


def indent(line):
    return len(line) - len(line.lstrip())


def enclosing_openers(lines, idx):
    """Block-opening lines that enclose lines[idx], innermost first.

    Fortran here is consistently indented, so walking upwards and keeping each
    line whose indentation is strictly smaller than everything seen so far
    yields the chain of IF / DO / SELECT statements around the target.
    """
    chain = []
    limit = indent(lines[idx])
    for i in range(idx - 1, -1, -1):
        text = lines[i]
        if not text.strip() or text.lstrip().startswith("!"):
            continue
        if indent(text) < limit:
            chain.append(text.strip())
            limit = indent(text)
    return chain


def calls_with_guard(body, call_text, guard):
    lines = body.splitlines()
    hits = [i for i, l in enumerate(lines) if call_text in l]
    assert hits, f"no call to {call_text}"
    return [any(guard in opener for opener in enclosing_openers(lines, i)) for i in hits]


# ---------------------------------------------------------------- river history

def test_route_writers_keep_the_master_out_of_the_group_collective():
    """The master is neither IO nor worker: entering route_shard_write_vector
    stops the run on 'layout not built' (route_hist_write_ucat/resv are called
    by every rank)."""
    src = read("main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90")
    for name in ("route_hist_write_ucat", "route_hist_write_resv"):
        body = routine(src, name)
        guarded = calls_with_guard(body, "CALL route_shard_write_vector",
                                   "IF (p_is_io .or. p_is_worker) THEN")
        assert guarded and all(guarded), f"{name}: unguarded shard write"


def test_bif_shard_dimensions_are_defined_per_file_not_per_run():
    """rh_bif_layout_built is only reset in route_hist_final, so dimensions tied
    to it exist only in the first history file of a run."""
    src = read("main/HYDRO/MOD_Grid_RiverLakeHistRoute.F90")
    body = routine(src, "route_hist_write_bif_matrix")
    lines = body.splitlines()
    calls = [i for i, l in enumerate(lines) if "CALL define_bif_shard_dims" in l]
    assert len(calls) == 1
    chain = enclosing_openers(lines, calls[0])
    assert not any("rh_bif_layout_built" in c for c in chain)
    assert any("ncio_var_exist" in c and "pth_global_id" in c for c in chain)
    assert "rh_bif_dims_file = rh_file_shard" in body
    assert "rh_bif_dims_file = ''" in routine(src, "route_hist_final")


def test_spinup_early_return_flushes_river_history_accumulators():
    """Routing accumulates every step; spinup skips hist_grid_riverlake_out,
    which is the only other place that clears them."""
    src = read("main/MOD_Hist.F90")
    early = src.split("IF (itstamp <= ptstamp) THEN", 1)[1].split("RETURN", 1)[0]
    assert "#ifdef GridRiverLakeFlow" in early
    assert "CALL flush_acc_fluxes_riverlake ()" in early
    hist = read("main/HYDRO/MOD_Grid_RiverLakeHist.F90")
    assert "PUBLIC :: flush_acc_fluxes_riverlake" in hist


def test_start_of_run_restart_is_not_reread_by_later_hist_init():
    flow = read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90")
    init = routine(flow, "grid_riverlake_flow_init")
    assert init.rstrip().endswith("END SUBROUTINE grid_riverlake_flow_init")
    assert init.index("gridriver_restart_file = ''") > init.index("CALL read_tracer_restart")


# ---------------------------------------------------- empty workers (numucat=0)

def test_reader_gives_empty_workers_zero_length_arrays():
    src = read("main/HYDRO/MOD_Grid_RiverLakeNetwork.F90")
    body = routine(src, "readin_riverlake_parameter")
    worker = body.split("ELSEIF (p_is_worker) THEN", 1)[1]
    assert "ELSE" in worker
    for line in ("allocate (rdata1d (0))", "allocate (rdata2d (ndim1,0))",
                 "allocate (idata1d (0))"):
        assert line in worker


def test_sediment_restart_write_does_not_touch_reader_arrays_on_empty_workers():
    src = read("main/TRACER/MOD_Tracer_Particle_Sediment.F90")
    body = routine(src, "write_sediment_restart")
    lines = body.splitlines()
    for needle in ("vector_gather_and_write(sed_frc(", "vector_gather_and_write(sed_slope(",
                   "vector_gather_and_write(topo_rivwth", "vector_gather_and_write(topo_rivlen"):
        idxs = [i for i, l in enumerate(lines) if needle in l]
        assert idxs, needle
        for i in idxs:
            assert any("p_is_worker .and. numucat > 0" in c
                       for c in enclosing_openers(lines, i)), needle


def test_single_point_crown_arrays_are_guarded_and_released():
    src = read("mksrfdata/MOD_SingleSrfdata.F90")
    for name in ("SITE_ncd_pfts", "SITE_ncw_pfts", "SITE_bcw_pfts"):
        assert f"IF (allocated({name})) deallocate({name})" in src
        assert re.search(rf"IF \(allocated\({name}\s*\)\) deallocate\({name}\s*\)", src)


# ------------------------------------------------------------------------ CaMa

def test_pthflw_is_mapped_to_river_cells_before_the_cell_writer():
    src = read("extends/CaMa/src/MOD_CaMa_Vars.F90")
    body = routine(src, "hist_out_cama")
    assert "real(D1PTHFLW_oAVG)" not in body
    assert "PTH_UPST(ipth)" in body and "sum(D1PTHFLW_oAVG(ipth,:))" in body
    assert "allocate (pthflw_cell(NSEQMAX,1))" in body
    # a corrupt pathway id must stop the run, not write outside the array
    assert "PTH_UPST(ipth) < 1 .or. PTH_UPST(ipth) > NSEQMAX" in body


def test_cama_history_attributes_follow_the_variable_not_record_one():
    """A record whose data was skipped still advances the time axis, so the
    first record a variable receives need not be record 1."""
    src = read("extends/CaMa/src/MOD_CaMa_Vars.F90")
    assert src.count("var_is_new = .not. ncio_var_exist (file_hist, varname, readflag = .false.)") == 3
    assert src.count("IF (itime_in_file == 1 .or. var_is_new) THEN") == 3
    assert "IF (itime_in_file == 1) THEN" not in src


def test_single_precision_cama_gets_conversion_specifics():
    src = read("extends/CaMa/src/cmf_coupling_budget_mod.F90")
    for name in ("budget_publish", "budget_runoff", "budget_debit"):
        assert f"INTERFACE {name}" in src
        assert f"SUBROUTINE {name}_r8(" in src
        assert f"SUBROUTINE {name}_rb(" in src


def test_budget_module_works_in_both_cama_precisions(tmp_path):
    fc = require_runnable_fortran_compiler(tmp_path)
    driver = tmp_path / "drv.f90"
    driver.write_text(
        """
program d
  use cmf_coupling_budget_mod
  use parkind1, only: jprb
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  integer :: x(1,1), y(1,1)
  real(real64) :: w(1,1), area(1,1), depth(1,1), frac(1,1), e(1,1), inf(1,1), er(1), ir(1)
  real(jprb) :: vol(1), fa(1), st(1), flow(1)
  real(real64) :: gflow(1,1)
  x=1; y=1; w=1; area=100
  call budget_init(x,y,w,area)
  vol=10; fa=50
  call budget_publish(vol,fa,depth,frac)
  if (abs(depth(1,1)-0.1_real64) > 1.e-6_real64) stop 1
  gflow=3
  call budget_runoff(gflow,flow)
  if (abs(flow(1)-3) > 1.e-5) stop 2
  e=1; inf=0; st=10
  call budget_debit(e,inf,st,er,ir)
  if (abs(st(1)-9) > 1.e-5) stop 3
  print *, 'BUDGET_OK'
end program
"""
    )
    src = ROOT / "extends/CaMa/src"
    for define in ([], ["-DSinglePrec_CMF"]):
        work = tmp_path / ("single" if define else "double")
        work.mkdir()
        cmd = [fc, "-cpp", *define, str(src / "parkind1.F90"),
               str(src / "cmf_coupling_budget_mod.F90"), str(driver), "-o", "drv"]
        built = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=120)
        assert built.returncode == 0, built.stdout + built.stderr
        ran = subprocess.run(["./drv"], cwd=work, capture_output=True, text=True, timeout=30)
        assert ran.returncode == 0 and "BUDGET_OK" in ran.stdout, ran.stdout + ran.stderr


# ------------------------------------------------------------------- interception

def _xsc_lines(src):
    body = re.search(r"SUBROUTINE LEAF_interception_CoLM2014 \(.*?END SUBROUTINE LEAF_interception_CoLM2014",
                     src, re.S).group()
    return body


def test_extended_interception_reports_the_same_old_pool_release_as_main():
    main = _xsc_lines(read("main/MOD_LeafInterception.F90"))
    ext = _xsc_lines(read("extends/interception/MOD_LeafInterception_Extended.F90"))
    for needle in ("xsc_rain_out          = xsc_rain / deltim",
                   "xsc_snow_out          = xsc_snow / deltim",
                   "xsc_rain = 0._r8", "xsc_snow = 0._r8"):
        assert needle in main and needle in ext, needle
    # module-level state in the extended file must be reset before first use
    assert ext.index("xsc_rain = 0._r8") < ext.index("xsc_rain  = max(0., ldew_rain-satcap_rain)")


# ----------------------------------------------------------------------- tracer

def test_tracer_history_numerators_are_persisted_with_their_denominator():
    """acctime_ucat is restored from the water restart; without the numerators a
    mid-window restart divides a partial sum by a full-window time."""
    src = read("main/TRACER/MOD_Tracer_RiverLake.F90")
    write = routine(src, "write_tracer_restart")
    read_ = routine(src, "read_tracer_restart")
    for prefix in ("trc_hist_stor_", "trc_hist_levsto_", "trc_hist_out_", "trc_hist_bifout_"):
        assert f"'{prefix}'" in write and f"'{prefix}'" in read_, prefix
    for name in ("trc_hist_water_storage", "trc_hist_levsto_water"):
        assert f"'{name}'" in write and f"'{name}'" in read_, name
    # optional on read (older restarts), collective on write (every rank enters)
    assert "CALL probe_riverlake_restart_vector (file_restart, varname, .false., on_disk)" in src
    # written before the commit marker, so a torn write stays uncommitted
    assert write.index("trc_hist_water_storage") < write.rindex("'trc_river_restart_complete', 1)")


def test_dead_tracer_code_stays_removed():
    src = read("main/TRACER/MOD_Tracer_RiverLake.F90")
    flow = read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90")
    assert "tracer_refresh_state" not in src and "tracer_refresh_state" not in flow
    assert "--- 10. Final concentration" not in src
    # only CoLMDEBUG resets and reads this accumulator, so only it may write it
    assert re.search(
        r"#ifdef CoLMDEBUG\n(?:\s*!.*\n)*\s*IF \(allocated\(trc_reactive_source\)\) THEN\n"
        r"\s*trc_reactive_source\(itrc, i\) = trc_reactive_source\(itrc, i\) \+ reactive_src\n"
        r"\s*ENDIF\n#endif", src)


def test_reservoir_history_stays_current_on_inland_depression_cells():
    src = read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90")
    block = src.split("IF (ucat_next(i) == -10) THEN\n                     ! Inland depression", 1)[1]
    block = block.split("CYCLE", 1)[0]
    assert "qresv_in(irsv)  = - sum_hflux_riv(i)" in block
    assert "qresv_out(irsv) = 0._r8" in block


def test_sediment_substep_diagnostic_scans_are_debug_only():
    src = read("main/TRACER/MOD_Tracer_Particle_Sediment.F90")
    assert re.search(
        r"#ifdef CoLMDEBUG\n(?:\s*!.*\n)+\s*IF \(numucat > 0\) THEN\n"
        r"\s*max_sedcon_local = max\(max_sedcon_local, maxval\(sedcon\)\)", src)


def _preprocessor_depth_at(text, needle):
    depth = 0
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith(("#if", "#ifdef", "#ifndef")):
            depth += 1
        elif stripped.startswith("#endif"):
            depth -= 1
        elif needle in line and not stripped.startswith("!"):
            return depth
    raise AssertionError(f"{needle} not found")


def test_sediment_diagnostic_packing_reads_only_defined_variables():
    """diag_max_global / diag_count_global are packed without CoLMDEBUG, so the
    variables they read must be set without it as well."""
    src = read("main/TRACER/MOD_Tracer_Particle_Sediment.F90")
    body = routine(src, "grid_sediment_calc")
    assert "diag_max_global = (/ max_sedcon_local" in body
    for name in ("max_sedcon_local", "max_sedout_local", "max_bedout_local",
                 "max_sedinp_local", "max_netflw_local", "max_shearvel_local",
                 "max_es_raw_local", "max_d_raw_local", "max_es_eff_local",
                 "max_d_eff_local"):
        assert _preprocessor_depth_at(body, f"{name} = 0._r8") == 0, name
    assert _preprocessor_depth_at(body, "n_flow_cancel_local = 0") == 0


def test_sediment_phase_timing_is_debug_only():
    """clk_rate is only set under CoLMDEBUG, so every read of it (and the system
    clock calls in the hot substep loop) must be under CoLMDEBUG too."""
    src = read("main/TRACER/MOD_Tracer_Particle_Sediment.F90")
    body = routine(src, "grid_sediment_calc")
    depth = 0
    seen = 0
    for line in body.splitlines():
        stripped = line.strip()
        if stripped.startswith(("#if", "#ifdef", "#ifndef")):
            depth += 1
        elif stripped.startswith("#endif"):
            depth -= 1
        elif ("clk_rate" in line or "clk_phase" in line) \
                and not stripped.startswith("!") and not stripped.startswith("integer"):
            seen += 1
            assert depth > 0, f"unconditional timing statement: {stripped}"
    assert seen >= 18


# ------------------------------------------- CaMa-aligned sediment / routing choices

def test_particle_inputs_are_sized_per_period_and_reuse_the_substep_flood_area():
    src = read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90")
    # allocated next to total_floodarea (once per routing period), freed at its end
    assert "allocate (particle_floodarea (numucat))" in src
    assert "IF (allocated(particle_floodarea)) deallocate(particle_floodarea)" in src
    assert "allocate(particle_floodarea(" not in src
    block = src.split("particle_floodarea(i) = total_floodarea(i)", 1)
    assert len(block) == 2, "particle flood area must reuse total_floodarea"
    assert "floodplain_curve(i)%floodarea" not in src.split("CALL tracer_lifecycle_route_diag_accumulate", 1)[0].rsplit(
        "IF (tracer_lifecycle_route_has_active()) THEN", 1)[1]


def test_sediment_advection_uses_interval_donor_limit_not_dry_cell_cfl():
    src = read("main/TRACER/MOD_Tracer_Particle_Sediment.F90")
    assert "dt_cell = sed_cfl_adv" not in src
    assert "CALL calc_sediment_advection(dt_morph" in src
    assert "sedout(ised,i) = min(sedout(ised,i), avail_sto(ised,i) / dt)" in src

def test_sediment_yield_threshold_is_applied_per_forcing_step():
    src = read("main/TRACER/MOD_Tracer_Particle_Sediment.F90")
    put = routine(src, "sediment_forcing_put")
    assert "IF (precip(i) * 86400._r8 > SED_PRECIP_THRESHOLD_MM_DAY) THEN" in put
    # exposure time accrues for every valid step, thresholded or not
    assert put.index("sed_precip_time(i) = sed_precip_time(i) + weight") > put.index("ENDIF")
    assert "SED_PRECIP_THRESHOLD_MM_DAY" not in routine(src, "calc_sediment_yield")


def test_routing_clip_is_covered_by_the_always_on_closure_check():
    """The clip creates water only into totalvol_aft; the closure check that
    reports it must not be compiled out."""
    body = routine(read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90"), "grid_riverlake_flow")
    head, tail = body.split("water_balance_err = totalvol_aft - totalvol_bef - totalrnof + totaldis", 1)
    assert "WARNING grid_riverlake_flow: water balance residual=" in tail
    # depth 0 inside the routine: not under CoLMDEBUG (or any other switch)
    assert _preprocessor_depth_at(body, "water_balance_err = totalvol_aft - totalvol_bef") == 0


# --------------------------------------------------------------------- namelist

def test_routing_max_dt_is_validated():
    src = read("share/MOD_Namelist.F90")
    assert "ieee_is_finite(DEF_GRIDBASED_ROUTING_MAX_DT)" in src
    assert "DEF_GRIDBASED_ROUTING_MAX_DT <= 0._r8" in src


# ------------------------------------------------- TRACER built with no tracers

def test_period_end_reset_skips_tracer_arrays_that_were_never_allocated():
    """TRACER compiled in but DEF_TRACER_NUM = 0: river_lake_tracer_init returns
    early, so acc_trc_inp / acc_rnof_ref / trc_dry_drain stay unallocated and
    the unconditional whole-array reset segfaulted at the first period end."""
    src = read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90")
    init = read("main/TRACER/MOD_Tracer_RiverLake.F90")
    assert re.search(r"IF \(ntracers <= 0\) RETURN", init), "premise: early return without tracers"
    for name in ("acc_trc_inp", "acc_rnof_ref", "trc_dry_drain"):
        assert not re.search(rf"^\s*{name} = 0\._r8", src, re.M), f"{name} is reset unguarded"
        assert re.search(rf"IF \(allocated\({name}\)\) {name} = 0\._r8", src), name


def test_tracer_history_vector_writers_import_elm_patch_for_unstructured():
    """MOD_Tracer_Hist did not compile with -DUNSTRUCTURED (no elm_patch)."""
    src = read("main/TRACER/MOD_Tracer_Hist.F90")
    for name in ("write_history_tracer_vector_2d", "write_history_tracer_ratio_vector_3d"):
        body = routine(src, name)
        branch = body.split("#else", 1)[1].split("#endif", 1)[0]
        assert "USE MOD_LandPatch, only: elm_patch" in branch, name


# ------------------------------------------------------------ push hot loops

def test_push_mapping_mode_is_compared_once_not_per_cell():
    """`trim(mode) == 'average'` inside the per-cell loop heap-allocates a
    temporary string for every cell of every substep (about 10% of a worker's
    time in the river-routing profile)."""
    src = read("share/MOD_WorkerPushData.F90")
    assert len(re.findall(r"do_average = \(trim\((?:batch_)?mode\) == 'average'\)", src)) == 3
    for line in src.splitlines():
        if line.strip().startswith("!"):
            continue
        if "trim(" in line and "'average'" in line:
            assert line.strip().startswith("do_average ="), line


# -------------------------------------------------------------- runtime: readers

def test_empty_worker_arrays_survive_bounds_checking(tmp_path):
    """The failure the reader fix prevents: a whole-array WHERE on an array that
    was never allocated aborts under -fbounds-check (Makeoptions.Mac-arm) and
    -fcheck=all (Makeoptions.github)."""
    fc = require_runnable_fortran_compiler(tmp_path)
    program = tmp_path / "t.f90"
    program.write_text(
        """
program t
  implicit none
  real(8), allocatable :: topo_area(:), prcp_area_uc(:)
  allocate(prcp_area_uc(0))
  allocate(topo_area(0))          ! what the reader now does on an empty worker
  where (topo_area > 0._8)
     prcp_area_uc = prcp_area_uc / topo_area
  elsewhere
     prcp_area_uc = 0._8
  end where
  print *, 'EMPTY_OK'
end program
"""
    )
    for flags in (["-fbounds-check"], ["-fcheck=all"]):
        built = subprocess.run([fc, *flags, str(program), "-o", "t"], cwd=tmp_path,
                               capture_output=True, text=True, timeout=60)
        assert built.returncode == 0, built.stderr
        ran = subprocess.run(["./t"], cwd=tmp_path, capture_output=True, text=True, timeout=30)
        assert ran.returncode == 0 and "EMPTY_OK" in ran.stdout, ran.stdout + ran.stderr
