from pathlib import Path
import re
import subprocess

import pytest

from fortran_test_support import SMOKE_TIMEOUT, require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]


def _flat(path: Path) -> str:
    return re.sub(r"\s+", " ", re.sub(r"&\s*", " ", path.read_text().lower()))


def test_cama_feedback_uses_fraction_once_and_scalar_extraction() -> None:
    source = _flat(ROOT / "extends" / "CaMa" / "src" / "cmf_calc_stonxt_mod.F90")

    assert "dwevapex = min(real(p2fldsto(iseq,1),kind=jprb),dt*d2wevap(iseq,1))" in source
    assert "dwinfiltex = min(real(p2fldsto(iseq,1),kind=jprb),dt*d2winfilt(iseq,1))" in source
    assert "d2fldfrc(iseq,1)*dt*d2wevap" not in source
    assert "d2fldfrc(iseq,1)*dt*d2winfilt" not in source
    assert "d2winfiltex = min(" not in source


def test_cama_infiltration_arrays_are_initialized_and_averaged() -> None:
    init = _flat(ROOT / "extends" / "CaMa" / "src" / "cmf_ctrl_vars_mod.F90")
    diag = _flat(ROOT / "extends" / "CaMa" / "src" / "cmf_calc_diag_mod.F90")

    assert "allocate( d2winfilt(nseqmax,1) ) d2winfilt(:,:)=0._jprb" in init
    assert "allocate(d2winfiltex_aavg(nseqmax,1))" in init
    assert "d2winfiltex_aavg(iseq,1) = 0._jprb" in diag
    assert "d2winfiltex_aavg(iseq,1) = d2winfiltex_aavg(iseq,1) / real(nadd_adp,kind=jprb)" in diag


def test_cama_to_colm_mapping_preserves_flood_water_volume() -> None:
    source = _flat(ROOT / "extends" / "CaMa" / "src" / "MOD_CaMa_colmCaMa.F90")

    assert "p2fldsto" in source
    assert "d2grarea" in source
    assert "budget_publish(p2fldsto(:,1),d2fldfrc(:,1)*d2grarea(:,1),flddepth_tmp,fldfrc_tmp)" in source
    assert "flood_credit=max(0._r8,flddepth_cama)" in source
    assert "flddepth_cama=flood_credit/fldfrc_cama*1000._r8" in source
    assert "fldfrc_cama=fldfrc_cama/100" not in source

    fractions = (0.0, 0.2)
    depths = (0.0, 1.0)
    remapped_fraction = sum(fractions) / 2.0
    remapped_water = sum(f * d for f, d in zip(fractions, depths)) / 2.0
    assert remapped_water / remapped_fraction == 1.0


def test_flood_sinks_are_limited_by_available_water_without_five_percent_cutoff() -> None:
    colm = _flat(ROOT / "main" / "CoLMMAIN.F90")
    hydro = _flat(ROOT / "main" / "MOD_SoilSnowHydrology.F90")
    thermal = _flat(ROOT / "extends" / "interception" / "MOD_Thermal_CanopyPhase_Extended.F90")

    assert "fldfrc .gt. 0.05" not in colm
    assert "fldfrc .gt. 0.05" not in hydro
    assert "fldfrc .gt. 0.05" not in thermal
    assert "fevpg_fld_local = min(max(fevpg_fld_local,0._r8), flddepth/deltim)" in thermal
    assert "condensation over" in thermal and "not credited" in thermal
    assert hydro.count("qinfl_fld_subgrid = min(gfld,max(0._r8,gfld-rsur_fld))") == 2
    assert hydro.count("flddepth = max(0._r8,flddepth-deltim*qinfl_fld_subgrid)") == 2


def test_cama_extraction_history_units_are_volumetric_rates() -> None:
    source = _flat(ROOT / "extends" / "CaMa" / "src" / "MOD_CaMa_Vars.F90")

    assert "'wevap', itime_in_file,'inundation water evaporation','m3/s'" in source
    assert "'winfilt', itime_in_file,'inundation water infiltration','m3/s'" in source


def test_bidirectional_extraction_conserves_water_volume_numerically() -> None:
    fractions = (0.25, 0.75)
    conditional_depths_m = (0.4, 1.2)
    areas_m2 = (2.0, 6.0)
    extraction_m_per_s = (0.01, 0.02)
    dt = 10.0

    initial = sum(f * d * a for f, d, a in zip(fractions, conditional_depths_m, areas_m2))
    requested = sum(f * q * a * dt for f, q, a in zip(fractions, extraction_m_per_s, areas_m2))
    remaining = sum(
        f * max(0.0, d - q * dt) * a
        for f, d, q, a in zip(fractions, conditional_depths_m, extraction_m_per_s, areas_m2)
    )

    assert remaining == pytest.approx(initial - requested)


def test_flood_evaporation_is_mixed_before_ground_temperature_and_excluded_from_land_water() -> None:
    colm = _flat(ROOT / "main" / "CoLMMAIN.F90")
    thermal = _flat(ROOT / "extends" / "interception" / "MOD_Thermal_CanopyPhase_Extended.F90")

    assert thermal.index("call get_fldevp") < thermal.index("call groundtemperature")
    assert "cgrnd = cgrnd_land*(1._r8-fldfrc_eff)" in thermal
    assert "fevpg_wat = fevpg - (hvap/htvp)*fevpg_fld" in thermal
    assert "qseva = min(wliq_soisno(lb)/deltim, fevpg_wat)" in thermal
    assert "fevpg = fevpg_fld + fevpg_wat" in thermal
    assert "lfevpg_ground = htvp*fevpg_wat + hvap*fevpg_fld" in thermal
    assert "lfevpa = lfevpl + lfevpg_ground" in thermal
    assert "fevpa_wb = fevpa - fevpg_fld" in colm
    assert "flood_input_wb = qinfl_fld" in colm
    assert "forc_prc+forc_prl+flood_input_wb-fevpa_wb-rnof" in colm
    assert "water_input_in = (forc_prc + forc_prl + flood_input_wb) * deltim" in colm
    assert "water_evap_in = fevpa_wb * deltim" in colm
    assert "cama flood infiltration has no composition for active water-borne tracers" in colm
    assert "call get_fldevp" not in colm


def test_flood_wet_dry_partition_numerically_conserves_water_and_energy() -> None:
    # Mirrors THERMAL's algebra: wet flood flux is explicit (zero derivative),
    # dry-land latent derivative is scaled by the unflooded fraction, and only
    # the dry component is exposed to WATER.
    fldfrc = 0.30
    dry = 1.0 - fldfrc
    dt = 1800.0
    tinc = 2.0
    hvap = 2.501e6
    hsub = 2.834e6

    dry_fevpg = 2.0e-5
    dry_cgrndl = 1.0e-6
    flood_raw = -1.0e-5  # dew over floodwater: policy is no CaMa credit.
    flood_sink = min(max(flood_raw, 0.0), 0.010 / dt) * fldfrc

    boundary_equiv = dry_fevpg * dry + (hvap / hsub) * flood_sink
    scaled_deriv = dry_cgrndl * dry
    corrected_equiv = boundary_equiv + tinc * scaled_deriv
    land_for_water = corrected_equiv - (hvap / hsub) * flood_sink
    physical_ground_evap = land_for_water + flood_sink
    latent_ground = hsub * land_for_water + hvap * flood_sink

    assert flood_sink == 0.0
    assert land_for_water == corrected_equiv
    assert physical_ground_evap == pytest.approx(dry * (dry_fevpg + tinc * dry_cgrndl))
    assert latent_ground == pytest.approx(hsub * dry * (dry_fevpg + tinc * dry_cgrndl))

    flood_raw = 4.0e-5
    flood_sink = min(max(flood_raw, 0.0), 0.010 / dt) * fldfrc
    boundary_equiv = dry_fevpg * dry + (hvap / hsub) * flood_sink
    corrected_equiv = boundary_equiv + tinc * scaled_deriv
    land_for_water = corrected_equiv - (hvap / hsub) * flood_sink
    physical_ground_evap = land_for_water + flood_sink
    latent_ground = hsub * land_for_water + hvap * flood_sink

    assert flood_sink > 0.0
    assert land_for_water == pytest.approx(dry * (dry_fevpg + tinc * dry_cgrndl))
    assert physical_ground_evap == land_for_water + flood_sink
    assert latent_ground == hsub * land_for_water + hvap * flood_sink


def test_flood_wet_dry_partition_fortran_block_regression(tmp_path: Path) -> None:
    """Compile/run a Fortran harness for the THERMAL CaMa flux formulas."""
    thermal = _flat(ROOT / "extends" / "interception" / "MOD_Thermal_CanopyPhase_Extended.F90")
    for formula in (
        "fevpg = (hvap/htvp)*fevpg_fld + fevpg_land*(1._r8-fldfrc_eff)",
        "cgrndl = cgrndl_land*(1._r8-fldfrc_eff)",
        "fevpg_wat = fevpg - (hvap/htvp)*fevpg_fld",
        "fevpg = fevpg_fld + fevpg_wat",
        "lfevpg_ground = htvp*fevpg_wat + hvap*fevpg_fld",
    ):
        assert formula in thermal

    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "flood_flux_block.f90"
    exe = tmp_path / "flood_flux_block"
    src.write_text(
        '''
program flood_flux_block
  implicit none
  integer, parameter :: r8 = selected_real_kind(12, 300)
  real(r8), parameter :: tol = 1.0e-12_r8
  real(r8), parameter :: hvap = 2.501e6_r8, htvp = 2.834e6_r8
  real(r8) :: fldfrc_eff, fevpg_fld, fevpg_fld_local, flddepth, deltim
  real(r8) :: fevpg_land, cgrndl_land, cgrnds_land, tinc_l, tinc_s
  real(r8) :: fevpg, fevpg_wat, cgrndl, cgrnds, lfevpg_ground
  real(r8) :: expected_land, expected_latent

  deltim = 1800._r8
  flddepth = 0.010_r8
  fldfrc_eff = 0.30_r8
  fevpg_land = 2.0e-5_r8
  cgrndl_land = 1.0e-6_r8
  cgrnds_land = 2.0e-6_r8
  tinc_l = 2.0_r8
  tinc_s = -1.0_r8

  fevpg_fld_local = -1.0e-5_r8
  fevpg_fld_local = min(max(fevpg_fld_local,0._r8), flddepth/deltim)
  fevpg_fld = fevpg_fld_local*fldfrc_eff
  if (abs(fevpg_fld) > tol) error stop 1
  call check_case()

  fevpg_fld_local = 4.0e-5_r8
  fevpg_fld_local = min(max(fevpg_fld_local,0._r8), flddepth/deltim)
  fevpg_fld = fevpg_fld_local*fldfrc_eff
  if (fevpg_fld <= 0._r8) error stop 2
  call check_case()

contains
  subroutine check_case()
    fevpg = fevpg_land*(1._r8-fldfrc_eff) + (hvap/htvp)*fevpg_fld
    cgrndl = cgrndl_land*(1._r8-fldfrc_eff)
    cgrnds = cgrnds_land*(1._r8-fldfrc_eff)
    fevpg = fevpg + tinc_l*cgrndl + tinc_s*cgrnds
    fevpg_wat = fevpg - (hvap/htvp)*fevpg_fld
    fevpg = fevpg_fld + fevpg_wat
    lfevpg_ground = htvp*fevpg_wat + hvap*fevpg_fld

    expected_land = (1._r8-fldfrc_eff) * (fevpg_land + tinc_l*cgrndl_land + tinc_s*cgrnds_land)
    expected_latent = htvp*expected_land + hvap*fevpg_fld
    if (abs(fevpg_wat - expected_land) > tol) error stop 3
    if (abs(fevpg - (expected_land + fevpg_fld)) > tol) error stop 4
    if (abs(lfevpg_ground - expected_latent) > 1.0e-7_r8) error stop 5
  end subroutine check_case
end program flood_flux_block
''',
        encoding="utf-8",
    )
    compiled = subprocess.run(
        [compiler, "-ffree-line-length-none", str(src), "-o", str(exe)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    ran = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True, timeout=SMOKE_TIMEOUT)
    assert ran.returncode == 0, ran.stdout + ran.stderr


def _fevpg_runoff_blocks() -> list[str]:
    source = (ROOT / "main/MOD_SoilSnowHydrology.F90").read_text()
    blocks = []
    start = 0
    while True:
        start = source.find("fevpg_runoff = fevpg(ipatch)", start)
        if start < 0:
            break
        end = source.find("qinfl_fld = 0._r8", start)
        assert end > start
        block = source[start:end]
        block = re.sub(r"^\s*#.*\n", "", block, flags=re.MULTILINE)
        block = re.sub(r"!.*", "", block)
        blocks.append(block.strip())
        start = end
    assert len(blocks) == 2
    return blocks


def test_vic_runoff_trial_uses_split_soil_snow_fluxes(tmp_path: Path) -> None:
    hydro = _flat(ROOT / "main/MOD_SoilSnowHydrology.F90")
    assert hydro.count("wliq_soisno(1:nl_soil), fevpg_runoff, rootflux") == 4

    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "vic_split_fevpg.f90"
    exe = tmp_path / "vic_split_fevpg"
    block1, block2 = _fevpg_runoff_blocks()
    src.write_text(
        """
program vic_split_fevpg
  implicit none
  integer, parameter :: r8 = selected_real_kind(12, 300)
  real(r8), parameter :: tol = 1.e-12_r8
  logical :: DEF_SPLIT_SOILSNOW
  integer :: patchtype, ipatch
  real(r8) :: fevpg_runoff, qseva, qsubl, qsdew, qfros
  real(r8) :: qseva_soil, qsubl_soil, qsdew_soil, qfros_soil
  real(r8) :: qseva_snow, qsubl_snow, qsdew_snow, qfros_snow
  real(r8) :: fevpg(1)

  ipatch = 1
  call run_cases_2014()
  call run_cases_vsf()

contains
  subroutine set_evap_case()
    fevpg = 99._r8
    qseva = 11._r8; qsubl = 13._r8; qsdew = 2._r8; qfros = 3._r8
    qseva_soil = 2._r8; qsubl_soil = 3._r8; qsdew_soil = 0.5_r8; qfros_soil = 0.25_r8
    qseva_snow = 5._r8; qsubl_snow = 7._r8; qsdew_snow = 1._r8; qfros_snow = 2._r8
  end subroutine set_evap_case

  subroutine set_condense_case()
    fevpg = 88._r8
    qseva = 1._r8; qsubl = 2._r8; qsdew = 5._r8; qfros = 7._r8
    qseva_soil = 0.5_r8; qsubl_soil = 0.25_r8; qsdew_soil = 2._r8; qfros_soil = 3._r8
    qseva_snow = 0.75_r8; qsubl_snow = 0.5_r8; qsdew_snow = 4._r8; qfros_snow = 5._r8
  end subroutine set_condense_case

  subroutine assert_close(got, want, code)
    real(r8), intent(in) :: got, want
    integer, intent(in) :: code
    if (abs(got - want) > tol) error stop code
  end subroutine assert_close

  subroutine run_cases_2014()
    call set_evap_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .true.
    call assign_block_2014()
    call assert_close(fevpg_runoff, 13.25_r8, 1)

    call set_condense_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .true.
    call assign_block_2014()
    call assert_close(fevpg_runoff, -12._r8, 2)

    call set_evap_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .false.
    call assign_block_2014()
    call assert_close(fevpg_runoff, 19._r8, 3)

    call set_condense_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .false.
    call assign_block_2014()
    call assert_close(fevpg_runoff, -9._r8, 4)

    call set_evap_case()
    patchtype = 1; DEF_SPLIT_SOILSNOW = .true.
    call assign_block_2014()
    call assert_close(fevpg_runoff, 99._r8, 5)
  end subroutine run_cases_2014

  subroutine run_cases_vsf()
    call set_evap_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .true.
    call assign_block_vsf()
    call assert_close(fevpg_runoff, 13.25_r8, 11)

    call set_condense_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .true.
    call assign_block_vsf()
    call assert_close(fevpg_runoff, -12._r8, 12)

    call set_evap_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .false.
    call assign_block_vsf()
    call assert_close(fevpg_runoff, 19._r8, 13)

    call set_condense_case()
    patchtype = 0; DEF_SPLIT_SOILSNOW = .false.
    call assign_block_vsf()
    call assert_close(fevpg_runoff, -9._r8, 14)

    call set_evap_case()
    patchtype = 1; DEF_SPLIT_SOILSNOW = .true.
    call assign_block_vsf()
    call assert_close(fevpg_runoff, 99._r8, 15)
  end subroutine run_cases_vsf

  subroutine assign_block_2014()
{block1}
  end subroutine assign_block_2014

  subroutine assign_block_vsf()
{block2}
  end subroutine assign_block_vsf
end program vic_split_fevpg
""".format(
            block1="\n".join(f"    {line}" for line in block1.splitlines()),
            block2="\n".join(f"    {line}" for line in block2.splitlines()),
        ),
        encoding="utf-8",
    )
    compiled = subprocess.run(
        [compiler, "-ffree-line-length-none", str(src), "-o", str(exe)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    ran = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True, timeout=SMOKE_TIMEOUT)
    assert ran.returncode == 0, ran.stdout + ran.stderr


def test_deprecated_closure_flags_do_not_advertise_unused_physics():
    source = (ROOT / 'extends/CaMa/src/cmf_ctrl_nmlist_mod.F90').read_text()
    assert 'LWEVAPFIX/LWINFILTFIX are deprecated and have no effect.' in source
    for path in (ROOT / 'run/CaMa').glob('*.nml'):
        for line in path.read_text().splitlines():
            if line.startswith(('LWEVAPFIX ', 'LWINFILTFIX ')):
                assert '.FALSE.' in line and 'no effect' in line
