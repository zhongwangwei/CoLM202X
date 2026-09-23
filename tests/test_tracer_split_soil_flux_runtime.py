"""Run the production soil-surface flux selector for all snow/split combinations."""

from pathlib import Path
import re
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]


def test_soil_surface_fluxes_follow_thermal_and_hydrology(tmp_path):
    source = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
    start = source.index("         IF (snl < 0 .and. .not. split_soilsnow) THEN")
    end = source.index("\n            ENDIF", start) + len("\n            ENDIF")
    selector = source[start:end]

    program = f"""
program check_split_soil_flux
  implicit none
  integer, parameter :: r8 = kind(1.d0)
  integer :: snl, isplit, isnows
  logical :: split_soilsnow
  real(r8) :: qseva_in, qsdew_in, qsubl_in, qfros_in
  real(r8) :: qseva_soil, qsdew_soil, qsubl_soil, qfros_soil
  real(r8) :: eff_qseva, eff_qsdew_topliq, eff_qsubl_top, eff_qfros_top
  real(r8) :: expected(4), actual(4)
  qseva_in = 0.1_r8
  qsdew_in = 0.2_r8
  qsubl_in = 0.3_r8
  qfros_in = 0.4_r8
  qseva_soil = 1.1_r8
  qsdew_soil = 1.2_r8
  qsubl_soil = 1.3_r8
  qfros_soil = 1.4_r8
  do isnows = 0, 1
    snl = -isnows
    do isplit = 0, 1
      split_soilsnow = isplit == 1
{selector}
      actual = [eff_qseva, eff_qsdew_topliq, eff_qsubl_top, eff_qfros_top]
      if (isnows == 1 .and. isplit == 0) then
        expected = 0._r8
      else if (isplit == 1) then
        expected = [qseva_soil, qsdew_soil, qsubl_soil, qfros_soil]
      else
        expected = [qseva_in, qsdew_in, qsubl_in, qfros_in]
      end if
      if (any(abs(actual - expected) > 1.e-12_r8)) stop 1
    end do
  end do
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_split_soil_flux.f90"
    exe = tmp_path / "check_split_soil_flux"
    src.write_text(program)
    subprocess.run([compiler, str(src), "-o", str(exe)], check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_vsf_soil_dew_respects_liquid_and_ice_pore_volume(tmp_path):
    """Execute the production post-solver dew allocation, not a copied formula."""
    source = (ROOT / "main/MOD_SoilSnowHydrology.F90").read_text()
    start = source.index("      IF (dew_input > 0._r8) THEN")
    end = source.index("      ENDIF", start) + len("      ENDIF")
    allocation = source[start:end]
    program = f"""
program check_vsf_dew
  implicit none
  integer, parameter :: r8 = kind(1.d0)
  real(r8), parameter :: denh2o = 1000._r8, denice = 917._r8
  real(r8) :: porsl(1), dz_soisno(1), wice_soisno(1), wliq_soisno(1)
  real(r8) :: dew_input, dew_capacity, dew_retained, dew_overflow, wdsrf
  real(r8) :: prior_water, prior_surface
  integer :: icase
  porsl = 0.4_r8
  dz_soisno = 0.03_r8
  do icase = 1, 3
    select case (icase)
    case (1) ! The failure mode: saturated soil plus condensation.
      wice_soisno = 0._r8
      wliq_soisno = 12._r8
      dew_input = 0.1_r8
    case (2) ! Ice excludes liquid by volume, not by ice mass.
      wice_soisno = 4._r8
      wliq_soisno = 7.6_r8
      dew_input = 0.1_r8
    case (3) ! Unsaturated top layer retains all condensation.
      wice_soisno = 0._r8
      wliq_soisno = 5._r8
      dew_input = 0.1_r8
    end select
    wdsrf = 0._r8
    dew_overflow = 0._r8
    prior_water = wliq_soisno(1)
    prior_surface = wdsrf
{allocation}
    if (abs((wliq_soisno(1)+wdsrf)-(prior_water+prior_surface+dew_input)) > 1.e-12_r8) stop 1
    if (wliq_soisno(1)/denh2o+wice_soisno(1)/denice > porsl(1)*dz_soisno(1)+1.e-14_r8) stop 2
    if (abs(dew_overflow-wdsrf) > 1.e-12_r8) stop 3
    if (icase == 1 .and. abs(dew_overflow-dew_input) > 1.e-12_r8) stop 4
    if (icase == 2 .and. dew_overflow <= 0._r8) stop 5
    if (icase == 3 .and. dew_overflow /= 0._r8) stop 6
  end do
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_vsf_dew.F90"
    exe = tmp_path / "check_vsf_dew"
    src.write_text(program)
    subprocess.run(
        [compiler, "-cpp", "-DTRACER", "-ffree-line-length-none", str(src), "-o", str(exe)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_frost_displaces_old_liquid_in_both_water_solvers(tmp_path):
    """Execute both production allocations, including frost without liquid dew."""
    source = (ROOT / "main/MOD_SoilSnowHydrology.F90").read_text()
    old_start = source.index("      dew_capacity = max((porsl(1)*dz_soisno(1)-wice_soisno(1)/denice)*denh2o")
    old_end = source.index("      IF (wdsrf > pondmx", old_start)
    vsf_start = source.index("      dew_capacity = max((porsl(1)*dz_soisno(1) - wice_soisno(1)/denice)*denh2o")
    vsf_end = source.index("      ! water imbalance mainly", vsf_start)
    allocations = (source[old_start:old_end], source[vsf_start:vsf_end])
    compiler = require_runnable_fortran_compiler(tmp_path)
    for name, allocation in zip(("2014", "vsf"), allocations):
        # WATER_2014 spells the dew diagnostic dew_overflow_trc; WATER_VSF
        # uses dew_overflow. Both distinguish it from frost_displaced.
        old = name == "2014"
        dew_name = "dew_overflow_trc" if old else "dew_overflow"
        frost_name = "frost_displaced_trc" if old else "frost_displaced"
        program = f"""
program check_frost
  implicit none
  integer, parameter :: r8 = kind(1.d0)
  real(r8), parameter :: denh2o = 1000._r8, denice = 917._r8
  real(r8) :: porsl(1), dz_soisno(1), wice_soisno(1), wliq_soisno(1)
  real(r8) :: dew_input, dew_capacity, dew_retained, dew_excess, frost_excess, wdsrf, ice_before_frost
  real(r8) :: dew_overflow, frost_displaced, dew_overflow_trc, frost_displaced_trc
  integer :: icase
  porsl = 0.4_r8
  dz_soisno = 0.03_r8
  do icase = 1, 3
    wdsrf = 0._r8
    ice_before_frost = 0._r8
    {dew_name} = 0._r8
    {frost_name} = 0._r8
    select case(icase)
    case(1) ! New frost expels old liquid, then dew arrives separately.
      wice_soisno = 2._r8
      wliq_soisno = 11._r8
      dew_input = 0.3_r8
    case(2) ! New frost without liquid dew must still expel old liquid.
      wice_soisno = 2._r8
      wliq_soisno = 11._r8
      dew_input = 0._r8
    case(3) ! With spare pores, liquid dew stays in soil.
      wice_soisno = 0._r8
      wliq_soisno = 10._r8
      dew_input = 0.3_r8
    end select
{allocation}
    if (abs(wliq_soisno(1)+wdsrf-(merge(10._r8,11._r8,icase==3)+dew_input)) > 1.e-12_r8) stop 1
    if (wliq_soisno(1)/denh2o+wice_soisno(1)/denice > porsl(1)*dz_soisno(1)+1.e-14_r8) stop 2
    if (icase <= 2 .and. {frost_name} <= 1._r8) stop 3
    if (icase == 1 .and. abs({dew_name}-dew_input) > 1.e-12_r8) stop 4
    if (icase == 2 .and. abs({dew_name}) > 1.e-12_r8) stop 5
    if (icase == 3 .and. ({frost_name} /= 0._r8 .or. {dew_name} /= 0._r8)) stop 6
  end do
end program
"""
        src = tmp_path / f"check_frost_{name}.F90"
        exe = tmp_path / f"check_frost_{name}"
        src.write_text(program)
        subprocess.run([compiler, "-cpp", "-DTRACER", "-ffree-line-length-none", str(src), "-o", str(exe)],
                       check=True, capture_output=True, text=True)
        subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_late_frost_exports_post_transport_soil_composition(tmp_path):
    """Same-step infiltration/layer transport precedes frost's old-water export."""
    source = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
    start = source.index("         IF (late_surface_water > trc_tiny) THEN\n            ! At the post-solver instant")
    marker = source.index("         ! 4b. WATER-phase", start)
    end = source.rfind("         ENDIF", start, marker) + len("         ENDIF")
    late_transfer = source[start:end]
    program = f"""
program check_late_frost
  implicit none
  integer, parameter :: r8 = kind(1.d0), ipatch=1, itrc=1
  real(r8), parameter :: trc_tiny=1.e-14_r8, trc_water_min_for_ratio=1.e-12_r8
  real(r8) :: trc_wliq_soisno(1,1,1), trc_wdsrf(1,1), trc_surface_residue(1,1)
  real(r8) :: trc_surface_solid(1,1)
  real(r8) :: a_trc_precip(1,1), a_trc_rsur(1,1), a_trc_rnof(1,1), trc_rnof_step(1,1)
  real(r8) :: water_shadow(1), layer_temp(1), frost_surface_water, dew_surface_water
  real(r8) :: late_surface_water, pending_surface_tracer, late_surface_ratio
  real(r8) :: late_water, late_runoff_water, trc_flux, wdsrf, rsur, deltim, expected_donor, before_total
  integer :: icase
  logical :: solute
  deltim = 1800._r8
  layer_temp = 270._r8
  do icase=1,2
    solute = icase==2
    ! Earlier in this same step, infiltration adds 2 mm and qlayer moves
    ! 3 mm downward. Top water is now 9 mm with a different composition
    ! from its pre-infiltration 10 mm / 2.0 tracer state.
    water_shadow = 9._r8
    trc_wliq_soisno = merge(0.09_r8,1.65_r8,solute)
    pending_surface_tracer = merge(0.5_r8,0.05_r8,solute)
    frost_surface_water = 1._r8
    dew_surface_water = 0.5_r8
    late_surface_water = frost_surface_water+dew_surface_water
    wdsrf=2._r8
    rsur=0.5_r8/deltim
    trc_wdsrf=0._r8
    trc_surface_residue=0._r8
    trc_surface_solid=0._r8
    a_trc_precip=0._r8
    a_trc_rsur=0._r8
    a_trc_rnof=0._r8
    trc_rnof_step=0._r8
    late_runoff_water=0.5_r8
    expected_donor = trc_wliq_soisno(1,1,1)/water_shadow(1)
    before_total = trc_wliq_soisno(1,1,1)+pending_surface_tracer
{late_transfer}
    if (abs(trc_wliq_soisno(1,1,1) - &
        (merge(0.09_r8,1.65_r8,solute)-expected_donor)) > 1.e-12_r8) stop 1
    if (abs(trc_wdsrf(1,1)+a_trc_rsur(1,1)+trc_wliq_soisno(1,1,1) - &
        (before_total+dew_surface_water*merge(0._r8,0.001_r8,solute))) > 1.e-12_r8) stop 2
    if (abs(a_trc_rsur(1,1)-rsur*deltim*late_surface_ratio) > 1.e-12_r8) stop 3
    if (abs(water_shadow(1)-8._r8) > 1.e-12_r8) stop 4
  end do
contains
  real(r8) function current_liq_ratio(j)
    integer, intent(in) :: j
    current_liq_ratio=trc_wliq_soisno(1,j,1)/water_shadow(j)
  end function
  real(r8) function deposition_ratio_for(temp, ice)
    real(r8), intent(in) :: temp
    logical, intent(in) :: ice
    deposition_ratio_for=merge(0._r8,0.001_r8,solute)
  end function
  logical function tracer_is_nonvolatile_solute(i)
    integer, intent(in) :: i
    tracer_is_nonvolatile_solute=solute
  end function
  subroutine CoLM_stop(message)
    character(*), intent(in) :: message
    print *, message
    stop 9
  end subroutine
  subroutine tracer_equilibrate_dissolved(i,water,dissolved,solid)
    integer, intent(in) :: i
    real(r8), intent(in) :: water
    real(r8), intent(inout) :: dissolved,solid
  end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_late_frost.F90"
    exe = tmp_path / "check_late_frost"
    src.write_text(program)
    subprocess.run([compiler, "-ffree-line-length-none", str(src), "-o", str(exe)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_water_2014_frost_pond_is_consumed_exactly_once_next_step(tmp_path):
    """The production carry + qinfl + frost snippets close across two steps."""
    source = (ROOT / "main/MOD_SoilSnowHydrology.F90").read_text()
    carry_start = source.index("      ! Carry the previous step's pond into this step's surface input once.")
    carry_end = source.index("#ifdef CROP", carry_start)
    carry = source[carry_start:carry_end]
    qinfl = "      qinfl = gwat - rsur - wdsrf/deltim"
    assert qinfl in source
    frost_start = source.index("      dew_capacity = max((porsl(1)*dz_soisno(1)-wice_soisno(1)/denice)*denh2o")
    frost_end = source.index("      err_solver =", frost_start)
    frost = source[frost_start:frost_end]
    program = f"""
program check_two_step_pond
  implicit none
  integer, parameter :: r8=kind(1.d0)
  real(r8), parameter :: denh2o=1000._r8, denice=917._r8
  integer :: patchtype, istep
  real(r8) :: porsl(1), dz_soisno(1), wice_soisno(1), wliq_soisno(1)
  real(r8) :: dew_input, dew_capacity, dew_retained, dew_excess, frost_excess, ice_before_frost
  real(r8) :: wdsrf, gwat, deltim, rsur, rnof, qinfl, pondmx
  real(r8) :: dew_overflow_trc, frost_displaced_trc, late_runoff_trc, rsur_before_late
  deltim=1800._r8; pondmx=5._r8; patchtype=0
  porsl=0.4_r8; dz_soisno=0.03_r8
  wliq_soisno=11._r8; wice_soisno=0._r8; wdsrf=0._r8
  do istep=1,2
    gwat=0._r8; rsur=0._r8; rnof=0._r8
{carry}
    if (istep==2) rsur=gwat ! frozen saturated soil sheds carried pond
{qinfl}
    if (abs(qinfl)>1.e-12_r8) stop 1
    rsur_before_late=rsur
    dew_input=0._r8
    ice_before_frost=wice_soisno(1)
    if (istep==1) wice_soisno(1)=2._r8 ! atmospheric frost
    dew_overflow_trc=0._r8
    frost_displaced_trc=0._r8
    late_runoff_trc=0._r8
{frost}
    if (istep==1) then
      if (frost_displaced_trc<=1._r8 .or. abs(wdsrf-frost_displaced_trc)>1.e-12_r8) stop 2
      if (abs(wliq_soisno(1)+wice_soisno(1)+wdsrf-13._r8)>1.e-12_r8) stop 3
    else
      if (frost_displaced_trc/=0._r8) stop 4
      if (wdsrf/=0._r8) stop 5
      if (abs(wliq_soisno(1)+wice_soisno(1)+rsur*deltim-13._r8)>1.e-12_r8) stop 6
    end if
  end do
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_two_step_pond.F90"
    exe = tmp_path / "check_two_step_pond"
    src.write_text(program)
    subprocess.run([compiler, "-cpp", "-DTRACER", "-ffree-line-length-none", str(src), "-o", str(exe)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_phs_root_return_books_net_transpiration_from_actual_donors(tmp_path):
    """Run the production return and vapor-booking blocks at U=R and R<U."""
    source = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
    donor_start = source.index("                     trc_flux = etroot_actual(j) * ratio_layer(j)")
    donor_end = source.index("                     root_gross_tracer = root_gross_tracer + trc_flux", donor_start) + len(
        "                     root_gross_tracer = root_gross_tracer + trc_flux"
    )
    donor_block = source[donor_start:donor_end]
    ratio_start = source.index("         return_ratio = xylem_ratio")
    ratio_end = source.index("               IF (transp_frac_active .and.", ratio_start)
    ratio_block = source[ratio_start:ratio_end]
    return_start = source.index("         IF (root_return_water > trc_tiny) THEN", ratio_end)
    return_end = source.index("         ! Baseflow leaves", return_start)
    return_block = source[return_start:return_end]
    book_start = source.index("                  IF (transp_frac_active) THEN", return_end)
    book_end = source.index("            ! ============================================================", book_start)
    book_block = source[book_start:book_end]
    program = f"""
program check_root_return
  implicit none
  integer, parameter :: r8=kind(1.d0), itrc=1, ipatch=1, nl_soil=2
  integer, parameter :: TRC_EVAP_KIND_TRANSP=1
  real(r8), parameter :: trc_tiny=1.e-14_r8
  real(r8) :: root_gross_water, root_gross_tracer, root_return_water, wa_bef
  real(r8) :: trc_flux, ratio_layer(nl_soil)
  real(r8) :: xylem_ratio, return_ratio, transp_ratio, transp_water_total
  real(r8) :: transp_source_tracer_total, transp_output_tracer, root_return_tracer_total
  real(r8) :: root_return_tracer, etroot_actual(nl_soil), etroot_aquifer, surface_root_return
  real(r8) :: water_shadow(nl_soil), trc_wliq_soisno(1,nl_soil,1), trc_wa(1,1), trc_wdsrf(1,1)
  real(r8) :: a_trc_evap(1,1), a_trc_transp(1,1), a_trc_transp_src(1,1)
  real(r8) :: a_water_transp(1,1), a_water_evap_gross(1,1)
  real(r8), allocatable :: trc_leaf_iso_storage(:,:)
  real(r8) :: u, r, m, expected_return, expected_net
  integer :: icase,j
  logical :: transp_frac_active
  allocate(trc_leaf_iso_storage(1,1))
  do icase=1,4
    u=0.1_r8
    r=merge(u,0.04_r8,icase==1)
    m=merge(0.0001_r8,0.0002_r8,icase==3) ! actual finite donor cap
    root_gross_water=u
    root_gross_tracer=0._r8
    root_return_water=r
    transp_water_total=u-r
    xylem_ratio=0.002_r8 ! theoretical source mixture; cap must override it
    transp_ratio=xylem_ratio
    transp_frac_active=icase==4
    trc_wliq_soisno=0._r8
    trc_wliq_soisno(1,2,1)=m
    ratio_layer=0.002_r8
    etroot_actual=[-r,u]
    j=2
{donor_block}
    if (abs(root_gross_tracer-m)>1.e-14_r8) stop 6
    transp_source_tracer_total=merge(root_gross_tracer,0._r8,transp_frac_active)
    trc_wa=0._r8; trc_wdsrf=0._r8; water_shadow=0._r8
    etroot_aquifer=0._r8; surface_root_return=0._r8; wa_bef=0._r8
    a_trc_evap=merge(0._r8,m,transp_frac_active)
    a_trc_transp=a_trc_evap; a_trc_transp_src=a_trc_evap
    a_water_transp=merge(0._r8,u,transp_frac_active)
    a_water_evap_gross=a_water_transp
    trc_leaf_iso_storage=0._r8
{ratio_block}
{return_block}
{book_block}
    expected_return=r*m/u
    expected_net=m-expected_return
    if (abs(trc_wliq_soisno(1,1,1)-expected_return)>1.e-14_r8) stop 1
    if (abs(transp_source_tracer_total-merge(expected_net,0._r8,transp_frac_active))>1.e-14_r8) stop 2
    if (abs(a_trc_evap(1,1)-expected_net)>1.e-14_r8) stop 3
    if (abs(a_water_transp(1,1)-(u-r))>1.e-14_r8) stop 4
    if (abs(trc_leaf_iso_storage(1,1))>1.e-14_r8) stop 5
  end do
contains
  logical function tracer_is_nonvolatile_solute(i)
    integer, intent(in) :: i
    tracer_is_nonvolatile_solute=.false.
  end function
  subroutine CoLM_stop(message)
    character(*), intent(in) :: message
    print *,message
    stop 9
  end subroutine
  subroutine check_isotope_aquifer(i,p,water,mass,stage)
    integer, intent(in) :: i,p
    real(r8), intent(in) :: water,mass
    character(*), intent(in) :: stage
  end subroutine
  subroutine tracer_book_evap_loss(i,p,mass,water,kind)
    integer, intent(in) :: i,p,kind
    real(r8), intent(in) :: mass,water
    a_trc_evap(i,p)=a_trc_evap(i,p)+mass
    a_trc_transp(i,p)=a_trc_transp(i,p)+mass
    a_water_transp(i,p)=a_water_transp(i,p)+water
    a_water_evap_gross(i,p)=a_water_evap_gross(i,p)+water
  end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_root_return.F90"
    exe = tmp_path / "check_root_return"
    src.write_text(program)
    subprocess.run([compiler, "-ffree-line-length-none", str(src), "-o", str(exe)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_flood_and_late_frost_keep_carrier_and_sequence(tmp_path):
    """Execute flood mixing with a later frost arrival and an early runoff."""
    source = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
    start = source.index("         late_ratio = ratio\n         IF (flood_water > 0._r8) THEN")
    end = source.index("         trc_wdsrf(itrc, ipatch) =", start)
    # The harness always supplies flood_tracer_input, unlike the optional
    # production dummy; keep the transport algebra verbatim.
    flood_mix = source[start:end].replace("present(flood_tracer_input)", ".true.")
    program = f"""
program check_flood_frost
  implicit none
  integer, parameter :: r8=kind(1.d0), itrc=1, ipatch=1
  real(r8), parameter :: trc_tiny=1.e-14_r8, trc_water_min_for_ratio=1.e-12_r8
  real(r8), parameter :: trc_delta_sanity_max=1000._r8
  type tracer_type
    real(r8) :: ref_ratio
  end type
  type(tracer_type) :: tracers(1)
  real(r8) :: ratio, late_ratio, late_tracer, late_water, trc_pool_total, flood_destination_water
  real(r8) :: flood_water, flood_ground_evap_tracer, flood_tracer_input(1)
  real(r8) :: trc_soil_upflow, wdsrf, late_runoff_water, late_surface_water
  real(r8) :: qinfl, deltim, trc_surface_residue(1,1), early_runoff_water
  integer :: icase
  deltim=1800._r8; tracers(1)%ref_ratio=0.002_r8
  flood_water=2._r8; flood_tracer_input=0.004_r8
  trc_soil_upflow=0._r8; trc_surface_residue=0._r8
  do icase=1,3
    select case(icase)
    case(1) ! Flood-only early carrier, 0.1 mm frost arrives after infiltration.
      ratio=0._r8; trc_pool_total=0._r8
      early_runoff_water=0._r8; late_runoff_water=0._r8
      wdsrf=1.1_r8; qinfl=1._r8/deltim; late_surface_water=0.1_r8
      flood_ground_evap_tracer=0._r8
    case(2) ! Old ordinary runoff leaves before flood/frost mixing.
      ratio=0.001_r8; trc_pool_total=0.0006_r8
      early_runoff_water=0.3_r8; late_runoff_water=0.2_r8
      wdsrf=1.8_r8; qinfl=0.7_r8/deltim; late_surface_water=0.1_r8
      flood_ground_evap_tracer=0._r8
    case(3) ! One mm flood pays pre-flood evaporation; only one mm carries on.
      ratio=0._r8; trc_pool_total=0._r8
      early_runoff_water=0._r8; late_runoff_water=0._r8
      wdsrf=0.1_r8; qinfl=1._r8/deltim; late_surface_water=0.1_r8
      flood_ground_evap_tracer=0.002_r8
    end select
    flood_destination_water=max(wdsrf,0._r8)+late_runoff_water-late_surface_water &
       + max(qinfl,0._r8)*deltim
    flood_destination_water=max(flood_destination_water,0._r8)
{flood_mix}
    if (icase==1 .and. abs(late_water-2._r8)>1.e-12_r8) stop 1
    if (icase==2 .and. abs(late_water-2.6_r8)>1.e-12_r8) stop 2
    if (icase==3 .and. abs(late_water-1._r8)>1.e-12_r8) stop 3
    if (abs(late_ratio-0.002_r8)>1.e-12_r8 .and. icase/=2) stop 4
    if (icase==2 .and. abs(late_ratio-0.0046_r8/2.6_r8)>1.e-12_r8) stop 5
  end do
contains
  logical function tracer_is_isotope(i)
    integer, intent(in) :: i
    tracer_is_isotope=.true.
  end function
  logical function tracer_is_nonvolatile_solute(i)
    integer, intent(in) :: i
    tracer_is_nonvolatile_solute=.false.
  end function
  subroutine CoLM_stop(message)
    character(*), intent(in) :: message
    print *,message
    stop 9
  end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_flood_frost.F90"
    exe = tmp_path / "check_flood_frost"
    src.write_text(program)
    subprocess.run([compiler, "-ffree-line-length-none", str(src), "-o", str(exe)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_late_pond_overflow_does_not_turn_exfiltration_into_evaporation(tmp_path):
    source = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
    start = source.index("            qgtop_est = qinfl +")
    end = source.index("\n\n", start)
    estimate = source[start:end]
    program = f"""
program check_exfil
  implicit none
  integer, parameter :: r8=kind(1.d0)
  real(r8), parameter :: trc_tiny=1.e-14_r8
  real(r8) :: qgtop_est, qinfl, wdsrf, wdsrf_bef, deltim
  real(r8) :: late_runoff_water, flood_water, late_surface_water
  deltim=1._r8; wdsrf_bef=0._r8
  qinfl=-1._r8; wdsrf=1._r8
  late_runoff_water=2._r8; flood_water=0._r8; late_surface_water=1._r8
{estimate}
  if (abs(qgtop_est-1._r8)>1.e-12_r8) stop 1
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_exfil.F90"
    exe = tmp_path / "check_exfil"
    src.write_text(program)
    subprocess.run([compiler, "-ffree-line-length-none", str(src), "-o", str(exe)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_new_frost_pure_ice_overcapacity_fails_closed(tmp_path):
    source = (ROOT / "main/MOD_SoilSnowHydrology.F90").read_text()
    compiler = require_runnable_fortran_compiler(tmp_path)
    for solver in ("WATER_2014", "WATER_VSF"):
        marker = f"CALL CoLM_stop('{solver}: new frost ice exceeds entire soil pore volume')"
        end = source.index(marker) + len(marker)
        start = source.rfind("      IF (wice_soisno(1) > ice_before_frost", 0, end)
        end = source.index("\n      ENDIF", end) + len("\n      ENDIF")
        guard = source[start:end]
        program = f"""
program check_frost_guard
  implicit none
  integer, parameter :: r8=kind(1.d0)
  real(r8), parameter :: denice=917._r8
  real(r8) :: wice_soisno(1), porsl(1), dz_soisno(1), ice_before_frost
  integer :: patchtype
  character(2) :: arg
  call get_command_argument(1,arg)
  patchtype=0
  porsl=0.4_r8; dz_soisno=0.03_r8; ice_before_frost=0._r8
  wice_soisno=merge(10._r8,12._r8,arg=='0')
  if (arg=='4') patchtype=1
  if (arg=='2' .or. arg=='4') then
    call check_guard(.true.)
  else if (arg=='3') then
    call check_guard(.false.)
  else
    call check_guard()
  endif
contains
  subroutine check_guard(defer_surface_ice_overflow)
    logical, optional, intent(in) :: defer_surface_ice_overflow
{guard}
  end subroutine
  subroutine CoLM_stop(message)
    character(*), intent(in) :: message
    print *,message
    stop 9
  end subroutine
end program
"""
        src = tmp_path / f"check_guard_{solver}.F90"
        exe = tmp_path / f"check_guard_{solver}"
        src.write_text(program)
        subprocess.run([compiler, str(src), "-o", str(exe)], check=True, capture_output=True, text=True)
        subprocess.run([str(exe), "0"], check=True, capture_output=True, text=True)
        subprocess.run([str(exe), "2"], check=True, capture_output=True, text=True)
        for case in ("1", "3", "4"):
            failed = subprocess.run([str(exe), case], capture_output=True, text=True)
            assert failed.returncode != 0 and "new frost ice exceeds" in failed.stdout


def test_vsf_pond_withdrawal_is_not_reported_as_infiltration():
    """The old pond-delta shortcut created tracer water with no soil-water gain."""
    hydro = (ROOT / "main/HYDRO/MOD_Hydro_SoilWater.F90").read_text()
    assert "pond_exchange = max(exchange_dp_before-ss_dp, 0._r8)" in hydro
    assert "qinfl = qgtop - (ss_dp - dp_m1 + pond_exchange)/dt" in hydro
    # Pure ET, pure subsurface runoff and their mixture all consume the same
    # pond, but none is an infiltration flux.
    for et, rsub in ((0.04, 0.0), (0.0, 0.04), (0.01, 0.03)):
        dp_before = 0.1
        dp_after = dp_before - et - rsub
        pond_exchange = max(dp_before - dp_after, 0.0)
        qinfl = 0.0 - (dp_after - dp_before + pond_exchange) / 1800.0
        assert abs(qinfl) < 1e-18


def test_vsf_exchange_reports_actual_surface_soil_and_aquifer_donors(tmp_path):
    """Execute the production exchange diagnostic around a controlled solver stub.

    The actual water solver chooses the donor pools; this test pins the
    pre/post-state attribution. Hydro's sp_zi/sp_dz geometry is already mm.
    """
    source = (ROOT / "main/HYDRO/MOD_Hydro_SoilWater.F90").read_text()
    start = source.index("      wexchange = rsubst * dt + deficit")
    end = source.index("      ! water table location", start)
    diagnostic = source[start:end]
    def mock_call(match):
        actual = match.group()
        amount = "deficit" if "nlev, deficit," in actual else ("rsubst*dt" if "nlev, rsubst*dt," in actual else "wexchange")
        return f"      call mock_exchange({amount}, ss_dp, ss_vliq, zwt, wa)"

    diagnostic, replacements = re.subn(
        r"      CALL soilwater_aquifer_exchange \(\s*&.*?\bwa, izwt\)",
        mock_call, diagnostic, flags=re.DOTALL,
    )
    assert replacements == 3
    program = f"""
program check_vsf_donors
  implicit none
  integer, parameter :: r8 = kind(1.d0), nlev = 3
  real(r8), parameter :: dt = 1800._r8
  real(r8) :: deficit, rsubst, wexchange, ss_dp, zwt, wa, dp_before
  real(r8) :: sp_zi(0:nlev), sp_dz(nlev), ss_vliq(nlev), porsl(nlev)
  real(r8) :: exchange_dp_before, pond_exchange, exchange_layer_before(nlev)
  real(r8) :: exchange_layer_after, exchange_wa_before, exchange_zwt_before
  real(r8) :: et_fraction, rsub_fraction, layer_debit
  real(r8) :: etroot_actual_out(nlev), etroot_aquifer_out, etroot_surface_out
  real(r8) :: rsub_layer_out(nlev), rsub_surface_out, rsub_aquifer_out
  real(r8) :: withdraw_surface, withdraw_layer, withdraw_aquifer
  logical :: DEF_USE_PLANTHYDRAULICS
  integer :: icase, ilev
  sp_zi = [0._r8, 100._r8, 200._r8, 300._r8]
  sp_dz = 100._r8
  porsl = 0.4_r8
  do icase = 1, 7
    ss_dp = 0.1_r8
    ss_vliq = 0.2_r8
    zwt = 1000._r8
    wa = -100._r8
    etroot_actual_out = 0._r8
    etroot_aquifer_out = 0._r8
    etroot_surface_out = 0._r8
    rsub_layer_out = 0._r8
    rsub_surface_out = 0._r8
    rsub_aquifer_out = 0._r8
    DEF_USE_PLANTHYDRAULICS = icase >= 5
    dp_before = ss_dp
    select case (icase)
    case (1) ! ET solely from the pond, not the aquifer.
      deficit = 0.03_r8
      rsubst = 0._r8
      withdraw_surface = 0.03_r8
      withdraw_layer = 0._r8
      withdraw_aquifer = 0._r8
    case (2) ! Baseflow solely from the shallow soil, in millimetres.
      deficit = 0._r8
      rsubst = 0.03_r8/dt
      withdraw_surface = 0._r8
      withdraw_layer = 0.03_r8
      withdraw_aquifer = 0._r8
    case (3) ! Shared aquifer withdrawal, with correct ET/rsub split.
      deficit = 0.01_r8
      rsubst = 0.02_r8/dt
      withdraw_surface = 0._r8
      withdraw_layer = 0._r8
      withdraw_aquifer = 0.03_r8
    case (4) ! Mixed pond, layer, aquifer donors.
      deficit = 0.01_r8
      rsubst = 0.02_r8/dt
      withdraw_surface = 0.01_r8
      withdraw_layer = 0.01_r8
      withdraw_aquifer = 0.01_r8
    case (5) ! Root return to pond must not become false soil exfiltration.
      deficit = -0.02_r8
      rsubst = 0._r8
      withdraw_surface = -0.02_r8
      withdraw_layer = 0._r8
      withdraw_aquifer = 0._r8
    case (6) ! Return first, then baseflow from the same mixed pond.
      deficit = -0.02_r8
      rsubst = 0.03_r8/dt
      withdraw_surface = 0.01_r8
      withdraw_layer = 0._r8
      withdraw_aquifer = 0._r8
    case (7) ! Gross opposing exchanges with exactly zero net pond change.
      deficit = -0.03_r8
      rsubst = 0.03_r8/dt
      withdraw_surface = 0._r8
      withdraw_layer = 0._r8
      withdraw_aquifer = 0._r8
    end select
{diagnostic}
    if (abs(etroot_surface_out+rsub_surface_out-withdraw_surface) > 1.e-11_r8) stop 1
    if (abs(sum(etroot_actual_out)+sum(rsub_layer_out)-withdraw_layer) > 1.e-11_r8) stop 2
    if (abs(etroot_aquifer_out+rsub_aquifer_out-withdraw_aquifer) > 1.e-11_r8) stop 3
    if (abs(etroot_surface_out+sum(etroot_actual_out)+etroot_aquifer_out-deficit) > 1.e-11_r8) stop 4
    if (abs(rsub_surface_out+sum(rsub_layer_out)+rsub_aquifer_out-rsubst*dt) > 1.e-11_r8) stop 5
    if (icase >= 5 .and. abs((ss_dp-dp_before+pond_exchange)/dt) > 1.e-13_r8) stop 7
  end do
contains
  subroutine mock_exchange(wex, dp, vliq, wt, aqu)
    real(r8), intent(in) :: wex
    real(r8), intent(inout) :: dp, vliq(nlev), wt, aqu
    if (icase >= 5) then
      dp = dp-wex
      return
    end if
    if (abs(wex-withdraw_surface-withdraw_layer-withdraw_aquifer) > 1.e-11_r8) stop 6
    dp = dp-withdraw_surface
    vliq(1) = vliq(1)-withdraw_layer/sp_dz(1)
    aqu = aqu-withdraw_aquifer
  end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = tmp_path / "check_vsf_donors.F90"
    exe = tmp_path / "check_vsf_donors"
    src.write_text(program)
    subprocess.run(
        [compiler, "-cpp", "-DTRACER", "-ffree-line-length-none", str(src), "-o", str(exe)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)
