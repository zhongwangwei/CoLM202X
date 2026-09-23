"""The production aquifer exchange must not erase isotope debt at wa=0."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
DEFS = (ROOT / "main/TRACER/MOD_Tracer_Defs.F90").read_text()
SOIL = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
REST = (ROOT / "main/TRACER/MOD_Tracer_Rest.F90").read_text()
SPECIAL = (ROOT / "main/TRACER/MOD_Tracer_SpecialPatches.F90").read_text()


def _routine(source: str, signature: str, end: str) -> str:
    return signature + source.split(signature, 1)[1].split(end, 1)[0] + end


def test_empty_tracer_special_patches_return_before_reference_index(tmp_path):
    """Both production entry guards must be bounds-safe for zero tracers."""
    snippets = []
    for name in ("tracer_glacier_patch", "tracer_waterbody_patch"):
        routine = _routine(SPECIAL, f"   SUBROUTINE {name}", f"   END SUBROUTINE {name}")
        start = routine.index("      IF (ntracers <= 0) RETURN")
        end = routine.index("      DO itrc = 1, ntracers", start)
        snippets.append(routine[start:end])
    program = f"""
program check_empty_special
 implicit none
 integer, parameter :: r8=kind(1.d0)
 integer :: ntracers, counter
 real(8), allocatable :: trc_aquifer_ref_water(:)
 ntracers=0; counter=0
 allocate(trc_aquifer_ref_water(0))
 call glacier(1)
 call waterbody(1)
 if(counter/=0) error stop 1
 deallocate(trc_aquifer_ref_water)
 ntracers=1
 allocate(trc_aquifer_ref_water(1));trc_aquifer_ref_water=0.d0
 call glacier(1)
 call waterbody(1)
 if(counter/=2) error stop 2
contains
 subroutine glacier(ipatch)
  integer,intent(in)::ipatch
{snippets[0]}
  counter=counter+1
 end subroutine
 subroutine waterbody(ipatch)
  integer,intent(in)::ipatch
{snippets[1]}
  counter=counter+1
 end subroutine
 subroutine CoLM_stop(message)
  character(*),intent(in)::message
  error stop 9
 end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    source = tmp_path / "empty_special.f90"
    exe = tmp_path / "empty_special"
    source.write_text(program)
    built = subprocess.run([compiler, "-ffree-line-length-0", "-fcheck=all", str(source), "-o", str(exe)],
                           capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(exe)], capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr


def test_aquifer_repayment_executes_production_exchange_and_guard(tmp_path):
    """Matched repayment passes; orphan/negative-ratio states stop before export."""
    isotope_guard = _routine(
        DEFS,
        "   pure logical FUNCTION tracer_aquifer_isotope_state_valid",
        "   END FUNCTION tracer_aquifer_isotope_state_valid",
    )
    actual_water = _routine(DEFS, "   pure real(r8) FUNCTION tracer_aquifer_actual_water", "   END FUNCTION tracer_aquifer_actual_water")
    actual_mass = _routine(DEFS, "   pure real(r8) FUNCTION tracer_aquifer_actual_mass", "   END FUNCTION tracer_aquifer_actual_mass")
    ratio = _routine(DEFS, "   pure real(r8) FUNCTION tracer_aquifer_isotope_ratio", "   END FUNCTION tracer_aquifer_isotope_ratio")
    stop_guard = _routine(
        SOIL, "   SUBROUTINE check_isotope_aquifer", "   END SUBROUTINE check_isotope_aquifer"
    )
    qcharge = "         aquifer_water_pre_qcharge = wa_bef - etroot_aquifer" + SOIL.split(
        "         aquifer_water_pre_qcharge = wa_bef - etroot_aquifer", 1
    )[1].split("         ! A previously dry aquifer", 1)[0]
    assert "'after aquifer root exchange'" in SOIL
    assert "'after aquifer baseflow'" in SOIL
    assert "CALL check_isotope_aquifer(itrc, ipatch, wa, trc_wa(itrc, ipatch), 'after qcharge')" in SOIL
    assert SOIL.index("'after qcharge'") < SOIL.index("'soil end'") < SOIL.index(
        "trc_wa(itrc, ipatch) = 0._r8", SOIL.index("'soil end'")
    )

    program = f"""
module production_aquifer_probe
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer, parameter :: r8=kind(1.d0), nl_soil=1
 real(r8), parameter :: trc_tiny=1.e-30_r8, trc_water_min_for_ratio=1.e-12_r8
 type tracer_info
   real(r8) :: ref_ratio=0.002_r8
 end type
 type(tracer_info) :: tracers(2)
 real(r8) :: soil_ratio
 real(r8) :: trc_aquifer_ref_water(1),trc_aquifer_ref_mass(2,1)
 real(r8) :: trc_wa(2,1),trc_wliq_soisno(2,1,1),trc_subsurface_residue(2,1)
 real(r8) :: trc_subsurface_solid(2,1),a_trc_qcharge(2,1),water_shadow(1)
contains
{actual_water}
{actual_mass}
{ratio}
{isotope_guard}
{stop_guard}
 logical function tracer_is_isotope(idx)
   integer, intent(in) :: idx
   tracer_is_isotope=idx==1
 end function
 logical function tracer_is_nonvolatile_solute(idx)
   integer, intent(in) :: idx
   tracer_is_nonvolatile_solute=idx==2
 end function
 real(r8) function layer_transport_ratio(idx)
   integer, intent(in) :: idx
   layer_transport_ratio=soil_ratio
 end function
 subroutine CoLM_stop(message)
   character(*), intent(in) :: message
   error stop 9
 end subroutine
end module
program check_aquifer
 use production_aquifer_probe
 implicit none
 integer :: itrc,ipatch,j,case_id
 character(16) :: arg
 real(r8) :: wa,wa_bef,qcharge_eff,ratio_src,trc_flux,deltim,source_fallback_ratio
 real(r8) :: etroot_aquifer,aquifer_water_pre_qcharge,rsub_source_aquifer
 real(r8) :: aquifer_ref_water,aquifer_ref_mass,aquifer_actual_mass
 logical :: resolved_rsub
 call get_command_argument(1,arg)
 read(arg,*) case_id
 itrc=1; ipatch=1; deltim=1._r8
 trc_wliq_soisno=1._r8; trc_subsurface_residue=0._r8
 trc_subsurface_solid=0._r8; a_trc_qcharge=0._r8
 water_shadow=1._r8; source_fallback_ratio=0.002_r8
 trc_aquifer_ref_water=0._r8;trc_aquifer_ref_mass=0._r8
 aquifer_ref_water=0._r8;aquifer_ref_mass=0._r8
 etroot_aquifer=0._r8;rsub_source_aquifer=0._r8;resolved_rsub=.false.
 select case(case_id)
 case(1) ! Exact repayment at the same isotope ratio.
   wa_bef=-1._r8; trc_wa=-.002_r8; qcharge_eff=1._r8; soil_ratio=.002_r8
 case(2) ! Different recharge ratio: no carrier for the leftover isotope.
   wa_bef=-1._r8; trc_wa=-.002_r8; qcharge_eff=1._r8; soil_ratio=.003_r8
 case(3) ! Opposite signs appear before the water debt reaches zero.
   wa_bef=-1._r8; trc_wa=-.002_r8; qcharge_eff=.9_r8; soil_ratio=.004_r8
 case(4) ! Ordinary positive aquifer remains valid.
   wa_bef=1._r8; trc_wa=.002_r8; qcharge_eff=.1_r8; soil_ratio=.003_r8
 case(5) ! A valid negative signed debt survives a restart unchanged.
   wa_bef=-1._r8; trc_wa=-.002_r8; qcharge_eff=0._r8; soil_ratio=.002_r8
 case(6) ! Finite but overflowing ratio must fail without an FP trap.
   wa_bef=2.e-12_r8; trc_wa=huge(1._r8); qcharge_eff=0._r8; soil_ratio=.002_r8
 case(7) ! Isotope-free water (delta=-1000 permil) is an allowed limit.
   wa_bef=1._r8; trc_wa=0._r8; qcharge_eff=0._r8; soil_ratio=.002_r8
 case(8) ! Root return changed sign before qcharge export: use stage donor.
   wa_bef=-.1_r8;etroot_aquifer=-.2_r8;trc_wa=.0002_r8
   qcharge_eff=-.05_r8;soil_ratio=.002_r8
 case(9) ! Root return exactly repays debt, then qcharge builds new debt.
   wa_bef=-.1_r8;etroot_aquifer=-.1_r8;trc_wa=0._r8
   qcharge_eff=-.05_r8;soil_ratio=.002_r8
 case(10) ! Aquifer baseflow also changes the qcharge donor.
   wa_bef=1._r8;etroot_aquifer=.2_r8;rsub_source_aquifer=.3_r8
   resolved_rsub=.true.;trc_wa=.001_r8;qcharge_eff=-.1_r8;soil_ratio=.002_r8
 case(11) ! Large matching debt repayment must cancel at round-off.
   wa_bef=-12345.6789_r8;trc_wa=wa_bef*.002_r8
   qcharge_eff=-wa_bef;soil_ratio=.002_r8;trc_wliq_soisno=1000._r8
 case(12) ! Heterogeneous recharge repays wa=-1 without orphaning M.
   wa_bef=-1._r8;trc_wa(1,1)=-.002_r8;qcharge_eff=1._r8;soil_ratio=.003_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(13) ! Heterogeneous recharge crosses wa to positive.
   wa_bef=-1._r8;trc_wa(1,1)=-.002_r8;qcharge_eff=2._r8;soil_ratio=.003_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(14) ! At wa=0 the finite reference remains a real donor.
   wa_bef=0._r8;trc_wa(1,1)=0._r8;qcharge_eff=-1._r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(15) ! A returned root litre joins the reference before qcharge.
   wa_bef=-1._r8;etroot_aquifer=-1._r8;trc_wa(1,1)=0._r8
   qcharge_eff=-1._r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(16) ! Baseflow's donor was already debited before qcharge.
   wa_bef=0._r8;rsub_source_aquifer=2._r8;resolved_rsub=.true.
   trc_wa(1,1)=-.004_r8;qcharge_eff=-1._r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(17) ! Extraction beyond the real reference carrier is fatal.
   wa_bef=-9.9_r8;trc_wa(1,1)=-.0198_r8;qcharge_eff=-.2_r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(18) ! Same patch's finite solute must NOT see isotope-only Vref.
   itrc=2;wa_bef=1._r8;trc_wa(2,1)=.5_r8
   qcharge_eff=-.2_r8;soil_ratio=.5_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(19) ! A genuinely empty aquifer cannot supply another qcharge export.
   wa_bef=-10._r8;trc_wa(1,1)=-.02_r8;qcharge_eff=-.1_r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(20) ! Nor can it supply an additional root withdrawal.
   wa_bef=-10._r8;trc_wa(1,1)=-.02_r8;etroot_aquifer=.1_r8
   qcharge_eff=0._r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(21) ! Nor can it supply an additional resolved baseflow withdrawal.
   wa_bef=-10._r8;trc_wa(1,1)=-.02_r8;rsub_source_aquifer=.1_r8
   resolved_rsub=.true.;qcharge_eff=0._r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 case(22) ! Exact exhaustion is a valid terminal state, not an orphan.
   wa_bef=0._r8;trc_wa(1,1)=0._r8;qcharge_eff=-10._r8;soil_ratio=.002_r8
   trc_aquifer_ref_water=10._r8;trc_aquifer_ref_mass(1,1)=.02_r8
 end select
 if(tracer_is_isotope(itrc)) then
   aquifer_ref_water=trc_aquifer_ref_water(ipatch)
   aquifer_ref_mass=trc_aquifer_ref_mass(itrc,ipatch)
 endif
 wa=wa_bef-etroot_aquifer-rsub_source_aquifer+qcharge_eff*deltim
 if(case_id<8 .or. case_id>=12 .and. case_id/=15 .and. case_id/=16) &
   call check_isotope_aquifer(itrc,ipatch,wa_bef,trc_wa(itrc,1),'loaded state')
{qcharge}
 call check_isotope_aquifer(itrc,ipatch,wa,trc_wa(itrc,1),'after qcharge')
 print '(2(ES24.15,1X))',wa,trc_wa(itrc,1)
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    source = tmp_path / "aquifer_debt.F90"
    exe = tmp_path / "aquifer_debt"
    source.write_text(program)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", str(source), "-o", str(exe)],
        capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    for case_id, expected in ((1, (0.0, 0.0)), (4, (1.1, 0.0023)),
                              (5, (-1.0, -0.002)), (7, (1.0, 0.0)),
                              (8, (0.05, 0.0001)), (9, (-0.05, -0.0001)),
                              (10, (0.4, 0.0008)), (11, (0.0, 0.0)),
                              (12, (0.0, 0.001)), (13, (1.0, 0.004)),
                              (14, (-1.0, -0.002)), (15, (-1.0, -0.002)),
                              (16, (-3.0, -0.006)), (18, (0.8, 0.4)),
                              (22, (-10.0, -0.02))):
        ran = subprocess.run([str(exe), str(case_id)], capture_output=True, text=True, timeout=30)
        assert ran.returncode == 0, (case_id, ran.stdout, ran.stderr)
        actual = tuple(float(value) for value in ran.stdout.split()[-2:])
        assert all(abs(a - e) < 1e-13 for a, e in zip(actual, expected))
    for case_id in (2, 3, 6, 17, 19, 20, 21):
        ran = subprocess.run([str(exe), str(case_id)], capture_output=True, text=True, timeout=30)
        assert ran.returncode != 0, (case_id, ran.stdout, ran.stderr)
        assert "Unresolved isotope aquifer debt" in ran.stdout


def test_wetland_and_restart_reuse_signed_aquifer_contract():
    wetland = _routine(SOIL, "   SUBROUTINE tracer_wetland", "   END SUBROUTINE tracer_wetland")
    assert "'wetland start'" in wetland
    assert "'wetland mixed pool'" in wetland
    assert "'wetland after loss'" in wetland
    assert "'wetland end'" in wetland
    reader = _routine(REST, "   SUBROUTINE read_land_tracer_restart", "   END SUBROUTINE read_land_tracer_restart")
    writer = _routine(REST, "   SUBROUTINE write_land_tracer_restart", "   END SUBROUTINE write_land_tracer_restart")
    assert "CALL validate_land_tracer_restart_state(wa)" in reader
    assert "CALL validate_land_tracer_restart_state(wa)" in writer
    assert reader.index("'trc_wa', trc_wa") < reader.index("CALL validate_land_tracer_restart_state(wa)")
    assert "trc_wa(itrc, ip) = wa(ip) * R_init" not in reader  # valid restart is not re-seeded


def test_wetland_two_step_debt_repayment_uses_production_mixing(tmp_path):
    actual_water = _routine(DEFS, "   pure real(r8) FUNCTION tracer_aquifer_actual_water", "   END FUNCTION tracer_aquifer_actual_water")
    actual_mass = _routine(DEFS, "   pure real(r8) FUNCTION tracer_aquifer_actual_mass", "   END FUNCTION tracer_aquifer_actual_mass")
    guard = _routine(
        DEFS,
        "   pure logical FUNCTION tracer_aquifer_isotope_state_valid",
        "   END FUNCTION tracer_aquifer_isotope_state_valid",
    )
    pool_start = SOIL.index("         pool_water = wdsrf_bef + wa_bef + wetwat_bef")
    pool_end = SOIL.index("         ! Wetland hydrology mixes the aquifer", pool_start)
    mixed_pool = SOIL[pool_start:pool_end]
    assert "pool_tracer = trc_wdsrf(itrc, ipatch) + trc_wa(itrc, ipatch)" in mixed_pool
    program = f"""
module wetland_probe
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer, parameter :: r8=kind(1.d0)
 real(r8), parameter :: trc_water_min_for_ratio=1.e-12_r8
 real(r8) :: trc_wdsrf(1,1),trc_wa(1,1),trc_wetwat(1,1)
contains
{actual_water}
{actual_mass}
{guard}
end module
program check_wetland
 use wetland_probe
 implicit none
 integer :: case_id,step,itrc,ipatch
 character(16) :: arg
 real(r8) :: wa_bef,wdsrf_bef,wetwat_bef,wresi_sum,q_rain_in,q_sm_in,q_dew_in,q_frost_in
 real(r8) :: trc_wresi_sum,pool_water,pool_tracer,rain_ratio
 real(r8) :: aquifer_ref_water,aquifer_ref_mass
 call get_command_argument(1,arg)
 read(arg,*) case_id
 itrc=1;ipatch=1;wa_bef=-1._r8;trc_wa=-.002_r8
 trc_wdsrf=0._r8;trc_wetwat=0._r8
 wdsrf_bef=0._r8;wetwat_bef=0._r8
 wresi_sum=0._r8;q_sm_in=0._r8;q_dew_in=0._r8;q_frost_in=0._r8
 trc_wresi_sum=0._r8
 aquifer_ref_water=0._r8;aquifer_ref_mass=0._r8
 if(case_id>=4) then
   aquifer_ref_water=10._r8;aquifer_ref_mass=.02_r8
 end if
 if(case_id==5) trc_wa(1,1)=(wa_bef+aquifer_ref_water)*.003_r8-aquifer_ref_mass
 do step=1,2
   q_rain_in=.5_r8
   rain_ratio=.002_r8
   if(step==2 .and. case_id==2) rain_ratio=.003_r8
   if(step==2 .and. case_id==3) then
      q_rain_in=.4_r8;rain_ratio=.004_r8
   endif
   if(step==2 .and. case_id>=4) rain_ratio=.003_r8
   if(case_id==5) rain_ratio=.003_r8
{mixed_pool}
   pool_tracer=pool_tracer+q_rain_in*rain_ratio
   if(.not.tracer_aquifer_isotope_state_valid(pool_water,pool_tracer,.002_r8)) then
      print '(A,I0,2(1X,ES18.9))','wetland orphan step ',step,pool_water,pool_tracer
      error stop 9
   endif
   wa_bef=pool_water-aquifer_ref_water
   trc_wa(1,1)=pool_tracer-aquifer_ref_mass
 enddo
 print '(2(ES24.15,1X))',wa_bef,trc_wa(1,1)
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    source = tmp_path / "wetland_debt.F90"
    exe = tmp_path / "wetland_debt"
    source.write_text(program)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", str(source), "-o", str(exe)],
        capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    matched = subprocess.run([str(exe), "1"], capture_output=True, text=True, timeout=30)
    assert matched.returncode == 0, matched.stdout + matched.stderr
    assert all(abs(float(value)) < 1e-13 for value in matched.stdout.split()[-2:])
    buffered = subprocess.run([str(exe), "4"], capture_output=True, text=True, timeout=30)
    assert buffered.returncode == 0, buffered.stdout + buffered.stderr
    assert all(abs(float(a) - e) < 1e-13 for a, e in zip(buffered.stdout.split()[-2:], (0.0, 0.0005)))
    nonzero_init_delta = subprocess.run([str(exe), "5"], capture_output=True, text=True, timeout=30)
    assert nonzero_init_delta.returncode == 0, nonzero_init_delta.stdout + nonzero_init_delta.stderr
    assert all(abs(float(a) - e) < 1e-13 for a, e in zip(nonzero_init_delta.stdout.split()[-2:], (0.0, 0.01)))
    assert "(wa(ip) + trc_aquifer_ref_water(ip)) * R_init" in REST
    for case_id in (2, 3):
        mismatched = subprocess.run([str(exe), str(case_id)], capture_output=True, text=True, timeout=30)
        assert mismatched.returncode != 0
        assert "wetland orphan step 2" in mismatched.stdout
