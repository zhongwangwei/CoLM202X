"""A tiny but finite flood carrier keeps a physically bounded isotope mass."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_SoilWater.F90"


def test_late_flood_microcarrier_uses_production_partition(tmp_path):
    source = SOURCE.read_text()
    carrier_start = source.index("      flood_destination_water = max(wdsrf,0._r8)+late_runoff_water-late_surface_water &")
    carrier_end = source.index("      flood_destination_water = max(flood_destination_water,0._r8)", carrier_start)
    carrier_end = source.index("\n", carrier_end) + 1
    carrier = source[carrier_start:carrier_end].replace("flood_destination_water", "late_water")
    start = source.index("            IF (late_water > trc_water_min_for_ratio .or. &")
    end = source.index("\n            ENDIF\n            trc_pool_total = late_tracer", start) + len("\n            ENDIF")
    partition = source[start:end]
    program = f"""
program check_microcarrier
  implicit none
  integer, parameter :: r8=kind(1.d0)
  real(r8), parameter :: trc_tiny=1.e-30_r8, trc_water_min_for_ratio=1.e-12_r8
  real(r8), parameter :: trc_delta_sanity_max=2.e3_r8
  type tracer_info
     real(r8) :: ref_ratio=2.0052e-3_r8
  end type
  type(tracer_info) :: tracers(1)
  real(r8) :: late_water,late_tracer,late_ratio,trc_surface_residue(1,1)
  real(r8) :: wdsrf,late_runoff_water,late_surface_water,qinfl,deltim
  integer :: itrc,ipatch,case_id
  character(16) :: arg
  itrc=1; ipatch=1; trc_surface_residue=0._r8
  deltim=1800._r8; qinfl=0._r8
  late_runoff_water=0._r8; late_surface_water=0._r8
  call get_command_argument(1,arg)
  read(arg,*) case_id
  select case(case_id)
  case(1) ! Measured IsoGSM flood boundary: water below generic floor.
     wdsrf=2.8422454e-14_r8; late_tracer=5.62818889e-17_r8
  case(2) ! With no late water the carrier equals the old pond+infiltration.
     wdsrf=1._r8; qinfl=1._r8/deltim; late_tracer=0.004_r8
  case(3) ! No carrier cannot transport positive isotope mass.
     wdsrf=0._r8; late_tracer=5.e-17_r8
  case(4) ! An excessive micro-pool ratio cannot bypass the dry guard.
     wdsrf=2.e-14_r8; late_tracer=2.e-16_r8
  case(5) ! A solute still deposits to dry surface residue.
     wdsrf=2.e-14_r8; late_tracer=5.e-17_r8
  case(6) ! Late frost is excluded from this step's flood infiltration.
     wdsrf=1.1_r8; late_surface_water=0.1_r8
     qinfl=1._r8/deltim; late_tracer=0.004_r8
  end select
  late_ratio=-1._r8
{carrier}
{partition}
  select case(case_id)
  case(1)
     if(abs(late_ratio-5.62818889e-17_r8/2.8422454e-14_r8)>1.e-14_r8) stop 1
     if(trc_surface_residue(1,1)/=0._r8) stop 2
  case(2)
     if(abs(late_ratio-0.002_r8)>1.e-14_r8) stop 3
  case(6)
     if(abs(late_water-2._r8)>1.e-14_r8) stop 6
     if(abs(late_ratio-0.002_r8)>1.e-14_r8) stop 7
  case(5)
     if(late_ratio/=0._r8) stop 4
     if(abs(trc_surface_residue(1,1)-5.e-17_r8)>1.e-29_r8) stop 5
  end select
contains
  logical function tracer_is_isotope(idx)
     integer, intent(in) :: idx
     tracer_is_isotope=case_id/=5
  end function
  logical function tracer_is_nonvolatile_solute(idx)
     integer, intent(in) :: idx
     tracer_is_nonvolatile_solute=case_id==5
  end function
  subroutine CoLM_stop(message)
     character(*), intent(in) :: message
     error stop 9
  end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    path = tmp_path / "microcarrier.f90"
    exe = tmp_path / "microcarrier"
    path.write_text(program)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", str(path), "-o", str(exe)],
        capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    for case_id in (1, 2, 5, 6):
        ran = subprocess.run([str(exe), str(case_id)], capture_output=True, text=True, timeout=30)
        assert ran.returncode == 0, (case_id, ran.stdout, ran.stderr)
    for case_id in (3, 4):
        ran = subprocess.run([str(exe), str(case_id)], capture_output=True, text=True, timeout=30)
        assert ran.returncode != 0, (case_id, ran.stdout, ran.stderr)


def test_all_flood_tracer_evaporates_when_solver_has_no_water_destination(tmp_path):
    source = SOURCE.read_text()
    start = source.index("         flood_ground_evap_water = 0._r8")
    end = source.index("         flood_ground_evap_tracer = 0._r8", start)
    allocation = source[start:end]
    program = f"""
program check_flood_evap_allocation
  implicit none
  integer, parameter :: r8=kind(1.d0)
  real(r8) :: flood_ground_evap_water, flood_water, surface_base_balance
  real(r8) :: gwat_evap, flood_destination_water
  integer :: icase
  do icase=1,3
     select case(icase)
     case(1) ! Measured failure: cancellation left isotope without a water destination.
        flood_water=2.828055473327107e-15_r8
        surface_base_balance=-3.740607399434487e-3_r8
        gwat_evap=2.828032946711190e-15_r8
        flood_destination_water=0._r8
     case(2) ! Finite destination: only the diagnosed one mm evaporates.
        flood_water=2._r8
        surface_base_balance=-1._r8
        gwat_evap=1.5_r8
        flood_destination_water=1._r8
     case(3) ! No pre-flood deficit: ordinary evaporation remains untouched.
        flood_water=2._r8
        surface_base_balance=0.5_r8
        gwat_evap=0.25_r8
        flood_destination_water=2.5_r8
     end select
{allocation}
     select case(icase)
     case(1)
        if(flood_ground_evap_water/=flood_water) stop 1
        if(gwat_evap/=0._r8) stop 2
     case(2)
        if(abs(flood_ground_evap_water-1._r8)>1.e-14_r8) stop 3
        if(abs(gwat_evap-0.5_r8)>1.e-14_r8) stop 4
     case(3)
        if(flood_ground_evap_water/=0._r8) stop 5
        if(abs(gwat_evap-0.25_r8)>1.e-14_r8) stop 6
     end select
  end do
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    path = tmp_path / "flood_evap_allocation.f90"
    exe = tmp_path / "flood_evap_allocation"
    path.write_text(program)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", str(path), "-o", str(exe)],
        capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(exe)], capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
