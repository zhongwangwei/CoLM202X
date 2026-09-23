"""A saturated PHS column must not turn negative ET residue into pond water."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]


def test_saturated_phs_net_return_needs_a_donor(tmp_path):
    source = (ROOT / "main/HYDRO/MOD_Hydro_SoilWater.F90").read_text()
    start = source.index("      ! With a saturated column there is no resolved root donor.")
    end = source.index("      ! Exchange water with aquifer", start)
    guard = source[start:end]
    driver = f"""
program check_phs_return
  implicit none
  integer, parameter :: r8=kind(1.d0)
  integer :: izwt
  real(r8) :: etr, deficit, dt, etroot(10)
  logical :: DEF_USE_PLANTHYDRAULICS
  character(len=16) :: mode
  call get_command_argument(1, mode)
  DEF_USE_PLANTHYDRAULICS=.true.
  dt=1800._r8
  etroot=[-5.0898000475722737e-8_r8,-1.3094824883084741e-7_r8, &
          -1.6566546632793867e-7_r8,-7.4660785671791052e-8_r8, &
           6.9832438520470777e-8_r8, 1.2786694409613638e-7_r8, &
           1.0394593669609849e-7_r8, 6.2205263081654699e-8_r8, &
           3.1813829860943006e-8_r8, 2.6508088342650939e-8_r8]
  izwt=1
  etr=-7.0830538362942706e-16_r8
  deficit=-1.2750219445783334e-12_r8
  select case (trim(mode))
  case ('resolved')
    izwt=2
  case ('positive_et')
    etr=1.e-5_r8
  case ('invalid')
    deficit=-1.e-5_r8
  end select
{guard}
  if (trim(mode)=='tiny' .and. deficit/=0._r8) error stop 1
  if (trim(mode)=='resolved' .and. deficit/=-1.2750219445783334e-12_r8) error stop 2
  if (trim(mode)=='positive_et' .and. deficit/=-1.2750219445783334e-12_r8) error stop 3
contains
  subroutine CoLM_stop(message)
    character(*),intent(in) :: message
    print *, trim(message)
    error stop 4
  end subroutine
end program
"""
    fc = require_runnable_fortran_compiler(tmp_path)
    src, exe = tmp_path / "phs_return.f90", tmp_path / "phs_return"
    src.write_text(driver)
    built = subprocess.run([fc, str(src), "-o", str(exe)], capture_output=True, text=True)
    assert built.returncode == 0, built.stdout + built.stderr
    for mode in ("tiny", "resolved", "positive_et"):
        ran = subprocess.run([str(exe), mode], capture_output=True, text=True)
        assert ran.returncode == 0, (mode, ran.stdout, ran.stderr)
    invalid = subprocess.run([str(exe), "invalid"], capture_output=True, text=True)
    assert invalid.returncode != 0
    assert "negative plant hydraulic transpiration without resolved root donor" in invalid.stdout


def test_roundoff_return_does_not_erase_resolved_reverse_roots(tmp_path):
    source = (ROOT / "main/HYDRO/MOD_Hydro_SoilWater.F90").read_text()
    start = source.index("      ! The aquifer exchange can subtract nearly equal layer storages")
    end = source.index("   END SUBROUTINE soil_water_vertical_movement", start)
    guard = source[start:end].replace("#ifdef TRACER", "").replace("#endif", "")
    driver = f"""
program check_return_decomposition
  implicit none
  integer, parameter :: r8=kind(1.d0)
  real(r8) :: dt, etroot(2), etroot_actual_out(2)
  real(r8) :: etroot_aquifer_out, etroot_surface_out
  logical :: DEF_USE_PLANTHYDRAULICS
  character(len=16) :: mode
  call get_command_argument(1, mode)
  dt=1800._r8
  etroot=[3.e-6_r8,-3.e-6_r8]
  DEF_USE_PLANTHYDRAULICS=.true.
  etroot_actual_out=[-5.7376325912628090e-12_r8,0._r8]
  etroot_aquifer_out=0._r8
  etroot_surface_out=0._r8
  select case (trim(mode))
  case ('reverse')
    etroot_actual_out=[1.e-3_r8,-1.e-3_r8]
  case ('orphan')
    etroot_actual_out=[-1.e-5_r8,0._r8]
  end select
{guard}
  select case (trim(mode))
  case ('residue')
    if (any(etroot_actual_out/=0._r8)) error stop 1
  case ('reverse')
    if (any(etroot_actual_out/=[1.e-3_r8,-1.e-3_r8])) error stop 2
  case ('orphan')
    if (etroot_actual_out(1)/=-1.e-5_r8) error stop 3
  end select
end program
"""
    fc = require_runnable_fortran_compiler(tmp_path)
    src, exe = tmp_path / "return_decomposition.f90", tmp_path / "return_decomposition"
    src.write_text(driver)
    built = subprocess.run([fc, str(src), "-o", str(exe)], capture_output=True, text=True)
    assert built.returncode == 0, built.stdout + built.stderr
    for mode in ("residue", "reverse", "orphan"):
        ran = subprocess.run([str(exe), mode], capture_output=True, text=True)
        assert ran.returncode == 0, (mode, ran.stdout, ran.stderr)
