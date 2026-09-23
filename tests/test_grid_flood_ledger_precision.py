"""The flood donor ledger uses transfers, not subtraction of large stores."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


FLOW = Path(__file__).resolve().parents[1] / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90"


def test_donor_ledger_stays_precise_with_large_pools_gain_and_clipping(tmp_path):
    source = FLOW.read_text()
    start = source.index("            DO j = 1, numucat\n               tracer_after = trc_mass(itrc,j)")
    end = source.index("\n            ENDDO", start) + len("\n            ENDDO")
    transfer = source[start:end]
    program = f"""
program check_flood_ledger
   implicit none
   integer, parameter :: r8=kind(1.d0), numucat=1
   integer :: case_id, itrc, j
   real(r8) :: trc_mass(1,1), trc_levsto(1,1), coefficient_uc(1), gain_uc(1)
   real(r8) :: flood_visible_tracer_uc(1,1), flood_protected_tracer_uc(1,1)
   real(r8) :: flood_visible_uc(1), flood_protected_uc(1)
   real(r8) :: tracer_exchange_ledger(3,1), tracer_after, tracer_cell_debit
   real(r8) :: expected, old_subtraction, before_visible, before_protected
   real(r8) :: expected_visible, expected_protected
   character(16) :: mode
   itrc=1
   call get_command_argument(1,mode)
   do case_id=1,5
      trc_mass=0._r8; trc_levsto=0._r8
      coefficient_uc=0._r8; gain_uc=0._r8
      flood_visible_tracer_uc=0._r8; flood_protected_tracer_uc=0._r8
      flood_visible_uc=0._r8; flood_protected_uc=0._r8
      tracer_exchange_ledger=0._r8
      select case(case_id)
      case(1) ! Tiny debit from two very large pools; old subtraction rounds it.
         trc_mass(1,1)=1.e10_r8; trc_levsto(1,1)=2.e10_r8
         flood_visible_tracer_uc(1,1)=1.e10_r8
         flood_protected_tracer_uc(1,1)=2.e10_r8
         coefficient_uc(1)=1.e-16_r8
         expected=3.e-6_r8
      case(2) ! Zero transfer must remain exactly zero at large storage.
         trc_mass(1,1)=1.e16_r8; trc_levsto(1,1)=1._r8
         expected=0._r8
      case(3) ! Distinct visible/protected debits and condensation gains.
         trc_mass(1,1)=1000._r8; trc_levsto(1,1)=500._r8
         flood_visible_tracer_uc(1,1)=400._r8
         flood_protected_tracer_uc(1,1)=100._r8
         flood_visible_uc(1)=100._r8; flood_protected_uc(1)=50._r8
         coefficient_uc(1)=0.25_r8; gain_uc(1)=0.01_r8
         expected=123.5_r8
      case(4) ! Vapour uptake alone is a negative river debit.
         trc_mass(1,1)=1._r8; trc_levsto(1,1)=2._r8
         flood_visible_uc(1)=0.5_r8; flood_protected_uc(1)=0.25_r8
         gain_uc(1)=0.1_r8
         expected=-0.075_r8
      case(5) ! Clamp only visible; protected storage remains untouched.
         trc_levsto(1,1)=1.e-3_r8
         flood_visible_tracer_uc(1,1)=5.e-11_r8
         coefficient_uc(1)=1._r8
         expected=0._r8
      end select
      before_visible=trc_mass(1,1)
      before_protected=trc_levsto(1,1)
      expected_visible=max(0._r8,before_visible-coefficient_uc(1)*flood_visible_tracer_uc(1,1) &
         +gain_uc(1)*flood_visible_uc(1))
      expected_protected=max(0._r8,before_protected-coefficient_uc(1)*flood_protected_tracer_uc(1,1) &
         +gain_uc(1)*flood_protected_uc(1))
{transfer}
      if(abs(tracer_exchange_ledger(3,1)-expected)>1.e-12_r8) stop 1
      if(trc_mass(1,1)/=expected_visible .or. trc_levsto(1,1)/=expected_protected) stop 5
      old_subtraction=(before_visible+before_protected)-trc_mass(1,1)-trc_levsto(1,1)
      if(case_id==1 .and. abs(old_subtraction-expected)<1.e-8_r8) stop 2
      if(case_id==2 .and. abs(old_subtraction+1._r8)>1.e-12_r8) stop 3
      if(case_id==5 .and. (trc_mass(1,1)/=0._r8 .or. trc_levsto(1,1)/=1.e-3_r8)) stop 6
   end do
   if(mode=='overdraft') then
      trc_mass=0._r8; trc_levsto=0._r8
      coefficient_uc=1._r8; gain_uc=0._r8
      flood_visible_tracer_uc=1.e-5_r8; flood_protected_tracer_uc=0._r8
      flood_visible_uc=0._r8; flood_protected_uc=0._r8
      tracer_exchange_ledger=0._r8
{transfer}
      stop 7
   endif
contains
   subroutine CoLM_stop(message)
      character(*), intent(in) :: message
      print *, message
      error stop 4
   end subroutine
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    path = tmp_path / "check_flood_ledger.f90"
    exe = tmp_path / "check_flood_ledger"
    path.write_text(program)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", str(path), "-o", str(exe)],
        capture_output=True,
        text=True,
    )
    assert built.returncode == 0, built.stderr
    ran = subprocess.run([str(exe)], capture_output=True, text=True)
    assert ran.returncode == 0, ran.stderr + ran.stdout
    overdraft = subprocess.run([str(exe), "overdraft"], capture_output=True, text=True)
    assert overdraft.returncode != 0
    assert "negative visible tracer" in overdraft.stdout
