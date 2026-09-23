"""Compile the production levee transfer kernel with tiny module stubs."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_gross_exchange_uses_old_protected_inventory_and_conserves_classes(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    kernel = "   SUBROUTINE transfer_levee_sediment(" + source.split(
        "   SUBROUTINE transfer_levee_sediment(", 1
    )[1].split("   END SUBROUTINE transfer_levee_sediment", 1)[0] + (
        "   END SUBROUTINE transfer_levee_sediment\n"
    )
    code = """
module MOD_Grid_RiverLakeNetwork
 integer, parameter :: numucat=1
end module
module probe
 use MOD_Grid_RiverLakeNetwork
 implicit none
 integer, parameter :: r8=kind(1.d0), nsed=2
 real(r8), parameter :: lambda=0.4_r8, MAX_SED_CONC=1._r8
 logical :: p_is_worker=.true.
 real(r8) :: sedsto(nsed,1), sedsto_protected(nsed,1), layer(nsed,1)
 real(r8) :: sedbed_protected(nsed,1), sedcon(nsed,1), setvel(nsed)
 real(r8) :: sed_acc_to_protected(1), sed_acc_from_protected(1)
 real(r8) :: sed_acc_protected_area(1), sed_acc_time(1)
 real(r8) :: netflw_adv_step(nsed,1),exch_d_adv_step(nsed,1)
 contains
 subroutine assert_sediment_mass_balance(context,index,before,after,expected_change)
  character(len=*), intent(in) :: context
  integer, intent(in) :: index
  real(r8), intent(in) :: before,after,expected_change
  if (abs(after-before-expected_change)>1.e-12_r8) error stop 'class budget'
 end subroutine
""" + kernel + """
end module
program main
 use probe
 implicit none
 real(r8) :: old_protected(nsed,1)
 sedsto(:,1)=[2._r8,1._r8]
 sedsto_protected(:,1)=[0.4_r8,0.2_r8]
 old_protected=sedsto_protected
 layer=0._r8; sedbed_protected=0._r8; setvel=0._r8
 sed_acc_to_protected=2._r8; sed_acc_from_protected=2._r8
 sed_acc_protected_area=1._r8; sed_acc_time=10._r8
 netflw_adv_step=0._r8;exch_d_adv_step=0._r8
 call transfer_levee_sediment(10._r8,10._r8,[10._r8],[5._r8],[5._r8], &
  old_protected,.true.)
 if (maxval(abs(sedsto(:,1)-[1.6_r8,0.8_r8]))>1.e-12_r8) error stop 'outbound'
 ! Represent an ordinary-face export. Retreat may not use the fresh protected arrival.
 sedsto(:,1)=sedsto(:,1)-[0.4_r8,0.2_r8]
 call transfer_levee_sediment(10._r8,10._r8,[10._r8],[5._r8],[5._r8], &
  old_protected,.false.)
 if (maxval(abs(sedsto(:,1)-[1.36_r8,0.68_r8]))>1.e-12_r8) error stop 'return'
 if (maxval(abs(sedsto_protected(:,1)-[0.64_r8,0.32_r8]))>1.e-12_r8) error stop 'fresh export'
 setvel=0.1_r8; sed_acc_from_protected=0._r8
 old_protected=sedsto_protected
 call transfer_levee_sediment(10._r8,10._r8,[10._r8],[5._r8],[5._r8], &
  old_protected,.false.)
 if (minval(sedbed_protected)<=0._r8) error stop 'protected settling'
 if (maxval(abs(exch_d_adv_step+netflw_adv_step))>1.e-12_r8) error stop 'deposition diagnostic sign'
 if (minval(exch_d_adv_step)<=0._r8) error stop 'protected D diagnostic missing'
 ! BIF already consumed all old protected suspended mass before this transfer.
 ! New overtopping in the same morphology interval cannot retreat immediately.
 sedsto(:,1)=[1._r8,0._r8];sedsto_protected=0._r8
 old_protected=sedsto_protected;sedbed_protected=0._r8;setvel=0._r8
 sed_acc_to_protected=5._r8;sed_acc_from_protected=5._r8
 call transfer_levee_sediment(10._r8,10._r8,[10._r8],[5._r8],[5._r8], &
  old_protected,.true.)
 call transfer_levee_sediment(10._r8,10._r8,[10._r8],[5._r8],[5._r8], &
  old_protected,.false.)
 if (abs(sedsto(1,1)-.5_r8)>1.e-12_r8) error stop 'new overtopping re-exported'
 ! BIF has already debited .2 of an initial M=1,V=10. Overtopping Q=2
 ! must still carry .2 from the shared initial concentration, not .16.
 sedsto(:,1)=[.8_r8,0._r8];sedsto_protected=0._r8
 sed_acc_to_protected=2._r8;sed_acc_from_protected=0._r8
 call transfer_levee_sediment(10._r8,10._r8,[10._r8],[0._r8],[2._r8], &
  old_protected,.true.,reshape([1._r8,0._r8],[2,1]),old_protected, &
  reshape([1._r8,1._r8],[2,1]),reshape([1._r8,1._r8],[2,1]))
 if (abs(sedsto_protected(1,1)-.2_r8)>1.e-12_r8) error stop 'levee diluted by BIF debit'
end program
"""
    path = tmp_path / "probe.f90"
    path.write_text(code)
    compiler = require_runnable_fortran_compiler(tmp_path)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "probe"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(tmp_path / "probe")], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
