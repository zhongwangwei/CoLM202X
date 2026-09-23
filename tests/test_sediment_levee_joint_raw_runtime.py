"""Gross levee demand must compete with ordinary faces before clipping."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_gross_exceeding_carrier_shares_donor_with_main_face(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    demand = "            levee_to_demand = 0._r8; levee_from_demand = 0._r8" + source.split(
        "            levee_to_demand = 0._r8; levee_from_demand = 0._r8", 1
    )[1].split("            IF (DEF_USE_BIFURCATION) THEN", 1)[0]
    scale = "   real(r8) FUNCTION joint_sediment_scale(" + source.split(
        "   real(r8) FUNCTION joint_sediment_scale(", 1
    )[1].split("   END FUNCTION joint_sediment_scale", 1)[0] + (
        "   END FUNCTION joint_sediment_scale\n"
    )
    transfer = "   SUBROUTINE transfer_levee_sediment(" + source.split(
        "   SUBROUTINE transfer_levee_sediment(", 1
    )[1].split("   END SUBROUTINE transfer_levee_sediment", 1)[0] + (
        "   END SUBROUTINE transfer_levee_sediment\n"
    )
    code = r"""
module MOD_Grid_RiverLakeNetwork
 integer,parameter :: numucat=1
end module
module probe
 use MOD_Grid_RiverLakeNetwork
 implicit none
 integer,parameter :: r8=kind(1.d0),nsed=1
 real(r8),parameter :: lambda=.4_r8,MAX_SED_CONC=1._r8
 logical :: p_is_worker=.true.,DEF_USE_LEVEE=.true.
 real(r8) :: donor_visible(1,1),donor_protected(1,1),rivsto_donor(1),protected_start(1)
 real(r8) :: levee_to_demand(1,1),levee_from_demand(1,1),dt_morph,deltime
 real(r8) :: sedsto(1,1),sedsto_protected(1,1),layer(1,1),sedbed_protected(1,1)
 real(r8) :: sed_acc_to_protected(1),sed_acc_from_protected(1)
 real(r8) :: sed_acc_protected_area(1),sed_acc_time(1),sedcon(1,1),setvel(1)
 real(r8) :: netflw_adv_step(1,1),exch_d_adv_step(1,1)
 contains
 subroutine assert_sediment_mass_balance(context,index,before,after,expected_change)
  character(len=*),intent(in) :: context
  integer,intent(in) :: index
  real(r8),intent(in) :: before,after,expected_change
  if(abs(after-before-expected_change)>1.e-12_r8) error stop 'class budget'
 end subroutine
 subroutine compute_levee_demand()
  integer :: i
""" + demand + r"""
 end subroutine
""" + scale + transfer + r"""
end module
program main
 use probe
 implicit none
 real(r8) :: common_scale(1,1),old_protected(1,1)
 dt_morph=1._r8;deltime=1._r8
 donor_visible=1._r8;donor_protected=0._r8
 rivsto_donor=1._r8;protected_start=0._r8
 sed_acc_to_protected=2._r8;sed_acc_from_protected=0._r8
 sedsto=1._r8;sedsto_protected=0._r8;old_protected=0._r8
 layer=0._r8;sedbed_protected=0._r8;setvel=0._r8;sedcon=1._r8
 sed_acc_protected_area=0._r8;sed_acc_time=1._r8
 netflw_adv_step=0._r8;exch_d_adv_step=0._r8
 call compute_levee_demand()
 if(abs(levee_to_demand(1,1)-2._r8)>1.e-12_r8) error stop 'raw levee demand clipped'
 common_scale=joint_sediment_scale(1._r8,1._r8+levee_to_demand(1,1))
 call transfer_levee_sediment(1._r8,1._r8,[1._r8],[0._r8],[0._r8], &
  old_protected,.true.,donor_visible,donor_protected,common_scale,common_scale)
 if(abs(sedsto_protected(1,1)-2._r8/3._r8)>1.e-12_r8) error stop 'levee unfair share'
 sedsto(1,1)=sedsto(1,1)-min(sedsto(1,1),1._r8*common_scale(1,1))
 if(abs(sedsto(1,1))>1.e-12_r8) error stop 'main unfair share'
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
