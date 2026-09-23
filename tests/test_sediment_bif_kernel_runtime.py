"""Exercise production BIF donor kernel with serial two-cell push maps."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_bif_layers_share_class_donor_and_keep_directions(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    kernels = "   real(r8) FUNCTION joint_sediment_scale(" + source.split(
        "   real(r8) FUNCTION joint_sediment_scale(", 1
    )[1].split("   END FUNCTION joint_sediment_scale", 1)[0] + (
        "   END FUNCTION joint_sediment_scale\n"
    )
    for name in ("prepare_bif_sediment", "apply_bif_sediment_credits"):
        kernels += "   SUBROUTINE " + name + "(" + source.split(
            "   SUBROUTINE " + name + "(", 1
        )[1].split("   END SUBROUTINE " + name, 1)[0] + (
            "   END SUBROUTINE " + name + "\n"
        )
    code = """
module MOD_Precision
 integer, parameter :: r8=kind(1.d0)
end module
module MOD_Const_Physical
 use MOD_Precision
 real(r8), parameter :: grav=9.81_r8
end module
module MOD_WorkerPushData
 use MOD_Precision
 type worker_pushdata_type
  integer :: mapping=0
 end type
 contains
 subroutine worker_push_data(mapping,vec_send,vec_recv,fillvalue,mode)
  type(worker_pushdata_type), intent(in) :: mapping
  real(r8), intent(in) :: vec_send(:)
  real(r8), intent(out) :: vec_recv(:)
  real(r8), optional, intent(in) :: fillvalue
  character(len=*), optional, intent(in) :: mode
  vec_recv=0._r8
  if (mapping%mapping==1) vec_recv(1)=vec_send(2)
  if (mapping%mapping==2) vec_recv(2)=vec_send(1)
 end subroutine
end module
module MOD_Grid_RiverLakeNetwork
 use MOD_Precision
 use MOD_WorkerPushData
 integer, parameter :: numucat=2, npthout_local=1, npthlev_bif=2
 integer :: pth_upst_local(1)=[1]
 real(r8) :: pth_wth(2,1)=1._r8
 type(worker_pushdata_type) :: push_bif_dn2pth=worker_pushdata_type(1)
 type(worker_pushdata_type) :: push_bif_influx=worker_pushdata_type(2)
end module
module MOD_Grid_RiverLakeLevee
 logical :: has_levee(2)=.false.
end module
module probe
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork
 implicit none
 integer, parameter :: nsed=1
 logical :: p_is_worker=.true., DEF_USE_LEVEE=.false.
 real(r8), parameter :: lambda=.4_r8, MAX_SED_CONC=1._r8
 real(r8), parameter :: SED_BEDLOAD_COEFF=17._r8
 real(r8), parameter :: SED_BALANCE_ABS_TOL=1.e-10_r8
 real(r8) :: psedD=2.65_r8,pwatD=1._r8
 real(r8) :: sedsto(1,2), sedsto_protected(1,2), layer(1,2), sedbed_protected(1,2)
 real(r8) :: sedcon(1,2), shearvel(2), critshearvel(1,2)
 real(r8) :: netflw_adv_step(1,2),exch_d_adv_step(1,2)
 real(r8) :: sed_acc_bif_forward(2,1),sed_acc_bif_reverse(2,1)
 real(r8) :: sed_acc_bif_forward_time(2,1),sed_acc_bif_reverse_time(2,1)
 contains
 subroutine CoLM_stop(message)
  character(len=*), intent(in) :: message
  print *, message
  error stop 1
 end subroutine
 subroutine assert_sediment_mass_balance(context,index,before,after,expected_change)
  character(len=*), intent(in) :: context
  integer, intent(in) :: index
  real(r8), intent(in) :: before,after,expected_change
  if (abs(after-before-expected_change)>1.e-12_r8) error stop 'BIF class budget'
 end subroutine
""" + kernels + """
end module
program main
 use probe
 use MOD_Grid_RiverLakeLevee
 implicit none
 real(r8) :: cv(1,2),cp(1,2),bv(1,2),bp(1,2)
 real(r8) :: ms(1,2),mb(1,2),lt(1,2),lf(1,2),sv(1,2),sp(1,2),bs(1,2)
 sedsto=0._r8;sedsto_protected=0._r8;layer=0._r8;sedbed_protected=0._r8
 netflw_adv_step=0._r8;exch_d_adv_step=0._r8
 shearvel=0._r8; critshearvel=1._r8
 sed_acc_bif_forward=0._r8;sed_acc_bif_reverse=0._r8
 sed_acc_bif_forward_time=0._r8;sed_acc_bif_reverse_time=0._r8
 ms=0._r8;mb=0._r8;lt=0._r8;lf=0._r8
 sedsto(1,1)=1._r8; sed_acc_bif_forward(:,1)=.8_r8
 call prepare_bif_sediment(1._r8,1._r8,[1._r8,2._r8],[0._r8,0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 if (abs(sedsto(1,1))>1.e-12_r8 .or. abs(cv(1,2)-1._r8)>1.e-12_r8) &
  error stop 'two forward layers overdrew visible donor'
 call apply_bif_sediment_credits(1._r8,[1._r8,2._r8],cv,cp,bv,bp)
 if (abs(sedsto(1,2)-1._r8)>1.e-12_r8) error stop 'forward credit'
 sedsto=0._r8;sedsto(1,2)=1._r8
 sed_acc_bif_forward=0._r8;sed_acc_bif_reverse(:,1)=.8_r8
 call prepare_bif_sediment(1._r8,1._r8,[2._r8,1._r8],[0._r8,0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 if (abs(sedsto(1,2))>1.e-12_r8 .or. abs(cv(1,1)-1._r8)>1.e-12_r8) &
  error stop 'two reverse layers overdrew downstream donor'
 ! The two BIF layers also compete for one physical bedload inventory.
 sedsto=0._r8;layer=0._r8;layer(1,1)=1._r8/(1._r8-lambda)
 shearvel(1)=2._r8;critshearvel(1,1)=.1_r8
 sed_acc_bif_reverse=0._r8;sed_acc_bif_forward(:,1)=.8_r8
 sed_acc_bif_forward_time(:,1)=1._r8
 call prepare_bif_sediment(1._r8,1._r8,[1._r8,2._r8],[0._r8,0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 if (abs(layer(1,1))>1.e-12_r8 .or. abs(bv(1,2)-1._r8)>1.e-12_r8) &
  error stop 'bedload donor capacity/stock'
 layer=0._r8;shearvel=0._r8
 has_levee=.true.; DEF_USE_LEVEE=.true.
 sedsto=0._r8;sedsto_protected=0._r8
 sedsto(1,1)=.5_r8;sedsto_protected(1,1)=.5_r8
 sed_acc_bif_forward(:,1)=.5_r8;sed_acc_bif_reverse=0._r8
 call prepare_bif_sediment(1._r8,1._r8,[1._r8,1._r8],[1._r8,1._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 if (abs(cv(1,2)-.25_r8)>1.e-12_r8 .or. abs(cp(1,2)-.25_r8)>1.e-12_r8) &
  error stop 'visible/protected destination split'
 ! Same start-of-interval M/V drives ordinary Q=2 and BIF Q=2.
 DEF_USE_LEVEE=.false.;has_levee=.false.
 sedsto=0._r8;sedsto(1,1)=1._r8;sedsto_protected=0._r8
 sed_acc_bif_forward=0._r8;sed_acc_bif_forward(1,1)=2._r8
 sed_acc_bif_reverse=0._r8;ms=0._r8;ms(1,1)=.2_r8
 call prepare_bif_sediment(1._r8,1._r8,[10._r8,10._r8],[0._r8,0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 if (abs(cv(1,2)-.2_r8)>1.e-12_r8 .or. abs(sv(1,1)-1._r8)>1.e-12_r8) &
  error stop 'main/BIF common concentration'
 if (abs(sedsto(1,1)-.8_r8)>1.e-12_r8) error stop 'BIF donor debit'
 ! Both Q=8 demands 0.8 but M=1: each must get half the stock.
 sedsto=0._r8;sedsto(1,1)=1._r8;sed_acc_bif_forward(1,1)=8._r8
 ms(1,1)=.8_r8
 call prepare_bif_sediment(1._r8,1._r8,[10._r8,10._r8],[0._r8,0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 if (abs(cv(1,2)-.5_r8)>1.e-12_r8 .or. abs(sv(1,1)-.625_r8)>1.e-12_r8) &
  error stop 'main/BIF common donor limiter'
end program
"""
    path = tmp_path / "probe.f90"
    path.write_text(code)
    compiler = require_runnable_fortran_compiler(tmp_path)
    built = subprocess.run(
        [compiler, "-cpp", "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "probe"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(tmp_path / "probe")], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
