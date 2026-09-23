"""Exercise the production morphology BIF kernel with cross-rank path maps."""

from pathlib import Path
import os
import shutil
import subprocess

import pytest


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_cross_rank_shared_forward_and_reverse_donors(tmp_path: Path) -> None:
    compiler = shutil.which("mpif90")
    launcher = shutil.which("mpirun") or shutil.which("mpiexec")
    if not compiler or not launcher:
        pytest.skip("MPI Fortran toolchain unavailable")
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
    code = r"""
module MOD_Precision
 integer, parameter :: r8=kind(1.d0)
end module
module MOD_Const_Physical
 use MOD_Precision
 real(r8), parameter :: grav=9.81_r8
end module
module MOD_WorkerPushData
 use MOD_Precision
 use mpi
 type worker_pushdata_type
  integer :: mapping=0
 end type
 integer :: test_rank, scenario=1, mpierr
 contains
 subroutine worker_push_data(mapping,vec_send,vec_recv,fillvalue,mode)
  type(worker_pushdata_type), intent(in) :: mapping
  real(r8), intent(in) :: vec_send(:)
  real(r8), intent(out) :: vec_recv(:)
  real(r8), optional, intent(in) :: fillvalue
  character(len=*), optional, intent(in) :: mode
  real(r8) :: one(3), one_all(3), two(2), two_all(2)
  integer :: dst2
  vec_recv=0._r8
  if (mapping%mapping==1) then
   one=0._r8
   if (size(vec_send)>0) one(test_rank+1)=vec_send(1)
   call MPI_Allreduce(one,one_all,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
   if (test_rank==0) then
    dst2=merge(2,1,scenario==1)
    vec_recv(1)=one_all(2);vec_recv(2)=one_all(dst2+1)
   endif
  else
   two=0._r8
   if (test_rank==0) two(:)=vec_send(:)
   call MPI_Allreduce(two,two_all,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
   if (test_rank==1) vec_recv(1)=two_all(1)+merge(two_all(2),0._r8,scenario==2)
   if (test_rank==2 .and. scenario==1) vec_recv(1)=two_all(2)
  endif
 end subroutine
end module
module MOD_Grid_RiverLakeNetwork
 use MOD_Precision
 use MOD_WorkerPushData
 integer :: numucat=1,npthout_local,npthlev_bif=1
 integer, allocatable :: pth_upst_local(:)
 real(r8), allocatable :: pth_wth(:,:)
 type(worker_pushdata_type) :: push_bif_dn2pth=worker_pushdata_type(1)
 type(worker_pushdata_type) :: push_bif_influx=worker_pushdata_type(2)
end module
module MOD_Grid_RiverLakeLevee
 logical :: has_levee(1)=.false.
end module
module probe
 use mpi
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork
 implicit none
 integer, parameter :: nsed=1
 integer :: p_comm_worker=MPI_COMM_WORLD,p_err
 logical :: p_is_worker=.true., DEF_USE_LEVEE=.false.
 real(r8), parameter :: lambda=.4_r8, MAX_SED_CONC=1._r8
 real(r8), parameter :: SED_BEDLOAD_COEFF=17._r8,SED_BALANCE_ABS_TOL=1.e-10_r8
 real(r8) :: psedD=2.65_r8,pwatD=1._r8
 real(r8) :: sedsto(1,1),sedsto_protected(1,1),layer(1,1),sedbed_protected(1,1)
 real(r8) :: sedcon(1,1),shearvel(1),critshearvel(1,1)
 real(r8) :: netflw_adv_step(1,1),exch_d_adv_step(1,1)
 real(r8), allocatable :: sed_acc_bif_forward(:,:),sed_acc_bif_reverse(:,:)
 real(r8), allocatable :: sed_acc_bif_forward_time(:,:),sed_acc_bif_reverse_time(:,:)
 contains
 subroutine CoLM_stop(message)
  character(len=*), intent(in) :: message
  print *,message
  call MPI_Abort(MPI_COMM_WORLD,9,p_err)
 end subroutine
 subroutine assert_sediment_mass_balance(context,index,before,after,expected_change)
  character(len=*), intent(in) :: context
  integer, intent(in) :: index
  real(r8), intent(in) :: before,after,expected_change
  if (abs(after-before-expected_change)>1.e-10_r8) call CoLM_stop('MPI BIF class budget')
 end subroutine
""" + kernels + r"""
end module
program main
 use probe
 use MOD_WorkerPushData, only: test_rank,scenario
 implicit none
 integer :: nproc
 real(r8) :: cv(1,1),cp(1,1),bv(1,1),bp(1,1)
 real(r8) :: ms(1,1),mb(1,1),lt(1,1),lf(1,1),sv(1,1),sp(1,1),bs(1,1)
 call MPI_Init(p_err)
 call MPI_Comm_rank(MPI_COMM_WORLD,test_rank,p_err)
 call MPI_Comm_size(MPI_COMM_WORLD,nproc,p_err)
 if (nproc/=3) call CoLM_stop('need three ranks')
 npthout_local=merge(2,0,test_rank==0)
 allocate(pth_upst_local(npthout_local),pth_wth(1,npthout_local))
 allocate(sed_acc_bif_forward(1,npthout_local),sed_acc_bif_reverse(1,npthout_local))
 allocate(sed_acc_bif_forward_time(1,npthout_local),sed_acc_bif_reverse_time(1,npthout_local))
 pth_upst_local=1;pth_wth=1._r8
 sedsto_protected=0._r8;layer=0._r8;sedbed_protected=0._r8
 netflw_adv_step=0._r8;exch_d_adv_step=0._r8
 shearvel=0._r8;critshearvel=1._r8
 sed_acc_bif_forward_time=0._r8;sed_acc_bif_reverse_time=0._r8
 ms=0._r8;mb=0._r8;lt=0._r8;lf=0._r8
 scenario=1;sedsto=0._r8;sed_acc_bif_reverse=0._r8
 if (test_rank==0) then
  sedsto(1,1)=1._r8;sed_acc_bif_forward=.8_r8
 endif
 call prepare_bif_sediment(1._r8,1._r8,[1._r8],[0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 call apply_bif_sediment_credits(1._r8,[1._r8],cv,cp,bv,bp)
 if (test_rank==0 .and. abs(sedsto(1,1))>1.e-12_r8) call CoLM_stop('forward donor overdraw')
 if (test_rank/=0 .and. abs(sedsto(1,1)-.5_r8)>1.e-12_r8) call CoLM_stop('forward fanout')
 scenario=2;sedsto=0._r8;sed_acc_bif_forward=0._r8
 if (test_rank==1) sedsto(1,1)=1._r8
 if (test_rank==0) sed_acc_bif_reverse=.8_r8
 call prepare_bif_sediment(1._r8,1._r8,[1._r8],[0._r8],ms,mb,lt,lf,sv,sp,bs,cv,cp,bv,bp)
 call apply_bif_sediment_credits(1._r8,[1._r8],cv,cp,bv,bp)
 if (test_rank==0 .and. abs(sedsto(1,1)-1._r8)>1.e-12_r8) call CoLM_stop('reverse fanin')
 if (test_rank==1 .and. abs(sedsto(1,1))>1.e-12_r8) call CoLM_stop('reverse donor overdraw')
 if (test_rank==2 .and. abs(sedsto(1,1))>1.e-12_r8) call CoLM_stop('empty worker')
 if (test_rank==0) print *,'SEDIMENT_BIF_MPI_OK'
 call MPI_Finalize(p_err)
end program
"""
    path = tmp_path / "probe.F90"
    path.write_text(code)
    built = subprocess.run(
        [compiler, "-cpp", "-DUSEMPI", "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "probe"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    env = os.environ.copy()
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    ran = subprocess.run(
        [launcher, "-n", "3", str(tmp_path / "probe")], cwd=tmp_path, env=env,
        capture_output=True, text=True, timeout=60,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "SEDIMENT_BIF_MPI_OK" in ran.stdout
