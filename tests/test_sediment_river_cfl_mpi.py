"""Sediment push-map validation remains collective, even without CFL reductions."""

from pathlib import Path
import os
import re
import shutil
import subprocess

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def routine(source: str, name: str) -> str:
    match = re.search(rf"^\s*SUBROUTINE {name}\b.*?^\s*END SUBROUTINE {name}\b", source, re.M | re.S | re.I)
    assert match, name
    return match.group()


def test_split_river_singleton_and_empty_worker_push_groups(tmp_path: Path) -> None:
    compiler, launcher = shutil.which("mpif90"), shutil.which("mpirun") or shutil.which("mpiexec")
    if not compiler or not launcher:
        pytest.skip("MPI Fortran toolchain unavailable")
    source = SOURCE.read_text()
    calc = routine(source, "grid_sediment_calc")
    assert "dt_cfl_global" not in calc
    assert "sed_push_groups_checked = .false." in routine(source, "grid_sediment_init")
    assert "sed_push_groups_checked = .false." in routine(source, "grid_sediment_final")
    guard = "      IF (.not. sed_push_groups_checked) THEN\n" + calc.split(
        "      IF (.not. sed_push_groups_checked) THEN\n", 1
    )[1].split("#endif", 1)[0]
    code = """
module probe
 use mpi
 implicit none
 integer :: p_comm_worker, p_comm_rivsys=MPI_COMM_NULL, p_np_worker, p_iam_worker, p_err
 integer :: numucat
 logical :: rivsys_by_multiple_procs, sed_push_groups_checked=.false.
 type push_type
  integer, allocatable :: n_to_other(:), n_from_other(:)
 end type
 type(push_type) :: push_next2ucat, push_ups2ucat
 contains
 subroutine CoLM_stop(message)
  character(len=*), intent(in) :: message
  print *, message
  call MPI_Abort(MPI_COMM_WORLD, 9, p_err)
 end subroutine
 subroutine check_mapping()
  integer :: iworker,river_color
  integer, allocatable :: worker_river_color(:)
  logical :: invalid_push_group,remote_peer
""" + guard + """
 end subroutine
end module
program probe_main
 use probe
 implicit none
 character(len=8) :: mode
 call MPI_Init(p_err)
 p_comm_worker=MPI_COMM_WORLD
 call MPI_Comm_rank(p_comm_worker,p_iam_worker,p_err)
 call MPI_Comm_size(p_comm_worker,p_np_worker,p_err)
 if (p_np_worker/=4) call CoLM_stop('need four workers')
 allocate(push_next2ucat%n_to_other(0:3),push_next2ucat%n_from_other(0:3))
 allocate(push_ups2ucat%n_to_other(0:3),push_ups2ucat%n_from_other(0:3))
 push_next2ucat%n_to_other=0; push_next2ucat%n_from_other=0
 push_ups2ucat%n_to_other=0; push_ups2ucat%n_from_other=0
 call get_command_argument(1,mode)
 if (p_iam_worker<2) then
  rivsys_by_multiple_procs=.true.; numucat=1
  call MPI_Comm_split(p_comm_worker,100,p_iam_worker,p_comm_rivsys,p_err)
  if (p_iam_worker==0) then
   push_next2ucat%n_to_other(1)=1; push_ups2ucat%n_to_other(1)=1
  else
   push_next2ucat%n_from_other(0)=1; push_ups2ucat%n_from_other(0)=1
  endif
 else
  rivsys_by_multiple_procs=.false.
  numucat=merge(1,0,p_iam_worker==2)
  call MPI_Comm_split(p_comm_worker,MPI_UNDEFINED,p_iam_worker,p_comm_rivsys,p_err)
 endif
 if (trim(mode)=='bad' .and. p_iam_worker==0) push_next2ucat%n_to_other(2)=1
 call check_mapping()
 if (p_iam_worker==0) call check_mapping() ! Cached: no unmatched collective.
 call MPI_Barrier(p_comm_worker,p_err)
 sed_push_groups_checked=.false.
 call check_mapping()
 if (.not. sed_push_groups_checked) call CoLM_stop('mapping not cached')
 if (p_iam_worker==0) print *, 'SEDIMENT_PUSH_GROUP_OK'
 call MPI_Finalize(p_err)
end program
"""
    path = tmp_path / "probe.F90"
    path.write_text(code)
    built = subprocess.run([compiler, "-cpp", "-DUSEMPI", "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "probe"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stdout + built.stderr
    env = os.environ.copy()
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    ran = subprocess.run([launcher, "-n", "4", str(tmp_path / "probe")], cwd=tmp_path, env=env,
                         capture_output=True, text=True, timeout=60)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "SEDIMENT_PUSH_GROUP_OK" in ran.stdout
    bad = subprocess.run([launcher, "-n", "4", str(tmp_path / "probe"), "bad"], cwd=tmp_path, env=env,
                         capture_output=True, text=True, timeout=30)
    assert bad.returncode != 0 and "sediment river push crosses" in bad.stdout + bad.stderr


def test_advection_snapshot_and_shared_donor_limiter_across_workers(tmp_path: Path) -> None:
    compiler, launcher = shutil.which("mpif90"), shutil.which("mpirun") or shutil.which("mpiexec")
    if not compiler or not launcher:
        pytest.skip("MPI Fortran toolchain unavailable")
    source = SOURCE.read_text()
    bodies = "\n".join(routine(source, name) for name in (
        "calc_sediment_advection", "calc_sediment_advection_one_direction", "limit_reverse_flux"
    ))
    code = """
module MOD_Precision
 integer, parameter :: r8=kind(1.d0)
end module
module MOD_Const_Physical
 use MOD_Precision
 real(r8), parameter :: grav=9.81_r8
end module
module MOD_Grid_RiverLakeNetwork
 use MOD_Precision
 integer :: numucat=1, push_next2ucat=1, push_ups2ucat=2
 integer :: all_next(4), ucat_next(1)
 real(r8) :: topo_rivwth(1)=1._r8
end module
module MOD_WorkerPushData
 use mpi
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork, only: all_next, ucat_next
 contains
 subroutine worker_push_data(mapping,source,dest,fillvalue,mode)
  integer, intent(in) :: mapping
  real(r8), intent(in) :: source(:),fillvalue
  real(r8), intent(out) :: dest(:)
  character(len=*), optional, intent(in) :: mode
  real(r8) :: gathered(4)
  integer :: rank,i,ierr
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Allgather(source(1),1,MPI_DOUBLE_PRECISION,gathered,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  dest=fillvalue
  if (mapping==1) then
   if (ucat_next(1)>0) dest(1)=gathered(ucat_next(1))
  else
   do i=1,4
    if (all_next(i)==rank+1) dest(1)=dest(1)+gathered(i)
   enddo
  endif
 end subroutine
end module
module sediment_probe
 use MOD_Precision
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer :: nsed=1
 logical :: p_is_worker=.true.
 real(r8) :: lambda=0.4_r8,psedD=2.65_r8,pwatD=1._r8
 real(r8), parameter :: MAX_SED_CONC=0.01_r8,SED_BEDLOAD_COEFF=17._r8
 real(r8), parameter :: SED_BALANCE_ABS_TOL=1.e-10_r8
 real(r8) :: sedcon(1,1),sedsto(1,1),layer(1,1),sedout(1,1),bedout(1,1)
 real(r8) :: critshearvel(1,1),shearvel(1),netflw_adv_step(1,1),exch_d_adv_step(1,1)
 contains
 subroutine CoLM_stop(message)
  use mpi
  character(len=*),intent(in) :: message
  integer :: ierr
  print *, message
  call MPI_Abort(MPI_COMM_WORLD,99,ierr)
 end subroutine
 subroutine assert_sediment_mass_balance(name,i,before,after,source)
  character(len=*),intent(in) :: name
  integer,intent(in) :: i
  real(r8),intent(in) :: before,after,source
  if (abs(after-before-source)>1.e-10_r8+1.e-10_r8*max(abs(before),abs(after),abs(source))) &
   call CoLM_stop('local mass mismatch')
 end subroutine
""" + bodies + """
end module
program main
 use mpi
 use sediment_probe
 use MOD_Grid_RiverLakeNetwork
 implicit none
 integer :: rank,ierr,scenario
 real(r8) :: flow(1),water(1),global_mass,local_mass,expected_mass,dt
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 do scenario=1,5
  expected_mass=0.1_r8; dt=1._r8; water=100._r8
  layer=0._r8; shearvel=0._r8; critshearvel=1._r8
  select case(scenario)
  case(1)
   all_next=[2,3,0,0]
   sedsto=merge(0.1_r8,0._r8,rank==0)
   flow=merge(1000._r8,0._r8,rank<2)
  case(2)
   all_next=[3,3,0,0]
   sedsto=merge(0.1_r8,0._r8,rank==2)
   flow=merge(-1000._r8,0._r8,rank<2)
  case(3)
   all_next=[3,3,4,0]
   sedsto=merge(0.1_r8,0._r8,rank==2)
   flow=0._r8
   if (rank<2) flow=-1000._r8
   if (rank==2) flow=50._r8
  case(4) ! x-(x/dt)*dt is negative by 9.31e-10, but only roundoff.
   all_next=[2,3,0,0]
   expected_mass=7661368.727868479_r8; dt=3600._r8; water=1.e9_r8
   sedsto=merge(expected_mass,0._r8,rank==0)
   flow=merge(1.e12_r8,0._r8,rank==0)
  case(5) ! Bulk bed storage 1e8 similarly undershoots after / (1-lambda).
   all_next=[2,3,0,0]
   expected_mass=(1._r8-lambda)*1.e8_r8; dt=3600._r8; water=1.e9_r8
   sedsto=0._r8
   layer=merge(1.e8_r8,0._r8,rank==0)
   shearvel=merge(1000._r8,0._r8,rank==0)
   critshearvel=0._r8
   flow=merge(1.e12_r8,0._r8,rank==0)
  end select
  ucat_next=all_next(rank+1)
  sedcon(1,1)=sedsto(1,1)/water(1)
  call calc_sediment_advection(dt,flow,abs(flow),water,water)
  local_mass=sum(sedsto)+(1._r8-lambda)*sum(layer)
  call MPI_Allreduce(local_mass,global_mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (abs(global_mass-expected_mass)>1.e-12_r8*expected_mass) call CoLM_stop('global mass mismatch')
  if (scenario==1) then
   if (rank==1 .and. abs(sedsto(1,1)-0.1_r8)>1.e-12_r8) call CoLM_stop('hop one missing')
   if (rank==2 .and. abs(sedsto(1,1))>1.e-12_r8) call CoLM_stop('same-step cascade')
  elseif (scenario==2) then
   if (rank<2 .and. abs(sedsto(1,1)-0.05_r8)>1.e-12_r8) call CoLM_stop('reverse donor cap')
  elseif (scenario==3) then
   if (rank<2 .and. abs(sedsto(1,1)-0.025_r8)>1.e-12_r8) call CoLM_stop('mixed reverse cap')
   if (rank==3 .and. abs(sedsto(1,1)-0.05_r8)>1.e-12_r8) call CoLM_stop('mixed forward flux')
  elseif (scenario==4 .or. scenario==5) then
   if (rank==0 .and. (sedsto(1,1)<0._r8 .or. layer(1,1)<0._r8)) &
    call CoLM_stop('negative donor after roundoff clamp')
  endif
 enddo
 if (rank==0) print *, 'SEDIMENT_MPI_SNAPSHOT_OK'
 call MPI_Finalize(ierr)
end program
"""
    path = tmp_path / "probe.F90"
    path.write_text(code)
    built = subprocess.run([compiler, "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "probe"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stdout + built.stderr
    env = os.environ.copy()
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    ran = subprocess.run([launcher, "-n", "4", str(tmp_path / "probe")], cwd=tmp_path, env=env,
                         capture_output=True, text=True, timeout=60)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "SEDIMENT_MPI_SNAPSHOT_OK" in ran.stdout


def test_donor_roundoff_clamps_only_machine_precision_errors(tmp_path: Path) -> None:
    compiler = shutil.which("gfortran")
    if not compiler:
        pytest.skip("Fortran compiler unavailable")
    source = SOURCE.read_text()
    update = source.split("         sed_roundoff =", 1)[1].split(
        "         layer(:,i) = max(layer(:,i), 0._r8)", 1
    )[0]
    update = "         sed_roundoff =" + update + "         layer(:,i) = max(layer(:,i), 0._r8)\n"
    code = """
program probe
 implicit none
 integer, parameter :: r8=kind(1.d0), nsed=1, i=1
 real(r8), parameter :: SED_BALANCE_ABS_TOL=1.e-10_r8
 real(r8) :: sedsto(1,1),layer(1,1),sedout(1,1),bedout(1,1)
 real(r8) :: sed_ups(1,1),bed_ups(1,1),sed_roundoff(1),bed_roundoff(1)
 real(r8) :: dt,lambda
 character(len=16) :: mode
 call get_command_argument(1,mode)
 lambda=0.4_r8; dt=3600._r8
 sedsto=0._r8; layer=0._r8; sedout=0._r8; bedout=0._r8
 sed_ups=0._r8; bed_ups=0._r8
 select case(trim(mode))
 case('suspended')
  sedsto=7661368.727868479_r8; sedout=sedsto/dt
 case('bed')
  layer=1.e8_r8; bedout=(1._r8-lambda)*layer/dt
 case('overdraft')
  sedsto=1._r8; sedout=1.01_r8/dt
 end select
""" + update + """
 if (minval(sedsto)<0._r8 .or. minval(layer)<0._r8) error stop 3
 print *, 'ROUNDING_ONLY_OK'
contains
 subroutine CoLM_stop(message)
  character(len=*),intent(in) :: message
  print *, message
  error stop 9
 end subroutine
end program
"""
    path = tmp_path / "roundoff.F90"
    path.write_text(code)
    built = subprocess.run([compiler, "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "roundoff"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stdout + built.stderr
    for mode in ("suspended", "bed"):
        ran = subprocess.run([str(tmp_path / "roundoff"), mode], capture_output=True, text=True, timeout=30)
        assert ran.returncode == 0 and "ROUNDING_ONLY_OK" in ran.stdout, ran.stdout + ran.stderr
    bad = subprocess.run([str(tmp_path / "roundoff"), "overdraft"], capture_output=True, text=True, timeout=30)
    assert bad.returncode != 0 and "sediment donor limiter produced negative inventory" in bad.stdout + bad.stderr
