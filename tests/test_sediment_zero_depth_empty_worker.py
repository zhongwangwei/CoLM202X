"""Compile the shipped sediment branches that used to fail on dry/empty workers."""

from pathlib import Path
import re
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def routine(source: str, name: str) -> str:
    match = re.search(
        rf"^\s*SUBROUTINE {name}\b.*?^\s*END SUBROUTINE {name}\b",
        source, re.M | re.S | re.I,
    )
    assert match is not None, name
    return match.group()


def run_probe(tmp_path: Path, code: str) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    source = tmp_path / "probe.f90"
    source.write_text(code)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", str(source), "-o", "probe"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run(["./probe"], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr


def test_legacy_restart_with_pending_water_window_requires_endpoints(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    gate = re.search(
        r"      IF \(schema_version < 2\) THEN.*?      ENDIF\n\n"
        r"      CALL require_sediment_restart_var\(file_restart, ncid, 'sed_acc_rivout'",
        source, re.S,
    )
    assert gate is not None
    body = gate.group().split("      CALL require_sediment_restart_var", 1)[0]
    body = re.sub(r"#ifdef USEMPI\n.*?#endif\n", "", body, flags=re.S)
    code = """
module probe
 implicit none
 integer, parameter :: r8=kind(1.d0)
 integer :: schema_version, ncid=0, ierr
 real(r8) :: sed_acc_time(1)
 logical :: p_is_worker=.true., p_is_io=.false., p_is_master=.false.
 logical :: legacy_nonzero, rejected
 contains
 subroutine CoLM_stop()
  rejected=.true.
 end subroutine
 integer function nf90_close(handle)
  integer, intent(in) :: handle
  nf90_close=0
 end function
 subroutine check(version, pending, must_reject)
  integer, intent(in) :: version
  logical, intent(in) :: pending, must_reject
  schema_version=version; rejected=.false.
  sed_acc_time=merge(1._r8,0._r8,pending)
""" + body + """
  if (rejected .neqv. must_reject) error stop 1
 end subroutine
end module
program main
 use probe
 call check(0,.false.,.false.)
 call check(1,.false.,.false.)
 call check(0,.true.,.true.)
 call check(1,.true.,.true.)
 call check(2,.true.,.false.)
end program
"""
    run_probe(tmp_path, code)


def test_zero_ignore_depth_is_finite_and_conserves_input(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    code = """
module MOD_Precision
 integer, parameter :: r8=kind(1.d0)
end module
module MOD_Grid_RiverLakeNetwork
 integer :: numucat=1
end module
module sediment_probe
 use MOD_Precision
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer :: nsed=1
 logical :: p_is_worker=.true.
 real(r8) :: lambda=0.4_r8, sed_ignore_dph=0._r8
 real(r8), parameter :: MAX_SED_CONC=0.01_r8, EXCH_SHEARVEL_MIN=1.e-4_r8
 real(r8), parameter :: EXCH_SHEARVEL_BLEND=2.e-4_r8, EXCH_ZD_MAX=100._r8
 real(r8) :: vonKar=0.4_r8
 real(r8) :: sedinp(1,1), sedsto(1,1), layer(1,1), sedcon(1,1)
 real(r8) :: netflw(1,1), exch_d_eff(1,1), exch_es_eff(1,1)
 real(r8) :: exch_d_raw(1,1), exch_es_raw(1,1), susvel(1,1)
 real(r8) :: setvel(1), shearvel(1)
 contains
 subroutine assert_sediment_mass_balance(name,i,before,after,source)
  character(len=*), intent(in) :: name
  integer, intent(in) :: i
  real(r8), intent(in) :: before,after,source
  if (abs(after-before-source)>1.e-12_r8) error stop 9
 end subroutine
""" + routine(source, "calc_sediment_exchange") + "\n" + routine(source, "apply_sediment_input") + """
end module
program probe
 use sediment_probe
 implicit none
 real(r8) :: water(1), area(1)
 water=0._r8; area=1._r8
 sedsto=0._r8; layer=0._r8; sedcon=0._r8; netflw=0._r8
 susvel=0._r8; setvel=0._r8; shearvel=0._r8; sedinp=0.12_r8
 call calc_sediment_exchange(1._r8,water,area)
 call apply_sediment_input(1._r8,water,area)
 if (.not.ieee_is_finite(sedcon(1,1))) error stop 1
 if (sedsto(1,1)/=0._r8 .or. abs((1._r8-lambda)*layer(1,1)-0.12_r8)>1.e-12_r8) error stop 2
 if (abs(netflw(1,1)+0.12_r8)>1.e-12_r8) error stop 3
 water=10._r8; sedinp=0.05_r8; netflw=0._r8
 call apply_sediment_input(1._r8,water,area)
 if (abs(sedsto(1,1)-0.05_r8)>1.e-12_r8) error stop 4
 if (abs(sedcon(1,1)-0.005_r8)>1.e-12_r8) error stop 5
end program
"""
    run_probe(tmp_path, code)


def test_empty_worker_writes_static_restart_metadata_without_arrays(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    writer = routine(source, "write_sediment_restart")
    start = writer.index("      DO ised = 1, nsed")
    end = writer.index("      DO ised = 1, nsed", start + 1)
    static_writes = writer[start:end]
    code = """
module probe_module
 implicit none
 integer, parameter :: r8=kind(1.d0)
 integer :: nsed=2, nlfp_sed=1, numucat=0, totalnumucat=1
 integer :: ucat_data_address(1)=0, n_calls=0
 logical :: p_is_worker=.true.
 real(r8) :: sDiam(2)=1._r8, setvel(2)=1._r8
 real(r8), allocatable :: sed_frc(:,:), sed_slope(:,:), topo_rivwth(:), topo_rivlen(:)
 contains
 subroutine write_sediment_scalar_meta(file_restart,vname,value)
  character(len=*), intent(in) :: file_restart,vname
  real(r8), intent(in) :: value
 end subroutine
 subroutine vector_gather_and_write(data,vlen,total,address,file_restart,vname,dimname)
  real(r8), intent(in) :: data(:)
  integer, intent(in) :: vlen,total,address(:)
  character(len=*), intent(in) :: file_restart,vname,dimname
  if (vlen/=numucat .or. size(data)<vlen) error stop 1
  n_calls=n_calls+1
 end subroutine
 subroutine metadata()
  integer :: ised,ilyr
  character(len=16) :: cised,cilyr
  character(len=16) :: file_restart='unused'
  real(r8) :: dummy_sed(1)
""" + static_writes + """
 end subroutine
end module
program probe
 use probe_module
 call metadata()
 if (n_calls/=5) error stop 2
 allocate(sed_frc(2,1),sed_slope(1,1),topo_rivwth(1),topo_rivlen(1))
 numucat=1; n_calls=0
 call metadata()
 if (n_calls/=5) error stop 3
end program
"""
    run_probe(tmp_path, code)


def test_advection_snapshot_limiter_and_empty_worker(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    bodies = "\n".join(
        routine(source, name) for name in (
            "calc_sediment_advection", "calc_sediment_advection_one_direction", "limit_reverse_flux"
        )
    )
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
 integer :: numucat=0, push_next2ucat=1, push_ups2ucat=2
 integer, allocatable :: ucat_next(:)
 real(r8), allocatable :: topo_rivwth(:)
end module
module MOD_WorkerPushData
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next
 contains
 subroutine worker_push_data(mapping,source,dest,fillvalue,mode)
  integer, intent(in) :: mapping
  real(r8), intent(in) :: source(:), fillvalue
  real(r8), intent(out) :: dest(:)
  character(len=*), optional, intent(in) :: mode
  integer :: i,j
  dest=fillvalue
  do i=1,numucat
   j=ucat_next(i)
   if (j<=0) cycle
   if (mapping==1) then
    dest(i)=source(j)
   else
    dest(j)=dest(j)+source(i)
   endif
  enddo
 end subroutine
end module
module sediment_probe
 use MOD_Precision
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer :: nsed=1
 logical :: p_is_worker=.true.
 real(r8) :: lambda=0.4_r8, psedD=2.65_r8, pwatD=1._r8
 real(r8), parameter :: MAX_SED_CONC=0.01_r8, SED_BEDLOAD_COEFF=17._r8
 real(r8), parameter :: SED_BALANCE_ABS_TOL=1.e-10_r8
 real(r8), allocatable :: sedcon(:,:), sedsto(:,:), layer(:,:), sedout(:,:), bedout(:,:)
 real(r8), allocatable :: critshearvel(:,:), shearvel(:), netflw_adv_step(:,:), exch_d_adv_step(:,:)
 contains
 subroutine CoLM_stop(message)
  character(len=*), intent(in) :: message
  print *, message
  error stop 99
 end subroutine
 subroutine assert_sediment_mass_balance(name,i,before,after,source)
  character(len=*), intent(in) :: name
  integer, intent(in) :: i
  real(r8), intent(in) :: before,after,source
  if (abs(after-before-source)>1.e-10_r8) error stop 9
 end subroutine
""" + bodies + """
end module
program probe
 use sediment_probe
 use MOD_Grid_RiverLakeNetwork
 implicit none
 integer :: n
 real(r8), allocatable :: flow(:), water(:), water_start(:)
 ! Empty and resized domains must still participate in all push calls.
 do n=0,2
  call setup(n)
  sedsto=0.05_r8; sedcon=0.005_r8; layer=0._r8
  flow=0._r8; water=10._r8
  call calc_sediment_advection(1._r8,flow,flow,water,water)
  if (any(abs(sedsto-0.05_r8)>1.e-12_r8)) error stop 1
 enddo
 ! A two-hop chain may advance only one hop in a morphology interval.
 call setup(3)
 ucat_next=[2,3,0]; water=100._r8
 sedsto(1,:)=[0.1_r8,0._r8,0._r8]; sedcon(1,:)=[0.001_r8,0._r8,0._r8]
 flow=[1000._r8,1000._r8,0._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(sedout(1,1)-0.1_r8)>1.e-12_r8 .or. abs(sedout(1,2))>1.e-12_r8) error stop 2
 if (abs(sedsto(1,2)-0.1_r8)>1.e-12_r8 .or. abs(sedsto(1,3))>1.e-12_r8) error stop 3
 ! A nearly dry high-throughflow donor is inventory-limited, not CFL-aborted.
 call setup(2)
 ucat_next=[2,0]; water=[1.e-9_r8,100._r8]
 sedsto(1,:)=[0.001_r8,0._r8]; sedcon(1,:)=[1.e6_r8,0._r8]
 flow=[1000._r8,0._r8]
 call calc_sediment_advection(3600._r8,flow,abs(flow),water,water)
 if (any(sedsto<0._r8) .or. abs(sum(sedsto)+(1._r8-lambda)*sum(layer)-0.001_r8)>1.e-12_r8) error stop 4
 ! Two reverse edges share the same donor, and their receipt is not re-exported.
 call setup(3)
 ucat_next=[3,3,0]; water=100._r8
 sedsto(1,:)=[0._r8,0._r8,0.1_r8]; sedcon(1,:)=[0._r8,0._r8,0.001_r8]
 flow=[-1000._r8,-1000._r8,0._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(sedout(1,1)+0.05_r8)>1.e-12_r8 .or. abs(sedout(1,2)+0.05_r8)>1.e-12_r8) error stop 5
 if (abs(sum(sedsto)-0.1_r8)>1.e-12_r8) error stop 6
 ! A forward edge and two reverse edges compete for the same initial stock.
 call setup(4)
 ucat_next=[3,3,4,0]; water=100._r8
 sedsto(1,:)=[0._r8,0._r8,0.1_r8,0._r8]
 sedcon(1,:)=[0._r8,0._r8,0.001_r8,0._r8]
 flow=[-1000._r8,-1000._r8,50._r8,0._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(sedout(1,1)+0.025_r8)>1.e-12_r8 .or. &
     abs(sedout(1,2)+0.025_r8) >1.e-12_r8 .or. &
     abs(sedout(1,3)-0.05_r8)>1.e-12_r8) error stop 7
 if (abs(sum(sedsto)-0.1_r8)>1.e-12_r8 .or. any(sedsto<0._r8)) error stop 8
 ! A shrinking carrier deposits excess before computing the donor flux.
 call setup(1)
 ucat_next=[-9]; water=[1._r8]
 sedsto(1,1)=0.1_r8; sedcon(1,1)=0.1_r8; flow=[1000._r8]
 layer(1,1)=0.1_r8/(1._r8-lambda); shearvel(1)=10._r8; critshearvel(1,1)=0.01_r8
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(sedout(1,1)-0.01_r8)>1.e-12_r8) error stop 16
 if (abs(bedout(1,1)-0.1_r8)>1.e-12_r8) error stop 19
 if (abs((1._r8-lambda)*layer(1,1)-0.09_r8)>1.e-12_r8) error stop 17
 if (abs(netflw_adv_step(1,1)+0.09_r8)>1.e-12_r8 .or. &
     abs(exch_d_adv_step(1,1)-0.09_r8)>1.e-12_r8) error stop 18
 ! Bedload from two reverse faces must share the donor's solid bed stock.
 call setup(3)
 ucat_next=[3,3,0]; water=100._r8
 sedsto=0._r8; sedcon=0._r8; layer(1,3)=1._r8/(1._r8-lambda)
 shearvel(3)=10._r8; critshearvel(1,3)=0.01_r8
 flow=[-1000._r8,-1000._r8,0._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(bedout(1,1)+0.5_r8)>1.e-12_r8 .or. abs(bedout(1,2)+0.5_r8)>1.e-12_r8) error stop 13
 if (abs((1._r8-lambda)*sum(layer)-1._r8)>1.e-12_r8 .or. any(layer<0._r8)) error stop 14
 ! Both outlet markers export exactly their limited face flux.
 call setup(2)
 ucat_next=[-9,-10]; water=100._r8
 sedsto=0.1_r8; sedcon=0.001_r8; flow=1000._r8
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(sum(sedsto)+sum(sedout)-0.2_r8)>1.e-12_r8) error stop 15
 ! Received mass beyond the current carrier's cap deposits once and is credited.
 call setup(2)
 ucat_next=[2,0]; water=[100._r8,1._r8]
 sedsto(1,:)=[0.1_r8,0._r8]; sedcon(1,:)=[0.001_r8,0._r8]
 flow=[1000._r8,0._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water,water)
 if (abs(sedsto(1,2)-0.01_r8)>1.e-12_r8) error stop 10
 if (abs((1._r8-lambda)*layer(1,2)-0.09_r8)>1.e-12_r8) error stop 11
 if (abs(netflw_adv_step(1,2)+0.09_r8)>1.e-12_r8 .or. &
     abs(exch_d_adv_step(1,2)-0.09_r8)>1.e-12_r8) error stop 12
 ! Falling donor 100->10 exports 0.9 at its initial 1% concentration;
 ! rising receiver 10->100 keeps the 0.9 arrival below its final 1% cap.
 call setup(2)
 ucat_next=[2,0]; water_start=[100._r8,10._r8]; water=[10._r8,100._r8]
 sedsto(1,:)=[1._r8,0._r8]; sedcon=0._r8; flow=[90._r8,0._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water_start,water)
 if (abs(sedout(1,1)-0.9_r8)>1.e-12_r8) error stop 20
 if (abs(sedsto(1,1)-0.1_r8)>1.e-12_r8 .or. abs(sedsto(1,2)-0.9_r8)>1.e-12_r8) error stop 21
 if (abs(sum(layer))>1.e-12_r8 .or. abs(sum(sedsto)-1._r8)>1.e-12_r8) error stop 22
 call setup(1)
 ucat_next=[-9]; water_start=[100._r8]; water=[10._r8]
 sedsto(1,1)=1._r8; sedcon=0._r8; flow=[90._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water_start,water)
 if (abs(sedout(1,1)-0.9_r8)>1.e-12_r8 .or. &
     abs(sedsto(1,1)-0.1_r8)>1.e-12_r8 .or. abs(sum(layer))>1.e-12_r8) error stop 23
 ! Splitting the same volume trajectory into two morphology intervals keeps
 ! the donor/end time levels paired and the integrated export unchanged.
 sedsto(1,1)=1._r8; layer=0._r8; water_start=[100._r8]; water=[55._r8]
 call calc_sediment_advection(0.5_r8,flow,abs(flow),water_start,water)
 if (abs(sedsto(1,1)-0.55_r8)>1.e-12_r8) error stop 24
 water_start=[55._r8]; water=[10._r8]
 call calc_sediment_advection(0.5_r8,flow,abs(flow),water_start,water)
 if (abs(sedsto(1,1)-0.1_r8)>1.e-12_r8 .or. abs(sum(layer))>1.e-12_r8) error stop 25
 ! M=1,V=10; BIF already debited Q=2 times C=.1. The ordinary
 ! Q=2 face still exports .2 from the shared initial concentration.
 call setup(1)
 ucat_next=[-9];water_start=[10._r8];water=[10._r8]
 sedsto(1,1)=.8_r8;sedcon(1,1)=.08_r8;flow=[2._r8]
 call calc_sediment_advection(1._r8,flow,abs(flow),water_start,water, &
  reshape([.1_r8],[1,1]),reshape([0._r8],[1,1]), &
  reshape([1._r8],[1,1]),reshape([1._r8],[1,1]))
 if (abs(sedout(1,1)-.2_r8)>1.e-12_r8) error stop 26
 contains
 subroutine setup(ncell)
  integer, intent(in) :: ncell
  numucat=ncell
  if (allocated(sedcon)) deallocate(sedcon,sedsto,layer,sedout,bedout,critshearvel, &
    shearvel,netflw_adv_step,exch_d_adv_step,topo_rivwth,ucat_next,flow,water,water_start)
  allocate(sedcon(1,ncell),sedsto(1,ncell),layer(1,ncell),sedout(1,ncell),bedout(1,ncell),critshearvel(1,ncell), &
    shearvel(ncell),netflw_adv_step(1,ncell),exch_d_adv_step(1,ncell),topo_rivwth(ncell),ucat_next(ncell), &
    flow(ncell),water(ncell),water_start(ncell))
  ucat_next=0; topo_rivwth=1._r8; shearvel=0._r8; critshearvel=1._r8; layer=0._r8
 end subroutine
end program
"""
    run_probe(tmp_path, code)
