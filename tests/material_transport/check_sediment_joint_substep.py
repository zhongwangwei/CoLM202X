#!/usr/bin/env python3
"""Extract current sediment advection; reproduce loss of constant concentration.
Usage: python3 repro_sediment_bidirectional.py REPO_ROOT [gfortran]
Exit zero means both coarse-step defects were reproduced and the fine-step
control passed. This is a bug reproduction, not an acceptance test.
All transport and donor-limit routines are extracted unchanged from REPO_ROOT;
the fixture only replaces distributed push maps with a serial three-cell chain.
"""
from pathlib import Path
import re,sys,subprocess,tempfile
import pathlib as _pl
default_root=_pl.Path(__file__).resolve().parents[2]
repo=Path(sys.argv[1]) if len(sys.argv)>1 else default_root
fc=sys.argv[2] if len(sys.argv)>2 else 'gfortran'
s=(repo/'main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90').read_text()
def routine(name):
 return re.search(r'(?ims)^   (?:RECURSIVE )?SUBROUTINE '+name+r'\s*\(.*?^   END SUBROUTINE '+name+r'\s*$',s).group(0)
stub='''module MOD_Precision
 integer,parameter::r8=kind(1.d0)
end module
module MOD_SPMD_Task
 logical::p_is_worker=.true.,p_is_io=.true.
 contains
 subroutine CoLM_stop(message)
 character(*),optional,intent(in)::message
 if(present(message)) print *,message
 error stop 2
 end subroutine
end module
module MOD_Const_Physical
 use MOD_Precision
 real(r8),parameter::grav=9.80616_r8
end module
module MOD_Grid_RiverLakeNetwork
 use MOD_Precision
 integer,parameter::numucat=3,push_next2ucat=1,push_ups2ucat=2
 integer::next(3)=[2,3,0]
 real(r8)::topo_rivwth(3)=1._r8
end module
module MOD_WorkerPushData
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork
 contains
 subroutine worker_push_data(which,src,dst,fillvalue,mode)
 integer,intent(in)::which
 real(r8),intent(in)::src(:)
 real(r8),intent(out)::dst(:)
 real(r8),intent(in)::fillvalue
 character(*),optional,intent(in)::mode
 integer::i
 dst=fillvalue
 if(which==1) then
  do i=1,numucat
   if(next(i)>0) dst(i)=src(next(i))
  enddo
 else
  do i=1,numucat
   if(next(i)>0) dst(next(i))=dst(next(i))+src(i)
  enddo
 endif
 end subroutine
end module
module sediment_extracted
 use MOD_Precision
 use MOD_SPMD_Task
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 integer::nsed=3
 real(r8)::lambda=.4_r8,psedD=2.65_r8,pwatD=1._r8
 real(r8),parameter::MAX_SED_CONC=.01_r8,SED_BEDLOAD_COEFF=17._r8
 real(r8),parameter::SED_BALANCE_ABS_TOL=1.e-10_r8,SED_BALANCE_REL_TOL=1.e-10_r8
 real(r8),allocatable::sedcon(:,:),sedsto(:,:),layer(:,:),sedout(:,:),bedout(:,:)
 real(r8),allocatable::netflw_adv_step(:,:),exch_d_adv_step(:,:),critshearvel(:,:),shearvel(:)
 contains
'''
for name in ['calc_sediment_advection','calc_sediment_advection_one_direction','ordinary_sediment_donor_demand','limit_reverse_flux','assert_sediment_mass_balance']:
 stub+=routine(name)+'\n'
stub+='''end module
program test
use sediment_extracted
use MOD_Grid_RiverLakeNetwork
implicit none
real(r8)::water(3),qsigned(3),qabs(3),before,dt
real(r8),allocatable::demand(:,:),bed_demand(:,:),joint_scale(:,:),conc_initial(:,:),bed_initial(:,:)
integer::idx,j,nsub
allocate(sedcon(nsed,3),sedsto(nsed,3),layer(nsed,3),sedout(nsed,3),bedout(nsed,3))
allocate(netflw_adv_step(nsed,3),exch_d_adv_step(nsed,3),critshearvel(nsed,3),shearvel(3))
allocate(demand(nsed,3),bed_demand(nsed,3),joint_scale(nsed,3),conc_initial(nsed,3),bed_initial(nsed,3))
water=10._r8
qsigned=0._r8
qabs=[2._r8,2._r8,0._r8]
critshearvel=1.e20_r8
shearvel=0._r8
! Gross water volume is 10m3 each way at both interfaces over a 10s period;
! net water divergence=0 in each cell. Each of the 3 grain classes starts at .001.
! Total .003 concentration is below MAX_SED_CONC=.01, even after the wrong doubling.
do idx=1,4
 sedsto=.01_r8
 sedcon=.001_r8
 layer=0._r8
 before=sum(sedsto)
 nsub=1
 if(idx==2) nsub=10
 if(idx==4) nsub=2
 dt=10._r8/nsub
 do j=1,nsub
  if(idx==3 .or. idx==4) then
   conc_initial=sedcon;bed_initial=layer
   call ordinary_sediment_donor_demand(dt,qsigned,qabs,conc_initial,bed_initial,demand,bed_demand)
   joint_scale=1._r8
   where(demand>0._r8) joint_scale=min(sedsto/max(demand,1.e-30_r8),1._r8)
   call calc_sediment_advection(dt,qsigned,qabs,water,water,conc_initial,bed_initial,joint_scale)
  else
   call calc_sediment_advection(dt,qsigned,qabs,water,water)
  endif
 enddo
 print '(A,I0,A,3ES18.9)', 'case ',idx,' grain-1 concentrations: ',sedcon(1,:)
 print '(A,3ES18.9)', 'grain-1 suspended stock: ',sedsto(1,:)
 print '(A,ES18.9)', 'solid volume residual: ',sum(sedsto)+(1._r8-lambda)*sum(layer)-before
 if(idx==3) then
  print '(A,I0,A,3ES18.9)','N01 case 3 (core fixed-scale alternative, not provider): ',idx,' C=',sedcon(1,:)
 else
  if(maxval(abs(sedcon-.001_r8))<=1.e-12_r8 .and. &
     abs(sum(sedsto)+(1._r8-lambda)*sum(layer)-before)<=1.e-12_r8) then
   print '(A,I0,A)','N01 PASS: case ',idx,' keeps uniform concentration and mass'
  else
   error stop 'N01 invariant failed'
  endif
 endif
enddo
end program
'''
with tempfile.TemporaryDirectory() as td:
 p=Path(td);(p/'test.f90').write_text(stub)
 subprocess.run([fc,'-O0','-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow',str(p/'test.f90'),'-o',str(p/'test')],check=True,cwd=p)
 subprocess.run([str(p/'test')],check=True,cwd=p)
