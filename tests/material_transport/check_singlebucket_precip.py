#!/usr/bin/env python3
"""Exercise production canopy transport with lightweight host stubs (gfortran)."""
import os
import pathlib
import re
import subprocess
import sys
import tempfile

root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[2]
def routine(path, name):
    text = (root / path).read_text()
    return re.search(rf'^   SUBROUTINE {name} \(.*?^   END SUBROUTINE {name}', text, re.M | re.S).group()

source = '''module fixture
implicit none
integer, parameter :: r8=kind(1.d0), ntracers=1
real(r8), parameter :: trc_tiny=1.e-20_r8, trc_water_min_for_ratio=1.e-12_r8
logical :: DEF_VEG_SNOW=.false., limited=.false.
type tracer
real(r8) :: max_dissolved_conc=2._r8
end type
type(tracer) :: tracers(1)
real(r8) :: trc_ldew_rain(1,1), trc_ldew_snow(1,1), a_trc_precip(1,1)
real(r8) :: trc_pg_rain_ground(1,1), trc_pg_snow_ground(1,1)
real(r8), allocatable :: trc_waterstorage(:,:)
real(r8) :: trc_canopy_solid(1,1), trc_surface_solid(1,1), trc_waterstorage_solid(1,1)
contains
logical function mt_uses_generic_water(i)
integer :: i
mt_uses_generic_water=.true.
end function
logical function mt_has_dissolved_limit(i)
integer :: i
mt_has_dissolved_limit=limited
end function
real(r8) function tracer_forcing_rain_value(i,p)
integer :: i,p
tracer_forcing_rain_value=5._r8
end function
real(r8) function tracer_forcing_snow_value(i,p)
integer :: i,p
tracer_forcing_snow_value=7._r8
end function
'''
source += routine('main/MATERIAL/FAMILY_SOLUTE/MOD_MT_Solute.F90', 'mt_equilibrate_dissolved')+'\n'
source += routine('main/MATERIAL/GENERIC_WATER/MOD_MT_Precip.F90', 'mt_precip')+'\n'
source += '''end module
program check
use fixture
implicit none
integer :: failures=0
allocate(trc_waterstorage(1,1))
! Old water is 10 with nonuniform rain/snow inventories 12 and 8 (R=2 combined).
call run('partial snow, warm mean leaf',10.d0,0.d0,0.d0,2.d0,0.d0,0.d0,0.d0,0.d0,0.d0,4.d0)
call run('partial rain, cold mean leaf',0.d0,10.d0,2.d0,0.d0,0.d0,0.d0,0.d0,0.d0,4.d0,0.d0)
call run('matching rain release',10.d0,0.d0,2.d0,0.d0,0.d0,0.d0,0.d0,0.d0,4.d0,0.d0)
call run('matching snow release',0.d0,10.d0,0.d0,2.d0,0.d0,0.d0,0.d0,0.d0,0.d0,4.d0)
call run('matching rain with input',10.d0,0.d0,1.d0,0.d0,4.d0,0.d0,2.d0,0.d0,12.d0,0.d0)
call run('simultaneous release',10.d0,0.d0,2.d0,3.d0,0.d0,0.d0,0.d0,0.d0,4.d0,6.d0)
call run('full mixed release',10.d0,0.d0,4.d0,6.d0,0.d0,0.d0,0.d0,0.d0,8.d0,12.d0)
call run('full rain release',0.d0,10.d0,10.d0,0.d0,0.d0,0.d0,0.d0,0.d0,20.d0,0.d0)
call run('full snow release',10.d0,0.d0,0.d0,10.d0,0.d0,0.d0,0.d0,0.d0,0.d0,20.d0)
call run('no release',10.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
! New input has different signatures; xsc must retain OLD signature.
call run('new precipitation and premix release',10.d0,0.d0,1.d0,2.d0,4.d0,3.d0,2.d0,1.d0,12.d0,18.d0)
call run('premix then mixed drip',10.d0,0.d0,1.d0,2.d0,4.d0,3.d0,2.d0,1.d0,12.d0+8.d0/3.d0,21.5d0,1.d0,0.5d0)
call run('nonuniform old phase inventories',8.d0,2.d0,0.d0,3.d0,0.d0,0.d0,0.d0,0.d0,0.d0,6.d0,separate_old=.true.)
call run('tiny snow release',10.d0,0.d0,0.d0,1.d-10,0.d0,0.d0,0.d0,0.d0,0.d0,2.d-10)
DEF_VEG_SNOW=.true.
call run('real phases retain separate signatures',8.d0,2.d0,1.d0,2.d0,0.d0,0.d0,0.d0,0.d0,1.5d0,8.d0)
DEF_VEG_SNOW=.false.
limited=.true.
call run('finite capacity partial snow with solid',10.d0,0.d0,0.d0,2.d0,0.d0,0.d0,0.d0,0.d0,0.d0,4.d0)
call run('finite capacity full mixed with solid',10.d0,0.d0,4.d0,6.d0,0.d0,0.d0,0.d0,0.d0,8.8d0,13.2d0)
if(failures>0) stop 1
print *, 'PASS single-bucket precipitation regression'
contains
subroutine run(label,wr,ws,xr,xs,pr,ps,ir,is,er,es,dr,ds,separate_old)
character(*) :: label
real(r8) :: wr,ws,xr,xs,pr,ps,ir,is,er,es,beg,remaining,actual_r,drip_r,drip_s
real(r8), optional :: dr,ds
logical, optional :: separate_old
logical :: lump
lump=.not.DEF_VEG_SNOW
if(present(separate_old)) lump=.not.separate_old
trc_ldew_rain=12.d0
trc_ldew_snow=8.d0
! Simulate caller lump: old phase inventories differ, then join into host bucket.
if (lump) then
 if(wr==0.d0) then
  trc_ldew_snow=20.d0
  trc_ldew_rain=0.d0
 else
  trc_ldew_rain=20.d0
  trc_ldew_snow=0.d0
 endif
endif
trc_canopy_solid=0.d0
if(limited) trc_canopy_solid=2.d0
trc_surface_solid=0.d0
trc_waterstorage=0.d0
trc_waterstorage_solid=0.d0
a_trc_precip=0.d0
beg=sum(trc_ldew_rain+trc_ldew_snow+trc_canopy_solid)
drip_r=0.d0
drip_s=0.d0
if(present(dr)) drip_r=dr
if(present(ds)) drip_s=ds
remaining=wr+ws-xr-xs+ir+is-drip_r-drip_s
call mt_precip(1,1.d0,pr,ps,ir+is,ir,is,pr-ir+xr+drip_r,ps-is+xs+drip_s,remaining,0.d0,wr,ws, &
 0.d0,ir,is,xr,xs,0.d0,0.d0)
actual_r=sum(trc_pg_rain_ground+trc_surface_solid)
write(*,'(A,6(1X,Z16.16))') trim(label), transfer(trc_ldew_rain(1,1),0_8), &
 transfer(trc_ldew_snow(1,1),0_8), transfer(trc_canopy_solid(1,1),0_8), &
 transfer(trc_pg_rain_ground(1,1),0_8), transfer(trc_pg_snow_ground(1,1),0_8), &
 transfer(trc_surface_solid(1,1),0_8)
call expect(abs(actual_r-er)<1.d-12,label//' rain')
call expect(abs(sum(trc_pg_snow_ground)-es)<1.d-12,label//' snow')
call expect(abs(sum(trc_ldew_rain+trc_ldew_snow+trc_canopy_solid)+actual_r+sum(trc_pg_snow_ground) &
 -beg-pr*5.d0-ps*7.d0)<1.d-12,label//' conservation')
end subroutine
subroutine expect(ok,label)
logical :: ok
character(*) :: label
if(.not.ok) then
 print *, 'FAIL ',label
 failures=failures+1
endif
end subroutine
end program
'''
with tempfile.TemporaryDirectory() as tmp:
    f = pathlib.Path(tmp) / 'check.f90'
    f.write_text(source)
    subprocess.run([os.environ.get('FC', 'gfortran'), '-ffree-line-length-none', '-fcheck=all', '-O0', str(f), '-o', str(f.with_suffix(''))], cwd=tmp, check=True)
    subprocess.run([str(f.with_suffix(''))], check=True)
