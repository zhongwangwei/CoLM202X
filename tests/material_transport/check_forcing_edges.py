#!/usr/bin/env python3
"""Compile production forcing routines with minimal host stubs; no MPI or data files."""
import pathlib
import re
import subprocess
import sys
import tempfile

root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[2]
forcing = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Forcing.F90').read_text()
namelist = (root / 'main/MATERIAL/MOD_MT_Namelist.F90').read_text()

def routine(source, name):
    pattern = rf'^ *[^\n]*\b(SUBROUTINE|FUNCTION) +{name}\s*\([^\n]*(?:\n.*?)*?^ *END +(SUBROUTINE|FUNCTION) +{name}\b[^\n]*'
    match = re.search(pattern, source, re.I | re.M)
    assert match, name
    return match.group(0)

names = ['tracer_forcing_allocate_state', 'tracer_forcing_deallocate_state',
         'tracer_forcing_prepare_step', 'tracer_forcing_precip_value',
         'tracer_forcing_rain_value', 'tracer_forcing_snow_value',
         'tracer_forcing_phase_fill_at', 'tracer_forcing_valid_value',
         'mt_forcing_lulcc_save', 'mt_forcing_lulcc_remap',
         'tracer_forcing_lulcc_map', 'mt_forcing_lulcc_restore']
state = forcing.split('   real(r8), allocatable :: trc_forc_precip_value', 1)[1].split('   PUBLIC ::', 1)[0]
settings = namelist.split('   IMPLICIT NONE', 1)[1].split('   PUBLIC ::', 1)[0]
source = '''module MOD_SPMD_Task
integer :: stop_count = 0
contains
subroutine CoLM_Stop(message)
character(*), optional :: message
stop_count = stop_count + 1
end subroutine
end module
module fixture
use MOD_SPMD_Task
use, intrinsic :: ieee_arithmetic
implicit none
integer, parameter :: r8 = kind(1.d0), ntracers = 2, N_land_classification = 2, CROPLAND = 2
integer, parameter :: DECODE_INVALID=0, DECODE_NEAR_DRY=1, DECODE_VALID=2
real(r8), parameter :: trc_tiny=1.e-20_r8, trc_forc_max_abs=1.e10_r8
logical :: p_is_worker=.true., trc_runtime_forcing_enabled=.true.
logical :: DEF_USE_PFT=.false., DEF_SOLO_PFT=.false., DEF_FAST_PC=.false.
logical :: isotope(2)=[.true.,.false.], precip_on(2)=.true., vapor_on(2)=.true.
integer :: patchtypes(2)=1
real(r8) :: forc_t(2)=273._r8
type patch_type
real(r8), allocatable :: subfrc(:)
end type
type(patch_type) :: elm_patch
''' + settings + '\nreal(r8), allocatable :: trc_forc_precip_value' + state + '\ncontains\n'
for name in names:
    source += routine(forcing, name) + '\n'
source += routine(namelist, 'check_settings') + '\n' + routine(namelist, 'mt_file_given') + '\n'
source += '''
logical function tracer_is_isotope(i)
integer :: i
tracer_is_isotope=isotope(i)
end function
logical function tracer_forcing_precip_configured(i)
integer :: i
tracer_forcing_precip_configured=precip_on(i)
end function
logical function tracer_forcing_vapor_configured(i)
integer :: i
tracer_forcing_vapor_configured=vapor_on(i)
end function
real(r8) function tracer_precip_default_ratio(i)
integer :: i
tracer_precip_default_ratio=10._r8
end function
real(r8) function tracer_vapor_default_ratio(i)
integer :: i
tracer_vapor_default_ratio=2._r8
end function
real(r8) function tracer_alpha_liq_vap(i,t)
integer :: i
real(r8) :: t
tracer_alpha_liq_vap=1.5_r8
end function
end module
program check_edges
use fixture
implicit none
integer :: i, scheme, failures=0
integer*8 :: new_element(1)
real(r8) :: weights(6)=[0.6_r8,0.49_r8,0._r8,1._r8,0.6_r8,0.6_r8], w, trace(1,0:2), expected
! Existing configurations with no phase fill may set MIN above the new MEMORY default.
DEF_MT_FORCING_MIN_PRECIP=5.e-5_r8
DEF_MT_FORCING_MEMORY_PRECIP=2.e-5_r8
call check_settings()
call expect(stop_count==0, 'namelist accepts unrelated threshold ordering')
vapor_on=.false.
call tracer_forcing_allocate_state(2)
call expect(stop_count==0, 'no phase fill accepts MIN > MEMORY')
vapor_on=.true.
precip_on=.false.
call tracer_forcing_allocate_state(2)
call expect(stop_count==0, 'no precipitation forcing accepts MIN > MEMORY')
precip_on=.true.
isotope=.false.
call tracer_forcing_allocate_state(2)
call expect(stop_count==0, 'no isotopes accepts MIN > MEMORY')
isotope(1)=.true.
call tracer_forcing_allocate_state(0)
call expect(stop_count==1, 'phase fill rejects MIN > MEMORY even on zero-patch rank')
stop_count=0
DEF_MT_FORCING_MIN_PRECIP=5.e-6_r8
call tracer_forcing_allocate_state(2)
call expect(stop_count==0, 'valid phase-fill thresholds')
trc_forc_precip_value=4._r8
trc_forc_precip_step=5._r8
trc_forc_vapor_value=2._r8
trc_forc_precip_status=DECODE_NEAR_DRY
call expect(tracer_forcing_rain_value(1,1)==3._r8, 'near-dry rain uses vapour equilibrium')
call expect(tracer_forcing_snow_value(1,1)==4._r8, 'near-dry snow retains memory')
call expect(tracer_forcing_rain_value(2,1)==4._r8, 'non-isotope fallback unchanged')
trc_forc_precip_status=DECODE_INVALID
call expect(tracer_forcing_rain_value(1,1)==4._r8, 'invalid wet rain retains cache')
trc_forc_has_precip=.true.
call expect(tracer_forcing_rain_value(1,1)==5._r8, 'valid rain uses current record')
call expect(tracer_forcing_snow_value(1,1)==5._r8, 'valid snow uses current record')
! Production save/remap/restore, including forcing-state teardown between save and remap.
allocate(elm_patch%subfrc(2)); elm_patch%subfrc=0.5_r8
do scheme=1,2
 do i=1,size(weights)
  w=weights(i); trace(1,:)=[0._r8,w,1._r8-w]
  elm_patch%subfrc=0.5_r8
  if(scheme==2) elm_patch%subfrc=[w,1._r8-w]
  call tracer_forcing_allocate_state(2)
  trc_forc_precip_value(1,:)=[2._r8,8._r8]
  trc_forc_precip_value(2,:)=[3._r8,9._r8]
  trc_forc_precip_reliable(:,1)=1._r8
  trc_forc_precip_reliable(:,2)=0._r8
  if(i==5) trc_forc_precip_reliable(:,2)=1._r8
  new_element=1_8
  if(i==6) new_element=2_8
  call mt_forcing_lulcc_save()
  call tracer_forcing_deallocate_state()
  if(scheme==1) then
   call mt_forcing_lulcc_remap([1],new_element,[1,2],[1_8,1_8],trace)
  else
   call mt_forcing_lulcc_remap([3],new_element,[1,2],[1_8,1_8])
  endif
  call tracer_forcing_allocate_state(1)
  call mt_forcing_lulcc_restore()
  expected=10._r8
  if(w>0._r8) expected=2._r8
  if(i==5) expected=w*2._r8+(1._r8-w)*8._r8
  if(i==6) expected=10._r8
  call expect(trc_forc_precip_value(1,1)==expected, 'LULCC only reliable sources, weight case')
  call expect(trc_forc_precip_reliable(1,1)==merge(1._r8,0._r8,w>0._r8.and.i/=6), 'LULCC reliable flag')
  expected=w*3._r8+(1._r8-w)*9._r8
  if(i==6) expected=10._r8
  call expect(trc_forc_precip_value(2,1)==expected, 'non-phase-fill area mean unchanged')
 enddo
enddo
call expect(stop_count==0, 'no unexpected CoLM_Stop')
DEF_MT_FORCING_MEMORY_PRECIP=-1._r8
call check_settings()
call expect(stop_count==1, 'negative memory threshold rejected')
if(failures>0) stop 1
print *, 'PASS: forcing fallback, threshold compatibility, LULCC reliable-source remap'
contains
subroutine expect(ok,label)
logical :: ok
character(*) :: label
if(.not.ok) then
 print *, 'FAIL: ',label
 failures=failures+1
endif
end subroutine
end program
'''
with tempfile.TemporaryDirectory(prefix='colm-forcing-check-') as tmp:
    path = pathlib.Path(tmp)
    (path / 'check.F90').write_text(source)
    subprocess.run(['gfortran', '-cpp', '-ffree-line-length-none', '-fcheck=all',
                    '-ffpe-trap=invalid,zero,overflow', '-O0', 'check.F90', '-o', 'check'], cwd=tmp, check=True)
    subprocess.run([str(path / 'check')], cwd=tmp, check=True)
