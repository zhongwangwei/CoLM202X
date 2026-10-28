#!/usr/bin/env python3
"""Compile actual step-entry validity, history call, lifecycle and acc1d snippets."""
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[2]
def read(name):
    return (root / 'main' / name).read_text()

hist = read('MOD_Hist.F90')
main, flux, acc = map(read, ['CoLMMAIN.F90', 'MOD_Vars_1DFluxes.F90', 'MOD_Vars_1DAccFluxes.F90'])
def extract(pattern, text, label):
    match = re.search(pattern, text, re.M | re.S)
    assert match, f'production {label} not found'
    return match.group(1)

entry = extract(r'(      is_dry_lake = .*?)(?=      IF \(\(patchtype <= 2\))', main, 'step-entry predicate')
call = extract(r'CALL acc1d \(fsno[^\n]*\n(.*?)\s*CALL acc1d \(sigf', acc, 'history call')
outcall = extract(r'(         ! fraction of saturation area on ground.*?)(?=         ! fraction of veg cover)', hist, 'frcsat output call')
normalization = extract(r'(      IF \(\.not\. is_hist\) RETURN.*?)(?=      select CASE \(HistForm\))', hist, 'actual output normalization')
acc1d = extract(r'(^   SUBROUTINE acc1d \(.*?^   END SUBROUTINE acc1d)', acc, 'acc1d')
def lifecycle(pattern, label, text=flux):
    # Preserve the compile guard from the production block, not a synthetic guard.
    return extract(r'(^#ifdef MaterialTransport\n(?:(?!^#endif).)*?' + pattern +
                   r'(?:(?!^#endif).)*?^#endif)', text, label)

decl = lifecycle(r'logical, allocatable :: frcsat_hist_valid\(:\)', 'flag declaration')
allocate = lifecycle(r'allocate \(frcsat_hist_valid\(numpatch\)\)', 'flag allocation')
deallocate = lifecycle(r'deallocate \(frcsat_hist_valid\)', 'flag deallocation')
count_decl = lifecycle(r'real\(r8\), allocatable :: nac_frcsat', 'count declaration', acc)
count_alloc = lifecycle(r'allocate \(nac_frcsat', 'count allocation', acc)
count_free = lifecycle(r'deallocate \(nac_frcsat', 'count deallocation', acc)
count_reset = extract(r'(            a_frcsat\s*\(:\) = spval.*?)(?=            a_sigf)', acc, 'sum and count reset')
# The production optional filter must be checked before the raw undefined value.
assert acc1d.index('IF (.not. filter(i)) CYCLE') < acc1d.index('IF (var(i) /= spval)')
source = '''module MOD_Precision
integer, parameter :: r8=kind(1.d0)
end module
module MOD_Vars_Global
use MOD_Precision
real(r8), parameter :: spval=-1.d36
end module
module fixture
use MOD_Precision
use MOD_Vars_Global
implicit none
integer, parameter :: numpatch=2
real(r8) :: frcsat(2),a_frcsat(2),nac=0.d0
logical :: p_is_worker=.true.,filter(2)=.true.
integer :: sumarea=0,itime_in_file=1
character(1) :: file_hist='x'
type hist_flags
logical :: frcsat=.true.
end type
type(hist_flags) :: DEF_hist_vars
logical :: enabled=.true.
''' + decl + '\n' + count_decl + '''
contains
logical function mt_any_enabled()
mt_any_enabled=enabled
end function
subroutine alloc()
''' + allocate + '\n' + count_alloc + '\n' + count_reset + '''
end subroutine
subroutine release()
''' + deallocate + '\n' + count_free + '''
end subroutine
subroutine capture(ipatch,patchtype,DEF_USE_Dynamic_Lake,wdsrf,zwt)
integer :: ipatch,patchtype
logical :: DEF_USE_Dynamic_Lake,is_dry_lake
real(r8) :: wdsrf,zwt
''' + entry + '''
end subroutine
subroutine history()
integer :: i
''' + call + '''
end subroutine
subroutine reset()
''' + count_reset + '''
end subroutine
subroutine output()
''' + outcall + '''
end subroutine
subroutine write_history_variable_2d(is_hist,acc_vec,file_hist,varname,itime_in_file,sumarea,filter,longname,units,acc_num,input_mode)
logical :: is_hist,filter(:)
real(r8) :: acc_vec(:)
real(r8), optional :: acc_num(:)
character(*) :: file_hist,varname,longname,units
character(*),optional :: input_mode
integer :: itime_in_file,sumarea
''' + normalization + '''
end subroutine
''' + acc1d + '''
end module
program check
use fixture
use, intrinsic :: ieee_arithmetic
implicit none
integer :: failures=0, cases=0
real(r8) :: nan
#ifdef MaterialTransport
logical, parameter :: compiled_mt=.true.
#else
logical, parameter :: compiled_mt=.false.
#endif
nan=ieee_value(0.d0,ieee_quiet_nan)
call alloc()
#ifdef MaterialTransport
call expect(all(frcsat_hist_valid),'allocation default valid')
#endif
call run('dry stale one',.true.,.true.,4,0.d0,1.d0,0.d0,1.d0,1.d0,spval,spval)
call run('dry initial missing',.true.,.true.,4,0.d0,1.d0,0.d0,1.d0,spval,spval,spval)
call run('invalid NaN never accumulated',.true.,.true.,4,0.d0,1.d0,0.d0,1.d0,nan,2.d0,2.d0)
call run('dry to wet uses entry',.true.,.true.,4,0.d0,1.d0,200.d0,0.d0,1.d0,spval,spval)
call run('wet to dry uses entry',.true.,.true.,4,200.d0,0.d0,0.d0,1.d0,1.d0,spval,1.d0)
call run('100mm and zero zwt valid',.true.,.true.,4,100.d0,0.d0,100.d0,0.d0,.5d0,spval,.5d0)
call run('below100 invalid',.true.,.true.,4,99.d0,0.d0,99.d0,0.d0,.5d0,spval,spval)
call run('positive zwt invalid',.true.,.true.,4,200.d0,1.d-9,200.d0,1.d-9,.5d0,spval,spval)
call run('static lake unchanged',.true.,.false.,4,0.d0,1.d0,0.d0,1.d0,.5d0,spval,.5d0)
call run('soil unchanged',.true.,.true.,0,0.d0,1.d0,0.d0,1.d0,.5d0,spval,.5d0)
call run('runtime materials off unchanged',.false.,.true.,4,0.d0,1.d0,0.d0,1.d0,1.d0,spval,1.d0)
call run('partial valid sum retained',.true.,.true.,4,0.d0,1.d0,0.d0,1.d0,1.d0,2.d0,2.d0)
call run('valid accumulated sum',.true.,.true.,0,0.d0,1.d0,0.d0,1.d0,.25d0,2.d0,2.25d0)
call run('valid signed zero bits',.true.,.true.,0,0.d0,1.d0,0.d0,1.d0,-0.d0,spval,-0.d0)
call mean_check()
call runtimeoff_mean()
call count_checks()
call release()
#ifdef MaterialTransport
call expect(.not.allocated(frcsat_hist_valid),'deallocated')
#endif
call alloc()
#ifdef MaterialTransport
call expect(all(frcsat_hist_valid),'reallocation resets validity')
#endif
call release()
if(failures>0) stop 1
print *, 'PASS production frcsat history scenarios=',cases,' compiled_MT=',compiled_mt
contains
subroutine expect(ok,label)
logical :: ok
character(*) :: label
if(.not.ok) then
 print *, 'FAIL ',label
 failures=failures+1
endif
end subroutine
subroutine runtimeoff_mean()
enabled=.false.
call reset()
frcsat=1.d0
call history()
frcsat=spval
call history()
nac=2.d0
call output()
call expect(all(a_frcsat==.5d0),'FINAL runtimeMT-off original nac')
end subroutine
subroutine count_checks()
#ifdef MaterialTransport
enabled=.true.
call reset()
call expect(all(nac_frcsat==0.d0),'flush count zero')
call expect(all(a_frcsat==spval),'flush sums missing')
nac_frcsat=-1.d0
a_frcsat=2.d0
frcsat=1.d0
call capture(1,0,.false.,0.d0,0.d0)
call capture(2,0,.false.,0.d0,0.d0)
call history()
call expect(all(nac_frcsat==-1.d0),'legacy unknown stays unknown after valid input')
call expect(all(a_frcsat==3.d0),'legacy unknown retains sum until output')
nac=3.d0
call output()
call expect(all(a_frcsat==spval),'unknown count final missing, not unnormalized sum')
call reset()
a_frcsat=4.d0
nac_frcsat=0.d0
call output()
call expect(all(a_frcsat==spval),'zero count final missing even with stale sum')
call reset()
call history()
nac=2.d0
call output()
call expect(all(a_frcsat==1.d0),'flush restores valid-only mean next window')
#endif
end subroutine
subroutine mean_check()
enabled=.true.
call reset()
frcsat=1.d0
call capture(1,4,.true.,200.d0,0.d0)
call capture(2,4,.true.,0.d0,1.d0)
call history()
frcsat=spval
call capture(1,4,.true.,0.d0,1.d0)
call history()
nac=2.d0
call output()
if(compiled_mt) then
 call expect(a_frcsat(1)==1.d0,'FINAL MEAN wet1 plus drymissing must equal1 not0.5')
 call expect(a_frcsat(2)==spval,'FINAL no valid samples missing')
else
 call expect(a_frcsat(1)==.5d0,'FINAL compiledMT-off original nac')
endif
end subroutine
subroutine run(label,on,dynamic,ptype,entryw,entryz,endw,endz,raw,oldsum,want)
character(*) :: label
logical :: on,dynamic
integer :: ptype
real(r8) :: entryw,entryz,endw,endz,raw,oldsum,want,w,z,expected
integer(kind=8) :: before(2)
cases=cases+1
enabled=on
frcsat=[raw,.25d0]
a_frcsat=[oldsum,spval]
#ifdef MaterialTransport
nac_frcsat=0.d0
#endif
before=transfer(frcsat,before)
w=entryw;z=entryz
call capture(1,ptype,dynamic,w,z)
call capture(2,0,.false.,w,z)
! Water-state changes during the step must not reclassify the captured branch.
w=endw;z=endz
call history()
expected=want
if(.not.compiled_mt .or. .not.on) then
 expected=oldsum
 if(raw/=spval) then
  if(oldsum/=spval) then
   expected=oldsum+raw
  else
   expected=raw
  endif
 endif
endif
if(ieee_is_nan(expected)) then
 call expect(ieee_is_nan(a_frcsat(1)),label//' unchanged NaN path')
else
 call expect(transfer(a_frcsat(1),0_8)==transfer(expected,0_8),label//' bits')
endif
call expect(a_frcsat(2)==.25d0,label//' other patch')
call expect(all(transfer(frcsat,before)==before),label//' raw unchanged')
end subroutine
end program
'''
compiler = shutil.which('gfortran')
assert compiler, 'gfortran is required'
with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory)
    (path / 'check.f90').write_text(source)
    for macros in [[], ['-DMaterialTransport']]:
        subprocess.run([compiler, '-cpp', *macros, '-ffree-line-length-none', '-fcheck=all',
                        str(path / 'check.f90'), '-o', str(path / 'check')], cwd=path, check=True)
        subprocess.run([str(path / 'check')], check=True)

# Exercise the actual compatibility/validation routine. Only NetCDF/storage and
# MPI metadata are replaced with deterministic one-block in-memory test doubles.
restart = read('MOD_HistAccRestart.F90')
sidecar = extract(r'(^   SUBROUTINE hacc_frcsat_count \(.*?^   END SUBROUTINE hacc_frcsat_count)',
                  restart, 'sidecar valid-count policy')
# A completed legacy window must not enter the partial-window field reader.
assert re.search(r'IF \(vmin > 0\._r8\) THEN\s+IF \(p_is_worker\) nac = vmin\s+CALL hacc_fields', restart)
assert "CALL hacc_frcsat_count (file, mode)" in restart
sum_field = extract(r"(      CALL hacc_field1d \(file, mode, 'a_frcsat'[^\n]*)", restart, 'actual sum restart field')
reader = extract(r'(^   SUBROUTINE read_marker \(.*?^   END SUBROUTINE read_marker)', restart, 'actual marker vector reader')
io_source = '''module MOD_Precision
integer,parameter::r8=kind(1.d0)
end module
module MOD_Vars_Global
use MOD_Precision
real(r8),parameter::spval=-1.d36
end module
module netcdf
integer,parameter::NF90_NOWRITE=0,NF90_NOERR=0,NF90_ENOTVAR=-1
logical::marker=.false.,field=.false.,checking_orphan=.false.
integer::marker_bad=0
real(8)::version=1.d0,saved(2),saved_sum(2)
contains
integer function nf90_open(file,mode,ncid)
character(*)::file
integer::mode,ncid
ncid=1
nf90_open=0
end function
integer function nf90_inq_varid(ncid,name,varid)
integer::ncid,varid
character(*)::name
nf90_inq_varid=NF90_ENOTVAR
select case(name)
case('hacc_frcsat_count_version')
 if(marker) nf90_inq_varid=NF90_NOERR
case('nac_frcsat')
 if(field) nf90_inq_varid=NF90_NOERR
end select
end function
integer function nf90_close(ncid)
integer::ncid
nf90_close=0
end function
end module
module MOD_Block
contains
subroutine get_filename_block(file,x,y,out)
character(*)::file,out
integer::x,y
out=file
end subroutine
end module
module MOD_LandPatch
integer,parameter::numpatch=3
type patchtype
integer::nblkgrp=1,xblkgrp(1)=[1],yblkgrp(1)=[1]
end type
type(patchtype)::landpatch
end module
module MOD_NetCDFVector
use MOD_LandPatch
use netcdf
use,intrinsic::ieee_arithmetic
contains
subroutine ncio_read_vector(file,name,patch,values)
character(*)::file,name
type(patchtype)::patch
real(8),allocatable::values(:)
allocate(values(3))
values=version
if(marker_bad==1) values(2)=ieee_value(0.d0,ieee_quiet_nan)
if(marker_bad==2) values(2)=ieee_value(0.d0,ieee_positive_inf)
end subroutine
end module
module MOD_MT_Namelist
logical::enabled=.true.
contains
logical function mt_any_enabled()
mt_any_enabled=enabled
end function
end module
module fixture_io
use MOD_Precision
use MOD_Vars_Global
use netcdf
use MOD_LandPatch
implicit none
integer,parameter::HACC_READ=1,HACC_WRITE=2
logical::p_is_io=.true.,p_is_master=.true.,p_is_worker=.true.
real(r8),allocatable::nac_frcsat(:)
real(r8)::nac=2.d0,a_frcsat(2)
contains
subroutine global_max(v)
integer::v(:)
end subroutine
subroutine CoLM_stop(message)
character(*)::message
if(checking_orphan) then
 if(any(nac_frcsat/=[1.d0,0.d0])) stop 98
 if(any(a_frcsat/=[1.d0,spval])) stop 98
 print *, 'ORPHAN_STATE_UNCHANGED'
endif
print *, 'EXPECTED_REJECTION ',message
stop 27
end subroutine
subroutine write_marker(file,name,value)
character(*)::file,name
real(r8)::value
marker=.true.
version=value
end subroutine
''' + reader + '''
subroutine require_var(file,name,patch,dimension,rank,a,b)
character(*)::file,name,dimension
type(patchtype)::patch
integer::rank,a,b
if(.not.marker) call CoLM_stop('missing version marker')
end subroutine
subroutine hacc_field1d(file,mode,name,values,patch,dimension)
character(*)::file,name,dimension
integer::mode
type(patchtype)::patch
real(r8)::values(:)
if(name=='a_frcsat') then
 if(mode==HACC_WRITE) then
  saved_sum=values
 else
  values=saved_sum
 endif
 return
endif
if(mode==HACC_WRITE) then
 saved=values;field=.true.
else
 if(.not.field) call CoLM_stop('missing required count field')
 values=saved
endif
end subroutine
subroutine restart_pair(file,mode)
character(*)::file
integer::mode
''' + sum_field + '''
call hacc_frcsat_count(file,mode)
end subroutine
''' + sidecar + '''
end module
program check_io
use fixture_io
use MOD_MT_Namelist
use,intrinsic::ieee_arithmetic
implicit none
character(32)::mode
allocate(nac_frcsat(2))
call get_command_argument(1,mode)
a_frcsat=[1.d0,spval]
nac_frcsat=[1.d0,0.d0]
call restart_pair('x',HACC_WRITE)
select case(trim(mode))
case('orphan_count')
marker=.false.;field=.true.;checking_orphan=.true.
case('mixed_nan_marker')
marker_bad=1
case('inf_marker')
marker_bad=2
case('missing')
field=.false.
case('version')
version=2.d0
case('negative')
saved=[-2.d0,0.d0]
case('too_many')
saved=[3.d0,0.d0]
case('fractional')
saved=[.5d0,0.d0]
case('nonfinite')
saved=[ieee_value(0.d0,ieee_quiet_nan),0.d0]
case('sum_mismatch')
saved=[0.d0,0.d0]
end select
if(len_trim(mode)>0) then
 call hacc_frcsat_count('x',HACC_READ)
 stop 99
endif
! New partial-window actual registered sum and count restoration.
nac_frcsat=0.d0
a_frcsat=spval
call restart_pair('x',HACC_READ)
if(any(nac_frcsat/=[1.d0,0.d0])) stop 1
if(any(a_frcsat/=[1.d0,spval])) stop 2
! Legacy partial window: no marker => unknown, preserving physical/sum data.
marker=.false.;field=.false.
call hacc_frcsat_count('x',HACC_READ)
if(any(nac_frcsat/=-1.d0)) stop 3
if(any(a_frcsat/=[1.d0,spval])) stop 4
! Unknown survives writing and rereading a newer sidecar.
call hacc_frcsat_count('x',HACC_WRITE)
nac_frcsat=0.d0
call hacc_frcsat_count('x',HACC_READ)
if(any(nac_frcsat/=-1.d0)) stop 5
! Runtime-off leaves legacy state alone.
enabled=.false.
nac_frcsat=0.d0
marker=.false.
call hacc_frcsat_count('x',HACC_READ)
if(any(nac_frcsat/=0.d0)) stop 6
print *, 'PASS production sidecar new/legacy/unknown/runtimeoff'
end program
'''
with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory)
    (path / 'sidecar.f90').write_text(io_source)
    for optimization in ['-O0', '-O2']:
        subprocess.run([compiler, '-cpp', optimization, '-ffree-line-length-none', '-fcheck=all',
                        str(path / 'sidecar.f90'), '-o', str(path / 'sidecar')], cwd=path, check=True)
        subprocess.run([str(path / 'sidecar')], check=True)
        for mode in ['missing', 'version', 'negative', 'too_many', 'fractional', 'nonfinite',
                     'sum_mismatch', 'mixed_nan_marker', 'inf_marker', 'orphan_count']:
            result = subprocess.run([str(path / 'sidecar'), mode], text=True, capture_output=True)
            assert result.returncode == 27 and 'EXPECTED_REJECTION' in result.stdout, (mode, result)
            if mode == 'orphan_count':
                assert 'ORPHAN_STATE_UNCHANGED' in result.stdout, result
            print('PASS production sidecar rejection:', optimization, mode)
