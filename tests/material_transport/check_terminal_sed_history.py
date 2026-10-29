#!/usr/bin/env python3
"""Directed audit reproduction, not a full CoLM/NetCDF integration test.

Compile the unchanged production sediment accumulator, output and flush routines,
and the history payload sections of the production restart reader/writer. Only
the IO/network dependencies are stubbed. Compare continuous, normal mid-window
restart, and terminal partial-window restart call orders. The terminal order
comes from CoLM -> MOD_Hist -> MOD_Grid_RiverLakeHist -> WRITE_TimeVariables.
The source audit below confirms the pre-flush sidecar has no sediment payload.

Usage: python3 repro_terminal_sed_history.py REPO_ROOT [gfortran]
REPO_ROOT defaults to the current directory; the compiler defaults to FC, then
gfortran from PATH. Generated Fortran, executable and module files are confined
to a temporary directory that is removed on exit. Results are printed to stdout.
"""
from pathlib import Path
import argparse
import os
import re
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('repo_root', nargs='?', default='.', metavar='REPO_ROOT',
                    help='CoLM source checkout (default: current directory)')
parser.add_argument('compiler', nargs='?', default=os.environ.get('FC', 'gfortran'),
                    metavar='gfortran', help='GNU Fortran executable (default: FC or gfortran)')
args = parser.parse_args()
ROOT = Path(args.repo_root).expanduser().resolve()
fc = shutil.which(os.path.expanduser(args.compiler))
if fc is None:
    parser.error(f'compiler not found: {args.compiler}')
SED = ROOT / 'main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90'
if not SED.is_file():
    parser.error(f'not a CoLM source checkout: {ROOT}')
text = SED.read_text()

def routine(name):
    match = re.search(r'^   SUBROUTINE ' + name + r'\b.*?^   END SUBROUTINE ' + name + r'\b',
                      text, flags=re.M | re.S | re.I)
    assert match, name
    return match.group(0)

writer = routine('write_sediment_restart')
start = writer.index('      IF (sedhistsnap_valid .and. p_is_worker')
end = writer.index('      CALL write_sediment_scalar_meta', start)
writer = ('real(r8), allocatable :: t2a(:,:), t2b(:,:), t2c(:,:), t2d(:,:), t2e(:,:), t2f(:,:), t1a(:)\n'
          'real(r8) :: t_acctime\n') + writer[start:end]
snapshot = routine('mt_sediment_snapshot_hist_acc')
reader = routine('read_sediment_restart')
start = reader.index('      DO ised = 1, nsed', reader.index("'sed_precip_time_vec'"))
end = reader.index('      IF (.not. has_schema)', start)
reader = reader[start:end]

sidecar = (ROOT / 'main/MOD_HistAccRestart.F90').read_text()
assert not re.search(r'CALL\s+.*(?:sediment|mt_riverlake_|lifecycle)', sidecar, re.I)
assert not re.search(r'CALL\s+mt_riverlake_(?:write|read)_hist_acc_restart',
                     '\n'.join(p.read_text(errors='replace') for p in (ROOT / 'main').rglob('*.F90')), re.I)

code = '''module MOD_Precision
  implicit none
  integer, parameter :: r8 = selected_real_kind(15)
end module
module MOD_Grid_RiverLakeNetwork
  implicit none
  integer, parameter :: numucat=1, totalnumucat=1
  integer :: ucat_data_address(1)=1
end module
module MOD_Grid_RiverLakeHistRoute
  use MOD_Precision
  implicit none
  real(r8) :: reported_conc=0, reported_out=0
contains
  subroutine route_hist_write_ucat(values, name, longname, units)
    real(r8), intent(in) :: values(:)
    character(*), intent(in) :: name, longname, units
    if (name == 'f_sedcon_1') reported_conc=values(1)
    if (name == 'f_sedout_1') reported_out=values(1)
  end subroutine
end module
module audit_memory_io
  use MOD_Precision
  implicit none
  character(64) :: names(8)=''
  real(r8) :: values(8)=0
  integer :: nwritten=0
contains
  subroutine vector_gather_and_write(data, nlocal, nglobal, address, file, name, dim)
    real(r8), intent(in) :: data(:)
    integer, intent(in) :: nlocal, nglobal, address(:)
    character(*), intent(in) :: file, name, dim
    integer :: k
    do k=1,nwritten
      if (names(k)==name) exit
    end do
    if (k>nwritten) nwritten=k
    if (k>8) error stop 'test IO field capacity'
    names(k)=name
    values(k)=data(1)
  end subroutine
  subroutine require_sediment_restart_var(file, ncid, name, data, nlocal, address)
    character(*), intent(in) :: file, name
    integer, intent(in) :: ncid, nlocal, address(:)
    real(r8), allocatable, intent(out) :: data(:)
    integer :: k
    do k=1,nwritten
      if (names(k)==name) exit
    end do
    if (k>nwritten) error stop 'missing test restart field'
    allocate(data(nlocal))
    data=values(k)
  end subroutine
end module
module audit_sediment
  use MOD_Precision
  use MOD_Grid_RiverLakeNetwork
  use MOD_Grid_RiverLakeHistRoute
  use audit_memory_io
  implicit none
  integer, parameter :: nsed=1
  logical :: p_is_worker=.true.
  type hist_flags
    logical :: sedcon=.true., sedout=.true., bedout=.true., sedinp=.true.
    logical :: netflw=.true., sedlayer=.true., shearvel=.true.
  end type
  type(hist_flags) :: DEF_hist_vars
  real(r8) :: sed_hist_acctime=0
  real(r8) :: sedcon(1,1),sedout(1,1),bedout(1,1),sedinp(1,1)
  real(r8) :: netflw(1,1),netflw_adv_step(1,1),layer(1,1),shearvel(1)
  real(r8), allocatable :: a_sedcon(:,:),a_sedout(:,:),a_bedout(:,:),a_sedinp(:,:)
  real(r8), allocatable :: a_netflw(:,:),a_layer(:,:),a_shearvel(:)
contains
  logical function sediment_particle_enabled()
    sediment_particle_enabled=.true.
  end function
  subroutine setup()
    allocate(a_sedcon(1,1),a_sedout(1,1),a_bedout(1,1),a_sedinp(1,1))
    allocate(a_netflw(1,1),a_layer(1,1),a_shearvel(1))
    call flush_sediment_history()
  end subroutine
  subroutine half_day(conc, flux)
    real(r8), intent(in) :: conc,flux
    sedcon=conc; sedout=flux; bedout=0; sedinp=0
    netflw=0; netflw_adv_step=0; layer=1; shearvel=0
    call accumulate_sediment_output(43200._r8)
    sed_hist_acctime=sed_hist_acctime+43200._r8
  end subroutine
  subroutine history_payload_write()
    character(*), parameter :: file_restart='test'
    integer :: ised
    character(16) :: cised
    real(r8) :: dummy_sed(1)
    real(r8), allocatable :: scalar_vec(:)
''' + writer + '''
  end subroutine
  subroutine history_payload_read()
    character(*), parameter :: file_restart='test'
    integer :: ised, ncid=0
    character(16) :: cised
    character(64) :: vname
    real(r8), allocatable :: buf(:)
''' + reader + '''
  end subroutine
''' + routine('accumulate_sediment_output') + '\n' + routine('write_sediment_history') + '\n' + routine('flush_sediment_history') + '''
end module
program reproduce
  use audit_sediment
  implicit none
  real(r8) :: full_conc, full_flux, split_conc, split_flux, ordinary_conc, ordinary_flux
  call setup()
  call half_day(0.001_r8,1._r8)
  call half_day(0.003_r8,3._r8)
  call write_sediment_history('test',1)
  full_conc=reported_conc; full_flux=reported_out
  print '(a,es16.8,a,f8.4)', 'continuous day concentration=',full_conc,' flux=',full_flux

  ! Ordinary checkpoint does not trigger terminal partial history output.
  call flush_sediment_history()
  call half_day(0.001_r8,1._r8)
  call history_payload_write()
  call flush_sediment_history()
  call history_payload_read()
  call half_day(0.003_r8,3._r8)
  call write_sediment_history('test',1)
  ordinary_conc=reported_conc; ordinary_flux=reported_out
  print '(a,es16.8,a,f8.4)', 'regular checkpoint concentration=',ordinary_conc,' flux=',ordinary_flux

  ! Terminal partial window: the sidecar omits sediment fields (source audited).
  ! MOD_Hist writes the partial window and the routing history flushes it.
  call flush_sediment_history()
  call half_day(0.001_r8,1._r8)
  call mt_sediment_snapshot_hist_acc()
  call write_sediment_history('test',1)
  call flush_sediment_history()
  ! WRITE_TimeVariables follows hist_out, serializing already-cleared fields.
  call history_payload_write()
  call history_payload_read()
  print '(a,f10.1)', 'terminal checkpoint restored sediment seconds=',sed_hist_acctime
  call half_day(0.003_r8,3._r8)
  call write_sediment_history('test',2)
  split_conc=reported_conc; split_flux=reported_out
  print '(a,es16.8,a,f8.4)', 'terminal split day concentration=',split_conc,' flux=',split_flux
  if (abs(full_conc-0.002_r8)>1.e-14_r8) error stop 'unexpected continuous average'
  if (ordinary_conc/=full_conc .or. ordinary_flux/=full_flux) error stop 'regular restart failed'
  if (abs(split_conc-full_conc)>1.e-12_r8 .or. abs(split_flux-full_flux)>1.e-8_r8) &
   error stop 'N05 terminal split does not match the continuous window'
  print *, 'N05 PASS: terminal partial-window restart matches continuous history window'
  if (abs(split_conc-full_conc)>1.e-12_r8 .or. abs(split_flux-full_flux)>1.e-8_r8) error stop 'N05 split mismatch'
  print *, 'N05 PASS: terminal partial-window restart matches continuous history window'
end program
'''
# N05 fix: make the terminal snapshot + pre-flush rows part of the audited module
code = code.replace(
 'real(r8), allocatable :: a_netflw(:,:),a_layer(:,:),a_shearvel(:)',
 'real(r8), allocatable :: a_netflw(:,:),a_layer(:,:),a_shearvel(:)\n'
 '  real(r8), allocatable, save :: sedhistsnap_con(:,:), sedhistsnap_out(:,:), sedhistsnap_bed(:,:)\n'
 '  real(r8), allocatable, save :: sedhistsnap_inp(:,:), sedhistsnap_net(:,:), sedhistsnap_lay(:,:)\n'
 '  real(r8), allocatable, save :: sedhistsnap_shear(:)\n'
 '  real(r8), save :: sedhistsnap_acctime=0\n'
 '  logical, save :: sedhistsnap_valid=.false.')
code = code.replace('end module' + chr(10) + 'program reproduce', snapshot + chr(10) + 'end module' + chr(10) + 'program reproduce')

with tempfile.TemporaryDirectory(prefix='colm-terminal-sed-history-') as temp:
    build = Path(temp)
    src = build / 'repro_terminal_sed_history.f90'
    exe = build / 'repro_terminal_sed_history'
    src.write_text(code)
    subprocess.run([fc, '-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                    '-ffree-line-length-none', '-J', str(build), str(src), '-o', str(exe)],
                   check=True, cwd=build)
    result = subprocess.run([str(exe)], text=True, capture_output=True, check=False, cwd=build)
    print(result.stdout, end='')
    if result.returncode != 0: raise SystemExit(result.returncode)
    print(result.stdout, end='')
    if result.stderr:
        print(result.stderr, end='')
