#!/usr/bin/env python3
"""F01 regression: main-step runoff hand-off to the river routing.

Production routines under test (extracted verbatim from
``main/MATERIAL/MOD_MT_Vars.F90`` with a minimal stub host):
``mt_add_rnof_step`` and ``mt_reset_rnof_step``.

Invariants
----------
* every lake/waterbody/glacier substep ADDS its tracer runoff mass to the
  main-step hand-off buffer (no per-substep overwrite or clear);
* reset clears the generic-water tracers once per main step only; provider-
  owned materials are untouched;
* with constant runoff the received mass equals the land history sum for
  1, 2 and 3 substeps; with variable runoff and with a zero last substep the
  received mass still equals the exact substep sum;
* negative/out-of-range input is ignored.

The test also asserts the structural wiring in the production sources:
the per-substep clears in MOD_MT_SpecialPatches / MOD_MT_Conservation are
gone, the special patches call ``mt_add_rnof_step``, and CoLM.F90 resets the
buffer once before CoLMDRIVER.
"""
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

root = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    pathlib.Path(__file__).resolve().parents[2]
vars_src = (root / 'main/MATERIAL/MOD_MT_Vars.F90').read_text()
spec_src = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_SpecialPatches.F90').read_text()
cons_src = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Conservation.F90').read_text()
colm_src = (root / 'main/CoLM.F90').read_text()


def routine(source, name):
    pattern = (r'^ *[^\n]*\b(SUBROUTINE|FUNCTION) +' + name +
               r'\s*\([^\n]*(?:\n.*?)*?^ *END +(SUBROUTINE|FUNCTION) +' + name + r'\b[^\n]*')
    match = re.search(pattern, source, re.I | re.M)
    assert match, f'production routine not found: {name}'
    return match.group(0)


fixture = r'''
module fixture
use, intrinsic :: iso_fortran_env, only: real64
implicit none
integer, parameter :: r8 = real64
integer, parameter :: ntracers = 3
integer, parameter :: numpatch = 2
real(r8), allocatable :: trc_rnof_step(:,:)
logical :: generic(ntracers) = [.true., .true., .false.]
contains
logical function mt_uses_generic_water(itrc)
  integer, intent(in) :: itrc
  mt_uses_generic_water = generic(itrc)
end function
'''

program = r'''
end module fixture

program check_rnof_handoff
use fixture
implicit none
integer :: nfail, npass, k
real(r8), parameter :: tol = 1.e-15_r8
real(r8) :: land_sum, river_recv, masses(3)
logical :: ok

nfail = 0
npass = 0
allocate(trc_rnof_step(ntracers, numpatch))
trc_rnof_step = 0._r8

! ---- 1 substep, constant runoff R = .002 * 1 mm (audit example) ----
masses = [0.00072_r8, 0._r8, 0._r8]
call run_step(masses, 1, land_sum, river_recv)
call report('F01 1 substep: river received == land accumulated', &
   abs(river_recv - land_sum) < tol, nfail, npass)

! ---- 2 substeps, constant runoff ----
masses = [0.00072_r8, 0.00072_r8, 0._r8]
call run_step(masses, 2, land_sum, river_recv)
call report('F01 2 substeps constant: river received == land accumulated', &
   abs(river_recv - land_sum) < tol, nfail, npass)

! ---- 3 substeps, constant runoff ----
masses = [0.00072_r8, 0.00072_r8, 0.00072_r8]
call run_step(masses, 3, land_sum, river_recv)
call report('F01 3 substeps constant: river received == land accumulated', &
   abs(river_recv - land_sum) < tol, nfail, npass)

! ---- 3 substeps, variable runoff, last substep zero ----
masses = [0.00091_r8, 0.00037_r8, 0.0_r8]
call run_step(masses, 3, land_sum, river_recv)
call report('F01 variable runoff + zero last substep: river == land', &
   abs(river_recv - land_sum) < tol, nfail, npass)
call report('F01 zero last substep keeps earlier substep mass', &
   river_recv > 0.0012_r8 .and. abs((river_recv - (0.00091_r8 + 0.00037_r8))) < tol, nfail, npass)

! ---- provider-owned material (itrc=3) is never touched ----
call mt_reset_rnof_step()
call mt_add_rnof_step(3, 1, 5.0_r8)
call report('F01 provider-owned tracer is not accumulated', &
   trc_rnof_step(3, 1) == 0._r8, nfail, npass)

! ---- negative and out-of-range input ignored ----
call mt_reset_rnof_step()
call mt_add_rnof_step(1, 1, -0.5_r8)
call mt_add_rnof_step(1, 3, 0.5_r8)
call mt_add_rnof_step(0, 1, 0.5_r8)
call report('F01 negative/out-of-range additions ignored', &
   trc_rnof_step(1, 1) == 0._r8, nfail, npass)

! ---- reset zeroes generic tracers only ----
trc_rnof_step(1, 1) = 1._r8
trc_rnof_step(3, 1) = 1._r8
call mt_reset_rnof_step()
call report('F01 reset zeroes generic-water tracers', &
   trc_rnof_step(1, 1) == 0._r8 .and. all(trc_rnof_step(2, :) == 0._r8), nfail, npass)
call report('F01 reset leaves provider-owned tracer alone', &
   trc_rnof_step(3, 1) == 1._r8, nfail, npass)

write(*,'(A,I0,A,I0,A)') 'PASS: ', npass, ' checks; FAIL: ', nfail, ' checks'
if (nfail > 0) stop 1

contains

subroutine run_step(masses, nsub, land_sum, river_recv)
real(r8), intent(in) :: masses(:)
integer, intent(in) :: nsub
real(r8), intent(out) :: land_sum, river_recv
integer :: s
land_sum = 0._r8
call mt_reset_rnof_step()
do s = 1, nsub
   ! land history accumulates every substep; hand-off buffer likewise
   land_sum = land_sum + masses(s)
   call mt_add_rnof_step(1, 1, masses(s))
end do
river_recv = trc_rnof_step(1, 1)
end subroutine

subroutine report(label, ok, nfail, npass)
character(len=*), intent(in) :: label
logical, intent(in) :: ok
integer, intent(inout) :: nfail, npass
if (ok) then
   write(*,'(A,A)') 'PASS: ', label
   npass = npass + 1
else
   write(*,'(A,A)') 'FAIL: ', label
   nfail = nfail + 1
end if
end subroutine

end program check_rnof_handoff
'''

# Structural checks on the production sources.
struct_ok = True
def check(label, cond):
    global struct_ok
    print(('PASS: ' if cond else 'FAIL: ') + label)
    struct_ok = struct_ok and cond

check('F01 SpecialPatches no longer clears trc_rnof_step per substep',
      'trc_rnof_step(itrc, ipatch) = 0._r8' not in spec_src)
check('F01 SpecialPatches accumulates through mt_add_rnof_step',
      spec_src.count('CALL mt_add_rnof_step') == 2)
check('F01 Conservation no longer clears trc_rnof_step',
      'trc_rnof_step(itrc, ipatch) = 0._r8' not in cons_src)
check('F01 CoLM.F90 resets the buffer once before CoLMDRIVER',
      'CALL mt_reset_rnof_step ()' in colm_src and
      colm_src.index('CALL mt_reset_rnof_step ()') < colm_src.index('CALL CoLMDRIVER'))
check('F01 MOD_MT_Vars helper is public',
      'PUBLIC :: mt_add_rnof_step, mt_reset_rnof_step' in vars_src)

gfortran = shutil.which('gfortran')
if not gfortran:
    print('SKIP: gfortran not available; this check needs a Fortran compiler')
    sys.exit(0 if struct_ok else 1)

with tempfile.TemporaryDirectory() as tmp:
    tmp = pathlib.Path(tmp)
    src = tmp / 'check_rnof_handoff.F90'
    exe = tmp / 'check_rnof_handoff'
    src.write_text(fixture +
                   routine(vars_src, 'mt_add_rnof_step') + '\n' +
                   routine(vars_src, 'mt_reset_rnof_step') + '\n' +
                   program.lstrip('\n'))
    proc = subprocess.run([gfortran, '-std=f2008', '-O0', '-ffree-line-length-none',
                           '-fcheck=bounds', '-o', str(exe), str(src)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        print('FAIL: compile of extracted production routines failed')
        print(proc.stdout)
        print(proc.stderr)
        sys.exit(1)
    run = subprocess.run([str(exe)], capture_output=True, text=True)
    print(run.stdout, end='')
    if run.returncode != 0:
        print(run.stderr, end='')
        sys.exit(1)
sys.exit(0 if struct_ok else 1)
