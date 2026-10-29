#!/usr/bin/env python3
"""N06 real-MPI regression for the production point-to-point push protocol.

Compiles the *extracted production* `worker_push_data_uniq_real8` from
share/MOD_WorkerPushData.F90 (USEMPI branch) with a minimal fixture that
provides only the pushdata type / module globals, then runs it with a real
2-rank OpenMPI:

* both ranks exchange a non-zero value -> each receives the peer value;
* both ranks exchange a *zero* value (the "local inventory is empty" case)
  -> the calls still match and complete, i.e. the unconditional fixed-order
  discipline used by the river forwarding cannot deadlock on zero data.

A data-dependent skip on one rank is *not* run as a test because it hangs by
design; the exact structural rule is asserted in check_audit_p1_fixes.py.
"""
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

root = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    pathlib.Path(__file__).resolve().parents[2]
src = (root / 'share/MOD_WorkerPushData.F90').read_text()

m = re.search(r'^\s*SUBROUTINE worker_push_data_uniq_real8\b.*?^\s*END SUBROUTINE '
              r'worker_push_data_uniq_real8\b[^\n]*', src, re.I | re.M | re.S)
if not m:
    print('FAIL: production worker_push_data_uniq_real8 not found')
    sys.exit(1)
routine = m.group(0)

mpirun = shutil.which('mpirun')
mpif90 = shutil.which('mpif90')
if not mpirun or not mpif90:
    print('SKIP: mpirun/mpif90 not available; real MPI check not run')
    sys.exit(0)

fixture = r'''
module fixture
use, intrinsic :: iso_fortran_env, only: real64
implicit none
integer, parameter :: r8 = real64
type :: pointer_int32_1d
   integer, allocatable :: val(:)
end type
type :: worker_pushdata_type
   integer :: num_req_uniq = 0
   integer :: nself = 0
   integer, allocatable :: self_from(:), self_to(:)
   integer, allocatable :: n_to_other(:), n_from_other(:)
   type(pointer_int32_1d), allocatable :: to_other(:), other_to(:)
end type
logical :: p_is_worker = .true.
integer :: p_np_worker = 2
integer :: p_comm_worker
integer :: p_err
end module
'''

program = r'''
program mpi_push_check
use mpi
use fixture
use pushmod, only: worker_push_data_uniq_real8
implicit none
type(worker_pushdata_type) :: push
integer :: rank, nprocs, other, nfail, ierr
real(r8) :: sendv(1), recv(1), expected

call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
if (nprocs /= 2) then
   write(*,'(A)') 'FAIL: this check needs exactly 2 ranks'
   call mpi_abort(MPI_COMM_WORLD, 1, ierr)
endif
p_comm_worker = MPI_COMM_WORLD
p_np_worker = nprocs
other = 1 - rank

push%num_req_uniq = 1
push%nself = 0
allocate(push%n_to_other(0:nprocs-1), push%n_from_other(0:nprocs-1))
allocate(push%to_other(0:nprocs-1), push%other_to(0:nprocs-1))
push%n_to_other = 0
push%n_from_other = 0
push%n_to_other(other) = 1
push%n_from_other(other) = 1
allocate(push%to_other(other)%val(1));   push%to_other(other)%val(1) = 1
allocate(push%other_to(other)%val(1));   push%other_to(other)%val(1) = 1

nfail = 0
! non-zero exchange
sendv(1) = real(rank + 1, r8)
recv = -1._r8
call worker_push_data_uniq_real8(push, sendv, recv, 0._r8)
expected = real(other + 1, r8)
if (abs(recv(1) - expected) > 1e-14_r8) then
   write(*,'(A,I0,A,ES14.6,A,ES14.6)') 'FAIL rank ', rank, ': recv=', recv(1), ' expected=', expected
   nfail = nfail + 1
endif

! zero-data exchange (local inventory empty) must still complete on both ranks
sendv(1) = 0._r8
recv = -1._r8
call worker_push_data_uniq_real8(push, sendv, recv, 0._r8)
if (abs(recv(1)) > 1e-14_r8) then
   write(*,'(A,I0,A,ES14.6)') 'FAIL rank ', rank, ': zero exchange recv=', recv(1)
   nfail = nfail + 1
endif

if (nfail == 0) write(*,'(A,I0,A)') 'PASS rank ', rank, ': fixed-order push matched with non-zero and zero data'
call mpi_finalize(ierr)
if (nfail /= 0) error stop 'N06 MPI push protocol check failed'
end program
'''

with tempfile.TemporaryDirectory() as td:
    td = pathlib.Path(td)
    (td / 'define.h').write_text('#define USEMPI\n')
    f90 = td / 'mpi_push_check.F90'
    wrapper = 'module pushmod\nuse fixture\nuse mpi\nimplicit none\ncontains\n'
    f90.write_text(fixture + wrapper + routine + '\nend module pushmod\n' + program)
    exe = td / 'mpi_push_check'
    comp = subprocess.run([mpif90, '-cpp', '-DUSEMPI', '-I' + str(td), '-ffree-line-length-none',
                           '-fcheck=all', '-O0', '-o', str(exe), str(f90)],
                          capture_output=True, text=True)
    if comp.returncode != 0:
        print('FAIL: compile of extracted production routine failed')
        print(comp.stdout)
        print(comp.stderr)
        sys.exit(1)
    try:
        run = subprocess.run([mpirun, '-np', '2', str(exe)], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        print('FAIL: real 2-rank MPI run hung (unmatched push)')
        sys.exit(1)
    print(run.stdout, end='')
    if run.returncode != 0:
        print(run.stderr)
        sys.exit(1)
    if run.stdout.count('PASS rank') != 2:
        print('FAIL: expected PASS from both ranks')
        sys.exit(1)
print('PASS: N06 real 2-rank MPI zero-data push exchange completed')
