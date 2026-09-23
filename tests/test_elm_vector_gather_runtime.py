"""Use the production gather statement with multiple and empty workers."""
from pathlib import Path
import os
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_element_gather_has_valid_receive_buffer_on_every_worker(tmp_path):
    compiler, launcher = shutil.which('mpif90'), shutil.which('mpirun')
    if not compiler or not launcher:
        pytest.skip('MPI Fortran toolchain unavailable')
    if os.environ.get('COLM_RUN_MPI_TESTS') != '1':
        pytest.skip('set COLM_RUN_MPI_TESTS=1 to run real MPI processes')
    source = (ROOT / 'mksrfdata/MOD_ElmVector.F90').read_text()
    start = source.index('\n', source.index('indexelm = landelm%eindex'))
    block = source[start:].split('ENDIF', 1)[1].split('IF (p_iam_worker == p_root) THEN', 1)[0]
    program = '''program probe
implicit none
include 'mpif.h'
integer :: p_iam_worker,p_np_worker,p_err,numelm,p_comm_worker
integer,parameter :: p_root=0
integer,allocatable :: numelm_worker(:)
call MPI_Init(p_err)
p_comm_worker=MPI_COMM_WORLD
call MPI_Comm_rank(p_comm_worker,p_iam_worker,p_err)
call MPI_Comm_size(p_comm_worker,p_np_worker,p_err)
! Root owns no elements, the other workers own distinct positive counts.
numelm=p_iam_worker
''' + block + '''
if (.not. allocated(numelm_worker)) stop 1
if (p_iam_worker==p_root) then
 if(any(numelm_worker/=[(numelm,numelm=0,p_np_worker-1)])) stop 2
endif
deallocate(numelm_worker)
call MPI_Finalize(p_err)
end program
'''
    path = tmp_path / 'gather.f90'
    path.write_text(program)
    exe = tmp_path / 'gather'
    subprocess.run([compiler, '-fcheck=all', str(path), '-o', str(exe)], check=True, capture_output=True)
    env = {**os.environ, 'OMPI_ALLOW_RUN_AS_ROOT': '1', 'OMPI_ALLOW_RUN_AS_ROOT_CONFIRM': '1',
           'OMPI_MCA_rmaps_base_oversubscribe': '1'}
    subprocess.run([launcher, '-np', '3', str(exe)], env=env, check=True, capture_output=True, timeout=60)
