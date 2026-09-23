"""Exercise the production namelist read and MPI broadcast for long tracer paths."""
from pathlib import Path
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'share/MOD_Namelist.F90'


@pytest.fixture(scope='module')
def probe(tmp_path_factory):
    compiler = shutil.which('mpifort')
    launcher = shutil.which('mpiexec')
    if not compiler or not launcher:
        pytest.skip('MPI Fortran compiler/launcher unavailable')
    tmp = tmp_path_factory.mktemp('tracer_param_capacity')
    (tmp / 'define.h').write_text(
        '#define USEMPI\n#define TRACER\n#define GridRiverLakeFlow\n'
        '#define vanGenuchten_Mualem_SOIL_MODEL\n'
    )
    (tmp / 'stub.f90').write_text('''
module MOD_Precision
  integer, parameter :: r8 = selected_real_kind(12)
end module
module MOD_SPMD_Task
  use mpi
  implicit none
  logical :: p_is_master
  integer :: p_address_master=0, p_comm_glb=MPI_COMM_WORLD, p_err
contains
  subroutine CoLM_stop(message)
    character(len=*), optional :: message
    if (present(message)) print *, trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 9, p_err)
  end subroutine
end module
''')
    (tmp / 'driver.f90').write_text('''
program driver
  use mpi
  use MOD_SPMD_Task
  use MOD_Namelist
  implicit none
  integer :: rank, n
  character(len=4096) :: wanted
  character(len=512) :: path, mode
  call MPI_Init(p_err)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, p_err)
  p_is_master = rank == 0
  call get_command_argument(1, path)
  call get_command_argument(2, mode)
  call read_namelist(trim(path))
  if (trim(mode) == 'overflow') then
    print *, 'UNEXPECTED_OVERFLOW_ACCEPTED'
    call MPI_Finalize(p_err)
    stop
  endif
  open(30, file='wanted.txt', status='old')
  read(30, '(A)') wanted
  close(30)
  n = len_trim(DEF_TRACER_PARAM_FILES)
  write(*,'(A,I0,A,I0)') 'RANK=', rank, ' LENGTH=', n
  if (n /= len_trim(wanted)) call MPI_Abort(MPI_COMM_WORLD, 21, p_err)
  if (trim(DEF_TRACER_PARAM_FILES) /= trim(wanted)) call MPI_Abort(MPI_COMM_WORLD, 22, p_err)
  call MPI_Finalize(p_err)
end program
''')
    def compile_one(*args):
        subprocess.run([compiler, '-cpp', '-ffree-line-length-none', '-I', str(tmp), *args],
                       cwd=tmp, check=True, capture_output=True, text=True, timeout=90)
    compile_one('-c', str(tmp / 'stub.f90'))
    compile_one('-c', str(SOURCE))
    compile_one('stub.o', 'MOD_Namelist.o', str(tmp / 'driver.f90'), '-o', 'probe.x')
    return tmp, launcher


def run_case(probe, mapping, mode):
    tmp, launcher = probe
    # The actual production module reads both namelists before broadcasting.
    (tmp / 'forcing.nml').write_text('&nl_colm_forcing /\n')
    (tmp / 'case.nml').write_text(
        "&nl_colm\n"
        f" DEF_dir_output = '{tmp}/output/'\n"
        f" DEF_forcing_namelist = '{tmp}/forcing.nml'\n"
        f" DEF_TRACER_PARAM_FILES = '{mapping}'\n"
        '/\n'
    )
    (tmp / 'wanted.txt').write_text(mapping + '\n')
    return subprocess.run([launcher, '-n', '3', str(tmp / 'probe.x'), 'case.nml', mode],
                          cwd=tmp, capture_output=True, text=True, timeout=90)


def test_six_species_mapping_over_512_bytes_survives_namelist_and_mpi(probe):
    tmp, _ = probe
    # Keep this over 512 bytes even when pytest uses a very short temp root.
    folder = tmp / ('p' * 80)
    folder.mkdir(exist_ok=True)
    names = ('H2_18O', 'HDO', 'CH4', 'SEDIMENT', 'CL', 'FINITE')
    mapping = ','.join(f'{name}:{folder / (name.lower() + ".nml")}' for name in names)
    assert 512 < len(mapping) < 2048
    for name in names:
        (folder / f'{name.lower()}.nml').write_text('&nl_colm_tracer_parameter /\n')
    result = run_case(probe, mapping, 'roundtrip')
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.count(f'LENGTH={len(mapping)}') == 3


def test_namelist_rejects_saturated_parameter_list_instead_of_truncating(probe):
    capacity = int(re.search(r'character\(len=(\d+)\)\s*::\s*DEF_TRACER_PARAM_FILES',
                             SOURCE.read_text(), re.I).group(1))
    mapping = 'CL:' + 'x' * (capacity + 10)
    result = run_case(probe, mapping, 'overflow')
    assert result.returncode != 0
    assert 'DEF_TRACER_PARAM_FILES exceeds' in result.stdout + result.stderr
    assert 'UNEXPECTED_OVERFLOW_ACCEPTED' not in result.stdout
