from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
FORCING_INPUT = ROOT / "main" / "TRACER" / "MOD_Tracer_ForcingInput.F90"


@pytest.fixture(scope="module")
def forcing_input_driver():
    compiler = shutil.which("mpif90") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("Fortran compiler is not available")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        (tmp / "define.h").write_text("#define TRACER\n", encoding="utf-8")
        (tmp / "precision.f90").write_text(
            """
module MOD_Precision
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
end module MOD_Precision
""",
            encoding="utf-8",
        )
        (tmp / "tracer_defs.f90").write_text(
            """
module MOD_Tracer_Defs
  implicit none
  integer :: ntracers = 1
  type :: tracer_info_type
    character(len=32) :: name = 'TEST'
  end type tracer_info_type
  type(tracer_info_type), allocatable :: tracers(:)
  character(len=256) :: parameter_file = ''
contains
  subroutine tracer_param_file_for_index(itrc, aliases, file_param, found)
    integer, intent(in) :: itrc
    character(len=*), intent(in) :: aliases
    character(len=*), intent(out) :: file_param
    logical, intent(out) :: found
    file_param = parameter_file
    found = len_trim(parameter_file) > 0
  end subroutine tracer_param_file_for_index

  function tracer_lower(raw) result(lower)
    character(len=*), intent(in) :: raw
    character(len=len(raw)) :: lower
    integer :: i, code
    lower = raw
    do i = 1, len(raw)
      code = iachar(lower(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    enddo
  end function tracer_lower
end module MOD_Tracer_Defs
""",
            encoding="utf-8",
        )
        (tmp / "spmd.f90").write_text(
            """
module MOD_SPMD_Task
  implicit none
  logical :: p_is_master = .true.
contains
  subroutine CoLM_stop(message)
    character(len=*), intent(in), optional :: message
    if (present(message)) write(*,'(A)') trim(message)
    error stop 1
  end subroutine CoLM_stop
end module MOD_SPMD_Task
""",
            encoding="utf-8",
        )
        (tmp / "driver.f90").write_text(
            """
program forcing_input_driver
  use MOD_Tracer_Defs
  use MOD_Tracer_ForcingInput
  implicit none
  type(tracer_forcing_spec_type) :: spec
  allocate(tracers(1))
  call get_command_argument(1, parameter_file)
  call tracer_forcing_input_load()
  write(*,'(I0)') tracer_forcing_input_count(1)
  if (tracer_forcing_input_count(1) > 0) then
    spec = tracer_forcing_input_get(1, 1)
    write(*,'(I0)') spec%dtime
  endif
end program forcing_input_driver
""",
            encoding="utf-8",
        )

        executable = tmp / "forcing_input_driver"
        subprocess.run(
            [
                compiler,
                "-cpp",
                "-ffree-line-length-0",
                "-I",
                str(tmp),
                str(tmp / "precision.f90"),
                str(tmp / "tracer_defs.f90"),
                str(tmp / "spmd.f90"),
                str(FORCING_INPUT),
                str(tmp / "driver.f90"),
                "-o",
                str(executable),
            ],
            cwd=tmp,
            check=True,
            capture_output=True,
            text=True,
        )
        yield executable, tmp


def run_forcing_driver(forcing_input_driver, text):
    executable, tmp = forcing_input_driver
    parameter_file = tmp / "forcing_parameter.nml"
    parameter_file.write_text(text, encoding="utf-8")
    return subprocess.run(
        [str(executable), str(parameter_file)], cwd=tmp, capture_output=True, text=True
    )


def test_absent_forcing_group_is_optional(forcing_input_driver):
    result = run_forcing_driver(
        forcing_input_driver,
        """
&nl_colm_tracer_parameter
/
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines()[-1] == "0"


def test_malformed_present_forcing_group_fails_fast(forcing_input_driver):
    result = run_forcing_driver(
        forcing_input_driver,
        """
&nl_colm_tracer_forcing
  forcing_num = 1
  unknown_field = 1
/
""",
    )
    assert result.returncode != 0
    assert "invalid &nl_colm_tracer_forcing" in result.stdout
    assert "unknown_field" in result.stdout.lower()


def test_nonpositive_forcing_dtime_fails_fast(forcing_input_driver):
    result = run_forcing_driver(
        forcing_input_driver,
        """
&nl_colm_tracer_forcing
  forcing_num = 1
  forcing_role = 'precip'
  forcing_fprefix = 'test'
  forcing_vname = 'test'
  forcing_dtime = 0
/
""",
    )
    assert result.returncode != 0
    assert "forcing_dtime(1)" in result.stdout


@pytest.mark.parametrize("forcing_num", [-1, 9])
def test_forcing_count_out_of_range_fails_fast(forcing_input_driver, forcing_num):
    result = run_forcing_driver(
        forcing_input_driver,
        f"""
&nl_colm_tracer_forcing
  forcing_num = {forcing_num}
/
""",
    )
    assert result.returncode != 0
    assert "forcing_num" in result.stdout
    assert "must be between 0 and 8" in result.stdout
