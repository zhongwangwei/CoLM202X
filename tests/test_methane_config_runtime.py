from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
SUBPROCESS_TIMEOUT = 60


def test_standard_methane_namelist_passes_the_real_validator(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)

    (tmp_path / "define.h").write_text(
        "#define TRACER\n#define BGC\n", encoding="utf-8"
    )
    (tmp_path / "precision.f90").write_text(
        """
module MOD_Precision
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
end module MOD_Precision
""",
        encoding="utf-8",
    )
    (tmp_path / "spmd.f90").write_text(
        """
module MOD_SPMD_Task
  implicit none
  logical :: p_is_master = .true.
contains
  subroutine CoLM_Stop(message)
    character(len=*), intent(in), optional :: message
    if (present(message)) write(*,'(A)') trim(message)
    error stop 1
  end subroutine CoLM_Stop
end module MOD_SPMD_Task
""",
        encoding="utf-8",
    )
    (tmp_path / "tracer_defs.f90").write_text(
        """
module MOD_Tracer_Defs
  implicit none
contains
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
    (tmp_path / "namelist.f90").write_text(
        """
module MOD_Namelist
  implicit none
  integer :: DEF_wetland_finundation_scheme = 0
  logical :: DEF_USE_Dynamic_Wetland = .true.
end module MOD_Namelist
""",
        encoding="utf-8",
    )
    (tmp_path / "driver.f90").write_text(
        """
program methane_config_driver
  use MOD_Tracer_Reactive_Methane_Const, only: read_methane_namelist
  implicit none
  character(len=1024) :: path
  call get_command_argument(1, path)
  call read_methane_namelist(trim(path))
  write(*,'(A)') 'VALID'
end program methane_config_driver
""",
        encoding="utf-8",
    )

    executable = tmp_path / "methane_config_driver"
    compile_result = subprocess.run(
        [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(tmp_path),
            str(tmp_path / "precision.f90"),
            str(tmp_path / "spmd.f90"),
            str(tmp_path / "tracer_defs.f90"),
            str(tmp_path / "namelist.f90"),
            str(ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90"),
            str(tmp_path / "driver.f90"),
            "-o",
            str(executable),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert compile_result.returncode == 0, compile_result.stdout + compile_result.stderr

    result = subprocess.run(
        [str(executable), str(ROOT / "run/standard_ch4_parameter.nml")],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines()[-1] == "VALID"

    invalid_text = (ROOT / "run/standard_ch4_parameter.nml").read_text(
        encoding="utf-8"
    ).replace(
        "DEF_METHANE%use_microbial_dormancy      = .false.",
        "DEF_METHANE%use_microbial_dormancy      = .true.",
    )
    invalid_file = tmp_path / "invalid_standard.nml"
    invalid_file.write_text(invalid_text, encoding="utf-8")
    invalid = subprocess.run(
        [str(executable), str(invalid_file)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert invalid.returncode != 0
    assert "use_microbial_dormancy requires use_microbial_pools" in invalid.stdout
