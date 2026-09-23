"""The single sediment interval rejects non-finite routing inputs without a CFL loop."""

from pathlib import Path
import re
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler

SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


@pytest.mark.parametrize("bad", ["storage", "discharge", "absolute", "timestep"])
def test_nonfinite_interval_input_fails_fast(tmp_path: Path, bad: str) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    calc = re.search(r"SUBROUTINE grid_sediment_calc\b.*?END SUBROUTINE grid_sediment_calc", SOURCE.read_text(), re.S).group()
    guard = calc.split("         IF (.not. ieee_is_finite(dt_morph)", 1)[1].split("         call_dt_min =", 1)[0]
    guard = "         IF (.not. ieee_is_finite(dt_morph)" + guard
    code = """
module probe
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer, parameter :: r8=kind(1.d0)
 real(r8) :: dt_morph=3600._r8, rivsto_donor(1)=[1._r8], &
  rivsto(1)=[1._r8], rivout(1)=[1._r8], rivout_abs(1)=[1._r8]
 contains
 subroutine CoLM_stop(message)
  character(len=*), intent(in) :: message
  print *, 'REJECTED_INTERVAL_INPUT: ', message
  error stop 1
 end subroutine
 subroutine check()
""" + guard + """
 end subroutine
end module
program main
 use probe
 select case ('BAD')
 case ('storage'); rivsto(1)=ieee_value(0._r8,ieee_quiet_nan)
 case ('discharge'); rivout(1)=ieee_value(0._r8,ieee_positive_inf)
 case ('absolute'); rivout_abs(1)=ieee_value(0._r8,ieee_quiet_nan)
 case ('timestep'); dt_morph=ieee_value(0._r8,ieee_quiet_nan)
 end select
 call check()
 error stop 2
end program
""".replace("BAD", bad)
    path = tmp_path / "probe.f90"
    path.write_text(code)
    built = subprocess.run([compiler, "-ffree-line-length-0", str(path), "-o", "probe"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run(["./probe"], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode != 0 and "REJECTED_INTERVAL_INPUT" in ran.stdout, ran.stdout + ran.stderr
