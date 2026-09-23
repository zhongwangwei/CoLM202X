"""The production evaporation limiter permits atmospheric isotope uptake."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_EvapLimit.F90"


def test_zero_pool_can_gain_isotope_from_vapor(tmp_path: Path):
    compiler = require_runnable_fortran_compiler(tmp_path)
    (tmp_path / "precision.f90").write_text(
        "module MOD_Precision\n integer, parameter :: r8=kind(1.d0)\nend module\n"
    )
    (tmp_path / "probe.f90").write_text(
        """program probe
 use MOD_Precision
 use MOD_Tracer_EvapLimit
 implicit none
 real(r8) :: loss
 loss = tracer_atmospheric_tracer_loss(0._r8, 10._r8, 1._r8, &
    298._r8, .false., vapor_ratio, 1.e-30_r8, 0._r8, .false.)
 if (abs(loss + 0.002_r8) > 1.e-14_r8) error stop 'signed uptake lost'
 contains
 real(r8) function vapor_ratio(source_ratio, temp_k, from_ice)
  real(r8), intent(in) :: source_ratio, temp_k
  logical, intent(in) :: from_ice
  vapor_ratio = -0.002_r8
 end function
end program
"""
    )
    exe = tmp_path / "probe"
    subprocess.run(
        [compiler, "-cpp", "-DTRACER", "-I", str(SOURCE.parents[2] / "include"),
         str(tmp_path / "precision.f90"), str(SOURCE), str(tmp_path / "probe.f90"),
         "-o", str(exe)],
        cwd=tmp_path, check=True, capture_output=True, text=True,
    )
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)
