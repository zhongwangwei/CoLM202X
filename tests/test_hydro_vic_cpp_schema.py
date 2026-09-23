from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
VIC_VARIABLES = ROOT / "main" / "HYDRO" / "MOD_Hydro_VIC_Variables.F90"


def test_cpp_preserves_both_water_table_curve_components():
    compiler = shutil.which("mpif90") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("Fortran compiler is not available")

    result = subprocess.run(
        [compiler, "-cpp", "-E", "-P", str(VIC_VARIABLES)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "real(r8) :: zwtvmoist_zwt" in result.stdout
    assert "real(r8) :: zwtvmoist_moist" in result.stdout


def test_cpp_generated_soil_type_exposes_zwt_and_moist_components():
    compiler = shutil.which("mpif90") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("Fortran compiler is not available")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        (tmp / "precision.f90").write_text(
            """
module MOD_Precision
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
end module MOD_Precision
""",
            encoding="utf-8",
        )
        (tmp / "vars_global.f90").write_text(
            """
module MOD_Vars_Global
  use MOD_Precision, only: r8
  implicit none
  integer, parameter :: nl_soil = 10
  real(r8) :: dz_soi(nl_soil) = 0.1_r8
  type :: simulation_time_type
    real(r8) :: timestep = 1800._r8
  end type simulation_time_type
  type(simulation_time_type) :: DEF_simulation_time
end module MOD_Vars_Global
""",
            encoding="utf-8",
        )
        (tmp / "consumer.f90").write_text(
            """
program vic_soil_type_consumer
  use MOD_Precision, only: r8
  use MOD_Hydro_VIC_Variables, only: soil_con_struct
  implicit none
  type(soil_con_struct) :: soil
  soil%zwtvmoist_zwt = 0._r8
  soil%zwtvmoist_moist = 0._r8
end program vic_soil_type_consumer
""",
            encoding="utf-8",
        )

        result = subprocess.run(
            [
                compiler,
                "-cpp",
                "-fdefault-real-8",
                "-ffree-line-length-0",
                "-I",
                str(tmp),
                "-J",
                str(tmp),
                str(tmp / "precision.f90"),
                str(tmp / "vars_global.f90"),
                str(VIC_VARIABLES),
                str(tmp / "consumer.f90"),
                "-o",
                str(tmp / "vic_soil_type_consumer"),
            ],
            cwd=tmp,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
