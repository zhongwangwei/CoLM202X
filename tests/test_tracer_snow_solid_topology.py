from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
SNOW_TOPOLOGY = ROOT / "main" / "MOD_SnowLayersCombineDivide.F90"


def test_snow_layer_combine_and_divide_conserve_solid_tracer_inventory():
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
        (tmp / "physical.f90").write_text(
            """
module MOD_Const_Physical
  use MOD_Precision, only: r8
  implicit none
  real(r8), parameter :: denice = 917.0_r8
  real(r8), parameter :: denh2o = 1000.0_r8
  real(r8), parameter :: tfrz = 273.15_r8
  real(r8), parameter :: cpice = 2117.27_r8
  real(r8), parameter :: cpliq = 4188.0_r8
  real(r8), parameter :: hfus = 333600.0_r8
end module MOD_Const_Physical
""",
            encoding="utf-8",
        )
        (tmp / "driver.f90").write_text(
            """
program snow_solid_topology_driver
  use MOD_Precision, only: r8
  use MOD_SnowLayersCombineDivide, only: snowlayerscombine, snowlayersdivide
  implicit none
  integer, parameter :: lb = -4
  integer :: snl
  real(r8) :: z(lb:1), dz(lb:1), zi(lb-1:1)
  real(r8) :: wliq(lb:1), wice(lb:1), temp(lb:1)
  real(r8) :: trc_wliq(1,lb:1), trc_wice(1,lb:1), trc_solid(1,lb:1)
  real(r8) :: trc_scv(1), scv, snowdp

  z = 0.0_r8; dz = 0.0_r8; zi = 0.0_r8
  wliq = 0.0_r8; wice = 0.0_r8; temp = 273.15_r8
  trc_wliq = 0.0_r8; trc_wice = 0.0_r8; trc_solid = 0.0_r8
  trc_scv = 0.0_r8
  snl = -2
  dz(-1) = 0.02_r8; dz(0) = 0.02_r8
  wice(-1) = 0.05_r8; wice(0) = 1.0_r8
  trc_solid(1,-1) = 0.01_r8; trc_solid(1,0) = 0.02_r8
  scv = sum(wliq(-1:0) + wice(-1:0))
  snowdp = sum(dz(-1:0))
  call snowlayerscombine(lb, snl, z, dz, zi, wliq, wice, temp, scv, snowdp, &
    trc_wliq=trc_wliq, trc_wice=trc_wice, trc_solid=trc_solid, trc_scv=trc_scv)
  write(*,'(A,I0,3(1X,ES24.16))') 'COMBINE=', snl, sum(trc_solid), &
    trc_solid(1,0), trc_solid(1,1)

  z = 0.0_r8; dz = 0.0_r8; zi = 0.0_r8
  wliq = 0.0_r8; wice = 0.0_r8; temp = 273.15_r8
  trc_wliq = 0.0_r8; trc_wice = 0.0_r8; trc_solid = 0.0_r8
  snl = -1
  dz(0) = 0.04_r8; wliq(0) = 0.2_r8; wice(0) = 1.0_r8
  trc_solid(1,0) = 0.02_r8
  call snowlayersdivide(lb, snl, z, dz, zi, wliq, wice, temp, &
    trc_wliq=trc_wliq, trc_wice=trc_wice, trc_solid=trc_solid)
  write(*,'(A,I0,3(1X,ES24.16))') 'DIVIDE=', snl, sum(trc_solid), &
    trc_solid(1,-1), trc_solid(1,0)
end program snow_solid_topology_driver
""",
            encoding="utf-8",
        )

        executable = tmp / "snow_solid_topology_driver"
        command = [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(tmp),
            str(tmp / "precision.f90"),
            str(tmp / "physical.f90"),
            str(SNOW_TOPOLOGY),
            str(tmp / "driver.f90"),
            "-o",
            str(executable),
        ]
        subprocess.run(command, cwd=tmp, check=True, capture_output=True, text=True)
        result = subprocess.run(
            [str(executable)], cwd=tmp, check=True, capture_output=True, text=True
        )

    lines = result.stdout.splitlines()
    combine = next(line for line in lines if line.startswith("COMBINE="))
    divide = next(line for line in lines if line.startswith("DIVIDE="))
    combine_values = combine.removeprefix("COMBINE=").split()
    divide_values = divide.removeprefix("DIVIDE=").split()

    assert int(combine_values[0]) == -1
    assert [float(value) for value in combine_values[1:]] == pytest.approx(
        [0.03, 0.03, 0.0]
    )
    assert int(divide_values[0]) == -2
    assert [float(value) for value in divide_values[1:]] == pytest.approx(
        [0.02, 0.01, 0.01]
    )
