"""Check the production flood/ocean evaporation viscosity polynomial numerically."""
from pathlib import Path
import re
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("relative_path", [
    "extends/CaMa/src/MOD_CaMa_colmCaMa.F90", "main/MOD_SimpleOcean.F90",
])
def test_surface_viscosity_uses_celsius(tmp_path, relative_path):
    compiler = require_runnable_fortran_compiler(tmp_path)
    source = (ROOT / relative_path).read_text()
    expression = re.search(r"(?im)^\s*visa\s*=([^\n]+)", source).group(1)
    probe = tmp_path / "viscosity.f90"
    probe.write_text(f"""program viscosity
use MOD_Precision
use MOD_Const_Physical, only: tfrz
implicit none
real(r8) :: tm, visa
integer :: i
do i = -20, 40, 20
  tm = tfrz + real(i, r8)
  visa = {expression}
  print '(ES24.16)', visa
end do
end program viscosity
""")
    executable = tmp_path / "viscosity"
    subprocess.run([
        compiler, "-fdefault-real-8", "-fcheck=all",
        "-ffpe-trap=invalid,zero,overflow", "-J", str(tmp_path),
        str(ROOT / "share/MOD_Precision.F90"),
        str(ROOT / "main/MOD_Const_Physical.F90"), str(probe),
        "-o", str(executable),
    ], cwd=tmp_path, check=True, capture_output=True, text=True)
    actual = subprocess.check_output([str(executable)], text=True)
    expected = [1.326e-5 * (1 + 6.542e-3 * c + 8.301e-6 * c**2 - 4.84e-9 * c**3)
                for c in (-20, 0, 20, 40)]
    assert [float(value) for value in actual.split()] == pytest.approx(expected, rel=1e-12, abs=0)
