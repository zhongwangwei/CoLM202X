"""One sediment pass per morphology interval is reported without CFL alarms."""

from pathlib import Path
import re
import subprocess

from fortran_test_support import require_runnable_fortran_compiler

SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def routine(source: str, name: str) -> str:
    match = re.search(rf"^\s*SUBROUTINE {name}\b.*?^\s*END SUBROUTINE {name}\b", source, re.M | re.S | re.I)
    assert match, name
    return match.group()


def test_interval_diagnostics_count_one_pass_per_interval(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    src = SOURCE.read_text()
    calc = routine(src, "grid_sediment_calc")
    assert calc.count("iter_adv = iter_adv + 1") == 1
    assert "dt_adv_remaining" not in calc
    assert "CALL report_sediment_advection_substeps(iter_sed, iter_adv, call_dt_min)" in calc
    code = """
module probe
 implicit none
 integer, parameter :: r8=kind(1.d0), SED_ADV_DIAG_PERIODIC=5
 integer :: sed_st_calls=0,sed_st_morph=0,sed_st_adv=0,p_iam_worker=0
 real(r8) :: sed_st_dt_min=huge(1._r8)
 contains
""" + routine(src, "report_sediment_advection_substeps") + """
end module
program main
 use probe
 integer :: i
 do i=1,6
  call report_sediment_advection_substeps(1,1,3600._r8)
 enddo
 if (sed_st_calls/=6 .or. sed_st_morph/=6 .or. sed_st_adv/=6) error stop 1
 if (abs(sed_st_dt_min-3600._r8)>1.e-12_r8) error stop 2
end program
"""
    path = tmp_path / "probe.F90"
    path.write_text(code)
    built = subprocess.run([compiler, "-cpp", "-ffree-line-length-0", str(path), "-o", "probe"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run(["./probe"], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert ran.stdout.count("Sediment advection: morph_intervals=1 worker_passes=1") == 5
