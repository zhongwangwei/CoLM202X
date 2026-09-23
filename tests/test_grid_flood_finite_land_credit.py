"""Exercise the production patch-credit expression with finite-solubility flood water."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


MAIN = Path(__file__).resolve().parents[1] / "main/CoLMMAIN.F90"


def test_finite_flood_infiltration_caps_actual_land_credit(tmp_path: Path):
    source = MAIN.read_text()
    expression = source.split("               flood_input_tracer = (flood_tracer_credit_patch(:,ipatch)", 1)[1]
    expression = "               flood_input_tracer = (flood_tracer_credit_patch(:,ipatch)" + expression.split(
        "               ENDDO", 1
    )[0] + "               ENDDO\n"
    assert "tracer_has_dissolved_limit(itrc_loc)" in expression
    code = """
module probe
 implicit none
 integer, parameter :: r8=kind(1.d0), ntracers=3
 type tracer_info
  real(r8) :: max_dissolved_conc
 end type
 type(tracer_info) :: tracers(ntracers)
 real(r8) :: flood_input_tracer(ntracers), flood_tracer_credit_patch(ntracers,1)
 real(r8) :: flood_tracer_evap_patch(ntracers,1), flood_credit_patch(1)
 real(r8) :: qinfl_fld, deltim, fevpg_fld
 integer :: ipatch=1, itrc_loc
 contains
 logical function tracer_has_dissolved_limit(itrc)
  integer, intent(in) :: itrc
  tracer_has_dissolved_limit = itrc == 1
 end function
 subroutine transfer()
""" + expression + """
 end subroutine
end module
program main
 use probe
 implicit none
 ! W=10 mm, E=4 mm, I=5 mm. Finite source is supersaturated
 ! after evaporation; only Cmax*I reaches land, not all source mass.
 deltim=1._r8; flood_credit_patch=0.01_r8; fevpg_fld=4._r8
 qinfl_fld=5._r8; tracers(1)%max_dissolved_conc=0.1_r8
 flood_tracer_credit_patch(:,1)=[1.2_r8,1.2_r8,1.2_r8]
 flood_tracer_evap_patch(:,1)=[0._r8,0._r8,0.3_r8]
 call transfer()
 if(abs(flood_input_tracer(1)-0.5_r8)>1.e-12_r8) error stop 'finite cap'
 if(abs(flood_input_tracer(2)-1._r8)>1.e-12_r8) error stop 'unlimited changed'
 if(abs(flood_input_tracer(3)-0.75_r8)>1.e-12_r8) error stop 'isotope changed'
 if(abs((1.2_r8-flood_input_tracer(1))-0.7_r8)>1.e-12_r8) &
  error stop 'donor residual not conserved'
 flood_tracer_credit_patch(1,1)=0._r8
 call transfer()
 if(flood_input_tracer(1)/=0._r8) error stop 'empty source fabricated mass'
end program
"""
    path = tmp_path / "probe.f90"
    path.write_text(code)
    compiler = require_runnable_fortran_compiler(tmp_path)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", "-fcheck=all", str(path), "-o", "probe"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(tmp_path / "probe")], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "flood_tracer_land_patch(:,ipatch) = flood_input_tracer" in source
