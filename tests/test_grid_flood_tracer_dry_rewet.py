"""A dry protected-side nonvolatile residue survives and redissolves on rewet."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_RiverLake.F90"


def test_production_levee_repartition_keeps_dry_residue(tmp_path: Path):
    source = SOURCE.read_text()
    kernel = "   SUBROUTINE levee_tracer_repartition (" + source.split(
        "   SUBROUTINE levee_tracer_repartition (", 1
    )[1].split("   END SUBROUTINE levee_tracer_repartition", 1)[0] + (
        "   END SUBROUTINE levee_tracer_repartition\n"
    )
    code = """
module MOD_Tracer_Defs
 integer, parameter :: r8=kind(1.d0)
 real(r8), parameter :: trc_tiny=1.e-30_r8
 type tracer_desc
  real(r8) :: max_dissolved_conc=huge(1.d0)
 end type
 type(tracer_desc) :: tracers(1)
 contains
 logical function tracer_uses_land_water_transport(itrc)
  integer, intent(in) :: itrc
  tracer_uses_land_water_transport=.true.
 end function
 logical function tracer_has_dissolved_limit(itrc)
  integer, intent(in) :: itrc
  tracer_has_dissolved_limit=.false.
 end function
end module
module probe
 use MOD_Tracer_Defs
 implicit none
 integer, parameter :: ntracers=1
 real(r8), allocatable :: trc_mass(:,:),trc_levsto(:,:)
 contains
 subroutine equilibrate_river_tracer_cell(icell,visible_water,protected_water)
  integer, intent(in) :: icell
  real(r8), intent(in) :: visible_water,protected_water
 end subroutine
""" + kernel + """
end module
program main
 use probe
 implicit none
 allocate(trc_mass(1,1),trc_levsto(1,1))
 ! Evaporation removed all protected water but not the nonvolatile solute.
 trc_mass=10._r8;trc_levsto=2._r8
 call levee_tracer_repartition(1,10._r8,0._r8,9._r8,1._r8)
 if(abs(trc_mass(1,1)-9._r8)>1.e-12_r8) error stop 'visible rewet debit'
 if(abs(trc_levsto(1,1)-3._r8)>1.e-12_r8) error stop 'dry residue lost'
 call levee_tracer_repartition(1,9._r8,1._r8,10._r8,0._r8)
 if(abs(trc_mass(1,1)-12._r8)>1.e-12_r8) error stop 'redissolution missing'
 if(abs(trc_levsto(1,1))>1.e-12_r8) error stop 'protected return'
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
