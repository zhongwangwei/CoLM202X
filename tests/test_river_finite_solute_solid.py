"""Compile the production river equilibrium and levee kernels for finite solute."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_RiverLake.F90"
FLOW = Path(__file__).resolve().parents[1] / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90"


def _subroutine(source: str, name: str) -> str:
    start = source.index(f"   SUBROUTINE {name} (")
    end = source.index(f"   END SUBROUTINE {name}", start)
    return source[start : source.index("\n", end) + 1]


def test_restart_and_legacy_fast_path_contract():
    source = SOURCE.read_text()
    restart_reader = _subroutine(source, "read_tracer_restart")
    restart_writer = _subroutine(source, "write_tracer_restart")
    initializer = _subroutine(source, "river_lake_tracer_init")
    assert "RIVER_TRACER_RESTART_SCHEMA_VERSION = 2" in source
    assert "schema /= 1 .and. schema /= RIVER_TRACER_RESTART_SCHEMA_VERSION" in source
    assert "river_restart_schema_loaded >= 2 .and. tracer_has_dissolved_limit(itrc)" in restart_reader
    assert "CALL probe_riverlake_restart_vector(file_restart, trim(varname), .true." in restart_reader
    assert "'trc_solid_'" in restart_writer and "'trc_levsto_solid_'" in restart_writer
    assert "IF (has_finite_solute) THEN" in initializer
    assert "IF (allocated(trc_solid))" in _subroutine(source, "tracer_substep")
    flow = FLOW.read_text()
    assert "sum(trc_solid(itrc,:)) + sum(trc_levsto_solid(itrc,:))" in flow
    assert "finite-solubility river solute needs a solid residue pool" not in flow


def test_finite_solute_dry_rewet_cap_and_levee_split(tmp_path: Path):
    source = SOURCE.read_text()
    code = """
module MOD_Tracer_Defs
 integer, parameter :: r8=kind(1.d0)
 real(r8), parameter :: trc_tiny=1.e-30_r8, trc_water_min_for_ratio=1.e-12_r8
 type tracer_desc
  real(r8) :: max_dissolved_conc=2._r8
 end type
 type(tracer_desc) :: tracers(1)
 contains
 logical function tracer_uses_land_water_transport(itrc)
  integer, intent(in) :: itrc
  tracer_uses_land_water_transport=.true.
 end function
 logical function tracer_has_dissolved_limit(itrc)
  integer, intent(in) :: itrc
  tracer_has_dissolved_limit=.true.
 end function
 subroutine tracer_equilibrate_dissolved(itrc,water_mass,dissolved_mass,solid_mass)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: water_mass
  real(r8), intent(inout) :: dissolved_mass,solid_mass
  real(r8) :: capacity,total_mass
  total_mass=dissolved_mass+solid_mass
  capacity=0._r8
  if(water_mass>trc_water_min_for_ratio) capacity=tracers(itrc)%max_dissolved_conc*water_mass
  dissolved_mass=min(total_mass,capacity)
  solid_mass=total_mass-dissolved_mass
 end subroutine
end module
module probe
 use MOD_Tracer_Defs
 implicit none
 integer, parameter :: ntracers=1
 real(r8), parameter :: TRC_RESTART_NEGATIVE_DUST=1.e-12_r8
 real(r8), allocatable :: trc_mass(:,:),trc_levsto(:,:),trc_solid(:,:),trc_levsto_solid(:,:)
 contains
 subroutine CoLM_stop(message)
  character(len=*), intent(in) :: message
  print *, message
  error stop
 end subroutine
""" + _subroutine(source, "equilibrate_river_tracer_cell") + _subroutine(
        source, "levee_tracer_repartition"
    ) + """
end module
program main
 use probe
 implicit none
 allocate(trc_mass(1,1),trc_levsto(1,1),trc_solid(1,1),trc_levsto_solid(1,1))
 trc_mass=12._r8; trc_levsto=3._r8; trc_solid=0._r8; trc_levsto_solid=0._r8
 call equilibrate_river_tracer_cell(1,4._r8,1._r8)
 if(abs(trc_mass(1,1)-8._r8)>1.e-12_r8) error stop 'visible cap'
 if(abs(trc_solid(1,1)-4._r8)>1.e-12_r8) error stop 'visible precipitate'
 if(abs(trc_levsto(1,1)-2._r8)>1.e-12_r8) error stop 'protected cap'
 if(abs(trc_levsto_solid(1,1)-1._r8)>1.e-12_r8) error stop 'protected precipitate'
 call equilibrate_river_tracer_cell(1,0._r8,0._r8)
 if(abs(trc_mass(1,1))+abs(trc_levsto(1,1))>1.e-12_r8) error stop 'dry mobile mass'
 if(abs(trc_solid(1,1)-12._r8)+abs(trc_levsto_solid(1,1)-3._r8)>1.e-12_r8) error stop 'dry conservation'
 call equilibrate_river_tracer_cell(1,5._r8,2._r8)
 if(abs(trc_mass(1,1)-10._r8)+abs(trc_solid(1,1)-2._r8)>1.e-12_r8) error stop 'visible rewet'
 if(abs(trc_levsto(1,1)-3._r8)+abs(trc_levsto_solid(1,1))>1.e-12_r8) error stop 'protected rewet'
 call levee_tracer_repartition(1,5._r8,2._r8,4._r8,3._r8)
 if(abs(trc_mass(1,1)-8._r8)+abs(trc_solid(1,1)-2._r8)>1.e-12_r8) error stop 'levee visible'
 if(abs(trc_levsto(1,1)-5._r8)+abs(trc_levsto_solid(1,1))>1.e-12_r8) error stop 'levee protected'
 if(abs(sum(trc_mass)+sum(trc_levsto)+sum(trc_solid)+sum(trc_levsto_solid)-15._r8)>1.e-12_r8) &
  error stop 'levee conservation'
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
