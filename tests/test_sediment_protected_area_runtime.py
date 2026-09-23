"""Protected settling area excludes the unprotected levee-side strip."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
FLOW = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text()


def test_protected_wet_area_at_crest_and_full_flood(tmp_path: Path) -> None:
    assignment = FLOW.split("particle_protected_area(i) = min(", 1)[1].split(") * topo_area(i))", 1)[0]
    expression = "min(" + assignment + ") * topo_area(i))"
    code = """
program probe
 implicit none
 integer,parameter::r8=kind(1.d0)
 integer::i
 real(r8)::particle_protected_area(1),levee_floodarea(1),levee_frc_data(1),topo_area(1)
 i=1;topo_area=100._r8;levee_frc_data=.3_r8
 levee_floodarea=30._r8
 particle_protected_area(i) = """ + expression + """
 if(abs(particle_protected_area(1))>1.e-12_r8) error stop 'crest protection area'
 levee_floodarea=100._r8
 particle_protected_area(i) = """ + expression + """
 if(abs(particle_protected_area(1)-70._r8)>1.e-12_r8) error stop 'full protection area'
end program
"""
    path = tmp_path / "probe.f90"
    path.write_text(code)
    compiler = require_runnable_fortran_compiler(tmp_path)
    built = subprocess.run(
        [compiler, "-ffree-line-length-0", str(path), "-o", "probe"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(tmp_path / "probe")], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
