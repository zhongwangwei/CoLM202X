"""A morphology donor's initial concentration feeds main and BIF faces equally."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_two_outlets_share_snapshot_concentration_and_limited_stock(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    helper = "   real(r8) FUNCTION joint_sediment_scale(" + source.split(
        "   real(r8) FUNCTION joint_sediment_scale(", 1
    )[1].split("   END FUNCTION joint_sediment_scale", 1)[0] + (
        "   END FUNCTION joint_sediment_scale\n"
    )
    code = """
module probe
 implicit none
 integer,parameter::r8=kind(1.d0)
 contains
""" + helper + """
end module
program main
 use probe
 implicit none
 real(r8)::stock,carrier,main_demand,bif_demand,scale
 stock=1._r8;carrier=10._r8
 main_demand=2._r8*stock/carrier
 bif_demand=2._r8*stock/carrier
 scale=joint_sediment_scale(stock,main_demand+bif_demand)
 if(abs(main_demand*scale-.2_r8)>1.e-14_r8) error stop 'main diluted by BIF debit'
 if(abs(bif_demand*scale-.2_r8)>1.e-14_r8) error stop 'BIF demand'
 stock=.3_r8
 scale=joint_sediment_scale(stock,main_demand+bif_demand)
 if(abs(main_demand*scale-.15_r8)>1.e-14_r8) error stop 'main not proportionally limited'
 if(abs(bif_demand*scale-.15_r8)>1.e-14_r8) error stop 'BIF not proportionally limited'
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
