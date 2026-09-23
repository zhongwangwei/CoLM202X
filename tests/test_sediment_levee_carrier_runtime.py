"""The first gross levee exchange and sediment donor share a pre-split carrier."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_first_repartition_preserves_pre_split_carrier(tmp_path: Path) -> None:
    source = SOURCE.read_text()
    kernel = "   SUBROUTINE sediment_levee_repartition(" + source.split(
        "   SUBROUTINE sediment_levee_repartition(", 1
    )[1].split("   END SUBROUTINE sediment_levee_repartition", 1)[0] + (
        "   END SUBROUTINE sediment_levee_repartition\n"
    )
    diag = source.split("   SUBROUTINE sediment_diag_accumulate(", 1)[1].split(
        "   END SUBROUTINE sediment_diag_accumulate", 1
    )[0]
    assert "IF (sed_acc_time(i) == 0._r8 .and. .not. sed_acc_pre_repartition_start(i))" in diag
    code = r"""
module MOD_Grid_RiverLakeNetwork
 integer, parameter :: numucat=1
end module
module probe
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer,parameter :: r8=kind(1.d0)
 real(r8),parameter :: SED_BALANCE_ABS_TOL=1.e-10_r8,SED_BALANCE_REL_TOL=1.e-10_r8
 logical :: p_is_worker=.true.
 real(r8),allocatable :: sed_acc_to_protected(:),sed_acc_from_protected(:)
 real(r8) :: sed_acc_time(1),sed_acc_rivsto_start(1),sed_acc_protected_start(1)
 logical :: sed_acc_pre_repartition_start(1)
 contains
 logical function sediment_particle_enabled()
  sediment_particle_enabled=.true.
 end function
 subroutine CoLM_stop(message)
  character(len=*),intent(in) :: message
  print *,message
  error stop
 end subroutine
""" + kernel + r"""
end module
program main
 use probe
 implicit none
 allocate(sed_acc_to_protected(1),sed_acc_from_protected(1))
 sed_acc_time=0;sed_acc_to_protected=0;sed_acc_from_protected=0
 sed_acc_pre_repartition_start=.false.
 call sediment_levee_repartition(1,100._r8,0._r8,50._r8,50._r8)
 if (abs(sed_acc_rivsto_start(1)-100._r8)>1.e-12_r8) error stop 'pre-split donor'
 if (abs(sed_acc_protected_start(1))>1.e-12_r8) error stop 'pre-split protected'
 if (abs(sed_acc_to_protected(1)-50._r8)>1.e-12_r8) error stop 'gross outward'
 ! A later repartition within this morphology window cannot overwrite start.
 call sediment_levee_repartition(1,50._r8,50._r8,60._r8,40._r8)
 if (abs(sed_acc_rivsto_start(1)-100._r8)>1.e-12_r8) error stop 'carrier overwritten'
 if (abs(sed_acc_from_protected(1)-10._r8)>1.e-12_r8) error stop 'gross return'
 ! M=1, V=100, Q=50 means 0.5 transported, not 1.0 using post-split V=50.
 if (abs(1._r8*sed_acc_to_protected(1)/sed_acc_rivsto_start(1)-0.5_r8)>1.e-12_r8) &
  error stop 'initial concentration'
 ! A half-routing checkpoint may begin with all water on the protected side.
 sed_acc_time=0;sed_acc_to_protected=0;sed_acc_from_protected=0
 sed_acc_pre_repartition_start=.false.
 call sediment_levee_repartition(1,0._r8,10._r8,5._r8,5._r8)
 if (abs(sed_acc_rivsto_start(1)-10._r8)>1.e-12_r8) error stop 'protected-only total'
 if (abs(sed_acc_protected_start(1)-10._r8)>1.e-12_r8) error stop 'protected-only donor'
 if (abs(sed_acc_from_protected(1)-5._r8)>1.e-12_r8) error stop 'protected-only gross'
 if (.not.sed_acc_pre_repartition_start(1)) error stop 'protected-only checkpoint flag'
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
