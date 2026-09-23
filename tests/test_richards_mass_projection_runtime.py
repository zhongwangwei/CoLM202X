"""Execute the production Richards finite-volume projection on small columns."""

from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "main/HYDRO/MOD_Hydro_SoilWater.F90"


def test_implicit_projection_is_atomic_bounded_and_conservative():
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran unavailable")

    source = SOURCE.read_text()
    start = source.index("   SUBROUTINE project_richards_liquid_water (")
    end = source.index("   END SUBROUTINE project_richards_liquid_water", start)
    helper = source[start : end + len("   END SUBROUTINE project_richards_liquid_water")]
    program = f"""
module probe_module
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer, parameter :: r8=selected_real_kind(15,307)
 integer, parameter :: BC_RAINFALL=2, BC_DRAINAGE=4
 contains
{helper}
end module
program probe
 use probe_module
 implicit none
 real(r8), parameter :: budget=64._r8*epsilon(1._r8)*500._r8
 real(r8) :: dz(10), vs(10), vr(10), q(0:10), wf(10), vl(10), wt(10)
 real(r8) :: wfm(10), vlm(10), wtm(10), before(10)
 logical :: ok
 integer :: j
 dz=100._r8; vs=.5_r8; vr=.1_r8
 q=0._r8; wf=0._r8; wt=0._r8; wfm=0._r8; wtm=0._r8
 vlm=.3_r8; vl=vlm

 ! A tiny real interface inflow must become soil storage, not vanish.
 q(0)=1.e-5_r8
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (.not.ok .or. abs(vl(1)-(.3_r8+1.e-7_r8))>1.e-14_r8) error stop 1
 if (any(abs(vl(2:)-.3_r8)>1.e-14_r8)) error stop 2

 ! Multiple substep fronts: previous water includes both wetting and water table.
 q=0._r8; vlm=.3_r8; vl=vlm; wfm=0._r8; wtm=0._r8
 wfm(1)=10._r8; wtm(1)=20._r8; wf(1)=12._r8; wt(1)=19._r8
 q(0)=.01_r8
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (.not.ok) error stop 3
 if (abs(vl(1)*69._r8+vs(1)*31._r8-(vlm(1)*70._r8+vs(1)*30._r8+.01_r8)) > budget) error stop 4

 ! Several individually tiny errors cannot evade the total-column budget.
 q=0._r8; wf=0._r8; wt=0._r8; wfm=0._r8; wtm=0._r8
 vlm=.3_r8; vl=vlm+9.e-13_r8
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (.not.ok .or. abs(sum((vl-vlm)*dz))>budget) error stop 5

 ! A fully saturated layer between unsaturated ones remains untouched;
 ! water_balance merges its row, whereas the projection checks raw layers.
 q=0._r8; q(0)=1.e-5_r8; vlm=.3_r8; vl=vlm
 wt=0._r8; wtm=0._r8; wt(5)=100._r8; wtm(5)=100._r8
 vl(5)=vs(5); vlm(5)=vs(5)
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (.not.ok .or. vl(5)/=vs(5)) error stop 9
 if (abs(vl(1)-(.3_r8+1.e-7_r8))>1.e-14_r8) error stop 10

 ! A nearly saturated layer cannot absorb a larger accepted influx.
 q=0._r8; vlm=.3_r8; vl=vlm; wt=0._r8
 wt(2)=99.9_r8; wtm=wt; vl(2)=.49999_r8; vlm(2)=vl(2)
 q(1)=1.e-4_r8; q(2)=0._r8
 before=vl
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (ok .or. any(vl/=before)) error stop 6

 ! A bad second layer must not partly commit the valid first-layer change.
 q=0._r8; q(0)=1.e-5_r8; q(1)=0._r8; q(2)=-1.e-4_r8
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (ok .or. any(vl/=before)) error stop 7

 ! Pond boundary inconsistency is not reattributed to an interior layer.
 q=0._r8
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,BC_RAINFALL,1.e-5_r8,3, &
      0._r8,0._r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (ok .or. any(vl/=before)) error stop 8

 ! Drainage/aquifer residual cannot be masked by soil-liquid projection.
 call project_richards_liquid_water(1,10,dz,1._r8,vs,vr,q,3,0._r8,BC_DRAINAGE, &
      0._r8,1.e-5_r8,wf,vl,wt,0._r8,0._r8,wfm,vlm,wtm,budget,ok)
 if (ok .or. any(vl/=before)) error stop 11
 print *, 'projection probe PASS'
end program
"""
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "probe.f90"
        binary = Path(tmp) / "probe"
        src.write_text(program)
        build = subprocess.run(
            [compiler, "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-ffree-line-length-none", str(src), "-o", str(binary)],
            capture_output=True,
            text=True,
        )
        assert build.returncode == 0, build.stderr
        result = subprocess.run([str(binary)], capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "projection probe PASS" in result.stdout
