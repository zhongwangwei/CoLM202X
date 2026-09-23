from pathlib import Path
import subprocess
import textwrap

ROOT = Path(__file__).resolve().parents[1]
MOD = ROOT / "extends" / "CaMa" / "src" / "cmf_coupling_budget_mod.F90"


def _compile_and_run(tmp_path: Path, source: str) -> subprocess.CompletedProcess[str]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    src = tmp_path / "probe.f90"
    exe = tmp_path / "probe"
    src.write_text(textwrap.dedent(source))
    subprocess.run(
        ["gfortran", "-ffree-line-length-none", "-fcheck=all", str(MOD), str(src), "-o", str(exe)],
        check=True,
        cwd=tmp_path,
        text=True,
        capture_output=True,
    )
    return subprocess.run([str(exe)], cwd=ROOT, text=True, capture_output=True)


def test_budget_rejects_donor_overdraw_nonfinite_and_invalid_indices(tmp_path: Path) -> None:
    checks = [
        ("overdraw", "call budget_debit(reshape([11d0],[1,1]),reshape([0d0],[1,1]),storage,e,i)", "exceeded its published grid credit"),
        ("nonfinite", "w=0d0; w=w/w; call budget_init(x,y,w,area)", "invalid exchange weights"),
        ("bad_index", "call budget_init(reshape([2],[1,1]),y,w,area)", "exchange index out of range"),
    ]
    for name, statement, message in checks:
        result = _compile_and_run(
            tmp_path / name,
            f"""
            program probe
              use, intrinsic :: iso_fortran_env, only: real64
              use cmf_coupling_budget_mod
              implicit none
              integer :: x(1,1), y(1,1)
              real(real64) :: w(1,1), area(1,1), storage(1), e(1), i(1), d(1,1), f(1,1)
              x=1; y=1; w=1; area=1; storage=10
              call budget_init(x,y,w,area)
              call budget_publish(storage,storage,d,f)
              {statement}
              error stop 'probe should have failed: {name}'
            end program
            """,
        )
        assert result.returncode != 0, result.stdout + result.stderr
        assert message.lower() in (result.stdout + result.stderr).lower()


def test_partial_coverage_conserves_covered_flux_only(tmp_path: Path) -> None:
    result = _compile_and_run(
        tmp_path / "partial",
        """
        program probe
          use, intrinsic :: iso_fortran_env, only: real64
          use cmf_coupling_budget_mod
          implicit none
          integer :: x(2,1), y(2,1)
          real(real64) :: w(2,1), area(1,1), gridflow(1,1), flow(2)
          x=1; y=1; w=reshape([2d0,6d0],[2,1]); area=reshape([100d0],[1,1])
          gridflow=40d0
          call budget_init(x,y,w,area)
          call budget_runoff(gridflow,flow)
          if (abs(sum(flow)-40d0)>1d-10) error stop 'covered flux not conserved'
          if (abs(flow(1)-10d0)>1d-10 .or. abs(flow(2)-30d0)>1d-10) error stop 'wrong covered-flux split'
        end program
        """,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_overlapping_donors_debit_by_published_provenance(tmp_path: Path) -> None:
    result = _compile_and_run(
        tmp_path / "provenance",
        """
        program probe
          use, intrinsic :: iso_fortran_env, only: real64
          use cmf_coupling_budget_mod
          implicit none
          integer :: x(2,1), y(2,1)
          real(real64) :: w(2,1), area(1,1), volume(2), floodarea(2), depth(1,1), fraction(1,1)
          real(real64) :: storage(2), evap(1,1), infil(1,1), eused(2), iused(2)
          x=1; y=1; w=1; area=100d0
          volume=[20d0,80d0]; floodarea=[10d0,20d0]; storage=volume; evap=10d0; infil=15d0
          call budget_init(x,y,w,area)
          call budget_publish(volume,floodarea,depth,fraction)
          if (abs(depth(1,1)-1d0)>1d-12) error stop 'wrong published depth'
          if (abs(fraction(1,1)-0.3d0)>1d-12) error stop 'wrong published fraction'
          call budget_debit(evap,infil,storage,eused,iused)
          if (abs(eused(1)-2d0)>1d-12 .or. abs(eused(2)-8d0)>1d-12) error stop 'evap provenance lost'
          if (abs(iused(1)-3d0)>1d-12 .or. abs(iused(2)-12d0)>1d-12) error stop 'infil provenance lost'
          if (abs(storage(1)-15d0)>1d-12 .or. abs(storage(2)-60d0)>1d-12) error stop 'storage debit wrong'
        end program
        """,
    )
    assert result.returncode == 0, result.stdout + result.stderr
