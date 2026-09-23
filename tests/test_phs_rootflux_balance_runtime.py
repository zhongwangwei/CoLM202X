from pathlib import Path
import subprocess

import pytest

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
SUBPROCESS_TIMEOUT = 60


@pytest.fixture(scope="module")
def phs_balance_driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("phs_rootflux_balance")
    compiler = require_runnable_fortran_compiler(workdir)
    (workdir / "define.h").write_text("", encoding="utf-8")
    (workdir / "precision.f90").write_text(
        "module MOD_Precision\n implicit none\n integer,parameter::r8=selected_real_kind(12)\nend module\n",
        encoding="utf-8",
    )
    (workdir / "driver.f90").write_text(
        """
program phs_rootflux_balance_driver
  use MOD_PHSRootfluxBalance, only: balance_phs_rootflux
  use MOD_Precision, only: r8
  implicit none
  character(len=32) :: case_name

  call get_command_argument(1, case_name)
  select case (trim(case_name))
  case ('balanced')
    call check_balanced_is_unchanged()
  case ('positive')
    call check_positive_fluxes_are_clipped_and_scaled()
  case ('fallback')
    call check_fallback_weights_are_used_without_positive_flux()
  case ('uniform')
    call check_uniform_split_is_used_without_positive_weight()
  case default
    error stop 'unknown case'
  end select
  write(*,'(A)') 'PASS'

contains
  subroutine assert_close(actual, expected, label)
    real(r8), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > 1.e-10_r8) then
      write(*,'(A,2ES24.16)') trim(label), actual, expected
      error stop 1
    endif
  end subroutine assert_close

  subroutine assert_sum_matches_etr(rootflux, etr)
    real(r8), intent(in) :: rootflux(:), etr
    call assert_close(sum(rootflux), etr, 'sum(rootflux) /= etr')
  end subroutine assert_sum_matches_etr

  subroutine check_balanced_is_unchanged()
    real(r8) :: etr, rootflux(2), original(2), fallback(2)
    etr = 2._r8
    rootflux = [0.5_r8, 1.5_r8]
    original = rootflux
    fallback = [1._r8, 1._r8]
    call balance_phs_rootflux(1, 0, etr, rootflux, fallback, 'balanced')
    call assert_close(rootflux(1), original(1), 'balanced layer 1 changed')
    call assert_close(rootflux(2), original(2), 'balanced layer 2 changed')
    call assert_sum_matches_etr(rootflux, etr)
  end subroutine check_balanced_is_unchanged

  subroutine check_positive_fluxes_are_clipped_and_scaled()
    real(r8) :: etr, rootflux(3), fallback(3)
    etr = 10._r8
    rootflux = [-1._r8, 2._r8, 3._r8]
    fallback = [1._r8, 1._r8, 1._r8]
    call balance_phs_rootflux(2, 0, etr, rootflux, fallback, 'positive')
    call assert_close(rootflux(1), 0._r8, 'negative rootflux was not clipped')
    call assert_close(rootflux(2), 4._r8, 'positive layer 2 not scaled')
    call assert_close(rootflux(3), 6._r8, 'positive layer 3 not scaled')
    call assert_sum_matches_etr(rootflux, etr)
  end subroutine check_positive_fluxes_are_clipped_and_scaled

  subroutine check_fallback_weights_are_used_without_positive_flux()
    real(r8) :: etr, rootflux(3), fallback(3)
    etr = 8._r8
    rootflux = [-1._r8, 0._r8, -2._r8]
    fallback = [1._r8, 2._r8, 1._r8]
    call balance_phs_rootflux(3, 0, etr, rootflux, fallback, 'fallback')
    call assert_close(rootflux(1), 2._r8, 'fallback layer 1')
    call assert_close(rootflux(2), 4._r8, 'fallback layer 2')
    call assert_close(rootflux(3), 2._r8, 'fallback layer 3')
    call assert_sum_matches_etr(rootflux, etr)
  end subroutine check_fallback_weights_are_used_without_positive_flux

  subroutine check_uniform_split_is_used_without_positive_weight()
    real(r8) :: etr, rootflux(3), fallback(3)
    etr = 9._r8
    rootflux = [-1._r8, 0._r8, -2._r8]
    fallback = [0._r8, 0._r8, 0._r8]
    call balance_phs_rootflux(4, 0, etr, rootflux, fallback, 'uniform')
    call assert_close(rootflux(1), 3._r8, 'uniform layer 1')
    call assert_close(rootflux(2), 3._r8, 'uniform layer 2')
    call assert_close(rootflux(3), 3._r8, 'uniform layer 3')
    call assert_sum_matches_etr(rootflux, etr)
  end subroutine check_uniform_split_is_used_without_positive_weight
end program phs_rootflux_balance_driver
""",
        encoding="utf-8",
    )
    executable = workdir / "phs_rootflux_balance_driver"
    compiled = subprocess.run(
        [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(workdir),
            str(workdir / "precision.f90"),
            str(ROOT / "extends/interception/MOD_PHSRootfluxBalance.F90"),
            str(workdir / "driver.f90"),
            "-o",
            str(executable),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    return executable


def run_case(executable: Path, case_name: str) -> None:
    result = subprocess.run(
        [str(executable), case_name],
        cwd=executable.parent,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines()[-1] == "PASS"


def test_balanced_rootflux_is_not_changed(phs_balance_driver: Path) -> None:
    run_case(phs_balance_driver, "balanced")


def test_positive_rootfluxes_are_clipped_and_scaled_to_etr(phs_balance_driver: Path) -> None:
    run_case(phs_balance_driver, "positive")


def test_fallback_weights_are_used_when_no_positive_rootflux_exists(phs_balance_driver: Path) -> None:
    run_case(phs_balance_driver, "fallback")


def test_uniform_split_is_used_when_fallback_weights_sum_to_zero(phs_balance_driver: Path) -> None:
    run_case(phs_balance_driver, "uniform")
