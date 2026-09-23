from pathlib import Path
import re
import shutil
import subprocess
import textwrap

import pytest

from fortran_test_support import netcdf_fortran_flags

ROOT = Path(__file__).resolve().parents[1]

def flat(path: str) -> str:
    text = (ROOT / path).read_text()
    text = re.sub(r"&\s*", " ", text)
    return re.sub(r"\s+", " ", text.lower())


def test_output_average_zero_window_is_explicitly_cleared() -> None:
    src = flat("extends/CaMa/src/cmf_calc_diag_mod.F90")
    routine = src[src.index("subroutine cmf_diag_getave_output"):src.index("end subroutine cmf_diag_getave_output")]
    assert "if (nadd_out <= 0._jprb) then" in routine
    assert "d2rivout_oavg(:,:) = 0._jprb" in routine
    assert "d2winfiltex_oavg(:,:) = 0._jprb" in routine
    assert "return" in routine
    assert routine.index("if (nadd_out <= 0._jprb) then") < routine.index("/ real(nadd_out,kind=jprb)")


def test_adaptive_average_zero_window_is_explicitly_cleared() -> None:
    src = flat("extends/CaMa/src/cmf_calc_diag_mod.F90")
    routine = src[src.index("subroutine cmf_diag_getave_adpstp"):src.index("end subroutine cmf_diag_getave_adpstp")]
    assert "if (nadd_adp <= 0._jprb) then" in routine
    assert "d2rivout_aavg(:,:) = 0._jprb" in routine
    assert "d1pthflwsum_aavg(:) = 0._jprb" in routine
    assert routine.index("if (nadd_adp <= 0._jprb) then") < routine.index("/ real(nadd_adp,kind=jprb)")


def test_nsteps_uses_wide_elapsed_seconds_and_rejects_bad_windows() -> None:
    src = flat("extends/CaMa/src/cmf_ctrl_time_mod.F90")
    routine = src[src.index("subroutine cmf_time_init"):src.index("end subroutine cmf_time_init")]
    assert "elapsed_min_8" in routine
    assert "elapsed_sec_8=elapsed_min_8*60_jpib" in routine
    assert "kminend <= kminstart" in routine
    assert "elapsed_sec_8 > int(huge(nsteps),kind=jpib) * int(dt,kind=jpib)" in routine
    assert "( (kminend-kminstart)*60_jpim )" not in routine


def test_restart_always_reconstructs_flood_stage_after_read() -> None:
    src = flat("extends/CaMa/src/cmf_drv_control_mod.F90")
    restart = src.index("if( lrestart )then call cmf_restart_init")
    stage_after = src.index("call cmf_physics_fldstg", restart)
    storage_only = src.index("if( lrestart .and. lstoonly )then")
    assert restart < stage_after < storage_only
    storage_block = src[storage_only:src.index("if ( loutini .and. loutput ) then")]
    assert "call cmf_physics_fldstg" not in storage_block
    assert "call cmf_calc_outpre" in storage_block


def test_default_flood_namelist_uses_adaptive_timestep() -> None:
    src = flat("run/CaMa/cama_flood.nml")
    assert "ladpstp = .true." in src


def test_adpstp_diag_keeps_simple_global_rate_contract() -> None:
    src = flat("extends/CaMa/src/cmf_calc_diag_mod.F90")
    routine = src[src.index("subroutine cmf_diag_avemax_adpstp"):src.index("end subroutine cmf_diag_avemax_adpstp")]
    assert "subroutine cmf_diag_avemax_adpstp " in routine
    assert "optional,intent(in)" not in routine
    assert "d2wevapex_aavg(iseq,1)= d2wevapex_aavg(iseq,1) +d2wevapex(iseq,1)*dt" in routine
    assert "d2winfiltex_aavg(iseq,1)= d2winfiltex_aavg(iseq,1) +d2winfiltex(iseq,1)*dt" in routine


def test_time_control_preserves_start_and_end_minutes() -> None:
    src = flat("extends/CaMa/src/cmf_ctrl_time_mod.F90")
    assert "integer(kind=jpim) :: smin" in src
    assert "namelist/nsimtime/ syear,smon,sday,shour,smin, eyear,emon,eday,ehour,emin" in src
    assert "ishhmm=shour*100_jpim+smin" in src
    assert "iehhmm=ehour*100_jpim+emin" in src
    assert "subroutine cmf_time_set_start_seconds(start_sec)" in src
    assert "mod(start_sec,60_jpim) /= 0" in src


def test_example_namelist_requests_double_precision_restart() -> None:
    src = flat("run/CaMa/cama_flood.nml")
    assert "lrestdbl = .true." in src
    assert "smin = 00" in src
    assert "emin = 00" in src


def test_partial_time_window_and_closed_basin_guards():
    time = flat("extends/CaMa/src/cmf_ctrl_time_mod.F90")
    driver = flat("extends/CaMa/src/cmf_drv_advance_mod.F90")
    couple = flat("extends/CaMa/src/MOD_CaMa_colmCaMa.F90")
    assert "nsteps=ceiling(" in time
    assert "dt=min(dt,real(kminend-kmin,kind=jprb)*60._jprb)" in driver
    assert "(coverage(c)>0._r8).neqv.(coverage(nextc)>0._r8)" in couple


def test_end_second_86400_is_accepted_by_actual_time_module(tmp_path: Path) -> None:
    gfortran = shutil.which("gfortran") or "/opt/homebrew/bin/gfortran"
    if not Path(gfortran).exists():
        pytest.skip("gfortran is not available")

    build = tmp_path / "build"
    build.mkdir()
    netcdf_include, netcdf_libs = netcdf_fortran_flags()
    flags = ["-cpp", "-free", "-fimplicit-none", "-ffree-line-length-none", f"-I{build}", *netcdf_include]
    sources = [
        "extends/CaMa/src/parkind1.F90",
        "extends/CaMa/src/yos_cmf_input.F90",
        "extends/CaMa/src/yos_cmf_time.F90",
        "extends/CaMa/src/yos_cmf_map.F90",
        "extends/CaMa/src/cmf_utils_mod.F90",
        "extends/CaMa/src/cmf_ctrl_time_mod.F90",
    ]
    objects: list[Path] = []
    for source in sources:
        src = ROOT / source
        obj = build / f"{src.stem}.o"
        subprocess.run(
            [gfortran, *flags, "-c", str(src), "-J", str(build), "-o", str(obj)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        objects.append(obj)

    probe = build / "time_end_86400_probe.F90"
    probe.write_text(
        textwrap.dedent(
            """
            PROGRAM time_end_86400_probe
              USE PARKIND1, ONLY: JPIM, JPRB
              USE YOS_CMF_INPUT, ONLY: DT, LOGNAM
              USE YOS_CMF_TIME, ONLY: YYYY0, ISHHMM, IEHHMM, KMINSTART, KMINEND, NSTEPS
              USE CMF_CTRL_TIME_MOD, ONLY: SYEAR, SMON, SDAY, EYEAR, EMON, EDAY, EHOUR, EMIN, &
                                         & CMF_TIME_SET_START_SECONDS, CMF_TIME_SET_END_SECONDS, CMF_TIME_INIT
              IMPLICIT NONE
              LOGNAM = 6
              DT = 3600._JPRB
              YYYY0 = 1900_JPIM
              SYEAR=1900_JPIM; SMON=1_JPIM; SDAY=1_JPIM
              EYEAR=1900_JPIM; EMON=1_JPIM; EDAY=1_JPIM
              CALL CMF_TIME_SET_START_SECONDS(0_JPIM)
              CALL CMF_TIME_SET_END_SECONDS(86400_JPIM)
              IF (EHOUR /= 24_JPIM .OR. EMIN /= 0_JPIM) ERROR STOP 'END_SEC=86400 did not map to 24:00'
              CALL CMF_TIME_INIT
              IF (ISHHMM /= 0_JPIM .OR. IEHHMM /= 2400_JPIM) ERROR STOP 'HHMM bounds wrong for END_SEC=86400'
              IF (KMINSTART /= 0_JPIM .OR. KMINEND /= 1440_JPIM) ERROR STOP 'minute bounds wrong for END_SEC=86400'
              IF (NSTEPS /= 24_JPIM) ERROR STOP 'NSTEPS wrong for one-day hourly window'
              PRINT *, 'time end 86400 ok', IEHHMM, KMINEND, NSTEPS
            END PROGRAM time_end_86400_probe
            """
        )
    )
    exe = build / "time_end_86400_probe"
    subprocess.run(
        [gfortran, *flags, str(probe), *(str(o) for o in objects), *netcdf_libs, "-o", str(exe)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    result = subprocess.run([str(exe)], cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    assert "time end 86400 ok" in result.stdout.lower()
