from __future__ import annotations

import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from fortran_test_support import netcdf_fortran_flags

ROOT = Path(__file__).resolve().parents[1]
ROUTING_NC = Path(os.environ.get("COLM_CAMA_ROUTING_NC",
    "/Volumes/Data01/Data/CoLMruntime/unitcatchment/grid_routing_data_15min.nc"))
GFORTRAN = shutil.which("gfortran") or "/opt/homebrew/bin/gfortran"
def _netcdf_flags() -> tuple[list[str], list[str]]:
    return netcdf_fortran_flags()


def _compile_runtime_probe(tmp_path: Path, driver_source: str | None = None) -> Path:
    if not Path(GFORTRAN).exists():
        pytest.skip("gfortran is not available")
    build = tmp_path / "build"
    build.mkdir()
    includes, libs = _netcdf_flags()
    flags = [
        "-fopenmp",
        "-O0",
        "-g",
        "-Wall",
        "-cpp",
        "-free",
        "-fimplicit-none",
        "-fbounds-check",
        "-fbacktrace",
        "-ffree-line-length-none",
        "-DUseCDF",
        "-DUseCDF_CMF",
        f"-I{ROOT / 'include'}",
        f"-I{build}",
        *includes,
    ]
    sources = [
        "share/MOD_Precision.F90",
        "extends/CaMa/src/parkind1.F90",
        "extends/CaMa/src/yos_cmf_input.F90",
        "extends/CaMa/src/yos_cmf_time.F90",
        "extends/CaMa/src/yos_cmf_map.F90",
        "extends/CaMa/src/yos_cmf_prog.F90",
        "extends/CaMa/src/yos_cmf_diag.F90",
        "extends/CaMa/src/cmf_utils_mod.F90",
        "extends/CaMa/src/cmf_coupling_budget_mod.F90",
        "extends/CaMa/src/cmf_calc_outflw_mod.F90",
        "extends/CaMa/src/cmf_calc_pthout_mod.F90",
        "extends/CaMa/src/cmf_calc_fldstg_mod.F90",
        "extends/CaMa/src/cmf_calc_stonxt_mod.F90",
        "extends/CaMa/src/cmf_opt_outflw_mod.F90",
        "extends/CaMa/src/cmf_ctrl_mpi_mod.F90",
        "extends/CaMa/src/cmf_ctrl_damout_mod.F90",
        "extends/CaMa/src/cmf_ctrl_tracer_mod.F90",
        "extends/CaMa/src/cmf_ctrl_levee_mod.F90",
        "extends/CaMa/src/cmf_ctrl_forcing_mod.F90",
        "extends/CaMa/src/cmf_ctrl_boundary_mod.F90",
        "extends/CaMa/src/cmf_ctrl_output_mod.F90",
        "extends/CaMa/src/cmf_ctrl_restart_mod.F90",
        "extends/CaMa/src/cmf_ctrl_sed_mod.F90",
        "extends/CaMa/src/cmf_calc_diag_mod.F90",
        "extends/CaMa/src/cmf_ctrl_physics_mod.F90",
        "extends/CaMa/src/cmf_ctrl_time_mod.F90",
        "extends/CaMa/src/cmf_ctrl_maps_mod.F90",
        "extends/CaMa/src/cmf_ctrl_vars_mod.F90",
        "extends/CaMa/src/cmf_ctrl_nmlist_mod.F90",
        "extends/CaMa/src/cmf_drv_control_mod.F90",
        "extends/CaMa/src/cmf_drv_advance_mod.F90",
    ]
    objects: list[Path] = []
    for source in sources:
        src = ROOT / source
        obj = build / (src.stem + ".o")
        subprocess.run(
            [GFORTRAN, *flags, "-c", str(src), "-J", str(build), "-o", str(obj)],
            cwd=build,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        objects.append(obj)

    probe = build / "cama_runtime_probe.F90"
    probe.write_text(
        textwrap.dedent(
            r'''
            PROGRAM cama_runtime_probe
              USE PARKIND1, ONLY: JPIM, JPRB, JPRD
              USE YOS_CMF_INPUT, ONLY: CSETFILE, CLOGOUT, LLOGOUT, DT, LWEVAP, LWINFILT, CSUFBIN
              USE YOS_CMF_TIME, ONLY: ISHHMM, IHHMM, JHHMM, JHOUR, JYYYYMMDD, KMIN, KMINSTART, KMINEND, KMINNEXT, NSTEPS
              USE YOS_CMF_INPUT, ONLY: LWEVAP, LWINFILT
              USE YOS_CMF_MAP, ONLY: NSEQALL
              USE YOS_CMF_PROG, ONLY: P2RIVSTO, P2FLDSTO, D2RIVOUT, D2FLDOUT, D2RIVOUT_PRE, D2FLDOUT_PRE, &
                                    & D2RIVDPH_PRE, D2FLDSTO_PRE, D2RUNOFF, D2ROFSUB, D2GDWRTN, D2WEVAP, D2WINFILT
              USE YOS_CMF_DIAG, ONLY: D2RIVINF, D2FLDINF, D2PTHOUT, D2FLDFRC, D2FLDDPH, D2STORGE, &
                                    & D2RIVOUT_oAVG, D2FLDOUT_oAVG, D2WINFILTEX, D2WEVAPEX, NADD_out
              USE CMF_DRV_CONTROL_MOD, ONLY: CMF_DRV_INPUT, CMF_DRV_INIT
              USE CMF_DRV_ADVANCE_MOD, ONLY: CMF_DRV_ADVANCE
              USE CMF_CTRL_TIME_MOD, ONLY: CMF_TIME_SET_START_SECONDS, CMF_TIME_SET_END_SECONDS
              USE CMF_CALC_DIAG_MOD, ONLY: CMF_DIAG_RESET_OUTPUT, CMF_DIAG_GETAVE_OUTPUT
              USE CMF_CALC_STONXT_MOD, ONLY: CMF_CALC_STONXT
              USE CMF_CTRL_PHYSICS_MOD, ONLY: CMF_PHYSICS_FLDSTG
              USE CMF_CTRL_RESTART_MOD, ONLY: CMF_RESTART_WRITE, CMF_RESTART_INIT, CRESTDIR, CVNREST, CRESTSTO, LRESTCDF, LRESTDBL
              USE CMF_COUPLING_BUDGET_MOD, ONLY: sinks_prepaid, budget_init, budget_publish, budget_runoff, budget_debit
              IMPLICIT NONE
              CHARACTER(LEN=512) :: nml, rstdir, cdate
              INTEGER(KIND=JPIM) :: narg, nwet
              REAL(KIND=JPRD) :: before, after
              REAL(KIND=JPRB) :: old_dt
              INTEGER :: bx(2,2), by(2,2)
              REAL(KIND=JPRD) :: bw(2,2), barea(2,2), bvol(2), bfare(2), bdepth(2,2), bfrac(2,2)
              REAL(KIND=JPRD) :: bgrid(2,2), bflow(2), bstore(2), bevap(2,2), binfil(2,2), bevap_used(2), binfil_used(2)

              narg = COMMAND_ARGUMENT_COUNT()
              IF (narg < 2) ERROR STOP 'usage: cama_runtime_probe namelist restart_dir'
              CALL GET_COMMAND_ARGUMENT(1, nml)
              CALL GET_COMMAND_ARGUMENT(2, rstdir)
              CSETFILE = TRIM(nml)
              CLOGOUT = TRIM(rstdir)//'/log_CaMa_probe.txt'
              LLOGOUT = .FALSE.

              CALL CMF_DRV_INPUT
              CALL CMF_TIME_SET_START_SECONDS(1800_JPIM)
              CALL CMF_TIME_SET_END_SECONDS(4500_JPIM)
              CALL CMF_DRV_INIT

              IF (ISHHMM /= 30 .OR. IHHMM /= 30) ERROR STOP 'start minute was not preserved'
              IF (KMINSTART /= 30 .OR. KMINEND /= 75 .OR. NSTEPS /= 2) ERROR STOP 'minute clock bounds wrong'
              old_dt = DT
              CALL CMF_DRV_ADVANCE(1_JPIM)
              IF (KMIN /= 60 .OR. IHHMM /= 100) ERROR STOP '30min CMF_DRV_ADVANCE clock wrong'
              CALL CMF_DRV_ADVANCE(1_JPIM)
              IF (KMIN /= 75 .OR. IHHMM /= 115) ERROR STOP 'partial last-DT CMF_DRV_ADVANCE clock wrong'
              IF(DT /= old_dt) ERROR STOP 'last-step clamp changed configured DT'

              CALL CMF_DIAG_RESET_OUTPUT
              CALL CMF_DIAG_GETAVE_OUTPUT
              IF (NADD_out /= 0._JPRB) ERROR STOP 'zero-window diagnostic changed sample counter'
              IF (MAXVAL(ABS(D2RIVOUT_oAVG)) /= 0._JPRB) ERROR STOP 'zero-window rivout average not cleared'
              IF (MAXVAL(ABS(D2FLDOUT_oAVG)) /= 0._JPRB) ERROR STOP 'zero-window fldout average not cleared'

              old_dt = DT
              DT = 60._JPRB
              P2RIVSTO = 0._JPRD; P2FLDSTO = 0._JPRD
              D2RIVOUT = 0._JPRB; D2FLDOUT = 0._JPRB; D2RIVINF = 0._JPRB; D2FLDINF = 0._JPRB; D2PTHOUT = 0._JPRB
              D2RUNOFF = 0._JPRB; D2ROFSUB = 0._JPRB; D2GDWRTN = 0._JPRB; D2WEVAP = 0._JPRB; D2WINFILT = 0._JPRB
              D2WEVAPEX = 0._JPRB; D2WINFILTEX = 0._JPRB; sinks_prepaid = .FALSE.; LWEVAP = .TRUE.; LWINFILT = .TRUE.
              P2FLDSTO(1,1) = 1000._JPRD; D2WINFILT(1,1) = 2._JPRB; D2WEVAP(1,1) = 3._JPRB
              CALL CMF_CALC_STONXT
              after = SUM(P2RIVSTO) + SUM(P2FLDSTO)
              IF (ABS(after - 700._JPRD) > 1.e-6_JPRD) ERROR STOP 'sink conservation smoke failed'
              IF (ABS(D2WINFILTEX(1,1)-2._JPRB) > 1.e-6_JPRB) ERROR STOP 'infil extraction rate changed'
              IF (ABS(D2WEVAPEX(1,1)-3._JPRB) > 1.e-6_JPRB) ERROR STOP 'evap extraction rate changed'

              P2RIVSTO = 0._JPRD; P2FLDSTO = 0._JPRD
              D2RIVOUT = 0._JPRB; D2FLDOUT = 0._JPRB; D2RIVINF = 0._JPRB; D2FLDINF = 0._JPRB; D2PTHOUT = 0._JPRB
              D2RUNOFF = 0._JPRB; D2ROFSUB = 0._JPRB; D2GDWRTN = 0._JPRB; D2WEVAP = 0._JPRB; D2WINFILT = 0._JPRB
              D2WEVAPEX = 0._JPRB; D2WINFILTEX = 0._JPRB; sinks_prepaid = .TRUE.; LWEVAP = .TRUE.; LWINFILT = .TRUE.
              P2FLDSTO(1,1) = 700._JPRD; D2FLDOUT(1,1) = 1._JPRB
              D2WINFILT(1,1) = 2._JPRB; D2WEVAP(1,1) = 3._JPRB; D2WINFILTEX(1,1) = 2._JPRB; D2WEVAPEX(1,1) = 3._JPRB
              CALL CMF_CALC_STONXT
              after = SUM(P2RIVSTO) + SUM(P2FLDSTO)
              IF (ABS(after - 640._JPRD) > 1.e-6_JPRD) ERROR STOP 'prepaid sinks were debited twice or routing not applied'
              IF (ABS(D2WINFILTEX(1,1)-2._JPRB) > 1.e-6_JPRB) ERROR STOP 'prepaid infil diagnostic was not preserved'
              IF (ABS(D2WEVAPEX(1,1)-3._JPRB) > 1.e-6_JPRB) ERROR STOP 'prepaid evap diagnostic was not preserved'

              P2RIVSTO = 0._JPRD; P2FLDSTO = 0._JPRD
              D2RIVOUT = 0._JPRB; D2FLDOUT = 0._JPRB; D2RIVINF = 0._JPRB; D2FLDINF = 0._JPRB; D2PTHOUT = 0._JPRB
              D2RUNOFF = 0._JPRB; D2ROFSUB = 0._JPRB; D2GDWRTN = 0._JPRB; D2WEVAP = 0._JPRB; D2WINFILT = 0._JPRB
              D2WEVAPEX = 4._JPRB; D2WINFILTEX = 0._JPRB; sinks_prepaid = .TRUE.; LWEVAP = .TRUE.; LWINFILT = .FALSE.
              P2FLDSTO(1,1) = 760._JPRD; D2FLDOUT(1,1) = 1._JPRB; D2WEVAP(1,1) = 4._JPRB
              CALL CMF_CALC_STONXT
              after = SUM(P2RIVSTO) + SUM(P2FLDSTO)
              IF (ABS(after - 700._JPRD) > 1.e-6_JPRD) ERROR STOP 'prepaid evap-only branch double debited'
              IF (ABS(D2WEVAPEX(1,1)-4._JPRB) > 1.e-6_JPRB) ERROR STOP 'prepaid evap-only diagnostic changed'
              sinks_prepaid = .FALSE.; LWEVAP = .TRUE.; LWINFILT = .TRUE.; DT = old_dt

              bx = RESHAPE([1,2,1,2], SHAPE(bx)); by = RESHAPE([1,1,2,2], SHAPE(by))
              bw = RESHAPE([60._JPRD,40._JPRD,40._JPRD,60._JPRD], SHAPE(bw)); barea = 100._JPRD
              bvol = [100._JPRD,50._JPRD]; bfare = [20._JPRD,10._JPRD]
              CALL budget_init(bx,by,bw,barea)
              CALL budget_publish(bvol,bfare,bdepth,bfrac)
              IF (ABS(SUM(bdepth*barea)-SUM(bvol)) > 1.e-10_JPRD) ERROR STOP 'budget publish volume not conservative'
              bgrid = RESHAPE([6._JPRD,4._JPRD,2._JPRD,8._JPRD], SHAPE(bgrid))
              CALL budget_runoff(bgrid,bflow)
              IF (ABS(SUM(bflow)-SUM(bgrid)) > 1.e-10_JPRD) ERROR STOP 'budget runoff mapping not conservative'
              bstore = bvol; bevap = 0._JPRD; binfil = 0._JPRD; bevap(1,1)=5._JPRD; binfil(2,2)=7._JPRD
              CALL budget_debit(bevap,binfil,bstore,bevap_used,binfil_used)
              IF (ABS(SUM(bevap_used)+SUM(binfil_used)-12._JPRD) > 1.e-10_JPRD) ERROR STOP 'budget debit not conservative'
              IF (ABS(SUM(bstore)-(SUM(bvol)-12._JPRD)) > 1.e-10_JPRD) ERROR STOP 'budget storage debit mismatch'

              P2RIVSTO = 0._JPRD; P2FLDSTO = 0._JPRD
              D2RIVOUT_PRE = 0._JPRB; D2FLDOUT_PRE = 0._JPRB; D2RIVDPH_PRE = 0._JPRB; D2FLDSTO_PRE = 0._JPRB
              nwet = MIN(100_JPIM, NSEQALL)
              P2FLDSTO(1:nwet,1) = 1.0e9_JPRD
              CALL CMF_PHYSICS_FLDSTG
              IF (MAXVAL(D2FLDFRC(1:nwet,1)) <= 0._JPRB .OR. MAXVAL(D2FLDDPH(1:nwet,1)) <= 0._JPRB) &
                ERROR STOP 'flood stage not published before restart'
              before = SUM(P2FLDSTO)

              CRESTDIR = TRIM(rstdir)//'/'
              CVNREST = 'runtime_restart'
              LRESTCDF = .FALSE.
              LRESTDBL = .TRUE.
              CALL CMF_RESTART_WRITE
              WRITE(cdate,'(I8.8,I2.2)') JYYYYMMDD,JHOUR
              CRESTSTO = TRIM(CRESTDIR)//TRIM(CVNREST)//TRIM(cdate)//TRIM(CSUFBIN)

              IF (narg > 2) THEN
                ! Legacy binary has no daily-history payload; this 01:15 file
                ! must be rejected rather than silently losing that history.
                CALL CMF_RESTART_INIT
                ERROR STOP 'mid-day binary restart was accepted'
              ENDIF

              PRINT *, 'cama runtime probe ok', KMINSTART, KMINEND, NSTEPS, before
            END PROGRAM cama_runtime_probe
            '''
        )
    )
    if driver_source is not None:
        probe.write_text(textwrap.dedent(driver_source))
    exe = build / "cama_runtime_probe"
    subprocess.run(
        [GFORTRAN, *flags, str(probe), *(str(o) for o in objects), *libs, "-o", str(exe)],
        cwd=build,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return exe


def _write_probe_namelist(path: Path) -> None:
    path.write_text(
        textwrap.dedent(
            f'''
            &NRUNVER
            LADPSTP=.TRUE., LFPLAIN=.TRUE., LKINE=.FALSE., LFLDOUT=.TRUE., LPTHOUT=.FALSE., LDAMOUT=.FALSE., LLEVEE=.FALSE.,
            LROSPLIT=.FALSE., LWEVAP=.TRUE., LWEVAPFIX=.TRUE., LWINFILT=.TRUE., LWINFILTFIX=.FALSE., LWEXTRACTRIV=.FALSE.,
            LSLOPEMOUTH=.FALSE., LGDWDLY=.FALSE., LSLPMIX=.FALSE., LMEANSL=.FALSE., LSEALEV=.FALSE., LOUTINS=.FALSE.,
            LRESTART=.FALSE., LSTOONLY=.FALSE., LOUTPUT=.FALSE., LOUTINI=.FALSE., LGRIDMAP=.TRUE., LLEAPYR=.TRUE.,
            LMAPEND=.FALSE., LSTG_ES=.FALSE. /
            &NDIMTIME
            CDIMINFO="NONE", DT=1800, IFRQ_INP=1 /
            &NPARAM
            PMANRIV=0.03D0, PMANFLD=0.10D0, PGRV=9.8D0, PDSTMTH=10000.D0, PCADP=0.7, PMINSLP=1.D-5,
            IMIS=-9999, RMIS=1.E36, DMIS=1.E36, CSUFBIN='.bin', CSUFVEC='.vec', CSUFPTH='.pth', CSUFCDF='.nc' /
            &NSIMTIME
            SYEAR=1900, SMON=1, SDAY=1, SHOUR=0, SMIN=30, EYEAR=1900, EMON=1, EDAY=1, EHOUR=1, EMIN=30 /
            &NMAP
            CROUTINGNC="{ROUTING_NC}", LMAPCDF=.FALSE., CNEXTXY="NONE", CGRAREA="NONE", CELEVTN="NONE", CNXTDST="NONE",
            CRIVLEN="NONE", CFLDHGT="NONE", CRIVWTH="NONE", CRIVHGT="NONE", CRIVMAN="NONE", CPTHOUT="NONE",
            CGDWDLY="", CMEANSL="", CRIVCLINC="", CRIVPARNC="", CMEANSLNC="" /
            &NRESTART
            CRESTSTO="", CRESTDIR="./", CVNREST="restart", LRESTCDF=.FALSE., LRESTDBL=.TRUE., IFRQ_RST=0 /
            &NFORCE
            LINTERP=.TRUE., LINPEND=.FALSE., LINPDAY=.FALSE., LINPCDF=.FALSE., LITRPCDF=.FALSE., CINPMAT="NONE", DROFUNIT=1.0,
            CROFDIR="./runoff/", CROFPRE="Roff____", CROFSUF=".one", CSUBDIR="./runoff/", CSUBPRE="Rsub____", CSUBSUF=".one",
            CROFCDF="NONE", CVNTIME="time", CVNROF="runoff", CVNSUB="NONE", SYEARIN=0, SMONIN=0, SDAYIN=0, SHOURIN=0 /
            '''
        )
    )


@pytest.mark.skipif(not ROUTING_NC.exists(), reason="bundled routing netCDF is unavailable")
@pytest.mark.parametrize("threads", [1, 4])
def test_legacy_cama_core_runtime_smoke_restart_history_and_budget(tmp_path: Path, threads: int) -> None:
    exe = _compile_runtime_probe(tmp_path)
    nml = tmp_path / "cama_runtime_probe.nml"
    restart_dir = tmp_path / "restart"
    restart_dir.mkdir()
    _write_probe_namelist(nml)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    result = subprocess.run(
        [str(exe), str(nml), str(restart_dir)],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout
    assert "cama runtime probe ok" in result.stdout
    rejected = subprocess.run(
        [str(exe), str(nml), str(restart_dir), "reject_midday"],
        cwd=tmp_path,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=120,
    )
    assert rejected.returncode != 0
    assert "legacy binary restart lacks a safe daily history boundary" in rejected.stdout
