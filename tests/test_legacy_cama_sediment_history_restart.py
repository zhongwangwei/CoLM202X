"""Sediment history survives the NC checkpoint without changing sediment physics."""

from pathlib import Path
import subprocess

from test_legacy_cama_runtime import _compile_runtime_probe


def test_sediment_history_is_wired_only_to_restart_and_partial_flush() -> None:
    source = (Path(__file__).resolve().parents[1] / "extends/CaMa/src/cmf_ctrl_sed_mod.F90").read_text().lower()
    history = (Path(__file__).resolve().parents[1] / "extends/CaMa/src/MOD_CaMa_Vars.F90").read_text().lower()
    assert "if (lrestart .and. lrestcdf) call cmf_sed_history_read_cdf" in source
    assert "call cmf_sed_history_write_cdf" in source
    for field in ("history_sedout", "history_sedinp", "history_bedout", "history_netflw"):
        assert field in source
    assert "raw_sed_history=d2sedv_avg" in history
    assert "d2sedv_avg=raw_sed_history" in history


def test_sediment_history_nc_roundtrip_and_missing_guard(tmp_path: Path) -> None:
    exe = _compile_runtime_probe(tmp_path, r"""
        program sediment_history_probe
          use netcdf
          use parkind1, only: JPRB, JPIM
          use yos_cmf_input, only: NX, NY, LOGNAM, CSUFCDF
          use yos_cmf_map, only: NSEQMAX, I1SEQX, I1SEQY, REGIONTHIS
          use yos_cmf_time, only: JYYYYMMDD, JHOUR, KMINSTART
          use yos_cmf_diag, only: NADD_out
          use cmf_ctrl_restart_mod, only: CRESTDIR, CVNREST, CRESTSTO, LRESTCDF, LLEGACY_DAILY_HISTORY
          use cmf_ctrl_sed_mod, only: nsed, d2sedv_avg, CMF_SED_HISTORY_WRITE_CDF, CMF_SED_HISTORY_READ_CDF
          implicit none
          integer :: ncid, lonid, latid, timeid, status, i, j, k
          real(JPRB) :: expected(2,2,4)
          character(len=20) :: mode
          call get_command_argument(1,mode)
          LOGNAM=6; NX=2; NY=2; NSEQMAX=2; nsed=2; REGIONTHIS=1
          allocate(I1SEQX(2),I1SEQY(2),d2sedv_avg(2,2,4))
          I1SEQX=[1,2]; I1SEQY=[1,2]
          CRESTDIR='./'; CVNREST='sed_probe'; CSUFCDF='.nc'
          CRESTSTO='sed_probe2000010100.nc'; LRESTCDF=.true.
          JYYYYMMDD=20000101; JHOUR=0
          status=nf90_create(trim(CRESTSTO),NF90_NETCDF4,ncid)
          if(status/=nf90_noerr) error stop 'create failed'
          status=nf90_def_dim(ncid,'lon',NX,lonid)
          status=nf90_def_dim(ncid,'lat',NY,latid)
          status=nf90_def_dim(ncid,'time',NF90_UNLIMITED,timeid)
          status=nf90_close(ncid)
          if(status/=nf90_noerr) error stop 'close failed'
          do k=1,4
            do j=1,2
              do i=1,2
                expected(i,j,k)=real(100*k+10*j+i,JPRB)
              enddo
            enddo
          enddo
          if(mode=='missing'.or.mode=='missing_midday_zero'.or.mode=='missing_midnight') then
            NADD_out=3600._JPRB; LLEGACY_DAILY_HISTORY=.false.; KMINSTART=720
            if(mode=='missing_midday_zero') then
              NADD_out=0._JPRB; LLEGACY_DAILY_HISTORY=.true.
            endif
            if(mode=='missing_midnight') then
              NADD_out=0._JPRB; LLEGACY_DAILY_HISTORY=.true.; KMINSTART=1440
            endif
            call CMF_SED_HISTORY_READ_CDF
            if(mode=='missing_midnight') then
              print *, 'legacy midnight sediment history accepted'
              stop
            endif
            error stop 'missing sediment history was accepted'
          endif
          d2sedv_avg=expected
          call CMF_SED_HISTORY_WRITE_CDF
          d2sedv_avg=0._JPRB
          call CMF_SED_HISTORY_READ_CDF
          if(any(d2sedv_avg/=expected)) error stop 'sediment history changed on restart'
          print *, 'sediment history roundtrip ok'
        end program
    """)
    good = subprocess.run([str(exe)], cwd=tmp_path, text=True, capture_output=True)
    assert good.returncode == 0, good.stdout + good.stderr
    assert "sediment history roundtrip ok" in good.stdout
    for mode in ("missing", "missing_midday_zero"):
        rejected = subprocess.run([str(exe), mode], cwd=tmp_path, text=True, capture_output=True)
        assert rejected.returncode != 0
        assert "checkpoint lacks sediment daily history" in rejected.stdout
    accepted = subprocess.run([str(exe), "missing_midnight"], cwd=tmp_path, text=True, capture_output=True)
    assert accepted.returncode == 0, accepted.stdout + accepted.stderr
    assert "legacy midnight sediment history accepted" in accepted.stdout
