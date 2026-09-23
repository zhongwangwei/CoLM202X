"""Exercise real CaMa namelist/map/forcing readers, not a Python path resolver."""
from pathlib import Path
import subprocess
import textwrap

import pytest

from test_legacy_cama_routing_nc import (
    ROOT, GFORTRAN, _netcdf_flags, _write_synthetic_routing_nc,
)
from test_legacy_cama_runtime import ROUTING_NC


@pytest.fixture(scope="module")
def input_probe(tmp_path_factory):
    if not Path(GFORTRAN).exists():
        pytest.skip("gfortran unavailable")
    build = tmp_path_factory.mktemp("cama_input_config")
    includes, libs = _netcdf_flags()
    flags = ["-cpp", "-DUseCDF_CMF", "-ffree-line-length-none", "-fcheck=all",
             "-ffpe-trap=invalid,zero,overflow", f"-I{build}", *includes]
    sources = [ROOT / "share/MOD_Precision.F90"] + [ROOT / "extends/CaMa/src" / name for name in (
        "parkind1.F90", "yos_cmf_input.F90", "yos_cmf_time.F90", "yos_cmf_map.F90",
        "yos_cmf_prog.F90", "yos_cmf_diag.F90", "cmf_utils_mod.F90",
        "cmf_ctrl_nmlist_mod.F90", "cmf_ctrl_maps_mod.F90", "cmf_ctrl_forcing_mod.F90",
    )]
    probe = build / "probe.F90"
    probe.write_text(textwrap.dedent("""
        program input_probe
          use parkind1
          use yos_cmf_input
          use yos_cmf_map, only: INPX, INPA, D2RIVELV
          use cmf_ctrl_nmlist_mod, only: CMF_CONFIG_NMLIST
          use cmf_ctrl_maps_mod, only: CMF_MAPS_NMLIST, CMF_RIVMAP_INIT, CMF_TOPO_INIT
          use cmf_ctrl_forcing_mod, only: CMF_FORCING_NMLIST, CMF_FORCING_INIT
          implicit none
          character(len=1024) :: routing_file, mode
          LOGNAM=6
          call get_command_argument(1,CSETFILE)
          call get_command_argument(2,routing_file)
          call get_command_argument(3,mode)
          call CMF_CONFIG_NMLIST
          if (mode=='standalone') then
            call CMF_MAPS_NMLIST
          else
            call CMF_MAPS_NMLIST(trim(routing_file))
          endif
          if (CROUTINGNC=='NONE') then
            if (NX/=7 .or. NY/=3 .or. NXIN/=2) error stop 'legacy dims lost'
          else
            if (NX/=1 .or. NY/=1 .or. NXIN/=1) error stop 'NC dims lost'
            call CMF_RIVMAP_INIT
            call CMF_TOPO_INIT
            call CMF_FORCING_NMLIST
            call CMF_FORCING_INIT
            if (INPX(1,1)/=1 .or. abs(INPA(1,1)-1._JPRB)>1.e-12_JPRB) &
                error stop 'NC matrix not read'
            print *, 'BED=',D2RIVELV(1,1)
          endif
          if (.not.LFPLAIN .or. .not.LWEVAP .or. LWINFILT .or. LPTHOUT) &
              error stop 'physics flags changed'
          if (DT/=1800 .or. IFRQ_INP/=2) error stop 'timing overridden'
          print *, 'SELECTED=',trim(CROUTINGNC)
          print *, 'INPUT_CONFIG_OK'
        end program
    """))
    exe = build / "probe"
    result = subprocess.run([GFORTRAN, *flags, *(str(p) for p in sources), str(probe),
                             *libs, "-o", str(exe)], cwd=build, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    return exe


def routing_nc(path, *, bed_schema=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_synthetic_routing_nc(path, bed_schema=bed_schema, include_inpmat=True)
    return path


def run_probe(exe, directory, *, nc="NONE", primary="null", legacy=False,
              mode="coupled", mapcdf=True):
    dims = directory / "legacy_dims.txt"
    if legacy:
        dims.write_text("7\n3\n2\n2\n1\n1\n\n-180\n180\n90\n-90\n")
    nml = directory / "input.nml"
    nml.write_text(f"""
&NRUNVER LFPLAIN=.TRUE., LWEVAP=.TRUE., LWINFILT=.FALSE., LPTHOUT=.FALSE., LMAPEND=.TRUE. /
&NDIMTIME CDIMINFO='{dims}', DT=1800, IFRQ_INP=2 /
&NPARAM /
&NMAP CROUTINGNC='{nc}', LMAPCDF={'.TRUE.' if mapcdf else '.FALSE.'},
 CNEXTXY='missing.bin', CRIVCLINC='missing.nc', CRIVPARNC='missing.nc',
 CGRAREA='missing.bin', CELEVTN='missing.bin', CNXTDST='missing.bin',
 CRIVLEN='missing.bin', CRIVHGT='missing.bin', CRIVWTH='missing.bin',
 CRIVMAN='missing.bin', CFLDHGT='missing.bin', CPTHOUT='missing.bin' /
&NFORCE LINTERP=.TRUE., LITRPCDF=.TRUE., CINPMAT='missing.nc', LINPEND=.TRUE. /
""")
    return subprocess.run([str(exe), str(nml), str(primary), mode], cwd=directory,
                          text=True, capture_output=True, timeout=30)


@pytest.mark.parametrize("mapcdf", [False, True])
def test_nc_bypasses_invalid_legacy_files_but_keeps_physics(input_probe, tmp_path, mapcdf):
    nc = routing_nc(tmp_path / "routing.nc")
    result = run_probe(input_probe, tmp_path, nc=nc, mode="standalone", mapcdf=mapcdf)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "INPUT_CONFIG_OK" in result.stdout
    assert "ignored" in result.stdout.lower()


def test_primary_file_overrides_old_cama_path(input_probe, tmp_path):
    primary = routing_nc(tmp_path / "primary.nc", bed_schema=True)
    result = run_probe(input_probe, tmp_path, primary=primary, nc="obsolete_missing.nc")
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"SELECTED={primary}" in result.stdout
    assert "10.000" in result.stdout


def test_bad_explicit_path_does_not_fall_back(input_probe, tmp_path):
    good = routing_nc(tmp_path / "legacy.nc")
    missing = tmp_path / "absent.nc"
    result = run_probe(input_probe, tmp_path, nc=good, primary=missing)
    assert result.returncode != 0
    assert "Routing netCDF not found" in result.stdout
    assert str(missing) in result.stdout


def test_none_preserves_legacy_dimensions(input_probe, tmp_path):
    result = run_probe(input_probe, tmp_path, nc="NONE", legacy=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "SELECTED=NONE" in result.stdout


def test_colm_passes_existing_namelist_fields_to_cama():
    source = (ROOT / "extends/CaMa/src/MOD_CaMa_colmCaMa.F90").read_text()
    assert "CALL CMF_DRV_INPUT(DEF_UnitCatchment_file, DEF_CaMa_Restart_file)" in source
    driver = (ROOT / "extends/CaMa/src/cmf_drv_control_mod.F90").read_text()
    assert "CALL CMF_MAPS_NMLIST(ROUTING_FILE)" in driver
    namelist = (ROOT / "share/MOD_Namelist.F90").read_text()
    assert "DEF_CaMa_Restart_file," in namelist
    assert "CALL mpi_bcast (DEF_CaMa_Restart_file" in namelist
    assert "IF(LPTHOUT) CVARSOUT=TRIM(CVARSOUT)//',pthflw,pthout'" in source


@pytest.mark.parametrize("real_routing,bifurcation", [
    (False, False),
    pytest.param(True, False, marks=pytest.mark.skipif(not ROUTING_NC.exists(), reason="routing dataset unavailable")),
    pytest.param(True, True, marks=pytest.mark.skipif(not ROUTING_NC.exists(), reason="routing dataset unavailable")),
], ids=["synthetic", "real-routing", "real-bifurcation"])
def test_coupled_nc_initializes_and_steps_without_any_cama_namelist(tmp_path, real_routing, bifurcation):
    from test_legacy_cama_runtime import _compile_runtime_probe

    nc = ROUTING_NC if real_routing else routing_nc(tmp_path / "routing with spaces.nc", bed_schema=True)
    exe = _compile_runtime_probe(tmp_path, r"""
        program nc_only
          use, intrinsic :: ieee_arithmetic
          use parkind1
          use yos_cmf_input
          use yos_cmf_prog, only: D2RUNOFF, P2RIVSTO, P2FLDSTO
          use cmf_drv_control_mod, only: CMF_DRV_INPUT, CMF_DRV_INIT, CMF_DRV_END
          use cmf_drv_advance_mod, only: CMF_DRV_ADVANCE
          use cmf_ctrl_forcing_mod, only: LINTERP, LINPCDF, LITRPCDF
          use cmf_ctrl_maps_mod, only: LMAPCDF
          use cmf_ctrl_output_mod, only: NVARSOUT
          use cmf_ctrl_restart_mod, only: CMF_RESTART_WRITE, CRESTDIR, LRESTCDF
          use cmf_ctrl_time_mod, only: CMF_TIME_SET_START_SECONDS
          implicit none
          character(len=512) :: path, restart_path, arg
          real(JPRD) :: expected
          call get_command_argument(1,path)
          restart_path='null'
          if(command_argument_count()>1) call get_command_argument(2,restart_path)
          CSETFILE='does_not_exist.nml'
          call CMF_DRV_INPUT(trim(path),trim(restart_path))
          if (trim(restart_path)/='null') call CMF_TIME_SET_START_SECONDS(3600_JPIM)
          if(CSETFILE/='NONE') error stop 'second namelist used'
          if(.not.LWEVAP.or..not.LWINFILT.or..not.LFPLAIN.or..not.LADPSTP) error stop 'coupling defaults'
          if(.not.LINTERP.or.LINPCDF.or.LITRPCDF.or.LMAPCDF.or.LMAPEND) error stop 'format defaults'
          if(LDAMOUT.or.LLEVEE.or.LSEDIMENT.or.LTRACE) error stop 'unsupported physics enabled'
          if(DT/=3600.or.IFRQ_INP/=1) error stop 'not hourly'
          ! The CoLM coupler sets LPTHOUT from DEF_USE_BIFURCATION after CMF_DRV_INPUT.
          call get_command_argument(4,arg)
          LPTHOUT=arg=='bifurcation'
          call CMF_DRV_INIT
          if(NVARSOUT<=0) error stop 'no history diagnostics'
          if(.not.LRESTCDF) error stop 'restart format not NC'
          if(LRESTART) then
            call get_command_argument(3,arg); read(arg,*) expected
            ! Python and Fortran may reduce large global arrays in different orders.
            if(abs(sum(P2RIVSTO+P2FLDSTO)-expected)>1.e-10_JPRD*max(1._JPRD,expected)) &
                error stop 'restart water lost'
            call CMF_DRV_END
            print *, 'NC_ONLY_RESUME_OK'
            stop
          endif
          D2RUNOFF=1._JPRB
          call CMF_DRV_ADVANCE(1_JPIM)
          if(any(.not.ieee_is_finite(P2RIVSTO)).or.any(.not.ieee_is_finite(P2FLDSTO))) &
              error stop 'nonfinite storage'
          if(sum(P2RIVSTO+P2FLDSTO)<=0._JPRD) error stop 'runoff was not routed'
          CRESTDIR='./'
          call CMF_RESTART_WRITE
          call CMF_DRV_END
          print *, 'NC_ONLY_STEP_OK'
        end program
    """)
    mode = "bifurcation" if bifurcation else "ordinary"
    result = subprocess.run([str(exe), str(nc), "null", "0", mode], cwd=tmp_path,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "NC_ONLY_STEP_OK" in result.stdout
    assert not (tmp_path / "does_not_exist.nml").exists()

    import netCDF4
    checkpoint, = tmp_path.glob("restart*.nc")
    with netCDF4.Dataset(checkpoint) as dataset:
        expected = float(dataset["rivsto"][:].sum() + dataset["fldsto"][:].sum())
    resumed = subprocess.run([str(exe), str(nc), str(checkpoint), str(expected), mode], cwd=tmp_path,
                             capture_output=True, text=True, timeout=60)
    assert resumed.returncode == 0, resumed.stdout + resumed.stderr
    assert "NC_ONLY_RESUME_OK" in resumed.stdout
    missing = subprocess.run([str(exe), str(nc), str(tmp_path / "missing_restart.nc")], cwd=tmp_path,
                             capture_output=True, text=True, timeout=60)
    assert missing.returncode != 0
    assert "DEF_CaMa_Restart_file not found" in missing.stdout
