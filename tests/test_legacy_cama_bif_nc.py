from pathlib import Path
import subprocess
import textwrap

import pytest

from test_legacy_cama_routing_nc import ROOT, GFORTRAN, _netcdf_flags, _write_synthetic_routing_nc


def _append_bif_nc(path: Path, *, level_major: bool, down: int = 1) -> None:
    netcdf4 = pytest.importorskip("netCDF4")
    with netcdf4.Dataset(path, "a") as nc:
        nc.createDimension("path", 1)
        nc.createDimension("lev", 3)
        nc.createVariable("bifurcation_upst", "i4", ("path",))[:] = [2]
        nc.createVariable("bifurcation_down", "i4", ("path",))[:] = [down]
        nc.createVariable("bifurcation_distance", "f8", ("path",))[:] = [500.0]
        if level_major:
            dims = ("lev", "path")
            elevation = [[1.0], [2.0], [3.0]]
            width = [[10.0], [20.0], [30.0]]
        else:
            dims = ("path", "lev")
            elevation = [[1.0, 2.0, 3.0]]
            width = [[10.0, 20.0, 30.0]]
        nc.createVariable("bifurcation_elevation", "f8", dims)[:, :] = elevation
        nc.createVariable("bifurcation_width", "f8", dims)[:, :] = width
        nc.createVariable("bifurcation_manning", "f8", ("lev",))[:] = [0.03, 0.04, 0.05]


def _compile_bif_probe(tmp_path: Path) -> Path:
    if not Path(GFORTRAN).exists():
        pytest.skip("gfortran is not available")
    build = tmp_path / "build"
    build.mkdir()
    includes, libs = _netcdf_flags()
    flags = [
        "-O0", "-g", "-cpp", "-free", "-fimplicit-none", "-ffree-line-length-none",
        "-DUseCDF", "-DUseCDF_CMF", f"-I{build}", *includes,
    ]
    sources = [
        "share/MOD_Precision.F90",
        "extends/CaMa/src/parkind1.F90",
        "extends/CaMa/src/yos_cmf_input.F90",
        "extends/CaMa/src/yos_cmf_time.F90",
        "extends/CaMa/src/yos_cmf_map.F90",
        "extends/CaMa/src/cmf_utils_mod.F90",
        "extends/CaMa/src/cmf_ctrl_maps_mod.F90",
    ]
    objects: list[Path] = []
    for source in sources:
        src = ROOT / source
        obj = build / f"{src.stem}.o"
        subprocess.run([GFORTRAN, *flags, "-c", str(src), "-J", str(build), "-o", str(obj)],
                       check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        objects.append(obj)
    probe = build / "bif_probe.F90"
    probe.write_text(textwrap.dedent("""
        PROGRAM bif_probe
          USE PARKIND1, ONLY: JPIM, JPRB
          USE YOS_CMF_INPUT, ONLY: CSETFILE, LOGNAM, LLOGOUT, LFPLAIN, LPTHOUT, LSLOPEMOUTH, &
                                 & LGDWDLY, LSLPMIX, LMEANSL, IMIS
          USE YOS_CMF_MAP, ONLY: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_WTH, PTH_ELV, PTH_MAN
          USE CMF_CTRL_MAPS_MOD, ONLY: CMF_MAPS_NMLIST, CMF_RIVMAP_INIT
          IMPLICIT NONE
          CHARACTER(LEN=512) :: nml, mode
          CALL GET_COMMAND_ARGUMENT(1, nml)
          mode = ''
          IF (COMMAND_ARGUMENT_COUNT()>=2) CALL GET_COMMAND_ARGUMENT(2, mode)
          CSETFILE = TRIM(nml)
          LOGNAM = 6
          LLOGOUT = .FALSE.
          LFPLAIN = .TRUE.
          LPTHOUT = .TRUE.
          LSLOPEMOUTH = .FALSE.
          LGDWDLY = .FALSE.
          LSLPMIX = .FALSE.
          LMEANSL = .FALSE.
          IMIS = -9999_JPIM
          CALL CMF_MAPS_NMLIST
          CALL CMF_RIVMAP_INIT
          IF (TRIM(mode)=='square') THEN
            IF (NPTHOUT/=2 .OR. NPTHLEV/=2) ERROR STOP 'square bif dimensions mismatch'
            IF (PTH_UPST(1)/=1 .OR. PTH_DOWN(1)/=2 .OR. PTH_UPST(2)/=2 .OR. PTH_DOWN(2)/=1) &
                ERROR STOP 'square bif endpoint remap mismatch'
            IF (ABS(PTH_ELV(1,2)-12._JPRB)>1.e-12_JPRB .OR. ABS(PTH_ELV(2,1)-21._JPRB)>1.e-12_JPRB) &
                ERROR STOP 'square bif elevation order mismatch'
            IF (ABS(PTH_WTH(1,2)-120._JPRB)>1.e-12_JPRB .OR. ABS(PTH_WTH(2,1)-210._JPRB)>1.e-12_JPRB) &
                ERROR STOP 'square bif width order mismatch'
            IF (ABS(PTH_MAN(2)-0.04_JPRB)>1.e-12_JPRB) ERROR STOP 'square bif manning mismatch'
            PRINT *, 'BIF_NC_SQUARE_OK'
          ELSE
            IF (NPTHOUT/=1 .OR. NPTHLEV/=3) ERROR STOP 'bif dimensions mismatch'
            IF (PTH_UPST(1)/=1 .OR. PTH_DOWN(1)/=2) ERROR STOP 'bif endpoint remap mismatch'
            IF (ABS(PTH_WTH(1,3)-30._JPRB)>1.e-12_JPRB) ERROR STOP 'bif width order mismatch'
            IF (ABS(PTH_ELV(1,2)-2._JPRB)>1.e-12_JPRB) ERROR STOP 'bif elevation order mismatch'
            IF (ABS(PTH_MAN(3)-0.05_JPRB)>1.e-12_JPRB) ERROR STOP 'bif manning mismatch'
            PRINT *, 'BIF_NC_OK'
          ENDIF
        END PROGRAM bif_probe
    """))
    exe = build / "bif_probe"
    subprocess.run([GFORTRAN, *flags, str(probe), *(str(o) for o in objects), *libs, "-o", str(exe)],
                   check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return exe


@pytest.mark.parametrize("level_major", [False, True])
def test_bif_nc_infers_dims_and_remaps_unsorted_routing(tmp_path: Path, level_major: bool) -> None:
    exe = _compile_bif_probe(tmp_path)
    nc = tmp_path / ("bif_level_major.nc" if level_major else "bif_path_major.nc")
    nml = tmp_path / "bif.nml"
    _write_synthetic_routing_nc(nc, bed_schema=True, legacy_sequence=False, unsorted_records=True)
    _append_bif_nc(nc, level_major=level_major)
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run([str(exe), str(nml)], cwd=tmp_path, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    assert result.returncode == 0, result.stdout
    assert "BIF_NC_OK" in result.stdout


def test_bif_nc_rejects_external_down_zero(tmp_path: Path) -> None:
    exe = _compile_bif_probe(tmp_path)
    nc = tmp_path / "bif_down_zero.nc"
    nml = tmp_path / "bif_down_zero.nml"
    _write_synthetic_routing_nc(nc, bed_schema=True, legacy_sequence=False, unsorted_records=True)
    _append_bif_nc(nc, level_major=False, down=0)
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run([str(exe), str(nml)], cwd=tmp_path, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    assert result.returncode != 0
    assert "invalid bifurcation sequence index" in result.stdout.lower()


def _append_square_bif_nc(path: Path, *, level_major: bool) -> None:
    netcdf4 = pytest.importorskip("netCDF4")
    with netcdf4.Dataset(path, "a") as nc:
        nc.createDimension("path2", 2)
        nc.createDimension("lev2", 2)
        nc.createVariable("bifurcation_upst", "i4", ("path2",))[:] = [2, 1]
        nc.createVariable("bifurcation_down", "i4", ("path2",))[:] = [1, 2]
        nc.createVariable("bifurcation_distance", "f8", ("path2",))[:] = [500.0, 600.0]
        if level_major:
            dims = ("lev2", "path2")
            elevation = [[11.0, 21.0], [12.0, 22.0]]
            width = [[110.0, 210.0], [120.0, 220.0]]
        else:
            dims = ("path2", "lev2")
            elevation = [[11.0, 12.0], [21.0, 22.0]]
            width = [[110.0, 120.0], [210.0, 220.0]]
        nc.createVariable("bifurcation_elevation", "f8", dims)[:, :] = elevation
        nc.createVariable("bifurcation_width", "f8", dims)[:, :] = width
        nc.createVariable("bifurcation_manning", "f8", ("lev2",))[:] = [0.03, 0.04]


@pytest.mark.parametrize("level_major", [False, True])
def test_square_bif_dims_use_dimids_not_lengths(tmp_path: Path, level_major: bool) -> None:
    exe = _compile_bif_probe(tmp_path)
    nc = tmp_path / ("square_level_major.nc" if level_major else "square_path_major.nc")
    nml = tmp_path / "square_bif.nml"
    _write_synthetic_routing_nc(nc, bed_schema=True, legacy_sequence=False, unsorted_records=True)
    _append_square_bif_nc(nc, level_major=level_major)
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run([str(exe), str(nml), "square"], cwd=tmp_path, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    assert result.returncode == 0, result.stdout
    assert "BIF_NC_SQUARE_OK" in result.stdout
