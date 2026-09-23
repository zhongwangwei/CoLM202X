from pathlib import Path
import os
import re
import shutil
import subprocess
import textwrap

import pytest

from fortran_test_support import netcdf_fortran_flags


ROOT = Path(__file__).resolve().parents[1]
ROUTING_NC = Path(os.environ.get(
    "COLM_CAMA_ROUTING_NC",
    "/Volumes/Data01/Data/CoLMruntime/unitcatchment/grid_routing_data_15min.nc",
))
GFORTRAN = shutil.which("gfortran") or "/opt/homebrew/bin/gfortran"


def _flat(path: Path) -> str:
    return re.sub(r"\s+", " ", re.sub(r"&\s*", " ", path.read_text().lower()))


def test_legacy_cama_accepts_one_vector_routing_netcdf() -> None:
    input_vars = _flat(ROOT / 'extends' / 'CaMa' / 'src' / 'yos_cmf_input.F90')
    maps = _flat(ROOT / 'extends' / 'CaMa' / 'src' / 'cmf_ctrl_maps_mod.F90')

    assert 'croutingnc' in input_vars
    assert 'namelist/nmap/' in maps and 'croutingnc' in maps
    assert 'call read_routing_header_cdf' in maps
    assert 'call read_routing_map_cdf' in maps
    assert 'call read_routing_topo_cdf' in maps


def test_vector_routing_path_uses_precomputed_topology_and_all_static_fields() -> None:
    maps = _flat(ROOT / 'extends' / 'CaMa' / 'src' / 'cmf_ctrl_maps_mod.F90')

    for varname in (
        'seq_x', 'seq_y', 'seq_next',
        'topo_area', 'topo_elevation', 'topo_distance', 'topo_rivlen',
        'topo_rivwth', 'topo_rivhgt', 'topo_rivman', 'topo_fldhgt',
        'bifurcation_upst', 'bifurcation_down', 'bifurcation_distance',
        'bifurcation_elevation', 'bifurcation_width', 'bifurcation_manning',
    ):
        assert f"'{varname}'" in maps

    assert 'if( croutingnc/="none" )then #ifdef usecdf_cmf call read_routing_map_cdf' in maps
    assert 'if( croutingnc=="none" )then write(lognam,*)' in maps
    assert 'call calc_1d_seq' in maps


def test_vector_routing_path_reads_bundled_input_matrix_without_bin() -> None:
    forcing = _flat(ROOT / 'extends' / 'CaMa' / 'src' / 'cmf_ctrl_forcing_mod.F90')

    assert 'call cmf_inpmat_init_routing_cdf' in forcing
    assert "'inpmat_x'" in forcing
    assert "'inpmat_y'" in forcing
    assert "'inpmat_area'" in forcing


def test_default_cama_namelist_uses_main_namelist_routing_nc() -> None:
    nml = _flat(ROOT / 'run' / 'CaMa' / 'cama_flood.nml')

    assert 'croutingnc = "none"' in nml
    assert "def_unitcatchment_file" in nml
    assert 'cdiminfo = "none"' in nml
    assert 'cinpmat = "none"' in nml
    assert 'lwinfilt = .true.' in nml


@pytest.mark.skipif(not ROUTING_NC.exists(), reason='local routing dataset is unavailable')
def test_requested_routing_nc_has_closed_valid_global_topology() -> None:
    netcdf4 = pytest.importorskip('netCDF4')
    import numpy as np

    with netcdf4.Dataset(ROUTING_NC) as dataset:
        seq = dataset['seq'][:]
        next_seq = dataset['seq_next'][:]
        x = dataset['seq_x'][:]
        y = dataset['seq_y'][:]
        nseqriv = int(dataset.nseqriv)
        nseqmax = int(dataset.nseqmax)

        assert np.array_equal(seq, np.arange(1, nseqmax + 1))
        assert np.all(next_seq[:nseqriv] > seq[:nseqriv])
        assert np.all(next_seq[nseqriv:] < 0)
        assert np.all((x >= 1) & (x <= int(dataset.nx)))
        assert np.all((y >= 1) & (y <= int(dataset.ny)))
        assert len(np.unique(np.column_stack((x, y)), axis=0)) == nseqmax
        assert np.all(np.diff(dataset['topo_fldhgt'][:], axis=1) >= 0.0)
        assert np.all(dataset['topo_area'][:] > 0.0)
        assert np.all(dataset['topo_rivlen'][:] > 0.0)
        assert np.all(dataset['topo_rivwth'][:] > 0.0)
        assert np.all(dataset['topo_rivhgt'][:] > 0.0)

        area = dataset['topo_area'][:]
        rivlen = dataset['topo_rivlen'][:]
        rivwth = dataset['topo_rivwth'][:]
        rivhgt = dataset['topo_rivhgt'][:]
        fldhgt = dataset['topo_fldhgt'][:]
        expected_storage = rivlen * rivwth * rivhgt
        previous_height = np.zeros(nseqmax)
        width_increment = area / rivlen / fldhgt.shape[1]
        for level in range(fldhgt.shape[1]):
            expected_storage += rivlen * (
                rivwth + width_increment * (level + 0.5)
            ) * (fldhgt[:, level] - previous_height)
            assert np.allclose(
                dataset['topo_fldstomax'][:, level], expected_storage, rtol=1.0e-12
            )
            previous_height = fldhgt[:, level]


def test_default_cama_coupling_is_hourly() -> None:
    defaults = _flat(ROOT / "extends/CaMa/src/cmf_ctrl_nmlist_mod.F90")
    assert re.search(r"ifrq_inp\s*=\s*1\s*!", defaults)
    assert re.search(r"dt\s*=\s*60\s*\*\s*60\s*!", defaults)
    for path in (ROOT / "run/CaMa").glob("*.nml"):
        nml = _flat(path)
        assert re.search(r"ifrq_inp\s*=\s*1\s*!", nml), path
        assert re.search(r"dt\s*=\s*3600\s*!", nml), path


def _netcdf_flags() -> tuple[list[str], list[str]]:
    return netcdf_fortran_flags()


def _storage_values(area: list[float], rivlen: list[float], rivwth: list[float], rivhgt: list[float], fldhgt: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    rivstomax = [l * w * h for l, w, h in zip(rivlen, rivwth, rivhgt)]
    nlfp = len(fldhgt[0])
    fldstomax: list[list[float]] = []
    for iseq, base in enumerate(rivstomax):
        row: list[float] = []
        prev_h = 0.0
        prev_sto = base
        width_increment = area[iseq] / rivlen[iseq] / nlfp
        for ilev in range(nlfp):
            now = prev_sto + rivlen[iseq] * (
                rivwth[iseq] + width_increment * (ilev + 0.5)
            ) * (fldhgt[iseq][ilev] - prev_h)
            row.append(now)
            prev_sto = now
            prev_h = fldhgt[iseq][ilev]
        fldstomax.append(row)
    return rivstomax, fldstomax


def _write_synthetic_routing_nc(
    path: Path,
    *,
    bed_schema: bool,
    distance: bool = True,
    distance_name: str = "topo_distance",
    legacy_sequence: bool = True,
    unsorted_records: bool = False,
    include_storage: bool = False,
    corrupt_storage: bool = False,
    cycle: bool = False,
    include_inpmat: bool = False,
) -> None:
    netcdf4 = pytest.importorskip("netCDF4")
    nseq = 2 if (unsorted_records or cycle) else 1
    seq_dim = "nseqmax" if legacy_sequence else "nseq"
    seq_ids = [1]
    seq_x = [1]
    seq_y = [1]
    seq_next = [-9]
    lon = [0.5, 1.5] if nseq == 2 else [0.5]
    topo_area = [1000.0]
    topo_dist = [1234.0]
    rivlen = [100.0]
    rivwth = [10.0]
    rivhgt = [3.0]
    rivman = [0.03]
    fldhgt = [[1.0, 2.0]]
    elev = [10.0 if bed_schema else 20.0]
    if nseq == 2:
        seq_ids = [2, 1] if unsorted_records else [1, 2]
        seq_x = [2, 1] if unsorted_records else [1, 2]
        seq_y = [1, 1]
        seq_next = [-9, 1] if unsorted_records else [2, -9]
        if cycle:
            seq_ids = [1, 2]
            seq_x = [1, 2]
            seq_next = [2, 1]
        topo_area = [1000.0, 2000.0]
        topo_dist = [1234.0, 2345.0]
        rivlen = [100.0, 200.0]
        rivwth = [10.0, 20.0]
        rivhgt = [3.0, 4.0]
        rivman = [0.03, 0.04]
        fldhgt = [[1.0, 2.0], [1.5, 3.0]]
        elev = [10.0, 30.0] if bed_schema else [20.0, 40.0]
    rivstomax, fldstomax = _storage_values(topo_area, rivlen, rivwth, rivhgt, fldhgt)
    if corrupt_storage:
        fldstomax[0][0] += 1.0

    with netcdf4.Dataset(path, "w") as nc:
        nc.createDimension("nx", 2 if nseq == 2 else 1)
        nc.createDimension("ny", 1)
        nc.createDimension("nlfp", 2)
        nc.createDimension("inpn", 1)
        nc.createDimension(seq_dim, nseq)
        if legacy_sequence:
            nc.createDimension("upnmax", 1)
            nc.nseqriv = 1 if nseq == 2 else 0
            nc.nseqall = nseq
        nc.west = 0.0
        nc.east = 2.0 if nseq == 2 else 1.0
        nc.north = 1.0
        nc.south = 0.0
        nc.createVariable("lon", "f8", ("nx",))[:] = lon
        nc.createVariable("lat", "f8", ("ny",))[:] = [0.5]
        if legacy_sequence or unsorted_records:
            nc.createVariable("seq", "i4", (seq_dim,))[:] = seq_ids
        nc.createVariable("seq_x", "i4", (seq_dim,))[:] = seq_x
        nc.createVariable("seq_y", "i4", (seq_dim,))[:] = seq_y
        nc.createVariable("seq_next", "i4", (seq_dim,))[:] = seq_next
        if legacy_sequence:
            upn = [0] if nseq == 1 else [0, 1]
            upst = [[0]] if nseq == 1 else [[0, 1]]
            nc.createVariable("seq_upst", "i4", ("upnmax", seq_dim))[:, :] = upst
            nc.createVariable("seq_upn", "i4", (seq_dim,))[:] = upn
        nc.createVariable("topo_area", "f8", (seq_dim,))[:] = topo_area
        if distance:
            nc.createVariable(distance_name, "f8", (seq_dim,))[:] = topo_dist
        nc.createVariable("topo_rivlen", "f8", (seq_dim,))[:] = rivlen
        nc.createVariable("topo_rivwth", "f8", (seq_dim,))[:] = rivwth
        nc.createVariable("topo_rivhgt", "f8", (seq_dim,))[:] = rivhgt
        nc.createVariable("topo_rivman", "f8", (seq_dim,))[:] = rivman
        nc.createVariable("topo_fldhgt", "f8", (seq_dim, "nlfp"))[:, :] = fldhgt
        if include_storage:
            nc.createVariable("topo_rivstomax", "f8", (seq_dim,))[:] = rivstomax
            nc.createVariable("topo_fldstomax", "f8", (seq_dim, "nlfp"))[:, :] = fldstomax
        if bed_schema:
            nc.createVariable("topo_rivelv", "f8", (seq_dim,))[:] = elev
        else:
            nc.createVariable("topo_elevation", "f8", (seq_dim,))[:] = elev
        if include_inpmat:
            nc.createVariable("inpmat_x", "i4", ("inpn", seq_dim))[:, :] = [[1] * nseq]
            nc.createVariable("inpmat_y", "i4", ("inpn", seq_dim))[:, :] = [[1] * nseq]
            nc.createVariable("inpmat_area", "f8", ("inpn", seq_dim))[:, :] = [[1.0] * nseq]

def _compile_routing_nc_probe(tmp_path: Path) -> Path:
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
        subprocess.run(
            [GFORTRAN, *flags, "-c", str(src), "-J", str(build), "-o", str(obj)],
            check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        objects.append(obj)
    probe = build / "routing_nc_probe.F90"
    probe.write_text(textwrap.dedent("""
        PROGRAM routing_nc_probe
          USE PARKIND1, ONLY: JPIM, JPRB
          USE YOS_CMF_INPUT, ONLY: CSETFILE, LOGNAM, LLOGOUT, LFPLAIN, LPTHOUT, LSLOPEMOUTH, &
                                 & LGDWDLY, LSLPMIX, LMEANSL, IMIS
          USE YOS_CMF_MAP, ONLY: D2ELEVTN, D2RIVELV, D2NXTDST, D2FLDHGT, &
                                 I1SEQX, I1SEQY, I2VECTOR, NSEQMAX
          USE CMF_CTRL_MAPS_MOD, ONLY: CMF_MAPS_NMLIST, CMF_RIVMAP_INIT, CMF_TOPO_INIT
          IMPLICIT NONE
          CHARACTER(LEN=512) :: nml, arg
          REAL(KIND=JPRB) :: bank, bed, dst
          INTEGER :: i
          IF (COMMAND_ARGUMENT_COUNT() < 4) ERROR STOP 'usage: probe nml bank bed dst'
          CALL GET_COMMAND_ARGUMENT(1, nml)
          CALL GET_COMMAND_ARGUMENT(2, arg); READ(arg,*) bank
          CALL GET_COMMAND_ARGUMENT(3, arg); READ(arg,*) bed
          CALL GET_COMMAND_ARGUMENT(4, arg); READ(arg,*) dst
          CSETFILE = TRIM(nml)
          LOGNAM = 6
          LLOGOUT = .FALSE.
          LFPLAIN = .TRUE.
          LPTHOUT = .FALSE.
          LSLOPEMOUTH = .FALSE.
          LGDWDLY = .FALSE.
          LSLPMIX = .FALSE.
          LMEANSL = .FALSE.
          IMIS = -9999_JPIM
          CALL CMF_MAPS_NMLIST
          CALL CMF_RIVMAP_INIT
          CALL CMF_TOPO_INIT
          IF (ABS(D2ELEVTN(1,1)-bank) > 1.e-8_JPRB) ERROR STOP 'bank elevation mismatch'
          IF (ABS(D2RIVELV(1,1)-bed) > 1.e-8_JPRB) ERROR STOP 'river-bed elevation mismatch'
          IF (ABS(D2NXTDST(1,1)-dst) > 1.e-8_JPRB) ERROR STOP 'downstream distance mismatch'
          DO i=1,NSEQMAX
            IF(I2VECTOR(I1SEQX(i),I1SEQY(i))/=i) ERROR STOP 'inverse grid mapping stale'
          ENDDO
          IF (bank==34._JPRB) THEN
            IF (ABS(D2FLDHGT(1,1,1)-1.5_JPRB)>1.e-8_JPRB) ERROR STOP 'square flood profile transposed'
          ENDIF
          PRINT *, 'routing nc probe ok', D2ELEVTN(1,1), D2RIVELV(1,1), D2NXTDST(1,1)
        END PROGRAM routing_nc_probe
    """))
    exe = build / "routing_nc_probe"
    subprocess.run(
        [GFORTRAN, *flags, str(probe), *(str(o) for o in objects), *libs, "-o", str(exe)],
        check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    return exe


def test_synthetic_routing_nc_supports_bank_and_bed_elevation_schemas(tmp_path: Path) -> None:
    exe = _compile_routing_nc_probe(tmp_path)
    for name, bed_schema, bank, bed, distance_name in (
        ("bank.nc", False, 20.0, 17.0, "topo_distance"),
        ("bed.nc", True, 13.0, 10.0, "topo_nxtdst"),
    ):
        nc = tmp_path / name
        nml = tmp_path / f"{name}.nml"
        _write_synthetic_routing_nc(nc, bed_schema=bed_schema, distance_name=distance_name)
        nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
        result = subprocess.run(
            [str(exe), str(nml), str(bank), str(bed), "1234.0"],
            cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
        )
        assert result.returncode == 0, result.stdout
        assert "routing nc probe ok" in result.stdout.lower()


def test_synthetic_routing_nc_rejects_missing_downstream_distance(tmp_path: Path) -> None:
    exe = _compile_routing_nc_probe(tmp_path)
    nc = tmp_path / "missing_distance.nc"
    nml = tmp_path / "missing_distance.nml"
    _write_synthetic_routing_nc(nc, bed_schema=True, distance=False)
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run(
        [str(exe), str(nml), "13.0", "10.0", "1234.0"],
        cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
    )
    assert result.returncode != 0
    assert "topo_distance/topo_nxtdst not present" in result.stdout.lower()
    assert "cannot derive downstream distance from topo_rivlen" in result.stdout.lower()



def test_canonical_routing_nc_works_without_sequence_attrs_or_upstream_vars(tmp_path: Path) -> None:
    exe = _compile_routing_nc_probe(tmp_path)
    nc = tmp_path / "canonical.nc"
    nml = tmp_path / "canonical.nml"
    _write_synthetic_routing_nc(
        nc,
        bed_schema=True,
        distance_name="topo_nxtdst",
        legacy_sequence=False,
        include_storage=True,
    )
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run(
        [str(exe), str(nml), "13.0", "10.0", "1234.0"],
        cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
    )
    assert result.returncode == 0, result.stdout
    assert "routing nc probe ok" in result.stdout.lower()


def test_canonical_routing_nc_reorders_unsorted_seq_records(tmp_path: Path) -> None:
    exe = _compile_routing_nc_probe(tmp_path)
    nc = tmp_path / "unsorted.nc"
    nml = tmp_path / "unsorted.nml"
    _write_synthetic_routing_nc(
        nc,
        bed_schema=True,
        legacy_sequence=False,
        unsorted_records=True,
    )
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run(
        [str(exe), str(nml), "34.0", "30.0", "2345.0"],
        cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
    )
    assert result.returncode == 0, result.stdout
    assert "routing nc probe ok" in result.stdout.lower()


def test_synthetic_routing_nc_rejects_cycle(tmp_path: Path) -> None:
    exe = _compile_routing_nc_probe(tmp_path)
    nc = tmp_path / "cycle.nc"
    nml = tmp_path / "cycle.nml"
    _write_synthetic_routing_nc(nc, bed_schema=True, legacy_sequence=False, cycle=True)
    nml.write_text(f'&NMAP CROUTINGNC="{nc}", LMAPCDF=.FALSE. /\n')
    result = subprocess.run(
        [str(exe), str(nml), "13.0", "10.0", "1234.0"],
        cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
    )
    assert result.returncode != 0
    assert "cycle" in result.stdout.lower()


def _compile_inpmat_probe(tmp_path: Path) -> Path:
    if not Path(GFORTRAN).exists():
        pytest.skip("gfortran is not available")
    build = tmp_path / "inpmat_build"
    build.mkdir()
    includes, libs = _netcdf_flags()
    flags = ["-O0", "-g", "-cpp", "-free", "-fimplicit-none", "-ffree-line-length-none", "-DUseCDF_CMF", f"-I{build}", *includes]
    sources = [
        "share/MOD_Precision.F90",
        "extends/CaMa/src/parkind1.F90",
        "extends/CaMa/src/yos_cmf_input.F90",
        "extends/CaMa/src/yos_cmf_time.F90",
        "extends/CaMa/src/yos_cmf_map.F90",
        "extends/CaMa/src/yos_cmf_prog.F90",
        "extends/CaMa/src/yos_cmf_diag.F90",
        "extends/CaMa/src/cmf_utils_mod.F90",
        "extends/CaMa/src/cmf_ctrl_maps_mod.F90",
        "extends/CaMa/src/cmf_ctrl_forcing_mod.F90",
    ]
    objects: list[Path] = []
    for source in sources:
        src = ROOT / source
        obj = build / f"{src.stem}.o"
        subprocess.run([GFORTRAN, *flags, "-c", str(src), "-J", str(build), "-o", str(obj)],
                       check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        objects.append(obj)
    probe = build / "inpmat_probe.F90"
    probe.write_text(textwrap.dedent("""
        PROGRAM inpmat_probe
          USE PARKIND1, ONLY: JPIM, JPRB
          USE YOS_CMF_INPUT, ONLY: CROUTINGNC, LOGNAM, NX, NY, INPN
          USE YOS_CMF_MAP, ONLY: NSEQMAX, INPX, INPY, INPA, ROUTING_OLD_TO_NEW
          USE CMF_CTRL_FORCING_MOD, ONLY: CMF_FORCING_INIT, LINTERP, LINPCDF
          IMPLICIT NONE
          INTEGER :: r1,r2
          CHARACTER(LEN=20) :: mode, shape
          CALL GET_COMMAND_ARGUMENT(1, CROUTINGNC)
          CALL GET_COMMAND_ARGUMENT(2, mode)
          CALL GET_COMMAND_ARGUMENT(3, shape)
          LOGNAM=6
          NX=2; NY=1; INPN=3; NSEQMAX=2
          IF(shape=='square') NSEQMAX=3
          r1=1; r2=2
          IF(mode=='reorder') THEN
            ALLOCATE(ROUTING_OLD_TO_NEW(NSEQMAX))
            ROUTING_OLD_TO_NEW(1:2)=[2,1]
            IF(NSEQMAX==3) ROUTING_OLD_TO_NEW(3)=3
            r1=2; r2=1
          ENDIF
          LINTERP=.TRUE.
          LINPCDF=.FALSE.
          CALL CMF_FORCING_INIT
          IF (INPX(r1,1)/=1 .OR. INPX(r1,2)/=2 .OR. INPX(r1,3)/=0 .OR. &
              INPX(r2,1)/=2 .OR. INPX(r2,2)/=1 .OR. INPX(r2,3)/=1) &
              ERROR STOP 'inpmat x order mismatch'
          IF (INPY(r1,1)/=1 .OR. INPY(r1,2)/=1 .OR. INPY(r1,3)/=0 .OR. &
              INPY(r2,1)/=1 .OR. INPY(r2,2)/=1 .OR. INPY(r2,3)/=1) &
              ERROR STOP 'inpmat y mismatch'
          IF (ABS(INPA(r1,1)-1._JPRB)>1.e-12_JPRB .OR. ABS(INPA(r1,2)-2._JPRB)>1.e-12_JPRB .OR. &
              ABS(INPA(r1,3)-0._JPRB)>1.e-12_JPRB .OR. ABS(INPA(r2,1)-3._JPRB)>1.e-12_JPRB .OR. &
              ABS(INPA(r2,2)-4._JPRB)>1.e-12_JPRB .OR. ABS(INPA(r2,3)-5._JPRB)>1.e-12_JPRB) &
              ERROR STOP 'inpmat area order mismatch'
          PRINT *, 'INPMAT_OK'
        END PROGRAM inpmat_probe
    """))
    exe = build / "inpmat_probe"
    subprocess.run([GFORTRAN, *flags, str(probe), *(str(o) for o in objects), *libs, "-o", str(exe)],
                   check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return exe


def _write_inpmat_nc(path: Path, *, seq_major: bool, square: bool = False) -> None:
    netcdf4 = pytest.importorskip("netCDF4")
    import numpy as np
    nseq = 3 if square else 2
    with netcdf4.Dataset(path, "w") as nc:
        nc.createDimension("nseq", nseq)
        nc.createDimension("inpn", 3)
        nc.createVariable("seq_x", "i4", ("nseq",))[:] = range(1, nseq + 1)
        x = [[1, 2, 0], [2, 1, 1]]
        y = [[1, 1, 1], [1, 1, 1]]
        a = [[1.0, 2.0, 0.0], [3.0, 4.0, 5.0]]
        if square:
            x.append([2, 2, 1]); y.append([1, 1, 1]); a.append([6., 7., 8.])
        dims = ("nseq", "inpn") if seq_major else ("inpn", "nseq")
        for name, dtype, values in (("x", "i4", x), ("y", "i4", y), ("area", "f8", a)):
            values = np.asarray(values)
            nc.createVariable("inpmat_" + name, dtype, dims)[:] = values if seq_major else values.T


def test_input_matrix_reader_preserves_dim_order_and_permutation(tmp_path: Path) -> None:
    exe = _compile_inpmat_probe(tmp_path)
    for square in (False, True):
        for seq_major in (True, False):
            nc = tmp_path / f"matrix_{square}_{seq_major}.nc"
            _write_inpmat_nc(nc, seq_major=seq_major, square=square)
            for mode in ("identity", "reorder"):
                result = subprocess.run([str(exe), str(nc), mode, "square" if square else "rectangle"],
                                        cwd=tmp_path, text=True, capture_output=True, timeout=60)
                assert result.returncode == 0, result.stdout + result.stderr
                assert "INPMAT_OK" in result.stdout


def test_float32_channel_geometry_accepts_roundoff_but_rejects_inconsistency(tmp_path):
    import numpy as np
    import netCDF4
    exe = _compile_routing_nc_probe(tmp_path)
    nc = tmp_path / "rounded.nc"
    nml = tmp_path / "rounded.nml"
    _write_synthetic_routing_nc(nc, bed_schema=True, legacy_sequence=False, include_storage=True)
    with netCDF4.Dataset(nc, "a") as dataset:
        dataset["topo_rivlen"][:] = np.float32(100.12345)
        dataset["topo_rivstomax"][:] = np.float32(100.12345 * 10 * 3)
    nml.write_text(f'&NMAP CROUTINGNC="{nc}" /\n')
    command = [str(exe), str(nml), "13.0", "10.0", "1234.0"]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    with netCDF4.Dataset(nc, "a") as dataset:
        dataset["topo_rivstomax"][:] *= 1.1
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0
    assert "inconsistent with CaMa river geometry" in result.stdout
