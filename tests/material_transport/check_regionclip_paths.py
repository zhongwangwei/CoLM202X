#!/usr/bin/env python3
"""RegionClip contract guard; optionally verify a real source/clip NetCDF pair.

No-argument checks use only the standard library. With SOURCE OUTPUT, requires
numpy/netCDF4 and compares every source block variable against its element,
pixel, patch or PFT subset (including attributes and mesh header metadata).
The output file creation timestamp is intentionally allowed to change.
"""
from pathlib import Path
import sys


def source_check():
    root = Path(__file__).resolve().parents[2]
    text = (root / 'mksrfdata/MOD_RegionClip.F90').read_text()
    assert 'USE MOD_MT_Methane_Preprocessing' not in text
    caller = (root / 'mksrfdata/MKSRFDATA.F90').read_text()
    assert 'requires_lake_soilc, requires_spatial_ph)' in caller
    assert 'DEF_dir_landdata, .false., .false.)' in caller
    assert 'MAX(1985, (lc_year / 5) * 5)' in text
    assert "input_year = ''" in text  # one layout decision, not field fallback
    assert "input_year = lc_cyear" in text
    for line in text.splitlines():
        if 'file_out' in line and 'dir_landdata_out' in line:
            if any('/' + d + '/' in line for d in ('mesh','landelm','landpatch','landpft','landhru','soil','htop','topography','lakedepth','pctpft')):
                assert 'lc_cyear' in line, line
    assert "'elmpixels', pixelmask" in text
    assert "'pctshared'" in text
    assert 'nf90_copy_att' in text
    assert 'Missing required lake carbon' in text
    assert "'/lulcc/'" in text and 'N_land_classification' in text
    assert 'IF (nelm_out == 0)' in text
    assert 'DEF_Runoff_SCHEME == 0 .and. DEF_TOPMOD_method /= 0' in text
    assert 'all(patchtypes_igbp(patchtypes_source) /= 0)' in text
    assert 'IF (count(vecmask) == 0) RETURN' in text
    assert 'Missing pixelset landpft for PFT-bearing source block' in text
    print('PASS RegionClip year/layout/mask structural contracts')


def compare(source, output, source_year='2005', output_year='2005'):
    import json
    import numpy as np
    from netCDF4 import Dataset
    source, output = Path(source), Path(output)
    def read(path, var):
        with Dataset(path) as ds:
            return np.asarray(ds[var][:])
    def path_for(p):
        rel = p.relative_to(source)
        parts = list(rel.parts)
        if len(parts) > 2 and parts[1] == source_year:
            parts[1] = output_year
        elif len(parts) == 2 and parts[0] not in ('LAI',):
            parts.insert(1,output_year)
        return output.joinpath(*parts)
    checked = 0; patches = lakes = 0
    for mesh in (source/'mesh'/source_year).glob('mesh_*.nc'):
        clipped = path_for(mesh)
        if not clipped.exists():
            continue  # excluded whole block; checked via header counts below
        suffix = mesh.name[len('mesh'):]
        ei = read(mesh,'elmindex'); eo = read(clipped,'elmindex')
        em = np.isin(ei, eo)
        assert len(eo) > 0 and em.sum() == len(eo)
        npxl=read(mesh,'elmnpxl'); pm=np.repeat(em,npxl)
        masks={'element':em, 'pixel':pm}
        for pset in ('landelm','landpatch','landpft','landhru'):
            p=source/pset/source_year/(pset+suffix)
            if p.exists(): masks[pset]=np.isin(read(p,'eindex'),eo)
        masks.update(patch=masks['landpatch'],pft=masks.get('landpft'))
        for p in source.rglob('*'+suffix):
            if len(p.relative_to(source).parts)>2 and p.relative_to(source).parts[1] != source_year:
                continue
            q=path_for(p)
            assert q.exists(), f'missing clipped field {q}'
            with Dataset(p) as a, Dataset(q) as b:
                assert a.ncattrs()==b.ncattrs(), (p,a.ncattrs(),b.ncattrs())
                for attr in a.ncattrs():
                    if attr == "create_time":
                        continue  # newly created output keeps its own creation timestamp
                    np.testing.assert_array_equal(a.getncattr(attr),b.getncattr(attr))
                assert set(a.variables)==set(b.variables), (p,set(a.variables)-set(b.variables))
                for name,v in a.variables.items():
                    w=b[name]; axis=0  # Fortran vector last dimension -> NetCDF Python first
                    mask=masks.get(v.dimensions[axis])
                    assert mask is not None, (p,name,v.dimensions)
                    np.testing.assert_array_equal(w[:],np.asarray(v[:])[mask],err_msg=str(p)+':'+name)
                    assert v.ncattrs()==w.ncattrs(), (p,name,v.ncattrs(),w.ncattrs())
                    for attr in v.ncattrs():
                        np.testing.assert_array_equal(v.getncattr(attr),w.getncattr(attr))
                    checked += 1
        t=read(output/'landpatch'/output_year/('landpatch'+suffix),'settyp'); lake=t==17
        patches+=len(t); lakes+=int(lake.sum())
        cfile=output/'soil'/output_year/('lake_soilc_patches'+suffix)
        if cfile.exists():
            c=read(cfile,'lake_soilc_patches')
            assert c.shape[0]==len(t) and (c[~lake]==0).all()
            for j in range(c.shape[1]):
                var=f'OM_density_s_l{min(8,max(1,j))}_patches'
                om=read(output/'soil'/output_year/(var+suffix),var)
                np.testing.assert_array_equal(c[lake,j],580.*om[lake])
    for name in ('block.nc','pixel.nc'):
        assert (source/name).read_bytes()==(output/name).read_bytes()
    with Dataset(source/'mesh'/source_year/'mesh.nc') as a, Dataset(output/'mesh'/output_year/'mesh.nc') as b:
        assert set(a.variables)==set(b.variables)
        for name in a.variables:
            if name != 'nelm_blk': np.testing.assert_array_equal(a[name][:],b[name][:])
        assert int(np.asarray(b['nelm_blk'][:]).sum()) == sum(len(read(p,'elmindex')) for p in (output/'mesh'/output_year).glob('mesh_*.nc'))
    assert checked and patches and lakes
    result=dict(variables_exact=checked,patches=patches,lake_patches=lakes,metadata_exact=True,lake_OM_exact=True)
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':
    if len(sys.argv)==1: source_check()
    else: compare(*sys.argv[1:])
