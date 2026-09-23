from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
MAPS = ROOT / "extends" / "CaMa" / "src" / "cmf_ctrl_maps_mod.F90"
FORCING = ROOT / "extends" / "CaMa" / "src" / "cmf_ctrl_forcing_mod.F90"


def _flat(path: Path) -> str:
    return re.sub(r"\s+", " ", re.sub(r"&\s*", " ", path.read_text().lower()))


def test_bundled_nc_rejects_unsupported_parallel_and_optional_physics() -> None:
    maps = _flat(MAPS)

    assert "croutingnc does not support usempi_cmf" in maps
    for flag, field in (
        ("lslpmix", "topo_mask_slope"),
        ("lslopemouth", "topo_elevslope"),
        ("lgdwdly", "topo_gdwdly"),
        ("lmeansl", "topo_meansl"),
    ):
        assert flag in maps
        assert field in maps
        assert "not present in routing netcdf" in maps

    assert "setting to zero" not in maps
    assert "nf90_inq_varid(ncid,'mask_slope',varid)" in maps
    assert "call ncerror(nf90_open(trim(crivparnc),nf90_nowrite,ncid)" in maps


def test_bundled_nc_validates_finite_topology_and_matrix_indices() -> None:
    maps = _flat(MAPS)
    forcing = _flat(FORCING)

    assert "ieee_is_finite" in maps
    assert "non-finite routing topography parameter" in maps
    assert "nupcalc" in maps and "routing sequence contains a cycle" in maps

    assert "ieee_is_finite" in forcing
    assert "non-finite input-matrix area" in forcing
    assert "where( inpa<=0._jprb )" in forcing
    assert "inpx=0" in forcing and "inpy=0" in forcing


def test_bundled_nc_has_safe_schema_fallbacks_only() -> None:
    maps = _flat(MAPS)

    assert "'topo_rivelv'" in maps
    assert "d2elevtn(:,1)=rivelv(:)+d2rivhgt(:,1)" in maps
    assert "'topo_nxtdst'" in maps
    assert "cannot derive downstream distance from topo_rivlen safely" in maps
