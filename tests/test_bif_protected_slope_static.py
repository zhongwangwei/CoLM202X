from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIF_SOURCE = ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90"


def test_bif_slope_uses_same_effective_surface_as_layer_depth() -> None:
    source = BIF_SOURCE.read_text()

    assert "zsurf_up = wdsrf_up_eff + rivelv_up" in source
    assert "zsurf_dn = wdsrf_dn_eff + rivelv_dn" in source
    assert "zsurf_up = wdsrf_up + rivelv_up" not in source
    assert "zsurf_dn = wdsrf_dn + rivelv_dn" not in source


def test_bif_effective_surface_keeps_channel_and_no_levee_fallbacks() -> None:
    source = BIF_SOURCE.read_text()

    assert "IF (ilev == 1) THEN\n               wdsrf_up_eff = wdsrf_up" in source
    assert "IF (ilev == 1) THEN\n               wdsrf_dn_eff = wdsrf_dn" in source
    assert "IF (DEF_USE_LEVEE .and. upstream_has_levee) THEN" in source
    assert "wdsrf_up_eff = protected_wdsrf_ucat(i_up)" in source
    assert "wdsrf_dn_eff = protected_wdsrf_ucat(i_dn)" in source
    assert "wdsrf_dn_eff = protected_wdsrf_dn_pth(ipth)" in source
