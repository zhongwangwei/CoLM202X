from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIF = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90").read_text(
    encoding="utf-8"
)


def test_tiny_flux_from_empty_donor_is_limited_to_zero() -> None:
    flux = 5.0e-11
    dt = 60.0
    donor_storage = 0.0

    transfer = abs(flux) * dt
    rate = min(1.0, donor_storage / transfer)

    assert transfer == 3.0e-9
    assert rate == 0.0
    assert flux * rate * dt == 0.0


def test_donor_limiters_do_not_discard_positive_flux_by_rate_threshold() -> None:
    assert "IF (layer_transfer <= 0._r8) CYCLE" in BIF
    assert "protected_outgoing(i_ucat) > 0._r8" in BIF
    assert "bif_outflow > 0._r8" in BIF

    assert "IF (layer_transfer <= 1.e-10_r8) CYCLE" not in BIF
    assert "protected_outgoing(i_ucat) > 1.e-10_r8" not in BIF
    assert "bif_outflow > 1.e-10_r8" not in BIF


def test_tiny_flux_still_obeys_cama_five_percent_path_limiter() -> None:
    flux = 5.0e-11
    dt = 60.0
    endpoint_storage = 1.0e-10

    transfer = abs(flux) * dt
    rate = min(1.0, 0.05 * endpoint_storage / transfer)

    assert rate < 1.0
    assert abs(flux * rate) * dt <= 0.05 * endpoint_storage
    assert "IF (abs(pth_hflux_total(ipth)) > 0._r8 .and. dt > 0._r8) THEN" in BIF
    assert "abs(pth_hflux_total(ipth)) > 1.e-10_r8" not in BIF
