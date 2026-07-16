from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIF = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeBifurcation.F90").read_text(
    encoding="utf-8"
)


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_dynamic_state_uses_one_fixed_order_batch() -> None:
    calc = routine(BIF, "bifurcation_calc")

    assert "type(worker_push_real8_field_type) :: dynamic_state_fields(7)" in calc
    assert "n_dynamic_state_fields = 4" in calc
    assert "dynamic_state_fields(5)%send => visible_storage_ucat" in calc
    assert "n_dynamic_state_fields = 5" in calc
    assert "dynamic_state_fields(6)%send => protected_wdsrf_ucat" in calc
    assert "n_dynamic_state_fields = 7" in calc
    assert calc.count(
        "CALL worker_push_data (push_bif_dn2pth, &\n"
        "         dynamic_state_fields(1:n_dynamic_state_fields))"
    ) == 1

    for scalar_call in (
        "push_bif_dn2pth, wdsrf_ucat",
        "push_bif_dn2pth, protected_wdsrf_ucat",
        "push_bif_dn2pth, storage_ucat",
        "push_bif_dn2pth, visible_storage_ucat",
        "push_bif_dn2pth, protected_storage_ucat",
        "push_bif_dn2pth, ucatfilter_r8",
    ):
        assert scalar_call not in calc


def test_reverse_rates_and_final_influx_each_use_one_batch() -> None:
    calc = routine(BIF, "bifurcation_calc")

    assert calc.count(
        "CALL worker_push_data (push_bif_influx, &\n"
        "         donor_reverse_fields(1:n_coupled_fields), mode = 'sum')"
    ) == 1
    assert calc.count(
        "CALL worker_push_data (push_bif_dn2pth, &\n"
        "         donor_rate_fields(1:n_coupled_fields))"
    ) == 1
    assert calc.count(
        "CALL worker_push_data (push_bif_influx, &\n"
        "         final_influx_fields(1:n_final_influx_fields), mode = 'sum')"
    ) == 1

    assert "n_coupled_fields = merge(2, 1, use_protected_fields)" in calc
    assert "n_final_influx_fields = merge(2, 1, npthlev_bif >= 2)" in calc

    for scalar_call in (
        "push_bif_influx, protected_neg_pth",
        "push_bif_influx, neg_pth",
        "push_bif_dn2pth, protected_out_rate",
        "push_bif_dn2pth, limiter_out_rate",
        "push_bif_influx, pth_hflux_total",
        "push_bif_influx, pth_lev_hflux_total",
    ):
        assert scalar_call not in calc


def test_protected_only_fields_are_omitted_without_multilevel_levee() -> None:
    calc = routine(BIF, "bifurcation_calc")

    assert "use_protected_fields = DEF_USE_LEVEE .and. npthlev_bif >= 2" in calc
    assert "IF (use_protected_fields) THEN" in calc
    assert "n_dynamic_state_fields = 7" in calc
    assert "n_coupled_fields = merge(2, 1, use_protected_fields)" in calc
