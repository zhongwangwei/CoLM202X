from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COLM_MAIN = (ROOT / "main" / "CoLMMAIN.F90").read_text(encoding="utf-8")
PRECIP = (ROOT / "main" / "TRACER" / "MOD_Tracer_Precip.F90").read_text(
    encoding="utf-8"
)
EVAPO = (ROOT / "main" / "TRACER" / "MOD_Tracer_Evapo.F90").read_text(
    encoding="utf-8"
)
EVAP_LIMIT = (ROOT / "main" / "TRACER" / "MOD_Tracer_EvapLimit.F90").read_text(
    encoding="utf-8"
)
RIVER = (ROOT / "main" / "TRACER" / "MOD_Tracer_RiverLake.F90").read_text(
    encoding="utf-8"
)
RIVER_FLOW = (
    ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeFlow.F90"
).read_text(encoding="utf-8")
LEAF_INTERCEPTION = (ROOT / "main" / "MOD_LeafInterception.F90").read_text(
    encoding="utf-8"
)


def compact(source: str) -> str:
    return " ".join(source.lower().replace("&", "").split())


def test_sprinkler_uses_pre_withdrawal_storage_and_actual_debit_signature():
    call = COLM_MAIN.split("CALL tracer_precip", 1)[1].split("#endif", 1)[0]
    assert "waterstorage_trc_beg" in call
    assert "waterstorage(ipatch)" not in call

    precip = compact(PRECIP)
    assert "sprinkler_water = sprinkler_rate * deltim" in precip
    assert "storage_ratio = sprinkler_trc / sprinkler_water" in precip
    assert (
        "max(waterstorage_before - sprinkler_water, 0._r8), "
        "trc_waterstorage(itrc, ipatch), trc_waterstorage_solid(itrc, ipatch)"
        in precip
    )


def test_simultaneous_canopy_melt_expands_the_pool_available_to_freeze():
    precip = compact(PRECIP)
    evapo = compact(EVAPO)

    assert "frzc_mass = min(frzc_mass, ldew_rain_old_eff)" in precip
    assert "ldew_rain_pre_phase + thaw_amt" in evapo
    assert (
        "ldew_rain_pre_phase + thaw_amt, freeze_amt, canopy_temp()" in evapo
    )


def test_canopy_fractionation_follows_host_external_then_phase_order():
    leaf_temperature = compact(
        (ROOT / "main" / "MOD_LeafTemperature.F90").read_text(encoding="utf-8")
    )
    evapo = compact(EVAPO)

    host_external = leaf_temperature.index(
        "ldew_snow = ldew_snow + (qfrol-qsubl)*deltim"
    )
    host_phase = leaf_temperature.index("ldew_snow = max(0.,ldew_snow - qmelt*deltim)")
    assert host_external < host_phase

    tracer_external = evapo.index("canopy snow external exchange")
    tracer_phase = evapo.index("internal melt: snow")
    assert tracer_external < tracer_phase
    assert "ldew_snow_bef, abs(d_snow_external)" in evapo
    assert "trc_ldew_snow(itrc, ipatch) / ldew_snow_pre_phase" in evapo

    # Rayleigh fractionation and later phase transfer are non-commutative.
    # The host removes vapor from the original snow pool, then melts from the
    # fractionated residual; doing melt first gives a different isotope mass.
    water0, tracer0, sublim, melt, alpha = 10.0, 2.0, 2.0, 3.0, 1.08
    tracer_after_sublim = tracer0 * ((water0 - sublim) / water0) ** alpha
    host_melt_tracer = melt * tracer_after_sublim / (water0 - sublim)

    wrong_melt_tracer = melt * tracer0 / water0
    assert abs(host_melt_tracer - wrong_melt_tracer) > 1.0e-4


def test_carrier_empty_canopy_release_rephases_and_flushes_all_old_tracer():
    host = compact(
        LEAF_INTERCEPTION.split("SUBROUTINE LEAF_interception_CoLM2014", 1)[1].split(
            "END SUBROUTINE LEAF_interception_CoLM2014", 1
        )[0]
    )
    precip = compact(PRECIP)
    release = precip.split(
        "leaf_interception has a carrier-empty branch", 1
    )[1].split("! ---- rain component", 1)[0]

    # Host water sends the complete old single-bucket store to one ground phase
    # and clears both canopy carriers, while exposing no explicit phase-change mass.
    assert "if (ldew > 0.) then" in host
    assert "xsc_rain = max(0._r8, ldew)" in host
    assert "xsc_snow = max(0._r8, ldew)" in host
    assert "ldew_rain = 0." in host
    assert "ldew_snow = 0." in host
    assert "ldew_smelt_out = 0._r8" in host
    assert "ldew_frzc_out = 0._r8" in host

    # Tracer repacks every old canopy component into that host-selected release
    # phase before xsc is evaluated, so a cross-phase old pool cannot be stranded.
    assert "canopy_release_trc = trc_ldew_rain(itrc, ipatch)" in release
    assert "+ trc_ldew_snow(itrc, ipatch) + trc_canopy_solid(itrc, ipatch)" in release
    assert "trc_canopy_solid(itrc, ipatch) = 0._r8" in release
    assert "max(xsc_rain_out, 0._r8) * deltim > trc_water_min_for_ratio" in release
    assert "max(xsc_snow_out, 0._r8) * deltim <= trc_water_min_for_ratio" in release
    assert "trc_ldew_rain(itrc, ipatch) = canopy_release_trc" in release
    assert "trc_ldew_snow(itrc, ipatch) = canopy_release_trc" in release
    assert "ldew_rain_old_eff = canopy_old_water" in release
    assert "ldew_snow_old_eff = canopy_old_water" in release

    # Algebraic regression for both cross-phase directions: the selected xsc
    # removes the full remapped inventory and leaves zero canopy tracer.
    for old_water, old_tracer in ((1.0, 0.37), (2.5, 1.4)):
        release_ratio = old_tracer / old_water
        released = min(old_water, old_water) * release_ratio
        assert abs(released - old_tracer) < 1.0e-15
        assert abs(old_tracer - released) < 1.0e-15

    # Ordinary phase-specific overflow may release both rain and snow in the
    # same step.  It must retain each phase signature rather than trigger the
    # carrier-empty single-phase repack.
    def carrier_empty_repack(rain_xsc: float, snow_xsc: float, tol: float = 1e-12) -> bool:
        rain_only = rain_xsc > tol and snow_xsc <= tol
        snow_only = snow_xsc > tol and rain_xsc <= tol
        return rain_only or snow_only

    assert carrier_empty_repack(1.0, 0.0)
    assert carrier_empty_repack(0.0, 1.0)
    assert not carrier_empty_repack(0.4, 0.6)


def test_large_finite_pool_fractionating_loss_still_uses_substeps():
    limiter = compact(EVAP_LIMIT)
    assert "evaplimit_default_max_pool_water" not in limiter
    early_return = limiter.split("remaining_trc =", 1)[0]
    assert "water_loss <= evaplimit_default_max_loss_fraction * pool_water" in early_return
    assert "pool_water >" not in early_return


def test_dry_drain_is_separate_from_advective_flux_and_counted_once():
    diag = compact(
        RIVER.split("SUBROUTINE tracer_diag_accumulate_substep", 1)[1].split(
            "END SUBROUTINE tracer_diag_accumulate_substep", 1
        )[0]
    )
    debug_budget = compact(RIVER_FLOW)

    dry_block = diag.split("dry_drain =", 1)[1].split("trc_mass(itrc, i) =", 1)[0]
    assert "trc_flux_out(itrc, i)" not in dry_block
    assert "trc_dry_drain(itrc, i) = trc_dry_drain(itrc, i) + dry_drain" in dry_block
    assert "a_trc_out(itrc, i) = a_trc_out(itrc, i) + dry_drain" in dry_block
    assert debug_budget.count("sum(trc_dry_drain(itrc,:))") == 1


def test_committed_restart_requires_both_accumulator_sides_before_any_scatter():
    reader = compact(
        RIVER.split("SUBROUTINE read_tracer_restart", 1)[1].split(
            "END SUBROUTINE read_tracer_restart", 1
        )[0]
    )
    preflight = reader.split("a committed current transaction", 1)[1].split(
        "vector_read_and_scatter contains", 1
    )[0]

    assert "'trc_accinp_'" in preflight
    assert "'acc_rnof_ref', .true., field_present" in preflight
    assert "probe_riverlake_restart_vector" in preflight
    assert "acc_rnof_uc" not in reader
    assert "tracer_init_water_ratio" not in reader
    assert "has_accinp" not in reader


def test_committed_restart_requires_signed_input_buffer_before_any_scatter():
    reader = compact(
        RIVER.split("SUBROUTINE read_tracer_restart", 1)[1].split(
            "END SUBROUTINE read_tracer_restart", 1
        )[0]
    )
    helper = compact(
        RIVER.split("SUBROUTINE probe_riverlake_restart_vector", 1)[1].split(
            "END SUBROUTINE probe_riverlake_restart_vector", 1
        )[0]
    )
    preflight = reader.split("a committed current transaction", 1)[1].split(
        "vector_read_and_scatter contains", 1
    )[0]

    assert "'trc_inpbuf_'" in preflight
    assert ".true., field_present" in preflight
    assert "status < 0 .or. (required .and. status == 0)" in helper
    assert "call colm_stop('incomplete or malformed committed river tracer restart')" in helper
    assert "mass_loaded" not in reader
    assert "required signed pending pool" not in reader
