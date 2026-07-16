from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
BGC = ROOT / "main" / "BGC"
TRACER = ROOT / "main" / "TRACER"


def source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def routine_body(text: str, name: str) -> str:
    match = re.search(
        rf"SUBROUTINE\s+{name}\b(?P<body>.*?)END\s+SUBROUTINE\s+{name}",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    assert match is not None, f"subroutine {name} not found"
    return match.group("body")


def without_continuations(text: str) -> str:
    return re.sub(r"&\s*\n\s*", " ", text)


def test_wetland_decomposition_runs_n_limited_fluxes_before_pool_update() -> None:
    shim = source(TRACER / "MOD_Tracer_Reactive_BgcShim.F90")
    body = routine_body(shim, "reactive_bgc_run_wetland_decomp")

    potential = body.index("CALL SoilBiogeochemPotential")
    competition = body.index("CALL SoilBiogeochemCompetitionNoPlant")
    decomposition = body.index("CALL SoilBiogeochemDecomp")
    assert potential < competition < decomposition
    assert "decomp_cpools_sourcesink  (1:nl_soil,:,ipatch) = 0._r8" in body
    assert "decomp_npools_sourcesink  (1:nl_soil,:,ipatch) = 0._r8" in body

    competition_source = source(BGC / "MOD_BGC_Soil_BiogeochemCompetition.F90")
    no_plant = routine_body(competition_source, "SoilBiogeochemCompetitionNoPlant")
    assert "patch_pft_s" not in no_plant
    assert "actual_immob_vr(j,i) = min" in no_plant
    assert "max(sminn_vr(j,i), 0._r8) / deltim" in no_plant


def test_wetland_finalizer_applies_only_extracted_decomposition_state_paths() -> None:
    carbon = source(BGC / "MOD_BGC_CNCStateUpdate1.F90")
    nitrogen = source(BGC / "MOD_BGC_Soil_BiogeochemNStateUpdate1.F90")
    driver = source(BGC / "MOD_BGC_driver.F90")
    link = source(TRACER / "MOD_Tracer_Reactive_Methane_BgcLink.F90")
    finalizer = without_continuations(
        routine_body(link, "tracer_ch4_bgc_finalize_step")
    )
    normalized_finalizer = " ".join(finalizer.split())

    assert "CALL CDecompStateUpdate(i, deltim, nl_soil, ndecomp_transitions)" in carbon
    assert (
        "CALL SoilBiogeochemNDecompStateUpdate(i, deltim, nl_soil, "
        "ndecomp_transitions)"
        in without_continuations(nitrogen)
    )
    assert (
        "CALL CDecompStateUpdate(ipatch, deltim, nl_soil, "
        "size(decomp_hr_vr,2), .true.)"
        in normalized_finalizer
    )
    assert (
        "CALL SoilBiogeochemNDecompStateUpdate(ipatch, deltim, nl_soil, "
        "size(decomp_hr_vr,2), .true.)"
        in normalized_finalizer
    )
    assert "decomp_cpools_vr(1:nl_soil,:,i) =" in carbon
    assert "decomp_npools_vr(1:nl_soil,:,i) =" in nitrogen
    assert "update_pools = .false." in carbon
    assert "update_pools = .false." in nitrogen
    assert "IF (update_pools) THEN" in carbon
    assert "IF (update_pools) THEN" in nitrogen
    assert "CDecompStateUpdate" not in driver
    assert "SoilBiogeochemNDecompStateUpdate" not in driver
    assert driver.index("CALL CStateUpdate1") < driver.index(
        "CALL SoilBiogeochemLittVertTransp"
    )
    assert driver.index("CALL SoilBiogeochemNStateUpdate1") < driver.index(
        "CALL SoilBiogeochemLittVertTransp"
    )
    assert "IF (patchtype == 2) THEN" in finalizer


def test_co2_history_partition_preserves_total_decomposed_carbon() -> None:
    link = source(TRACER / "MOD_Tracer_Reactive_Methane_BgcLink.F90")
    physics = source(TRACER / "MOD_Tracer_Reactive_Methane_Physics.F90")
    constants = source(TRACER / "MOD_Tracer_Reactive_Methane_Const.F90")
    history = source(TRACER / "MOD_Tracer_Reactive_Methane_Hist.F90")
    finalizer = without_continuations(
        routine_body(link, "tracer_ch4_bgc_finalize_step")
    )

    assert "net_methane_unsat = -methane_prod_tot_unsat + methane_oxid_tot_unsat" in physics
    assert "net_methane_sat   = -methane_prod_tot_sat   + methane_oxid_tot_sat" in physics
    assert "catomw = 12.011_r8" in constants
    assert "co2_hr = total_hr + catomw * net_methane" in finalizer
    assert ".not. ieee_is_finite(total_hr)" in finalizer
    assert ".not. ieee_is_finite(net_methane)" in finalizer
    assert ".not. ieee_is_finite(co2_hr)" in finalizer
    assert "decomp_hr(ipatch) = max(co2_hr, 0._r8)" in finalizer
    assert "er(ipatch) = ar(ipatch) + decomp_hr(ipatch)" in finalizer
    assert "CH4 oxidation minus production; applied to soil/wetland BGC CO2 respiration" in history
    assert "average net methane correction to CO2 flux" not in history
    assert "decomp_hr_vr" not in finalizer.split("total_hr =", 1)[1].split(
        "co2_hr =", 1
    )[1]

    # net_methane = oxidation - production [mol C m-2 s-1].  Therefore the
    # BGC pool loss is recovered from the published CO2-C flux by subtracting
    # the already-applied methane correction.
    catomw = 12.011
    total_hr = 0.36
    net_methane = -0.01
    co2_hr = total_hr + catomw * net_methane
    assert abs((co2_hr - catomw * net_methane) - total_hr) < 1.0e-14


def test_decomposition_source_sink_algebra_conserves_carbon_and_nitrogen() -> None:
    # Representative nonterminal and terminal transitions exercise the same
    # signs used by CDecompStateUpdate and SoilBiogeochemNDecompStateUpdate.
    donors = (0, 1, 2)
    receivers = (1, 2, None)
    hr = (0.08, 0.03, 0.02)
    ctransfer = (0.12, 0.05, 0.0)
    ntransfer = (0.006, 0.003, 0.0)
    sminn_flux = (0.001, -0.0005, 0.002)
    deltim = 1800.0

    carbon_delta = [0.0, 0.0, 0.0]
    organic_n_delta = [0.0, 0.0, 0.0]
    mineral_n_delta = 0.0
    for donor, receiver, hr_flux, c_flux, n_flux, mineral_flux in zip(
        donors, receivers, hr, ctransfer, ntransfer, sminn_flux
    ):
        carbon_delta[donor] -= (hr_flux + c_flux) * deltim
        organic_n_delta[donor] -= n_flux * deltim
        if receiver is not None:
            carbon_delta[receiver] += c_flux * deltim
            organic_n_delta[receiver] += (n_flux + mineral_flux) * deltim
            mineral_n_delta -= mineral_flux * deltim
        else:
            organic_n_delta[donor] -= mineral_flux * deltim
            mineral_n_delta += mineral_flux * deltim

    assert abs(sum(carbon_delta) + sum(hr) * deltim) < 1.0e-12
    assert abs(sum(organic_n_delta) + mineral_n_delta) < 1.0e-12


def test_bgc_finalization_precedes_patch_diagnostics() -> None:
    driver = source(TRACER / "MOD_Tracer_Reactive_Methane_Driver.F90")
    methane_driver = routine_body(driver, "methane_driver")

    # Physics may run twice for a mixed soil/rice patch. Finalization belongs
    # to the outer driver and occurs once after the completed columns merge.
    aggregate = methane_driver.index("CALL aggregate_methane_columns")
    finalize = methane_driver.index("CALL tracer_ch4_bgc_finalize_step", aggregate)
    assert aggregate < finalize
    assert "tracer_ch4_bgc_finalize_step" not in source(
        TRACER / "MOD_Tracer_Reactive_Methane_Physics.F90"
    )
