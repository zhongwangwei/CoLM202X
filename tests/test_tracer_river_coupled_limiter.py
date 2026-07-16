from math import isclose
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RIVER = ROOT / "main" / "TRACER" / "MOD_Tracer_RiverLake.F90"
NETWORK = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeNetwork.F90"
FLOW = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeFlow.F90"
BIF = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeBifurcation.F90"


def _safe_fixed_point(mass, outgoing, edges, max_iter=256, tol=1.0e-12):
    """Mirror the monotone donor-rate equations used by tracer_substep."""
    rates = [min(max(m, 0.0) / out, 1.0) if out > 1.0e-30 else 1.0 for m, out in zip(mass, outgoing)]
    for _ in range(max_iter):
        incoming = [0.0] * len(mass)
        for donor, receiver, raw_amount in edges:
            incoming[receiver] += raw_amount * rates[donor]
        new = [
            max(old, min((max(m, 0.0) + inc) / out, 1.0)) if out > 1.0e-30 else 1.0
            for old, m, inc, out in zip(rates, mass, incoming, outgoing)
        ]
        if max(abs(a - b) for a, b in zip(new, rates)) <= tol:
            return new
        rates = new
    raise AssertionError("test fixed point did not converge")


def test_three_cell_reverse_chain_no_longer_spends_unscaled_incoming_credit():
    # Nominal network C -> B -> A -> sea; all three face fluxes are negative,
    # so actual water moves sea -> A -> B -> C during this one-second substep.
    # The water solver accepts both cells because it limits NET, not gross, loss.
    dt = 1.0
    volume_a, volume_b = 0.25, 1.0
    q_sea_a, q_a_b, q_b_c = 0.25, 0.45, 1.40
    assert isclose(volume_a - (q_a_b - q_sea_a) * dt, 0.05, abs_tol=1.0e-15)
    assert isclose(volume_b - (q_b_c - q_a_b) * dt, 0.05, abs_tol=1.0e-15)

    # volflux=max(volume, abs(own-face flux)*dt) gives unit raw concentration
    # in A and B. Sea/boundary concentration is zero for this solute example.
    mass = [0.25, 1.0, 0.0]  # A, B, C
    assert mass[0] / max(volume_a, q_sea_a * dt) == 1.0
    assert mass[1] / max(volume_b, q_a_b * dt) == 1.0
    outgoing = [q_a_b, q_b_c, 0.0]
    edges = [(0, 1, q_a_b), (1, 2, q_b_c)]

    # Former one-shot limiter used raw incoming at B, so it set r_B=1 even
    # though A subsequently cut that incoming to its available 0.25 mass.
    old_rate_a = mass[0] / outgoing[0]
    old_rate_b = min((mass[1] + q_a_b) / outgoing[1], 1.0)
    old_mass_b = mass[1] + q_a_b * old_rate_a - outgoing[1] * old_rate_b
    assert isclose(old_mass_b, -0.15, rel_tol=0.0, abs_tol=1.0e-14)

    rates = _safe_fixed_point(mass, outgoing, edges)
    final = [mass[0], mass[1], mass[2]]
    for donor, receiver, raw_amount in edges:
        moved = raw_amount * rates[donor]
        final[donor] -= moved
        final[receiver] += moved

    assert isclose(rates[0], 5.0 / 9.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert isclose(rates[1], 1.25 / 1.40, rel_tol=0.0, abs_tol=1.0e-14)
    assert min(final) >= -1.0e-15
    assert isclose(final[1], 0.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert isclose(sum(final), sum(mass), rel_tol=0.0, abs_tol=1.0e-14)


def test_source_iterates_actual_scaled_main_and_bif_incoming_monotonically():
    source = RIVER.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]

    assert "Start at T(0)" in substep
    initial = substep.split("Start at T(0)", 1)[1].split(
        "DO limiter_iter = 1, limiter_max_iter", 1
    )[0]
    assert "trc_in_mass(i)" not in initial
    assert "max(trc_mass(itrc, i), 0._r8) / trc_out_mass(i)" in initial
    assert "max(trc_levsto(itrc, i), 0._r8) / trc_out_mass_lev(i)" in initial

    iteration = substep.split("DO limiter_iter = 1, limiter_max_iter", 1)[1].split(
        "IF (.not. limiter_converged) THEN", 1
    )[0]
    assert "trc_inp_step(i) = max(trc_flux(i) * rate_cell(i), 0._r8)" in iteration
    assert "max(-trc_flux(i) * rate_next(i), 0._r8)" in iteration
    assert "trc_pth_fl = trc_pth_fl * trc_rate" in iteration
    assert "trc_rate = rate_cell_lev(i_up)" in iteration
    assert "trc_rate = rate_dn_pth_lev(ipth)" in iteration
    assert "limiter_rate_new = max(rate_cell(i), limiter_rate_new)" in iteration
    assert "limiter_rate_new = max(rate_cell_lev(i), limiter_rate_new)" in iteration


def test_limiter_uses_only_collectives_whose_water_loops_are_lockstepped():
    source = RIVER.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]

    assert "IF (bif_workspace_active) THEN" in substep
    assert "p_comm_worker, p_err" in substep
    assert "ELSEIF (rivsys_by_multiple_procs) THEN" in substep
    assert "p_comm_rivsys, p_err" in substep
    assert "river tracer donor limiter did not converge" in substep
    assert "TRC_LIMITER_RATE_TOL" in substep
    assert "ieee_is_finite(limiter_delta_global)" in substep
    assert "negative river tracer mass after coupled donor limiter" in substep
    assert "negative protected tracer mass after coupled donor limiter" in substep


def test_iteration_cap_covers_long_acyclic_information_paths():
    # Jacobi communication moves a newly authorised rate one edge per round.
    # A hard 256 cap leaves the end of this valid 300-edge path unresolved.
    reach = 300
    rates = [1.0] + [0.0] * reach
    for _ in range(256):
        rates = [1.0] + [max(rates[i], rates[i - 1]) for i in range(1, len(rates))]
    assert rates[-1] == 0.0

    dynamic_cap = 2 * len(rates) + 1
    for _ in range(256, dynamic_cap):
        new = [1.0] + [max(rates[i], rates[i - 1]) for i in range(1, len(rates))]
        if new == rates:
            break
        rates = new
    assert rates[-1] == 1.0

    source = RIVER.read_text(encoding="utf-8")
    assert "TRC_LIMITER_MAX_ITER" not in source
    assert "TRC_LIMITER_CYCLE_MIN_ITER" not in source
    assert "totalnumucat > (huge(limiter_max_iter) - 1) / 2" in source
    assert "limiter_max_iter = 2 * totalnumucat + 1" in source


def test_bif_water_limiter_keeps_subunit_tracer_dependencies_on_main_forest():
    flow = FLOW.read_text(encoding="utf-8")
    bif = BIF.read_text(encoding="utf-8")
    source = RIVER.read_text(encoding="utf-8")

    # Flow supplies gross ordinary donor outflow, including reverse faces.
    assert "normal_outgoing_rate(i) = hflux_fc(i)" in flow
    assert "normal_outgoing_rate = normal_outgoing_rate + mflux_sumups" in flow
    # BIF can use only storage left after that gross ordinary demand; the rate
    # is then applied to every outgoing BIF layer before tracer_substep sees it.
    assert "remaining_capacity = max(storage_ref / dt_cell - normal_outflow, 0._r8)" in bif
    assert "limiter_out_rate(i_ucat) = min(1._r8, remaining_capacity / bif_outflow)" in bif
    assert "bif_hflux_lev(ilev, ipth) = bif_hflux_lev(ilev, ipth) * rate" in bif
    assert "sub-unit-rate dependencies propagate only along the ordinary river" in source


def test_negative_roundoff_dust_never_forms_a_negative_transport_flux():
    source = RIVER.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]
    concentration = substep.split("! --- 1. Concentration", 1)[1].split(
        "! --- 3. Get downstream concentration", 1
    )[0]
    update = substep.split("! --- 8. Update tracer mass ---", 1)[1].split(
        "! --- 9. Save flux", 1
    )[0]

    assert "trc_mass(itrc, i) < -TRC_RESTART_NEGATIVE_DUST" in concentration
    assert "trc_mass(itrc, i) = 0._r8" in concentration
    assert "trc_levsto(itrc, i) < -TRC_RESTART_NEGATIVE_DUST" in concentration
    assert "trc_levsto(itrc, i) = 0._r8" in concentration
    assert concentration.index("trc_mass(itrc, i) = 0._r8") < concentration.index(
        "trc_conc_flux(i) = trc_mass(itrc, i) / volwater"
    )
    assert "trc_mass(itrc, i) = max(trc_mass_new, 0._r8)" in update
    assert "trc_levsto(itrc, i) = max(trc_levsto(itrc, i), 0._r8)" in update


def test_dry_protected_pool_never_borrows_visible_concentration():
    dryoff = 1.0e-6
    protected_volume = 1.0e-9
    protected_mass = 0.4e-9
    protected_water_out = protected_volume  # water limiter maximum for dt=1
    visible_concentration = 1.0

    # The former fallback to visible concentration overdrew the protected
    # tracer pool even though the protected WATER limiter was satisfied.
    old_outgoing = visible_concentration * protected_water_out
    assert old_outgoing > protected_mass
    protected_concentration = protected_mass / protected_volume
    assert isclose(
        protected_concentration * protected_water_out,
        protected_mass,
        rel_tol=0.0,
        abs_tol=1.0e-30,
    )

    source = RIVER.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]
    concentration = substep.split("! --- 1. Concentration", 1)[1].split(
        "! --- 3. Get downstream concentration", 1
    )[0]
    protected = concentration.split("trc_prot_conc_flux(i) = trc_conc_flux(i)", 1)[1]

    assert "IF (has_levee(i)) THEN" in protected
    assert "has_levee(i) .and. levsto(i) > trc_v_dry_off" not in protected
    assert "IF (levsto(i) > 0._r8) THEN" in protected
    assert "trc_levsto(itrc, i) / levsto(i)" in protected
    assert "trc_prot_conc_flux(i) = 0._r8" in protected
    assert "max(levsto(i), trc_v_dry_off)" not in protected
    assert "ieee_is_finite(levsto(i))" in protected
    assert "levsto(i) < 0._r8" in protected


def test_reverse_main_edge_uses_one_river_system_timestep_contract():
    network = NETWORK.read_text(encoding="utf-8")
    flow = FLOW.read_text(encoding="utf-8")
    source = RIVER.read_text(encoding="utf-8")

    # Every ordinary upstream/downstream edge inherits one rivermouth label;
    # split systems use irivsys=1 and reduce that scalar dt on their communicator.
    assert "rivermouth(uc_up2down(i)) = rivermouth(j)" in network
    assert "numrivsys  = 1" in network
    assert "irivsys(:) = 1" in network
    assert "MPI_MIN" in flow
    assert "p_comm_rivsys" in flow
    assert "this edge/receiver dt is also the downstream donor dt" in source


def test_wet_flux_uses_true_pool_concentration_and_preserves_constant_state():
    volume = 1.0
    mass = 1.0
    incoming_water = 0.5
    outgoing_water = 1.4
    volume_next = volume + incoming_water - outgoing_water
    assert isclose(volume_next, 0.1, abs_tol=1.0e-15)

    # Uniform concentration is one on both incoming and local water. The true
    # pool concentration plus coupled credit leaves exactly C*V_next.
    true_concentration = mass / volume
    final_mass = mass + incoming_water - true_concentration * outgoing_water
    assert isclose(final_mass, volume_next, abs_tol=1.0e-15)

    # The former inflated denominator would under-export tracer and destroy C=1.
    old_concentration = mass / max(volume, outgoing_water)
    old_final_mass = mass + incoming_water - old_concentration * outgoing_water
    assert not isclose(old_final_mass, volume_next, abs_tol=1.0e-15)

    source = RIVER.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]
    assert "IF (volwater > trc_v_dry_off) THEN" in substep
    assert "trc_conc_flux(i) = trc_mass(itrc, i) / volwater" in substep
    assert "max(volwater, abs(hflux_fc(i)) * dt_i)" not in substep
    assert "volflux = max(abs(hflux_fc(i)) * dt_i, trc_v_dry_off)" in substep


def test_post_advection_mass_is_never_resynced_to_initial_signature():
    # The former 1e-6 relative snap silently removed this small but legitimate
    # transported/restarted gradient without booking a source or sink.
    r_fill = 2.0e-3
    volume_next = 100.0
    conservative_mass = r_fill * volume_next + 1.0e-7
    old_snapped_mass = r_fill * volume_next
    assert conservative_mass != old_snapped_mass

    source = RIVER.read_text(encoding="utf-8")
    substep = source.split("SUBROUTINE tracer_substep", 1)[1].split(
        "END SUBROUTINE tracer_substep", 1
    )[0]
    update = substep.split("! --- 8. Update tracer mass ---", 1)[1].split(
        "! --- 9. Save flux", 1
    )[0]
    assert "trc_mass(itrc, i) = max(trc_mass_new, 0._r8)" in update
    assert "trc_mass(itrc, i) = R_fill" not in substep
    assert "fixed_sig_rel_tol" not in substep
    assert "tracer_can_use_fixed_signature" not in substep
    assert "tracer_fractionation_active" not in substep
    assert "trc_runtime_forced" not in substep
