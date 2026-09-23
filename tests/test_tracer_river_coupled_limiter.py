import inspect
import os
import re
import shutil
import subprocess
import tempfile
from math import isclose
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
RIVER = ROOT / "main" / "TRACER" / "MOD_Tracer_RiverLake.F90"
NETWORK = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeNetwork.F90"
FLOW = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeFlow.F90"
BIF = ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeBifurcation.F90"


@pytest.mark.skipif(os.environ.get("COLM_RUN_MPI_TESTS") != "1", reason="opt-in MPI integration")
def test_production_tracer_substep_mpi_budget_and_fallback():
    compiler = shutil.which("mpif90")
    launcher = shutil.which("mpiexec") or shutil.which("mpirun")
    make = shutil.which("gmake") or shutil.which("make")
    if not (ROOT / ".bld").is_dir() or not all(
        (compiler, launcher, make, shutil.which("nf-config"))
    ):
        pytest.skip("built CoLM modules and MPI/NetCDF toolchain required")
    with tempfile.TemporaryDirectory(prefix="colm-tracer-bif-") as directory:
        work = Path(directory)
        flags = ["-O0", "-g", "-fopenmp", "-fdefault-real-8", "-ffree-form",
                 "-fcheck=all", "-fallow-argument-mismatch", "-cpp", "-ffree-line-length-0", "-w",
                 "-Iinclude", f"-I{work}", "-I.bld", f"-J{work}"]
        tracer_object = work / "MOD_Tracer_RiverLake.o"
        subprocess.run([compiler, *flags, "-c", str(RIVER), "-o",
                        str(tracer_object)], cwd=ROOT, check=True)
        harness = work / "harness.o"
        subprocess.run([compiler, *flags, "-c", "tests/river_bif_tracer_harness.F90",
                        "-o", str(harness)], cwd=ROOT, check=True)
        stub = work / "providers.f90"
        stub.write_text("SUBROUTINE register_all_tracer_providers()\nEND SUBROUTINE\n")
        subprocess.run([compiler, "-c", str(stub), "-o", str(work / "providers.o")],
                       cwd=ROOT, check=True)
        objects = subprocess.check_output([
            make, "--no-print-directory",
            "--eval=p:;@echo $(OBJS_SHARED_T) $(TRACER_RUNTIME_CONFIG_OBJS_T) $(OBJS_BASIC_T)",
            "p"], cwd=ROOT, text=True).split()
        objects = [str(tracer_object) if obj == ".bld/MOD_Tracer_RiverLake.o"
                   else obj for obj in objects]
        netcdf = subprocess.check_output(["nf-config", "--flibs"], text=True).split()
        brew = shutil.which("brew")
        if brew:
            netcdf.insert(0, "-L" + subprocess.check_output(
                [brew, "--prefix", "netcdf"], text=True).strip() + "/lib")
        binary = work / "harness"
        subprocess.run([compiler, "-fopenmp", "-o", str(binary), *objects,
                        str(harness), str(work / "providers.o"), *netcdf,
                        "-llapack", "-lblas"], cwd=ROOT, check=True)
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS="1", OMPI_MCA_rmaps_base_oversubscribe="1")
        for ranks in (2, 3, 4):
            result = subprocess.run([launcher, "-n", str(ranks), str(binary)],
                                    cwd=ROOT, env=env, text=True, capture_output=True,
                                    timeout=60, check=True)
            assert f"TRACER_DYNAMIC ranks={ranks} failures=0" in result.stdout
            assert result.stdout.count("TRACER_DYNAMIC scenario=") == 4
            assert "TRACER_BUDGET initial=" in result.stdout
            assert "TRACER_FINITE_BIF moved=" in result.stdout
            assert "TRACER_UNIFORM cold_placeholder error=" in result.stdout
            assert result.stdout.count("TRACER_UNIFORM transport=") == 3


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


def test_long_chain_requires_topology_bound_and_is_never_truncated():
    # Jacobi communication moves a newly authorised donor rate exactly one edge
    # per round, so a legitimate long chain needs more rounds than any small
    # fast-path threshold.  Truncating there would under-export tracer, and the
    # dry-cell cleanup books the leftovers to the trc_dry_drain SINK instead of
    # delivering them downstream, so the loop must run to the 2*N+1 topology
    # bound and the small threshold may only be reported.
    reach = 300
    rates = [1.0] + [0.0] * reach
    rounds = 0
    for _ in range(2 * len(rates) + 1):
        new = [1.0] + [max(rates[i], rates[i - 1]) for i in range(1, len(rates))]
        if new == rates:
            break
        rates = new
        rounds += 1
    assert rounds > 32          # a 32-round budget is NOT enough
    assert rates[-1] == 1.0     # the topology bound does resolve the chain

    source = RIVER.read_text(encoding="utf-8")
    # The loop bound is the topology bound; the small threshold is diagnostic.
    assert "limiter_max_iter = 2 * totalnumucat + 1" in source
    assert "min(TRC_LIMITER_SOFT_ITER" not in source
    assert "TRC_LIMITER_SOFT_ITER" in source
    assert "IF (limiter_iter == TRC_LIMITER_SOFT_ITER)" in source
    # No residual-based truncation / no "accept the approximate iterate" path.
    assert "TRC_LIMITER_CAPPED_FATAL_RESID" not in source
    assert "TRC_LIMITER_MAX_ITER" not in source
    assert "conservative iterate" not in source
    # Exhausting the topology bound is still fatal (a broken limiter).
    assert "CALL CoLM_stop('river tracer donor limiter did not converge')" in source
    assert "totalnumucat > (huge(limiter_max_iter) - 1) / 2" in source
    # Iteration accounting is published, and the caller reports it globally.
    stats = source.split("SUBROUTINE tracer_limiter_stats", 1)[1].split("END SUBROUTINE", 1)[0]
    assert "over_soft" in stats and "capped" not in stats
    assert "trc_limiter_iter_sum = trc_limiter_iter_sum + limiter_iters_used" in source
    flow = FLOW.read_text(encoding="utf-8")
    assert "CALL tracer_limiter_stats" in flow
    assert "mpi_allreduce (MPI_IN_PLACE, lim_calls, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)" in flow
    assert "mpi_allreduce (MPI_IN_PLACE, lim_iter_peak, 1, MPI_INTEGER, MPI_MAX, p_comm_worker, p_err)" in flow


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


def test_bif_generic_tracers_are_not_rejected_before_transport():
    namelist = (ROOT / "share/MOD_Namelist.F90").read_text(encoding="utf-8")
    sediment = (ROOT / "main/TRACER/MOD_Tracer_Particle_Sediment.F90").read_text(encoding="utf-8")
    assert "IF (DEF_USE_BIFURCATION .and. DEF_TRACER_NUM > 0)" not in namelist
    # Particle sediment now has its own BIF transport path.
    assert "sediment bifurcation transport is not yet implemented" not in sediment
    assert "CALL apply_bif_sediment_credits" in sediment


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


def test_limiter_tolerance_is_the_one_the_replica_verifies():
    # The Python replica above only proves anything about the shipped loop if it
    # stops at the same tolerance.  The tolerance is absolute: rates live in
    # [0,1], so a relative scale would always be 1 and only hide the loosening.
    source = RIVER.read_text(encoding="utf-8")
    shipped = float(
        re.search(r"TRC_LIMITER_RATE_TOL\s*=\s*([0-9.]+e-?[0-9]+)_r8", source).group(1)
    )
    assert shipped == inspect.signature(_safe_fixed_point).parameters["tol"].default
    assert "TRC_LIMITER_RATE_TOL_REL" not in source
    assert "limiter_rate_scale" not in source
    assert source.count("limiter_delta_global <= TRC_LIMITER_RATE_TOL") == 2


def test_fast_path_is_exact_for_heterogeneous_visible_and_protected_pools():
    # Pools 0/1 are visible, 2 is behind a levee. The first edge reverses
    # the nominal main channel; the others are opposed BIF layer transfers.
    # Every donor's gross WATER outflow is below its own initial storage.
    volumes = [2.0, 1.5, 0.4]
    edges_water = [(1, 0, 0.3), (0, 2, 0.5), (2, 1, 0.2)]
    outgoing_water = [sum(q for src, _, q in edges_water if src == i)
                      for i in range(len(volumes))]
    assert all(out <= vol for out, vol in zip(outgoing_water, volumes))

    # The two dissolved tracers deliberately have different, nonuniform
    # signatures, including a low-concentration protected compartment.
    for concentrations in ([0.2, 0.9, 0.03], [1.7, 0.01, 0.8]):
        mass = [vol * conc for vol, conc in zip(volumes, concentrations)]
        edges = [(src, dst, q * concentrations[src])
                 for src, dst, q in edges_water]
        outgoing = [sum(amount for src, _, amount in edges if src == i)
                    for i in range(len(volumes))]
        assert all(out <= held for out, held in zip(outgoing, mass))
        assert _safe_fixed_point(mass, outgoing, edges) == [1.0] * len(mass)
        final = mass.copy()
        for src, dst, amount in edges:
            final[src] -= amount
            final[dst] += amount
        assert min(final) >= 0.0
        assert isclose(sum(final), sum(mass), rel_tol=0.0, abs_tol=1e-15)


def test_dry_multi_face_donor_requires_global_fallback():
    # A dry cell's queued tracer uses an edge-volume concentration, not
    # mass/actual-water-storage; two outgoing faces can each request it.
    # Rank 0 sees only the safe donor, while rank 1 must enter Jacobi.
    mass = [1.0, 0.2, 0.0, 0.0]
    edges = [(0, 1, 0.4), (1, 2, 0.2), (1, 3, 0.2)]
    outgoing = [0.4, 0.4, 0.0, 0.0]
    initial_rates = [min(m / out, 1.0) if out else 1.0
                     for m, out in zip(mass, outgoing)]
    assert initial_rates == [1.0, 0.5, 1.0, 1.0]
    assert max(1.0 - initial_rates[i] for i in (0, 2, 3)) == 0.0
    assert max(1.0 - initial_rates[i] for i in (1,)) == 0.5
    rates = _safe_fixed_point(mass, outgoing, edges)
    assert rates == [1.0] * len(mass)  # same-step incoming credit resolves it
    # Without the all-rank branch decision, rank 0 would skip the rate push
    # and rank 1 would deadlock, even though the final fixed point is unity.
    source = RIVER.read_text(encoding="utf-8")
    initial = source.split("Skip communication entirely when T(0)==1 everywhere", 1)[1].split(
        "DO limiter_iter = 1, limiter_max_iter", 1
    )[0]
    assert "mpi_allreduce(limiter_delta, limiter_delta_global" in initial
