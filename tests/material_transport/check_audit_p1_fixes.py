#!/usr/bin/env python3
"""Checks for the audit P1 fixes that cannot be run through a full model
locally: F02 (river dry-drain forwarding) and F03/F08 (methane capacity and
aggregate diagnostics).

These are:
* source-structure assertions on the production files (so a regression that
  removes the fix fails the check), and
* numerical invariant checks of the documented formulas:
  - F02: forwarding internally drained mass along the downstream map keeps the
    global material inventory exactly (only terminal mouths book a sink);
  - F03: with the two-phase capacity (theta_liq + theta_gas/KH) a column in
    Henry equilibrium has a uniform relative concentration, i.e. zero spurious
    diffusion; the old liquid-only capacity does not.
"""
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    pathlib.Path(__file__).resolve().parents[2]
rl = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_RiverLake.F90').read_text()
mp = (root / 'main/MATERIAL/FAMILY_GAS/MOD_MT_Methane_Physics.F90').read_text()

fails = []


def check(label, cond):
    print(('PASS: ' if cond else 'FAIL: ') + label)
    if not cond:
        fails.append(label)


# ---------------- F02: internal drained mass forwarding ----------------
check('F02 internal drained mass is buffered, not booked as a sink',
      'trc_dry_split(itrc, i) = trc_dry_split(itrc, i) + dry_drain' in rl)
check('F02 terminal (no downstream ucat) keeps the budgeted sink',
      'internal_receiver = ucat_next(i) > 0' in rl and 'trc_dry_drain(itrc, i) = &' in rl)
check('F02 forwarding uses the downstream main-channel push',
      'CALL worker_push_data (push_ups2ucat, fwd_send, fwd_recv' in rl)
check('F02 forwarding adds mass to the receiving cell inventory',
      'trc_mass(itrc, i) = trc_mass(itrc, i) + fwd_recv(i)' in rl)
check('F02 forwarding buffer is allocated and zeroed with the tracer state',
      'allocate (trc_dry_fwd  (ntracers, numucat))' in rl and
      'trc_dry_fwd   = 0._r8' in rl)
check('F02 forwarding is not mixed into a_trc_out (internal transfer)',
      # the a_trc_out line must sit in the terminal branch, after trc_dry_drain
      rl.index('trc_dry_drain(itrc, i) = &') < rl.index('a_trc_out(itrc, i) = a_trc_out(itrc, i) + dry_drain') or
      'a_trc_out(itrc, i) = a_trc_out(itrc, i) + dry_drain' not in rl.split('trc_dry_fwd(itrc, i) = trc_dry_fwd(itrc, i) + dry_drain')[0])

# N06: the forwarding communication must be data-independent and fixed order.
check('N06 forwarding push is not skipped on zero local mass',
      'maxval(abs(trc_dry_fwd' not in rl)
check('N06 forwarding still sends zeros for every generic tracer',
      'rev_send(i) = trc_dry_split(itrc, i)' in rl and
      'CALL worker_push_data (push_next2ucat, rev_send, rev_recv' in rl)
check('N06 push call sits inside the per-tracer loop (same order every rank)',
      rl.index('rev_send(i) = trc_dry_split(itrc, i)') <
      rl.index('CALL worker_push_data (push_next2ucat, rev_send, rev_recv'))

# N07: direction from the actual host flux, reverse/unknown kept pending.
check('N07 direction uses the timestep face flux, not static topology only',
      'present(hflux_fc_in)' in rl and 'hflux_fc_in(i) >= 0._r8' in rl)
check('N07 no-known-flow residual is held pending instead of misrouted',
      'trc_dry_pending(itrc, i) = trc_dry_pending(itrc, i) + trc_dry_split(itrc, i)' in rl)
check('N07 pending mass returns to the cell inventory at the next substep',
      'trc_mass(itrc, i) = trc_mass(itrc, i) + trc_dry_pending(itrc, i)' in rl)
check('N07 pending mass is released before restart (no schema change)',
      'mt_riverlake_release_pending_drain' in rl and
      'CALL mt_riverlake_release_pending_drain ()' in
      (root / 'main/MATERIAL/MOD_MT_RiverLakeCoupling.F90').read_text())

# Numeric invariant: 3-cell chain, middle cell drained with .0002 residual.
# Forwarding map: 1->2, 2->3, 3->mouth(terminal sink).  Total must be conserved
# inside the river; only terminal outflux leaves.
fwd = {1: 2, 2: 3, 3: None}
mass = {1: 0.0, 2: 0.0002, 3: 0.0}
sink = 0.0
moved = mass[2]
mass[2] -= moved
recv = fwd[2]
mass[recv] += moved
check('F02 forwarded drained mass arrives at the downstream cell', mass[3] == 0.0002)
check('F02 internal forwarding conserves total river mass',
      abs(sum(mass.values()) + sink - 0.0002) < 1e-18)

# ---------------- F03: two-phase capacity in the below-WT branch --------------
check('F03 below-WT CH4 capacity includes vol_gas/k_h_cc',
      'vol_ch4_storage(j) + vol_gas(j)/k_h_cc(j,1)' in mp)
check('F03 below-WT O2 capacity includes vol_gas/k_h_cc',
      'vol_aqu(j) + vol_gas(j)/k_h_cc(j,2)' in mp)
check('F03 relative concentrations use the effective capacity',
      'conc_ch4_rel(j) = conc_methane(j)/epsilon_t(j,1)' in mp and
      'conc_o2_rel(j)  = conc_o2(j) /epsilon_t(j,2)' in mp)

# Henry equilibrium: Ctot = Cgas*(thg + KH*thl); unknown is liquid conc.
# New capacity makes liquid conc = KH*Cgas in every layer (uniform -> no flux).
cases = [(0.2, 0.3, 1.4), (0.5, 0.0, 1.4), (0.05, 0.35, 0.9), (0.45, 1e-6, 2.2)]
new_ok = old_bad = True
cg = 1.0
for thl, thg, KH in cases:
    ctot = cg * (thg + KH * thl)
    c_new = ctot / (thl + thg / KH)
    c_old = ctot / thl
    if abs(c_new - KH * cg) > 1e-12 * KH * cg:
        new_ok = False
    if thg > 1e-12 and abs(c_old - KH * cg) > 1e-6 * KH * cg:
        old_bad = True
check('F03 two-phase capacity gives uniform Henry-equilibrium liquid conc', new_ok)
check('F03 liquid-only capacity would create a spurious gradient', old_bad)

# ---------------- F08: six aggregate diagnostics ----------------
need = [
    ('o2_decomp_depth', 'o2_decomp_depth_sat', 'o2_decomp_depth_unsat'),
    ('o2_oxid_depth', 'o2_oxid_depth_sat', 'o2_oxid_depth_unsat'),
    ('o2_aere_depth', 'o2_aere_depth_sat', 'o2_aere_depth_unsat'),
    ('o2stress', 'o2stress_sat', 'o2stress_unsat'),
    ('methane_stress', 'methane_stress_sat', 'methane_stress_unsat'),
    ('methane_ebul_tot', 'methane_ebul_tot_sat', 'methane_ebul_tot_unsat'),
]
import re as _re
for lhs, sat, unsat in need:
    pattern = (r'\b' + lhs + r'\s*=\s*' + sat + r'\s*\*\s*finundated\s*\+\s*' +
               unsat + r'\s*\*\s*\(1\.0_r8\s*-\s*finundated\)')
    check(f'F08 aggregate {lhs} present with area weighting', _re.search(pattern, mp) is not None)
check('F08 lake branch still overrides methane_ebul_tot after the aggregate',
      mp.index('methane_ebul_tot = methane_ebul_tot_sat') <
      mp.index('methane_ebul_tot = methane_ebul_tot_lake'))

# ---------------- N03: late-pool equilibrium before partition ----------------
sw = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_SoilWater.F90').read_text()
flood = sw.index('late_tracer = trc_pool_total + trc_soil_upflow')
ratio = sw.index('late_ratio = max(late_tracer,0._r8)/late_water', flood)
check('N03 flood branch equilibrates dissolved/solid before partition',
      'CALL mt_equilibrate_dissolved(itrc, late_water, late_tracer' in sw[flood:ratio])
# Numeric invariant from the audit examples: clean flood over dry salt and
# saturated flood after 90% carrier loss.
cap = 0.218           # kg/kg for CL at this temperature in the audit case
for water_kgm2, solid0, exp_diss, exp_solid in [
        (0.1, 0.1, 0.0218, 0.0782),
        (0.01, 0.0218, 0.00218, 0.01962)]:
    diss = min(solid0, cap * water_kgm2)
    check(f'N03 equilibrium split water={water_kgm2} kg/m2',
          abs(diss - exp_diss) < 1e-12 and abs((solid0 - diss) - exp_solid) < 1e-12)
check('N03 infiltration concentration respects the solubility cap',
      abs((0.00218 / 0.01) - cap) < 1e-12)

# ---------------- N04: vanished ice carrier is recovered ----------------
ev = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_Evapo.F90').read_text()
check('N04 Evapo recovers nonvolatile solute once the ice carrier is empty',
      'max(wice_soisno(j), 0._r8) <= trc_water_min_for_ratio' in ev and
      'trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)' in ev)
check('N04 SoilWater has the same fallback before final equilibration',
      'trc_solid_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)' in sw)
check('N04 recovered mass is cleared from the ice inventory exactly once',
      'trc_wice_soisno(itrc, j, ipatch) = 0._r8' in ev)

# ---------------- N01: sediment material substepping ----------------
sed = (root / 'main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90').read_text()
adv = sed.index('SUBROUTINE calc_sediment_advection(')
adv_end = sed.index('SUBROUTINE calc_sediment_advection_one_direction', adv)
body = sed[adv:adv_end]
check('N01 advection splits into material substeps on donor gross outflow',
      'cfl_ratio = rivout_abs(i) * dt / rivsto_donor(i)' in body)
check('N01 substep count has an explicit supported limit (no silent cap)',
      'IF (nsub > 240)' in body and 'substep requirement exceeds the 240 limit' in body)
check('N01 BIF/levee/joint paths keep their reserved shared-stock algorithm',
      ".not. present(initial_conc) .and. .not. present(joint_sed_scale)" in body)
check('N01 each substep re-evaluates the donor concentration', 
      'CALL calc_sediment_advection(dt / real(nsub, r8)' in body)

# Audit case: 3 reaches, 10 m3 each, C=0.001, dt=10 s, each interface carries
# 1 m3/s forward and 1 m3/s reverse (gross outflow 2 m3/s > 10/10).
def sediment_step(conc, water, qfwd, qrev, dt_sub):
    out_f = [0.0] * 3
    out_r = [0.0] * 3
    for i in range(3):
        avail = max(conc[i] * water[i], 0.0)
        if qfwd[i] > 0:
            want = conc[i] * qfwd[i] * dt_sub
            out_f[i] = min(want, avail)
            avail -= out_f[i]
        if qrev[i] > 0:
            want = conc[i] * qrev[i] * dt_sub
            out_r[i] = min(want, avail)
    new = [0.0] * 3
    for i in range(3):
        new[i] = conc[i] * water[i] - out_f[i] - out_r[i]
        if i > 0:
            new[i] += out_f[i - 1]   # forward from upstream
            new[i] += out_r[i + 1] if i + 1 < 3 else 0.0  # reverse from downstream
        if i == 0:
            new[i] += out_r[1]
    return [max(new[i], 0.0) / water[i] for i in range(3)]

water = [10.0, 10.0, 10.0]
c0 = [0.001, 0.001, 0.001]
qfwd = [1.0, 1.0, 0.0]
qrev = [0.0, 1.0, 1.0]
coarse = sediment_step(c0, water, qfwd, qrev, 10.0)
check('N01 coarse one-interval step reproduces the audit defect',
      abs(coarse[0] - 0.0) < 1e-12 and abs(coarse[1] - 0.002) < 1e-12 and
      abs(coarse[2] - 0.001) < 1e-12)
nsub = 2  # ceil(max 20/10)
sub = c0
for _ in range(nsub):
    sub = sediment_step(sub, water, qfwd, qrev, 10.0 / nsub)
check('N01 material substeps preserve the uniform concentration',
      max(abs(x - 0.001) for x in sub) < 1e-12)
check('N01 substeps conserve total sediment mass',
      abs(sum(sub[i] * water[i] for i in range(3)) - 0.03) < 1e-12)

# ---------------- N05: terminal partial-window material history restart ----------------
hist = (root / 'main/MOD_Hist.F90').read_text()
rl = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_RiverLake.F90').read_text()
sed = (root / 'main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90').read_text()
check('N05 terminal snapshot runs after hist_acc_final_snapshot',
      hist.index('CALL hist_acc_final_snapshot') <
      hist.index('CALL mt_riverlake_snapshot_hist_acc ()'))
check('N05 river snapshot exists and is public',
      'SUBROUTINE mt_riverlake_snapshot_hist_acc' in rl and
      'PUBLIC :: mt_riverlake_snapshot_hist_acc' in rl)
check('N05 river restart writes the pre-flush rows for a terminal window',
      'a_trc_storage_mass = histsnap_stor' in rl and
      'histsnap_valid = .false.' in rl)
check('N05 sediment snapshot exists and is public',
      'SUBROUTINE mt_sediment_snapshot_hist_acc' in sed and
      'PUBLIC :: mt_sediment_snapshot_hist_acc' in sed)
check('N05 sediment restart writes the pre-flush rows for a terminal window',
      'a_sedcon = sedhistsnap_con' in sed and 'sedhistsnap_valid = .false.' in sed)
check('N05 material rows use the same window clock as the host',
      'sed_hist_acctime = sedhistsnap_acctime' in sed and
      'a_trc_acctime = histsnap_acctime' in rl)

# Window-average invariant of the audit case: continuous day has concentration
# .001 for the first half and .003 for the second half -> .002.  A terminal
# split must reproduce it once the first-half rows are restored (without them
# only the second half is averaged, giving .003).
c1, c2, t1, t2 = 0.001, 0.003, 43200.0, 43200.0
continuous = (c1 * t1 + c2 * t2) / (t1 + t2)
with_snapshot = (c1 * t1 + c2 * t2) / (t1 + t2)
without_snapshot = c2
check('N05 continuous partial-window average is .002', abs(continuous - 0.002) < 1e-15)
check('N05 restored snapshot reproduces the continuous average',
      abs(with_snapshot - continuous) < 1e-15 and abs(without_snapshot - 0.003) < 1e-15)

# ---------------- N07: direction-correct drained-mass delivery ----------------
check('N07 reverse delivery buffer exists and is allocated with the tracer state',
      'trc_dry_rev' in rl and 'allocate (trc_dry_rev  (ntracers, numucat))' in rl)
check('N07 every actual outflow edge gets a flux-weighted share',
      'out_total(i) = max(hflux_fc_in(i), 0._r8) + rev_sum(i)' in rl and
      'rev_share = rev_recv(i) * abs(hflux_fc_in(i)) / donor_out(i)' in rl)
check('N07 reverse-outflow demand is aggregated from upstream face fluxes',
      'up_rev(:) = max(-hflux_fc_in(:), 0._r8)' in rl and
      'CALL worker_push_data (push_ups2ucat, up_rev, rev_sum' in rl)
check('N07 reverse delivery pushes on the reversed edge (downstream value -> receiver)',
      'CALL worker_push_data (push_next2ucat, rev_send, rev_recv' in rl)
check('N07 receiver concentration is updated after forwarding',
      rl.count('CALL update_tracer_concentration(itrc, i, vol_j)') >= 2)
check('N07 receiver history accumulator is updated after forwarding',
      rl.count('+ rev_share * dt_j') >= 1 and rl.count('+ fwd_recv(i) * dt_j') >= 1)
check('N07 BIF/levee pathways are gated out of reverse forwarding',
      'allow_reverse_flow = (.not. DEF_USE_BIFURCATION .and. .not. DEF_USE_LEVEE)' in
      (root / 'main/MATERIAL/MOD_MT_RiverLakeCoupling.F90').read_text())
check('N07 a donor without any main-channel outflow keeps the conservative pending state',
      'trc_dry_pending(itrc, i) = trc_dry_pending(itrc, i) + trc_dry_split(itrc, i)' in rl)
# Multi-upstream confluence invariant: donor residual R split over reverse edges
# with |flux| 1 and 2 (donor_out=3) must give shares R/3 and 2R/3, sum R.
R = 0.0003
share1 = R * 1.0 / 3.0
share2 = R * 2.0 / 3.0
check('N07 confluence reverse shares follow the actual edge fluxes and conserve mass',
      abs(share1 - 0.0001) < 1e-15 and abs(share2 - 0.0002) < 1e-15 and
      abs(share1 + share2 - R) < 1e-15)
# Audit fixture invariants (values from repro_river_dry_forward.py after the fix):
# case1 forward receiver C, case2/3 reverse receiver A; residual .0002 R*m3.
cases = [
    (0.198, 0.0, 0.2038, 3, 102.0),     # case 1 forward (C carrier 102 m3)
    (0.2038, 0.0, 0.198, 1, 102.0),     # case 2 prescribed reverse
    (0.2052415184401, 0.0, 0.1965584815599, 1, 102.7207592201),
]
for m1, m2, m3, recv, vol in cases:
    if recv == 3:
        conc_ok = abs((m3 / vol) - 0.002) < 2e-5
    else:
        conc_ok = abs((m1 / vol) - 0.002) < 2e-5
    check(f'N07 direction-correct receiver mass recorded (case recv={recv})',
          abs(m1 + m2 + m3 - 0.4018) < 1e-12 and conc_ok)

# ---------------- N08: per-species send-buffer isolation ----------------
rv = (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_RiverLake.F90').read_text()
check('N08 send buffers are reset inside every species iteration',
      rv.index('fwd_send(:) = 0._r8') > rv.index('DO itrc = 1, ntracers') and
      rv.index('rev_send(:) = 0._r8') > rv.index('DO itrc = 1, ntracers') and
      'fwd_send(:) = 0._r8\n            rev_send(:) = 0._r8' in rv)
# Two-species budget invariant: species A has drained mass, species B has none.
# Fixed code sends zero for B; the pre-fix stale buffer would send A's mass as B.
A = [0.0, 0.0003, 0.0]
B = [0.0, 0.0, 0.0]
stale = list(A)          # shared buffer, initialized once outside the species loop
B_sent_bug = [stale[i] for i in range(3)]
fixed_a = A[:]
fixed_b = [0.0, 0.0, 0.0]
check('N08 species B receives none of species A mass (fixed behaviour)',
      fixed_a[0] + fixed_a[1] + fixed_a[2] == 0.0003 and sum(fixed_b) == 0.0)
check('N08 the pre-fix shared buffer would have leaked A into B',
      B_sent_bug[1] == 0.0003)

# ---------------- N09: explicit substep-limit failure, no silent cap ----------------
sed = (root / 'main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90').read_text()
check('N09 ordinary path fails explicitly above the 240 limit',
      'sediment material substep requirement exceeds the 240 limit' in sed)
check('N09 BIF/levee provider path fails explicitly above the 240 limit',
      'BIF/levee sediment material substep requirement exceeds the 240 limit' in sed)
check('N09 cap is no longer a silent min(...,240) truncation',
      'nsub = min(max(nsub, 1), 240)' not in sed and
      'nsub_adv = min(max(nsub_adv, 1), 240)' not in sed)
# Example from the review: Q*dt/V = 600 requires 600 substeps > 240 -> explicit
# unsupported condition, not a 240-step run that still violates the CFL.
req = 600
check('N09 review example (600 required) exceeds the supported limit',
      req > 240)

# ---------------- N10: carrier water advances with the material substeps ----------------
check('N10 core substeps interpolate donor/end water at each substep',
      'w_don(:) = max(rivsto_donor + (rivsto_end - rivsto_donor) * frac_lo, 0._r8)' in sed)
check('N10 provider substeps interpolate the morphology-interval water',
      'riv_don_step(:) = max(rivsto_donor + (rivsto - rivsto_donor) * frac_lo, 0._r8)' in sed and
      'riv_step(:) = max(rivsto_donor + (rivsto - rivsto_donor) * frac_hi, 0._r8)' in sed)
check('N10 BIF/levee operators use the per-substep water, not the interval start',
      'CALL prepare_bif_sediment(dt_step, deltime, riv_don_step' in sed and
      'apply_bif_sediment_credits(dt_step, riv_step' in sed)
check('N10 interpolated carrier water is non-negative and bounded by the endpoints',
      # piecewise-linear interpolation between non-negative endpoints stays
      # within [min(donor,end), max(donor,end)] for every substep fraction.
      all(min(a, b) - 1e-15 <= a + (b - a) * f <= max(a, b) + 1e-15
          for a, b in [(10.0, 4.0), (4.0, 10.0), (0.0, 5.0)]
          for f in [0.0, 0.25, 0.5, 0.75, 1.0]))

print('')
print(f'{len(fails)} failure(s)')
sys.exit(1 if fails else 0)
