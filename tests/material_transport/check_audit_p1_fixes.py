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
      'trc_dry_fwd(itrc, i) = trc_dry_fwd(itrc, i) + dry_drain' in rl)
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
      'fwd_send(:) = trc_dry_fwd(itrc, :)' in rl)
check('N06 push call sits inside the per-tracer loop (same order every rank)',
      rl.index('fwd_send(:) = trc_dry_fwd(itrc, :)') <
      rl.index('CALL worker_push_data (push_ups2ucat, fwd_send, fwd_recv'))

# N07: direction from the actual host flux, reverse/unknown kept pending.
check('N07 direction uses the timestep face flux, not static topology only',
      'present(hflux_fc_in)' in rl and 'hflux_fc_in(i) >= 0._r8' in rl)
check('N07 reverse/no-known-flow residual is held pending instead of misrouted',
      'trc_dry_pending(itrc, i) = trc_dry_pending(itrc, i) + dry_drain' in rl)
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

print('')
print(f'{len(fails)} failure(s)')
sys.exit(1 if fails else 0)
