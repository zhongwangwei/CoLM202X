#!/usr/bin/env python3
"""Trace δ18O drift through canopy → soil/snow column over the 5-day tracer test.

Under the Phase-1 invariant (R_atm = R_init = init_delta everywhere), all
pools should sit at exactly -10‰ for δ18O for all time. Any drift is a
bookkeeping bug. This script reports day-by-day δ statistics for each pool
so we can localise the source.

Usage: python3 trace_tracer_chain.py
"""
import os, glob
import numpy as np
import netCDF4

HIST = "/media/zhwei/data02/zhwei/tracer_test/Amazon_levee_test60/history"
DAYS = [1, 2, 3, 4, 5]
INIT_DELTA = -10.0  # δ18O init
Rsmow_18O = 2.0052e-3
TRACER = "H2_18O"

# For numerical stability: only count cells with finite, positive R, and
# water-mass weight large enough that R noise stays bounded. Without this
# the tail of nearly-empty cells dominates.
TINY_R = 1e-12


def load_R(ds, name):
    if name not in ds.variables:
        return None
    a = ds.variables[name][:]
    return a


def stats(R):
    """Return (n, median, mean_abs_drift_‰, q05, q95) for finite R>0."""
    R = np.asarray(R, dtype=np.float64)
    if hasattr(R, "compressed"):
        R = R.compressed()
    finite = np.isfinite(R) & (R > TINY_R)
    R = R[finite]
    if R.size == 0:
        return 0, np.nan, np.nan, np.nan, np.nan
    delta = (R / Rsmow_18O - 1.0) * 1000.0
    drift = delta - INIT_DELTA
    return R.size, float(np.median(delta)), float(np.mean(drift)), float(np.percentile(delta, 5)), float(np.percentile(delta, 95))


print(f"=== δ18O drift vs init {INIT_DELTA}‰ — pool-by-pool, day-by-day ===")
print("(median > init means soil enriching in heavy isotope = light leak)\n")

pools_2d = ["ldew", "wa", "wdsrf", "wetwat", "scv"]
header = f"{'pool':<12s} {'day':<5s} {'n':>6s} {'median(‰)':>11s} {'mean_drift':>11s} {'q05(‰)':>9s} {'q95(‰)':>9s}"

for pool in pools_2d:
    print(header)
    for day in DAYS:
        fn = f"{HIST}/Amazon_levee_test60_hist_1995-01-{day:02d}.nc"
        if not os.path.exists(fn):
            continue
        ds = netCDF4.Dataset(fn)
        R = load_R(ds, f"f_trc_conc_{pool}_{TRACER}")
        ds.close()
        if R is None:
            continue
        n, med, mn, q05, q95 = stats(R)
        print(f"{pool:<12s} {day:<5d} {n:>6d} {med:>+11.4f} {mn:>+11.4f} {q05:>+9.3f} {q95:>+9.3f}")
    print()

# Soil/snow column: split by layer to localise vertical drift
print("=== soisno per-layer δ18O (layer 1=top soil, ..., 10=deepest) ===")
print(header)
for day in DAYS:
    fn = f"{HIST}/Amazon_levee_test60_hist_1995-01-{day:02d}.nc"
    if not os.path.exists(fn):
        continue
    ds = netCDF4.Dataset(fn)
    R3d = load_R(ds, f"f_trc_conc_soisno_{TRACER}")
    ds.close()
    # shape (1, 15, 80, 96) — 15 = 5 snow + 10 soil; layer index has snow first?
    # check first day
    if R3d is None:
        continue
    shape = R3d.shape
    nlayer = shape[-3]
    for k in range(nlayer):
        Rk = R3d[0, k, :, :]
        n, med, mn, q05, q95 = stats(Rk)
        if n == 0:
            continue
        print(f"L{k+1:<11d} {day:<5d} {n:>6d} {med:>+11.4f} {mn:>+11.4f} {q05:>+9.3f} {q95:>+9.3f}")
    print()

# River-side
print("=== river-side (unitcat) δ18O ===")
print(header)
for day in DAYS:
    fn = f"{HIST}/Amazon_levee_test60_hist_unitcat_1995-01-{day:02d}.nc"
    if not os.path.exists(fn):
        continue
    ds = netCDF4.Dataset(fn)
    R = load_R(ds, f"f_trc_conc_{TRACER}")
    ds.close()
    if R is None:
        continue
    n, med, mn, q05, q95 = stats(R)
    print(f"{'river':<12s} {day:<5d} {n:>6d} {med:>+11.4f} {mn:>+11.4f} {q05:>+9.3f} {q95:>+9.3f}")
