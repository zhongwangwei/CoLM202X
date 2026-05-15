#!/usr/bin/env python3
"""Compare levee-on vs levee-off Amazon 2-yr unitcat history.

Focuses on flood fraction (f_floodfrc), flooded area (f_floodarea),
floodplain storage (f_fldsto), floodplain depth (f_flddph), river storage
(f_rivsto), and (when available) levee-protected storage (f_levsto).

Outputs:
  - per-month domain-mean time series (CSV + plot)
  - 2-yr mean spatial difference map (PNG)
  - summary stats to stdout
"""
import glob, os, sys
import numpy as np
import netCDF4

ON_DIR  = "/media/zhwei/data02/zhwei/levee/history"
OFF_DIR = "/media/zhwei/data02/zhwei/levee_off/history"
OUT_DIR = "/media/zhwei/data02/zhwei/levee/compare_levee_onoff"
os.makedirs(OUT_DIR, exist_ok=True)

# Amazon study window from test.nml
LAT_MIN, LAT_MAX = -15.0, 5.0
LON_MIN, LON_MAX = -72.0, -48.0

VARS = ["f_floodfrc", "f_floodarea", "f_flddph", "f_fldsto", "f_rivsto",
        "f_wdpth_ucat", "f_levsto", "f_levdph"]


def load_pair(month_tag):
    on_path  = f"{ON_DIR}/Amazon_levee_test60_hist_unitcat_{month_tag}.nc"
    off_path = f"{OFF_DIR}/Amazon_levee_test60_hist_unitcat_{month_tag}.nc"
    if not (os.path.exists(on_path) and os.path.exists(off_path)):
        return None, None
    return netCDF4.Dataset(on_path), netCDF4.Dataset(off_path)


def domain_slice(ds):
    lat = ds.variables["lat_ucat"][:]
    lon = ds.variables["lon_ucat"][:]
    j = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
    i = np.where((lon >= LON_MIN) & (lon <= LON_MAX))[0]
    return j, i, lat[j], lon[i]


def domain_mean(ds, var, j, i):
    if var not in ds.variables:
        return np.nan, 0
    a = ds.variables[var][0, :, :]
    sub = a[np.ix_(j, i)]
    if hasattr(sub, "mask"):
        valid = sub.compressed()
    else:
        valid = sub[np.isfinite(sub)]
    if valid.size == 0:
        return np.nan, 0
    return float(valid.mean()), int(valid.size)


# ---------- 1. monthly time series ----------
months = []
for year in (1995, 1996):
    for m in range(1, 13):
        months.append(f"{year}-{m:02d}")

rows = []
for mt in months:
    on, off = load_pair(mt)
    if on is None:
        rows.append((mt, *([np.nan]*len(VARS)*2)))
        continue
    j, i, _, _ = domain_slice(on)
    rec = [mt]
    for v in VARS:
        on_m, n_on  = domain_mean(on,  v, j, i)
        off_m, n_off = domain_mean(off, v, j, i)
        rec.extend([on_m, off_m])
    rows.append(tuple(rec))
    on.close(); off.close()

# csv
csv_path = f"{OUT_DIR}/monthly_domain_mean.csv"
header = "month," + ",".join(f"{v}_on,{v}_off" for v in VARS)
with open(csv_path, "w") as fh:
    fh.write(header + "\n")
    for r in rows:
        fh.write(",".join(str(x) for x in r) + "\n")
print(f"wrote {csv_path}")

# ---------- 2. 2-yr mean spatial difference ----------
def two_year_mean(srcdir, var):
    acc = None; cnt = 0
    for mt in months:
        p = f"{srcdir}/Amazon_levee_test60_hist_unitcat_{mt}.nc"
        if not os.path.exists(p):
            continue
        ds = netCDF4.Dataset(p)
        if var not in ds.variables:
            ds.close(); continue
        a = ds.variables[var][0, :, :]
        if acc is None:
            acc = np.ma.zeros_like(a, dtype=np.float64)
            mask = np.ma.getmaskarray(a).copy()
        acc = acc + a.astype(np.float64)
        mask &= np.ma.getmaskarray(a)
        cnt += 1
        ds.close()
    if acc is None:
        return None, None, None
    mean = acc / cnt
    if isinstance(mean, np.ma.MaskedArray):
        mean.mask = mask
    return mean, mask, cnt


# Just print summary numbers, not plots (avoids matplotlib dep noise).
print()
print("=== 2-year domain-mean summary (Amazon: lat[-15,5], lon[-72,-48]) ===")
print(f"{'var':<14} {'levee_on':>12} {'levee_off':>12} {'on-off':>12} {'rel_pct':>10}")
arr = np.array([row[1:] for row in rows], dtype=float)
for k, v in enumerate(VARS):
    on_col  = arr[:, 2*k]
    off_col = arr[:, 2*k + 1]
    on_mean = np.nanmean(on_col)
    off_mean = np.nanmean(off_col)
    diff = on_mean - off_mean
    rel = 100.0 * diff / off_mean if (off_mean and np.isfinite(off_mean) and off_mean != 0) else np.nan
    print(f"{v:<14} {on_mean:>12.5g} {off_mean:>12.5g} {diff:>+12.5g} {rel:>+9.3f}%")

# ---------- 3. peak flood months (look where ON-OFF is largest) ----------
print()
print("=== floodfrc by month (domain mean): on / off / on-off ===")
k_ff = VARS.index("f_floodfrc")
for r, mt in zip(rows, months):
    on_m  = r[1 + 2*k_ff]
    off_m = r[1 + 2*k_ff + 1]
    if np.isfinite(on_m) and np.isfinite(off_m):
        print(f"  {mt}  on={on_m:.5f}  off={off_m:.5f}  diff={on_m-off_m:+.5f}")

# ---------- 4. spatial diff for f_floodfrc, save as NetCDF ----------
mean_on,  _, _ = two_year_mean(ON_DIR,  "f_floodfrc")
mean_off, _, _ = two_year_mean(OFF_DIR, "f_floodfrc")
if mean_on is not None and mean_off is not None:
    diff = mean_on - mean_off
    out_nc = f"{OUT_DIR}/floodfrc_2yr_diff.nc"
    # need lat/lon coords; pull from any month file
    src = netCDF4.Dataset(f"{ON_DIR}/Amazon_levee_test60_hist_unitcat_1995-01.nc")
    lat = src.variables["lat_ucat"][:]; lon = src.variables["lon_ucat"][:]
    src.close()
    with netCDF4.Dataset(out_nc, "w") as ds:
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        v_lat = ds.createVariable("lat", "f8", ("lat",)); v_lat[:] = lat; v_lat.units="degrees_north"
        v_lon = ds.createVariable("lon", "f8", ("lon",)); v_lon[:] = lon; v_lon.units="degrees_east"
        v_on  = ds.createVariable("floodfrc_on",  "f4", ("lat","lon")); v_on[:]  = mean_on
        v_off = ds.createVariable("floodfrc_off", "f4", ("lat","lon")); v_off[:] = mean_off
        v_d   = ds.createVariable("floodfrc_diff","f4", ("lat","lon")); v_d[:]   = diff
        v_d.long_name = "2-yr-mean floodfrc difference (LEVEE_on minus LEVEE_off)"
    print(f"wrote {out_nc}")

    # Amazon-window stats on diff
    j = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
    i = np.where((lon >= LON_MIN) & (lon <= LON_MAX))[0]
    sub = diff[np.ix_(j, i)]
    if hasattr(sub, "compressed"):
        s = sub.compressed()
    else:
        s = sub[np.isfinite(sub)]
    if s.size > 0:
        print()
        print("=== Amazon-window 2-yr mean floodfrc difference (on-off) ===")
        print(f"  count of valid cells : {s.size}")
        print(f"  mean diff            : {s.mean():+.6f}")
        print(f"  median diff          : {np.median(s):+.6f}")
        print(f"  min  diff            : {s.min():+.6f}")
        print(f"  max  diff            : {s.max():+.6f}")
        print(f"  RMS  diff            : {np.sqrt(np.mean(s**2)):+.6f}")
        print(f"  cells |diff|>1e-3    : {np.sum(np.abs(s)>1e-3)}")
        print(f"  cells |diff|>1e-2    : {np.sum(np.abs(s)>1e-2)}")
