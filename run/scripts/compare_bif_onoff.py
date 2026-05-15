#!/usr/bin/env python3
"""Compare bif-on (BIFURCATION=true, LEVEE=off) vs baseline (both off).

Uses the same Amazon-window subset as compare_levee_onoff.py."""
import os
import numpy as np
import netCDF4

ON_DIR  = "/media/zhwei/data02/zhwei/bif_on/history"
OFF_DIR = "/media/zhwei/data02/zhwei/levee_off/history"
OUT_DIR = "/media/zhwei/data02/zhwei/bif_on/compare_bif_onoff"
os.makedirs(OUT_DIR, exist_ok=True)

LAT_MIN, LAT_MAX = -15.0, 5.0
LON_MIN, LON_MAX = -72.0, -48.0

VARS = ["f_floodfrc", "f_floodarea", "f_flddph", "f_fldsto", "f_rivsto",
        "f_wdpth_ucat", "f_discharge", "f_veloc_riv"]


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


months = []
for year in (1995, 1996):
    for m in range(1, 13):
        months.append(f"{year}-{m:02d}")

rows = []
for mt in months:
    on_p = f"{ON_DIR}/Amazon_levee_test60_hist_unitcat_{mt}.nc"
    off_p = f"{OFF_DIR}/Amazon_levee_test60_hist_unitcat_{mt}.nc"
    if not (os.path.exists(on_p) and os.path.exists(off_p)):
        rows.append((mt, *([np.nan]*len(VARS)*2)))
        continue
    on  = netCDF4.Dataset(on_p);  off = netCDF4.Dataset(off_p)
    j, i, _, _ = domain_slice(on)
    rec = [mt]
    for v in VARS:
        on_m, _  = domain_mean(on,  v, j, i)
        off_m, _ = domain_mean(off, v, j, i)
        rec.extend([on_m, off_m])
    rows.append(tuple(rec))
    on.close(); off.close()

csv_path = f"{OUT_DIR}/monthly_domain_mean.csv"
header = "month," + ",".join(f"{v}_bifon,{v}_bifoff" for v in VARS)
with open(csv_path, "w") as fh:
    fh.write(header + "\n")
    for r in rows:
        fh.write(",".join(str(x) for x in r) + "\n")
print(f"wrote {csv_path}")

print()
print("=== 2-year domain-mean summary (Amazon) ===")
print(f"{'var':<14} {'bif_on':>12} {'bif_off':>12} {'on-off':>14} {'rel_pct':>10}")
arr = np.array([row[1:] for row in rows], dtype=float)
for k, v in enumerate(VARS):
    on_col  = arr[:, 2*k]
    off_col = arr[:, 2*k + 1]
    on_mean = np.nanmean(on_col)
    off_mean = np.nanmean(off_col)
    diff = on_mean - off_mean
    rel = 100.0 * diff / off_mean if (off_mean and np.isfinite(off_mean) and off_mean != 0) else np.nan
    print(f"{v:<14} {on_mean:>12.5g} {off_mean:>12.5g} {diff:>+14.5g} {rel:>+9.3f}%")

print()
print("=== floodfrc by month: bifon / bifoff / diff ===")
k_ff = VARS.index("f_floodfrc")
for r, mt in zip(rows, months):
    on_m  = r[1 + 2*k_ff]
    off_m = r[1 + 2*k_ff + 1]
    if np.isfinite(on_m) and np.isfinite(off_m):
        print(f"  {mt}  bif_on={on_m:.5f}  bif_off={off_m:.5f}  diff={on_m-off_m:+.5f}")


def two_year_mean(srcdir, var):
    acc = None; cnt = 0; mask = None
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
        acc += a.astype(np.float64)
        mask &= np.ma.getmaskarray(a)
        cnt += 1
        ds.close()
    if acc is None:
        return None
    mean = acc / cnt
    if isinstance(mean, np.ma.MaskedArray):
        mean.mask = mask
    return mean


for var_to_save in ("f_floodfrc", "f_discharge"):
    mean_on  = two_year_mean(ON_DIR,  var_to_save)
    mean_off = two_year_mean(OFF_DIR, var_to_save)
    if mean_on is None or mean_off is None:
        continue
    diff = mean_on - mean_off
    out_nc = f"{OUT_DIR}/{var_to_save}_2yr_diff.nc"
    src = netCDF4.Dataset(f"{ON_DIR}/Amazon_levee_test60_hist_unitcat_1995-01.nc")
    lat = src.variables["lat_ucat"][:]; lon = src.variables["lon_ucat"][:]
    src.close()
    with netCDF4.Dataset(out_nc, "w") as ds:
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        v_lat = ds.createVariable("lat", "f8", ("lat",)); v_lat[:] = lat; v_lat.units="degrees_north"
        v_lon = ds.createVariable("lon", "f8", ("lon",)); v_lon[:] = lon; v_lon.units="degrees_east"
        v_on  = ds.createVariable(f"{var_to_save}_bifon",  "f4", ("lat","lon")); v_on[:]  = mean_on
        v_off = ds.createVariable(f"{var_to_save}_bifoff", "f4", ("lat","lon")); v_off[:] = mean_off
        v_d   = ds.createVariable(f"{var_to_save}_diff",   "f4", ("lat","lon")); v_d[:]   = diff
        v_d.long_name = f"2-yr-mean {var_to_save} difference (BIF_on minus BIF_off)"
    print(f"wrote {out_nc}")

    # Amazon-window stats
    j = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
    i = np.where((lon >= LON_MIN) & (lon <= LON_MAX))[0]
    sub = diff[np.ix_(j, i)]
    if hasattr(sub, "compressed"):
        s = sub.compressed()
    else:
        s = sub[np.isfinite(sub)]
    if s.size > 0:
        print(f"  Amazon-box {var_to_save} 2-yr-mean diff: mean={s.mean():+.5g}  med={np.median(s):+.5g}  min={s.min():+.5g}  max={s.max():+.5g}  RMS={np.sqrt(np.mean(s**2)):+.5g}")
        print(f"  cells |diff|>0.001: {np.sum(np.abs(s)>1e-3)}, |diff|>0.01: {np.sum(np.abs(s)>1e-2)}")

# Top cells: where bifurcation has biggest discharge / floodfrc effect
ds = netCDF4.Dataset(f"{OUT_DIR}/f_floodfrc_2yr_diff.nc")
lat = ds.variables["lat"][:]; lon = ds.variables["lon"][:]
diff = ds.variables["f_floodfrc_diff"][:]
on  = ds.variables["f_floodfrc_bifon"][:]
off = ds.variables["f_floodfrc_bifoff"][:]
ds.close()
j = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
i = np.where((lon >= LON_MIN) & (lon <= LON_MAX))[0]
sub_diff = diff[np.ix_(j,i)]
sub_on = on[np.ix_(j,i)]
sub_off = off[np.ix_(j,i)]
sub_lat = lat[j]; sub_lon = lon[i]

flat = sub_diff.filled(np.nan) if hasattr(sub_diff,'filled') else sub_diff
flat_on = sub_on.filled(np.nan) if hasattr(sub_on,'filled') else sub_on
flat_off = sub_off.filled(np.nan) if hasattr(sub_off,'filled') else sub_off

print()
print("=== Top 10 cells where BIF DECREASED floodfrc ===")
order = np.argsort(np.where(np.isfinite(flat), flat, np.inf), axis=None)
shown = 0
for k in order:
    if shown >= 10: break
    if not np.isfinite(flat.flat[k]): continue
    jj, ii = np.unravel_index(k, sub_diff.shape)
    print(f"  lat={sub_lat[jj]:>+6.2f}  lon={sub_lon[ii]:>+7.2f}  on={flat_on[jj,ii]:>8.4f}  off={flat_off[jj,ii]:>8.4f}  diff={flat[jj,ii]:>+9.4f}")
    shown += 1
print()
print("=== Top 10 cells where BIF INCREASED floodfrc ===")
order = np.argsort(np.where(np.isfinite(flat), -flat, np.inf), axis=None)
shown = 0
for k in order:
    if shown >= 10: break
    if not np.isfinite(flat.flat[k]): continue
    jj, ii = np.unravel_index(k, sub_diff.shape)
    print(f"  lat={sub_lat[jj]:>+6.2f}  lon={sub_lon[ii]:>+7.2f}  on={flat_on[jj,ii]:>8.4f}  off={flat_off[jj,ii]:>8.4f}  diff={flat[jj,ii]:>+9.4f}")
    shown += 1
