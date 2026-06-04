#!/usr/bin/env python3
"""Build CoLM TRACER chloride precipitation-concentration forcing files.

The TRACER runtime precip forcing expects a value with units tracer-mass per
water-mass.  For chloride wet deposition this script converts

    Cl wet deposition amount / precipitation water amount -> kg Cl kg-1 H2O

and writes files named like MOD_Tracer_Forcing expects for tracer forcing:

    <fprefix>_<YYYY>_<MM>.nc   (default groupby=month)
    <fprefix>_<YYYY>.nc        (groupby=year)
    <fprefix>_<YYYY>_<MM>_<DD>.nc (groupby=day)

Typical use with monthly wet-deposition totals and monthly precipitation totals:

    python make_chloride_precip_conc.py \
      --wetdep nadp_cldep_monthly.nc --wetdep-var cl_wetdep \
      --wetdep-unit mg_m-2 \
      --precip forcing_prcp_monthly.nc --precip-var precip \
      --precip-unit mm \
      --out-dir /path/to/forcing/tracer --fprefix CL_CONC --output-dtime-seconds 21600

The output variable defaults to cl_precip_conc and should be referenced from
&nl_colm_tracer_forcing with forcing_input_mode='value'.
"""

from __future__ import annotations

import argparse
import calendar
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import xarray as xr


WETDEP_UNITS = {
    "kg_m-2": (1.0, False),
    "kg/m2": (1.0, False),
    "kg_m-2_s-1": (1.0, True),
    "kg/m2/s": (1.0, True),
    "g_m-2": (1.0e-3, False),
    "g/m2": (1.0e-3, False),
    "g_m-2_s-1": (1.0e-3, True),
    "g/m2/s": (1.0e-3, True),
    "mg_m-2": (1.0e-6, False),
    "mg/m2": (1.0e-6, False),
    "mg_m-2_s-1": (1.0e-6, True),
    "mg/m2/s": (1.0e-6, True),
}

PRECIP_UNITS = {
    # 1 mm water over 1 m2 is approximately 1 kg water.
    "kg_m-2": (1.0, False),
    "kg/m2": (1.0, False),
    "kg_m-2_s-1": (1.0, True),
    "kg/m2/s": (1.0, True),
    "mm": (1.0, False),
    "mm_s-1": (1.0, True),
    "mm/s": (1.0, True),
    "mm_day-1": (1.0 / 86400.0, True),
    "mm/day": (1.0 / 86400.0, True),
    "m": (1000.0, False),
    "m_s-1": (1000.0, True),
    "m/s": (1000.0, True),
}


def _norm_unit(unit: str) -> str:
    return unit.strip().replace(" ", "_").replace("−", "-")


def _time_name(da: xr.DataArray) -> str:
    for name in da.dims:
        if name.lower() == "time" or "time" in name.lower():
            return name
    raise ValueError(f"Cannot find a time dimension in {da.name!r}: dims={da.dims}")


def _coord_name(da: xr.DataArray, candidates: Iterable[str]) -> str | None:
    lower = {name.lower(): name for name in list(da.coords) + list(da.dims)}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _period_seconds(time: xr.DataArray, fallback_days: float | None) -> xr.DataArray:
    values = time.values
    if values.size == 0:
        raise ValueError("empty time coordinate")
    if values.size == 1:
        if fallback_days is not None:
            seconds = np.array([fallback_days * 86400.0], dtype="float64")
        else:
            # Monthly wet-deposition products are often processed one month at
            # a time.  Use that calendar month as the record length when no
            # neighboring time stamp is available.
            start = np.datetime64(values[0], "ns")
            stamp = np.datetime_as_string(start, unit="D")
            year = int(stamp[0:4])
            month = int(stamp[5:7])
            seconds = np.array([calendar.monthrange(year, month)[1] * 86400.0], dtype="float64")
    else:
        t = values.astype("datetime64[ns]")
        dt = np.diff(t).astype("timedelta64[ns]").astype("float64") / 1.0e9
        last = dt[-1]
        seconds = np.concatenate([dt, [last]])
    return xr.DataArray(seconds, coords={time.dims[0]: time}, dims=(time.dims[0],))


def _amount_from_unit(da: xr.DataArray, unit: str, table: dict[str, tuple[float, bool]],
                      seconds: xr.DataArray) -> xr.DataArray:
    key = _norm_unit(unit)
    if key not in table:
        allowed = ", ".join(sorted(table))
        raise ValueError(f"Unsupported unit {unit!r}; allowed: {allowed}")
    factor, is_rate = table[key]
    out = da.astype("float64") * factor
    if is_rate:
        out = out * seconds
    return out


def _maybe_interp_to_target(ds: xr.Dataset, target: xr.Dataset, varname: str) -> xr.Dataset:
    da = ds[varname]
    lat = _coord_name(da, ["lat", "latitude"])
    lon = _coord_name(da, ["lon", "longitude"])
    target_lat = _coord_name(next(iter(target.data_vars.values())), ["lat", "latitude"])
    target_lon = _coord_name(next(iter(target.data_vars.values())), ["lon", "longitude"])
    if not lat or not lon or not target_lat or not target_lon:
        raise ValueError("Cannot identify lat/lon coordinates for --target-grid interpolation")
    return ds.interp({lat: target[target_lat], lon: target[target_lon]}, method="linear")


def _canonical_time_coord(ds: xr.Dataset, varname: str) -> xr.Dataset:
    da = ds[varname]
    tname = _time_name(da)
    if tname != "time":
        ds = ds.rename({tname: "time"})
    return ds



def _expand_records(ds: xr.Dataset, dtime_seconds: int) -> xr.Dataset:
    """Repeat each source-record value on a fixed CoLM forcing time step."""
    if dtime_seconds <= 0:
        raise ValueError("--output-dtime-seconds must be positive")
    times = ds.time.values.astype("datetime64[ns]")
    if times.size == 0:
        raise ValueError("cannot expand an empty time coordinate")
    if times.size == 1:
        # Prefer one calendar month for monthly products when no next stamp is
        # available. This is intentionally simple; users with single-record
        # nonmonthly inputs should provide an already expanded file instead.
        start = np.datetime64(times[0], "ns")
        stamp = np.datetime_as_string(start, unit="D")
        year = int(stamp[0:4])
        month = int(stamp[5:7])
        ndays = calendar.monthrange(year, month)[1]
        deltas = np.array([ndays * 86400.0], dtype="float64")
    else:
        deltas = np.diff(times).astype("timedelta64[ns]").astype("float64") / 1.0e9
        deltas = np.concatenate([deltas, [deltas[-1]]])

    pieces = []
    for i, (start, seconds) in enumerate(zip(times, deltas)):
        nstep = max(1, int(np.ceil(seconds / float(dtime_seconds))))
        offsets = (np.arange(nstep, dtype="int64") * int(dtime_seconds)).astype("timedelta64[s]")
        new_time = start + offsets.astype("timedelta64[ns]")
        piece = ds.isel(time=i).expand_dims(time=new_time)
        pieces.append(piece)
    return xr.concat(pieces, dim="time")

def _write_grouped(ds: xr.Dataset, out_dir: Path, fprefix: str, groupby: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    groupby = groupby.lower()
    if groupby == "month":
        keys = sorted({(int(t.dt.year), int(t.dt.month)) for t in ds.time})
        for year, month in keys:
            sub = ds.where((ds.time.dt.year == year) & (ds.time.dt.month == month), drop=True)
            path = out_dir / f"{fprefix}_{year:04d}_{month:02d}.nc"
            sub.to_netcdf(path)
            written.append(path)
    elif groupby == "year":
        keys = sorted({int(t.dt.year) for t in ds.time})
        for year in keys:
            sub = ds.where(ds.time.dt.year == year, drop=True)
            path = out_dir / f"{fprefix}_{year:04d}.nc"
            sub.to_netcdf(path)
            written.append(path)
    elif groupby == "day":
        keys = sorted({(int(t.dt.year), int(t.dt.month), int(t.dt.day)) for t in ds.time})
        for year, month, day in keys:
            sub = ds.where(
                (ds.time.dt.year == year) & (ds.time.dt.month == month) & (ds.time.dt.day == day),
                drop=True,
            )
            path = out_dir / f"{fprefix}_{year:04d}_{month:02d}_{day:02d}.nc"
            sub.to_netcdf(path)
            written.append(path)
    else:
        raise ValueError("--groupby must be month, year, or day")
    return written


def build(args: argparse.Namespace) -> list[Path]:
    wet = xr.open_dataset(args.wetdep)
    pr = xr.open_dataset(args.precip)
    if args.wetdep_var not in wet:
        raise KeyError(f"{args.wetdep_var!r} not in {args.wetdep}")
    if args.precip_var not in pr:
        raise KeyError(f"{args.precip_var!r} not in {args.precip}")

    wet = _canonical_time_coord(wet, args.wetdep_var)
    pr = _canonical_time_coord(pr, args.precip_var)

    if args.target_grid:
        target = xr.open_dataset(args.target_grid)
        wet = _maybe_interp_to_target(wet, target, args.wetdep_var)
        pr = _maybe_interp_to_target(pr, target, args.precip_var)

    wet_da, pr_da = xr.align(wet[args.wetdep_var], pr[args.precip_var], join="inner")
    if wet_da.sizes.get("time", 0) == 0:
        raise ValueError("wetdep and precip have no overlapping time records")

    seconds = _period_seconds(wet_da["time"], args.fallback_period_days)
    wet_amount = _amount_from_unit(wet_da, args.wetdep_unit, WETDEP_UNITS, seconds)
    precip_amount = _amount_from_unit(pr_da, args.precip_unit, PRECIP_UNITS, seconds)

    dry = precip_amount <= args.precip_min_kg_m2
    conc = xr.where(dry, 0.0, wet_amount / precip_amount)
    conc = conc.clip(min=0.0, max=args.max_conc_kg_kg)
    conc.name = args.output_var
    conc.attrs.update({
        "long_name": "chloride concentration in precipitation from wet deposition",
        "units": "kg kg-1",
        "description": "Computed as chloride wet-deposition amount divided by precipitation water amount for each record.",
        "wetdep_source": str(args.wetdep),
        "wetdep_variable": args.wetdep_var,
        "wetdep_unit_input": args.wetdep_unit,
        "precip_source": str(args.precip),
        "precip_variable": args.precip_var,
        "precip_unit_input": args.precip_unit,
        "precip_min_kg_m2": args.precip_min_kg_m2,
    })

    out = conc.to_dataset()
    if args.write_diagnostics:
        out["cl_wetdep_amount"] = wet_amount
        out["cl_wetdep_amount"].attrs.update({"units": "kg m-2", "long_name": "chloride wet deposition amount per record"})
        out["precip_amount"] = precip_amount
        out["precip_amount"].attrs.update({"units": "kg m-2", "long_name": "precipitation water amount per record"})
        out["cl_precip_dry_mask"] = dry.astype("i1")
        out["cl_precip_dry_mask"].attrs.update({"units": "1", "long_name": "1 where precipitation amount is at or below threshold"})

    # Keep common CoLM coordinate names/attrs if present; xarray preserves them.
    out.attrs.update({
        "title": "CoLM TRACER chloride precipitation-concentration forcing",
        "Conventions": "CF-1.8",
        "history": "created by preprocess/Forcings/ChlorideWetDep/make_chloride_precip_conc.py",
    })

    if args.output_dtime_seconds is not None:
        out = _expand_records(out, args.output_dtime_seconds)

    return _write_grouped(out, Path(args.out_dir), args.fprefix, args.groupby)


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wetdep", required=True, help="NetCDF file containing chloride wet deposition")
    p.add_argument("--wetdep-var", required=True, help="chloride wet-deposition variable name")
    p.add_argument("--wetdep-unit", required=True, help="input wetdep unit, e.g. mg_m-2 or kg_m-2_s-1")
    p.add_argument("--precip", required=True, help="NetCDF file containing precipitation")
    p.add_argument("--precip-var", required=True, help="precipitation variable name")
    p.add_argument("--precip-unit", required=True, help="input precipitation unit, e.g. mm or kg_m-2_s-1")
    p.add_argument("--out-dir", required=True, help="output directory under DEF_dir_forcing/tracer or similar")
    p.add_argument("--fprefix", default="CL_CONC", help="CoLM tracer forcing fprefix; default CL_CONC")
    p.add_argument("--output-var", default="cl_precip_conc", help="output variable name; default cl_precip_conc")
    p.add_argument("--groupby", choices=["month", "year", "day"], default="month", help="CoLM tracer forcing file grouping")
    p.add_argument(
        "--output-dtime-seconds",
        type=int,
        help="repeat each source value to this fixed CoLM forcing interval, e.g. 21600 for 6-hourly records",
    )
    p.add_argument("--target-grid", help="optional NetCDF whose rectilinear lat/lon grid is used for linear interpolation")
    p.add_argument("--fallback-period-days", type=float, help="record length for single-record rate inputs")
    p.add_argument("--precip-min-kg-m2", type=float, default=1.0e-9, help="dry threshold for precipitation amount per record")
    p.add_argument("--max-conc-kg-kg", type=float, default=1.0e-2, help="safety cap for concentration; default 1e-2 kg/kg")
    p.add_argument("--write-diagnostics", action="store_true", help="also write wetdep amount, precip amount, and dry mask")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    paths = build(args)
    print(f"wrote {len(paths)} file(s):")
    for path in paths:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
