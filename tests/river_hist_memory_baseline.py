#!/usr/bin/env python3
"""Root-rank memory baseline for the river-history `one` path.

Step 1.5 of ``.omx/plans/river-history-sharded-output.md``.

Quantifies the global-master allocations that the sharded design must remove,
using the formulas the plan fixes:

    unitcat fields : 8 * (totalnumucat + nlon*nlat + max_local_ucat)   bytes
    BIF fields     : 8 * npthlev_bif * (totalnpthout + max_local_path) bytes

These are lower bounds on the *peak per-variable* root footprint: they cover
``wdata(totalvlen)`` and ``wdata2d(nlon,nlat)`` in
``vector_gather_map2grid_and_write`` plus the ``rcache`` staging buffer in
``vector_gather_to_master``, and exclude NetCDF's own internal copies and any
compression buffer. Measured RSS/HWM must therefore be recorded alongside.

Parameters come either from an explicit ``--param`` set or from a CoLM river
network file, so the number quoted in the go/no-go is the target case's, not a
guess.

Usage
-----
    # explicit
    tests/river_hist_memory_baseline.py --totalnumucat 2200000 \
        --nlon 4320 --nlat 2160 --totalnpthout 180000 --npthlev 3 \
        --workers 512 --node-memory-gb 256

    # from a river network / history file
    tests/river_hist_memory_baseline.py --from-netcdf /path/hist_unitcat.nc \
        --workers 512 --node-memory-gb 256
"""

from __future__ import annotations

import argparse
import json
import sys

BYTES_PER_REAL8 = 8


def unitcat_root_bytes(totalnumucat: int, nlon: int, nlat: int, max_local_ucat: int) -> int:
    """Peak root bytes for one gathered+regridded unit-catchment field."""
    return BYTES_PER_REAL8 * (totalnumucat + nlon * nlat + max_local_ucat)


def bif_root_bytes(npthlev: int, totalnpthout: int, max_local_path: int) -> int:
    """Peak root bytes for the bifurcation pathway matrix."""
    return BYTES_PER_REAL8 * npthlev * (totalnpthout + max_local_path)


def even_split(total: int, workers: int, imbalance: float = 1.0) -> int:
    """Largest per-worker share.

    ``imbalance`` scales the even share to reflect that the river decomposition
    follows basins, not equal counts, so the busiest worker holds more than
    ``total/workers``. It only affects the ``max_local`` staging buffer, which
    is the smallest of the three terms; the dominant ``totalnumucat`` and
    ``nlon*nlat`` terms are independent of the decomposition. Pass the observed
    ratio (max_local * workers / total) from a real run when it is known.
    """
    if workers <= 0:
        return total
    return min(total, int(-(-total // workers) * max(imbalance, 1.0)))


def read_from_netcdf(path: str) -> dict[str, int]:
    from netCDF4 import Dataset

    with Dataset(path) as ds:
        dims = {k: len(v) for k, v in ds.dimensions.items() if not v.isunlimited()}
        out: dict[str, int] = {}
        for key, names in (
            ("nlon", ("lon_ucat", "lon", "longitude")),
            ("nlat", ("lat_ucat", "lat", "latitude")),
            ("totalnpthout", ("bifurcation_pathway", "npthout")),
            ("npthlev", ("bifurcation_level", "npthlev")),
            ("totalnumucat", ("unitcat", "numucat", "ucat")),
        ):
            for name in names:
                if name in dims:
                    out[key] = dims[name]
                    break
    if "totalnumucat" not in out and {"nlon", "nlat"} <= out.keys():
        # Regridded history carries no unitcat dimension; the grid is the
        # upper bound on how many unit catchments can be present.
        out["totalnumucat"] = out["nlon"] * out["nlat"]
        out["totalnumucat_inferred"] = 1
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--from-netcdf", help="river history/network file to read dimensions from")
    p.add_argument("--totalnumucat", type=int)
    p.add_argument("--nlon", type=int)
    p.add_argument("--nlat", type=int)
    p.add_argument("--totalnpthout", type=int, default=0)
    p.add_argument("--npthlev", type=int, default=3)
    p.add_argument("--workers", type=int, default=0,
                   help="worker rank count; used to derive max local lengths")
    p.add_argument("--max-local-ucat", type=int,
                   help="override the derived largest per-worker unitcat count")
    p.add_argument("--max-local-path", type=int)
    p.add_argument("--imbalance", type=float, default=1.0,
                   help="basin decomposition imbalance factor for max-local "
                        "staging buffers (1.0 = perfectly even)")
    p.add_argument("--node-memory-gb", type=float, default=0.0,
                   help="per-node memory budget, for the go/no-go 25%% threshold")
    p.add_argument("--json", action="store_true", help="emit machine-readable output")
    args = p.parse_args(argv)

    params: dict[str, int] = {}
    if args.from_netcdf:
        params.update(read_from_netcdf(args.from_netcdf))
    for key in ("totalnumucat", "nlon", "nlat", "totalnpthout", "npthlev"):
        value = getattr(args, key)
        if value:
            params[key] = value

    missing = [k for k in ("totalnumucat", "nlon", "nlat") if not params.get(k)]
    if missing:
        p.error(f"missing required parameter(s): {', '.join(missing)} "
                "(pass explicitly or via --from-netcdf)")

    inferred = params.pop("totalnumucat_inferred", 0)

    max_local_ucat = args.max_local_ucat or even_split(
        params["totalnumucat"], args.workers, args.imbalance)
    max_local_path = args.max_local_path or even_split(
        params.get("totalnpthout", 0), args.workers, args.imbalance)

    ucat = unitcat_root_bytes(params["totalnumucat"], params["nlon"], params["nlat"],
                              max_local_ucat)
    bif = bif_root_bytes(params.get("npthlev", 3), params.get("totalnpthout", 0),
                         max_local_path)
    peak = max(ucat, bif)

    result = {
        **params,
        "workers": args.workers,
        "max_local_ucat": max_local_ucat,
        "max_local_path": max_local_path,
        "unitcat_root_bytes": ucat,
        "bif_root_bytes": bif,
        "peak_root_bytes": peak,
        "totalnumucat_inferred_from_grid": bool(inferred),
    }
    if args.node_memory_gb > 0:
        budget = args.node_memory_gb * 1024**3
        share = peak / budget
        result["node_memory_bytes"] = int(budget)
        result["peak_fraction_of_node"] = share
        result["exceeds_go_nogo_25pct"] = share >= 0.25

    if args.json:
        print(json.dumps(result, indent=2))
        return 0

    print("river-history root memory baseline (plan step 1.5)")
    print(f"  totalnumucat        {params['totalnumucat']:,}"
          + ("   [inferred from grid size]" if inferred else ""))
    print(f"  grid  nlon x nlat   {params['nlon']:,} x {params['nlat']:,}"
          f" = {params['nlon']*params['nlat']:,}")
    print(f"  totalnpthout        {params.get('totalnpthout', 0):,}"
          f"   npthlev {params.get('npthlev', 3)}")
    print(f"  workers             {args.workers or 'n/a'}")
    print(f"  max local ucat      {max_local_ucat:,}")
    print(f"  max local pathway   {max_local_path:,}"
          + (f"   [imbalance x{args.imbalance:g}]" if args.imbalance != 1.0 else ""))
    print()
    print(f"  unitcat field root peak   {ucat/1e6:12.1f} MB")
    print(f"  BIF matrix   root peak    {bif/1e6:12.1f} MB")
    print(f"  ---> peak root lower bound {peak/1e6:11.1f} MB"
          f"  ({peak/1024**3:.2f} GiB)")
    if args.node_memory_gb > 0:
        print(f"  node budget               {args.node_memory_gb:12.1f} GB"
              f"   -> {result['peak_fraction_of_node']*100:.1f}% of one node")
        verdict = "GO" if result["exceeds_go_nogo_25pct"] else "does not by itself justify"
        print(f"  go/no-go memory criterion (>=25%): {verdict}")
    print()
    print("  NOTE: lower bound only -- excludes NetCDF internal copies and the")
    print("        compression buffer. Record measured RSS/HWM next to this.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
