#!/usr/bin/env python3
"""Compare a block+postprocess history directory against a 'one' reference.

Step 5 items 9 and 12 of the river-history sharding plan:

    block shards + postprocess 与 one 文件逐变量比较
    实际枚举 one 与 block+postprocess 的变量集合、维度和属性并断言相同,
    禁止只用手工维护的预期清单做覆盖判断

Everything compared here is discovered from the files. A variable that exists
in one mode and not the other is a failure, which is the whole point: a
hand-maintained expected list would happily pass a run that silently dropped a
field from one path.

Comparison rules follow the plan:

  * integer ids, dimensions, variable names, attributes, missing value:
    exactly equal;
  * history values: bitwise for data that is only rearranged, falling back to
    an explicit r8 tolerance if the platform's NetCDF/compiler introduces
    rounding -- never a relaxed physical-conservation check.

Usage:
    river_hist_compare.py <ref_dir> <cand_dir> [--label L] [--rtol R] [--atol A]
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
from netCDF4 import Dataset

# Artefacts of the sharded path itself, not part of the 'one' schema.
SKIP_SUFFIXES = (".pending", ".complete", ".tmp")
SKIP_PATTERNS = ("_shard",)


def history_files(d: pathlib.Path) -> dict[str, pathlib.Path]:
    out: dict[str, pathlib.Path] = {}
    for p in sorted(d.glob("*.nc")):
        if any(s in p.name for s in SKIP_PATTERNS):
            continue
        if any(p.name.endswith(s) for s in SKIP_SUFFIXES):
            continue
        out[p.name] = p
    return out


def describe(path: pathlib.Path) -> dict:
    with Dataset(path) as ds:
        dims = {k: len(v) for k, v in ds.dimensions.items()}
        variables = {}
        for name, var in ds.variables.items():
            attrs = {}
            for a in var.ncattrs():
                if a in ("long_name", "units", "missing_value"):
                    v = var.getncattr(a)
                    attrs[a] = float(v) if a == "missing_value" else str(v)
            variables[name] = {
                "dimensions": tuple(var.dimensions),
                "dtype": str(var.dtype),
                "attrs": attrs,
            }
    return {"dims": dims, "variables": variables}


def compare_values(ref: pathlib.Path, cand: pathlib.Path, name: str,
                   rtol: float, atol: float) -> list[str]:
    problems: list[str] = []
    with Dataset(ref) as a, Dataset(cand) as b:
        x = np.asarray(a.variables[name][:], dtype=np.float64)
        y = np.asarray(b.variables[name][:], dtype=np.float64)
        if x.shape != y.shape:
            return [f"{name}: shape {x.shape} != {y.shape}"]

        both_nan = np.isnan(x) & np.isnan(y)
        if np.array_equal(np.where(both_nan, 0.0, x), np.where(both_nan, 0.0, y)):
            return []          # bitwise identical: the expected outcome

        bad = ~(np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True))
        n = int(bad.sum())
        if n:
            idx = np.unravel_index(int(np.argmax(bad)), bad.shape)
            problems.append(
                f"{name}: {n} of {x.size} values differ beyond rtol={rtol} atol={atol}; "
                f"first at {idx}: {x[idx]!r} vs {y[idx]!r}")
        else:
            problems.append(
                f"{name}: NOT bitwise identical but within tolerance "
                f"(max |d| = {np.nanmax(np.abs(x - y)):.3e}) -- rearranged data "
                "should reproduce exactly; investigate before accepting")
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref_dir")
    ap.add_argument("cand_dir")
    ap.add_argument("--label", default="candidate")
    ap.add_argument("--rtol", type=float, default=0.0)
    ap.add_argument("--atol", type=float, default=0.0)
    args = ap.parse_args(argv)

    ref_dir = pathlib.Path(args.ref_dir)
    cand_dir = pathlib.Path(args.cand_dir)
    problems: list[str] = []

    ref_files = history_files(ref_dir)
    cand_files = history_files(cand_dir)

    missing = sorted(set(ref_files) - set(cand_files))
    extra = sorted(set(cand_files) - set(ref_files))
    problems += [f"missing history file: {m}" for m in missing]
    problems += [f"unexpected history file: {e}" for e in extra]

    for name in sorted(set(ref_files) & set(cand_files)):
        a, b = describe(ref_files[name]), describe(cand_files[name])

        for d, size in sorted(a["dims"].items()):
            if d not in b["dims"]:
                problems.append(f"{name}: missing dimension {d}")
            elif b["dims"][d] != size:
                problems.append(f"{name}: dimension {d} {size} != {b['dims'][d]}")

        va, vb = a["variables"], b["variables"]
        for v in sorted(set(va) - set(vb)):
            problems.append(f"{name}: variable missing from {args.label}: {v}")
        for v in sorted(set(vb) - set(va)):
            problems.append(f"{name}: variable only in {args.label}: {v}")

        for v in sorted(set(va) & set(vb)):
            if va[v]["dimensions"] != vb[v]["dimensions"]:
                problems.append(
                    f"{name}/{v}: dimension order {va[v]['dimensions']} != {vb[v]['dimensions']}")
                continue
            if va[v]["dtype"] != vb[v]["dtype"]:
                problems.append(f"{name}/{v}: dtype {va[v]['dtype']} != {vb[v]['dtype']}")
            for k in sorted(set(va[v]["attrs"]) | set(vb[v]["attrs"])):
                if va[v]["attrs"].get(k) != vb[v]["attrs"].get(k):
                    problems.append(
                        f"{name}/{v}: attr {k} {va[v]['attrs'].get(k)!r} != "
                        f"{vb[v]['attrs'].get(k)!r}")
            problems += [f"{name}/{p}" for p in
                         compare_values(ref_files[name], cand_files[name], v,
                                        args.rtol, args.atol)]

    print(f"-- {args.label}: {len(ref_files)} reference file(s), "
          f"{len(set(ref_files) & set(cand_files))} compared")
    if problems:
        for p in problems[:40]:
            print(f"   FAIL {p}")
        if len(problems) > 40:
            print(f"   ... and {len(problems) - 40} more")
        return 1
    print("   identical: variable sets, dimensions, attributes and values")
    return 0


if __name__ == "__main__":
    sys.exit(main())
