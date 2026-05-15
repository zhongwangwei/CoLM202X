#!/usr/bin/env python3
"""Reset CoLM land tracer restart pools to fixed initial isotope ratios.

This intentionally rewrites only tracer variables. The corresponding water
state variables are left unchanged.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


RESET_MAP = {
    "trc_ldew_rain": ("ldew_rain", "nonnegative"),
    "trc_ldew_snow": ("ldew_snow", "nonnegative"),
    "trc_wliq_soisno": ("wliq_soisno", "nonnegative"),
    "trc_wice_soisno": ("wice_soisno", "nonnegative"),
    "trc_wa": ("wa", "signed"),
    "trc_wdsrf": ("wdsrf", "nonnegative"),
    "trc_wetwat": ("wetwat", "nonnegative"),
    "trc_scv": ("scv", "nonnegative"),
    "trc_waterstorage": ("waterstorage", "nonnegative"),
}

DELTA_RESET_NAMES = (
    "trc_leaf_delta_e",
    "trc_leaf_delta_b",
)

ZERO_RESET_NAMES = (
    "trc_leaf_water_moles",
    "trc_leaf_iso_storage",
)

ONE_RESET_NAMES = (
    "trc_leaf_peclet",
)


def parse_csv_float(value: str) -> np.ndarray:
    return np.array([float(x.strip()) for x in value.split(",") if x.strip()], dtype=np.float64)


def iter_nc_files(path: Path):
    if path.is_file():
        yield path
        return
    yield from sorted(path.rglob("*.nc"))


def tracer_axis(var) -> int | None:
    try:
        return var.dimensions.index("tracer")
    except ValueError:
        return None


def reset_array(ds: Dataset, tracer_name: str, water_name: str, policy: str, ratios: np.ndarray):
    tracer = ds.variables[tracer_name]
    axis = tracer_axis(tracer)
    if axis is None:
        return None

    if len(ds.dimensions["tracer"]) != len(ratios):
        raise ValueError(
            f"{tracer_name}: file tracer dimension is {len(ds.dimensions['tracer'])}, "
            f"but {len(ratios)} ratios were provided"
        )

    water = np.asarray(ds.variables[water_name][...], dtype=np.float64)
    if policy == "nonnegative":
        water = np.maximum(water, 0.0)

    expected_water_shape = tuple(
        len(ds.dimensions[dim]) for dim in tracer.dimensions if dim != "tracer"
    )
    if water.shape != expected_water_shape:
        raise ValueError(
            f"{tracer_name}: water shape {water.shape} does not match tracer "
            f"non-tracer shape {expected_water_shape}"
        )

    ratio_shape = [1] * len(tracer.dimensions)
    ratio_shape[axis] = len(ratios)
    return np.expand_dims(water, axis) * ratios.reshape(ratio_shape)


def reset_constant_by_tracer(ds: Dataset, tracer_name: str, values: np.ndarray):
    tracer = ds.variables[tracer_name]
    axis = tracer_axis(tracer)
    if axis is None:
        return None

    if len(ds.dimensions["tracer"]) != len(values):
        raise ValueError(
            f"{tracer_name}: file tracer dimension is {len(ds.dimensions['tracer'])}, "
            f"but {len(values)} values were provided"
        )

    value_shape = [1] * len(tracer.dimensions)
    value_shape[axis] = len(values)
    return np.broadcast_to(values.reshape(value_shape), tracer.shape)


def reset_file(path: Path, ratios: np.ndarray, init_delta: np.ndarray, apply: bool) -> int:
    mode = "r+" if apply else "r"
    changed = 0
    with Dataset(path, mode) as ds:
        if "tracer" not in ds.dimensions:
            return 0

        for tracer_name, (water_name, policy) in RESET_MAP.items():
            if tracer_name not in ds.variables or water_name not in ds.variables:
                continue
            arr = reset_array(ds, tracer_name, water_name, policy, ratios)
            if arr is None:
                continue
            changed += 1
            if apply:
                ds.variables[tracer_name][...] = arr

        for tracer_name in DELTA_RESET_NAMES:
            if tracer_name not in ds.variables:
                continue
            arr = reset_constant_by_tracer(ds, tracer_name, init_delta)
            if arr is None:
                continue
            changed += 1
            if apply:
                ds.variables[tracer_name][...] = arr

        for tracer_name in ZERO_RESET_NAMES:
            if tracer_name not in ds.variables:
                continue
            changed += 1
            if apply:
                ds.variables[tracer_name][...] = 0.0

        for tracer_name in ONE_RESET_NAMES:
            if tracer_name not in ds.variables:
                continue
            changed += 1
            if apply:
                ds.variables[tracer_name][...] = 1.0
    return changed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path, help="Restart NetCDF file or directory to scan recursively")
    parser.add_argument("--ref-ratio", default="2.0052e-3,1.5576e-4")
    parser.add_argument("--init-delta", default="-10.0,-70.0")
    parser.add_argument("--apply", action="store_true", help="Rewrite files in place")
    parser.add_argument("--backup-dir", type=Path, help="Copy each changed file here before rewriting")
    args = parser.parse_args()

    ref_ratio = parse_csv_float(args.ref_ratio)
    init_delta = parse_csv_float(args.init_delta)
    if ref_ratio.shape != init_delta.shape:
        raise SystemExit("--ref-ratio and --init-delta must have the same length")
    ratios = ref_ratio * (1.0 + init_delta / 1000.0)

    files_seen = 0
    files_changed = 0
    vars_changed = 0

    for path in iter_nc_files(args.path):
        files_seen += 1
        try:
            nvars = reset_file(path, ratios, init_delta, apply=False)
            if nvars == 0:
                continue

            files_changed += 1
            vars_changed += nvars
            print(f"{'RESET' if args.apply else 'WOULD_RESET'} {path} vars={nvars}")

            if args.apply:
                if args.backup_dir is not None:
                    backup_path = args.backup_dir / path.name
                    backup_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(path, backup_path)
                reset_file(path, ratios, init_delta, apply=True)
        except Exception as exc:
            print(f"SKIP {path}: {exc}")

    action = "reset" if args.apply else "would reset"
    print(f"Scanned {files_seen} files; {action} {vars_changed} variables in {files_changed} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
