"""Schema/value lock for the river-history reference output.

Step 1.2 of ``.omx/plans/river-history-sharded-output.md``:

    记录变量名、维度顺序、属性、missing value、时间轴和数值。
    验收：测试能识别变量缺失、维度转置和 ID 错位。

This module is deliberately independent of how the file was produced, so the
same checks apply to

* the current ``one``-mode reference written by
  ``tests/river_hist_baseline_harness.F90``, and
* the postprocess-aggregated file produced later from IO-group shards
  (plan step 4), which must be schema- and value-identical.

The ID-misplacement check does not merely diff two files: it recomputes the
expected grid from the unit-catchment id mapping, so a candidate that is
self-consistent but reassembled in rank order is still rejected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from netCDF4 import Dataset

SPVAL = -1.0e36
# Values in the reference encode value = gid + 1e-3*ivar + 1e-6*itime, so a
# tolerance well below 1e-6 still separates neighbouring ids unambiguously.
VALUE_ATOL = 1.0e-9


@dataclass
class VarSpec:
    name: str
    dimensions: tuple[str, ...]
    dtype: str
    attrs: dict[str, Any] = field(default_factory=dict)


@dataclass
class FileSpec:
    dimensions: dict[str, int]
    variables: dict[str, VarSpec]
    time: list[float]


def describe(path: str) -> FileSpec:
    """Capture the parts of the schema the plan requires to stay stable."""
    with Dataset(path) as ds:
        dims = {name: (0 if d.isunlimited() else len(d)) for name, d in ds.dimensions.items()}
        variables = {}
        for name, var in ds.variables.items():
            attrs = {}
            for key in ("long_name", "units", "missing_value"):
                if key in var.ncattrs():
                    value = var.getncattr(key)
                    attrs[key] = float(value) if key == "missing_value" else str(value)
            variables[name] = VarSpec(
                name=name,
                dimensions=tuple(var.dimensions),
                dtype=str(var.dtype),
                attrs=attrs,
            )
        time = ds.variables["time"][:].tolist() if "time" in ds.variables else []
    return FileSpec(dimensions=dims, variables=variables, time=time)


def compare_schema(ref: FileSpec, cand: FileSpec) -> list[str]:
    """Return every schema difference; empty list means identical."""
    problems: list[str] = []

    missing = sorted(set(ref.variables) - set(cand.variables))
    for name in missing:
        problems.append(f"missing variable: {name}")

    extra = sorted(set(cand.variables) - set(ref.variables))
    for name in extra:
        problems.append(f"unexpected variable: {name}")

    for name in sorted(set(ref.variables) & set(cand.variables)):
        rv, cv = ref.variables[name], cand.variables[name]
        if rv.dimensions != cv.dimensions:
            problems.append(
                f"dimension order differs for {name}: {rv.dimensions} != {cv.dimensions}"
            )
        if rv.dtype != cv.dtype:
            problems.append(f"dtype differs for {name}: {rv.dtype} != {cv.dtype}")
        for key in sorted(set(rv.attrs) | set(cv.attrs)):
            a, b = rv.attrs.get(key), cv.attrs.get(key)
            if isinstance(a, float) and isinstance(b, float):
                if not math.isclose(a, b, rel_tol=0.0, abs_tol=abs(a) * 1e-12):
                    problems.append(f"attr {key} differs for {name}: {a} != {b}")
            elif a != b:
                problems.append(f"attr {key} differs for {name}: {a!r} != {b!r}")

    for name, size in sorted(ref.dimensions.items()):
        if name not in cand.dimensions:
            problems.append(f"missing dimension: {name}")
        elif cand.dimensions[name] != size:
            problems.append(
                f"dimension {name} size differs: {size} != {cand.dimensions[name]}"
            )

    if len(ref.time) != len(cand.time):
        problems.append(f"time record count differs: {len(ref.time)} != {len(cand.time)}")
    else:
        for i, (a, b) in enumerate(zip(ref.time, cand.time)):
            if a != b:
                problems.append(f"time[{i}] differs: {a} != {b}")

    return problems


def expected_ucat_grid(
    total_numucat: int, nlon: int, nlat: int, ivar: int, itime: int
) -> np.ndarray:
    """Rebuild the grid the harness must have written, from the id mapping alone.

    Mirrors ``build_coordinates`` and ``fill_local_vector`` in
    ``tests/river_hist_baseline_harness.F90``.
    """
    grid = np.full((nlat, nlon), SPVAL, dtype=np.float64)
    for gid in range(1, total_numucat + 1):
        x = (gid - 1) % nlon
        y = ((gid - 1) // nlon) % nlat
        grid[y, x] = gid + 1.0e-3 * ivar + 1.0e-6 * itime
    return grid


def check_ucat_values(
    path: str, varname: str, total_numucat: int, nlon: int, nlat: int, ivar: int, itime: int
) -> list[str]:
    """Reject any unit-catchment id landing in the wrong grid cell."""
    expected = expected_ucat_grid(total_numucat, nlon, nlat, ivar, itime)
    with Dataset(path) as ds:
        if varname not in ds.variables:
            return [f"missing variable: {varname}"]
        actual = np.asarray(ds.variables[varname][itime - 1, :, :], dtype=np.float64)

    if actual.shape != expected.shape:
        return [f"shape differs for {varname}: {expected.shape} != {actual.shape}"]

    bad = ~np.isclose(actual, expected, rtol=0.0, atol=VALUE_ATOL, equal_nan=True)
    if not bad.any():
        return []

    count = int(bad.sum())
    ys, xs = np.nonzero(bad)
    first = f"(lat={ys[0]}, lon={xs[0]}) expected {expected[ys[0], xs[0]]!r} got {actual[ys[0], xs[0]]!r}"
    return [f"{varname}: {count} cell(s) hold the wrong unit-catchment value; first {first}"]


def check_bif_values(
    path: str,
    varname: str,
    npthlev: int,
    total_npthout: int,
    itime: int,
    level_dim: str = "bifurcation_level",
    pathway_dim: str = "bifurcation_pathway",
) -> list[str]:
    """Pathway columns must follow pth_global_id, never rank order.

    The axis order is resolved from the variable's dimension names rather than
    assumed: CoLM's ``ncio_write_serial_time`` takes dimension names in Fortran
    (fastest-first) order, so a Fortran ``(level, pathway)`` array lands on disk
    as ``(time, pathway, level)``. Reading by name keeps this check valid for
    both that file and any aggregated file that stores the other order.
    """
    with Dataset(path) as ds:
        if varname not in ds.variables:
            return [f"missing variable: {varname}"]
        var = ds.variables[varname]
        dims = list(var.dimensions)
        if level_dim not in dims or pathway_dim not in dims:
            return [f"{varname}: expected dimensions {level_dim!r} and {pathway_dim!r}, got {dims}"]
        actual = np.asarray(var[itime - 1, :, :], dtype=np.float64)
        # positions within the post-time slice
        sliced = [d for d in dims if d != "time"]
        if sliced.index(level_dim) > sliced.index(pathway_dim):
            actual = actual.T

    expected = np.empty((npthlev, total_npthout), dtype=np.float64)
    for gid in range(1, total_npthout + 1):
        for lev in range(1, npthlev + 1):
            expected[lev - 1, gid - 1] = gid + 1.0e-3 * lev + 1.0e-6 * itime

    if actual.shape != expected.shape:
        return [f"shape differs for {varname}: {expected.shape} != {actual.shape}"]

    bad = ~np.isclose(actual, expected, rtol=0.0, atol=VALUE_ATOL, equal_nan=True)
    if not bad.any():
        return []
    return [f"{varname}: {int(bad.sum())} pathway cell(s) reassembled out of global-id order"]
