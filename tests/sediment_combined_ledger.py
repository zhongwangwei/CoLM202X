"""Independent per-class inventory ledger for two gridriver restart files.

The optional flux log is produced by the *temporary test build*'s SED_LEDGER
instrumentation. Only its input and outlet columns are used; the inventory
and residual columns printed by that build are deliberately ignored.
"""

import argparse
import json
import math
import re
from pathlib import Path

from netCDF4 import Dataset


FLUX = re.compile(
    r"SED_LEDGER\s+(\d+)\s+(\d+)\s+"
    r"[\d.Ee+-]+\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+"
    r"[\d.Ee+-]+\s+[\d.Ee+-]+"
)


def inventory(path: Path) -> dict[int, dict[str, float]]:
    with Dataset(path) as ds:
        nclass = int(ds["sed_n_meta"][0])
        nlayer = int(ds["sed_totlyrnum_meta"][0])
        porosity = ds["sed_lambda_meta"][:]
        if not (1 <= nclass <= 100 and 1 <= nlayer <= 100):
            raise ValueError("invalid sediment metadata")
        out = {}
        for cls in range(1, nclass + 1):
            def sum_field(name: str) -> float:
                return math.fsum(ds[name][:].filled(0).tolist())

            visible_suspended = sum_field(f"sedsto_{cls}")
            protected_suspended = sum_field(f"sedsto_protected_{cls}")
            protected_bed = sum_field(f"sedbed_protected_{cls}")
            # The ordinary layer and deposits store BULK bed volume; the
            # protected bed stores solid volume already.
            ordinary_bed = math.fsum(
                ((1.0 - porosity) * ds[f"layer_{cls}"][:]).filled(0).tolist()
            )
            ordinary_bed += math.fsum(
                math.fsum(
                    ((1.0 - porosity) * ds[f"seddep_{cls}_{layer}"][:]).filled(0).tolist()
                )
                for layer in range(1, nlayer + 1)
            )
            parts = dict(
                visible_suspended=visible_suspended,
                protected_suspended=protected_suspended,
                ordinary_bed_solid=ordinary_bed,
                protected_bed_solid=protected_bed,
            )
            parts["total"] = math.fsum(parts.values())
            out[cls] = parts
        return out


def read_fluxes(path: Path) -> dict[int, dict[str, float]]:
    entries: dict[int, dict[str, list[float]]] = {}
    for line in path.read_text().splitlines():
        match = FLUX.search(line)
        if match:
            _, cls, inp, outlet = match.groups()
            entry = entries.setdefault(int(cls), {"land_input": [], "outlet": []})
            entry["land_input"].append(float(inp))
            entry["outlet"].append(float(outlet))
    if not entries:
        raise ValueError("no SED_LEDGER flux rows")
    return {
        cls: {key: math.fsum(values) for key, values in entry.items()}
        for cls, entry in entries.items()
    }


def audit(start: Path, end: Path, flux_log: Path, tolerance: float) -> dict:
    before, after, fluxes = inventory(start), inventory(end), read_fluxes(flux_log)
    if before.keys() != after.keys() or before.keys() != fluxes.keys():
        raise ValueError("sediment class mismatch between restart and flux data")
    results = {}
    for cls in before:
        inp, outlet = fluxes[cls]["land_input"], fluxes[cls]["outlet"]
        delta = after[cls]["total"] - before[cls]["total"]
        residual = delta - inp + outlet
        results[str(cls)] = dict(
            start=before[cls], end=after[cls],
            land_input=inp, outlet=outlet,
            inventory_delta=delta, residual=residual,
            internal_conversion=dict(
                visible_suspended=after[cls]["visible_suspended"] - before[cls]["visible_suspended"],
                protected_suspended=after[cls]["protected_suspended"] - before[cls]["protected_suspended"],
                ordinary_bed_solid=after[cls]["ordinary_bed_solid"] - before[cls]["ordinary_bed_solid"],
                protected_bed_solid=after[cls]["protected_bed_solid"] - before[cls]["protected_bed_solid"],
            ),
        )
        if not math.isfinite(residual) or abs(residual) > tolerance:
            raise AssertionError(f"class {cls} budget residual {residual} exceeds {tolerance}")
    return {
        "formula": "M_end-M_start-land_input+outlet; M=visible_suspended+protected_suspended+(1-lambda)*(layer+seddep)+protected_bed_solid",
        "internal_flux_convention": "Dry/cap deposition and settling transfer suspended solid to the corresponding bed inventory; they are included in component deltas and cancel in M, not external sinks.",
        "classes": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", type=Path)
    parser.add_argument("end", type=Path)
    parser.add_argument("flux_log", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-3)
    args = parser.parse_args()
    print(json.dumps(audit(args.start, args.end, args.flux_log, args.tolerance), indent=2))
