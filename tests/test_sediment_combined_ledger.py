"""The offline ledger counts transfers among all four class inventories."""

import json
from pathlib import Path
import subprocess
import sys

from netCDF4 import Dataset


SCRIPT = Path(__file__).with_name("sediment_combined_ledger.py")


def _write_state(path: Path, visible: float, protected: float, bed: float) -> None:
    with Dataset(path, "w") as ds:
        ds.createDimension("ucatch", 1)
        for name, value in {
            "sed_n_meta": 1,
            "sed_totlyrnum_meta": 1,
            "sed_lambda_meta": 0.4,
            "sedsto_1": visible,
            "sedsto_protected_1": protected,
            "sedbed_protected_1": bed,
            "layer_1": 0,
            "seddep_1_1": 0,
        }.items():
            ds.createVariable(name, "f8", ("ucatch",))[:] = [value]


def test_offline_ledger_balances_visible_and_protected_stocks(tmp_path: Path) -> None:
    start, end, flux = (tmp_path / name for name in ("start.nc", "end.nc", "flux.log"))
    _write_state(start, 1.0, 0.0, 0.0)
    _write_state(end, 0.6, 0.2, 0.3)
    flux.write_text("SED_LEDGER 0 1 1.0 0.2 0.1 1.1 0.0\n")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(start), str(end), str(flux)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    cls = json.loads(result.stdout)["classes"]["1"]
    assert abs(cls["residual"]) < 1e-12
    assert cls["end"]["protected_suspended"] == 0.2
    assert cls["end"]["protected_bed_solid"] == 0.3
