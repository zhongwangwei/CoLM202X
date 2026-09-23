"""DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT: optional momentum limit on the routing substep.

The routing substep is the smallest of three per-cell limits: CFL, storage
(don't drain a cell) and momentum (don't let a cell's momentum pass through zero).
On the 15-min network the momentum limit sets the substep in most substeps.  CaMa-Flood
has no such limit.  The switch drops only the third limit and defaults to keeping it,
so existing runs do not change.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


def test_switch_defaults_to_the_existing_behaviour_and_reaches_every_rank():
    src = read("share/MOD_Namelist.F90")
    assert re.search(r"logical\s+::\s+DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT\s*=\s*\.true\.", src)
    assert re.search(r"^\s*DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT,\s*&", src, re.M)
    assert re.search(r"mpi_bcast \(DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT,1\s+,mpi_logical", src)


def test_only_the_momentum_limit_is_switched_and_its_formula_is_unchanged():
    src = read("main/HYDRO/MOD_Grid_RiverLakeFlow.F90")
    assert "DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT" in src.split("SUBROUTINE grid_riverlake_flow (", 1)[1].split("IMPLICIT NONE", 1)[0]
    block = src.split("! constraint 3: avoid change of flow direction", 1)[1].split("dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))", 1)[0]
    assert re.search(r"IF \(DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT\) THEN\s+IF \(\.not\. is_built_resv\(i\)\) THEN", block)
    assert "abs(momen_riv(i) * topo_rivare(i) / (sum_mflux_riv(i)-sum_zgrad_riv(i)))" in block
    assert "(abs(veloc_riv(i)) > 0.1_r8)" in block
    # CFL and storage limits stay unconditional
    before = src.split("! constraint 3: avoid change of flow direction", 1)[0]
    tail = before.split("! constraint 1: CFL condition", 1)[1]
    assert "DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT" not in tail
