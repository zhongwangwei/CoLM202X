"""CaMa main-NC coupling: DEF_CaMa_StrictDomain and DEF_CaMa_FloodFeedback.

StrictDomain = .false. (default) lets a regional domain run without complete
drainage basins: runoff over grid cells that no routing cell receives is dropped
and counted instead of stopping the run.  FloodFeedback = .true. (default) keeps the
two-way flood evaporation/infiltration exchange; .false. makes the coupling one-way.
"""
import pathlib
import re
import shutil
import subprocess
import tempfile

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
BUDGET = ROOT / "extends/CaMa/src/cmf_coupling_budget_mod.F90"


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


# ------------------------------------------------------------------ runtime probe

DRIVER = '''program check
use cmf_coupling_budget_mod
implicit none
integer :: x(2,1), y(2,1)
real(8) :: w(2,1), area(3,1), gridflow(3,1), flow(2)
x(:,1)=[1,2]; y=1; w=1d0; area(:,1)=[10d0,10d0,10d0]
call budget_init(x,y,w,area)
gridflow(:,1)=[1d0,2d0,4d0]        ! grid cell 3 has no routing cell
if (%s) call budget_set_strict(.false.)
call budget_runoff(gridflow,flow)
if (abs(sum(flow)-3d0)>1d-12) stop 2                 ! routed part only
if (abs(budget_unrouted()-4d0)>1d-12) stop 3         ! dropped part reported
gridflow(3,1)=0d0
call budget_runoff(gridflow,flow)
if (abs(budget_unrouted())>1d-12) stop 4             ! nothing dropped when every cell has a recipient
print '(A)', 'BUDGET_OK'
end program
'''


def _build(tmp, relaxed):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran not available")
    driver = pathlib.Path(tmp) / f"check_{int(relaxed)}.f90"
    driver.write_text(DRIVER % (".true." if relaxed else ".false."))
    exe = pathlib.Path(tmp) / f"check_{int(relaxed)}"
    subprocess.run([compiler, "-fcheck=all", str(BUDGET), str(driver), "-o", str(exe)],
                   cwd=tmp, check=True, capture_output=True)
    return exe


def test_runoff_without_a_recipient_is_dropped_and_reported_unless_strict():
    with tempfile.TemporaryDirectory() as tmp:
        ran = subprocess.run([str(_build(tmp, relaxed=True))], capture_output=True, text=True)
        assert ran.returncode == 0 and "BUDGET_OK" in ran.stdout, ran.stdout + ran.stderr


def test_the_budget_module_itself_is_strict_until_told_otherwise():
    with tempfile.TemporaryDirectory() as tmp:
        ran = subprocess.run([str(_build(tmp, relaxed=False))], capture_output=True, text=True)
        assert ran.returncode != 0
        assert "runoff has no routing recipient" in ran.stdout + ran.stderr


# ----------------------------------------------------------------------- wiring

def test_switches_are_declared_with_the_documented_defaults_and_reach_every_rank():
    src = read("share/MOD_Namelist.F90")
    assert re.search(r"logical\s+::\s+DEF_CaMa_FloodFeedback\s*=\s*\.true\.", src)
    assert re.search(r"logical\s+::\s+DEF_CaMa_StrictDomain\s*=\s*\.false\.", src)
    for name in ("DEF_CaMa_FloodFeedback", "DEF_CaMa_StrictDomain"):
        assert re.search(rf"^\s*{name},\s*&", src, re.M), name
        assert re.search(rf"mpi_bcast \({name}\s+,1\s+,mpi_logical", src), name


def test_domain_checks_stop_only_in_strict_mode():
    src = read("extends/CaMa/src/MOD_CaMa_colmCaMa.F90")
    assert "CALL budget_set_strict(DEF_CaMa_StrictDomain)" in src
    for message in ("regional domain cuts a river link", "regional domain cuts a bidirectional bifurcation"):
        for match in re.finditer(re.escape(message), src):
            before = src[max(0, match.start() - 160):match.start()]
            assert "IF(DEF_CaMa_StrictDomain)" in before, message
    assert "runoff_unrouted" in src.split("SUBROUTINE colm_cama_exit", 1)[1]
    budget = read("extends/CaMa/src/cmf_coupling_budget_mod.F90")
    assert "IF(strict_recipients .AND. ANY(" in budget
    assert re.search(r"LOGICAL :: strict_recipients = \.TRUE\.", budget)


def test_feedback_switch_reaches_cama_before_it_reads_its_input():
    colm = read("extends/CaMa/src/MOD_CaMa_colmCaMa.F90")
    assert colm.index("LCOLMFEEDBACK = DEF_CaMa_FloodFeedback") < \
        colm.index("CALL CMF_DRV_INPUT(DEF_UnitCatchment_file, DEF_CaMa_Restart_file)")
    nml = read("extends/CaMa/src/cmf_ctrl_nmlist_mod.F90")
    assert "LWEVAP=LCOLMFEEDBACK" in nml and "LWINFILT=LCOLMFEEDBACK" in nml
    assert re.search(r"LCOLMFEEDBACK\s*=\s*\.TRUE\.", read("extends/CaMa/src/yos_cmf_input.F90"))
