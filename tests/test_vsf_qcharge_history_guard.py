from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_unsupported_vsf_qcharge_is_neither_accumulated_nor_output():
    acc = (ROOT / "main/MOD_Vars_1DAccFluxes.F90").read_text()
    hist = (ROOT / "main/MOD_Hist.F90").read_text()
    assert "IF (.not. DEF_USE_VariablySaturatedFlow) CALL acc1d (qcharge, a_qcharge)" in acc
    assert "DEF_hist_vars%qcharge &\n            .and. (.not.DEF_USE_VariablySaturatedFlow)" in hist
