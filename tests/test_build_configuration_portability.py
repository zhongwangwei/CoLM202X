from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_gnu_makeoptions_discovers_libraries_and_gates_x86_flags() -> None:
    source = (ROOT / "include/Makeoptions.github").read_text(encoding="utf-8")

    assert "/usr/lib/x86_64-linux-gnu" not in source
    assert "-xcheck=stkovf" not in source
    assert "nf-config" in source
    assert "nc-config" in source
    assert "pkg-config" in source
    assert "--exists netcdf-fortran netcdf" in source
    assert "ifeq ($(shell uname -m),x86_64)" in source
    assert "-mcmodel=medium" in source
    assert "GFORTRAN_MAJOR" in source
    gate = source.index("GFORTRAN_MAJOR")
    assert gate < source.index("FOPTS += -fallow-argument-mismatch")


def test_ci_compiles_supported_ch4_pft_and_pc_configurations() -> None:
    cases = (ROOT / ".github/workflows/TestCaseLists").read_text(encoding="utf-8")

    assert "GRID_PFT_TRACER_CH4" in cases
    assert "GRID_PC_TRACER_CH4" in cases
    assert "SITE_PFT_TRACER_CH4" not in cases
