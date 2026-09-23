from pathlib import Path
import re
import subprocess

from fortran_test_support import SMOKE_TIMEOUT, require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
NAMELIST = ROOT / "share" / "MOD_Namelist.F90"
AGGREGATION = ROOT / "mksrfdata" / "Aggregation_CanopyStructure.F90"
HTOP_READIN = ROOT / "mkinidata" / "MOD_HtopReadin.F90"
TIME_INVARIANTS = ROOT / "main" / "MOD_Vars_TimeInvariants.F90"
INTERCEPTION_SOURCES = (
    ROOT / "main" / "MOD_LeafInterception.F90",
    ROOT / "extends" / "interception" / "MOD_LeafInterception_Extended.F90",
)


def _flat(path: Path) -> str:
    return re.sub(r"\s+", " ", re.sub(r"&\s*", " ", path.read_text(encoding="utf-8").lower()))


def _extract_capacity(source: str) -> str:
    start = source.index("PURE REAL(r8) FUNCTION canopy_storage_capacity_colm2024")
    end_marker = "END FUNCTION canopy_storage_capacity_colm2024"
    end = source.index(end_marker, start) + len(end_marker)
    return source[start:end]


def test_usgs_colm2024_is_not_unconditionally_downgraded_in_namelist() -> None:
    source = NAMELIST.read_text(encoding="utf-8")
    assert "CoLM2024 interception capacity is not available for LULC_USGS" not in source
    assert "DEF_Interception_scheme is set to 1 (CoLM2014) for LULC_USGS" not in source


def test_preprocessed_usgs_namelist_preserves_scheme8_support(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    inc = tmp_path / "include"
    inc.mkdir()
    inc.joinpath("define.h").write_text("#define LULC_USGS\n", encoding="utf-8")
    result = subprocess.run(
        [compiler, "-cpp", "-E", "-P", "-I", str(inc), str(NAMELIST)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CoLM2024 interception capacity is not available for LULC_USGS" not in result.stdout
    assert "DEF_Interception_scheme = 8" in result.stdout


def test_usgs_canopy_structure_aggregation_is_compiled() -> None:
    source = AGGREGATION.read_text(encoding="utf-8")
    active = source[source.index("#if (defined LULC_USGS") :]
    assert "defined LULC_USGS" in active.splitlines()[0]
    for variable in (
        "NEEDLELEAF_CROWN_DEPTH",
        "NEEDLELEAF_CROWN_WIDTH",
        "BROADLEAF_CROWN_WIDTH",
        "ncd_patches",
        "ncw_patches",
        "bcw_patches",
    ):
        assert variable in active


def test_compiled_usgs_class_mapping_uses_only_supported_formulas(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    support = tmp_path / "support.F90"
    support.write_text(
        """
module MOD_Precision
  integer, parameter :: r8 = selected_real_kind(12)
end module MOD_Precision
""",
        encoding="utf-8",
    )

    driver = tmp_path / "driver.F90"
    driver.write_text(
        """
program driver
  use MOD_Precision
  use capacity_under_test, only: canopy_storage_capacity_colm2024
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  integer, parameter :: cls(12) = [8,11,12,13,14,15,6,9,10,18,21,22]
  real(r8), parameter :: expected(12) = [2._r8/3._r8,0.4_r8,0.75_r8,0.4_r8,0.75_r8,0.575_r8, &
                                         0.3_r8,0.3_r8,0.3_r8,0.3_r8,0.3_r8,0.3_r8]
  real(r8) :: got, nan
  integer :: i
  do i=1,size(cls)
    got = canopy_storage_capacity_colm2024(0.1_r8,2._r8,1._r8,2._r8,0._r8, &
                                            12._r8,5._r8,4._r8,4._r8,cls(i),.false.)
    if (abs(got-expected(i)) > 1.e-12_r8) error stop i
  enddo
  nan = ieee_value(0._r8, ieee_quiet_nan)
  got = canopy_storage_capacity_colm2024(0.1_r8,2._r8,1._r8,2._r8,0._r8, &
                                          12._r8,nan,4._r8,4._r8,13,.false.)
  if (abs(got-0.4_r8) > 1.e-12_r8) error stop 101
  got = canopy_storage_capacity_colm2024(0.1_r8,2._r8,1._r8,2._r8,0._r8, &
                                          12._r8,5._r8,4._r8,nan,12,.false.)
  if (abs(got-0.75_r8) > 1.e-12_r8) error stop 102
  print *, 'USGS_CAPACITY_OK'
end program driver
""",
        encoding="utf-8",
    )

    for index, source_path in enumerate(INTERCEPTION_SOURCES):
        module = tmp_path / f"capacity_{index}.F90"
        module.write_text(
            "module capacity_under_test\n"
            "  use MOD_Precision\n"
            "  implicit none\n"
            "contains\n"
            + _extract_capacity(source_path.read_text(encoding="utf-8"))
            + "\nend module capacity_under_test\n",
            encoding="utf-8",
        )
        exe = tmp_path / f"capacity_{index}"
        built = subprocess.run(
            [
                compiler,
                "-cpp",
                "-DLULC_USGS",
                "-fcheck=all",
                "-ffpe-trap=invalid,zero,overflow",
                "-ffree-line-length-none",
                str(support),
                str(module),
                str(driver),
                "-o",
                str(exe),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=SMOKE_TIMEOUT,
        )
        assert built.returncode == 0, built.stdout + built.stderr
        ran = subprocess.run(
            [str(exe)], cwd=tmp_path, capture_output=True, text=True, timeout=SMOKE_TIMEOUT
        )
        assert ran.returncode == 0, ran.stdout + ran.stderr
        assert "USGS_CAPACITY_OK" in ran.stdout


def test_interception_readme_documents_usgs_scope_and_global_fallback() -> None:
    readme = _flat(ROOT / "extends" / "interception" / "README.md")
    assert "shrubland (8)" in readme
    assert "broadleaf forest (11/13)" in readme
    assert "needleleaf forest (12/14)" in readme
    assert "mixed forest (15)" in readme
    assert "select scheme 1 for the whole run" in readme
    assert "must never mix invented per-patch structure" in readme
    assert "classification mean `htop0`" in readme
    assert "not an observed crown height" in readme


def test_usgs_initialization_and_restart_apply_the_same_global_structure_gate() -> None:
    for path in (HTOP_READIN, TIME_INVARIANTS):
        source = _flat(path)
        assert "case (11,13)" in source
        assert "case (12,14)" in source
        assert "case (15)" in source
        assert "case (8)" in source
        assert (
            "call mpi_allreduce(mpi_in_place, usgs_counts, 4, mpi_integer, "
            "mpi_sum, p_comm_glb, p_err)"
        ) in source
        assert "if (usgs_counts(2) < usgs_counts(1)) then" in source
        assert "downgrading the whole run to colm2014 interception" in source
        assert "def_interception_scheme = 1" in source
        # A shrub contributes to the supported count, so a shrub-only domain
        # skips the no-supported-vegetation downgrade without invented crown data.
        assert "usgs_counts(3) = usgs_counts(3) + 1" in source
        assert "elseif (usgs_counts(1)+usgs_counts(3) == 0) then" in source
