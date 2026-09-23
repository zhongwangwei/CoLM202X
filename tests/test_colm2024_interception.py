from pathlib import Path

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
NORMAL = (ROOT / "main" / "MOD_LeafInterception.F90").read_text()
EXTENDED = (
    ROOT / "extends" / "interception" / "MOD_LeafInterception_Extended.F90"
).read_text()
INVARIANTS = (ROOT / "main" / "MOD_Vars_TimeInvariants.F90").read_text()
READIN = (ROOT / "mkinidata" / "MOD_HtopReadin.F90").read_text()
NAMELIST = (ROOT / "share" / "MOD_Namelist.F90").read_text()


def compact(source: str) -> str:
    return " ".join(source.lower().replace("&", "").split())


def capacity(
    veg_class: int,
    is_pft: bool,
    wind: float,
    htop: float,
    ncd: float,
    ncw: float,
    bcw: float,
    legacy: float,
) -> float:
    if is_pft:
        canopy_type = 1 if 1 <= veg_class <= 3 else 2 if 4 <= veg_class <= 8 else 3 if 9 <= veg_class <= 11 else 0
    else:
        canopy_type = 1 if veg_class in (1, 3) else 2 if veg_class in (2, 4) else 4 if veg_class == 5 else 3 if veg_class in (6, 7) else 0

    needle_valid = 0 < ncd < 1000 and 0 < ncw < 1000
    broad_valid = 0 < bcw < 1000 and 0 < htop < 1000
    needle = (min(11, max(3, ncd)) + min(7, max(2.9, ncw))) / (
        4 * (1 + min(3.6, max(1, wind)))
    )
    ratio = min(7, max(1, htop / bcw)) if broad_valid else 1
    broad = min(8, max(2, bcw)) / (2 * (min(4, max(1.5, wind)) + ratio))

    if canopy_type == 1 and needle_valid:
        return needle
    if canopy_type == 2 and broad_valid:
        return broad
    if canopy_type == 3:
        return 0.5 * (1 + 1 / (1 + min(4, max(1, wind))))
    if canopy_type == 4 and needle_valid and broad_valid:
        return 0.5 * (needle + broad)
    return legacy


def test_colm2024_capacity_formulas_and_class_tables():
    assert capacity(1, False, 2, 12, 5, 4, 4, 0.3) == 0.75
    assert capacity(4, True, 2, 12, 5, 4, 4, 0.3) == 0.4
    assert capacity(5, False, 2, 12, 5, 4, 4, 0.3) == 0.575
    assert capacity(9, True, 2, 12, -1e36, -1e36, -1e36, 0.3) == 2 / 3
    assert capacity(1, False, 2, 12, -1e36, -1e36, 4, 0.3) == 0.3


def test_both_build_paths_expose_colm2024_without_mutating_structure_inputs():
    for source in (NORMAL, EXTENDED):
        text = compact(source)
        assert "subroutine leaf_interception_colm2024" in text
        assert "subroutine leaf_interception_colm202x" not in text
        assert "real(r8), intent(in) :: htop, ncd, ncw, bcw" in text
        assert "canopy_storage_capacity_colm2024" in text

    assert "p,.true.,ncd_p(i),ncw_p(i),bcw_p(i),htop_p(i)" in compact(EXTENDED)


def test_colm2014_snow_unloading_has_consistent_rate_units():
    for source in (NORMAL, EXTENDED):
        text = compact(source)
        assert "ft = max(0._r8, (tleaf - tfrz) / 1.87e5_r8)" in text
        assert "tex_snow = max(0._r8, ldew_snow) * (fv+ft)" in text
        assert "ldew_snow/deltim) * (fv+ft)" not in text


def test_canopy_structure_is_persisted_and_2024_is_default():
    for name in ("ncd", "ncw", "bcw"):
        assert f"allocatable :: {name}" in INVARIANTS
        assert f"'{name}_patches'" in READIN
        assert f"'{name}_pfts'" in READIN
    assert "DEF_Interception_scheme = 8" in NAMELIST


if __name__ == "__main__":
    test_colm2024_capacity_formulas_and_class_tables()
    test_both_build_paths_expose_colm2024_without_mutating_structure_inputs()
    test_colm2014_snow_unloading_has_consistent_rate_units()
    test_canopy_structure_is_persisted_and_2024_is_default()



def compile_and_run_fortran(tmp_path: Path, program: str) -> str:
    import subprocess
    compiler = require_runnable_fortran_compiler(tmp_path)
    support = tmp_path / "support.F90"
    support.write_text(
        """
module MOD_Precision
  integer, parameter :: r8 = selected_real_kind(12)
end module
module MOD_Const_Physical
  use MOD_Precision
  real(r8), parameter :: tfrz=273.15_r8, denh2o=1000._r8, denice=917._r8
  real(r8), parameter :: cpliq=4188._r8, cpice=2117._r8, hfus=3.337e5_r8
end module
module MOD_Namelist
  use MOD_Precision
  integer :: DEF_Interception_scheme = 8
  logical :: DEF_VEG_SNOW = .true.
end module
"""
    )
    interception = (ROOT / "main" / "MOD_LeafInterception.F90").read_text()
    interception = interception.split("#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)")[0] + "\nEND MODULE MOD_LeafInterception\n"
    interception_src = tmp_path / "MOD_LeafInterception_subset.F90"
    interception_src.write_text(interception)
    main_src = tmp_path / "driver.F90"
    main_src.write_text(program)
    exe = tmp_path / "a.out"
    cmd = [
        compiler, "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-cpp", "-ffree-form", "-ffree-line-length-0",
        "-I", str(ROOT / "include"), "-I", str(tmp_path),
        str(support), str(interception_src), str(main_src),
        "-o", str(exe),
    ]
    subprocess.run(cmd, check=True, cwd=tmp_path, capture_output=True, text=True)
    return subprocess.run([str(exe)], check=True, cwd=tmp_path, capture_output=True, text=True).stdout


def test_fortran_colm2024_capacity_distinguishes_same_class_same_height_pfts(tmp_path):
    out = compile_and_run_fortran(
        tmp_path,
        """
program driver
  use MOD_Precision
  use MOD_LeafInterception, only: canopy_storage_capacity_colm2024
  implicit none
  real(r8) :: a, b
  a = canopy_storage_capacity_colm2024(0.1_r8,2._r8,1._r8,2._r8,0._r8,12._r8,3._r8,3._r8,4._r8,1,.true.)
  b = canopy_storage_capacity_colm2024(0.1_r8,2._r8,1._r8,2._r8,0._r8,12._r8,9._r8,6._r8,4._r8,1,.true.)
  if (abs(a-b) < 1.e-12_r8) error stop 1
  print *, a, b
end program
""",
    )
    vals = [float(x) for x in out.split()]
    assert vals[0] != vals[1]


def test_fortran_colm2024_snow_matches_unchanged_colm2014_capacity(tmp_path):
    program = """
program driver
  use MOD_Precision
  use MOD_LeafInterception, only: LEAF_interception_CoLM2014, LEAF_interception_CoLM2024
  implicit none
  real(r8) :: tleaf, ldew1, rain1, snow1, ldew2, rain2, snow2
  real(r8) :: pg_rain, pg_snow, qintr, qintr_rain, qintr_snow
  real(r8) :: gross_rain, gross_snow, xsc_rain, xsc_snow, smelt, frzc, heat
  tleaf = 268._r8
  ldew1 = 20._r8; rain1 = 0._r8; snow1 = 20._r8
  ldew2 = ldew1; rain2 = rain1; snow2 = snow1
  call LEAF_interception_CoLM2014(1800._r8,0.1_r8,0._r8,0._r8,0.01_r8,1._r8,2._r8,1._r8,268._r8,tleaf, &
       0._r8,0._r8,0._r8,0._r8,0._r8,68._r8,ldew1,rain1,snow1,0.1_r8,10._r8, &
       pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,gross_rain,gross_snow,xsc_rain,xsc_snow,smelt,frzc,heat)
  tleaf = 268._r8
  call LEAF_interception_CoLM2024(1800._r8,0.1_r8,0._r8,0._r8,0.01_r8,1._r8,2._r8,1._r8,268._r8,tleaf, &
       0._r8,0._r8,0._r8,0._r8,0._r8,68._r8,1,.true.,9._r8,6._r8,4._r8,12._r8, &
       ldew2,rain2,snow2,0.1_r8,10._r8,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow, &
       gross_rain,gross_snow,xsc_rain,xsc_snow,smelt,frzc,heat)
  if (abs(snow1-snow2) > 1.e-10_r8) error stop 2
  if (abs(ldew1-ldew2) > 1.e-10_r8) error stop 3
  if (abs(snow1-14.4_r8) > 1.e-10_r8) error stop 4
  print *, ldew1, snow1, ldew2, snow2
end program
"""
    out = compile_and_run_fortran(tmp_path, program)
    vals = [float(x) for x in out.split()]
    assert vals[0] == vals[2]
    assert vals[1] == vals[3]


def test_fortran_colm2014_and_2024_recover_legacy_unsplit_restart_storage(tmp_path):
    out = compile_and_run_fortran(
        tmp_path,
        """
program driver
  use MOD_Precision
  use MOD_LeafInterception, only: LEAF_interception_CoLM2014, LEAF_interception_CoLM2024
  implicit none
  integer :: scheme
  real(r8) :: tleaf, ldew, rain, snow
  real(r8) :: pg_rain, pg_snow, qintr, qintr_rain, qintr_snow
  real(r8) :: gross_rain, gross_snow, xsc_rain, xsc_snow, smelt, frzc, heat
  do scheme=1,8,7
    tleaf=280._r8; ldew=0.3_r8; rain=-1.e36_r8; snow=-1.e36_r8
    if (scheme == 1) then
      call LEAF_interception_CoLM2014(1800._r8,0.1_r8,0._r8,0._r8,0.01_r8,1._r8,2._r8,1._r8,280._r8,tleaf, &
           0._r8,0._r8,0._r8,0._r8,0._r8,68._r8,ldew,rain,snow,0.1_r8,10._r8, &
           pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,gross_rain,gross_snow,xsc_rain,xsc_snow,smelt,frzc,heat)
    else
      call LEAF_interception_CoLM2024(1800._r8,0.1_r8,0._r8,0._r8,0.01_r8,1._r8,2._r8,1._r8,280._r8,tleaf, &
           0._r8,0._r8,0._r8,0._r8,0._r8,68._r8,1,.true.,9._r8,6._r8,4._r8,12._r8, &
           ldew,rain,snow,0.1_r8,10._r8,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow, &
           gross_rain,gross_snow,xsc_rain,xsc_snow,smelt,frzc,heat)
    endif
    if (abs(ldew-0.3_r8) > 1.e-12_r8 .or. abs(rain-0.3_r8) > 1.e-12_r8 .or. abs(snow) > 1.e-12_r8) error stop 1
  enddo
  print *, 'RESTART_STORAGE_OK'
end program
""",
    )
    assert "RESTART_STORAGE_OK" in out


def test_leaf_temperature_uses_exact_ipft_index_not_global_class_height_scan():
    for path in [ROOT / "main" / "MOD_LeafTemperature.F90", ROOT / "extends" / "interception" / "MOD_LeafTemperature_Extended.F90"]:
        text = compact(path.read_text())
        assert "integer, intent(in), optional :: ipft_index" in text
        assert "ncd_p(ipft)" in text
        assert "pftclass(j)" not in text
        assert "htop_p(j)" not in text

    for path in [ROOT / "main" / "MOD_Thermal.F90", ROOT / "extends" / "interception" / "MOD_Thermal_CanopyPhase_Extended.F90"]:
        assert "ipft_index=i" in compact(path.read_text())


def test_sources_thread_2024_liquid_capacity_without_changing_snow_capacity():
    for source in (NORMAL, EXTENDED):
        text = compact(source)
        colm2014 = compact(extract_block(source, "SUBROUTINE LEAF_interception_CoLM2014", "END SUBROUTINE LEAF_interception_CoLM2014"))
        assert "satcap_rain_override=satcap_2024" in text
        assert "call leaf_interception_colm2014 (deltim,dewmx," in text
        assert "dewmx_2024" not in text
        assert "if (.not. def_veg_snow) satcap = satcap_rain" in text
        assert "satcap_snow = 48._r8*satcap" in colm2014
        assert "46._r8/bifall" not in colm2014

    for path in [
        ROOT / "main" / "MOD_LeafTemperature.F90",
        ROOT / "main" / "MOD_LeafTemperaturePC.F90",
        ROOT / "extends" / "interception" / "MOD_LeafTemperature_Extended.F90",
        ROOT / "extends" / "interception" / "MOD_LeafTemperaturePC_Extended.F90",
    ]:
        text = compact(path.read_text())
        assert "use mod_leafinterception, only: canopy_storage_capacity_colm2024" in text
        assert "colm2024_capacity_formula" not in text
        assert "satcap_rain_override" in text


def extract_block(source: str, start: str, end: str) -> str:
    a = source.lower().index(start.lower())
    b = source.lower().index(end.lower(), a) + len(end)
    return source[a:b]


def compile_and_run_dewfraction(tmp_path: Path, files: list[Path]) -> str:
    import subprocess
    compiler = require_runnable_fortran_compiler(tmp_path)
    support = tmp_path / "support_dew.F90"
    support.write_text(
        """
module MOD_Precision
  integer, parameter :: r8 = selected_real_kind(12)
end module
module MOD_Const_Physical
  use MOD_Precision
  real(r8), parameter :: tfrz=273.15_r8
end module
module MOD_Namelist
  use MOD_Precision
  integer :: DEF_Interception_scheme = 8
  logical :: DEF_VEG_SNOW = .true.
  real(r8) :: DEF_MATSIRO_CWCAP_SCALE = 1._r8
end module
"""
    )
    modules = []
    calls = []
    for idx, path in enumerate(files, 1):
        src = path.read_text()
        is_ext = "Extended" in path.name
        body = extract_block(src, "SUBROUTINE dewfraction", "END SUBROUTINE dewfraction")
        extra = ""
        if is_ext:
            extra = "\n".join(
                extract_block(src, start, end)
                for start, end in [
                    ("FUNCTION canopy_rain_capacity_for_fwet", "END FUNCTION canopy_rain_capacity_for_fwet"),
                    ("FUNCTION canopy_snow_capacity_for_fwet", "END FUNCTION canopy_snow_capacity_for_fwet"),
                    ("FUNCTION canopy_snow_wetfrac", "END FUNCTION canopy_snow_wetfrac"),
                ]
            )
        mod = f"dew_mod_{idx}"
        mod_src = tmp_path / f"{mod}.F90"
        mod_src.write_text(
            f"module {mod}\n"
            "  use MOD_Precision\n"
            "  use MOD_Namelist\n"
            "  implicit none\n"
            "contains\n"
            f"{body}\n{extra}\n"
            f"end module {mod}\n"
        )
        modules.append(mod_src)
        call_args = "1._r8,2._r8,1._r8,0.1_r8,ldew,ldew_rain,ldew_snow,fwet,fdry,cap" if not is_ext else \
            "1._r8,2._r8,1._r8,0.1_r8,280._r8,ldew,ldew_rain,ldew_snow,fwet,fdry,cap"
        snow_args = "1._r8,2._r8,1._r8,0.1_r8,ldew,ldew_rain,ldew_snow,fwet,fdry,cap" if not is_ext else \
            "1._r8,2._r8,1._r8,0.1_r8,268._r8,ldew,ldew_rain,ldew_snow,fwet,fdry,cap"
        calls.append(f"""
  subroutine run_{idx}()
    use {mod}, only: dewfraction
    real(r8) :: ldew, ldew_rain, ldew_snow, fwet, fdry, cap, expected
    cap = 0.15_r8
    DEF_Interception_scheme = 8
    DEF_VEG_SNOW = .false.
    ldew = cap; ldew_rain = cap; ldew_snow = 0._r8
    call dewfraction({call_args})
    if (abs(fwet-1._r8) > 1.e-12_r8) error stop {100+idx}
    ldew = 0.5_r8*cap; ldew_rain = ldew
    call dewfraction({call_args})
    if (fwet < 0.62_r8 .or. fwet > 0.64_r8) error stop {110+idx}
    DEF_VEG_SNOW = .true.
    ldew = cap; ldew_rain = cap; ldew_snow = 0._r8
    call dewfraction({call_args})
    if (abs(fwet-1._r8) > 1.e-12_r8) error stop {120+idx}
    ldew = 0.5_r8*cap; ldew_rain = ldew
    call dewfraction({call_args})
    if (fwet < 0.62_r8 .or. fwet > 0.64_r8) error stop {130+idx}
    DEF_Interception_scheme = 1
    DEF_VEG_SNOW = .true.
    cap = 99._r8
    ldew_rain = 0._r8; ldew_snow = 48._r8*0.1_r8*3._r8; ldew = ldew_snow
    call dewfraction({snow_args})
    if (abs(fwet-1._r8) > 1.e-12_r8) error stop {140+idx}
    ldew_snow = 0.5_r8*ldew_snow; ldew = ldew_snow
    call dewfraction({snow_args})
    expected = fwet
    DEF_Interception_scheme = 8
    cap = 0.15_r8
    call dewfraction({snow_args})
    if (abs(fwet-expected) > 1.e-12_r8) error stop {150+idx}
  end subroutine
""")
    program = tmp_path / "driver_dew.F90"
    program.write_text(
        "program driver_dew\n  use MOD_Precision\n  use MOD_Namelist\n  implicit none\n"
        + "\n".join(f"  call run_{i}()" for i in range(1, len(files)+1))
        + "\ncontains\n"
        + "\n".join(calls)
        + "\nend program\n"
    )
    exe = tmp_path / "dew.out"
    cmd = [compiler, "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-cpp", "-ffree-form", "-ffree-line-length-0", "-I", str(tmp_path), str(support), *map(str, modules), str(program), "-o", str(exe)]
    subprocess.run(cmd, check=True, cwd=tmp_path, capture_output=True, text=True)
    return subprocess.run([str(exe)], check=True, cwd=tmp_path, capture_output=True, text=True).stdout


def test_compiled_dewfraction_consumers_use_2024_capacity_and_legacy_snow(tmp_path):
    compile_and_run_dewfraction(
        tmp_path,
        [
            ROOT / "main" / "MOD_LeafTemperature.F90",
            ROOT / "main" / "MOD_LeafTemperaturePC.F90",
            ROOT / "extends" / "interception" / "MOD_LeafTemperature_Extended.F90",
            ROOT / "extends" / "interception" / "MOD_LeafTemperaturePC_Extended.F90",
        ],
    )
