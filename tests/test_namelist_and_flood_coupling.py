from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def _flat(path: Path) -> str:
    return re.sub(r"\s+", " ", re.sub(r"&\s*", " ", path.read_text().lower()))


def test_timestep_rejects_nonfinite_and_nonpositive_values() -> None:
    source = _flat(ROOT / "share" / "MOD_Namelist.F90")
    validation = source.split("error: timestep must be finite", 1)[0][-300:]

    assert "ieee_is_finite(def_simulation_time%timestep)" in validation
    assert "def_simulation_time%timestep <= 0._r8" in source


def test_flood_depth_is_recovered_from_jointly_remapped_water() -> None:
    source = _flat(ROOT / "main" / "HYDRO" / "MOD_Grid_RiverLakeFlow.F90")
    publish = source.split("subroutine publish_fldfrc_to_patches", 1)[1].split(
        "end subroutine publish_fldfrc_to_patches", 1
    )[0]

    assert "fldwat_uc(i) = fldfrc_uc(i) * max(0._r8, total_flooddepth_in(i))" in publish
    assert "flddph_patch(i) = max(0._r8, fldwat_patch(i)) / fldfrc_patch(i)" in publish

    fractions = (0.0, 0.2)
    depths = (0.0, 1.0)
    remapped_fraction = sum(fractions) / 2.0
    remapped_water = sum(f * d for f, d in zip(fractions, depths)) / 2.0
    conditional_depth = remapped_water / remapped_fraction

    assert remapped_water == 0.1
    assert conditional_depth == 1.0


def test_cama_selection_disables_competing_routing_macro(tmp_path):
    import shutil
    import subprocess
    import pytest

    cpp = shutil.which('cpp')
    if cpp is None:
        pytest.skip('C preprocessor unavailable')
    header = (ROOT / 'include/define.h').read_text()
    # Flip only the user-facing selector, not later safety undefs.
    header = header.replace('#undef CaMa_Flood', '#define CaMa_Flood', 1)
    check = tmp_path / 'check.h'
    check.write_text(header + '\n#ifdef CaMa_Flood\nCAMA_SELECTED\n#endif\n'
                     '#ifdef GridRiverLakeFlow\nGRID_SELECTED\n#endif\n')
    result = subprocess.run([cpp, '-P', str(check)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'CAMA_SELECTED' in result.stdout
    assert 'GRID_SELECTED' not in result.stdout


def test_regression_sources_are_not_ignored():
    import subprocess

    files = ['tests/test_colm2024_interception.py',
             'tests/test_sediment_precip_normalization.py', 'tests/new_regression.py']
    result = subprocess.run(['git', 'check-ignore', '--no-index', *files], cwd=ROOT,
                            capture_output=True, text=True)
    assert result.returncode == 1, result.stdout + result.stderr


def test_vegetation_snow_default_is_on_and_explicit_false_still_parseable():
    compact = _flat(ROOT / "share" / "MOD_Namelist.F90")
    assert "logical :: def_veg_snow = .true." in compact
    assert "def_veg_snow" in compact.split("namelist", 1)[1]

    # Minimal namelist parse: default remains true when omitted, explicit false overrides.
    import subprocess
    import textwrap
    import tempfile
    import shutil
    import pytest

    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran unavailable")
    with tempfile.TemporaryDirectory() as td:
        tdir = Path(td)
        program = tdir / "check_veg_snow.f90"
        program.write_text(textwrap.dedent("""
            program check_veg_snow
              implicit none
              logical :: DEF_VEG_SNOW = .true.
              namelist /nl_colm/ DEF_VEG_SNOW
              open(10, file='empty.nml', status='replace')
              write(10,'(a)') '&nl_colm /'
              close(10)
              open(10, file='empty.nml', status='old')
              read(10, nml=nl_colm)
              close(10)
              if (.not. DEF_VEG_SNOW) error stop 1
              open(11, file='false.nml', status='replace')
              write(11,'(a)') '&nl_colm DEF_VEG_SNOW = .false. /'
              close(11)
              open(11, file='false.nml', status='old')
              read(11, nml=nl_colm)
              close(11)
              if (DEF_VEG_SNOW) error stop 2
            end program
        """))
        exe = tdir / "check_veg_snow"
        result = subprocess.run([compiler, str(program), "-o", str(exe)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        result = subprocess.run([str(exe)], cwd=tdir, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
