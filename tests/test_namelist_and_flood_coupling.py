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
