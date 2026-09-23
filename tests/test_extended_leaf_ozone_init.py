from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_extended_leaf_ozone_off_coefficients_are_set_before_first_use():
    for name, suffix in (
        ("MOD_LeafTemperature_Extended.F90", ""),
        ("MOD_LeafTemperaturePC_Extended.F90", "(i)"),
    ):
        source = (ROOT / "extends/interception" / name).read_text()
        anchor = "fevpl_bef = 0." if not suffix else "fevpl_bef(:) = 0."
        begin = source.index(anchor)
        init = source.index("IF (.not. DEF_USE_OZONESTRESS) THEN", begin)
        first_use = source.index(f"gs0sun{suffix} =", begin)
        assert begin < init < first_use
        early = source[init : source.index("ENDIF", init)]
        for key in ("o3coefv_sun", "o3coefg_sun", "o3coefv_sha", "o3coefg_sha"):
            assert f"{key}{suffix} = 1.0_r8" in early
