from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OZONE = (ROOT / "main/MOD_Ozone.F90").read_text(encoding="utf-8")
TIME_VARIABLES = (ROOT / "main/MOD_Vars_TimeVariables.F90").read_text(
    encoding="utf-8"
)


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_initial_ozone_grid_is_mapped_before_first_timestep() -> None:
    initializer = routine(OZONE, "init_ozone_data")
    read_position = initializer.index("CALL ncio_read_block_time")
    map_position = initializer.index("CALL mg2p_ozone%grid2pset")

    assert read_position < map_position
    assert initializer.count("CALL mg2p_ozone%grid2pset") == 1


def test_unallocated_pft_coefficient_where_statements_stay_disabled() -> None:
    coefficient_names = (
        "o3coefv_sun_p",
        "o3coefv_sha_p",
        "o3coefg_sun_p",
        "o3coefg_sha_p",
    )
    active_lines = [
        line.strip()
        for line in TIME_VARIABLES.splitlines()
        if line.strip() and not line.lstrip().startswith("!")
    ]

    for name in coefficient_names:
        assert not any(line.startswith(f"WHERE ({name}") for line in active_lines)
        assert f"!      WHERE ({name} <= -1.e30) {name} = 1.0_r8" in TIME_VARIABLES
