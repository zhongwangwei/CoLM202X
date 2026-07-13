from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def source(relative_path: str) -> str:
    return (ROOT / relative_path).read_text()


def test_effective_conductance_does_not_feed_next_timestep() -> None:
    physics = source("main/TRACER/MOD_Tracer_Reactive_Methane_Physics.F90")

    assert "grnd_methane_cond_base, grnd_methane_cond_effective" in physics
    assert "grnd_methane_cond_effective = spec_grnd_cond(1)" in physics
    assert "grnd_cond_base = max(grnd_methane_cond" not in physics
    assert physics.count("DEF_METHANE%grnd_methane_cond_default, grnd_methane_cond_") == 2


def test_lulcc_collective_reload_is_outside_worker_remap() -> None:
    driver = source("main/LULCC/MOD_Lulcc_Driver.F90")
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")

    worker_remap = driver.index("IF (p_is_worker .and. allocated(patchclass)")
    reload_call = driver.index("CALL tracer_reactive_reload_lulcc_inputs ()")
    assert worker_remap < reload_call
    assert "ENDIF\n      ! GIEMS broadcasts" in driver[worker_remap:reload_call]

    remap_start = methane.index("SUBROUTINE ch4_reactive_remap_lulcc_state")
    remap_end = methane.index("END SUBROUTINE ch4_reactive_remap_lulcc_state", remap_start)
    reload_start = methane.index("SUBROUTINE ch4_reactive_reload_lulcc_inputs")
    reload_end = methane.index("END SUBROUTINE ch4_reactive_reload_lulcc_inputs", reload_start)
    remap_body = methane[remap_start:remap_end]
    reload_body = methane[reload_start:reload_end]
    assert "read_methane_giems" not in remap_body
    assert "read_methane_ph_patch" not in remap_body
    assert "read_methane_giems" in reload_body
    assert "read_methane_ph_patch" in reload_body


def test_reactive_lulcc_area_order_is_new_then_old() -> None:
    files = (
        "main/TRACER/MOD_Tracer_Reactive.F90",
        "main/TRACER/MOD_Tracer_Reactive_Methane.F90",
        "main/TRACER/MOD_Tracer_Reactive_Methane_State.F90",
        "main/TRACER/MOD_Tracer_Reactive_Methane_Microbes.F90",
    )

    for relative_path in files:
        text = source(relative_path)
        assert "old_patch_area, new_patch_area" not in text
    assert "landpatch%pctshared, landpatch_%pctshared" in source(
        "main/LULCC/MOD_Lulcc_Driver.F90"
    )
