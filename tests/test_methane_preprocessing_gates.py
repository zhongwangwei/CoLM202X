from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def source(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def test_methane_surface_inputs_and_bgc_off_state_are_fail_closed() -> None:
    mksrfdata = source("mksrfdata/MKSRFDATA.F90")
    helper_source = source(
        "main/TRACER/MOD_Tracer_Reactive_Methane_Preprocessing.F90"
    )
    defs = source("main/TRACER/MOD_Tracer_Defs.F90")

    assert "CALL methane_preprocessing_requirements" in mksrfdata
    assert "IF (requires_lake_soilc)" in mksrfdata
    assert "IF (requires_spatial_ph)" in mksrfdata

    helper = helper_source.split(
        "SUBROUTINE methane_preprocessing_requirements", 1
    )[1].split("END SUBROUTINE methane_preprocessing_requirements", 1)[0]
    assert "CALL tracer_defs_init" in helper
    assert "igas_ch4 = tracer_index_for_name('CH4', 'METHANE')" in helper
    assert "IF (.not. methane_is_active()) RETURN" in helper
    assert "IF (.not. tracer_is_gas(igas_ch4))" in helper
    assert "CH4/METHANE preprocessing descriptor must use family=gas" in helper
    assert "CALL tracer_param_file_for_index (igas_ch4, 'CH4,METHANE'" in helper
    assert "CALL read_methane_namelist" in helper
    assert "requires_lake_soilc = DEF_METHANE%allowlakeprod" in helper
    assert "requires_spatial_ph = DEF_METHANE%use_spatial_ph" in helper

    bgc_off_guard = defs.split("#ifndef BGC", 1)[1].split("#endif", 1)[0]
    assert "'CH4'" in bgc_off_guard and "'METHANE'" in bgc_off_guard
    assert "tracers(itrc)%state_owner == STATE_OWNER_PROVIDER" in bgc_off_guard
    assert "species-owned CH4 requires compiling with BGC" in bgc_off_guard
