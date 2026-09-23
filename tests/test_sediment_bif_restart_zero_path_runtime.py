"""BIF gross matrices are ephemeral and need no zero-path restart I/O."""

from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"


def test_bif_gross_restart_uses_zero_init_and_boundary_assertion() -> None:
    source = SOURCE.read_text()
    reader = source.split("SUBROUTINE read_sediment_restart(", 1)[1].split(
        "END SUBROUTINE read_sediment_restart", 1
    )[0]
    writer = source.split("SUBROUTINE write_sediment_restart(", 1)[1].split(
        "END SUBROUTINE write_sediment_restart", 1
    )[0]
    fields = (
        "sed_acc_bif_forward", "sed_acc_bif_reverse",
        "sed_acc_bif_forward_time", "sed_acc_bif_reverse_time",
    )
    for field in fields:
        assert f"{field} = 0._r8" in reader
        assert f"any({field} /= 0._r8)" in writer
        assert f"'{field}'" not in reader + writer
    assert "vector_read_matrix_and_scatter" not in reader
    assert "vector_gather_matrix_to_master" not in writer
    assert "SED_RESTART_SCHEMA_VERSION = 5" in source
