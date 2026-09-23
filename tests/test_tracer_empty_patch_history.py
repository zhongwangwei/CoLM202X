"""History ranks without patches still need safe zero-length tracer rows."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def test_allocator_keeps_tracer_rows_on_empty_patch_ranks():
    vars_src = (ROOT / "main/TRACER/MOD_Tracer_Vars.F90").read_text()
    body = vars_src.split("SUBROUTINE allocate_Tracer_Vars", 1)[1].split(
        "END SUBROUTINE allocate_Tracer_Vars", 1
    )[0]
    assert "IF (ntracers <= 0) RETURN" in body
    assert "numpatch <= 0" not in body
    for row in ("a_trc_precip", "a_water_precip", "a_trc_rnof", "a_water_rnof"):
        assert re.search(rf"allocate\(\s*{row}\s*\(ntracers, numpatch\)\)", body)
