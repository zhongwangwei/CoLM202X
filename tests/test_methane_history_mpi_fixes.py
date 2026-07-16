from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "main/MOD_HistVector.F90").read_text(encoding="utf-8")


def _routine(name: str) -> str:
    match = re.search(
        rf"SUBROUTINE\s+{name}\b(.*?)END\s+SUBROUTINE\s+{name}",
        SOURCE,
        re.IGNORECASE | re.DOTALL,
    )
    assert match
    return match.group(1)


def test_vector_history_pins_each_reused_tag_to_one_worker_without_global_barriers():
    for rank in ("2d", "3d", "4d"):
        body = _routine(f"aggregate_to_vector_and_write_{rank}")
        assert not re.search(r"\bcall\s+mpi_barrier\b", body, re.IGNORECASE)
        assert re.search(r"\bcall\s+mpi_send\b", body, re.IGNORECASE)
        assert re.search(r"\bcall\s+mpi_recv\b", body, re.IGNORECASE)
        assert "isrc = p_address_worker(iwork)" in body
        assert re.search(
            r"mpi_recv\s*\(\s*mesg\s*,\s*2\s*,\s*MPI_INTEGER\s*,\s*isrc\s*,",
            body,
            re.IGNORECASE,
        )
        assert "mesg(1) /= isrc" in body
        assert re.search(
            r"mpi_recv\s*\(\s*rcache.*?MPI_REAL8\s*,\s*isrc\s*,",
            body,
            re.IGNORECASE | re.DOTALL,
        )


def test_writeback_cleanup_distinguishes_header_and_data_nodes():
    source = (ROOT / "main/MOD_HistWriteBack.F90").read_text(encoding="utf-8")
    header = _routine_from(source, "hist_writeback_var_header")
    data = _routine_from(source, "hist_writeback_var")

    for body in (header, data):
        assert "allocated(HistSendBuffer%senddata)" in body
        assert re.search(r"MPI_(?:Test|Wait)all\s*\(\s*2", body, re.IGNORECASE)
        assert re.search(r"MPI_(?:Test|Wait)all\s*\(\s*3", body, re.IGNORECASE)

    assert "size(TempSendBuffer%senddata)" in header
    assert "size(TempSendBuffer%senddata)" in data
    assert re.search(r"deallocate\s*\(TempSendBuffer%senddata\)", header)
    assert re.search(r"deallocate\s*\(TempSendBuffer%senddata\)", data)
    for pointer in ("HistSendBuffer", "LastSendBuffer", "timenodes", "lasttime"):
        assert re.search(rf"pointer\s*::\s*{pointer}\s*=>\s*null\(\)", source)
    exit_body = _routine_from(source, "hist_writeback_exit")
    assert "TotalMemSize = 0" in exit_body
    time_body = _routine_from(source, "hist_writeback_latlon_time")
    assert "integer :: req (5)" in source
    assert "req_zero" not in source
    assert "dataid_zero" not in source
    assert "integer :: dataid = 0" in source
    assert "mpi_isend (lasttime%dataid" in time_body
    assert "lasttime%req(1)" in time_body
    assert re.search(r"MPI_TestAll\s*\(\s*5\s*,\s*timenodes%req", time_body)
    assert re.search(r"MPI_WaitAll\s*\(\s*5\s*,\s*timenodes%req", exit_body)


def _routine_from(source: str, name: str) -> str:
    match = re.search(
        rf"SUBROUTINE\s+{name}\b(.*?)END\s+SUBROUTINE\s+{name}",
        source,
        re.IGNORECASE | re.DOTALL,
    )
    assert match
    return match.group(1)
