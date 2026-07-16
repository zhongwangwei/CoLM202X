from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def test_gridded_history_tags_fit_the_mpi_guaranteed_range() -> None:
    gridded = (ROOT / "main/MOD_HistGridded.F90").read_text(encoding="utf-8")
    writeback = (ROOT / "main/MOD_HistWriteBack.F90").read_text(encoding="utf-8")

    first = int(re.search(r"hist_data_id_first\s*=\s*(\d+)", gridded).group(1))
    count = int(re.search(r"hist_data_id_count\s*=\s*(\d+)", gridded).group(1))
    last = first + count - 1

    assert first > 0
    assert last * 10 + 1 <= 32767
    assert gridded.count("hist_data_id = mod(hist_data_id-hist_data_id_first+1") == 3
    assert writeback.count("dataid*10") == 4
    assert writeback.count("dataid*10+1") == 2
