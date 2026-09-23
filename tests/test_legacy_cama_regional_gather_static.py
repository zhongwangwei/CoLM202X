"""CaMa exchange arrays must use indices of the full routing grid on a regional domain.

grid_concat_type concatenates only the blocks a run covers, so the segment
displacements (gdsp) count from the first covered column / row, not from column /
row 1 of the grid.  The coupler keeps its master arrays on the full NX x NY routing
grid (the input matrix and the flood credits are indexed that way), so placing a
gathered region at gdsp+1 fed a regional run's runoff to the wrong routing cells:
on a 52 x 22-cell domain every covered cell landed at (1..52, 1..22), no cell had a
routing recipient, and no water ever reached the river.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8")


def test_grid_concat_records_where_the_region_starts_in_the_full_grid():
    src = read("share/MOD_Grid.F90")
    assert re.search(r"integer\s+::\s+ilon0\s*=\s*1,\s*ilat0\s*=\s*1", src)
    assert "this%ilat0 = ilat_l" in src
    assert "this%ilon0 = ilon_w" in src


def test_every_master_array_access_of_the_gathered_region_is_offset_to_the_full_grid():
    src = read("extends/CaMa/src/MOD_CaMa_Vars.F90")
    # no access with the region-local displacement alone remains
    assert "MasterVar (xdsp+1:xdsp+xcnt" not in src
    gathers = re.findall(r"MasterVar \(cama_master_cols\(xdsp,xcnt,size\(MasterVar,1\)\), "
                         r"ydsp\+cama_gather%ilat0:ydsp\+cama_gather%ilat0\+ycnt-1\)", src)
    assert len(gathers) == 4, "two gathers to the master and two scatters from it"
    body = src.split("FUNCTION cama_master_cols", 1)[1].split("END FUNCTION cama_master_cols", 1)[0]
    assert "mod(dsp + cama_gather%ilon0 - 1 + k - 1, n) + 1" in body, "columns wrap at the dateline"
