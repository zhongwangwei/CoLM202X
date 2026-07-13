from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
RIVER = ROOT / "main" / "TRACER" / "MOD_Tracer_RiverLake.F90"


def test_all_generic_river_loops_use_the_transport_contract():
    source = RIVER.read_text(encoding="utf-8")

    assert "tracer_uses_land_water_transport" in source.split("IMPLICIT NONE", 1)[0]
    for helper in (
        "riverlake_legacy_tracer_count_meta",
        "riverlake_legacy_tracer_namehash_meta",
    ):
        source = source.replace(_function_body(source, helper), "")

    tracer_loops = re.findall(
        r"DO\s+itrc\w*\s*=\s*1\s*,\s*ntracers",
        source,
        flags=re.IGNORECASE,
    )
    guarded_loops = re.findall(
        r"DO\s+(?P<index>itrc\w*)\s*=\s*1\s*,\s*ntracers\s*"
        r"IF\s*\(\.not\.\s*tracer_uses_land_water_transport\("
        r"(?P=index)\)\)\s*CYCLE",
        source,
        flags=re.IGNORECASE,
    )
    assert tracer_loops
    assert len(guarded_loops) == len(tracer_loops)


def _function_body(source: str, name: str) -> str:
    match = re.search(
        rf"(?:integer|real\(r8\))\s+FUNCTION\s+{name}\s*\(\s*\)(.*?)"
        rf"END\s+FUNCTION\s+{name}",
        source,
        flags=re.IGNORECASE | re.DOTALL,
    )
    assert match, f"missing function {name}"
    return match.group(1)


def test_restart_writer_keeps_the_new_transport_capability_signature():
    source = RIVER.read_text(encoding="utf-8")
    writer = source.split("SUBROUTINE write_tracer_restart", 1)[1]
    writer = writer.split("END SUBROUTINE write_tracer_restart", 1)[0]

    assert "riverlake_tracer_count_meta()" in writer
    assert "riverlake_tracer_namehash_meta()" in writer
    assert "riverlake_legacy_tracer_count_meta()" not in writer
    assert "riverlake_legacy_tracer_namehash_meta()" not in writer


def test_restart_reader_accepts_only_complete_new_or_legacy_signatures():
    source = RIVER.read_text(encoding="utf-8")
    current_count = _function_body(source, "riverlake_tracer_count_meta")
    current_hash = _function_body(source, "riverlake_tracer_namehash_meta")
    legacy_count = _function_body(source, "riverlake_legacy_tracer_count_meta")
    legacy_hash = _function_body(source, "riverlake_legacy_tracer_namehash_meta")
    reader = source.split("SUBROUTINE read_tracer_restart", 1)[1]
    reader = reader.split("END SUBROUTINE read_tracer_restart", 1)[0]

    assert "tracer_uses_land_water_transport(itrc)" in current_count
    assert "tracer_uses_land_water_transport(itrc)" in current_hash
    assert "tracer_is_particle(itrc)" in legacy_count
    assert "tracer_is_particle(itrc)" in legacy_hash

    assert "meta_matches_current" in reader
    assert "meta_matches_legacy" in reader
    assert re.search(
        r"legacy_meta_complete\s*=\s*has_trc_n_meta\s*\.and\.\s*has_namehash_meta",
        reader,
        flags=re.IGNORECASE,
    )
    reader_without_continuations = re.sub(r"\s*&\s*", " ", reader)
    assert re.search(
        r"meta_bad\s*=\s*meta_bad\s*\.or\.\s*\(\.not\.\s*"
        r"meta_matches_current\s*\.and\.\s*\.not\.\s*"
        r"\(legacy_meta_complete\s*\.and\.\s*meta_matches_legacy\)\)",
        reader_without_continuations,
        flags=re.IGNORECASE,
    )
    assert re.search(
        r"mpi_allreduce\s*\(\s*MPI_IN_PLACE\s*,\s*meta_bad\s*,\s*1\s*,\s*"
        r"MPI_LOGICAL\s*,\s*MPI_LOR",
        reader,
        flags=re.IGNORECASE,
    )

    accepted = lambda current, legacy, complete: current or (complete and legacy)
    assert accepted(current=True, legacy=False, complete=True)
    assert accepted(current=False, legacy=True, complete=True)
    assert not accepted(current=False, legacy=True, complete=False)
    assert not accepted(current=False, legacy=False, complete=True)
