from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LAI_READIN = (ROOT / "main/MOD_LAIReadin.F90").read_text(encoding="utf-8")
VECTOR_IO = (ROOT / "share/MOD_NetCDFVector.F90").read_text(encoding="utf-8")


def test_lai_uses_the_official_warn_and_continue_policy() -> None:
    assert LAI_READIN.count("IF (ncio_vector_exist") == 7
    assert LAI_READIN.count("Warning: LAI_") == 4
    assert LAI_READIN.count("Warning: SAI_") == 3
    assert "ncio_read_vector_complete" not in LAI_READIN

    probe = VECTOR_IO.split("logical FUNCTION ncio_vector_exist", 1)[1].split(
        "END FUNCTION ncio_vector_exist", 1
    )[0]
    assert "any_data_exists = .true." in probe
    assert "MPI_LOR" in probe
