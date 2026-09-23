"""Land tracer restart vectors require their dimension before the first write."""

from pathlib import Path


def test_transport_dimension_precedes_restart_vectors():
    source = (
        Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Rest.F90"
    ).read_text()
    body = source.split("SUBROUTINE write_land_tracer_restart", 1)[1].split(
        "END SUBROUTINE write_land_tracer_restart", 1
    )[0].split("CALL validate_land_tracer_restart_state(wa)", 1)[1]
    assert body.index("CALL write_land_tracer_transaction_marker(file_restart, 0)") < body.index(
        "CALL write_land_tracer_descriptor_metadata(file_restart)"
    ) < body.index("CALL ncio_write_vector(file_restart, 'trc_ldew_rain'")
    assert body.rindex("CALL write_land_tracer_transaction_marker(file_restart, 1)") > body.index(
        "CALL ncio_write_vector(file_restart, 'trc_leaf_iso_storage'"
    )
