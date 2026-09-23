from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_pixelset_gather_receives_into_contiguous_rank_one_buffer():
    source = (ROOT / "share/MOD_Pixelset.F90").read_text()
    routine = source.split("SUBROUTINE vec_gather_scatter_set", 1)[1].split(
        "END SUBROUTINE vec_gather_scatter_set", 1
    )[0]
    assert "allocate (gathered_counts(0:p_np_group-1))" in routine
    assert "gathered_counts, 1, MPI_INTEGER" in routine
    assert "this%vecgs%vcnt(:,iblk,jblk) = gathered_counts" in routine
    assert "this%vecgs%vcnt(:,iblk,jblk), 1, MPI_INTEGER" not in routine


def test_tracer_restart_nonworker_passes_allocated_zero_extent_arrays():
    source = (ROOT / "main/MOD_Vars_TimeVariables.F90").read_text()
    routine = source.split("SUBROUTINE WRITE_TimeVariables", 1)[1].split(
        "END SUBROUTINE WRITE_TimeVariables", 1
    )[0]
    assert "empty_patch(0), empty_soilsnow(maxsnl+1:nl_soil,0)" in routine
    assert "IF (allocated(ldew_rain)) THEN" in routine
    assert "IF (numpatch < 0) ERROR STOP 'invalid tracer restart patch count'" in routine
    assert "IF (p_is_worker) THEN\n               IF (numpatch /= 0) ERROR STOP" in routine
    assert "incomplete tracer restart water state" in routine
    assert "tracer restart water shape mismatch" in routine
    assert "CALL write_tracer_restart_all(file_restart, maxsnl, nl_soil, 0," in routine
    assert "empty_patch, empty_patch, empty_soilsnow, empty_soilsnow" in routine
