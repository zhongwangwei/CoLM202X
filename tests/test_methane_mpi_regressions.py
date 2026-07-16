from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def source(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def routine(text: str, name: str) -> str:
    return text.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


def test_giems_rejects_wrong_rank_and_dimension_order_collectively() -> None:
    giems = source("main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90")
    read = routine(giems, "read_methane_giems")

    assert "ndims /= 3" in read
    for index, name, length in (
        (1, "longitude", "nlon"),
        (2, "latitude", "nlat"),
        (3, "time", "ntime"),
    ):
        assert f"vdims({index})" in read
        assert re.search(
            rf"trim\(dname\)\s*/=\s*'{name}'.*?dlen\s*/=\s*{length}",
            read,
            re.IGNORECASE | re.DOTALL,
        )

    error_bcast = read.index("MPI_Bcast(giems_metadata_error")
    stop = read.index("CALL CoLM_stop", error_bcast)
    assert error_bcast < stop


def test_giems_fails_closed_and_verifies_the_hardcoded_month_axis() -> None:
    giems = source("main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90")
    read = routine(giems, "read_methane_giems")
    time_check = routine(giems, "validate_giems_time_axis")

    assert "giems_expected_months = 348" in read
    assert "ntime /= giems_expected_months" in read
    assert "CALL validate_giems_time_axis" in read
    assert "GIEMS file not found" in read
    assert "GIEMS open failed" in read
    for failure in ("GIEMS file not found", "GIEMS open failed"):
        branch = read.split(failure, 1)[1][:400]
        assert "giems_metadata_error = 1" in branch

    assert "nf90_inq_varid(ncid, 'time', time_vid)" in time_check
    assert "nf90_get_var(ncid, time_vid, time_values)" in time_check
    assert "nf90_get_att(ncid, time_vid, 'units', time_units)" in time_check
    assert "days since 1992-01-01" in time_check
    assert "nf90_get_att(ncid, time_vid, 'calendar', time_calendar)" in time_check
    assert "proleptic_gregorian" in time_check
    assert "nf90_get_att(ncid, time_vid, 'time_resolution', time_resolution)" in time_check
    assert "monthly_average" in time_check
    assert "month_days(12)" in time_check
    assert "expected_day = expected_day + days_this_month" in time_check
    assert "ieee_is_finite(time_values(i))" in time_check


def test_giems_counts_documented_zero_flags_in_climatology_and_rejects_corruption() -> None:
    giems = source("main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90")
    read = routine(giems, "read_methane_giems")

    value_branch = read.split("v = patch_values(month_in_chunk, ipatch)", 1)[1].split(
        "deallocate(patch_values, requested_values)", 1
    )[0]
    assert "ieee_is_nan(v)" in value_branch
    for flag in ("-999._r4", "-998._r4", "-997._r4"):
        assert flag in value_branch
    assert "giems_value_error = 0" in read
    assert "giems_value_error = 1" in value_branch
    assert "ccnt(mo, ipatch) = ccnt(mo, ipatch) + 1" in value_branch
    assert value_branch.index("ccnt(mo, ipatch) = ccnt(mo, ipatch) + 1") > value_branch.index(
        "ENDIF"
    )
    assert "MPI_Allreduce(giems_value_error" in value_branch
    assert "CALL CoLM_stop" in value_branch


def test_giems_chunks_monthly_reads_and_distributions_by_twelve() -> None:
    giems = source("main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90")
    read = routine(giems, "read_methane_giems")
    chunk = read.split("DO t = 1, ntime, chunk_max", 1)[1].split(
        "! Finalize climatology", 1
    )[0]

    assert "giems_chunk_months = 12" in read
    assert "giems_max_packed_values = 16 * 1024 * 1024" in read
    assert "chunk_n = min(chunk_max, ntime - t + 1)" in chunk
    assert re.search(
        r"allocate\(requested_values\(chunk_n,\s*max\(1,\s*total_requests\)\)\)",
        chunk,
    )
    assert re.search(
        r"allocate\(patch_values\(chunk_n,\s*max\(1,\s*numpatch\)\)\)", chunk
    )
    assert "chunk_counts(:) = request_counts(:) * chunk_n" in chunk
    assert "chunk_displs(:) = request_displs(:) * chunk_n" in chunk
    assert "chunk_n * numpatch" in chunk

    slab_read = chunk.index("nf90_get_var")
    error_bcast = chunk.index("MPI_Bcast(giems_block_error")
    scatter = chunk.index("MPI_Scatterv(requested_values")
    assert slab_read < error_bcast < scatter
    assert chunk.count("MPI_Bcast(giems_block_error") == 1
    assert chunk.count("MPI_Scatterv(requested_values") == 1

    assert re.search(
        r"chunk_max\s*=\s*min\(chunk_max,\s*max\(1,\s*"
        r"giems_max_packed_values\s*/\s*max\(1,\s*total_requests\)\)\)",
        read,
    )
    assert "MPI_Bcast(chunk_max" in read

    chunks = [min(12, 26 - start + 1) for start in range(1, 27, 12)]
    assert chunks == [12, 12, 2]


def test_giems_checks_chunked_mpi_counts_before_distribution() -> None:
    giems = source("main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90")
    read = routine(giems, "read_methane_giems")

    assert "total_requests > huge(total_requests) - request_counts(irequest)" in read
    assert "total_requests > huge(total_requests) / chunk_max" in read
    assert "numpatch > huge(numpatch) / chunk_max" in read
    allreduce = read.index("MPI_Allreduce(giems_count_error")
    stop = read.index("CALL CoLM_stop", allreduce)
    gather = read.index("MPI_Gatherv(pixel_index")
    assert allreduce < stop < gather


def test_giems_rejects_unrepresentable_pixels_and_invalid_patch_mapping() -> None:
    giems = source("main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90")
    read = routine(giems, "read_methane_giems")

    overflow_guard = "int(nlat, i8) * int(nlon, i8) > int(huge(0), i8)"
    assert overflow_guard in read
    assert read.index(overflow_guard) < read.index(
        "pixel_index(ipatch) = (best_iy(ipatch) - 1) * nlon + best_ix(ipatch)"
    )

    assert "ieee_is_finite(patchlatr_in(ipatch))" in read
    assert "ieee_is_finite(patchlonr_in(ipatch))" in read
    mapping_reduce = read.index("MPI_Allreduce(giems_mapping_error")
    mapping_stop = read.index("CALL CoLM_stop", mapping_reduce)
    pixel_allocation = read.index("allocate(pixel_index", mapping_stop)
    assert mapping_reduce < mapping_stop < pixel_allocation

    # The current API has centres, not patch footprints: keep this limitation
    # explicit so a future change does not claim a conservative remap.
    declaration = read.split("integer, intent(in)  :: numpatch", 1)[0]
    assert "patchlatr_in(:), patchlonr_in(:)" in declaration
    assert "underlying mesh-pixel bounds and area weights" in giems


def test_giems_chunk_layout_preserves_month_patch_order() -> None:
    ntime = 26
    rank_counts = [2, 0, 3]
    rank_displs = [0, 2, 2]
    monthly = [[1000 * month + patch for patch in range(5)] for month in range(ntime)]

    for chunk_max in (12, 3):
        received = [[[] for _ in range(count)] for count in rank_counts]
        for start in range(0, ntime, chunk_max):
            chunk_n = min(chunk_max, ntime - start)
            # Fortran requested_values(chunk_n,total_requests) is contiguous by patch.
            packed = [
                monthly[start + month][patch]
                for patch in range(5)
                for month in range(chunk_n)
            ]
            for rank, count in enumerate(rank_counts):
                begin = rank_displs[rank] * chunk_n
                local = packed[begin : begin + count * chunk_n]
                for patch in range(count):
                    received[rank][patch].extend(
                        local[patch * chunk_n : (patch + 1) * chunk_n]
                    )

        assert received[0] == [
            [row[0] for row in monthly],
            [row[1] for row in monthly],
        ]
        assert received[1] == []
        assert received[2] == [
            [row[2] for row in monthly],
            [row[3] for row in monthly],
            [row[4] for row in monthly],
        ]


def test_giems_requires_worker_patch_coordinates_before_collective_read() -> None:
    methane = source("main/TRACER/MOD_Tracer_Reactive_Methane.F90")
    require = routine(methane, "require_methane_giems_patch_coords")

    assert methane.count("CALL require_methane_giems_patch_coords") == 2
    assert "p_is_worker .and. local_numpatch > 0" in require
    assert ".not. allocated(patchlatr)" in require
    assert ".not. allocated(patchlonr)" in require
    assert "size(patchlatr) < local_numpatch" in require
    assert "size(patchlonr) < local_numpatch" in require
    assert "MPI_Allreduce" in require
    assert "MPI_LOR" in require
    assert "CALL CoLM_stop" in require

    init = routine(methane, "ch4_reactive_init")
    reload = routine(methane, "ch4_reactive_reload_lulcc_inputs")
    for body, count_name in ((init, "numpatch"), (reload, "nnew")):
        guarded_read = body.split(
            f"CALL require_methane_giems_patch_coords ({count_name})", 1
        )[1]
        assert f"IF (p_is_worker .and. {count_name} > 0) THEN" in guarded_read
        assert "allocated(patchlatr) .and. allocated(patchlonr)" not in guarded_read


def test_missing_bgc_warning_is_emitted_by_one_reachable_worker() -> None:
    driver = source("main/TRACER/MOD_Tracer_Reactive_Methane_Driver.F90")
    warning = driver.split("IF (patchtype == 0 .and. .not. bgc_inputs_ready) THEN", 1)[
        1
    ].split("ENDIF", 2)[0]

    assert "p_iam_worker == 0" in warning
    assert "p_is_master" not in warning
    assert "warned_missing_bgc_inputs" in warning


def test_methane_scoped_runtime_checks_use_colm_stop_not_local_abort() -> None:
    files = (
        "main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90",
        "main/TRACER/MOD_Tracer_Reactive_Methane_Driver.F90",
        "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90",
        "main/TRACER/MOD_Tracer_Reactive_Methane_Hist.F90",
        "main/TRACER/MOD_Tracer_Reactive.F90",
    )
    combined = "\n".join(source(path) for path in files)
    assert not re.search(r"(?i)\bcall\s+abort\b", combined)
    assert "CALL CoLM_stop" in source(
        "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90"
    )


def test_standard_methane_run_uses_production_history_and_fails_on_bad_lake_depth() -> None:
    standard = source("run/standard_ch4_parameter.nml")

    assert "DEF_METHANE%ch4_history_vars  = 'core'" in standard
    assert "DEF_METHANE%lake_zero_depth_fatal = .true." in standard


def test_core_history_skips_unselected_sumarea_phases_and_temp_vectors() -> None:
    history = source("main/TRACER/MOD_Tracer_Reactive_Methane_Hist.F90")

    for gate in ("need_lake_history", "need_extra_history", "need_soil_history"):
        assert gate in history
    assert re.search(
        r"IF \(need_lake_history.*?CALL mp2g_hist%get_sumarea \(sumarea, filter\)",
        history,
        re.DOTALL,
    )
    assert re.search(
        r"IF \(need_extra_history.*?CALL mp2g_hist%get_sumarea \(sumarea, filter\)",
        history,
        re.DOTALL,
    )
    assert re.search(
        r"IF \(need_soil_history.*?CALL mp2g_hist%get_sumarea \(sumarea, filter\)",
        history,
        re.DOTALL,
    )

    derived = history.split("IF (need_derived_ch4) THEN", 1)[1].split(
        "ENDIF", 1
    )[0]
    assert "IF (need_active_without_lake) allocate" in derived
    assert "IF (need_lake_intensive) allocate" in derived
    assert "IF (need_rice_intensive) THEN" in derived


def test_global_lake_history_uses_patch_specific_accumulation_counts() -> None:
    history = source("main/TRACER/MOD_Tracer_Reactive_Methane_Hist.F90")

    assert re.search(
        r"hist_ch4_global_with_lake\(ipatch\)\s*=\s*&?\s*"
        r"a_methane_surf_flux_tot\(ipatch\)\s*/\s*a_methane_acc_num\(ipatch\)",
        history,
        re.DOTALL,
    )
    assert re.search(
        r"hist_ch4_global_with_lake\(ipatch\)\s*=\s*&?\s*"
        r"a_methane_surf_flux_tot_lake\(ipatch\)\s*/\s*"
        r"a_methane_acc_num_lake\(ipatch\)",
        history,
        re.DOTALL,
    )
    assert re.search(
        r"CALL\s+write_history_variable_2d\s*\(\.true\.,\s*&?\s*"
        r"hist_ch4_global_with_lake.*?"
        r"'f_methane_surf_flux_global_total_with_lake'.*?"
        r"acc_num\s*=\s*hist_ch4_acc_one",
        history,
        re.DOTALL | re.IGNORECASE,
    )

    diagnostics = {
        "phys": (
            "hist_ch4_global_phys_with_lake",
            "a_methane_surf_flux_tot_phys",
            "f_methane_surf_flux_global_phys_with_lake",
        ),
        "balance": (
            "hist_ch4_global_balance_with_lake",
            "a_methane_balance_residual",
            "f_methane_balance_residual_global_with_lake",
        ),
        "clip": (
            "hist_ch4_global_clip_credit_with_lake",
            "a_methane_ch4_clip_credit",
            "f_methane_ch4_clip_credit_global_with_lake",
        ),
    }
    for _, (temporary, accumulator, field) in diagnostics.items():
        assert f"allocate ({temporary}(numpatch))" in history
        assert re.search(
            rf"{temporary}\(ipatch\)\s*=\s*&?\s*"
            rf"{accumulator}\(ipatch\)\s*/\s*a_methane_acc_num\(ipatch\)",
            history,
            re.DOTALL,
        )
        assert re.search(
            rf"{temporary}\(ipatch\)\s*=\s*&?\s*"
            rf"{accumulator}\(ipatch\)\s*/\s*a_methane_acc_num_lake\(ipatch\)",
            history,
            re.DOTALL,
        )
        assert re.search(
            rf"CALL\s+write_history_variable_2d\s*\(\.true\.,\s*&?\s*"
            rf"{temporary}.*?'{field}'.*?acc_num\s*=\s*hist_ch4_acc_one",
            history,
            re.DOTALL | re.IGNORECASE,
        )


def test_methane_history_returns_before_collectives_when_disabled() -> None:
    history = source("main/TRACER/MOD_Tracer_Reactive_Methane_Hist.F90")
    body = routine(history, "methane_reactive_history")
    first_sumarea = body.index("get_sumarea")

    assert body.index(".not. DEF_METHANE%write_ch4_history") < first_sumarea
    assert body.index("methane_history_selector_is_off") < first_sumarea
    assert re.search(
        r"IF \(need_active_history \.and\. HistForm == 'Gridded'\) THEN\s+"
        r"CALL mp2g_hist%get_sumarea \(sumarea, filter\)",
        body,
    )
    selector = history.split("FUNCTION methane_history_all_land_only", 1)[1].split(
        "END FUNCTION methane_history_all_land_only", 1
    )[0]
    for name in (
        "f_methane_surf_flux_tot",
        "f_methane_surf_flux_wetland",
        "f_methane_surf_flux_soil",
        "f_methane_surf_flux_lake",
        "f_methane_surf_flux_rice",
        "f_methane_surf_flux_active_total_without_lake",
        "f_methane_surf_flux_global_total_with_lake",
        "f_methane_surf_flux_global_phys_with_lake",
        "f_methane_balance_residual_global_with_lake",
        "f_methane_ch4_clip_credit_global_with_lake",
    ):
        assert name in selector


def test_history_writeback_memory_limit_counts_real8_bytes_prospectively() -> None:
    writeback = source("main/MOD_HistWriteBack.F90")
    body = writeback.split("SUBROUTINE hist_writeback_var (", 1)[1].split(
        "END SUBROUTINE hist_writeback_var", 1
    )[0]

    assert re.search(
        r"HistRealBytes\s*=\s*storage_size\(0\._r8\)\s*/\s*8", writeback
    )
    assert "TotalMemSize + buffer_mem_size > MaxHistMemSize" in body
    assert "int(size(TempSendBuffer%senddata), 8) * HistRealBytes" in body
    assert "TotalMemSize = TotalMemSize + buffer_mem_size" in body


def test_single_reactive_history_callback_avoids_defensive_deep_copies() -> None:
    reactive = source("main/TRACER/MOD_Tracer_Reactive.F90")
    history = routine(reactive, "tracer_reactive_history")

    assert "active_history_callbacks == 1" in history
    direct = history.split("active_history_callbacks == 1", 1)[1].split("ELSE", 1)[0]
    assert "%history (file_hist, itime_in_file, sumarea, filter" in direct
