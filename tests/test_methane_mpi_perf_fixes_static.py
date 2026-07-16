from pathlib import Path
import math
import re


ROOT = Path(__file__).resolve().parents[1]
GIEMS = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_GIEMS.F90").read_text()
ACC = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_AccFlux.F90").read_text()
HIST = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_Hist.F90").read_text()


def _compact(source: str) -> str:
    source = re.sub(r"\s*&\s*", " ", source.lower())
    return re.sub(r"\s+", " ", source)


def test_giems_rejects_every_wrong_variable_dimension_collectively():
    source = _compact(GIEMS)
    assert "ndims /= 3" in source
    assert "trim(dname) /= 'longitude'" in source
    assert "trim(dname) /= 'latitude'" in source
    assert "trim(dname) /= 'time'" in source
    assert "metadata = [giems_metadata_error, ntime, nlat, nlon]" in source
    assert "call mpi_bcast(metadata, 4, mpi_integer" in source
    assert "giems_active = giems_metadata_error == 0" in source
    assert "call colm_stop ('error: invalid giems metadata" in source
    rank_query = source.index("nf90_inquire_variable(ncid, vid, ndims=ndims)")
    dimid_query = source.index("nf90_inquire_variable(ncid, vid, dimids=vdims)")
    assert rank_query < source.index("ndims /= 3") < dimid_query


def test_giems_deduplicates_rank_local_pixels_before_network_distribution():
    source = _compact(GIEMS)
    assert "pixel_to_unique" in source
    assert "patch_to_unique" in source
    assert "pixel_index_unique" in source
    assert "call mpi_gather(n_unique" in source
    assert "call mpi_gatherv(pixel_index_unique, n_unique" in source
    assert "call mpi_scatterv(requested_values, chunk_counts, chunk_displs" in source
    assert "unique_values, chunk_n * n_unique" in source
    assert "v = unique_values(month_in_chunk, patch_to_unique(ipatch))" in source


def test_giems_mpi_buffers_cover_zero_patch_ranks_and_match_communicator():
    source = _compact(GIEMS)
    assert "allocate(pixel_index(max(1, numpatch)))" not in source
    assert "call mpi_comm_size(p_comm_glb, comm_size" in source
    assert "comm_size /= p_np_glb" in source
    assert "allocate(request_counts(comm_size), request_displs(comm_size))" in source
    assert "allocate(pixel_index_unique(max(1, numpatch)))" in source
    assert "allocate(unique_values(chunk_n, max(1, n_unique)))" in source
    assert "allocate(all_pixel_index(max(1, total_requests)))" in source
    assert "allocate(requested_values(chunk_n, max(1, total_requests)))" in source
    assert source.count("call check_giems_mpi(") >= 10
    assert "subroutine check_giems_mpi(ierr, operation)" in source


def test_giems_uses_source_precision_for_time_series_and_fewer_metadata_collectives():
    source = GIEMS.lower()
    assert re.search(
        r"real\(r4\),\s*allocatable,\s*public\s*::\s*giems_ts_wetland_frac",
        source,
    )
    assert re.search(
        r"real\(r8\),\s*allocatable,\s*public\s*::\s*giems_clim_wetland_frac",
        source,
    )
    assert "giems_chunk_months = 12" in source
    assert "do t = 1, ntime, chunk_max" in source
    assert "chunk_counts(:) = request_counts(:) * chunk_n" in source


def test_giems_nearest_pixel_keeps_the_two_dimensional_five_degree_limit():
    source = _compact(GIEMS)
    assert "dlat_min = dmin" in source
    assert "dlon_min = dmin" in source
    assert "sqrt(dlat_min**2 + dlon_min**2) <= 5._r8" in source
    # Axis-wise checks would accept this point, while the intended 2-D radius
    # correctly rejects it.
    assert 4.0 <= 5.0 and 4.0 <= 5.0 and math.hypot(4.0, 4.0) > 5.0


def test_accumulation_is_gated_by_selected_history_groups():
    source = _compact(ACC)
    assert "mhist_on => methane_history_enabled" in source
    for flag in (
        "need_main_layers",
        "need_main_surface",
        "need_main_totals",
        "need_main_concentrations",
        "need_unsat_sat",
        "need_lake",
        "need_extras",
        "need_microbes",
    ):
        assert flag in source
        assert f"if ({flag}) then" in source
    assert "if (.not. need_any) return" in source
    assert "logical, save :: history_groups_initialized = .false." in source
    assert "trim(cached_history_vars) /= trim(def_methane%ch4_history_vars)" in source


def test_every_history_selector_has_an_accumulation_dependency():
    history_names = set(re.findall(r"mhist_on\('([^']+)'\)", HIST, re.I))
    start = ACC.index("need_main_layers = any([")
    end = ACC.index("cached_write_history = DEF_METHANE%write_ch4_history", start)
    dependency_block = ACC[start:end]
    dependency_names = set(
        re.findall(r"mhist_on\('([^']+)'\)", dependency_block, re.I)
    )
    assert history_names <= dependency_names


def test_history_restart_fields_are_symmetric_and_include_selector_marker():
    # Runtime selection may skip arithmetic.  Every persisted accumulator field
    # must still have a symmetric reader; the selector marker is intentionally
    # part of that schema and avoids relying on a brittle fixed field count.
    write_names = set(
        re.findall(r"ncio_write_vector\s*\([^,]+,\s*'([^']+)'", ACC, re.I)
    )
    read_names = set(
        re.findall(r"ncio_read_vector\s*\([^,]+,\s*'([^']+)'", ACC, re.I)
    )
    assert write_names == read_names
    assert "ch4_acc_history_selector_hash" in write_names


def test_restart_selector_change_starts_a_clean_history_window():
    source = _compact(ACC)
    assert "integer(i8) function methane_history_selector_fingerprint()" in source
    assert "ishftc(methane_history_selector_fingerprint, 7, 52)" in source
    assert "ncio_vector_var_present( file_restart, 'ch4_acc_history_selector_hash'" in source
    assert "history_selector_changed = any(selector_marker /= selector_fingerprint)" in source
    assert "call mpi_allreduce(mpi_in_place, history_selector_changed" in source
    changed = source.index(
        "if (history_selector_changed .or. invalid_selector_marker) then"
    )
    assert source.index("call flush_methane_acc_fluxes ()", changed) < source.index("return", changed)
    first_accumulator_read = source.index("'ch4_a_net_methane'", changed)
    assert source.index("return", changed) < first_accumulator_read
    assert "merge(257, 256, def_methane%write_ch4_history)" in source
    assert "merge(513, 512, def_methane%use_microbial_pools)" in source


def test_runtime_selector_change_also_resets_before_new_group_accumulation():
    source = _compact(ACC)
    changed = source.index("selector_changed = history_groups_initialized")
    reset = source.index("if (selector_changed) call flush_methane_acc_fluxes ()", changed)
    first_accumulation = source.index("if (need_main_layers) then", reset)
    assert changed < reset < first_accumulation


def test_legacy_restart_discards_incomplete_pre_schema3_history_window():
    source = _compact(ACC)
    assert "older committed schemas are readable" in source
    assert "window must be discarded" in source
    assert "if (has_history_selector_marker) then" in source
    reset = source.index("if (restart_schema_active < 3) then")
    assert source.index("call flush_methane_acc_fluxes ()", reset) > reset


def test_global_physical_flux_and_balance_include_lakes_with_sample_counts():
    source = _compact(HIST)
    for name in (
        "f_methane_surf_flux_global_phys_with_lake",
        "f_methane_balance_residual_global_with_lake",
        "f_methane_ch4_clip_credit_global_with_lake",
    ):
        assert name in source
    assert "hist_ch4_global_phys_with_lake" in source
    assert "hist_ch4_global_balance_with_lake" in source
    assert "hist_ch4_global_clip_credit_with_lake" in source
    assert "a_methane_surf_flux_tot_phys(ipatch) / a_methane_acc_num(ipatch)" in source
    assert "a_methane_balance_residual(ipatch) / a_methane_acc_num(ipatch)" in source
    assert "a_methane_surf_flux_tot_phys(ipatch) / a_methane_acc_num_lake(ipatch)" in source
    assert "a_methane_balance_residual(ipatch) / a_methane_acc_num_lake(ipatch)" in source
    assert "a_methane_ch4_clip_credit(ipatch) / a_methane_acc_num(ipatch)" in source
    assert "a_methane_ch4_clip_credit(ipatch) / a_methane_acc_num_lake(ipatch)" in source
    for array_name in (
        "hist_ch4_global_with_lake",
        "hist_ch4_global_phys_with_lake",
        "hist_ch4_global_balance_with_lake",
        "hist_ch4_global_clip_credit_with_lake",
    ):
        assert re.search(
            rf"call write_history_variable_2d\s*\(\.true\.,\s*&?\s*{array_name},.*?acc_num=hist_ch4_acc_one\)",
            source,
        )
