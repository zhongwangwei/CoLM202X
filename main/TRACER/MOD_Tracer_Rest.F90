#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Rest

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracer_init_water_ratio, trc_tiny, &
      trc_water_min_for_ratio, tracers, &
      tracer_uses_delta_diagnostics, tracer_uses_land_water_transport, tracer_is_nonvolatile_solute, &
      tracer_equilibrate_dissolved
   USE MOD_Tracer_Vars
   USE MOD_LandPatch, only: landpatch
   USE MOD_Block, only: get_filename_block
   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize
   USE MOD_NetCDFVector, only: ncio_read_vector, ncio_write_vector
   USE MOD_Tracer_Reactive, only: tracer_reactive_write_restart

   IMPLICIT NONE

   PRIVATE :: tracer_dim_matches

CONTAINS

   !-------------------------------------------------------------------
   ! Verify the fastest-varying extent of a restart variable matches the
   ! currently-compiled ntracers. The land restart writer stores each
   ! pool as a single (tracer, [soilsnow,] patch) variable, and
   ! ncio_read_vector → ncio_read_serial auto-reallocates the per-block
   ! scatter buffer to the on-disk shape. If the file was written with a
   ! different DEF_TRACER_NUM, the downstream mpi_scatterv still uses the
   ! current ntracers as its element multiplier and will misalign or
   ! read past the buffer. Detect the mismatch up front so the caller
   ! can fall back to tracer_init_from_water instead of crashing.
   !-------------------------------------------------------------------
   logical FUNCTION tracer_dim_matches (file_restart, varname, expect_soilsnow)
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_io, p_comm_io, p_comm_group, p_err, p_root, &
         MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
#else
      USE MOD_SPMD_Task, only: p_is_io
#endif
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      ! Optional second-dimension size (soilsnow extent = nl_soil-maxsnl).
      ! When present, tracer_dim_matches additionally requires
      ! varsize(2) == expect_soilsnow; this guards against a restart
      ! written with different maxsnl / nl_soil values, which would
      ! otherwise pass the ntracers-only check and then silently misalign
      ! per-layer data via ncio_read_vector's reshape.
      integer, intent(in), optional :: expect_soilsnow
      integer, allocatable :: varsize(:)
      integer :: iblkgrp, iblk, jblk, expected_rank
      integer :: counts(3)
      character(len=256) :: fileblock
      logical :: block_shape_ok

      counts(:) = 0
      expected_rank = merge(3, 2, present(expect_soilsnow))

      ! Vector restart files are split by block via get_filename_block().
      ! Checking the unsuffixed base filename makes hot starts look like
      ! old/missing tracer restarts and silently reinitializes NSS state.
      IF (p_is_io) THEN
         counts(1) = landpatch%nblkgrp
         DO iblkgrp = 1, landpatch%nblkgrp
            iblk = landpatch%xblkgrp(iblkgrp)
            jblk = landpatch%yblkgrp(iblkgrp)
            CALL get_filename_block(file_restart, iblk, jblk, fileblock)

            IF (.not. ncio_var_exist(fileblock, varname, readflag = .false.)) THEN
               CYCLE
            ENDIF
            counts(2) = counts(2) + 1

            CALL ncio_inquire_varsize(fileblock, varname, varsize)
            block_shape_ok = .false.
            IF (allocated(varsize)) THEN
               IF (size(varsize) == expected_rank) THEN
                  block_shape_ok = varsize(1) == ntracers .and. &
                     varsize(expected_rank) == landpatch%vecgs%vlen(iblk,jblk)
                  IF (block_shape_ok .and. present(expect_soilsnow)) THEN
                     block_shape_ok = varsize(2) == expect_soilsnow
                  ENDIF
               ENDIF
            ENDIF
            IF (allocated(varsize)) deallocate(varsize)
            IF (.not. block_shape_ok) counts(3) = counts(3) + 1
         ENDDO
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL mpi_allreduce(MPI_IN_PLACE, counts, 3, MPI_INTEGER, MPI_SUM, p_comm_io, p_err)
      ENDIF
      CALL mpi_bcast(counts, 3, MPI_INTEGER, p_root, p_comm_group, p_err)
#endif

      tracer_dim_matches = counts(1) > 0 .and. counts(2) == counts(1) .and. counts(3) == 0
   END FUNCTION tracer_dim_matches

   SUBROUTINE tracer_init_from_water (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, waterstorage)

      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      ! Optional: per-patch irrigation reservoir from
      ! MOD_Vars_TimeVariables. Present when DEF_USE_IRRIGATION is on
      ! under CROP. When passed, trc_waterstorage is cold-started to
      ! waterstorage*R_init so the first tracer_save_storage sees a
      ! correct starting inventory; without this the pool would be 0
      ! and the first step's balance check would under-count the
      ! reservoir mass that the water side already tracks.
      real(r8), intent(in), optional :: waterstorage(numpatch)

      integer  :: itrc, ip, j, snl_local
      real(r8) :: R_init

      IF (allocated(trc_surface_residue)) trc_surface_residue = 0._r8
      IF (allocated(trc_subsurface_residue)) trc_subsurface_residue = 0._r8
      IF (allocated(trc_solid_soisno)) trc_solid_soisno = 0._r8
      IF (allocated(trc_canopy_solid)) trc_canopy_solid = 0._r8
      IF (allocated(trc_surface_solid)) trc_surface_solid = 0._r8
      IF (allocated(trc_subsurface_solid)) trc_subsurface_solid = 0._r8
      IF (allocated(trc_waterstorage_solid)) trc_waterstorage_solid = 0._r8
      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         R_init = tracer_init_water_ratio(itrc)
         DO ip = 1, numpatch
            trc_ldew_rain(itrc, ip) = ldew_rain(ip) * R_init
            trc_ldew_snow(itrc, ip) = ldew_snow(ip) * R_init
            CALL tracer_equilibrate_dissolved(itrc, max(ldew_rain(ip), 0._r8), &
               trc_ldew_rain(itrc, ip), trc_canopy_solid(itrc, ip))
            DO j = maxsnl + 1, nl_soil
               trc_wliq_soisno(itrc, j, ip) = max(wliq_soisno(j, ip), 0._r8) * R_init
               trc_wice_soisno(itrc, j, ip) = max(wice_soisno(j, ip), 0._r8) * R_init
               CALL tracer_equilibrate_dissolved(itrc, max(wliq_soisno(j, ip), 0._r8), &
                  trc_wliq_soisno(itrc, j, ip), trc_solid_soisno(itrc, j, ip))
            ENDDO
            ! wa uses the SIGNED water value so a hydrology restart with an
            ! aquifer-debt state (wa<0, recorded by the wetland branch at
            ! MOD_SoilSnowHydrology.F90:1200-1203) starts with a matching
            ! trc_wa<0. Without this the first wetland recovery step after a
            ! cold tracer-start would see pool_water = wa_bef + inputs but
            ! pool_tracer = 0 + inputs*R_atm (trc_wa_bef was clamped to 0),
            ! producing a 2× over-concentration in the recovered wetwat.
            ! wdsrf/wetwat stay with max(.,0) since WATER_VSF never leaves
            ! them negative by construction (overflow case L1196-1207).
            trc_wa    (itrc, ip) = wa(ip) * R_init
            trc_wdsrf (itrc, ip) = max(wdsrf(ip),  0._r8) * R_init
            trc_wetwat(itrc, ip) = max(wetwat(ip), 0._r8) * R_init
            CALL tracer_equilibrate_dissolved(itrc, max(wdsrf(ip), 0._r8), &
               trc_wdsrf(itrc, ip), trc_surface_solid(itrc, ip))
            CALL tracer_equilibrate_dissolved(itrc, max(wetwat(ip), 0._r8), &
               trc_wetwat(itrc, ip), trc_surface_solid(itrc, ip))
            CALL tracer_equilibrate_dissolved(itrc, wa(ip), &
               trc_wa(itrc, ip), trc_subsurface_solid(itrc, ip))
            ! Reconstruct snl from the snow layer water content (same
            ! recipe as CoLMMAIN L770-774): snl counts non-empty snow
            ! layers from the top. snl is a runtime scalar, not a
            ! persistent per-patch array, so we infer it here.
            ! Relies on CoLM's snow-column contiguity convention: active
            ! snow layers occupy indices [snl+1 : 0] with no gaps, so
            ! iterating from j=0 downward and EXIT-ing at the first empty
            ! slot correctly recovers snl. If a future snow path violates
            ! this (gaps in the column), replace EXIT with a count of
            ! non-empty slots instead.
            snl_local = 0
            DO j = 0, maxsnl + 1, -1
               IF (wliq_soisno(j, ip) + wice_soisno(j, ip) > 0._r8) THEN
                  snl_local = snl_local - 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            ! trc_scv holds the pre-layer snow tracer pool. It is only
            ! populated when snl==0 (thin snow, no layer yet); once a
            ! snow layer is created, tracer lives in trc_wice/wliq and
            ! trc_scv stays zero.
            IF (snl_local == 0) THEN
               trc_scv(itrc, ip) = max(scv(ip), 0._r8) * R_init
            ELSE
               trc_scv(itrc, ip) = 0._r8
            ENDIF
            IF (present(waterstorage) .and. allocated(trc_waterstorage)) THEN
               trc_waterstorage(itrc, ip) = max(waterstorage(ip), 0._r8) * R_init
               CALL tracer_equilibrate_dissolved(itrc, max(waterstorage(ip), 0._r8), &
                  trc_waterstorage(itrc, ip), trc_waterstorage_solid(itrc, ip))
            ENDIF
         ENDDO
      ENDDO
      CALL zero_particle_land_tracer_state()
   END SUBROUTINE tracer_init_from_water

   SUBROUTINE read_land_tracer_restart (file_restart, maxsnl, nl_soil, found_restart, &
      scv_missing, waterstorage_missing)
      USE MOD_SPMD_Task, only: p_is_io, p_is_worker, p_iam_io, p_root
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil
      logical, intent(out) :: found_restart
      ! .true. iff the restart file is otherwise complete but the trc_scv
      ! field is absent (old-format restart). The caller uses this to
      ! rebuild trc_scv from the concurrent scv water state instead of
      ! leaving it at the post-allocation zero — without this, a hot
      ! start with non-zero scv but no layered snow silently loses the
      ! thin-snow tracer once it melts.
      logical, optional, intent(out) :: scv_missing
      ! .true. iff the restart is otherwise usable but lacks the optional
      ! irrigation-reservoir tracer pool.
      logical, optional, intent(out) :: waterstorage_missing
      integer :: itrc
      logical :: has_leaf_delta_e, has_leaf_delta_b
      logical :: has_leaf_peclet, has_leaf_water_moles, has_leaf_iso_storage
      logical :: reset_legacy_leaf_e, reset_legacy_leaf_b
      logical :: has_nonvolatile_solute
      logical :: restart_complete, field_matches

      found_restart = .false.
      IF (present(scv_missing)) scv_missing = .false.
      IF (present(waterstorage_missing)) waterstorage_missing = .false.
      IF (ntracers <= 0) RETURN
      has_nonvolatile_solute = .false.
      DO itrc = 1, ntracers
         has_nonvolatile_solute = has_nonvolatile_solute .or. tracer_is_nonvolatile_solute(itrc)
      ENDDO

      ! The master/control rank does not belong to the landpatch vector
      ! IO/worker group under MPI. Its landpatch%nblkgrp is therefore 0,
      ! so probing restart files there would falsely report a cold start.
      IF (.not. (p_is_io .or. p_is_worker)) THEN
         found_restart = .true.
         RETURN
      ENDIF

      ! Reject the restart if required tracer variables are absent from
      ! the block-split vector files, if ntracers changed, or if the
      ! per-layer extent no longer matches nl_soil-maxsnl.
      restart_complete = .true.
      field_matches = tracer_dim_matches(file_restart, 'trc_ldew_rain')
      restart_complete = restart_complete .and. field_matches
      field_matches = tracer_dim_matches(file_restart, 'trc_ldew_snow')
      restart_complete = restart_complete .and. field_matches
      field_matches = tracer_dim_matches(file_restart, 'trc_wliq_soisno', &
                                         expect_soilsnow = nl_soil - maxsnl)
      restart_complete = restart_complete .and. field_matches
      field_matches = tracer_dim_matches(file_restart, 'trc_wice_soisno', &
                                         expect_soilsnow = nl_soil - maxsnl)
      restart_complete = restart_complete .and. field_matches
      field_matches = tracer_dim_matches(file_restart, 'trc_wa')
      restart_complete = restart_complete .and. field_matches
      field_matches = tracer_dim_matches(file_restart, 'trc_wdsrf')
      restart_complete = restart_complete .and. field_matches
      field_matches = tracer_dim_matches(file_restart, 'trc_wetwat')
      restart_complete = restart_complete .and. field_matches
      IF (.not. restart_complete) THEN
         IF (p_is_io .and. p_iam_io == p_root) WRITE(*,*) &
            'Tracer restart ntracers/soilsnow mismatch in ', &
            TRIM(file_restart), '. Using cold-start.'
         RETURN
      ENDIF

      CALL ncio_read_vector(file_restart, 'trc_ldew_rain', ntracers, landpatch, trc_ldew_rain)
      CALL ncio_read_vector(file_restart, 'trc_ldew_snow', ntracers, landpatch, trc_ldew_snow)
      CALL ncio_read_vector(file_restart, 'trc_wliq_soisno', ntracers, nl_soil-maxsnl, landpatch, trc_wliq_soisno)
      CALL ncio_read_vector(file_restart, 'trc_wice_soisno', ntracers, nl_soil-maxsnl, landpatch, trc_wice_soisno)
      CALL ncio_read_vector(file_restart, 'trc_wa', ntracers, landpatch, trc_wa)
      CALL ncio_read_vector(file_restart, 'trc_wdsrf', ntracers, landpatch, trc_wdsrf)
      CALL ncio_read_vector(file_restart, 'trc_wetwat', ntracers, landpatch, trc_wetwat)
      IF (allocated(trc_surface_residue)) trc_surface_residue = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_surface_residue')) THEN
         CALL ncio_read_vector(file_restart, 'trc_surface_residue', ntracers, landpatch, trc_surface_residue)
      ELSE
         IF (has_nonvolatile_solute .and. p_is_io .and. p_iam_io == p_root) WRITE(*,*) &
            'Tracer restart has no trc_surface_residue; initializing the optional pool to zero.'
      ENDIF
      IF (allocated(trc_subsurface_residue)) trc_subsurface_residue = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_subsurface_residue')) THEN
         CALL ncio_read_vector(file_restart, 'trc_subsurface_residue', ntracers, landpatch, trc_subsurface_residue)
      ELSE
         IF (has_nonvolatile_solute .and. p_is_io .and. p_iam_io == p_root) WRITE(*,*) &
            'Tracer restart has no trc_subsurface_residue; initializing the optional pool to zero.'
      ENDIF
      IF (allocated(trc_solid_soisno)) trc_solid_soisno = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_solid_soisno', &
                             expect_soilsnow = nl_soil - maxsnl)) THEN
         CALL ncio_read_vector(file_restart, 'trc_solid_soisno', ntracers, nl_soil-maxsnl, &
            landpatch, trc_solid_soisno)
      ENDIF
      IF (allocated(trc_canopy_solid)) trc_canopy_solid = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_canopy_solid')) THEN
         CALL ncio_read_vector(file_restart, 'trc_canopy_solid', ntracers, landpatch, trc_canopy_solid)
      ENDIF
      IF (allocated(trc_surface_solid)) trc_surface_solid = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_surface_solid')) THEN
         CALL ncio_read_vector(file_restart, 'trc_surface_solid', ntracers, landpatch, trc_surface_solid)
      ENDIF
      IF (allocated(trc_subsurface_solid)) trc_subsurface_solid = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_subsurface_solid')) THEN
         CALL ncio_read_vector(file_restart, 'trc_subsurface_solid', ntracers, landpatch, trc_subsurface_solid)
      ENDIF
      IF (allocated(trc_waterstorage_solid)) trc_waterstorage_solid = 0._r8
      IF (tracer_dim_matches(file_restart, 'trc_waterstorage_solid')) THEN
         CALL ncio_read_vector(file_restart, 'trc_waterstorage_solid', ntracers, landpatch, trc_waterstorage_solid)
      ENDIF
      IF (tracer_dim_matches(file_restart, 'trc_scv')) THEN
         CALL ncio_read_vector(file_restart, 'trc_scv', ntracers, landpatch, trc_scv)
      ELSE
         IF (present(scv_missing)) scv_missing = .true.
         IF (p_is_io .and. p_iam_io == p_root) WRITE(*,*) &
            'Tracer restart has no trc_scv; will rebuild from scv state.'
      ENDIF
      ! trc_waterstorage is optional in restart for backward compat
      ! with files written before the irrigation bookkeeping fix.
      ! The ncio_read_vector collective must be entered consistently by
      ! IO and worker ranks. Guarding it with ALLOCATED(...) is unsafe
      ! here because IO/master can legitimately have numpatch<=0 while
      ! workers hold the actual patch array, which would make only a
      ! subset of ranks enter the MPI scatter/gather path. Follow the
      ! same rank-uniform pattern as the other tracer restart pools:
      ! when the variable exists and dimensions match, everybody calls
      ! ncio_read_vector; when it is absent, let the caller rebuild it
      ! from the loaded patch tracer ratio.
      IF (tracer_dim_matches(file_restart, 'trc_waterstorage')) THEN
         CALL ncio_read_vector(file_restart, 'trc_waterstorage', ntracers, landpatch, trc_waterstorage)
      ELSE
         IF (present(waterstorage_missing)) waterstorage_missing = .true.
         IF (p_is_io .and. p_iam_io == p_root) WRITE(*,*) &
            'Tracer restart has no compatible trc_waterstorage; will rebuild from patch tracer ratio.'
      ENDIF

      IF (allocated(trc_leaf_delta_e)) THEN
         DO itrc = 1, ntracers
            IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
            IF (tracer_uses_delta_diagnostics(itrc)) THEN
               trc_leaf_delta_e(itrc, :) = tracers(itrc)%init_delta
            ELSE
               trc_leaf_delta_e(itrc, :) = 0._r8
            ENDIF
         ENDDO
      ENDIF
      IF (allocated(trc_leaf_delta_b)) THEN
         DO itrc = 1, ntracers
            IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
            IF (tracer_uses_delta_diagnostics(itrc)) THEN
               trc_leaf_delta_b(itrc, :) = tracers(itrc)%init_delta
            ELSE
               trc_leaf_delta_b(itrc, :) = 0._r8
            ENDIF
         ENDDO
      ENDIF
      IF (allocated(trc_leaf_peclet)) trc_leaf_peclet = 1._r8
      IF (allocated(trc_leaf_water_moles)) trc_leaf_water_moles = 0._r8
      IF (allocated(trc_leaf_iso_storage)) trc_leaf_iso_storage = 0._r8
      has_leaf_delta_e = tracer_dim_matches(file_restart, 'trc_leaf_delta_e')
      has_leaf_delta_b = tracer_dim_matches(file_restart, 'trc_leaf_delta_b')
      has_leaf_peclet = tracer_dim_matches(file_restart, 'trc_leaf_peclet')
      has_leaf_water_moles = tracer_dim_matches(file_restart, 'trc_leaf_water_moles')
      has_leaf_iso_storage = tracer_dim_matches(file_restart, 'trc_leaf_iso_storage')
      IF (has_leaf_delta_e) THEN
         CALL ncio_read_vector(file_restart, 'trc_leaf_delta_e', ntracers, landpatch, trc_leaf_delta_e)
      ENDIF
      IF (has_leaf_delta_b) THEN
         CALL ncio_read_vector(file_restart, 'trc_leaf_delta_b', ntracers, landpatch, trc_leaf_delta_b)
      ENDIF
      IF (has_leaf_peclet) THEN
         CALL ncio_read_vector(file_restart, 'trc_leaf_peclet', ntracers, landpatch, trc_leaf_peclet)
      ENDIF
      IF (has_leaf_water_moles) THEN
         CALL ncio_read_vector(file_restart, 'trc_leaf_water_moles', ntracers, landpatch, trc_leaf_water_moles)
      ENDIF
      IF (has_leaf_iso_storage) THEN
         CALL ncio_read_vector(file_restart, 'trc_leaf_iso_storage', ntracers, landpatch, trc_leaf_iso_storage)
      ENDIF
      IF (p_is_io .and. p_iam_io == p_root .and. &
          (.not. has_leaf_delta_e .or. .not. has_leaf_delta_b .or. &
           .not. has_leaf_peclet .or. .not. has_leaf_water_moles .or. &
           .not. has_leaf_iso_storage)) THEN
         WRITE(*,'(A)') 'Tracer restart lacks one or more leaf NSS diagnostic fields; initialized defaults.'
      ENDIF

      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         IF (.not. tracer_uses_delta_diagnostics(itrc)) THEN
            IF (allocated(trc_leaf_delta_e)) trc_leaf_delta_e(itrc, :) = 0._r8
            IF (allocated(trc_leaf_delta_b)) trc_leaf_delta_b(itrc, :) = 0._r8
            IF (allocated(trc_leaf_peclet)) trc_leaf_peclet(itrc, :) = 1._r8
            IF (allocated(trc_leaf_water_moles)) trc_leaf_water_moles(itrc, :) = 0._r8
            IF (allocated(trc_leaf_iso_storage)) trc_leaf_iso_storage(itrc, :) = 0._r8
            CYCLE
         ENDIF
         reset_legacy_leaf_e = allocated(trc_leaf_delta_e) .and. &
            abs(tracers(itrc)%init_delta) > trc_tiny .and. &
            all(abs(trc_leaf_delta_e(itrc, :)) <= trc_tiny)
         reset_legacy_leaf_b = allocated(trc_leaf_delta_b) .and. &
            abs(tracers(itrc)%init_delta) > trc_tiny .and. &
            all(abs(trc_leaf_delta_b(itrc, :)) <= trc_tiny)
         IF (reset_legacy_leaf_e) trc_leaf_delta_e(itrc, :) = tracers(itrc)%init_delta
         IF (reset_legacy_leaf_b) trc_leaf_delta_b(itrc, :) = tracers(itrc)%init_delta
         IF ((reset_legacy_leaf_e .or. reset_legacy_leaf_b) .and. &
             p_is_io .and. p_iam_io == p_root) THEN
            WRITE(*,'(A,I0,A,F10.3)') &
               'Tracer restart legacy all-zero leaf NSS delta reset for tracer ', &
               itrc, ' to init_delta=', tracers(itrc)%init_delta
         ENDIF
      ENDDO
      CALL zero_particle_land_tracer_state()
      found_restart = .true.
   END SUBROUTINE read_land_tracer_restart

   !-------------------------------------------------------------------
   ! Rebuild only trc_scv from the current scv state, leaving every
   ! other tracer pool as-is. Used after read_land_tracer_restart when
   ! the restart file was otherwise complete but lacked trc_scv (old
   ! format). Uses the same snl-from-water heuristic as
   ! tracer_init_from_water to decide when trc_scv is meaningful.
   !-------------------------------------------------------------------
   SUBROUTINE tracer_init_scv_from_water (numpatch, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno, scv)
      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: scv(numpatch)

      integer  :: itrc, ip, j, snl_local
      real(r8) :: R_init

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         R_init = tracer_init_water_ratio(itrc)
         DO ip = 1, numpatch
            snl_local = 0
            DO j = 0, maxsnl + 1, -1
               IF (wliq_soisno(j, ip) + wice_soisno(j, ip) > 0._r8) THEN
                  snl_local = snl_local - 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            IF (snl_local == 0) THEN
               trc_scv(itrc, ip) = max(scv(ip), 0._r8) * R_init
            ELSE
               trc_scv(itrc, ip) = 0._r8
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE tracer_init_scv_from_water

   SUBROUTINE tracer_init_waterstorage_from_ratio (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, waterstorage)
      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      real(r8), intent(in) :: waterstorage(numpatch)

      integer  :: itrc, ip, j, snl_local
      real(r8) :: R_init, ratio, water_ref, tracer_ref

      IF (ntracers <= 0) RETURN
      IF (.not. allocated(trc_waterstorage)) RETURN

      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         R_init = tracer_init_water_ratio(itrc)
         DO ip = 1, numpatch
            snl_local = 0
            DO j = 0, maxsnl + 1, -1
               IF (wliq_soisno(j, ip) + wice_soisno(j, ip) > 0._r8) THEN
                  snl_local = snl_local - 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
            water_ref = max(ldew_rain(ip), 0._r8) + max(ldew_snow(ip), 0._r8) &
                      + max(wa(ip), 0._r8) + max(wdsrf(ip), 0._r8) &
                      + max(wetwat(ip), 0._r8) + max(scv(ip), 0._r8)
            tracer_ref = max(trc_ldew_rain(itrc, ip), 0._r8) &
                       + max(trc_ldew_snow(itrc, ip), 0._r8) &
                       + max(trc_wa(itrc, ip), 0._r8) &
                       + max(trc_wdsrf(itrc, ip), 0._r8) &
                       + max(trc_wetwat(itrc, ip), 0._r8) &
                       + max(trc_scv(itrc, ip), 0._r8)
            DO j = snl_local + 1, nl_soil
               water_ref = water_ref + max(wliq_soisno(j, ip), 0._r8)
               tracer_ref = tracer_ref + max(trc_wliq_soisno(itrc, j, ip), 0._r8)
               water_ref = water_ref + max(wice_soisno(j, ip), 0._r8)
               tracer_ref = tracer_ref + max(trc_wice_soisno(itrc, j, ip), 0._r8)
            ENDDO
            IF (water_ref > trc_water_min_for_ratio) THEN
               ratio = tracer_ref / water_ref
            ELSE
               ratio = R_init
            ENDIF
            trc_waterstorage(itrc, ip) = max(waterstorage(ip), 0._r8) * ratio
            CALL tracer_equilibrate_dissolved(itrc, max(waterstorage(ip), 0._r8), &
               trc_waterstorage(itrc, ip), trc_waterstorage_solid(itrc, ip))
         ENDDO
      ENDDO
   END SUBROUTINE tracer_init_waterstorage_from_ratio

   SUBROUTINE tracer_enforce_solubility_from_water (numpatch, maxsnl, nl_soil, &
      ldew_rain, wliq_soisno, wa, wdsrf, wetwat, waterstorage)
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch), wdsrf(numpatch), wetwat(numpatch)
      real(r8), intent(in), optional :: waterstorage(numpatch)

      integer :: itrc, ip, j
      real(r8) :: surface_water, surface_tracer

      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         DO ip = 1, numpatch
            CALL tracer_equilibrate_dissolved(itrc, max(ldew_rain(ip), 0._r8), &
               trc_ldew_rain(itrc, ip), trc_canopy_solid(itrc, ip))
            DO j = maxsnl + 1, nl_soil
               CALL tracer_equilibrate_dissolved(itrc, max(wliq_soisno(j, ip), 0._r8), &
                  trc_wliq_soisno(itrc, j, ip), trc_solid_soisno(itrc, j, ip))
            ENDDO
            surface_water = max(wdsrf(ip), 0._r8) + max(wetwat(ip), 0._r8)
            surface_tracer = trc_wdsrf(itrc, ip) + trc_wetwat(itrc, ip)
            CALL tracer_equilibrate_dissolved(itrc, surface_water, surface_tracer, &
               trc_surface_solid(itrc, ip))
            IF (surface_water > trc_water_min_for_ratio) THEN
               trc_wdsrf(itrc, ip) = surface_tracer * max(wdsrf(ip), 0._r8) / surface_water
               trc_wetwat(itrc, ip) = surface_tracer * max(wetwat(ip), 0._r8) / surface_water
            ELSE
               trc_wdsrf(itrc, ip) = 0._r8
               trc_wetwat(itrc, ip) = 0._r8
            ENDIF
            CALL tracer_equilibrate_dissolved(itrc, wa(ip), trc_wa(itrc, ip), &
               trc_subsurface_solid(itrc, ip))
            IF (present(waterstorage) .and. allocated(trc_waterstorage)) THEN
               CALL tracer_equilibrate_dissolved(itrc, max(waterstorage(ip), 0._r8), &
                  trc_waterstorage(itrc, ip), trc_waterstorage_solid(itrc, ip))
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE tracer_enforce_solubility_from_water

   SUBROUTINE write_land_tracer_restart (file_restart, maxsnl, nl_soil, numpatch, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, waterstorage)
      USE MOD_SPMD_Task, only: p_is_worker
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil, numpatch
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      real(r8), intent(in), optional :: waterstorage(numpatch)

      integer :: itrc, j
      real(r8) :: R_init
      real(r8), allocatable :: restart_patch(:,:)
      real(r8), allocatable :: restart_soilsnow(:,:,:)
      logical :: have_patch_data

      IF (ntracers <= 0) RETURN

      CALL zero_particle_land_tracer_state()

      have_patch_data = p_is_worker .and. numpatch > 0

      allocate(restart_patch(ntracers, numpatch))

      restart_patch(:, :) = 0._r8
      IF (have_patch_data) restart_patch(:, :) = trc_ldew_rain(:, :)
      CALL ncio_write_vector(file_restart, 'trc_ldew_rain', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (have_patch_data) restart_patch(:, :) = trc_ldew_snow(:, :)
      CALL ncio_write_vector(file_restart, 'trc_ldew_snow', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      allocate(restart_soilsnow(ntracers, maxsnl+1:nl_soil, numpatch))

      restart_soilsnow(:, :, :) = 0._r8
      IF (have_patch_data) restart_soilsnow(:, :, :) = trc_wliq_soisno(:, maxsnl+1:nl_soil, :)
      CALL ncio_write_vector(file_restart, 'trc_wliq_soisno', 'tracer', ntracers, 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, &
         restart_soilsnow, DEF_REST_CompressLevel)

      restart_soilsnow(:, :, :) = 0._r8
      IF (have_patch_data) restart_soilsnow(:, :, :) = trc_wice_soisno(:, maxsnl+1:nl_soil, :)
      CALL ncio_write_vector(file_restart, 'trc_wice_soisno', 'tracer', ntracers, 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, &
         restart_soilsnow, DEF_REST_CompressLevel)

      restart_soilsnow(:, :, :) = 0._r8
      IF (allocated(trc_solid_soisno) .and. have_patch_data) &
         restart_soilsnow(:, :, :) = trc_solid_soisno(:, maxsnl+1:nl_soil, :)
      CALL ncio_write_vector(file_restart, 'trc_solid_soisno', 'tracer', ntracers, &
         'soilsnow', nl_soil-maxsnl, 'patch', landpatch, restart_soilsnow, DEF_REST_CompressLevel)

      deallocate(restart_soilsnow)

      restart_patch(:, :) = 0._r8
      IF (have_patch_data) restart_patch(:, :) = trc_wa(:, :)
      CALL ncio_write_vector(file_restart, 'trc_wa', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (have_patch_data) restart_patch(:, :) = trc_wdsrf(:, :)
      CALL ncio_write_vector(file_restart, 'trc_wdsrf', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (have_patch_data) restart_patch(:, :) = trc_wetwat(:, :)
      CALL ncio_write_vector(file_restart, 'trc_wetwat', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_surface_residue) .and. have_patch_data) &
         restart_patch(:, :) = trc_surface_residue(:, :)
      CALL ncio_write_vector(file_restart, 'trc_surface_residue', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_subsurface_residue) .and. have_patch_data) &
         restart_patch(:, :) = trc_subsurface_residue(:, :)
      CALL ncio_write_vector(file_restart, 'trc_subsurface_residue', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_canopy_solid) .and. have_patch_data) restart_patch(:, :) = trc_canopy_solid(:, :)
      CALL ncio_write_vector(file_restart, 'trc_canopy_solid', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_surface_solid) .and. have_patch_data) restart_patch(:, :) = trc_surface_solid(:, :)
      CALL ncio_write_vector(file_restart, 'trc_surface_solid', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_subsurface_solid) .and. have_patch_data) restart_patch(:, :) = trc_subsurface_solid(:, :)
      CALL ncio_write_vector(file_restart, 'trc_subsurface_solid', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_waterstorage_solid) .and. have_patch_data) restart_patch(:, :) = trc_waterstorage_solid(:, :)
      CALL ncio_write_vector(file_restart, 'trc_waterstorage_solid', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (have_patch_data) restart_patch(:, :) = trc_scv(:, :)
      CALL ncio_write_vector(file_restart, 'trc_scv', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)
      ! Persist the irrigation-reservoir tracer. This write must be
      ! entered by all ranks, not just the workers that currently own a
      ! non-empty patch vector: otherwise workers can advance one
      ! collective further than the IO root and the next restart write
      ! (e.g. PFT thermk_p) will see a patch-count/PFT-count mismatch.
      ! The pool itself is valid even when irrigation is off; it stays
      ! at 0 and readers that do not use it simply ignore the variable.
      restart_patch(:, :) = 0._r8
      IF (allocated(trc_waterstorage) .and. have_patch_data) restart_patch(:, :) = trc_waterstorage(:, :)
      CALL ncio_write_vector(file_restart, 'trc_waterstorage', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_leaf_delta_e) .and. have_patch_data) restart_patch(:, :) = trc_leaf_delta_e(:, :)
      CALL ncio_write_vector(file_restart, 'trc_leaf_delta_e', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_leaf_delta_b) .and. have_patch_data) restart_patch(:, :) = trc_leaf_delta_b(:, :)
      CALL ncio_write_vector(file_restart, 'trc_leaf_delta_b', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_leaf_peclet) .and. have_patch_data) restart_patch(:, :) = trc_leaf_peclet(:, :)
      CALL ncio_write_vector(file_restart, 'trc_leaf_peclet', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_leaf_water_moles) .and. have_patch_data) restart_patch(:, :) = trc_leaf_water_moles(:, :)
      CALL ncio_write_vector(file_restart, 'trc_leaf_water_moles', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      restart_patch(:, :) = 0._r8
      IF (allocated(trc_leaf_iso_storage) .and. have_patch_data) restart_patch(:, :) = trc_leaf_iso_storage(:, :)
      CALL ncio_write_vector(file_restart, 'trc_leaf_iso_storage', 'tracer', ntracers, 'patch', landpatch, &
         restart_patch, DEF_REST_CompressLevel)

      deallocate(restart_patch)
   END SUBROUTINE write_land_tracer_restart

   SUBROUTINE write_tracer_restart_all (file_restart, maxsnl, nl_soil, numpatch, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, &
      compress, waterstorage)
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil, numpatch, compress
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)
      real(r8), intent(in) :: scv(numpatch)
      real(r8), intent(in), optional :: waterstorage(numpatch)

      IF (present(waterstorage)) THEN
         CALL write_land_tracer_restart(file_restart, maxsnl, nl_soil, numpatch, &
            ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv, waterstorage)
      ELSE
         CALL write_land_tracer_restart(file_restart, maxsnl, nl_soil, numpatch, &
            ldew_rain, ldew_snow, wliq_soisno, wice_soisno, wa, wdsrf, wetwat, scv)
      ENDIF

      CALL tracer_reactive_write_restart(file_restart, compress)

   END SUBROUTINE write_tracer_restart_all

END MODULE MOD_Tracer_Rest
#endif
