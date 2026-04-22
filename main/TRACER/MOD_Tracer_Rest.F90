#include <define.h>

MODULE MOD_Tracer_Rest

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, delta_to_R, trc_tiny
   USE MOD_Tracer_Vars
   USE MOD_LandPatch, only: landpatch
   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_var_exist, ncio_inquire_varsize
   USE MOD_NetCDFVector, only: ncio_read_vector, ncio_write_vector

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

      tracer_dim_matches = .false.
      CALL ncio_inquire_varsize(file_restart, varname, varsize)
      IF (allocated(varsize)) THEN
         IF (size(varsize) >= 1) THEN
            tracer_dim_matches = (varsize(1) == ntracers)
            IF (tracer_dim_matches .and. present(expect_soilsnow)) THEN
               IF (size(varsize) >= 2) THEN
                  tracer_dim_matches = (varsize(2) == expect_soilsnow)
               ELSE
                  tracer_dim_matches = .false.
               ENDIF
            ENDIF
         ENDIF
         deallocate(varsize)
      ENDIF
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

      DO itrc = 1, ntracers
         R_init = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
         DO ip = 1, numpatch
            trc_ldew_rain(itrc, ip) = ldew_rain(ip) * R_init
            trc_ldew_snow(itrc, ip) = ldew_snow(ip) * R_init
            DO j = maxsnl + 1, nl_soil
               trc_wliq_soisno(itrc, j, ip) = max(wliq_soisno(j, ip), 0._r8) * R_init
               trc_wice_soisno(itrc, j, ip) = max(wice_soisno(j, ip), 0._r8) * R_init
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
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE tracer_init_from_water

   SUBROUTINE read_land_tracer_restart (file_restart, maxsnl, nl_soil, found_restart, scv_missing)
      USE MOD_SPMD_Task, only: p_is_master
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
      logical :: file_ok

      found_restart = .false.
      IF (present(scv_missing)) scv_missing = .false.
      IF (ntracers <= 0) RETURN

      inquire(file = trim(file_restart), exist = file_ok)
      IF (.not. file_ok) RETURN

      ! Check all required variables exist before reading
      IF (.not. ncio_var_exist(file_restart, 'trc_ldew_rain',   readflag = .false.) .or. &
          .not. ncio_var_exist(file_restart, 'trc_ldew_snow',   readflag = .false.) .or. &
          .not. ncio_var_exist(file_restart, 'trc_wliq_soisno', readflag = .false.) .or. &
          .not. ncio_var_exist(file_restart, 'trc_wice_soisno', readflag = .false.) .or. &
          .not. ncio_var_exist(file_restart, 'trc_wa',          readflag = .false.) .or. &
          .not. ncio_var_exist(file_restart, 'trc_wdsrf',       readflag = .false.) .or. &
          .not. ncio_var_exist(file_restart, 'trc_wetwat',      readflag = .false.)) THEN
         IF (p_is_master) WRITE(*,*) 'Tracer restart variables not found in ', &
            TRIM(file_restart), '. Using cold-start.'
         RETURN
      ENDIF

      ! Reject the restart if any variable's tracer dimension no longer
      ! matches ntracers (DEF_TRACER_NUM changed between runs), or the
      ! per-layer extent no longer matches nl_soil-maxsnl (maxsnl or
      ! nl_soil changed). Without these guards ncio_read_vector's
      ! scatter misaligns element-wise.
      IF (.not. tracer_dim_matches(file_restart, 'trc_ldew_rain'  ) .or. &
          .not. tracer_dim_matches(file_restart, 'trc_ldew_snow'  ) .or. &
          .not. tracer_dim_matches(file_restart, 'trc_wliq_soisno', &
                                   expect_soilsnow = nl_soil - maxsnl) .or. &
          .not. tracer_dim_matches(file_restart, 'trc_wice_soisno', &
                                   expect_soilsnow = nl_soil - maxsnl) .or. &
          .not. tracer_dim_matches(file_restart, 'trc_wa'         ) .or. &
          .not. tracer_dim_matches(file_restart, 'trc_wdsrf'      ) .or. &
          .not. tracer_dim_matches(file_restart, 'trc_wetwat'     )) THEN
         IF (p_is_master) WRITE(*,*) 'Tracer restart ntracers/soilsnow mismatch in ', &
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
      IF (ncio_var_exist(file_restart, 'trc_scv', readflag = .false.)) THEN
         IF (tracer_dim_matches(file_restart, 'trc_scv')) THEN
            CALL ncio_read_vector(file_restart, 'trc_scv', ntracers, landpatch, trc_scv)
         ELSE
            IF (present(scv_missing)) scv_missing = .true.
            IF (p_is_master) WRITE(*,*) 'Tracer restart trc_scv ntracers mismatch; will rebuild from scv state.'
         ENDIF
      ELSE
         IF (present(scv_missing)) scv_missing = .true.
         IF (p_is_master) WRITE(*,*) 'Tracer restart has no trc_scv; will rebuild from scv state.'
      ENDIF
      ! trc_waterstorage is optional in restart for backward compat
      ! with files written before the irrigation bookkeeping fix.
      ! When missing, leave at 0: the first tracer_save_storage call
      ! will re-sync it from hydrology's waterstorage under the
      ! Phase-1 invariant. When present, read it so Phase-2 (future
      ! time-varying refill R) round-trips the reservoir ratio.
      IF (allocated(trc_waterstorage)) THEN
         IF (ncio_var_exist(file_restart, 'trc_waterstorage', readflag = .false.)) THEN
            IF (tracer_dim_matches(file_restart, 'trc_waterstorage')) THEN
               CALL ncio_read_vector(file_restart, 'trc_waterstorage', ntracers, landpatch, trc_waterstorage)
            ELSE
               IF (p_is_master) WRITE(*,*) &
                  'Tracer restart trc_waterstorage ntracers mismatch; will rebuild via tracer_save_storage sync.'
            ENDIF
         ENDIF
      ENDIF
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
         R_init = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
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

   SUBROUTINE write_land_tracer_restart (file_restart, maxsnl, nl_soil)
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil

      IF (ntracers <= 0) RETURN

      CALL ncio_write_vector(file_restart, 'trc_ldew_rain', 'tracer', ntracers, 'patch', landpatch, &
         trc_ldew_rain, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_ldew_snow', 'tracer', ntracers, 'patch', landpatch, &
         trc_ldew_snow, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_wliq_soisno', 'tracer', ntracers, 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, &
         trc_wliq_soisno, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_wice_soisno', 'tracer', ntracers, 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, &
         trc_wice_soisno, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_wa', 'tracer', ntracers, 'patch', landpatch, &
         trc_wa, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_wdsrf', 'tracer', ntracers, 'patch', landpatch, &
         trc_wdsrf, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_wetwat', 'tracer', ntracers, 'patch', landpatch, &
         trc_wetwat, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'trc_scv', 'tracer', ntracers, 'patch', landpatch, &
         trc_scv, DEF_REST_CompressLevel)
      ! Persist the irrigation-reservoir tracer. Allocated unconditionally
      ! by allocate_Tracer_Vars (stays at 0 when DEF_USE_IRRIGATION is
      ! off), so writing it always produces a valid variable — readers
      ! without irrigation simply ignore it. Without this write, a
      ! mid-run restart loses any reservoir tracer mass that has been
      ! debited by irrigation within the current step.
      IF (allocated(trc_waterstorage)) THEN
         CALL ncio_write_vector(file_restart, 'trc_waterstorage', 'tracer', ntracers, 'patch', landpatch, &
            trc_waterstorage, DEF_REST_CompressLevel)
      ENDIF
   END SUBROUTINE write_land_tracer_restart

END MODULE MOD_Tracer_Rest
