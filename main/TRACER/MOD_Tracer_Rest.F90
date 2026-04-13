#include <define.h>

MODULE MOD_Tracer_Rest

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, delta_to_R, trc_tiny
   USE MOD_Tracer_Vars
   USE MOD_LandPatch, only: landpatch
   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_NetCDFVector, only: ncio_read_vector, ncio_write_vector

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_init_from_water (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat)

      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)

      integer  :: itrc, ip, j
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
            trc_wa    (itrc, ip) = max(wa(ip),     0._r8) * R_init
            trc_wdsrf (itrc, ip) = max(wdsrf(ip),  0._r8) * R_init
            trc_wetwat(itrc, ip) = max(wetwat(ip), 0._r8) * R_init
         ENDDO
      ENDDO
   END SUBROUTINE tracer_init_from_water

   SUBROUTINE read_land_tracer_restart (file_restart, maxsnl, nl_soil, found_restart)
      USE MOD_SPMD_Task, only: p_is_master
      USE MOD_NetCDFSerial, only: ncio_read_serial
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil
      logical, intent(out) :: found_restart
      logical :: file_ok
      integer :: ntracers_file
      character(len=32), allocatable :: names_file(:)
      integer :: i

      found_restart = .false.
      IF (ntracers <= 0) RETURN

      inquire(file = trim(file_restart), exist = file_ok)
      IF (.not. file_ok) RETURN
      IF (.not. ncio_var_exist(file_restart, 'trc_ldew_rain', readflag = .false.)) RETURN

      ! Check tracer count matches (dimension 'tracer' in the file)
      IF (ncio_var_exist(file_restart, 'tracer_names', readflag = .false.)) THEN
         ! Read tracer names and verify consistency
         CALL ncio_read_serial(file_restart, 'tracer_names', names_file)
         ntracers_file = size(names_file)
         IF (ntracers_file /= ntracers) THEN
            IF (p_is_master) WRITE(*,*) 'WARNING: tracer restart has ', &
               ntracers_file, ' tracers but current config has ', ntracers, &
               '. Falling back to cold-start.'
            deallocate(names_file)
            RETURN
         ENDIF
         ! Check names match
         DO i = 1, ntracers
            IF (TRIM(names_file(i)) /= TRIM(tracers(i)%name)) THEN
               IF (p_is_master) WRITE(*,*) 'WARNING: tracer restart name mismatch at index ', i, &
                  ': file=', TRIM(names_file(i)), ' config=', TRIM(tracers(i)%name), &
                  '. Falling back to cold-start.'
               deallocate(names_file)
               RETURN
            ENDIF
         ENDDO
         deallocate(names_file)
      ENDIF

      CALL ncio_read_vector(file_restart, 'trc_ldew_rain', ntracers, landpatch, trc_ldew_rain)
      CALL ncio_read_vector(file_restart, 'trc_ldew_snow', ntracers, landpatch, trc_ldew_snow)
      CALL ncio_read_vector(file_restart, 'trc_wliq_soisno', ntracers, nl_soil-maxsnl, landpatch, trc_wliq_soisno)
      CALL ncio_read_vector(file_restart, 'trc_wice_soisno', ntracers, nl_soil-maxsnl, landpatch, trc_wice_soisno)
      CALL ncio_read_vector(file_restart, 'trc_wa', ntracers, landpatch, trc_wa)
      CALL ncio_read_vector(file_restart, 'trc_wdsrf', ntracers, landpatch, trc_wdsrf)
      CALL ncio_read_vector(file_restart, 'trc_wetwat', ntracers, landpatch, trc_wetwat)
      IF (ncio_var_exist(file_restart, 'trc_scv', readflag = .false.)) THEN
         CALL ncio_read_vector(file_restart, 'trc_scv', ntracers, landpatch, trc_scv)
      ENDIF
      found_restart = .true.
   END SUBROUTINE read_land_tracer_restart

   SUBROUTINE write_land_tracer_restart (file_restart, maxsnl, nl_soil)
      USE MOD_NetCDFSerial, only: ncio_write_serial
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: maxsnl, nl_soil
      character(len=32) :: trc_names(ntracers)
      integer :: i

      IF (ntracers <= 0) RETURN

      ! Save tracer names for restart validation
      DO i = 1, ntracers
         trc_names(i) = tracers(i)%name
      ENDDO
      CALL ncio_write_serial(file_restart, 'tracer_names', trc_names, 'tracer', ntracers)

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
   END SUBROUTINE write_land_tracer_restart

END MODULE MOD_Tracer_Rest
