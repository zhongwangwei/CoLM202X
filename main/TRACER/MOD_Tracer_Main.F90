#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Main

   USE MOD_Precision
   USE MOD_Tracer_Defs
   USE MOD_Tracer_Vars
   USE MOD_Tracer_Precip
   USE MOD_Tracer_Evapo
   USE MOD_Tracer_SoilWater
   USE MOD_Tracer_Snow
   USE MOD_Tracer_Conservation
   USE MOD_Tracer_Hist
   USE MOD_Tracer_Rest
   USE MOD_Tracer_SoilInit
   USE MOD_SPMD_Task

   IMPLICIT NONE

   PUBLIC :: tracer_init, tracer_final

CONTAINS

   SUBROUTINE tracer_init (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, file_restart, waterstorage, init_month)

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
      character(len=*), intent(in), optional :: file_restart
      ! Optional: irrigation reservoir from MOD_Vars_TimeVariables.
      ! Passed by CoLM.F90 only under CROP+DEF_USE_IRRIGATION. Used for
      ! cold-start initialisation of trc_waterstorage; the restart path
      ! already handles persistence via read_land_tracer_restart.
      real(r8), intent(in), optional :: waterstorage(numpatch)
      integer,  intent(in), optional :: init_month
      logical :: found_restart, scv_missing, waterstorage_missing
      integer :: soil_init_month

      CALL tracer_defs_init()
      IF (ntracers <= 0) RETURN
      ! Defensive: if a previous tracer_init left state behind (e.g. a
      ! LULCC transition that rebuilt numpatch without an explicit
      ! tracer_final call), release it so the re-allocation below uses
      ! the new numpatch. Without this the allocate() call would trip
      ! Fortran's "already allocated" runtime check, OR the save-level
      ! snap_* arrays in MOD_Tracer_Conservation would keep the old
      ! numpatch and then misalign with the new tracer arrays.
      CALL deallocate_Tracer_Vars()
      CALL deallocate_tracer_conservation()
      CALL allocate_Tracer_Vars(numpatch, maxsnl, nl_soil)
      found_restart = .false.
      scv_missing   = .false.
      waterstorage_missing = .false.
      IF (present(file_restart)) THEN
         CALL read_land_tracer_restart(file_restart, maxsnl, nl_soil, &
            found_restart, scv_missing, waterstorage_missing)
      ENDIF
#ifdef USEMPI
      ! The control/master rank does not own vector restart blocks, so
      ! read_land_tracer_restart can only make a local decision there.
      ! Make restart usability global before any cold-start path that
      ! contains MPI collectives, including tracer soil initialisation.
      CALL mpi_allreduce(MPI_IN_PLACE, found_restart, 1, MPI_LOGICAL, &
         MPI_LAND, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, scv_missing, 1, MPI_LOGICAL, &
         MPI_LOR, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, waterstorage_missing, 1, MPI_LOGICAL, &
         MPI_LOR, p_comm_glb, p_err)
      scv_missing = scv_missing .and. found_restart
      waterstorage_missing = waterstorage_missing .and. found_restart
#endif
      IF (.not. found_restart) THEN
         soil_init_month = 1
         IF (present(init_month)) soil_init_month = init_month
         IF (present(waterstorage)) THEN
            CALL tracer_init_from_water(numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv, waterstorage)
         ELSE
            CALL tracer_init_from_water(numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv)
         ENDIF
         CALL tracer_soil_init_from_file(soil_init_month, numpatch, maxsnl, nl_soil, &
            wliq_soisno, wice_soisno)
      ELSE
         IF (waterstorage_missing .and. present(waterstorage)) THEN
            CALL tracer_init_waterstorage_from_ratio(numpatch, maxsnl, nl_soil, &
               ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
               wa, wdsrf, wetwat, scv, waterstorage)
         ENDIF
         IF (scv_missing) THEN
            ! Old-format restart: every other pool loaded from file, only
            ! trc_scv needs to be rebuilt from the hydrology scv.
            CALL tracer_init_scv_from_water(numpatch, maxsnl, nl_soil, &
               wliq_soisno, wice_soisno, scv)
         ENDIF
      ENDIF
   END SUBROUTINE tracer_init

   SUBROUTINE tracer_final ()
      CALL deallocate_Tracer_Vars()
      CALL deallocate_tracer_conservation()
      CALL tracer_defs_final()
   END SUBROUTINE tracer_final

END MODULE MOD_Tracer_Main
#endif
