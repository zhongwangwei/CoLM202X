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
   USE MOD_Tracer_Rest
   USE MOD_Tracer_SoilInit
   USE MOD_SPMD_Task
#ifdef BGC
   USE MOD_Tracer_Methane_Registry, only: methane_registry_init, igas_ch4
   USE MOD_Tracer_Methane_State,    only: allocate_methane_state, &
      init_methane_wetland_fraction_cache, deallocate_methane_state, &
      read_methane_restart, initialize_methane_lake_soilc_from_surface
   USE MOD_Tracer_Methane_AccFlux,  only: allocate_methane_acc_fluxes, &
      deallocate_methane_acc_fluxes
   USE MOD_Tracer_Methane_Microbes, only: allocate_methane_microbes_state, &
      deallocate_methane_microbes_state, read_methane_microbes_restart
   USE MOD_Tracer_Methane_Const,    only: read_methane_namelist, &
      configure_methane_inundation_mode, DEF_METHANE
   USE MOD_Tracer_Methane_GIEMS,    only: allocate_methane_giems, &
      deallocate_methane_giems, read_methane_giems, giems_active
   USE MOD_Tracer_Methane_pH,       only: allocate_methane_ph, &
      deallocate_methane_ph, read_methane_ph_patch
   USE MOD_Tracer_Methane_VegOverride, only: allocate_wetland_aere_overrides, &
      deallocate_wetland_aere_overrides
   USE MOD_Namelist, only: DEF_USE_METHANE_para, DEF_file_METHANE_para, &
      DEF_file_GIEMS, DEF_wetland_finundation_scheme, DEF_USE_SpatialpH
   USE MOD_Vars_TimeInvariants, only: patchtype, lake_soilc_srf, patchlatr, patchlonr
#endif

   IMPLICIT NONE

   PUBLIC :: tracer_init, tracer_final

CONTAINS

#ifdef BGC
   SUBROUTINE tracer_methane_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)

      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer, intent(in) :: lc_year
      integer, intent(in) :: jdate(3)
      character(len=*), intent(in) :: casename
      character(len=*), intent(in) :: dir_restart
      character(len=*), intent(in) :: dir_landdata

      character(len=256) :: file_restart_trc
      character(len=14)  :: cdate_restart
      character(len=256) :: cyear_restart
      real(r8), allocatable :: giems_dummy_patch(:)

      ! Resolve igas_ch4 / igas_o2 / igas_co2 from the registered tracer
      ! list. methane_registry_init is a no-op when none of CH4/O2/CO2 appear
      ! in DEF_TRACER_NAMES.
      CALL methane_registry_init()
      IF (igas_ch4 <= 0) RETURN

      IF (DEF_USE_METHANE_para) THEN
         CALL read_methane_namelist (DEF_file_METHANE_para)
      END IF
      CALL configure_methane_inundation_mode ()
      CALL allocate_methane_state (numpatch)
      CALL init_methane_wetland_fraction_cache (numpatch)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL allocate_methane_microbes_state (numpatch)
      ENDIF
      CALL allocate_methane_acc_fluxes (numpatch)

      write(cyear_restart,'(i4.4)') lc_year
      write(cdate_restart,'(i4.4,"-",i3.3,"-",i5.5)') jdate(1), jdate(2), jdate(3)
      file_restart_trc = trim(dir_restart)//'/'//trim(cdate_restart)//'/'//trim(casename)// &
                         '_restart_'//trim(cdate_restart)//'_lc'//trim(cyear_restart)//'.nc'
      CALL read_methane_restart (file_restart_trc)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL read_methane_microbes_restart (file_restart_trc)
      ENDIF
      CALL initialize_methane_lake_soilc_from_surface (patchtype, lake_soilc_srf, DEF_METHANE%allowlakeprod)

      ! Optional GIEMS-MC monthly wetland climatology for scheme 5.
      ! Scheme 5 is observation-driven and must not silently degrade to
      ! finundated=0 when its file is absent or unreadable.
      CALL allocate_methane_giems (numpatch)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         IF (p_is_worker .and. allocated(patchlatr) .and. allocated(patchlonr) .and. &
             size(patchlatr) >= numpatch .and. size(patchlonr) >= numpatch) THEN
            CALL read_methane_giems (DEF_file_GIEMS, patchlatr, patchlonr, numpatch)
         ELSE
            allocate(giems_dummy_patch(0))
            CALL read_methane_giems (DEF_file_GIEMS, giems_dummy_patch, giems_dummy_patch, 0)
            deallocate(giems_dummy_patch)
         ENDIF
         IF (.not. giems_active) THEN
            CALL CoLM_stop (' ***** ERROR: GIEMS file could not be loaded for methane scheme 5.')
         ENDIF
      ENDIF

      ! Optional spatial soil pH for the Dunfield 1993 pH factor in
      ! methane production. mksrfdata precomputes this from native PHH2O1
      ! into landdata/soil/<lc_year>/methane_ph_patches.nc; runtime reads
      ! only the patch vector to avoid global PHH2O broadcasts.
      CALL allocate_methane_ph (numpatch)
      IF (DEF_USE_SpatialpH) THEN
         CALL read_methane_ph_patch (trim(dir_landdata)//'/soil/'//trim(cyear_restart)// &
            '/methane_ph_patches.nc', numpatch)
      ENDIF

      ! Per-patch wetland aerenchyma overrides (5-zone climate proxy:
      ! tropical reed/papyrus, tropical swamp, temperate marsh, boreal
      ! fen, Sphagnum bog).  Populated by get_wetland_veg_proxy in
      ! methane_driver, consumed by SiteOxAere in Physics.
      CALL allocate_wetland_aere_overrides (numpatch)

   END SUBROUTINE tracer_methane_init

   SUBROUTINE tracer_methane_final ()

      IMPLICIT NONE

      IF (igas_ch4 <= 0) RETURN

      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()

   END SUBROUTINE tracer_methane_final
#endif

   SUBROUTINE tracer_init (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv, file_restart, waterstorage, init_month, &
      lc_year, jdate, casename, dir_restart, dir_landdata)

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
      integer, intent(in), optional :: lc_year
      integer, intent(in), optional :: jdate(3)
      character(len=*), intent(in), optional :: casename
      character(len=*), intent(in), optional :: dir_restart
      character(len=*), intent(in), optional :: dir_landdata
      logical :: found_restart, scv_missing, waterstorage_missing
      integer :: soil_init_month

      ! Defensive: if a previous tracer_init left state behind (e.g. a
      ! LULCC transition that rebuilt numpatch without an explicit
      ! tracer_final call), release it so the re-allocation below uses
      ! the new numpatch. Without this the allocate() call would trip
      ! Fortran's "already allocated" runtime check, OR the save-level
      ! snap_* arrays in MOD_Tracer_Conservation would keep the old
      ! numpatch and then misalign with the new tracer arrays.
      CALL deallocate_Tracer_Vars()
      CALL deallocate_tracer_conservation()
      CALL tracer_defs_init()
#ifdef BGC
      IF (present(lc_year) .and. present(jdate) .and. present(casename) .and. &
          present(dir_restart) .and. present(dir_landdata)) THEN
         CALL tracer_methane_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)
      ENDIF
#endif
      IF (ntracers <= 0) RETURN
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
#ifdef BGC
      CALL tracer_methane_final ()
#endif
      CALL deallocate_Tracer_Vars()
      CALL deallocate_tracer_conservation()
      CALL tracer_defs_final()
   END SUBROUTINE tracer_final

END MODULE MOD_Tracer_Main
#endif
