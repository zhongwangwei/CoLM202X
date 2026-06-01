#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane
!=======================================================================
! CH4 reactive-tracer implementation behind the generic reactive dispatch.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: p_is_worker, p_is_master, CoLM_stop
   USE MOD_Tracer_Defs, only: tracers, tracer_param_file_for_index, tracer_lower, tracer_upper
   USE MOD_Namelist, only: DEF_file_GIEMS, DEF_wetland_finundation_scheme
   USE MOD_Vars_TimeInvariants, only: patchtype, lake_soilc_srf, patchlatr, patchlonr
   USE MOD_Tracer_Reactive_Methane_Registry, only: methane_registry_init, methane_registry_refresh, igas_ch4
   USE MOD_Tracer_Reactive_Methane_State,    only: allocate_methane_state, &
      init_methane_wetland_fraction_cache, deallocate_methane_state, &
      read_methane_restart, write_methane_restart, initialize_methane_lake_soilc_from_surface, &
      save_methane_lulcc_state, remap_methane_lulcc_state, &
      publish_methane_levee_flood_patch, publish_methane_flood_patch
   USE MOD_Tracer_Reactive_Methane_AccFlux,  only: allocate_methane_acc_fluxes, &
      deallocate_methane_acc_fluxes, flush_methane_acc_fluxes, accumulate_methane_fluxes, &
      read_methane_accflux_restart, write_methane_accflux_restart
   USE MOD_Tracer_Reactive_Methane_Microbes, only: allocate_methane_microbes_state, &
      deallocate_methane_microbes_state, read_methane_microbes_restart, &
      write_methane_microbes_restart, save_methane_microbes_lulcc_state, &
      remap_methane_microbes_lulcc_state
   USE MOD_Tracer_Reactive_Methane_Const,    only: read_methane_namelist, &
      configure_methane_inundation_mode, DEF_METHANE
   USE MOD_Tracer_Reactive_Methane_GIEMS,    only: allocate_methane_giems, &
      deallocate_methane_giems, read_methane_giems, giems_active
   USE MOD_Tracer_Reactive_Methane_pH,       only: allocate_methane_ph, &
      deallocate_methane_ph, read_methane_ph_patch
   USE MOD_Tracer_Reactive_Methane_VegOverride, only: allocate_wetland_aere_overrides, &
      deallocate_wetland_aere_overrides
   USE MOD_Tracer_Reactive_Methane_Impl, only: ch4_impl_lake_step, &
      ch4_impl_wetland_decomp, ch4_impl_soil_step, ch4_impl_report
   USE MOD_Tracer_Reactive_Methane_Hist, only: methane_reactive_history

   IMPLICIT NONE
   PRIVATE

   logical, save :: registry_init_reported = .false.
   character(len=512), save :: last_methane_ph_patch_file = ''

   PUBLIC :: ch4_reactive_name, ch4_reactive_aliases
   PUBLIC :: ch4_reactive_has_name
   PUBLIC :: ch4_reactive_refresh_registry
   PUBLIC :: ch4_reactive_init, ch4_reactive_final
   PUBLIC :: ch4_reactive_lake_step, ch4_reactive_wetland_decomp
   PUBLIC :: ch4_reactive_soil_step, ch4_reactive_report
   PUBLIC :: ch4_reactive_write_restart, ch4_reactive_read_restart
   PUBLIC :: ch4_reactive_flush_acc_fluxes, ch4_reactive_accumulate_fluxes
   PUBLIC :: ch4_reactive_save_lulcc_state, ch4_reactive_remap_lulcc_state
   PUBLIC :: ch4_reactive_publish_levee_flood, ch4_reactive_publish_flood
   PUBLIC :: ch4_register_reactive_callbacks

CONTAINS

   SUBROUTINE ch4_register_reactive_callbacks ()

      USE MOD_Tracer_Reactive, only: register_reactive_callbacks
      IMPLICIT NONE

      CALL register_reactive_callbacks (ch4_reactive_name(), ch4_reactive_aliases(), &
         has_fn=ch4_reactive_has_name, refresh_fn=ch4_reactive_refresh_registry, init_fn=ch4_reactive_init, &
         read_restart_fn=ch4_reactive_read_restart, write_restart_fn=ch4_reactive_write_restart, &
         lake_step_fn=ch4_reactive_lake_step, wetland_decomp_fn=ch4_reactive_wetland_decomp, &
         soil_step_fn=ch4_reactive_soil_step, report_fn=ch4_reactive_report, &
         flush_acc_fluxes_fn=ch4_reactive_flush_acc_fluxes, &
         accumulate_fluxes_fn=ch4_reactive_accumulate_fluxes, &
         history_fn=methane_reactive_history, save_lulcc_fn=ch4_reactive_save_lulcc_state, &
         remap_lulcc_fn=ch4_reactive_remap_lulcc_state, &
         publish_levee_flood_fn=ch4_reactive_publish_levee_flood, &
         publish_flood_fn=ch4_reactive_publish_flood, final_fn=ch4_reactive_final)

   END SUBROUTINE ch4_register_reactive_callbacks

   character(len=32) FUNCTION ch4_reactive_name ()

      IMPLICIT NONE

      ch4_reactive_name = 'CH4'

   END FUNCTION ch4_reactive_name

   character(len=128) FUNCTION ch4_reactive_aliases ()

      IMPLICIT NONE

      ch4_reactive_aliases = 'METHANE'

   END FUNCTION ch4_reactive_aliases

   logical FUNCTION ch4_reactive_has ()

      IMPLICIT NONE

      ch4_reactive_has = (igas_ch4 > 0)

   END FUNCTION ch4_reactive_has

   logical FUNCTION ch4_reactive_has_name (tracer_name)

      IMPLICIT NONE
      character(len=*), intent(in) :: tracer_name

      ch4_reactive_has_name = ch4_reactive_is_alias(tracer_name) .and. ch4_reactive_has()

   END FUNCTION ch4_reactive_has_name

   SUBROUTINE ch4_reactive_refresh_registry ()

      IMPLICIT NONE

      IF (.not. registry_init_reported .and. allocated(tracers)) THEN
         CALL methane_registry_init ()
         registry_init_reported = .true.
      ELSE
         CALL methane_registry_refresh ()
      ENDIF

   END SUBROUTINE ch4_reactive_refresh_registry

   SUBROUTINE ch4_reactive_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)

      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer, intent(in) :: lc_year
      integer, intent(in) :: jdate(3)
      character(len=*), intent(in) :: casename
      character(len=*), intent(in) :: dir_restart
      character(len=*), intent(in) :: dir_landdata

      character(len=256) :: cyear_restart
      character(len=256) :: file_param
      logical :: use_param
      real(r8), allocatable :: giems_dummy_patch(:)

      IF (.not. ch4_reactive_has()) RETURN

      CALL ch4_reactive_resolve_param_file (use_param, file_param)
      IF (use_param) THEN
         CALL read_methane_namelist (file_param)
      END IF
      CALL configure_methane_inundation_mode ()

#ifndef CROP
      IF (DEF_METHANE%enable_rice_paddy) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: DEF_METHANE%enable_rice_paddy requires compiling with CROP.'
         CALL CoLM_stop (' ***** ERROR: rice-paddy methane requested without CROP')
      ENDIF
#endif

      ! ch4_reactive_init can be re-entered by LULCC/restart workflows.
      ! Drop any previous CH4 allocatables before allocating the current
      ! numpatch layout; the deallocators are no-ops when not allocated.
      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()

      CALL allocate_methane_state (numpatch)
      CALL init_methane_wetland_fraction_cache (numpatch)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL allocate_methane_microbes_state (numpatch)
      ENDIF
      CALL allocate_methane_acc_fluxes (numpatch)

      write(cyear_restart,'(i4.4)') lc_year
      IF (p_is_worker .and. allocated(patchtype) .and. allocated(lake_soilc_srf)) THEN
         CALL initialize_methane_lake_soilc_from_surface (patchtype, lake_soilc_srf, DEF_METHANE%allowlakeprod)
      ENDIF

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

      last_methane_ph_patch_file = ''
      CALL allocate_methane_ph (numpatch)
      IF (DEF_METHANE%use_spatial_ph) THEN
         last_methane_ph_patch_file = trim(dir_landdata)//'/soil/'//trim(cyear_restart)// &
            '/methane_ph_patches.nc'
         CALL read_methane_ph_patch (trim(last_methane_ph_patch_file), numpatch)
      ENDIF

      CALL allocate_wetland_aere_overrides (numpatch)

   END SUBROUTINE ch4_reactive_init


   SUBROUTINE ch4_reactive_resolve_param_file (use_param, file_param)

      IMPLICIT NONE
      logical, intent(out) :: use_param
      character(len=*), intent(out) :: file_param

      logical :: found

      CALL tracer_param_file_for_index (igas_ch4, ch4_reactive_aliases(), file_param, found)
      use_param = found .and. len_trim(file_param) > 0 .and. trim(tracer_lower(file_param)) /= 'null'
      IF (.not. use_param) THEN
         CALL CoLM_stop (' ***** ERROR: CH4 requires DEF_TRACER_PARAM_FILES to include a CH4 parameter file')
      ENDIF

   END SUBROUTINE ch4_reactive_resolve_param_file

   logical FUNCTION ch4_reactive_is_alias (raw_name)

      IMPLICIT NONE
      character(len=*), intent(in) :: raw_name

      character(len=32) :: name

      name = tracer_upper(raw_name)
      ch4_reactive_is_alias = trim(name) == 'CH4' .or. trim(name) == 'METHANE'

   END FUNCTION ch4_reactive_is_alias

   SUBROUTINE ch4_reactive_read_restart (file_restart)

      USE MOD_NetCDFSerial, only: ncio_var_exist
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      logical :: file_has_pools

      IF (.not. ch4_reactive_has()) RETURN
      CALL read_methane_restart (file_restart)
      CALL read_methane_accflux_restart (file_restart)

      ! C4: warn on a microbial-pool restart<->runtime-flag mismatch. The read
      ! below is flag-gated, so without this a restart written with pools ON but
      ! resumed with use_microbial_pools=.false. silently drops the prognostic
      ! biomass pools (and the reverse silently cold-starts them from B_init),
      ! breaking restart reproducibility with no message. Master-only probe.
      IF (p_is_master) THEN
         file_has_pools = ncio_var_exist(file_restart, 'ch4_B_methanogen', readflag = .false.)
         IF (file_has_pools .and. (.not. DEF_METHANE%use_microbial_pools)) THEN
            WRITE(*,'(A)') 'WARNING: methane microbial-pool fields are present in the '// &
               'restart but use_microbial_pools=.false.; those pools are being ignored.'
         ELSEIF ((.not. file_has_pools) .and. DEF_METHANE%use_microbial_pools) THEN
            WRITE(*,'(A)') 'WARNING: use_microbial_pools=.true. but the restart has no '// &
               'microbial-pool fields; biomass cold-started from B_init (not reproducible).'
         ENDIF
      ENDIF

      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL read_methane_microbes_restart (file_restart)
      ENDIF

   END SUBROUTINE ch4_reactive_read_restart

   SUBROUTINE ch4_reactive_write_restart (file_restart, compress)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress

      IF (.not. ch4_reactive_has()) RETURN
      CALL write_methane_restart(file_restart, compress)
      CALL write_methane_accflux_restart(file_restart, compress)
      CALL write_methane_microbes_restart(file_restart, compress)

   END SUBROUTINE ch4_reactive_write_restart

   SUBROUTINE ch4_reactive_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim_phy
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub

      IF (ch4_reactive_has()) THEN
         CALL ch4_impl_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)
      ENDIF

   END SUBROUTINE ch4_reactive_lake_step

   SUBROUTINE ch4_reactive_wetland_decomp (ipatch)

      IMPLICIT NONE
      integer, intent(in) :: ipatch

      IF (ch4_reactive_has()) CALL ch4_impl_wetland_decomp (ipatch)

   END SUBROUTINE ch4_reactive_wetland_decomp

   SUBROUTINE ch4_reactive_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

      IF (ch4_reactive_has()) CALL ch4_impl_soil_step (istep_local, ipatch, idate, deltim)

   END SUBROUTINE ch4_reactive_soil_step

   SUBROUTINE ch4_reactive_report ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL ch4_impl_report ()

   END SUBROUTINE ch4_reactive_report

   SUBROUTINE ch4_reactive_flush_acc_fluxes ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL flush_methane_acc_fluxes ()

   END SUBROUTINE ch4_reactive_flush_acc_fluxes

   SUBROUTINE ch4_reactive_accumulate_fluxes ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL accumulate_methane_fluxes ()

   END SUBROUTINE ch4_reactive_accumulate_fluxes

   SUBROUTINE ch4_reactive_save_lulcc_state ()

      IMPLICIT NONE

      IF (.not. ch4_reactive_has()) RETURN
      CALL save_methane_lulcc_state ()
      IF (DEF_METHANE%use_microbial_pools) CALL save_methane_microbes_lulcc_state ()

   END SUBROUTINE ch4_reactive_save_lulcc_state

   SUBROUTINE ch4_reactive_remap_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
      lccpct_patches, old_patch_area)

      IMPLICIT NONE
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,:)
      real(r8), intent(in), optional :: old_patch_area(:)
      integer :: nnew
      real(r8), allocatable :: giems_dummy_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      CALL remap_methane_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
         lccpct_patches, old_patch_area)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL remap_methane_microbes_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
            lccpct_patches, old_patch_area)
      ENDIF

      ! LULCC can change the worker-local patch count.  Methane state and
      ! microbial pools are remapped above; all other patch-sized methane
      ! buffers must be resized in the same transaction.  Leaving the old
      ! accumulator size would make acc1d(state, accumulator) write past the
      ! accumulator when numpatch grows, or keep stale/mis-aligned entries
      ! when it shrinks.
      nnew = size(patchclass_new)
      CALL flush_methane_acc_fluxes ()
      CALL deallocate_methane_acc_fluxes ()
      CALL allocate_methane_acc_fluxes (nnew)

      CALL deallocate_methane_giems ()
      CALL allocate_methane_giems (nnew)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         IF (p_is_worker .and. allocated(patchlatr) .and. allocated(patchlonr) .and. &
             size(patchlatr) >= nnew .and. size(patchlonr) >= nnew) THEN
            CALL read_methane_giems (DEF_file_GIEMS, patchlatr, patchlonr, nnew)
         ELSE
            allocate(giems_dummy_patch(0))
            CALL read_methane_giems (DEF_file_GIEMS, giems_dummy_patch, giems_dummy_patch, 0)
            deallocate(giems_dummy_patch)
         ENDIF
         IF (.not. giems_active) THEN
            CALL CoLM_stop (' ***** ERROR: GIEMS file could not be loaded for methane scheme 5.')
         ENDIF
      ENDIF

      CALL deallocate_methane_ph ()
      CALL allocate_methane_ph (nnew)
      IF (DEF_METHANE%use_spatial_ph .and. len_trim(last_methane_ph_patch_file) > 0) THEN
         CALL read_methane_ph_patch (trim(last_methane_ph_patch_file), nnew)
      ENDIF

      CALL deallocate_wetland_aere_overrides ()
      CALL allocate_wetland_aere_overrides (nnew)

   END SUBROUTINE ch4_reactive_remap_lulcc_state

   SUBROUTINE ch4_reactive_publish_levee_flood (fldfrc_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      CALL publish_methane_levee_flood_patch (fldfrc_patch)

   END SUBROUTINE ch4_reactive_publish_levee_flood

   SUBROUTINE ch4_reactive_publish_flood (fldfrc_patch, flddph_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)
      real(r8), intent(in) :: flddph_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      CALL publish_methane_flood_patch (fldfrc_patch, flddph_patch)

   END SUBROUTINE ch4_reactive_publish_flood

   SUBROUTINE ch4_reactive_final ()

      IMPLICIT NONE

      registry_init_reported = .false.
      IF (.not. ch4_reactive_has()) RETURN

      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()

   END SUBROUTINE ch4_reactive_final


END MODULE MOD_Tracer_Reactive_Methane
#endif
