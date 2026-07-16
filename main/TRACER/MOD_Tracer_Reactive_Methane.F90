#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane
!=======================================================================
! CH4 reactive-tracer implementation behind the generic reactive dispatch.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: p_is_worker, p_is_master, CoLM_stop
#ifdef USEMPI
   USE MPI
   USE MOD_SPMD_Task, only: p_comm_glb, p_err
#endif
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
      validate_methane_microbes_restart_values, &
      write_methane_microbes_restart, save_methane_microbes_lulcc_state, &
      remap_methane_microbes_lulcc_state
   USE MOD_Tracer_Reactive_Methane_Const,    only: read_methane_namelist, &
      configure_methane_inundation_mode, methane_history_accumulation_mode, DEF_METHANE
   USE MOD_Tracer_Reactive_Methane_GIEMS,    only: allocate_methane_giems, &
      deallocate_methane_giems, read_methane_giems, giems_active
   USE MOD_Tracer_Reactive_Methane_pH,       only: allocate_methane_ph, &
      deallocate_methane_ph, read_methane_ph_patch
   USE MOD_Tracer_Reactive_Methane_VegOverride, only: allocate_wetland_aere_overrides, &
      deallocate_wetland_aere_overrides
   USE MOD_Tracer_Reactive_Methane_Impl, only: ch4_impl_lake_step, &
      ch4_impl_wetland_decomp, ch4_impl_soil_step
   USE MOD_Tracer_Reactive_Methane_Hist, only: methane_reactive_history

   IMPLICIT NONE
   PRIVATE

   integer, parameter :: METHANE_RESTART_SCHEMA_VERSION = 4
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
   PUBLIC :: ch4_reactive_reload_lulcc_inputs
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
         reload_lulcc_fn=ch4_reactive_reload_lulcc_inputs, &
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
      IF (p_is_worker .and. allocated(patchtype)) THEN
         IF (DEF_METHANE%allowlakeprod .and. count(patchtype == 4) > 0 .and. &
             .not. allocated(lake_soilc_srf)) THEN
            CALL CoLM_stop (' ***** ERROR: lake CH4 production requires allocated lake_soilc surface data.')
         ENDIF
         IF (allocated(lake_soilc_srf)) THEN
            CALL initialize_methane_lake_soilc_from_surface (patchtype, lake_soilc_srf, DEF_METHANE%allowlakeprod)
         ENDIF
      ENDIF

      CALL allocate_methane_giems (numpatch)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         CALL require_methane_giems_patch_coords (numpatch)
         IF (p_is_worker .and. numpatch > 0) THEN
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

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_vector_group_presence, ncio_set_complete_require_present
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, parameter :: nmicrobe_restart_fields = 25
      character(len=40), parameter :: microbe_restart_fields(nmicrobe_restart_fields) = &
         [character(len=40) :: &
          'ch4_B_methanogen', 'ch4_B_methanotroph', &
          'ch4_B_methanogen_dormant', 'ch4_B_methanotroph_dormant', &
          'ch4_B_methanogen_soil', 'ch4_B_methanogen_rice', &
          'ch4_B_methanotroph_soil', 'ch4_B_methanotroph_rice', &
          'ch4_B_methanogen_dormant_soil', 'ch4_B_methanogen_dormant_rice', &
          'ch4_B_methanotroph_dormant_soil', 'ch4_B_methanotroph_dormant_rice', &
          'ch4_a_B_methanogen', 'ch4_a_B_methanotroph', &
          'ch4_a_B_methanogen_dormant', 'ch4_a_B_methanotroph_dormant', &
          'ch4_a_f_T_methanogen', 'ch4_a_f_S_methanogen', &
          'ch4_a_f_O2_methanogen', 'ch4_a_f_T_methanotroph', &
          'ch4_a_methanogen_growth_rate', 'ch4_a_methanotroph_growth_rate', &
          'ch4_a_microbial_prod_potential', 'ch4_a_microbial_oxid_potential', &
          'ch4_a_methane_acc_num_microbe']
      logical :: file_has_pools, file_has_component_pools
      logical :: file_has_microbe_accumulators, strict_restart
      logical :: history_mode_changed
      logical :: microbe_field_present(nmicrobe_restart_fields)
      integer :: n_microbe_state_fields, n_microbe_component_fields
      integer :: n_microbe_accumulator_fields
      integer :: restart_schema

      IF (.not. ch4_reactive_has()) RETURN
      CALL validate_methane_restart_transaction (file_restart, strict_restart, restart_schema, &
         history_mode_changed)
      ! C4: warn on a microbial-pool restart<->runtime-flag mismatch. The read
      ! below is flag-gated, so without this a restart written with pools ON but
      ! resumed with use_microbial_pools=.false. silently drops the prognostic
      ! biomass pools (and the reverse silently cold-starts them from B_init),
      ! breaking restart reproducibility with no message. Master-only probe.
      CALL ncio_vector_group_presence(file_restart, microbe_restart_fields, &
         landpatch, microbe_field_present)
      n_microbe_state_fields = count(microbe_field_present(1:4))
      n_microbe_component_fields = count(microbe_field_present(5:12))
      n_microbe_accumulator_fields = count(microbe_field_present(13:nmicrobe_restart_fields))
      file_has_pools = n_microbe_state_fields == 4
      file_has_component_pools = n_microbe_component_fields == 8
      file_has_microbe_accumulators = &
         n_microbe_accumulator_fields == nmicrobe_restart_fields - 12
      IF (strict_restart .and. &
          ((n_microbe_state_fields /= 0 .and. n_microbe_state_fields /= 4) .or. &
           (n_microbe_component_fields /= 0 .and. n_microbe_component_fields /= 8) .or. &
           (n_microbe_accumulator_fields /= 0 .and. &
            n_microbe_accumulator_fields /= nmicrobe_restart_fields - 12) .or. &
           (file_has_pools .neqv. file_has_microbe_accumulators) .or. &
           (file_has_component_pools .and. .not. file_has_pools) .or. &
           (restart_schema >= 3 .and. (file_has_pools .neqv. file_has_component_pools)))) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: committed methane restart contains a partial/inconsistent microbial feature group.'
         CALL CoLM_stop()
      ENDIF
      IF (.not. strict_restart .and. &
          ((n_microbe_state_fields /= 0 .and. n_microbe_state_fields /= 4) .or. &
           (n_microbe_component_fields /= 0 .and. n_microbe_component_fields /= 8) .or. &
           (n_microbe_accumulator_fields /= 0 .and. &
            n_microbe_accumulator_fields /= nmicrobe_restart_fields - 12))) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'WARNING: legacy methane restart has a partial microbial feature group; missing fields cold-start.'
      ENDIF
      IF (p_is_master) THEN
         IF (file_has_pools .and. (.not. DEF_METHANE%use_microbial_pools)) THEN
            WRITE(*,'(A)') 'WARNING: methane microbial-pool fields are present in the '// &
               'restart but use_microbial_pools=.false.; those pools are being ignored.'
         ELSEIF ((.not. file_has_pools) .and. DEF_METHANE%use_microbial_pools) THEN
            WRITE(*,'(A)') 'WARNING: use_microbial_pools=.true. but the restart has no '// &
               'microbial-pool fields; biomass cold-started from B_init (not reproducible).'
         ENDIF
      ENDIF

      ! Core schema fields are mandatory.  Microbial pools are a feature-
      ! conditional group: if the checkpoint contains the group, strict reads
      ! require every member; if it does not, enabling pools cold-starts from
      ! B_init as documented above rather than reclassifying a valid checkpoint
      ! as corrupt.
      CALL ncio_set_complete_require_present (strict_restart)
      ! Pass the validated schema so migration defaults are restricted to
      ! fields introduced after schema 1.
      CALL read_methane_restart (file_restart, strict_restart, restart_schema)
      ! The accumulator reader uses the same schema boundary for the six new
      ! lake accumulators.
      CALL read_methane_accflux_restart (file_restart, file_has_pools, &
         file_has_microbe_accumulators, strict_restart, restart_schema)
      IF (restart_schema >= 2 .and. history_mode_changed) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'WARNING: methane history accumulation mode changed across restart; resetting the partial window.'
         CALL flush_methane_acc_fluxes ()
      ENDIF
      IF (file_has_pools) THEN
         IF (DEF_METHANE%use_microbial_pools) THEN
            CALL read_methane_microbes_restart (file_restart, strict_restart, restart_schema >= 3, &
               file_has_component_pools)
         ELSE
            CALL validate_methane_microbes_restart_values(file_restart, strict_restart, &
               file_has_component_pools)
         ENDIF
      ENDIF
      CALL ncio_set_complete_require_present (.false.)

   END SUBROUTINE ch4_reactive_read_restart

   SUBROUTINE ch4_reactive_write_restart (file_restart, compress)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress

      IF (.not. ch4_reactive_has()) RETURN
      ! Invalidate an older commit before overwriting any member field.  The
      ! final marker is written only after state, accumulators and optional
      ! microbial pools have all completed on every vector block.
      CALL write_methane_restart_marker(file_restart, 'ch4_restart_complete', 0._r8, compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_restart_schema', &
         real(METHANE_RESTART_SCHEMA_VERSION, r8), compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_history_accumulation_mode', &
         real(methane_history_accumulation_mode(), r8), compress)
      CALL write_methane_restart(file_restart, compress)
      CALL write_methane_accflux_restart(file_restart, compress)
      CALL write_methane_microbes_restart(file_restart, compress)
      CALL write_methane_restart_marker(file_restart, 'ch4_restart_complete', 1._r8, compress)

   END SUBROUTINE ch4_reactive_write_restart

   SUBROUTINE write_methane_restart_marker (file_restart, varname, value, compress)

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart, varname
      real(r8), intent(in) :: value
      integer, intent(in) :: compress
      real(r8), allocatable :: marker(:)

      IF (p_is_worker) THEN
         allocate(marker(landpatch%nset))
         marker(:) = value
      ELSE
         allocate(marker(0))
      ENDIF
      CALL ncio_write_vector(file_restart, varname, 'patch', landpatch, marker, compress)
      deallocate(marker)

   END SUBROUTINE write_methane_restart_marker

   SUBROUTINE validate_methane_restart_transaction (file_restart, strict_restart, restart_schema, &
         history_mode_changed)

      USE MOD_LandPatch, only: landpatch
      USE MOD_NetCDFVector, only: ncio_read_vector_complete, ncio_vector_group_presence
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master, p_comm_glb, p_err, &
         MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR
#else
      USE MOD_SPMD_Task, only: p_is_worker, p_is_master
#endif
      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      character(len=40), parameter :: restart_metadata_fields(3) = &
         [character(len=40) :: 'ch4_restart_schema', 'ch4_restart_complete', &
          'ch4_history_accumulation_mode']
      logical, intent(out) :: strict_restart
      integer, intent(out) :: restart_schema
      logical, intent(out) :: history_mode_changed
      real(r8), allocatable :: schema(:), committed(:), history_mode(:)
      logical :: has_schema, has_commit, has_history_mode, invalid_marker
      logical :: restart_metadata_present(3)
      integer :: local_schema

      strict_restart = .false.
      restart_schema = 0
      history_mode_changed = .false.
      ! Both markers absent identifies a legacy restart.  Any one-sided marker
      ! is an interrupted new-format write and must never fall back to legacy.
      CALL ncio_vector_group_presence(file_restart, restart_metadata_fields, &
         landpatch, restart_metadata_present)
      has_schema = restart_metadata_present(1)
      has_commit = restart_metadata_present(2)
      IF (.not. has_schema .and. .not. has_commit) RETURN
      IF (.not. has_schema .or. .not. has_commit) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: incomplete methane restart transaction metadata; refusing mixed checkpoint state.'
         CALL CoLM_stop()
      ENDIF

      CALL ncio_read_vector_complete(file_restart, 'ch4_restart_schema', landpatch, schema)
      CALL ncio_read_vector_complete(file_restart, 'ch4_restart_complete', landpatch, committed)
      local_schema = 0
      IF (p_is_worker .and. allocated(schema)) THEN
         IF (any(schema == real(METHANE_RESTART_SCHEMA_VERSION, r8))) THEN
            local_schema = METHANE_RESTART_SCHEMA_VERSION
         ELSEIF (any(schema == 3._r8)) THEN
            local_schema = 3
         ELSEIF (any(schema == 2._r8)) THEN
            local_schema = 2
         ELSEIF (any(schema == 1._r8)) THEN
            local_schema = 1
         ENDIF
      ENDIF
      restart_schema = local_schema
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, restart_schema, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
#endif
      invalid_marker = .false.
      IF (p_is_worker) THEN
         IF (allocated(schema)) invalid_marker = invalid_marker .or. &
            any(schema /= real(restart_schema, r8))
         IF (allocated(committed)) invalid_marker = invalid_marker .or. any(committed /= 1._r8)
      ENDIF
      invalid_marker = invalid_marker .or. &
         (restart_schema /= 1 .and. restart_schema /= 2 .and. restart_schema /= 3 .and. &
          restart_schema /= METHANE_RESTART_SCHEMA_VERSION)
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, invalid_marker, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)
#endif
      IF (invalid_marker) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: methane restart transaction is uncommitted or has an unsupported schema.'
         CALL CoLM_stop()
      ENDIF

      ! The coarse mode catches off/core/custom changes.  AccFlux separately
      ! validates the exact custom selector fingerprint because custom lists
      ! can accumulate different array families.
      has_history_mode = restart_metadata_present(3)
      IF (restart_schema >= 3 .and. .not. has_history_mode) THEN
         IF (p_is_master) WRITE(*,'(A)') &
            'ERROR: schema-3+ methane restart is missing the history accumulation mode.'
         CALL CoLM_stop()
      ENDIF
      history_mode_changed = .not. has_history_mode
      IF (has_history_mode) THEN
         CALL ncio_read_vector_complete(file_restart, 'ch4_history_accumulation_mode', &
            landpatch, history_mode)
         IF (p_is_worker .and. allocated(history_mode)) THEN
            history_mode_changed = any(history_mode /= &
               real(methane_history_accumulation_mode(), r8))
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, history_mode_changed, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_glb, p_err)
#endif
         IF (allocated(history_mode)) deallocate(history_mode)
      ENDIF
      IF (allocated(schema)) deallocate(schema)
      IF (allocated(committed)) deallocate(committed)
      strict_restart = .true.

   END SUBROUTINE validate_methane_restart_transaction

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

   SUBROUTINE ch4_reactive_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      IF (ch4_reactive_has()) CALL ch4_impl_wetland_decomp (ipatch, deltim)

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

      ! Generic tracer conservation is reported unconditionally from
      ! MOD_Tracer_LandPhase::tracer_report so non-CH4 tracer runs are covered too.

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
      lccpct_patches, new_patch_area, old_patch_area)

      IMPLICIT NONE
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,:)
      real(r8), intent(in), optional :: new_patch_area(:)
      real(r8), intent(in), optional :: old_patch_area(:)
      integer :: nnew

      IF (.not. ch4_reactive_has()) RETURN
      CALL remap_methane_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
         lccpct_patches, new_patch_area, old_patch_area)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL remap_methane_microbes_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
            lccpct_patches, new_patch_area, old_patch_area)
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

      CALL deallocate_wetland_aere_overrides ()
      CALL allocate_wetland_aere_overrides (nnew)

   END SUBROUTINE ch4_reactive_remap_lulcc_state

   SUBROUTINE ch4_reactive_reload_lulcc_inputs (lc_year, dir_landdata)

      IMPLICIT NONE
      integer, intent(in) :: lc_year
      character(len=*), intent(in) :: dir_landdata
      integer :: nnew
      character(len=4) :: cyear_lulcc
      real(r8), allocatable :: giems_dummy_patch(:)

      IF (.not. ch4_reactive_has()) RETURN
      IF (DEF_METHANE%use_spatial_ph .and. &
          (lc_year <= 0 .or. len_trim(dir_landdata) == 0)) THEN
         CALL CoLM_stop (' ***** ERROR: methane LULCC pH reload requires current landdata directory and year.')
      ENDIF

      ! GIEMS and spatial-pH loading are collective.  Derive the worker-local
      ! size here so non-worker ranks participate with zero-length vectors.
      nnew = 0
      IF (p_is_worker .and. allocated(patchtype)) nnew = size(patchtype)

      CALL deallocate_methane_giems ()
      CALL allocate_methane_giems (nnew)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         CALL require_methane_giems_patch_coords (nnew)
         IF (p_is_worker .and. nnew > 0) THEN
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
      IF (DEF_METHANE%use_spatial_ph) THEN
         write(cyear_lulcc,'(i4.4)') lc_year
         last_methane_ph_patch_file = trim(dir_landdata)//'/soil/'//trim(cyear_lulcc)// &
            '/methane_ph_patches.nc'
         CALL read_methane_ph_patch (trim(last_methane_ph_patch_file), nnew)
      ENDIF

   END SUBROUTINE ch4_reactive_reload_lulcc_inputs

   SUBROUTINE require_methane_giems_patch_coords (local_numpatch)

      IMPLICIT NONE
      integer, intent(in) :: local_numpatch
      logical :: invalid_coords

      invalid_coords = .false.
      IF (p_is_worker .and. local_numpatch > 0) THEN
         IF (.not. allocated(patchlatr) .or. .not. allocated(patchlonr)) THEN
            invalid_coords = .true.
         ELSE
            invalid_coords = size(patchlatr) < local_numpatch .or. &
                             size(patchlonr) < local_numpatch
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL MPI_Allreduce (MPI_IN_PLACE, invalid_coords, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)
#endif
      IF (invalid_coords) THEN
         CALL CoLM_stop (' ***** ERROR: GIEMS methane input requires coordinates for every worker-local patch.')
      ENDIF

   END SUBROUTINE require_methane_giems_patch_coords

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
      last_methane_ph_patch_file = ''

      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()

   END SUBROUTINE ch4_reactive_final


END MODULE MOD_Tracer_Reactive_Methane
#endif
