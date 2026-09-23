#include <define.h>

#ifdef GridRiverLakeFlow
#if defined(RIVERLAKE_PERF_TRACE) && !defined(RIVERLAKE_PERF_DIAG)
#error "RIVERLAKE_PERF_TRACE requires RIVERLAKE_PERF_DIAG"
#endif
MODULE MOD_Grid_RiverLakeFlow
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   River Lake flow.
!
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_RiverLakeTimeVars
   USE MOD_Grid_Reservoir
   USE MOD_Grid_RiverLakeHistState
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levsto, levdph, &
      levee_init, read_levee_restart, levee_apply_protected_flux, &
      levee_repartition_storage, levee_fldstg, &
      levee_visible_volume_from_stage, levee_final
   USE MOD_Grid_RiverLakeBifurcation, only: bifurcation_init, bifurcation_calc, &
      read_bifurcation_restart, bifurcation_final, bifurcation_invalidate_static_dn, &
      bif_hflux_sum, bif_hflux_lev, bif_lev_hflux_sum, bif_path_active
#ifdef TRACER
   USE MOD_Tracer_Lifecycle, only: tracer_lifecycle_route_has_active, tracer_lifecycle_route_init, &
      tracer_lifecycle_route_calc, tracer_lifecycle_route_final, tracer_lifecycle_route_diag_accumulate, &
      tracer_lifecycle_route_forcing_put, tracer_lifecycle_route_read_restart, &
      tracer_lifecycle_route_sediment_bif_accumulate, tracer_lifecycle_route_sediment_levee_repartition
#endif
#ifdef TRACER
   USE MOD_Tracer_RiverLake, only: river_lake_tracer_init, tracer_init_from_water, &
      tracer_input_from_runoff, &
      tracer_substep, tracer_flush_acc, tracer_limiter_stats, &
      read_tracer_restart, river_lake_tracer_final, acc_trc_inp, acc_rnof_ref, trc_mass, trc_inp_buf, trc_flux_out, &
      tracer_diag_accumulate_substep, &
      trc_levsto, trc_solid, trc_levsto_solid, trc_dry_drain, trc_reactive_source, &
      levee_tracer_repartition, equilibrate_river_tracer_cell, &
      get_cell_volume_dep => get_cell_volume, trc_conc_dep => trc_conc
      USE MOD_Tracer_Lifecycle, only: tracer_lifecycle_publish_levee_flood_patch, &
         tracer_lifecycle_publish_flood_patch, &
         tracer_lifecycle_has_levee_flood_publisher, tracer_lifecycle_has_flood_publisher
#endif
   IMPLICIT NONE

   real(r8), parameter :: RIVERMIN  = RIVERLAKE_DRY_DEPTH
   real(r8), parameter :: RIVERLAKE_FLOOD_MISSING_VALUE = -1.e30_r8
   real(r8), parameter :: ROUTING_PATHOLOGICAL_DT_FALLBACK = 10._r8

   real(r8), save :: acctime_rnof_max
   integer, save :: routing_zero_dt_warn_count = 0
   integer, save :: routing_mass_warn_count = 0

   ! acctime_rnof (scalar) and acc_rnof_uc (:) are owned by
   ! MOD_Grid_RiverLakeTimeVars (imported via the module-wide USE above)
   ! so their mid-period state is serialised by WRITE/READ_GridRiverLakeTimeVars.
   logical,  allocatable :: filter_rnof (:)
   ! Patch-side flood exchange for the optional land/routing feedback.  These
   ! are derived credits/fluxes, not additional hydrodynamic state.
   real(r8), allocatable :: flood_depth_patch(:), flood_fraction_patch(:)
   real(r8), allocatable :: flood_evap_patch(:), flood_infil_patch(:)
   ! A publication is a finite credit, retained until the next routing call.
   ! The two donor pools must not be conflated by the reverse land mapping.
   real(r8), allocatable :: flood_visible_uc(:), flood_protected_uc(:)
   logical, allocatable :: flood_reservoir_uc(:)
   real(r8), allocatable :: flood_credit_patch(:)
   real(r8), allocatable :: flood_grid_area(:)
   real(r8), allocatable :: flood_evap_acc(:), flood_infil_acc(:)
#ifdef TRACER
   ! Derived publication credits; the prognostic pools remain in MOD_Tracer_RiverLake.
   real(r8), allocatable :: flood_tracer_credit_patch(:,:)
   real(r8), allocatable :: flood_tracer_evap_patch(:,:)
   real(r8), allocatable :: flood_tracer_land_patch(:,:)
   real(r8), allocatable :: flood_visible_tracer_uc(:,:), flood_protected_tracer_uc(:,:)
   integer, save :: flood_tracer_ledger_reports = 0
#endif
   real(r8), save :: flood_evap_period = 0._r8, flood_infil_period = 0._r8
CONTAINS

   ! ---------
   SUBROUTINE grid_riverlake_flow_init (start_year)

   USE MOD_LandPatch,           only: numpatch
   USE MOD_Forcing,             only: forcmask_pch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Vars_Global,         only: spval
#ifdef TRACER
   USE MOD_Tracer_Defs,          only: ntracers, tracer_uses_land_water_transport, &
      tracer_is_nonvolatile_solute, tracer_has_dissolved_limit
#endif
   IMPLICIT NONE

      integer, intent(in) :: start_year

#ifdef TRACER
      logical :: trc_restart_found
      logical :: has_flood_tracer
      logical, allocatable :: trc_missing(:)
      real(r8), allocatable :: wdsrf_safe(:), volresv_safe(:)
      integer,  allocatable :: ucat2resv_safe(:)
      logical, allocatable :: is_built_resv_init(:)
#endif
      logical :: bif_restart_loaded
      integer :: i, irsv
      real(r8), allocatable :: grid_area_local(:)

      acctime_rnof_max = DEF_GRIDBASED_ROUTING_MAX_DT
      ! acctime_rnof / acc_rnof_uc are allocated + zero-initialised in
      ! allocate_GridRiverLakeTimeVars and may then be overwritten by
      ! READ_GridRiverLakeTimeVars when a restart holds persisted values.
      ! Do NOT zero them here — that would clobber a mid-period recovery.

      ! excluding (patchtype >= 99), virtual patches and those forcing missed
      IF (p_is_worker) THEN
         allocate (filter_rnof (numpatch))
         IF (numpatch > 0) THEN
            filter_rnof = patchtype < 99
            filter_rnof = filter_rnof .and. patchmask
            IF (DEF_forcing%has_missing_value) THEN
               filter_rnof = filter_rnof .and. forcmask_pch
            ENDIF
         ENDIF
      ENDIF

#ifdef TRACER
      CALL tracer_lifecycle_route_init()
      IF (len_trim(gridriver_restart_file) > 0) THEN
         CALL tracer_lifecycle_route_read_restart(gridriver_restart_file)
      ENDIF
#endif

      ! Always call levee_init: when DEF_USE_LEVEE=.false. it allocates the
      ! levee arrays in inert state (has_levee=.false. everywhere), so guards
      ! like `IF (DEF_USE_LEVEE .and. has_levee(i) ...)` can safely evaluate
      ! both operands under ifx -check bounds.
      CALL levee_init()
      IF (len_trim(gridriver_restart_file) > 0) THEN
         CALL read_levee_restart(gridriver_restart_file, &
            restart_transaction_validated, restart_feature_manifest_present, &
            restart_levee_enabled, &
            fold_protected_to_visible = .not. DEF_USE_LEVEE, &
            volwater_ucat_io = volwater_ucat, &
            volwater_ucat_valid_io = volwater_ucat_valid, &
            wdsrf_ucat_in = wdsrf_ucat)
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         CALL bifurcation_init()
         IF (len_trim(gridriver_restart_file) > 0) THEN
            CALL read_bifurcation_restart(gridriver_restart_file, &
               wdsrf_ucat_prev_restart_found, bif_restart_loaded, &
               restart_transaction_validated, restart_feature_manifest_present, &
               restart_bifurcation_enabled, restart_levee_enabled)
            IF (.not. bif_restart_loaded) THEN
               ! Previous depth and pathway momentum form one numerical state
               ! unit. If either half is absent, cold-start both together.
               wdsrf_ucat_prev = wdsrf_ucat
               wdsrf_ucat_prev_valid = .true.
            ENDIF
         ENDIF
      ENDIF

      ! A legacy restart can mark tracked river storage valid while leaving a
      ! zero placeholder at a wet stage. Materialize the carrier before the
      ! first flood publication/debit or cold tracer seeding; otherwise a
      ! feedback debit can replace the wet stage with depth(0) while retaining
      ! the stage-derived tracer mass.
      IF (DEF_GridRiverLake_FloodFeedback .and. p_is_worker) THEN
         DO i = 1, numucat
            IF ((.not. volwater_ucat_valid) .or. &
                (volwater_ucat(i) <= 0._r8 .and. wdsrf_ucat(i) > RIVERMIN)) THEN
               IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  volwater_ucat(i) = levee_visible_volume_from_stage(i, wdsrf_ucat(i), levsto(i))
               ELSE
                  volwater_ucat(i) = floodplain_curve(i)%volume(wdsrf_ucat(i))
               ENDIF
            ENDIF
         ENDDO
         volwater_ucat_valid = .true.
         ! A built reservoir uses volresv rather than volwater_ucat for the
         ! first publication. Recover its legacy missing-value sentinel from
         ! stage before feedback can turn that stage into depth(0) as well.
         IF (allocated(volresv) .and. allocated(ucat2resv) .and. allocated(dam_build_year)) THEN
            DO i = 1, numucat
               IF (lake_type(i) /= 2) CYCLE
               irsv = ucat2resv(i)
               IF (irsv < 1 .or. irsv > size(volresv) .or. irsv > size(dam_build_year)) CYCLE
               IF (start_year < dam_build_year(irsv)) CYCLE
               IF (volresv(irsv) == spval) &
                  volresv(irsv) = floodplain_curve(i)%volume(wdsrf_ucat(i))
            ENDDO
         ENDIF
      ENDIF

#ifdef TRACER
         trc_restart_found = .false.
         CALL river_lake_tracer_init()
         IF (ntracers > 0) THEN
            allocate(trc_missing(ntracers))
            trc_missing = .true.   ! assume every tracer missing if no restart
            IF (len_trim(gridriver_restart_file) > 0) THEN
               CALL read_tracer_restart(gridriver_restart_file, trc_restart_found, trc_missing)
            ENDIF
            ! Cold start per-tracer: restart failures (no file, or specific
            ! tracer variables absent) fall back to init-from-water only for
            ! the missing tracers, preserving any that loaded successfully.
            ! wdsrf_ucat / volresv are populated by READ_GridRiverLakeTimeVars
            ! before this routine runs. volresv and ucat2resv are unallocated
            ! on reservoir-free workers, so wrap them into size-0 proxies.
            IF (any(trc_missing)) THEN
               IF (allocated(wdsrf_ucat)) THEN
                  allocate(wdsrf_safe(size(wdsrf_ucat)));   wdsrf_safe = wdsrf_ucat
               ELSE
                  allocate(wdsrf_safe(0))
               ENDIF
               IF (allocated(volresv)) THEN
                  allocate(volresv_safe(size(volresv)));    volresv_safe = volresv
               ELSE
                  allocate(volresv_safe(0))
               ENDIF
               IF (allocated(ucat2resv)) THEN
                  allocate(ucat2resv_safe(size(ucat2resv))); ucat2resv_safe = ucat2resv
               ELSE
                  allocate(ucat2resv_safe(0))
               ENDIF
               allocate(is_built_resv_init(size(wdsrf_safe)))
               is_built_resv_init = .false.
               IF (allocated(lake_type) .and. allocated(dam_build_year)) THEN
                  DO i = 1, min(size(is_built_resv_init), size(lake_type), size(ucat2resv_safe))
                     IF (lake_type(i) /= 2) CYCLE
                     irsv = ucat2resv_safe(i)
                     IF (irsv < 1 .or. irsv > size(dam_build_year)) CYCLE
                     is_built_resv_init(i) = start_year >= dam_build_year(irsv)
                  ENDDO
               ENDIF
               CALL tracer_init_from_water(wdsrf_safe, volresv_safe, ucat2resv_safe, trc_missing, is_built_resv_init)
               deallocate(wdsrf_safe, volresv_safe, ucat2resv_safe, is_built_resv_init)
            ENDIF
            deallocate(trc_missing)
         ENDIF
#endif

      ! Every reader has consumed the start-of-run restart.  Clear the path so a
      ! later hist_init (LULCC re-initialisation) cannot reload start-of-run
      ! history accumulators over the current ones.
      gridriver_restart_file = ''

      IF (p_is_worker) THEN
         allocate(flood_depth_patch(numpatch), flood_fraction_patch(numpatch), &
            flood_evap_patch(numpatch), flood_infil_patch(numpatch))
         flood_depth_patch = 0._r8
         flood_fraction_patch = 0._r8
         flood_evap_patch = 0._r8
         flood_infil_patch = 0._r8
         IF (DEF_GridRiverLake_FloodFeedback) THEN
#ifdef TRACER
            has_flood_tracer = .false.
            DO i = 1, ntracers
               IF (.not. tracer_uses_land_water_transport(i)) CYCLE
               has_flood_tracer = .true.
            ENDDO
            IF (has_flood_tracer) THEN
               allocate(flood_tracer_credit_patch(ntracers,numpatch), flood_tracer_evap_patch(ntracers,numpatch), &
                  flood_tracer_land_patch(ntracers,numpatch), &
                  flood_visible_tracer_uc(ntracers,numucat), flood_protected_tracer_uc(ntracers,numucat))
               flood_tracer_credit_patch = 0._r8
               flood_tracer_evap_patch = 0._r8
               flood_tracer_land_patch = 0._r8
            ENDIF
#endif
            allocate(flood_visible_uc(numucat), flood_protected_uc(numucat), &
               flood_credit_patch(numpatch), &
               flood_evap_acc(numpatch), flood_infil_acc(numpatch), flood_reservoir_uc(numucat), &
               flood_grid_area(numinpm))
            flood_evap_acc = 0._r8
            flood_infil_acc = 0._r8
            flood_evap_period = 0._r8
            flood_infil_period = 0._r8
         ENDIF
      ENDIF
      IF (DEF_GridRiverLake_FloodFeedback .and. p_is_worker) THEN
         allocate(grid_area_local(numinpm))
         grid_area_local = 0._r8
         flood_credit_patch = 1._r8
         IF (numpatch > 0) CALL worker_remap_data_pset2grid(remap_patch2inpm, &
            flood_credit_patch, grid_area_local, fillvalue=0._r8, filter=flood_credit_patch>0._r8)
         CALL worker_push_data(allreduce_inpm, grid_area_local, flood_grid_area, fillvalue=0._r8)
         flood_grid_area = max(flood_grid_area, push_ucat2inpm%sum_area)
      ENDIF
      IF (DEF_GridRiverLake_FloodFeedback) CALL publish_flood_feedback(start_year)

   END SUBROUTINE grid_riverlake_flow_init

   ! ---------
   SUBROUTINE grid_riverlake_flow (year, deltime)

   USE MOD_Utils
   USE MOD_Namelist,       only: DEF_Reservoir_Method, DEF_USE_LEVEE, DEF_USE_BIFURCATION, &
      DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT
   USE MOD_Vars_1DFluxes,  only: rnof
   USE MOD_Forcing,        only: forcmask_pch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_LandPatch,      only: elm_patch, numpatch
   USE MOD_Const_Physical, only: grav
   USE MOD_Vars_Global,    only: spval
   USE MOD_WorkerPushData, only: worker_push_real8_field_type
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
#ifdef TRACER
   USE MOD_Tracer_Defs,    only: ntracers, tracer_uses_land_water_transport, &
      tracer_has_dissolved_limit, tracer_equilibrate_dissolved
#endif
#ifdef TRACER
   USE MOD_Tracer_Vars,    only: trc_rnof_step
#endif
#ifdef TRACER
   USE MOD_Vars_1DForcing, only: forc_prc, forc_prl
#endif
   IMPLICIT NONE

   integer,  intent(in) :: year
   real(r8), intent(in) :: deltime

   ! Local Variables
   integer  :: i, irsv, ntimestep, ipth, i_up, itrc
   integer  :: lim_calls, lim_iter_sum, lim_iter_peak, lim_over_soft
   integer, save :: lim_diag_printed = 0
   real(r8) :: dt_this

   real(r8), allocatable :: rnof_gd(:)
   real(r8), allocatable :: rnof_uc(:)
   real(r8), allocatable :: trc_rnof_gd(:,:)
   real(r8), allocatable :: trc_rnof_uc(:,:)

#ifdef TRACER
   real(r8), allocatable, target :: prcp_gd(:)
   real(r8), allocatable, target :: prcp_uc(:)
   real(r8), allocatable, target :: prcp_area_gd(:), prcp_area_uc(:)
   real(r8), allocatable :: prcp_pch(:)
   logical, allocatable :: filter_prcp(:)
   type(worker_push_real8_field_type) :: prcp_push_fields(2)
   real(r8), allocatable :: particle_floodarea(:)
   real(r8), allocatable :: particle_protected_area(:)
   real(r8), allocatable :: particle_water_storage_start(:)
   real(r8), allocatable :: particle_water_storage(:)
   real(r8), allocatable :: particle_protected_start(:), particle_protected_end(:)
#endif

   logical,  allocatable :: is_built_resv(:)

   real(r8), allocatable, target :: wdsrf_next(:)
   real(r8), allocatable, target :: veloc_next(:)
   type(worker_push_real8_field_type) :: downstream_state_fields(2)

   real(r8), allocatable, target :: hflux_fc(:)
   real(r8), allocatable, target :: mflux_fc(:)
   real(r8), allocatable, target :: zgrad_dn(:)

   real(r8), allocatable, target :: hflux_resv(:)
   real(r8), allocatable, target :: mflux_resv(:)
   type(worker_push_real8_field_type) :: reservoir_flux_fields(2)

   real(r8), allocatable, target :: hflux_sumups(:)
   real(r8), allocatable, target :: mflux_sumups(:)
   real(r8), allocatable, target :: zgrad_sumups(:)
   type(worker_push_real8_field_type) :: upstream_flux_fields(3)

   real(r8), allocatable :: sum_hflux_riv(:)
   real(r8), allocatable :: sum_hflux_base(:)
   real(r8), allocatable :: normal_outgoing_rate(:)
   real(r8), allocatable :: ordinary_scale(:), ordinary_scale_next(:)
   real(r8), allocatable :: sum_mflux_riv(:)
   real(r8), allocatable :: sum_zgrad_riv(:)

   real(r8) :: veloct_fc, height_fc
   real(r8) :: bedelv_fc, height_up, height_dn
   real(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: volwater, friction, floodarea
   real(r8) :: rivsto_hist
   real(r8) :: visible_hflux, protected_hflux, protected_clip
   real(r8) :: fldfrc_levee
   real(r8) :: vis_vol_bef_lv, levsto_bef_lv
   real(r8) :: vis_vol_bef_lv2, levsto_bef_lv2
   real(r8), allocatable :: volresv_safe(:)
   integer,  allocatable :: ucat2resv_safe(:)
   integer :: itrc_dep
   real(r8) :: frac_remove, trc_removed, vol_post
   real(r8), allocatable :: levee_floodarea(:)
   real(r8), allocatable :: total_floodarea(:)   ! general floodarea (levee+floodplain), per ucat
   real(r8), allocatable :: total_flooddepth(:)  ! floodplain water depth flddph [m], per ucat
   real(r8),  allocatable :: dt_res(:), dt_all(:)
   logical,   allocatable :: ucatfilter(:)
   logical :: loop_active, next_loop_active
   real(r8) :: totalvol_bef, totalvol_aft, totalrnof, totaldis
   real(r8) :: totalflood_evap, totalflood_infil
   real(r8) :: water_balance_err, water_balance_tol
   real(r8) :: water_balance_vec(6)
#ifdef CoLMDEBUG
   real(r8) :: totalclip
   real(r8), allocatable :: trc_mass_bef(:), trc_mass_aft(:)
   real(r8), allocatable :: trc_mass_inp(:), trc_mass_dis(:), trc_mass_reactive(:)
   real(r8) :: bif_flux_sum_total, bif_flux_sum_max, bif_protected_clip_sum
   integer  :: itrc_dbg
#endif

#ifdef CoLMDEBUG
#ifdef TRACER
      IF (ntracers > 0) THEN
         allocate (trc_mass_bef(ntracers), trc_mass_aft(ntracers), &
                   trc_mass_inp(ntracers), trc_mass_dis(ntracers), &
                   trc_mass_reactive(ntracers))
      ELSE
         allocate (trc_mass_bef(0), trc_mass_aft(0), &
                   trc_mass_inp(0), trc_mass_dis(0), trc_mass_reactive(0))
      END IF
#else
      allocate (trc_mass_bef(0), trc_mass_aft(0), &
                trc_mass_inp(0), trc_mass_dis(0), trc_mass_reactive(0))
#endif
      trc_mass_bef = 0._r8
      trc_mass_aft = 0._r8
      trc_mass_inp = 0._r8
      trc_mass_dis = 0._r8
      trc_mass_reactive = 0._r8
#endif

      totalvol_bef = 0._r8
      totalvol_aft = 0._r8
      totalrnof = 0._r8
      totaldis = 0._r8
      totalflood_evap = 0._r8
      totalflood_infil = 0._r8
      water_balance_err = 0._r8
      water_balance_tol = 0._r8

      IF (p_is_worker) THEN
         allocate (rnof_gd (numinpm))
         allocate (rnof_uc (numucat))
         rnof_gd = 0._r8
#ifdef TRACER
         IF (ntracers > 0) THEN
            allocate (trc_rnof_gd (ntracers, numinpm))
            allocate (trc_rnof_uc (ntracers, numucat))
            trc_rnof_gd = 0._r8
            trc_rnof_uc = 0._r8
         END IF
#endif

         ! Forcing coverage can change after initialization.  Water and tracer
         ! runoff must be remapped from exactly the same current patch support.
         IF (numpatch > 0) THEN
            filter_rnof = patchtype < 99 .and. patchmask
            IF (DEF_forcing%has_missing_value) filter_rnof = filter_rnof .and. forcmask_pch
            CALL worker_remap_data_pset2grid (remap_patch2inpm, rnof, rnof_gd, &
               fillvalue = 0., filter = filter_rnof)
         ENDIF

         IF (numinpm > 0) THEN
            WHERE (push_ucat2inpm%sum_area > 0)
               rnof_gd = rnof_gd / push_ucat2inpm%sum_area
            END WHERE
         ENDIF

         CALL worker_push_data (push_inpm2ucat, rnof_gd, rnof_uc, &
            fillvalue = 0., mode = 'sum')

#ifdef TRACER
         IF (ntracers > 0) THEN
            DO itrc = 1, ntracers
               IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
               IF (numpatch > 0) &
                  CALL worker_remap_data_pset2grid(remap_patch2inpm, trc_rnof_step(itrc, :), trc_rnof_gd(itrc, :), &
                     fillvalue = 0._r8, filter = filter_rnof)
               IF (numinpm > 0) THEN
                  WHERE (push_ucat2inpm%sum_area > 0._r8)
                     trc_rnof_gd(itrc, :) = trc_rnof_gd(itrc, :) / push_ucat2inpm%sum_area
                  END WHERE
               ENDIF
               CALL worker_push_data(push_inpm2ucat, trc_rnof_gd(itrc, :), trc_rnof_uc(itrc, :), &
                  fillvalue = 0._r8, mode = 'sum')
            ENDDO
         END IF
#endif

         IF (numucat > 0) THEN
            acc_rnof_uc = acc_rnof_uc + rnof_uc*1.e-3*deltime

            ! Accumulate tracer input associated with this runoff increment
#ifdef TRACER
                  IF (ntracers > 0) THEN
                     ! trc_rnof_uc is in R×mm (from land: rsur*ratio*dt, rsur in mm/s)
                     ! Area-weighted rnof_uc * 1.e-3 * deltime is in m3.
                     ! Convert area-weighted tracer runoff from mm to m3 too.
                     CALL tracer_input_from_runoff(rnof_uc*1.e-3*deltime, numucat, trc_rnof_uc*1.e-3)
                  ENDIF
#endif
         ENDIF

         deallocate(rnof_gd)
         deallocate(rnof_uc)
         IF (allocated(trc_rnof_gd)) deallocate(trc_rnof_gd)
         IF (allocated(trc_rnof_uc)) deallocate(trc_rnof_uc)

#ifdef TRACER
         IF (tracer_lifecycle_route_has_active()) THEN
            ! Intensive rain forcing uses the same valid support in both mappings.
            ! Runoff above remains an extensive volume: never renormalize it here.
            allocate(prcp_pch(numpatch), filter_prcp(numpatch))
            allocate(prcp_gd(numinpm), prcp_area_gd(numinpm))
            allocate(prcp_uc(numucat), prcp_area_uc(numucat))
            prcp_pch = 0._r8
            prcp_gd = 0._r8
            prcp_area_gd = 0._r8
            IF (numpatch > 0) THEN
               filter_prcp = patchtype < 99 .and. patchmask
               IF (DEF_forcing%has_missing_value) filter_prcp = filter_prcp .and. forcmask_pch
               filter_prcp = filter_prcp .and. ieee_is_finite(forc_prc) .and. ieee_is_finite(forc_prl)
               ! Do not compare NaNs under trapping builds, or add missing sentinels.
               DO i = 1, numpatch
                  IF (.not. filter_prcp(i)) CYCLE
                  filter_prcp(i) = forc_prc(i) /= spval .and. forc_prl(i) /= spval
                  IF (filter_prcp(i)) filter_prcp(i) = forc_prc(i) >= 0._r8 .and. forc_prl(i) >= 0._r8
               ENDDO
               WHERE (filter_prcp) prcp_pch = forc_prc + forc_prl
               CALL worker_remap_data_pset2grid(remap_patch2inpm, prcp_pch, prcp_gd, &
                  fillvalue = 0._r8, filter = filter_prcp)
            ENDIF
            ! A dry but valid patch contributes area even when its rain is zero.
            prcp_pch = 1._r8
            IF (numpatch > 0) &
               CALL worker_remap_data_pset2grid(remap_patch2inpm, prcp_pch, prcp_area_gd, &
                  fillvalue = 0._r8, filter = filter_prcp)
            WHERE (push_ucat2inpm%sum_area > 0._r8)
               prcp_gd = prcp_gd / push_ucat2inpm%sum_area
               prcp_area_gd = prcp_area_gd / push_ucat2inpm%sum_area
            ELSEWHERE
               prcp_gd = 0._r8
               prcp_area_gd = 0._r8
            END WHERE
            ! Rain and its valid area share one mapping and one mode, so they go
            ! in a single batched exchange rather than two.
            prcp_push_fields(1)%send => prcp_gd
            prcp_push_fields(1)%recv => prcp_uc
            prcp_push_fields(1)%fillvalue = 0._r8
            prcp_push_fields(2)%send => prcp_area_gd
            prcp_push_fields(2)%recv => prcp_area_uc
            prcp_push_fields(2)%fillvalue = 0._r8
            CALL worker_push_data(push_inpm2ucat, prcp_push_fields, mode = 'sum')
            WHERE (prcp_area_uc > 0._r8)
               prcp_uc = prcp_uc / prcp_area_uc
            ELSEWHERE
               prcp_uc = 0._r8
            END WHERE
            ! Existing erosion uses uniform rain over the full catchment. Its
            ! sampled mean is weighted by valid area * time, not missing zeros.
            ! Zero coverage means no observation, distinct from observed zero rain.
            WHERE (topo_area > 0._r8)
               prcp_area_uc = prcp_area_uc / topo_area
            ELSEWHERE
               prcp_area_uc = 0._r8
            END WHERE
            CALL tracer_lifecycle_route_forcing_put(prcp_uc, deltime, prcp_area_uc)
            deallocate(prcp_pch, filter_prcp, prcp_gd, prcp_area_gd, prcp_uc, prcp_area_uc)
         ENDIF
#endif

      ENDIF


      IF (DEF_GridRiverLake_FloodFeedback .and. p_is_worker) THEN
         flood_evap_acc = flood_evap_acc + flood_evap_patch*deltime
         flood_infil_acc = flood_infil_acc + flood_infil_patch*deltime
         flood_evap_patch = 0._r8
         flood_infil_patch = 0._r8
      ENDIF
      IF (DEF_GridRiverLake_FloodFeedback) THEN
         CALL debit_flood_feedback(totalflood_evap, totalflood_infil)
         flood_evap_period = flood_evap_period + totalflood_evap
         flood_infil_period = flood_infil_period + totalflood_infil
         CALL publish_flood_feedback(year)
      ENDIF
      acctime_rnof = acctime_rnof + deltime

      IF (acctime_rnof+0.01 < acctime_rnof_max) THEN
         RETURN
      ENDIF


         IF (p_is_worker) THEN

            ! ROUTING_DT_BUFFERS_ARE_ALLOCATED_FOR_EMPTY_WORKERS:
            ! allocate all routing work arrays even on worker ranks with
            ! numucat==0.  Those ranks still participate in MPI collectives and
            ! worker_push_data calls while other workers are active; zero-length
            ! assumed-shape arrays are safe, unallocated arrays are not.
            allocate (is_built_resv (numucat))
            allocate (wdsrf_next    (numucat))
            allocate (veloc_next    (numucat))
            downstream_state_fields(1)%send => wdsrf_ucat
            downstream_state_fields(1)%recv => wdsrf_next
            downstream_state_fields(1)%fillvalue = spval
            downstream_state_fields(2)%send => veloc_riv
            downstream_state_fields(2)%recv => veloc_next
            downstream_state_fields(2)%fillvalue = 0._r8
            allocate (hflux_fc      (numucat))
            allocate (mflux_fc      (numucat))
            allocate (zgrad_dn      (numucat))
            allocate (sum_hflux_riv (numucat))
            IF (DEF_USE_BIFURCATION) allocate (sum_hflux_base(numucat))
            allocate (normal_outgoing_rate(numucat), ordinary_scale(numucat), ordinary_scale_next(numucat))
            allocate (sum_mflux_riv (numucat))
            allocate (sum_zgrad_riv (numucat))
            allocate (ucatfilter    (numucat))

            ! Always allocate levee_floodarea so guards
            ! `IF (DEF_USE_LEVEE .and. levee_floodarea(i) > 0.)` survive
            ! ifx -check bounds when LEVEE=off (entries stay zero, guard
            ! short-circuits logically).
            allocate (levee_floodarea (numucat))
            levee_floodarea = 0.

            ! General flood area (levee+floodplain) for methane scheme 7.
            allocate (total_floodarea (numucat))
            total_floodarea = 0.
            allocate (total_flooddepth (numucat))
            total_flooddepth = 0.
#ifdef TRACER
            ! Per-substep inputs of the particle (sediment) diagnostics; sized
            ! once per routing period instead of once per substep.
            allocate (particle_floodarea (numucat))
            allocate (particle_protected_area(numucat))
            allocate (particle_water_storage_start (numucat))
            allocate (particle_water_storage (numucat))
            allocate (particle_protected_start(numucat), particle_protected_end(numucat))
            particle_floodarea = 0._r8
            particle_protected_area = 0._r8
            particle_water_storage_start = 0._r8
            particle_water_storage = 0._r8
            particle_protected_start = 0._r8
            particle_protected_end = 0._r8
#endif

            allocate (hflux_sumups  (numucat))
            allocate (mflux_sumups  (numucat))
            allocate (zgrad_sumups  (numucat))

            upstream_flux_fields(1)%send => hflux_fc
            upstream_flux_fields(1)%recv => hflux_sumups
            upstream_flux_fields(1)%fillvalue = 0._r8
            upstream_flux_fields(2)%send => mflux_fc
            upstream_flux_fields(2)%recv => mflux_sumups
            upstream_flux_fields(2)%fillvalue = 0._r8
            upstream_flux_fields(3)%send => zgrad_dn
            upstream_flux_fields(3)%recv => zgrad_sumups
            upstream_flux_fields(3)%fillvalue = 0._r8

            IF (DEF_Reservoir_Method > 0) THEN
               allocate (hflux_resv (numucat))
               allocate (mflux_resv (numucat))
               reservoir_flux_fields(1)%send => hflux_resv
               reservoir_flux_fields(1)%recv => hflux_sumups
               reservoir_flux_fields(1)%fillvalue = 0._r8
               reservoir_flux_fields(2)%send => mflux_resv
               reservoir_flux_fields(2)%recv => mflux_sumups
               reservoir_flux_fields(2)%fillvalue = 0._r8
            ENDIF

            allocate (dt_res (numrivsys))
            allocate (dt_all (numrivsys))

            ! ROUTING_SAFE_ARRAYS_REUSED_PER_ROUTING_CALL:
            ! volresv and ucat2resv are unallocated on reservoir-free
            ! workers, but downstream routines use assumed-shape dummies.
            ! Allocate zero-length proxies once per routing call and refresh
            ! the reservoir volumes before each consumer instead of doing
            ! alloc/dealloc churn inside every substep/bif iteration.
            IF (allocated(volresv)) THEN
               allocate (volresv_safe(size(volresv)))
               volresv_safe = volresv
            ELSE
               allocate (volresv_safe(0))
            ENDIF
            IF (allocated(ucat2resv)) THEN
               allocate (ucat2resv_safe(size(ucat2resv)))
               ucat2resv_safe = ucat2resv
            ELSE
               allocate (ucat2resv_safe(0))
            ENDIF

         ! Tracer conservation: snapshot old state and the queued input
         ! before merging that input into the pending pool below.
#ifdef TRACER
         ! Reset this period diagnostic even without CoLMDEBUG; reactive
         ! chemistry can still produce increments in production builds.
         IF (allocated(trc_reactive_source)) trc_reactive_source = 0._r8
#endif
#ifdef CoLMDEBUG
#ifdef TRACER
         IF (numucat > 0) THEN
            DO itrc = 1, ntracers
               IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
               trc_mass_bef(itrc) = sum(trc_mass(itrc,:)) + sum(trc_inp_buf(itrc,:))
               IF (allocated(trc_levsto)) &
                  trc_mass_bef(itrc) = trc_mass_bef(itrc) + sum(trc_levsto(itrc,:))
               IF (allocated(trc_solid)) trc_mass_bef(itrc) = trc_mass_bef(itrc) &
                  + sum(trc_solid(itrc,:)) + sum(trc_levsto_solid(itrc,:))
               trc_mass_inp(itrc) = sum(acc_trc_inp(itrc,:))
               trc_mass_dis(itrc) = 0._r8
               trc_mass_reactive(itrc) = 0._r8
            ENDDO
         ELSE
            trc_mass_bef = 0._r8
            trc_mass_inp = 0._r8
            trc_mass_dis = 0._r8
            trc_mass_reactive = 0._r8
         ENDIF
#else
         trc_mass_bef = 0._r8
         trc_mass_inp = 0._r8
         trc_mass_dis = 0._r8
         trc_mass_reactive = 0._r8
#endif
#endif

#ifdef TRACER
         ! Water is added in full before adaptive routing. Merge its tracer
         ! counterpart exactly once at the same boundary. Positive pending
         ! mass is released by tracer_substep; signed debt remains queued.
         IF (numucat > 0) THEN
            DO itrc = 1, ntracers
               IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
               trc_inp_buf(itrc, :) = trc_inp_buf(itrc, :) + acc_trc_inp(itrc, :)
               acc_trc_inp(itrc, :) = 0._r8
            ENDDO
         ENDIF
#endif

         totalrnof = sum(acc_rnof_uc)
         totalvol_bef = 0._r8

         DO i = 1, numucat

            is_built_resv(i) = .false.
            IF (lake_type(i) == 2) THEN
               irsv = ucat2resv(i)
               IF (year >= dam_build_year(irsv)) THEN
                  is_built_resv(i) = .true.
                  IF (volresv(irsv) == spval) THEN
                     volresv(irsv) = floodplain_curve(i)%volume (wdsrf_ucat(i))
                  ELSE
                     wdsrf_ucat(i) = floodplain_curve(i)%depth (volresv(irsv))
                  ENDIF
               ENDIF
            ENDIF

            IF (.not. is_built_resv(i)) THEN
               momen_riv(i) = wdsrf_ucat(i) * veloc_riv(i)
               IF (volwater_ucat_valid) THEN
                  ! Persistent tracked volume (restored from restart or
                  ! previous call). Old restarts can contain volwater_ucat
                  ! as an all-zero placeholder even when stage is wet; rebuild
                  ! those cells from stage once, including levee visible water.
                  IF (volwater_ucat(i) > 0._r8 .or. wdsrf_ucat(i) <= RIVERMIN) THEN
                     volwater = volwater_ucat(i)
                  ELSEIF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                     volwater = levee_visible_volume_from_stage(i, wdsrf_ucat(i), levsto(i))
                  ELSE
                     volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
                  ENDIF
               ELSEIF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  volwater = levee_visible_volume_from_stage(i, wdsrf_ucat(i), levsto(i))
               ELSE
                  volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
               ENDIF
            ELSE
               ! water in reservoirs is assumued to be stationary.
               momen_riv(i) = 0
               veloc_riv(i) = 0
               volwater = volresv(ucat2resv(i))
            ENDIF

            IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
               totalvol_bef = totalvol_bef + volwater + levsto(i)
            ELSE
               totalvol_bef = totalvol_bef + volwater
            ENDIF

            volwater = volwater + acc_rnof_uc(i)

            IF (.not. is_built_resv(i)) THEN
               IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  ! volwater already includes runoff, while the matching
                  ! tracer is still pending in trc_inp_buf. Include that
                  ! net positive pool in the visible-side concentration
                  ! used for the first levee repartition.
                  CALL levee_repartition_storage(i, volwater, wdsrf_ucat(i), &
                     fldfrc_levee, vis_vol_bef_lv, levsto_bef_lv)
                  volwater_ucat(i) = volwater
#ifdef TRACER
                  IF (ntracers > 0) CALL tracer_lifecycle_route_sediment_levee_repartition(i, &
                     vis_vol_bef_lv, levsto_bef_lv, volwater_ucat(i), levsto(i))
                  ! The optional pending pool is unallocated when no tracer is
                  ! configured; do not form an array section before the callee.
                  IF (ntracers > 0) CALL levee_tracer_repartition(i, &
                     vis_vol_bef_lv, levsto_bef_lv, &
                     volwater_ucat(i), levsto(i), &
                     pending_trc_pool = trc_inp_buf(:, i))
#endif
               ELSE
                  volwater_ucat(i) = volwater
                  wdsrf_ucat(i) = floodplain_curve(i)%depth (volwater)
               ENDIF
               IF (wdsrf_ucat(i) > RIVERMIN) THEN
                   veloc_riv(i) = momen_riv(i) / wdsrf_ucat(i)
                ELSE
                   veloc_riv(i) = 0.
               ENDIF
            ELSE
               volresv(ucat2resv(i)) = volwater
            ENDIF

         ENDDO

         ! From here on volwater_ucat is defined for every non-reservoir cell
         ! in this routing call, even when no compatible restart volume existed.
         volwater_ucat_valid = .true.

         ntimestep = 0
            totaldis  = 0._r8
#ifdef CoLMDEBUG
            totalclip = 0._r8
            bif_flux_sum_total = 0._r8
            bif_flux_sum_max   = 0._r8
            bif_protected_clip_sum = 0._r8
#endif

         dt_res(:) = acctime_rnof

         ! CaMa-style routing uses one adaptive substep DT for the whole
         ! worker domain. Keep every worker in this loop until all river
         ! systems have exhausted dt_res, otherwise a worker that finishes
         ! early may skip collectives while others still need push_data or
         ! global-DT synchronization.

         ! Static bifurcation downstream path fields (riverbed elevation, levee
         ! mask, reservoir mask) are constant within this routing call; flag them
         ! stale so bifurcation_calc pushes them once on the first sub-step and
         ! reuses the cached path buffers thereafter.  Runs on every worker (sets
         ! a module flag only) to keep the one-per-call push collective-matched.
         IF (DEF_USE_BIFURCATION) CALL bifurcation_invalidate_static_dn ()

         ! Legacy restarts and cold starts have no previous BIF depth. Use the
         ! current stage for the first semi-implicit channel update; subsequent
         ! substeps advance this state at the same point as CaMa CALC_VARS_PRE.
         IF (DEF_USE_BIFURCATION .and. .not. wdsrf_ucat_prev_valid) THEN
            wdsrf_ucat_prev = wdsrf_ucat
            wdsrf_ucat_prev_valid = .true.
         ENDIF

         IF (DEF_USE_BIFURCATION) THEN
            ! acctime_rnof is broadcast on restart and advanced identically on
            ! every rank. Use that global scalar for loop entry so even an
            ! empty worker enters the collective-bearing BIF loop, without a
            ! separate MPI_LOR just to recover the same global decision.
            IF (.not. ieee_is_finite(acctime_rnof)) THEN
               loop_active = .true.
            ELSE
               loop_active = acctime_rnof > 0._r8
            ENDIF
         ELSE
            loop_active = any(dt_res > 0._r8)
         ENDIF

         DO WHILE (loop_active)

            ntimestep = ntimestep + 1

            ! River systems can finish at different substeps. Inactive cells
            ! are skipped below but still participate in the unfiltered MPI
            ! push, so clear last-substep fluxes before rebuilding active ones.
            hflux_fc = 0._r8
            mflux_fc = 0._r8
            zgrad_dn = 0._r8

            ! Water depth and velocity use the same one-to-one downstream
            ! mapping.  Pack them into one peer message per routing substep.
            CALL worker_push_data (push_next2ucat, downstream_state_fields)

            ! Preserve the original 60 s upper bound; the local CFL, storage,
            ! and momentum constraints below may shorten this further.
            dt_all(:) = min(dt_res(:), 60._r8)
            ! In BIF mode every active system starts with the same residual
            ! time: the first substep starts from acctime_rnof and every later
            ! one subtracts the globally reduced dt.  Synchronizing this
            ! initial residual is therefore redundant.  Reduce only once,
            ! after the local CFL/storage/momentum constraints below, before
            ! any cross-system BIF flux is evaluated.

            DO i = 1, numucat

               ! Bounds-guard irivsys: a corrupt/partial restart or
               ! malformed network metadata can leave irivsys(i)==0,
               ! negative, or > size(dt_all). Unlike the tracer_substep
               ! guards at MOD_Tracer_RiverLake.F90:482,569,786, this
               ! is the gatekeeper of ucatfilter — every downstream
               ! dt_all(irivsys(i)) use in this loop trusts ucatfilter.
               ! Treat an invalid entry as inactive instead of segfaulting.
               IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) THEN
                  ucatfilter(i) = dt_all(irivsys(i)) > 0
               ELSE
                  ucatfilter(i) = .false.
               ENDIF

               IF (.not. ucatfilter(i)) CYCLE

               sum_hflux_riv(i) = 0.
               sum_mflux_riv(i) = 0.
               sum_zgrad_riv(i) = 0.

               ! reservoir
               IF (is_built_resv(i)) THEN
                  hflux_fc(i) = 0.
                  mflux_fc(i) = 0.
                  zgrad_dn(i) = 0.
                  CYCLE
               ENDIF

               IF ((ucat_next(i) > 0) .or. (ucat_next(i) == -9)) THEN

                  IF (ucat_next(i) > 0) THEN
                     ! both rivers are dry.
                     IF ((wdsrf_ucat(i) < RIVERMIN) .and. (wdsrf_next(i) < RIVERMIN)) THEN
                        hflux_fc(i) = 0
                        mflux_fc(i) = 0
                        zgrad_dn(i) = 0
                        CYCLE
                     ENDIF
                  ENDIF

                  ! reconstruction of height of water near interface
                  IF (ucat_next(i) > 0) THEN
                     bedelv_fc = max(topo_rivelv(i), bedelv_next(i))
                     height_up = max(0., wdsrf_ucat(i)+topo_rivelv(i)-bedelv_fc)
                     height_dn = max(0., wdsrf_next(i)+bedelv_next(i)-bedelv_fc)
                  ELSEIF (ucat_next(i) == -9) THEN ! for river mouth
                     bedelv_fc = topo_rivelv(i)
                     height_up = wdsrf_ucat (i)
                     ! sea level is assumed to be 0. and sea bed is assumed to be negative infinity.
                     height_dn = max(0., - bedelv_fc)
                  ENDIF

                  ! velocity at river downstream face (middle region in Riemann problem)
                  veloct_fc = 0.5 * (veloc_riv(i) + veloc_next(i)) &
                     + sqrt(grav * height_up) - sqrt(grav * height_dn)

                  ! height of water at downstream face (middle region in Riemann problem)
                  height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn)) &
                     + 0.25 * (veloc_riv(i) - veloc_next(i))) ** 2

                  IF (height_up > 0) THEN
                     vwave_up = min(veloc_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                  ELSE
                     vwave_up = veloc_next(i) - 2.0 * sqrt(grav*height_dn)
                  ENDIF

                  IF (height_dn > 0) THEN
                     vwave_dn = max(veloc_next(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))
                  ELSE
                     vwave_dn = veloc_riv(i) + 2.0 * sqrt(grav*height_up)
                  ENDIF

                  hflux_up = veloc_riv(i)  * height_up
                  hflux_dn = veloc_next(i) * height_dn
                  mflux_up = veloc_riv(i)**2  * height_up + 0.5*grav * height_up**2
                  mflux_dn = veloc_next(i)**2 * height_dn + 0.5*grav * height_dn**2

                  IF (vwave_up >= 0.) THEN
                     hflux_fc(i) = outletwth(i) * hflux_up
                     mflux_fc(i) = outletwth(i) * mflux_up
                  ELSEIF (vwave_dn <= 0.) THEN
                     hflux_fc(i) = outletwth(i) * hflux_dn
                     mflux_fc(i) = outletwth(i) * mflux_dn
                  ELSE
                     hflux_fc(i) = outletwth(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                        + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                     mflux_fc(i) = outletwth(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                        + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                  ENDIF

                  sum_zgrad_riv(i) = sum_zgrad_riv(i) + outletwth(i) * 0.5*grav * height_up**2

                  zgrad_dn(i) = outletwth(i) * 0.5*grav * height_dn**2

               ELSEIF (ucat_next(i) == -99) THEN
                  ! downstream is not in model region.
                  ! assume: 1. downstream river bed is equal to this river bed.
                  !         2. downstream water surface is equal to this river depth.
                  !         3. downstream water velocity is equal to this velocity.

                  veloc_riv(i) = max(veloc_riv(i), 0.)

                  IF (wdsrf_ucat(i) > topo_rivhgt(i)) THEN

                     ! reconstruction of height of water near interface
                     height_up = wdsrf_ucat (i)
                     height_dn = topo_rivhgt(i)

                     veloct_fc = veloc_riv(i) + sqrt(grav * height_up) - sqrt(grav * height_dn)
                     height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn))) ** 2

                     vwave_up = min(veloc_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                     vwave_dn = max(veloc_riv(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))

                     hflux_up = veloc_riv(i) * height_up
                     hflux_dn = veloc_riv(i) * height_dn
                     mflux_up = veloc_riv(i)**2 * height_up + 0.5*grav * height_up**2
                     mflux_dn = veloc_riv(i)**2 * height_dn + 0.5*grav * height_dn**2

                     IF (vwave_up >= 0.) THEN
                        hflux_fc(i) = outletwth(i) * hflux_up
                        mflux_fc(i) = outletwth(i) * mflux_up
                     ELSEIF (vwave_dn <= 0.) THEN
                        hflux_fc(i) = outletwth(i) * hflux_dn
                        mflux_fc(i) = outletwth(i) * mflux_dn
                     ELSE
                        hflux_fc(i) = outletwth(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                           + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                        mflux_fc(i) = outletwth(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                           + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                     ENDIF

                     sum_zgrad_riv(i) = sum_zgrad_riv(i) + outletwth(i) * 0.5*grav * height_up**2

                  ELSE
                     hflux_fc(i) = 0
                     mflux_fc(i) = 0
                  ENDIF

               ELSEIF (ucat_next(i) == -10) THEN ! inland depression
                  hflux_fc(i) = 0
                  mflux_fc(i) = 0
               ENDIF

               sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
               sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

            ENDDO

            CALL worker_push_data (push_ups2ucat, upstream_flux_fields, mode = 'sum')

            IF (numucat > 0) THEN
               WHERE (ucatfilter)
                  sum_hflux_riv = sum_hflux_riv - hflux_sumups
                  sum_mflux_riv = sum_mflux_riv - mflux_sumups
                  sum_zgrad_riv = sum_zgrad_riv - zgrad_sumups
               END WHERE
            ENDIF

            ! reservoir operation.
            IF (DEF_Reservoir_Method > 0) THEN

               hflux_resv = 0._r8
               mflux_resv = 0._r8

               DO i = 1, numucat

                  IF (.not. ucatfilter(i)) CYCLE

                  IF (ucat_next(i) == -10) THEN
                     ! Inland depression: hflux_fc is fixed at zero, so nothing is
                     ! released, but the reservoir history must still report the
                     ! current inflow instead of the last value it was given.
                     IF (is_built_resv(i)) THEN
                        irsv = ucat2resv(i)
                        qresv_in(irsv)  = - sum_hflux_riv(i)
                        qresv_out(irsv) = 0._r8
                     ENDIF
                     CYCLE
                  ENDIF

                  IF (is_built_resv(i)) THEN

                     irsv = ucat2resv(i)
                     qresv_in(irsv) = - sum_hflux_riv(i)

                     IF (volresv(irsv) > 1.e-4 * volresv_total(irsv)) THEN
                        CALL reservoir_operation (DEF_Reservoir_Method, &
                           irsv, qresv_in(irsv), volresv(irsv), qresv_out(irsv))
                     ELSE
                        qresv_out (irsv) = 0.
                     ENDIF

                     hflux_fc(i) = qresv_out(irsv)
                     mflux_fc(i) = qresv_out(irsv) * sqrt(2*grav*wdsrf_ucat(i))

                     sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
                     sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

                     hflux_resv(i) = hflux_fc(i)
                     mflux_resv(i) = mflux_fc(i)
                  ENDIF

               ENDDO

               CALL worker_push_data (push_ups2ucat, reservoir_flux_fields, mode = 'sum')

               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_riv - hflux_sumups
                     sum_mflux_riv = sum_mflux_riv - mflux_sumups
                  END WHERE
               ENDIF

            ENDIF

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

               dt_this = dt_all(irivsys(i))

               ! constraint 1: CFL condition (only for rivers)
               IF (.not. is_built_resv(i)) THEN
                  IF ((veloc_riv(i) /= 0._r8) .or. (wdsrf_ucat(i) > 0._r8)) THEN
                     dt_this = min(dt_this, topo_rivlen(i) / &
                        (abs(veloc_riv(i)) + sqrt(grav * wdsrf_ucat(i))) * 0.8_r8)
                  ENDIF
               ENDIF

               ! constraint 2: avoid negative visible/reservoir water storage
               IF (sum_hflux_riv(i) > 0._r8) THEN
                  IF (.not. is_built_resv(i)) THEN
                     volwater = volwater_ucat(i)
                  ELSE
                     volwater = volresv(ucat2resv(i))
                  ENDIF
                  ! Tiny storage is protected by the gross outgoing-face
                  ! limiter below, rather than forcing a pathological dt.
                  IF (volwater > 1.e-6_r8) &
                     dt_this = min(dt_this, volwater / sum_hflux_riv(i))
               ENDIF

               ! constraint 3: avoid change of flow direction (only for rivers).
               ! Optional: DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT = .false. leaves the
               ! momentum update to the semi-implicit friction term, as CaMa-Flood does.
               IF (DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT) THEN
                  IF (.not. is_built_resv(i)) THEN
                     IF ((abs(veloc_riv(i)) > 0.1_r8) &
                        .and. (veloc_riv(i) * (sum_mflux_riv(i)-sum_zgrad_riv(i)) > 0._r8)) THEN
                        dt_this = min(dt_this, &
                           abs(momen_riv(i) * topo_rivare(i) / (sum_mflux_riv(i)-sum_zgrad_riv(i))))
                     ENDIF
                  ENDIF
               ENDIF

               dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))

            ENDDO

            ! Bifurcation fluxes are applied once below.  The old predictive
            ! dt-feedback loop is intentionally gone: ordinary routing still
            ! uses CFL/storage/momentum adaptive dt; BIF uses storage limiters.

            ! A net-flux dt bound cannot protect a nearly dry donor with two
            ! outgoing faces: upstream inflow can mask one outgoing face, and
            ! the old tiny-storage exemption skipped the bound altogether.
            ! Limit each donor's GROSS ordinary outflow to its current storage.
            ! Positive faces belong to this cell; negative faces belong to the
            ! downstream cell, so exchange both gross reverse outflow and the
            ! resulting donor scale before changing either shared face.
            normal_outgoing_rate = 0._r8
            hflux_sumups = 0._r8
            DO i = 1, numucat
               IF (.not. ucatfilter(i)) CYCLE
               IF (hflux_fc(i) >= 0._r8) THEN
                  normal_outgoing_rate(i) = hflux_fc(i)
               ELSE
                  hflux_sumups(i) = -hflux_fc(i)
               ENDIF
            ENDDO
            CALL worker_push_data(push_ups2ucat, hflux_sumups, mflux_sumups, fillvalue = 0._r8, mode = 'sum')
            WHERE (ucatfilter)
               normal_outgoing_rate = normal_outgoing_rate + mflux_sumups
            END WHERE
            ordinary_scale = 1._r8
            DO i = 1, numucat
               IF (.not. ucatfilter(i)) CYCLE
               IF (normal_outgoing_rate(i) <= 0._r8) CYCLE
               IF (is_built_resv(i)) THEN
                  volwater = volresv(ucat2resv(i))
               ELSE
                  volwater = volwater_ucat(i)
               ENDIF
               ! Final synchronization can only shorten this provisional dt,
               ! so capping against it is conservative even across workers.
               ! A pathological provisional dt is handled by the synchronizer;
               ! transfer no ordinary water from that donor in the meantime.
               IF (.not. ieee_is_finite(dt_all(irivsys(i)))) THEN
                  ordinary_scale(i) = 0._r8
               ELSEIF (dt_all(irivsys(i)) <= 0._r8) THEN
                  ordinary_scale(i) = 0._r8
               ELSE
                  ordinary_scale(i) = min(1._r8, max(volwater, 0._r8) / &
                     (normal_outgoing_rate(i) * dt_all(irivsys(i))))
               ENDIF
               normal_outgoing_rate(i) = normal_outgoing_rate(i) * ordinary_scale(i)
            ENDDO
            CALL worker_push_data(push_next2ucat, ordinary_scale, ordinary_scale_next, fillvalue = 1._r8)
            DO i = 1, numucat
               IF (.not. ucatfilter(i)) CYCLE
               IF (hflux_fc(i) >= 0._r8) THEN
                  dt_this = ordinary_scale(i)
               ELSE
                  dt_this = ordinary_scale_next(i)
               ENDIF
               hflux_fc(i) = hflux_fc(i) * dt_this
               mflux_fc(i) = mflux_fc(i) * dt_this
               IF (is_built_resv(i)) qresv_out(ucat2resv(i)) = hflux_fc(i)
            ENDDO
            ! Rebuild net fluxes from the capped faces.  Sender and receiver
            ! then use one identical volume/momentum transfer, including MPI
            ! boundaries, reverse flow, reservoirs and the tracer substep.
            CALL worker_push_data(push_ups2ucat, upstream_flux_fields, mode = 'sum')
            WHERE (ucatfilter)
               sum_hflux_riv = hflux_fc - hflux_sumups
               sum_mflux_riv = mflux_fc - mflux_sumups
            END WHERE
            DO i = 1, numucat
               IF (.not. ucatfilter(i)) CYCLE
               IF (is_built_resv(i)) qresv_in(ucat2resv(i)) = hflux_sumups(i)
            ENDDO
            ! The face limiter also scales momentum fluxes.  Recheck the
            ! optional non-reversal bound against their FINAL net force: an
            ! upstream donor may have been capped more than this cell.
            IF (DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT) THEN
               DO i = 1, numucat
                  IF (.not. ucatfilter(i) .or. is_built_resv(i)) CYCLE
                  IF ((abs(veloc_riv(i)) > 0.1_r8) .and. &
                     (veloc_riv(i) * (sum_mflux_riv(i)-sum_zgrad_riv(i)) > 0._r8)) THEN
                     dt_this = abs(momen_riv(i) * topo_rivare(i) / &
                        (sum_mflux_riv(i)-sum_zgrad_riv(i)))
                     ! Never let the synchronizer's invalid-dt fallback grow
                     ! a step after faces were capped against a smaller one.
                     IF (.not. ieee_is_finite(dt_this)) &
                        CALL CoLM_stop('grid_riverlake_flow: invalid final momentum dt')
                     IF (dt_this <= 0._r8) &
                        CALL CoLM_stop('grid_riverlake_flow: non-positive final momentum dt')
                     dt_all(irivsys(i)) = min(dt_all(irivsys(i)), dt_this)
                  ENDIF
               ENDDO
            ENDIF
            IF (DEF_USE_BIFURCATION) THEN
               ! One global step keeps cross-system BIF transfers paired.
               CALL sync_global_routing_dt(dt_res, dt_all, next_loop_active)
#ifdef USEMPI
            ELSE IF (rivsys_by_multiple_procs) THEN
               CALL mpi_allreduce (MPI_IN_PLACE, dt_all, 1, MPI_REAL8, MPI_MIN, &
                  p_comm_rivsys, p_err)
#endif
            ENDIF
            IF (DEF_USE_BIFURCATION) THEN
               IF (numucat > 0) sum_hflux_base = sum_hflux_riv
               ! Production BIF transport: storage limiters inside
               ! bifurcation_calc prevent donor overdraft without a predictive
               ! dt-feedback loop.
               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_base
                  END WHERE
               ENDIF

                     IF (allocated(volresv)) volresv_safe = volresv
                     CALL bifurcation_calc(wdsrf_ucat, wdsrf_ucat_prev, &
                              volwater_ucat, volwater_ucat_valid, &
                              volresv_safe, is_built_resv, dt_all, &
                        irivsys, ucatfilter, normal_outgoing_rate)

                     ! CaMa saves D2RIVDPH_PRE after calculating pathway flow
                     ! and before advancing storage/stage. Preserve that
                     ! ordering so the next adaptive substep sees the depth
                     ! paired with the carried pathway momentum.
                     wdsrf_ucat_prev = wdsrf_ucat
                     wdsrf_ucat_prev_valid = .true.

               IF (allocated(a_bifflw_lev) .and. allocated(a_bifflw_acctime)) THEN
                  DO ipth = 1, npthout_local
                     i_up = pth_upst_local(ipth)
                     IF (i_up < 1 .or. i_up > numucat) CYCLE
                     IF (.not. ucatfilter(i_up)) CYCLE
                     IF (allocated(bif_path_active)) THEN
                        IF (ipth <= size(bif_path_active)) THEN
                           IF (.not. bif_path_active(ipth)) CYCLE
                        ENDIF
                     ENDIF
                     a_bifflw_lev(:, ipth) = a_bifflw_lev(:, ipth) &
                        + bif_hflux_lev(:, ipth) * dt_all(irivsys(i_up))
                     a_bifflw_acctime(ipth) = a_bifflw_acctime(ipth) + dt_all(irivsys(i_up))
                  ENDDO
               ENDIF

               ! Add final bifurcation volume flux to main routing (volume
               ! only, not momentum).
               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_base + bif_hflux_sum
                  END WHERE
               ENDIF
#ifdef CoLMDEBUG
               ! Bifurcation conservation: the GLOBAL sum of bif_hflux_sum
               ! across all workers should be ~0 each step (what leaves
               ! one cell enters another). Track cumulative and worst-case.
               IF (numucat > 0) THEN
                  dt_this = sum(bif_hflux_sum)
               ELSE
                  dt_this = 0._r8
               ENDIF
#ifdef USEMPI
               CALL mpi_allreduce (MPI_IN_PLACE, dt_this, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
               bif_flux_sum_total = bif_flux_sum_total + dt_this
               bif_flux_sum_max = max(bif_flux_sum_max, abs(dt_this))
#endif
            ENDIF

            ! Per-sub-step tracer transport: advance tracer in lockstep
            ! with water BEFORE the water state update so concentration
            ! is computed from the pre-update volume.
#ifdef TRACER
            IF (ntracers > 0 .and. DEF_USE_BIFURCATION) &
               CALL tracer_lifecycle_route_sediment_bif_accumulate(dt_all, irivsys, ucatfilter, bif_hflux_lev)
            IF (allocated(volresv)) volresv_safe = volresv
            IF (DEF_USE_BIFURCATION) THEN
               CALL tracer_substep (acctime_rnof, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
                  wdsrf_ucat, ucatfilter, volresv_safe, ucat2resv_safe, is_built_resv, &
                  do_bif = .true., bif_hflux_lev_in = bif_hflux_lev, npthout_local_in = npthout_local)
            ELSE
               CALL tracer_substep (acctime_rnof, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
                  wdsrf_ucat, ucatfilter, volresv_safe, ucat2resv_safe, is_built_resv)
            ENDIF
#endif

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

                  IF (.not. is_built_resv(i)) THEN
                     volwater = volwater_ucat(i)
                  ELSE
                  volwater = volresv(ucat2resv(i))
               ENDIF

#ifdef TRACER
               ! The particle provider needs the carrier before this HYDRO
               ! substep, including protected levee water when present.
               particle_water_storage_start(i) = max(volwater, 0._r8)
               particle_protected_start(i) = 0._r8
               IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) &
                  particle_water_storage_start(i) = particle_water_storage_start(i) + max(levsto(i), 0._r8)
               IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) &
                  particle_protected_start(i) = max(levsto(i), 0._r8)
#endif

                  visible_hflux = sum_hflux_riv(i)
                  protected_hflux = 0._r8
                  IF (DEF_USE_LEVEE .and. DEF_USE_BIFURCATION .and. allocated(bif_lev_hflux_sum)) THEN
                     IF ((.not. is_built_resv(i)) .and. has_levee(i) .and. i <= size(bif_lev_hflux_sum)) THEN
                        protected_hflux = bif_lev_hflux_sum(i)
                        visible_hflux = sum_hflux_riv(i) - protected_hflux
                     ENDIF
                  ENDIF
                  volwater = volwater - visible_hflux * dt_all(irivsys(i))
                  IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                        CALL levee_apply_protected_flux(i, protected_hflux, &
                           dt_all(irivsys(i)), protected_clip)
                        IF (protected_clip > 0._r8) THEN
                           write(*,'(A,I0,A,ES12.4)') &
                              'ERROR bifurcation protected limiter failed: ucat=', i, &
                              ' clipped=', protected_clip
                           CALL CoLM_stop('BIF protected-side limiter failed')
                        ENDIF
#ifdef CoLMDEBUG
                        bif_protected_clip_sum = bif_protected_clip_sum + protected_clip
#endif
                     ENDIF
#ifdef CoLMDEBUG
               IF (volwater < 0._r8) totalclip = totalclip - volwater
#endif
               ! Water created by this clip is not silent in any build: it enters
               ! totalvol_aft, so it appears in the always-on closure check
               ! (water_balance_err vs water_balance_tol, the analogue of CaMa's
               ! CALC_WATBAL) at the end of the routing period.  CoLMDEBUG only
               ! breaks the amount out as totalclip.
               volwater = max(volwater, 0.)

               ! Inland depression overflow is a post-transport water correction.
               ! ordinary routing tracer transport is handled earlier in tracer_substep.
               ! only this explicit correction updates tracer locally so the
               ! removed tracer follows the same sink.
               IF (ucat_next(i) == -10) THEN
                  IF (volwater > topo_rivstomax(i)) THEN
                     ! Remove excess water after transport has already run.
                     hflux_fc(i) = (volwater - topo_rivstomax(i)) / dt_all(irivsys(i))
                     IF (is_built_resv(i)) qresv_out(ucat2resv(i)) = hflux_fc(i)
                     ! Remove corresponding tracer proportionally.
                     ! Update trc_flux_out so the unified discharge diagnostic
                     ! at line ~895 (ucat_next <= 0) picks up the correct value.
                     !
#ifdef TRACER
                        IF (volwater > 1.e-6_r8) THEN
                              frac_remove = (volwater - topo_rivstomax(i)) / volwater
                              IF (allocated(volresv)) volresv_safe = volresv
                                 CALL get_cell_volume_dep(i, &
                                    floodplain_curve(i)%depth(topo_rivstomax(i)), &
                                    volresv_safe, ucat2resv_safe, vol_post)
                           vol_post = max(vol_post, 1.e-6_r8)
                        DO itrc_dep = 1, ntracers
                           IF (.not. tracer_uses_land_water_transport(itrc_dep)) CYCLE
                           IF (allocated(trc_solid)) THEN
                              IF (tracer_has_dissolved_limit(itrc_dep)) &
                                 CALL tracer_equilibrate_dissolved(itrc_dep, volwater, &
                                    trc_mass(itrc_dep, i), trc_solid(itrc_dep, i))
                           ENDIF
                           trc_removed = trc_mass(itrc_dep, i) * frac_remove
                           trc_mass(itrc_dep, i) = trc_mass(itrc_dep, i) - trc_removed
                           trc_flux_out(itrc_dep, i) = trc_removed / dt_all(irivsys(i))
                           trc_conc_dep(itrc_dep, i) = trc_mass(itrc_dep, i) / vol_post
                        ENDDO
                     END IF
#endif
                     volwater = topo_rivstomax(i)
                  ENDIF
               ENDIF

               IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                  ! CaMa's simplified levee scheme applies all pathway fluxes
                  ! to total storage, then restores the static visible/protected
                  ! partition.  The split pools above are transport/limiter
                  ! bookkeeping only; retaining that transient split here made
                  ! river stage inconsistent with visible storage.
                  CALL levee_repartition_storage(i, volwater, wdsrf_ucat(i), &
                     fldfrc_levee, vis_vol_bef_lv2, levsto_bef_lv2)
                  volwater_ucat(i) = volwater
                  levee_floodarea(i) = fldfrc_levee * topo_area(i)
#ifdef TRACER
                     IF (ntracers > 0) CALL tracer_lifecycle_route_sediment_levee_repartition(i, &
                        vis_vol_bef_lv2, levsto_bef_lv2, volwater_ucat(i), levsto(i))
                     CALL levee_tracer_repartition(i, &
                        vis_vol_bef_lv2, levsto_bef_lv2, &
                        volwater_ucat(i), levsto(i))
#endif
                  ELSE
                     volwater_ucat(i) = volwater
                     wdsrf_ucat(i) = floodplain_curve(i)%depth (volwater)
                     IF (DEF_USE_LEVEE) levee_floodarea(i) = 0.
               ENDIF

               IF (is_built_resv(i)) THEN
                  volresv(ucat2resv(i)) = volwater
               ENDIF

               IF ((.not. is_built_resv(i)) .and. (wdsrf_ucat(i) >= RIVERMIN)) THEN
                  friction = grav * topo_rivman(i)**2 / wdsrf_ucat(i)**(7.0/3.0) * abs(momen_riv(i))
                  momen_riv(i) = (momen_riv(i) &
                     - (sum_mflux_riv(i) - sum_zgrad_riv(i)) / topo_rivare(i) * dt_all(irivsys(i))) &
                     / (1 + friction * dt_all(irivsys(i)))
                  veloc_riv(i) = momen_riv(i) / wdsrf_ucat(i)
               ELSE
                  momen_riv(i) = 0
                  veloc_riv(i) = 0
               ENDIF

               ! inland depression river
               IF ((.not. is_built_resv(i)) .and. (ucat_next(i) == -10)) THEN
                  momen_riv(i) = min(0., momen_riv(i))
                  veloc_riv(i) = min(0., veloc_riv(i))
               ENDIF

               veloc_riv(i) = min(veloc_riv(i),  20.)
               veloc_riv(i) = max(veloc_riv(i), -20.)
               ! Keep momen_riv consistent with the clamped velocity; next
               ! substep reads friction from abs(momen_riv) and flux from
               ! veloc_riv, so they must describe the same physical state.
               IF ((.not. is_built_resv(i)) .and. (wdsrf_ucat(i) >= RIVERMIN)) THEN
                  momen_riv(i) = veloc_riv(i) * wdsrf_ucat(i)
               ENDIF

            ENDDO

#ifdef TRACER
                  IF (p_is_worker .and. numucat > 0) THEN
                     IF (allocated(volresv)) volresv_safe = volresv
                     CALL tracer_diag_accumulate_substep (dt_all, irivsys, ucatfilter, wdsrf_ucat, &
                        volresv_safe, ucat2resv_safe, is_built_resv)
                  END IF
#endif

            DO i = 1, numucat
               IF (ucatfilter(i)) THEN

                  IF (ucat_next(i) <= 0) THEN
                     totaldis = totaldis + hflux_fc(i)*dt_all(irivsys(i))
#ifdef CoLMDEBUG
                     ! Accumulate tracer discharge at river mouth
#ifdef TRACER
                     IF (ntracers > 0) THEN
                        DO itrc = 1, ntracers
                           IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
                           trc_mass_dis(itrc) = trc_mass_dis(itrc) &
                              + trc_flux_out(itrc,i)*dt_all(irivsys(i))
                        ENDDO
                     END IF
#endif
#endif
                  ENDIF

                  acctime_ucat(i) = acctime_ucat(i) + dt_all(irivsys(i))

                  a_wdsrf_ucat(i) = a_wdsrf_ucat(i) + wdsrf_ucat(i) * dt_all(irivsys(i))
                  a_veloc_riv (i) = a_veloc_riv (i) + veloc_riv (i) * dt_all(irivsys(i))
                  a_discharge (i) = a_discharge (i) + hflux_fc  (i) * dt_all(irivsys(i))

                  IF (DEF_USE_LEVEE .and. levee_floodarea(i) > 0.) THEN
                     floodarea = levee_floodarea(i)
                  ELSE
                     floodarea = floodplain_curve(i)%floodarea (wdsrf_ucat(i))
                  ENDIF
                  a_floodarea (i) = a_floodarea (i) + floodarea * dt_all(irivsys(i))
                  total_floodarea (i) = floodarea
                  IF (DEF_USE_LEVEE .and. levee_floodarea(i) > 0._r8) THEN
                     total_flooddepth(i) = max(levdph(i), max(wdsrf_ucat(i) - floodplain_curve(i)%rivhgt, 0._r8))
                  ELSE
                     total_flooddepth(i) = max(wdsrf_ucat(i) - floodplain_curve(i)%rivhgt, 0._r8)
                  ENDIF

                  ! River/floodplain storage separation, total storage, surface elevation
                  ! Partition the actual visible volume at bankfull capacity.
                  ! Above bankfull, the volume-depth curve need not contain
                  ! rivare*wdsrf; using that rectangle can exceed total storage.
                  ! rivsto = below-bank storage; fldsto = visible overbank storage
                  ! flddph = max(wdsrf - rivhgt, 0) (depth above channel banks)
                  ! storge = total_volume (+ levsto if levee enabled)
                  ! sfcelv = bed_elevation + wdsrf (matches CaMa: D2RIVELV + D2RIVDPH)
                     IF (is_built_resv(i)) THEN
                        volwater = volresv(ucat2resv(i))
                     ELSE
                        volwater = volwater_ucat(i)
                     ENDIF
                  rivsto_hist = min(volwater, floodplain_curve(i)%rivstomax)
                  a_rivsto(i) = a_rivsto(i) + rivsto_hist * dt_all(irivsys(i))
                  a_fldsto(i) = a_fldsto(i) &
                     + (volwater - rivsto_hist) * dt_all(irivsys(i))
                  a_flddph(i) = a_flddph(i) &
                     + max(wdsrf_ucat(i) - floodplain_curve(i)%rivhgt, 0._r8) * dt_all(irivsys(i))
                  IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                     a_storge(i) = a_storge(i) + (volwater + levsto(i)) * dt_all(irivsys(i))
                  ELSE
                     a_storge(i) = a_storge(i) + volwater * dt_all(irivsys(i))
                  ENDIF
                  a_sfcelv(i) = a_sfcelv(i) + (topo_rivelv(i) + wdsrf_ucat(i)) * dt_all(irivsys(i))

                  IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                     a_levsto(i) = a_levsto(i) + levsto(i) * dt_all(irivsys(i))
                     a_levdph(i) = a_levdph(i) + levdph(i) * dt_all(irivsys(i))
                  ENDIF

                  IF (DEF_USE_BIFURCATION) THEN
                     a_bifout(i) = a_bifout(i) + bif_hflux_sum(i) * dt_all(irivsys(i))
                  ENDIF

                  IF (is_built_resv(i)) THEN
                     irsv = ucat2resv(i)
                     acctime_resv(irsv) = acctime_resv(irsv) + dt_all(irivsys(i))
                     a_volresv   (irsv) = a_volresv  (irsv) + volresv  (irsv) * dt_all(irivsys(i))
                     a_qresv_in  (irsv) = a_qresv_in (irsv) + qresv_in (irsv) * dt_all(irivsys(i))
                     a_qresv_out (irsv) = a_qresv_out(irsv) + qresv_out(irsv) * dt_all(irivsys(i))
                  ENDIF

               ENDIF
            ENDDO

            dt_res = dt_res - dt_all

#ifdef TRACER
            IF (tracer_lifecycle_route_has_active()) THEN
               DO i = 1, numucat
                  IF (ucatfilter(i)) THEN
                     particle_protected_end(i) = 0._r8
                     ! Same levee/floodplain flood area the history block above
                     ! stored for this substep; do not evaluate the curve again.
                     particle_floodarea(i) = total_floodarea(i)
                     particle_protected_area(i) = 0._r8
                     IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) &
                        particle_protected_area(i) = min(max(levee_floodarea(i) - &
                           levee_frc_data(i) * topo_area(i), 0._r8), &
                           (1._r8 - levee_frc_data(i)) * topo_area(i))
                     ! Sediment concentration uses the same full water volume
                     ! HYDRO just advanced, not a channel-width rectangle that
                     ! omits floodplain (and protected levee) storage.
                     IF (is_built_resv(i) .and. allocated(volresv) .and. allocated(ucat2resv)) THEN
                        irsv = ucat2resv(i)
                        IF (irsv >= 1 .and. irsv <= size(volresv) .and. volresv(irsv) /= spval) THEN
                           particle_water_storage(i) = max(volresv(irsv), 0._r8)
                        ELSE
                           particle_water_storage(i) = max(volwater_ucat(i), 0._r8)
                        ENDIF
                     ELSE
                        particle_water_storage(i) = max(volwater_ucat(i), 0._r8)
                        IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                           particle_water_storage(i) = particle_water_storage(i) + max(levsto(i), 0._r8)
                           particle_protected_end(i) = max(levsto(i), 0._r8)
                        ENDIF
                     ENDIF
                  ELSE
                     particle_floodarea(i) = 0._r8
                     particle_protected_area(i) = 0._r8
                     particle_water_storage(i) = 0._r8
                  ENDIF
               ENDDO
               CALL tracer_lifecycle_route_diag_accumulate(dt_all, irivsys, ucatfilter, &
                  veloc_riv, wdsrf_ucat, particle_water_storage_start, particle_water_storage, &
                  hflux_fc, particle_floodarea, particle_protected_start, particle_protected_end, &
                  particle_protected_area)
            ENDIF
#endif

            IF (DEF_USE_BIFURCATION) THEN
               ! sync_global_routing_dt derives this from the same global
               ! reduction that selected the just-completed substep.  Avoid a
               ! second collective solely to determine loop continuation.
               loop_active = next_loop_active
            ELSE
               loop_active = any(dt_res > 0)
            ENDIF

         ENDDO

               ! Keep restart-visible state consistent after the substep loop.
               volwater_ucat_valid = .true.

         totalvol_aft = 0._r8
         DO i = 1, numucat
               IF (.not. is_built_resv(i)) THEN
                  IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                     volwater = volwater_ucat(i) + levsto(i)
                  ELSE
                     volwater = volwater_ucat(i)
                  ENDIF
               totalvol_aft = totalvol_aft + volwater
            ELSE
               volwater = volresv(ucat2resv(i))
               totalvol_aft = totalvol_aft + volwater
            ENDIF
         ENDDO
#ifdef CoLMDEBUG
         ! Tracer conservation: compute total mass after routing.
         ! Include protected-side pool so the levee repartition
         ! stays internally closed in the before/after accounting.
#ifdef TRACER
         IF (numucat > 0) THEN
            DO itrc = 1, ntracers
               IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
               trc_mass_aft(itrc) = sum(trc_mass(itrc,:)) + sum(trc_inp_buf(itrc,:))
               IF (allocated(trc_levsto)) &
                  trc_mass_aft(itrc) = trc_mass_aft(itrc) + sum(trc_levsto(itrc,:))
               IF (allocated(trc_solid)) trc_mass_aft(itrc) = trc_mass_aft(itrc) &
                  + sum(trc_solid(itrc,:)) + sum(trc_levsto_solid(itrc,:))
            ENDDO
         ELSE
            trc_mass_aft = 0._r8
         END IF
#else
         trc_mass_aft = 0._r8
#endif
#endif
      ENDIF

#ifdef USEMPI
      water_balance_vec = (/ totalvol_bef, totalvol_aft, totalrnof, totaldis, &
         flood_evap_period, flood_infil_period /)
      IF (.not. p_is_worker) water_balance_vec = 0._r8
      IF (DEF_GridRiverLake_FloodFeedback) THEN
         CALL mpi_allreduce (MPI_IN_PLACE, water_balance_vec, 6, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      ELSE
         CALL mpi_allreduce (MPI_IN_PLACE, water_balance_vec, 4, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      ENDIF
      totalvol_bef = water_balance_vec(1)
      totalvol_aft = water_balance_vec(2)
      totalrnof    = water_balance_vec(3)
      totaldis     = water_balance_vec(4)
#endif
      IF (p_is_master .and. (water_balance_vec(5) > 0._r8 .or. water_balance_vec(6) > 0._r8)) &
         write(*,'(A,2(1X,ES24.16))') 'Grid flood feedback evap/infil [m3]:', &
            water_balance_vec(5), water_balance_vec(6)
      IF (p_is_master .and. DEF_GridRiverLake_FloodFeedback) &
         write(*,'(A,4(1X,ES24.16))') 'Grid flood route S0/S1/R/Q [m3]:', &
            totalvol_bef, totalvol_aft, totalrnof, totaldis
      flood_evap_period = 0._r8
      flood_infil_period = 0._r8

      water_balance_err = totalvol_aft - totalvol_bef - totalrnof + totaldis
      water_balance_tol = max(1.e-6_r8, 1.e-10_r8 * max(abs(totalvol_bef), abs(totalvol_aft), &
         abs(totalrnof), abs(totaldis), 1._r8))
      IF (p_is_master .and. (water_balance_err /= water_balance_err .or. &
          abs(water_balance_err) > water_balance_tol)) THEN
         IF (routing_mass_warn_count < 5) THEN
            write(*,'(A,ES12.4,A,ES12.4,A)') &
               'WARNING grid_riverlake_flow: water balance residual=', water_balance_err, &
               ' m3 exceeds tolerance=', water_balance_tol, ' m3'
         ENDIF
         routing_mass_warn_count = routing_mass_warn_count + 1
      ENDIF

#ifdef CoLMDEBUG
#ifdef USEMPI
      IF (.not. p_is_worker) ntimestep = 0
      CALL mpi_allreduce (MPI_IN_PLACE, ntimestep, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)

      IF (.not. p_is_worker) totalclip = 0._r8
      IF (.not. p_is_worker) bif_protected_clip_sum = 0._r8

      CALL mpi_allreduce (MPI_IN_PLACE, totalclip, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, bif_protected_clip_sum, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
#ifdef TRACER
      IF (p_is_worker .and. numucat > 0) THEN
         DO itrc = 1, ntracers
            IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
            trc_mass_dis(itrc) = trc_mass_dis(itrc) + sum(trc_dry_drain(itrc,:))
            IF (allocated(trc_reactive_source)) THEN
               trc_mass_reactive(itrc) = sum(trc_reactive_source(itrc,:))
            ENDIF
         ENDDO
      END IF
#endif
         IF (.not. p_is_worker) THEN
            bif_flux_sum_total = 0._r8
            bif_flux_sum_max   = 0._r8
            bif_protected_clip_sum = 0._r8
            trc_mass_bef = 0._r8
         trc_mass_aft = 0._r8
         trc_mass_inp = 0._r8
         trc_mass_dis = 0._r8
         trc_mass_reactive = 0._r8
      ENDIF
         CALL mpi_allreduce (MPI_IN_PLACE, bif_flux_sum_total, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, bif_flux_sum_max,   1, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
#ifdef TRACER
      IF (ntracers > 0) THEN
         CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_bef, ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_aft, ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_inp, ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_dis, ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_reactive, ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      END IF
#endif
#endif

#ifdef CoLMDEBUG
      IF (p_is_master) THEN
         write(*,'(/,A)') 'Checking River Routing Flow ...'
         write(*,'(A,F12.5,A)') 'River Lake Flow minimum average timestep: ', acctime_rnof/ntimestep, ' seconds'
         write(*,'(A,ES8.1,A)') 'Total water before :  ', totalvol_bef,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total runoff :        ', totalrnof, ' m^3'
         write(*,'(A,ES8.1,A)') 'Total discharge :     ', totaldis,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total water change :  ', totalvol_aft-totalvol_bef,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total water balance : ', totalvol_aft-totalvol_bef-totalrnof+totaldis,  ' m^3'
         write(*,'(A,ES8.1,A)') 'Total water clipping :', totalclip, ' m^3'
            IF (DEF_USE_BIFURCATION) THEN
               write(*,'(A,ES10.2,A)') 'Bif global net flux (cumul)     : ', bif_flux_sum_total, ' m^3/s (should be ~0)'
               write(*,'(A,ES10.2,A)') 'Bif global net flux (max/step)  : ', bif_flux_sum_max,   ' m^3/s (should be ~0)'
               write(*,'(A,ES10.2,A)') 'Bif protected-side limiter residual clip : ', bif_protected_clip_sum, ' m^3'
            ENDIF
#ifdef TRACER
            DO itrc_dbg = 1, ntracers
               IF (.not. tracer_uses_land_water_transport(itrc_dbg)) CYCLE
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass before  : ', &
                  trc_mass_bef(itrc_dbg), ' R*m3'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass input   : ', &
                  trc_mass_inp(itrc_dbg), ' R*m3'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass discharge: ', &
                  trc_mass_dis(itrc_dbg), ' R*m3'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass reactive : ', &
                  trc_mass_reactive(itrc_dbg), ' R*m3'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass after   : ', &
                  trc_mass_aft(itrc_dbg), ' R*m3'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass change  : ', &
                  trc_mass_aft(itrc_dbg) - trc_mass_bef(itrc_dbg), ' R*m3'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass balance : ', &
                  trc_mass_aft(itrc_dbg) - trc_mass_bef(itrc_dbg) &
                  - trc_mass_inp(itrc_dbg) + trc_mass_dis(itrc_dbg) &
                  - trc_mass_reactive(itrc_dbg), ' R*m3 (should be ~0)'
            ENDDO
#endif
      ENDIF  ! p_is_master

#endif

#ifdef TRACER
      ! Donor-limiter work for the routing period just finished.  Every worker
      ! must take part: each rank reads and clears its own counters, then the
      ! totals are reduced (SUM for work, MAX for the peak) so the report covers
      ! the whole domain instead of one rank's share.  Printed for the first few
      ! periods (the repo's one-shot convention) and whenever some rank needed
      ! the long-chain fast path.
      IF (p_is_worker) THEN
         CALL tracer_limiter_stats (lim_calls, lim_iter_sum, lim_iter_peak, lim_over_soft, &
            reset = .true.)
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, lim_calls, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_iter_sum, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_over_soft, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_iter_peak, 1, MPI_INTEGER, MPI_MAX, p_comm_worker, p_err)
#endif
         IF (p_iam_worker == p_root) THEN
            IF (lim_calls > 0 .and. (lim_diag_printed < 5 .or. lim_over_soft > 0)) THEN
               IF (lim_diag_printed < 5) lim_diag_printed = lim_diag_printed + 1
               write(*,'(A,I0,A,F9.3,A,I0,A,I0)') 'River tracer donor limiter: calls=', lim_calls, &
                  ' mean iter=', real(lim_iter_sum, r8) / real(lim_calls, r8), &
                  ' peak=', lim_iter_peak, ' over_soft=', lim_over_soft
            ENDIF
         ENDIF
      ENDIF
#endif

#ifdef TRACER
      IF (tracer_lifecycle_route_has_active() .and. p_is_worker) THEN
         ! All workers must participate (MPI point-to-point inside push_data).
         ! Particle tracers compute their own flood-exposure diagnostics from
         ! per-routing-period accumulators, not history-period averages.
         CALL tracer_lifecycle_route_calc(acctime_rnof)
      ENDIF
#endif

#ifdef TRACER
      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            ! With TRACER compiled in but no tracer configured
            ! (DEF_TRACER_NUM = 0) river_lake_tracer_init returns early and
            ! leaves these arrays unallocated.
            IF (allocated(acc_trc_inp)) acc_trc_inp = 0._r8
            IF (allocated(acc_rnof_ref)) acc_rnof_ref = 0._r8
            IF (allocated(trc_dry_drain)) trc_dry_drain = 0._r8
         ENDIF
      END IF
#endif

      acctime_rnof = 0.

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            acc_rnof_uc = 0.
         ENDIF
      ENDIF

            ! ---- Publish per-ucat levee/general flood diagnostics to
         !      reactive tracers.  Because this publish occurs at the end
         !      of routing, reactive land tracers normally consume it on
         !      the next land step.
         !      Flow: ucat -> inpm grid (average) -> landpatch (remap).
         !      Inactive when GridRiverLakeFlow is undef (this whole file
         !      is gated by that macro).
#ifdef TRACER
         IF (allocated(levee_floodarea)) THEN
            CALL publish_levee_fldfrc_to_patches (levee_floodarea)
         ENDIF
         IF (allocated(total_floodarea)) THEN
            CALL publish_fldfrc_to_patches (total_floodarea, total_flooddepth)
         ENDIF
#endif
         IF (DEF_GridRiverLake_FloodFeedback) CALL publish_flood_feedback(year)

      IF (allocated(is_built_resv)) deallocate(is_built_resv)
      IF (allocated(wdsrf_next   )) deallocate(wdsrf_next   )
      IF (allocated(veloc_next   )) deallocate(veloc_next   )
      IF (allocated(hflux_fc     )) deallocate(hflux_fc     )
      IF (allocated(mflux_fc     )) deallocate(mflux_fc     )
      IF (allocated(zgrad_dn     )) deallocate(zgrad_dn     )
      IF (allocated(hflux_resv   )) deallocate(hflux_resv   )
      IF (allocated(mflux_resv   )) deallocate(mflux_resv   )
      IF (allocated(hflux_sumups )) deallocate(hflux_sumups )
      IF (allocated(mflux_sumups )) deallocate(mflux_sumups )
      IF (allocated(zgrad_sumups )) deallocate(zgrad_sumups )
      IF (allocated(sum_hflux_riv)) deallocate(sum_hflux_riv)
      IF (allocated(sum_hflux_base)) deallocate(sum_hflux_base)
      IF (allocated(normal_outgoing_rate)) deallocate(normal_outgoing_rate)
      IF (allocated(ordinary_scale)) deallocate(ordinary_scale)
      IF (allocated(ordinary_scale_next)) deallocate(ordinary_scale_next)
      IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
      IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)
      IF (allocated(ucatfilter      )) deallocate(ucatfilter      )
      IF (allocated(levee_floodarea)) deallocate(levee_floodarea)
         IF (allocated(total_floodarea)) deallocate(total_floodarea)
         IF (allocated(total_flooddepth)) deallocate(total_flooddepth)
#ifdef TRACER
         IF (allocated(particle_floodarea)) deallocate(particle_floodarea)
         IF (allocated(particle_protected_area)) deallocate(particle_protected_area)
         IF (allocated(particle_water_storage_start)) deallocate(particle_water_storage_start)
         IF (allocated(particle_water_storage)) deallocate(particle_water_storage)
         IF (allocated(particle_protected_start)) deallocate(particle_protected_start, particle_protected_end)
#endif
            IF (allocated(dt_res         )) deallocate(dt_res         )
            IF (allocated(dt_all       )) deallocate(dt_all       )
         IF (allocated(volresv_safe)) deallocate(volresv_safe)
         IF (allocated(ucat2resv_safe)) deallocate(ucat2resv_safe)
#ifdef CoLMDEBUG
      IF (allocated(trc_mass_bef)) deallocate(trc_mass_bef, trc_mass_aft, &
                                              trc_mass_inp, trc_mass_dis, trc_mass_reactive)
#endif

   END SUBROUTINE grid_riverlake_flow

#ifdef TRACER
      SUBROUTINE publish_levee_fldfrc_to_patches (levee_floodarea_in)
      !-------------------------------------------------------------------
      ! Push the per-ucat levee floodplain fraction through the inpm grid
      ! back to reactive tracers as a per-patch levee flood fraction.
      !-------------------------------------------------------------------
      USE MOD_Grid_RiverLakeNetwork, only: numucat, numinpm, topo_area, &
                                           push_ucat2inpm, remap_patch2inpm
      USE MOD_WorkerPushData, only: worker_push_data, worker_remap_data_grid2pset
      USE MOD_LandPatch, only: numpatch
      USE MOD_SPMD_Task

      real(r8), intent(in) :: levee_floodarea_in(:)
      real(r8), allocatable :: fldfrc_uc(:), fldfrc_gd(:), fldfrc_patch(:)
      integer :: i

            IF (.not. p_is_worker) RETURN
            IF (.not. tracer_lifecycle_has_levee_flood_publisher()) RETURN

      ! Build per-ucat fldfrc
      IF (numucat > 0) THEN
         allocate (fldfrc_uc(numucat))
         DO i = 1, numucat
            IF (topo_area(i) > 0._r8 .and. i <= size(levee_floodarea_in)) THEN
               fldfrc_uc(i) = min(1._r8, max(0._r8, &
                  levee_floodarea_in(i) / topo_area(i)))
            ELSE
               fldfrc_uc(i) = 0._r8
            ENDIF
         ENDDO
      ELSE
         allocate (fldfrc_uc(0))
      ENDIF

      ! ucat -> inpm grid (area-weighted average)
      IF (numinpm > 0) THEN
         allocate (fldfrc_gd(numinpm))
      ELSE
         allocate (fldfrc_gd(0))
      ENDIF
            CALL worker_push_data (push_ucat2inpm, fldfrc_uc, fldfrc_gd, &
               fillvalue = RIVERLAKE_FLOOD_MISSING_VALUE, mode = 'average')

            ! inpm grid -> landpatch (area-weighted remap)
            allocate(fldfrc_patch(max(0,numpatch)))
            CALL worker_remap_data_grid2pset (remap_patch2inpm, fldfrc_gd, &
               fldfrc_patch, fillvalue = RIVERLAKE_FLOOD_MISSING_VALUE, mode = 'average')
            WHERE (fldfrc_patch == RIVERLAKE_FLOOD_MISSING_VALUE) fldfrc_patch = 0._r8
            CALL tracer_lifecycle_publish_levee_flood_patch (fldfrc_patch)

         deallocate(fldfrc_uc, fldfrc_gd, fldfrc_patch)

   END SUBROUTINE publish_levee_fldfrc_to_patches

      SUBROUTINE publish_fldfrc_to_patches (total_floodarea_in, total_flooddepth_in)
      !-------------------------------------------------------------------
      ! Publish general flood area and depth to reactive tracers.
      !-------------------------------------------------------------------
      USE MOD_Grid_RiverLakeNetwork, only: numucat, numinpm, topo_area, &
                                           push_ucat2inpm, remap_patch2inpm
      USE MOD_WorkerPushData, only: worker_push_data, worker_remap_data_grid2pset
      USE MOD_LandPatch, only: numpatch
      USE MOD_SPMD_Task

      real(r8), intent(in) :: total_floodarea_in(:)
      real(r8), intent(in) :: total_flooddepth_in(:)
      real(r8), allocatable :: fldfrc_uc(:), fldfrc_gd(:)
      real(r8), allocatable :: fldwat_uc(:), fldwat_gd(:)
      real(r8), allocatable :: fldfrc_patch(:), fldwat_patch(:), flddph_patch(:)
      integer :: i

            IF (.not. p_is_worker) RETURN
            IF (.not. tracer_lifecycle_has_flood_publisher()) RETURN

      IF (numucat > 0) THEN
         allocate (fldfrc_uc(numucat), fldwat_uc(numucat))
         fldwat_uc(:) = 0._r8
         DO i = 1, numucat
            IF (topo_area(i) > 0._r8 .and. i <= size(total_floodarea_in)) THEN
               fldfrc_uc(i) = min(1._r8, max(0._r8, &
                  total_floodarea_in(i) / topo_area(i)))
               IF (i <= size(total_flooddepth_in)) THEN
                  ! Remap patch-mean floodwater depth f*d, not conditional
                  ! depth d independently from its flooded-area fraction f.
                  fldwat_uc(i) = fldfrc_uc(i) * max(0._r8, total_flooddepth_in(i))
               ENDIF
            ELSE
               fldfrc_uc(i) = 0._r8
               fldwat_uc(i) = 0._r8
            ENDIF
         ENDDO
      ELSE
         allocate (fldfrc_uc(0), fldwat_uc(0))
      ENDIF

      IF (numinpm > 0) THEN
         allocate (fldfrc_gd(numinpm), fldwat_gd(numinpm))
      ELSE
         allocate (fldfrc_gd(0), fldwat_gd(0))
      ENDIF
            CALL worker_push_data (push_ucat2inpm, fldfrc_uc, fldfrc_gd, &
               fillvalue = RIVERLAKE_FLOOD_MISSING_VALUE, mode = 'average')

            allocate(fldfrc_patch(max(0,numpatch)), fldwat_patch(max(0,numpatch)), &
               flddph_patch(max(0,numpatch)))
            CALL worker_remap_data_grid2pset (remap_patch2inpm, fldfrc_gd, &
               fldfrc_patch, fillvalue = RIVERLAKE_FLOOD_MISSING_VALUE, mode = 'average')
            WHERE (fldfrc_patch == RIVERLAKE_FLOOD_MISSING_VALUE) fldfrc_patch = 0._r8

            CALL worker_push_data (push_ucat2inpm, fldwat_uc, fldwat_gd, &
               fillvalue = RIVERLAKE_FLOOD_MISSING_VALUE, mode = 'average')
            CALL worker_remap_data_grid2pset (remap_patch2inpm, fldwat_gd, &
               fldwat_patch, fillvalue = RIVERLAKE_FLOOD_MISSING_VALUE, mode = 'average')
            WHERE (fldwat_patch == RIVERLAKE_FLOOD_MISSING_VALUE) fldwat_patch = 0._r8
            flddph_patch = 0._r8
            DO i = 1, numpatch
               IF (fldfrc_patch(i) > 0._r8) THEN
                  flddph_patch(i) = max(0._r8, fldwat_patch(i)) / fldfrc_patch(i)
               ENDIF
            ENDDO
            CALL tracer_lifecycle_publish_flood_patch (fldfrc_patch, flddph_patch)

         deallocate(fldfrc_uc, fldfrc_gd, fldwat_uc, fldwat_gd, &
            fldfrc_patch, fldwat_patch, flddph_patch)

   END SUBROUTINE publish_fldfrc_to_patches
#endif


   SUBROUTINE sync_global_routing_dt(dt_res, dt_all, next_loop_active)

      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

      real(r8), intent(in)    :: dt_res(:)
      real(r8), intent(inout) :: dt_all(:)
      logical,  intent(out)   :: next_loop_active

      real(r8) :: dt_reduce(2), dt_global, max_residual, remaining_after
      logical  :: local_pathological
      logical  :: pathological_mask(size(dt_res))
      integer  :: i

      pathological_mask = .false.
      DO i = 1, size(dt_res)
         ! Do not compare a NaN with zero: Fortran does not guarantee
         ! short-circuit evaluation and production builds may trap invalid FP.
         IF (.not. ieee_is_finite(dt_res(i))) CYCLE
         IF (dt_res(i) > 0._r8) THEN
            IF (.not. ieee_is_finite(dt_all(i))) THEN
               pathological_mask(i) = .true.
            ELSEIF (dt_all(i) <= 0._r8) THEN
               pathological_mask(i) = .true.
            ENDIF
         ENDIF
      ENDDO
      local_pathological = any(pathological_mask)

      IF (local_pathological) THEN
         DO i = 1, size(dt_res)
            IF (pathological_mask(i)) &
               dt_all(i) = min(ROUTING_PATHOLOGICAL_DT_FALLBACK, dt_res(i))
         ENDDO

         IF (p_is_worker .and. p_iam_worker == p_root .and. routing_zero_dt_warn_count < 5) THEN
            routing_zero_dt_warn_count = routing_zero_dt_warn_count + 1
            write(*,'(A)') 'WARNING grid_riverlake_flow: non-positive or non-finite adaptive dt; ' // &
               'using a bounded fallback to avoid a stalled routing loop.'
         ENDIF
      ENDIF

      ! One collective returns both the global adaptive step D and the global
      ! maximum residual R.  The second MPI_MIN lane carries -R.  After every
      ! active system subtracts the same D, another substep is needed exactly
      ! when R-D remains positive.
      dt_reduce = [huge(1._r8), 0._r8]
      IF (any(.not. ieee_is_finite(dt_res))) THEN
         ! A negative-HUGE sentinel makes every rank fail together after the
         ! collective instead of allowing a NaN residual to desynchronize the
         ! collective-bearing BIF loop.
         dt_reduce(1) = -huge(1._r8)
      ELSEIF (any(dt_res > 0._r8)) THEN
         dt_reduce(1) = minval(dt_all, mask = dt_res > 0._r8)
         dt_reduce(2) = -maxval(dt_res, mask = dt_res > 0._r8)
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, dt_reduce, 2, MPI_REAL8, MPI_MIN, &
         p_comm_worker, p_err)
#endif

      IF (.not. ieee_is_finite(dt_reduce(1))) THEN
         CALL CoLM_stop('grid_riverlake_flow: invalid synchronized routing dt')
      ELSEIF (dt_reduce(1) <= -0.5_r8 * huge(1._r8)) THEN
         CALL CoLM_stop('grid_riverlake_flow: non-finite routing residual')
      ENDIF

      dt_global = dt_reduce(1)
      max_residual = -dt_reduce(2)
      next_loop_active = .false.

      IF (dt_global < 0.5_r8 * huge(1._r8)) THEN
         IF (.not. ieee_is_finite(max_residual)) THEN
            CALL CoLM_stop('grid_riverlake_flow: invalid synchronized routing dt')
         ELSEIF (dt_global <= 0._r8) THEN
            CALL CoLM_stop('grid_riverlake_flow: invalid synchronized routing dt')
         ENDIF

         remaining_after = max_residual - dt_global
         IF (max_residual > 0._r8 .and. remaining_after >= max_residual) THEN
            CALL CoLM_stop('grid_riverlake_flow: routing dt makes no numerical progress')
         ENDIF
         next_loop_active = remaining_after > 0._r8

         WHERE (dt_res > 0._r8)
            dt_all = dt_global
         ELSEWHERE
            dt_all = 0._r8
         END WHERE
      ELSE
         dt_all = 0._r8
      ENDIF

   END SUBROUTINE sync_global_routing_dt

   SUBROUTINE publish_flood_feedback(year)
   ! Publish only exposed overbank/protected storage.  The same ucat/grid
   ! intersection areas are used in reverse by debit_flood_feedback.
   USE MOD_LandPatch, only: numpatch
#ifdef TRACER
   USE MOD_Tracer_Defs, only: ntracers, tracer_uses_land_water_transport
#endif
   IMPLICIT NONE
   integer, intent(in) :: year
   integer :: i
   real(r8) :: visible, protected, stage, protected_stage, protected_new, fraction
   real(r8), allocatable :: density(:), grid_visible(:), grid_protected(:)
   real(r8), allocatable :: frac_uc(:), frac_grid(:)
#ifdef TRACER
   integer :: itrc
   real(r8), allocatable :: grid_tracer(:), patch_tracer(:)
#endif

      IF (.not. p_is_worker) RETURN
      allocate(density(numucat), grid_visible(numinpm), grid_protected(numinpm), &
         frac_uc(numucat), frac_grid(numinpm))
      flood_visible_uc = 0._r8
      flood_protected_uc = 0._r8
      flood_reservoir_uc = .false.
#ifdef TRACER
      IF (allocated(flood_visible_tracer_uc)) THEN
         flood_visible_tracer_uc = 0._r8
         flood_protected_tracer_uc = 0._r8
         flood_tracer_credit_patch = 0._r8
         flood_tracer_evap_patch = 0._r8
         flood_tracer_land_patch = 0._r8
      ENDIF
#endif
      frac_uc = 0._r8
      DO i = 1, numucat
         IF (lake_type(i) == 2 .and. allocated(ucat2resv)) THEN
            IF (ucat2resv(i) > 0 .and. allocated(dam_build_year)) THEN
               flood_reservoir_uc(i) = year >= dam_build_year(ucat2resv(i))
            ENDIF
         ENDIF
         visible = max(0._r8, volwater_ucat(i))
         IF (flood_reservoir_uc(i)) visible = max(0._r8, volresv(ucat2resv(i)))
         IF (.not. volwater_ucat_valid) visible = max(0._r8, floodplain_curve(i)%volume(wdsrf_ucat(i)))
         IF (flood_reservoir_uc(i)) visible = max(0._r8, volresv(ucat2resv(i)))
         protected = 0._r8
         IF (DEF_USE_LEVEE .and. has_levee(i)) protected = max(0._r8, levsto(i))
#ifdef TRACER
         IF (allocated(trc_solid)) CALL equilibrate_river_tracer_cell(i, visible, protected)
#endif
         IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
            CALL levee_fldstg(i, visible+protected, stage, protected_new, protected_stage, fraction)
         ELSE
            stage = floodplain_curve(i)%depth(visible)
            fraction = floodplain_curve(i)%floodarea(stage) / max(topo_area(i), 1._r8)
         ENDIF
         IF (fraction <= 0._r8) CYCLE
         flood_visible_uc(i) = max(0._r8, visible - topo_rivstomax(i))
         flood_protected_uc(i) = protected
         frac_uc(i) = min(1._r8, fraction)
#ifdef TRACER
         IF (allocated(flood_visible_tracer_uc)) THEN
            IF (visible > 0._r8) flood_visible_tracer_uc(:,i) = &
               max(trc_mass(:,i), 0._r8) * (flood_visible_uc(i) / visible)
            IF (protected > 0._r8) flood_protected_tracer_uc(:,i) = max(trc_levsto(:,i), 0._r8)
         ENDIF
#endif
      ENDDO

      density = 0._r8
      WHERE (push_inpm2ucat%sum_area > 0._r8) &
         density = flood_visible_uc / max(push_inpm2ucat%sum_area, tiny(1._r8))
      CALL worker_push_data(push_ucat2inpm, density, grid_visible, fillvalue=0._r8, mode='sum')
      density = 0._r8
      WHERE (push_inpm2ucat%sum_area > 0._r8) &
         density = flood_protected_uc / max(push_inpm2ucat%sum_area, tiny(1._r8))
      CALL worker_push_data(push_ucat2inpm, density, grid_protected, fillvalue=0._r8, mode='sum')
      CALL worker_push_data(push_ucat2inpm, frac_uc, frac_grid, fillvalue=0._r8, mode='average')
      WHERE (flood_grid_area > 0._r8)
         frac_grid = min(1._r8, max(0._r8, frac_grid))
         grid_visible = grid_visible / max(flood_grid_area, tiny(1._r8))
         grid_protected = grid_protected / max(flood_grid_area, tiny(1._r8))
      ELSEWHERE
         frac_grid = 0._r8
         grid_visible = 0._r8
         grid_protected = 0._r8
      END WHERE
      IF (numpatch > 0) THEN
         CALL worker_remap_data_grid2pset(remap_patch2inpm, grid_visible, &
            flood_credit_patch, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
         WHERE (flood_credit_patch == RIVERLAKE_FLOOD_MISSING_VALUE) flood_credit_patch = 0._r8
         CALL worker_remap_data_grid2pset(remap_patch2inpm, grid_protected, &
            flood_depth_patch, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
         WHERE (flood_depth_patch == RIVERLAKE_FLOOD_MISSING_VALUE) flood_depth_patch = 0._r8
         flood_credit_patch = max(0._r8, flood_credit_patch) + max(0._r8, flood_depth_patch)
         CALL worker_remap_data_grid2pset(remap_patch2inpm, frac_grid, &
            flood_fraction_patch, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
         WHERE (flood_fraction_patch == RIVERLAKE_FLOOD_MISSING_VALUE) flood_fraction_patch = 0._r8
         flood_fraction_patch = min(1._r8, max(0._r8, flood_fraction_patch))
         flood_depth_patch = 0._r8
         WHERE (flood_fraction_patch > epsilon(1._r8)) &
            flood_depth_patch = 1000._r8 * flood_credit_patch / max(flood_fraction_patch, tiny(1._r8))
      ENDIF
#ifdef TRACER
      IF (allocated(flood_visible_tracer_uc)) THEN
         allocate(grid_tracer(numinpm), patch_tracer(numpatch))
         DO itrc = 1, ntracers
            IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
            density = 0._r8
            WHERE (push_inpm2ucat%sum_area > 0._r8) &
               density = flood_visible_tracer_uc(itrc,:) / &
                  max(push_inpm2ucat%sum_area, tiny(1._r8))
            CALL worker_push_data(push_ucat2inpm, density, grid_tracer, fillvalue=0._r8, mode='sum')
            WHERE (flood_grid_area > 0._r8)
               grid_tracer = grid_tracer / max(flood_grid_area, tiny(1._r8))
            ELSEWHERE
               grid_tracer = 0._r8
            END WHERE
            IF (numpatch > 0) THEN
               CALL worker_remap_data_grid2pset(remap_patch2inpm, grid_tracer, &
                  flood_tracer_credit_patch(itrc,:), fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
               WHERE (flood_tracer_credit_patch(itrc,:) == RIVERLAKE_FLOOD_MISSING_VALUE) &
                  flood_tracer_credit_patch(itrc,:) = 0._r8
            ENDIF
            density = 0._r8
            WHERE (push_inpm2ucat%sum_area > 0._r8) &
               density = flood_protected_tracer_uc(itrc,:) / &
                  max(push_inpm2ucat%sum_area, tiny(1._r8))
            CALL worker_push_data(push_ucat2inpm, density, grid_tracer, fillvalue=0._r8, mode='sum')
            WHERE (flood_grid_area > 0._r8)
               grid_tracer = grid_tracer / max(flood_grid_area, tiny(1._r8))
            ELSEWHERE
               grid_tracer = 0._r8
            END WHERE
            IF (numpatch > 0) THEN
               CALL worker_remap_data_grid2pset(remap_patch2inpm, grid_tracer, &
                  patch_tracer, fillvalue=RIVERLAKE_FLOOD_MISSING_VALUE, mode='average')
               WHERE (patch_tracer == RIVERLAKE_FLOOD_MISSING_VALUE) patch_tracer = 0._r8
               flood_tracer_credit_patch(itrc,:) = 1000._r8 * &
                  (flood_tracer_credit_patch(itrc,:) + patch_tracer)
            ENDIF
         ENDDO
      ENDIF
#endif
   END SUBROUTINE publish_flood_feedback

   SUBROUTINE debit_flood_feedback(evap_volume, infil_volume)
   USE MOD_LandPatch, only: numpatch
#ifdef TRACER
   USE MOD_Tracer_Defs, only: ntracers, tracer_uses_land_water_transport, &
      tracer_is_nonvolatile_solute, tracer_has_dissolved_limit
#endif
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE
   real(r8), intent(out) :: evap_volume, infil_volume
   real(r8), allocatable :: ratio_patch(:), ratio_grid(:), ratio_uc(:), fraction_sum(:)
#ifdef TRACER
   real(r8), allocatable :: gain_patch(:), coefficient_uc(:), gain_uc(:)
   real(r8), allocatable :: tracer_exchange_ledger(:,:)
   real(r8) :: patch_area, tracer_cell_debit, ledger_tol
   real(r8) :: water_credit, tracer_credit, vapor_loss, positive_loss, vapor_gain
   real(r8) :: infiltrated_fraction, evaporated_fraction, tracer_after
   integer :: itrc
#endif
   real(r8) :: debit_fraction, sink, credit, tol, fraction
   real(r8) :: visible_before, protected_before
   integer :: i, j

      evap_volume = 0._r8
      infil_volume = 0._r8
      IF (.not. p_is_worker) RETURN
      allocate(ratio_patch(numpatch), ratio_grid(numinpm), ratio_uc(numucat), fraction_sum(numucat))
      fraction_sum = 0._r8
#ifdef TRACER
      IF (allocated(flood_visible_tracer_uc)) THEN
         allocate(gain_patch(numpatch), coefficient_uc(numucat), gain_uc(numucat), &
            tracer_exchange_ledger(3,ntracers))
         tracer_exchange_ledger = 0._r8
      ENDIF
#endif
      DO i = 1, numpatch
         IF (.not. ieee_is_finite(flood_evap_acc(i)) .or. &
             .not. ieee_is_finite(flood_infil_acc(i))) &
            CALL CoLM_stop('grid flood feedback: nonfinite land sink')
         sink = (flood_evap_acc(i) + flood_infil_acc(i)) * 1.e-3_r8
         credit = flood_credit_patch(i)
         tol = 1.e-10_r8 * max(credit, 1.e-6_r8)
         IF (sink < -tol .or. sink > credit + tol) &
            CALL CoLM_stop('grid flood feedback: land exceeded published patch credit')
      ENDDO
      ! Land flux is patch-mean mm.  The ratio to published patch-mean
      ! metres reconstructs each patch/grid contribution without assuming
      ! that a patch belongs to only one grid.
      DO i = 1, 2
         ratio_patch = 0._r8
         IF (i == 1) THEN
            WHERE (flood_credit_patch > 0._r8) &
               ratio_patch = max(0._r8, flood_evap_acc*1.e-3_r8) / max(flood_credit_patch, tiny(1._r8))
         ELSE
            WHERE (flood_credit_patch > 0._r8) &
               ratio_patch = max(0._r8, flood_infil_acc*1.e-3_r8) / max(flood_credit_patch, tiny(1._r8))
         ENDIF
         CALL worker_remap_data_pset2grid(remap_patch2inpm, ratio_patch, ratio_grid, &
            fillvalue=0._r8, filter=flood_credit_patch>0._r8)
         WHERE (flood_grid_area > 0._r8)
            ratio_grid = ratio_grid / max(flood_grid_area, tiny(1._r8))
         ELSEWHERE
            ratio_grid = 0._r8
         END WHERE
         IF (any(ratio_grid < -1.e-10_r8) .or. any(ratio_grid > 1._r8+1.e-10_r8)) &
            CALL CoLM_stop('grid flood feedback: grid sink exceeded donor credit')
         CALL worker_push_data(push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
         DO j = 1, numucat
            IF (push_inpm2ucat%sum_area(j) > 0._r8) THEN
               debit_fraction = min(1._r8, max(0._r8, ratio_uc(j) / push_inpm2ucat%sum_area(j)))
            ELSE
               debit_fraction = 0._r8
            ENDIF
            IF (i == 1) THEN
               evap_volume = evap_volume + debit_fraction * &
                  (flood_visible_uc(j) + flood_protected_uc(j))
            ELSE
               infil_volume = infil_volume + debit_fraction * &
                  (flood_visible_uc(j) + flood_protected_uc(j))
            ENDIF
            fraction_sum(j) = fraction_sum(j) + debit_fraction
         ENDDO
      ENDDO
      IF (any(fraction_sum > 1._r8+1.e-10_r8)) &
         CALL CoLM_stop('grid flood feedback: donor overdraft')
#ifdef TRACER
      IF (allocated(flood_visible_tracer_uc)) THEN
         ! Patch-scale open-water evaporation (signed for isotopic vapour
         ! uptake) precedes infiltration. Adjoint coefficients preserve each
         ! source pool's concentration without homogenising ucatchments.
         DO itrc = 1, ntracers
            IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
            ratio_patch = 0._r8
            gain_patch = 0._r8
            DO i = 1, numpatch
               water_credit = flood_credit_patch(i)*1000._r8
               IF (water_credit <= 0._r8) CYCLE
               tracer_credit = flood_tracer_credit_patch(itrc,i)
               vapor_loss = flood_tracer_evap_patch(itrc,i)
               IF (.not. ieee_is_finite(vapor_loss)) &
                  CALL CoLM_stop('grid flood feedback: nonfinite isotope vapor exchange')
               positive_loss = max(vapor_loss,0._r8)
               vapor_gain = max(-vapor_loss,0._r8)
               ledger_tol = max(1.e-12_r8, 1.e-10_r8*abs(tracer_credit))
               IF (positive_loss > tracer_credit+ledger_tol) &
                  CALL CoLM_stop('grid flood feedback: isotope evaporation overdrew tracer')
               evaporated_fraction = 0._r8
               IF (tracer_credit > 0._r8) evaporated_fraction = &
                  min(1._r8,positive_loss/tracer_credit)
               infiltrated_fraction = 0._r8
               IF (flood_infil_acc(i) > 0._r8) THEN
                  IF (water_credit <= flood_evap_acc(i)) &
                     CALL CoLM_stop('grid flood feedback: infiltration after full evaporation')
                  infiltrated_fraction = flood_infil_acc(i)/(water_credit-flood_evap_acc(i))
               ENDIF
               IF (infiltrated_fraction > 1._r8+1.e-10_r8) &
                  CALL CoLM_stop('grid flood feedback: tracer infiltration overdraft')
               infiltrated_fraction = min(1._r8,max(0._r8,infiltrated_fraction))
               ratio_patch(i) = evaporated_fraction + &
                  (1._r8-evaporated_fraction)*infiltrated_fraction
               IF (tracer_has_dissolved_limit(itrc) .and. tracer_credit > 0._r8) &
                  ratio_patch(i) = min(1._r8,max(0._r8,flood_tracer_land_patch(itrc,i)/tracer_credit))
               gain_patch(i) = vapor_gain*(1._r8-infiltrated_fraction)/water_credit
               patch_area = sum(remap_patch2inpm%areapart(i)%val)*1.e-3_r8
               IF (tracer_has_dissolved_limit(itrc)) THEN
                  IF (flood_tracer_land_patch(itrc,i) < -ledger_tol .or. &
                      flood_tracer_land_patch(itrc,i) > &
                      (tracer_credit-vapor_loss)*infiltrated_fraction + ledger_tol) &
                     CALL CoLM_stop('grid flood feedback: finite solute input exceeds donor exchange')
               ELSE
                  IF (abs(flood_tracer_land_patch(itrc,i) - &
                      (tracer_credit-vapor_loss)*infiltrated_fraction) > &
                      max(1.e-12_r8,1.e-10_r8*max(abs(tracer_credit),abs(flood_tracer_land_patch(itrc,i))))) &
                     CALL CoLM_stop('grid flood feedback: land tracer input differs from donor exchange')
               ENDIF
               tracer_exchange_ledger(1,itrc) = tracer_exchange_ledger(1,itrc) + &
                  flood_tracer_land_patch(itrc,i)*patch_area
               tracer_exchange_ledger(2,itrc) = tracer_exchange_ledger(2,itrc) + &
                  vapor_loss*patch_area
            ENDDO
            CALL worker_remap_data_pset2grid(remap_patch2inpm, ratio_patch, ratio_grid, &
               fillvalue=0._r8, filter=flood_credit_patch>0._r8)
            WHERE (flood_grid_area > 0._r8)
               ratio_grid = ratio_grid/max(flood_grid_area,tiny(1._r8))
            ELSEWHERE
               ratio_grid = 0._r8
            END WHERE
            CALL worker_push_data(push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
            coefficient_uc = 0._r8
            WHERE (push_inpm2ucat%sum_area > 0._r8) &
               coefficient_uc = min(1._r8,max(0._r8,ratio_uc/push_inpm2ucat%sum_area))
            gain_uc = 0._r8
            IF (.not. tracer_is_nonvolatile_solute(itrc)) THEN
               CALL worker_remap_data_pset2grid(remap_patch2inpm, gain_patch, ratio_grid, &
                  fillvalue=0._r8, filter=flood_credit_patch>0._r8)
               WHERE (flood_grid_area > 0._r8)
                  ratio_grid = ratio_grid/max(flood_grid_area,tiny(1._r8))
               ELSEWHERE
                  ratio_grid = 0._r8
               END WHERE
               CALL worker_push_data(push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
               WHERE (push_inpm2ucat%sum_area > 0._r8) &
                  gain_uc = ratio_uc/push_inpm2ucat%sum_area
            ENDIF
            DO j = 1, numucat
               tracer_after = trc_mass(itrc,j) - coefficient_uc(j)*flood_visible_tracer_uc(itrc,j) &
                  + gain_uc(j)*flood_visible_uc(j)
               IF (tracer_after < -1.e-10_r8*max(1._r8,trc_mass(itrc,j))) &
                  CALL CoLM_stop('grid flood feedback: negative visible tracer')
               ! Accumulate explicit transfer instead of subtracting large
               ! stored masses; include the effect of the existing zero clamp.
               tracer_cell_debit = coefficient_uc(j)*flood_visible_tracer_uc(itrc,j) &
                  - gain_uc(j)*flood_visible_uc(j) + min(tracer_after,0._r8)
               trc_mass(itrc,j) = max(0._r8,tracer_after)
               tracer_after = trc_levsto(itrc,j) - coefficient_uc(j)*flood_protected_tracer_uc(itrc,j) &
                  + gain_uc(j)*flood_protected_uc(j)
               IF (tracer_after < -1.e-10_r8*max(1._r8,trc_levsto(itrc,j))) &
                  CALL CoLM_stop('grid flood feedback: negative protected tracer')
               tracer_cell_debit = tracer_cell_debit + coefficient_uc(j)*flood_protected_tracer_uc(itrc,j) &
                  - gain_uc(j)*flood_protected_uc(j) + min(tracer_after,0._r8)
               ! Unlimited nonvolatile solute retains the legacy protected
               ! pool behavior; finite solute is partitioned into solid below.
               trc_levsto(itrc,j) = max(0._r8,tracer_after)
               tracer_exchange_ledger(3,itrc) = tracer_exchange_ledger(3,itrc) + &
                  tracer_cell_debit
            ENDDO
         ENDDO
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, tracer_exchange_ledger, 3*ntracers, MPI_REAL8, &
            MPI_SUM, p_comm_worker, p_err)
#endif
         DO itrc = 1, ntracers
            IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
            ledger_tol = max(1.e-8_r8, 1.e-10_r8*maxval(abs(tracer_exchange_ledger(:,itrc))))
            IF (abs(tracer_exchange_ledger(1,itrc)+tracer_exchange_ledger(2,itrc) - &
                   tracer_exchange_ledger(3,itrc)) > ledger_tol) &
               CALL CoLM_stop('grid flood feedback: land/river tracer exchange ledger mismatch')
            IF (p_iam_worker == p_root .and. flood_tracer_ledger_reports < 5 .and. &
                tracer_exchange_ledger(1,itrc) > 0._r8) THEN
               write(*,'(A,I0,4(1X,ES16.8))') 'Grid flood tracer land/vapor/donor/residual:', itrc, &
                  tracer_exchange_ledger(1,itrc), tracer_exchange_ledger(2,itrc), &
                  tracer_exchange_ledger(3,itrc), &
                  tracer_exchange_ledger(1,itrc)+tracer_exchange_ledger(2,itrc)-tracer_exchange_ledger(3,itrc)
               flood_tracer_ledger_reports = flood_tracer_ledger_reports + 1
            ENDIF
         ENDDO
      ENDIF
#endif
      DO j = 1, numucat
         IF (flood_reservoir_uc(j)) THEN
            volresv(ucat2resv(j)) = max(0._r8, volresv(ucat2resv(j)) - fraction_sum(j)*flood_visible_uc(j))
            volwater_ucat(j) = volresv(ucat2resv(j))
         ELSE
            volwater_ucat(j) = max(0._r8, volwater_ucat(j) - fraction_sum(j)*flood_visible_uc(j))
         ENDIF
         IF (DEF_USE_LEVEE .and. has_levee(j)) &
            levsto(j) = max(0._r8, levsto(j) - fraction_sum(j)*flood_protected_uc(j))
#ifdef TRACER
         IF (allocated(trc_solid)) THEN
            protected_before = 0._r8
            IF (DEF_USE_LEVEE .and. has_levee(j)) protected_before = levsto(j)
            CALL equilibrate_river_tracer_cell(j, volwater_ucat(j), protected_before)
         ENDIF
#endif
         IF (fraction_sum(j) > 0._r8) THEN
            IF (DEF_USE_LEVEE .and. has_levee(j) .and. .not. flood_reservoir_uc(j)) THEN
               CALL levee_repartition_storage(j, volwater_ucat(j), wdsrf_ucat(j), &
                  fraction, visible_before, protected_before)
#ifdef TRACER
               IF (ntracers > 0) CALL tracer_lifecycle_route_sediment_levee_repartition(j, &
                  visible_before, protected_before, volwater_ucat(j), levsto(j))
               IF (ntracers > 0) CALL levee_tracer_repartition(j, &
                  visible_before, protected_before, volwater_ucat(j), levsto(j))
#endif
            ELSEIF (flood_reservoir_uc(j)) THEN
               wdsrf_ucat(j) = floodplain_curve(j)%depth(volresv(ucat2resv(j)))
            ELSE
               wdsrf_ucat(j) = floodplain_curve(j)%depth(volwater_ucat(j))
            ENDIF
         ENDIF
      ENDDO
      flood_evap_acc = 0._r8
      flood_infil_acc = 0._r8
   END SUBROUTINE debit_flood_feedback

   ! ---------
   SUBROUTINE grid_riverlake_flow_final ()

      CALL riverlake_network_final ()

      IF (DEF_Reservoir_Method > 0) THEN
         CALL reservoir_final ()
      ENDIF

#ifdef TRACER
      CALL tracer_lifecycle_route_final()
#endif

         CALL levee_final()
      CALL bifurcation_final()
#ifdef TRACER
      CALL river_lake_tracer_final()
#endif

      ! acc_rnof_uc is owned by MOD_Grid_RiverLakeTimeVars and freed by
      ! deallocate_GridRiverLakeTimeVars; don't deallocate it here.
      IF (allocated(filter_rnof)) deallocate(filter_rnof)
      IF (allocated(flood_depth_patch)) deallocate(flood_depth_patch, flood_fraction_patch, &
         flood_evap_patch, flood_infil_patch)
      IF (allocated(flood_visible_uc)) deallocate(flood_visible_uc, flood_protected_uc, &
         flood_credit_patch, flood_evap_acc, flood_infil_acc, flood_reservoir_uc, flood_grid_area)
#ifdef TRACER
      IF (allocated(flood_tracer_credit_patch)) deallocate(flood_tracer_credit_patch, &
         flood_tracer_evap_patch, flood_tracer_land_patch, &
         flood_visible_tracer_uc, flood_protected_tracer_uc)
#endif

   END SUBROUTINE grid_riverlake_flow_final

END MODULE MOD_Grid_RiverLakeFlow
#endif
