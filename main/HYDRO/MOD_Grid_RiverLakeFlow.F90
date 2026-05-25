#include <define.h>

#ifdef GridRiverLakeFlow
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
   USE MOD_Grid_RiverLakeHist
   USE MOD_Grid_RiverLakeLevee
   USE MOD_Grid_RiverLakeBifurcation
#ifdef GridRiverLakeSediment
   USE MOD_Grid_RiverLakeSediment, only: grid_sediment_init, grid_sediment_calc, &
      grid_sediment_final, sediment_diag_accumulate, sediment_forcing_put, &
      read_sediment_restart
#endif
#ifdef TRACER
   USE MOD_Grid_RiverLakeTracer, only: tracer_init, tracer_init_from_water, &
      tracer_input_from_runoff, &
      tracer_substep, tracer_flush_acc, &
      read_tracer_restart, tracer_final, acc_trc_inp, acc_rnof_ref, trc_mass, trc_inp_buf, trc_flux_out, &
      tracer_refresh_state, tracer_diag_accumulate_substep, &
      trc_levsto, trc_dry_drain, trc_reactive_source, levee_tracer_repartition, &
      get_cell_volume_dep => get_cell_volume, trc_conc_dep => trc_conc
#endif
#if (defined TRACER) && (defined BGC)
   USE MOD_Tracer_Methane_Registry, only: igas_ch4
#endif
   IMPLICIT NONE

   real(r8), parameter :: RIVERMIN  = 1.e-5_r8

   real(r8), save :: acctime_rnof_max

   ! acctime_rnof (scalar) and acc_rnof_uc (:) are owned by
   ! MOD_Grid_RiverLakeTimeVars (imported via the module-wide USE above)
   ! so their mid-period state is serialised by WRITE/READ_GridRiverLakeTimeVars.
   logical,  allocatable :: filter_rnof (:)
CONTAINS

   ! ---------
   SUBROUTINE grid_riverlake_flow_init ()

   USE MOD_LandPatch,           only: numpatch
   USE MOD_Forcing,             only: forcmask_pch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
#ifdef TRACER
   USE MOD_Tracer_Defs,          only: ntracers
#endif
   IMPLICIT NONE

#ifdef TRACER
      logical :: trc_restart_found
      logical, allocatable :: trc_missing(:)
      real(r8), allocatable :: wdsrf_safe(:), volresv_safe(:)
      integer,  allocatable :: ucat2resv_safe(:)
#endif

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

#ifdef GridRiverLakeSediment
      CALL grid_sediment_init()
      IF (len_trim(gridriver_restart_file) > 0) THEN
         CALL read_sediment_restart(gridriver_restart_file)
      ENDIF
#endif

      ! Always call levee_init: when DEF_USE_LEVEE=.false. it allocates the
      ! levee arrays in inert state (has_levee=.false. everywhere), so guards
      ! like `IF (DEF_USE_LEVEE .and. has_levee(i) ...)` can safely evaluate
      ! both operands under ifx -check bounds.
      CALL levee_init()
      IF (DEF_USE_LEVEE .and. len_trim(gridriver_restart_file) > 0) THEN
         CALL read_levee_restart(gridriver_restart_file)
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         CALL bifurcation_init()
         IF (len_trim(gridriver_restart_file) > 0) THEN
            CALL read_bifurcation_restart(gridriver_restart_file)
         ENDIF
      ENDIF

#ifdef TRACER
         trc_restart_found = .false.
         CALL tracer_init()
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
               CALL tracer_init_from_water(wdsrf_safe, volresv_safe, ucat2resv_safe, trc_missing)
               deallocate(wdsrf_safe, volresv_safe, ucat2resv_safe)
            ENDIF
            deallocate(trc_missing)
         ENDIF
#endif

   END SUBROUTINE grid_riverlake_flow_init

   ! ---------
   SUBROUTINE grid_riverlake_flow (year, deltime)

   USE MOD_Utils
   USE MOD_Namelist,       only: DEF_Reservoir_Method, DEF_USE_SEDIMENT, DEF_USE_LEVEE, DEF_USE_BIFURCATION
   USE MOD_Vars_1DFluxes,  only: rnof
   USE MOD_Mesh,           only: numelm
   USE MOD_LandPatch,      only: elm_patch, numpatch
   USE MOD_Const_Physical, only: grav
   USE MOD_Vars_Global,    only: spval
#ifdef TRACER
   USE MOD_Tracer_Defs,    only: ntracers
#endif
#ifdef TRACER
   USE MOD_Tracer_Vars,    only: trc_rnof_step
#endif
#ifdef GridRiverLakeSediment
   USE MOD_Vars_1DForcing, only: forc_prc, forc_prl
#endif
   IMPLICIT NONE

   integer,  intent(in) :: year
   real(r8), intent(in) :: deltime

   ! Local Variables
   integer  :: i, j, irsv, ntimestep, ipth, i_up, itrc
   real(r8) :: dt_this
   integer  :: sed_clock_start, sed_clock_end, sed_clock_rate
   real(r8) :: sed_elapsed

   real(r8), allocatable :: rnof_gd(:)
   real(r8), allocatable :: rnof_uc(:)
   real(r8), allocatable :: trc_rnof_gd(:,:)
   real(r8), allocatable :: trc_rnof_uc(:,:)

#ifdef GridRiverLakeSediment
   real(r8), allocatable :: prcp_gd(:)
   real(r8), allocatable :: prcp_uc(:)
   real(r8), allocatable :: prcp_pch(:)
   real(r8), allocatable :: floodarea_sed(:)
#endif

   logical,  allocatable :: is_built_resv(:)

   real(r8), allocatable :: wdsrf_next(:)
   real(r8), allocatable :: veloc_next(:)

   real(r8), allocatable :: hflux_fc(:)
   real(r8), allocatable :: mflux_fc(:)
   real(r8), allocatable :: zgrad_dn(:)

   real(r8), allocatable :: hflux_resv(:)
   real(r8), allocatable :: mflux_resv(:)

   real(r8), allocatable :: hflux_sumups(:)
   real(r8), allocatable :: mflux_sumups(:)
   real(r8), allocatable :: zgrad_sumups(:)

   real(r8), allocatable :: sum_hflux_riv(:)
   real(r8), allocatable :: sum_hflux_base(:)
   real(r8), allocatable :: sum_mflux_riv(:)
   real(r8), allocatable :: sum_zgrad_riv(:)

   real(r8) :: veloct_fc, height_fc, momen_fc, zsurf_fc
   real(r8) :: bedelv_fc, height_up, height_dn
   real(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: volwater, friction, floodarea
   real(r8) :: visible_hflux, protected_hflux, protected_avail
   real(r8) :: vol_total_levee, fldfrc_levee
   real(r8) :: vis_vol_bef_lv, levsto_bef_lv
   real(r8) :: vis_vol_bef_lv2, levsto_bef_lv2
   real(r8), allocatable :: volresv_safe(:)
   real(r8), allocatable :: volresv_safe_trc(:)
   real(r8), allocatable :: volresv_safe_dep(:)
   real(r8), allocatable :: volresv_safe_diag(:)
   integer,  allocatable :: ucat2resv_safe_trc(:)
   integer,  allocatable :: ucat2resv_safe_dep(:)
   integer,  allocatable :: ucat2resv_safe_diag(:)
   integer :: itrc_dep
   real(r8) :: frac_remove, trc_removed, vol_post
   real(r8), allocatable :: levee_floodarea(:)
   real(r8), allocatable :: total_floodarea(:)   ! general floodarea (levee+floodplain), per ucat
   real(r8), allocatable :: total_flooddepth(:)  ! floodplain water depth flddph [m], per ucat
   real(r8),  allocatable :: dt_res(:), dt_all(:)
   logical,   allocatable :: ucatfilter(:)
   logical :: loop_active, dt_changed
   integer :: ibif_iter
   real(r8) :: dt_before_sync, dt_after_sync
#ifdef CoLMDEBUG
   real(r8) :: totalvol_bef, totalvol_aft, totalrnof, totaldis, totalclip
   real(r8), allocatable :: trc_mass_bef(:), trc_mass_aft(:)
   real(r8), allocatable :: trc_mass_inp(:), trc_mass_dis(:), trc_mass_reactive(:)
   real(r8) :: bif_flux_sum_total, bif_flux_sum_max
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


      IF (p_is_worker) THEN
         allocate (rnof_gd (numinpm))
         allocate (rnof_uc (numucat))
#ifdef TRACER
         IF (ntracers > 0) THEN
            allocate (trc_rnof_gd (ntracers, numinpm))
            allocate (trc_rnof_uc (ntracers, numucat))
            trc_rnof_gd = 0._r8
            trc_rnof_uc = 0._r8
         END IF
#endif

         CALL worker_remap_data_pset2grid (remap_patch2inpm, rnof, rnof_gd, &
            fillvalue = 0., filter = filter_rnof)

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
	                  ! rnof_uc * 1.e-3 * deltime is in m (depth, not volume)
	                  ! Need to convert trc_rnof_uc from mm to m to match water units
	                  CALL tracer_input_from_runoff(rnof_uc*1.e-3*deltime, numucat, trc_rnof_uc*1.e-3)
	               ENDIF
#endif
         ENDIF

         deallocate(rnof_gd)
         deallocate(rnof_uc)
         IF (allocated(trc_rnof_gd)) deallocate(trc_rnof_gd)
         IF (allocated(trc_rnof_uc)) deallocate(trc_rnof_uc)

#ifdef GridRiverLakeSediment
         IF (DEF_USE_SEDIMENT) THEN
            ! Allocate zero-length arrays on empty workers to avoid passing unallocated
            ! arrays to assumed-shape dummy arguments in MPI communication routines.
            IF (numpatch > 0) THEN
               allocate (prcp_pch (numpatch))
               prcp_pch = forc_prc + forc_prl
            ELSE
               allocate (prcp_pch (0))
            ENDIF
            IF (numinpm > 0) THEN
               allocate (prcp_gd (numinpm))
            ELSE
               allocate (prcp_gd (0))
            ENDIF
            IF (numucat > 0) THEN
               allocate (prcp_uc (numucat))
            ELSE
               allocate (prcp_uc (0))
            ENDIF

            CALL worker_remap_data_pset2grid (remap_patch2inpm, prcp_pch, prcp_gd, &
               fillvalue = 0., filter = filter_rnof)

            IF (numinpm > 0) THEN
               WHERE (push_ucat2inpm%sum_area > 0)
                  prcp_gd = prcp_gd / push_ucat2inpm%sum_area
               END WHERE
            ENDIF

            CALL worker_push_data (push_inpm2ucat, prcp_gd, prcp_uc, &
               fillvalue = 0., mode = 'sum')

            ! Convert from area-integrated [mm/s * m²] back to flux density [mm/s].
            ! push_data(mode='sum') produces area-integrated values (like rnof_uc),
            ! but the sediment yield formula expects a rate and multiplies by area internally.
            IF (numucat > 0) THEN
               WHERE (topo_area > 0._r8)
                  prcp_uc = prcp_uc / topo_area
               END WHERE
            ENDIF

            CALL sediment_forcing_put(prcp_uc, deltime)

            deallocate(prcp_pch)
            deallocate(prcp_gd)
            deallocate(prcp_uc)
         ENDIF
#endif

      ENDIF


      acctime_rnof = acctime_rnof + deltime

      IF (acctime_rnof+0.01 < acctime_rnof_max) THEN
         RETURN
      ENDIF


      IF (p_is_worker) THEN

         IF (numucat > 0) THEN

            allocate (is_built_resv (numucat))
            allocate (wdsrf_next    (numucat))
            allocate (veloc_next    (numucat))
            allocate (hflux_fc      (numucat))
            allocate (mflux_fc      (numucat))
            allocate (zgrad_dn      (numucat))
            allocate (sum_hflux_riv (numucat))
            allocate (sum_hflux_base(numucat))
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

            allocate (hflux_sumups  (numucat))
            allocate (mflux_sumups  (numucat))
            allocate (zgrad_sumups  (numucat))

            IF (DEF_Reservoir_Method > 0) THEN
               allocate (hflux_resv (numucat))
               allocate (mflux_resv (numucat))
            ENDIF

            allocate (dt_res (numrivsys))
            allocate (dt_all (numrivsys))

         ENDIF

#ifdef CoLMDEBUG
         totalrnof = sum(acc_rnof_uc)
         totalvol_bef = 0.
#endif

         DO i = 1, numucat

            is_built_resv(i) = .false.
            IF (DEF_USE_LEVEE .and. has_levee(i) .and. lake_type(i) == 2) THEN
               write(*,'(A,I0,A,I0)') 'ERROR: grid river-lake cell cannot be both reservoir and levee; ucat=', &
                  ucat_ucid(i), ' local=', i
               CALL CoLM_stop ('grid river-lake cell cannot be both reservoir and levee')
            ENDIF
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

#ifdef CoLMDEBUG
            IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
               totalvol_bef = totalvol_bef + volwater + levsto(i)
            ELSE
               totalvol_bef = totalvol_bef + volwater
            ENDIF
#endif

            volwater = volwater + acc_rnof_uc(i)

            IF (.not. is_built_resv(i)) THEN
               IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  ! Pass acc_trc_inp(:, i) as pending pool. At
                  ! this point volwater already includes acc_rnof_uc(i)
                  ! (line above) but trc_mass has not yet absorbed
                  ! acc_trc_inp via the trc_inp_buf release; without the
                  ! pending arg, the ratio (trc_mass / vis_vol_bef) is
                  ! diluted and the transfer leaks Phase-1 R_init.
                  vis_vol_bef_lv = volwater
                  levsto_bef_lv  = levsto(i)
                  vol_total_levee = volwater + levsto(i)
                  CALL levee_fldstg(i, vol_total_levee, wdsrf_ucat(i), &
                     levsto(i), levdph(i), fldfrc_levee)
                  volwater_ucat(i) = vol_total_levee - levsto(i)
#ifdef TRACER
                     ! Also pass acc_rnof_ref(i) so the routine can
                     ! debit the runoff WATER reference proportionally
                     ! when the pending tracer pool spills (overflow
                     ! path). Keeps acc_trc_inp/acc_rnof_ref
                     ! in lockstep instead of drifting apart.
                     CALL levee_tracer_repartition(i, &
                        vis_vol_bef_lv, levsto_bef_lv, &
                        volwater_ucat(i), levsto(i), &
                        pending_trc_pool = acc_trc_inp(:, i), &
                        pending_water_ref = acc_rnof_ref(i))
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
#ifdef CoLMDEBUG
         totaldis  = 0.
         totalclip = 0._r8
         bif_flux_sum_total = 0._r8
         bif_flux_sum_max   = 0._r8
         ! Tracer conservation: save total mass BEFORE input addition.
         ! Include protected-side pool so mass trapped behind a
         ! levee is not treated as "missing" by the before/after check.
#ifdef TRACER
         IF (numucat > 0) THEN
            DO itrc = 1, ntracers
               trc_mass_bef(itrc) = sum(trc_mass(itrc,:)) + sum(trc_inp_buf(itrc,:))
               IF (allocated(trc_levsto)) &
                  trc_mass_bef(itrc) = trc_mass_bef(itrc) + sum(trc_levsto(itrc,:))
               trc_mass_inp(itrc) = sum(acc_trc_inp(itrc,:))
               trc_mass_dis(itrc) = 0._r8
               trc_mass_reactive(itrc) = 0._r8
            ENDDO
            IF (allocated(trc_reactive_source)) trc_reactive_source = 0._r8
         ELSE
            trc_mass_bef = 0._r8
            trc_mass_inp = 0._r8
            trc_mass_dis = 0._r8
            trc_mass_reactive = 0._r8
         END IF
#else
         trc_mass_bef = 0._r8
         trc_mass_inp = 0._r8
         trc_mass_dis = 0._r8
         trc_mass_reactive = 0._r8
#endif
#endif

         dt_res(:) = acctime_rnof

         ! CaMa-style routing uses one adaptive substep DT for the whole
         ! worker domain. Keep every worker in this loop until all river
         ! systems have exhausted dt_res, otherwise a worker that finishes
         ! early may skip collectives while others still need push_data or
         ! global-DT synchronization.
         loop_active = any(dt_res > 0)
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, loop_active, 1, MPI_LOGICAL, &
            MPI_LOR, p_comm_worker, p_err)
#endif

         DO WHILE (loop_active)

            ntimestep = ntimestep + 1

            CALL worker_push_data (push_next2ucat, wdsrf_ucat, wdsrf_next, fillvalue = spval)
            ! velocity in ocean or inland depression is assumed to be 0.
            CALL worker_push_data (push_next2ucat, veloc_riv,  veloc_next, fillvalue = 0.)

            dt_all(:) = min(dt_res(:), 60.)
            CALL sync_global_routing_dt(dt_res, dt_all)

            DO i = 1, numucat

               ! Bounds-guard irivsys: a corrupt/partial restart or
               ! malformed network metadata can leave irivsys(i)==0,
               ! negative, or > size(dt_all). Unlike the tracer_substep
               ! guards at MOD_Grid_RiverLakeTracer.F90:482,569,786, this
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

            CALL worker_push_data (push_ups2ucat, hflux_fc, hflux_sumups, fillvalue = 0., mode = 'sum')
            CALL worker_push_data (push_ups2ucat, mflux_fc, mflux_sumups, fillvalue = 0., mode = 'sum')
            CALL worker_push_data (push_ups2ucat, zgrad_dn, zgrad_sumups, fillvalue = 0., mode = 'sum')

            IF (numucat > 0) THEN
               WHERE (ucatfilter)
                  sum_hflux_riv = sum_hflux_riv - hflux_sumups
                  sum_mflux_riv = sum_mflux_riv - mflux_sumups
                  sum_zgrad_riv = sum_zgrad_riv - zgrad_sumups
               END WHERE
            ENDIF

            ! reservoir operation.
            IF (DEF_Reservoir_Method > 0) THEN

               DO i = 1, numucat

                  IF ((.not. ucatfilter(i)) .or. (ucat_next(i) == -10)) CYCLE

                  hflux_resv(i) = 0.
                  mflux_resv(i) = 0.

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

               CALL worker_push_data (push_ups2ucat, hflux_resv, hflux_sumups, fillvalue = 0., mode = 'sum')
               CALL worker_push_data (push_ups2ucat, mflux_resv, mflux_sumups, fillvalue = 0., mode = 'sum')

               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_riv - hflux_sumups
                     sum_mflux_riv = sum_mflux_riv - mflux_sumups
                  END WHERE
               ENDIF

            ENDIF

            IF (numucat > 0) sum_hflux_base = sum_hflux_riv

	            ! ----- Bifurcation pathways -----
	            IF (DEF_USE_BIFURCATION) THEN
	               IF (allocated(volresv)) THEN
	                  allocate(volresv_safe(size(volresv))); volresv_safe = volresv
	               ELSE
	                  allocate(volresv_safe(0))
	               ENDIF
		               CALL bifurcation_calc(wdsrf_ucat, volresv_safe, is_built_resv, dt_all, &
		                  irivsys, ucatfilter, update_state = .false.)
	               deallocate(volresv_safe)
               ! Add predicted bifurcation volume flux to the timestep
               ! constraints (volume only, not momentum). The prediction is
               ! removed and recomputed after the final synchronized dt_all
               ! is known.
               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_riv + bif_hflux_sum
                  END WHERE
               ENDIF
            ENDIF

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

               dt_this = dt_all(irivsys(i))

               ! constraint 1: CFL condition (only for rivers)
               IF (.not. is_built_resv(i)) THEN
                  IF ((veloc_riv(i) /= 0.) .or. (wdsrf_ucat(i) > 0.)) THEN
                     dt_this = min(dt_this, topo_rivlen(i)/(abs(veloc_riv(i))+sqrt(grav*wdsrf_ucat(i)))*0.8)
                  ENDIF
               ENDIF

               ! constraint 2: Avoid negative values of water.  On leveed
               ! cells, bifurcation layer-2+ flux belongs to the protected
               ! pool (`levsto`), not to the visible routing volume.  Limit
               ! each donor compartment against its own available storage so
               ! protected-side outflow is not artificially throttled by a
               ! dry visible side (and vice versa).
               visible_hflux = sum_hflux_riv(i)
               protected_hflux = 0._r8
               IF (DEF_USE_LEVEE .and. DEF_USE_BIFURCATION .and. allocated(bif_lev_hflux_sum)) THEN
                  IF ((.not. is_built_resv(i)) .and. has_levee(i) .and. i <= size(bif_lev_hflux_sum)) THEN
                     protected_hflux = bif_lev_hflux_sum(i)
                     visible_hflux = sum_hflux_riv(i) - protected_hflux
                  ENDIF
               ENDIF
               IF (.not. is_built_resv(i)) THEN
                  volwater = volwater_ucat(i)
                  IF (visible_hflux > 0._r8) dt_this = min(dt_this, volwater / visible_hflux)
                  IF (protected_hflux > 0._r8) THEN
                     protected_avail = 0._r8
                     IF (DEF_USE_LEVEE .and. has_levee(i)) protected_avail = max(levsto(i), 0._r8)
                     dt_this = min(dt_this, protected_avail / protected_hflux)
                  ENDIF
               ELSE
                  IF (sum_hflux_riv(i) > 0._r8) THEN
                     volwater = volresv(ucat2resv(i))
                     dt_this = min(dt_this, volwater / sum_hflux_riv(i))
                  ENDIF
               ENDIF

               ! constraint 3: Avoid change of flow direction (only for rivers)
               ! IF (.not. is_built_resv(i)) THEN
               !    IF ((abs(veloc_riv(i)) > 0.1) &
               !       .and. (veloc_riv(i) * (sum_mflux_riv(i)-sum_zgrad_riv(i)) > 0)) THEN
               !       dt_this = min(dt_this, &
               !          abs(momen_riv(i) * topo_rivare(i) / (sum_mflux_riv(i)-sum_zgrad_riv(i))))
               !    ENDIF
               ! ENDIF

               dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))

            ENDDO

            ! No dedicated CFL constraint for bifurcation pathways.
            ! This is an engineering trade-off (not a strict stability proof):
            ! the momentum driving term (mflux_lev) is explicit, so the scheme
            ! is not unconditionally stable in theory. In practice, stability
            ! is maintained by: (1) 5% storage limiter per pathway,
            ! (2) volume constraint (constraint 2) via sum_hflux_riv,
            ! (3) semi-implicit friction damping in momentum update,
            ! (4) velocity clamp at +/-20 m/s.
            ! This follows CaMa-Flood's approach and avoids short pathways
            ! (small pth_dst) forcing excessively small timesteps.
            ! If bifurcation instability is observed, re-adding a CFL
            ! constraint based on pth_dst is the first thing to try.

            CALL sync_global_routing_dt(dt_res, dt_all)

            IF (DEF_USE_BIFURCATION .and. numucat > 0) THEN
               WHERE (ucatfilter)
                  sum_hflux_riv = sum_hflux_riv - bif_hflux_sum
               END WHERE
            ENDIF

            IF (DEF_USE_BIFURCATION) THEN
               ! Bifurcation's storage limiter depends on dt. The first
               ! prediction above may therefore become inconsistent after the
               ! synchronized global dt is reduced by CFL/storage constraints.
               ! Iterate with update_state=.false. until the final bif fluxes
               ! also satisfy the no-negative-storage constraint, then advance
               ! bifurcation momentum exactly once with that final dt.
               DO ibif_iter = 1, 8
                  IF (numucat > 0) THEN
                     WHERE (ucatfilter)
                        sum_hflux_riv = sum_hflux_base
                     END WHERE
                  ENDIF

	                  IF (allocated(volresv)) THEN
	                     allocate(volresv_safe(size(volresv))); volresv_safe = volresv
	                  ELSE
	                     allocate(volresv_safe(0))
	                  ENDIF
		                  CALL bifurcation_calc(wdsrf_ucat, volresv_safe, is_built_resv, dt_all, &
		                     irivsys, ucatfilter, update_state = .false.)
	                  deallocate(volresv_safe)

                  IF (numucat > 0) THEN
                     WHERE (ucatfilter)
                        sum_hflux_riv = sum_hflux_base + bif_hflux_sum
                     END WHERE
                  ENDIF

                  dt_changed = .false.
                  IF (any(dt_res > 0._r8)) THEN
                     dt_before_sync = maxval(dt_all, mask = dt_res > 0._r8)
                  ELSE
                     dt_before_sync = 0._r8
                  ENDIF

                  DO i = 1, numucat
                     IF (.not. ucatfilter(i)) CYCLE
	                     visible_hflux = sum_hflux_riv(i)
                     protected_hflux = 0._r8
                     IF (DEF_USE_LEVEE .and. DEF_USE_BIFURCATION .and. allocated(bif_lev_hflux_sum)) THEN
                        IF ((.not. is_built_resv(i)) .and. has_levee(i) .and. i <= size(bif_lev_hflux_sum)) THEN
                           protected_hflux = bif_lev_hflux_sum(i)
                           visible_hflux = sum_hflux_riv(i) - protected_hflux
                        ENDIF
                     ENDIF
                     dt_this = dt_all(irivsys(i))
                     IF (.not. is_built_resv(i)) THEN
                        IF (visible_hflux > 0._r8) dt_this = min(dt_this, volwater_ucat(i) / visible_hflux)
                        IF (protected_hflux > 0._r8) THEN
                           protected_avail = 0._r8
                           IF (DEF_USE_LEVEE .and. has_levee(i)) protected_avail = max(levsto(i), 0._r8)
                           dt_this = min(dt_this, protected_avail / protected_hflux)
                        ENDIF
	                     ELSE
	                        IF (sum_hflux_riv(i) > 0._r8) THEN
	                           volwater = volresv(ucat2resv(i))
	                           dt_this = min(dt_this, volwater / sum_hflux_riv(i))
	                        ENDIF
	                     ENDIF
                     IF (dt_this < dt_all(irivsys(i)) - max(1.e-9_r8, 1.e-12_r8 * dt_all(irivsys(i)))) THEN
                        dt_all(irivsys(i)) = dt_this
                        dt_changed = .true.
                     ENDIF
                  ENDDO

                  CALL sync_global_routing_dt(dt_res, dt_all)

                  IF (any(dt_res > 0._r8)) THEN
                     dt_after_sync = maxval(dt_all, mask = dt_res > 0._r8)
                  ELSE
                     dt_after_sync = 0._r8
                  ENDIF
                  IF (abs(dt_after_sync - dt_before_sync) > &
                      max(1.e-9_r8, 1.e-12_r8 * max(abs(dt_before_sync), abs(dt_after_sync)))) THEN
                     dt_changed = .true.
                  ENDIF
                  IF (.not. dt_changed) EXIT
               ENDDO

               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_base
                  END WHERE
               ENDIF

	               IF (allocated(volresv)) THEN
	                  allocate(volresv_safe(size(volresv))); volresv_safe = volresv
	               ELSE
	                  allocate(volresv_safe(0))
	               ENDIF
	               CALL bifurcation_calc(wdsrf_ucat, volresv_safe, is_built_resv, dt_all, &
	                  irivsys, ucatfilter, update_state = .true.)
	               deallocate(volresv_safe)

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
	               IF (allocated(volresv)) THEN
	                  allocate(volresv_safe_trc(size(volresv))); volresv_safe_trc = volresv
	               ELSE
	                  allocate(volresv_safe_trc(0))
	               ENDIF
	               IF (allocated(ucat2resv)) THEN
	                  allocate(ucat2resv_safe_trc(size(ucat2resv))); ucat2resv_safe_trc = ucat2resv
	               ELSE
	                  allocate(ucat2resv_safe_trc(0))
	               ENDIF
	               IF (DEF_USE_BIFURCATION .and. allocated(bif_hflux_lev)) THEN
	                  CALL tracer_substep (acctime_rnof, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
	                     wdsrf_ucat, ucatfilter, &
	                     volresv_safe_trc, ucat2resv_safe_trc, is_built_resv, &
	                     .true., bif_hflux_lev, npthout_local)
	               ELSE
	                  CALL tracer_substep (acctime_rnof, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
	                     wdsrf_ucat, ucatfilter, &
	                     volresv_safe_trc, ucat2resv_safe_trc, is_built_resv, &
	                     .false., reshape((/0._r8/), (/1,1/)), 0)
	               ENDIF
	               deallocate(volresv_safe_trc, ucat2resv_safe_trc)
#endif

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

	               IF (.not. is_built_resv(i)) THEN
	                  volwater = volwater_ucat(i)
	               ELSE
                  volwater = volresv(ucat2resv(i))
               ENDIF

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
	                  levsto(i) = max(levsto(i) - protected_hflux * dt_all(irivsys(i)), 0._r8)
	               ENDIF
#ifdef CoLMDEBUG
               IF (volwater < 0._r8) totalclip = totalclip - volwater
#endif
               volwater = max(volwater, 0.)

               ! Inland depression overflow is a post-transport water correction.
               ! ordinary routing tracer transport is handled earlier in tracer_substep.
               ! only this explicit correction updates tracer locally so the
               ! removed tracer follows the same sink.
               IF (ucat_next(i) == -10) THEN
                  IF (volwater > topo_rivstomax(i)) THEN
                     ! Remove excess water after transport has already run.
                     hflux_fc(i) = (volwater - topo_rivstomax(i)) / dt_all(irivsys(i))
                     ! Remove corresponding tracer proportionally.
                     ! Update trc_flux_out so the unified discharge diagnostic
                     ! at line ~895 (ucat_next <= 0) picks up the correct value.
                     !
#ifdef TRACER
	                     IF (volwater > 1.e-6_r8) THEN
		                        frac_remove = (volwater - topo_rivstomax(i)) / volwater
                              IF (allocated(volresv)) THEN
                                 allocate(volresv_safe_dep(size(volresv))); volresv_safe_dep = volresv
                              ELSE
                                 allocate(volresv_safe_dep(0))
                              ENDIF
                              IF (allocated(ucat2resv)) THEN
                                 allocate(ucat2resv_safe_dep(size(ucat2resv))); ucat2resv_safe_dep = ucat2resv
                              ELSE
                                 allocate(ucat2resv_safe_dep(0))
                              ENDIF
	                           CALL get_cell_volume_dep(i, &
	                              floodplain_curve(i)%depth(topo_rivstomax(i)), &
	                              volresv_safe_dep, ucat2resv_safe_dep, vol_post)
                              deallocate(volresv_safe_dep, ucat2resv_safe_dep)
                           vol_post = max(vol_post, 1.e-6_r8)
                        DO itrc_dep = 1, ntracers
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
                  ! Static levee repartition after transport.  Main-channel
                  ! and bif layer-1 fluxes have already updated the visible
                  ! routing volume; bif layer-2+ flux has already updated
                  ! `levsto` above.  These post-transport compartment
                  ! volumes are the BEFORE state for levee_fldstg's internal
                  ! visible/protected redistribution.
                  vis_vol_bef_lv2 = max(volwater, 0._r8)
                  levsto_bef_lv2  = max(levsto(i), 0._r8)
                  vol_total_levee = volwater + levsto(i)
                  CALL levee_fldstg(i, vol_total_levee, wdsrf_ucat(i), &
                     levsto(i), levdph(i), fldfrc_levee)
                  volwater_ucat(i) = vol_total_levee - levsto(i)
                  levee_floodarea(i) = fldfrc_levee * topo_area(i)
#ifdef TRACER
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
	               IF (allocated(volresv)) THEN
	                  allocate(volresv_safe_diag(size(volresv))); volresv_safe_diag = volresv
	               ELSE
	                  allocate(volresv_safe_diag(0))
	               ENDIF
	               IF (allocated(ucat2resv)) THEN
	                  allocate(ucat2resv_safe_diag(size(ucat2resv))); ucat2resv_safe_diag = ucat2resv
	               ELSE
	                  allocate(ucat2resv_safe_diag(0))
	               ENDIF
	               CALL tracer_diag_accumulate_substep (dt_all, irivsys, ucatfilter, wdsrf_ucat, &
	                  volresv_safe_diag, ucat2resv_safe_diag)
	               deallocate(volresv_safe_diag, ucat2resv_safe_diag)
	            END IF
#endif

            DO i = 1, numucat
               IF (ucatfilter(i)) THEN

#ifdef CoLMDEBUG
                  IF (ucat_next(i) <= 0) THEN
                     totaldis = totaldis + hflux_fc(i)*dt_all(irivsys(i))
                     ! Accumulate tracer discharge at river mouth
#ifdef TRACER
                     IF (ntracers > 0) THEN
                        DO itrc = 1, ntracers
                           trc_mass_dis(itrc) = trc_mass_dis(itrc) &
                              + trc_flux_out(itrc,i)*dt_all(irivsys(i))
                        ENDDO
                     END IF
#endif
                  ENDIF
#endif

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
                  ! rivsto = rivare * wdsrf (matches CaMa-Flood: P2RIVSTO = RIVLEN*RIVWTH*RIVDPH)
                  ! fldsto = total_volume - rivsto
                  ! flddph = max(wdsrf - rivhgt, 0) (depth above channel banks)
                  ! storge = total_volume (+ levsto if levee enabled)
                  ! sfcelv = bed_elevation + wdsrf (matches CaMa: D2RIVELV + D2RIVDPH)
	                  IF (is_built_resv(i)) THEN
	                     volwater = volresv(ucat2resv(i))
	                  ELSE
	                     volwater = volwater_ucat(i)
	                  ENDIF
                  a_rivsto(i) = a_rivsto(i) + floodplain_curve(i)%rivare * wdsrf_ucat(i) * dt_all(irivsys(i))
                  a_fldsto(i) = a_fldsto(i) &
                     + max(volwater - floodplain_curve(i)%rivare * wdsrf_ucat(i), 0._r8) * dt_all(irivsys(i))
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

#ifdef GridRiverLakeSediment
            IF (DEF_USE_SEDIMENT) THEN
               IF (numucat > 0) THEN
                  allocate(floodarea_sed(numucat))
                  DO i = 1, numucat
                     IF (ucatfilter(i)) THEN
                        IF (DEF_USE_LEVEE .and. levee_floodarea(i) > 0.) THEN
                           floodarea_sed(i) = levee_floodarea(i)
                        ELSE
                           floodarea_sed(i) = floodplain_curve(i)%floodarea(wdsrf_ucat(i))
                        ENDIF
                     ELSE
                        floodarea_sed(i) = 0._r8
                     ENDIF
                  ENDDO
               ELSE
                  allocate(floodarea_sed(0))
               ENDIF
               CALL sediment_diag_accumulate(dt_all, irivsys, ucatfilter, &
                  veloc_riv, wdsrf_ucat, hflux_fc, floodarea_sed)
               deallocate(floodarea_sed)
            ENDIF
#endif

            loop_active = any(dt_res > 0)
#ifdef USEMPI
            CALL mpi_allreduce (MPI_IN_PLACE, loop_active, 1, MPI_LOGICAL, &
               MPI_LOR, p_comm_worker, p_err)
#endif

         ENDDO

		         ! Keep restart-visible state consistent after the substep loop.
		         volwater_ucat_valid = .true.

#ifdef CoLMDEBUG
         totalvol_aft = 0.
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
         ! Tracer conservation: compute total mass after routing.
         ! Include protected-side pool so the levee repartition
         ! stays internally closed in the before/after accounting.
#ifdef TRACER
         IF (numucat > 0) THEN
            DO itrc = 1, ntracers
               trc_mass_aft(itrc) = sum(trc_mass(itrc,:)) + sum(trc_inp_buf(itrc,:))
               IF (allocated(trc_levsto)) &
                  trc_mass_aft(itrc) = trc_mass_aft(itrc) + sum(trc_levsto(itrc,:))
            ENDDO
         ELSE
            trc_mass_aft = 0._r8
         END IF
#else
         trc_mass_aft = 0._r8
#endif
#endif
      ENDIF

#ifdef CoLMDEBUG
#ifdef USEMPI
      IF (.not. p_is_worker) ntimestep = 0
      CALL mpi_allreduce (MPI_IN_PLACE, ntimestep, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)

      IF (.not. p_is_worker) totalvol_bef = 0.
      IF (.not. p_is_worker) totalvol_aft = 0.
      IF (.not. p_is_worker) totalrnof    = 0.
      IF (.not. p_is_worker) totaldis     = 0.
      IF (.not. p_is_worker) totalclip    = 0.

      CALL mpi_allreduce (MPI_IN_PLACE, totalvol_bef, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalvol_aft, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalrnof,    1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totaldis,     1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalclip,    1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
#ifdef TRACER
      IF (p_is_worker .and. numucat > 0) THEN
         DO itrc = 1, ntracers
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
         ENDIF
#ifdef TRACER
            DO itrc_dbg = 1, ntracers
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass before  : ', &
                  trc_mass_bef(itrc_dbg), ' kg'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass input   : ', &
                  trc_mass_inp(itrc_dbg), ' kg'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass discharge: ', &
                  trc_mass_dis(itrc_dbg), ' kg'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass reactive : ', &
                  trc_mass_reactive(itrc_dbg), ' kg'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass after   : ', &
                  trc_mass_aft(itrc_dbg), ' kg'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass change  : ', &
                  trc_mass_aft(itrc_dbg) - trc_mass_bef(itrc_dbg), ' kg'
               write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc_dbg, ') mass balance : ', &
                  trc_mass_aft(itrc_dbg) - trc_mass_bef(itrc_dbg) &
                  - trc_mass_inp(itrc_dbg) + trc_mass_dis(itrc_dbg) &
                  - trc_mass_reactive(itrc_dbg), ' kg (should be ~0)'
            ENDDO
#endif
      ENDIF  ! p_is_master

#endif

#ifdef GridRiverLakeSediment
      IF (DEF_USE_SEDIMENT .and. p_is_worker) THEN
         ! All workers must participate (MPI point-to-point inside push_data).
         ! fldfrc is now computed inside grid_sediment_calc from per-routing-period
         ! accumulators (sed_acc_floodarea), not from history-period averages.
         CALL grid_sediment_calc(acctime_rnof)
      ENDIF
#endif

#ifdef TRACER
      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            acc_trc_inp = 0._r8
            acc_rnof_ref = 0._r8
            trc_dry_drain = 0._r8
         ENDIF
      END IF
#endif

      acctime_rnof = 0.

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            acc_rnof_uc = 0.
         ENDIF
      ENDIF

	      ! ---- Publish per-ucat levee floodplain fraction to per-patch
	      !      f_inund_levee_patch so methane scheme 7 can read the latest
	      !      completed routing state.  Because this publish occurs at the
	      !      end of routing, CH4 normally consumes it on the next land step.
	      !      Flow: ucat -> inpm grid (average) -> landpatch (remap).
	      !      Inactive when GridRiverLakeFlow is undef (this whole file
	      !      is gated by that macro).
#if (defined TRACER) && (defined BGC)
      IF (igas_ch4 > 0 .and. allocated(levee_floodarea)) THEN
         CALL publish_levee_fldfrc_to_patches (levee_floodarea)
      ENDIF
      IF (igas_ch4 > 0 .and. allocated(total_floodarea)) THEN
         CALL publish_fldfrc_to_patches (total_floodarea, total_flooddepth)
      ENDIF
#endif

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
      IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
      IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)
      IF (allocated(ucatfilter      )) deallocate(ucatfilter      )
      IF (allocated(levee_floodarea)) deallocate(levee_floodarea)
      IF (allocated(total_floodarea)) deallocate(total_floodarea)
      IF (allocated(total_flooddepth)) deallocate(total_flooddepth)
      IF (allocated(dt_res         )) deallocate(dt_res         )
      IF (allocated(dt_all       )) deallocate(dt_all       )
#ifdef CoLMDEBUG
      IF (allocated(trc_mass_bef)) deallocate(trc_mass_bef, trc_mass_aft, &
                                              trc_mass_inp, trc_mass_dis, trc_mass_reactive)
#endif

   END SUBROUTINE grid_riverlake_flow

#if (defined TRACER) && (defined BGC)
   SUBROUTINE publish_levee_fldfrc_to_patches (levee_floodarea_in)
   !-------------------------------------------------------------------
   ! Push the per-ucat levee floodplain fraction through the inpm grid
   ! back to per-patch f_inund_levee_patch for methane scheme 7.
   !-------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, numinpm, topo_area, &
                                        push_ucat2inpm, remap_patch2inpm
   USE MOD_WorkerPushData, only: worker_push_data, worker_remap_data_grid2pset
   USE MOD_Tracer_Methane_State, only: f_inund_levee_patch
   USE MOD_SPMD_Task

   real(r8), intent(in) :: levee_floodarea_in(:)
   real(r8), allocatable :: fldfrc_uc(:), fldfrc_gd(:)
   integer :: i

      IF (.not. allocated(f_inund_levee_patch)) RETURN
      IF (.not. p_is_worker) RETURN

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
         fillvalue = 0._r8, mode = 'average')

      ! inpm grid -> landpatch (area-weighted remap)
      CALL worker_remap_data_grid2pset (remap_patch2inpm, fldfrc_gd, &
         f_inund_levee_patch, fillvalue = 0._r8, mode = 'average')

      deallocate(fldfrc_uc, fldfrc_gd)

   END SUBROUTINE publish_levee_fldfrc_to_patches

   SUBROUTINE publish_fldfrc_to_patches (total_floodarea_in, total_flooddepth_in)
   !-------------------------------------------------------------------
   ! Publish general flood area and depth for methane scheme 7.
   ! Pushed to per-patch f_inund_flood_patch and f_inund_flood_depth_patch.
   !-------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, numinpm, topo_area, &
                                        push_ucat2inpm, remap_patch2inpm
   USE MOD_WorkerPushData, only: worker_push_data, worker_remap_data_grid2pset
   USE MOD_Tracer_Methane_State, only: f_inund_flood_patch, f_inund_flood_depth_patch
   USE MOD_SPMD_Task

   real(r8), intent(in) :: total_floodarea_in(:)
   real(r8), intent(in) :: total_flooddepth_in(:)
   real(r8), allocatable :: fldfrc_uc(:), fldfrc_gd(:)
   real(r8), allocatable :: flddph_uc(:), flddph_gd(:)
   integer :: i

      IF (.not. allocated(f_inund_flood_patch)) RETURN
      IF (.not. p_is_worker) RETURN

      IF (numucat > 0) THEN
         allocate (fldfrc_uc(numucat), flddph_uc(numucat))
         flddph_uc(:) = 0._r8
         DO i = 1, numucat
            IF (topo_area(i) > 0._r8 .and. i <= size(total_floodarea_in)) THEN
               fldfrc_uc(i) = min(1._r8, max(0._r8, &
                  total_floodarea_in(i) / topo_area(i)))
               IF (i <= size(total_flooddepth_in)) THEN
                  flddph_uc(i) = max(0._r8, total_flooddepth_in(i))
               ELSE
                  flddph_uc(i) = 0._r8
               ENDIF
            ELSE
               fldfrc_uc(i) = 0._r8
               flddph_uc(i) = 0._r8
            ENDIF
         ENDDO
      ELSE
         allocate (fldfrc_uc(0), flddph_uc(0))
      ENDIF

      IF (numinpm > 0) THEN
         allocate (fldfrc_gd(numinpm), flddph_gd(numinpm))
      ELSE
         allocate (fldfrc_gd(0), flddph_gd(0))
      ENDIF
      CALL worker_push_data (push_ucat2inpm, fldfrc_uc, fldfrc_gd, &
         fillvalue = 0._r8, mode = 'average')

      CALL worker_remap_data_grid2pset (remap_patch2inpm, fldfrc_gd, &
         f_inund_flood_patch, fillvalue = 0._r8, mode = 'average')

      IF (allocated(f_inund_flood_depth_patch)) THEN
         CALL worker_push_data (push_ucat2inpm, flddph_uc, flddph_gd, &
            fillvalue = 0._r8, mode = 'average')
         CALL worker_remap_data_grid2pset (remap_patch2inpm, flddph_gd, &
            f_inund_flood_depth_patch, fillvalue = 0._r8, mode = 'average')
      ENDIF

      deallocate(fldfrc_uc, fldfrc_gd, flddph_uc, flddph_gd)

   END SUBROUTINE publish_fldfrc_to_patches
#endif


   SUBROUTINE sync_global_routing_dt(dt_res, dt_all)

      real(r8), intent(in)    :: dt_res(:)
      real(r8), intent(inout) :: dt_all(:)

      real(r8) :: dt_global

      dt_global = huge(1._r8)
      IF (any(dt_res > 0._r8)) THEN
         dt_global = minval(dt_all, mask = dt_res > 0._r8)
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, dt_global, 1, MPI_REAL8, MPI_MIN, &
         p_comm_worker, p_err)
#endif

      IF (dt_global < 0.5_r8 * huge(1._r8)) THEN
         WHERE (dt_res > 0._r8)
            dt_all = dt_global
         ELSEWHERE
            dt_all = 0._r8
         END WHERE
      ELSE
         dt_all = 0._r8
      ENDIF

   END SUBROUTINE sync_global_routing_dt

   ! ---------
   SUBROUTINE grid_riverlake_flow_final ()

      CALL riverlake_network_final ()

      IF (DEF_Reservoir_Method > 0) THEN
         CALL reservoir_final ()
      ENDIF

#ifdef GridRiverLakeSediment
      CALL grid_sediment_final()
#endif

	      CALL levee_final()
      IF (DEF_USE_BIFURCATION) CALL bifurcation_final()
#ifdef TRACER
      CALL tracer_final()
#endif

      ! acc_rnof_uc is owned by MOD_Grid_RiverLakeTimeVars and freed by
      ! deallocate_GridRiverLakeTimeVars; don't deallocate it here.
      IF (allocated(filter_rnof)) deallocate(filter_rnof)

   END SUBROUTINE grid_riverlake_flow_final

END MODULE MOD_Grid_RiverLakeFlow
#endif
