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
   USE MOD_Grid_RiverLakeTracer, only: tracer_init, tracer_init_from_water, &
      tracer_input_from_runoff, &
      tracer_substep, tracer_flush_acc, &
      read_tracer_restart, tracer_final, acc_trc_inp, acc_rnof_ref, trc_mass, trc_inp_buf, trc_flux_out, &
      tracer_refresh_state, tracer_diag_accumulate_substep, &
      trc_levsto, trc_dry_drain, levee_tracer_repartition
   IMPLICIT NONE

   real(r8), parameter :: RIVERMIN  = 1.e-5_r8
   logical, parameter :: dbg_skip_bif_init = .false.
   logical, parameter :: dbg_skip_bif_restart = .false.
   logical, parameter :: dbg_skip_bif_calc = .false.

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
   IMPLICIT NONE

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

      IF (DEF_USE_LEVEE) THEN
         CALL levee_init()
         IF (len_trim(gridriver_restart_file) > 0) THEN
            CALL read_levee_restart(gridriver_restart_file)
         ENDIF
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         IF (dbg_skip_bif_init) THEN
            IF (p_is_master) THEN
               write(*,'(A)') 'DBG bifskip master bifurcation_init skipped'
               call flush(6)
            ENDIF
            IF (p_is_worker .and. p_iam_worker == 0) THEN
               write(*,'(A)') 'DBG bifskip worker0 bifurcation_init skipped'
               call flush(6)
            ENDIF
         ELSE
            CALL bifurcation_init()
         ENDIF
         IF (len_trim(gridriver_restart_file) > 0) THEN
            IF (dbg_skip_bif_restart) THEN
               IF (p_is_master) THEN
                  write(*,'(A)') 'DBG bifskip master read_bifurcation_restart skipped'
                  call flush(6)
               ENDIF
               IF (p_is_worker .and. p_iam_worker == 0) THEN
                  write(*,'(A)') 'DBG bifskip worker0 read_bifurcation_restart skipped'
                  call flush(6)
               ENDIF
            ELSE
               CALL read_bifurcation_restart(gridriver_restart_file)
            ENDIF
         ENDIF
      ENDIF

      IF (DEF_USE_TRACER) THEN
         BLOCK
         USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, volresv
         USE MOD_Grid_Reservoir,         only: ucat2resv
         USE MOD_Tracer_Defs,            only: ntracers
         logical :: trc_restart_found
         logical, allocatable :: trc_missing(:)
         real(r8), allocatable :: wdsrf_safe(:), volresv_safe(:)
         integer,  allocatable :: ucat2resv_safe(:)
         trc_restart_found = .false.
         CALL tracer_init()
         allocate(trc_missing(max(ntracers, 1)))
         trc_missing = .true.   ! assume every tracer missing if no restart
         IF (len_trim(gridriver_restart_file) > 0) THEN
            CALL read_tracer_restart(gridriver_restart_file, trc_restart_found, trc_missing)
         ENDIF
         ! Cold start per-tracer: restart failures (no file, or specific
         ! tracer variables absent) fall back to init-from-water only for
         ! the missing tracers, preserving any that loaded successfully.
         ! wdsrf_ucat / volresv are populated by READ_GridRiverLakeTimeVars
         ! before this routine runs. volresv and ucat2resv are unallocated
         ! on reservoir-free workers, so wrap them into size-0 proxies
         ! (same pattern as tracer_substep).
         IF (any(trc_missing(1:ntracers))) THEN
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
            CALL tracer_init_from_water(wdsrf_safe, volresv_safe, ucat2resv_safe, &
               trc_missing(1:ntracers))
            deallocate(wdsrf_safe, volresv_safe, ucat2resv_safe)
         ENDIF
         deallocate(trc_missing)
         END BLOCK
      ENDIF

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
   USE MOD_Tracer_Defs,    only: ntracers
   USE MOD_Tracer_Vars,    only: trc_rnof_step
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
   real(r8), allocatable :: sum_mflux_riv(:)
   real(r8), allocatable :: sum_zgrad_riv(:)

   real(r8) :: veloct_fc, height_fc, momen_fc, zsurf_fc
   real(r8) :: bedelv_fc, height_up, height_dn
   real(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: volwater, friction, floodarea
   real(r8) :: vol_total_levee, fldfrc_levee
   real(r8), allocatable :: levee_floodarea(:)
   real(r8),  allocatable :: dt_res(:), dt_all(:)
   logical,   allocatable :: ucatfilter(:)
   real(r8),  allocatable :: hflux_avg_trc(:)
   logical,   allocatable :: trc_filter(:)
   real(r8),  allocatable :: volresv_trc(:)
   integer,   allocatable :: ucat2resv_trc(:)
   real(r8),  allocatable :: trc_acc_discharge(:) ! Per-routing-period discharge accumulator for tracer
   real(r8),  allocatable :: trc_acctime(:)       ! Per-routing-period time accumulator for tracer
   real(r8),  allocatable :: trc_acc_bifflw_lev(:,:) ! Per-routing-period bif pathway flux for tracer
   real(r8),  allocatable :: trc_acc_bifflw_time(:)  ! Per-routing-period bif pathway time for tracer
   logical :: loop_active
#ifdef CoLMDEBUG
   real(r8) :: totalvol_bef, totalvol_aft, totalrnof, totaldis
   real(r8) :: trc_mass_bef = 0, trc_mass_aft = 0, trc_mass_inp = 0, trc_mass_dis = 0
   real(r8) :: bif_flux_sum_total, bif_flux_sum_max
   integer  :: itrc_dbg
#endif


      IF (p_is_worker) THEN
         allocate (rnof_gd (numinpm))
         allocate (rnof_uc (numucat))
         IF (DEF_USE_TRACER .and. ntracers > 0) THEN
            allocate (trc_rnof_gd (ntracers, numinpm))
            allocate (trc_rnof_uc (ntracers, numucat))
            trc_rnof_gd = 0._r8
            trc_rnof_uc = 0._r8
         ENDIF

         CALL worker_remap_data_pset2grid (remap_patch2inpm, rnof, rnof_gd, &
            fillvalue = 0., filter = filter_rnof)

         IF (numinpm > 0) THEN
            WHERE (push_ucat2inpm%sum_area > 0)
               rnof_gd = rnof_gd / push_ucat2inpm%sum_area
            END WHERE
         ENDIF

         CALL worker_push_data (push_inpm2ucat, rnof_gd, rnof_uc, &
            fillvalue = 0., mode = 'sum')

         IF (DEF_USE_TRACER .and. ntracers > 0) THEN
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
#ifdef CoLMDEBUG
            BLOCK
            integer :: i, imax_gd, imax_uc
            real(r8) :: ratio_gd, ratio_uc, ratio_gd_max, ratio_uc_max

            ratio_gd_max = 0._r8
            ratio_uc_max = 0._r8
            imax_gd = 0
            imax_uc = 0

            DO i = 1, numinpm
               ratio_gd = trc_rnof_gd(1, i) / max(rnof_gd(i) * deltime, 1.e-30_r8)
               IF (ratio_gd > ratio_gd_max) THEN
                  ratio_gd_max = ratio_gd
                  imax_gd = i
               ENDIF
            ENDDO

            DO i = 1, numucat
               ratio_uc = trc_rnof_uc(1, i) / max(rnof_uc(i) * deltime, 1.e-30_r8)
               IF (ratio_uc > ratio_uc_max) THEN
                  ratio_uc_max = ratio_uc
                  imax_uc = i
               ENDIF
            ENDDO

            IF (ratio_gd_max > 5.e-3_r8 .and. imax_gd > 0) THEN
               WRITE(*,'(A,I8,A,E10.3,A,E10.3,A,E10.3)') &
                  ' DBG_TRCRMAP_GD: gdid=', inpm_gdid(imax_gd), &
                  ' rnof=', rnof_gd(imax_gd), &
                  ' trc_rnof=', trc_rnof_gd(1, imax_gd), &
                  ' ratio=', ratio_gd_max
            ENDIF

            IF (ratio_uc_max > 5.e-3_r8 .and. imax_uc > 0) THEN
               WRITE(*,'(A,I8,A,I8,A,E10.3,A,E10.3,A,E10.3)') &
                  ' DBG_TRCRMAP_UC: ucat=', imax_uc, &
                  ' next=', ucat_next(imax_uc), &
                  ' rnof=', rnof_uc(imax_uc), &
                  ' trc_rnof=', trc_rnof_uc(1, imax_uc), &
                  ' ratio=', ratio_uc_max
            ENDIF
            END BLOCK
#endif
         ENDIF

         IF (numucat > 0) THEN
            acc_rnof_uc = acc_rnof_uc + rnof_uc*1.e-3*deltime

            ! Accumulate tracer input associated with this runoff increment
            IF (DEF_USE_TRACER) THEN
               IF (ntracers > 0) THEN
                  ! trc_rnof_uc is in R×mm (from land: rsur*ratio*dt, rsur in mm/s)
                  ! rnof_uc * 1.e-3 * deltime is in m (depth, not volume)
                  ! Need to convert trc_rnof_uc from mm to m to match water units
                  CALL tracer_input_from_runoff(rnof_uc*1.e-3*deltime, numucat, trc_rnof_uc*1.e-3)
               ELSE
                  CALL tracer_input_from_runoff(rnof_uc*1.e-3*deltime, numucat)
               ENDIF
            ENDIF
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
            allocate (sum_mflux_riv (numucat))
            allocate (sum_zgrad_riv (numucat))
            allocate (ucatfilter    (numucat))

            IF (DEF_USE_LEVEE) THEN
               allocate (levee_floodarea (numucat))
               levee_floodarea = 0.
            ENDIF

            allocate (hflux_sumups  (numucat))
            allocate (mflux_sumups  (numucat))
            allocate (zgrad_sumups  (numucat))

            IF (DEF_Reservoir_Method > 0) THEN
               allocate (hflux_resv (numucat))
               allocate (mflux_resv (numucat))
            ENDIF

            allocate (dt_res (numrivsys))
            allocate (dt_all (numrivsys))

            IF (DEF_USE_TRACER) THEN
               allocate (trc_acc_discharge (numucat))
               allocate (trc_acctime       (numucat))
               trc_acc_discharge = 0._r8
               trc_acctime       = 0._r8
               ! Always allocate (zero-length if no bifurcation) for safe passing
               IF (DEF_USE_BIFURCATION) THEN
                  allocate (trc_acc_bifflw_lev  (npthlev_bif, npthout_local))
                  allocate (trc_acc_bifflw_time (npthout_local))
               ELSE
                  allocate (trc_acc_bifflw_lev  (0, 0))
                  allocate (trc_acc_bifflw_time (0))
               ENDIF
               trc_acc_bifflw_lev  = 0._r8
               trc_acc_bifflw_time = 0._r8
            ENDIF

         ENDIF

#ifdef CoLMDEBUG
         totalrnof = sum(acc_rnof_uc)
         totalvol_bef = 0.
#endif

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
               IF (DEF_USE_LEVEE .and. has_levee(i) .and. volwater_ucat_valid) THEN
                  ! Persistent tracked volume (restored from restart or previous call)
                  volwater = volwater_ucat(i)
               ELSEIF (DEF_USE_LEVEE .and. has_levee(i) .and. levsto(i) > 0._r8) THEN
                  ! Old restart without volwater_ucat: reconstruct per case.
                  ! Case-3: wdsrf pinned at levee crest, volume(wdsrf) includes protected side
                  !   → non-protected volume = levee_topsto (exact)
                  ! Case-4: volume(wdsrf) = vol_total → volwater = volume(wdsrf) - levsto
                  IF (wdsrf_ucat(i) <= floodplain_curve(i)%rivhgt + levee_hgt(i) + 1.e-6_r8) THEN
                     volwater = levee_topsto(i)
                  ELSE
                     volwater = floodplain_curve(i)%volume(wdsrf_ucat(i)) - levsto(i)
                  ENDIF
                  volwater = max(volwater, 0._r8)
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
                  ! H1 fix: levee repartition now also moves tracer mass
                  ! between visible (trc_mass) and protected (trc_levsto)
                  ! pools in lockstep with the water transfer, so mass
                  ! no longer stays pinned to the visible side while
                  ! water crosses the levee.
                  BLOCK
                  real(r8) :: vis_vol_bef_lv, levsto_bef_lv
                  vis_vol_bef_lv = volwater
                  levsto_bef_lv  = levsto(i)
                  vol_total_levee = volwater + levsto(i)
                  CALL levee_fldstg(i, vol_total_levee, wdsrf_ucat(i), &
                     levsto(i), levdph(i), fldfrc_levee)
                  volwater_ucat(i) = vol_total_levee - levsto(i)
                  IF (DEF_USE_TRACER) THEN
                     CALL levee_tracer_repartition(i, &
                        vis_vol_bef_lv, levsto_bef_lv, &
                        volwater_ucat(i), levsto(i))
                  ENDIF
                  END BLOCK
               ELSE
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

         ntimestep = 0
#ifdef CoLMDEBUG
         totaldis  = 0.
         bif_flux_sum_total = 0._r8
         bif_flux_sum_max   = 0._r8
         ! Tracer conservation: save total mass BEFORE input addition.
         ! H1 fix: include protected-side pool so mass trapped behind a
         ! levee is not treated as "missing" by the before/after check.
         IF (DEF_USE_TRACER .and. numucat > 0) THEN
            trc_mass_bef = sum(trc_mass(1,:)) + sum(trc_inp_buf(1,:))
            IF (allocated(trc_levsto)) trc_mass_bef = trc_mass_bef + sum(trc_levsto(1,:))
            trc_mass_inp = sum(acc_trc_inp(1,:))
            trc_mass_dis = 0._r8
         ELSE
            trc_mass_bef = 0._r8
            trc_mass_inp = 0._r8
            trc_mass_dis = 0._r8
         ENDIF
#endif

         dt_res(:) = acctime_rnof

         ! When bifurcation is on, pathways can cross river systems. Each
         ! worker's local dt_res is per-system and can reach 0 at different
         ! iterations. Without synchronisation, a worker that finishes early
         ! exits this loop and stops calling worker_push_data while other
         ! workers still need that exchange — deadlocking the isend/irecv
         ! pairs. Fix: all workers keep looping as long as ANY worker has
         ! dt_res > 0. "Done" workers have dt_all = 0, ucatfilter = false,
         ! and their loop bodies are no-ops except for the collective pushes.
         loop_active = any(dt_res > 0)
#ifdef USEMPI
         IF (DEF_USE_BIFURCATION) THEN
            CALL mpi_allreduce (MPI_IN_PLACE, loop_active, 1, MPI_LOGICAL, &
               MPI_LOR, p_comm_worker, p_err)
         ENDIF
#endif

         DO WHILE (loop_active)

            ntimestep = ntimestep + 1

            CALL worker_push_data (push_next2ucat, wdsrf_ucat, wdsrf_next, fillvalue = spval)
            ! velocity in ocean or inland depression is assumed to be 0.
            CALL worker_push_data (push_next2ucat, veloc_riv,  veloc_next, fillvalue = 0.)

            dt_all(:) = min(dt_res(:), 60.)

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

            ! ----- Bifurcation pathways -----
            IF (DEF_USE_BIFURCATION .and. .not. dbg_skip_bif_calc) THEN
               CALL bifurcation_calc(wdsrf_ucat, volresv, is_built_resv, dt_all, irivsys, ucatfilter)
               IF (allocated(a_bifflw_lev) .and. allocated(a_bifflw_acctime)) THEN
                  DO ipth = 1, npthout_local
                     i_up = pth_upst_local(ipth)
                     IF (i_up < 1 .or. i_up > numucat) CYCLE
                     IF (.not. ucatfilter(i_up)) CYCLE
                     a_bifflw_lev(:, ipth) = a_bifflw_lev(:, ipth) &
                        + bif_hflux_lev(:, ipth) * dt_all(irivsys(i_up))
                     a_bifflw_acctime(ipth) = a_bifflw_acctime(ipth) + dt_all(irivsys(i_up))
                     IF (DEF_USE_TRACER) THEN
                        trc_acc_bifflw_lev(:, ipth) = trc_acc_bifflw_lev(:, ipth) &
                           + bif_hflux_lev(:, ipth) * dt_all(irivsys(i_up))
                        trc_acc_bifflw_time(ipth) = trc_acc_bifflw_time(ipth) + dt_all(irivsys(i_up))
                     ENDIF
                  ENDDO
               ENDIF
               ! Add bifurcation volume flux to main routing (volume only, not momentum)
               IF (numucat > 0) THEN
                  WHERE (ucatfilter)
                     sum_hflux_riv = sum_hflux_riv + bif_hflux_sum
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

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

               dt_this = dt_all(irivsys(i))

               ! constraint 1: CFL condition (only for rivers)
               IF (.not. is_built_resv(i)) THEN
                  IF ((veloc_riv(i) /= 0.) .or. (wdsrf_ucat(i) > 0.)) THEN
                     dt_this = min(dt_this, topo_rivlen(i)/(abs(veloc_riv(i))+sqrt(grav*wdsrf_ucat(i)))*0.8)
                  ENDIF
               ENDIF

               ! constraint 2: Avoid negative values of water
               IF (sum_hflux_riv(i) > 0) THEN
                  IF (.not. is_built_resv(i)) THEN
                     ! for river or lake catchment
                     IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                        volwater = volwater_ucat(i)
                     ELSE
                        volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
                     ENDIF
                  ELSE
                     ! for reservoir
                     volwater = volresv(ucat2resv(i))
                  ENDIF

                  dt_this = min(dt_this, volwater / sum_hflux_riv(i))

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

#ifdef USEMPI
            IF (rivsys_by_multiple_procs) THEN
               CALL mpi_allreduce (MPI_IN_PLACE, dt_all, 1, MPI_REAL8, MPI_MIN, p_comm_rivsys, p_err)
            ENDIF
#endif

            ! Per-sub-step tracer transport: advance tracer in lockstep
            ! with water BEFORE the water state update so concentration
            ! is computed from the pre-update volume.
            IF (DEF_USE_TRACER) THEN
               IF (DEF_USE_BIFURCATION .and. .not. dbg_skip_bif_calc &
                  .and. allocated(bif_hflux_lev)) THEN
                  CALL tracer_substep (acctime_rnof, dt_all, irivsys, hflux_fc, wdsrf_ucat, ucatfilter, &
                     volresv, ucat2resv, is_built_resv, &
                     .true., bif_hflux_lev, npthout_local)
               ELSE
                  CALL tracer_substep (acctime_rnof, dt_all, irivsys, hflux_fc, wdsrf_ucat, ucatfilter, &
                     volresv, ucat2resv, is_built_resv, &
                     .false., reshape((/0._r8/), (/1,1/)), 0)
               ENDIF
            ENDIF

            DO i = 1, numucat

               IF (.not. ucatfilter(i)) CYCLE

               IF (.not. is_built_resv(i)) THEN
                  IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                     volwater = volwater_ucat(i)
                  ELSE
                     volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
                  ENDIF
               ELSE
                  volwater = volresv(ucat2resv(i))
               ENDIF

               volwater = volwater - sum_hflux_riv(i) * dt_all(irivsys(i))
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
                     ! H3: compute the post-clamp cell volume via the same
                     ! helper the rest of tracer uses (reservoir / leveed /
                     ! plain floodplain-curve branches). topo_rivstomax only
                     ! matches get_cell_volume in the plain depression case;
                     ! leveed depressions carry water inside volwater_ucat,
                     ! so using topo_rivstomax there would make trc_conc_dep
                     ! transiently disagree with the immediate follow-up in
                     ! tracer_refresh_state (MOD_Grid_RiverLakeTracer:457).
                     IF (DEF_USE_TRACER .and. volwater > 1.e-6_r8) THEN
                        BLOCK
                        USE MOD_Grid_RiverLakeTracer, only: &
                           trc_conc_dep => trc_conc, get_cell_volume_dep => get_cell_volume
                        integer :: itrc_dep
                        real(r8) :: frac_remove, trc_removed
                        real(r8) :: vol_post
                        frac_remove = (volwater - topo_rivstomax(i)) / volwater
                        ! Post-clamp wdsrf_ucat is set at L950 from
                        ! floodplain_curve(...)%depth(topo_rivstomax). Mirror
                        ! that projection to recover the volume the tracer
                        ! conc denominator should use.
                        CALL get_cell_volume_dep(i, &
                           floodplain_curve(i)%depth(topo_rivstomax(i)), &
                           volresv, ucat2resv, vol_post)
                        vol_post = max(vol_post, 1.e-6_r8)
                        DO itrc_dep = 1, ntracers
                           trc_removed = trc_mass(itrc_dep, i) * frac_remove
                           trc_mass(itrc_dep, i) = trc_mass(itrc_dep, i) - trc_removed
                           trc_flux_out(itrc_dep, i) = trc_removed / dt_all(irivsys(i))
                           trc_conc_dep(itrc_dep, i) = trc_mass(itrc_dep, i) / vol_post
                        ENDDO
                        END BLOCK
                     ENDIF
                     volwater = topo_rivstomax(i)
                  ENDIF
               ENDIF

               IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i))) THEN
                  ! H1 fix: repartition tracer mass between visible and
                  ! protected-side pools alongside the water transfer.
                  BLOCK
                  real(r8) :: vis_vol_bef_lv2, levsto_bef_lv2
                  vis_vol_bef_lv2 = volwater
                  levsto_bef_lv2  = levsto(i)
                  vol_total_levee = volwater + levsto(i)
                  CALL levee_fldstg(i, vol_total_levee, wdsrf_ucat(i), &
                     levsto(i), levdph(i), fldfrc_levee)
                  volwater_ucat(i) = vol_total_levee - levsto(i)
                  levee_floodarea(i) = fldfrc_levee * topo_area(i)
                  IF (DEF_USE_TRACER) THEN
                     CALL levee_tracer_repartition(i, &
                        vis_vol_bef_lv2, levsto_bef_lv2, &
                        volwater_ucat(i), levsto(i))
                  ENDIF
                  END BLOCK
               ELSE
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

            IF (DEF_USE_TRACER .and. p_is_worker .and. numucat > 0) THEN
               CALL tracer_diag_accumulate_substep (dt_all, irivsys, ucatfilter, wdsrf_ucat, &
                  volresv, ucat2resv)
            ENDIF

            DO i = 1, numucat
               IF (ucatfilter(i)) THEN

#ifdef CoLMDEBUG
                  IF (ucat_next(i) <= 0) THEN
                     totaldis = totaldis + hflux_fc(i)*dt_all(irivsys(i))
                     ! Accumulate tracer discharge at river mouth
                     IF (DEF_USE_TRACER .and. ntracers > 0) THEN
                        trc_mass_dis = trc_mass_dis + trc_flux_out(1,i)*dt_all(irivsys(i))
                     ENDIF
                  ENDIF
#endif

                  acctime_ucat(i) = acctime_ucat(i) + dt_all(irivsys(i))

                  a_wdsrf_ucat(i) = a_wdsrf_ucat(i) + wdsrf_ucat(i) * dt_all(irivsys(i))
                  a_veloc_riv (i) = a_veloc_riv (i) + veloc_riv (i) * dt_all(irivsys(i))
                  a_discharge (i) = a_discharge (i) + hflux_fc  (i) * dt_all(irivsys(i))

                  IF (DEF_USE_TRACER) THEN
                     trc_acc_discharge(i) = trc_acc_discharge(i) + hflux_fc(i) * dt_all(irivsys(i))
                     trc_acctime(i) = trc_acctime(i) + dt_all(irivsys(i))
                  ENDIF

                  IF (DEF_USE_LEVEE .and. levee_floodarea(i) > 0.) THEN
                     floodarea = levee_floodarea(i)
                  ELSE
                     floodarea = floodplain_curve(i)%floodarea (wdsrf_ucat(i))
                  ENDIF
                  a_floodarea (i) = a_floodarea (i) + floodarea * dt_all(irivsys(i))

                  ! River/floodplain storage separation, total storage, surface elevation
                  ! rivsto = rivare * wdsrf (matches CaMa-Flood: P2RIVSTO = RIVLEN*RIVWTH*RIVDPH)
                  ! fldsto = total_volume - rivsto
                  ! flddph = max(wdsrf - rivhgt, 0) (depth above channel banks)
                  ! storge = total_volume (+ levsto if levee enabled)
                  ! sfcelv = bed_elevation + wdsrf (matches CaMa: D2RIVELV + D2RIVDPH)
                  IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                     volwater = volwater_ucat(i)
                  ELSE
                     volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
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

                  IF (DEF_USE_BIFURCATION .and. .not. dbg_skip_bif_calc) THEN
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
            IF (DEF_USE_BIFURCATION .or. DEF_USE_TRACER) THEN
               CALL mpi_allreduce (MPI_IN_PLACE, loop_active, 1, MPI_LOGICAL, &
                  MPI_LOR, p_comm_worker, p_err)
            ENDIF
#endif

         ENDDO

         ! After first routing call, volwater_ucat is valid for all levee cells
         IF (DEF_USE_LEVEE) volwater_ucat_valid = .true.

#ifdef CoLMDEBUG
         totalvol_aft = 0.
         DO i = 1, numucat
            IF (.not. is_built_resv(i)) THEN
               IF (DEF_USE_LEVEE .and. has_levee(i)) THEN
                  volwater = volwater_ucat(i) + levsto(i)
               ELSE
                  volwater = floodplain_curve(i)%volume (wdsrf_ucat(i))
               ENDIF
               totalvol_aft = totalvol_aft + volwater
            ELSE
               totalvol_aft = totalvol_aft + volresv(ucat2resv(i))
            ENDIF
         ENDDO
         ! Tracer conservation: compute total mass after routing.
         ! H1 fix: include protected-side pool so the levee repartition
         ! stays internally closed in the before/after accounting.
         IF (DEF_USE_TRACER .and. numucat > 0) THEN
            trc_mass_aft = sum(trc_mass(1,:)) + sum(trc_inp_buf(1,:))
            IF (allocated(trc_levsto)) trc_mass_aft = trc_mass_aft + sum(trc_levsto(1,:))
         ELSE
            trc_mass_aft = 0._r8
         ENDIF
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

      CALL mpi_allreduce (MPI_IN_PLACE, totalvol_bef, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalvol_aft, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totalrnof,    1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, totaldis,     1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (DEF_USE_TRACER .and. p_is_worker .and. numucat > 0) THEN
         trc_mass_dis = trc_mass_dis + sum(trc_dry_drain(1,:))
      ENDIF
      IF (.not. p_is_worker) THEN
         bif_flux_sum_total = 0._r8
         bif_flux_sum_max   = 0._r8
         trc_mass_bef = 0._r8
         trc_mass_aft = 0._r8
         trc_mass_inp = 0._r8
         trc_mass_dis = 0._r8
      ENDIF
      CALL mpi_allreduce (MPI_IN_PLACE, bif_flux_sum_total, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, bif_flux_sum_max,   1, MPI_REAL8, MPI_MAX, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_bef, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_aft, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_inp, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_dis, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
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
         IF (DEF_USE_BIFURCATION) THEN
            write(*,'(A,ES10.2,A)') 'Bif global net flux (cumul)     : ', bif_flux_sum_total, ' m^3/s (should be ~0)'
            write(*,'(A,ES10.2,A)') 'Bif global net flux (max/step)  : ', bif_flux_sum_max,   ' m^3/s (should be ~0)'
         ENDIF
         IF (DEF_USE_TRACER) THEN
            write(*,'(A,ES12.4,A)') 'Tracer(1) mass before  : ', trc_mass_bef, ' kg'
            write(*,'(A,ES12.4,A)') 'Tracer(1) mass input   : ', trc_mass_inp, ' kg'
            write(*,'(A,ES12.4,A)') 'Tracer(1) mass discharge:', trc_mass_dis, ' kg'
            write(*,'(A,ES12.4,A)') 'Tracer(1) mass after   : ', trc_mass_aft, ' kg'
            write(*,'(A,ES12.4,A)') 'Tracer(1) mass change  : ', trc_mass_aft - trc_mass_bef, ' kg'
            write(*,'(A,ES12.4,A)') 'Tracer(1) mass balance : ', &
               trc_mass_aft - trc_mass_bef - trc_mass_inp + trc_mass_dis, ' kg (should be ~0)'
         ENDIF
      ENDIF  ! p_is_master

      ! Max concentration diagnostic — runs on each worker, prints locally
      IF (DEF_USE_TRACER .and. p_is_worker .and. numucat > 0) THEN
         BLOCK
         USE MOD_Grid_RiverLakeTracer, only: trc_conc_dbg => trc_conc, trc_v_dry_off_dbg => trc_v_dry_off
         integer :: i, imax, num_dry_with_mass
         real(r8) :: volw_i, volw_max, conc_max, inp_mass_ratio, inp_out_ratio, inp_rnof_ratio
         real(r8) :: dry_mass_total, dry_mass_max
         conc_max = 0._r8
         imax = 0
         dry_mass_total = 0._r8
         dry_mass_max = 0._r8
         num_dry_with_mass = 0
         DO i = 1, numucat
            IF (lake_type(i) == 2 .and. is_built_resv(i) .and. size(volresv) > 0) THEN
               volw_i = volresv(ucat2resv(i))
            ELSEIF (DEF_USE_LEVEE .and. has_levee(i)) THEN
               volw_i = volwater_ucat(i)
            ELSE
               volw_i = floodplain_curve(i)%volume(wdsrf_ucat(i))
            ENDIF

            IF (DEF_USE_LEVEE .and. has_levee(i) .and. (.not. is_built_resv(i)) .and. &
                volwater_ucat(i) <= trc_v_dry_off_dbg .and. levsto(i) > trc_v_dry_off_dbg) THEN
               ! H1 fix: trc_levsto now holds the protected-side pool, so
               ! visible-dry + levsto>0 is no longer a mass-orphan case.
               ! Keep the guard but warn only when the visible-side still
               ! carries residual mass (would indicate levee_tracer_repartition
               ! missed the transfer — a regression indicator).
               IF (abs(trc_mass(1, i)) > 1.e-9_r8) THEN
                  WRITE(*,'(A,I8,A,E10.3,A,E10.3,A,E10.3,A,E10.3)') &
                     ' DBG_LEVTRC: ucat=', i, &
                     ' vis_vol=', volwater_ucat(i), &
                     ' levsto=', levsto(i), &
                     ' trc_mass=', trc_mass(1, i), &
                     ' trc_levsto=', trc_levsto(1, i)
                  WRITE(*,'(A)') &
                     ' DBG_LEVTRC: visible-side tracer residual on dry+levsto cell — check repartition'
               ENDIF
            ENDIF

            IF (volw_i <= trc_v_dry_off_dbg) THEN
               IF (abs(trc_mass(1, i)) > 1.e-30_r8) THEN
                  dry_mass_total = dry_mass_total + trc_mass(1, i)
                  dry_mass_max = max(dry_mass_max, abs(trc_mass(1, i)))
                  num_dry_with_mass = num_dry_with_mass + 1
               ENDIF
               CYCLE
            ENDIF

            IF (trc_conc_dbg(1, i) > conc_max) THEN
               conc_max = trc_conc_dbg(1, i)
               imax = i
            ENDIF
         ENDDO
         IF (num_dry_with_mass > 0) THEN
            WRITE(*,'(A,E10.3,A,I8,A,E10.3)') &
               ' DBG_DRYTRC: dry_mass_total=', dry_mass_total, &
               ' num_dry_with_mass=', num_dry_with_mass, &
               ' dry_mass_max=', dry_mass_max
         ENDIF
         IF (imax > 0 .and. conc_max > 5.0e-3_r8) THEN  ! only print if notably above R_input (~2e-3)
            IF (lake_type(imax) == 2 .and. is_built_resv(imax) .and. size(volresv) > 0) THEN
               volw_max = volresv(ucat2resv(imax))
            ELSEIF (DEF_USE_LEVEE .and. has_levee(imax)) THEN
               volw_max = volwater_ucat(imax)
            ELSE
               volw_max = floodplain_curve(imax)%volume(wdsrf_ucat(imax))
            ENDIF
            inp_rnof_ratio = acc_trc_inp(1,imax) / max(acc_rnof_ref(imax), 1.e-30_r8)
            inp_mass_ratio = acc_trc_inp(1,imax) / max(abs(trc_mass(1,imax)), 1.e-30_r8)
            inp_out_ratio  = acc_trc_inp(1,imax) / max(abs(trc_flux_out(1,imax)), 1.e-30_r8)
            WRITE(*,'(A,I8,A,I8,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3)') &
               ' DBG_MAXCONC: ucat=', imax, &
               ' next=', ucat_next(imax), &
               ' conc=', trc_conc_dbg(1,imax), &
               ' mass=', trc_mass(1,imax), &
               ' vol=', volw_max, &
               ' inp=', acc_trc_inp(1,imax), &
               ' inp_rnof=', inp_rnof_ratio, &
               ' inp_mass=', inp_mass_ratio, &
               ' inp_out=', inp_out_ratio, &
               ' out=', trc_flux_out(1,imax), &
               ' hflux=', hflux_fc(imax), &
               ' wdsrf=', wdsrf_ucat(imax)
         ENDIF
         END BLOCK
      ENDIF
#endif

#ifdef GridRiverLakeSediment
      IF (DEF_USE_SEDIMENT .and. p_is_worker) THEN
         ! All workers must participate (MPI point-to-point inside push_data).
         ! fldfrc is now computed inside grid_sediment_calc from per-routing-period
         ! accumulators (sed_acc_floodarea), not from history-period averages.
         CALL grid_sediment_calc(acctime_rnof)
      ENDIF
#endif

      ! ----- Tracer transport -----
      ! All workers must participate (MPI point-to-point inside worker_push_data).
      IF (DEF_USE_TRACER .and. p_is_worker) THEN
         ! Allocate on all workers (zero-length if numucat=0)
         allocate (hflux_avg_trc (numucat))
         allocate (trc_filter    (numucat))

         IF (numucat > 0) THEN
            ! Use per-routing-period accumulators (not history-window ones)
            WHERE (trc_acctime > 0._r8)
               hflux_avg_trc = trc_acc_discharge / trc_acctime
            ELSE WHERE
               hflux_avg_trc = 0._r8
            END WHERE
            trc_filter = trc_acctime > 0._r8
         ENDIF

         ! Safe wrappers: volresv/ucat2resv may not be allocated on reservoir-free workers
         IF (allocated(volresv)) THEN
            allocate (volresv_trc(size(volresv)));  volresv_trc = volresv
         ELSE
            allocate (volresv_trc(0))
         ENDIF
         IF (allocated(ucat2resv)) THEN
            allocate (ucat2resv_trc(size(ucat2resv)));  ucat2resv_trc = ucat2resv
         ELSE
            allocate (ucat2resv_trc(0))
         ENDIF

         ! tracer_calc is now replaced by per-sub-step tracer_substep
         ! calls inside the DO WHILE loop above. The old single-step
         ! tracer_calc used time-averaged flux which decoupled tracer
         ! from water, causing unphysical behaviour for signed tracers.

         IF (numucat > 0) THEN
            acc_trc_inp = 0._r8
            acc_rnof_ref = 0._r8
            trc_dry_drain = 0._r8
         ENDIF

         deallocate (hflux_avg_trc)
         deallocate (trc_filter)
         deallocate (volresv_trc)
         deallocate (ucat2resv_trc)
         IF (allocated(trc_acc_discharge  )) deallocate(trc_acc_discharge  )
         IF (allocated(trc_acctime        )) deallocate(trc_acctime        )
         IF (allocated(trc_acc_bifflw_lev )) deallocate(trc_acc_bifflw_lev )
         IF (allocated(trc_acc_bifflw_time)) deallocate(trc_acc_bifflw_time)
      ENDIF

      acctime_rnof = 0.

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            acc_rnof_uc = 0.
         ENDIF
      ENDIF

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
      IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
      IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)
      IF (allocated(ucatfilter      )) deallocate(ucatfilter      )
      IF (allocated(levee_floodarea)) deallocate(levee_floodarea)
      IF (allocated(dt_res         )) deallocate(dt_res         )
      IF (allocated(dt_all       )) deallocate(dt_all       )

   END SUBROUTINE grid_riverlake_flow

   ! ---------
   SUBROUTINE grid_riverlake_flow_final ()

      CALL riverlake_network_final ()

      IF (DEF_Reservoir_Method > 0) THEN
         CALL reservoir_final ()
      ENDIF

#ifdef GridRiverLakeSediment
      CALL grid_sediment_final()
#endif

      IF (DEF_USE_LEVEE) CALL levee_final()
      IF (DEF_USE_BIFURCATION) CALL bifurcation_final()
      IF (DEF_USE_TRACER) CALL tracer_final()

      ! acc_rnof_uc is owned by MOD_Grid_RiverLakeTimeVars and freed by
      ! deallocate_GridRiverLakeTimeVars; don't deallocate it here.
      IF (allocated(filter_rnof)) deallocate(filter_rnof)

   END SUBROUTINE grid_riverlake_flow_final

END MODULE MOD_Grid_RiverLakeFlow
#endif
