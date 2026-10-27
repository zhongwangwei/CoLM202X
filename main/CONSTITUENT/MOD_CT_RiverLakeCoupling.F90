#include <define.h>

#if (defined GridRiverLakeFlow) && (defined ConstituentTransport)
MODULE MOD_CT_RiverLakeCoupling
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Constituents in the grid river-lake flow: the calls made from
!   MOD_Grid_RiverLakeFlow (each one marked "#ifdef ConstituentTransport" there).
!
!   The routines read the host water state (passed in, or from the river
!   modules) and change only the constituent pools of MOD_CT_RiverLake and the
!   flood-feedback constituent arrays below.  No host variable is changed, so
!   the host solution is the same as without them.
!
!   Units: river pools are R*m^3 (ratio times water volume); land-side flood
!   arrays are R*mm over the patch, as the land constituent pools.
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_CT_Defs, only: ntracers, ct_uses_land_water_transport, ct_has_dissolved_limit, &
      ct_equilibrate_dissolved, tracer_is_nonvolatile_solute
   USE MOD_CT_RiverLake
   IMPLICIT NONE

   ! ----- grid flood feedback -----
   ! Land side (ntracers, numpatch), read and written in CoLMMAIN:
   !   flood_tracer_credit_patch : constituent of the flood water published to the patch
   !   flood_tracer_evap_patch   : constituent lost (>0) or gained (<0) by flood evaporation
   !   flood_tracer_land_patch   : constituent that entered the soil with flood infiltration
   real(r8), allocatable :: flood_tracer_credit_patch(:,:)
   real(r8), allocatable :: flood_tracer_evap_patch  (:,:)
   real(r8), allocatable :: flood_tracer_land_patch  (:,:)
   ! River side (ntracers, numucat): constituent of the published visible and
   ! protected (levee) flood water of each unit catchment
   real(r8), allocatable :: flood_visible_tracer_uc  (:,:)
   real(r8), allocatable :: flood_protected_tracer_uc(:,:)
   integer, save :: flood_tracer_ledger_reports = 0

   ! reservoir state passed to the transport (empty when there are no reservoirs)
   real(r8), allocatable :: volresv_ct  (:)
   integer,  allocatable :: ucat2resv_ct(:)

   integer, save :: lim_diag_printed = 0

#ifdef CoLMDEBUG
   ! constituent balance of one routing call (cf. the water balance printed there)
   real(r8), allocatable :: trc_mass_bef(:), trc_mass_inp(:), trc_mass_dis(:)
#endif

   PUBLIC :: flood_tracer_credit_patch, flood_tracer_evap_patch, flood_tracer_land_patch

CONTAINS

   ! ---------
   SUBROUTINE ct_rlc_init (start_year, restart_file)

   ! River constituent pools: from the river restart file, else (or for the
   ! constituents missing in it) from the water state at the default ratio.
   ! Called after the river water state is complete (rebuild_volwater_ucat).

   USE MOD_Grid_RiverLakeNetwork,  only: lake_type
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, volresv
   USE MOD_Grid_Reservoir,         only: ucat2resv, dam_build_year
   IMPLICIT NONE

   integer,          intent(in) :: start_year
   character(len=*), intent(in) :: restart_file

   integer :: i, irsv
   logical :: trc_restart_found
   logical,  allocatable :: trc_missing(:), is_built_resv_init(:)
   real(r8), allocatable :: wdsrf_safe(:), volresv_safe(:)
   integer,  allocatable :: ucat2resv_safe(:)

      CALL ct_riverlake_init ()
      IF (ntracers <= 0) RETURN

      allocate (trc_missing(ntracers))
      trc_missing = .true.
      trc_restart_found = .false.
      IF (len_trim(restart_file) > 0) THEN
         CALL ct_riverlake_read_restart (restart_file, trc_restart_found, trc_missing)
      ENDIF

      IF (any(trc_missing)) THEN
         IF (allocated(wdsrf_ucat)) THEN
            allocate (wdsrf_safe(size(wdsrf_ucat)));    wdsrf_safe = wdsrf_ucat
         ELSE
            allocate (wdsrf_safe(0))
         ENDIF
         IF (allocated(volresv)) THEN
            allocate (volresv_safe(size(volresv)));     volresv_safe = volresv
         ELSE
            allocate (volresv_safe(0))
         ENDIF
         IF (allocated(ucat2resv)) THEN
            allocate (ucat2resv_safe(size(ucat2resv))); ucat2resv_safe = ucat2resv
         ELSE
            allocate (ucat2resv_safe(0))
         ENDIF
         ! reservoirs already built in the start year hold their own volume
         allocate (is_built_resv_init(size(wdsrf_safe)))
         is_built_resv_init = .false.
         IF (allocated(lake_type) .and. allocated(dam_build_year)) THEN
            DO i = 1, min(size(is_built_resv_init), size(lake_type), size(ucat2resv_safe))
               IF (lake_type(i) /= 2) CYCLE
               irsv = ucat2resv_safe(i)
               IF (irsv < 1 .or. irsv > size(dam_build_year)) CYCLE
               is_built_resv_init(i) = start_year >= dam_build_year(irsv)
            ENDDO
         ENDIF
         CALL ct_riverlake_init_from_water (wdsrf_safe, volresv_safe, ucat2resv_safe, &
            trc_missing, is_built_resv_init)
         deallocate (wdsrf_safe, volresv_safe, ucat2resv_safe, is_built_resv_init)
      ENDIF

      deallocate (trc_missing)

   END SUBROUTINE ct_rlc_init

   ! ---------
   SUBROUTINE ct_rlc_runoff_input (filter_rnof, rnof_uc, deltime)

   ! Constituent of the land runoff of this land step, mapped to the unit
   ! catchments as the runoff water is (patch -> input grid -> unit catchment),
   ! into the accumulation that is merged at the next routing call.
   ! Called on all workers (the pushes are collective).

   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, remap_patch2inpm, &
      push_inpm2ucat, push_ucat2inpm
   USE MOD_WorkerPushData
   USE MOD_CT_Vars,               only: trc_rnof_step
   IMPLICIT NONE

   logical,  intent(in) :: filter_rnof(:)   ! patches whose runoff enters the rivers
   real(r8), intent(in) :: rnof_uc(:)       ! runoff water [mm/s] per unit catchment
   real(r8), intent(in) :: deltime          ! land time step [s]

   integer :: itrc
   real(r8), allocatable :: trc_rnof_gd(:,:), trc_rnof_uc(:,:)

      IF (.not. p_is_worker) RETURN
      IF (ntracers <= 0) RETURN

      allocate (trc_rnof_gd (ntracers, numinpm))
      allocate (trc_rnof_uc (ntracers, numucat))
      trc_rnof_gd = 0._r8
      trc_rnof_uc = 0._r8

      DO itrc = 1, ntracers
         IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
         IF (numpatch > 0) &
            CALL worker_remap_data_pset2grid (remap_patch2inpm, trc_rnof_step(itrc,:), &
               trc_rnof_gd(itrc,:), fillvalue = 0._r8, filter = filter_rnof)
         IF (numinpm > 0) THEN
            WHERE (push_ucat2inpm%sum_area > 0._r8)
               trc_rnof_gd(itrc,:) = trc_rnof_gd(itrc,:) / push_ucat2inpm%sum_area
            END WHERE
         ENDIF
         CALL worker_push_data (push_inpm2ucat, trc_rnof_gd(itrc,:), trc_rnof_uc(itrc,:), &
            fillvalue = 0._r8, mode = 'sum')
      ENDDO

      ! water reference [m^3] as acc_rnof_uc; constituent [R*mm*m^2 -> R*m^3]
      IF (numucat > 0) &
         CALL ct_riverlake_input_from_runoff (rnof_uc*1.e-3*deltime, numucat, trc_rnof_uc*1.e-3)

      deallocate (trc_rnof_gd, trc_rnof_uc)

   END SUBROUTINE ct_rlc_runoff_input

   ! ---------
   SUBROUTINE ct_rlc_begin_routing ()

   ! Start of a routing call (workers): the runoff constituent accumulated over
   ! the land steps joins the input buffer, as acc_rnof_uc joins the water.

   USE MOD_Grid_RiverLakeNetwork,  only: numucat
   USE MOD_Grid_RiverLakeTimeVars, only: volresv
   USE MOD_Grid_Reservoir,         only: ucat2resv
   IMPLICIT NONE

   integer :: itrc

      IF (.not. p_is_worker) RETURN
      IF (ntracers <= 0) RETURN

      IF (allocated(volresv_ct)) deallocate (volresv_ct)
      IF (allocated(ucat2resv_ct)) deallocate (ucat2resv_ct)
      IF (allocated(volresv)) THEN
         allocate (volresv_ct(size(volresv)))
         volresv_ct = volresv
      ELSE
         allocate (volresv_ct(0))
      ENDIF
      IF (allocated(ucat2resv)) THEN
         allocate (ucat2resv_ct(size(ucat2resv)))
         ucat2resv_ct = ucat2resv
      ELSE
         allocate (ucat2resv_ct(0))
      ENDIF

      IF (allocated(trc_decay_source)) trc_decay_source = 0._r8

#ifdef CoLMDEBUG
      CALL ensure_debug_arrays ()
      IF (numucat > 0) THEN
         DO itrc = 1, ntracers
            IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
            trc_mass_bef(itrc) = pool_total(itrc)
            trc_mass_inp(itrc) = sum(acc_trc_inp(itrc,:))
         ENDDO
      ENDIF
#endif

      IF (numucat > 0) THEN
         DO itrc = 1, ntracers
            IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
            trc_inp_buf(itrc,:) = trc_inp_buf(itrc,:) + acc_trc_inp(itrc,:)
            acc_trc_inp(itrc,:) = 0._r8
         ENDDO
      ENDIF

   END SUBROUTINE ct_rlc_begin_routing

   ! ---------
   SUBROUTINE ct_rlc_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft, with_input)

   ! The host has moved water between the visible and the protected (levee)
   ! side of unit catchment i; move the constituent with it.  with_input: the
   ! water of this call's runoff is already in vis_bef, its constituent still
   ! in the input buffer (before the first substep).

   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: vis_bef, lev_bef, vis_aft, lev_aft   ! [m^3]
   logical,  intent(in), optional :: with_input

      IF (ntracers <= 0) RETURN
      IF (present(with_input)) THEN
         IF (with_input) THEN
            CALL ct_riverlake_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft, &
               pending_trc_pool = trc_inp_buf(:,i))
            RETURN
         ENDIF
      ENDIF
      CALL ct_riverlake_levee_repartition (i, vis_bef, lev_bef, vis_aft, lev_aft)

   END SUBROUTINE ct_rlc_levee_repartition

   ! ---------
   SUBROUTINE ct_rlc_equilibrate_cell (i, visible, with_levee)

   ! Dissolution equilibrium of unit catchment i for the given visible water and
   ! the protected water of its levee (if any), before the host takes water out.

   USE MOD_Grid_RiverLakeLevee, only: levsto
   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: visible      ! [m^3]
   logical,  intent(in) :: with_levee

   real(r8) :: protected

      IF (ntracers <= 0) RETURN
      protected = 0._r8
      IF (with_levee) protected = levsto(i)
      CALL ct_riverlake_equilibrate_cell (i, visible, protected)

   END SUBROUTINE ct_rlc_equilibrate_cell

   ! ---------
   SUBROUTINE ct_rlc_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, ucatfilter, is_built_resv)

   ! Constituent transport of one routing substep, with the water fluxes of the
   ! substep, before the host updates the storages with them.

   USE MOD_Grid_RiverLakeNetwork,     only: npthout_local
   USE MOD_Grid_RiverLakeTimeVars,    only: wdsrf_ucat, volresv
   USE MOD_Grid_RiverLakeBifurcation, only: bif_hflux_lev
   IMPLICIT NONE

   real(r8), intent(in) :: dt_ref
   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   real(r8), intent(in) :: hflux_fc(:)
   real(r8), intent(in) :: sum_hflux_riv(:)
   logical,  intent(in) :: ucatfilter(:)
   logical,  intent(in) :: is_built_resv(:)

      IF (.not. p_is_worker) RETURN
      IF (ntracers <= 0) RETURN

      IF (allocated(volresv)) volresv_ct = volresv
      IF (DEF_USE_BIFURCATION) THEN
         CALL ct_riverlake_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
            wdsrf_ucat, ucatfilter, volresv_ct, ucat2resv_ct, is_built_resv, &
            do_bif = .true., bif_hflux_lev_in = bif_hflux_lev, npthout_local_in = npthout_local)
      ELSE
         CALL ct_riverlake_substep (dt_ref, dt_all, irivsys, hflux_fc, sum_hflux_riv, &
            wdsrf_ucat, ucatfilter, volresv_ct, ucat2resv_ct, is_built_resv)
      ENDIF

   END SUBROUTINE ct_rlc_substep

   ! ---------
   SUBROUTINE ct_rlc_inland_depression (i, volwater, volmax, dt)

   ! The host removes the water above volmax from an inland depression; the
   ! constituent leaves with it at the cell's dissolved ratio.

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: volwater   ! water before the removal [m^3]
   real(r8), intent(in) :: volmax     ! water kept [m^3]
   real(r8), intent(in) :: dt         ! substep [s]

   integer  :: itrc
   real(r8) :: frac_remove, trc_removed

      IF (ntracers <= 0) RETURN
      IF (volwater <= 1.e-6_r8) RETURN

      frac_remove = (volwater - volmax) / volwater
      DO itrc = 1, ntracers
         IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
         IF (allocated(trc_solid)) THEN
            IF (ct_has_dissolved_limit(itrc)) &
               CALL ct_equilibrate_dissolved (itrc, volwater, trc_mass(itrc,i), trc_solid(itrc,i))
         ENDIF
         trc_removed = trc_mass(itrc,i) * frac_remove
         trc_mass(itrc,i) = trc_mass(itrc,i) - trc_removed
         trc_flux_out(itrc,i) = trc_removed / dt
      ENDDO

   END SUBROUTINE ct_rlc_inland_depression

   ! ---------
   SUBROUTINE ct_rlc_diag_accumulate (dt_all, irivsys, ucatfilter, is_built_resv)

   ! History accumulation of the substep, after the host storage update.

   USE MOD_Grid_RiverLakeNetwork,  only: numucat
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, volresv
   IMPLICIT NONE

   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   logical,  intent(in) :: ucatfilter(:)
   logical,  intent(in) :: is_built_resv(:)

      IF (.not. p_is_worker) RETURN
      IF (ntracers <= 0 .or. numucat <= 0) RETURN

      IF (allocated(volresv)) volresv_ct = volresv
      CALL ct_riverlake_diag_accumulate (dt_all, irivsys, ucatfilter, wdsrf_ucat, &
         volresv_ct, ucat2resv_ct, is_built_resv)

   END SUBROUTINE ct_rlc_diag_accumulate

   ! ---------
   SUBROUTINE ct_rlc_debug_discharge (i, dt)

   ! CoLMDEBUG: constituent leaving the network at outlet i in this substep.

   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: dt
#ifdef CoLMDEBUG
   integer :: itrc

      IF (ntracers <= 0) RETURN
      CALL ensure_debug_arrays ()
      DO itrc = 1, ntracers
         IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
         trc_mass_dis(itrc) = trc_mass_dis(itrc) + trc_flux_out(itrc,i) * dt
      ENDDO
#endif

   END SUBROUTINE ct_rlc_debug_discharge

   ! ---------
   SUBROUTINE ct_rlc_end_routing ()

   ! End of a routing call, on all processes: CoLMDEBUG constituent balance,
   ! donor-limiter statistics, and reset of the per-call accumulations.

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   integer :: lim_calls, lim_iter_sum, lim_iter_peak, lim_over_soft
#ifdef CoLMDEBUG
   integer :: itrc
   real(r8), allocatable :: trc_mass_aft(:), trc_mass_decay(:)
#endif

      IF (ntracers <= 0) RETURN

#ifdef CoLMDEBUG
      CALL ensure_debug_arrays ()
      allocate (trc_mass_aft(ntracers), trc_mass_decay(ntracers))
      trc_mass_aft = 0._r8
      trc_mass_decay = 0._r8
      IF (p_is_worker .and. numucat > 0) THEN
         DO itrc = 1, ntracers
            IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
            trc_mass_aft(itrc) = pool_total(itrc)
            trc_mass_dis(itrc) = trc_mass_dis(itrc) + sum(trc_dry_drain(itrc,:))
            IF (allocated(trc_decay_source)) trc_mass_decay(itrc) = sum(trc_decay_source(itrc,:))
         ENDDO
      ELSE
         trc_mass_bef = 0._r8
         trc_mass_inp = 0._r8
         trc_mass_dis = 0._r8
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_bef,   ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_aft,   ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_inp,   ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_dis,   ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, trc_mass_decay, ntracers, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         DO itrc = 1, ntracers
            IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
            write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc, ') mass before  : ', trc_mass_bef(itrc), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc, ') mass input   : ', trc_mass_inp(itrc), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc, ') mass discharge: ', trc_mass_dis(itrc), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc, ') mass decay   : ', trc_mass_decay(itrc), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc, ') mass after   : ', trc_mass_aft(itrc), ' R*m3'
            write(*,'(A,I0,A,ES12.4,A)') 'Tracer(', itrc, ') mass balance : ', &
               trc_mass_aft(itrc) - trc_mass_bef(itrc) - trc_mass_inp(itrc) + trc_mass_dis(itrc) &
               - trc_mass_decay(itrc), ' R*m3'
         ENDDO
      ENDIF
      trc_mass_bef = 0._r8
      trc_mass_inp = 0._r8
      trc_mass_dis = 0._r8
      deallocate (trc_mass_aft, trc_mass_decay)
#endif

      IF (p_is_worker) THEN
         CALL ct_riverlake_limiter_stats (lim_calls, lim_iter_sum, lim_iter_peak, lim_over_soft, &
            reset = .true.)
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, lim_calls,     1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, lim_iter_sum,  1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
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

         IF (numucat > 0) THEN
            IF (allocated(acc_trc_inp))   acc_trc_inp   = 0._r8
            IF (allocated(acc_rnof_ref))  acc_rnof_ref  = 0._r8
            IF (allocated(trc_dry_drain)) trc_dry_drain = 0._r8
         ENDIF
      ENDIF

   END SUBROUTINE ct_rlc_end_routing

   ! ---------
   SUBROUTINE ct_rlc_flood_init ()

   ! Flood-feedback constituent arrays on the current land patches: at start-up
   ! and after LULCC (called from init_flood_feedback, workers).

   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   integer :: itrc
   logical :: has_flood_tracer

      IF (.not. p_is_worker) RETURN

      IF (allocated(flood_tracer_credit_patch)) deallocate (flood_tracer_credit_patch, &
         flood_tracer_evap_patch, flood_tracer_land_patch, &
         flood_visible_tracer_uc, flood_protected_tracer_uc)

      has_flood_tracer = .false.
      DO itrc = 1, ntracers
         IF (ct_uses_land_water_transport(itrc)) has_flood_tracer = .true.
      ENDDO
      IF (.not. has_flood_tracer) RETURN

      allocate (flood_tracer_credit_patch(ntracers,numpatch), flood_tracer_evap_patch(ntracers,numpatch), &
         flood_tracer_land_patch(ntracers,numpatch), &
         flood_visible_tracer_uc(ntracers,numucat), flood_protected_tracer_uc(ntracers,numucat))
      flood_tracer_credit_patch = 0._r8
      flood_tracer_evap_patch   = 0._r8
      flood_tracer_land_patch   = 0._r8
      flood_visible_tracer_uc   = 0._r8
      flood_protected_tracer_uc = 0._r8

   END SUBROUTINE ct_rlc_flood_init

   ! ---------
   SUBROUTINE ct_rlc_flood_withhold ()

   ! No flood is published to the land: no constituent either.

      IF (.not. p_is_worker) RETURN
      IF (allocated(flood_tracer_credit_patch)) flood_tracer_credit_patch = 0._r8

   END SUBROUTINE ct_rlc_flood_withhold

   ! ---------
   SUBROUTINE ct_rlc_flood_publish_reset ()

   ! Start of publish_flood_feedback: a new flood state replaces the old one.

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(flood_visible_tracer_uc)) RETURN
      flood_visible_tracer_uc   = 0._r8
      flood_protected_tracer_uc = 0._r8
      flood_tracer_credit_patch = 0._r8
      flood_tracer_evap_patch   = 0._r8
      flood_tracer_land_patch   = 0._r8

   END SUBROUTINE ct_rlc_flood_publish_reset

   ! ---------
   SUBROUTINE ct_rlc_flood_publish_cell (i, visible, protected, visible_flood)

   ! Constituent of the flood water published from unit catchment i:
   ! the visible part above the river channel at the visible dissolved ratio,
   ! and all protected (levee) water.

   IMPLICIT NONE
   integer,  intent(in) :: i
   real(r8), intent(in) :: visible        ! visible water [m^3]
   real(r8), intent(in) :: protected      ! protected water [m^3]
   real(r8), intent(in) :: visible_flood  ! published visible flood water [m^3]

      IF (.not. allocated(flood_visible_tracer_uc)) RETURN
      IF (visible > 0._r8) flood_visible_tracer_uc(:,i) = &
         max(trc_mass(:,i), 0._r8) * (visible_flood / visible)
      IF (protected > 0._r8) flood_protected_tracer_uc(:,i) = max(trc_levsto(:,i), 0._r8)

   END SUBROUTINE ct_rlc_flood_publish_cell

   ! ---------
   SUBROUTINE ct_rlc_flood_publish_push (flood_grid_area)

   ! Map the published flood constituent to the patches, as the flood water is
   ! (unit catchment -> input grid -> patch), into flood_tracer_credit_patch.
   ! Called on all workers (the pushes are collective).

   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, remap_patch2inpm, &
      push_inpm2ucat, push_ucat2inpm
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: flood_grid_area(:)   ! land area of each input grid cell

   real(r8), parameter :: MISSING = -1.e30_r8
   integer :: itrc, iside
   real(r8), allocatable :: density(:), grid_tracer(:), patch_tracer(:)

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(flood_visible_tracer_uc)) RETURN

      allocate (density(numucat), grid_tracer(numinpm), patch_tracer(numpatch))
      DO itrc = 1, ntracers
         IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
         DO iside = 1, 2
            ! sum_area is allocated only on workers with requests (numucat > 0)
            density = 0._r8
            IF (numucat > 0) THEN
               IF (iside == 1) THEN
                  WHERE (push_inpm2ucat%sum_area > 0._r8) &
                     density = flood_visible_tracer_uc(itrc,:) / max(push_inpm2ucat%sum_area, tiny(1._r8))
               ELSE
                  WHERE (push_inpm2ucat%sum_area > 0._r8) &
                     density = flood_protected_tracer_uc(itrc,:) / max(push_inpm2ucat%sum_area, tiny(1._r8))
               ENDIF
            ENDIF
            CALL worker_push_data (push_ucat2inpm, density, grid_tracer, fillvalue=0._r8, mode='sum')
            WHERE (flood_grid_area > 0._r8)
               grid_tracer = grid_tracer / max(flood_grid_area, tiny(1._r8))
            ELSEWHERE
               grid_tracer = 0._r8
            END WHERE
            IF (numpatch > 0) THEN
               CALL worker_remap_data_grid2pset (remap_patch2inpm, grid_tracer, &
                  patch_tracer, fillvalue=MISSING, mode='average')
               WHERE (patch_tracer == MISSING) patch_tracer = 0._r8
               IF (iside == 1) THEN
                  flood_tracer_credit_patch(itrc,:) = patch_tracer
               ELSE
                  ! [m] -> [mm], as flood_credit_patch*1000 is the water credit [mm]
                  flood_tracer_credit_patch(itrc,:) = 1000._r8 * &
                     (flood_tracer_credit_patch(itrc,:) + patch_tracer)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      deallocate (density, grid_tracer, patch_tracer)

   END SUBROUTINE ct_rlc_flood_publish_push

   ! ---------
   SUBROUTINE ct_rlc_flood_debit (flood_credit_patch, flood_evap_acc, flood_infil_acc, &
         flood_grid_area, flood_visible_uc, flood_protected_uc)

   ! Constituent side of debit_flood_feedback (workers, collective): take from
   ! the river pools what the land took with flood evaporation (vapour
   ! exchange, flood_tracer_evap_patch) and infiltration (flood_tracer_land_patch),
   ! and stop if the land-side and river-side totals do not close.

   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_LandPatch,             only: numpatch
   USE MOD_Grid_RiverLakeNetwork, only: numinpm, numucat, remap_patch2inpm, push_inpm2ucat
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: flood_credit_patch(:)   ! published water credit [m]
   real(r8), intent(in) :: flood_evap_acc(:)       ! land flood evaporation [mm]
   real(r8), intent(in) :: flood_infil_acc(:)      ! land flood infiltration [mm]
   real(r8), intent(in) :: flood_grid_area(:)
   real(r8), intent(in) :: flood_visible_uc(:)     ! published visible flood water [m^3]
   real(r8), intent(in) :: flood_protected_uc(:)   ! published protected flood water [m^3]

   integer  :: i, j, itrc
   real(r8) :: patch_area, tracer_cell_debit, ledger_tol
   real(r8) :: water_credit, tracer_credit, vapor_loss, positive_loss, vapor_gain
   real(r8) :: infiltrated_fraction, evaporated_fraction, tracer_after
   real(r8), allocatable :: ratio_patch(:), ratio_grid(:), ratio_uc(:)
   real(r8), allocatable :: gain_patch(:), coefficient_uc(:), gain_uc(:)
   real(r8), allocatable :: tracer_exchange_ledger(:,:)   ! (land, vapour, river donor)

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(flood_visible_tracer_uc)) RETURN

      allocate (ratio_patch(numpatch), ratio_grid(numinpm), ratio_uc(numucat))
      allocate (gain_patch(numpatch), coefficient_uc(numucat), gain_uc(numucat))
      allocate (tracer_exchange_ledger(3,ntracers))
      tracer_exchange_ledger = 0._r8

      DO itrc = 1, ntracers
         IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
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
            IF (tracer_credit > 0._r8) evaporated_fraction = min(1._r8,positive_loss/tracer_credit)
            infiltrated_fraction = 0._r8
            IF (flood_infil_acc(i) > 0._r8) THEN
               IF (water_credit <= flood_evap_acc(i)) THEN
                  IF (flood_infil_acc(i) > 1.e-12_r8*max(water_credit, 1._r8)) &
                     CALL CoLM_stop('grid flood feedback: infiltration after full evaporation')
               ELSE
                  infiltrated_fraction = flood_infil_acc(i)/(water_credit-flood_evap_acc(i))
               ENDIF
            ENDIF
            IF (infiltrated_fraction > 1._r8+1.e-10_r8) &
               CALL CoLM_stop('grid flood feedback: tracer infiltration overdraft')
            infiltrated_fraction = min(1._r8,max(0._r8,infiltrated_fraction))
            ratio_patch(i) = evaporated_fraction + (1._r8-evaporated_fraction)*infiltrated_fraction
            IF (ct_has_dissolved_limit(itrc) .and. tracer_credit > 0._r8) &
               ratio_patch(i) = min(1._r8,max(0._r8,flood_tracer_land_patch(itrc,i)/tracer_credit))
            gain_patch(i) = vapor_gain*(1._r8-infiltrated_fraction)/water_credit
            patch_area = sum(remap_patch2inpm%areapart(i)%val)*1.e-3_r8
            IF (ct_has_dissolved_limit(itrc)) THEN
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
            tracer_exchange_ledger(2,itrc) = tracer_exchange_ledger(2,itrc) + vapor_loss*patch_area
         ENDDO

         ! fraction of the published constituent taken, per unit catchment
         ratio_grid = 0._r8
         IF (numpatch > 0) CALL worker_remap_data_pset2grid (remap_patch2inpm, ratio_patch, ratio_grid, &
            fillvalue=0._r8, filter=flood_credit_patch>0._r8)
         WHERE (flood_grid_area > 0._r8)
            ratio_grid = ratio_grid/max(flood_grid_area,tiny(1._r8))
         ELSEWHERE
            ratio_grid = 0._r8
         END WHERE
         CALL worker_push_data (push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
         coefficient_uc = 0._r8
         IF (numucat > 0) THEN
            WHERE (push_inpm2ucat%sum_area > 0._r8) &
               coefficient_uc = min(1._r8,max(0._r8,ratio_uc/push_inpm2ucat%sum_area))
         ENDIF

         ! vapour gained by the flood water still in the river (isotopes only)
         gain_uc = 0._r8
         IF (.not. tracer_is_nonvolatile_solute(itrc)) THEN
            ratio_grid = 0._r8
            IF (numpatch > 0) CALL worker_remap_data_pset2grid (remap_patch2inpm, gain_patch, ratio_grid, &
               fillvalue=0._r8, filter=flood_credit_patch>0._r8)
            WHERE (flood_grid_area > 0._r8)
               ratio_grid = ratio_grid/max(flood_grid_area,tiny(1._r8))
            ELSEWHERE
               ratio_grid = 0._r8
            END WHERE
            CALL worker_push_data (push_inpm2ucat, ratio_grid, ratio_uc, fillvalue=0._r8, mode='sum')
            IF (numucat > 0) THEN
               WHERE (push_inpm2ucat%sum_area > 0._r8) gain_uc = ratio_uc/push_inpm2ucat%sum_area
            ENDIF
         ENDIF

         DO j = 1, numucat
            tracer_after = trc_mass(itrc,j) - coefficient_uc(j)*flood_visible_tracer_uc(itrc,j) &
               + gain_uc(j)*flood_visible_uc(j)
            IF (tracer_after < -1.e-10_r8*max(1._r8,trc_mass(itrc,j))) &
               CALL CoLM_stop('grid flood feedback: negative visible tracer')
            tracer_cell_debit = coefficient_uc(j)*flood_visible_tracer_uc(itrc,j) &
               - gain_uc(j)*flood_visible_uc(j) + min(tracer_after,0._r8)
            trc_mass(itrc,j) = max(0._r8,tracer_after)
            tracer_after = trc_levsto(itrc,j) - coefficient_uc(j)*flood_protected_tracer_uc(itrc,j) &
               + gain_uc(j)*flood_protected_uc(j)
            IF (tracer_after < -1.e-10_r8*max(1._r8,trc_levsto(itrc,j))) &
               CALL CoLM_stop('grid flood feedback: negative protected tracer')
            tracer_cell_debit = tracer_cell_debit + coefficient_uc(j)*flood_protected_tracer_uc(itrc,j) &
               - gain_uc(j)*flood_protected_uc(j) + min(tracer_after,0._r8)
            trc_levsto(itrc,j) = max(0._r8,tracer_after)
            tracer_exchange_ledger(3,itrc) = tracer_exchange_ledger(3,itrc) + tracer_cell_debit
         ENDDO
      ENDDO

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, tracer_exchange_ledger, 3*ntracers, MPI_REAL8, &
         MPI_SUM, p_comm_worker, p_err)
#endif
      DO itrc = 1, ntracers
         IF (.not. ct_uses_land_water_transport(itrc)) CYCLE
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

      deallocate (ratio_patch, ratio_grid, ratio_uc, gain_patch, coefficient_uc, gain_uc, &
         tracer_exchange_ledger)

   END SUBROUTINE ct_rlc_flood_debit

   ! ---------
   SUBROUTINE ct_rlc_final ()

      CALL ct_riverlake_final ()

      IF (allocated(flood_tracer_credit_patch)) deallocate (flood_tracer_credit_patch, &
         flood_tracer_evap_patch, flood_tracer_land_patch, &
         flood_visible_tracer_uc, flood_protected_tracer_uc)
      IF (allocated(volresv_ct))   deallocate (volresv_ct)
      IF (allocated(ucat2resv_ct)) deallocate (ucat2resv_ct)
#ifdef CoLMDEBUG
      IF (allocated(trc_mass_bef)) deallocate (trc_mass_bef, trc_mass_inp, trc_mass_dis)
#endif

   END SUBROUTINE ct_rlc_final

#ifdef CoLMDEBUG
   ! ---------
   SUBROUTINE ensure_debug_arrays ()

      IF (allocated(trc_mass_bef)) RETURN
      allocate (trc_mass_bef(ntracers), trc_mass_inp(ntracers), trc_mass_dis(ntracers))
      trc_mass_bef = 0._r8
      trc_mass_inp = 0._r8
      trc_mass_dis = 0._r8

   END SUBROUTINE ensure_debug_arrays

   ! ---------
   real(r8) FUNCTION pool_total (itrc)

   ! all river pools of constituent itrc on this worker [R*m^3]

   integer, intent(in) :: itrc

      pool_total = sum(trc_mass(itrc,:)) + sum(trc_inp_buf(itrc,:))
      IF (allocated(trc_levsto)) pool_total = pool_total + sum(trc_levsto(itrc,:))
      IF (allocated(trc_solid))  pool_total = pool_total + sum(trc_solid(itrc,:)) &
         + sum(trc_levsto_solid(itrc,:))

   END FUNCTION pool_total
#endif

END MODULE MOD_CT_RiverLakeCoupling
#endif
