#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeBifurcation
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Bifurcation (multi-channel flow) module for grid-based river-lake routing.
!   Computes water exchange through bifurcation pathways using a Riemann-solver-
!   based shallow water formulation consistent with the main channel solver in
!   MOD_Grid_RiverLakeFlow.
!
! Created by CoLM team, April 2026
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_LEVEE
   USE MOD_SPMD_Task
   USE MOD_Grid_Reservoir, only: ucat2resv
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levdph, levsto, &
      volwater_ucat, volwater_ucat_valid, levee_visible_volume_from_stage
   USE MOD_Grid_RiverLakeNetwork, only: floodplain_curve
   IMPLICIT NONE

   real(r8), parameter :: BIFMIN = 1.e-5_r8
   real(r8), parameter :: BIF_STORAGE_EPS = 1.e-6_r8

   ! ----- State variables -----
   real(r8), allocatable :: pth_veloc     (:,:)  ! velocity (npthlev, npthout_local) [m/s]
   real(r8), allocatable :: pth_momen     (:,:)  ! momentum (npthlev, npthout_local) [m^2/s]
   real(r8), allocatable :: bif_hflux_lev (:,:)  ! effective volume flux per pathway layer [m^3/s]
   real(r8), allocatable :: bif_hflux_sum (:)    ! net volume flux per ucat [m^3/s]
   logical, allocatable :: bif_path_active (:)    ! true when pathway was eligible this call
   ! Per-ucat layer-2+ ("overland" / above-bankfull) bif
   ! water flux. tracer_substep already subtracts the matching layer-2+
   ! tracer flux from trc_levsto, but the water-side accumulator
   ! `sum_hflux_riv` (MOD_Grid_RiverLakeFlow.F90:820) lumps every layer
   ! onto the visible side. The second levee_tracer_repartition uses
   ! `levsto_bef = levsto(i)` as the ratio denominator while
   ! `trc_levsto` already reflects the layer-2+ subtraction — denominator
   ! and numerator are out of phase, leaking the Phase-1 R_init invariant
   ! at the levee/bif interface. Exposing this per-cell layer-2+ flux
   ! lets the caller compute the post-bif protected water and align both.
   real(r8), allocatable :: bif_lev_hflux_sum(:) ! net layer-2+ vol flux per ucat [m^3/s]
   logical, save :: dbg_bif_restart_checked = .false.
   integer, save :: n_dt_mismatch_warn_count = 0

   ! ----- Persistent scratch buffers (allocated once in bifurcation_init,
   !       reused every bifurcation_calc call to avoid alloc/dealloc churn).
   !       TARGET attribute lets bifurcation_calc alias them via pointers.  -----
   real(r8), allocatable, target, save :: wdsrf_dn_pth_buf     (:)
   real(r8), allocatable, target, save :: rivelv_dn_pth_buf    (:)
   real(r8), allocatable, target, save :: protected_wdsrf_ucat_buf  (:)
   real(r8), allocatable, target, save :: protected_wdsrf_dn_pth_buf (:)
   real(r8), allocatable, target, save :: has_levee_dn_pth_buf (:)
   real(r8), allocatable, target, save :: storage_dn_pth_buf   (:)
   real(r8), allocatable, target, save :: pth_hflux_total_buf  (:)
   real(r8), allocatable, target, save :: pth_limiter_rate_buf (:)
   real(r8), allocatable, target, save :: bif_influx_buf       (:)
   real(r8), allocatable, target, save :: has_levee_r8_buf     (:)
   real(r8), allocatable, target, save :: storage_ucat_buf     (:)
   real(r8), allocatable, target, save :: limiter_outgoing_buf (:)
   real(r8), allocatable, target, save :: limiter_incoming_buf (:)
   real(r8), allocatable, target, save :: limiter_out_rate_buf (:)
   real(r8), allocatable, target, save :: limiter_in_rate_buf  (:)

   ! ----- cross-rank limiter buffers -----
   real(r8), allocatable, target, save :: pos_pth_buf              (:)
   real(r8), allocatable, target, save :: neg_pth_buf              (:)
   real(r8), allocatable, target, save :: bif_pos_recv_buf         (:)
   real(r8), allocatable, target, save :: bif_neg_recv_buf         (:)
   real(r8), allocatable, target, save :: limiter_in_rate_pth_buf  (:)
   real(r8), allocatable, target, save :: limiter_out_rate_pth_buf (:)

   ! ----- per-pathway layer-2+ flux + remote influx buffers -----
   real(r8), allocatable, target, save :: pth_lev_hflux_total_buf  (:)
   real(r8), allocatable, target, save :: bif_lev_influx_buf       (:)

   ! ----- cross-river-system sync buffers -----
   ! ucatfilter materialised as real8 so we can push it to pathway space
   ! and let source-rank pathways see whether their downstream cell is
   ! active in the current adaptive sub-step. Real-valued mask follows the
   ! same pattern as has_levee_r8.
   real(r8), allocatable, target, save :: ucatfilter_r8_buf        (:)
   real(r8), allocatable, target, save :: ucatfilter_dn_pth_buf    (:)
   real(r8), allocatable, target, save :: is_resv_r8_buf           (:)
   real(r8), allocatable, target, save :: is_resv_dn_pth_buf       (:)
   real(r8), allocatable, target, save :: dt_ucat_buf              (:)
   real(r8), allocatable, target, save :: dt_dn_pth_buf            (:)

   PUBLIC :: bifurcation_init
   PUBLIC :: bifurcation_calc
   PUBLIC :: read_bifurcation_restart
   PUBLIC :: write_bifurcation_restart
   PUBLIC :: bifurcation_final
   PUBLIC :: bif_lev_hflux_sum
   PUBLIC :: bif_path_active

CONTAINS

   ! =========================================================================
   SUBROUTINE bifurcation_init ()
   ! =========================================================================
   !
   ! Allocate and initialize bifurcation state arrays.
   !
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: numucat, npthout_local, npthlev_bif

   IMPLICIT NONE

      ! bifurcation_init is contractually one-shot per run, but the
      ! allocates below will blow up with a runtime error if the module is
      ! ever re-initialised (e.g. during a future restart-in-place refactor).
      ! A single allocated() check on the first state array is enough to
      ! make re-entry a no-op instead of a hard stop, since all buffers
      ! below are allocated as a group.
      IF (allocated(pth_veloc)) RETURN

      IF (.not. p_is_worker) THEN
         CALL allocate_bifurcation_arrays (0, 0, 0)
         RETURN
      ENDIF

      CALL allocate_bifurcation_arrays (numucat, npthlev_bif, npthout_local)

   END SUBROUTINE bifurcation_init


   SUBROUTINE allocate_bifurcation_arrays (nucat, nlev, npth)
   ! =========================================================================
   !
   ! Allocate all persistent bifurcation buffers. Zero-length dimensions are
   ! intentional on non-worker ranks so module state is always allocated safely.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer, intent(in) :: nucat
   integer, intent(in) :: nlev
   integer, intent(in) :: npth

      allocate (pth_veloc     (nlev, npth))
      allocate (pth_momen     (nlev, npth))
      allocate (bif_hflux_lev (nlev, npth))
      allocate (bif_hflux_sum (nucat))
      allocate (bif_path_active (npth))
      allocate (bif_lev_hflux_sum (nucat))

      pth_veloc     (:,:) = 0._r8
      pth_momen     (:,:) = 0._r8
      bif_hflux_lev (:,:) = 0._r8
      bif_hflux_sum (:)   = 0._r8
      bif_path_active(:)  = .false.
      bif_lev_hflux_sum(:) = 0._r8

      ! Persistent scratch buffers used by bifurcation_calc. Sized once here
      ! using the network-fixed dimensions; zero-length is valid.
      allocate (wdsrf_dn_pth_buf     (npth))
      allocate (rivelv_dn_pth_buf    (npth))
      allocate (protected_wdsrf_ucat_buf   (nucat))
      allocate (protected_wdsrf_dn_pth_buf (npth))
      allocate (has_levee_dn_pth_buf (npth))
      allocate (storage_dn_pth_buf   (npth))
      allocate (pth_hflux_total_buf  (npth))
      allocate (pth_limiter_rate_buf (npth))
      allocate (bif_influx_buf       (nucat))
      allocate (has_levee_r8_buf     (nucat))
      allocate (storage_ucat_buf     (nucat))
      allocate (limiter_outgoing_buf (nucat))
      allocate (limiter_incoming_buf (nucat))
      allocate (limiter_out_rate_buf (nucat))
      allocate (limiter_in_rate_buf  (nucat))

      allocate (pos_pth_buf              (npth))
      allocate (neg_pth_buf              (npth))
      allocate (bif_pos_recv_buf         (nucat))
      allocate (bif_neg_recv_buf         (nucat))
      allocate (limiter_in_rate_pth_buf  (npth))
      allocate (limiter_out_rate_pth_buf (npth))

      allocate (ucatfilter_r8_buf        (nucat))
      allocate (ucatfilter_dn_pth_buf    (npth))
      allocate (is_resv_r8_buf           (nucat))
      allocate (is_resv_dn_pth_buf       (npth))
      allocate (dt_ucat_buf              (nucat))
      allocate (dt_dn_pth_buf            (npth))

      allocate (pth_lev_hflux_total_buf  (npth))
      allocate (bif_lev_influx_buf       (nucat))

   END SUBROUTINE allocate_bifurcation_arrays


   ! =========================================================================
   SUBROUTINE bifurcation_calc (wdsrf_ucat, volresv, is_built_resv, dt_all, irivsys, ucatfilter, update_state)
   ! =========================================================================
   !
   ! Compute bifurcation fluxes for one sub-timestep using a two-rarefaction
   ! Riemann solver. The net volume flux per unit catchment is accumulated
   ! into bif_hflux_sum(:).
   !
   ! =========================================================================

   USE MOD_Const_Physical, only: grav
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakeNetwork, only: &
      numucat, npthout_local, npthlev_bif, &
      pth_upst_local, pth_down_local, &
      pth_dst, pth_elv, pth_wth, pth_man, &
      push_bif_dn2pth, push_bif_influx, &
      topo_rivelv

   IMPLICIT NONE

   real(r8), intent(in) :: wdsrf_ucat (:)   ! water depth above riverbed [m]
   real(r8), intent(in) :: volresv    (:)   ! reservoir water volume [m^3]
   logical,  intent(in) :: is_built_resv (:) ! reservoir-built mask
   real(r8), intent(in) :: dt_all     (:)   ! timestep per ucat [s]
   integer,  intent(in) :: irivsys    (:)   ! river system id per ucat
   logical,  intent(in) :: ucatfilter (:)   ! active ucat filter
   logical,  intent(in), optional :: update_state ! true: advance pathway momentum state

   ! Pointer aliases to persistent module-level buffers (no per-call alloc).
   real(r8), pointer :: wdsrf_dn_pth     (:) => null()
   real(r8), pointer :: rivelv_dn_pth    (:) => null()
   real(r8), pointer :: protected_wdsrf_ucat  (:) => null()
   real(r8), pointer :: protected_wdsrf_dn_pth (:) => null()
   real(r8), pointer :: has_levee_dn_pth (:) => null()
   real(r8), pointer :: storage_dn_pth   (:) => null()
   real(r8), pointer :: storage_ucat     (:) => null()
   real(r8), pointer :: pth_hflux_total  (:) => null()
   real(r8), pointer :: bif_influx       (:) => null()
   real(r8), pointer :: has_levee_r8     (:) => null()
   real(r8), pointer :: limiter_outgoing (:) => null()
   real(r8), pointer :: limiter_incoming (:) => null()
   real(r8), pointer :: limiter_out_rate (:) => null()
   real(r8), pointer :: limiter_in_rate  (:) => null()
   real(r8), pointer :: pth_limiter_rate (:) => null()
   real(r8), pointer :: pos_pth              (:) => null()
   real(r8), pointer :: neg_pth              (:) => null()
   real(r8), pointer :: bif_pos_recv         (:) => null()
   real(r8), pointer :: bif_neg_recv         (:) => null()
   real(r8), pointer :: limiter_in_rate_pth  (:) => null()
   real(r8), pointer :: limiter_out_rate_pth (:) => null()
   real(r8), pointer :: ucatfilter_r8        (:) => null()
   real(r8), pointer :: ucatfilter_dn_pth    (:) => null()
   real(r8), pointer :: is_resv_r8           (:) => null()
   real(r8), pointer :: is_resv_dn_pth       (:) => null()
   real(r8), pointer :: dt_ucat              (:) => null()
   real(r8), pointer :: dt_dn_pth            (:) => null()
   real(r8), pointer :: pth_lev_hflux_total  (:) => null()
   real(r8), pointer :: bif_lev_influx       (:) => null()

   real(r8) :: rivelv_up, rivelv_dn, wdsrf_up, wdsrf_dn, wdsrf_up_eff, wdsrf_dn_eff
   real(r8) :: bedelv_pth, height_up, height_dn
   real(r8) :: v_up, v_dn
   real(r8) :: veloct_fc, height_fc
   real(r8) :: vwave_up, vwave_dn
   real(r8) :: hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: hflux_lev, mflux_lev  ! flux for this layer
   real(r8) :: width_pth, pth_are
   real(r8) :: friction
   real(r8) :: dt, dt_dn, dt_cell, storage_up, storage_dn, storage_ref, rate
   real(r8) :: path_transfer
   integer  :: ipth, ilev, i_up, i_dn, i_ucat
   integer  :: n_dt_mismatch_skip, n_dt_mismatch_skip_glb
   logical  :: do_update_state

      IF (.not. p_is_worker) RETURN

      ! Contract: bifurcation_init must have been called before bifurcation_calc.
      ! Guarding so a missing init shows up as an explicit message instead of a
      ! silent segfault on the zero-assign below.
      IF (.not. allocated(bif_hflux_sum) .or. .not. allocated(bif_hflux_lev)) THEN
         write(*,'(A,I0)') 'ERROR: bifurcation_calc called before bifurcation_init on glb=', p_iam_glb
         call flush(6)
         CALL CoLM_stop ()
      ENDIF

      do_update_state = .true.
      IF (present(update_state)) do_update_state = update_state
      n_dt_mismatch_skip = 0
      n_dt_mismatch_skip_glb = 0

      ! Reset net flux accumulator
      bif_hflux_sum(:) = 0._r8
      bif_hflux_lev(:,:) = 0._r8
      IF (allocated(bif_path_active)) bif_path_active(:) = .false.
      IF (allocated(bif_lev_hflux_sum)) bif_lev_hflux_sum(:) = 0._r8

      ! Alias persistent module-level buffers (allocated once in
      ! bifurcation_init). Zero-length buffers are still valid targets.
      wdsrf_dn_pth     => wdsrf_dn_pth_buf
      rivelv_dn_pth    => rivelv_dn_pth_buf
      protected_wdsrf_ucat   => protected_wdsrf_ucat_buf
      protected_wdsrf_dn_pth => protected_wdsrf_dn_pth_buf
      has_levee_dn_pth => has_levee_dn_pth_buf
      storage_dn_pth   => storage_dn_pth_buf
      pth_hflux_total  => pth_hflux_total_buf
      pth_limiter_rate => pth_limiter_rate_buf
      bif_influx       => bif_influx_buf
      has_levee_r8     => has_levee_r8_buf
      storage_ucat     => storage_ucat_buf
      limiter_outgoing => limiter_outgoing_buf
      limiter_incoming => limiter_incoming_buf
      limiter_out_rate => limiter_out_rate_buf
      limiter_in_rate  => limiter_in_rate_buf
      pos_pth              => pos_pth_buf
      neg_pth              => neg_pth_buf
      bif_pos_recv         => bif_pos_recv_buf
      bif_neg_recv         => bif_neg_recv_buf
      limiter_in_rate_pth  => limiter_in_rate_pth_buf
      limiter_out_rate_pth => limiter_out_rate_pth_buf
      ucatfilter_r8        => ucatfilter_r8_buf
      ucatfilter_dn_pth    => ucatfilter_dn_pth_buf
      is_resv_r8           => is_resv_r8_buf
      is_resv_dn_pth       => is_resv_dn_pth_buf
      dt_ucat              => dt_ucat_buf
      dt_dn_pth            => dt_dn_pth_buf
      pth_lev_hflux_total  => pth_lev_hflux_total_buf
      bif_lev_influx       => bif_lev_influx_buf

      wdsrf_dn_pth     (:) = 0._r8
      rivelv_dn_pth    (:) = 0._r8
      protected_wdsrf_ucat  (:) = 0._r8
      protected_wdsrf_dn_pth(:) = 0._r8
      has_levee_dn_pth (:) = 0._r8
      storage_dn_pth   (:) = 0._r8
      pth_hflux_total(:) = 0._r8
      ! storage_ucat is fully overwritten by the available_storage_ucat
      ! loop below, so there is no need to zero it here.
      limiter_outgoing(:) = 0._r8
      limiter_incoming(:) = 0._r8
      limiter_out_rate(:) = 1._r8
      limiter_in_rate (:) = 1._r8
      pth_limiter_rate(:) = 1._r8
      pth_lev_hflux_total(:) = 0._r8
      bif_lev_influx     (:) = 0._r8
      ! has_levee_r8 is referenced only under DEF_USE_LEVEE guards, but
      ! zero it unconditionally so a future refactor that decouples the
      ! guards can't read uninitialised values.
      has_levee_r8(:) = 0._r8
      pos_pth             (:) = 0._r8
      neg_pth             (:) = 0._r8
      bif_pos_recv        (:) = 0._r8
      bif_neg_recv        (:) = 0._r8
      limiter_in_rate_pth (:) = 1._r8
      limiter_out_rate_pth(:) = 1._r8
      ! materialise ucatfilter as r8 so the push machinery can ship
      ! it to each pathway's destination cell. The dn_pth view is filled by
      ! the Step 2 push below.
      ucatfilter_r8    (:) = 0._r8
      ucatfilter_dn_pth(:) = 0._r8
      is_resv_r8       (:) = 0._r8
      is_resv_dn_pth   (:) = 0._r8
      dt_ucat          (:) = 0._r8
      dt_dn_pth        (:) = 0._r8
      IF (numucat > 0) THEN
         WHERE (ucatfilter) ucatfilter_r8 = 1._r8
         WHERE (is_built_resv) is_resv_r8 = 1._r8
         DO i_ucat = 1, numucat
            IF (irivsys(i_ucat) > 0 .and. irivsys(i_ucat) <= size(dt_all)) THEN
               dt_ucat(i_ucat) = dt_all(irivsys(i_ucat))
            ENDIF
         ENDDO
      ENDIF

      IF (DEF_USE_LEVEE) THEN
         IF (allocated(has_levee)) THEN
            WHERE (has_levee)
               has_levee_r8 = 1._r8
            END WHERE
            DO i_ucat = 1, numucat
               IF (has_levee(i_ucat)) THEN
                  protected_wdsrf_ucat(i_ucat) = floodplain_curve(i_ucat)%rivhgt + levdph(i_ucat)
               ELSE
                  protected_wdsrf_ucat(i_ucat) = wdsrf_ucat(i_ucat)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF (numucat > 0) THEN
         DO i_ucat = 1, numucat
            storage_ucat(i_ucat) = available_storage_ucat(i_ucat, wdsrf_ucat(i_ucat), volresv, is_built_resv)
         ENDDO
      ENDIF

      ! ----- Step 2: Get downstream cell state via push objects -----
      ! Use -9999 as fillvalue to mark pathways with no valid downstream cell
      ! (domain boundary or unresolved remote). These are skipped in Step 3.
      CALL worker_push_data (push_bif_dn2pth, wdsrf_ucat,  wdsrf_dn_pth,  fillvalue = -9999._r8)
      CALL worker_push_data (push_bif_dn2pth, topo_rivelv, rivelv_dn_pth, fillvalue = -9999._r8)
      IF (DEF_USE_LEVEE) THEN
         CALL worker_push_data (push_bif_dn2pth, protected_wdsrf_ucat, protected_wdsrf_dn_pth, fillvalue = 0._r8)
         CALL worker_push_data (push_bif_dn2pth, has_levee_r8,         has_levee_dn_pth,       fillvalue = 0._r8)
      ENDIF
      CALL worker_push_data (push_bif_dn2pth, storage_ucat, storage_dn_pth, fillvalue = -9999._r8)
      CALL worker_push_data (push_bif_dn2pth, is_resv_r8, is_resv_dn_pth, fillvalue = 0._r8)
      CALL worker_push_data (push_bif_dn2pth, dt_ucat, dt_dn_pth, fillvalue = 0._r8)
      ! push destination ucat active-mask so the source rank can
      ! skip pathways whose downstream cell lives in a river system that
      ! has already exhausted its adaptive sub-step. Fillvalue = 0 means
      ! unresolved remote/boundary destinations are treated as inactive
      ! (those are already handled by the wdsrf/rivelv boundary check
      ! above; this fill is only a defensive default).
      CALL worker_push_data (push_bif_dn2pth, ucatfilter_r8, ucatfilter_dn_pth, fillvalue = 0._r8)

      ! ----- Step 3: Riemann solver for each pathway / layer -----
      DO ipth = 1, npthout_local

         i_up = pth_upst_local(ipth)

         ! Skip if upstream ucat is not active
         IF (i_up < 1 .or. i_up > numucat) CYCLE
         IF (.not. ucatfilter(i_up)) CYCLE
         IF (i_up <= size(is_built_resv) .and. is_built_resv(i_up)) CYCLE

         ! Skip if downstream cell is outside domain (CaMa: CYCLE when JSEQP<=0)
         IF (wdsrf_dn_pth(ipth) < -9000._r8 .or. rivelv_dn_pth(ipth) < -9000._r8) CYCLE

         rivelv_up = topo_rivelv(i_up)
         wdsrf_up  = wdsrf_ucat(i_up)
         rivelv_dn = rivelv_dn_pth(ipth)
         wdsrf_dn  = wdsrf_dn_pth(ipth)
         dt        = dt_ucat(i_up)

         pth_hflux_total(ipth) = 0._r8

         i_dn = pth_down_local(ipth)

         ! skip pathway when destination cell is in a river system
         ! that has exhausted its adaptive sub-step (dt = 0). Applying the
         ! flux would update bif_hflux_sum(i_up) on the source side, but
         ! the destination contribution is dropped by WHERE(ucatfilter) in
         ! MOD_Grid_RiverLakeFlow, producing a cross-system conservation
         ! violation. Deferring the transfer to the next top-level routing
         ! step is safe: pathway momentum state is carried forward
         ! unchanged, and the two systems resynchronise at the next outer
         ! routing step.
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            IF (.not. ucatfilter(i_dn)) CYCLE
            IF (i_dn <= size(is_built_resv) .and. is_built_resv(i_dn)) CYCLE
            dt_dn = dt_ucat(i_dn)
         ELSE
            IF (ucatfilter_dn_pth(ipth) < 0.5_r8) CYCLE
            IF (is_resv_dn_pth(ipth) > 0.5_r8) CYCLE
            dt_dn = dt_dn_pth(ipth)
         ENDIF

         IF (dt <= 0._r8 .or. dt_dn <= 0._r8) CYCLE
         IF (abs(dt - dt_dn) > max(1.e-9_r8, 1.e-12_r8 * max(abs(dt), abs(dt_dn)))) THEN
            n_dt_mismatch_skip = n_dt_mismatch_skip + 1
            CYCLE
         ENDIF

         ! Zero-length pathways have no physical length scale for the
         ! momentum update; clear all layer states and skip this pathway.
         IF (pth_dst(ipth) <= 0._r8) THEN
            IF (do_update_state) THEN
               pth_momen(:, ipth) = 0._r8
               pth_veloc(:, ipth) = 0._r8
            ENDIF
            CYCLE
         ENDIF

         IF (allocated(bif_path_active)) bif_path_active(ipth) = .true.

         DO ilev = 1, npthlev_bif

            width_pth = pth_wth(ilev, ipth)
            IF (width_pth <= 0._r8) THEN
               IF (do_update_state) THEN
                  pth_momen(ilev, ipth) = 0._r8
                  pth_veloc(ilev, ipth) = 0._r8
               ENDIF
               CYCLE
            ENDIF

            IF (ilev == 1) THEN
               wdsrf_up_eff = wdsrf_up
            ELSE
               IF (DEF_USE_LEVEE .and. allocated(has_levee) .and. has_levee(i_up)) THEN
                  wdsrf_up_eff = protected_wdsrf_ucat(i_up)
               ELSE
                  wdsrf_up_eff = wdsrf_up
               ENDIF
            ENDIF

            IF (ilev == 1) THEN
               wdsrf_dn_eff = wdsrf_dn
            ELSE
               IF (i_dn > 0 .and. i_dn <= numucat) THEN
                  IF (DEF_USE_LEVEE .and. allocated(has_levee) .and. has_levee(i_dn)) THEN
                     wdsrf_dn_eff = protected_wdsrf_ucat(i_dn)
                  ELSE
                     wdsrf_dn_eff = wdsrf_dn
                  ENDIF
               ELSEIF (DEF_USE_LEVEE .and. has_levee_dn_pth(ipth) > 0.5_r8) THEN
                  wdsrf_dn_eff = protected_wdsrf_dn_pth(ipth)
               ELSE
                  wdsrf_dn_eff = wdsrf_dn
               ENDIF
            ENDIF

            ! Effective bed elevation for this pathway layer
            ! Use max of upstream bed, downstream bed, and pathway layer elevation
            ! (same logic as main channel: bedelv_fc = max(rivelv_i, bedelv_next_i))
            bedelv_pth = max(rivelv_up, rivelv_dn, pth_elv(ilev, ipth))

            ! Flow depths at both ends relative to pathway bed
            height_up = max(0._r8, wdsrf_up_eff + rivelv_up - bedelv_pth)
            height_dn = max(0._r8, wdsrf_dn_eff + rivelv_dn - bedelv_pth)

            ! Skip if both ends are dry
            IF (height_up < BIFMIN .and. height_dn < BIFMIN) THEN
               hflux_lev = 0._r8
               IF (do_update_state) THEN
                  pth_momen(ilev, ipth) = 0._r8
                  pth_veloc(ilev, ipth) = 0._r8
               ENDIF
               CYCLE
            ENDIF

            ! Upstream velocity from pathway state; downstream velocity = 0
            v_up = pth_veloc(ilev, ipth)
            v_dn = 0._r8

            ! --- Two-rarefaction Riemann solver ---

            ! Middle-state velocity and height
            veloct_fc = 0.5_r8 * (v_up + v_dn) &
               + sqrt(grav * height_up) - sqrt(grav * height_dn)

            height_fc = 1._r8/grav * (0.5_r8*(sqrt(grav*height_up) + sqrt(grav*height_dn)) &
               + 0.25_r8 * (v_up - v_dn)) ** 2

            ! Wave speeds
            IF (height_up > 0._r8) THEN
               vwave_up = min(v_up - sqrt(grav*height_up), veloct_fc - sqrt(grav*height_fc))
            ELSE
               vwave_up = v_dn - 2.0_r8 * sqrt(grav*height_dn)
            ENDIF

            IF (height_dn > 0._r8) THEN
               vwave_dn = max(v_dn + sqrt(grav*height_dn), veloct_fc + sqrt(grav*height_fc))
            ELSE
               vwave_dn = v_up + 2.0_r8 * sqrt(grav*height_up)
            ENDIF

            ! Fluxes at left and right states
            hflux_up = v_up  * height_up
            hflux_dn = v_dn  * height_dn
            mflux_up = v_up**2  * height_up + 0.5_r8*grav * height_up**2
            mflux_dn = v_dn**2  * height_dn + 0.5_r8*grav * height_dn**2

            ! Select flux based on wave structure, scaled by pathway width
            IF (vwave_up >= 0._r8) THEN
               hflux_lev = width_pth * hflux_up
               mflux_lev = width_pth * mflux_up
            ELSEIF (vwave_dn <= 0._r8) THEN
               hflux_lev = width_pth * hflux_dn
               mflux_lev = width_pth * mflux_dn
            ELSE
               hflux_lev = width_pth * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                  + vwave_up*vwave_dn*(height_dn - height_up)) / (vwave_dn - vwave_up)
               mflux_lev = width_pth * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                  + vwave_up*vwave_dn*(hflux_dn - hflux_up)) / (vwave_dn - vwave_up)
            ENDIF

            IF (do_update_state) THEN
               ! --- Update pathway momentum (semi-implicit friction, same as main channel) ---
               ! pth_are = pathway top area (width * pathway length)
               ! Momentum equation: d(h*v)/dt = -d(h*v^2 + gh^2/2)/dx - friction
               ! Discrete: momen_new = (momen_old - mflux/are*dt) / (1 + friction*dt)
               pth_are = width_pth * pth_dst(ipth)

               ! Both ends < BIFMIN was already CYCLE'd above, so max(...) >= BIFMIN
               ! holds here unconditionally.
               friction = grav * pth_man(ilev)**2 &
                  / max(height_up, height_dn)**(7._r8/3._r8) * abs(pth_momen(ilev, ipth))
               pth_momen(ilev, ipth) = (pth_momen(ilev, ipth) &
                  - mflux_lev / pth_are * dt) &
                  / (1._r8 + friction * dt)
               pth_veloc(ilev, ipth) = pth_momen(ilev, ipth) / max(height_up, height_dn)

               ! Clamp velocity and keep momentum in sync so the next substep
               ! reads a consistent (v_up, pth_momen) pair: v_up comes from
               ! pth_veloc and friction reads abs(pth_momen). Without the
               ! re-sync, a clamp-only update leaves |pth_momen|/max_h larger
               ! than the clamped |pth_veloc|, over-estimating friction on the
               ! next step.
               pth_veloc(ilev, ipth) = max(-20._r8, min(20._r8, pth_veloc(ilev, ipth)))
               pth_momen(ilev, ipth) = pth_veloc(ilev, ipth) * max(height_up, height_dn)
            ENDIF

            ! Accumulate total flux for this pathway [m^3/s]
            pth_hflux_total(ipth) = pth_hflux_total(ipth) + hflux_lev
            bif_hflux_lev(ilev, ipth) = hflux_lev

         ENDDO  ! ilev

         ! ----- Step 4: 5% storage limiter -----
         IF (abs(pth_hflux_total(ipth)) > 1.e-10_r8 .and. dt > 0._r8) THEN
            storage_up = storage_ucat(i_up)
            storage_dn = storage_dn_pth(ipth)
            ! Do not promote nearly dry pathway endpoints to an artificial
            ! 1 m3 storage. With a globally synchronized routing DT, that
            ! floor lets dry bifurcation paths force pathological tiny DTs.
            storage_ref = min(storage_up, storage_dn)
            storage_ref = max(storage_ref, BIF_STORAGE_EPS)
            pth_limiter_rate(ipth) = min(1._r8, 0.05_r8 * storage_ref / (abs(pth_hflux_total(ipth)) * dt))
         ENDIF

         ! Source-side accumulation only (i_up is always local).
         ! Destination-side contributions are added below via cross-rank push,
         ! covering both local and remote pathways feeding this cell.
         path_transfer = abs(pth_hflux_total(ipth))
         IF (pth_hflux_total(ipth) >= 0._r8) THEN
            limiter_outgoing(i_up) = limiter_outgoing(i_up) + path_transfer
            pos_pth(ipth) = pth_hflux_total(ipth)
         ELSE
            limiter_incoming(i_up) = limiter_incoming(i_up) + path_transfer
            neg_pth(ipth) = -pth_hflux_total(ipth)
         ENDIF

      ENDDO  ! ipth

      n_dt_mismatch_skip_glb = n_dt_mismatch_skip
#ifdef USEMPI
      CALL mpi_allreduce (n_dt_mismatch_skip, n_dt_mismatch_skip_glb, 1, MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
#endif
      IF (p_iam_worker == p_root .and. n_dt_mismatch_skip_glb > 0 .and. n_dt_mismatch_warn_count < 5) THEN
         write(*,'(A,I0,A)') 'WARNING bifurcation_calc: skipped ', n_dt_mismatch_skip_glb, &
            ' pathway(s) because source/destination routing dt differ; cross-river-system bifurcation may be inactive.'
         n_dt_mismatch_warn_count = n_dt_mismatch_warn_count + 1
      ENDIF
      ! ----- cross-rank destination accumulation -----
      ! Push pos/neg pathway flux components to their destination cells.
      ! sum mode collects contributions from all ranks owning the source.
      CALL worker_push_data (push_bif_influx, pos_pth, bif_pos_recv, &
         fillvalue = 0._r8, mode = 'sum')
      CALL worker_push_data (push_bif_influx, neg_pth, bif_neg_recv, &
         fillvalue = 0._r8, mode = 'sum')

      ! Each cell now sees the full destination-side contribution from
      ! every pathway pointing AT it (local and remote combined).
      DO i_ucat = 1, numucat
         limiter_incoming(i_ucat) = limiter_incoming(i_ucat) + bif_pos_recv(i_ucat)
         limiter_outgoing(i_ucat) = limiter_outgoing(i_ucat) + bif_neg_recv(i_ucat)
      ENDDO

      DO i_ucat = 1, numucat
         dt_cell = dt_ucat(i_ucat)
         storage_ref = max(storage_ucat(i_ucat), BIF_STORAGE_EPS)
         IF (limiter_outgoing(i_ucat) > 1.e-10_r8 .and. dt_cell > 0._r8) THEN
            limiter_out_rate(i_ucat) = min(1._r8, 0.05_r8 * storage_ref / (limiter_outgoing(i_ucat) * dt_cell))
         ENDIF
         IF (limiter_incoming(i_ucat) > 1.e-10_r8 .and. dt_cell > 0._r8) THEN
            limiter_in_rate(i_ucat) = min(1._r8, 0.05_r8 * storage_ref / (limiter_incoming(i_ucat) * dt_cell))
         ENDIF
      ENDDO

      ! ----- push destination cell rates back to each pathway -----
      ! For pathways whose downstream cell is on a remote rank, the cell rate
      ! is computed there; we pull it via push_bif_dn2pth. For local-downstream
      ! pathways the push delivers the same value as a direct lookup.
      ! Fillvalue=1 keeps boundary pathways unconstrained from this side.
      CALL worker_push_data (push_bif_dn2pth, limiter_in_rate,  limiter_in_rate_pth,  fillvalue = 1._r8)
      CALL worker_push_data (push_bif_dn2pth, limiter_out_rate, limiter_out_rate_pth, fillvalue = 1._r8)

      DO ipth = 1, npthout_local

         i_up = pth_upst_local(ipth)
         IF (i_up < 1 .or. i_up > numucat) CYCLE
         IF (.not. ucatfilter(i_up)) CYCLE
         IF (i_up <= size(is_built_resv) .and. is_built_resv(i_up)) CYCLE
         IF (wdsrf_dn_pth(ipth) < -9000._r8 .or. rivelv_dn_pth(ipth) < -9000._r8) CYCLE
         ! Zero-length pathways are skipped up in the first DO ipth loop, so
         ! pth_hflux_total(ipth) stays 0 here and this loop's scaling is a
         ! no-op. No need to re-check pth_dst.

         i_dn = pth_down_local(ipth)
         ! mirror the first DO loop's cross-system sync guard here
         ! so we never scale pth_momen / pth_veloc by a rate that the
         ! upstream system's other pathways produced for a pathway whose
         ! destination is not participating in this sub-step.
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            IF (.not. ucatfilter(i_dn)) CYCLE
            IF (i_dn <= size(is_built_resv) .and. is_built_resv(i_dn)) CYCLE
            dt_dn = dt_ucat(i_dn)
         ELSE
            IF (ucatfilter_dn_pth(ipth) < 0.5_r8) CYCLE
            IF (is_resv_dn_pth(ipth) > 0.5_r8) CYCLE
            dt_dn = dt_dn_pth(ipth)
         ENDIF
         dt = dt_ucat(i_up)
         IF (dt <= 0._r8 .or. dt_dn <= 0._r8) CYCLE
         IF (abs(dt - dt_dn) > max(1.e-9_r8, 1.e-12_r8 * max(abs(dt), abs(dt_dn)))) CYCLE
         ! Source rate is always local; destination rate comes from push above.
         IF (pth_hflux_total(ipth) >= 0._r8) THEN
            pth_limiter_rate(ipth) = min(pth_limiter_rate(ipth), limiter_out_rate(i_up))
            pth_limiter_rate(ipth) = min(pth_limiter_rate(ipth), limiter_in_rate_pth(ipth))
         ELSE
            pth_limiter_rate(ipth) = min(pth_limiter_rate(ipth), limiter_in_rate(i_up))
            pth_limiter_rate(ipth) = min(pth_limiter_rate(ipth), limiter_out_rate_pth(ipth))
         ENDIF

         rate = pth_limiter_rate(ipth)
         IF (rate < 1._r8) THEN
            pth_hflux_total(ipth) = pth_hflux_total(ipth) * rate
            bif_hflux_lev(:, ipth) = bif_hflux_lev(:, ipth) * rate
            IF (do_update_state) THEN
               pth_momen(:, ipth) = pth_momen(:, ipth) * rate
               pth_veloc(:, ipth) = pth_veloc(:, ipth) * rate
            ENDIF
         ENDIF

         ! ----- Step 5: Accumulate to upstream ucat (local) -----
         bif_hflux_sum(i_up) = bif_hflux_sum(i_up) + pth_hflux_total(ipth)

         IF (npthlev_bif >= 2) THEN
            pth_lev_hflux_total(ipth) = sum(bif_hflux_lev(2:npthlev_bif, ipth))
            bif_lev_hflux_sum(i_up) = bif_lev_hflux_sum(i_up) + pth_lev_hflux_total(ipth)
         ENDIF

         ! Also accumulate directly to downstream ucat if local
         ! (i_dn was set above)
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            bif_hflux_sum(i_dn) = bif_hflux_sum(i_dn) - pth_hflux_total(ipth)
            IF (npthlev_bif >= 2) THEN
               bif_lev_hflux_sum(i_dn) = bif_lev_hflux_sum(i_dn) - pth_lev_hflux_total(ipth)
            ENDIF
         ENDIF

      ENDDO  ! ipth

      ! ----- Step 6: Scatter flux to remote downstream ucats -----
      IF (numucat > 0) THEN
         bif_influx(:) = 0._r8
         bif_lev_influx(:) = 0._r8
      ENDIF
      CALL worker_push_data (push_bif_influx, pth_hflux_total, bif_influx, &
         fillvalue = 0._r8, mode = 'sum')
      IF (npthlev_bif >= 2) THEN
         CALL worker_push_data (push_bif_influx, pth_lev_hflux_total, bif_lev_influx, &
            fillvalue = 0._r8, mode = 'sum')
      ENDIF

      ! ----- Step 7: Subtract remote inflow, avoiding double-counting -----
      ! Remove contributions from local pathways (already counted in step 5)
      DO ipth = 1, npthout_local
         i_dn = pth_down_local(ipth)
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            bif_influx(i_dn) = bif_influx(i_dn) - pth_hflux_total(ipth)
            IF (npthlev_bif >= 2) THEN
               bif_lev_influx(i_dn) = bif_lev_influx(i_dn) - pth_lev_hflux_total(ipth)
            ENDIF
         ENDIF
      ENDDO

      ! Apply remote inflow to bif_hflux_sum
      DO i_ucat = 1, numucat
         bif_hflux_sum(i_ucat) = bif_hflux_sum(i_ucat) - bif_influx(i_ucat)
      ENDDO
      IF (npthlev_bif >= 2) THEN
         DO i_ucat = 1, numucat
            bif_lev_hflux_sum(i_ucat) = bif_lev_hflux_sum(i_ucat) - bif_lev_influx(i_ucat)
         ENDDO
      ENDIF

      ! Disassociate pointer aliases (targets remain allocated, owned by module)
      nullify (wdsrf_dn_pth, rivelv_dn_pth, protected_wdsrf_ucat, protected_wdsrf_dn_pth, has_levee_dn_pth, &
               storage_dn_pth, pth_hflux_total, pth_limiter_rate, bif_influx, &
               has_levee_r8, storage_ucat, limiter_outgoing, limiter_incoming, &
               limiter_out_rate, limiter_in_rate, &
               pos_pth, neg_pth, bif_pos_recv, bif_neg_recv, &
               limiter_in_rate_pth, limiter_out_rate_pth, &
               ucatfilter_r8, ucatfilter_dn_pth, is_resv_r8, is_resv_dn_pth, &
               dt_ucat, dt_dn_pth, &
               pth_lev_hflux_total, bif_lev_influx)

   END SUBROUTINE bifurcation_calc


   ! =========================================================================
   FUNCTION available_storage_ucat (i, wdsrf, volresv_in, is_built_resv_in) RESULT(storage)
   ! =========================================================================
   !
   ! Current host-model available storage semantics for limiter reference.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
	   real(r8), intent(in) :: volresv_in(:)
	   logical,  intent(in) :: is_built_resv_in(:)
	   real(r8)             :: storage
	   logical              :: has_levee_cell
	   real(r8), parameter  :: stage_restart_tol = 1.e-5_r8

      storage = 0._r8

      IF (i < 1) RETURN

      ! Size-check is_built_resv_in before dereferencing it; the dummy is
      ! caller-provided and we cannot rely on a fixed length contract.
      IF (allocated(ucat2resv) .and. i <= size(is_built_resv_in)) THEN
         IF (is_built_resv_in(i)) THEN
            IF (i <= size(ucat2resv)) THEN
               IF (ucat2resv(i) > 0 .and. ucat2resv(i) <= size(volresv_in)) THEN
                  storage = volresv_in(ucat2resv(i))
                  RETURN
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      has_levee_cell = .false.
      IF (DEF_USE_LEVEE .and. allocated(has_levee)) THEN
         IF (i <= size(has_levee)) has_levee_cell = has_levee(i)
      ENDIF

	      IF (volwater_ucat_valid .and. allocated(volwater_ucat)) THEN
	         IF (i <= size(volwater_ucat)) THEN
	            IF (volwater_ucat(i) > 0._r8 .or. wdsrf <= stage_restart_tol) THEN
	               storage = volwater_ucat(i)
	               IF (has_levee_cell) storage = storage + MAX(levsto(i), 0._r8)
	            ELSEIF (has_levee_cell) THEN
	               ! Old restart compatibility: some restarts store
	               ! volwater_ucat as zero placeholders for wet levee cells.
	               storage = levee_visible_volume_from_stage(i, wdsrf, levsto(i)) + MAX(levsto(i), 0._r8)
	            ELSE
	               storage = storage_with_river_prism(i, wdsrf)
	            ENDIF
	         ELSE
	            storage = storage_with_river_prism(i, wdsrf)
	         ENDIF
	      ELSEIF (has_levee_cell .and. allocated(levsto)) THEN
	         storage = levee_visible_volume_from_stage(i, wdsrf, levsto(i)) + MAX(levsto(i), 0._r8)
	      ELSE
	         storage = storage_with_river_prism(i, wdsrf)
	      ENDIF

   END FUNCTION available_storage_ucat


   ! =========================================================================
   FUNCTION storage_with_river_prism (i, wdsrf) RESULT(storage)
   ! =========================================================================
   !
   ! Limiter storage reference using the same above-bank river-prism convention
   ! as levee visible/total storage.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
   real(r8)             :: storage

      storage = floodplain_curve(i)%volume(wdsrf)
      IF (wdsrf > floodplain_curve(i)%rivhgt) THEN
         storage = storage + floodplain_curve(i)%rivare * (wdsrf - floodplain_curve(i)%rivhgt)
      ENDIF
      storage = MAX(storage, 0._r8)

   END FUNCTION storage_with_river_prism


   ! =========================================================================
   SUBROUTINE read_bifurcation_restart (file_restart)
   ! =========================================================================
   !
   ! Read bifurcation pathway state from restart in global pathway order.
   !
   ! =========================================================================

   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_Vector_ReadWrite, only: vector_read_matrix_and_scatter
   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical :: has_pth_veloc, has_pth_momen
   integer, allocatable :: global_id_read(:)
   real(r8), allocatable :: matrix_dummy(:,:)
   integer :: has_flags(2)

      ! vector_read_matrix_and_scatter contains mpi_barrier(p_comm_glb)
      ! and mpi_recv/mpi_send on p_comm_glb, so every rank in that
      ! communicator (master, workers AND IO ranks) must enter this
      ! routine or the barrier hangs forever. The IO-only ranks take the
      ! "ELSE" branch below along with the master-when-non-worker case,
      ! passing zero-length dummy arrays and participating only in the
      ! collectives. Do NOT short-circuit with `p_is_master .or.
      ! p_is_worker` here.
      IF (npthlev_bif <= 0 .or. totalnpthout <= 0) RETURN

      ! Master probes restart variable existence then broadcasts the
      ! result. Having all ranks call ncio_var_exist directly would mean
      ! p_np_glb concurrent nf90_open calls on the same file; on an
      ! HDF5-backed NetCDF4 build the default file locking serialises or
      ! deadlocks those opens. One master-only probe + broadcast avoids
      ! the storm entirely and keeps all ranks on the same boolean.
      IF (p_is_master) THEN
         has_pth_veloc = ncio_var_exist(file_restart, 'pth_veloc', readflag = .false.)
         has_pth_momen = ncio_var_exist(file_restart, 'pth_momen', readflag = .false.)
         has_flags(1) = merge(1, 0, has_pth_veloc)
         has_flags(2) = merge(1, 0, has_pth_momen)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (has_flags, 2, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      has_pth_veloc = (has_flags(1) /= 0)
      has_pth_momen = (has_flags(2) /= 0)

      if (has_pth_veloc .and. has_pth_momen) then
         IF (p_is_worker) THEN
            CALL vector_read_matrix_and_scatter ( &
               file_restart, pth_veloc, npthlev_bif, npthout_local, 'pth_veloc', pth_global_id, totalnpthout)
            CALL vector_read_matrix_and_scatter ( &
               file_restart, pth_momen, npthlev_bif, npthout_local, 'pth_momen', pth_global_id, totalnpthout)
         ELSE
            allocate (global_id_read(0))
            CALL vector_read_matrix_and_scatter ( &
               file_restart, matrix_dummy, npthlev_bif, 0, 'pth_veloc', global_id_read, totalnpthout)
            CALL vector_read_matrix_and_scatter ( &
               file_restart, matrix_dummy, npthlev_bif, 0, 'pth_momen', global_id_read, totalnpthout)
            deallocate (global_id_read)
         ENDIF
      else
         ! Missing pathway variables mean the restart file predates bifurcation
         ! support, so keep the zero-initialized cold-start state for compatibility.
         IF (p_is_worker) THEN
            pth_veloc(:,:) = 0._r8
            pth_momen(:,:) = 0._r8
         ENDIF
      endif

      ! Sanitize loaded state: zero out NaN, Inf, and unphysically extreme
      ! values so a corrupt restart cannot poison subsequent bifurcation_calc
      ! steps. Velocity is clamped to +/-20 m/s every substep, so anything
      ! above 50 m/s on load is almost certainly garbage. Momentum = velocity
      ! * depth, and plausible depths are O(10 m), so 1e4 m^2/s is a generous
      ! cap. Sanitize (velocity, momentum) as a pair: if either field at a
      ! given (ilev, ipth) is bad, clear BOTH so the solver never sees an
      ! inconsistent (v_up, pth_momen) state on the first substep (line 353
      ! reads v_up = pth_veloc, line 406 reads friction from pth_momen --
      ! they must describe the same physical state).
      IF (p_is_worker .and. npthout_local > 0 .and. allocated(pth_veloc) .and. allocated(pth_momen)) THEN
         WHERE (pth_veloc /= pth_veloc .or. abs(pth_veloc) > 50._r8 .or. &
                pth_momen /= pth_momen .or. abs(pth_momen) > 1.e4_r8)
            pth_veloc = 0._r8
            pth_momen = 0._r8
         END WHERE
      ENDIF

      CALL debug_check_bifurcation_restart_state (has_pth_veloc, has_pth_momen)

   END SUBROUTINE read_bifurcation_restart

   SUBROUTINE debug_check_bifurcation_restart_state (has_pth_veloc, has_pth_momen)

   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   logical, intent(in) :: has_pth_veloc
   logical, intent(in) :: has_pth_momen

   integer :: nbad_gid
   integer :: nnan_veloc
   integer :: nnan_momen
   real(r8) :: maxabs_veloc
   real(r8) :: maxabs_momen

      IF (dbg_bif_restart_checked) RETURN
      dbg_bif_restart_checked = .true.

      IF (.not. p_is_worker) RETURN

      nbad_gid = 0
      nnan_veloc = 0
      nnan_momen = 0
      maxabs_veloc = 0._r8
      maxabs_momen = 0._r8

      IF (npthout_local > 0) THEN
         nbad_gid = count(pth_global_id < 1 .or. pth_global_id > totalnpthout)
         nnan_veloc = count(pth_veloc /= pth_veloc)
         nnan_momen = count(pth_momen /= pth_momen)
         maxabs_veloc = maxval(abs(pth_veloc))
         maxabs_momen = maxval(abs(pth_momen))
      ENDIF

      IF (nbad_gid > 0 .or. nnan_veloc > 0 .or. nnan_momen > 0) THEN
         write(*,'(A,I0,A,L1,A,L1,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)') &
            'WARNING bifurcation restart state glb=', p_iam_glb, &
            ' has_veloc=', has_pth_veloc, &
            ' has_momen=', has_pth_momen, &
            ' bad_gid=', nbad_gid, &
            ' nan_vel=', nnan_veloc, &
            ' nan_mom=', nnan_momen, &
            ' maxabs_vel=', maxabs_veloc, &
            ' maxabs_mom=', maxabs_momen
         call flush(6)
      ENDIF

   END SUBROUTINE debug_check_bifurcation_restart_state


   ! =========================================================================
   SUBROUTINE write_bifurcation_restart (file_restart)
   ! =========================================================================
   !
   ! Write bifurcation pathway state in global pathway order.
   !
   ! =========================================================================

   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_define_dimension, ncio_write_serial
   USE MOD_Vector_ReadWrite, only: vector_gather_matrix_to_master
   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: ncol_local_write
   integer, allocatable :: pth_global_id_write(:)
   real(r8), allocatable :: wdata(:,:)
   real(r8), allocatable :: pth_veloc_write(:,:), pth_momen_write(:,:)

      IF (totalnpthout <= 0 .or. npthlev_bif <= 0) RETURN

      ! Guard: bifurcation state may not be initialised (e.g. mkinidata).
      ! Write zero state so the restart file has valid variables.
      IF (p_is_worker) THEN
         IF (allocated(pth_veloc) .and. allocated(pth_momen)) THEN
            ncol_local_write = npthout_local
            pth_veloc_write = pth_veloc
            pth_momen_write = pth_momen
            allocate (pth_global_id_write (ncol_local_write))
            pth_global_id_write(:) = pth_global_id(:)
         ELSE
            ! Not initialised: contribute zero-state columns
            ncol_local_write = npthout_local
            allocate (pth_veloc_write     (npthlev_bif, ncol_local_write))
            allocate (pth_momen_write     (npthlev_bif, ncol_local_write))
            allocate (pth_global_id_write (ncol_local_write))
            pth_veloc_write = 0._r8
            pth_momen_write = 0._r8
            pth_global_id_write(:) = pth_global_id(:)
         ENDIF
      ELSE
         ncol_local_write = 0
         allocate (pth_veloc_write     (npthlev_bif, 0))
         allocate (pth_momen_write     (npthlev_bif, 0))
         allocate (pth_global_id_write (0))
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_define_dimension(file_restart, 'bifurcation_level',   npthlev_bif)
         CALL ncio_define_dimension(file_restart, 'bifurcation_pathway', totalnpthout)
      ENDIF

      ! Gather local (npthlev_bif, npthout_local) state into global
      ! (npthlev_bif, totalnpthout) pathway order using pth_global_id.
      CALL vector_gather_matrix_to_master ( &
         pth_veloc_write, npthlev_bif, ncol_local_write, totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'pth_veloc', wdata, &
            'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF

      CALL vector_gather_matrix_to_master ( &
         pth_momen_write, npthlev_bif, ncol_local_write, totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'pth_momen', wdata, &
            'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF

      deallocate (pth_veloc_write)
      deallocate (pth_momen_write)
      deallocate (pth_global_id_write)

   END SUBROUTINE write_bifurcation_restart


   ! =========================================================================
   SUBROUTINE bifurcation_final ()
   ! =========================================================================
   !
   ! Deallocate all bifurcation arrays.
   !
   ! =========================================================================

   IMPLICIT NONE

      IF (allocated(pth_veloc))     deallocate (pth_veloc)
      IF (allocated(pth_momen))     deallocate (pth_momen)
      IF (allocated(bif_hflux_lev)) deallocate (bif_hflux_lev)
      IF (allocated(bif_hflux_sum)) deallocate (bif_hflux_sum)
      IF (allocated(bif_path_active)) deallocate (bif_path_active)
      IF (allocated(bif_lev_hflux_sum)) deallocate (bif_lev_hflux_sum)

      IF (allocated(wdsrf_dn_pth_buf    )) deallocate (wdsrf_dn_pth_buf    )
      IF (allocated(rivelv_dn_pth_buf   )) deallocate (rivelv_dn_pth_buf   )
      IF (allocated(protected_wdsrf_ucat_buf  )) deallocate (protected_wdsrf_ucat_buf  )
      IF (allocated(protected_wdsrf_dn_pth_buf)) deallocate (protected_wdsrf_dn_pth_buf)
      IF (allocated(has_levee_dn_pth_buf)) deallocate (has_levee_dn_pth_buf)
      IF (allocated(storage_dn_pth_buf  )) deallocate (storage_dn_pth_buf  )
      IF (allocated(pth_hflux_total_buf )) deallocate (pth_hflux_total_buf )
      IF (allocated(pth_limiter_rate_buf)) deallocate (pth_limiter_rate_buf)
      IF (allocated(bif_influx_buf      )) deallocate (bif_influx_buf      )
      IF (allocated(has_levee_r8_buf    )) deallocate (has_levee_r8_buf    )
      IF (allocated(storage_ucat_buf    )) deallocate (storage_ucat_buf    )
      IF (allocated(limiter_outgoing_buf)) deallocate (limiter_outgoing_buf)
      IF (allocated(limiter_incoming_buf)) deallocate (limiter_incoming_buf)
      IF (allocated(limiter_out_rate_buf)) deallocate (limiter_out_rate_buf)
      IF (allocated(limiter_in_rate_buf )) deallocate (limiter_in_rate_buf )

      IF (allocated(pos_pth_buf             )) deallocate (pos_pth_buf             )
      IF (allocated(neg_pth_buf             )) deallocate (neg_pth_buf             )
      IF (allocated(bif_pos_recv_buf        )) deallocate (bif_pos_recv_buf        )
      IF (allocated(bif_neg_recv_buf        )) deallocate (bif_neg_recv_buf        )
      IF (allocated(limiter_in_rate_pth_buf )) deallocate (limiter_in_rate_pth_buf )
      IF (allocated(limiter_out_rate_pth_buf)) deallocate (limiter_out_rate_pth_buf)

      IF (allocated(ucatfilter_r8_buf       )) deallocate (ucatfilter_r8_buf       )
      IF (allocated(ucatfilter_dn_pth_buf   )) deallocate (ucatfilter_dn_pth_buf   )
      IF (allocated(is_resv_r8_buf          )) deallocate (is_resv_r8_buf          )
      IF (allocated(is_resv_dn_pth_buf      )) deallocate (is_resv_dn_pth_buf      )
      IF (allocated(dt_ucat_buf             )) deallocate (dt_ucat_buf             )
      IF (allocated(dt_dn_pth_buf           )) deallocate (dt_dn_pth_buf           )

      IF (allocated(pth_lev_hflux_total_buf )) deallocate (pth_lev_hflux_total_buf )
      IF (allocated(bif_lev_influx_buf      )) deallocate (bif_lev_influx_buf      )

   END SUBROUTINE bifurcation_final

END MODULE MOD_Grid_RiverLakeBifurcation
#endif
