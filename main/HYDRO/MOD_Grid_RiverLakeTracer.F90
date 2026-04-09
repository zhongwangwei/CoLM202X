#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeTracer
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Passive tracer transport module for GridRiverLakeFlow.
!   Supports arbitrary number of tracers (e.g., delta18O, deltaD).
!
!   Physics:
!     - Upwind advection: tracer flux = concentration × water flux
!     - Flux limiter for mass conservation (prevents negative storage)
!     - Fully-mixed assumption within each unit catchment
!
!   Reference: CaMa-Flood tracer module (cmf_ctrl_tracer_mod.F90)
!
! Created by CoLM, Apr 2026
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE

   !-------------------------------------------------------------------------------------
   ! Module variables
   !-------------------------------------------------------------------------------------
   integer,  save :: ntracers             ! Number of tracers
   character(len=32), allocatable :: tracer_names(:)  ! Tracer names

   ! State variables (prognostic)
   real(r8), allocatable :: trc_mass  (:,:)   ! Tracer mass storage [mass] (ntracers, numucat)
   real(r8), allocatable :: trc_conc  (:,:)   ! Tracer concentration [mass/m3] (ntracers, numucat)

   ! Flux variables (diagnostic, per routing period)
   real(r8), allocatable :: trc_flux_out (:,:) ! Tracer outflux [mass/s] (ntracers, numucat)

   ! Input: tracer mass flux from runoff [mass/s] (ntracers, numucat)
   real(r8), allocatable :: acc_trc_inp  (:,:) ! Accumulated tracer input [mass] (ntracers, numucat)

   ! Bifurcation net flux (saved from last tracer_calc for diagnostics)
   real(r8), allocatable :: trc_bif_net_saved (:,:) ! Per-cell net bif flux [mass/s] (ntracers, numucat)

   ! History accumulators
   real(r8), allocatable :: a_trc_conc   (:,:) ! Accumulated tracer conc [mass/m3 * s] (ntracers, numucat)
   real(r8), allocatable :: a_trc_out    (:,:) ! Accumulated tracer outflux [mass/s * s] (ntracers, numucat)
   real(r8), allocatable :: a_trc_bifout (:,:) ! Accumulated tracer bif net flux [mass/s * s] (ntracers, numucat)

   PUBLIC :: tracer_init
   PUBLIC :: tracer_input_from_runoff
   PUBLIC :: tracer_calc
   PUBLIC :: tracer_diag_accumulate
   PUBLIC :: tracer_flush_acc
   PUBLIC :: read_tracer_restart
   PUBLIC :: write_tracer_restart
   PUBLIC :: write_tracer_history
   PUBLIC :: tracer_final

CONTAINS

   !-------------------------------------------------------------------------------------
   ! Parse comma-separated tracer names
   !-------------------------------------------------------------------------------------
   SUBROUTINE parse_tracer_names(namestr, names, n)

   IMPLICIT NONE
   character(len=*), intent(in) :: namestr
   character(len=32), intent(out) :: names(:)
   integer, intent(in) :: n

   integer :: i, istart, itrc
   character(len=256) :: str

      str = trim(namestr)
      istart = 1
      itrc = 0
      names(:) = ''

      DO i = 1, len_trim(str)
         IF (str(i:i) == ',') THEN
            itrc = itrc + 1
            IF (itrc <= n) names(itrc) = sanitize_ncname(str(istart:i-1))
            istart = i + 1
         ENDIF
      ENDDO
      ! last name
      itrc = itrc + 1
      IF (itrc <= n) names(itrc) = sanitize_ncname(str(istart:len_trim(str)))

      ! Validate: parsed count must match requested count
      IF (itrc /= n) THEN
         IF (p_is_master) THEN
            write(*,'(A,I0,A,I0,A)') ' WARNING: DEF_TRACER_NAMES has ', itrc, &
               ' tokens but DEF_TRACER_NUM = ', n, '. Filling missing with tracer_N.'
         ENDIF
         DO i = itrc+1, n
            write(names(i), '(A,I0)') 'tracer_', i
         ENDDO
      ENDIF

   END SUBROUTINE parse_tracer_names


   !-------------------------------------------------------------------------------------
   ! Sanitize a string for use as a NetCDF variable name component.
   ! Keep only alphanumeric characters, underscore, and period.
   !-------------------------------------------------------------------------------------
   FUNCTION sanitize_ncname(raw) RESULT(clean)
   IMPLICIT NONE
   character(len=*), intent(in) :: raw
   character(len=32) :: clean
   integer :: i, j
   character(len=1) :: c

      clean = ''
      j = 0
      DO i = 1, len_trim(raw)
         c = raw(i:i)
         IF ((c >= 'A' .and. c <= 'Z') .or. (c >= 'a' .and. c <= 'z') .or. &
             (c >= '0' .and. c <= '9') .or. c == '_' .or. c == '.') THEN
            j = j + 1
            IF (j <= 32) clean(j:j) = c
         ENDIF
      ENDDO

      ! Guard against empty result
      IF (len_trim(clean) == 0) clean = 'unnamed'

   END FUNCTION sanitize_ncname


   !-------------------------------------------------------------------------------------
   ! Initialize tracer module
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_init ()

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   integer :: i, j

      ntracers = DEF_TRACER_NUM

      allocate (tracer_names(ntracers))
      CALL parse_tracer_names(DEF_TRACER_NAMES, tracer_names, ntracers)

      ! Ensure unique names: append index if duplicates found after sanitization
      DO i = 1, ntracers
         DO j = 1, i-1
            IF (trim(tracer_names(i)) == trim(tracer_names(j))) THEN
               write(tracer_names(i), '(A,A,I0)') trim(tracer_names(i)), '_', i
               EXIT
            ENDIF
         ENDDO
      ENDDO

      IF (p_is_worker) THEN
         ! Allocate on ALL workers (zero-length if numucat=0) for MPI safety.
         allocate (trc_mass     (ntracers, numucat))
         allocate (trc_conc     (ntracers, numucat))
         allocate (trc_flux_out (ntracers, numucat))
         allocate (acc_trc_inp        (ntracers, numucat))
         allocate (trc_bif_net_saved  (ntracers, numucat))
         allocate (a_trc_conc         (ntracers, numucat))
         allocate (a_trc_out          (ntracers, numucat))
         allocate (a_trc_bifout       (ntracers, numucat))

         trc_mass     = 0._r8
         trc_conc     = 0._r8
         trc_flux_out = 0._r8
         acc_trc_inp  = 0._r8
         trc_bif_net_saved = 0._r8
         a_trc_conc   = 0._r8
         a_trc_out    = 0._r8
         a_trc_bifout = 0._r8
      ENDIF

      IF (p_is_master) THEN
         write(*,'(A,I4,A)') ' Tracer module initialised with ', ntracers, ' tracers:'
         write(*,'(A,*(A,:,", "))') '   Names: ', (trim(tracer_names(i)), i=1, ntracers)
      ENDIF

   END SUBROUTINE tracer_init


   !-------------------------------------------------------------------------------------
   ! Accumulate tracer input associated with runoff.
   ! Called each land-model timestep (before routing accumulation threshold).
   !
   ! trc_conc_runoff(itrc, i): concentration of tracer itrc in runoff at ucat i [mass/m3]
   ! rnof_uc(i):               runoff volume accumulated this step [m3]
   !   (= rnof_uc_raw * 1.e-3 * deltime, same as acc_rnof_uc increment)
   !
   ! Tracer mass input = conc_in_runoff × runoff_volume
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_input_from_runoff (rnof_uc_vol, numucat_in)

   IMPLICIT NONE
   real(r8), intent(in) :: rnof_uc_vol(:)  ! Runoff volume this step [m3]
   integer,  intent(in) :: numucat_in

   integer :: i, itrc
   real(r8) :: conc_d18O, conc_dD

      ! For testing: assign synthetic isotope concentrations to runoff.
      ! delta18O ~ -10 permil => R/Rsmow = 0.990 => conc ~ 0.990 * Rsmow
      ! We store delta values directly in permil for simplicity.
      ! In production, this would come from the land surface model.
      conc_d18O = -10.0_r8   ! permil (VSMOW)
      conc_dD   = -70.0_r8   ! permil (VSMOW)

      DO i = 1, numucat_in
         IF (rnof_uc_vol(i) > 0._r8) THEN
            ! Tracer "mass" = delta_value × water_volume
            ! This preserves mass-weighted mixing:
            !   mixed_delta = sum(delta_i * vol_i) / sum(vol_i)
            DO itrc = 1, ntracers
               IF (itrc == 1) THEN
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + conc_d18O * rnof_uc_vol(i)
               ELSEIF (itrc == 2) THEN
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + conc_dD * rnof_uc_vol(i)
               ELSE
                  ! Additional tracers: zero input by default
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + 0._r8
               ENDIF
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE tracer_input_from_runoff


   !-------------------------------------------------------------------------------------
   ! Main tracer transport calculation.
   ! Called once per routing period (after all sub-timesteps of water routing).
   !
   ! Uses time-averaged discharge (hflux_fc) from the routing sub-steps.
   ! This is a single-step explicit advection over the full routing period,
   ! using the effective (time-averaged) water fluxes.
   !
   ! Algorithm (following CaMa-Flood):
   !   1. Add runoff tracer input to storage
   !   2. Compute concentration = mass / volume
   !   3. Compute upwind tracer flux = conc × water_flux
   !   4. Flux-limit to prevent negative mass
   !   5. Update mass: mass += (inflow - outflow) * dt
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_calc (acctime, wdsrf, veloc, hflux_avg, &
      numucat_in, ucatfilter_in, volresv_in, ucat2resv_in, &
      do_bifurcation, prd_bifflw_lev, prd_bifflw_time)
   !
   ! Joint main-channel + bifurcation tracer transport (CaMa-Flood semantics).
   ! All fluxes are computed from the SAME concentration snapshot, then a
   ! shared flux limiter is applied, and mass is updated in one pass.
   !
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, &
      floodplain_curve, topo_area, lake_type, &
      push_ups2ucat, push_next2ucat, &
      npthout_local, npthlev_bif, pth_upst_local, pth_down_local, &
      push_bif_influx, push_bif_dn2pth
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, volwater_ucat
   USE MOD_WorkerPushData
   USE MOD_Vars_Global,           only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: acctime             ! Routing period length [s]
   real(r8), intent(in) :: wdsrf(:)            ! Final water depth [m]
   real(r8), intent(in) :: veloc(:)            ! Final velocity [m/s]
   real(r8), intent(in) :: hflux_avg(:)        ! Time-averaged discharge [m3/s]
   integer,  intent(in) :: numucat_in
   logical,  intent(in) :: ucatfilter_in(:)
   real(r8), intent(in) :: volresv_in(:)       ! Reservoir volumes [m3]
   integer,  intent(in) :: ucat2resv_in(:)     ! ucat → reservoir index mapping
   logical,  intent(in) :: do_bifurcation      ! Whether to include bifurcation
   real(r8), intent(in) :: prd_bifflw_lev(:,:) ! Per-period bif pathway flux [m3/s * s]
   real(r8), intent(in) :: prd_bifflw_time(:)  ! Per-period bif pathway time [s]

   ! Local variables
   integer :: i, itrc, ipth, i_up, i_dn
   real(r8) :: volwater, pth_flux_avg, conc_up, conc_dn
   real(r8), allocatable :: trc_conc_next(:,:)   ! Downstream concentration (main)
   real(r8), allocatable :: trc_flux_ups(:,:)    ! Upstream tracer flux sum (main)
   real(r8), allocatable :: trc_out_mass(:)      ! Total outgoing tracer mass per cell
   real(r8), allocatable :: rate_cell(:)         ! Limiter rate per cell (reused per tracer)
   real(r8), allocatable :: rate_next(:)         ! Downstream cell's rate (via push_next2ucat)
   ! Bifurcation locals
   real(r8), allocatable :: trc_pth_flux(:,:)    ! Per-pathway tracer flux
   real(r8), allocatable :: trc_bif_influx(:,:)  ! Remote bif tracer inflow
   real(r8), allocatable :: trc_pth_1trc(:)      ! Single-tracer pathway flux for push
   real(r8), allocatable :: conc_dn_pth(:)       ! Downstream concentration per pathway
   real(r8), allocatable :: rate_dn_pth(:)       ! Downstream rate per pathway (via push_bif_dn2pth)
   real(r8), allocatable :: trc_bif_net(:,:)     ! Net bif tracer flux per cell

      IF (acctime <= 0._r8) RETURN

      ! Allocate on all workers (zero-length if numucat=0) for MPI safety.
      allocate (trc_conc_next (ntracers, numucat_in))
      allocate (trc_flux_ups  (ntracers, numucat_in))
      allocate (trc_out_mass  (numucat_in))
      allocate (rate_cell     (numucat_in))
      allocate (rate_next     (numucat_in))
      allocate (trc_bif_net   (ntracers, numucat_in))
      trc_bif_net = 0._r8

      IF (do_bifurcation .and. npthlev_bif > 0) THEN
         allocate (trc_pth_flux   (ntracers, npthout_local))
         allocate (trc_bif_influx (ntracers, numucat_in))
         allocate (trc_pth_1trc   (npthout_local))
         allocate (conc_dn_pth    (npthout_local))
         allocate (rate_dn_pth    (npthout_local))
         trc_pth_flux   = 0._r8
         trc_bif_influx = 0._r8
      ENDIF

      ! ===== Step 1: Add accumulated runoff tracer input to storage =====
      DO i = 1, numucat_in
         DO itrc = 1, ntracers
            trc_mass(itrc, i) = trc_mass(itrc, i) + acc_trc_inp(itrc, i)
         ENDDO
      ENDDO

      ! ===== Step 2: Concentration snapshot (frozen for all flux calcs) =====
      DO i = 1, numucat_in
         CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
         DO itrc = 1, ntracers
            trc_conc(itrc, i) = trc_mass(itrc, i) / volwater
         ENDDO
      ENDDO

      ! ===== Step 3: Main-channel flux (from snapshot) =====
      ! 3a: Get downstream concentration for reverse flow
      DO itrc = 1, ntracers
         CALL worker_push_data (push_next2ucat, trc_conc(itrc,:), trc_conc_next(itrc,:), fillvalue = 0._r8)
      ENDDO

      ! 3b: Upwind tracer outflux
      DO i = 1, numucat_in
         DO itrc = 1, ntracers
            IF (hflux_avg(i) >= 0._r8) THEN
               trc_flux_out(itrc, i) = trc_conc(itrc, i) * hflux_avg(i)
            ELSE
               trc_flux_out(itrc, i) = trc_conc_next(itrc, i) * hflux_avg(i)
            ENDIF
         ENDDO
      ENDDO

      ! ===== Step 4: Bifurcation pathway flux (from same snapshot) =====
      IF (do_bifurcation .and. npthlev_bif > 0) THEN
         DO itrc = 1, ntracers
            CALL worker_push_data (push_bif_dn2pth, trc_conc(itrc,:), conc_dn_pth, fillvalue = 0._r8)

            DO ipth = 1, npthout_local
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat_in) CYCLE
               IF (prd_bifflw_time(ipth) <= 0._r8) CYCLE

               pth_flux_avg = sum(prd_bifflw_lev(:, ipth)) / prd_bifflw_time(ipth)
               conc_up = trc_conc(itrc, i_up)
               conc_dn = conc_dn_pth(ipth)

               IF (pth_flux_avg >= 0._r8) THEN
                  trc_pth_flux(itrc, ipth) = conc_up * pth_flux_avg
               ELSE
                  trc_pth_flux(itrc, ipth) = conc_dn * pth_flux_avg
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      ! ===== Step 5: CaMa-style shared limiter (pathway-level) =====
      !
      ! CaMa P2STOOUT/D2RATE semantics:
      !   P2STOOUT(i) = max(main_out(i), 0)          -- positive main outflow
      !               + sum_upstream_j max(-main_out(j), 0)  -- reverse inflow from upstream = outflow from i
      !               + sum_pathways max(pth_flux, 0)  -- positive bif outflow from i
      !               + sum_pathways max(-pth_flux, 0) -- reverse bif inflow = outflow from i
      !   D2RATE(i) = min(storage(i) / P2STOOUT(i), 1)
      !
      ! Then scale:
      !   positive main_out(i) *= D2RATE(i)
      !   negative main_out(i) *= D2RATE(downstream_of_i)  -- sender is downstream
      !   positive pth_flux(ipth) *= D2RATE(upstream_cell)
      !   negative pth_flux(ipth) *= D2RATE(downstream_cell)  -- including remote
      !
      DO itrc = 1, ntracers

         ! 5a: Compute P2STOOUT — total sender-side outflow per cell (no cancellation)
         DO i = 1, numucat_in
            ! Main channel: positive outflux from this cell
            trc_out_mass(i) = max(trc_flux_out(itrc, i), 0._r8) * acctime
         ENDDO

         ! Main channel: reverse inflow from upstream cells counts as outflow from this cell
         ! (CaMa line 591-595: P2STOOUT(ISEQ) += max(-D2TRCOUT(JSEQ), 0) for upstream JSEQ)
         ! This is equivalent to: for each cell i, if its upstream j has reverse flow
         ! (trc_flux_out(j) < 0), then mass leaves i toward j, counted in P2STOOUT(i).
         ! We handle this by using push_ups2ucat to aggregate max(-flux, 0) from upstream.
         ! But simpler: for cell i with trc_flux_out(i) < 0, the donor is ucat_next(i).
         ! So we add max(-trc_flux_out(i), 0) to P2STOOUT(ucat_next(i)).
         ! Since we can't directly index ucat_next here, we use push_ups2ucat.
         ! Actually simpler: just accumulate on a temporary and push.
         ! For now, use the direct CaMa approach: iterate upstream neighbors.
         ! In our push-based model, we need a helper push. Let's use push_next2ucat
         ! to send the reverse-flow contribution to the downstream cell.
         ! Scratch: use trc_flux_ups temporarily for this.
         DO i = 1, numucat_in
            ! How much mass leaves this cell's downstream neighbor due to reverse flow at i?
            ! If trc_flux_out(i) < 0: downstream cell donates |flux|*dt of tracer.
            rate_cell(i) = max(-trc_flux_out(itrc, i), 0._r8) * acctime
         ENDDO
         ! Push reverse-out contribution to the downstream cell (which is the actual sender)
         CALL worker_push_data (push_ups2ucat, rate_cell, rate_next, fillvalue = 0._r8, mode = 'sum')
         DO i = 1, numucat_in
            trc_out_mass(i) = trc_out_mass(i) + rate_next(i)
         ENDDO

         ! Bifurcation pathways
         IF (do_bifurcation .and. npthlev_bif > 0) THEN
            DO ipth = 1, npthout_local
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat_in) CYCLE
               ! Positive flux: upstream is sender
               trc_out_mass(i_up) = trc_out_mass(i_up) &
                  + max(trc_pth_flux(itrc, ipth), 0._r8) * acctime
               ! Negative flux: downstream is sender
               i_dn = pth_down_local(ipth)
               IF (i_dn > 0 .and. i_dn <= numucat_in) THEN
                  trc_out_mass(i_dn) = trc_out_mass(i_dn) &
                     + max(-trc_pth_flux(itrc, ipth), 0._r8) * acctime
               ENDIF
            ENDDO
            ! Note: reverse bif outflow from remote downstream cells is not counted
            ! in their P2STOOUT here (they compute their own rate on their own worker).
         ENDIF

         ! 5b: Compute D2RATE per cell
         DO i = 1, numucat_in
            IF (trc_out_mass(i) > 1.e-30_r8) THEN
               rate_cell(i) = min(abs(trc_mass(itrc, i)) / trc_out_mass(i), 1._r8)
            ELSE
               rate_cell(i) = 1._r8
            ENDIF
         ENDDO

         ! 5c: Get downstream cell's rate (for reverse main flow scaling)
         CALL worker_push_data (push_next2ucat, rate_cell, rate_next, fillvalue = 1._r8)

         ! 5d: Scale main-channel flux
         DO i = 1, numucat_in
            IF (trc_flux_out(itrc, i) >= 0._r8) THEN
               trc_flux_out(itrc, i) = trc_flux_out(itrc, i) * rate_cell(i)
            ELSE
               ! Reverse: sender is downstream → use downstream rate
               trc_flux_out(itrc, i) = trc_flux_out(itrc, i) * rate_next(i)
            ENDIF
         ENDDO

         ! 5e: Get downstream pathway rate (for reverse bif flow, including remote)
         IF (do_bifurcation .and. npthlev_bif > 0) THEN
            CALL worker_push_data (push_bif_dn2pth, rate_cell, rate_dn_pth, fillvalue = 1._r8)

            DO ipth = 1, npthout_local
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat_in) CYCLE

               IF (trc_pth_flux(itrc, ipth) >= 0._r8) THEN
                  ! Normal: sender is upstream → use upstream rate
                  trc_pth_flux(itrc, ipth) = trc_pth_flux(itrc, ipth) * rate_cell(i_up)
               ELSE
                  ! Reverse: sender is downstream → use downstream rate (possibly remote)
                  trc_pth_flux(itrc, ipth) = trc_pth_flux(itrc, ipth) * rate_dn_pth(ipth)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      ! ===== Step 6: Aggregate fluxes after limiting =====
      ! 6a: Aggregate upstream main-channel tracer inflow
      DO itrc = 1, ntracers
         CALL worker_push_data (push_ups2ucat, trc_flux_out(itrc,:), trc_flux_ups(itrc,:), &
            fillvalue = 0._r8, mode = 'sum')
      ENDDO

      ! 6b: Rebuild per-cell bifurcation net from scaled pathway fluxes
      trc_bif_net = 0._r8
      IF (do_bifurcation .and. npthlev_bif > 0) THEN
         DO ipth = 1, npthout_local
            i_up = pth_upst_local(ipth)
            IF (i_up < 1 .or. i_up > numucat_in) CYCLE
            DO itrc = 1, ntracers
               trc_bif_net(itrc, i_up) = trc_bif_net(itrc, i_up) + trc_pth_flux(itrc, ipth)
            ENDDO
            i_dn = pth_down_local(ipth)
            IF (i_dn > 0 .and. i_dn <= numucat_in) THEN
               DO itrc = 1, ntracers
                  trc_bif_net(itrc, i_dn) = trc_bif_net(itrc, i_dn) - trc_pth_flux(itrc, ipth)
               ENDDO
            ENDIF
         ENDDO

         ! Scatter scaled flux to remote downstream
         DO itrc = 1, ntracers
            trc_pth_1trc(:) = trc_pth_flux(itrc, :)
            CALL worker_push_data (push_bif_influx, trc_pth_1trc, trc_bif_influx(itrc,:), &
               fillvalue = 0._r8, mode = 'sum')
         ENDDO
         ! Remove local double-count
         DO ipth = 1, npthout_local
            i_dn = pth_down_local(ipth)
            IF (i_dn > 0 .and. i_dn <= numucat_in) THEN
               DO itrc = 1, ntracers
                  trc_bif_influx(itrc, i_dn) = trc_bif_influx(itrc, i_dn) - trc_pth_flux(itrc, ipth)
               ENDDO
            ENDIF
         ENDDO
         ! Apply remote inflow
         DO i = 1, numucat_in
            DO itrc = 1, ntracers
               trc_bif_net(itrc, i) = trc_bif_net(itrc, i) - trc_bif_influx(itrc, i)
            ENDDO
         ENDDO
      ENDIF

      ! ===== Step 7: Single mass update (main + bif) =====
      DO i = 1, numucat_in
         DO itrc = 1, ntracers
            trc_mass(itrc, i) = trc_mass(itrc, i) &
               + (- trc_flux_out(itrc, i) + trc_flux_ups(itrc, i) &
                  - trc_bif_net(itrc, i)) * acctime
         ENDDO
      ENDDO

      ! ===== Step 8: Save bifurcation net flux for diagnostics =====
      DO i = 1, numucat_in
         DO itrc = 1, ntracers
            trc_bif_net_saved(itrc, i) = trc_bif_net(itrc, i)
         ENDDO
      ENDDO

      ! ===== Step 9: Update final concentration =====
      DO i = 1, numucat_in
         CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
         DO itrc = 1, ntracers
            trc_conc(itrc, i) = trc_mass(itrc, i) / volwater
         ENDDO
      ENDDO

      deallocate (trc_conc_next)
      deallocate (trc_flux_ups)
      deallocate (trc_out_mass)
      deallocate (rate_cell)
      deallocate (rate_next)
      deallocate (trc_bif_net)
      IF (allocated(trc_pth_flux  )) deallocate (trc_pth_flux  )
      IF (allocated(trc_bif_influx)) deallocate (trc_bif_influx)
      IF (allocated(trc_pth_1trc  )) deallocate (trc_pth_1trc  )
      IF (allocated(conc_dn_pth   )) deallocate (conc_dn_pth   )
      IF (allocated(rate_dn_pth   )) deallocate (rate_dn_pth   )

   END SUBROUTINE tracer_calc


   !-------------------------------------------------------------------------------------
   ! Helper: get water volume for a unit catchment cell, respecting reservoir/levee state.
   !-------------------------------------------------------------------------------------
   SUBROUTINE get_cell_volume (icell, wdsrf_cell, volresv_in, ucat2resv_in, volwater)

   USE MOD_Grid_RiverLakeNetwork, only: floodplain_curve, lake_type
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, volwater_ucat
   IMPLICIT NONE

   integer,  intent(in)  :: icell
   real(r8), intent(in)  :: wdsrf_cell
   real(r8), intent(in)  :: volresv_in(:)
   integer,  intent(in)  :: ucat2resv_in(:)
   real(r8), intent(out) :: volwater

      IF (lake_type(icell) == 2 .and. size(volresv_in) > 0) THEN
         volwater = volresv_in(ucat2resv_in(icell))
      ELSEIF (DEF_USE_LEVEE .and. has_levee(icell)) THEN
         volwater = volwater_ucat(icell)
      ELSE
         volwater = floodplain_curve(icell)%volume(wdsrf_cell)
      ENDIF
      volwater = max(volwater, 1.e-6_r8)

   END SUBROUTINE get_cell_volume


   !-------------------------------------------------------------------------------------
   ! Accumulate tracer diagnostics for history output
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_diag_accumulate (acctime)

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: acctime  ! Routing period [s]

   integer :: i, itrc

      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         DO itrc = 1, ntracers
            a_trc_conc  (itrc, i) = a_trc_conc  (itrc, i) + trc_conc(itrc, i) * acctime
            a_trc_out   (itrc, i) = a_trc_out   (itrc, i) + trc_flux_out(itrc, i) * acctime
            a_trc_bifout(itrc, i) = a_trc_bifout(itrc, i) + trc_bif_net_saved(itrc, i) * acctime
         ENDDO
      ENDDO

   END SUBROUTINE tracer_diag_accumulate


   !-------------------------------------------------------------------------------------
   ! Flush accumulated tracer diagnostics
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_flush_acc ()

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

      IF (numucat > 0) THEN
         IF (allocated(a_trc_conc  )) a_trc_conc   = 0._r8
         IF (allocated(a_trc_out   )) a_trc_out    = 0._r8
         IF (allocated(a_trc_bifout)) a_trc_bifout = 0._r8
      ENDIF

   END SUBROUTINE tracer_flush_acc


   !-------------------------------------------------------------------------------------
   ! Read tracer restart
   !-------------------------------------------------------------------------------------
   SUBROUTINE read_tracer_restart (file_restart)

   USE MOD_NetCDFSerial,          only: ncio_var_exist
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: itrc
   logical :: has_var
   character(len=64) :: varname
   real(r8), allocatable :: tmpvec(:)

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      DO itrc = 1, ntracers
         write(varname, '(A,A)') 'trc_mass_', trim(tracer_names(itrc))

         has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)

         IF (has_var) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) THEN
               trc_mass(itrc, :) = tmpvec(:)
            ENDIF
         ELSE
            ! Cold start: tracer variable not in restart file, initialise to zero
            IF (p_is_master) THEN
               write(*,'(A,A,A)') '  Tracer restart variable "', trim(varname), '" not found, cold start (zero).'
            ENDIF
            IF (p_is_worker .and. numucat > 0) THEN
               trc_mass(itrc, :) = 0._r8
            ENDIF
         ENDIF
      ENDDO

      deallocate (tmpvec)

   END SUBROUTINE read_tracer_restart


   !-------------------------------------------------------------------------------------
   ! Write tracer restart
   !-------------------------------------------------------------------------------------
   SUBROUTINE write_tracer_restart (file_restart)

   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: itrc
   character(len=64) :: varname
   real(r8), allocatable :: tmpvec(:)

      ! Guard: tracer module may not be initialised (e.g. mkinidata)
      IF (.not. allocated(tracer_names)) RETURN
      IF (.not. allocated(trc_mass))     RETURN

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      DO itrc = 1, ntracers
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = trc_mass(itrc, :)
         ENDIF

         write(varname, '(A,A)') 'trc_mass_', trim(tracer_names(itrc))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')
      ENDDO

      deallocate (tmpvec)

   END SUBROUTINE write_tracer_restart


   !-------------------------------------------------------------------------------------
   ! Write tracer history output
   !-------------------------------------------------------------------------------------
   SUBROUTINE write_tracer_history (file_hist_ucat, itime_in_file_ucat, acctime_hist)

   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address, &
      x_ucat, y_ucat, griducat
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist_ucat
   integer,  intent(in) :: itime_in_file_ucat
   real(r8), intent(in) :: acctime_hist  ! Total accumulated history time [s]

   integer :: itrc
   character(len=64) :: varname
   character(len=128) :: longname
   real(r8), allocatable :: tmpvec(:)

      IF (acctime_hist <= 0._r8) RETURN

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      DO itrc = 1, ntracers

         ! --- Tracer concentration ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = a_trc_conc(itrc, :) / acctime_hist
         ENDIF

         write(varname, '(A,A)') 'f_trc_conc_', trim(tracer_names(itrc))
         write(longname, '(A,A,A)') 'tracer concentration (', trim(tracer_names(itrc)), ')'

         CALL vector_gather_map2grid_and_write ( tmpvec, numucat,                        &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, trim(varname), 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            trim(longname), 'permil*m3/m3')

         ! --- Tracer outflux ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = a_trc_out(itrc, :) / acctime_hist
         ENDIF

         write(varname, '(A,A)') 'f_trc_flux_', trim(tracer_names(itrc))
         write(longname, '(A,A,A)') 'tracer outflux (', trim(tracer_names(itrc)), ')'

         CALL vector_gather_map2grid_and_write ( tmpvec, numucat,                        &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, trim(varname), 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            trim(longname), 'permil*m3/s')

         ! --- Tracer bifurcation net flux ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = a_trc_bifout(itrc, :) / acctime_hist
         ENDIF

         write(varname, '(A,A)') 'f_trc_bifout_', trim(tracer_names(itrc))
         write(longname, '(A,A,A)') 'tracer net bifurcation outflux (', trim(tracer_names(itrc)), ')'

         CALL vector_gather_map2grid_and_write ( tmpvec, numucat,                        &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, trim(varname), 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            trim(longname), 'permil*m3/s')

      ENDDO

      deallocate (tmpvec)

   END SUBROUTINE write_tracer_history


   !-------------------------------------------------------------------------------------
   ! Deallocate tracer module
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_final ()

   IMPLICIT NONE

      IF (allocated(tracer_names )) deallocate (tracer_names )
      IF (allocated(trc_mass     )) deallocate (trc_mass     )
      IF (allocated(trc_conc     )) deallocate (trc_conc     )
      IF (allocated(trc_flux_out )) deallocate (trc_flux_out )
      IF (allocated(acc_trc_inp       )) deallocate (acc_trc_inp       )
      IF (allocated(trc_bif_net_saved)) deallocate (trc_bif_net_saved)
      IF (allocated(a_trc_conc       )) deallocate (a_trc_conc       )
      IF (allocated(a_trc_out        )) deallocate (a_trc_out        )
      IF (allocated(a_trc_bifout     )) deallocate (a_trc_bifout     )

   END SUBROUTINE tracer_final

END MODULE MOD_Grid_RiverLakeTracer
#endif
