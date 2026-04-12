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
   real(r8), allocatable :: trc_mass  (:,:)   ! Heavy-water mass storage [kg] (ntracers, numucat)
   real(r8), allocatable :: trc_conc  (:,:)   ! Heavy-water concentration [kg/m3] (ntracers, numucat)

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
   PUBLIC :: check_tracer_state
   PUBLIC :: tracer_substep

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
   ! Accumulate heavy-water mass input from runoff.
   ! Called each land-model timestep (before routing accumulation threshold).
   !
   ! rnof_uc_vol(i):             runoff volume this step [m3]
   ! trc_rnof_ext(itrc, i):      (optional) heavy-water mass from land tracer [kg]
   !
   ! If trc_rnof_ext is present, use it directly (coupled to land tracer system).
   ! Otherwise, compute default heavy-water mass from ref_ratio and init_delta.
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_input_from_runoff (rnof_uc_vol, numucat_in, trc_rnof_ext)

      USE MOD_Tracer_Defs, only: tracers

      IMPLICIT NONE
      real(r8), intent(in) :: rnof_uc_vol(:)
      integer,  intent(in) :: numucat_in
      real(r8), intent(in), optional :: trc_rnof_ext(:,:)

      integer :: i, itrc
      real(r8) :: R_default

      DO i = 1, numucat_in
         IF (rnof_uc_vol(i) > 0._r8) THEN
            DO itrc = 1, ntracers
               IF (present(trc_rnof_ext)) THEN
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + trc_rnof_ext(itrc, i)
               ELSE
                  R_default = tracers(itrc)%ref_ratio * &
                     (1._r8 + tracers(itrc)%init_delta / 1000._r8)
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + R_default * rnof_uc_vol(i)
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
   real(r8), allocatable :: pth_wflux_avg(:)     ! Water flux direction per bif pathway [m3/s]

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
         allocate (pth_wflux_avg  (npthout_local))
         trc_pth_flux   = 0._r8
         trc_bif_influx = 0._r8
         pth_wflux_avg  = 0._r8
      ENDIF

      ! ===== Step 1: Add accumulated runoff tracer input to storage =====
      DO i = 1, numucat_in
         DO itrc = 1, ntracers
            trc_mass(itrc, i) = trc_mass(itrc, i) + acc_trc_inp(itrc, i)
         ENDDO
      ENDDO

      ! ===== Step 2: Concentration snapshot (frozen for all flux calcs) =====
      ! For cells with negligible water volume, set concentration to zero
      ! instead of dividing by the 1e-6 m³ floor, which would produce
      ! unphysical extreme values (e.g. -1e14 permil).
      DO i = 1, numucat_in
         CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
         IF (volwater > 1._r8) THEN
            DO itrc = 1, ntracers
               trc_conc(itrc, i) = trc_mass(itrc, i) / volwater
            ENDDO
         ELSE
            DO itrc = 1, ntracers
               trc_conc(itrc, i) = 0._r8
            ENDDO
         ENDIF
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
         ! Pre-compute water flux direction per pathway (tracer-independent).
         DO ipth = 1, npthout_local
            i_up = pth_upst_local(ipth)
            IF (i_up < 1 .or. i_up > numucat_in) CYCLE
            IF (prd_bifflw_time(ipth) <= 0._r8) CYCLE
            pth_wflux_avg(ipth) = sum(prd_bifflw_lev(:, ipth)) / prd_bifflw_time(ipth)
         ENDDO

         DO itrc = 1, ntracers
            CALL worker_push_data (push_bif_dn2pth, trc_conc(itrc,:), conc_dn_pth, fillvalue = 0._r8)

            DO ipth = 1, npthout_local
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat_in) CYCLE
               IF (prd_bifflw_time(ipth) <= 0._r8) CYCLE

               conc_up = trc_conc(itrc, i_up)
               conc_dn = conc_dn_pth(ipth)

               IF (pth_wflux_avg(ipth) >= 0._r8) THEN
                  trc_pth_flux(itrc, ipth) = conc_up * pth_wflux_avg(ipth)
               ELSE
                  trc_pth_flux(itrc, ipth) = conc_dn * pth_wflux_avg(ipth)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      ! ===== Step 5: Water-direction-based limiter =====
      !
      ! The donor (sender) of tracer is determined by WATER FLUX DIRECTION,
      ! not by tracer flux sign. This is critical for signed tracers like
      ! isotope delta values whose "mass" (delta × volume) can be negative.
      ! Using tracer flux sign would misidentify normal downstream flow
      ! carrying negative delta as reverse flow, breaking conservation.
      !
      ! P2STOOUT(i) = total |tracer| leaving cell i in all directions:
      !   + |flux_out(i)|  when hflux(i) >= 0   (normal: i is sender)
      !   + pushed contributions from upstream cells j where hflux(j) < 0
      !     (reverse: downstream cell = i is sender, mass flows toward j)
      !   + bif pathway contributions (using water flux direction)
      !
      ! D2RATE(i) = min(|trc_mass(i)| / P2STOOUT(i), 1)
      !
      DO itrc = 1, ntracers

         ! 5a: P2STOOUT — total outgoing |tracer| per cell, using water direction
         DO i = 1, numucat_in
            IF (hflux_avg(i) >= 0._r8) THEN
               trc_out_mass(i) = abs(trc_flux_out(itrc, i)) * acctime
            ELSE
               trc_out_mass(i) = 0._r8
            ENDIF
         ENDDO

         ! Reverse flow (hflux < 0): the downstream cell is the donor.
         ! Push |flux| to the downstream cell's P2STOOUT via push_ups2ucat.
         DO i = 1, numucat_in
            IF (hflux_avg(i) < 0._r8) THEN
               rate_cell(i) = abs(trc_flux_out(itrc, i)) * acctime
            ELSE
               rate_cell(i) = 0._r8
            ENDIF
         ENDDO
         CALL worker_push_data (push_ups2ucat, rate_cell, rate_next, fillvalue = 0._r8, mode = 'sum')
         DO i = 1, numucat_in
            trc_out_mass(i) = trc_out_mass(i) + rate_next(i)
         ENDDO

         ! Bifurcation pathways: use pth_wflux_avg to determine sender
         IF (do_bifurcation .and. npthlev_bif > 0) THEN
            DO ipth = 1, npthout_local
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat_in) CYCLE
               IF (pth_wflux_avg(ipth) >= 0._r8) THEN
                  trc_out_mass(i_up) = trc_out_mass(i_up) &
                     + abs(trc_pth_flux(itrc, ipth)) * acctime
               ELSE
                  i_dn = pth_down_local(ipth)
                  IF (i_dn > 0 .and. i_dn <= numucat_in) THEN
                     trc_out_mass(i_dn) = trc_out_mass(i_dn) &
                        + abs(trc_pth_flux(itrc, ipth)) * acctime
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         ! 5b: D2RATE per cell
         DO i = 1, numucat_in
            IF (trc_out_mass(i) > 1.e-30_r8) THEN
               rate_cell(i) = min(abs(trc_mass(itrc, i)) / trc_out_mass(i), 1._r8)
            ELSE
               rate_cell(i) = 1._r8
            ENDIF
         ENDDO

         ! 5c: Get downstream cell's rate (for reverse main flow)
         CALL worker_push_data (push_next2ucat, rate_cell, rate_next, fillvalue = 1._r8)

         ! 5d: Scale main-channel flux by sender's rate (water direction)
         DO i = 1, numucat_in
            IF (hflux_avg(i) >= 0._r8) THEN
               trc_flux_out(itrc, i) = trc_flux_out(itrc, i) * rate_cell(i)
            ELSE
               trc_flux_out(itrc, i) = trc_flux_out(itrc, i) * rate_next(i)
            ENDIF
         ENDDO

         ! 5e: Scale bif pathway flux by sender's rate (water direction)
         IF (do_bifurcation .and. npthlev_bif > 0) THEN
            CALL worker_push_data (push_bif_dn2pth, rate_cell, rate_dn_pth, fillvalue = 1._r8)

            DO ipth = 1, npthout_local
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat_in) CYCLE
               IF (pth_wflux_avg(ipth) >= 0._r8) THEN
                  trc_pth_flux(itrc, ipth) = trc_pth_flux(itrc, ipth) * rate_cell(i_up)
               ELSE
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
         IF (volwater > 1._r8) THEN
            DO itrc = 1, ntracers
               trc_conc(itrc, i) = trc_mass(itrc, i) / volwater
            ENDDO
         ELSE
            DO itrc = 1, ntracers
               trc_conc(itrc, i) = 0._r8
            ENDDO
         ENDIF
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
      IF (allocated(pth_wflux_avg )) deallocate (pth_wflux_avg )

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
   ! Per-sub-step tracer transport (called inside the DO WHILE routing loop).
   ! Uses instantaneous water fluxes, not time-averaged, so tracer and water
   ! advance in lockstep.  This avoids the "water left but tracer stayed"
   ! artefact of the old single-step-per-period approach.
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_substep (dt_all, irivsys, hflux_fc, wdsrf, ucatfilter, &
      volresv, ucat2resv, is_built_resv, &
      do_bif, bif_hflux_lev_in, npthout_local_in)

   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, &
      floodplain_curve, lake_type, push_ups2ucat, push_next2ucat, &
      npthlev_bif, pth_upst_local, pth_down_local, &
      push_bif_influx, push_bif_dn2pth
   USE MOD_Grid_RiverLakeLevee, only: has_levee, volwater_ucat
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   real(r8), intent(in) :: hflux_fc(:)
   real(r8), intent(in) :: wdsrf(:)
   logical,  intent(in) :: ucatfilter(:)
   real(r8), intent(in) :: volresv(:)
   integer,  intent(in) :: ucat2resv(:)
   logical,  intent(in) :: is_built_resv(:)
   logical,  intent(in) :: do_bif
   real(r8), intent(in) :: bif_hflux_lev_in(:,:)
   integer,  intent(in) :: npthout_local_in

   integer  :: i, itrc, ipth, i_up, i_dn
   real(r8) :: volwater, dt_i, pth_wflux, trc_pth_fl
   real(r8), allocatable :: conc_next(:)
   real(r8), allocatable :: flux_ups(:)
   real(r8), allocatable :: trc_flux(:)
   real(r8), allocatable :: conc_dn_pth(:)
   real(r8), allocatable :: trc_pth_1trc(:)
   real(r8), allocatable :: bif_recv(:)
   real(r8), allocatable :: bif_net(:)

      IF (numucat <= 0 .and. npthout_local_in <= 0) RETURN

      allocate (conc_next    (numucat))
      allocate (flux_ups     (numucat))
      allocate (trc_flux     (numucat))
      allocate (bif_net      (numucat))
      IF (do_bif .and. npthout_local_in > 0) THEN
         allocate (conc_dn_pth  (npthout_local_in))
         allocate (trc_pth_1trc (npthout_local_in))
         allocate (bif_recv     (numucat))
      ENDIF

      DO itrc = 1, ntracers

         ! --- 1. Concentration from pre-update state ---
         ! Use the TRUE cell volume (same as water routing, no 1e-6
         ! floor) to keep tracer exactly proportional to water. For
         ! cells with volume < 1e-6 m³ (effectively dry), zero BOTH
         ! concentration AND mass to prevent floating-point residual
         ! mass from producing extreme concentrations or NaN. This
         ! loses a negligible amount of tracer (<< 1e-10 of total)
         ! but maintains exact proportionality for all wet cells.
         DO i = 1, numucat
            IF (lake_type(i) == 2 .and. size(volresv) > 0) THEN
               volwater = volresv(ucat2resv(i))
            ELSEIF (DEF_USE_LEVEE .and. has_levee(i)) THEN
               volwater = volwater_ucat(i)
            ELSE
               volwater = floodplain_curve(i)%volume(wdsrf(i))
            ENDIF
            IF (volwater > 1.e-6_r8) THEN
               trc_conc(itrc, i) = trc_mass(itrc, i) / volwater
            ELSE
               trc_conc(itrc, i) = 0._r8
               trc_mass(itrc, i) = 0._r8
            ENDIF
         ENDDO

         ! --- 2. Get downstream concentration (main channel) ---
         CALL worker_push_data (push_next2ucat, trc_conc(itrc,:), conc_next, fillvalue = 0._r8)

         ! --- 3. Upwind main-channel tracer flux ---
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) THEN
               trc_flux(i) = 0._r8
               CYCLE
            ENDIF
            IF (hflux_fc(i) >= 0._r8) THEN
               trc_flux(i) = trc_conc(itrc, i) * hflux_fc(i)
            ELSE
               trc_flux(i) = conc_next(i) * hflux_fc(i)
            ENDIF
         ENDDO

         ! --- 4. Aggregate upstream tracer flux ---
         CALL worker_push_data (push_ups2ucat, trc_flux, flux_ups, fillvalue = 0._r8, mode = 'sum')

         ! --- 5. Bifurcation pathway tracer flux ---
         bif_net(:) = 0._r8
         IF (do_bif .and. npthout_local_in > 0) THEN
            CALL worker_push_data (push_bif_dn2pth, trc_conc(itrc,:), conc_dn_pth, fillvalue = 0._r8)

            trc_pth_1trc(:) = 0._r8
            DO ipth = 1, npthout_local_in
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               IF (.not. ucatfilter(i_up)) CYCLE

               pth_wflux = sum(bif_hflux_lev_in(:, ipth))
               IF (pth_wflux >= 0._r8) THEN
                  trc_pth_fl = trc_conc(itrc, i_up) * pth_wflux
               ELSE
                  trc_pth_fl = conc_dn_pth(ipth) * pth_wflux
               ENDIF

               bif_net(i_up) = bif_net(i_up) + trc_pth_fl
               i_dn = pth_down_local(ipth)
               IF (i_dn > 0 .and. i_dn <= numucat) THEN
                  bif_net(i_dn) = bif_net(i_dn) - trc_pth_fl
               ENDIF
               trc_pth_1trc(ipth) = trc_pth_fl
            ENDDO

            ! Remote bif scatter + de-duplicate
            CALL worker_push_data (push_bif_influx, trc_pth_1trc, bif_recv, &
               fillvalue = 0._r8, mode = 'sum')
            DO ipth = 1, npthout_local_in
               i_dn = pth_down_local(ipth)
               IF (i_dn > 0 .and. i_dn <= numucat) THEN
                  bif_recv(i_dn) = bif_recv(i_dn) - trc_pth_1trc(ipth)
               ENDIF
            ENDDO
            DO i = 1, numucat
               bif_net(i) = bif_net(i) - bif_recv(i)
            ENDDO
         ENDIF

         ! --- 6. Limiter ---
         ! No separate tracer limiter needed. The water routing already
         ! constrains dt via CFL + volume constraints such that
         ! hflux_fc*dt <= volwater. Since tracer flux = conc * hflux_fc
         ! (upwind), tracer_out*dt = conc * hflux_fc * dt <= |conc*vol|
         ! = |mass|. Tracer mass stays bounded by construction.
         !
         ! A mass-based limiter would reduce tracer outflow WITHOUT
         ! reducing water outflow, causing the cell to over-concentrate
         ! (tracer stays, water leaves). Removing it ensures tracer and
         ! water leave in the same proportion, preserving concentration.

         ! --- 7. Update tracer mass ---
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) CYCLE
            dt_i = dt_all(irivsys(i))
            IF (dt_i <= 0._r8) CYCLE

            trc_mass(itrc, i) = trc_mass(itrc, i) &
               + (- trc_flux(i) + flux_ups(i) - bif_net(i)) * dt_i
         ENDDO

         ! --- 8. Save flux for diagnostics ---
         DO i = 1, numucat
            trc_flux_out(itrc, i) = trc_flux(i)
            trc_bif_net_saved(itrc, i) = bif_net(i)
         ENDDO

      ENDDO  ! itrc

      ! --- 9. Final concentration (post-update) ---
      DO i = 1, numucat
         IF (lake_type(i) == 2 .and. size(volresv) > 0) THEN
            volwater = volresv(ucat2resv(i))
         ELSEIF (DEF_USE_LEVEE .and. has_levee(i)) THEN
            volwater = volwater_ucat(i)
         ELSE
            volwater = floodplain_curve(i)%volume(wdsrf(i))
         ENDIF
         IF (volwater > 1.e-6_r8) THEN
            DO itrc = 1, ntracers
               trc_conc(itrc, i) = trc_mass(itrc, i) / volwater
            ENDDO
         ELSE
            DO itrc = 1, ntracers
               trc_conc(itrc, i) = 0._r8
               trc_mass(itrc, i) = 0._r8
            ENDDO
         ENDIF
      ENDDO

      deallocate (conc_next, flux_ups, trc_flux, bif_net)
      IF (allocated(conc_dn_pth )) deallocate (conc_dn_pth )
      IF (allocated(trc_pth_1trc)) deallocate (trc_pth_1trc)
      IF (allocated(bif_recv    )) deallocate (bif_recv    )

   END SUBROUTINE tracer_substep


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

   integer :: itrc, has_flag
   logical :: has_var
   character(len=64) :: varname
   real(r8), allocatable :: tmpvec(:)

      ! vector_read_and_scatter contains mpi_barrier(p_comm_glb), so
      ! ALL ranks (master, workers, IO) must enter this routine.
      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      DO itrc = 1, ntracers
         write(varname, '(A,A)') 'trc_mass_', trim(tracer_names(itrc))

         ! Master-only file probe + broadcast to avoid concurrent opens
         IF (p_is_master) THEN
            has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
            has_flag = merge(1, 0, has_var)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
         has_var = (has_flag /= 0)

         IF (has_var) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) THEN
               trc_mass(itrc, :) = tmpvec(:)
            ENDIF
         ELSE
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

      ! Guard: tracer module may not be initialised (e.g. mkinidata).
      ! tracer_names is allocated on ALL ranks by tracer_init (before the
      ! p_is_worker guard), so this check is safe for master/IO.
      ! Do NOT check allocated(trc_mass) here — trc_mass is worker-only,
      ! but vector_gather_and_write has mpi_barrier(p_comm_glb) that
      ! requires ALL ranks to participate.
      IF (.not. allocated(tracer_names)) RETURN
      IF (p_is_worker .and. .not. allocated(trc_mass)) RETURN

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
   ! Diagnostic check: print min/max of tracer state variables
   !-------------------------------------------------------------------------------------
   SUBROUTINE check_tracer_state ()

   USE MOD_SPMD_Task
   USE MOD_RangeCheck
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   integer :: itrc
   character(len=64) :: label
   real(r8), allocatable :: tmp(:)

      ! Workers that never ran tracer_init have no data to check.
      ! Master and IO ranks don't allocate trc_mass but MUST still enter
      ! this routine: master participates in check_vector_data's recv/print
      ! side, and all ranks share the same call count to avoid orphaning
      ! worker-to-master mpi_send messages.
      IF (p_is_worker .and. (.not. allocated(trc_mass) .or. .not. allocated(trc_conc))) RETURN

      IF (p_is_master) THEN
         write(*,'(/,A)') 'Checking Tracer Variables ...'
      ENDIF

      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmp(numucat))
      ELSE
         allocate (tmp(0))
      ENDIF

      DO itrc = 1, ntracers
         IF (p_is_worker .and. numucat > 0) tmp = trc_mass(itrc,:)
         write(label,'(A,A,A)') 'trc_mass_', trim(tracer_names(itrc)), ' [mass]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) tmp = trc_conc(itrc,:)
         write(label,'(A,A,A)') 'trc_conc_', trim(tracer_names(itrc)), ' [m/m3]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) tmp = trc_flux_out(itrc,:)
         write(label,'(A,A,A)') 'trc_flux_', trim(tracer_names(itrc)), ' [m/s]'
         CALL check_vector_data (label, tmp)
      ENDDO

      deallocate (tmp)

   END SUBROUTINE check_tracer_state


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
