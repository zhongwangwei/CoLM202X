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
   ! L1: share the single authoritative `ntracers` from MOD_Tracer_Defs to
   ! avoid the land/river tracer modules holding two independent copies
   ! that could silently diverge on re-init.
   USE MOD_Tracer_Defs, only: ntracers
   IMPLICIT NONE

   !-------------------------------------------------------------------------------------
   ! Module variables
   !-------------------------------------------------------------------------------------
   character(len=32), allocatable :: tracer_names(:)  ! Tracer names

   ! State variables (prognostic)
   ! trc_mass is built from R_default * water_volume, so its units are
   ! [R*m3] (isotope ratio times cell water volume), NOT kg of heavy water.
   ! trc_conc = trc_mass / volwater is therefore a dimensionless ratio [R].
   real(r8), allocatable :: trc_mass  (:,:)   ! Tracer amount [R*m3] (ntracers, numucat)
   real(r8), allocatable :: trc_conc  (:,:)   ! Tracer ratio  [R]    (ntracers, numucat)

   ! H1 fix: protected-side (behind-levee) tracer pool paired with
   ! `levsto` in MOD_Grid_RiverLakeLevee. Levee repartition moves water
   ! between visible-side (volwater_ucat, driving trc_mass) and
   ! protected-side storage; without its own tracer pool the mass would
   ! stay pinned to the visible side while water moves, systematically
   ! inflating visible-side concentration (and orphaning mass when
   ! visible-side goes dry — previously only tracked via DBG_LEVTRC).
   real(r8), allocatable :: trc_levsto(:,:)   ! Protected-side tracer amount [R*m3] (ntracers, numucat)

   ! Flux variables (diagnostic, per routing period)
   real(r8), allocatable :: trc_flux_out  (:,:) ! Tracer outflux [R*m3/s] (ntracers, numucat)
   real(r8), allocatable :: trc_dry_drain(:,:) ! Internal sink from dry-cell forced drain [R*m3]

   ! Input: tracer amount flux from runoff
   real(r8), allocatable :: acc_trc_inp  (:,:) ! Accumulated tracer input [R*m3] (ntracers, numucat)
   real(r8), allocatable :: trc_inp_buf  (:,:) ! Buffered runoff tracer awaiting release [R*m3]
   real(r8), allocatable :: acc_rnof_ref (:)   ! Accumulated runoff water volume paired with acc_trc_inp [m3]

   ! Bifurcation net flux (saved from last tracer_substep for diagnostics)
   real(r8), allocatable :: trc_bif_net_saved (:,:) ! Per-cell net bif flux [mass/s] (ntracers, numucat)

   ! Near-dry volume floor used only for transport concentration
   ! stabilization; prognostic tracer mass always remains single-pool.
   real(r8), parameter :: trc_v_dry_off = 1.e-6_r8
   real(r8), parameter :: inp_cap_factor = 3._r8

   ! History accumulators
   real(r8), allocatable :: a_trc_conc   (:,:) ! Accumulated tracer conc [mass/m3 * s] (ntracers, numucat)
   real(r8), allocatable :: a_trc_out    (:,:) ! Accumulated tracer outflux [mass/s * s] (ntracers, numucat)
   real(r8), allocatable :: a_trc_bifout (:,:) ! Accumulated tracer bif net flux [mass/s * s] (ntracers, numucat)

   PUBLIC :: tracer_init
   PUBLIC :: tracer_init_from_water
   PUBLIC :: tracer_input_from_runoff
   PUBLIC :: tracer_flush_acc
   PUBLIC :: read_tracer_restart
   PUBLIC :: write_tracer_restart
   PUBLIC :: write_tracer_history
   PUBLIC :: tracer_final
   PUBLIC :: check_tracer_state
   PUBLIC :: tracer_substep
   PUBLIC :: tracer_refresh_state
   PUBLIC :: tracer_diag_accumulate_substep
   ! H3: exposed so the inland-depression overflow fix in
   ! MOD_Grid_RiverLakeFlow can use the same reservoir/levee/floodplain
   ! volume selection as tracer_refresh_state instead of hard-wiring
   ! topo_rivstomax (which disagrees on leveed cells).
   PUBLIC :: get_cell_volume
   PUBLIC :: trc_inp_buf
   PUBLIC :: acc_rnof_ref
   PUBLIC :: trc_levsto
   PUBLIC :: trc_dry_drain
   PUBLIC :: levee_tracer_repartition

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

      ! H4: empty string yields zero tokens, matching land-side parse_csv
      ! (MOD_Tracer_Defs.F90:158 where `IF (j <= slen)` skips the
      ! trailing-segment add). Without this the river side would invent a
      ! single `'unnamed'` tracer (from sanitize_ncname of a zero-length
      ! substring), which then diverged from the land-side count.
      IF (len_trim(str) == 0) THEN
         IF (p_is_master) THEN
            write(*,'(A,I0,A)') ' WARNING: DEF_TRACER_NAMES is empty but DEF_TRACER_NUM = ', &
               n, '. Filling all names with tracer_N.'
         ENDIF
         DO i = 1, n
            write(names(i), '(A,I0)') 'tracer_', i
         ENDDO
         RETURN
      ENDIF

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

      ! ntracers now populated by tracer_defs_init() from CoLM.F90 before
      ! any HYDRO tracer_init runs — no local assignment needed.
      allocate (tracer_names(ntracers))
      CALL parse_tracer_names(DEF_TRACER_NAMES, tracer_names, ntracers)

      ! Ensure unique names: append index if duplicates found after sanitization.
      ! Budget the base name so `base // '_' // i` does not re-truncate to a
      ! duplicate at len=32 (mirrors the dedup guard in
      ! MOD_Tracer_Defs:tracer_defs_init).
      BLOCK
         integer :: sfx_len, max_base
         character(len=32) :: sfx, base
         DO i = 1, ntracers
            DO j = 1, i-1
               IF (trim(tracer_names(i)) == trim(tracer_names(j))) THEN
                  write(sfx, '(A,I0)') '_', i
                  sfx_len = len_trim(sfx)
                  max_base = max(1, 32 - sfx_len)
                  base = tracer_names(i)
                  IF (len_trim(base) > max_base) base = base(1:max_base)
                  write(tracer_names(i), '(A,A)') trim(base), trim(sfx)
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      END BLOCK

      IF (p_is_worker) THEN
         ! Allocate on ALL workers (zero-length if numucat=0) for MPI safety.
         allocate (trc_mass     (ntracers, numucat))
         allocate (trc_conc     (ntracers, numucat))
         allocate (trc_flux_out (ntracers, numucat))
         allocate (trc_dry_drain(ntracers, numucat))
         allocate (acc_trc_inp        (ntracers, numucat))
         allocate (trc_inp_buf        (ntracers, numucat))
         allocate (acc_rnof_ref       (numucat))
         allocate (trc_bif_net_saved  (ntracers, numucat))
         allocate (a_trc_conc         (ntracers, numucat))
         allocate (a_trc_out          (ntracers, numucat))
         allocate (a_trc_bifout       (ntracers, numucat))
         allocate (trc_levsto         (ntracers, numucat))

         trc_mass     = 0._r8
         trc_conc     = 0._r8
         trc_flux_out  = 0._r8
         trc_dry_drain = 0._r8
         acc_trc_inp  = 0._r8
         trc_inp_buf  = 0._r8
         acc_rnof_ref = 0._r8
         trc_bif_net_saved = 0._r8
         a_trc_conc   = 0._r8
         a_trc_out    = 0._r8
         a_trc_bifout = 0._r8
         trc_levsto   = 0._r8
      ENDIF

      IF (p_is_master) THEN
         write(*,'(A,I4,A)') ' Tracer module initialised with ', ntracers, ' tracers:'
         write(*,'(A,*(A,:,", "))') '   Names: ', (trim(tracer_names(i)), i=1, ntracers)
      ENDIF

   END SUBROUTINE tracer_init


   !-------------------------------------------------------------------------------------
   ! Cold-start river tracer mass from the current water volume.
   ! Mirrors the land-side tracer_init_from_water logic: when no tracer
   ! restart variable is present in the gridriver restart file but the
   ! hydrology restart provides a non-zero wdsrf_ucat/volresv, seed
   ! trc_mass = volwater * R_init. Without this, any hot-start of water
   ! starts with tracer=0 and bleeds spurious dilution into downstream
   ! concentrations until runoff "re-colours" the reservoir.
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_init_from_water (wdsrf, volresv_in, ucat2resv_in, missing_mask)

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, levsto
   USE MOD_Tracer_Defs,           only: tracers, delta_to_R
   IMPLICIT NONE

   real(r8), intent(in) :: wdsrf(:)
   real(r8), intent(in) :: volresv_in(:)
   integer,  intent(in) :: ucat2resv_in(:)
   ! Optional per-tracer selector. When present, only tracers with
   ! missing_mask(itrc) == .true. are re-seeded from the water volume;
   ! this avoids overwriting tracers that loaded successfully from a
   ! partially-complete restart file (e.g. user renamed one tracer or
   ! added a new tracer to DEF_TRACER_NAMES without re-writing the file).
   ! Without the mask every tracer is cold-started (equivalent to
   ! passing an all-.true. mask).
   logical,  optional, intent(in) :: missing_mask(:)

   integer  :: i, itrc
   real(r8) :: volwater, R_init
   logical  :: do_init
   integer  :: n_init

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN
      IF (.not. allocated(trc_mass)) RETURN

      n_init = 0
      DO itrc = 1, ntracers
         do_init = .true.
         IF (present(missing_mask)) THEN
            IF (itrc <= size(missing_mask)) do_init = missing_mask(itrc)
         ENDIF
         IF (.not. do_init) CYCLE
         R_init = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
         DO i = 1, numucat
            CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
            trc_mass(itrc, i) = max(volwater, 0._r8) * R_init
            ! H1 fix: seed protected-side pool from levsto so a cold
            ! start with water already behind the levee doesn't produce
            ! a phantom conservation jump on the first repartition.
            IF (DEF_USE_LEVEE .and. allocated(trc_levsto)) THEN
               IF (has_levee(i)) THEN
                  trc_levsto(itrc, i) = max(levsto(i), 0._r8) * R_init
               ELSE
                  trc_levsto(itrc, i) = 0._r8
               ENDIF
            ENDIF
         ENDDO
         ! Recompute concentration so diagnostics consistent before first step.
         DO i = 1, numucat
            CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
            CALL update_tracer_concentration(itrc, i, volwater)
         ENDDO
         n_init = n_init + 1
      ENDDO

      IF (p_is_master .and. n_init > 0) THEN
         write(*,'(A,I0,A)') '  River tracer cold-started from water volume for ', n_init, ' tracer(s).'
      ENDIF

   END SUBROUTINE tracer_init_from_water


   !-------------------------------------------------------------------------------------
   ! Accumulate heavy-water mass input from runoff.
   ! Called each land-model timestep (before routing accumulation threshold).
   !
   ! rnof_uc_depth(i):           runoff water column [m] (runoff flux * dt,
   !                             area-weighted per ucat cell — NOT a volume).
   ! trc_rnof_ext(itrc, i):      (optional) tracer amount from land tracer
   !                             in matching [R*m] units (depth-based, not mass).
   !
   ! If trc_rnof_ext is present, use it directly (coupled to land tracer system).
   ! Otherwise, compute default heavy-water mass from ref_ratio and init_delta.
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_input_from_runoff (rnof_uc_depth, numucat_in, trc_rnof_ext)

      ! Re-uses the per-tracer ref_ratio / init_delta cached by
      ! tracer_defs_init (MOD_Tracer_Defs.F90:82-84). Previously this
      ! routine re-parsed the CSV namelist strings via a local helper on
      ! every land timestep — wasteful and a divergence hazard (two
      ! parse code paths for the same inputs).
      USE MOD_Tracer_Defs, only: tracers, delta_to_R
      IMPLICIT NONE
      real(r8), intent(in) :: rnof_uc_depth(:)
      integer,  intent(in) :: numucat_in
      real(r8), intent(in), optional :: trc_rnof_ext(:,:)

      integer :: i, itrc
      real(r8), allocatable :: R_default(:)

      ! Cache per-tracer R_default once; absent trc_rnof_ext it is used
      ! for every cell and was previously recomputed inside the inner
      ! loop (delta_to_R is pure but costs trig/div per call).
      IF (.not. present(trc_rnof_ext)) THEN
         allocate(R_default(ntracers))
         DO itrc = 1, ntracers
            R_default(itrc) = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
         ENDDO
      ENDIF

      DO i = 1, numucat_in
         IF (rnof_uc_depth(i) > 0._r8) THEN
            acc_rnof_ref(i) = acc_rnof_ref(i) + rnof_uc_depth(i)
            DO itrc = 1, ntracers
               IF (present(trc_rnof_ext)) THEN
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + trc_rnof_ext(itrc, i)
               ELSE
                  acc_trc_inp(itrc, i) = acc_trc_inp(itrc, i) + R_default(itrc) * rnof_uc_depth(i)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (allocated(R_default)) deallocate(R_default)

   END SUBROUTINE tracer_input_from_runoff


   !-------------------------------------------------------------------------------------
   ! Helper: get water volume for a unit catchment cell, respecting reservoir/levee state.
   !-------------------------------------------------------------------------------------
   SUBROUTINE get_cell_volume (icell, wdsrf_cell, volresv_in, ucat2resv_in, volwater)

   USE MOD_Grid_RiverLakeNetwork, only: floodplain_curve, lake_type
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, volwater_ucat, &
      volwater_ucat_valid, levee_hgt, levee_topsto, levsto
   USE MOD_Vars_Global,           only: spval
   IMPLICIT NONE

   integer,  intent(in)  :: icell
   real(r8), intent(in)  :: wdsrf_cell
   real(r8), intent(in)  :: volresv_in(:)
   integer,  intent(in)  :: ucat2resv_in(:)
   real(r8), intent(out) :: volwater

      IF (lake_type(icell) == 2 .and. size(volresv_in) > 0 .and. size(ucat2resv_in) > 0) THEN
         IF (ucat2resv_in(icell) > 0 .and. ucat2resv_in(icell) <= size(volresv_in) &
            .and. volresv_in(ucat2resv_in(icell)) /= spval) THEN
            volwater = volresv_in(ucat2resv_in(icell))
         ELSE
            volwater = floodplain_curve(icell)%volume(wdsrf_cell)
         ENDIF
      ELSEIF (DEF_USE_LEVEE .and. has_levee(icell)) THEN
         ! Mirror MOD_Grid_RiverLakeFlow.F90:481 fallback chain so the
         ! tracer path does not read volwater_ucat=0 when the cached
         ! per-cell volume is still invalid (fresh levee_init without
         ! a levee restart, or old restart that lacked volwater_ucat).
         ! Without this, tracer_init_from_water / refresh reads 0 for
         ! levee cells that physically hold water, yielding trc_mass=0
         ! and diluting downstream concentrations until the first
         ! routing step flips volwater_ucat_valid to .true.
         IF (volwater_ucat_valid) THEN
            volwater = volwater_ucat(icell)
         ELSEIF (levsto(icell) > 0._r8) THEN
            IF (wdsrf_cell <= floodplain_curve(icell)%rivhgt + levee_hgt(icell) + 1.e-6_r8) THEN
               volwater = levee_topsto(icell)
            ELSE
               volwater = floodplain_curve(icell)%volume(wdsrf_cell) - levsto(icell)
            ENDIF
         ELSE
            volwater = floodplain_curve(icell)%volume(wdsrf_cell)
         ENDIF
         volwater = max(volwater, 0._r8)
      ELSE
         volwater = floodplain_curve(icell)%volume(wdsrf_cell)
      ENDIF

   END SUBROUTINE get_cell_volume


   !-------------------------------------------------------------------------------------
   ! H1 fix: repartition tracer mass between visible-side (`trc_mass`)
   ! and protected-side (`trc_levsto`) pools in lockstep with the water
   ! redistribution that `levee_fldstg` just performed.
   !
   ! Inputs are the visible-side and protected-side water volumes BEFORE
   ! and AFTER the levee flood-staging call. Two conventions:
   !   d_lev > 0 : water moved visible -> protected. Tracer transferred
   !               at the visible-side's current ratio.
   !   d_lev < 0 : water moved protected -> visible. Tracer transferred
   !               at the protected-side's current ratio.
   ! This preserves asymmetric concentrations across the levee (unlike
   ! a uniform pool-ratio split which homogenises the two compartments).
   !-------------------------------------------------------------------------------------
   SUBROUTINE levee_tracer_repartition (icell, vis_vol_bef, levsto_bef, vis_vol_aft, levsto_aft)

   USE MOD_Tracer_Defs, only: trc_tiny, tracers, delta_to_R
   IMPLICIT NONE
   integer,  intent(in) :: icell
   real(r8), intent(in) :: vis_vol_bef, levsto_bef
   real(r8), intent(in) :: vis_vol_aft, levsto_aft

   integer  :: itrc
   real(r8) :: d_lev, ratio, trc_move

      IF (.not. allocated(trc_mass) .or. .not. allocated(trc_levsto)) RETURN
      IF (icell < 1 .or. icell > size(trc_mass, 2)) RETURN

      d_lev = levsto_aft - levsto_bef

      IF (abs(d_lev) < trc_tiny) RETURN

      DO itrc = 1, ntracers
         IF (d_lev > 0._r8) THEN
            ! Visible -> protected
            IF (vis_vol_bef > trc_tiny) THEN
               ratio = trc_mass(itrc, icell) / vis_vol_bef
            ELSE
               ratio = 0._r8
            ENDIF
            trc_move = d_lev * ratio
            ! Don't overdraw the donor; for signed-negative mass, the
            ! ratio sign carries through so both sides move in tandem.
            IF (trc_move > 0._r8) THEN
               trc_move = min(trc_move, max(trc_mass(itrc, icell), 0._r8))
            ELSE
               trc_move = max(trc_move, min(trc_mass(itrc, icell), 0._r8))
            ENDIF
            trc_mass  (itrc, icell) = trc_mass  (itrc, icell) - trc_move
            trc_levsto(itrc, icell) = trc_levsto(itrc, icell) + trc_move
         ELSE
            ! Protected -> visible
            IF (levsto_bef > trc_tiny) THEN
               ratio = trc_levsto(itrc, icell) / levsto_bef
            ELSE
               ratio = 0._r8
            ENDIF
            trc_move = abs(d_lev) * ratio
            IF (trc_move > 0._r8) THEN
               trc_move = min(trc_move, max(trc_levsto(itrc, icell), 0._r8))
            ELSE
               trc_move = max(trc_move, min(trc_levsto(itrc, icell), 0._r8))
            ENDIF
            trc_levsto(itrc, icell) = trc_levsto(itrc, icell) - trc_move
            trc_mass  (itrc, icell) = trc_mass  (itrc, icell) + trc_move
         ENDIF
      ENDDO

   END SUBROUTINE levee_tracer_repartition


   !-------------------------------------------------------------------------------------
   ! Refresh concentration from single-pool mass and current water volume.
   !-------------------------------------------------------------------------------------
   SUBROUTINE update_tracer_concentration (itrc, icell, volwater)

   IMPLICIT NONE

   integer,  intent(in) :: itrc, icell
   real(r8), intent(in) :: volwater

      ! Intentionally does NOT clamp trc_mass here, even for tiny
      ! negatives (e.g. -1e-15 R*m³). Clamping state silently would
      ! swallow any future transport / flux-limiter bug that drifts
      ! trc_mass negative progressively — each substep would re-zero
      ! the residual below an absolute threshold and hide the drift.
      ! Instead, display-side reports (check_tracer_state) dust-clamp
      ! for log cleanliness, and report a WARNING when magnitudes
      ! exceed the FP-noise floor. See tracer_mass_fp_dust in
      ! check_tracer_state below for the cosmetic threshold.
      IF (volwater <= trc_v_dry_off) THEN
         trc_conc(itrc, icell) = 0._r8
      ELSE
         trc_conc(itrc, icell) = trc_mass(itrc, icell) / volwater
      ENDIF

   END SUBROUTINE update_tracer_concentration


   !-------------------------------------------------------------------------------------
   ! Refresh tracer concentration from the final hydrologic state.
   ! This is called after the routing loop has updated water volumes so that
   ! diagnostics/history read a concentration consistent with the final step state.
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_refresh_state (wdsrf, volresv_in, ucat2resv_in)

   USE MOD_Grid_RiverLakeNetwork, only: numucat, floodplain_curve, lake_type
   USE MOD_Grid_RiverLakeLevee,   only: has_levee, volwater_ucat
   IMPLICIT NONE

   real(r8), intent(in) :: wdsrf(:)
   real(r8), intent(in) :: volresv_in(:)
   integer,  intent(in) :: ucat2resv_in(:)

   integer :: i, itrc
   real(r8) :: volwater

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         CALL get_cell_volume(i, wdsrf(i), volresv_in, ucat2resv_in, volwater)
         DO itrc = 1, ntracers
            CALL update_tracer_concentration(itrc, i, volwater)
         ENDDO
      ENDDO

   END SUBROUTINE tracer_refresh_state


   !-------------------------------------------------------------------------------------
   ! Accumulate tracer history diagnostics for one routing substep using the final
   ! post-update water state and any post-transport flux corrections already applied
   ! by the flow solver (for example inland-depression overflow adjustments).
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_diag_accumulate_substep (dt_all, irivsys, ucatfilter, wdsrf, volresv_in, ucat2resv_in)

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: dt_all(:)
   integer,  intent(in) :: irivsys(:)
   logical,  intent(in) :: ucatfilter(:)
   real(r8), intent(in) :: wdsrf(:)
   real(r8), intent(in) :: volresv_in(:)
   integer,  intent(in) :: ucat2resv_in(:)

   integer :: i, itrc
   real(r8) :: dt_i

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      CALL tracer_refresh_state(wdsrf, volresv_in, ucat2resv_in)

      DO i = 1, numucat
         IF (.not. ucatfilter(i)) CYCLE
         IF (irivsys(i) <= 0 .or. irivsys(i) > size(dt_all)) CYCLE
         dt_i = dt_all(irivsys(i))
         IF (dt_i <= 0._r8) CYCLE

         DO itrc = 1, ntracers
            a_trc_conc  (itrc, i) = a_trc_conc  (itrc, i) + trc_conc(itrc, i) * dt_i
            a_trc_out   (itrc, i) = a_trc_out   (itrc, i) + trc_flux_out(itrc, i) * dt_i
            a_trc_bifout(itrc, i) = a_trc_bifout(itrc, i) + trc_bif_net_saved(itrc, i) * dt_i
         ENDDO
      ENDDO

   END SUBROUTINE tracer_diag_accumulate_substep


   !-------------------------------------------------------------------------------------
   ! Per-sub-step tracer transport (called inside the DO WHILE routing loop).
   ! Uses instantaneous water fluxes, not time-averaged, so tracer and water
   ! advance in lockstep.  This avoids the "water left but tracer stayed"
   ! artefact of the old single-step-per-period approach.
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_substep (dt_ref, dt_all, irivsys, hflux_fc, wdsrf, ucatfilter, &
      volresv, ucat2resv, is_built_resv, &
      do_bif, bif_hflux_lev_in, npthout_local_in)

   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, &
      floodplain_curve, lake_type, push_ups2ucat, push_next2ucat, &
      npthlev_bif, pth_upst_local, pth_down_local, &
      push_bif_influx, push_bif_dn2pth
   USE MOD_Grid_RiverLakeLevee, only: has_levee, volwater_ucat
   USE MOD_WorkerPushData
   USE MOD_Tracer_Defs, only: trc_tiny, tracers, delta_to_R
   IMPLICIT NONE

   real(r8), intent(in) :: dt_ref
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
   real(r8) :: volwater, dt_i, pth_wflux, trc_pth_fl, volflux, inj_frac
   real(r8) :: trc_inj_tau, m_cap, m_room, m_tau, release, R_cap
   real(r8), allocatable :: conc_next(:)
   real(r8), allocatable :: flux_ups(:)
   real(r8), allocatable :: trc_flux(:)
   real(r8), allocatable :: trc_conc_flux(:)
   real(r8), allocatable :: conc_dn_pth(:)
   real(r8), allocatable :: trc_pth_1trc(:)
   real(r8), allocatable :: bif_recv(:)
   real(r8), allocatable :: bif_net(:)
   ! Mass-based flux limiter variables
   real(r8), allocatable :: trc_out_mass(:)  ! Total outgoing |tracer| per cell
   real(r8), allocatable :: rate_cell(:)     ! Limiter rate per cell
   real(r8), allocatable :: rate_next(:)     ! Downstream cell rate (via push)
   real(r8), allocatable :: rate_dn_pth(:)   ! Downstream rate per bif pathway

      IF (numucat <= 0 .and. npthout_local_in <= 0) RETURN

      allocate (conc_next    (numucat))
      allocate (flux_ups     (numucat))
      allocate (trc_flux     (numucat))
      allocate (trc_conc_flux(numucat))
      allocate (bif_net      (numucat))
      allocate (trc_out_mass (numucat))
      allocate (rate_cell    (numucat))
      allocate (rate_next    (numucat))
      IF (do_bif .and. npthout_local_in > 0) THEN
         allocate (conc_dn_pth  (npthout_local_in))
         allocate (trc_pth_1trc (npthout_local_in))
         allocate (bif_recv     (numucat))
         allocate (rate_dn_pth  (npthout_local_in))
      ENDIF

      trc_inj_tau = dt_ref

      DO itrc = 1, ntracers

         ! --- 1. Inject this routing-period runoff tracer proportionally in time ---
         DO i = 1, numucat
            dt_i = 0._r8
            IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
            IF (dt_i <= 0._r8 .or. dt_ref <= 0._r8) CYCLE
            inj_frac = dt_i / dt_ref
            trc_inp_buf(itrc, i) = trc_inp_buf(itrc, i) + acc_trc_inp(itrc, i) * inj_frac
         ENDDO

         ! --- 2. Concentration from pre-update single-pool state ---
         DO i = 1, numucat
            CALL get_cell_volume(i, wdsrf(i), volresv, ucat2resv, volwater)
            dt_i = 0._r8
            IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
            ! M1 fix: do not release buffer into a dry cell. Without water,
            ! released mass has nowhere to host a concentration
            ! (update_tracer_concentration below would zero trc_conc and
            ! leave trc_mass as hidden inventory). The dry-drain pass in
            ! Section 10 below folds any orphan trc_mass + trc_inp_buf into
            ! trc_flux_out as an exit flux.
            IF (volwater > trc_v_dry_off) THEN
               ! m_cap caps release only while this period has fresh input; once
               ! acc_trc_inp has been consumed, any residual buffer must still be
               ! able to drain (otherwise mass freezes on cells that no longer
               ! receive runoff — acc_trc_inp/acc_rnof_ref are reset every period).
               IF (acc_trc_inp(itrc, i) > trc_tiny) THEN
                  m_cap = inp_cap_factor * (acc_trc_inp(itrc, i) / max(acc_rnof_ref(i), 1.e-30_r8)) * volwater
                  m_room = max(0._r8, m_cap - trc_mass(itrc, i))
               ELSE
                  ! Residual buffer from a previous routing period must still
                  ! respect the current host-cell water volume. Releasing the
                  ! whole buffer into a tiny-but-wet cell creates nonphysical
                  ! concentration spikes (DBG_MAXCONC) because acc_trc_inp /
                  ! acc_rnof_ref have already been reset. Cap the fallback
                  ! concentration by the larger of the current cell ratio and
                  ! the tracer's Phase-1 baseline R_init.
                  R_cap = max(delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio), &
                     max(trc_mass(itrc, i), 0._r8) / max(volwater, trc_v_dry_off))
                  m_cap = inp_cap_factor * R_cap * volwater
                  m_room = max(0._r8, m_cap - trc_mass(itrc, i))
               ENDIF
               m_tau = trc_inp_buf(itrc, i) * dt_i / max(trc_inj_tau, dt_i)
               ! H2: clamp release to non-negative. In the normal path
               ! trc_inp_buf is accumulated from non-negative rnof input and
               ! only drained by previous releases, so it never goes below zero;
               ! this guard makes sure an upstream corruption (partial restart,
               ! FP drift near trc_tiny) can't propagate negative mass into
               ! trc_mass via a negative release.
               release = max(0._r8, min(trc_inp_buf(itrc, i), min(m_room, m_tau)))
               trc_inp_buf(itrc, i) = trc_inp_buf(itrc, i) - release
               trc_mass(itrc, i) = trc_mass(itrc, i) + release
            ENDIF
            CALL update_tracer_concentration(itrc, i, volwater)
            volflux = max(volwater, max(hflux_fc(i), 0._r8) * dt_i)
            trc_conc_flux(i) = trc_mass(itrc, i) / max(volflux, trc_v_dry_off)
         ENDDO

         ! --- 3. Get downstream concentration (main channel) ---
         CALL worker_push_data (push_next2ucat, trc_conc_flux, conc_next, fillvalue = 0._r8)

         ! --- 4. Upwind main-channel tracer flux ---
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) THEN
               trc_flux(i) = 0._r8
               CYCLE
            ENDIF
            IF (hflux_fc(i) >= 0._r8) THEN
               trc_flux(i) = trc_conc_flux(i) * hflux_fc(i)
            ELSE
               trc_flux(i) = conc_next(i) * hflux_fc(i)
            ENDIF
         ENDDO

         ! --- 5. Aggregate upstream tracer flux ---
         CALL worker_push_data (push_ups2ucat, trc_flux, flux_ups, fillvalue = 0._r8, mode = 'sum')

         ! --- 6. Bifurcation pathway tracer flux ---
         bif_net(:) = 0._r8
         IF (do_bif .and. npthout_local_in > 0) THEN
            CALL worker_push_data (push_bif_dn2pth, trc_conc_flux, conc_dn_pth, fillvalue = 0._r8)

            trc_pth_1trc(:) = 0._r8
            DO ipth = 1, npthout_local_in
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               IF (.not. ucatfilter(i_up)) CYCLE

               pth_wflux = sum(bif_hflux_lev_in(:, ipth))
               IF (pth_wflux >= 0._r8) THEN
                  trc_pth_fl = trc_conc_flux(i_up) * pth_wflux
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

         ! --- 7. Mass-based flux limiter ---
         ! CFL constrains NET water flux (outflow - inflow), but tracer
         ! outflow = conc × GROSS outflow.  When upstream inflow brings
         ! water with lower concentration (e.g. cold start), total tracer
         ! export can exceed the cell's mass.  This limiter prevents
         ! negative mass while preserving upwind proportionality.

         ! 6a: P2STOOUT — total outgoing |tracer| per cell (water-direction)
         ! H1: guard dt_all(irivsys(i)) the same way as section 0/8 so a
         ! corrupt/uninitialised irivsys doesn't segfault here even if the
         ! water-flow main loop happens to skip it via ucatfilter. See
         ! MOD_Grid_RiverLakeFlow.F90:585 for the primary defense.
         DO i = 1, numucat
            IF (hflux_fc(i) >= 0._r8 .and. &
                irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) THEN
               trc_out_mass(i) = abs(trc_flux(i)) * dt_all(irivsys(i))
            ELSE
               trc_out_mass(i) = 0._r8
            ENDIF
         ENDDO
         ! Reverse flow: downstream cell is donor → push to its P2STOOUT
         DO i = 1, numucat
            IF (hflux_fc(i) < 0._r8 .and. &
                irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) THEN
               rate_cell(i) = abs(trc_flux(i)) * dt_all(irivsys(i))
            ELSE
               rate_cell(i) = 0._r8
            ENDIF
         ENDDO
         CALL worker_push_data (push_ups2ucat, rate_cell, rate_next, fillvalue = 0._r8, mode = 'sum')
         DO i = 1, numucat
            trc_out_mass(i) = trc_out_mass(i) + rate_next(i)
         ENDDO
         ! Bif pathways: add sender-side contribution
         IF (do_bif .and. npthout_local_in > 0) THEN
            DO ipth = 1, npthout_local_in
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               pth_wflux = sum(bif_hflux_lev_in(:, ipth))
               IF (pth_wflux >= 0._r8) THEN
                  IF (irivsys(i_up) > 0 .and. irivsys(i_up) <= size(dt_all)) THEN
                     trc_out_mass(i_up) = trc_out_mass(i_up) &
                        + abs(trc_pth_1trc(ipth)) * dt_all(irivsys(i_up))
                  ENDIF
               ELSE
                  i_dn = pth_down_local(ipth)
                  IF (i_dn > 0 .and. i_dn <= numucat) THEN
                     IF (irivsys(i_dn) > 0 .and. irivsys(i_dn) <= size(dt_all)) THEN
                        trc_out_mass(i_dn) = trc_out_mass(i_dn) &
                           + abs(trc_pth_1trc(ipth)) * dt_all(irivsys(i_dn))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         ! 6b: D2RATE per cell
         ! Clamp donor mass to [0, ...): a spuriously negative trc_mass
         ! must NOT authorise positive outflow scaled to its magnitude
         ! (the old `abs(trc_mass)` let a -5e-10 cell still export at
         ! half-rate of a 1e-9 cell, pushing trc_mass further negative
         ! rather than stopping the leak). True transport bugs that
         ! produce large negative mass are still surfaced via the
         ! check_tracer_state watchdog at trc_mass_neg_warn; this
         ! change just prevents the limiter from papering over them.
         DO i = 1, numucat
            IF (trc_out_mass(i) > 1.e-30_r8) THEN
               rate_cell(i) = min(max(trc_mass(itrc, i), 0._r8) / trc_out_mass(i), 1._r8)
            ELSE
               rate_cell(i) = 1._r8
            ENDIF
         ENDDO

         ! 6c: Get downstream cell rate (for reverse main flow)
         CALL worker_push_data (push_next2ucat, rate_cell, rate_next, fillvalue = 1._r8)

         ! 6d: Scale main-channel flux by sender rate (water direction)
         DO i = 1, numucat
            IF (hflux_fc(i) >= 0._r8) THEN
               trc_flux(i) = trc_flux(i) * rate_cell(i)
            ELSE
               trc_flux(i) = trc_flux(i) * rate_next(i)
            ENDIF
         ENDDO

         ! 6e: Scale bif pathway flux by sender rate (water direction)
         IF (do_bif .and. npthout_local_in > 0) THEN
            CALL worker_push_data (push_bif_dn2pth, rate_cell, rate_dn_pth, fillvalue = 1._r8)
            DO ipth = 1, npthout_local_in
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               pth_wflux = sum(bif_hflux_lev_in(:, ipth))
               IF (pth_wflux >= 0._r8) THEN
                  trc_pth_1trc(ipth) = trc_pth_1trc(ipth) * rate_cell(i_up)
               ELSE
                  trc_pth_1trc(ipth) = trc_pth_1trc(ipth) * rate_dn_pth(ipth)
               ENDIF
            ENDDO
            ! Rebuild bif_net from scaled pathway fluxes
            bif_net(:) = 0._r8
            DO ipth = 1, npthout_local_in
               i_up = pth_upst_local(ipth)
               IF (i_up < 1 .or. i_up > numucat) CYCLE
               bif_net(i_up) = bif_net(i_up) + trc_pth_1trc(ipth)
               i_dn = pth_down_local(ipth)
               IF (i_dn > 0 .and. i_dn <= numucat) THEN
                  bif_net(i_dn) = bif_net(i_dn) - trc_pth_1trc(ipth)
               ENDIF
            ENDDO
            ! Re-scatter scaled flux to remote downstream
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

         ! Re-aggregate upstream tracer flux after limiting
         CALL worker_push_data (push_ups2ucat, trc_flux, flux_ups, fillvalue = 0._r8, mode = 'sum')

         ! --- 8. Update tracer mass ---
         DO i = 1, numucat
            IF (.not. ucatfilter(i)) CYCLE
            dt_i = 0._r8
            IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
            IF (dt_i <= 0._r8) CYCLE

            trc_mass(itrc, i) = trc_mass(itrc, i) &
               + (- trc_flux(i) + flux_ups(i) - bif_net(i)) * dt_i
         ENDDO

         ! --- 9. Save flux for diagnostics (only for active cells) ---
         DO i = 1, numucat
            IF (ucatfilter(i)) THEN
               trc_flux_out(itrc, i) = trc_flux(i)
               trc_bif_net_saved(itrc, i) = bif_net(i)
            ENDIF
         ENDDO

      ENDDO  ! itrc

      ! --- 10. Final concentration (post-update) + dry-cell drain ---
      !
      ! M1 fix: when a cell has run dry (volwater <= trc_v_dry_off) but
      ! still holds trc_mass / trc_inp_buf, fold the orphan mass into
      ! trc_flux_out as an exit flux so it leaves the river system
      ! instead of accumulating as hidden inventory (previously tracked
      ! only by the DBG_DRYTRC counter at MOD_Grid_RiverLakeFlow).
      ! dt_i==0 falls back to simply clearing state; that residual mass
      ! then rejoins the system-wide diagnostic through the conservation
      ! summary's trc_mass_aft drop.
      BLOCK
      real(r8) :: dry_drain
      DO i = 1, numucat
         CALL get_cell_volume(i, wdsrf(i), volresv, ucat2resv, volwater)
         dt_i = 0._r8
         IF (irivsys(i) > 0 .and. irivsys(i) <= size(dt_all)) dt_i = dt_all(irivsys(i))
         DO itrc = 1, ntracers
            IF (volwater <= trc_v_dry_off) THEN
               dry_drain = trc_mass(itrc, i) + trc_inp_buf(itrc, i)
               IF (abs(dry_drain) > trc_tiny) THEN
                  IF (dt_i > 0._r8) THEN
                     trc_flux_out(itrc, i) = trc_flux_out(itrc, i) + dry_drain / dt_i
                  ENDIF
                  IF (allocated(trc_dry_drain)) trc_dry_drain(itrc, i) = trc_dry_drain(itrc, i) + dry_drain
                  trc_mass(itrc, i)    = 0._r8
                  trc_inp_buf(itrc, i) = 0._r8
               ENDIF
            ENDIF
            CALL update_tracer_concentration(itrc, i, volwater)
         ENDDO
      ENDDO
      END BLOCK

      deallocate (conc_next, flux_ups, trc_flux, trc_conc_flux, bif_net)
      deallocate (trc_out_mass, rate_cell, rate_next)
      IF (allocated(conc_dn_pth )) deallocate (conc_dn_pth )
      IF (allocated(trc_pth_1trc)) deallocate (trc_pth_1trc)
      IF (allocated(bif_recv    )) deallocate (bif_recv    )
      IF (allocated(rate_dn_pth )) deallocate (rate_dn_pth )

   END SUBROUTINE tracer_substep


   !-------------------------------------------------------------------------------------
   ! Flush accumulated tracer diagnostics
   !-------------------------------------------------------------------------------------
   SUBROUTINE tracer_flush_acc ()

   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

      ! Reset the *history* accumulators only. acc_trc_inp / acc_rnof_ref are
      ! routing-period state, not history state, and are reset at the end of
      ! each routing period in MOD_Grid_RiverLakeFlow.F90:1263-1264. Resetting
      ! acc_rnof_ref here too would leave acc_trc_inp / acc_rnof_ref with
      ! mismatched windows when history flushes happen between routing periods
      ! (acctime_rnof_max > deltim), inflating m_cap on the next substep.
      IF (numucat > 0) THEN
         IF (allocated(a_trc_conc  )) a_trc_conc   = 0._r8
         IF (allocated(a_trc_out   )) a_trc_out    = 0._r8
         IF (allocated(a_trc_bifout)) a_trc_bifout = 0._r8
      ENDIF

   END SUBROUTINE tracer_flush_acc

   !-------------------------------------------------------------------------------------
   ! Read tracer restart
   !-------------------------------------------------------------------------------------
   SUBROUTINE read_tracer_restart (file_restart, found_restart, missing_mask)

   USE MOD_NetCDFSerial,          only: ncio_var_exist
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   ! .true. iff every tracer's trc_mass* variable was loaded from the file.
   ! Missing variables fall back to zero here but the caller uses this flag
   ! to trigger an init-from-water cold start so non-zero river storage
   ! does not start at zero tracer mass.
   logical, optional, intent(out) :: found_restart
   ! Per-tracer selector, sized (ntracers). On exit, .true. for tracers
   ! whose trc_mass* variable was not present in the file — the caller
   ! then cold-starts only those, preserving successfully loaded ones.
   ! Previously a single missing tracer would clobber every tracer via
   ! the all_found global switch.
   logical, optional, intent(out) :: missing_mask(:)

   integer :: itrc, has_flag
   logical :: has_var, has_active, has_inactive
   logical :: all_found
   character(len=64) :: varname
   real(r8), allocatable :: tmpvec(:)

      all_found = .true.
      IF (present(missing_mask)) missing_mask = .false.

      ! vector_read_and_scatter contains mpi_barrier(p_comm_glb), so
      ! ALL ranks (master, workers, IO) must enter this routine.
      IF (p_is_worker .and. numucat > 0) THEN
         allocate (tmpvec(numucat))
      ELSE
         allocate (tmpvec(0))
      ENDIF

      DO itrc = 1, ntracers
         write(varname, '(A,A)') 'trc_mass_', trim(tracer_names(itrc))

         ! Master-only file probe + broadcast to avoid concurrent opens.
         IF (p_is_master) THEN
            has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
            has_flag = merge(1, 0, has_var)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
         IF (has_flag /= 0) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) trc_mass(itrc, :) = tmpvec(:)
         ELSE
            IF (p_is_master) THEN
               write(varname, '(A,A)') 'trc_mass_active_', trim(tracer_names(itrc))
               has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
               has_flag = merge(1, 0, has_var)
            ENDIF
#ifdef USEMPI
            CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
            has_active = (has_flag /= 0)

            IF (has_active) THEN
               write(varname, '(A,A)') 'trc_mass_inactive_', trim(tracer_names(itrc))
               IF (p_is_master) THEN
                  has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
                  has_flag = merge(1, 0, has_var)
               ENDIF
#ifdef USEMPI
               CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
               has_inactive = (has_flag /= 0)

               write(varname, '(A,A)') 'trc_mass_active_', trim(tracer_names(itrc))
               CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
               IF (p_is_worker .and. numucat > 0) trc_mass(itrc, :) = tmpvec(:)

               IF (has_inactive) THEN
                  write(varname, '(A,A)') 'trc_mass_inactive_', trim(tracer_names(itrc))
                  CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
                  IF (p_is_worker .and. numucat > 0) trc_mass(itrc, :) = trc_mass(itrc, :) + tmpvec(:)
               ENDIF
            ELSE
               IF (p_is_master) THEN
                  write(*,'(A,A,A)') '  Tracer restart variable "', trim(varname), '" not found, cold start from water.'
               ENDIF
               IF (p_is_worker .and. numucat > 0) trc_mass(itrc, :) = 0._r8
               all_found = .false.
               IF (present(missing_mask)) THEN
                  IF (itrc <= size(missing_mask)) missing_mask(itrc) = .true.
               ENDIF
            ENDIF
         ENDIF

         write(varname, '(A,A)') 'trc_inpbuf_', trim(tracer_names(itrc))
         IF (p_is_master) THEN
            has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
            has_flag = merge(1, 0, has_var)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
         IF (has_flag /= 0) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) trc_inp_buf(itrc, :) = tmpvec(:)
         ELSE
            ! M2 fix: trc_inpbuf_* missing on an otherwise-complete restart
            ! (typical when upgrading from a pre-persistence build) used to
            ! silently reset trc_inp_buf to zero, dropping whatever runoff
            ! mass was queued for release at write time. Warn so the user
            ! knows the per-period buffer was dropped.
            !
            ! Do NOT flip missing_mask here: trc_mass_* / trc_levsto_* may
            ! have loaded successfully, and the caller cold-starts the
            ! whole tracer when the mask fires — overwriting valid state.
            ! The worst case (lost pending runoff buffer) is contained:
            ! its mass re-enters via the next tracer_input_from_runoff.
            IF (p_is_master) THEN
               write(*,'(A,A,A)') &
                  ' WARNING read_tracer_restart: "', trim(varname), &
                  '" absent; trc_inp_buf reset to 0 (prognostic state kept).'
            ENDIF
            IF (p_is_worker .and. numucat > 0) trc_inp_buf(itrc, :) = 0._r8
         ENDIF

         ! H1 fix: protected-side tracer pool. On pre-H1 restarts the
         ! variable is absent; previously we left trc_levsto=0, which
         ! silently dropped the protected-side mass until some later
         ! levee_fldstg happened to shift water into/out of levsto.
         ! Now: if present, read directly; if absent but DEF_USE_LEVEE
         ! and levsto>0, inline-seed trc_levsto = levsto*R_init under
         ! the Phase-1 invariant — same recipe as tracer_init_from_water
         ! but without triggering the full cold-start (which would wipe
         ! trc_mass / trc_inp_buf that loaded successfully).
         write(varname, '(A,A)') 'trc_levsto_', trim(tracer_names(itrc))
         IF (p_is_master) THEN
            has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
            has_flag = merge(1, 0, has_var)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
         IF (has_flag /= 0) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) trc_levsto(itrc, :) = tmpvec(:)
         ELSE
            BLOCK
            USE MOD_Grid_RiverLakeLevee, only: has_levee_bf => has_levee, levsto_bf => levsto
            USE MOD_Tracer_Defs,         only: tracers_bf => tracers, delta_to_R_bf => delta_to_R
            integer  :: ii_bf
            real(r8) :: R_bf
            logical  :: reported_bf
            reported_bf = .false.
            IF (DEF_USE_LEVEE .and. p_is_worker .and. numucat > 0 &
                .and. allocated(levsto_bf) .and. allocated(has_levee_bf)) THEN
               R_bf = delta_to_R_bf(tracers_bf(itrc)%init_delta, tracers_bf(itrc)%ref_ratio)
               DO ii_bf = 1, numucat
                  IF (has_levee_bf(ii_bf) .and. levsto_bf(ii_bf) > 0._r8) THEN
                     trc_levsto(itrc, ii_bf) = levsto_bf(ii_bf) * R_bf
                     reported_bf = .true.
                  ENDIF
               ENDDO
            ENDIF
            IF (p_is_master) THEN
               write(*,'(A,A,A)') &
                  ' NOTE read_tracer_restart: "', trim(varname), &
                  '" absent; trc_levsto seeded from current levsto under Phase-1 invariant.'
            ENDIF
            END BLOCK
         ENDIF

         ! Per-tracer routing-period accumulator. Absent in old-format
         ! restarts: keep the zero initialisation so behaviour is unchanged
         ! when the file predates this persistence.
         write(varname, '(A,A)') 'trc_accinp_', trim(tracer_names(itrc))
         IF (p_is_master) THEN
            has_var = ncio_var_exist(file_restart, trim(varname), readflag = .false.)
            has_flag = merge(1, 0, has_var)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
         IF (has_flag /= 0) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, trim(varname), ucat_data_address)
            IF (p_is_worker .and. numucat > 0) acc_trc_inp(itrc, :) = tmpvec(:)
         ENDIF

      ENDDO

      ! Shared-across-tracers runoff reference for the routing period.
      ! Mirrors the per-tracer accinp recovery above.
      !
      ! Warn on desync with the water-side routing accumulator: if the file
      ! contains `acc_rnof_uc` (water accumulator, read by
      ! MOD_Grid_RiverLakeTimeVars.F90:114-120) but is missing
      ! `acc_rnof_ref` (tracer accumulator), the first routing period after
      ! restart will route water that was already accumulated last period
      ! with zero tracer input — producing a spurious dilution pulse
      ! downstream. Most often this happens with mixed-origin restarts
      ! (hand-edited file, tracer enabled mid-run, or a restart written by
      ! an older build without this persistence).
      IF (p_is_master) THEN
         has_var = ncio_var_exist(file_restart, 'acc_rnof_ref', readflag = .false.)
         has_flag = merge(1, 0, has_var)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (has_flag /= 0) THEN
         CALL vector_read_and_scatter (file_restart, tmpvec, numucat, 'acc_rnof_ref', ucat_data_address)
         IF (p_is_worker .and. numucat > 0) acc_rnof_ref(:) = tmpvec(:)
      ELSE
         ! acc_rnof_ref absent — check whether water side has acc_rnof_uc.
         ! If present, backfill: set acc_rnof_ref = acc_rnof_uc and
         ! approximate acc_trc_inp(itrc,:) = acc_rnof_uc(:) * R_default(itrc)
         ! so the already-accumulated runoff water carries its R_init tracer
         ! instead of zero. Under Phase 1 (constant R) this is exact; under
         ! Phase 2 the actual R during the missing period may have differed,
         ! but it's still far better than the zero that would otherwise
         ! produce a one-period dilution pulse downstream.
         IF (p_is_master) THEN
            has_var = ncio_var_exist(file_restart, 'acc_rnof_uc', readflag = .false.)
            has_flag = merge(1, 0, has_var)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (has_flag, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
         IF (has_flag /= 0) THEN
            CALL vector_read_and_scatter (file_restart, tmpvec, numucat, 'acc_rnof_uc', ucat_data_address)
            IF (p_is_worker .and. numucat > 0) THEN
               BLOCK
               USE MOD_Tracer_Defs, only: tracers, delta_to_R
               integer  :: itrc_bf
               real(r8) :: R_bf
               acc_rnof_ref(:) = tmpvec(:)
               DO itrc_bf = 1, ntracers
                  R_bf = delta_to_R(tracers(itrc_bf)%init_delta, tracers(itrc_bf)%ref_ratio)
                  acc_trc_inp(itrc_bf, :) = acc_rnof_ref(:) * R_bf
               ENDDO
               END BLOCK
            ENDIF
            IF (p_is_master) THEN
               write(*,'(A)') '  NOTE (read_tracer_restart): acc_rnof_ref missing but acc_rnof_uc present.'
               write(*,'(A)') '    Backfilled acc_rnof_ref from acc_rnof_uc and approximated acc_trc_inp'
               write(*,'(A)') '    using R_init. Exact under Phase 1; approximate under Phase 2.'
            ENDIF
         ENDIF
      ENDIF

      deallocate (tmpvec)

      IF (present(found_restart)) found_restart = all_found

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
      IF (p_is_worker .and. (.not. allocated(trc_mass))) RETURN

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

         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = trc_inp_buf(itrc, :)
         ENDIF
         write(varname, '(A,A)') 'trc_inpbuf_', trim(tracer_names(itrc))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')

         ! Per-tracer routing-period accumulator. Paired with acc_rnof_ref
         ! below; both must persist so a restart written mid-period
         ! recovers the queued land-tracer input instead of dropping it.
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = acc_trc_inp(itrc, :)
         ENDIF
         write(varname, '(A,A)') 'trc_accinp_', trim(tracer_names(itrc))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')

         ! H1 fix: protected-side tracer pool. Writing unconditionally
         ! so even runs that never touched a levee cell produce a
         ! zero-filled variable, letting restart consumers rely on its
         ! presence without a per-case branch.
         IF (p_is_worker .and. numucat > 0) THEN
            IF (allocated(trc_levsto)) THEN
               tmpvec(:) = trc_levsto(itrc, :)
            ELSE
               tmpvec(:) = 0._r8
            ENDIF
         ENDIF
         write(varname, '(A,A)') 'trc_levsto_', trim(tracer_names(itrc))
         CALL vector_gather_and_write ( &
            tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')
      ENDDO

      ! Shared runoff reference paired with acc_trc_inp for the active
      ! routing period. Written once (not per tracer).
      IF (p_is_worker .and. numucat > 0) THEN
         tmpvec(:) = acc_rnof_ref(:)
      ENDIF
      CALL vector_gather_and_write ( &
         tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, 'acc_rnof_ref', 'ucatch')

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

      ! Output-side FP-dust clamp threshold. Same rationale as
      ! check_tracer_state: trc_mass can carry sub-picomass alternating
      ! residue from upwind flux add/subtract, which propagates to
      ! trc_conc / a_trc_conc / a_trc_out / a_trc_bifout. The underlying
      ! state is preserved intact (so transport / limiter regressions
      ! remain visible via the WARNING watchdog); only this display
      ! slice hides the noise floor.
      BLOCK
      real(r8), parameter :: trc_hist_fp_dust = 1.0e-12_r8

      DO itrc = 1, ntracers

         ! --- Tracer concentration ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = a_trc_conc(itrc, :) / acctime_hist
            WHERE (abs(tmpvec) < trc_hist_fp_dust) tmpvec = 0._r8
         ENDIF

         write(varname, '(A,A)') 'f_trc_conc_', trim(tracer_names(itrc))
         write(longname, '(A,A,A)') 'tracer ratio (', trim(tracer_names(itrc)), ')'

         CALL vector_gather_map2grid_and_write ( tmpvec, numucat,                        &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, trim(varname), 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            trim(longname), 'R')

         ! --- Tracer outflux ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = a_trc_out(itrc, :) / acctime_hist
            WHERE (abs(tmpvec) < trc_hist_fp_dust) tmpvec = 0._r8
         ENDIF

         write(varname, '(A,A)') 'f_trc_flux_', trim(tracer_names(itrc))
         write(longname, '(A,A,A)') 'tracer outflux (', trim(tracer_names(itrc)), ')'

         CALL vector_gather_map2grid_and_write ( tmpvec, numucat,                        &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, trim(varname), 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            trim(longname), 'R*m3/s')

         ! --- Tracer bifurcation net flux ---
         IF (p_is_worker .and. numucat > 0) THEN
            tmpvec(:) = a_trc_bifout(itrc, :) / acctime_hist
            WHERE (abs(tmpvec) < trc_hist_fp_dust) tmpvec = 0._r8
         ENDIF

         write(varname, '(A,A)') 'f_trc_bifout_', trim(tracer_names(itrc))
         write(longname, '(A,A,A)') 'tracer net bifurcation outflux (', trim(tracer_names(itrc)), ')'

         CALL vector_gather_map2grid_and_write ( tmpvec, numucat,                        &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, trim(varname), 'lon_ucat', 'lat_ucat', itime_in_file_ucat,     &
            trim(longname), 'R*m3/s')

      ENDDO
      END BLOCK

      deallocate (tmpvec)

   END SUBROUTINE write_tracer_history


   !-------------------------------------------------------------------------------------
   ! Diagnostic check: print min/max of tracer state variables
   !-------------------------------------------------------------------------------------
   SUBROUTINE check_tracer_state ()

   USE MOD_SPMD_Task
   USE MOD_RangeCheck
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Tracer_Defs, only: tracers, trc_tiny, trc_delta_sanity_max
   USE MOD_Vars_Global, only: spval
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

      ! Display-only thresholds (state itself is never modified here):
      !   trc_mass_fp_dust: absolute trc_mass below this is FP residue
      !       from alternating substep add/subtract (order 1e-15 R*m³);
      !       dust-clamp for log cleanliness only. Well below any
      !       physical input (rnof*R_init is typically >= 1e-9).
      !   trc_mass_neg_warn: trc_mass below this signed threshold is
      !       NOT noise — indicates transport/limiter drift. Emit a
      !       diagnostic once so the root cause can be traced instead
      !       of silently swallowed by a state-side clamp.
      BLOCK
      real(r8), parameter :: trc_mass_fp_dust  = 1.0e-12_r8
      real(r8), parameter :: trc_mass_neg_warn = -1.0e-6_r8
      real(r8) :: worst_neg_mass
      integer  :: icell_worst, nbad

      DO itrc = 1, ntracers
         IF (p_is_worker .and. numucat > 0) THEN
            tmp = trc_mass(itrc,:)
            ! Significant-negative watchdog (state unchanged)
            worst_neg_mass = 0._r8
            icell_worst    = 0
            nbad           = 0
            IF (numucat > 0) THEN
               icell_worst = minloc(tmp, dim=1)
               IF (icell_worst >= 1 .and. icell_worst <= numucat) THEN
                  IF (tmp(icell_worst) < trc_mass_neg_warn) THEN
                     worst_neg_mass = tmp(icell_worst)
                     nbad = count(tmp < trc_mass_neg_warn)
                  ENDIF
               ENDIF
            ENDIF
            IF (worst_neg_mass < trc_mass_neg_warn) THEN
               write(*,'(A,A,A,I6,A,I8,A,E12.5)') &
                  ' WARNING check_tracer_state: ', trim(tracer_names(itrc)), &
                  ' trc_mass < -1e-6 in ', nbad, ' cell(s); worst @ucat=', &
                  icell_worst, ' mass=', worst_neg_mass
            ENDIF
            ! Cosmetic dust clamp for the min/max line only
            WHERE (abs(tmp) < trc_mass_fp_dust) tmp = 0._r8
         ENDIF
         write(label,'(A,A,A)') 'trc_mass_', trim(tracer_names(itrc)), ' [R*m3]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) THEN
            tmp = trc_conc(itrc,:)
            ! Same cosmetic dust clamp on the conc display path — a tiny
            ! negative trc_mass (below fp_dust threshold) propagates into
            ! trc_conc with opposite sign on division; we don't want that
            ! to show up in min/max or in downstream diagnostic readers.
            WHERE (abs(tmp) < trc_mass_fp_dust) tmp = 0._r8
         ENDIF
         write(label,'(A,A,A)') 'trc_conc_', trim(tracer_names(itrc)), ' [R]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) THEN
            tmp = spval
            WHERE (trc_conc(itrc,:) > trc_tiny)
               tmp = (trc_conc(itrc,:) / tracers(itrc)%ref_ratio - 1.0_r8) * 1000.0_r8
            END WHERE
            ! Hard sanity cap: tracer transport around near-dry cells
            ! can transiently spike trc_conc by 10^4+ due to rapid
            ! volwater drops (CFL on water flux vs tracer flux); those
            ! deltas would otherwise dominate the min/max line and
            ! hide the healthy pool range. Mirror the land-side
            ! trc_delta_sanity_max clamp.
            WHERE (tmp /= spval .and. abs(tmp) > trc_delta_sanity_max) &
               tmp = spval
         ENDIF
         write(label,'(A,A,A)') 'trc_delta_', trim(tracer_names(itrc)), ' [permil]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) THEN
            tmp = trc_flux_out(itrc,:)
            ! FP-dust clamp on the flux display path (same rationale
            ! as trc_mass / trc_conc above). Real flux magnitudes are
            ! many orders above trc_mass_fp_dust so this only hides
            ! the alternating +/- 1e-15 numerical residue.
            WHERE (abs(tmp) < trc_mass_fp_dust) tmp = 0._r8
         ENDIF
         write(label,'(A,A,A)') 'trc_outflux_', trim(tracer_names(itrc)), ' [R*m3/s]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) tmp = acc_trc_inp(itrc,:)
         write(label,'(A,A,A)') 'trc_inp_', trim(tracer_names(itrc)), ' [R*m3]'
         CALL check_vector_data (label, tmp)

         IF (p_is_worker .and. numucat > 0) tmp = trc_inp_buf(itrc,:)
         write(label,'(A,A,A)') 'trc_inpbuf_', trim(tracer_names(itrc)), ' [R*m3]'
         CALL check_vector_data (label, tmp)

         ! H1 fix: surface protected-side pool so the levee-on case has a
         ! min/max line alongside trc_mass. Same FP-dust clamp rationale.
         IF (p_is_worker .and. numucat > 0) THEN
            IF (allocated(trc_levsto)) THEN
               tmp = trc_levsto(itrc,:)
               WHERE (abs(tmp) < trc_mass_fp_dust) tmp = 0._r8
            ELSE
               tmp = 0._r8
            ENDIF
         ENDIF
         write(label,'(A,A,A)') 'trc_levsto_', trim(tracer_names(itrc)), ' [R*m3]'
         CALL check_vector_data (label, tmp)
      ENDDO
      END BLOCK

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
      IF (allocated(trc_inp_buf       )) deallocate (trc_inp_buf       )
      IF (allocated(acc_rnof_ref      )) deallocate (acc_rnof_ref      )
      IF (allocated(trc_bif_net_saved)) deallocate (trc_bif_net_saved)
      IF (allocated(a_trc_conc       )) deallocate (a_trc_conc       )
      IF (allocated(a_trc_out        )) deallocate (a_trc_out        )
      IF (allocated(a_trc_bifout     )) deallocate (a_trc_bifout     )
      IF (allocated(trc_dry_drain    )) deallocate (trc_dry_drain    )
      IF (allocated(trc_levsto       )) deallocate (trc_levsto       )

   END SUBROUTINE tracer_final

END MODULE MOD_Grid_RiverLakeTracer
#endif
