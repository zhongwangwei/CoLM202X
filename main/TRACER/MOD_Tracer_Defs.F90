#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Defs

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_TRACER_NUM, DEF_TRACER_NAMES, DEF_TRACER_TYPES, &
      DEF_TRACER_MRAT, DEF_TRACER_REF_RATIO, DEF_TRACER_INIT_DELTA, &
      DEF_TRACER_REACTIVE_DECAY_RATE, DEF_TRACER_USE_FRACTIONATION

   IMPLICIT NONE
   SAVE

   type :: tracer_info_type
      character(len=32) :: name
      character(len=16) :: category
      real(r8)          :: mol_weight
      real(r8)          :: ref_ratio
      real(r8)          :: init_delta
      real(r8)          :: reactive_decay_rate
      logical           :: has_fractionation
   end type tracer_info_type

   integer :: ntracers = 0
   type(tracer_info_type), allocatable :: tracers(:)

	   real(r8), parameter :: Rsmow_18O = 2.0052e-3_r8
	   real(r8), parameter :: Rsmow_D   = 1.5576e-4_r8
	   real(r8), parameter :: trc_tiny = 1.0e-30_r8
	   ! ------------------------------------------------------------------
	   ! DELTA DIAGNOSTIC INVARIANT
	   !
	   ! Any official isotope delta diagnostic must aggregate tracer material
	   ! and water first, then convert once:
	   !
	   !    delta = mass_to_delta(sum(tracer_mass), sum(water_mass), ref_ratio)
	   !
	   ! Never average per-step/per-cell delta values, and do not average R
	   ! unless that R is explicitly weighted by its matching water amount.
	   ! Arithmetic means of delta or unweighted R bias low-water states and can
	   ! create apparent isotope drift even when tracer mass is conserved.
	   ! ------------------------------------------------------------------
	   ! Physical minimum water mass (kg/m² ≈ mm) for delta diagnostics.
   ! Below this the water pool is floating-point noise (e.g. wa ≈ 1e-5
   ! mm in a patch that just went near zero), and R_sample = trc/water
   ! amplifies round-off into ±10000‰ deltas that are not physical.
   ! 1e-3 kg/m² = 1 μm water depth; this keeps legitimate trace water
   ! (dew/frost puddles, residual soil moisture) in the diagnostic while
   ! filtering pure FP dust. Previously 1e-6 (1 nm) was too permissive —
   ! flux deltas (qinfl/qcharge) with water_mass = rate*deltim near the
   ! FP floor produced arbitrarily large ratios that flowed into
   ! trc_delta_* history / check outputs.
   real(r8), parameter :: trc_water_min_for_delta = 1.0e-3_r8
   ! Flux isotope deltas need a stricter accumulated-water denominator than
   ! storage deltas. Below this history-window water loss, the flux signature is
   ! dominated by NSS/storage tendency or round-off rather than a resolved flux.
   real(r8), parameter :: trc_flux_water_min_for_delta = 1.0e-1_r8
   ! Hard sanity cap on |delta| (permil). Even extreme Rayleigh
   ! depletion in snowpack or evap residuals rarely exceeds a few
   ! hundred permil for meteoric isotope ratios; anything past this is
   ! a sign of tracer/water state divergence, not real physics. Used
   ! by delta diagnostics to clamp to spval so non-physical spikes do
   ! not poison min/max reports or downstream post-processing. 2000 ‰
   ! gives comfortable headroom for strongly fractionated pools.
   real(r8), parameter :: trc_delta_sanity_max = 2.0e3_r8

   PUBLIC :: tracer_defs_init, tracer_defs_final
   PUBLIC :: mass_to_delta, delta_to_R, R_to_mass
   PUBLIC :: tracer_is_isotope, tracer_is_conservative, tracer_is_reactive
   PUBLIC :: tracer_uses_delta_diagnostics, tracer_can_use_fixed_signature
   PUBLIC :: tracer_init_water_ratio, tracer_reactive_decay_fraction
   PUBLIC :: ntracers, tracers
   PUBLIC :: tracer_info_type
   PUBLIC :: Rsmow_18O, Rsmow_D, trc_tiny, trc_water_min_for_delta, &
             trc_flux_water_min_for_delta, trc_delta_sanity_max

CONTAINS

   SUBROUTINE tracer_defs_init ()
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer :: i
      character(len=32), allocatable :: tokens(:)
      integer :: ntokens, n_stored, tok_cap
      integer :: j, sfx_len, max_base
      character(len=32) :: sfx, base

      ntracers = DEF_TRACER_NUM
      IF (ntracers < 0) THEN
         IF (p_is_master) WRITE(*,'(A,I0)') 'ERROR tracer_defs_init: DEF_TRACER_NUM must be >= 0, got ', ntracers
         CALL CoLM_stop()
      ENDIF
      IF (ntracers > 1000) THEN
         IF (p_is_master) WRITE(*,'(A,I0,A)') 'ERROR tracer_defs_init: DEF_TRACER_NUM=', ntracers, &
            ' is unreasonably large; raise the guard intentionally if needed.'
         CALL CoLM_stop()
      ENDIF
      IF (allocated(tracers)) RETURN   ! already initialised (e.g. by earlier all-rank call)
      IF (ntracers == 0) RETURN
      allocate(tracers(ntracers))

      ! Size tokens() to cover ntracers exactly (at least 16) so a user with
      ! DEF_TRACER_NUM > 100 no longer silently drops excess entries to the
      ! `tracer_N` default. parse_csv still guards against a CSV longer than
      ! the array.
      tok_cap = max(ntracers, 16)
      allocate(tokens(tok_cap))

      CALL parse_csv(DEF_TRACER_NAMES, tokens, ntokens)
      ! parse_csv keeps counting past the tokens() capacity to
      ! report the raw comma count; only the first `size(tokens)` entries
      ! are actually populated. Use `n_stored` as the real upper bound so
      ! the fill-default loop below covers slots that ntokens "logically"
      ! claims but parse_csv had no room to write.
      n_stored = min(ntokens, size(tokens))
      DO i = 1, min(ntracers, n_stored)
         tracers(i)%name = sanitize_ncname(ADJUSTL(TRIM(tokens(i))))
      ENDDO
      ! Fill default names when DEF_TRACER_NAMES has fewer stored tokens
      ! than DEF_TRACER_NUM. Without this, tracers(i)%name past n_stored
      ! would be an uninitialised character that later gets spliced into
      ! NetCDF variable names and log labels.
      DO i = n_stored + 1, ntracers
         write(tracers(i)%name, '(A,I0)') 'tracer_', i
      ENDDO
      ! Ensure uniqueness after sanitisation. Two raw tokens that differ
      ! only by illegal characters (e.g. 'H2O-16' vs 'H2O_16') collapse
      ! to the same clean token; MOD_Hist would then emit duplicate
      ! NetCDF variable names and crash. Mirrors the river-side
      ! dedup loop in MOD_Grid_RiverLakeTracer:tracer_init.
         DO i = 1, ntracers
            DO j = 1, i - 1
               IF (trim(tracers(i)%name) == trim(tracers(j)%name)) THEN
                  ! Budget the base name so `base // '_' // i` survives the
                  ! len=32 truncation. Without this, two raw tokens that
                  ! already fill 32 chars collide again post-suffix, yielding
                  ! a duplicate NetCDF variable name from MOD_Hist.
                  write(sfx, '(A,I0)') '_', i
                  sfx_len = len_trim(sfx)
                  max_base = max(1, 32 - sfx_len)
                  base = tracers(i)%name
                  IF (len_trim(base) > max_base) base = base(1:max_base)
                  write(tracers(i)%name, '(A,A)') trim(base), trim(sfx)
                  EXIT
               ENDIF
            ENDDO
         ENDDO

      CALL parse_csv(DEF_TRACER_TYPES, tokens, ntokens)
      ! Same n_stored cap as the name loop above — avoids leaving
      ! tracers(size(tokens)+1 : ntokens)%category uninitialised when the
      ! TYPES CSV (or DEF_TRACER_NUM) exceeds the tokens() capacity.
      n_stored = min(ntokens, size(tokens))
      DO i = 1, min(ntracers, n_stored)
         tracers(i)%category = canonical_tracer_category(tokens(i))
      ENDDO
      DO i = n_stored + 1, ntracers
         tracers(i)%category = 'isotope'
      ENDDO

      CALL parse_csv_real(DEF_TRACER_MRAT, ntracers, tracers(:)%mol_weight, 18.0_r8)
      CALL parse_csv_real(DEF_TRACER_REF_RATIO, ntracers, tracers(:)%ref_ratio, 1.0_r8)
      CALL parse_csv_real(DEF_TRACER_INIT_DELTA, ntracers, tracers(:)%init_delta, 0.0_r8)
      CALL parse_csv_real(DEF_TRACER_REACTIVE_DECAY_RATE, ntracers, &
         tracers(:)%reactive_decay_rate, 0.0_r8)
      IF (DEF_TRACER_USE_FRACTIONATION .and. p_is_master) THEN
         WRITE(*,'(A)') 'WARNING tracer_defs_init: isotope fractionation is experimental; enabled paths include ' // &
            'precipitation/evaporation, phase change, transpiration, wetland, glacier, and waterbody approximations.'
      ENDIF

      DO i = 1, ntracers
         IF (.not. tracer_is_isotope(i) .and. .not. tracer_is_conservative(i) .and. &
             .not. tracer_is_reactive(i)) THEN
            IF (p_is_master) THEN
               write(*,'(A,A,A,A)') ' WARNING tracer_defs_init: unknown tracer category "', &
                  trim(tracers(i)%category), '" for ', trim(tracers(i)%name)
               write(*,'(A)') '   Falling back to isotope behavior.'
            ENDIF
            tracers(i)%category = 'isotope'
         ENDIF
         IF (tracers(i)%reactive_decay_rate < 0._r8) THEN
            IF (p_is_master) THEN
               write(*,'(A,A,A)') ' WARNING tracer_defs_init: negative reactive decay rate for ', &
                  trim(tracers(i)%name), '; reset to zero.'
            ENDIF
            tracers(i)%reactive_decay_rate = 0._r8
         ENDIF
         tracers(i)%has_fractionation = tracer_is_isotope(i) .and. DEF_TRACER_USE_FRACTIONATION
      ENDDO
      deallocate(tokens)
   END SUBROUTINE tracer_defs_init

   SUBROUTINE tracer_defs_final ()
      IF (allocated(tracers)) deallocate(tracers)
      ntracers = 0
   END SUBROUTINE tracer_defs_final

   !-------------------------------------------------------------------
   ! Tracer category framework.
   !
   ! The transport math is intentionally expressed as "water flux times
   ! tracer per water" for all categories. For isotopes this is R_sample;
   ! for conservative/reactive tracers it is concentration per unit water.
   ! Phase-1 still uses DEF_TRACER_INIT_DELTA for both concepts to keep
   ! the namelist stable; a later conservative/reactive implementation
   ! should split that into explicit INIT_CONC / forcing concentration.
   !-------------------------------------------------------------------
   FUNCTION canonical_tracer_category (raw) RESULT(category)
      character(len=*), intent(in) :: raw
      character(len=16) :: category
      integer :: i, ia

      category = ADJUSTL(TRIM(raw))
      DO i = 1, len_trim(category)
         ia = iachar(category(i:i))
         IF (ia >= iachar('A') .and. ia <= iachar('Z')) THEN
            category(i:i) = achar(ia + iachar('a') - iachar('A'))
         ENDIF
      ENDDO
   END FUNCTION canonical_tracer_category

   logical FUNCTION tracer_is_isotope (itrc)
      integer, intent(in) :: itrc
      tracer_is_isotope = allocated(tracers) .and. itrc >= 1 .and. itrc <= ntracers .and. &
         trim(tracers(itrc)%category) == 'isotope'
   END FUNCTION tracer_is_isotope

   logical FUNCTION tracer_is_conservative (itrc)
      integer, intent(in) :: itrc
      tracer_is_conservative = allocated(tracers) .and. itrc >= 1 .and. itrc <= ntracers .and. &
         trim(tracers(itrc)%category) == 'conservative'
   END FUNCTION tracer_is_conservative

   logical FUNCTION tracer_is_reactive (itrc)
      integer, intent(in) :: itrc
      tracer_is_reactive = allocated(tracers) .and. itrc >= 1 .and. itrc <= ntracers .and. &
         trim(tracers(itrc)%category) == 'reactive'
   END FUNCTION tracer_is_reactive

   logical FUNCTION tracer_uses_delta_diagnostics (itrc)
      integer, intent(in) :: itrc
      ! Delta/NSS diagnostics are isotope-specific. Conservative/reactive
      ! tracers store concentration per unit water and must not be converted
      ! with ref_ratio.
      tracer_uses_delta_diagnostics = tracer_is_isotope(itrc)
   END FUNCTION tracer_uses_delta_diagnostics

   logical FUNCTION tracer_can_use_fixed_signature (itrc)
      integer, intent(in) :: itrc
      ! Fixed-signature resync/snap is a Phase-1 shortcut for isotope
      ! tracers with inactive fractionation and no runtime forcing. Ordinary
      ! conservative/reactive tracers may carry spatial gradients, restart
      ! memory, or source/sink history, so they must preserve the current
      ! pool-mixed concentration instead of being snapped to init_delta.
      tracer_can_use_fixed_signature = tracer_is_isotope(itrc)
   END FUNCTION tracer_can_use_fixed_signature

   real(r8) FUNCTION tracer_init_water_ratio (itrc)
      integer, intent(in) :: itrc

      tracer_init_water_ratio = 0._r8
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      IF (tracer_is_isotope(itrc)) THEN
         tracer_init_water_ratio = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
      ELSE
         tracer_init_water_ratio = tracers(itrc)%init_delta
      ENDIF
   END FUNCTION tracer_init_water_ratio

   real(r8) FUNCTION tracer_reactive_decay_fraction (itrc, deltim)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: deltim
      real(r8) :: kdt

      tracer_reactive_decay_fraction = 0._r8
      IF (.not. tracer_is_reactive(itrc)) RETURN
      IF (deltim <= 0._r8) RETURN
      IF (tracers(itrc)%reactive_decay_rate <= 0._r8) RETURN

      kdt = min(tracers(itrc)%reactive_decay_rate * deltim, 700._r8)
      tracer_reactive_decay_fraction = 1._r8 - exp(-kdt)
   END FUNCTION tracer_reactive_decay_fraction

	   real(r8) FUNCTION mass_to_delta (trc_mass, water_mass, ref_ratio)
	      USE MOD_Vars_Global, only: spval
	      real(r8), intent(in) :: trc_mass, water_mass, ref_ratio
	      real(r8) :: R_sample
	      ! Public delta conversion point. Callers should pass already-aggregated
	      ! tracer and water amounts for the requested history/check window.
	      ! Do not call this on per-step values and then average the returned
	      ! delta; aggregate mass and water first, convert last.
	      ! abs(water_mass) guard: signed wa (aquifer-debt from wetland deficit)
      ! and its matching signed trc_wa still produce a well-defined ratio
      ! (both numerator and denominator share sign → R_sample positive).
      ! Physically non-negative pools (wliq/wice/wdsrf/wetwat/scv/ldew/flux)
      ! are unaffected since abs() equals value there.
      IF (abs(water_mass) > trc_tiny .and. ref_ratio > trc_tiny) THEN
         R_sample = trc_mass / water_mass
         IF (R_sample < -trc_tiny) THEN
            ! Sign mismatch between trc_mass and water_mass: upstream state
            ! machine is broken (e.g. trc_wa and wa drifted to opposite
            ! signs). Returning spval rather than silently clamping so
            ! history readers and assert tools can flag the cell instead
            ! of treating a bogus delta as physical.
#if (defined CoLMDEBUG)
            write(*,'(A,E12.5,A,E12.5)') &
               ' WARNING mass_to_delta: R_sample<0 trc=', trc_mass, ' water=', water_mass
#endif
            mass_to_delta = spval
         ELSE
            mass_to_delta = (R_sample / ref_ratio - 1.0_r8) * 1000.0_r8
         ENDIF
      ELSE
         mass_to_delta = 0.0_r8
      ENDIF
   END FUNCTION mass_to_delta

   pure real(r8) FUNCTION delta_to_R (delta, ref_ratio)
      real(r8), intent(in) :: delta, ref_ratio
      delta_to_R = ref_ratio * (1.0_r8 + delta / 1000.0_r8)
   END FUNCTION delta_to_R

   pure real(r8) FUNCTION R_to_mass (delta, water_mass, ref_ratio)
      real(r8), intent(in) :: delta, water_mass, ref_ratio
      R_to_mass = water_mass * delta_to_R(delta, ref_ratio)
   END FUNCTION R_to_mass

   !-------------------------------------------------------------------
   ! Sanitize a string for use as a NetCDF variable name component.
   ! Keeps alphanumerics, underscore and period; replaces anything else
   ! by dropping it. Mirrors MOD_Grid_RiverLakeTracer.sanitize_ncname so
   ! land-side history/diagnostics names (which splice tracers(i)%name)
   ! stay consistent with river-side output naming.
   !-------------------------------------------------------------------
   FUNCTION sanitize_ncname (raw) RESULT(clean)
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
      IF (len_trim(clean) == 0) clean = 'unnamed'
   END FUNCTION sanitize_ncname

   SUBROUTINE parse_csv (csvstr, tokens, ntokens)
      USE MOD_SPMD_Task, only: p_is_master
      character(len=*), intent(in)  :: csvstr
      character(len=32), intent(out) :: tokens(:)
      integer, intent(out) :: ntokens
      integer :: i, j, slen, cap
      logical :: warned
      ! deferred-length buffer absorbs whatever length the
      ! caller's csvstr carries. The previous fixed-256 buf silently
      ! truncated a longer DEF_TRACER_* namelist string (e.g. a future
      ! ntracers > 8 with full-length names) and parse_csv would then
      ! drop the trailing tokens without warning.
      character(len=:), allocatable :: buf
      buf = ADJUSTL(TRIM(csvstr))
      slen = LEN_TRIM(buf)
      ntokens = 0; j = 1
      ! Cap writes at size(tokens) so a pathological DEF_TRACER_* with
      ! more than size(tokens) entries can't write past the end of the
      ! caller's fixed array (tracer_defs_init / parse_csv_real both pass
      ! tokens(100)). We still keep counting ntokens so the caller sees the
      ! raw comma count, but the warning fires once per invocation.
      cap = size(tokens)
      warned = .false.
      DO i = 1, slen
         IF (buf(i:i) == ',') THEN
            ntokens = ntokens + 1
            IF (ntokens <= cap) THEN
               tokens(ntokens) = ADJUSTL(TRIM(buf(j:i-1)))
            ELSEIF (.not. warned) THEN
               IF (p_is_master) THEN
                  write(*,'(A,I0,A)') ' WARNING: parse_csv: token count exceeds capacity (', &
                     cap, '); excess entries ignored.'
               ENDIF
               warned = .true.
            ENDIF
            j = i + 1
         ENDIF
      ENDDO
      IF (j <= slen) THEN
         ntokens = ntokens + 1
         IF (ntokens <= cap) THEN
            tokens(ntokens) = ADJUSTL(TRIM(buf(j:slen)))
         ELSEIF (.not. warned) THEN
            IF (p_is_master) THEN
               write(*,'(A,I0,A)') ' WARNING: parse_csv: token count exceeds capacity (', &
                  cap, '); excess entries ignored.'
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE parse_csv

   SUBROUTINE parse_csv_real (csvstr, n, vals, default_val)
      character(len=*), intent(in)  :: csvstr
      integer, intent(in) :: n
      real(r8), intent(out) :: vals(n)
      real(r8), intent(in) :: default_val
      character(len=32), allocatable :: tokens(:)
      integer :: ntokens, i, ios, tok_cap
      ! Match tokens capacity to the caller's requested n so a DEF_TRACER_*
      ! CSV with >100 entries no longer gets silently truncated to the
      ! default_val fallback. parse_csv still reports a WARNING if the
      ! comma count exceeds the array.
      tok_cap = max(n, 16)
      allocate(tokens(tok_cap))
      CALL parse_csv(csvstr, tokens, ntokens)
      DO i = 1, n
         IF (i <= ntokens .and. i <= size(tokens)) THEN
            READ(tokens(i), *, iostat=ios) vals(i)
            IF (ios /= 0) vals(i) = default_val
         ELSE
            vals(i) = default_val
         ENDIF
      ENDDO
      deallocate(tokens)
   END SUBROUTINE parse_csv_real

END MODULE MOD_Tracer_Defs
#endif
