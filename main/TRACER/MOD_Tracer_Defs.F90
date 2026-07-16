#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Defs

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_TRACER_NUM, DEF_TRACER_NAMES, DEF_TRACER_TYPES, &
      DEF_TRACER_MRAT, DEF_TRACER_REF_RATIO, DEF_TRACER_INIT_DELTA, &
      DEF_TRACER_REACTIVE_DECAY_RATE, DEF_TRACER_PARAM_FILES, &
      DEF_TRACER_USE_FRACTIONATION
   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite

   IMPLICIT NONE
   SAVE

   type :: tracer_info_type
      character(len=32) :: name
      character(len=16) :: category
      ! unit_kind controls diagnostic labels only; it does not rescale forcing
      ! or state values. Numeric inputs must already use the declared ratio.
      character(len=32) :: unit_kind = 'tracer_per_water'
      real(r8)          :: mol_weight
      real(r8)          :: ref_ratio
      real(r8)          :: init_delta
      ! Explicit non-isotope concentrations.  init_delta remains the legacy
      ! fallback so existing namelists retain their former signatures.
      real(r8)          :: init_conc
      real(r8)          :: precip_default_conc
      real(r8)          :: vapor_default_conc
      real(r8)          :: reactive_decay_rate
      ! Ionic charge is registry metadata; no electroneutrality solver exists.
      integer           :: charge = 0
      logical           :: has_fractionation
      ! Species-owned reactive implementations may opt out of the generic
      ! canopy/soil/snow/aquifer tracer pools during registry refresh.
      logical           :: uses_land_water_transport = .true.
      logical           :: is_nonvolatile = .false.
      logical           :: uses_reaction = .false.
   end type tracer_info_type

   type :: tracer_parameter_type
      character(len=32) :: unit_kind = 'tracer_per_water'
      real(r8) :: mol_weight = 18.0_r8
      real(r8) :: ref_ratio = 1.0_r8
      real(r8) :: init_delta = 0.0_r8
      real(r8) :: init_conc = huge(1.0_r8)
      real(r8) :: precip_default_conc = huge(1.0_r8)
      real(r8) :: vapor_default_conc = huge(1.0_r8)
      real(r8) :: reactive_decay_rate = 0.0_r8
      integer  :: charge = 0
      logical  :: uses_generic_land_water_transport = .true.
      logical  :: is_nonvolatile = .false.
      logical  :: uses_reaction = .false.
   end type tracer_parameter_type

   integer :: ntracers = 0
   type(tracer_info_type), allocatable :: tracers(:)

      real(r8), parameter :: Rsmow_18O = 2.0052e-3_r8
      real(r8), parameter :: Rsmow_D   = 1.5576e-4_r8
      real(r8), parameter :: trc_tiny = 1.0e-30_r8
      ! Numerical zero and physical denominator floors are deliberately
      ! separate.  trc_tiny is only an arithmetic nonzero guard.  Transport
      ! ratios from finite water pools need a larger physically resolved water
      ! floor, otherwise near-dry pools can export arbitrarily large R =
      ! tracer/water before downstream sanity caps see the problem.
      ! Units follow land water storage conventions (kg/m2 ~= mm).
      real(r8), parameter :: trc_water_min_for_ratio = 1.0e-12_r8
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
   PUBLIC :: tracer_is_particle, tracer_is_nonvolatile_solute, tracer_uses_land_water_transport
   PUBLIC :: tracer_set_land_water_transport
   PUBLIC :: tracer_uses_delta_diagnostics, tracer_can_use_fixed_signature
   PUBLIC :: tracer_concentration_units
   PUBLIC :: tracer_init_water_ratio, tracer_precip_default_ratio, tracer_vapor_default_ratio
   PUBLIC :: tracer_reactive_decay_fraction
   PUBLIC :: tracer_param_file_for_index
   PUBLIC :: tracer_lower, tracer_upper
   PUBLIC :: ntracers, tracers
   PUBLIC :: tracer_info_type
      PUBLIC :: Rsmow_18O, Rsmow_D, trc_tiny, trc_water_min_for_ratio, &
                trc_water_min_for_delta, trc_flux_water_min_for_delta, &
                trc_delta_sanity_max

CONTAINS

   SUBROUTINE tracer_defs_init ()
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer :: i
      character(len=32), allocatable :: tokens(:)
      integer :: ntokens, n_stored, tok_cap
      integer :: j, sfx_len, max_base, suffix_index
      character(len=32) :: sfx, base
      logical :: duplicate_name

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
      CALL warn_csv_count('DEF_TRACER_NAMES', ntokens, ntracers)
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
      ! NetCDF variable names and crash. River/lake reuses these names.
         DO i = 1, ntracers
            base = tracers(i)%name
            suffix_index = i
            DO
               duplicate_name = .false.
               DO j = 1, i - 1
                  IF (trim(tracers(i)%name) == trim(tracers(j)%name)) THEN
                     duplicate_name = .true.
                     EXIT
                  ENDIF
               ENDDO
               IF (.not. duplicate_name) EXIT

               ! Re-scan after adding a suffix: the candidate itself may
               ! already be a user-provided name (A,A_3,A is the minimal
               ! example). Advance deterministically until the closure is
               ! unique across every prior tracer.
               write(sfx, '(A,I0)') '_', suffix_index
               suffix_index = suffix_index + 1
               sfx_len = len_trim(sfx)
               max_base = max(1, 32 - sfx_len)
               IF (len_trim(base) > max_base) THEN
                  write(tracers(i)%name, '(A,A)') base(1:max_base), trim(sfx)
               ELSE
                  write(tracers(i)%name, '(A,A)') trim(base), trim(sfx)
               ENDIF
            ENDDO
         ENDDO

      CALL parse_csv(DEF_TRACER_TYPES, tokens, ntokens)
      CALL warn_csv_count('DEF_TRACER_TYPES', ntokens, ntracers)
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

      DO i = 1, ntracers
         CALL set_tracer_category_defaults (i)
      ENDDO

      CALL parse_csv_real(DEF_TRACER_MRAT, ntracers, tracers(:)%mol_weight, 18.0_r8, &
         'DEF_TRACER_MRAT')
      CALL parse_csv_real(DEF_TRACER_REF_RATIO, ntracers, tracers(:)%ref_ratio, 1.0_r8, &
         'DEF_TRACER_REF_RATIO')
      CALL parse_csv_real(DEF_TRACER_INIT_DELTA, ntracers, tracers(:)%init_delta, 0.0_r8, &
         'DEF_TRACER_INIT_DELTA')
      ! Legacy main-namelist input used one scalar for initialization and both
      ! forcing fallbacks.  Parameter-file fields below may override each one
      ! independently without changing old configurations.
      tracers(:)%init_conc = tracers(:)%init_delta
      tracers(:)%precip_default_conc = tracers(:)%init_delta
      tracers(:)%vapor_default_conc = tracers(:)%init_delta
      CALL parse_csv_real(DEF_TRACER_REACTIVE_DECAY_RATE, ntracers, &
         tracers(:)%reactive_decay_rate, 0.0_r8, 'DEF_TRACER_REACTIVE_DECAY_RATE')

      CALL apply_tracer_param_files ()

      IF (DEF_TRACER_USE_FRACTIONATION .and. p_is_master) THEN
         WRITE(*,'(A)') 'WARNING tracer_defs_init: isotope fractionation is experimental; enabled paths include ' // &
            'precipitation/evaporation, phase change, transpiration, wetland, glacier, and waterbody approximations.'
      ENDIF

      DO i = 1, ntracers
         IF (.not. tracer_is_isotope(i) .and. .not. tracer_is_conservative(i) .and. &
             .not. tracer_is_reactive(i) .and. .not. tracer_is_particle(i)) THEN
            IF (p_is_master) THEN
               write(*,'(A,A,A,A)') 'ERROR tracer_defs_init: unknown tracer category "', &
                  trim(tracers(i)%category), '" for ', trim(tracers(i)%name)
               write(*,'(A)') '   Valid categories: isotope, conservative, reactive, particle.'
            ENDIF
            CALL CoLM_stop()
         ENDIF
         IF (.not. ieee_is_finite(tracers(i)%reactive_decay_rate)) THEN
            IF (p_is_master) THEN
               write(*,'(A,A,A)') 'ERROR tracer_defs_init: non-finite reactive_decay_rate for ', &
                  trim(tracers(i)%name), '.'
            ENDIF
            CALL CoLM_stop()
         ELSEIF (tracers(i)%reactive_decay_rate < 0._r8) THEN
            IF (p_is_master) THEN
               write(*,'(A,A,A)') ' WARNING tracer_defs_init: negative reactive decay rate for ', &
                  trim(tracers(i)%name), '; reset to zero.'
            ENDIF
            tracers(i)%reactive_decay_rate = 0._r8
         ENDIF
         CALL validate_tracer_descriptor (i)
         tracers(i)%has_fractionation = tracer_is_isotope(i) .and. DEF_TRACER_USE_FRACTIONATION
      ENDDO
      deallocate(tokens)
   END SUBROUTINE tracer_defs_init


   SUBROUTINE apply_tracer_param_files ()
      IMPLICIT NONE
      integer :: itrc
      character(len=512) :: file_param
      logical :: found

      IF (.not. allocated(tracers)) RETURN
      DO itrc = 1, ntracers
         CALL tracer_param_file_for_index (itrc, '', file_param, found)
         IF (found) CALL read_tracer_parameter_file (itrc, file_param)
      ENDDO
   END SUBROUTINE apply_tracer_param_files

   SUBROUTINE read_tracer_parameter_file (itrc, nlfile)
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer, intent(in) :: itrc
      character(len=*), intent(in) :: nlfile

      type(tracer_parameter_type) :: DEF_TRACER
      logical :: fexists
      integer :: ierr, unit_nml
      character(len=512) :: iomsg
      namelist /nl_colm_tracer_parameter/ DEF_TRACER

      IF (itrc < 1 .or. itrc > ntracers) RETURN
      IF (len_trim(nlfile) <= 0 .or. trim(tracer_lower(nlfile)) == 'null') RETURN

      INQUIRE (file=trim(nlfile), exist=fexists)
      IF (.not. fexists) THEN
         IF (p_is_master) write(*,'(A,A)') 'ERROR read_tracer_parameter_file: missing tracer parameter file: ', &
            trim(nlfile)
         CALL CoLM_stop()
      ENDIF
      IF (.not. tracer_parameter_group_present(nlfile)) RETURN

      DEF_TRACER%unit_kind          = tracers(itrc)%unit_kind
      DEF_TRACER%mol_weight          = tracers(itrc)%mol_weight
      DEF_TRACER%ref_ratio           = tracers(itrc)%ref_ratio
      DEF_TRACER%init_delta          = tracers(itrc)%init_delta
      DEF_TRACER%init_conc           = huge(1.0_r8)
      DEF_TRACER%precip_default_conc = huge(1.0_r8)
      DEF_TRACER%vapor_default_conc  = huge(1.0_r8)
      DEF_TRACER%reactive_decay_rate = tracers(itrc)%reactive_decay_rate
      DEF_TRACER%charge              = tracers(itrc)%charge
      DEF_TRACER%uses_generic_land_water_transport = tracers(itrc)%uses_land_water_transport
      DEF_TRACER%is_nonvolatile      = tracers(itrc)%is_nonvolatile
      DEF_TRACER%uses_reaction       = tracers(itrc)%uses_reaction

      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form='FORMATTED')
      iomsg = ''
      read(unit_nml, nml=nl_colm_tracer_parameter, iostat=ierr, iomsg=iomsg)
      close(unit_nml)
      IF (ierr /= 0) THEN
         IF (p_is_master) THEN
            write(*,'(3A)') 'ERROR read_tracer_parameter_file: invalid &nl_colm_tracer_parameter in ', &
               trim(nlfile), ': '
            write(*,'(A)') trim(iomsg)
         ENDIF
         CALL CoLM_stop()
      ENDIF

      tracers(itrc)%unit_kind          = tracer_lower(DEF_TRACER%unit_kind)
      tracers(itrc)%mol_weight          = DEF_TRACER%mol_weight
      tracers(itrc)%ref_ratio           = DEF_TRACER%ref_ratio
      tracers(itrc)%init_delta          = DEF_TRACER%init_delta
      tracers(itrc)%init_conc           = DEF_TRACER%init_delta
      tracers(itrc)%precip_default_conc = DEF_TRACER%init_delta
      tracers(itrc)%vapor_default_conc  = DEF_TRACER%init_delta
      IF (DEF_TRACER%init_conc /= huge(1.0_r8)) &
         tracers(itrc)%init_conc = DEF_TRACER%init_conc
      IF (DEF_TRACER%precip_default_conc /= huge(1.0_r8)) &
         tracers(itrc)%precip_default_conc = DEF_TRACER%precip_default_conc
      IF (DEF_TRACER%vapor_default_conc /= huge(1.0_r8)) &
         tracers(itrc)%vapor_default_conc = DEF_TRACER%vapor_default_conc
      tracers(itrc)%reactive_decay_rate = DEF_TRACER%reactive_decay_rate
      tracers(itrc)%charge              = DEF_TRACER%charge
      tracers(itrc)%uses_land_water_transport = DEF_TRACER%uses_generic_land_water_transport
      tracers(itrc)%is_nonvolatile      = DEF_TRACER%is_nonvolatile
      tracers(itrc)%uses_reaction       = DEF_TRACER%uses_reaction
   END SUBROUTINE read_tracer_parameter_file

   logical FUNCTION tracer_parameter_group_present (nlfile)
      IMPLICIT NONE
      character(len=*), intent(in) :: nlfile
      character(len=1024) :: line
      integer :: ierr, unit_nml

      tracer_parameter_group_present = .false.
      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form='FORMATTED')
      DO
         read(unit_nml, '(A)', iostat=ierr) line
         IF (ierr /= 0) EXIT
         line = adjustl(line)
         IF (line(1:1) == '!') CYCLE
         IF (index(tracer_lower(line), '&nl_colm_tracer_parameter') == 1 .or. &
             index(tracer_lower(line), '$nl_colm_tracer_parameter') == 1) THEN
            tracer_parameter_group_present = .true.
            EXIT
         ENDIF
      ENDDO
      close(unit_nml)
   END FUNCTION tracer_parameter_group_present

   SUBROUTINE tracer_param_file_for_index (itrc, aliases, file_param, found)
      IMPLICIT NONE
      integer, intent(in) :: itrc
      character(len=*), intent(in) :: aliases
      character(len=*), intent(out) :: file_param
      logical, intent(out) :: found

      integer :: start_pos, end_pos, list_len, colon_pos, positional_index
      character(len=512) :: entry, key, value

      file_param = ''
      found = .false.
      IF (itrc < 1 .or. .not. allocated(tracers) .or. itrc > ntracers) RETURN

      list_len = len_trim(DEF_TRACER_PARAM_FILES)
      IF (list_len <= 0) RETURN
      IF (trim(tracer_lower(DEF_TRACER_PARAM_FILES)) == 'null') RETURN

      positional_index = 0
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = tracer_param_next_entry_end(DEF_TRACER_PARAM_FILES, start_pos, list_len)
         entry = adjustl(trim(DEF_TRACER_PARAM_FILES(start_pos:end_pos)))
         IF (len_trim(entry) > 0) THEN
            colon_pos = index(entry, ':')
            IF (colon_pos > 0) THEN
               key = adjustl(trim(entry(:colon_pos-1)))
               value = adjustl(trim(entry(colon_pos+1:)))
               IF (tracer_param_key_matches(itrc, key, aliases)) THEN
                  file_param = trim(value)
                  found = len_trim(file_param) > 0 .and. trim(tracer_lower(file_param)) /= 'null'
                  RETURN
               ENDIF
            ELSE
               positional_index = positional_index + 1
               IF (positional_index == itrc) THEN
                  file_param = trim(entry)
                  found = len_trim(file_param) > 0 .and. trim(tracer_lower(file_param)) /= 'null'
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         start_pos = end_pos + 2
      ENDDO
   END SUBROUTINE tracer_param_file_for_index

   logical FUNCTION tracer_param_key_matches (itrc, raw_key, aliases)
      IMPLICIT NONE
      integer, intent(in) :: itrc
      character(len=*), intent(in) :: raw_key, aliases

      integer :: start_pos, end_pos, list_len
      character(len=128) :: one_alias

      tracer_param_key_matches = .false.
      IF (itrc < 1 .or. .not. allocated(tracers) .or. itrc > ntracers) RETURN
      IF (tracer_param_equal(raw_key, tracers(itrc)%name)) THEN
         tracer_param_key_matches = .true.
         RETURN
      ENDIF

      list_len = len_trim(aliases)
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = tracer_param_next_entry_end(aliases, start_pos, list_len)
         one_alias = adjustl(trim(aliases(start_pos:end_pos)))
         IF (tracer_param_equal(raw_key, one_alias)) THEN
            tracer_param_key_matches = .true.
            RETURN
         ENDIF
         start_pos = end_pos + 2
      ENDDO
   END FUNCTION tracer_param_key_matches

   logical FUNCTION tracer_param_equal (a, b)
      IMPLICIT NONE
      character(len=*), intent(in) :: a, b
      tracer_param_equal = trim(tracer_lower(a)) == trim(tracer_lower(b))
   END FUNCTION tracer_param_equal

   integer FUNCTION tracer_param_next_entry_end (raw_list, start_pos, list_len)
      IMPLICIT NONE
      character(len=*), intent(in) :: raw_list
      integer, intent(in) :: start_pos, list_len
      integer :: i

      tracer_param_next_entry_end = list_len
      DO i = start_pos, list_len
         IF (raw_list(i:i) == ',' .or. raw_list(i:i) == ';') THEN
            tracer_param_next_entry_end = i - 1
            RETURN
         ENDIF
      ENDDO
   END FUNCTION tracer_param_next_entry_end

   FUNCTION tracer_lower (raw) RESULT(out)
      IMPLICIT NONE
      character(len=*), intent(in) :: raw
      character(len=max(1,len_trim(raw))) :: out
      integer :: i, ia

      out = adjustl(trim(raw))
      DO i = 1, len_trim(out)
         ia = iachar(out(i:i))
         IF (ia >= iachar('A') .and. ia <= iachar('Z')) THEN
            out(i:i) = achar(ia + iachar('a') - iachar('A'))
         ENDIF
      ENDDO
   END FUNCTION tracer_lower

   FUNCTION tracer_upper (raw) RESULT(out)
      IMPLICIT NONE
      character(len=*), intent(in) :: raw
      character(len=max(1,len_trim(raw))) :: out
      integer :: i, ia

      out = adjustl(trim(raw))
      DO i = 1, len_trim(out)
         ia = iachar(out(i:i))
         IF (ia >= iachar('a') .and. ia <= iachar('z')) THEN
            out(i:i) = achar(ia - iachar('a') + iachar('A'))
         ENDIF
      ENDDO
   END FUNCTION tracer_upper

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
   ! DEF_TRACER_INIT_DELTA remains the compatibility fallback.  Non-isotope
   ! parameter files can independently set init_conc, precip_default_conc,
   ! and vapor_default_conc.
   !-------------------------------------------------------------------
   FUNCTION canonical_tracer_category (raw) RESULT(category)
      character(len=*), intent(in) :: raw
      character(len=16) :: category

      category = tracer_lower(raw)
   END FUNCTION canonical_tracer_category

   SUBROUTINE set_tracer_category_defaults (itrc)
      integer, intent(in) :: itrc

      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      tracers(itrc)%unit_kind = 'tracer_per_water'
      tracers(itrc)%charge = 0
      tracers(itrc)%uses_land_water_transport = .true.
      tracers(itrc)%is_nonvolatile = .false.
      tracers(itrc)%uses_reaction = .false.

      SELECT CASE (trim(tracers(itrc)%category))
      CASE ('isotope')
         tracers(itrc)%unit_kind = 'ratio'
      CASE ('conservative')
         tracers(itrc)%is_nonvolatile = .true.
      CASE ('reactive')
         tracers(itrc)%uses_reaction = .true.
      CASE ('particle')
         tracers(itrc)%unit_kind = 'volume_fraction'
         tracers(itrc)%uses_land_water_transport = .false.
      END SELECT
   END SUBROUTINE set_tracer_category_defaults

   SUBROUTINE validate_tracer_descriptor (itrc)
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      integer, intent(in) :: itrc
      logical :: supported_unit

      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      tracers(itrc)%unit_kind = tracer_lower(tracers(itrc)%unit_kind)

      supported_unit = trim(tracers(itrc)%unit_kind) == 'ratio' .or. &
         trim(tracers(itrc)%unit_kind) == 'tracer_per_water' .or. &
         trim(tracers(itrc)%unit_kind) == 'mass_fraction' .or. &
         trim(tracers(itrc)%unit_kind) == 'volume_fraction' .or. &
         trim(tracers(itrc)%unit_kind) == 'species_owned'
      IF (.not. supported_unit) THEN
         CALL tracer_descriptor_error(itrc, 'unit_kind', 'unsupported unit convention')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%mol_weight)) THEN
         CALL tracer_descriptor_error(itrc, 'mol_weight', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%ref_ratio)) THEN
         CALL tracer_descriptor_error(itrc, 'ref_ratio', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%init_delta)) THEN
         CALL tracer_descriptor_error(itrc, 'init_delta', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%init_conc)) THEN
         CALL tracer_descriptor_error(itrc, 'init_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%precip_default_conc)) THEN
         CALL tracer_descriptor_error(itrc, 'precip_default_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%vapor_default_conc)) THEN
         CALL tracer_descriptor_error(itrc, 'vapor_default_conc', 'must be finite')
      ENDIF
      IF (.not. ieee_is_finite(tracers(itrc)%reactive_decay_rate)) THEN
         CALL tracer_descriptor_error(itrc, 'reactive_decay_rate', 'must be finite')
      ENDIF
      IF (tracers(itrc)%mol_weight < 0._r8) THEN
         CALL tracer_descriptor_error(itrc, 'mol_weight', 'must be non-negative')
      ENDIF
      IF (tracers(itrc)%ref_ratio <= 0._r8) THEN
         CALL tracer_descriptor_error(itrc, 'ref_ratio', 'must be positive')
      ENDIF
      IF (tracer_is_isotope(itrc)) THEN
         IF (tracers(itrc)%init_delta < -1000._r8) THEN
            CALL tracer_descriptor_error(itrc, 'init_delta', 'must be at least -1000 permil')
         ENDIF
      ELSE
         IF (tracers(itrc)%init_conc < 0._r8) THEN
            CALL tracer_descriptor_error(itrc, 'init_conc', 'must be non-negative')
         ENDIF
         IF (tracers(itrc)%precip_default_conc < 0._r8) THEN
            CALL tracer_descriptor_error(itrc, 'precip_default_conc', 'must be non-negative')
         ENDIF
         IF (tracers(itrc)%vapor_default_conc < 0._r8) THEN
            CALL tracer_descriptor_error(itrc, 'vapor_default_conc', 'must be non-negative')
         ENDIF
      ENDIF
      IF (tracer_is_isotope(itrc) .and. trim(tracers(itrc)%unit_kind) /= 'ratio') THEN
         CALL tracer_descriptor_error(itrc, 'unit_kind', 'isotope tracers require ratio')
      ENDIF
      IF (tracer_is_conservative(itrc) .and. .not. tracers(itrc)%is_nonvolatile) THEN
         CALL tracer_descriptor_error(itrc, 'is_nonvolatile', &
            'generic conservative transport currently requires a nonvolatile solute')
      ENDIF
      IF (tracers(itrc)%charge /= 0 .and. .not. tracers(itrc)%is_nonvolatile) THEN
         CALL tracer_descriptor_error(itrc, 'charge', &
            'non-zero ionic charge requires a nonvolatile land-water solute')
      ENDIF
      IF ((tracer_is_isotope(itrc) .or. tracer_is_conservative(itrc)) .and. &
          .not. tracers(itrc)%uses_land_water_transport) THEN
         CALL tracer_descriptor_error(itrc, 'uses_generic_land_water_transport', &
            'isotope and conservative tracers require generic land-water transport')
      ENDIF
      IF (tracer_is_particle(itrc) .and. tracers(itrc)%uses_land_water_transport) THEN
         CALL tracer_descriptor_error(itrc, 'uses_generic_land_water_transport', &
            'particle tracers use species-owned suspended and bed pools')
      ENDIF
      IF (tracers(itrc)%is_nonvolatile .and. &
          .not. (tracer_is_conservative(itrc) .or. tracer_is_reactive(itrc))) THEN
         CALL tracer_descriptor_error(itrc, 'is_nonvolatile', &
            'is only valid for conservative or reactive solutes')
      ENDIF
      IF (tracers(itrc)%is_nonvolatile .and. &
          .not. tracers(itrc)%uses_land_water_transport) THEN
         CALL tracer_descriptor_error(itrc, 'is_nonvolatile', &
            'requires generic land-water transport')
      ENDIF
      IF (tracer_is_reactive(itrc) .neqv. tracers(itrc)%uses_reaction) THEN
         CALL tracer_descriptor_error(itrc, 'uses_reaction', &
            'must match category=reactive in the current registry')
      ENDIF
      IF (tracers(itrc)%reactive_decay_rate > 0._r8 .and. .not. tracers(itrc)%uses_reaction) THEN
         CALL tracer_descriptor_error(itrc, 'reactive_decay_rate', 'requires uses_reaction=.true.')
      ENDIF

      IF (.not. tracers(itrc)%uses_land_water_transport .and. &
          .not. tracer_is_particle(itrc) .and. &
          trim(tracers(itrc)%unit_kind) /= 'species_owned') THEN
         CALL tracer_descriptor_error(itrc, 'unit_kind', &
            'species-owned state requires unit_kind=species_owned')
      ENDIF

#ifndef BGC
      IF (tracer_is_reactive(itrc) .and. &
          (trim(tracer_upper(tracers(itrc)%name)) == 'CH4' .or. &
           trim(tracer_upper(tracers(itrc)%name)) == 'METHANE') .and. &
          trim(tracers(itrc)%unit_kind) == 'species_owned') THEN
         CALL tracer_descriptor_error(itrc, 'unit_kind', &
            'species-owned CH4 requires compiling with BGC')
      ENDIF
#endif

   CONTAINS

      SUBROUTINE tracer_descriptor_error (itrc_local, field_name, reason)
         integer, intent(in) :: itrc_local
         character(len=*), intent(in) :: field_name, reason

         IF (p_is_master) THEN
            write(*,'(5A)') 'ERROR tracer_defs_init: tracer "', trim(tracers(itrc_local)%name), &
               '" field ', trim(field_name), ': '//trim(reason)
         ENDIF
         CALL CoLM_stop()
      END SUBROUTINE tracer_descriptor_error

   END SUBROUTINE validate_tracer_descriptor

   logical FUNCTION tracer_is_isotope (itrc)
      integer, intent(in) :: itrc
      tracer_is_isotope = .false.
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_is_isotope = trim(tracers(itrc)%category) == 'isotope'
   END FUNCTION tracer_is_isotope

   logical FUNCTION tracer_is_conservative (itrc)
      integer, intent(in) :: itrc
      tracer_is_conservative = .false.
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_is_conservative = trim(tracers(itrc)%category) == 'conservative'
   END FUNCTION tracer_is_conservative

   logical FUNCTION tracer_is_reactive (itrc)
      integer, intent(in) :: itrc
      tracer_is_reactive = .false.
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_is_reactive = trim(tracers(itrc)%category) == 'reactive'
   END FUNCTION tracer_is_reactive

   logical FUNCTION tracer_is_nonvolatile_solute (itrc)
      ! Nonvolatile conservative or reactive solutes are dissolved material,
      ! not water isotopologues. They move with liquid runoff/infiltration but
      ! atmospheric phase changes remove/add water only; evap/subl/transp
      ! leave solute behind and dew/frost do not import it unless a separate
      ! deposition source is provided.
      integer, intent(in) :: itrc
      tracer_is_nonvolatile_solute = .false.
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_is_nonvolatile_solute = tracer_uses_land_water_transport(itrc) .and. &
         tracers(itrc)%is_nonvolatile
   END FUNCTION tracer_is_nonvolatile_solute

   logical FUNCTION tracer_is_particle (itrc)
      ! Particle tracers (e.g. suspended sediment) carry concentration per
      ! unit water like conservative tracers, but additionally settle and
      ! exchange mass with a non-water bed pool. Diagnostics are concentration
      ! / bed-change, never delta, so they share the non-isotope branches of
      ! tracer_uses_delta_diagnostics / tracer_can_use_fixed_signature.
      integer, intent(in) :: itrc
      tracer_is_particle = .false.
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_is_particle = trim(tracers(itrc)%category) == 'particle'
   END FUNCTION tracer_is_particle

   logical FUNCTION tracer_uses_land_water_transport (itrc)
      ! Generic land-water tracer modules transport scalar signatures tied to
      ! canopy/soil/snow/aquifer water pools. Particle tracers own separate
      ! species-specific state and must not be initialized, transported,
      ! diagnosed, or restarted through those land-water pools.
      integer, intent(in) :: itrc
      tracer_uses_land_water_transport = .false.
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_uses_land_water_transport = .not. tracer_is_particle(itrc) .and. &
         tracers(itrc)%uses_land_water_transport
   END FUNCTION tracer_uses_land_water_transport

   SUBROUTINE tracer_set_land_water_transport (itrc, enabled)
      integer, intent(in) :: itrc
      logical, intent(in) :: enabled

      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracers(itrc)%uses_land_water_transport = enabled
      IF (tracer_is_reactive(itrc)) THEN
         IF (enabled .and. trim(tracers(itrc)%unit_kind) == 'species_owned') THEN
            tracers(itrc)%unit_kind = 'tracer_per_water'
         ELSEIF (.not. enabled) THEN
            tracers(itrc)%unit_kind = 'species_owned'
         ENDIF
      ENDIF
   END SUBROUTINE tracer_set_land_water_transport

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

   FUNCTION tracer_concentration_units (itrc) RESULT(units)
      integer, intent(in) :: itrc
      character(len=32) :: units

      units = 'tracer/water'
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      SELECT CASE (trim(tracers(itrc)%unit_kind))
      CASE ('ratio')
         units = 'R'
      CASE ('mass_fraction')
         units = 'kg/kg water'
      CASE ('volume_fraction')
         units = 'm3/m3 water'
      CASE ('species_owned')
         units = 'species-owned'
      CASE DEFAULT
         units = 'tracer/water'
      END SELECT
   END FUNCTION tracer_concentration_units

   real(r8) FUNCTION tracer_init_water_ratio (itrc)
      integer, intent(in) :: itrc

      tracer_init_water_ratio = 0._r8
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      IF (tracer_is_isotope(itrc)) THEN
         tracer_init_water_ratio = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
      ELSE
         tracer_init_water_ratio = tracers(itrc)%init_conc
      ENDIF
   END FUNCTION tracer_init_water_ratio

   real(r8) FUNCTION tracer_precip_default_ratio (itrc)
      integer, intent(in) :: itrc

      tracer_precip_default_ratio = 0._r8
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      IF (tracer_is_isotope(itrc)) THEN
         tracer_precip_default_ratio = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
      ELSE
         tracer_precip_default_ratio = tracers(itrc)%precip_default_conc
      ENDIF
   END FUNCTION tracer_precip_default_ratio

   real(r8) FUNCTION tracer_vapor_default_ratio (itrc)
      integer, intent(in) :: itrc

      tracer_vapor_default_ratio = 0._r8
      IF (.not. allocated(tracers)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      IF (tracer_is_isotope(itrc)) THEN
         tracer_vapor_default_ratio = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
      ELSE
         tracer_vapor_default_ratio = tracers(itrc)%vapor_default_conc
      ENDIF
   END FUNCTION tracer_vapor_default_ratio

   real(r8) FUNCTION tracer_reactive_decay_fraction (itrc, deltim)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: deltim
      real(r8) :: kdt

      tracer_reactive_decay_fraction = 0._r8
      IF (.not. tracer_is_reactive(itrc)) RETURN
      IF (.not. tracers(itrc)%uses_reaction) RETURN
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
   ! by dropping it. Land and river/lake history both splice tracers(i)%name.
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

   SUBROUTINE warn_csv_count(field_name, ntokens, expected)
      USE MOD_SPMD_Task, only: p_is_master
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: ntokens, expected

      IF (ntokens == expected) RETURN
      IF (p_is_master) THEN
         WRITE(*,'(A,A,A,I0,A,I0,A)') ' WARNING: ', trim(field_name), ' has ', &
            ntokens, ' entries but DEF_TRACER_NUM=', expected, '; defaults/fallbacks will be used.'
      ENDIF
   END SUBROUTINE warn_csv_count

   SUBROUTINE parse_csv_real (csvstr, n, vals, default_val, field_name)
      USE MOD_SPMD_Task, only: p_is_master
      character(len=*), intent(in)  :: csvstr
      integer, intent(in) :: n
      real(r8), intent(out) :: vals(n)
      real(r8), intent(in) :: default_val
      character(len=*), intent(in), optional :: field_name
      character(len=32), allocatable :: tokens(:)
      integer :: ntokens, i, ios, tok_cap
      ! Match tokens capacity to the caller's requested n so a DEF_TRACER_*
      ! CSV with >100 entries no longer gets silently truncated to the
      ! default_val fallback. parse_csv still reports a WARNING if the
      ! comma count exceeds the array.
      tok_cap = max(n, 16)
      allocate(tokens(tok_cap))
      CALL parse_csv(csvstr, tokens, ntokens)
      IF (present(field_name)) CALL warn_csv_count(field_name, ntokens, n)
      DO i = 1, n
         IF (i <= ntokens .and. i <= size(tokens)) THEN
            READ(tokens(i), *, iostat=ios) vals(i)
            IF (ios /= 0) THEN
               vals(i) = default_val
               IF (present(field_name) .and. p_is_master) THEN
                  WRITE(*,'(A,A,A,I0,A,A,A,E12.5)') ' WARNING: ', trim(field_name), &
                     ' entry ', i, ' is not numeric ("', trim(tokens(i)), '"); using default ', &
                     default_val
               ENDIF
            ENDIF
         ELSE
            vals(i) = default_val
         ENDIF
      ENDDO
      deallocate(tokens)
   END SUBROUTINE parse_csv_real

END MODULE MOD_Tracer_Defs
#endif
