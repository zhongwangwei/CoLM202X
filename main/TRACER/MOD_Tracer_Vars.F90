#include <define.h>

MODULE MOD_Tracer_Vars

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers

   IMPLICIT NONE
   SAVE

   real(r8), allocatable :: trc_ldew_rain  (:,:)
   real(r8), allocatable :: trc_ldew_snow  (:,:)
   real(r8), allocatable :: trc_wliq_soisno(:,:,:)
   real(r8), allocatable :: trc_wice_soisno(:,:,:)
   real(r8), allocatable :: trc_wa         (:,:)
   real(r8), allocatable :: trc_wdsrf      (:,:)
   real(r8), allocatable :: trc_wetwat     (:,:)

   ! Tracer mirror of MOD_Vars_TimeVariables::waterstorage. Irrigation
   ! (drip/flood/paddy/sprinkler) is an internal transfer from the
   ! per-patch waterstorage reservoir to soil / canopy; water-side
   ! balance (CoLMMAIN:806) treats it that way, so the atmospheric-input
   ! term (forc_prc+forc_prl) does NOT include irrigation. Without a
   ! matching tracer reservoir the tracer side would have to record
   ! irrigation as external input, creating a systematic offset between
   ! the water and tracer conservation accountings. trc_waterstorage
   ! carries the isotope composition of the reservoir and is debited
   ! by exactly the same volume the water side withdraws each step.
   real(r8), allocatable :: trc_waterstorage(:,:)     ! (ntracers, numpatch)

   ! Snow water equivalent tracer for pre-layer accumulation [kg/m2]
   ! Tracks tracer in scv before a snow layer is created (snl=0).
   ! When the first snow layer forms, trc_scv is transferred to trc_wice_soisno.
   real(r8), allocatable :: trc_scv        (:,:)  ! (ntracers, numpatch)

   ! Throughfall tracer reaching ground [kg/m2 per step]
   ! Rain and snow kept separate: rain enters the surface mixed pool via
   ! tracer_soil_water; snow enters trc_scv via tracer_newsnow (mirrors
   ! water side where gwat only carries pg_rain). A previous "trc_pg_to_ground
   ! = rain+snow" combined var double-counted snow tracer when used in the
   ! surface pool — removed to prevent reintroduction.
   real(r8), allocatable :: trc_pg_rain_ground(:,:)    ! (ntracers, numpatch)
   real(r8), allocatable :: trc_pg_snow_ground(:,:)    ! (ntracers, numpatch)
   real(r8), allocatable :: trc_rnof_step(:,:)         ! (ntracers, numpatch), this land step only

   ! Tracer carried by snowmelt sm from thin-snow (snl==0) scv during THERMAL.
   ! Set in CoLMMAIN after THERMAL based on the actual trc_scv reduction;
   ! consumed by tracer_soil_water when snl==0 .and. sm>0. Replaces the old
   ! "sm * R_precip" injection which used a fixed init delta instead of the
   ! real scv ratio.
   real(r8), allocatable :: trc_sm_carry(:,:)          ! (ntracers, numpatch)

   real(r8), allocatable :: a_trc_precip   (:,:)
   real(r8), allocatable :: a_trc_evap     (:,:)
   real(r8), allocatable :: a_trc_rsur     (:,:)
   real(r8), allocatable :: a_trc_rsub     (:,:)
   real(r8), allocatable :: a_trc_rnof     (:,:)
   real(r8), allocatable :: a_trc_qinfl    (:,:)
   real(r8), allocatable :: a_trc_qcharge  (:,:)

   real(r8), allocatable :: trc_storage_beg(:,:)
   real(r8), allocatable :: trc_balance_err(:,:)

   real(r8), allocatable :: a_trc_ldew_mass (:,:)
   real(r8), allocatable :: a_water_ldew    (:)
   real(r8), allocatable :: a_trc_soil_mass (:,:,:)  ! (ntracers, nl_soil, numpatch)
   real(r8), allocatable :: a_water_soil    (:,:)     ! (nl_soil, numpatch)
   real(r8), allocatable :: a_trc_snow_mass (:,:,:)  ! (ntracers, abs(maxsnl), numpatch) snow layers
   real(r8), allocatable :: a_water_snow    (:,:)     ! (abs(maxsnl), numpatch) snow layers

   ! Surface / aquifer / wetland / thin-snow diagnostic accumulators.
   ! tracer_balance_check already includes these four pools in the conservation
   ! sum, but the history writer only had canopy + soil + snow, so users could
   ! not see the δ values of wa (aquifer), wdsrf (surface water), wetwat
   ! (wetland pool) or scv (pre-layer thin snow). These accumulators let
   ! MOD_Hist emit a per-pool tracer ratio diagnostic alongside the existing
   ! ldew / soisno outputs.
   real(r8), allocatable :: a_trc_wa_mass    (:,:)   ! (ntracers, numpatch)
   real(r8), allocatable :: a_water_wa       (:)     ! (numpatch)
   real(r8), allocatable :: a_trc_wdsrf_mass (:,:)
   real(r8), allocatable :: a_water_wdsrf    (:)
   real(r8), allocatable :: a_trc_wetwat_mass(:,:)
   real(r8), allocatable :: a_water_wetwat   (:)
   real(r8), allocatable :: a_trc_scv_mass   (:,:)
   real(r8), allocatable :: a_water_scv      (:)

   PUBLIC :: allocate_Tracer_Vars, deallocate_Tracer_Vars, flush_Tracer_Acc
   PUBLIC :: sync_tracer_patch_phase1

CONTAINS

   SUBROUTINE allocate_Tracer_Vars (numpatch, maxsnl, nl_soil)
      IMPLICIT NONE
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      IF (ntracers <= 0 .or. numpatch <= 0) RETURN

      allocate(trc_ldew_rain   (ntracers, numpatch));           trc_ldew_rain   = 0._r8
      allocate(trc_ldew_snow   (ntracers, numpatch));           trc_ldew_snow   = 0._r8
      allocate(trc_wliq_soisno (ntracers, maxsnl+1:nl_soil, numpatch)); trc_wliq_soisno = 0._r8
      allocate(trc_wice_soisno (ntracers, maxsnl+1:nl_soil, numpatch)); trc_wice_soisno = 0._r8
      allocate(trc_wa          (ntracers, numpatch));           trc_wa          = 0._r8
      allocate(trc_wdsrf       (ntracers, numpatch));           trc_wdsrf       = 0._r8
      allocate(trc_wetwat      (ntracers, numpatch));           trc_wetwat      = 0._r8
      allocate(trc_scv         (ntracers, numpatch));           trc_scv         = 0._r8
      allocate(trc_waterstorage(ntracers, numpatch));           trc_waterstorage = 0._r8
      allocate(trc_pg_rain_ground(ntracers, numpatch));      trc_pg_rain_ground = 0._r8
      allocate(trc_pg_snow_ground(ntracers, numpatch));      trc_pg_snow_ground = 0._r8
      allocate(trc_rnof_step(ntracers, numpatch));           trc_rnof_step = 0._r8
      allocate(trc_sm_carry (ntracers, numpatch));           trc_sm_carry  = 0._r8

      allocate(a_trc_precip    (ntracers, numpatch));           a_trc_precip    = 0._r8
      allocate(a_trc_evap      (ntracers, numpatch));           a_trc_evap      = 0._r8
      allocate(a_trc_rsur      (ntracers, numpatch));           a_trc_rsur      = 0._r8
      allocate(a_trc_rsub      (ntracers, numpatch));           a_trc_rsub      = 0._r8
      allocate(a_trc_rnof      (ntracers, numpatch));           a_trc_rnof      = 0._r8
      allocate(a_trc_qinfl     (ntracers, numpatch));           a_trc_qinfl     = 0._r8
      allocate(a_trc_qcharge   (ntracers, numpatch));           a_trc_qcharge   = 0._r8

      allocate(trc_storage_beg (ntracers, numpatch));           trc_storage_beg = 0._r8
      allocate(trc_balance_err (ntracers, numpatch));           trc_balance_err = 0._r8

      allocate(a_trc_ldew_mass (ntracers, numpatch));           a_trc_ldew_mass = 0._r8
      allocate(a_water_ldew    (numpatch));                     a_water_ldew    = 0._r8
      allocate(a_trc_soil_mass (ntracers, nl_soil, numpatch));  a_trc_soil_mass = 0._r8
      allocate(a_water_soil    (nl_soil, numpatch));             a_water_soil    = 0._r8
      allocate(a_trc_snow_mass (ntracers, abs(maxsnl), numpatch)); a_trc_snow_mass = 0._r8
      allocate(a_water_snow    (abs(maxsnl), numpatch));            a_water_snow    = 0._r8

      allocate(a_trc_wa_mass    (ntracers, numpatch));          a_trc_wa_mass     = 0._r8
      allocate(a_water_wa       (numpatch));                    a_water_wa        = 0._r8
      allocate(a_trc_wdsrf_mass (ntracers, numpatch));          a_trc_wdsrf_mass  = 0._r8
      allocate(a_water_wdsrf    (numpatch));                    a_water_wdsrf     = 0._r8
      allocate(a_trc_wetwat_mass(ntracers, numpatch));          a_trc_wetwat_mass = 0._r8
      allocate(a_water_wetwat   (numpatch));                    a_water_wetwat    = 0._r8
      allocate(a_trc_scv_mass   (ntracers, numpatch));          a_trc_scv_mass    = 0._r8
      allocate(a_water_scv      (numpatch));                    a_water_scv       = 0._r8
   END SUBROUTINE allocate_Tracer_Vars

   SUBROUTINE deallocate_Tracer_Vars ()
      IMPLICIT NONE
      IF (allocated(trc_ldew_rain  )) deallocate(trc_ldew_rain  )
      IF (allocated(trc_ldew_snow  )) deallocate(trc_ldew_snow  )
      IF (allocated(trc_wliq_soisno)) deallocate(trc_wliq_soisno)
      IF (allocated(trc_wice_soisno)) deallocate(trc_wice_soisno)
      IF (allocated(trc_wa         )) deallocate(trc_wa         )
      IF (allocated(trc_wdsrf      )) deallocate(trc_wdsrf      )
      IF (allocated(trc_wetwat     )) deallocate(trc_wetwat     )
      IF (allocated(trc_scv        )) deallocate(trc_scv        )
      IF (allocated(trc_waterstorage)) deallocate(trc_waterstorage)
      IF (allocated(trc_pg_rain_ground)) deallocate(trc_pg_rain_ground)
      IF (allocated(trc_pg_snow_ground)) deallocate(trc_pg_snow_ground)
      IF (allocated(trc_rnof_step    )) deallocate(trc_rnof_step    )
      IF (allocated(trc_sm_carry     )) deallocate(trc_sm_carry     )
      IF (allocated(a_trc_precip   )) deallocate(a_trc_precip   )
      IF (allocated(a_trc_evap     )) deallocate(a_trc_evap     )
      IF (allocated(a_trc_rsur     )) deallocate(a_trc_rsur     )
      IF (allocated(a_trc_rsub     )) deallocate(a_trc_rsub     )
      IF (allocated(a_trc_rnof     )) deallocate(a_trc_rnof     )
      IF (allocated(a_trc_qinfl    )) deallocate(a_trc_qinfl    )
      IF (allocated(a_trc_qcharge  )) deallocate(a_trc_qcharge  )
      IF (allocated(trc_storage_beg)) deallocate(trc_storage_beg)
      IF (allocated(trc_balance_err)) deallocate(trc_balance_err)
      IF (allocated(a_trc_ldew_mass)) deallocate(a_trc_ldew_mass)
      IF (allocated(a_water_ldew   )) deallocate(a_water_ldew   )
      IF (allocated(a_trc_soil_mass)) deallocate(a_trc_soil_mass)
      IF (allocated(a_water_soil   )) deallocate(a_water_soil   )
      IF (allocated(a_trc_snow_mass)) deallocate(a_trc_snow_mass)
      IF (allocated(a_water_snow   )) deallocate(a_water_snow   )
      IF (allocated(a_trc_wa_mass    )) deallocate(a_trc_wa_mass    )
      IF (allocated(a_water_wa       )) deallocate(a_water_wa       )
      IF (allocated(a_trc_wdsrf_mass )) deallocate(a_trc_wdsrf_mass )
      IF (allocated(a_water_wdsrf    )) deallocate(a_water_wdsrf    )
      IF (allocated(a_trc_wetwat_mass)) deallocate(a_trc_wetwat_mass)
      IF (allocated(a_water_wetwat   )) deallocate(a_water_wetwat   )
      IF (allocated(a_trc_scv_mass   )) deallocate(a_trc_scv_mass   )
      IF (allocated(a_water_scv      )) deallocate(a_water_scv      )
   END SUBROUTINE deallocate_Tracer_Vars

   SUBROUTINE flush_Tracer_Acc ()
      IMPLICIT NONE
      IF (allocated(a_trc_precip )) a_trc_precip  = 0._r8
      IF (allocated(a_trc_evap   )) a_trc_evap    = 0._r8
      IF (allocated(a_trc_rsur   )) a_trc_rsur    = 0._r8
      IF (allocated(a_trc_rsub   )) a_trc_rsub    = 0._r8
      IF (allocated(a_trc_rnof   )) a_trc_rnof    = 0._r8
      IF (allocated(a_trc_qinfl  )) a_trc_qinfl   = 0._r8
      IF (allocated(a_trc_qcharge)) a_trc_qcharge  = 0._r8
      IF (allocated(a_trc_ldew_mass)) a_trc_ldew_mass = 0._r8
      IF (allocated(a_water_ldew   )) a_water_ldew    = 0._r8
      IF (allocated(a_trc_soil_mass)) a_trc_soil_mass = 0._r8
      IF (allocated(a_water_soil   )) a_water_soil    = 0._r8
      IF (allocated(a_trc_snow_mass)) a_trc_snow_mass = 0._r8
      IF (allocated(a_water_snow   )) a_water_snow    = 0._r8
      IF (allocated(a_trc_wa_mass    )) a_trc_wa_mass     = 0._r8
      IF (allocated(a_water_wa       )) a_water_wa        = 0._r8
      IF (allocated(a_trc_wdsrf_mass )) a_trc_wdsrf_mass  = 0._r8
      IF (allocated(a_water_wdsrf    )) a_water_wdsrf     = 0._r8
      IF (allocated(a_trc_wetwat_mass)) a_trc_wetwat_mass = 0._r8
      IF (allocated(a_water_wetwat   )) a_water_wetwat    = 0._r8
      IF (allocated(a_trc_scv_mass   )) a_trc_scv_mass    = 0._r8
      IF (allocated(a_water_scv      )) a_water_scv       = 0._r8
      IF (allocated(trc_pg_rain_ground)) trc_pg_rain_ground = 0._r8
      IF (allocated(trc_pg_snow_ground)) trc_pg_snow_ground = 0._r8
      IF (allocated(trc_rnof_step    )) trc_rnof_step     = 0._r8
      IF (allocated(trc_sm_carry     )) trc_sm_carry      = 0._r8
   END SUBROUTINE flush_Tracer_Acc

   !---------------------------------------------------------------
   ! Phase-1 re-synchronisation of a patch's prognostic tracer pools
   ! to the current water state. Used by branches that do NOT run the
   ! full tracer_precip / tracer_evapo / tracer_soil_water pipeline
   ! (glacier, land-water-bodies, urban-model), where prognostic
   ! pools would otherwise stay at their initial values while the
   ! water pools evolve freely — yielding arbitrarily large
   ! `trc / water` ratios (the `δ_wa ≈ +1039‰` class of artefact
   ! observed in early runs).
   !
   ! Under Phase 1 (constant atmospheric R = R_init, no fractionation)
   ! the pool ratio must be R_init at every step by mass balance, so
   ! rebuilding trc_* = water_* * R_init is exact. Phase 2 (time-
   ! varying forcing R, fractionation) must replace this with a full
   ! mixed-box update per pool.
   !
   ! Optional `ldew_rain` / `ldew_snow`: when present, rebuild the
   ! canopy tracer from the provided values (urban patches have real
   ! canopy state). Glacier / waterbody callers omit them and rely on
   ! their caller-side "purge foreign pools" step to zero trc_ldew
   ! before save_storage.
   !
   ! This routine intentionally does NOT touch trc_wetwat or
   ! trc_waterstorage — those are class-specific and the caller is
   ! responsible for either purging them (glacier/lake/urban) or
   ! leaving them for the normal soil/wetland pipeline (soil-ground).
   !
   ! `snl < 0` means layered snow exists → keep trc_scv=0 (the snow
   ! tracer lives in trc_wice/wliq of the snow layers); otherwise
   ! trc_scv = scv * R_init (pre-layer accumulation).
   !---------------------------------------------------------------
   SUBROUTINE sync_tracer_patch_phase1 (ipatch, snl, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno, wa, wdsrf, scv, &
      ldew_rain, ldew_snow)
      USE MOD_Tracer_Defs, only: tracers, delta_to_R
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wa, wdsrf, scv
      real(r8), intent(in), optional :: ldew_rain, ldew_snow

      integer  :: itrc, j
      real(r8) :: R_init

      IF (ntracers <= 0) RETURN
      IF (.not. allocated(trc_wliq_soisno)) RETURN

      DO itrc = 1, ntracers
         R_init = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! Canopy: rebuild when the caller supplies ldew state
         ! (urban). Callers that omit these args are expected to
         ! have explicitly zeroed trc_ldew_rain/snow before
         ! save_storage (see glacier/lake blocks in CoLMMAIN.F90).
         IF (present(ldew_rain)) THEN
            trc_ldew_rain(itrc, ipatch) = max(ldew_rain, 0._r8) * R_init
         ENDIF
         IF (present(ldew_snow)) THEN
            trc_ldew_snow(itrc, ipatch) = max(ldew_snow, 0._r8) * R_init
         ENDIF

         DO j = maxsnl + 1, nl_soil
            trc_wliq_soisno(itrc, j, ipatch) = max(wliq_soisno(j), 0._r8) * R_init
            trc_wice_soisno(itrc, j, ipatch) = max(wice_soisno(j), 0._r8) * R_init
         ENDDO

         ! wa carried signed to mirror a possible aquifer-debt state on
         ! boundary patches (lake bed infiltration); wdsrf non-negative
         ! by hydrology construction, so a max(.,0) guard suffices.
         trc_wa   (itrc, ipatch) = wa * R_init
         trc_wdsrf(itrc, ipatch) = max(wdsrf, 0._r8) * R_init

         IF (snl < 0) THEN
            trc_scv(itrc, ipatch) = 0._r8
         ELSE
            trc_scv(itrc, ipatch) = max(scv, 0._r8) * R_init
         ENDIF
      ENDDO
   END SUBROUTINE sync_tracer_patch_phase1

END MODULE MOD_Tracer_Vars
