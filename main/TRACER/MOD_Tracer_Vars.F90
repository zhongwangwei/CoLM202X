#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Vars

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, tracer_is_particle, &
      tracer_uses_land_water_transport

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

   ! Non-steady-state leaf-water isotope memory for transpiration.
   ! These are diagnostic state variables for the isotope fractionation
   ! calculation; hydrology does not carry a prognostic leaf-tissue water pool.
   real(r8), allocatable :: trc_leaf_delta_e(:,:)      ! evaporation-site leaf water delta [permil]
   real(r8), allocatable :: trc_leaf_delta_b(:,:)      ! bulk leaf water delta [permil]
   real(r8), allocatable :: trc_leaf_peclet (:,:)      ! previous Peclet mixing factor [-]
   real(r8), allocatable :: trc_leaf_water_moles(:,:)  ! previous leaf water content [mol/m2]
   real(r8), allocatable :: trc_leaf_iso_storage(:,:)  ! NSS isotope storage anomaly [R*mm]
   logical,  allocatable :: trc_runtime_forced(:)      ! tracer has runtime atmospheric forcing

		   real(r8), allocatable :: a_trc_precip   (:,:)
	   integer, parameter :: TRC_EVAP_KIND_TOTAL       = 0
	   integer, parameter :: TRC_EVAP_KIND_TRANSP      = 1
	   integer, parameter :: TRC_EVAP_KIND_SOILEVAP    = 2
	   integer, parameter :: TRC_EVAP_KIND_CANOPYEVAP  = 3
	   integer, parameter :: TRC_EVAP_KIND_SUBL        = 4
	   integer, parameter :: TRC_EVAP_KIND_WETLAND     = 5
	   real(r8), allocatable :: a_trc_evap     (:,:)
	   real(r8), allocatable :: a_water_evap_gross(:,:)
	   real(r8), allocatable :: a_trc_transp   (:,:)
	   real(r8), allocatable :: a_trc_transp_src(:,:)
	   real(r8), allocatable :: a_water_transp (:,:)
	   real(r8), allocatable :: a_trc_soilevap (:,:)
	   real(r8), allocatable :: a_water_soilevap(:,:)
	   real(r8), allocatable :: a_trc_canopyevap(:,:)
	   real(r8), allocatable :: a_water_canopyevap(:,:)
	   real(r8), allocatable :: a_trc_subl(:,:)
	   real(r8), allocatable :: a_water_subl(:,:)
	   real(r8), allocatable :: a_trc_wetland_evap(:,:)
	   real(r8), allocatable :: a_water_wetland_evap(:,:)
   real(r8), allocatable :: a_trc_rsur     (:,:)
   real(r8), allocatable :: a_trc_rsub     (:,:)
   real(r8), allocatable :: a_trc_rnof     (:,:)
   real(r8), allocatable :: a_trc_qinfl    (:,:)
   real(r8), allocatable :: a_trc_qcharge  (:,:)

   real(r8), allocatable :: trc_storage_beg(:,:)
   real(r8), allocatable :: trc_balance_err(:,:)
	   ! Net reactive source/sink applied during the current land step.
	   ! Positive adds tracer to storage; negative removes tracer.
	   real(r8), allocatable :: trc_reactive_source_step(:,:)
	   ! Explicit numerical source/sink used only to reconcile solver residuals
	   ! that have no physical flux signature. Positive adds tracer to storage.
	   real(r8), allocatable :: trc_numerical_residual_step(:,:)

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
	   real(r8), allocatable :: a_trc_wa_debt_mass(:,:)  ! signed-debt magnitude diagnostics
	   real(r8), allocatable :: a_water_wa_debt   (:)    ! aquifer debt magnitude
   real(r8), allocatable :: a_trc_wdsrf_mass (:,:)
   real(r8), allocatable :: a_water_wdsrf    (:)
   real(r8), allocatable :: a_trc_wetwat_mass(:,:)
   real(r8), allocatable :: a_water_wetwat   (:)
   real(r8), allocatable :: a_trc_scv_mass   (:,:)
   real(r8), allocatable :: a_water_scv      (:)

	   PUBLIC :: allocate_Tracer_Vars, deallocate_Tracer_Vars, flush_Tracer_Acc
	   PUBLIC :: zero_particle_land_tracer_state
	   PUBLIC :: sync_tracer_patch_phase1, sync_tracer_patch_ratio
	   PUBLIC :: tracer_book_evap_loss

CONTAINS

   SUBROUTINE allocate_Tracer_Vars (numpatch, maxsnl, nl_soil)
      IMPLICIT NONE
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      integer :: itrc
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
      allocate(trc_leaf_delta_e(ntracers, numpatch));        trc_leaf_delta_e = 0._r8
      allocate(trc_leaf_delta_b(ntracers, numpatch));        trc_leaf_delta_b = 0._r8
      allocate(trc_leaf_peclet (ntracers, numpatch));        trc_leaf_peclet  = 1._r8
      allocate(trc_leaf_water_moles(ntracers, numpatch));    trc_leaf_water_moles = 0._r8
      allocate(trc_leaf_iso_storage(ntracers, numpatch));    trc_leaf_iso_storage = 0._r8
      allocate(trc_runtime_forced(ntracers));                trc_runtime_forced = .false.
      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         trc_leaf_delta_e(itrc, :) = tracers(itrc)%init_delta
         trc_leaf_delta_b(itrc, :) = tracers(itrc)%init_delta
      ENDDO

	      allocate(a_trc_precip    (ntracers, numpatch));           a_trc_precip    = 0._r8
		      allocate(a_trc_evap      (ntracers, numpatch));           a_trc_evap      = 0._r8
		      allocate(a_water_evap_gross(ntracers, numpatch));         a_water_evap_gross = 0._r8
		      allocate(a_trc_transp    (ntracers, numpatch));           a_trc_transp    = 0._r8
		      allocate(a_trc_transp_src(ntracers, numpatch));           a_trc_transp_src = 0._r8
		      allocate(a_water_transp  (ntracers, numpatch));           a_water_transp  = 0._r8
		      allocate(a_trc_soilevap  (ntracers, numpatch));           a_trc_soilevap  = 0._r8
		      allocate(a_water_soilevap(ntracers, numpatch));           a_water_soilevap = 0._r8
		      allocate(a_trc_canopyevap(ntracers, numpatch));           a_trc_canopyevap = 0._r8
		      allocate(a_water_canopyevap(ntracers, numpatch));         a_water_canopyevap = 0._r8
		      allocate(a_trc_subl      (ntracers, numpatch));           a_trc_subl      = 0._r8
		      allocate(a_water_subl    (ntracers, numpatch));           a_water_subl    = 0._r8
		      allocate(a_trc_wetland_evap(ntracers, numpatch));         a_trc_wetland_evap = 0._r8
		      allocate(a_water_wetland_evap(ntracers, numpatch));       a_water_wetland_evap = 0._r8
      allocate(a_trc_rsur      (ntracers, numpatch));           a_trc_rsur      = 0._r8
      allocate(a_trc_rsub      (ntracers, numpatch));           a_trc_rsub      = 0._r8
      allocate(a_trc_rnof      (ntracers, numpatch));           a_trc_rnof      = 0._r8
      allocate(a_trc_qinfl     (ntracers, numpatch));           a_trc_qinfl     = 0._r8
      allocate(a_trc_qcharge   (ntracers, numpatch));           a_trc_qcharge   = 0._r8

	      allocate(trc_storage_beg (ntracers, numpatch));           trc_storage_beg = 0._r8
	      allocate(trc_balance_err (ntracers, numpatch));           trc_balance_err = 0._r8
	      allocate(trc_reactive_source_step(ntracers, numpatch));   trc_reactive_source_step = 0._r8
	      allocate(trc_numerical_residual_step(ntracers, numpatch)); trc_numerical_residual_step = 0._r8

	      allocate(a_trc_ldew_mass (ntracers, numpatch));           a_trc_ldew_mass = 0._r8
      allocate(a_water_ldew    (numpatch));                     a_water_ldew    = 0._r8
      allocate(a_trc_soil_mass (ntracers, nl_soil, numpatch));  a_trc_soil_mass = 0._r8
      allocate(a_water_soil    (nl_soil, numpatch));             a_water_soil    = 0._r8
      allocate(a_trc_snow_mass (ntracers, abs(maxsnl), numpatch)); a_trc_snow_mass = 0._r8
      allocate(a_water_snow    (abs(maxsnl), numpatch));            a_water_snow    = 0._r8

	      allocate(a_trc_wa_mass    (ntracers, numpatch));          a_trc_wa_mass     = 0._r8
	      allocate(a_water_wa       (numpatch));                    a_water_wa        = 0._r8
	      allocate(a_trc_wa_debt_mass(ntracers, numpatch));         a_trc_wa_debt_mass = 0._r8
	      allocate(a_water_wa_debt   (numpatch));                   a_water_wa_debt    = 0._r8
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
      IF (allocated(trc_leaf_delta_e )) deallocate(trc_leaf_delta_e )
      IF (allocated(trc_leaf_delta_b )) deallocate(trc_leaf_delta_b )
      IF (allocated(trc_leaf_peclet  )) deallocate(trc_leaf_peclet  )
      IF (allocated(trc_leaf_water_moles)) deallocate(trc_leaf_water_moles)
      IF (allocated(trc_leaf_iso_storage)) deallocate(trc_leaf_iso_storage)
      IF (allocated(trc_runtime_forced)) deallocate(trc_runtime_forced)
		      IF (allocated(a_trc_precip   )) deallocate(a_trc_precip   )
		      IF (allocated(a_trc_evap     )) deallocate(a_trc_evap     )
		      IF (allocated(a_water_evap_gross)) deallocate(a_water_evap_gross)
		      IF (allocated(a_trc_transp   )) deallocate(a_trc_transp   )
		      IF (allocated(a_trc_transp_src)) deallocate(a_trc_transp_src)
		      IF (allocated(a_water_transp )) deallocate(a_water_transp )
		      IF (allocated(a_trc_soilevap )) deallocate(a_trc_soilevap )
		      IF (allocated(a_water_soilevap)) deallocate(a_water_soilevap)
		      IF (allocated(a_trc_canopyevap)) deallocate(a_trc_canopyevap)
		      IF (allocated(a_water_canopyevap)) deallocate(a_water_canopyevap)
		      IF (allocated(a_trc_subl     )) deallocate(a_trc_subl     )
		      IF (allocated(a_water_subl   )) deallocate(a_water_subl   )
		      IF (allocated(a_trc_wetland_evap)) deallocate(a_trc_wetland_evap)
		      IF (allocated(a_water_wetland_evap)) deallocate(a_water_wetland_evap)
      IF (allocated(a_trc_rsur     )) deallocate(a_trc_rsur     )
      IF (allocated(a_trc_rsub     )) deallocate(a_trc_rsub     )
      IF (allocated(a_trc_rnof     )) deallocate(a_trc_rnof     )
      IF (allocated(a_trc_qinfl    )) deallocate(a_trc_qinfl    )
      IF (allocated(a_trc_qcharge  )) deallocate(a_trc_qcharge  )
	      IF (allocated(trc_storage_beg)) deallocate(trc_storage_beg)
	      IF (allocated(trc_balance_err)) deallocate(trc_balance_err)
	      IF (allocated(trc_reactive_source_step)) deallocate(trc_reactive_source_step)
	      IF (allocated(trc_numerical_residual_step)) deallocate(trc_numerical_residual_step)
	      IF (allocated(a_trc_ldew_mass)) deallocate(a_trc_ldew_mass)
      IF (allocated(a_water_ldew   )) deallocate(a_water_ldew   )
      IF (allocated(a_trc_soil_mass)) deallocate(a_trc_soil_mass)
      IF (allocated(a_water_soil   )) deallocate(a_water_soil   )
      IF (allocated(a_trc_snow_mass)) deallocate(a_trc_snow_mass)
      IF (allocated(a_water_snow   )) deallocate(a_water_snow   )
	      IF (allocated(a_trc_wa_mass    )) deallocate(a_trc_wa_mass    )
	      IF (allocated(a_water_wa       )) deallocate(a_water_wa       )
	      IF (allocated(a_trc_wa_debt_mass)) deallocate(a_trc_wa_debt_mass)
	      IF (allocated(a_water_wa_debt   )) deallocate(a_water_wa_debt   )
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
		      IF (allocated(a_water_evap_gross)) a_water_evap_gross = 0._r8
		      IF (allocated(a_trc_transp )) a_trc_transp  = 0._r8
		      IF (allocated(a_trc_transp_src)) a_trc_transp_src = 0._r8
		      IF (allocated(a_water_transp)) a_water_transp = 0._r8
		      IF (allocated(a_trc_soilevap)) a_trc_soilevap = 0._r8
		      IF (allocated(a_water_soilevap)) a_water_soilevap = 0._r8
		      IF (allocated(a_trc_canopyevap)) a_trc_canopyevap = 0._r8
		      IF (allocated(a_water_canopyevap)) a_water_canopyevap = 0._r8
		      IF (allocated(a_trc_subl    )) a_trc_subl = 0._r8
		      IF (allocated(a_water_subl  )) a_water_subl = 0._r8
		      IF (allocated(a_trc_wetland_evap)) a_trc_wetland_evap = 0._r8
		      IF (allocated(a_water_wetland_evap)) a_water_wetland_evap = 0._r8
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
	      IF (allocated(a_trc_wa_debt_mass)) a_trc_wa_debt_mass = 0._r8
	      IF (allocated(a_water_wa_debt   )) a_water_wa_debt    = 0._r8
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
	      IF (allocated(trc_reactive_source_step)) trc_reactive_source_step = 0._r8
	      IF (allocated(trc_numerical_residual_step)) trc_numerical_residual_step = 0._r8
	   END SUBROUTINE flush_Tracer_Acc

   SUBROUTINE zero_particle_land_tracer_state ()
      IMPLICIT NONE
      integer :: itrc

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         IF (.not. tracer_is_particle(itrc)) CYCLE

         IF (allocated(trc_ldew_rain  )) trc_ldew_rain  (itrc, :)    = 0._r8
         IF (allocated(trc_ldew_snow  )) trc_ldew_snow  (itrc, :)    = 0._r8
         IF (allocated(trc_wliq_soisno)) trc_wliq_soisno(itrc, :, :) = 0._r8
         IF (allocated(trc_wice_soisno)) trc_wice_soisno(itrc, :, :) = 0._r8
         IF (allocated(trc_wa         )) trc_wa         (itrc, :)    = 0._r8
         IF (allocated(trc_wdsrf      )) trc_wdsrf      (itrc, :)    = 0._r8
         IF (allocated(trc_wetwat     )) trc_wetwat     (itrc, :)    = 0._r8
         IF (allocated(trc_scv        )) trc_scv        (itrc, :)    = 0._r8
         IF (allocated(trc_waterstorage)) trc_waterstorage(itrc, :)  = 0._r8
         IF (allocated(trc_pg_rain_ground)) trc_pg_rain_ground(itrc, :) = 0._r8
         IF (allocated(trc_pg_snow_ground)) trc_pg_snow_ground(itrc, :) = 0._r8
         IF (allocated(trc_rnof_step    )) trc_rnof_step    (itrc, :) = 0._r8
         IF (allocated(trc_sm_carry     )) trc_sm_carry     (itrc, :) = 0._r8
         IF (allocated(trc_leaf_delta_e )) trc_leaf_delta_e (itrc, :) = 0._r8
         IF (allocated(trc_leaf_delta_b )) trc_leaf_delta_b (itrc, :) = 0._r8
         IF (allocated(trc_leaf_peclet  )) trc_leaf_peclet  (itrc, :) = 1._r8
         IF (allocated(trc_leaf_water_moles)) trc_leaf_water_moles(itrc, :) = 0._r8
         IF (allocated(trc_leaf_iso_storage)) trc_leaf_iso_storage(itrc, :) = 0._r8

         IF (allocated(a_trc_precip   )) a_trc_precip   (itrc, :) = 0._r8
         IF (allocated(a_trc_evap     )) a_trc_evap     (itrc, :) = 0._r8
         IF (allocated(a_water_evap_gross)) a_water_evap_gross(itrc, :) = 0._r8
         IF (allocated(a_trc_transp   )) a_trc_transp   (itrc, :) = 0._r8
         IF (allocated(a_trc_transp_src)) a_trc_transp_src(itrc, :) = 0._r8
         IF (allocated(a_water_transp )) a_water_transp (itrc, :) = 0._r8
         IF (allocated(a_trc_soilevap )) a_trc_soilevap (itrc, :) = 0._r8
         IF (allocated(a_water_soilevap)) a_water_soilevap(itrc, :) = 0._r8
         IF (allocated(a_trc_canopyevap)) a_trc_canopyevap(itrc, :) = 0._r8
         IF (allocated(a_water_canopyevap)) a_water_canopyevap(itrc, :) = 0._r8
         IF (allocated(a_trc_subl     )) a_trc_subl     (itrc, :) = 0._r8
         IF (allocated(a_water_subl   )) a_water_subl   (itrc, :) = 0._r8
         IF (allocated(a_trc_wetland_evap)) a_trc_wetland_evap(itrc, :) = 0._r8
         IF (allocated(a_water_wetland_evap)) a_water_wetland_evap(itrc, :) = 0._r8
         IF (allocated(a_trc_rsur     )) a_trc_rsur     (itrc, :) = 0._r8
         IF (allocated(a_trc_rsub     )) a_trc_rsub     (itrc, :) = 0._r8
         IF (allocated(a_trc_rnof     )) a_trc_rnof     (itrc, :) = 0._r8
         IF (allocated(a_trc_qinfl    )) a_trc_qinfl    (itrc, :) = 0._r8
         IF (allocated(a_trc_qcharge  )) a_trc_qcharge  (itrc, :) = 0._r8
         IF (allocated(trc_storage_beg)) trc_storage_beg(itrc, :) = 0._r8
         IF (allocated(trc_balance_err)) trc_balance_err(itrc, :) = 0._r8
         IF (allocated(trc_reactive_source_step)) trc_reactive_source_step(itrc, :) = 0._r8
         IF (allocated(trc_numerical_residual_step)) trc_numerical_residual_step(itrc, :) = 0._r8
         IF (allocated(a_trc_ldew_mass)) a_trc_ldew_mass(itrc, :) = 0._r8
         IF (allocated(a_trc_soil_mass)) a_trc_soil_mass(itrc, :, :) = 0._r8
         IF (allocated(a_trc_snow_mass)) a_trc_snow_mass(itrc, :, :) = 0._r8
         IF (allocated(a_trc_wa_mass    )) a_trc_wa_mass    (itrc, :) = 0._r8
         IF (allocated(a_trc_wa_debt_mass)) a_trc_wa_debt_mass(itrc, :) = 0._r8
         IF (allocated(a_trc_wdsrf_mass )) a_trc_wdsrf_mass (itrc, :) = 0._r8
         IF (allocated(a_trc_wetwat_mass)) a_trc_wetwat_mass(itrc, :) = 0._r8
         IF (allocated(a_trc_scv_mass   )) a_trc_scv_mass   (itrc, :) = 0._r8
      ENDDO
   END SUBROUTINE zero_particle_land_tracer_state

	   SUBROUTINE tracer_book_evap_loss (itrc, ipatch, tracer_mass, water_mass, evap_kind)
	      IMPLICIT NONE
	      integer,  intent(in) :: itrc, ipatch
	      real(r8), intent(in) :: tracer_mass, water_mass
	      integer,  intent(in), optional :: evap_kind
	      integer :: kind

	      IF (ntracers <= 0) RETURN
	      IF (.not. allocated(a_trc_evap)) RETURN
	      kind = TRC_EVAP_KIND_TOTAL
	      IF (present(evap_kind)) kind = evap_kind

	      a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + tracer_mass
	      IF (allocated(a_water_evap_gross)) THEN
	         a_water_evap_gross(itrc, ipatch) = a_water_evap_gross(itrc, ipatch) &
	            + max(water_mass, 0._r8)
	      ENDIF

	      SELECT CASE (kind)
	      CASE (TRC_EVAP_KIND_TRANSP)
	         IF (allocated(a_trc_transp)) a_trc_transp(itrc, ipatch) = &
	            a_trc_transp(itrc, ipatch) + tracer_mass
	         IF (allocated(a_water_transp)) a_water_transp(itrc, ipatch) = &
	            a_water_transp(itrc, ipatch) + max(water_mass, 0._r8)
	      CASE (TRC_EVAP_KIND_SOILEVAP)
	         IF (allocated(a_trc_soilevap)) a_trc_soilevap(itrc, ipatch) = &
	            a_trc_soilevap(itrc, ipatch) + tracer_mass
	         IF (allocated(a_water_soilevap)) a_water_soilevap(itrc, ipatch) = &
	            a_water_soilevap(itrc, ipatch) + max(water_mass, 0._r8)
	      CASE (TRC_EVAP_KIND_CANOPYEVAP)
	         IF (allocated(a_trc_canopyevap)) a_trc_canopyevap(itrc, ipatch) = &
	            a_trc_canopyevap(itrc, ipatch) + tracer_mass
	         IF (allocated(a_water_canopyevap)) a_water_canopyevap(itrc, ipatch) = &
	            a_water_canopyevap(itrc, ipatch) + max(water_mass, 0._r8)
	      CASE (TRC_EVAP_KIND_SUBL)
	         IF (allocated(a_trc_subl)) a_trc_subl(itrc, ipatch) = &
	            a_trc_subl(itrc, ipatch) + tracer_mass
	         IF (allocated(a_water_subl)) a_water_subl(itrc, ipatch) = &
	            a_water_subl(itrc, ipatch) + max(water_mass, 0._r8)
	      CASE (TRC_EVAP_KIND_WETLAND)
	         IF (allocated(a_trc_wetland_evap)) a_trc_wetland_evap(itrc, ipatch) = &
	            a_trc_wetland_evap(itrc, ipatch) + tracer_mass
	         IF (allocated(a_water_wetland_evap)) a_water_wetland_evap(itrc, ipatch) = &
	            a_water_wetland_evap(itrc, ipatch) + max(water_mass, 0._r8)
	      END SELECT
	   END SUBROUTINE tracer_book_evap_loss

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
   ! rebuilding trc_* = water_* * R_init is exact. Reactive tracers are
   ! not fixed-signature: keep their current patch-mixed concentration
   ! and then apply the water-state sync, otherwise reaction history is
   ! erased on simplified branches.
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
   ! trc_scv = scv * R_sync (pre-layer accumulation).
   !---------------------------------------------------------------
   SUBROUTINE sync_tracer_patch_phase1 (ipatch, snl, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno, wa, wdsrf, scv, &
      ldew_rain, ldew_snow)
      USE MOD_Tracer_Defs, only: tracer_init_water_ratio, tracer_can_use_fixed_signature, trc_tiny
      USE MOD_Tracer_Frac, only: tracer_fractionation_active
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wa, wdsrf, scv
      real(r8), intent(in), optional :: ldew_rain, ldew_snow

      integer  :: itrc, j
      real(r8) :: R_init, R_sync, trc_sum, water_sum
      logical  :: fixed_signature_sync

      IF (ntracers <= 0) RETURN
      IF (.not. allocated(trc_wliq_soisno)) RETURN

      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         R_init = tracer_init_water_ratio(itrc)
         fixed_signature_sync = tracer_can_use_fixed_signature(itrc) .and. &
            .not. tracer_fractionation_active(itrc)
         IF (allocated(trc_runtime_forced)) THEN
            fixed_signature_sync = fixed_signature_sync .and. .not. trc_runtime_forced(itrc)
         ENDIF
         R_sync = R_init
         IF (.not. fixed_signature_sync) THEN
            trc_sum = 0._r8
            water_sum = 0._r8
            IF (present(ldew_rain)) THEN
               trc_sum = trc_sum + trc_ldew_rain(itrc, ipatch)
               water_sum = water_sum + max(ldew_rain, 0._r8)
            ENDIF
            IF (present(ldew_snow)) THEN
               trc_sum = trc_sum + trc_ldew_snow(itrc, ipatch)
               water_sum = water_sum + max(ldew_snow, 0._r8)
            ENDIF
            DO j = maxsnl + 1, nl_soil
               trc_sum = trc_sum + trc_wliq_soisno(itrc, j, ipatch) &
                  + trc_wice_soisno(itrc, j, ipatch)
               water_sum = water_sum + max(wliq_soisno(j), 0._r8) &
                  + max(wice_soisno(j), 0._r8)
            ENDDO
            trc_sum = trc_sum + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch)
            water_sum = water_sum + wa + max(wdsrf, 0._r8)
            IF (snl >= 0) THEN
               trc_sum = trc_sum + trc_scv(itrc, ipatch)
               water_sum = water_sum + max(scv, 0._r8)
            ENDIF
            IF (abs(water_sum) > trc_tiny) R_sync = trc_sum / water_sum
         ENDIF

         ! Canopy: rebuild when the caller supplies ldew state
         ! (urban). Callers that omit these args are expected to
         ! have explicitly zeroed trc_ldew_rain/snow before
         ! save_storage (see glacier/lake blocks in CoLMMAIN.F90).
         IF (present(ldew_rain)) THEN
            trc_ldew_rain(itrc, ipatch) = max(ldew_rain, 0._r8) * R_sync
         ENDIF
         IF (present(ldew_snow)) THEN
            trc_ldew_snow(itrc, ipatch) = max(ldew_snow, 0._r8) * R_sync
         ENDIF

         DO j = maxsnl + 1, nl_soil
            trc_wliq_soisno(itrc, j, ipatch) = max(wliq_soisno(j), 0._r8) * R_sync
            trc_wice_soisno(itrc, j, ipatch) = max(wice_soisno(j), 0._r8) * R_sync
         ENDDO

         ! wa carried signed to mirror a possible aquifer-debt state on
         ! boundary patches (lake bed infiltration); wdsrf non-negative
         ! by hydrology construction, so a max(.,0) guard suffices.
         trc_wa   (itrc, ipatch) = wa * R_sync
         trc_wdsrf(itrc, ipatch) = max(wdsrf, 0._r8) * R_sync

         IF (snl < 0) THEN
            trc_scv(itrc, ipatch) = 0._r8
         ELSE
            trc_scv(itrc, ipatch) = max(scv, 0._r8) * R_sync
         ENDIF
      ENDDO
   END SUBROUTINE sync_tracer_patch_phase1

   !---------------------------------------------------------------
   ! Rebuild one tracer's prognostic storage from a caller-owned
   ! mixed-box ratio. This is used by patch branches that do not run
   ! the full tracer hydrology pipeline (glacier / waterbody) but can
   ! no longer assume a fixed R_init signature when runtime forcing or
   ! fractionation is enabled.
   !
   ! `snl` is the accounting lower bound used by tracer_save_storage:
   ! layers below snl+1 are zeroed so hidden snow-layer tracer cannot
   ! leak outside the caller's conservation account.
   !---------------------------------------------------------------
   SUBROUTINE sync_tracer_patch_ratio (itrc, ipatch, snl, maxsnl, nl_soil, &
      wliq_soisno, wice_soisno, wa, wdsrf, scv, R_mix)
      IMPLICIT NONE
      integer,  intent(in) :: itrc, ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wa, wdsrf, scv, R_mix

      integer :: j

      IF (ntracers <= 0) RETURN
      IF (.not. allocated(trc_wliq_soisno)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN

      IF (allocated(trc_ldew_rain)) trc_ldew_rain(itrc, ipatch) = 0._r8
      IF (allocated(trc_ldew_snow)) trc_ldew_snow(itrc, ipatch) = 0._r8

      DO j = maxsnl + 1, nl_soil
         IF (j >= snl + 1) THEN
            trc_wliq_soisno(itrc, j, ipatch) = max(wliq_soisno(j), 0._r8) * R_mix
            trc_wice_soisno(itrc, j, ipatch) = max(wice_soisno(j), 0._r8) * R_mix
         ELSE
            trc_wliq_soisno(itrc, j, ipatch) = 0._r8
            trc_wice_soisno(itrc, j, ipatch) = 0._r8
         ENDIF
      ENDDO

      trc_wa   (itrc, ipatch) = wa * R_mix
      trc_wdsrf(itrc, ipatch) = max(wdsrf, 0._r8) * R_mix

      IF (snl < 0) THEN
         trc_scv(itrc, ipatch) = 0._r8
      ELSE
         trc_scv(itrc, ipatch) = max(scv, 0._r8) * R_mix
      ENDIF
   END SUBROUTINE sync_tracer_patch_ratio

END MODULE MOD_Tracer_Vars
#endif
