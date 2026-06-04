#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Evapo

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, trc_tiny, tracer_uses_land_water_transport, &
      tracer_init_water_ratio, tracer_is_nonvolatile_solute, tracers, trc_delta_sanity_max
   USE MOD_Tracer_Forcing, only: tracer_forcing_vapor_value
   USE MOD_Tracer_Frac, only: tracer_fractionation_active, tracer_alpha_kinetic_craig_gordon, &
      tracer_craig_gordon_evap_ratio, tracer_equilibrium_deposition_ratio, &
      tracer_rayleigh_freezing_loss, tracer_surface_relhum
   USE MOD_Tracer_EvapLimit, only: tracer_evaporative_tracer_loss
	   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, &
	      trc_wliq_soisno, trc_wice_soisno, &
	      trc_numerical_residual_step, &
	      a_trc_precip, tracer_book_evap_loss, &
	      TRC_EVAP_KIND_CANOPYEVAP, TRC_EVAP_KIND_SOILEVAP, TRC_EVAP_KIND_SUBL

   IMPLICIT NONE

CONTAINS

   !---------------------------------------------------------------
   ! Delta-based ET tracer update after THERMAL.
   !
   ! Compare post-THERMAL vs pre-THERMAL water states:
   !   Decrease → evaporation/sublimation/transpiration (output)
   !   Increase → dew/frost (input) or thaw/freeze (internal)
   !---------------------------------------------------------------
   SUBROUTINE tracer_evapo (ipatch, deltim, snl, nl_soil, &
      ldew_rain, ldew_snow, ldew_rain_bef, ldew_snow_bef, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      canopy_smelt_mass_th, canopy_frzc_mass_th, &
      soil_thaw_mass_th, soil_frzc_mass_th, &
      tleaf_frac, t_soisno_frac, forc_q_frac, forc_psrf_frac)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: ldew_rain_bef, ldew_snow_bef
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil)
      ! Explicit canopy snow→rain melt and rain→snow freeze
      ! masses (mm of water-equivalent over deltim) reported by THERMAL /
      ! LeafTemperature. When present, they replace the d_rain/d_snow
      ! sign-pattern heuristic so a simultaneous (qmelt + qevpl) or
      ! (qfrz + qsubl) episode no longer collapses one branch into the
      ! other. The heuristic is preserved as a fallback for callers that
      ! do not (yet) plumb these values — important for hot-start /
      ! third-party drivers that bind tracer_evapo directly.
      real(r8), intent(in), optional :: canopy_smelt_mass_th
      real(r8), intent(in), optional :: canopy_frzc_mass_th
      real(r8), intent(in), optional :: soil_thaw_mass_th(snl+1:nl_soil)
      real(r8), intent(in), optional :: soil_frzc_mass_th(snl+1:nl_soil)
      real(r8), intent(in), optional :: tleaf_frac
      real(r8), intent(in), optional :: t_soisno_frac(snl+1:nl_soil)
      real(r8), intent(in), optional :: forc_q_frac
      real(r8), intent(in), optional :: forc_psrf_frac

      integer  :: itrc, j, lb
      real(r8) :: ratio, trc_flux, R_atm, R_vapor
	      real(r8) :: d_rain, d_snow, d_wliq, d_wice, water_loss, trc_resid
      real(r8) :: thaw_amt, freeze_amt
      ! Post-internal-transfer pool sizes used as the denominator for
      ! evaporation / sublimation. Using *_soisno_bef would mix a post-thaw
      ! tracer numerator with a pre-thaw water denominator.
      real(r8) :: wliq_post_phase, wice_post_phase
      ! Canopy phase-change detection variables.
      ! Mirrors the soil-layer thaw/freeze handling at L94-119: separate
      ! INTERNAL phase transfer (snow↔rain) from EXTERNAL atm I/O
      ! (evap / dew / sublim / frost). Without this split, LeafTemperature's
      ! canopy melt (which moves mass from ldew_snow to ldew_rain) was
      ! misclassified as simultaneous "dew on canopy" + "sublimation from
      ! canopy snow", inflating both a_trc_precip and a_trc_evap and
      ! polluting trc_ldew_rain with R_atm signature instead of the
      ! original snow signature R_canopy_snow.
      real(r8) :: ldew_rain_post_phase, ldew_snow_post_phase
      real(r8) :: d_rain_external, d_snow_external
      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE
         R_atm = tracer_forcing_vapor_value(itrc, ipatch)
         R_vapor = R_atm

         ! --- Canopy rain & snow: detect INTERNAL phase transfer first ---
         d_rain = ldew_rain - ldew_rain_bef
         d_snow = ldew_snow - ldew_snow_bef

         ! Prefer the explicit qmelt / qfrz masses reported
         ! by THERMAL when the caller plumbs them. The d_rain/d_snow
         ! sign-pattern heuristic (legacy fallback below) silently
         ! collapses thaw_amt to min(|d_snow|, d_rain), which under-counts
         ! when canopy melt and rain-pool evaporation fire in the same
         ! step (e.g. tl ~ tfrz with qmelt > 0 AND qevpl > 0 → d_rain =
         ! qmelt - qevpl < qmelt; min picks d_rain instead of qmelt).
         ! Phase-1 invariant hides this because all R = R_atm = R_init,
         ! but Phase 2 (time-varying R_atm) would mis-charge the missing
         ! qmelt fraction to sublimation/dew. Clamp to the matching pool
         ! so caller-side rounding noise can't push the migration past
         ! the actual mass that crossed.
         thaw_amt   = 0._r8   ! snow→rain (canopy melt)
         freeze_amt = 0._r8   ! rain→snow (canopy refreeze)
         IF (present(canopy_smelt_mass_th) .or. present(canopy_frzc_mass_th)) THEN
            IF (present(canopy_smelt_mass_th)) THEN
               thaw_amt = max(canopy_smelt_mass_th, 0._r8)
               thaw_amt = min(thaw_amt, max(ldew_snow_bef, 0._r8))
            ENDIF
            IF (present(canopy_frzc_mass_th)) THEN
               freeze_amt = max(canopy_frzc_mass_th, 0._r8)
               freeze_amt = min(freeze_amt, max(ldew_rain_bef, 0._r8))
            ENDIF
         ELSE
            ! Legacy heuristic — used only when the caller has not
            ! yet been updated to forward THERMAL's explicit phase-change
            ! mass. Misattributes mass when melt+evap (or freeze+frost)
            ! coexist in one step; see comment above.
            IF (d_snow < -trc_tiny .and. d_rain > trc_tiny) THEN
               thaw_amt = min(abs(d_snow), d_rain)
            ENDIF
            IF (d_rain < -trc_tiny .and. d_snow > trc_tiny) THEN
               freeze_amt = min(abs(d_rain), d_snow)
            ENDIF
         ENDIF

         ! --- Internal MELT: trc_ldew_snow → trc_ldew_rain ---
         IF (thaw_amt > trc_tiny) THEN
            IF (ldew_snow_bef > trc_tiny) THEN
               ratio = trc_ldew_snow(itrc, ipatch) / ldew_snow_bef
               trc_flux = min(thaw_amt * ratio, max(trc_ldew_snow(itrc, ipatch), 0._r8))
               trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) - trc_flux
               trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! --- Internal FREEZE: trc_ldew_rain → trc_ldew_snow ---
         IF (freeze_amt > trc_tiny) THEN
            IF (ldew_rain_bef > trc_tiny) THEN
               trc_flux = tracer_rayleigh_freezing_loss(itrc, trc_ldew_rain(itrc, ipatch), &
                  ldew_rain_bef, freeze_amt, canopy_temp())
               trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) - trc_flux
               trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! Reconstruct the canopy pool sizes that participate in the
         ! external phase changes below: rain pool after melt-in /
         ! freeze-out is ldew_rain_bef + thaw - freeze.
         ldew_rain_post_phase = max(ldew_rain_bef + thaw_amt - freeze_amt, 0._r8)
         ldew_snow_post_phase = max(ldew_snow_bef - thaw_amt + freeze_amt, 0._r8)

         ! --- Canopy rain: net change beyond phase = evap/dew ---
         d_rain_external = d_rain - thaw_amt + freeze_amt
         IF (d_rain_external < -trc_tiny) THEN
            ! Net rain loss = canopy evaporation
            IF (ldew_rain_post_phase > trc_tiny) THEN
	               trc_flux = evaporative_tracer_loss(trc_ldew_rain(itrc, ipatch), &
	                  ldew_rain_post_phase, abs(d_rain_external), canopy_temp(), .false.)
	               trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) - trc_flux
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, abs(d_rain_external), &
	                  TRC_EVAP_KIND_CANOPYEVAP)
            ENDIF
         ELSEIF (d_rain_external > trc_tiny) THEN
            ! Net rain gain = dew deposition
            ratio = deposition_ratio_for(canopy_temp(), .false.)
            trc_flux = d_rain_external * ratio
            trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         trc_ldew_rain(itrc, ipatch) = max(trc_ldew_rain(itrc, ipatch), 0._r8)

         ! --- Canopy snow: net change beyond phase = sublim/frost ---
         d_snow_external = d_snow + thaw_amt - freeze_amt
         IF (d_snow_external < -trc_tiny) THEN
            ! Net snow loss = canopy sublimation
            IF (ldew_snow_post_phase > trc_tiny) THEN
	               trc_flux = evaporative_tracer_loss(trc_ldew_snow(itrc, ipatch), &
	                  ldew_snow_post_phase, abs(d_snow_external), canopy_temp(), .true.)
	               trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) - trc_flux
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, abs(d_snow_external), &
	                  TRC_EVAP_KIND_SUBL)
            ENDIF
         ELSEIF (d_snow_external > trc_tiny) THEN
            ! Net snow gain = frost deposition
            ratio = deposition_ratio_for(canopy_temp(), .true.)
            trc_flux = d_snow_external * ratio
            trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         trc_ldew_snow(itrc, ipatch) = max(trc_ldew_snow(itrc, ipatch), 0._r8)

         ! --- Soil+snow layers: combined liquid+ice ---
         DO j = lb, nl_soil
            d_wliq = wliq_soisno(j) - wliq_soisno_bef(j)
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            ! Determine freeze/thaw amounts (internal transfer, not I/O)
            thaw_amt   = 0._r8  ! ice→liquid
            freeze_amt = 0._r8  ! liquid→ice

            IF (present(soil_thaw_mass_th) .or. present(soil_frzc_mass_th)) THEN
               IF (present(soil_thaw_mass_th)) THEN
                  thaw_amt = min(max(soil_thaw_mass_th(j), 0._r8), max(wice_soisno_bef(j), 0._r8))
               ENDIF
               IF (present(soil_frzc_mass_th)) THEN
                  freeze_amt = min(max(soil_frzc_mass_th(j), 0._r8), max(wliq_soisno_bef(j), 0._r8))
               ENDIF
            ELSE
               IF (d_wice < -trc_tiny .and. d_wliq > trc_tiny) THEN
                  thaw_amt = min(abs(d_wice), d_wliq)
               ENDIF
               IF (d_wliq < -trc_tiny .and. d_wice > trc_tiny) THEN
                  freeze_amt = min(abs(d_wliq), d_wice)
               ENDIF
            ENDIF

            ! --- THAW: ice → liquid (internal transfer) ---
            IF (thaw_amt > trc_tiny) THEN
               IF (wice_soisno_bef(j) > trc_tiny) THEN
                  ratio = trc_wice_soisno(itrc, j, ipatch) / wice_soisno_bef(j)
                  trc_flux = min(thaw_amt * ratio, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ENDIF

            ! --- FREEZE: liquid → ice (internal transfer) ---
            IF (freeze_amt > trc_tiny) THEN
               IF (wliq_soisno_bef(j) > trc_tiny) THEN
                  trc_flux = tracer_rayleigh_freezing_loss(itrc, trc_wliq_soisno(itrc, j, ipatch), &
                     wliq_soisno_bef(j), freeze_amt, layer_temp(j))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ENDIF

            ! Reconstruct the water pool sizes that participate in the external
            ! phase changes below: pre-evap liquid = pre-THERMAL liquid + thaw
            ! - freeze, pre-sublimation ice = pre-THERMAL ice - thaw + freeze.
	            wliq_post_phase = max(wliq_soisno_bef(j) + thaw_amt - freeze_amt, 0._r8)
	            wice_post_phase = max(wice_soisno_bef(j) - thaw_amt + freeze_amt, 0._r8)

	            ! THERMAL reports only internal phase changes for snow layers.
	            ! Atmospheric snow exchange (qsdew/qseva/qfros/qsubl) is applied
	            ! later by snowwater and mirrored in tracer_soil_water. Treating
	            ! any remaining snow-layer residual here as dew/frost/evap/subl
	            ! would create tracer fluxes with no matching water-side flux.
	            IF (j < 1) CYCLE

	            ! --- Net liquid change beyond thaw/freeze = evap/dew ---
            ! Net external liquid change = d_wliq - thaw + freeze
            !   (thaw adds liquid internally, freeze removes liquid internally)
	            trc_flux = d_wliq - thaw_amt + freeze_amt
	            IF (trc_flux < -trc_tiny) THEN
	               IF (j == 1) THEN
	                  ! Net exposed-soil liquid loss = evaporation.
	                  IF (wliq_post_phase > trc_tiny) THEN
	                     water_loss = abs(trc_flux)
	                     trc_flux = evaporative_tracer_loss(trc_wliq_soisno(itrc, j, ipatch), &
	                        wliq_post_phase, water_loss, layer_temp(j), .false.)
	                     trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
	                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, water_loss, &
	                        TRC_EVAP_KIND_SOILEVAP)
	                  ENDIF
	               ELSE ! DEEP_STORAGE_RESIDUAL liquid loss: not atmospheric soil evaporation.
	                  IF (wliq_post_phase > trc_tiny) THEN
	                     ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_post_phase
	                  ELSE
	                     ratio = tracer_init_water_ratio(itrc)
	                  ENDIF
	                  trc_resid = min(abs(trc_flux) * max(ratio, 0._r8), &
	                     max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
	                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_resid
	                  IF (allocated(trc_numerical_residual_step)) THEN
	                     trc_numerical_residual_step(itrc, ipatch) = &
	                        trc_numerical_residual_step(itrc, ipatch) - trc_resid
	                  ENDIF
	               ENDIF
            ELSEIF (trc_flux > trc_tiny) THEN
               ! Net liquid gain = dew deposition, only on the exposed soil
               ! surface. Deeper positive residuals are internal storage
               ! redistribution and must not be recoloured as atmospheric dew.
               IF (j == 1) THEN
                  ratio = deposition_ratio_for(layer_temp(j), .false.)
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux * ratio
                  a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux * ratio
               ELSE ! DEEP_STORAGE_RESIDUAL liquid gain: preserve local/fallback signature.
                  IF (wliq_post_phase > trc_tiny) THEN
                     ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_post_phase
                  ELSE
                     ratio = tracer_init_water_ratio(itrc)
                  ENDIF
                  trc_resid = trc_flux * max(ratio, 0._r8)
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_resid
                  IF (allocated(trc_numerical_residual_step)) THEN
                     trc_numerical_residual_step(itrc, ipatch) = &
                        trc_numerical_residual_step(itrc, ipatch) + trc_resid
                  ENDIF
               ENDIF
            ENDIF

            ! --- Net ice change beyond thaw/freeze = sublimation/frost ---
	            trc_flux = d_wice + thaw_amt - freeze_amt
	            IF (trc_flux < -trc_tiny) THEN
	               IF (j == 1) THEN
	                  ! Net exposed-soil ice loss = sublimation.
	                  IF (wice_post_phase > trc_tiny) THEN
	                     water_loss = abs(trc_flux)
	                     trc_flux = evaporative_tracer_loss(trc_wice_soisno(itrc, j, ipatch), &
	                        wice_post_phase, water_loss, layer_temp(j), .true.)
	                     trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
	                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, water_loss, &
	                        TRC_EVAP_KIND_SUBL)
	                  ENDIF
	               ELSE ! DEEP_STORAGE_RESIDUAL ice loss: not atmospheric sublimation.
	                  IF (wice_post_phase > trc_tiny) THEN
	                     ratio = trc_wice_soisno(itrc, j, ipatch) / wice_post_phase
	                  ELSE
	                     ratio = tracer_init_water_ratio(itrc)
	                  ENDIF
	                  trc_resid = min(abs(trc_flux) * max(ratio, 0._r8), &
	                     max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
	                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_resid
	                  IF (allocated(trc_numerical_residual_step)) THEN
	                     trc_numerical_residual_step(itrc, ipatch) = &
	                        trc_numerical_residual_step(itrc, ipatch) - trc_resid
	                  ENDIF
	               ENDIF
            ELSEIF (trc_flux > trc_tiny) THEN
               ! Net ice gain = frost deposition, only on the exposed soil
               ! surface. Deeper positive residuals are internal storage
               ! redistribution and must not be recoloured as atmospheric frost.
               IF (j == 1) THEN
                  ratio = deposition_ratio_for(layer_temp(j), .true.)
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux * ratio
                  a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux * ratio
               ELSE ! DEEP_STORAGE_RESIDUAL ice gain: preserve local/fallback signature.
                  IF (wice_post_phase > trc_tiny) THEN
                     ratio = trc_wice_soisno(itrc, j, ipatch) / wice_post_phase
                  ELSE
                     ratio = tracer_init_water_ratio(itrc)
                  ENDIF
                  trc_resid = trc_flux * max(ratio, 0._r8)
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_resid
                  IF (allocated(trc_numerical_residual_step)) THEN
                     trc_numerical_residual_step(itrc, ipatch) = &
                        trc_numerical_residual_step(itrc, ipatch) + trc_resid
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      CONTAINS

      real(r8) FUNCTION evaporative_tracer_loss (pool_trc, pool_water, water_loss, temp_k, from_ice)
         real(r8), intent(in) :: pool_trc
         real(r8), intent(in) :: pool_water
         real(r8), intent(in) :: water_loss
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         IF (tracer_is_nonvolatile_solute(itrc)) THEN
            evaporative_tracer_loss = 0._r8
         ELSE
            ! r_max ceiling only for fractionating isotope tracers; above the
            ! ceiling, evaporation becomes non-fractionating so the pool R stops
            ! increasing.  Nonvolatile solutes are handled above as zero
            ! evaporative tracer loss.
            evaporative_tracer_loss = tracer_evaporative_tracer_loss(pool_trc, pool_water, &
               water_loss, temp_k, from_ice, evap_ratio_for, trc_tiny, &
               merge(tracers(itrc)%ref_ratio * (1._r8 + trc_delta_sanity_max / 1000._r8), &
                     0._r8, tracer_fractionation_active(itrc)))
         ENDIF
      END FUNCTION evaporative_tracer_loss

      real(r8) FUNCTION canopy_temp ()
         canopy_temp = 273.15_r8
         IF (present(tleaf_frac)) canopy_temp = tleaf_frac
      END FUNCTION canopy_temp

      real(r8) FUNCTION layer_temp (jlay)
         integer, intent(in) :: jlay

         layer_temp = 273.15_r8
         IF (present(t_soisno_frac)) layer_temp = t_soisno_frac(jlay)
      END FUNCTION layer_temp

      real(r8) FUNCTION evap_ratio_for (source_ratio, temp_k, from_ice)
         real(r8), intent(in) :: source_ratio
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8) :: relhum, alpha_k

         evap_ratio_for = source_ratio
         IF (.not. tracer_fractionation_active(itrc)) RETURN
         IF (.not. present(forc_q_frac) .or. .not. present(forc_psrf_frac)) RETURN

         relhum = tracer_surface_relhum(forc_q_frac, forc_psrf_frac, temp_k, from_ice)
         alpha_k = tracer_alpha_kinetic_craig_gordon(itrc, from_ice)
         evap_ratio_for = tracer_craig_gordon_evap_ratio(itrc, source_ratio, R_vapor, &
            temp_k, relhum, alpha_k, from_ice)
         ! Net evaporation/sublimation should not remove heavy isotope at
         ! a ratio higher than the finite source pool ratio in this
         ! one-way storage update. If the Craig-Gordon exchange term would
         ! do that, the residual pool is driven artificially toward
         ! delta=-1000 in a single step.
         evap_ratio_for = min(evap_ratio_for, max(source_ratio, 0._r8))
      END FUNCTION evap_ratio_for

      real(r8) FUNCTION deposition_ratio_for (temp_k, from_ice)
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         IF (tracer_is_nonvolatile_solute(itrc)) THEN
            deposition_ratio_for = 0._r8
            RETURN
         ENDIF
         deposition_ratio_for = R_atm
         IF (.not. tracer_fractionation_active(itrc)) RETURN
         deposition_ratio_for = tracer_equilibrium_deposition_ratio(itrc, R_vapor, temp_k, from_ice)
      END FUNCTION deposition_ratio_for

   END SUBROUTINE tracer_evapo

END MODULE MOD_Tracer_Evapo
#endif
