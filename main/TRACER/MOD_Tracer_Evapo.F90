#include <define.h>

MODULE MOD_Tracer_Evapo

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, &
      trc_wliq_soisno, trc_wice_soisno, &
      a_trc_evap, a_trc_precip

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
      wliq_soisno_bef, wice_soisno_bef)

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

      integer  :: itrc, j, lb
      real(r8) :: ratio, trc_flux, R_atm
      real(r8) :: d_rain, d_snow, d_wliq, d_wice
      real(r8) :: thaw_amt, freeze_amt
      ! H1: post-internal-transfer pool sizes used as the denominator for
      ! evaporation / sublimation. Using *_soisno_bef would mix a post-thaw
      ! tracer numerator with a pre-thaw water denominator.
      real(r8) :: wliq_post_phase, wice_post_phase
      ! Audit fix TR-1: canopy phase-change detection variables.
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
         R_atm = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! --- Canopy rain & snow: detect INTERNAL phase transfer first ---
         d_rain = ldew_rain - ldew_rain_bef
         d_snow = ldew_snow - ldew_snow_bef

         ! Audit fix TR-1: detect canopy melt (snow↓ + rain↑) and freeze
         ! (rain↓ + snow↑) using the same heuristic as the soil-layer
         ! block below (L94-99). Mass moved internally between the two
         ! canopy pools must NOT be charged to a_trc_evap or a_trc_precip.
         thaw_amt   = 0._r8   ! snow→rain (canopy melt)
         freeze_amt = 0._r8   ! rain→snow (canopy refreeze)
         IF (d_snow < -trc_tiny .and. d_rain > trc_tiny) THEN
            thaw_amt = min(abs(d_snow), d_rain)
         ENDIF
         IF (d_rain < -trc_tiny .and. d_snow > trc_tiny) THEN
            freeze_amt = min(abs(d_rain), d_snow)
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
               ratio = trc_ldew_rain(itrc, ipatch) / ldew_rain_bef
               trc_flux = min(freeze_amt * ratio, max(trc_ldew_rain(itrc, ipatch), 0._r8))
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
               ratio = trc_ldew_rain(itrc, ipatch) / ldew_rain_post_phase
               trc_flux = min(abs(d_rain_external) * ratio, max(trc_ldew_rain(itrc, ipatch), 0._r8))
               trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ELSEIF (d_rain_external > trc_tiny) THEN
            ! Net rain gain = dew deposition
            trc_flux = d_rain_external * R_atm
            trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         trc_ldew_rain(itrc, ipatch) = max(trc_ldew_rain(itrc, ipatch), 0._r8)

         ! --- Canopy snow: net change beyond phase = sublim/frost ---
         d_snow_external = d_snow + thaw_amt - freeze_amt
         IF (d_snow_external < -trc_tiny) THEN
            ! Net snow loss = canopy sublimation
            IF (ldew_snow_post_phase > trc_tiny) THEN
               ratio = trc_ldew_snow(itrc, ipatch) / ldew_snow_post_phase
               trc_flux = min(abs(d_snow_external) * ratio, max(trc_ldew_snow(itrc, ipatch), 0._r8))
               trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ELSEIF (d_snow_external > trc_tiny) THEN
            ! Net snow gain = frost deposition
            trc_flux = d_snow_external * R_atm
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

            IF (d_wice < -trc_tiny .and. d_wliq > trc_tiny) THEN
               thaw_amt = min(abs(d_wice), d_wliq)
            ENDIF
            IF (d_wliq < -trc_tiny .and. d_wice > trc_tiny) THEN
               freeze_amt = min(abs(d_wliq), d_wice)
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
                  ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno_bef(j)
                  trc_flux = min(freeze_amt * ratio, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ENDIF

            ! Reconstruct the water pool sizes that participate in the external
            ! phase changes below: pre-evap liquid = pre-THERMAL liquid + thaw
            ! - freeze, pre-sublimation ice = pre-THERMAL ice - thaw + freeze.
            wliq_post_phase = max(wliq_soisno_bef(j) + thaw_amt - freeze_amt, 0._r8)
            wice_post_phase = max(wice_soisno_bef(j) - thaw_amt + freeze_amt, 0._r8)

            ! --- Net liquid change beyond thaw/freeze = evap/dew ---
            ! Net external liquid change = d_wliq - thaw + freeze
            !   (thaw adds liquid internally, freeze removes liquid internally)
            trc_flux = d_wliq - thaw_amt + freeze_amt
            IF (trc_flux < -trc_tiny) THEN
               ! Net liquid loss = evaporation/transpiration
               IF (wliq_post_phase > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_post_phase
                  trc_flux = min(abs(trc_flux) * ratio, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
               ENDIF
            ELSEIF (trc_flux > trc_tiny) THEN
               ! Net liquid gain = dew deposition
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux * R_atm
               a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux * R_atm
            ENDIF

            ! --- Net ice change beyond thaw/freeze = sublimation/frost ---
            trc_flux = d_wice + thaw_amt - freeze_amt
            IF (trc_flux < -trc_tiny) THEN
               ! Net ice loss = sublimation
               IF (wice_post_phase > trc_tiny) THEN
                  ratio = trc_wice_soisno(itrc, j, ipatch) / wice_post_phase
                  trc_flux = min(abs(trc_flux) * ratio, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
               ENDIF
            ELSEIF (trc_flux > trc_tiny) THEN
               ! Net ice gain = frost deposition
               trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux * R_atm
               a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux * R_atm
            ENDIF

         ENDDO
      ENDDO

   END SUBROUTINE tracer_evapo

END MODULE MOD_Tracer_Evapo
