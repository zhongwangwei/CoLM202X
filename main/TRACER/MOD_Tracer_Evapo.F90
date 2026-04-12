#include <define.h>

MODULE MOD_Tracer_Evapo

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, &
      trc_wliq_soisno, trc_wice_soisno, &
      a_trc_evap, a_trc_trans, a_trc_precip

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

      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         R_atm = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! --- Canopy rain ---
         d_rain = ldew_rain - ldew_rain_bef
         IF (d_rain < -trc_tiny) THEN
            IF (ldew_rain_bef > trc_tiny) THEN
               ratio = trc_ldew_rain(itrc, ipatch) / ldew_rain_bef
               trc_flux = min(abs(d_rain) * ratio, max(trc_ldew_rain(itrc, ipatch), 0._r8))
               trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ELSEIF (d_rain > trc_tiny) THEN
            trc_flux = d_rain * R_atm
            trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         trc_ldew_rain(itrc, ipatch) = max(trc_ldew_rain(itrc, ipatch), 0._r8)

         ! --- Canopy snow ---
         d_snow = ldew_snow - ldew_snow_bef
         IF (d_snow < -trc_tiny) THEN
            IF (ldew_snow_bef > trc_tiny) THEN
               ratio = trc_ldew_snow(itrc, ipatch) / ldew_snow_bef
               trc_flux = min(abs(d_snow) * ratio, max(trc_ldew_snow(itrc, ipatch), 0._r8))
               trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ELSEIF (d_snow > trc_tiny) THEN
            trc_flux = d_snow * R_atm
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

            ! --- Net liquid change beyond thaw/freeze = evap/dew ---
            ! Net external liquid change = d_wliq - thaw + freeze
            !   (thaw adds liquid internally, freeze removes liquid internally)
            trc_flux = d_wliq - thaw_amt + freeze_amt
            IF (trc_flux < -trc_tiny) THEN
               ! Net liquid loss = evaporation/transpiration
               IF (wliq_soisno_bef(j) > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j, ipatch) / max(wliq_soisno_bef(j), trc_tiny)
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
               IF (wice_soisno_bef(j) > trc_tiny) THEN
                  ratio = trc_wice_soisno(itrc, j, ipatch) / max(wice_soisno_bef(j), trc_tiny)
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
