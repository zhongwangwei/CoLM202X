#include <define.h>

MODULE MOD_Tracer_Precip

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, delta_to_R, trc_tiny
   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, a_trc_precip

   IMPLICIT NONE

CONTAINS

   !---------------------------------------------------------------
   ! Delta-based canopy tracer update after LEAF_INTERCEPTION.
   !
   ! Logic:
   !   ldew changed from old to new due to interception/drip/overflow.
   !   Δldew_rain = ldew_rain_new - ldew_rain_old
   !   If Δ > 0: canopy gained water from precipitation → add tracer at precip R
   !   If Δ < 0: canopy lost water (drip/overflow) → remove tracer at canopy ratio
   !
   !   Precipitation input counted = (forc_rain + forc_snow) * R * dt
   !   This is the TOTAL entering the system (interception + throughfall).
   !---------------------------------------------------------------
   SUBROUTINE tracer_precip (ipatch, deltim, &
      forc_rain, forc_snow, qintr, qintr_rain, qintr_snow, &
      pg_rain, pg_snow, ldew_rain, ldew_snow, &
      ldew_rain_old, ldew_snow_old)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: forc_rain, forc_snow
      real(r8), intent(in) :: qintr, qintr_rain, qintr_snow
      real(r8), intent(in) :: pg_rain, pg_snow
      real(r8), intent(in) :: ldew_rain, ldew_snow         ! AFTER interception
      real(r8), intent(in) :: ldew_rain_old, ldew_snow_old ! BEFORE interception

      integer  :: itrc
      real(r8) :: R_input
      real(r8) :: d_rain, d_snow, ratio

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         R_input = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! Count total precipitation as system input
         a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) &
            + (forc_rain + forc_snow) * R_input * deltim

         ! --- Canopy rain water: delta-based ---
         d_rain = ldew_rain - ldew_rain_old
         IF (d_rain > trc_tiny) THEN
            ! Canopy gained rain water (interception): add at precip ratio
            trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + d_rain * R_input
         ELSEIF (d_rain < -trc_tiny) THEN
            ! Canopy lost rain water (drip/overflow): remove at canopy ratio
            IF (ldew_rain_old > trc_tiny) THEN
               ratio = trc_ldew_rain(itrc, ipatch) / ldew_rain_old
               trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + d_rain * ratio
            ENDIF
         ENDIF
         trc_ldew_rain(itrc, ipatch) = max(trc_ldew_rain(itrc, ipatch), 0._r8)

         ! --- Canopy snow water: delta-based ---
         d_snow = ldew_snow - ldew_snow_old
         IF (d_snow > trc_tiny) THEN
            trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) + d_snow * R_input
         ELSEIF (d_snow < -trc_tiny) THEN
            IF (ldew_snow_old > trc_tiny) THEN
               ratio = trc_ldew_snow(itrc, ipatch) / ldew_snow_old
               trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) + d_snow * ratio
            ENDIF
         ENDIF
         trc_ldew_snow(itrc, ipatch) = max(trc_ldew_snow(itrc, ipatch), 0._r8)

      ENDDO

   END SUBROUTINE tracer_precip

END MODULE MOD_Tracer_Precip
