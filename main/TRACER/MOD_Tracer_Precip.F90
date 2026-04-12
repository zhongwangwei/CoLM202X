#include <define.h>

MODULE MOD_Tracer_Precip

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, delta_to_R, trc_tiny
   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, a_trc_precip, &
      trc_pg_to_ground

   IMPLICIT NONE

CONTAINS

   !---------------------------------------------------------------
   ! Canopy tracer update after LEAF_INTERCEPTION.
   !
   ! Physical process:
   !   1. Precipitation arrives: forc_rain, forc_snow
   !   2. Part intercepted (qintr): mixes with EXISTING canopy water
   !   3. Part passes through (throughfall): keeps precipitation R
   !   4. If canopy saturates, excess drips: carries MIXED canopy R
   !
   ! Tracer to ground = throughfall*R_precip + drip*R_mixed
   !   (stored in trc_pg_to_ground for use by tracer_soil_water)
   !
   ! Canopy tracer: first mix interception with old, then remove drip.
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
      real(r8) :: intercepted, drip, throughfall
      real(r8) :: trc_mixed, water_mixed, R_mixed
      real(r8) :: trc_throughfall, trc_drip
      real(r8) :: d_ldew

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         R_input = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! Count total precipitation as system input
         a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) &
            + (forc_rain + forc_snow) * R_input * deltim

         ! ---- Rain component ----
         intercepted = qintr_rain * deltim      ! intercepted amount [mm]
         d_ldew = ldew_rain - ldew_rain_old     ! net canopy change [mm]

         ! Drip = intercepted that didn't stay on canopy
         ! d_ldew = intercepted - drip  =>  drip = intercepted - d_ldew
         drip = max(intercepted - d_ldew, 0._r8)

         ! Throughfall = precipitation that was never intercepted
         throughfall = max(forc_rain * deltim - intercepted, 0._r8)

         ! Mix interception with existing canopy water
         water_mixed = ldew_rain_old + intercepted
         trc_mixed = trc_ldew_rain(itrc, ipatch) + intercepted * R_input

         IF (water_mixed > trc_tiny) THEN
            R_mixed = trc_mixed / water_mixed
         ELSE
            R_mixed = R_input
         ENDIF

         ! Canopy tracer: mixed pool minus drip
         trc_ldew_rain(itrc, ipatch) = max(trc_mixed - drip * R_mixed, 0._r8)

         ! Tracer to ground: throughfall (at R_precip) + drip (at R_mixed)
         trc_throughfall = throughfall * R_input
         trc_drip = drip * R_mixed

         ! ---- Snow component (same logic) ----
         intercepted = qintr_snow * deltim
         d_ldew = ldew_snow - ldew_snow_old
         drip = max(intercepted - d_ldew, 0._r8)
         throughfall = max(forc_snow * deltim - intercepted, 0._r8)

         water_mixed = ldew_snow_old + intercepted
         trc_mixed = trc_ldew_snow(itrc, ipatch) + intercepted * R_input

         IF (water_mixed > trc_tiny) THEN
            R_mixed = trc_mixed / water_mixed
         ELSE
            R_mixed = R_input
         ENDIF

         trc_ldew_snow(itrc, ipatch) = max(trc_mixed - drip * R_mixed, 0._r8)

         trc_throughfall = trc_throughfall + throughfall * R_input
         trc_drip = trc_drip + drip * R_mixed

         ! Store total tracer reaching ground (for tracer_soil_water)
         trc_pg_to_ground(itrc, ipatch) = trc_throughfall + trc_drip

      ENDDO

   END SUBROUTINE tracer_precip

END MODULE MOD_Tracer_Precip
