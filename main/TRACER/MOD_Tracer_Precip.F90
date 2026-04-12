#include <define.h>

MODULE MOD_Tracer_Precip

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, delta_to_R, trc_tiny
   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, a_trc_precip

   IMPLICIT NONE

CONTAINS

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
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: ldew_rain_old, ldew_snow_old

      integer  :: itrc
      real(r8) :: R_input, trc_forc_rain, trc_forc_snow
      real(r8) :: trc_qintr_rain, trc_qintr_snow
      real(r8) :: frac_rain_intr, frac_snow_intr

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         R_input = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         trc_forc_rain = forc_rain * R_input
         trc_forc_snow = forc_snow * R_input

         a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) &
            + (trc_forc_rain + trc_forc_snow) * deltim

         IF (forc_rain > trc_tiny) THEN
            frac_rain_intr = qintr_rain / forc_rain
         ELSE
            frac_rain_intr = 0._r8
         ENDIF
         IF (forc_snow > trc_tiny) THEN
            frac_snow_intr = qintr_snow / forc_snow
         ELSE
            frac_snow_intr = 0._r8
         ENDIF

         frac_rain_intr = min(max(frac_rain_intr, 0._r8), 1._r8)
         frac_snow_intr = min(max(frac_snow_intr, 0._r8), 1._r8)

         trc_qintr_rain = trc_forc_rain * frac_rain_intr
         trc_qintr_snow = trc_forc_snow * frac_snow_intr

         trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) + trc_qintr_rain * deltim
         trc_ldew_snow(itrc, ipatch) = trc_ldew_snow(itrc, ipatch) + trc_qintr_snow * deltim
      ENDDO

   END SUBROUTINE tracer_precip

END MODULE MOD_Tracer_Precip
