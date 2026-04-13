#include <define.h>

MODULE MOD_Tracer_Hist

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, mass_to_delta, trc_tiny
   USE MOD_Tracer_Vars

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_hist_accumulate (ipatch, snl, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)

      integer :: itrc, j, jsnow

      IF (ntracers <= 0) RETURN

      ! Canopy water
      a_water_ldew(ipatch) = a_water_ldew(ipatch) + (ldew_rain + ldew_snow)
      DO itrc = 1, ntracers
         a_trc_ldew_mass(itrc, ipatch) = a_trc_ldew_mass(itrc, ipatch) &
            + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
      ENDDO

      ! Soil layers (1:nl_soil) — fixed dimension, safe for MPI gather
      DO j = 1, nl_soil
         a_water_soil(j, ipatch) = a_water_soil(j, ipatch) + wliq_soisno(j) + wice_soisno(j)
         DO itrc = 1, ntracers
            a_trc_soil_mass(itrc, j, ipatch) = a_trc_soil_mass(itrc, j, ipatch) &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
      ENDDO

      ! Snow layers — separate accumulator with fixed dimension abs(maxsnl).
      ! Snow layer index j ranges from snl+1 to 0 (negative indices).
      ! Map to positive index: jsnow = j - maxsnl  (1-based, from top to bottom)
      ! e.g., maxsnl=-5: j=-4 → jsnow=1, j=-3 → jsnow=2, ..., j=0 → jsnow=5
      IF (snl < 0) THEN
         DO j = snl + 1, 0
            jsnow = j - maxsnl  ! positive index (1 to abs(maxsnl))
            a_water_snow(jsnow, ipatch) = a_water_snow(jsnow, ipatch) &
               + wliq_soisno(j) + wice_soisno(j)
            DO itrc = 1, ntracers
               a_trc_snow_mass(itrc, jsnow, ipatch) = a_trc_snow_mass(itrc, jsnow, ipatch) &
                  + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
            ENDDO
         ENDDO
      ENDIF

      trc_hist_nac = trc_hist_nac + 1
   END SUBROUTINE tracer_hist_accumulate

END MODULE MOD_Tracer_Hist
