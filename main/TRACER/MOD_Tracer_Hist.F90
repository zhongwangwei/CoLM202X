#include <define.h>

MODULE MOD_Tracer_Hist

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, mass_to_delta, trc_tiny
   USE MOD_Tracer_Vars

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_hist_accumulate (ipatch, nl_soil, ldew_rain, ldew_snow, &
      wliq_soisno, wice_soisno)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, nl_soil
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: wliq_soisno(1:nl_soil)
      real(r8), intent(in) :: wice_soisno(1:nl_soil)

      integer :: itrc, j

      IF (ntracers <= 0) RETURN

      ! Canopy water
      a_water_ldew(ipatch) = a_water_ldew(ipatch) + (ldew_rain + ldew_snow)
      DO itrc = 1, ntracers
         a_trc_ldew_mass(itrc, ipatch) = a_trc_ldew_mass(itrc, ipatch) &
            + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
      ENDDO

      ! Soil layers only (1:nl_soil). Snow layer history deferred —
      ! extending to snl+1:nl_soil caused MPI buffer mismatch in history output.
      DO j = 1, nl_soil
         a_water_soil(j, ipatch) = a_water_soil(j, ipatch) + wliq_soisno(j) + wice_soisno(j)
         DO itrc = 1, ntracers
            a_trc_soil_mass(itrc, j, ipatch) = a_trc_soil_mass(itrc, j, ipatch) &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
      ENDDO

      trc_hist_nac = trc_hist_nac + 1
   END SUBROUTINE tracer_hist_accumulate

END MODULE MOD_Tracer_Hist
