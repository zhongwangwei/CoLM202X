#include <define.h>

MODULE MOD_Tracer_Main

   USE MOD_Precision
   USE MOD_Tracer_Defs
   USE MOD_Tracer_Vars
   USE MOD_Tracer_Precip
   USE MOD_Tracer_Evapo
   USE MOD_Tracer_SoilWater
   USE MOD_Tracer_Snow
   USE MOD_Tracer_Runoff
   USE MOD_Tracer_Conservation
   USE MOD_Tracer_Hist
   USE MOD_Tracer_Rest

   IMPLICIT NONE

   PUBLIC :: tracer_init, tracer_final

CONTAINS

   SUBROUTINE tracer_init (numpatch, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat)

      IMPLICIT NONE
      integer,  intent(in) :: numpatch, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain(numpatch)
      real(r8), intent(in) :: ldew_snow(numpatch)
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil, numpatch)
      real(r8), intent(in) :: wa(numpatch)
      real(r8), intent(in) :: wdsrf(numpatch)
      real(r8), intent(in) :: wetwat(numpatch)

      CALL tracer_defs_init()
      IF (ntracers <= 0) RETURN
      CALL allocate_Tracer_Vars(numpatch, maxsnl, nl_soil)
      CALL tracer_init_from_water(numpatch, maxsnl, nl_soil, &
         ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
         wa, wdsrf, wetwat)
   END SUBROUTINE tracer_init

   SUBROUTINE tracer_final ()
      CALL deallocate_Tracer_Vars()
      CALL tracer_defs_final()
   END SUBROUTINE tracer_final

END MODULE MOD_Tracer_Main
