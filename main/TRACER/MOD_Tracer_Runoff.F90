#include <define.h>

MODULE MOD_Tracer_Runoff

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, trc_tiny
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wdsrf, trc_wa, &
      a_trc_rsur, a_trc_rsub, a_trc_rnof

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_runoff (ipatch, deltim, nl_soil, &
      rsur, rsub, rnof, wdsrf, wliq_soisno)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: nl_soil
      real(r8), intent(in) :: rsur, rsub, rnof, wdsrf
      real(r8), intent(in) :: wliq_soisno(1:nl_soil)

      integer  :: itrc, j
      real(r8) :: ratio, trc_flux

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         IF (rsur > trc_tiny) THEN
            IF (wdsrf > trc_tiny) THEN
               ratio = trc_wdsrf(itrc, ipatch) / wdsrf
            ELSE
               ratio = 0._r8
            ENDIF
            trc_flux = rsur * ratio * deltim
            trc_flux = min(trc_flux, max(trc_wdsrf(itrc, ipatch), 0._r8))
            trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) - trc_flux
            a_trc_rsur(itrc, ipatch) = a_trc_rsur(itrc, ipatch) + trc_flux
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_flux
         ENDIF

         IF (rsub > trc_tiny) THEN
            j = nl_soil
            IF (wliq_soisno(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno(j)
               trc_flux = rsub * ratio * deltim
               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               a_trc_rsub(itrc, ipatch) = a_trc_rsub(itrc, ipatch) + trc_flux
               a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE tracer_runoff

END MODULE MOD_Tracer_Runoff
