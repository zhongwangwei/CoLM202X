#include <define.h>

MODULE MOD_Tracer_Snow

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, trc_tiny, delta_to_R, tracers
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_newsnow (ipatch, snl, snl_old, wliq_soisno, wice_soisno)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, snl_old
      real(r8), intent(in) :: wliq_soisno(snl+1:0)
      real(r8), intent(in) :: wice_soisno(snl+1:0)
      integer  :: itrc, j
      real(r8) :: R_precip

      IF (ntracers <= 0) RETURN
      IF (snl >= 0) RETURN

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
         IF (snl < snl_old .and. snl < 0) THEN
            j = snl + 1
            IF (wice_soisno(j) > trc_tiny) THEN
               trc_wice_soisno(itrc, j, ipatch) = wice_soisno(j) * R_precip
            ENDIF
            IF (wliq_soisno(j) > trc_tiny) THEN
               trc_wliq_soisno(itrc, j, ipatch) = wliq_soisno(j) * R_precip
            ENDIF
         ENDIF
      ENDDO
   END SUBROUTINE tracer_newsnow

   SUBROUTINE tracer_snow_layer_adj (ipatch, maxsnl, nl_soil, &
      snl, snl_old, wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch, maxsnl, nl_soil, snl, snl_old
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_bef(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_bef(maxsnl+1:nl_soil)

      integer  :: itrc, j
      real(r8) :: total_trc_wliq, total_trc_wice, total_wliq, total_wice

      IF (ntracers <= 0) RETURN
      IF (snl >= 0 .and. snl_old >= 0) RETURN

      DO itrc = 1, ntracers
         total_trc_wliq = 0._r8; total_trc_wice = 0._r8
         DO j = maxsnl + 1, 0
            total_trc_wliq = total_trc_wliq + max(trc_wliq_soisno(itrc, j, ipatch), 0._r8)
            total_trc_wice = total_trc_wice + max(trc_wice_soisno(itrc, j, ipatch), 0._r8)
         ENDDO

         total_wliq = 0._r8; total_wice = 0._r8
         DO j = snl + 1, 0
            total_wliq = total_wliq + max(wliq_soisno(j), 0._r8)
            total_wice = total_wice + max(wice_soisno(j), 0._r8)
         ENDDO

         DO j = maxsnl + 1, 0
            IF (j >= snl + 1 .and. j <= 0) THEN
               IF (total_wliq > trc_tiny) THEN
                  trc_wliq_soisno(itrc, j, ipatch) = total_trc_wliq * (wliq_soisno(j) / total_wliq)
               ELSE
                  trc_wliq_soisno(itrc, j, ipatch) = 0._r8
               ENDIF
               IF (total_wice > trc_tiny) THEN
                  trc_wice_soisno(itrc, j, ipatch) = total_trc_wice * (wice_soisno(j) / total_wice)
               ELSE
                  trc_wice_soisno(itrc, j, ipatch) = 0._r8
               ENDIF
            ELSE
               trc_wliq_soisno(itrc, j, ipatch) = 0._r8
               trc_wice_soisno(itrc, j, ipatch) = 0._r8
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE tracer_snow_layer_adj

END MODULE MOD_Tracer_Snow
