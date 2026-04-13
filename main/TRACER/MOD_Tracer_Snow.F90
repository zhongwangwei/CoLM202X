#include <define.h>

MODULE MOD_Tracer_Snow

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, trc_tiny, delta_to_R, tracers
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, trc_pg_snow_ground, trc_scv

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_newsnow (ipatch, snl, snl_old, pg_snow, deltim, &
      wliq_soisno, wice_soisno, wice_soisno_bef)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, snl_old
      real(r8), intent(in) :: pg_snow                    ! snow throughfall [mm/s]
      real(r8), intent(in) :: deltim
      ! These arrays are only valid when snl < 0
      real(r8), intent(in), optional :: wliq_soisno(:)
      real(r8), intent(in), optional :: wice_soisno(:)
      real(r8), intent(in), optional :: wice_soisno_bef(:)

      integer  :: itrc, j
      real(r8) :: R_precip, d_wice

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! Every step: accumulate snowfall tracer into trc_scv
         ! (mirrors water model's scv += pg_snow*dt at every step)
         trc_scv(itrc, ipatch) = trc_scv(itrc, ipatch) + trc_pg_snow_ground(itrc, ipatch)

         IF (snl_old == 0 .and. snl < 0) THEN
            ! Case A: First snow layer just created from accumulated scv.
            ! wice_soisno(0) = scv (accumulated over multiple steps).
            ! Transfer trc_scv (accumulated tracer) to the new layer.
            j = snl + 1  ! = 0
            trc_wice_soisno(itrc, j, ipatch) = trc_scv(itrc, ipatch)
            IF (present(wliq_soisno)) THEN
               IF (size(wliq_soisno) > 0) THEN
                  trc_wliq_soisno(itrc, j, ipatch) = wliq_soisno(1) * R_precip
               ENDIF
            ENDIF
            ! Clear trc_scv: all transferred to the snow layer
            trc_scv(itrc, ipatch) = 0._r8

         ELSEIF (snl < 0 .and. snl < snl_old) THEN
            ! Case B: New layer added to existing snow pack (deeper snow).
            ! This is a less common case (snow pack growing a new layer).
            j = snl + 1
            IF (present(wice_soisno)) THEN
               trc_wice_soisno(itrc, j, ipatch) = wice_soisno(1) * R_precip
            ENDIF
            IF (present(wliq_soisno)) THEN
               IF (size(wliq_soisno) > 0) THEN
                  trc_wliq_soisno(itrc, j, ipatch) = wliq_soisno(1) * R_precip
               ENDIF
            ENDIF
            ! trc_scv should be 0 when layers already exist
            trc_scv(itrc, ipatch) = 0._r8

         ELSEIF (snl < 0 .and. snl == snl_old) THEN
            ! Case C: Snow added to existing top layer (no new node).
            j = snl + 1
            IF (present(wice_soisno) .and. present(wice_soisno_bef)) THEN
               d_wice = wice_soisno(1) - wice_soisno_bef(1)
               IF (d_wice > trc_tiny) THEN
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) &
                     + trc_pg_snow_ground(itrc, ipatch)
               ENDIF
            ENDIF
            ! trc_scv should be 0 when layers exist
            trc_scv(itrc, ipatch) = 0._r8

         ELSE
            ! snl == 0: no snow layer yet. trc_scv accumulates (done above).
            ! Nothing else to do.
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
