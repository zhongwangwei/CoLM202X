#include <define.h>

MODULE MOD_Tracer_Snow

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, trc_tiny, delta_to_R, tracers
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, trc_pg_snow_ground, &
      trc_scv, trc_wetwat, a_trc_evap

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_newsnow (ipatch, patchtype, snl, snl_old, pg_snow, deltim, &
      scv, scv_bef, wetwat_val, &
      wliq_soisno, wice_soisno, wice_soisno_bef)
      USE MOD_Const_Physical, only: tfrz
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, patchtype, snl, snl_old
      real(r8), intent(in) :: pg_snow                    ! snow throughfall [mm/s]
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: scv                        ! post-newsnow scv [mm]
      real(r8), intent(in) :: scv_bef                    ! pre-newsnow scv [mm]
      real(r8), intent(in) :: wetwat_val                 ! current wetwat [mm]
      ! These arrays are only valid when snl < 0
      real(r8), intent(in), optional :: wliq_soisno(:)
      real(r8), intent(in), optional :: wice_soisno(:)
      real(r8), intent(in), optional :: wice_soisno_bef(:)

      integer  :: itrc, j
      real(r8) :: R_precip, d_wice
      real(r8) :: scv_before_snowfall, trc_lost, ratio_scv

      IF (ntracers <= 0) RETURN

      ! scv before this step's snowfall addition
      scv_before_snowfall = scv_bef

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

            ! Sync trc_scv with water-side scv sinks:
            ! 1. Warm wetland: scv → wetwat, scv=0 (MOD_NewSnow.F90:75-79)
            ! 2. Thin-snow melt: scv reduced (MOD_PhaseChange.F90:237-250)
            ! Handle by comparing trc_scv's implied scv with actual scv.
            ! If scv decreased, proportionally reduce trc_scv.
            IF (scv < trc_tiny) THEN
               ! scv was zeroed (warm wetland or full melt)
               ! Transfer remaining trc_scv to appropriate sink
               IF (patchtype == 2) THEN
                  ! Wetland: scv → wetwat. Transfer tracer too.
                  trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) + trc_scv(itrc, ipatch)
               ELSE
                  ! Melt: scv → sm → runoff/soil. Count as melt output.
                  ! In Phase 1, this is effectively lost to the system
                  ! as snowmelt water enters the gwat path.
                  ! For now, count as evaporation (conservative: tracer leaves system).
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_scv(itrc, ipatch)
               ENDIF
               trc_scv(itrc, ipatch) = 0._r8
            ELSEIF (scv < scv_bef + pg_snow*deltim - trc_tiny) THEN
               ! scv partially reduced (thin-snow melt took some)
               ! scv_after_snowfall = scv_bef + pg_snow*dt (what it would be without melt)
               ! actual scv < scv_after_snowfall means some was removed
               ! Proportionally reduce trc_scv
               ratio_scv = scv / max(scv_bef + pg_snow*deltim, trc_tiny)
               ratio_scv = max(min(ratio_scv, 1._r8), 0._r8)
               trc_lost = trc_scv(itrc, ipatch) * (1._r8 - ratio_scv)
               trc_scv(itrc, ipatch) = trc_scv(itrc, ipatch) * ratio_scv
               ! Melt tracer goes to sm → gwat → soil/runoff
               ! For conservation, count as output (enters water system via sm)
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_lost
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
