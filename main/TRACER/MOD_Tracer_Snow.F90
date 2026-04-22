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
      real(r8) :: scv_before_snowfall

      IF (ntracers <= 0) RETURN

      ! scv before this step's snowfall addition
      scv_before_snowfall = scv_bef

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

#if (defined CoLMDEBUG)
         ! Invariant: when the pre-newsnow state has layered snow
         ! (snl_old < 0), trc_scv must be zero — pre-layer accumulation
         ! stops once a real snow layer exists, and Cases B/C below will
         ! unconditionally reset trc_scv. A non-zero residual here means
         ! somebody else wrote trc_scv while snl_old<0 (a real bug
         ! upstream, not a simple rounding artefact). Flag it once per
         ! tracer per occurrence so the root cause can be traced before
         ! the reset silently drops the mass.
         IF (snl_old < 0 .and. abs(trc_scv(itrc, ipatch)) > trc_tiny) THEN
            write(*,'(A,I8,A,I3,A,E12.5)') &
               ' WARNING tracer_newsnow: trc_scv residual under layered snow ipatch=', &
               ipatch, ' itrc=', itrc, ' trc_scv_pre=', trc_scv(itrc, ipatch)
         ENDIF
#endif

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
            !
            ! scv sink handling (warm wetland melt, thin-snow melt) is done
            ! AFTER THERMAL in CoLMMAIN, not here. tracer_newsnow only
            ! accumulates and transfers on layer creation.
            !
            ! One exception: warm wetland zeroing scv happens INSIDE newsnow
            ! (before THERMAL), so handle it here.
            IF (patchtype == 2 .and. scv < trc_tiny) THEN
               ! Warm wetland: newsnow transferred scv → wetwat, scv=0
               trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) + trc_scv(itrc, ipatch)
               trc_scv(itrc, ipatch) = 0._r8
            ENDIF
         ENDIF
      ENDDO
   END SUBROUTINE tracer_newsnow

   ! tracer_snow_layer_adj was removed: the prior post-hoc redistribution
   ! (total-tracer / total-water share) homogenised the snow column every
   ! time snowlayerscombine / snowlayersdivide fired, erasing vertical
   ! isotope gradients that the model is otherwise supposed to preserve.
   ! The fix carries trc_wliq / trc_wice / trc_scv through the exact same
   ! per-layer topology as water via the new optional hooks on
   ! MOD_SnowLayersCombineDivide's combine/divide routines.

END MODULE MOD_Tracer_Snow
