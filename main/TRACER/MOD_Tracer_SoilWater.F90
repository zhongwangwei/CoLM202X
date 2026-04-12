#include <define.h>

MODULE MOD_Tracer_SoilWater

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, &
      trc_wa, trc_wdsrf, trc_wetwat, &
      a_trc_qinfl, a_trc_qcharge, trc_pg_to_ground

   IMPLICIT NONE

CONTAINS

   !---------------------------------------------------------------
   ! Delta-based soil/water tracer update after WATER.
   !
   ! WATER modifies wliq_soisno, wice_soisno, wdsrf, wa, wetwat.
   ! We compare post-WATER vs pre-WATER states.
   !
   ! Surface water (wdsrf): receives throughfall, loses to infiltration/runoff.
   !   Use mixed-pool approach: old_trc + throughfall_trc as source pool.
   !
   ! Soil layers: compare Δwliq. Positive = gained (from infiltration/above layer).
   !   Negative = lost (to below layer/groundwater).
   !
   ! Aquifer (wa): compare Δwa directly.
   !---------------------------------------------------------------
   SUBROUTINE tracer_soil_water (ipatch, deltim, snl, nl_soil, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      wa, wa_bef, wdsrf, wdsrf_bef, &
      wetwat, wetwat_bef, pg_rain, pg_snow, &
      rsur, rsub, qinfl, qcharge, &
      use_vsf)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wa, wa_bef, wdsrf, wdsrf_bef
      real(r8), intent(in) :: wetwat, wetwat_bef
      real(r8), intent(in) :: pg_rain, pg_snow
      real(r8), intent(in) :: rsur, rsub, qinfl, qcharge
      logical,  intent(in) :: use_vsf

      integer  :: itrc, j, lb
      real(r8) :: R_precip
      real(r8) :: trc_flux, ratio
      real(r8) :: d_wdsrf, d_wa, d_wliq, d_wetwat
      real(r8) :: trc_throughfall, trc_pool_total, water_pool_total

      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! ============================================================
         ! Surface water (wdsrf): mixed-pool delta approach
         !
         ! The surface water pool receives throughfall and loses to
         ! infiltration and runoff. WATER processes all simultaneously.
         !
         ! Conceptual pool before WATER processing:
         !   water_pool = wdsrf_bef + (pg_rain + pg_snow) * dt
         !   trc_pool   = trc_wdsrf_old + throughfall_trc
         !
         ! After WATER: wdsrf is the residual.
         ! The removed water (infiltration + runoff) carries the pool ratio.
         ! ============================================================

         ! Use pre-computed throughfall tracer (includes correct mixing
         ! of throughfall at R_precip + canopy drip at R_canopy_mixed)
         trc_throughfall = trc_pg_to_ground(itrc, ipatch)
         trc_pool_total  = trc_wdsrf(itrc, ipatch) + trc_throughfall
         water_pool_total = wdsrf_bef + (pg_rain + pg_snow) * deltim

         d_wdsrf = wdsrf - wdsrf_bef

         IF (water_pool_total > trc_tiny) THEN
            ratio = trc_pool_total / water_pool_total
            ! New surface water tracer = residual water * pool ratio
            trc_wdsrf(itrc, ipatch) = max(wdsrf, 0._r8) * ratio
         ELSE
            trc_wdsrf(itrc, ipatch) = 0._r8
         ENDIF

         ! Track infiltration diagnostic
         IF (qinfl > trc_tiny .and. water_pool_total > trc_tiny) THEN
            a_trc_qinfl(itrc, ipatch) = a_trc_qinfl(itrc, ipatch) + qinfl * ratio * deltim
         ENDIF

         ! ============================================================
         ! Soil layers: delta-based
         !
         ! Each layer's wliq changed due to infiltration, inter-layer
         ! flux, root uptake, etc. We already handled ET changes in
         ! tracer_evapo. The remaining changes here are from WATER:
         ! infiltration, redistribution, drainage.
         !
         ! For layers that gained water: source is the layer above
         !   (or surface for top layer). Use source pool ratio.
         ! For layers that lost water: remove at current ratio.
         ! ============================================================

         DO j = max(lb, 1), nl_soil
            d_wliq = wliq_soisno(j) - wliq_soisno_bef(j)

            IF (d_wliq > trc_tiny) THEN
               ! Layer gained water.
               ! Source: from above (infiltration for top layer, percolation for others)
               ! Phase 1: all water has same R, so use R_precip
               ! In general: would use ratio of source layer
               IF (j == max(lb, 1)) THEN
                  ! Top soil layer: gains from surface water pool
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) &
                     + d_wliq * ratio  ! ratio from surface pool (computed above)
               ELSE
                  ! Lower layer: gains from layer above
                  IF (wliq_soisno_bef(j-1) > trc_tiny) THEN
                     trc_flux = d_wliq * (trc_wliq_soisno(itrc, j-1, ipatch) / &
                        max(wliq_soisno(j-1), trc_tiny))
                  ELSE
                     trc_flux = d_wliq * R_precip
                  ENDIF
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               ENDIF

            ELSEIF (d_wliq < -trc_tiny) THEN
               ! Layer lost water: remove at current ratio
               IF (wliq_soisno_bef(j) > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno_bef(j)
                  trc_flux = min(abs(d_wliq) * ratio, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               ENDIF
            ENDIF
         ENDDO

         ! ============================================================
         ! Freeze/thaw during WATER (if not handled in THERMAL)
         ! ============================================================
         DO j = lb, nl_soil
            d_wliq = wliq_soisno(j) - wliq_soisno_bef(j)
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            ! Only handle freeze/thaw within WATER (not re-doing THERMAL)
            IF (d_wice > trc_tiny .and. d_wliq < -trc_tiny) THEN
               ! Freeze during WATER
               trc_flux = min(abs(d_wliq), d_wice)
               IF (wliq_soisno_bef(j) > trc_tiny .and. trc_flux > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j, ipatch) / max(wliq_soisno_bef(j), trc_tiny)
                  trc_flux = min(trc_flux * ratio, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ELSEIF (d_wice < -trc_tiny .and. d_wliq > trc_tiny) THEN
               ! Thaw during WATER
               trc_flux = min(abs(d_wice), d_wliq)
               IF (wice_soisno_bef(j) > trc_tiny .and. trc_flux > trc_tiny) THEN
                  ratio = trc_wice_soisno(itrc, j, ipatch) / max(wice_soisno_bef(j), trc_tiny)
                  trc_flux = min(trc_flux * ratio, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ENDIF
         ENDDO

         ! ============================================================
         ! Aquifer (wa): delta-based
         ! ============================================================
         d_wa = wa - wa_bef
         IF (d_wa > trc_tiny) THEN
            ! Aquifer gained water (recharge from bottom soil layer)
            j = nl_soil
            IF (wliq_soisno_bef(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / max(wliq_soisno(j), trc_tiny)
            ELSE
               ratio = R_precip
            ENDIF
            trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) + d_wa * ratio
            a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) + d_wa * ratio
         ELSEIF (d_wa < -trc_tiny) THEN
            ! Aquifer lost water (discharge to soil)
            IF (wa_bef > trc_tiny) THEN
               ratio = trc_wa(itrc, ipatch) / max(wa_bef, trc_tiny)
               trc_flux = min(abs(d_wa) * ratio, max(trc_wa(itrc, ipatch), 0._r8))
               trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - trc_flux
               a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) - trc_flux
            ENDIF
         ENDIF

         ! ============================================================
         ! Wetland (VSF mode): delta-based
         ! ============================================================
         IF (use_vsf) THEN
            d_wetwat = wetwat - wetwat_bef
            IF (d_wetwat > trc_tiny) THEN
               trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) + d_wetwat * R_precip
            ELSEIF (d_wetwat < -trc_tiny) THEN
               IF (wetwat_bef > trc_tiny) THEN
                  ratio = trc_wetwat(itrc, ipatch) / wetwat_bef
                  trc_flux = min(abs(d_wetwat) * ratio, max(trc_wetwat(itrc, ipatch), 0._r8))
                  trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) - trc_flux
               ENDIF
            ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE tracer_soil_water

END MODULE MOD_Tracer_SoilWater
