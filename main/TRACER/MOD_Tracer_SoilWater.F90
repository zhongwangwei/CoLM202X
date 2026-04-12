#include <define.h>

MODULE MOD_Tracer_SoilWater

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, &
      trc_wa, trc_wdsrf, trc_wetwat, &
      a_trc_qinfl, a_trc_qcharge, a_trc_rsur, a_trc_rsub, a_trc_rnof, &
      trc_pg_to_ground

   IMPLICIT NONE

CONTAINS

   !---------------------------------------------------------------
   ! Flux-driven tracer update after WATER_VSF.
   !
   ! Surface water: mixed-pool approach (throughfall + old wdsrf)
   ! Soil layers: qlayer-driven (exact inter-layer flux tracking)
   ! Aquifer: qcharge-driven
   ! Freeze/thaw: delta-based (same-layer internal transfer)
   !
   ! REQUIRES: DEF_USE_VariablySaturatedFlow = .true.
   !           (WATER_2014 does not output qlayer)
   !---------------------------------------------------------------
   SUBROUTINE tracer_soil_water (ipatch, deltim, snl, nl_soil, &
      qlayer, qcharge, rsur, rsub, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      wa, wa_bef, wdsrf, wdsrf_bef, &
      wetwat, wetwat_bef, pg_rain, pg_snow)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: qlayer(0:nl_soil)     ! inter-layer flux [mm/s]
      real(r8), intent(in) :: qcharge               ! groundwater recharge [mm/s]
      real(r8), intent(in) :: rsur, rsub            ! runoff [mm/s]
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)     ! post-WATER
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)     ! post-WATER
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil) ! pre-WATER (post-THERMAL)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil) ! pre-WATER
      real(r8), intent(in) :: wa, wa_bef, wdsrf, wdsrf_bef
      real(r8), intent(in) :: wetwat, wetwat_bef
      real(r8), intent(in) :: pg_rain, pg_snow

      integer  :: itrc, j, lb
      real(r8) :: R_precip
      real(r8) :: trc_flux, ratio, ratio_src
      real(r8) :: d_wice, d_wliq, d_wa, d_wetwat
      real(r8) :: trc_throughfall, trc_pool_total, water_pool_total
      real(r8) :: ratio_layer(1:nl_soil)  ! pre-WATER tracer ratio per layer

      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! ============================================================
         ! 0. Compute pre-WATER tracer ratios for all soil layers
         !    (after tracer_evapo, ratios should be exact R_input in Phase 1)
         ! ============================================================
         DO j = 1, nl_soil
            IF (wliq_soisno_bef(j) > trc_tiny) THEN
               ratio_layer(j) = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno_bef(j)
            ELSE
               ratio_layer(j) = R_precip
            ENDIF
         ENDDO

         ! ============================================================
         ! 1. Surface water: mixed-pool approach
         !    Pool = old_wdsrf_trc + throughfall_trc
         !    Distributed to: wdsrf_new, rsur, and soil (via qlayer(0))
         ! ============================================================
         trc_throughfall = trc_pg_to_ground(itrc, ipatch)
         trc_pool_total  = trc_wdsrf(itrc, ipatch) + trc_throughfall

         ! Surface pool ratio
         water_pool_total = max(wdsrf, 0._r8) + max(rsur, 0._r8) * deltim &
                          + max(qlayer(0), 0._r8) * deltim
         IF (water_pool_total > trc_tiny .and. trc_pool_total > trc_tiny) THEN
            ratio = trc_pool_total / water_pool_total
         ELSEIF (trc_pool_total <= trc_tiny) THEN
            ratio = 0._r8
         ELSE
            ratio = R_precip
         ENDIF

         trc_wdsrf(itrc, ipatch) = max(wdsrf, 0._r8) * ratio

         ! Surface runoff output
         IF (rsur > trc_tiny) THEN
            a_trc_rsur(itrc, ipatch) = a_trc_rsur(itrc, ipatch) + rsur * ratio * deltim
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + rsur * ratio * deltim
         ENDIF

         ! Infiltration diagnostic
         IF (qlayer(0) > trc_tiny) THEN
            a_trc_qinfl(itrc, ipatch) = a_trc_qinfl(itrc, ipatch) + qlayer(0) * ratio * deltim
         ENDIF

         ! ============================================================
         ! 2. Soil layers: qlayer-driven flux tracking
         !
         !    qlayer(j) > 0: downward from j to j+1, source = layer j
         !    qlayer(j) < 0: upward from j+1 to j, source = layer j+1
         !
         !    Each flux removes from source, adds to destination.
         !    Perfectly conservative: remove = add for every flux.
         ! ============================================================

         ! --- qlayer(0): surface → layer 1 ---
         IF (qlayer(0) > trc_tiny) THEN
            ! Downward: surface to soil. Tracer at surface pool ratio.
            trc_flux = qlayer(0) * ratio * deltim
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) + trc_flux
         ELSEIF (qlayer(0) < -trc_tiny) THEN
            ! Upward: soil to surface. Tracer at layer 1 ratio.
            trc_flux = abs(qlayer(0)) * ratio_layer(1) * deltim
            trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, 1, ipatch), 0._r8))
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_flux
            ! Goes to surface water (add to trc_wdsrf)
            trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) + trc_flux
         ENDIF

         ! --- qlayer(j): layer j ↔ layer j+1 ---
         DO j = 1, nl_soil - 1
            IF (qlayer(j) > trc_tiny) THEN
               ! Downward: j → j+1
               trc_flux = qlayer(j) * ratio_layer(j) * deltim
               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
               trc_wliq_soisno(itrc, j,   ipatch) = trc_wliq_soisno(itrc, j,   ipatch) - trc_flux
               trc_wliq_soisno(itrc, j+1, ipatch) = trc_wliq_soisno(itrc, j+1, ipatch) + trc_flux
            ELSEIF (qlayer(j) < -trc_tiny) THEN
               ! Upward: j+1 → j
               trc_flux = abs(qlayer(j)) * ratio_layer(j+1) * deltim
               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j+1, ipatch), 0._r8))
               trc_wliq_soisno(itrc, j+1, ipatch) = trc_wliq_soisno(itrc, j+1, ipatch) - trc_flux
               trc_wliq_soisno(itrc, j,   ipatch) = trc_wliq_soisno(itrc, j,   ipatch) + trc_flux
            ENDIF
         ENDDO

         ! ============================================================
         ! 3. Groundwater: qcharge-driven
         !    qcharge > 0: soil bottom → aquifer (recharge)
         !    qcharge < 0: aquifer → soil bottom (discharge)
         ! ============================================================
         IF (qcharge > trc_tiny) THEN
            j = nl_soil
            trc_flux = qcharge * ratio_layer(j) * deltim
            trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
            trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) + trc_flux
            a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) + trc_flux
         ELSEIF (qcharge < -trc_tiny) THEN
            IF (wa_bef > trc_tiny) THEN
               ratio_src = trc_wa(itrc, ipatch) / max(wa_bef, trc_tiny)
            ELSE
               ratio_src = R_precip
            ENDIF
            trc_flux = abs(qcharge) * ratio_src * deltim
            trc_flux = min(trc_flux, max(trc_wa(itrc, ipatch), 0._r8))
            trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - trc_flux
            j = nl_soil
            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
            a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) - trc_flux
         ENDIF

         ! ============================================================
         ! 4. Subsurface runoff: from bottom soil layer
         ! ============================================================
         IF (rsub > trc_tiny) THEN
            j = nl_soil
            IF (wliq_soisno_bef(j) > trc_tiny) THEN
               ratio_src = trc_wliq_soisno(itrc, j, ipatch) / &
                  max(wliq_soisno(j), trc_tiny)
               trc_flux = rsub * ratio_src * deltim
               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               a_trc_rsub(itrc, ipatch) = a_trc_rsub(itrc, ipatch) + trc_flux
               a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! ============================================================
         ! 5. Freeze/thaw during WATER: delta-based (same-layer only)
         ! ============================================================
         DO j = lb, nl_soil
            d_wliq = wliq_soisno(j) - wliq_soisno_bef(j)
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            IF (d_wice > trc_tiny .and. d_wliq < -trc_tiny) THEN
               ! Freeze: liquid → ice
               trc_flux = min(abs(d_wliq), d_wice)
               IF (wliq_soisno_bef(j) > trc_tiny .and. trc_flux > trc_tiny) THEN
                  ratio_src = trc_wliq_soisno(itrc, j, ipatch) / max(wliq_soisno_bef(j), trc_tiny)
                  trc_flux = min(trc_flux * ratio_src, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ELSEIF (d_wice < -trc_tiny .and. d_wliq > trc_tiny) THEN
               ! Thaw: ice → liquid
               trc_flux = min(abs(d_wice), d_wliq)
               IF (wice_soisno_bef(j) > trc_tiny .and. trc_flux > trc_tiny) THEN
                  ratio_src = trc_wice_soisno(itrc, j, ipatch) / max(wice_soisno_bef(j), trc_tiny)
                  trc_flux = min(trc_flux * ratio_src, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ENDIF
         ENDDO

         ! ============================================================
         ! 6. Wetland water: delta-based
         ! ============================================================
         d_wetwat = wetwat - wetwat_bef
         IF (d_wetwat > trc_tiny) THEN
            trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) + d_wetwat * R_precip
         ELSEIF (d_wetwat < -trc_tiny) THEN
            IF (wetwat_bef > trc_tiny) THEN
               ratio_src = trc_wetwat(itrc, ipatch) / wetwat_bef
               trc_flux = min(abs(d_wetwat) * ratio_src, max(trc_wetwat(itrc, ipatch), 0._r8))
               trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) - trc_flux
            ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE tracer_soil_water

END MODULE MOD_Tracer_SoilWater
