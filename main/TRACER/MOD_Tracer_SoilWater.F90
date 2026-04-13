#include <define.h>

MODULE MOD_Tracer_SoilWater

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, &
      trc_wa, trc_wdsrf, trc_wetwat, &
      a_trc_precip, a_trc_evap, &
      a_trc_qinfl, a_trc_qcharge, a_trc_rsur, a_trc_rsub, a_trc_rnof, &
      trc_pg_to_ground

   IMPLICIT NONE

CONTAINS

   !---------------------------------------------------------------
   ! Flux-driven tracer update after WATER_VSF.
   !
   ! Surface water: mixed-pool approach (throughfall + old wdsrf + soil upflow)
   ! Soil layers: qlayer-driven (exact inter-layer flux tracking)
   ! Aquifer: qcharge-driven
   ! Subsurface runoff: pre-WATER ratio
   ! Freeze/thaw: delta-based (same-layer internal transfer)
   !
   ! REQUIRES: DEF_USE_VariablySaturatedFlow = .true.
   !---------------------------------------------------------------
   SUBROUTINE tracer_soil_water (ipatch, deltim, snl, nl_soil, &
      qlayer, qcharge, rsur, rsub, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, sm, fsno, &
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
      real(r8), intent(in) :: qseva_soil            ! soil evaporation in WATER [mm/s]
      real(r8), intent(in) :: qsdew_soil            ! soil dew in WATER [mm/s]
      real(r8), intent(in) :: qsubl_soil            ! soil sublimation in WATER [mm/s]
      real(r8), intent(in) :: qfros_soil            ! soil frost in WATER [mm/s]
      real(r8), intent(in) :: sm                    ! snowmelt [mm/s]
      real(r8), intent(in) :: fsno                  ! snow fraction [-]
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
      real(r8) :: d_wetwat
      real(r8) :: trc_throughfall, trc_pool_total, water_pool_total
      real(r8) :: ratio_layer(1:nl_soil)  ! pre-WATER tracer ratio per layer
      real(r8) :: trc_soil_upflow         ! tracer from soil when qlayer(0)<0
      real(r8) :: gwat_evap              ! surface evaporation removed from gwat [mm]
      real(r8) :: trc_gwat_evap          ! corresponding tracer removed

      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! ============================================================
         ! 0. Compute pre-WATER tracer ratios for all soil layers
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
         !
         !    First handle soil→surface upflow (qlayer(0)<0) so that
         !    the upflow tracer participates in the same-step mixed pool.
         !    Then compute ratio and distribute to wdsrf/rsur/infiltration.
         ! ============================================================

         ! Start with old surface tracer + throughfall from canopy
         trc_pool_total  = trc_wdsrf(itrc, ipatch) + trc_pg_to_ground(itrc, ipatch)

         ! If soil pushes water to surface (qlayer(0)<0), add to mixed pool
         ! BEFORE computing ratio, so it participates in rsur allocation.
         trc_soil_upflow = 0._r8
         IF (qlayer(0) < -trc_tiny) THEN
            trc_soil_upflow = abs(qlayer(0)) * ratio_layer(1) * deltim
            trc_soil_upflow = min(trc_soil_upflow, max(trc_wliq_soisno(itrc, 1, ipatch), 0._r8))
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_soil_upflow
            trc_pool_total = trc_pool_total + trc_soil_upflow
         ENDIF

         ! WATER computes gwat = pg_rain + sm - qseva_soil (no snow case).
         ! The qseva_soil term was subtracted from gwat but not from trc_pool.
         ! Remove the corresponding tracer as evaporation output.
         ! (qsdew_soil adds water/tracer; qfros/qsubl handled later for ice)
         gwat_evap = max(qseva_soil, 0._r8) * deltim
         IF (gwat_evap > trc_tiny .and. trc_pool_total > trc_tiny) THEN
            ! Evaporation draws from the mixed surface pool at pool ratio
            ! Approximate ratio before evaporation
            trc_gwat_evap = gwat_evap * (trc_pool_total / &
               max(trc_pool_total / R_precip, gwat_evap + trc_tiny))
            trc_gwat_evap = min(trc_gwat_evap, trc_pool_total)
            trc_pool_total = trc_pool_total - trc_gwat_evap
            a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_gwat_evap
         ENDIF

         ! Add snowmelt tracer (sm) and soil dew tracer (qsdew_soil)
         ! to the surface pool (these increase gwat)
         IF (sm > trc_tiny) THEN
            ! Snowmelt: tracer from snow layers. Use top snow layer ratio or R_precip.
            ! Phase 1: all snow has R_precip.
            trc_pool_total = trc_pool_total + sm * R_precip * deltim
         ENDIF
         IF (qsdew_soil > trc_tiny) THEN
            ! Dew deposition: atmospheric input
            trc_pool_total = trc_pool_total + qsdew_soil * R_precip * deltim
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + qsdew_soil * R_precip * deltim
         ENDIF

         ! Use OUTPUT-side water pool for ratio computation.
         ! This guarantees: distributed = water_pool * ratio = trc_pool (exact).
         water_pool_total = max(wdsrf, 0._r8) + max(rsur, 0._r8) * deltim &
                          + max(qlayer(0), 0._r8) * deltim

         ! Compute mixed-pool ratio
         IF (water_pool_total > trc_tiny .and. trc_pool_total > trc_tiny) THEN
            ratio = trc_pool_total / water_pool_total
         ELSEIF (trc_pool_total <= trc_tiny) THEN
            ratio = 0._r8
         ELSE
            ratio = R_precip
         ENDIF

         ! Distribute: wdsrf residual, surface runoff, infiltration
         trc_wdsrf(itrc, ipatch) = max(wdsrf, 0._r8) * ratio

         IF (rsur > trc_tiny) THEN
            a_trc_rsur(itrc, ipatch) = a_trc_rsur(itrc, ipatch) + rsur * ratio * deltim
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + rsur * ratio * deltim
         ENDIF

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
         ! ============================================================

         ! --- qlayer(0): surface → layer 1 (downward only) ---
         ! Upward case was already handled above (added to mixed pool).
         IF (qlayer(0) > trc_tiny) THEN
            trc_flux = qlayer(0) * ratio * deltim
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) + trc_flux
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
         ! 4. Subsurface runoff: post qlayer/qcharge bottom-layer ratio
         !    rsub occurs after infiltration, inter-layer exchange, and
         !    groundwater exchange have updated the bottom layer's tracer.
         !    Use the current (post-update) bottom-layer concentration.
         ! ============================================================
         IF (rsub > trc_tiny) THEN
            j = nl_soil
            ! Reconstruct pre-rsub water: WATER's final wliq already has rsub
            ! removed, so add it back to get the state rsub drew from.
            ratio_src = wliq_soisno(j) + rsub * deltim
            IF (ratio_src > trc_tiny) THEN
               ratio_src = trc_wliq_soisno(itrc, j, ipatch) / ratio_src
            ELSE
               ratio_src = ratio_layer(j)
            ENDIF
            trc_flux = rsub * ratio_src * deltim
            trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
            a_trc_rsub(itrc, ipatch) = a_trc_rsub(itrc, ipatch) + trc_flux
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_flux
         ENDIF

         ! ============================================================
         ! 5. WATER-phase ice external fluxes: qfros_soil / qsubl_soil
         !    These are atmosphere↔soil ice exchanges applied by WATER
         !    directly to wice_soisno(1). They are NOT internal freeze/thaw.
         !
         !    Internal freeze/thaw during WATER is not yet tracked here
         !    (needs separating d_wice into external + internal components).
         !    THERMAL-phase freeze/thaw is handled in tracer_evapo.
         ! ============================================================
         IF (qfros_soil > trc_tiny) THEN
            ! Frost deposition: atmosphere → ice in top soil layer
            trc_flux = qfros_soil * R_precip * deltim
            trc_wice_soisno(itrc, 1, ipatch) = trc_wice_soisno(itrc, 1, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         IF (qsubl_soil > trc_tiny) THEN
            ! Sublimation: ice in top soil layer → atmosphere
            IF (wice_soisno_bef(1) > trc_tiny) THEN
               ratio_src = trc_wice_soisno(itrc, 1, ipatch) / max(wice_soisno_bef(1), trc_tiny)
               trc_flux = qsubl_soil * ratio_src * deltim
               trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, 1, ipatch), 0._r8))
               trc_wice_soisno(itrc, 1, ipatch) = trc_wice_soisno(itrc, 1, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

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
