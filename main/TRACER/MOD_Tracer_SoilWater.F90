#include <define.h>

MODULE MOD_Tracer_SoilWater

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, &
      trc_wa, trc_wdsrf, trc_wetwat, &
      a_trc_precip, a_trc_evap, &
      a_trc_qinfl, a_trc_qcharge, a_trc_rsur, a_trc_rsub, a_trc_rnof, &
      trc_pg_to_ground, trc_pg_rain_ground

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
      qseva_in, qsdew_in, qsubl_in, qfros_in, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
      sm, fsno, split_soilsnow, &
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
      real(r8), intent(in) :: qseva_in              ! total evaporation (no _soil suffix) [mm/s]
      real(r8), intent(in) :: qsdew_in              ! total dew [mm/s]
      real(r8), intent(in) :: qsubl_in              ! total sublimation [mm/s]
      real(r8), intent(in) :: qfros_in              ! total frost [mm/s]
      real(r8), intent(in) :: qseva_soil            ! soil-only evaporation [mm/s]
      real(r8), intent(in) :: qsdew_soil            ! soil-only dew [mm/s]
      real(r8), intent(in) :: qsubl_soil            ! soil-only sublimation [mm/s]
      real(r8), intent(in) :: qfros_soil            ! soil-only frost [mm/s]
      real(r8), intent(in) :: qseva_snow            ! snow-only evaporation [mm/s]
      real(r8), intent(in) :: qsdew_snow            ! snow-only dew [mm/s]
      real(r8), intent(in) :: qsubl_snow            ! snow-only sublimation [mm/s]
      real(r8), intent(in) :: qfros_snow            ! snow-only frost [mm/s]
      real(r8), intent(in) :: sm                    ! snowmelt [mm/s]
      real(r8), intent(in) :: fsno                  ! snow fraction [-]
      logical,  intent(in) :: split_soilsnow        ! DEF_SPLIT_SOILSNOW
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)     ! post-WATER
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)     ! post-WATER
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil) ! pre-WATER (post-THERMAL)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil) ! pre-WATER
      real(r8), intent(in) :: wa, wa_bef, wdsrf, wdsrf_bef
      real(r8), intent(in) :: wetwat, wetwat_bef
      real(r8), intent(in) :: pg_rain, pg_snow

      integer  :: itrc, j, lb, lb_snow
      real(r8) :: R_precip
      real(r8) :: trc_flux, ratio, ratio_src
      real(r8) :: d_wice, d_wetwat
      real(r8) :: trc_pool_total, water_pool_total
      real(r8) :: ratio_layer(1:nl_soil)  ! pre-WATER tracer ratio per layer
      real(r8) :: trc_soil_upflow         ! tracer from soil when qlayer(0)<0
      ! Snow top-layer effective flux selectors
      real(r8) :: eff_qseva_snow, eff_qsdew_snow, eff_qsubl_snow, eff_qfros_snow
      ! Cross-phase deficit variables
      real(r8) :: water_ice_pool, water_liq_pool  ! post-deposit pool sizes
      real(r8) :: subl_water, evap_water           ! water to remove
      real(r8) :: deficit_water                    ! excess beyond primary phase
      real(r8) :: gwat_evap              ! evaporation subtracted from gwat [mm]
      real(r8) :: trc_gwat_evap         ! corresponding tracer removed
      real(r8) :: eff_qseva             ! effective evap used in gwat
      real(r8) :: eff_qsdew_topliq      ! effective dew on soil layer 1 liquid
      real(r8) :: eff_qsubl_top         ! effective sublimation on soil layer 1 ice
      real(r8) :: eff_qfros_top         ! effective frost on soil layer 1 ice
      real(r8) :: pool_ratio            ! actual ratio of mixed pool
      real(r8) :: trc_gwat_snow         ! tracer leaving snow bottom [kg/m2 per step]
      real(r8) :: snow_rain_input       ! rain entering snowwater [mm]
      real(r8) :: snow_dew_input        ! dew entering top snow liquid [mm]
      real(r8) :: snow_evap_output      ! evap removed from top snow liquid [mm]
      real(r8) :: snow_frost_input      ! frost added to top snow ice [mm]
      real(r8) :: snow_subl_output      ! sublimation removed from top snow ice [mm]
      real(r8) :: snow_wliq_before_flow ! top snow liquid before percolation [mm]
      real(r8) :: qin_snow, qout_snow   ! reconstructed snowwater qin/qout [mm]
      real(r8) :: trc_qin_snow, trc_qout_snow
      real(r8) :: water_before_flow, trc_before_flow

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
         ! 0b. Snow layer external fluxes (when snow exists)
         !     snowwater applies qseva_snow/qsdew_snow/qsubl_snow/qfros_snow
         !     to the top snow layer. Track these as system I/O.
         ! ============================================================
         IF (snl < 0) THEN
            lb_snow = snl + 1  ! top snow layer index

            ! Select flux variables based on split path
            IF (split_soilsnow) THEN
               eff_qseva_snow = qseva_snow
               eff_qsdew_snow = qsdew_snow
               eff_qsubl_snow = qsubl_snow
               eff_qfros_snow = qfros_snow
            ELSE
               eff_qseva_snow = qseva_in
               eff_qsdew_snow = qsdew_in
               eff_qsubl_snow = qsubl_in
               eff_qfros_snow = qfros_in
            ENDIF

            ! Step 1: Ice phase — frost deposits, sublimation removes
            ! snowwater does: wice(lb) += (qfros - qsubl)*dt
            ! If wice goes negative, deficit removes from wliq.
            trc_flux = max(eff_qfros_snow, 0._r8) * R_precip * deltim
            trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux

            ! Post-frost ice pool (water and tracer)
            water_ice_pool = wice_soisno_bef(lb_snow) + max(eff_qfros_snow, 0._r8) * deltim
            subl_water = max(eff_qsubl_snow, 0._r8) * deltim

            IF (subl_water > trc_tiny) THEN
               IF (water_ice_pool > trc_tiny .and. subl_water <= water_ice_pool) THEN
                  ! Normal case: enough ice to cover sublimation
                  ratio_src = trc_wice_soisno(itrc, lb_snow, ipatch) / water_ice_pool
                  trc_flux = subl_water * ratio_src
                  trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
               ELSE
                  ! Ice insufficient (or empty): drain all ice, remainder from liquid
                  ! Part 1: drain ice completely
                  trc_flux = max(trc_wice_soisno(itrc, lb_snow, ipatch), 0._r8)
                  trc_wice_soisno(itrc, lb_snow, ipatch) = 0._r8
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
                  ! Part 2: excess from liquid at LIQUID ratio
                  deficit_water = subl_water - max(water_ice_pool, 0._r8)
                  IF (deficit_water > trc_tiny .and. wliq_soisno_bef(lb_snow) > trc_tiny) THEN
                     ratio_src = trc_wliq_soisno(itrc, lb_snow, ipatch) / wliq_soisno_bef(lb_snow)
                     trc_flux = deficit_water * ratio_src
                     trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, lb_snow, ipatch), 0._r8))
                     trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
                     a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
                  ENDIF
               ENDIF
            ENDIF
            ! Update water_ice_pool to post-subl value for use in Step 2's ice deficit
            water_ice_pool = max(water_ice_pool - subl_water, 0._r8)

            ! Step 2: Liquid phase — rain + dew deposits, evaporation removes
            ! snowwater does: wliq(lb) += (pg_rain + qsdew - qseva)*dt
            ! If wliq goes negative, deficit removes from wice.
            !
            ! Add RAIN tracer here (before evaporation), mirroring snowwater's
            ! order: wliq += (pg_rain + qsdew - qseva)*dt applied together.
            ! This ensures evaporation can draw from rain that arrived this step.
            ! The rain tracer is then NOT re-added in section 0c percolation.
            IF (split_soilsnow) THEN
               trc_flux = trc_pg_rain_ground(itrc, ipatch) * fsno
            ELSE
               trc_flux = trc_pg_rain_ground(itrc, ipatch)
            ENDIF
            trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) + trc_flux

            ! Add dew
            trc_flux = max(eff_qsdew_snow, 0._r8) * R_precip * deltim
            trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux

            ! Evaporate from the full pool (wliq_bef + qsubl_deficit + rain + dew)
            IF (eff_qseva_snow > trc_tiny) THEN
               ! Reconstruct pre-evap water pool (matches snowwater):
               ! wliq_bef + (pg_rain + qsdew)*dt + any ice deficit transferred earlier
               IF (split_soilsnow) THEN
                  water_liq_pool = wliq_soisno_bef(lb_snow) + (pg_rain*fsno + max(eff_qsdew_snow,0._r8)) * deltim
               ELSE
                  water_liq_pool = wliq_soisno_bef(lb_snow) + (pg_rain + max(eff_qsdew_snow,0._r8)) * deltim
               ENDIF
               ! Account for qsubl deficit transferred from ice (Step 1)
               IF (water_ice_pool < eff_qsubl_snow * deltim) THEN
                  water_liq_pool = water_liq_pool - (eff_qsubl_snow * deltim - water_ice_pool)
               ENDIF
               water_liq_pool = max(water_liq_pool, 0._r8)
               evap_water = eff_qseva_snow * deltim

               IF (water_liq_pool > trc_tiny .and. evap_water <= water_liq_pool) THEN
                  ! Normal case: enough liquid
                  ratio_src = trc_wliq_soisno(itrc, lb_snow, ipatch) / water_liq_pool
                  trc_flux = evap_water * ratio_src
                  trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
               ELSE
                  ! Liquid insufficient (or empty): drain liquid, remainder from ice
                  ! Part 1: drain liquid completely
                  trc_flux = max(trc_wliq_soisno(itrc, lb_snow, ipatch), 0._r8)
                  trc_wliq_soisno(itrc, lb_snow, ipatch) = 0._r8
                  a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
                  ! Part 2: excess from ice at post-frost-post-subl ICE ratio
                  deficit_water = evap_water - max(water_liq_pool, 0._r8)
                  IF (deficit_water > trc_tiny .and. water_ice_pool > trc_tiny) THEN
                     ! water_ice_pool is already post-frost, post-subl (updated at end of Step 1)
                     ratio_src = trc_wice_soisno(itrc, lb_snow, ipatch) / water_ice_pool
                     trc_flux = deficit_water * ratio_src
                     trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, lb_snow, ipatch), 0._r8))
                     trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
                     a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
                  ENDIF
               ENDIF
            ENDIF

         ENDIF  ! snl < 0

         ! ============================================================
         ! 0c. Snow percolation tracer tracking
         !
         !     Reconstruct the explicit snowwater loop:
         !       1. Apply top-layer external liquid/ice fluxes
         !       2. March qin/qout from top snow layer to bottom
         !       3. Bottom outflow becomes gwat tracer into surface pool
         !
         !     Use the separated rain tracer from tracer_precip so only
         !     liquid rain enters snowwater; snowfall tracer stays in snow.
         ! ============================================================
         trc_gwat_snow = 0._r8
         IF (snl < 0) THEN
            IF (split_soilsnow) THEN
               snow_rain_input  = max(pg_rain * fsno, 0._r8) * deltim
               snow_dew_input   = max(qsdew_snow, 0._r8) * deltim
               snow_evap_output = max(qseva_snow, 0._r8) * deltim
               snow_frost_input = max(qfros_snow, 0._r8) * deltim
               snow_subl_output = max(qsubl_snow, 0._r8) * deltim
            ELSE
               snow_rain_input  = max(pg_rain, 0._r8) * deltim
               snow_dew_input   = max(qsdew_in, 0._r8) * deltim
               snow_evap_output = max(qseva_in, 0._r8) * deltim
               snow_frost_input = max(qfros_in, 0._r8) * deltim
               snow_subl_output = max(qsubl_in, 0._r8) * deltim
            ENDIF
            ! Rain tracer was already added to trc_wliq in Step 2 (section 0b).
            ! Top snow layer starts percolation with trc_qin_snow = 0.
            trc_qin_snow = 0._r8

            snow_wliq_before_flow = wliq_soisno_bef(lb_snow)
            IF (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output < 0._r8) THEN
               snow_wliq_before_flow = snow_wliq_before_flow + &
                  (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output)
            ENDIF
            snow_wliq_before_flow = snow_wliq_before_flow + snow_rain_input + snow_dew_input - snow_evap_output
            IF (snow_wliq_before_flow < 0._r8) THEN
               snow_wliq_before_flow = 0._r8
            ENDIF

            qin_snow = 0._r8
            DO j = lb_snow, 0
               IF (j == lb_snow) THEN
                  water_before_flow = max(snow_wliq_before_flow, 0._r8)
                  trc_before_flow = max(trc_wliq_soisno(itrc, j, ipatch), 0._r8) + trc_qin_snow
               ELSE
                  water_before_flow = max(wliq_soisno_bef(j), 0._r8) + qin_snow
                  trc_before_flow = max(trc_wliq_soisno(itrc, j, ipatch), 0._r8) + trc_qin_snow
               ENDIF

               qout_snow = max(water_before_flow - max(wliq_soisno(j), 0._r8), 0._r8)
               qout_snow = min(qout_snow, water_before_flow)

               IF (qout_snow > trc_tiny .and. water_before_flow > trc_tiny .and. trc_before_flow > trc_tiny) THEN
                  ratio_src = trc_before_flow / water_before_flow
                  trc_qout_snow = qout_snow * ratio_src
                  trc_qout_snow = min(trc_qout_snow, trc_before_flow)
               ELSE
                  trc_qout_snow = 0._r8
               ENDIF

               trc_wliq_soisno(itrc, j, ipatch) = max(trc_before_flow - trc_qout_snow, 0._r8)

               qin_snow = qout_snow
               trc_qin_snow = trc_qout_snow
            ENDDO

            trc_gwat_snow = trc_qin_snow
         ENDIF

         ! ============================================================
         ! 1. Surface water: mixed-pool approach
         !
         !    First handle soil→surface upflow (qlayer(0)<0) so that
         !    the upflow tracer participates in the same-step mixed pool.
         !    Then compute ratio and distribute to wdsrf/rsur/infiltration.
         ! ============================================================

         ! If soil pushes water to surface (qlayer(0)<0), add to mixed pool
         ! BEFORE computing ratio, so it participates in rsur allocation.
         trc_soil_upflow = 0._r8
         IF (qlayer(0) < -trc_tiny) THEN
            trc_soil_upflow = abs(qlayer(0)) * ratio_layer(1) * deltim
            trc_soil_upflow = min(trc_soil_upflow, max(trc_wliq_soisno(itrc, 1, ipatch), 0._r8))
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_soil_upflow
         ENDIF

         ! ---- Determine effective fluxes applied to gwat / soil layer 1 ----
         !
         ! Three cases for how WATER handles evaporation:
         !
         ! Case A: no snow (lb>=1)
         !   gwat = pg_rain + sm - qseva  (= qseva_soil, same thing)
         !   ice fluxes on soil layer 1: qfros/qsubl (= qfros_soil/qsubl_soil)
         !   → surface pool must deduct qseva, add qsdew
         !
         ! Case B: snow + split_soilsnow
         !   snowwater handles snow portion with qseva_snow/qsdew_snow/etc.
         !   gwat from snow bottom → then gwat += pg_rain*(1-fsno) - qseva_soil
         !   ice fluxes on soil layer 1: qfros_soil/qsubl_soil
         !   → surface pool deducts qseva_soil (soil portion only)
         !
         ! Case C: snow + NOT split_soilsnow
         !   snowwater handles ALL qseva/qsdew/qsubl/qfros on snow layers
         !   gwat is snowwater's OUTPUT (already post-evaporation)
         !   → surface pool should NOT deduct qseva again (already in gwat)
         !   → ice fluxes on soil layer 1: none separate (snowwater handled them)
         !     but qfros/qsubl still may affect soil layer 1 via wice update
         !
         IF (snl < 0 .and. .not. split_soilsnow) THEN
            ! Case C: evap/dew handled inside snowwater, not in surface pool
            eff_qseva = 0._r8
            eff_qsdew_topliq = 0._r8
            eff_qsubl_top = 0._r8
            eff_qfros_top = 0._r8
         ELSEIF (snl < 0 .and. split_soilsnow) THEN
            ! Case B: only soil portion affects surface pool
            eff_qseva = qseva_soil
            eff_qsdew_topliq = qsdew_soil
            eff_qsubl_top = qsubl_soil
            eff_qfros_top = qfros_soil
         ELSE
            ! Case A: no snow, all fluxes apply
            eff_qseva = qseva_in
            eff_qsdew_topliq = qsdew_in
            eff_qsubl_top = qsubl_in
            eff_qfros_top = qfros_in
         ENDIF

         ! Build pre-evaporation surface tracer pool.
         IF (snl < 0) THEN
            trc_pool_total = trc_wdsrf(itrc, ipatch) + trc_soil_upflow + trc_gwat_snow
            IF (split_soilsnow) THEN
               trc_pool_total = trc_pool_total + trc_pg_rain_ground(itrc, ipatch) * (1._r8 - fsno)
            ENDIF
         ELSE
            trc_pool_total = trc_wdsrf(itrc, ipatch) + trc_pg_to_ground(itrc, ipatch) &
                           + trc_soil_upflow
         ENDIF

         ! ---- Remove evaporation tracer from mixed pool ----
         ! gwat -= eff_qseva, so trc_pool must also be reduced.
         ! Use the ACTUAL pool ratio: trc_pool / (trc_pool / pool_ratio)
         ! where pool_ratio = trc_pool / water_pool_before_evap.
         ! Simplification: evap tracer = evap_water * (trc_pool / water_pool_before_evap)
         ! water_pool_before_evap ≈ trc_pool / pool_ratio. We don't know the exact
         ! water pool size at this point, but we can use:
         !   pool_ratio = trc_pool / water_pool_before_evap
         !   trc_evap = eff_qseva * dt * pool_ratio = eff_qseva * dt * trc_pool / water_bef
         ! water_bef = water_pool_after_evap + eff_qseva*dt
         !           = (wdsrf + rsur*dt + qlayer(0)*dt) + eff_qseva*dt  [reconstruct pre-evap]
         gwat_evap = max(eff_qseva, 0._r8) * deltim
         IF (gwat_evap > trc_tiny .and. trc_pool_total > trc_tiny) THEN
            ! Reconstruct pre-evap water pool
            water_pool_total = max(wdsrf, 0._r8) + max(rsur, 0._r8) * deltim &
                             + max(qlayer(0), 0._r8) * deltim + gwat_evap
            pool_ratio = trc_pool_total / max(water_pool_total, trc_tiny)
            trc_gwat_evap = gwat_evap * pool_ratio
            trc_gwat_evap = min(trc_gwat_evap, trc_pool_total)
            trc_pool_total = trc_pool_total - trc_gwat_evap
            a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_gwat_evap
         ENDIF

         ! ---- Add no-snow snowmelt tracer to pool ----
         IF (snl >= 0 .and. sm > trc_tiny) THEN
            trc_pool_total = trc_pool_total + sm * R_precip * deltim
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
         ! 4b. WATER-phase liquid external flux on soil layer 1: dew
         !     qsdew / qsdew_soil are applied directly to wliq_soisno(1)
         !     by WATER, not to the surface mixed pool.
         ! ============================================================
         IF (eff_qsdew_topliq > trc_tiny) THEN
            trc_flux = eff_qsdew_topliq * R_precip * deltim
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
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
         ! Use effective fluxes (handles both split and non-split paths)
         IF (eff_qfros_top > trc_tiny) THEN
            trc_flux = eff_qfros_top * R_precip * deltim
            trc_wice_soisno(itrc, 1, ipatch) = trc_wice_soisno(itrc, 1, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         IF (eff_qsubl_top > trc_tiny) THEN
            IF (wice_soisno_bef(1) > trc_tiny) THEN
               ratio_src = trc_wice_soisno(itrc, 1, ipatch) / max(wice_soisno_bef(1), trc_tiny)
               trc_flux = eff_qsubl_top * ratio_src * deltim
               trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, 1, ipatch), 0._r8))
               trc_wice_soisno(itrc, 1, ipatch) = trc_wice_soisno(itrc, 1, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! ============================================================
         ! 5b. Internal freeze/thaw during WATER
         !
         !     d_wice_total = wice_post - wice_pre (total ice change)
         !     d_wice_external = (qfros_soil - qsubl_soil) * dt  (layer 1 only)
         !     d_wice_internal = d_wice_total - d_wice_external  (pure phase change)
         !
         !     d_wice_internal > 0 → freeze (liquid → ice)
         !     d_wice_internal < 0 → thaw  (ice → liquid)
         ! ============================================================
         DO j = lb, nl_soil
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            ! Subtract external ice fluxes to isolate internal freeze/thaw
            IF (j == 1) THEN
               ! Soil layer 1: subtract effective external ice flux
               d_wice = d_wice - (eff_qfros_top - eff_qsubl_top) * deltim
            ENDIF
            IF (j < 1 .and. j == lb .and. snl < 0) THEN
               ! Top snow layer: external ice flux is applied in snowwater.
               IF (split_soilsnow) THEN
                  d_wice = d_wice - (qfros_snow - qsubl_snow) * deltim
               ELSE
                  d_wice = d_wice - (qfros_in - qsubl_in) * deltim
               ENDIF
            ENDIF

            IF (d_wice > trc_tiny) THEN
               ! Internal freeze: liquid → ice
               IF (wliq_soisno_bef(j) > trc_tiny) THEN
                  ratio_src = trc_wliq_soisno(itrc, j, ipatch) / max(wliq_soisno_bef(j), trc_tiny)
                  trc_flux = d_wice * ratio_src
                  trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ELSEIF (d_wice < -trc_tiny) THEN
               ! Internal thaw: ice → liquid
               IF (wice_soisno_bef(j) > trc_tiny) THEN
                  ratio_src = trc_wice_soisno(itrc, j, ipatch) / max(wice_soisno_bef(j), trc_tiny)
                  trc_flux = abs(d_wice) * ratio_src
                  trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
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
