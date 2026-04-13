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
      real(r8) :: trc_throughfall, trc_pool_total, water_pool_total
      real(r8) :: ratio_layer(1:nl_soil)  ! pre-WATER tracer ratio per layer
      real(r8) :: trc_soil_upflow         ! tracer from soil when qlayer(0)<0
      real(r8) :: gwat_evap              ! evaporation subtracted from gwat [mm]
      real(r8) :: gwat_dew              ! dew added to gwat [mm]
      real(r8) :: trc_gwat_evap         ! corresponding tracer removed
      real(r8) :: eff_qseva             ! effective evap used in gwat
      real(r8) :: eff_qsdew             ! effective dew used in gwat
      real(r8) :: eff_qsubl_top         ! effective sublimation on top layer (ice)
      real(r8) :: eff_qfros_top         ! effective frost on top layer (ice)
      real(r8) :: pool_ratio            ! actual ratio of mixed pool
      ! Snow percolation reconstruction variables
      real(r8) :: trc_gwat_snow         ! tracer exiting snow bottom → surface pool
      real(r8) :: srf_input_water       ! liquid water input to top snow layer [mm]
      real(r8) :: trc_srf_input         ! tracer in surface input to snow
      real(r8) :: qin_water, trc_qin    ! water/tracer flowing into current snow layer
      real(r8) :: qout_water, trc_qout  ! water/tracer flowing out of current snow layer
      real(r8) :: d_wliq_j              ! per-layer wliq delta
      real(r8) :: trc_total_j           ! total tracer in layer after receiving qin
      real(r8) :: water_total_j         ! total water in layer after receiving qin

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

            IF (split_soilsnow) THEN
               ! DEF_SPLIT_SOILSNOW: snow fluxes use _snow variables
               IF (qsdew_snow > trc_tiny) THEN
                  trc_flux = qsdew_snow * R_precip * deltim
                  trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) + trc_flux
                  a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
               ENDIF
               IF (qseva_snow > trc_tiny) THEN
                  IF (wliq_soisno_bef(lb_snow) > trc_tiny) THEN
                     ratio_src = trc_wliq_soisno(itrc, lb_snow, ipatch) / max(wliq_soisno_bef(lb_snow), trc_tiny)
                     trc_flux = qseva_snow * ratio_src * deltim
                     trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, lb_snow, ipatch), 0._r8))
                     trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
                     a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
                  ENDIF
               ENDIF
               IF (qfros_snow > trc_tiny) THEN
                  trc_flux = qfros_snow * R_precip * deltim
                  trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) + trc_flux
                  a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
               ENDIF
               IF (qsubl_snow > trc_tiny) THEN
                  IF (wice_soisno_bef(lb_snow) > trc_tiny) THEN
                     ratio_src = trc_wice_soisno(itrc, lb_snow, ipatch) / max(wice_soisno_bef(lb_snow), trc_tiny)
                     trc_flux = qsubl_snow * ratio_src * deltim
                     trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, lb_snow, ipatch), 0._r8))
                     trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
                     a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
                  ENDIF
               ENDIF
            ENDIF
            ! Note: when NOT split_soilsnow, qseva/qsdew/qsubl/qfros (without _snow suffix)
            ! are applied to snow layers by snowwater. These are the same as qseva_soil etc.
            ! but we already handle qseva_soil in the surface pool section below.
            ! The snow layer internal changes (percolation etc.) are captured by the
            ! delta approach for snow layers in sections 5a/5b.
         ENDIF

         ! ============================================================
         ! 0c. Snow percolation tracer tracking (when snow exists)
         !
         !     Reconstruct layer-by-layer qin/qout from wliq deltas:
         !       qout(j) = qin(j) - Δwliq(j)
         !     At each layer, incoming tracer mixes with existing,
         !     outflow carries the mixed ratio.
         !
         !     External fluxes (qseva/qsdew/qfros/qsubl) were already
         !     handled in 0b, so wliq_bef already reflects pre-percolation
         !     state EXCEPT for the surface input.
         !
         !     Output: trc_gwat_snow = tracer exiting snow bottom → surface pool
         ! ============================================================
         trc_gwat_snow = 0._r8
         IF (snl < 0) THEN
            ! Surface liquid input to top snow layer (same as snowwater):
            IF (split_soilsnow) THEN
               srf_input_water = (pg_rain * fsno + max(qsdew_snow,0._r8) - max(qseva_snow,0._r8)) * deltim
            ELSE
               srf_input_water = (pg_rain + max(qsdew_in,0._r8) - max(qseva_in,0._r8)) * deltim
            ENDIF
            srf_input_water = max(srf_input_water, 0._r8)

            ! Tracer in surface input: rain portion at throughfall ratio,
            ! dew portion at R_precip (already handled in 0b for split path)
            ! For simplicity, use throughfall fraction of trc_pg_to_ground
            IF (split_soilsnow) THEN
               ! Only the fsno fraction of throughfall goes through snow
               trc_srf_input = trc_pg_to_ground(itrc, ipatch) * fsno
            ELSE
               trc_srf_input = trc_pg_to_ground(itrc, ipatch)
            ENDIF

            ! Sequential top→bottom percolation
            qin_water = srf_input_water
            trc_qin   = trc_srf_input
            DO j = snl + 1, 0
               d_wliq_j = wliq_soisno(j) - wliq_soisno_bef(j)
               ! Reconstruct outflow: qout = qin - Δwliq
               qout_water = max(qin_water - d_wliq_j, 0._r8)

               ! Mix incoming tracer with existing layer tracer
               trc_total_j = trc_wliq_soisno(itrc, j, ipatch) + trc_qin
               water_total_j = max(wliq_soisno_bef(j) + qin_water, trc_tiny)

               IF (qout_water > trc_tiny) THEN
                  ratio_src = trc_total_j / water_total_j
                  trc_qout = qout_water * ratio_src
                  trc_qout = min(trc_qout, max(trc_total_j, 0._r8))
               ELSE
                  trc_qout = 0._r8
               ENDIF

               ! Update layer tracer
               trc_wliq_soisno(itrc, j, ipatch) = max(trc_total_j - trc_qout, 0._r8)

               ! Set up for next layer
               qin_water = qout_water
               trc_qin = trc_qout
            ENDDO

            ! gwat tracer = what exits bottom snow layer
            trc_gwat_snow = trc_qout
         ENDIF

         ! ============================================================
         ! 1. Surface water: mixed-pool approach
         !
         !    First handle soil→surface upflow (qlayer(0)<0) so that
         !    the upflow tracer participates in the same-step mixed pool.
         !    Then compute ratio and distribute to wdsrf/rsur/infiltration.
         ! ============================================================

         ! Build surface pool tracer:
         !   If snow: gwat tracer from snow bottom + bare-ground throughfall
         !   If no snow: full throughfall from canopy
         IF (snl < 0 .and. split_soilsnow) THEN
            ! Split: snow portion went through snow → trc_gwat_snow
            !        bare portion = (1-fsno) fraction of throughfall → direct
            trc_pool_total = trc_wdsrf(itrc, ipatch) &
               + trc_gwat_snow + trc_pg_to_ground(itrc, ipatch) * (1._r8 - fsno)
         ELSEIF (snl < 0) THEN
            ! Not split: all rain went through snow → trc_gwat_snow
            trc_pool_total = trc_wdsrf(itrc, ipatch) + trc_gwat_snow
         ELSE
            ! No snow: all throughfall goes directly to surface
            trc_pool_total = trc_wdsrf(itrc, ipatch) + trc_pg_to_ground(itrc, ipatch)
         ENDIF

         ! If soil pushes water to surface (qlayer(0)<0), add to mixed pool
         ! BEFORE computing ratio, so it participates in rsur allocation.
         trc_soil_upflow = 0._r8
         IF (qlayer(0) < -trc_tiny) THEN
            trc_soil_upflow = abs(qlayer(0)) * ratio_layer(1) * deltim
            trc_soil_upflow = min(trc_soil_upflow, max(trc_wliq_soisno(itrc, 1, ipatch), 0._r8))
            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_soil_upflow
            trc_pool_total = trc_pool_total + trc_soil_upflow
         ENDIF

         ! ---- Determine effective fluxes used in gwat computation ----
         ! WATER_VSF uses different variables depending on path:
         !   split_soilsnow AND snow:  gwat includes qseva_soil (soil portion)
         !   NOT split_soilsnow AND snow: gwat includes qseva (total, via snowwater)
         !   no snow (lb>=1):  gwat = pg_rain + sm - qseva_soil  (or qseva, same thing)
         IF (split_soilsnow .and. snl < 0) THEN
            eff_qseva = qseva_soil     ! soil portion only
            eff_qsdew = qsdew_soil
            eff_qsubl_top = qsubl_soil  ! on soil layer 1
            eff_qfros_top = qfros_soil
         ELSE
            eff_qseva = qseva_in       ! total (includes snow portion when no split)
            eff_qsdew = qsdew_in
            eff_qsubl_top = qsubl_in
            eff_qfros_top = qfros_in
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

         ! ---- Add snowmelt and dew tracer to pool ----
         IF (sm > trc_tiny) THEN
            trc_pool_total = trc_pool_total + sm * R_precip * deltim
         ENDIF
         gwat_dew = max(eff_qsdew, 0._r8) * deltim
         IF (gwat_dew > trc_tiny) THEN
            trc_pool_total = trc_pool_total + gwat_dew * R_precip
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + gwat_dew * R_precip
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
            IF (j < 1 .and. j == lb .and. snl < 0 .and. split_soilsnow) THEN
               ! Top snow layer (split path): has qfros_snow/qsubl_snow
               d_wice = d_wice - (qfros_snow - qsubl_snow) * deltim
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
