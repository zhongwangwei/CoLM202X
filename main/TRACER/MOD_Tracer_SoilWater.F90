#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_SoilWater

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, trc_delta_sanity_max, &
      tracer_init_water_ratio
   USE MOD_Tracer_Forcing, only: tracer_forcing_precip_value, tracer_forcing_vapor_value, &
      tracer_forcing_has_vapor
   USE MOD_Tracer_Frac, only: tracer_fractionation_active, tracer_surface_relhum, &
      tracer_transpiration_nss_ratio, tracer_diffusivity_ratio_air, &
      tracer_craig_gordon_evap_ratio, tracer_alpha_liq_vap, &
      tracer_rayleigh_freezing_loss, tracer_equilibrium_deposition_ratio
	   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, &
	      trc_wa, trc_wdsrf, trc_wetwat, trc_waterstorage, &
		      a_trc_precip, a_trc_transp_src, tracer_book_evap_loss, &
		      TRC_EVAP_KIND_TRANSP, TRC_EVAP_KIND_SOILEVAP, TRC_EVAP_KIND_SUBL, &
		      TRC_EVAP_KIND_WETLAND, &
		      a_trc_qinfl, a_trc_qcharge, a_trc_rsur, a_trc_rsub, a_trc_rnof, &
	      trc_pg_rain_ground, trc_rnof_step, trc_sm_carry, &
	      trc_leaf_delta_e, trc_leaf_delta_b, trc_leaf_peclet, trc_leaf_water_moles, &
	      trc_leaf_iso_storage, trc_numerical_residual_step

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
      qlayer, qinfl, qcharge, rsur, rsub, &
      qseva_in, qsdew_in, qsubl_in, qfros_in, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
      sm, fsno, split_soilsnow, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      wa, wa_bef, wdsrf, wdsrf_bef, &
      wetwat, wetwat_bef, pg_rain, pg_snow, &
      etroot, wblc_ice_sink, &
      etroot_actual, etroot_aquifer, &
      qflx_irrig_ground, waterstorage_patch, &
      tleaf_frac, t_soisno_frac, forc_q_frac, forc_psrf_frac, lai_frac, rst_frac)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: qlayer(0:nl_soil)     ! inter-layer flux [mm/s]
      real(r8), intent(in) :: qinfl                 ! diagnostic infiltration flux [mm/s]
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
      ! Per-layer transpiration demand (mm/s) as applied by WATER_VSF's
      ! soil_water_vertical_movement (MOD_Hydro_SoilWater.F90:296). The
      ! tracer path subtracts etroot(j)*dt at each layer's pre-WATER ratio
      ! so transpired water carries its tracer to a_trc_evap — previously
      ! this was a silent sink that concentrated trc_wliq over time without
      ! breaking xerr_tracer (both Δstorage and output were dropped).
      real(r8), intent(in) :: etroot(1:nl_soil)
      ! Per-layer ice mass (kg/m2) withdrawn by WATER_VSF's wblc
      ! reconciliation (MOD_SoilSnowHydrology.F90:1069-1081). The water
      ! left the column as ET output, so the tracer path must remove the
      ! matching trc_wice mass (at the layer's ice ratio) and account it
      ! in a_trc_evap. Previously this appeared to Section 5b as a thaw
      ! (d_wice<0), silently moving tracer from ice to liquid instead of
      ! out of the column — producing over-concentration on cold-region /
      ! high-ET patches. Zero for wetland branches that skip wblc.
      real(r8), intent(in) :: wblc_ice_sink(1:nl_soil)
      ! Per-layer actual ET water (mm, not mm/s) after the deficit
      ! cascade in MOD_Hydro_SoilWater.F90:soil_water_vertical_movement,
      ! and the aquifer share absorbing any remaining deficit (mm).
      ! Tracer path uses these instead of the raw demand `etroot(j)*dt`
      ! so that a shallow dry layer reports 0 (its demand cascaded to
      ! a deeper layer, which now reports that extra uptake) and the
      ! aquifer deduction is applied via trc_wa. Without this, dry
      ! layers would get spurious tracer removal clipped by min() and
      ! deep-cascade / aquifer ET would silently leave tracer behind.
      real(r8), intent(in) :: etroot_actual(1:nl_soil)
      real(r8), intent(in) :: etroot_aquifer

      ! Combined drip + flood + paddy irrigation rate (mm/s). WATER_VSF
      ! adds this to gwat at MOD_SoilSnowHydrology.F90:828-832 AFTER the
      ! qseva subtraction, so post-WATER wdsrf + rsur*dt + qlayer(0)*dt
      ! implicitly contains the irrigation water. The tracer path must
      ! inject the matching R*dt into the mixed surface pool or the
      ! downstream wdsrf/rsur/qlayer(0) distribution is under-supplied
      ! and the tracer budget shows a spurious step-level dilution.
      ! Sprinkler is handled separately in tracer_precip (merged into
      ! the rain/canopy path by LEAF_INTERCEPTION). 0 when irrigation
      ! is disabled.
      real(r8), intent(in) :: qflx_irrig_ground
      real(r8), intent(in), optional :: waterstorage_patch
      real(r8), intent(in), optional :: tleaf_frac
      real(r8), intent(in), optional :: t_soisno_frac(snl+1:nl_soil)
      real(r8), intent(in), optional :: forc_q_frac
      real(r8), intent(in), optional :: forc_psrf_frac
      real(r8), intent(in), optional :: lai_frac
      real(r8), intent(in), optional :: rst_frac

      integer  :: itrc, j, lb, lb_snow
      real(r8) :: R_precip, R_atm
      real(r8) :: trc_flux, ratio, ratio_src, storage_ratio
      real(r8) :: qcharge_eff
      real(r8) :: d_wice, d_wetwat
      real(r8) :: trc_pool_total, water_pool_total
	      real(r8) :: ratio_layer(1:nl_soil)  ! pre-WATER tracer ratio per layer
	      real(r8) :: water_shadow(1:nl_soil)
	      real(r8) :: water_resid, water_shadow_ratio
      real(r8) :: trc_soil_upflow         ! tracer from true soil-to-surface upflow
      real(r8) :: trc_soil_evap           ! tracer removed by top-boundary evaporation deficit
      real(r8) :: qgtop_est
      ! Snow top-layer effective flux selectors
      real(r8) :: eff_qseva_snow, eff_qsdew_snow, eff_qsubl_snow, eff_qfros_snow
      ! Cross-phase deficit variables
      real(r8) :: water_ice_pool, water_liq_pool  ! post-deposit pool sizes
      real(r8) :: subl_water, evap_water           ! water to remove
      real(r8) :: deficit_water                    ! excess beyond primary phase
      real(r8) :: water_ice_pool_prefrost          ! pre-subl ice pool for deficit check
      real(r8) :: gwat_evap              ! evaporation subtracted from gwat [mm]
      real(r8) :: trc_gwat_evap         ! corresponding tracer removed
      real(r8) :: eff_qseva             ! effective evap used in gwat
      real(r8) :: eff_qsdew_topliq      ! effective dew on soil layer 1 liquid
      real(r8) :: eff_qsubl_top         ! effective sublimation on soil layer 1 ice
      real(r8) :: eff_qfros_top         ! effective frost on soil layer 1 ice
      real(r8) :: pool_ratio            ! actual ratio of mixed pool
      real(r8) :: top_infil_water       ! surface-to-soil infiltration water [mm/step]
      real(r8) :: top_exfil_water       ! soil-to-surface exfiltration water [mm/step]
      real(r8) :: top_boundary_out_water ! total negative qinfl water [mm/step]
      real(r8) :: top_soil_evap_water   ! negative qinfl caused by qseva deficit [mm/step]
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
      real(r8) :: wliq_pre_phase, wice_pre_phase   ! phase-change basis
      ! Actual external wice change on top snow layer after deficit handling:
      ! nominal (qfros - qsubl)*dt understates removal when qsubl
      ! exceeds post-frost ice, because the deficit is taken from liquid
      ! instead. Section 5b uses this to isolate internal freeze/thaw.
      real(r8) :: d_wice_ext_snow
      ! Actual external wice change on soil layer 1 after the
      ! max(0., wice_bef + (qfros-qsubl)*dt) clamp in
      ! MOD_SoilSnowHydrology.F90:1101. Symmetric to d_wice_ext_snow for
      ! top snow layer: when qsubl*dt exceeds wice_bef+qfros*dt,
      ! the water side floors wice at 0, so the true external change is
      ! -wice_bef rather than the nominal (qfros-qsubl)*dt. Without this
      ! correction, Section 5b's d_wice residual becomes spuriously
      ! positive and triggers a phantom internal freeze, moving liquid
      ! tracer into the already-empty ice slot.
      real(r8) :: d_wice_ext_soil1
      real(r8) :: transp_water_total, xylem_tracer_total, xylem_ratio
      real(r8) :: aquifer_ratio
      real(r8) :: transp_source_tracer_total, transp_output_tracer
      real(r8) :: transp_ratio, remove_ratio, R_vapor
      real(r8) :: source_fallback_ratio
      real(r8) :: relhum_leaf, leaf_area_use, rst_use, tleaf_use
      real(r8) :: leaf_delta_e_new, leaf_delta_b_new, leaf_peclet_new, leaf_moles_new
      real(r8) :: soil_resid_trc
      logical  :: transp_frac_active

      IF (ntracers <= 0) RETURN
      lb = snl + 1
      qcharge_eff = qcharge
      IF (abs(qcharge_eff) > 0.5_r8 * abs(spval)) qcharge_eff = 0._r8

      DO itrc = 1, ntracers
         R_precip = tracer_forcing_precip_value(itrc, ipatch)
         R_atm = tracer_forcing_vapor_value(itrc, ipatch)
         soil_resid_trc = 0._r8
         water_shadow(1:nl_soil) = wliq_soisno_bef(1:nl_soil)

         ! ============================================================
         ! 0. Compute pre-WATER tracer ratios for all soil layers.
         ! Dry-source fallbacks must not use current precipitation forcing:
         ! IsoGSM normalized ratios are noisy when precipitation is tiny,
         ! and root/aquifer uptake should inherit existing source water.
         ! ============================================================
         source_fallback_ratio = tracer_init_water_ratio(itrc)
         DO j = nl_soil, 1, -1
            IF (wliq_soisno_bef(j) > trc_tiny) THEN
               source_fallback_ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno_bef(j)
               EXIT
            ENDIF
         ENDDO
         DO j = 1, nl_soil
            IF (wliq_soisno_bef(j) > trc_tiny) THEN
               ratio_layer(j) = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno_bef(j)
            ELSE
	               ratio_layer(j) = source_fallback_ratio
	            ENDIF
	         ENDDO

	         ! ============================================================
	         ! 0a. Transpiration per soil layer.
         !
         !    WATER_VSF's soil_water_vertical_movement removes water for
         !    transpiration with a cascade: a dry shallow layer's unmet
         !    demand propagates to deeper layers, and any residual is
         !    absorbed by the aquifer via soilwater_aquifer_exchange
         !    (MOD_Hydro_SoilWater.F90:323-351). `etroot_actual(j)` is the
         !    mm of water actually drawn from layer j after that cascade;
         !    `etroot_aquifer` is the mm absorbed by the aquifer. Using
         !    these — instead of the raw demand etroot(j)*dt — matches the
         !    water side layer-by-layer and lets trc_wa absorb the deep
         !    ET tracer that would otherwise be silently orphaned.
         !    No-fractionation assumption: ratio_layer stays unchanged
         !    after this step, so section 2's qlayer transfer can still
         !    use the same ratio_layer.
	         ! ============================================================
		         aquifer_ratio = source_fallback_ratio
		         IF (abs(wa_bef) > trc_tiny) aquifer_ratio = trc_wa(itrc, ipatch) / wa_bef
		         transp_water_total = max(etroot_aquifer, 0._r8)
		         xylem_tracer_total = transp_water_total * aquifer_ratio
	         DO j = 1, nl_soil
	            IF (etroot_actual(j) > trc_tiny) THEN
	               transp_water_total = transp_water_total + etroot_actual(j)
	               xylem_tracer_total = xylem_tracer_total + etroot_actual(j) * ratio_layer(j)
	            ENDIF
	         ENDDO
		         xylem_ratio = source_fallback_ratio
		         IF (transp_water_total > trc_tiny) xylem_ratio = xylem_tracer_total / transp_water_total
			         transp_ratio = xylem_ratio
		         transp_source_tracer_total = 0._r8
		         transp_output_tracer = 0._r8
		         transp_frac_active = transp_water_total > trc_tiny .and. tracer_fractionation_active(itrc)
		         IF (transp_frac_active .and. &
		             present(tleaf_frac) .and. present(forc_q_frac) .and. present(forc_psrf_frac) .and. &
		             present(lai_frac) .and. present(rst_frac) .and. allocated(trc_leaf_delta_e)) THEN
	            tleaf_use = tleaf_frac
	            leaf_area_use = lai_frac
	            rst_use = rst_frac
	            IF (tracer_forcing_has_vapor(itrc, ipatch)) THEN
	               R_vapor = R_atm
	            ELSE
	               R_vapor = xylem_ratio / max(tracer_alpha_liq_vap(itrc, tleaf_use), trc_tiny)
	            ENDIF
	            relhum_leaf = tracer_surface_relhum(forc_q_frac, forc_psrf_frac, tleaf_use, .false.)
	            CALL tracer_transpiration_nss_ratio(itrc, xylem_ratio, R_vapor, &
	               tleaf_use, relhum_leaf, forc_psrf_frac, transp_water_total, deltim, leaf_area_use, rst_use, &
	               trc_leaf_delta_e(itrc, ipatch), trc_leaf_delta_b(itrc, ipatch), &
	               trc_leaf_peclet(itrc, ipatch), trc_leaf_water_moles(itrc, ipatch), &
	               transp_ratio, leaf_delta_e_new, leaf_delta_b_new, leaf_peclet_new, leaf_moles_new)
	            trc_leaf_delta_e(itrc, ipatch) = leaf_delta_e_new
	            trc_leaf_delta_b(itrc, ipatch) = leaf_delta_b_new
	            trc_leaf_peclet(itrc, ipatch) = leaf_peclet_new
	            trc_leaf_water_moles(itrc, ipatch) = leaf_moles_new
	         ENDIF

		         DO j = 1, nl_soil
		            IF (etroot_actual(j) > trc_tiny .and. wliq_soisno_bef(j) > trc_tiny) THEN
		               remove_ratio = ratio_layer(j)
		               trc_flux = etroot_actual(j) * remove_ratio
		               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
			               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
			               IF (transp_frac_active) THEN
			                  transp_source_tracer_total = transp_source_tracer_total + trc_flux
				               ELSE
					                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
					                     etroot_actual(j), TRC_EVAP_KIND_TRANSP)
					                  a_trc_transp_src(itrc, ipatch) = a_trc_transp_src(itrc, ipatch) + trc_flux
				               ENDIF
			               water_shadow(j) = water_shadow(j) - etroot_actual(j)
			            ENDIF
			         ENDDO

         ! Aquifer share of ET (deficit cascade residual). Use SIGNED
         ! wa_bef so a pre-existing aquifer debt (wa_bef<0) still yields
         ! a well-defined ratio; the min/max clamp prevents trc_wa from
         ! going more negative than its current signed value.
			         IF (etroot_aquifer > trc_tiny) THEN
			            remove_ratio = aquifer_ratio
			            trc_flux = etroot_aquifer * remove_ratio
			            IF (wa_bef > trc_tiny) trc_flux = min(trc_flux, max(trc_wa(itrc, ipatch), 0._r8))
			            trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - trc_flux
			            IF (transp_frac_active) THEN
			               transp_source_tracer_total = transp_source_tracer_total + trc_flux
			            ELSE
				               CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
				                  etroot_aquifer, TRC_EVAP_KIND_TRANSP)
				               a_trc_transp_src(itrc, ipatch) = a_trc_transp_src(itrc, ipatch) + trc_flux
			            ENDIF
		         ENDIF

			         IF (transp_frac_active) THEN
			            a_trc_transp_src(itrc, ipatch) = a_trc_transp_src(itrc, ipatch) &
			               + transp_source_tracer_total
			            ! NSS changes only the transpired vapor signature. Root uptake
		            ! still removes source/xylem water; the source-vapor difference
		            ! is held as leaf isotope isostorage to close tracer mass.
		            transp_output_tracer = transp_water_total * transp_ratio
			            CALL tracer_book_evap_loss(itrc, ipatch, transp_output_tracer, &
			               transp_water_total, TRC_EVAP_KIND_TRANSP)
		            IF (allocated(trc_leaf_iso_storage)) THEN
		               trc_leaf_iso_storage(itrc, ipatch) = trc_leaf_iso_storage(itrc, ipatch) &
		                  + transp_source_tracer_total - transp_output_tracer
		            ENDIF
		         ENDIF

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
            trc_flux = max(eff_qfros_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .true.) * deltim
            trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux

            ! Post-frost ice pool (water and tracer)
            water_ice_pool = wice_soisno_bef(lb_snow) + max(eff_qfros_snow, 0._r8) * deltim
            water_ice_pool_prefrost = water_ice_pool  ! save pre-subl value for Step 2 deficit check
            subl_water = max(eff_qsubl_snow, 0._r8) * deltim

            IF (subl_water > trc_tiny) THEN
               IF (water_ice_pool > trc_tiny .and. subl_water <= water_ice_pool) THEN
                  ! Normal case: enough ice to cover sublimation
	                  trc_flux = atmospheric_loss_tracer(trc_wice_soisno(itrc, lb_snow, ipatch), &
	                     water_ice_pool, subl_water, layer_temp(lb_snow), .true.)
	                  trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, subl_water, &
	                     TRC_EVAP_KIND_SUBL)
               ELSE
                  ! Ice insufficient (or empty): drain all ice, remainder from liquid
	                  ! Part 1: drain ice completely
	                  trc_flux = max(trc_wice_soisno(itrc, lb_snow, ipatch), 0._r8)
	                  trc_wice_soisno(itrc, lb_snow, ipatch) = 0._r8
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
	                     min(subl_water, max(water_ice_pool, 0._r8)), TRC_EVAP_KIND_SUBL)
                  ! Part 2: excess from liquid at LIQUID ratio
                  deficit_water = subl_water - max(water_ice_pool, 0._r8)
                  IF (deficit_water > trc_tiny .and. wliq_soisno_bef(lb_snow) > trc_tiny) THEN
	                     trc_flux = atmospheric_loss_tracer(trc_wliq_soisno(itrc, lb_snow, ipatch), &
	                        wliq_soisno_bef(lb_snow), deficit_water, layer_temp(lb_snow), .false.)
	                     trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
	                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, deficit_water, &
	                        TRC_EVAP_KIND_SUBL)
                  ENDIF
               ENDIF
            ENDIF
            ! Update water_ice_pool to post-subl value for use in Step 2's ice deficit
            water_ice_pool = max(water_ice_pool - subl_water, 0._r8)

            ! Record the ACTUAL external wice change for Section 5b.
            ! In no-deficit case this equals (qfros - qsubl)*dt as before;
            ! in deficit case the excess sublimation was drawn from liquid,
            ! so the net wice change is only -(wice_bef + qfros*dt).
            d_wice_ext_snow = max(eff_qfros_snow, 0._r8) * deltim &
                            - min(subl_water, water_ice_pool_prefrost)

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
            trc_flux = max(eff_qsdew_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .false.) * deltim
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
               ! Use PRE-subl ice pool to check if deficit actually occurred
               IF (water_ice_pool_prefrost < subl_water) THEN
                  water_liq_pool = water_liq_pool - (subl_water - water_ice_pool_prefrost)
               ENDIF
               water_liq_pool = max(water_liq_pool, 0._r8)
               evap_water = eff_qseva_snow * deltim

               IF (water_liq_pool > trc_tiny .and. evap_water <= water_liq_pool) THEN
                  ! Normal case: enough liquid
	                  trc_flux = atmospheric_loss_tracer(trc_wliq_soisno(itrc, lb_snow, ipatch), &
	                     water_liq_pool, evap_water, layer_temp(lb_snow), .false.)
	                  trc_wliq_soisno(itrc, lb_snow, ipatch) = trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, evap_water, &
	                     TRC_EVAP_KIND_SOILEVAP)
               ELSE
                  ! Liquid insufficient (or empty): drain liquid, remainder from ice
	                  ! Part 1: drain liquid completely
	                  trc_flux = max(trc_wliq_soisno(itrc, lb_snow, ipatch), 0._r8)
	                  trc_wliq_soisno(itrc, lb_snow, ipatch) = 0._r8
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
	                     min(evap_water, max(water_liq_pool, 0._r8)), TRC_EVAP_KIND_SOILEVAP)
	                  ! Part 2: excess from ice at post-frost-post-subl ICE ratio
	                  deficit_water = evap_water - max(water_liq_pool, 0._r8)
	                  IF (deficit_water > trc_tiny .and. water_ice_pool > trc_tiny) THEN
	                     deficit_water = min(deficit_water, water_ice_pool)
	                     d_wice_ext_snow = d_wice_ext_snow - deficit_water
	                     ! water_ice_pool is already post-frost, post-subl (updated at end of Step 1)
		                     trc_flux = atmospheric_loss_tracer(trc_wice_soisno(itrc, lb_snow, ipatch), &
		                        water_ice_pool, deficit_water, layer_temp(lb_snow), .true.)
		                     trc_wice_soisno(itrc, lb_snow, ipatch) = trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
		                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, deficit_water, &
		                        TRC_EVAP_KIND_SOILEVAP)
	                     water_ice_pool = water_ice_pool - deficit_water
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
         !    Negative qinfl has two physically different meanings:
         !      (a) qseva exceeds available surface input, so WATER extracts
         !          the deficit directly from layer 1 as evaporation output.
         !      (b) soil water truly exfiltrates into the surface pool.
         !    Only (b) should enter the surface mixed pool. Case (a) must be
         !    booked to a_trc_evap, otherwise active fractionation shows a
         !    storage loss equal to the negative-qinfl tracer diagnostic.
         ! ============================================================

         trc_soil_upflow = 0._r8
         trc_soil_evap = 0._r8
         top_exfil_water = 0._r8
         top_soil_evap_water = 0._r8

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
	         top_infil_water = max(qinfl, 0._r8) * deltim

         top_boundary_out_water = max(-qinfl, 0._r8) * deltim
         IF (top_boundary_out_water > trc_tiny) THEN
            ! qinfl = qgtop - d(wdsrf)/dt in WATER. Reconstruct qgtop from
            ! the final ponding change; for the evaporation-deficit cases
            ! that dominate active ptype=0 residuals, rsur is zero and this
            ! identifies negative qgtop exactly.
            qgtop_est = qinfl + (wdsrf - wdsrf_bef) / max(deltim, trc_tiny)

            IF (eff_qseva > trc_tiny .and. qgtop_est < -trc_tiny) THEN
               top_soil_evap_water = top_boundary_out_water
	               trc_soil_evap = atmospheric_loss_tracer(trc_wliq_soisno(itrc, 1, ipatch), &
	                  max(water_shadow(1), 0._r8), top_soil_evap_water, layer_temp(1), .false.)
	               trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_soil_evap
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_soil_evap, top_soil_evap_water, &
	                  TRC_EVAP_KIND_SOILEVAP)
               water_shadow(1) = water_shadow(1) - top_soil_evap_water
            ELSE
               top_exfil_water = top_boundary_out_water
               trc_soil_upflow = top_exfil_water * ratio_layer(1)
               trc_soil_upflow = min(trc_soil_upflow, max(trc_wliq_soisno(itrc, 1, ipatch), 0._r8))
               trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_soil_upflow
               water_shadow(1) = water_shadow(1) - top_exfil_water
            ENDIF
         ENDIF

	         ! Build pre-evaporation surface tracer pool.
         IF (snl < 0) THEN
            trc_pool_total = trc_wdsrf(itrc, ipatch) + trc_gwat_snow
            IF (split_soilsnow) THEN
               trc_pool_total = trc_pool_total + trc_pg_rain_ground(itrc, ipatch) * (1._r8 - fsno)
            ENDIF
         ELSE
            ! Only RAIN throughfall enters the surface mixed pool; snowfall
            ! throughfall is captured in trc_scv by tracer_newsnow (mirroring
            ! the water side where `gwat = pg_rain + sm - qseva` never
            ! contains pg_snow — see MOD_SoilSnowHydrology.F90:789,809).
            ! Using trc_pg_to_ground here would double-count pg_snow tracer
            ! (once in trc_scv, once in trc_pool_total), breaking
            ! conservation on snl==0 patches with simultaneous rain+snow.
            trc_pool_total = trc_wdsrf(itrc, ipatch) + trc_pg_rain_ground(itrc, ipatch)
         ENDIF

         ! WATER uses gwat = pg_rain + sm - qseva when snl>=0, so snowmelt
         ! participates in the same pre-evaporation surface pool as rain.
         IF (snl >= 0 .and. sm > trc_tiny) THEN
            trc_pool_total = trc_pool_total + trc_sm_carry(itrc, ipatch)
         ENDIF

         ! ---- Remove evaporation tracer from mixed pool ----
         ! gwat -= eff_qseva, so trc_pool must also be reduced.
         ! Use the ACTUAL pool ratio: trc_pool / (trc_pool / pool_ratio)
         ! where pool_ratio = trc_pool / water_pool_before_evap.
         ! Simplification: evap tracer = evap_water * (trc_pool / water_pool_before_evap)
         ! water_pool_before_evap ≈ trc_pool / pool_ratio. We don't know the exact
         ! water pool size at this point, but we can use:
         !   pool_ratio = trc_pool / water_pool_before_evap
         !   trc_evap = surface_evap_water * pool_ratio
         !            = surface_evap_water * trc_pool / water_bef
         ! water_bef = water_pool_after_evap + surface_evap_water
         ! If qseva exceeds the surface input/ponding supply, WATER_VSF
         ! draws the deficit directly from layer 1 as negative qgtop/qinfl;
         ! that deficit is already booked above as top_soil_evap_water.
         ! Only the remaining qseva can carry the surface mixed-pool
         ! signature. Using the full qseva here double-counts water and
         ! dilutes the surface tracer flux delta.
         gwat_evap = max(max(eff_qseva, 0._r8) * deltim - top_soil_evap_water, 0._r8)
         IF (gwat_evap > trc_tiny .and. trc_pool_total > trc_tiny) THEN
            ! Reconstruct pre-evap water pool
            water_pool_total = max(wdsrf, 0._r8) + max(rsur, 0._r8) * deltim &
                             + top_infil_water + gwat_evap
	            trc_gwat_evap = atmospheric_loss_tracer(trc_pool_total, water_pool_total, &
	               gwat_evap, layer_temp(1), .false.)
	            trc_pool_total = trc_pool_total - trc_gwat_evap
	            CALL tracer_book_evap_loss(itrc, ipatch, trc_gwat_evap, gwat_evap, &
	               TRC_EVAP_KIND_SOILEVAP)
         ENDIF

         ! Soil exfiltration is produced by the VSF solve after qseva has
         ! already been applied to gwat, so it should feed runoff/wdsrf but
         ! should not be available to this step's ground evaporation.
         trc_pool_total = trc_pool_total + trc_soil_upflow

         ! ---- Add drip/flood/paddy irrigation tracer ----
         ! Water side (MOD_SoilSnowHydrology.F90:828-832) adds the
         ! irrigation mass to gwat AFTER qseva, so none of it is
         ! evaporated at this step. Irrigation is an INTERNAL transfer
         ! from the per-patch waterstorage reservoir to the surface
         ! pool — water-side totwb/endwb include waterstorage
         ! (CoLMMAIN:806), so the atmospheric input term does not count
         ! irrigation. The tracer path mirrors this: debit
         ! trc_waterstorage at its current ratio (Phase-1 invariant
         ! keeps it at R_init, so we use R_precip as a proxy) and do
         ! NOT add to a_trc_precip. Without this mirroring the tracer
         ! budget claims more atmospheric input than the water budget.
         IF (qflx_irrig_ground > trc_tiny) THEN
            storage_ratio = R_precip
            IF (allocated(trc_waterstorage) .and. present(waterstorage_patch)) THEN
               IF (waterstorage_patch > trc_tiny) THEN
                  storage_ratio = trc_waterstorage(itrc, ipatch) / max(waterstorage_patch, trc_tiny)
               ENDIF
            ENDIF
            trc_flux = qflx_irrig_ground * deltim * storage_ratio
            IF (allocated(trc_waterstorage)) THEN
               trc_flux = min(max(trc_flux, 0._r8), max(trc_waterstorage(itrc, ipatch), 0._r8))
            ENDIF
            trc_pool_total = trc_pool_total + trc_flux
            IF (allocated(trc_waterstorage)) THEN
               trc_waterstorage(itrc, ipatch) = trc_waterstorage(itrc, ipatch) - trc_flux
            ENDIF
         ENDIF

         ! Use OUTPUT-side water pool for ratio computation.
         ! This guarantees: distributed = water_pool * ratio = trc_pool (exact).
         water_pool_total = max(wdsrf, 0._r8) + max(rsur, 0._r8) * deltim &
                          + top_infil_water

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
            trc_flux = rsur * ratio * deltim
            a_trc_rsur(itrc, ipatch) = a_trc_rsur(itrc, ipatch) + trc_flux
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_flux
            trc_rnof_step(itrc, ipatch) = trc_rnof_step(itrc, ipatch) + trc_flux
         ENDIF

         IF (qinfl > trc_tiny) THEN
            a_trc_qinfl(itrc, ipatch) = a_trc_qinfl(itrc, ipatch) + qinfl * ratio * deltim
         ELSEIF (qinfl < -trc_tiny) THEN
            a_trc_qinfl(itrc, ipatch) = a_trc_qinfl(itrc, ipatch) + qinfl * ratio_layer(1) * deltim
         ENDIF

         ! ============================================================
         ! 2. Soil layers: qlayer-driven flux tracking
         !
         !    qlayer(j) > 0: downward from j to j+1, source = layer j
         !    qlayer(j) < 0: upward from j+1 to j, source = layer j+1
         !
         !    Each flux removes from source, adds to destination.
         ! ============================================================

         ! --- surface → layer 1 infiltration (downward only) ---
         ! WATER_VSF reports this as qinfl = qgtop - d(wdsrf)/dt. qlayer(0)
         ! can be zero when the top-layer storage changes through the VSF
         ! wetting-front representation, so using qlayer(0) here misses the
         ! actual water entering layer 1.
	         IF (qinfl > trc_tiny) THEN
	            trc_flux = qinfl * ratio * deltim
	            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) + trc_flux
	            water_shadow(1) = water_shadow(1) + qinfl * deltim
	         ENDIF

         ! --- qlayer(j): layer j ↔ layer j+1 ---
         DO j = 1, nl_soil - 1
            IF (qlayer(j) > trc_tiny) THEN
               ! Downward: j → j+1
               trc_flux = qlayer(j) * ratio_layer(j) * deltim
	               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
	               trc_wliq_soisno(itrc, j,   ipatch) = trc_wliq_soisno(itrc, j,   ipatch) - trc_flux
	               trc_wliq_soisno(itrc, j+1, ipatch) = trc_wliq_soisno(itrc, j+1, ipatch) + trc_flux
	               water_shadow(j)   = water_shadow(j)   - qlayer(j) * deltim
	               water_shadow(j+1) = water_shadow(j+1) + qlayer(j) * deltim
	            ELSEIF (qlayer(j) < -trc_tiny) THEN
	               ! Upward: j+1 → j
	               trc_flux = abs(qlayer(j)) * ratio_layer(j+1) * deltim
	               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j+1, ipatch), 0._r8))
	               trc_wliq_soisno(itrc, j+1, ipatch) = trc_wliq_soisno(itrc, j+1, ipatch) - trc_flux
	               trc_wliq_soisno(itrc, j,   ipatch) = trc_wliq_soisno(itrc, j,   ipatch) + trc_flux
	               water_shadow(j+1) = water_shadow(j+1) + qlayer(j) * deltim
	               water_shadow(j)   = water_shadow(j)   - qlayer(j) * deltim
	            ENDIF
	         ENDDO

         ! ============================================================
         ! 3. Groundwater: qcharge-driven
         ! ============================================================
         IF (qcharge_eff > trc_tiny) THEN
            j = nl_soil
            trc_flux = qcharge_eff * ratio_layer(j) * deltim
	            trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
	            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
	            trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) + trc_flux
	            a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) + trc_flux
	            water_shadow(j) = water_shadow(j) - qcharge_eff * deltim
	         ELSEIF (qcharge_eff < -trc_tiny) THEN
            IF (abs(wa_bef) > trc_tiny) THEN
               ratio_src = trc_wa(itrc, ipatch) / wa_bef
            ELSE
               ratio_src = source_fallback_ratio
            ENDIF
            trc_flux = abs(qcharge_eff) * ratio_src * deltim
            IF (wa_bef > trc_tiny) THEN
               trc_flux = min(trc_flux, max(trc_wa(itrc, ipatch), 0._r8))
            ENDIF
            trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - trc_flux
	            j = nl_soil
	            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
	            a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) - trc_flux
	            water_shadow(j) = water_shadow(j) - qcharge_eff * deltim
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
            trc_rnof_step(itrc, ipatch) = trc_rnof_step(itrc, ipatch) + trc_flux
            water_shadow(j) = water_shadow(j) - rsub * deltim
         ENDIF

         ! ============================================================
         ! 4b. WATER-phase liquid external flux on soil layer 1: dew
         !     qsdew / qsdew_soil are applied directly to wliq_soisno(1)
         !     by WATER, not to the surface mixed pool.
         ! ============================================================
         IF (eff_qsdew_topliq > trc_tiny) THEN
	            trc_flux = eff_qsdew_topliq * deposition_ratio_for(layer_temp(1), .false.) * deltim
	            trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) + trc_flux
	            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
	            water_shadow(1) = water_shadow(1) + eff_qsdew_topliq * deltim
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
         ! Record the actual external wice change on soil layer 1 BEFORE
         ! mutating trc_wice, so Section 5b can subtract exactly what the
         ! water side did (including the max(0,...) clamp). Formula:
         !   qfros_applied = max(qfros, 0) * dt            (always adds)
         !   qsubl_actual  = min(max(qsubl,0)*dt, wice_bef + qfros_applied)
         !   d_wice_ext    = qfros_applied - qsubl_actual
         ! Normal case: qsubl*dt <= wice_bef + qfros*dt  →  d_wice_ext = (qfros-qsubl)*dt
         ! Clamped case: qsubl*dt  > wice_bef + qfros*dt  →  d_wice_ext = -wice_bef
         d_wice_ext_soil1 = max(eff_qfros_top, 0._r8) * deltim &
            - min(max(eff_qsubl_top, 0._r8) * deltim, &
                  max(wice_soisno_bef(1), 0._r8) + max(eff_qfros_top, 0._r8) * deltim)

         ! Use effective fluxes (handles both split and non-split paths)
         IF (eff_qfros_top > trc_tiny) THEN
            trc_flux = eff_qfros_top * deposition_ratio_for(layer_temp(1), .true.) * deltim
            trc_wice_soisno(itrc, 1, ipatch) = trc_wice_soisno(itrc, 1, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF
         IF (eff_qsubl_top > trc_tiny) THEN
            ! trc_wice was just updated by frost above, so the right
            ! denominator is the post-frost ice pool, not wice_soisno_bef(1).
            wice_pre_phase = max(wice_soisno_bef(1) + max(eff_qfros_top, 0._r8) * deltim, 0._r8)
            subl_water = eff_qsubl_top * deltim
            ! Mirror the snow-top deficit path. Water side
            ! (MOD_SoilSnowHydrology.F90:1101) clamps wice at 0 when
            ! qsubl*dt > wice_bef+qfros*dt and charges the excess
            ! against wliq. Previously the tracer side only min()'d
            ! against trc_wice, so the deficit vapour carried away no
            ! tracer and a_trc_evap was systematically short in
            ! sublimation-heavy patches.
            IF (wice_pre_phase > trc_tiny .and. subl_water <= wice_pre_phase) THEN
               ! Normal case: ice covers sublimation.
	               trc_flux = atmospheric_loss_tracer(trc_wice_soisno(itrc, 1, ipatch), &
	                  wice_pre_phase, subl_water, layer_temp(1), .true.)
	               trc_wice_soisno(itrc, 1, ipatch) = trc_wice_soisno(itrc, 1, ipatch) - trc_flux
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, subl_water, &
	                  TRC_EVAP_KIND_SUBL)
            ELSEIF (subl_water > trc_tiny) THEN
               ! Deficit: drain ice completely, then pull the remainder
               ! from liquid at pre-WATER wliq ratio (ratio_layer(1) was
               ! cached at the top of the step before any layer-1 mutation).
	               trc_flux = max(trc_wice_soisno(itrc, 1, ipatch), 0._r8)
	               trc_wice_soisno(itrc, 1, ipatch) = 0._r8
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
	                  min(subl_water, max(wice_pre_phase, 0._r8)), TRC_EVAP_KIND_SUBL)
               deficit_water = subl_water - max(wice_pre_phase, 0._r8)
	               IF (deficit_water > trc_tiny .and. wliq_soisno_bef(1) > trc_tiny) THEN
		                  trc_flux = atmospheric_loss_tracer(trc_wliq_soisno(itrc, 1, ipatch), &
	                        max(water_shadow(1), 0._r8), deficit_water, layer_temp(1), .false.)
		                  trc_wliq_soisno(itrc, 1, ipatch) = trc_wliq_soisno(itrc, 1, ipatch) - trc_flux
		                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, deficit_water, &
		                     TRC_EVAP_KIND_SUBL)
	                  water_shadow(1) = water_shadow(1) - deficit_water
	               ENDIF
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
         !
         ! ACCURACY LIMIT (H8 reconstruction):
         !     The denominators below are reconstructed from POST-WATER state:
         !         pre-freeze wliq = wliq_soisno(j) + d_wice
         !         pre-thaw  wice = wice_soisno(j) + |d_wice|
         !     This recovers the exact pre-phase pool size ONLY if freeze/thaw
         !     was the LAST operation in WATER's substep order. If WATER does
         !     percolation AFTER phase change in some layer, the reconstruction
         !     under/over-states the actual pool at phase-change time.
         !
         !     Non-fractionating tracers still use conservative ratio
         !     transfer. Fractionating isotope tracers now use a Rayleigh
         !     finite-pool freezing update with this reconstructed pre-freeze
         !     liquid pool. That is exact only if the reconstruction matches
         !     WATER's phase-change intermediate state. If a future WATER_VSF
         !     moves phase change before/among flow steps, pass explicit
         !     phase-change masses or move this tracer logic into WATER_VSF.
         !     Thaw remains conservative melting of existing ice.
         ! ============================================================
         DO j = lb, nl_soil
            d_wice = wice_soisno(j) - wice_soisno_bef(j)

            ! Subtract external ice fluxes to isolate internal freeze/thaw
            IF (j == 1) THEN
               ! Soil layer 1: use the ACTUAL external change recorded above
               ! in Section 5. Using the nominal (qfros-qsubl)*dt would
               ! misattribute a `max(0, ...)` clamp residual to internal
               ! freezing under extreme sublimation (qsubl*dt > wice_bef+qfros*dt).
               d_wice = d_wice - d_wice_ext_soil1
            ENDIF
            IF (j < 1 .and. j == lb .and. snl < 0) THEN
               ! Top snow layer: external ice flux is applied in snowwater.
               ! Use actual external change recorded in Section 0b, which
               ! correctly accounts for qsubl deficit drawn from liquid.
               d_wice = d_wice - d_wice_ext_snow
            ENDIF

            ! wblc reconciliation (WATER_VSF:1069-1081) pulls ice mass to
            ! cover the ET liquid-water deficit; that water leaves the
            ! column as ET OUTPUT, not as internal thaw. Remove the
            ! matching tracer from trc_wice at the pre-wblc ice ratio and
            ! book it into a_trc_evap. Then add wblc_ice_sink back into
            ! d_wice so the freeze/thaw classification below sees only
            ! the true internal phase-change residual (≈ 0 under the
            ! current WATER_VSF, which performs no in-place phase change
            ! on soil layers — this keeps the branch ready for future
            ! WATER versions that do).
            IF (j >= 1) THEN
               IF (wblc_ice_sink(j) > trc_tiny) THEN
                  wice_pre_phase = max(wice_soisno(j) + wblc_ice_sink(j), 0._r8)
                  IF (wice_pre_phase > trc_tiny) THEN
                     ratio_src = trc_wice_soisno(itrc, j, ipatch) / wice_pre_phase
	                     trc_flux = wblc_ice_sink(j) * ratio_src
	                     trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
	                     trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
	                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, wblc_ice_sink(j), &
	                        TRC_EVAP_KIND_SOILEVAP)
                  ENDIF
                  d_wice = d_wice + wblc_ice_sink(j)
               ENDIF
            ENDIF

            ! Use the pool size at freeze/thaw time, not pre-WATER. By the time
            ! we reach this block trc_wliq has been updated by percolation /
            ! qcharge / rsub / dew, so wliq_soisno_bef no longer matches the
            ! tracer pool that participates in the phase change. Reconstruct
            ! pre-freeze wliq = post-WATER wliq + d_wice (the part that froze
            ! out of liquid into ice); pre-thaw wice = post-WATER wice + |d_wice|.
            IF (d_wice > trc_tiny) THEN
               ! Internal freeze: liquid → ice
               wliq_pre_phase = max(wliq_soisno(j) + d_wice, 0._r8)
               IF (wliq_pre_phase > trc_tiny) THEN
                  trc_flux = tracer_rayleigh_freezing_loss(itrc, &
                     trc_wliq_soisno(itrc, j, ipatch), wliq_pre_phase, &
                     d_wice, layer_temp(j))
                  trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
                  IF (j >= 1) water_shadow(j) = water_shadow(j) - d_wice
               ENDIF
            ELSEIF (d_wice < -trc_tiny) THEN
               ! Internal thaw: ice → liquid
               wice_pre_phase = max(wice_soisno(j) - d_wice, 0._r8)   ! d_wice < 0, so -d_wice > 0
               IF (wice_pre_phase > trc_tiny) THEN
                  ratio_src = trc_wice_soisno(itrc, j, ipatch) / wice_pre_phase
                  trc_flux = abs(d_wice) * ratio_src
	                  trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
	                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
	                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
	                  IF (j >= 1) water_shadow(j) = water_shadow(j) - d_wice
	               ENDIF
	            ENDIF
	         ENDDO

	         ! The VSF solver closes each layer only to its nonlinear tolerance.
	         ! Carry that numerical storage residual at the layer's current
	         ! concentration so it cannot appear as isotope fractionation.
	         DO j = 1, nl_soil
	            water_resid = wliq_soisno(j) - water_shadow(j)
	            IF (abs(water_resid) > trc_tiny) THEN
	               IF (water_shadow(j) > trc_tiny) THEN
	                  water_shadow_ratio = trc_wliq_soisno(itrc, j, ipatch) / water_shadow(j)
	               ELSE
	                  water_shadow_ratio = source_fallback_ratio
	               ENDIF
	               IF (water_resid >= 0._r8) THEN
	                  trc_flux = water_resid * water_shadow_ratio
	                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
	                  soil_resid_trc = soil_resid_trc + trc_flux
	               ELSE
	                  trc_flux = min(-water_resid * water_shadow_ratio, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
	                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
	                  soil_resid_trc = soil_resid_trc - trc_flux
	               ENDIF
	            ENDIF
	         ENDDO

		         ! The VSF layer residual is a numerical closure term, not a
		         ! physical aquifer / runoff / evaporation flux. Keep it out of
		         ! trc_wa and expose it to the conservation checker as an
		         ! explicit numerical source/sink.
			         IF (abs(soil_resid_trc) > trc_tiny) THEN
			            IF (allocated(trc_numerical_residual_step)) THEN
			               trc_numerical_residual_step(itrc, ipatch) = &
			                  trc_numerical_residual_step(itrc, ipatch) + soil_resid_trc
			            ELSE
			               trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - soil_resid_trc
			            ENDIF
			         ENDIF
				         IF (abs(wa) <= trc_tiny) THEN
				            IF (abs(trc_wa(itrc, ipatch)) > trc_tiny) THEN
				               ! With exactly zero aquifer water there is no signed
				               ! storage/debt state to mirror. Keep conservation
				               ! explicit, but do not persist an orphan trc_wa.
				               IF (allocated(trc_numerical_residual_step)) THEN
				                  trc_numerical_residual_step(itrc, ipatch) = &
				                     trc_numerical_residual_step(itrc, ipatch) &
				                     - trc_wa(itrc, ipatch)
				               ENDIF
				            ENDIF
				            trc_wa(itrc, ipatch) = 0._r8
				         ENDIF

		         ! ============================================================
		         ! 6. Wetland water: standard soil path should not mutate it.
	         ! ============================================================
	         d_wetwat = wetwat - wetwat_bef
#if (defined CoLMDEBUG)
	         IF (abs(d_wetwat) > trc_tiny) THEN
	            write(*,'(A,I8,A,I3,A,E12.5,A,E12.5,A,E12.5)') &
	               ' WARNING tracer_soil_water: wetwat changed in standard soil path ipatch=', &
	               ipatch, ' itrc=', itrc, ' d_wetwat=', d_wetwat, &
	               ' wetwat_bef=', wetwat_bef, ' wetwat=', wetwat
	         ENDIF
#endif

	      ENDDO

   CONTAINS

      real(r8) FUNCTION layer_temp (jlay)
         integer, intent(in) :: jlay

         layer_temp = 273.15_r8
         IF (present(t_soisno_frac)) layer_temp = t_soisno_frac(jlay)
      END FUNCTION layer_temp

      real(r8) FUNCTION atmospheric_loss_tracer (pool_trc, pool_water, water_loss, temp_k, from_ice)
         real(r8), intent(in) :: pool_trc
         real(r8), intent(in) :: pool_water
         real(r8), intent(in) :: water_loss
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8), parameter :: max_evap_frac_step = 0.05_r8
         real(r8) :: source_ratio, flux_ratio
         real(r8) :: remaining_water, remaining_trc, loss_left, loss_step, trc_loss_step

         atmospheric_loss_tracer = 0._r8
         IF (pool_water <= trc_tiny .or. water_loss <= trc_tiny) RETURN
         IF (pool_trc <= trc_tiny) RETURN

         IF (water_loss >= pool_water * (1._r8 - 1.e-12_r8)) THEN
            atmospheric_loss_tracer = max(pool_trc, 0._r8)
            RETURN
         ENDIF

         source_ratio = max(pool_trc, 0._r8) / pool_water
         flux_ratio = evap_ratio_for(source_ratio, temp_k, from_ice)
         IF (water_loss <= max_evap_frac_step * pool_water .or. &
             abs(flux_ratio - source_ratio) <= 1.e-12_r8 * max(source_ratio, trc_tiny)) THEN
            atmospheric_loss_tracer = min(water_loss * flux_ratio, max(pool_trc, 0._r8))
            RETURN
         ENDIF

         ! Finite-pool evaporation can remove nearly all water in one land step.
         ! Integrating the loss prevents a single light vapor ratio from leaving
         ! an unrealistically enriched trace residual.
         remaining_water = pool_water
         remaining_trc = max(pool_trc, 0._r8)
         loss_left = water_loss
         DO WHILE (loss_left > trc_tiny .and. remaining_water > trc_tiny .and. &
                   remaining_trc > trc_tiny)
            loss_step = min(loss_left, max_evap_frac_step * remaining_water)
            source_ratio = remaining_trc / remaining_water
            flux_ratio = evap_ratio_for(source_ratio, temp_k, from_ice)
            trc_loss_step = min(loss_step * flux_ratio, remaining_trc)
            remaining_trc = remaining_trc - trc_loss_step
            remaining_water = remaining_water - loss_step
            loss_left = loss_left - loss_step
         ENDDO

         atmospheric_loss_tracer = min(max(pool_trc, 0._r8), &
            max(pool_trc, 0._r8) - max(remaining_trc, 0._r8))
      END FUNCTION atmospheric_loss_tracer

      real(r8) FUNCTION evap_ratio_for (source_ratio, temp_k, from_ice)
         real(r8), intent(in) :: source_ratio
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8) :: relhum, alpha_k

         evap_ratio_for = source_ratio
         IF (.not. tracer_fractionation_active(itrc)) RETURN
         IF (.not. present(forc_q_frac) .or. .not. present(forc_psrf_frac)) RETURN

         relhum = tracer_surface_relhum(forc_q_frac, forc_psrf_frac, temp_k, from_ice)
         alpha_k = tracer_diffusivity_ratio_air(itrc)
         evap_ratio_for = tracer_craig_gordon_evap_ratio(itrc, source_ratio, R_atm, &
            temp_k, relhum, alpha_k, from_ice)
         evap_ratio_for = min(evap_ratio_for, max(source_ratio, 0._r8))
      END FUNCTION evap_ratio_for

      real(r8) FUNCTION deposition_ratio_for (temp_k, from_ice)
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         deposition_ratio_for = R_atm
         IF (.not. tracer_fractionation_active(itrc)) RETURN
         deposition_ratio_for = tracer_equilibrium_deposition_ratio(itrc, R_atm, temp_k, from_ice)
      END FUNCTION deposition_ratio_for

   END SUBROUTINE tracer_soil_water

   !---------------------------------------------------------------
   ! Wetland-specific tracer update (patchtype==2 .and.
   ! .not. DEF_USE_Dynamic_Wetland).
   !
   ! WATER_VSF's wetland branch (MOD_SoilSnowHydrology.F90 L1170+) merges
   ! wdsrf + wa + wetwat + (gwat-etr+dew+frost-subl)*dt + sum(wresi) into a
   ! single wetwat pool, then overflows back to wdsrf / deficits into wa.
   ! The normal soil-water tracer path (tracer_soil_water) assumes an
   ! independent wdsrf/wa/wetwat accounting and silently loses the
   ! pool-to-wetwat mass transfer. This routine mirrors the water merge
   ! exactly:
   !   pool_water  = wdsrf_bef + wa_bef + wetwat_bef + wresi + inputs
   !   pool_tracer = trc_wdsrf + trc_wa + trc_wetwat + trc_wresi + trc_inputs
   !   loss        = (qseva + qsubl + etr)*dt  removed at pool ratio
   !   remaining   = redistributed to post-WATER wdsrf / wa / wetwat / rsur
   !
   ! Snow layers (lb<0) are handled up front with the same top-layer flux +
   ! percolation logic as tracer_soil_water sections 0b/0c, producing
   ! trc_gwat_snow which joins the wetland pool instead of the surface pool.
   !---------------------------------------------------------------
   SUBROUTINE tracer_wetland (ipatch, deltim, snl, nl_soil, &
      rsur, &
      qseva_in, qsdew_in, qsubl_in, qfros_in, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
      etr, sm, fsno, split_soilsnow, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_bef, wice_soisno_bef, &
      wa, wa_bef, wdsrf, wdsrf_bef, &
      wetwat, wetwat_bef, pg_rain, pg_snow, &
      t_soisno, porsl, dz_soisno, qflx_irrig_ground, waterstorage_patch, &
      forc_q_frac, forc_psrf_frac)

      USE MOD_Const_Physical, only: tfrz
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, nl_soil
      real(r8), intent(in) :: deltim, rsur
      real(r8), intent(in) :: qseva_in, qsdew_in, qsubl_in, qfros_in
      real(r8), intent(in) :: qseva_soil, qsdew_soil, qsubl_soil, qfros_soil
      real(r8), intent(in) :: qseva_snow, qsdew_snow, qsubl_snow, qfros_snow
      real(r8), intent(in) :: etr, sm, fsno
      logical,  intent(in) :: split_soilsnow
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_bef(snl+1:nl_soil)
      real(r8), intent(in) :: wa, wa_bef, wdsrf, wdsrf_bef
      real(r8), intent(in) :: wetwat, wetwat_bef, pg_rain, pg_snow
      real(r8), intent(in) :: t_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: porsl(1:nl_soil), dz_soisno(snl+1:nl_soil)
      ! Combined drip + flood + paddy irrigation rate (mm/s). Added to
      ! gwat in WATER_VSF (MOD_SoilSnowHydrology.F90:828-832) ahead of
      ! the wetland-merge branch, so wetwat_post already reflects this
      ! water. Without mirroring it here, pool_tracer is short by
      ! `qflx_irrig_ground*dt*R_precip` while pool_water is short by
      ! the same volume — the ratio stays roughly right but the
      ! absolute tracer mass shrinks and xerr_tracer grows by the
      ! irrigation input per step.
      real(r8), intent(in) :: qflx_irrig_ground
      real(r8), intent(in), optional :: waterstorage_patch
      real(r8), intent(in), optional :: forc_q_frac
      real(r8), intent(in), optional :: forc_psrf_frac

      integer  :: itrc, j, lb, lb_snow
      real(r8) :: R_atm, ratio_src, trc_flux, storage_ratio
      real(r8) :: wresi_j, wresi_sum, trc_wresi_j, trc_wresi_sum
      real(r8) :: q_rain_in, q_dew_in, q_frost_in
      real(r8) :: q_evap_out, q_subl_out, q_etr_out, q_sm_in
      real(r8) :: pool_water, pool_tracer, pool_ratio
      real(r8) :: trc_loss, loss_water, trc_rsur_local
      real(r8) :: trc_dew_input, trc_frost_input, trc_evap_loss, trc_subl_loss, trc_etr_loss
      real(r8) :: redist_storage, redist_target, redist_resid, redist_fix
      real(r8) :: trc_gwat_snow_local
      ! Water mass (mm) leaving the snow bottom during percolation — the
      ! hydrology side feeds this into wetwat via `gwat*deltim`. Tracer
      ! pool_water must mirror it so the mixed-pool ratio stays consistent
      ! with the water balance (otherwise trc_gwat_snow_local inflates
      ! pool_tracer without a matching pool_water contribution).
      real(r8) :: gwat_snow_local
      real(r8) :: eff_qseva_snow, eff_qsdew_snow, eff_qsubl_snow, eff_qfros_snow
      real(r8) :: water_ice_pool, water_ice_pool_prefrost, water_liq_pool
      real(r8) :: subl_water, evap_water, deficit_water
      real(r8) :: snow_rain_input, snow_dew_input, snow_evap_output
      real(r8) :: snow_frost_input, snow_subl_output, snow_wliq_before_flow
      real(r8) :: qin_snow, qout_snow, water_before_flow, trc_before_flow
      real(r8) :: trc_qin_snow, trc_qout_snow
      logical  :: fractionate_pool_loss
      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         R_atm = tracer_forcing_vapor_value(itrc, ipatch)

         !--------------------------------------------------------
         ! 1) Snow-top external fluxes + percolation → trc_gwat_snow
         !    (mirror of tracer_soil_water sections 0b / 0c)
         !--------------------------------------------------------
         trc_gwat_snow_local = 0._r8
         gwat_snow_local     = 0._r8
         IF (snl < 0) THEN
            lb_snow = lb

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

            ! Ice phase: frost in, sublim out
            trc_flux = max(eff_qfros_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .true.) * deltim
            trc_wice_soisno(itrc, lb_snow, ipatch) = &
               trc_wice_soisno(itrc, lb_snow, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux

            water_ice_pool = wice_soisno_bef(lb_snow) + max(eff_qfros_snow, 0._r8) * deltim
            water_ice_pool_prefrost = water_ice_pool
            subl_water = max(eff_qsubl_snow, 0._r8) * deltim

            IF (subl_water > trc_tiny) THEN
               IF (water_ice_pool > trc_tiny .and. subl_water <= water_ice_pool) THEN
	                  trc_flux = atmospheric_loss_tracer(trc_wice_soisno(itrc, lb_snow, ipatch), &
	                     water_ice_pool, subl_water, layer_temp(lb_snow), .true.)
	                  trc_wice_soisno(itrc, lb_snow, ipatch) = &
	                     trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, subl_water, &
	                     TRC_EVAP_KIND_SUBL)
               ELSE
	                  trc_flux = max(trc_wice_soisno(itrc, lb_snow, ipatch), 0._r8)
	                  trc_wice_soisno(itrc, lb_snow, ipatch) = 0._r8
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
	                     min(subl_water, max(water_ice_pool, 0._r8)), TRC_EVAP_KIND_SUBL)
                  deficit_water = subl_water - max(water_ice_pool, 0._r8)
                  IF (deficit_water > trc_tiny .and. wliq_soisno_bef(lb_snow) > trc_tiny) THEN
	                     trc_flux = atmospheric_loss_tracer(trc_wliq_soisno(itrc, lb_snow, ipatch), &
	                        wliq_soisno_bef(lb_snow), deficit_water, layer_temp(lb_snow), .false.)
	                     trc_wliq_soisno(itrc, lb_snow, ipatch) = &
	                        trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
	                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, deficit_water, &
	                        TRC_EVAP_KIND_SUBL)
                  ENDIF
               ENDIF
            ENDIF
            water_ice_pool = max(water_ice_pool - subl_water, 0._r8)

            ! Liquid phase: rain + dew in, evap out
            IF (split_soilsnow) THEN
               trc_flux = trc_pg_rain_ground(itrc, ipatch) * fsno
            ELSE
               trc_flux = trc_pg_rain_ground(itrc, ipatch)
            ENDIF
            trc_wliq_soisno(itrc, lb_snow, ipatch) = &
               trc_wliq_soisno(itrc, lb_snow, ipatch) + trc_flux

            trc_flux = max(eff_qsdew_snow, 0._r8) * &
               deposition_ratio_for(layer_temp(lb_snow), .false.) * deltim
            trc_wliq_soisno(itrc, lb_snow, ipatch) = &
               trc_wliq_soisno(itrc, lb_snow, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux

            IF (eff_qseva_snow > trc_tiny) THEN
               IF (split_soilsnow) THEN
                  water_liq_pool = wliq_soisno_bef(lb_snow) &
                     + (pg_rain*fsno + max(eff_qsdew_snow, 0._r8)) * deltim
               ELSE
                  water_liq_pool = wliq_soisno_bef(lb_snow) &
                     + (pg_rain + max(eff_qsdew_snow, 0._r8)) * deltim
               ENDIF
               IF (water_ice_pool_prefrost < subl_water) THEN
                  water_liq_pool = water_liq_pool - (subl_water - water_ice_pool_prefrost)
               ENDIF
               water_liq_pool = max(water_liq_pool, 0._r8)
               evap_water = eff_qseva_snow * deltim

               IF (water_liq_pool > trc_tiny .and. evap_water <= water_liq_pool) THEN
	                  trc_flux = atmospheric_loss_tracer(trc_wliq_soisno(itrc, lb_snow, ipatch), &
	                     water_liq_pool, evap_water, layer_temp(lb_snow), .false.)
	                  trc_wliq_soisno(itrc, lb_snow, ipatch) = &
	                     trc_wliq_soisno(itrc, lb_snow, ipatch) - trc_flux
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, evap_water, &
	                     TRC_EVAP_KIND_SOILEVAP)
               ELSE
	                  trc_flux = max(trc_wliq_soisno(itrc, lb_snow, ipatch), 0._r8)
	                  trc_wliq_soisno(itrc, lb_snow, ipatch) = 0._r8
	                  CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, &
	                     min(evap_water, max(water_liq_pool, 0._r8)), TRC_EVAP_KIND_SOILEVAP)
                  deficit_water = evap_water - max(water_liq_pool, 0._r8)
                  IF (deficit_water > trc_tiny .and. water_ice_pool > trc_tiny) THEN
	                     trc_flux = atmospheric_loss_tracer(trc_wice_soisno(itrc, lb_snow, ipatch), &
	                        water_ice_pool, deficit_water, layer_temp(lb_snow), .true.)
	                     trc_wice_soisno(itrc, lb_snow, ipatch) = &
	                        trc_wice_soisno(itrc, lb_snow, ipatch) - trc_flux
	                     CALL tracer_book_evap_loss(itrc, ipatch, trc_flux, deficit_water, &
	                        TRC_EVAP_KIND_SOILEVAP)
                  ENDIF
               ENDIF
            ENDIF

            ! Snow percolation → trc_gwat_snow_local
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

            trc_qin_snow = 0._r8
            snow_wliq_before_flow = wliq_soisno_bef(lb_snow)
            IF (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output < 0._r8) THEN
               snow_wliq_before_flow = snow_wliq_before_flow &
                  + (wice_soisno_bef(lb_snow) + snow_frost_input - snow_subl_output)
            ENDIF
            snow_wliq_before_flow = snow_wliq_before_flow &
               + snow_rain_input + snow_dew_input - snow_evap_output
            snow_wliq_before_flow = max(snow_wliq_before_flow, 0._r8)

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

               IF (qout_snow > trc_tiny .and. water_before_flow > trc_tiny &
                   .and. trc_before_flow > trc_tiny) THEN
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

            trc_gwat_snow_local = trc_qin_snow
            gwat_snow_local     = qin_snow
         ENDIF   ! snl < 0

         !--------------------------------------------------------
         ! 2) Strip wresi (soil layer excess over saturation) into pool
         !    Mirrors MOD_SoilSnowHydrology.F90 L1187-1192.
         !--------------------------------------------------------
         trc_wresi_sum = 0._r8
         wresi_sum     = 0._r8
         DO j = 1, nl_soil
            IF (t_soisno(j) > tfrz) THEN
               wresi_j = max(wliq_soisno_bef(j) - porsl(j)*dz_soisno(j)*1000._r8, 0._r8)
               IF (wresi_j > trc_tiny) THEN
                  wresi_sum = wresi_sum + wresi_j
                  IF (wliq_soisno_bef(j) > trc_tiny) THEN
                     ratio_src = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno_bef(j)
                     trc_wresi_j = min(wresi_j * ratio_src, &
                        max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                     trc_wliq_soisno(itrc, j, ipatch) = &
                        trc_wliq_soisno(itrc, j, ipatch) - trc_wresi_j
                     trc_wresi_sum = trc_wresi_sum + trc_wresi_j
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         !--------------------------------------------------------
         ! 3) Build mixed wetland pool (water + tracer) PRE loss
         !
         !    WATER_VSF wetwat formula variants (lb = snl+1):
         !      lb>=1, no split : gwat = pg_rain + sm - qseva
         !                        wetwat += (gwat - etr + qsdew + qfros - qsubl)*dt + wresi
         !      lb>=1, split    : gwat = pg_rain + sm - qseva_soil
         !                        wetwat += (gwat - etr + qsdew_s + qfros_s - qsubl_s)*dt + wresi
         !      lb<0 , no split : gwat = snowwater_output (rain/dew/frost/subl/evap all absorbed)
         !                        wetwat += (gwat - etr)*dt + wresi
         !      lb<0 , split    : gwat = snowwater_out + pg_rain*(1-fsno) - qseva_soil
         !                        wetwat += (gwat - etr + qsdew_s + qfros_s - qsubl_s)*dt + wresi
         !--------------------------------------------------------
         IF (snl < 0) THEN
            IF (split_soilsnow) THEN
               q_rain_in  = max(pg_rain * (1._r8 - fsno), 0._r8) * deltim
               q_dew_in   = max(qsdew_soil, 0._r8) * deltim
               q_frost_in = max(qfros_soil, 0._r8) * deltim
               q_evap_out = max(qseva_soil, 0._r8) * deltim
               q_subl_out = max(qsubl_soil, 0._r8) * deltim
            ELSE
               ! All atmospheric exchange absorbed by snowwater; no extra terms.
               q_rain_in  = 0._r8
               q_dew_in   = 0._r8
               q_frost_in = 0._r8
               q_evap_out = 0._r8
               q_subl_out = 0._r8
            ENDIF
            q_sm_in = 0._r8   ! snowmelt contributes via trc_gwat_snow_local
         ELSE
            IF (split_soilsnow) THEN
               q_rain_in  = max(pg_rain, 0._r8) * deltim
               q_dew_in   = max(qsdew_soil, 0._r8) * deltim
               q_frost_in = max(qfros_soil, 0._r8) * deltim
               q_evap_out = max(qseva_soil, 0._r8) * deltim
               q_subl_out = max(qsubl_soil, 0._r8) * deltim
            ELSE
               q_rain_in  = max(pg_rain, 0._r8) * deltim
               q_dew_in   = max(qsdew_in, 0._r8) * deltim
               q_frost_in = max(qfros_in, 0._r8) * deltim
               q_evap_out = max(qseva_in, 0._r8) * deltim
               q_subl_out = max(qsubl_in, 0._r8) * deltim
            ENDIF
            q_sm_in = max(sm, 0._r8) * deltim
         ENDIF
         q_etr_out = max(etr, 0._r8) * deltim

         ! Use SIGNED arithmetic for wdsrf_bef / wa_bef / wetwat_bef so the
         ! aquifer "debt" (wa<0 from a prior wetwat<0 deficit branch at
         ! MOD_SoilSnowHydrology.F90 L1200-1203) is mirrored exactly. Using
         ! max(.,0) here previously dropped wa_bef<0 from pool_water while
         ! keeping its positive inputs, which silently lost tracer when the
         ! debt was repaid (e.g. wa_bef=-10, input=+20 -> pool_water was 20
         ! instead of 10, so only half the 20R_atm input ended up in wetwat).
         ! wdsrf / wetwat are non-negative by WATER_VSF construction, but
         ! writing them signed keeps the formula uniform.
         pool_water = wdsrf_bef + wa_bef + wetwat_bef &
                    + wresi_sum + q_rain_in + q_sm_in + q_dew_in + q_frost_in

         pool_tracer = trc_wdsrf(itrc, ipatch) + trc_wa(itrc, ipatch) &
                     + trc_wetwat(itrc, ipatch) + trc_wresi_sum

         ! Drip / flood / paddy irrigation (if enabled) joins gwat in
         ! WATER_VSF ahead of the wetland merge, so it lands inside
         ! wetwat_post. Water side treats this as an internal transfer
         ! from the per-patch waterstorage reservoir (CoLMMAIN:806
         ! includes waterstorage in totwb/endwb), so the atmospheric
         ! input term ignores it. The tracer path mirrors that — debit
         ! trc_waterstorage at R_atm (Phase-1 invariant), add the
         ! matching mass to the wetland pool, and do NOT increment
         ! a_trc_precip. See MOD_Tracer_SoilWater.F90 main
         ! tracer_soil_water branch for the symmetric treatment.
         IF (qflx_irrig_ground > trc_tiny) THEN
            pool_water  = pool_water  + qflx_irrig_ground * deltim
            storage_ratio = R_atm
            IF (allocated(trc_waterstorage) .and. present(waterstorage_patch)) THEN
               IF (waterstorage_patch > trc_tiny) THEN
                  storage_ratio = trc_waterstorage(itrc, ipatch) / max(waterstorage_patch, trc_tiny)
               ENDIF
            ENDIF
            trc_flux    = qflx_irrig_ground * deltim * storage_ratio
            IF (allocated(trc_waterstorage)) THEN
               trc_flux = min(max(trc_flux, 0._r8), max(trc_waterstorage(itrc, ipatch), 0._r8))
            ENDIF
            pool_tracer = pool_tracer + trc_flux
            IF (allocated(trc_waterstorage)) THEN
               trc_waterstorage(itrc, ipatch) = trc_waterstorage(itrc, ipatch) - trc_flux
            ENDIF
         ENDIF

         IF (snl < 0) THEN
            ! Rain that reaches wetland pool directly is only (1-fsno) portion
            ! in split mode; in non-split mode all rain went through snow.
            IF (split_soilsnow) THEN
               pool_tracer = pool_tracer + trc_pg_rain_ground(itrc, ipatch) * (1._r8 - fsno)
            ENDIF
            pool_tracer = pool_tracer + trc_gwat_snow_local
            ! Mirror trc_gwat_snow_local with its water mass so
            ! pool_ratio = pool_tracer / pool_water stays consistent.
            ! Water side: wetwat += gwat*dt with gwat including the
            ! snow-bottom outflow (snowwater's qout at layer 0). Without
            ! this companion term, pool_water is systematically low and
            ! pool_ratio is overstated, over-draining trc_loss for snowy
            ! wetland patches (non-split mode is the worst case — no
            ! atmospheric term enters pool_water otherwise).
            pool_water  = pool_water  + gwat_snow_local
         ELSE
            pool_tracer = pool_tracer + trc_pg_rain_ground(itrc, ipatch) &
                                      + trc_sm_carry(itrc, ipatch)
         ENDIF

         ! Atmospheric input tracer (dew + frost)
         IF (q_dew_in + q_frost_in > trc_tiny) THEN
            trc_dew_input = q_dew_in * deposition_ratio_for(layer_temp(1), .false.)
            trc_frost_input = q_frost_in * deposition_ratio_for(layer_temp(1), .true.)
            pool_tracer = pool_tracer + trc_dew_input + trc_frost_input
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_dew_input + trc_frost_input
         ENDIF

         !--------------------------------------------------------
         ! 4) Remove evap/subl at atmospheric-exchange signatures and
         !    transpiration at mixed pool ratio.
         !    trc_loss is UNCAPPED — over-drawn losses push pool_tracer
         !    negative, which becomes tracer debt assigned to trc_wa at
         !    step 5 (matching the water-side debt wa<0). Capping here
         !    would break conservation across the deficit→recovery cycle.
         !--------------------------------------------------------
         loss_water = q_evap_out + q_subl_out + q_etr_out
         IF (abs(pool_water) > trc_tiny) THEN
            pool_ratio = pool_tracer / pool_water
         ELSE
            pool_ratio = R_atm
         ENDIF

         fractionate_pool_loss = pool_water > trc_tiny .and. pool_tracer > trc_tiny
         IF (fractionate_pool_loss) THEN
            trc_evap_loss = q_evap_out * evap_ratio_for(pool_ratio, layer_temp(1), .false.)
            trc_subl_loss = q_subl_out * evap_ratio_for(pool_ratio, layer_temp(1), .true.)
            trc_etr_loss  = q_etr_out  * pool_ratio
         ELSE
            trc_evap_loss = q_evap_out * pool_ratio
            trc_subl_loss = q_subl_out * pool_ratio
            trc_etr_loss  = q_etr_out  * pool_ratio
         ENDIF
         trc_loss = trc_evap_loss + trc_subl_loss + trc_etr_loss
	         IF (abs(trc_loss) > trc_tiny) THEN
	            ! Keep the total evaporation bucket identical, but split
	            ! process diagnostics: wetland q_etr is transpiration, not
	            ! open-water/soil evaporation.
	            IF (abs(trc_evap_loss + trc_subl_loss) > trc_tiny) THEN
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_evap_loss + trc_subl_loss, &
	                  q_evap_out + q_subl_out, TRC_EVAP_KIND_WETLAND)
	            ENDIF
	            IF (abs(trc_etr_loss) > trc_tiny) THEN
	               CALL tracer_book_evap_loss(itrc, ipatch, trc_etr_loss, q_etr_out, &
	                  TRC_EVAP_KIND_TRANSP)
	               a_trc_transp_src(itrc, ipatch) = a_trc_transp_src(itrc, ipatch) + trc_etr_loss
	            ENDIF
	            pool_tracer = pool_tracer - trc_loss
	         ENDIF
         pool_water = pool_water - loss_water

         !--------------------------------------------------------
         ! 5) Redistribute to post-WATER wdsrf / wa / wetwat / rsur.
         !    Water balance holds in signed form:
         !      pool_water_post_loss = wetwat + wdsrf + wa + rsur*dt
         !    with wa<0 in the deficit branch. wetwat and wdsrf are non-
         !    negative by WATER_VSF construction, so only wa carries the
         !    debt sign. pool_ratio is the mixed ratio of the post-loss
         !    pool and applies uniformly to all three post-WATER pools.
         !--------------------------------------------------------
         trc_rsur_local = 0._r8
         IF (pool_tracer < 0._r8) THEN
            ! Deficit branch (covers both wa<0 "true debt" and the rarer
            ! wa>=0 sign-mismatch state). The negative tracer belongs to the
            ! signed aquifer debt, not to positive water stores. Keep wetwat,
            ! wdsrf, and any exported rsur at a physically positive signature
            ! and absorb the compensating mass into trc_wa so patch-level
            ! conservation still closes. Setting positive stores to R=0
            ! produced apparent -999 permil wetland/surface-water deltas.
            pool_ratio = R_atm
            trc_wetwat(itrc, ipatch) = max(wetwat, 0._r8) * pool_ratio
            trc_wdsrf (itrc, ipatch) = max(wdsrf,  0._r8) * pool_ratio
            IF (rsur > trc_tiny) THEN
               trc_rsur_local = rsur * deltim * pool_ratio
               a_trc_rsur   (itrc, ipatch) = a_trc_rsur   (itrc, ipatch) + trc_rsur_local
               a_trc_rnof   (itrc, ipatch) = a_trc_rnof   (itrc, ipatch) + trc_rsur_local
               trc_rnof_step(itrc, ipatch) = trc_rnof_step(itrc, ipatch) + trc_rsur_local
            ENDIF
            trc_wa(itrc, ipatch) = pool_tracer - trc_rsur_local &
                                  - trc_wetwat(itrc, ipatch) - trc_wdsrf(itrc, ipatch)
         ELSE
            IF (abs(pool_water) > trc_tiny) THEN
               pool_ratio = pool_tracer / pool_water
            ELSE
               pool_ratio = 0._r8
            ENDIF

            trc_wetwat(itrc, ipatch) = wetwat * pool_ratio
            trc_wdsrf (itrc, ipatch) = wdsrf  * pool_ratio
            trc_wa    (itrc, ipatch) = wa     * pool_ratio

            IF (rsur > trc_tiny) THEN
               trc_rsur_local = rsur * deltim * pool_ratio
               a_trc_rsur   (itrc, ipatch) = a_trc_rsur   (itrc, ipatch) + trc_rsur_local
               a_trc_rnof   (itrc, ipatch) = a_trc_rnof   (itrc, ipatch) + trc_rsur_local
               trc_rnof_step(itrc, ipatch) = trc_rnof_step(itrc, ipatch) + trc_rsur_local
            ENDIF
         ENDIF

         ! The wetland water solver can leave a tiny signed mismatch between
         ! pool_water and the post-WATER destinations. Keep that internal
         ! redistribution residual in storage, not in the tracer budget. This
         ! applies to all tracer modes; dynamic forcing invalidates the old
         ! no-fractionation shortcut that reset wa to the current vapor ratio.
         redist_storage = trc_wetwat(itrc, ipatch) + trc_wdsrf(itrc, ipatch) &
                        + trc_wa(itrc, ipatch)
         redist_target = pool_tracer - trc_rsur_local
         redist_resid = redist_storage - redist_target
         IF (abs(redist_resid) > trc_tiny) THEN
            IF (redist_resid > 0._r8) THEN
               redist_fix = min(redist_resid, max(trc_wetwat(itrc, ipatch), 0._r8))
               trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) - redist_fix
               redist_resid = redist_resid - redist_fix

               redist_fix = min(redist_resid, max(trc_wdsrf(itrc, ipatch), 0._r8))
               trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) - redist_fix
               redist_resid = redist_resid - redist_fix

               trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - redist_resid
            ELSE
               IF (wetwat > trc_tiny) THEN
                  trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) - redist_resid
               ELSEIF (wdsrf > trc_tiny) THEN
                  trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) - redist_resid
               ELSE
                  trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - redist_resid
               ENDIF
            ENDIF
         ENDIF
      ENDDO

   CONTAINS

      real(r8) FUNCTION layer_temp (jlay)
         integer, intent(in) :: jlay

         layer_temp = 273.15_r8
         IF (jlay >= lbound(t_soisno, 1) .and. jlay <= ubound(t_soisno, 1)) THEN
            layer_temp = t_soisno(jlay)
         ENDIF
      END FUNCTION layer_temp

      real(r8) FUNCTION atmospheric_loss_tracer (pool_trc, pool_water, water_loss, temp_k, from_ice)
         real(r8), intent(in) :: pool_trc
         real(r8), intent(in) :: pool_water
         real(r8), intent(in) :: water_loss
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8), parameter :: max_evap_frac_step = 0.05_r8
         real(r8) :: source_ratio, flux_ratio
         real(r8) :: remaining_water, remaining_trc, loss_left, loss_step, trc_loss_step

         atmospheric_loss_tracer = 0._r8
         IF (pool_water <= trc_tiny .or. water_loss <= trc_tiny) RETURN
         IF (pool_trc <= trc_tiny) RETURN

         IF (water_loss >= pool_water * (1._r8 - 1.e-12_r8)) THEN
            atmospheric_loss_tracer = max(pool_trc, 0._r8)
            RETURN
         ENDIF

         source_ratio = max(pool_trc, 0._r8) / pool_water
         flux_ratio = evap_ratio_for(source_ratio, temp_k, from_ice)
         IF (water_loss <= max_evap_frac_step * pool_water .or. &
             abs(flux_ratio - source_ratio) <= 1.e-12_r8 * max(source_ratio, trc_tiny)) THEN
            atmospheric_loss_tracer = min(water_loss * flux_ratio, max(pool_trc, 0._r8))
            RETURN
         ENDIF

         ! Finite-pool evaporation can remove nearly all water in one land step.
         ! Integrating the loss prevents a single light vapor ratio from leaving
         ! an unrealistically enriched trace residual.
         remaining_water = pool_water
         remaining_trc = max(pool_trc, 0._r8)
         loss_left = water_loss
         DO WHILE (loss_left > trc_tiny .and. remaining_water > trc_tiny .and. &
                   remaining_trc > trc_tiny)
            loss_step = min(loss_left, max_evap_frac_step * remaining_water)
            source_ratio = remaining_trc / remaining_water
            flux_ratio = evap_ratio_for(source_ratio, temp_k, from_ice)
            trc_loss_step = min(loss_step * flux_ratio, remaining_trc)
            remaining_trc = remaining_trc - trc_loss_step
            remaining_water = remaining_water - loss_step
            loss_left = loss_left - loss_step
         ENDDO

         atmospheric_loss_tracer = min(max(pool_trc, 0._r8), &
            max(pool_trc, 0._r8) - max(remaining_trc, 0._r8))
      END FUNCTION atmospheric_loss_tracer

      real(r8) FUNCTION evap_ratio_for (source_ratio, temp_k, from_ice)
         real(r8), intent(in) :: source_ratio
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
         real(r8) :: relhum, alpha_k

         evap_ratio_for = source_ratio
         IF (.not. tracer_fractionation_active(itrc)) RETURN
         IF (.not. present(forc_q_frac) .or. .not. present(forc_psrf_frac)) RETURN

         relhum = tracer_surface_relhum(forc_q_frac, forc_psrf_frac, temp_k, from_ice)
         alpha_k = tracer_diffusivity_ratio_air(itrc)
         evap_ratio_for = tracer_craig_gordon_evap_ratio(itrc, source_ratio, R_atm, &
            temp_k, relhum, alpha_k, from_ice)
         evap_ratio_for = min(evap_ratio_for, max(source_ratio, 0._r8))
      END FUNCTION evap_ratio_for

      real(r8) FUNCTION deposition_ratio_for (temp_k, from_ice)
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice

         deposition_ratio_for = R_atm
         IF (.not. tracer_fractionation_active(itrc)) RETURN
         deposition_ratio_for = tracer_equilibrium_deposition_ratio(itrc, R_atm, temp_k, from_ice)
      END FUNCTION deposition_ratio_for

   END SUBROUTINE tracer_wetland

END MODULE MOD_Tracer_SoilWater
#endif
