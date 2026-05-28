#include <define.h>
#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_Driver
!=======================================================================
! Driver wrapper for Methane reactive-tracer physics.
!=======================================================================

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: methane_driver

CONTAINS

	SUBROUTINE methane_driver (istep,i,idate,patchclass,patchtype,deltim,lb,snl,dlon,dlat,&!input
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
			forc_t,forc_pbot,forc_po2m,forc_pco2m,forc_us,forc_vs,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai,rootr,fsatmax,fsatdcf,frcsat,f_h2osfc, &
		is_rice_paddy_in, rice_pft_frac_in)

		use MOD_Precision
		use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
		use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
		use MOD_Tracer_Reactive_Methane_Const
		use MOD_Namelist, only : DEF_USE_VariablySaturatedFlow
		use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad,z_soi,zi_soi,dz_soi
		use MOD_Tracer_Reactive_Methane_Physics
		use MOD_SPMD_Task
			USE MOD_Tracer_Reactive_Methane_BgcLink, only: tracer_ch4_bgc_patch_inputs, &
			     get_wetland_veg_proxy, get_rice_veg_proxy, is_paddy_rice_live, &
			     organic_max, get_biome_f_methane, get_biome_redoxlag
			USE MOD_Tracer_Reactive_Methane_VegOverride, only: wetland_aere_active
		USE MOD_Tracer_Reactive_Methane_Microbes, only: methane_microbes_step, &
		     microbial_prod_potential, microbial_oxid_potential
		! ----- Methane tracer state (module-level allocatables) -----
		USE MOD_Tracer_Reactive_Methane_State, only: net_methane, methane_prod_depth, o2_decomp_depth, &
		     co2_decomp_depth, methane_oxid_depth, o2_oxid_depth, co2_oxid_depth, &
		     methane_aere_depth, methane_tran_depth, o2_aere_depth, co2_aere_depth, methane_ebul_depth, &
		     o2stress, methane_stress, &
		     methane_surf_flux_tot, methane_surf_flux_tot_phys, methane_surf_aere, methane_surf_ebul, methane_surf_diff, &
		     methane_surf_diff_phys, &
		     methane_balance_residual, methane_ch4_clip_credit, o2_cap_loss, o2_cap_gain, &
		     methane_ebul_tot, methane_prod_tot, methane_oxid_tot, &
		     co2_decomp_tot, co2_oxid_tot, co2_aere_tot, co2_net_tot, &
		     totcol_methane, grnd_methane_cond, conc_o2, conc_methane, &
		     net_methane_unsat, net_methane_sat, &
		     methane_prod_depth_unsat, methane_prod_depth_sat, &
		     o2_decomp_depth_unsat, o2_decomp_depth_sat, &
		     co2_decomp_depth_unsat, co2_decomp_depth_sat, &
		     methane_oxid_depth_unsat, methane_oxid_depth_sat, &
		     o2_oxid_depth_unsat, o2_oxid_depth_sat, &
		     co2_oxid_depth_unsat, co2_oxid_depth_sat, &
		     methane_aere_depth_unsat, methane_aere_depth_sat, &
		     methane_tran_depth_unsat, methane_tran_depth_sat, &
		     o2_aere_depth_unsat, o2_aere_depth_sat, &
		     co2_aere_depth_unsat, co2_aere_depth_sat, &
		     methane_ebul_depth_unsat, methane_ebul_depth_sat, &
		     o2stress_unsat, o2stress_sat, methane_stress_unsat, methane_stress_sat, &
		     methane_surf_flux_tot_unsat, methane_surf_flux_tot_sat, &
		     methane_surf_aere_unsat, methane_surf_aere_sat, &
		     methane_surf_ebul_unsat, methane_surf_ebul_sat, &
		     methane_surf_diff_unsat, methane_surf_diff_sat, &
		     methane_surf_diff_phys_unsat, methane_surf_diff_phys_sat, &
		     methane_ebul_tot_unsat, methane_ebul_tot_sat, &
		     methane_prod_tot_unsat, methane_prod_tot_sat, &
		     methane_oxid_tot_unsat, methane_oxid_tot_sat, &
		     co2_decomp_tot_unsat, co2_decomp_tot_sat, &
		     co2_oxid_tot_unsat, co2_oxid_tot_sat, &
		     co2_net_tot_unsat, co2_net_tot_sat, &
		     totcol_methane_unsat, totcol_methane_sat, &
		     grnd_methane_cond_unsat, grnd_methane_cond_sat, &
		     conc_o2_unsat, conc_o2_sat, conc_methane_unsat, conc_methane_sat, &
		     methane_prod_depth_lake, methane_oxid_depth_lake, methane_ebul_depth_lake, &
		     co2_decomp_depth_lake, co2_oxid_depth_lake, &
		     methane_surf_ebul_lake, methane_surf_diff_lake, methane_surf_flux_tot_lake, &
		     methane_prod_tot_lake, methane_oxid_tot_lake, methane_ebul_tot_lake, &
		     co2_decomp_tot_lake, co2_oxid_tot_lake, co2_net_tot_lake, &
		     totcol_methane_lake, grnd_methane_cond_lake, conc_o2_lake, conc_methane_lake, &
		     c_atm, forc_pmethanem, layer_sat_lag, lake_soilc, &
		     annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
		     tempavg_agnpp, tempavg_bgnpp, annsum_counter, &
		     tempavg_somhr, tempavg_finrw, &
		     fsat_bef, finundated_lag, methane_dfsat_tot, &
		     biome_f_methane_patch, biome_redoxlag_patch, &
		     f_inund_flood_patch

		IMPLICIT NONE
		integer ,intent(in) :: istep

		integer ,intent(in) :: i         ! patch index
		integer ,intent(in) :: idate(1:3)  ! current date (year, day of the year, seconds of the day)
		integer, intent(in) :: &
							patchclass  ,&! land patch class of USGS classification or others
							patchtype     ! land patch type (0=soil, 1=urban and built-up,
								          ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
		real(r8),intent(in) :: deltim    ! time step in seconds
		integer ,intent(in) :: &
				lb,&
				snl
		real(r8),intent(in) :: dlon     ! longitude (degrees)
		real(r8),intent(in) :: dlat     ! latitude (degrees)

		real(r8),intent(in) :: &
				z_soisno   (maxsnl+1:nl_soil) , &! layer depth (m)
				dz_soisno  (maxsnl+1:nl_soil) , &! layer thickness (m)
				zi_soisno  (maxsnl:nl_soil)   , &! interface level below a "z" level (m)
				t_soisno   (maxsnl+1:nl_soil) , &! soil + snow layer temperature [K]
				t_grnd                 		  , &! ground surface temperature [k]
				wliq_soisno(maxsnl+1:nl_soil) , &! liquid water (kg/m2)
				wice_soisno(maxsnl+1:nl_soil) , &! ice lens (kg/m2)
				forc_t                        , &! temperature at agcm reference height [kelvin]
					forc_pbot                     , &! atmosphere pressure at the bottom of the atmos. model level [pa]
					forc_po2m                     , &! partial pressure of O2 at observational height [pa]
					forc_pco2m                    , &! partial pressure of CO2 at observational height [Pa]
					forc_us                       , &! eastward wind speed [m/s], used as U10 proxy for lake gas exchange
					forc_vs                       , &! northward wind speed [m/s], used as U10 proxy for lake gas exchange
					zwt                           , &! the depth to water table [m]
				rootfr     (1:nl_soil)        , &! fraction of roots in each soil layer
				snowdp                        , &! snow depth (m)
				wat                           , &! total water storage [mm] (reserved interface)
				rsur                          , &! surface runoff [mm H2O/s] (reserved interface)
				etr                           , &! transpiration rate [mm/s]
				lakedepth                     , &! lake depth (m), used by CTSM-aligned lake CH4 ebullition pressure
				lake_icefrac(1:nl_lake)       , &! lake frozen mass fraction for lake CH4 ebullition suppression
				wdsrf                         , &! depth of surface water [mm]
		wetwat                        , &! water storage in wetland [mm]
				bsw         (1:nl_soil)       , &! clapp and hornbereger "b" parameter [-]
				smp         (1:nl_soil)       , &! soil matrix potential [mm]
				porsl       (1:nl_soil)       , &! fraction of soil that is voids [-]
				lai                           , &! leaf area index
				rootr       (1:nl_soil)       , &! water exchange between soil and root. Positive: soil->root [?]

				fsatmax                       , &! maximum saturated area fraction [-]
				fsatdcf                       , &! decay factor in calculation of saturated area fraction [1/m]
		frcsat                        , &! fraction of saturation area
				f_h2osfc                        ! fraction of surface water [-], maintained by WATER_2014 / WATER_VSF

		! Optional R1 rice paddy inputs from CoLMDRIVER.  When absent the
		! routine behaves exactly as before (default-off).
		logical, intent(in), optional :: is_rice_paddy_in
		real(r8), intent(in), optional :: rice_pft_frac_in

		logical  :: is_rice_paddy
		logical  :: is_floodplain_active
		real(r8) :: rice_pft_frac

		real(r8):: &
				crootfr  (1:nl_soil)     , &! fraction of roots for carbon in each soil layer
				pH                       , &! soil water pH
				cellorg  (1:nl_soil)     , &! column 3D org (kg/m3 organic matter)
				t_h2osfc             	    ! surface water temperature

			real(r8) :: somhr_loc, lithr_loc, rr_loc, agnpp_loc, bgnpp_loc, annsum_npp_loc
			real(r8) :: hr_vr_loc(1:nl_soil), fphr_loc(1:nl_soil)
			real(r8) :: o_scalar_loc(1:nl_soil), pot_f_nit_vr_loc(1:nl_soil)
				real(r8) :: microbe_conc_o2(1:nl_soil), microbe_conc_ch4(1:nl_soil)
				real(r8) :: microbial_prod_potential_eff(1:nl_soil)
				real(r8) :: microbial_oxid_potential_eff(1:nl_soil)
			real(r8) :: lai_eff
			real(r8) :: rootfr_eff(1:nl_soil)
			real(r8) :: rootr_eff(1:nl_soil)
				real(r8) :: forc_t_eff, forc_pbot_eff, forc_po2m_eff, forc_pco2m_eff
				real(r8) :: forc_us_eff, forc_vs_eff
			real(r8) :: fprev
			logical  :: bgc_inputs_ready
			logical, save :: warned_missing_bgc_inputs = .false.

		! Use raw t_grnd as proxy for surface-water temperature.
		! Do NOT clamp to tfrz: the ponddiff formula handles sub-freezing
		! T (parabolic in t_c) without issue, and the downstream ice-block
		! path (Physics.F90 line ~2274) needs t_h2osfc < tfrz to fire.
		! Clamping was the cause of the dead ice-block branch.
		t_h2osfc = t_grnd

		! Sanitize atmospheric forcing used by methane gas exchange.  The
		! subroutine arguments are intent(in), so keep local copies rather
		! than mutating caller state.  This prevents coast/domain-edge spval
		! or NaN forcing from poisoning c_atm and lake ebullition pressure.

		! R1 rice paddy: unpack optional inputs, default off.
		IF (present(is_rice_paddy_in)) THEN
		   is_rice_paddy = is_rice_paddy_in
		ELSE
		   is_rice_paddy = .false.
		ENDIF
		IF (present(rice_pft_frac_in)) THEN
		   rice_pft_frac = rice_pft_frac_in
		ELSE
		   rice_pft_frac = 0._r8
		ENDIF

		forc_t_eff = forc_t
		if (ieee_is_nan(forc_t_eff) .or. abs(forc_t_eff) >= 0.5_r8*abs(spval) .or. &
		    forc_t_eff < 150._r8 .or. forc_t_eff > 350._r8) forc_t_eff = 288.15_r8

		forc_pbot_eff = forc_pbot
		if (ieee_is_nan(forc_pbot_eff) .or. abs(forc_pbot_eff) >= 0.5_r8*abs(spval) .or. &
		    forc_pbot_eff <= 0._r8) forc_pbot_eff = 101325._r8

		forc_po2m_eff = forc_po2m
		if (ieee_is_nan(forc_po2m_eff) .or. abs(forc_po2m_eff) >= 0.5_r8*abs(spval) .or. &
		    forc_po2m_eff <= 0._r8) forc_po2m_eff = 0.2095_r8 * forc_pbot_eff

			forc_pco2m_eff = forc_pco2m
			if (ieee_is_nan(forc_pco2m_eff) .or. abs(forc_pco2m_eff) >= 0.5_r8*abs(spval) .or. &
			    forc_pco2m_eff < 0._r8) forc_pco2m_eff = 415.e-6_r8 * forc_pbot_eff

			forc_us_eff = forc_us
			if (ieee_is_nan(forc_us_eff) .or. abs(forc_us_eff) >= 0.5_r8*abs(spval)) forc_us_eff = 0._r8

			forc_vs_eff = forc_vs
			if (ieee_is_nan(forc_vs_eff) .or. abs(forc_vs_eff) >= 0.5_r8*abs(spval)) forc_vs_eff = 0._r8

			CALL tracer_ch4_bgc_patch_inputs (i, rootfr, crootfr, pH, cellorg, &
			     somhr_loc, lithr_loc, hr_vr_loc, rr_loc, agnpp_loc, bgnpp_loc, &
			     annsum_npp_loc, fphr_loc, o_scalar_loc, pot_f_nit_vr_loc, bgc_inputs_ready)
			IF (patchtype == 0 .and. .not. bgc_inputs_ready) THEN
				IF (p_is_master .and. .not. warned_missing_bgc_inputs) THEN
					write(6,*) 'WARNING: Soil methane running with sanitized BGC defaults before BGC/PFT inputs are initialized.'
					warned_missing_bgc_inputs = .true.
				ENDIF
			ENDIF

			! ---- Wetland vegetation proxy (patchtype==2 only) -------------------
			! IGBP class 11 wetland patches carry no PFT, so NPP / rootfr arrive
			! as 0 and methane_aere (Wania 2010) cannot fire fully, dropping
			! the plant-mediated CH4 pathway that observations show carries
			! 50-90% of total wetland CH4 efflux (Bridgham 2013 GCB).
			!
			! LAI handling: CoLM mksrfdata aggregates Yuan+2011 LAI onto
			! patchtype==2 patches (fveg0_igbp(11)=1), so `lai` has real
			! data + DEF_LAI_MONTHLY seasonality.  get_wetland_veg_proxy
			! trusts it when valid (0.1 < lai < 20) and only falls back to
			! a climate-zone peak when missing.
			!
			! For patchtype /= 2: use the inputs as-is (vegetated patches
			! already have their own PFT-driven LAI/NPP/roots, which
			! correctly drives the aerenchyma path even under seasonal
			! inundation, e.g. floodplain pasture or seasonally flooded forest).
			lai_eff       = lai
			rootfr_eff(:) = rootfr(:)
			rootr_eff(:)  = rootr(:)
			! Clear any prior aerenchyma override (so non-wetland patches use
			! the DEF_METHANE%* defaults via the getter fallback path).
			IF (allocated(wetland_aere_active)) THEN
				IF (i >= 1 .and. i <= size(wetland_aere_active)) &
					wetland_aere_active(i) = .false.
			ENDIF
			IF (patchtype == 2) THEN
				CALL get_wetland_veg_proxy (dlat, cellorg(1), lai, i, &
				     lai_eff, annsum_npp_loc, agnpp_loc, bgnpp_loc, rootfr_eff)
				! BgcLink set crootfr from the original (zero) rootfr; replace it
				! with the proxy profile so methane_prod distributes root
				! respiration into the right layers.
				crootfr(1:nl_soil) = rootfr_eff(1:nl_soil)
				! Wetland patches have no PFT-level root water uptake or root
				! respiration diagnostics.  Use the same shallow proxy profile
				! for transpiration loss and a conservative belowground-respiration
				! proxy so root O2 demand is not silently zero.
				rootr_eff(1:nl_soil) = rootfr_eff(1:nl_soil)
				rr_loc = max(rr_loc, 0.5_r8 * bgnpp_loc)
			ENDIF
			! R4 rice paddy aerenchyma: rice patches (patchtype==0) get their own
			! Wania-style override (Zone 6) so methane_aere uses rice tiller
			! geometry instead of the natural-wetland default.  NPP / rootfr stay
			! whatever CN provides for the actual rice CFT — only the aerenchyma
			! geometry channel is overridden.
			!
			! Gate on croplive_p so winter stubble/fallow on rice patches drops
			! back to the wetland default (otherwise the LAI-scaled aere keeps
			! a residual 50% transport even after harvest).  rice_paddy_alive
			! flag is set by the same any_paddy_rice_live() helper used in
			! methane() for finundated override, keeping the two physics
			! decisions consistent.
			! All BGC state access goes through MOD_Tracer_Reactive_Methane_BgcLink
			! (architectural rule enforced by scripts/ci/check_methane_integration.sh).
			IF (is_rice_paddy .and. rice_pft_frac > 0._r8) THEN
			   IF (is_paddy_rice_live(i)) THEN
			      CALL get_rice_veg_proxy (lai, i)
			   ENDIF
			ENDIF

			! Biome-specific f_methane lookup (Bridgham 2013 GCB).  Sets
			! biome_f_methane_patch(i) which methane_prod consumes when
			! DEF_METHANE%use_biome_f_methane is true; otherwise returns
			! the legacy DEF_METHANE%f_methane scalar (backwards compatible).
			IF (allocated(biome_f_methane_patch) .and. &
			    i >= 1 .and. i <= size(biome_f_methane_patch)) THEN
				is_floodplain_active = patchtype == 0 .and. .not. is_rice_paddy .and. &
				                       DEF_METHANE%use_routing_for_soil .and. &
				                       allocated(f_inund_flood_patch) .and. &
				                       i >= 1 .and. i <= size(f_inund_flood_patch) .and. &
				                       f_inund_flood_patch(i) > DEF_METHANE%hybrid_soil_threshold
				biome_f_methane_patch(i) = get_biome_f_methane (patchtype, dlat, cellorg(1), &
				                                                is_rice_paddy, is_floodplain_active)
			ENDIF

			! Biome-specific redoxlag lookup (Pangala 2017, Whalen 1990).
			! Tropical warm wetlands respond faster (~5-15d); cold boreal
			! peat slower (~30-60d).  Consumed by methane when
			! DEF_METHANE%use_biome_redoxlag is true.
			IF (allocated(biome_redoxlag_patch) .and. &
			    i >= 1 .and. i <= size(biome_redoxlag_patch)) THEN
				biome_redoxlag_patch(i) = get_biome_redoxlag (patchtype, dlat, cellorg(1), &
				                                              is_rice_paddy)
			ENDIF
					! f_h2osfc is maintained by the water module before methane_driver.
					! Lake CH4 uses the CTSM-style sediment-carbon pathway, not the optional
					! soil microbial-pool override; keep lake microbial pools from evolving as
					! diagnostics-only state.
					IF (patchtype /= 4) THEN
						IF (abs(fsat_bef(i)) < 1.e30_r8 .and. fsat_bef(i) >= 0._r8 .and. &
						    fsat_bef(i) <= 1._r8) THEN
							fprev = fsat_bef(i)
							microbe_conc_o2(:) = conc_o2_sat(1:nl_soil,i) * fprev + &
								conc_o2_unsat(1:nl_soil,i) * (1._r8 - fprev)
							microbe_conc_ch4(:) = conc_methane_sat(1:nl_soil,i) * fprev + &
								conc_methane_unsat(1:nl_soil,i) * (1._r8 - fprev)
						ELSE
							microbe_conc_o2(:) = conc_o2(1:nl_soil,i)
							microbe_conc_ch4(:) = conc_methane(1:nl_soil,i)
						ENDIF
							CALL methane_microbes_step (i, deltim, t_soisno(1:nl_soil), &
							     microbe_conc_o2, microbe_conc_ch4, hr_vr_loc, cellorg)
						ELSE
							! Lake CH4 bypasses the optional soil microbial-pool override.
							! Clear per-layer potentials so accumulated diagnostics do not
							! retain spval/restart values on lake patches.
							IF (allocated(microbial_prod_potential) .and. allocated(microbial_oxid_potential)) THEN
								IF (i >= lbound(microbial_prod_potential,2) .and. i <= ubound(microbial_prod_potential,2) .and. &
								    i >= lbound(microbial_oxid_potential,2) .and. i <= ubound(microbial_oxid_potential,2)) THEN
									microbial_prod_potential(1:nl_soil,i) = 0._r8
									microbial_oxid_potential(1:nl_soil,i) = 0._r8
								ENDIF
							ENDIF
						ENDIF

						microbial_prod_potential_eff(:) = 0._r8
						microbial_oxid_potential_eff(:) = 0._r8
						IF (allocated(microbial_prod_potential) .and. allocated(microbial_oxid_potential)) THEN
							IF (i >= lbound(microbial_prod_potential,2) .and. i <= ubound(microbial_prod_potential,2) .and. &
							    i >= lbound(microbial_oxid_potential,2) .and. i <= ubound(microbial_oxid_potential,2)) THEN
								microbial_prod_potential_eff(:) = microbial_prod_potential(1:nl_soil,i)
								microbial_oxid_potential_eff(:) = microbial_oxid_potential(1:nl_soil,i)
							ENDIF
						ENDIF

			CALL methane (istep,i,idate(1:3),patchclass,patchtype,lb,snl,dlon,dlat,deltim,&
		z_soisno(maxsnl+1:),dz_soisno(maxsnl+1:),zi_soisno(maxsnl:),t_soisno(maxsnl+1:),&
		t_grnd,wliq_soisno(maxsnl+1:),wice_soisno(maxsnl+1:),&
			forc_t_eff,forc_pbot_eff,forc_po2m_eff,forc_pco2m_eff,forc_us_eff,forc_vs_eff,&
		zwt,rootfr_eff,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai_eff,rootr_eff,&
		annsum_npp_loc, rr_loc,&
		fsatmax,fsatdcf,frcsat,&
		agnpp_loc, bgnpp_loc, somhr_loc,&
		crootfr(1:nl_soil), lithr_loc, hr_vr_loc(1:nl_soil), o_scalar_loc(1:nl_soil), &
		fphr_loc(1:nl_soil), pot_f_nit_vr_loc(1:nl_soil), pH,&
			cellorg(1:nl_soil),t_h2osfc,organic_max,&
				microbial_prod_potential_eff, microbial_oxid_potential_eff, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane(i), &
		methane_prod_depth(1:nl_soil,i), o2_decomp_depth(1:nl_soil,i), co2_decomp_depth(1:nl_soil,i), &
		methane_oxid_depth(1:nl_soil,i), o2_oxid_depth(1:nl_soil,i), co2_oxid_depth(1:nl_soil,i), &
		methane_aere_depth(1:nl_soil,i), methane_tran_depth(1:nl_soil,i), &
		o2_aere_depth(1:nl_soil,i), co2_aere_depth(1:nl_soil,i), &
		methane_ebul_depth(1:nl_soil,i), &
		o2stress(1:nl_soil,i), methane_stress(1:nl_soil,i), &
		methane_surf_flux_tot(i), methane_surf_flux_tot_phys(i), methane_surf_aere(i), methane_surf_ebul(i), methane_surf_diff(i), &
		methane_surf_diff_phys(i), &
		methane_balance_residual(i), methane_ch4_clip_credit(i), o2_cap_loss(i), o2_cap_gain(i), &
		methane_ebul_tot(i), methane_prod_tot(i), methane_oxid_tot(i), &
		co2_decomp_tot(i), co2_oxid_tot(i), co2_aere_tot(i), co2_net_tot(i), &
		totcol_methane(i), grnd_methane_cond(i), conc_o2(1:nl_soil,i), conc_methane(1:nl_soil,i), &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane_unsat(i), net_methane_sat(i), &
		methane_prod_depth_unsat(1:nl_soil,i), methane_prod_depth_sat(1:nl_soil,i), &
		o2_decomp_depth_unsat(1:nl_soil,i), o2_decomp_depth_sat(1:nl_soil,i), &
		co2_decomp_depth_unsat(1:nl_soil,i), co2_decomp_depth_sat(1:nl_soil,i), &
		methane_oxid_depth_unsat(1:nl_soil,i), methane_oxid_depth_sat(1:nl_soil,i), &
		o2_oxid_depth_unsat(1:nl_soil,i), o2_oxid_depth_sat(1:nl_soil,i), &
		co2_oxid_depth_unsat(1:nl_soil,i), co2_oxid_depth_sat(1:nl_soil,i), &
		methane_aere_depth_unsat(1:nl_soil,i), methane_aere_depth_sat(1:nl_soil,i), &
		methane_tran_depth_unsat(1:nl_soil,i), methane_tran_depth_sat(1:nl_soil,i), &
		o2_aere_depth_unsat(1:nl_soil,i), o2_aere_depth_sat(1:nl_soil,i), &
		co2_aere_depth_unsat(1:nl_soil,i), co2_aere_depth_sat(1:nl_soil,i), &
		methane_ebul_depth_unsat(1:nl_soil,i), methane_ebul_depth_sat(1:nl_soil,i), &
		o2stress_unsat(1:nl_soil,i), o2stress_sat(1:nl_soil,i), &
		methane_stress_unsat(1:nl_soil,i), methane_stress_sat(1:nl_soil,i), &
		methane_surf_flux_tot_unsat(i), methane_surf_flux_tot_sat(i), &
		methane_surf_aere_unsat(i), methane_surf_aere_sat(i), &
		methane_surf_ebul_unsat(i), methane_surf_ebul_sat(i), &
		methane_surf_diff_unsat(i), methane_surf_diff_sat(i), &
		methane_surf_diff_phys_unsat(i), methane_surf_diff_phys_sat(i), &
		methane_ebul_tot_unsat(i), methane_ebul_tot_sat(i), &
		methane_prod_tot_unsat(i), methane_prod_tot_sat(i), &
		methane_oxid_tot_unsat(i), methane_oxid_tot_sat(i), &
		co2_decomp_tot_unsat(i), co2_decomp_tot_sat(i), &
		co2_oxid_tot_unsat(i), co2_oxid_tot_sat(i), &
		co2_net_tot_unsat(i), co2_net_tot_sat(i), &
		totcol_methane_unsat(i), totcol_methane_sat(i), &
		grnd_methane_cond_unsat(i), grnd_methane_cond_sat(i), &
		conc_o2_unsat(1:nl_soil,i), conc_o2_sat(1:nl_soil,i), &
		conc_methane_unsat(1:nl_soil,i), conc_methane_sat(1:nl_soil,i), &
		methane_prod_depth_lake(1:nl_soil,i), methane_oxid_depth_lake(1:nl_soil,i), &
		methane_ebul_depth_lake(1:nl_soil,i), &
		co2_decomp_depth_lake(1:nl_soil,i), co2_oxid_depth_lake(1:nl_soil,i), &
		methane_surf_ebul_lake(i), methane_surf_diff_lake(i), methane_surf_flux_tot_lake(i), &
		methane_prod_tot_lake(i), methane_oxid_tot_lake(i), methane_ebul_tot_lake(i), &
		co2_decomp_tot_lake(i), co2_oxid_tot_lake(i), co2_net_tot_lake(i), &
		totcol_methane_lake(i), grnd_methane_cond_lake(i), conc_o2_lake(1:nl_soil,i), conc_methane_lake(1:nl_soil,i), &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm(1:3,i), forc_pmethanem(i), layer_sat_lag(1:nl_soil,i), lake_soilc(1:nl_soil,i), &
		annavg_agnpp(i), annavg_bgnpp(i), annavg_somhr(i), annavg_finrw(i), &
		tempavg_agnpp(i), tempavg_bgnpp(i), annsum_counter(i), tempavg_somhr(i), &
		tempavg_finrw(i), fsat_bef(i), finundated_lag(i), methane_dfsat_tot(i), f_h2osfc, &
		is_rice_paddy_in=is_rice_paddy, rice_pft_frac_in=rice_pft_frac)

	END SUBROUTINE methane_driver

END MODULE MOD_Tracer_Reactive_Methane_Driver
#endif
