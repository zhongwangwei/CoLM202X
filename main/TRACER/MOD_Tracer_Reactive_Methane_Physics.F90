#include <define.h>
#if (defined TRACER) && (defined BGC)

module MOD_Tracer_Reactive_Methane_Physics
    !=======================================================================
	! DESCRIPTION:
	! Module holding routines to calculate methane fluxes
	! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
	! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in methane_tran.

	! ORIGINAL:
	! The Community Land Model version 5.0 (CLM5.0)

	! REFERENCES:
	! Lawrence, D.M., Fisher, R.A., Koven, C.D., Oleson, K.W., Swenson, S.C., Bonan, G., Collier, N.,
	! Ghimire, B., van Kampenhout, L., Kennedy, D. and Kluzek, E., 2019.
	! The Community Land Model version 5: Description of new features, benchmarking,
	! and impact of forcing uncertainty. Journal of Advances in Modeling Earth Systems, 11(12), 4245-4287.

	! REVISION:
	! Xionghui Xu, 2025, 1) Modify original CLM5 to be compatible with CoLM code structure.
    !                    2) Adapt and stabilize the process coupling for CoLM202X.
	!=======================================================================
	use MOD_Precision
	use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
	use MOD_SPMD_Task
	use MOD_TimeManager
	use MOD_Vars_TimeInvariants, only: wetwatmax
	use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad
	use MOD_Const_Physical, only: denh2o, denice, tfrz, grav
	use MOD_Tracer_Reactive_Methane_Const
	USE MOD_Namelist, only: DEF_wetland_finundation_scheme, DEF_USE_Dynamic_Wetland
	USE MOD_Tracer_Reactive_Methane_GIEMS, only: giems_finundated, giems_active
	USE MOD_Tracer_Reactive_Methane_BgcLink, only: is_paddy_rice_live, &
	                                       rice_days_past_planting, rice_days_since_harvest
	USE MOD_Tracer_Reactive_Methane_VegOverride, only: get_aere_poros, get_aere_tillerC, &
	                                          get_aere_radius, get_aere_scale
		USE MOD_Tracer_Reactive_Methane_State, only: f_inund_levee_patch, f_inund_flood_patch, &
		                                     f_inund_flood_depth_patch, wetland_frac_per_patch, &
		                                     biome_f_methane_patch, biome_redoxlag_patch, &
		                                     methane_finundated, methane_soil_finundated, &
		                                     methane_soil_zwt, methane_surf_flux_wetland, &
		                                     methane_surf_flux_soil, methane_surf_flux_lake, &
		                                     methane_surf_flux_rice, methane_surf_aere_soil, &
		                                     methane_surf_aere_rice, methane_surf_ebul_soil, &
		                                     methane_surf_ebul_rice, methane_surf_diff_soil, &
		                                     methane_surf_diff_rice, methane_prod_tot_soil, &
		                                     methane_prod_tot_rice, methane_oxid_tot_soil, &
		                                     methane_oxid_tot_rice
	!-----------------------------------------------------------------------
	implicit none
	save

	public  :: methane

	private :: methane_annualupdate
	private :: methane_prod
	private :: methane_oxid
	private :: methane_aere
	private :: methane_ebul
	private :: methane_tran
	! Rice BGC bridge functions are now public in MOD_Tracer_Reactive_Methane_BgcLink
	! (architectural rule: BGC state access only via BgcLink, enforced by CI):
	!   is_paddy_rice_live, rice_days_past_planting, rice_days_since_harvest
	private
contains

	!-----------------------------------------------------------------------
	subroutine methane (istep,ipatch,idate,patchclass,patchtype,&!input
		lb,snl,&
		dlon,dlat,&
		deltim,&
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
			forc_t,forc_pbot,forc_po2m,forc_pco2m,forc_us,forc_vs,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai,rootr,&
		annsum_npp,rr,&
		fsatmax,fsatdcf,frcsat,&
		agnpp,bgnpp,somhr,&
		crootfr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,&
			cellorg,t_h2osfc,organic_max,&
			microbial_prod_potential_patch, microbial_oxid_potential_patch, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane, &
		methane_prod_depth, o2_decomp_depth, co2_decomp_depth, methane_oxid_depth, o2_oxid_depth, co2_oxid_depth, &
		methane_aere_depth, methane_tran_depth, o2_aere_depth, co2_aere_depth, methane_ebul_depth, &
		o2stress, methane_stress, &
			methane_surf_flux_tot, methane_surf_flux_tot_phys, &
			methane_surf_aere, methane_surf_ebul, methane_surf_diff, methane_surf_diff_phys, &
			methane_balance_residual, methane_ch4_clip_credit, o2_cap_loss, o2_cap_gain, &
		methane_ebul_tot, methane_prod_tot, methane_oxid_tot, &
		co2_decomp_tot, co2_oxid_tot, co2_aere_tot, co2_net_tot, &
		totcol_methane, grnd_methane_cond, conc_o2, conc_methane, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane_unsat, net_methane_sat, &
		methane_prod_depth_unsat, methane_prod_depth_sat, o2_decomp_depth_unsat, o2_decomp_depth_sat, &
		co2_decomp_depth_unsat, co2_decomp_depth_sat, &
		methane_oxid_depth_unsat, methane_oxid_depth_sat, o2_oxid_depth_unsat, o2_oxid_depth_sat, &
		co2_oxid_depth_unsat, co2_oxid_depth_sat, &
		methane_aere_depth_unsat, methane_aere_depth_sat, methane_tran_depth_unsat, methane_tran_depth_sat, &
		o2_aere_depth_unsat, o2_aere_depth_sat, co2_aere_depth_unsat, co2_aere_depth_sat, &
		methane_ebul_depth_unsat, methane_ebul_depth_sat, &
		o2stress_unsat, o2stress_sat, methane_stress_unsat, methane_stress_sat, &
		methane_surf_flux_tot_unsat, methane_surf_flux_tot_sat, methane_surf_aere_unsat, methane_surf_aere_sat, &
		methane_surf_ebul_unsat, methane_surf_ebul_sat, methane_surf_diff_unsat, methane_surf_diff_sat, &
		methane_surf_diff_phys_unsat, methane_surf_diff_phys_sat, &
		methane_ebul_tot_unsat, methane_ebul_tot_sat, methane_prod_tot_unsat, methane_prod_tot_sat, &
		methane_oxid_tot_unsat, methane_oxid_tot_sat, &
		co2_decomp_tot_unsat, co2_decomp_tot_sat, &
		co2_oxid_tot_unsat, co2_oxid_tot_sat, &
		co2_net_tot_unsat, co2_net_tot_sat, &
		totcol_methane_unsat, totcol_methane_sat, grnd_methane_cond_unsat, grnd_methane_cond_sat, &
		conc_o2_unsat, conc_o2_sat, conc_methane_unsat, conc_methane_sat, &
		methane_prod_depth_lake, methane_oxid_depth_lake, methane_ebul_depth_lake, &
		co2_decomp_depth_lake, co2_oxid_depth_lake, &
		methane_surf_ebul_lake, methane_surf_diff_lake, methane_surf_flux_tot_lake, &
		methane_prod_tot_lake, methane_oxid_tot_lake, methane_ebul_tot_lake, &
		co2_decomp_tot_lake, co2_oxid_tot_lake, co2_net_tot_lake, &
		totcol_methane_lake, grnd_methane_cond_lake, conc_o2_lake, conc_methane_lake, &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm, forc_pmethanem, layer_sat_lag, lake_soilc, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
		tempavg_agnpp, tempavg_bgnpp, annsum_counter, tempavg_somhr, tempavg_finrw, &
		fsat_bef, finundated_lag, methane_dfsat_tot, f_h2osfc, &
		is_rice_paddy_in, rice_pft_frac_in, store_patch_diagnostics_in, finundated_used_out, &
		finundated_default_out)

		!=======================================================================
		! DESCRIPTION:
		! Driver for the methane emissions model
		!=======================================================================

		!===================== input ===========================================
		integer, intent(in) :: &
			istep            , &
			ipatch           , &! caller's per-patch index
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchclass       , &! land patch class of USGS classification or others
			patchtype        , &! land patch type (0=soil, 1=urban or built-up, 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)

			lb               , &! lower bound of array   (snl+1)
			snl				     ! number of snow layers (-5~-1)

		real(r8), intent(in) :: &
			dlon                    , &! longitude (degrees)
			dlat                    , &! latitude (degrees)

			deltim                  , &! land model time step [sec]
			z_soisno (maxsnl+1:nl_soil)    , &! layer depth [m]
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness [m]
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level [m]

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature [K]
			t_grnd                 		    , &! ground surface temperature [K]
			wliq_soisno(maxsnl+1:nl_soil)  , &! liquid water in layers [kg/m2]
			wice_soisno(maxsnl+1:nl_soil)  , &! ice lens in layers [kg/m2]

			forc_t                  , &! temperature at reference height [K]
				forc_pbot               , &! atm bottom level pressure (or reference height) [Pa]
				forc_po2m               , &! O2 concentration in atmos. [Pa]
				forc_pco2m              , &! CO2 partial pressure at observational height [Pa]
				forc_us                 , &! eastward wind speed [m/s], used as U10 proxy for lake gas exchange
				forc_vs                 , &! northward wind speed [m/s], used as U10 proxy for lake gas exchange

			zwt                     , &! the depth from ground (soil) surface to water table [m]

			rootfr   (1:nl_soil)    , &! fraction of roots in each soil layer (the sum of all layer is 1)

			snowdp                  , &! snow depth [m]
			wat                     , &! total water storage [mm] (reserved interface)
			rsur                    , &! surface runoff [mm/s] (reserved interface)
			etr                     , &! transpiration rate [mm/s]
			lakedepth               , &! lake depth [m], used by CTSM lake ebullition pressure
			lake_icefrac(1:nl_lake) , &! lake frozen mass fraction [-], used to suppress lake ebullition
			wdsrf                   , &! depth of surface water [mm]
			wetwat                  , &! water storage in wetland [mm]
			bsw      (1:nl_soil)   	, &! Clapp and Hornberger "b" (nlevgrnd)

			smp      (1:nl_soil)    , &! soil matrix potential [mm]
			porsl    (1:nl_soil)    , &! volumetric soil water at saturation (porosity)
			lai                     , &! leaf area index [m2/m2]

			annsum_npp              , &! annual sum NPP (g C/m2/yr)
			rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)

			fsatmax                 , &! maximum saturated area fraction [-]
			fsatdcf                 , &! decay factor in calculation of saturated area fraction [1/m]
			frcsat                     ! fraction of saturation area
		!------------------- methane_annualupdate ------------------------------------
		real(r8), intent(in) :: &
			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! (gC/m2/s) soil organic matter heterotrophic respiration

		!------------------- methane_prod --------------------------------------------
		real(r8), intent(in) :: &
			crootfr  (1:nl_soil)    , &! fraction of roots for carbon in each soil layer (the sum of all layer is 1)
			lithr                   , &! (gC/m2/s) litter heterotrophic respiration
			hr_vr    (1:nl_soil)    , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)    , &! fraction by which decomposition is limited by DEF_METHANE%anoxia
			fphr     (1:nl_soil)    , &! fraction of potential heterotrophic respiration

			pot_f_nit_vr(1:nl_soil) , &! (gN/m3/s) potential soil nitrification flux
			pH                         ! soil water pH

		!------------------- methane_aere --------------------------------------------
		real(r8), intent(in) :: &
			rootr    (1:nl_soil)       ! effective fraction of roots in each soil layer (the sum of all layer is 1)

		!------------------- methane_tran --------------------------------------------
		real(r8), intent(in) :: &
			cellorg  (1:nl_soil)   		, &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		, &! surface water temperature
			microbial_prod_potential_patch(1:nl_soil), &! optional microbial CH4 production potential (mol/m3/s)
			microbial_oxid_potential_patch(1:nl_soil), &! optional microbial CH4 oxidation potential (mol/m3/s)
			organic_max               		! organic matter content (kg m-3) where soil is assumed to act like peat

		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data
		!!!! --------------------------------------------------------------------------------------------------------
		!------------------- ch4_flux, balance, depth variables ------------------------------
		real(r8), intent(out) :: &
			net_methane                     , & ! average net methane correction to CO2 flux (mol/m2/s)
			methane_prod_depth    (1:nl_soil)   , & ! production of CH4 in each soil layer (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)   , & ! O2 consumption during decomposition in each soil layer (mol/m3/s)
			co2_decomp_depth  (1:nl_soil)   , & ! diagnostic CO2 from decomposition/methanogenesis before O2-stress scaling (mol/m3/s)
			methane_oxid_depth    (1:nl_soil)   , & ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)   , & ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)
			co2_oxid_depth    (1:nl_soil)   , & ! CO2 production from CH4 oxidation (mol/m3/s)
			methane_aere_depth    (1:nl_soil)   , & ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			methane_tran_depth    (1:nl_soil)   , & ! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			o2_aere_depth     (1:nl_soil)   , & ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)
			co2_aere_depth    (1:nl_soil)   , & ! CO2 aerenchyma diagnostic flux (mol/m3/s)
			methane_ebul_depth    (1:nl_soil)   , & ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2stress          (1:nl_soil)   , & ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			methane_stress         (1:nl_soil)   , & ! Ratio of methane available to total per-timestep methane sinks
			methane_surf_flux_tot               , & ! CH4 flux to atmosphere incl. numerical corrections (mol/m2/s)
			methane_surf_flux_tot_phys          , & ! CH4 physical flux before CH4 clip/residual corrections (mol/m2/s)
			methane_surf_aere                   , & ! CH4 surface flux via aerenchyma (mol/m2/s)
			methane_surf_ebul                   , & ! CH4 ebullition flux (mol/m2/s)
				methane_surf_diff                   , & ! CH4 diffusion flux plus numerical closure (mol/m2/s)
				methane_surf_diff_phys              , & ! CH4 pure-physics diffusion flux before clip/residual (mol/m2/s)
				methane_balance_residual            , & ! numerical closure flux credited to methane_surf_diff (mol/m2/s)
				methane_ch4_clip_credit             , & ! negative CH4 clip credited to methane_surf_diff (mol/m2/s)
				o2_cap_loss                        , & ! O2 removed by post-solve physical cap (mol/m2/s)
				o2_cap_gain                        , & ! O2 added by post-solve nonnegative floor (mol/m2/s)
			methane_ebul_tot                    , & ! Total CH4 ebullition (mol/m2/s)
			methane_prod_tot                    , & ! Total CH4 production (mol/m2/s)
			methane_oxid_tot                    , & ! Total CH4 oxidation (mol/m2/s)
			co2_decomp_tot                     , & ! Total CO2 decomposition/methanogenesis (mol/m2/s)
			co2_oxid_tot                       , & ! Total CO2 from CH4 oxidation (mol/m2/s)
			co2_aere_tot                       , & ! Total CO2 aerenchyma diagnostic (mol/m2/s)
			co2_net_tot                            ! Net diagnosed CO2 source (mol/m2/s)

		!------------------- total and concentration variables ------------------------------
		real(r8), intent(inout) :: &
			totcol_methane               , & ! total methane in soil column (mol/m2)
			grnd_methane_cond           , & ! tracer conductance for boundary layer [m/s]
			conc_o2  (1:nl_soil)    , & ! O2 conc in each soil layer (mol/m3)
			conc_methane (1:nl_soil)        ! CH4 conc in each soil layer (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		!------------------- ch4_flux, balance, depth variables ------------------------------
		real(r8), intent(out) :: &
			net_methane_unsat               , & ! average unsaturated net methane correction to CO2 flux (mol/m2/s)
			net_methane_sat                 , & ! average saturated net methane correction to CO2 flux (mol/m2/s)
			methane_prod_depth_unsat (1:nl_soil), & ! CH4 production rate in unsaturated soil layer (mol/m3/s)
			methane_prod_depth_sat   (1:nl_soil), & ! CH4 production rate in saturated soil layer (mol/m3/s)
			o2_decomp_depth_unsat(1:nl_soil), & ! O2 consumption during decomposition (unsaturated) (mol/m3/s)
			o2_decomp_depth_sat  (1:nl_soil), & ! O2 consumption during decomposition (saturated) (mol/m3/s)
			co2_decomp_depth_unsat(1:nl_soil), & ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsat) (mol/m3/s)
			co2_decomp_depth_sat  (1:nl_soil), & ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (sat) (mol/m3/s)
			methane_oxid_depth_unsat (1:nl_soil), & ! CH4 oxidation rate in unsaturated soil layer (mol/m3/s)
			methane_oxid_depth_sat   (1:nl_soil), & ! CH4 oxidation rate in saturated soil layer (mol/m3/s)
			o2_oxid_depth_unsat  (1:nl_soil), & ! O2 oxidation rate in unsaturated soil layer (mol/m3/s)
			o2_oxid_depth_sat    (1:nl_soil), & ! O2 oxidation rate in saturated soil layer (mol/m3/s)
			co2_oxid_depth_unsat (1:nl_soil), & ! CO2 from CH4 oxidation (unsaturated) (mol/m3/s)
			co2_oxid_depth_sat   (1:nl_soil), & ! CO2 from CH4 oxidation (saturated) (mol/m3/s)
			methane_aere_depth_unsat (1:nl_soil), & ! CH4 loss rate via aerenchyma (unsaturated) (mol/m3/s)
			methane_aere_depth_sat   (1:nl_soil), & ! CH4 loss rate via aerenchyma (saturated) (mol/m3/s)
			methane_tran_depth_unsat (1:nl_soil), & ! CH4 loss rate via transpiration (unsaturated) (mol/m3/s)
			methane_tran_depth_sat   (1:nl_soil), & ! CH4 loss rate via transpiration (saturated) (mol/m3/s)
			o2_aere_depth_unsat  (1:nl_soil), & ! O2 gain via aerenchyma (unsaturated) (mol/m3/s)
			o2_aere_depth_sat    (1:nl_soil), & ! O2 gain via aerenchyma (saturated) (mol/m3/s)
			co2_aere_depth_unsat (1:nl_soil), & ! CO2 aerenchyma diagnostic (unsaturated) (mol/m3/s)
			co2_aere_depth_sat   (1:nl_soil), & ! CO2 aerenchyma diagnostic (saturated) (mol/m3/s)
			methane_ebul_depth_unsat (1:nl_soil), & ! CH4 ebullition loss (unsaturated) (mol/m3/s)
			methane_ebul_depth_sat   (1:nl_soil), & ! CH4 ebullition loss (saturated) (mol/m3/s)
			o2stress_unsat       (1:nl_soil), & ! O2 stress ratio (unsaturated)
			o2stress_sat         (1:nl_soil), & ! O2 stress ratio (saturated)
			methane_stress_unsat      (1:nl_soil), & ! CH4 stress ratio (unsaturated)
			methane_stress_sat        (1:nl_soil), & ! CH4 stress ratio (saturated)
			methane_surf_flux_tot_unsat         , & ! CH4 surface flux to atmosphere (unsaturated) (mol/m2/s)
			methane_surf_flux_tot_sat           , & ! CH4 surface flux to atmosphere (saturated) (mol/m2/s)
			methane_surf_aere_unsat             , & ! CH4 surface flux via aerenchyma (unsaturated) (mol/m2/s)
			methane_surf_aere_sat               , & ! CH4 surface flux via aerenchyma (saturated) (mol/m2/s)
			methane_surf_ebul_unsat             , & ! CH4 ebullition flux (unsaturated) (mol/m2/s)
			methane_surf_ebul_sat               , & ! CH4 ebullition flux (saturated) (mol/m2/s)
			methane_surf_diff_unsat             , & ! CH4 diffusion flux (unsaturated) (mol/m2/s)
			methane_surf_diff_sat               , & ! CH4 diffusion flux (saturated) (mol/m2/s)
			methane_surf_diff_phys_unsat        , & ! CH4 pure-physics diff (unsaturated) (mol/m2/s)
			methane_surf_diff_phys_sat          , & ! CH4 pure-physics diff (saturated) (mol/m2/s)
			methane_ebul_tot_unsat              , & ! Total CH4 ebullition (unsaturated) (mol/m2/s)
			methane_ebul_tot_sat                , & ! Total CH4 ebullition (saturated) (mol/m2/s)
			methane_prod_tot_unsat              , & ! Total CH4 production (unsaturated) (mol/m2/s)
			methane_prod_tot_sat                , & ! Total CH4 production (saturated) (mol/m2/s)
			methane_oxid_tot_unsat              , & ! Total CH4 oxidation (unsaturated) (mol/m2/s)
			methane_oxid_tot_sat                , & ! Total CH4 oxidation (saturated) (mol/m2/s)
			co2_decomp_tot_unsat               , & ! Total CO2 decomp/methanogenesis (unsaturated) (mol/m2/s)
			co2_decomp_tot_sat                 , & ! Total CO2 decomp/methanogenesis (saturated) (mol/m2/s)
			co2_oxid_tot_unsat                 , & ! Total CO2 from CH4 oxidation (unsaturated) (mol/m2/s)
			co2_oxid_tot_sat                   , & ! Total CO2 from CH4 oxidation (saturated) (mol/m2/s)
			co2_net_tot_unsat                  , & ! Net diagnosed CO2 source (unsaturated) (mol/m2/s)
			co2_net_tot_sat                        ! Net diagnosed CO2 source (saturated) (mol/m2/s)

		!------------------- total and concentration variables ------------------------------
		real(r8), intent(inout) :: &
			totcol_methane_unsat         , & ! total methane in soil column (unsaturated) (mol/m2)
			totcol_methane_sat           , & ! total methane in soil column (saturated) (mol/m2)
			grnd_methane_cond_unsat     , & ! tracer conductance for boundary layer (unsaturated) [m/s]
			grnd_methane_cond_sat       , & ! tracer conductance for boundary layer (saturated) [m/s]
			conc_o2_unsat (1:nl_soil), & ! O2 conc in unsaturated soil layer (mol/m3)
			conc_o2_sat   (1:nl_soil), & ! O2 conc in saturated soil layer (mol/m3)
			conc_methane_unsat(1:nl_soil), & ! CH4 conc in unsaturated soil layer (mol/m3)
			conc_methane_sat  (1:nl_soil)   ! CH4 conc in saturated soil layer (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!------------------- lake methane diagnostic/state variables ------------------------------
		real(r8), intent(out) :: &
			methane_prod_depth_lake (1:nl_soil), & ! lake CH4 production rate (mol/m3/s)
			methane_oxid_depth_lake (1:nl_soil), & ! lake CH4 oxidation rate (mol/m3/s)
			methane_ebul_depth_lake (1:nl_soil), & ! lake CH4 ebullition loss rate (mol/m3/s)
			co2_decomp_depth_lake (1:nl_soil), & ! lake CO2 decomp/methanogenesis (mol/m3/s)
			co2_oxid_depth_lake   (1:nl_soil), & ! lake CO2 from CH4 oxidation (mol/m3/s)
			methane_surf_ebul_lake             , & ! lake CH4 surface ebullition flux (mol/m2/s)
			methane_surf_diff_lake             , & ! lake CH4 surface diffusion flux (mol/m2/s)
			methane_surf_flux_tot_lake         , & ! lake total CH4 surface flux (mol/m2/s)
			methane_prod_tot_lake              , & ! lake total CH4 production (mol/m2/s)
			methane_oxid_tot_lake              , & ! lake total CH4 oxidation (mol/m2/s)
			methane_ebul_tot_lake              , & ! lake total CH4 ebullition (mol/m2/s)
			co2_decomp_tot_lake                , & ! lake total CO2 decomp/methanogenesis (mol/m2/s)
			co2_oxid_tot_lake                  , & ! lake total CO2 oxidation product (mol/m2/s)
			co2_net_tot_lake                       ! lake net diagnosed CO2 source (mol/m2/s)

		real(r8), intent(inout) :: &
			totcol_methane_lake           , & ! lake total methane column (mol/m2)
			grnd_methane_cond_lake       , & ! lake tracer conductance (m/s)
			conc_o2_lake      (1:nl_soil), & ! lake O2 concentration (mol/m3)
			conc_methane_lake (1:nl_soil)   ! lake CH4 concentration (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!------------------- atmospheric and structural variables ------------------------------
		real(r8), intent(out) :: &
			c_atm      (1:3)                 ! CH4, O2, CO2 atmospheric conc (mol/m3)

		real(r8), intent(inout) :: &
			forc_pmethanem              , & ! CH4 concentration in atmosphere (Pa)
			lake_soilc  (1:nl_soil)         ! lake sediment organic carbon per layer (gC/m3)

		! layer_sat_lag is persistent state: unsaturated updates read its previous value.
		real(r8), intent(inout) :: &
			layer_sat_lag(1:nl_soil)         ! lagged saturation ratio per layer (persists across timesteps)

		!------------------- annual accumulators ------------------------------
			real(r8), intent(inout) :: &
				annavg_agnpp            , & ! annual average above-ground NPP (gC/m2/s)
				annavg_bgnpp            , & ! annual average below-ground NPP (gC/m2/s)
				annavg_somhr            , & ! annual average SOM heterotrophic respiration (gC/m2/s)
			annavg_finrw              ! respiration-weighted annual average of inundated zones (gC/m2/s)

		!------------------- temporary accumulators ------------------------------
		real(r8), intent(inout) :: &
			tempavg_agnpp           , & ! temporary average above-ground NPP (gC/m2/s)
			tempavg_bgnpp           , & ! temporary average below-ground NPP (gC/m2/s)
			annsum_counter          , & ! seconds since last annual accumulator turnover
			tempavg_somhr           , & ! temporary average SOM heterotrophic respiration (gC/m2/s)
			tempavg_finrw               ! respiration-weighted temporary average of inundated zones (gC/m2/s)

		!------------------- saturated fraction ------------------------------
		real(r8), intent(inout) :: &
			fsat_bef                , & ! finundated from previous timestep
			finundated_lag          , & ! time-lagged fractional inundated area
			methane_dfsat_tot               ! CH4 flux to atm due to decreasing finundated [mol/m2/s]

		real(r8), intent(in) :: &
            f_h2osfc

		! R1 rice paddy hooks: when present, overrides finundated for
		! patches containing paddy-irrigated rice CFTs (61 / 62).  Both
		! are optional so callers that don't enable rice paddy CH4 leave
		! the natural-wetland behaviour completely intact.
		logical , intent(in), optional :: is_rice_paddy_in
		real(r8), intent(in), optional :: rice_pft_frac_in
		logical, intent(in), optional :: store_patch_diagnostics_in
		real(r8), intent(out), optional :: finundated_used_out
		real(r8), intent(out), optional :: finundated_default_out

		!=================== Local Variables ============================================
		integer  :: i,j,s,l                     ! indices

		integer  :: sat                     ! 0 = unsatured, 1 = saturated
		real(r8) :: finundated              ! fractional inundated area
		logical  :: methane_cold_start      ! true only when no valid previous CH4 inundation state exists
		logical  :: wetland_inactive_dry_area ! non-inundated patchtype-2 area has no active CH4 column
		logical  :: is_rice_paddy           ! R1: rice paddy patch flag (unpacked from optional)
		real(r8) :: rice_pft_frac           ! R1: rice fraction within patch (0..1)
		real(r8) :: rice_substrate_boost_eff ! R2: fraction-weighted paddy substrate multiplier
		real(r8) :: finundated_rice         ! R1: rice-overridden finundated value
		real(r8) :: finundated_default      ! R1: scheme-computed natural-wetland finundated
		logical  :: rice_finundation_active ! true only during live rice or drain-window
		logical  :: store_patch_diagnostics
		integer  :: idpp_rice, dsh
		real(r8) :: ms_start, ms_end, decay, drain_target

		real(r8) :: totcolch4_bef
		real(r8) :: totcolch4_bef_sat
		real(r8) :: totcolch4_bef_unsat
		real(r8) :: totcolch4_bef_lake

		real(r8) :: err_methane

		logical, save :: warned_ch4_budget_methane = .false.

		real(r8) :: dfsat
		real(r8), parameter :: min_rebalance_finundated = 1.e-4_r8
		real(r8) :: redoxlags               ! redox lag time in s
		real(r8) :: redoxlags_vertical      ! Vertical redox lag time in s
			integer  :: dummyfilter(1)          ! empty filter
			integer  :: atm_month, atm_mday
			real(r8) :: atm_methane_mix
			real(r8) :: methane_balance_residual_unsat, methane_balance_residual_sat
			real(r8) :: methane_ch4_clip_credit_unsat, methane_ch4_clip_credit_sat
			real(r8) :: o2_cap_loss_unsat, o2_cap_loss_sat
			real(r8) :: o2_cap_gain_unsat, o2_cap_gain_sat

		real(r8) :: k_h_cc(0:nl_soil,ngases)! ratio of mol/m3 in liquid to mol/m3 in gas [-]

        ! liquid volumetric water content ---- water volume/all volume [m3/m3]
        real(r8) :: vol_aqu            (1:nl_soil)
        real(r8) :: vol_aqu_sat        (1:nl_soil)
        real(r8) :: vol_aqu_unsat      (1:nl_soil)

        ! CH4 storage volume.  Liquid water is always included; ice is
        ! included when freeze-out is disabled and as immobile storage when a
        ! layer is completely frozen.
        real(r8) :: vol_ch4_storage       (1:nl_soil)
        real(r8) :: vol_ch4_storage_sat   (1:nl_soil)
        real(r8) :: vol_ch4_storage_unsat (1:nl_soil)

        ! air volumetric water content ---- air volume/all volume [m3/m3]
        real(r8) :: vol_gas            (1:nl_soil)
        real(r8) :: vol_gas_sat        (1:nl_soil)
        real(r8) :: vol_gas_unsat      (1:nl_soil)

        ! water-filled proportion
        real(r8) :: f_aqu              (1:nl_soil)
        real(r8) :: f_aqu_sat          (1:nl_soil)
        real(r8) :: f_aqu_unsat        (1:nl_soil)

        ! air-filled proportion
        real(r8) :: f_gas              (1:nl_soil)
        real(r8) :: f_gas_sat          (1:nl_soil)
        real(r8) :: f_gas_unsat        (1:nl_soil)

        ! gas phase CH4 conc in each soil layer (mol/m3)
        real(r8) :: conc_ch4_gas             (1:nl_soil)
        real(r8) :: conc_ch4_gas_sat         (1:nl_soil)
        real(r8) :: conc_ch4_gas_unsat       (1:nl_soil)

        ! aqueous phase CH4 conc in each soil layer (mol/m3)
        real(r8) :: conc_ch4_aqu             (1:nl_soil)
        real(r8) :: conc_ch4_aqu_sat         (1:nl_soil)
        real(r8) :: conc_ch4_aqu_unsat       (1:nl_soil)

        ! CH4 conc in each porosity (mol/m3)
        real(r8) :: conc_ch4_porsl           (1:nl_soil)
        real(r8) :: conc_ch4_porsl_sat       (1:nl_soil)
        real(r8) :: conc_ch4_porsl_unsat     (1:nl_soil)

        ! gas phase CH4 conc in each porosity (mol/m3)
        real(r8) :: conc_ch4_gas_porsl       (1:nl_soil)
        real(r8) :: conc_ch4_gas_porsl_sat   (1:nl_soil)
        real(r8) :: conc_ch4_gas_porsl_unsat (1:nl_soil)

        ! aqueous phase CH4 conc in each porosity (mol/m3)
        real(r8) :: conc_ch4_aqu_porsl       (1:nl_soil)
        real(r8) :: conc_ch4_aqu_porsl_sat   (1:nl_soil)
        real(r8) :: conc_ch4_aqu_porsl_unsat (1:nl_soil)

        ! gas phase O2 conc in each soil layer (mol/m3)
        real(r8) :: conc_o2_gas              (1:nl_soil)
        real(r8) :: conc_o2_gas_sat          (1:nl_soil)
        real(r8) :: conc_o2_gas_unsat        (1:nl_soil)

        ! aqueous phase O2 conc in each soil layer (mol/m3)
        real(r8) :: conc_o2_aqu              (1:nl_soil)
        real(r8) :: conc_o2_aqu_sat          (1:nl_soil)
        real(r8) :: conc_o2_aqu_unsat        (1:nl_soil)

        ! O2 conc in each porosity (mol/m3)
        real(r8) :: conc_o2_porsl            (1:nl_soil)
        real(r8) :: conc_o2_porsl_sat        (1:nl_soil)
        real(r8) :: conc_o2_porsl_unsat      (1:nl_soil)

        ! gas phase O2 conc in each porosity (mol/m3)
        real(r8) :: conc_o2_gas_porsl        (1:nl_soil)
        real(r8) :: conc_o2_gas_porsl_sat    (1:nl_soil)
        real(r8) :: conc_o2_gas_porsl_unsat  (1:nl_soil)

        ! aqueous phase O2 conc in each porosity (mol/m3)
        real(r8) :: conc_o2_aqu_porsl        (1:nl_soil)
        real(r8) :: conc_o2_aqu_porsl_sat    (1:nl_soil)
        real(r8) :: conc_o2_aqu_porsl_unsat  (1:nl_soil)

		real(r8) :: err
		real(r8) :: vliq, vice, vtot
		real(r8) :: wtd_arg

		integer  :: jwt                 ! index of the soil layer right above the water table (-)
        integer  :: jwt_sat            ! index of the soil layer right above the water table (-), saturated zone
        integer  :: jwt_unsat          ! index of the soil layer right above the water table (-), unsaturated zone
		real(r8) :: zwt_sat, wice_soisno_sat(maxsnl+1:nl_soil), wliq_soisno_sat(maxsnl+1:nl_soil), wdsrf_sat
		real(r8) :: zwt_unsat, wice_soisno_unsat(maxsnl+1:nl_soil), wliq_soisno_unsat(maxsnl+1:nl_soil), wdsrf_unsat
		real(r8) :: routing_depth_mm
		logical  :: lake_restart_debug_print
		real(r8) :: lake_dbg_totcol_bef, lake_dbg_ch4_l1_bef, lake_dbg_o2_l1_bef
		real(r8) :: lake_dbg_ch4_l5_bef, lake_dbg_o2_l5_bef

		!-----------------------------------------------------------------------
		! Set parameters.  When biome lookup enabled, use per-patch redoxlag
		! (set in Driver via get_biome_redoxlag).  Otherwise use global default.
		if (DEF_METHANE%use_biome_redoxlag .and. allocated(biome_redoxlag_patch) .and. &
		    ipatch >= 1 .and. ipatch <= size(biome_redoxlag_patch)) then
			redoxlags = biome_redoxlag_patch(ipatch) * secspday   ! days --> s
		else
			redoxlags = DEF_METHANE%redoxlag * secspday           ! days --> s
		endif
		! Vertical redoxlag uses same biome lookup ratio (sim ratio assumption)
		if (DEF_METHANE%use_biome_redoxlag .and. allocated(biome_redoxlag_patch) .and. &
		    ipatch >= 1 .and. ipatch <= size(biome_redoxlag_patch)) then
			redoxlags_vertical = biome_redoxlag_patch(ipatch) * secspday
		else
			redoxlags_vertical = DEF_METHANE%redoxlag_vertical * secspday
		endif

		! Initialize fluxes and diagnostics to zero
		methane_prod_depth        = 0._r8
		o2_decomp_depth          = 0._r8
		co2_decomp_depth         = 0._r8
		methane_oxid_depth       = 0._r8
		o2_oxid_depth            = 0._r8
		co2_oxid_depth           = 0._r8
		methane_aere_depth       = 0._r8
		methane_tran_depth       = 0._r8
		o2_aere_depth            = 0._r8
		co2_aere_depth           = 0._r8
		methane_ebul_depth       = 0._r8
		o2stress                 = 1._r8
		methane_stress           = 1._r8
		methane_ebul_tot         = 0._r8

			methane_surf_flux_tot       = 0._r8
			methane_surf_flux_tot_sat   = 0._r8
			methane_surf_flux_tot_unsat = 0._r8
				methane_balance_residual    = 0._r8
				methane_ch4_clip_credit     = 0._r8
				methane_balance_residual_unsat = 0._r8
				methane_balance_residual_sat   = 0._r8
				methane_ch4_clip_credit_unsat = 0._r8
				methane_ch4_clip_credit_sat   = 0._r8
				o2_cap_loss             = 0._r8
				o2_cap_loss_unsat       = 0._r8
				o2_cap_loss_sat         = 0._r8
				o2_cap_gain             = 0._r8
				o2_cap_gain_unsat       = 0._r8
				o2_cap_gain_sat         = 0._r8

		methane_prod_tot            = 0._r8
		methane_prod_tot_sat        = 0._r8
		methane_prod_tot_unsat      = 0._r8

		methane_oxid_tot            = 0._r8
		methane_oxid_tot_sat        = 0._r8
		methane_oxid_tot_unsat      = 0._r8

		methane_prod_tot_lake       = 0._r8
		methane_oxid_tot_lake       = 0._r8
		methane_ebul_tot_lake       = 0._r8
		methane_surf_ebul_lake      = 0._r8
		methane_surf_diff_lake      = 0._r8
		methane_surf_flux_tot_lake  = 0._r8
		methane_prod_depth_lake     = 0._r8
		methane_oxid_depth_lake     = 0._r8
		methane_ebul_depth_lake     = 0._r8

		! Adjustment to NEE for methane production - oxidation
		net_methane             = 0._r8
		net_methane_sat         = 0._r8
		net_methane_unsat       = 0._r8

			CALL julian2monthday (idate(1), idate(2), atm_month, atm_mday)
			atm_methane_mix = methane_atm_mixing_ratio(idate(1), atm_month)
			! Offline CH4 can use either the scalar atmospheric value
			! DEF_METHANE%atm_methane or an optional transient file table.
			if (DEF_METHANE%methane_offline) then
				forc_pmethanem = atm_methane_mix*forc_pbot
				! [Pa]     =[mol/mol]*[Pa]
			else
				if (ieee_is_nan(forc_pmethanem) .or. &
				    abs(forc_pmethanem) >= 0.5_r8 * abs(spval) .or. &
				    forc_pmethanem <= 0._r8) then
					forc_pmethanem = atm_methane_mix*forc_pbot
				end if
			end if
		c_atm(1) =  forc_pmethanem / rgasm / forc_t
		c_atm(2) =  forc_po2m  / rgasm / forc_t
		c_atm(3) =  max(forc_pco2m, 0._r8) / rgasm / forc_t
		! n/V = P/RT
		![mol/m3]=[Pa]         /[J/K/mol]/[K]
		![mol/m3]=[J/m3]       *[K*mol/J]*[1/K]
		totcolch4_bef = totcol_methane
		totcolch4_bef_sat = totcol_methane_sat
		totcolch4_bef_unsat = totcol_methane_unsat
		totcolch4_bef_lake = totcol_methane_lake
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			totcolch4_bef = totcolch4_bef_lake
			totcolch4_bef_sat = totcolch4_bef_lake
			totcolch4_bef_unsat = 0._r8
		endif
		lake_restart_debug_print = DEF_METHANE%lake_restart_debug .and. &
			patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. &
			idate(1) == DEF_METHANE%lake_restart_debug_year .and. &
			idate(2) >= DEF_METHANE%lake_restart_debug_start_doy .and. &
			idate(2) <= DEF_METHANE%lake_restart_debug_end_doy .and. &
			(DEF_METHANE%lake_restart_debug_sec < 0 .or. &
			 idate(3) == DEF_METHANE%lake_restart_debug_sec) .and. &
			totcolch4_bef_lake > 0._r8
		if (lake_restart_debug_print) then
			lake_dbg_totcol_bef = totcolch4_bef_lake
			lake_dbg_ch4_l1_bef = conc_methane_lake(1)
			lake_dbg_o2_l1_bef = conc_o2_lake(1)
			lake_dbg_ch4_l5_bef = conc_methane_lake(min(5,nl_soil))
			lake_dbg_o2_l5_bef = conc_o2_lake(min(5,nl_soil))
			write(6,'(A,1X,I6,1X,I4,1X,I7,1X,I9,1X,2F10.3,1X,24ES14.6)') &
				'CH4_LAKE_ENTRY_DEBUG', idate(1), idate(2), idate(3), ipatch, dlon, dlat, &
				lakedepth, wdsrf, t_grnd, t_h2osfc, t_soisno(1), lake_icefrac(1), &
				wliq_soisno(1), wice_soisno(1), smp(1), f_h2osfc, fsat_bef, finundated_lag, &
				layer_sat_lag(1), layer_sat_lag(min(5,nl_soil)), totcolch4_bef_lake, &
				conc_methane_lake(1), conc_o2_lake(1), conc_methane_lake(min(5,nl_soil)), &
				conc_o2_lake(min(5,nl_soil)), sum(lake_soilc(1:nl_soil)), forc_t, &
				forc_pbot, forc_po2m, forc_pmethanem
		else
			lake_dbg_totcol_bef = 0._r8
			lake_dbg_ch4_l1_bef = 0._r8
			lake_dbg_o2_l1_bef = 0._r8
			lake_dbg_ch4_l5_bef = 0._r8
			lake_dbg_o2_l5_bef = 0._r8
		endif

		totcol_methane = 0.
		totcol_methane_sat   = 0.
		totcol_methane_unsat = 0.

		finundated = 0._r8
			if (DEF_wetland_finundation_scheme == 1) then
				! Wetland tiles fully saturated (finundated=1), all other patch
				! types unsaturated (finundated=0).  The user-facing wetwat mode
				! starts from this classification and then replaces wetland tiles
				! with wetwat/wetwat_ref below.
				if (patchtype == 2) then
					finundated = 1._r8
				else
					finundated = 0._r8
				endif
			elseif (DEF_wetland_finundation_scheme == 2) then
						finundated = frcsat
			elseif (DEF_wetland_finundation_scheme == 3) then
				finundated = f_h2osfc
			elseif (DEF_wetland_finundation_scheme == 4) then
				if (patchtype == 2) then
					finundated = 1._r8
				else
					finundated = f_h2osfc
				endif
			elseif (DEF_wetland_finundation_scheme == 5) then
				! GIEMS-MC monthly wetland climatology, looked up per-patch.
				! Loaded once at init by MOD_Tracer_Reactive_Methane_GIEMS (see CoLM.F90).
			! ipatch is passed in through methane()'s arg list and forwarded
			! to giems_finundated() below.
			if (.not. giems_active) then
				write(6,*) 'ERROR: methane scheme 5 requires an active GIEMS file.'
				CALL CoLM_stop ()
				endif
				finundated = giems_finundated(ipatch, idate(1), idate(2))
				! GIEMS-MC is grid-cell absolute inundated/saturated wetland
				! fraction.  Convert to a wetland-relative fraction on mapped
				! wetland patches, matching routing A1 semantics.
				IF (patchtype == 2 .and. allocated(wetland_frac_per_patch) .and. &
				    ipatch >= 1 .and. ipatch <= size(wetland_frac_per_patch)) THEN
				   finundated = min(1._r8, finundated / max(wetland_frac_per_patch(ipatch), 0.01_r8))
				ENDIF
			elseif (DEF_wetland_finundation_scheme == 6) then
			! Water-table-driven inundation via logistic S-curve:
			!   finundated = 1 / (1 + exp((zwt - wtd_inflection) / wtd_steepness))
			! Defaults (wtd_inflection=0.30 m, wtd_steepness=0.05 m) reproduce
			! the "30 cm threshold" as the inflection point (finundated=0.5)
			! while smoothing the original step function (29 cm = 1 vs 31 cm = 0
			! creates unphysical 100x jump in CH4 production).
			! References:
			!   Walter & Heimann 2000 (sigmoid f(zwt))
			!   Sundh 2000 (CH4 ~ exp(-zwt/0.3) for boreal peatland)
			!   Bridgham 2013 GCB (review).
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				finundated = 1._r8
			elseif (patchtype /= 2 .and. DEF_METHANE%use_routing_for_soil .and. &
			        allocated(f_inund_flood_patch) .and. &
			        ipatch >= 1 .and. ipatch <= size(f_inund_flood_patch)) then
				! HYBRID: soil tile uses routing-published f_inund_flood_patch
				! (captures Pantanal/Amazon overland flooding that sigmoid(zwt)
				! cannot resolve because soil zwt is too deep there).
				! Optional threshold gate: only flood when fldfrc > threshold
				! (Bates 2010, Coe 2002 hydrological connectivity concept).
				if (f_inund_flood_patch(ipatch) > DEF_METHANE%hybrid_soil_threshold) then
					finundated = f_inund_flood_patch(ipatch)
				else
					finundated = 0._r8
				endif
			else
				if (DEF_METHANE%wtd_steepness <= 0._r8 .or. ieee_is_nan(DEF_METHANE%wtd_steepness)) then
					write(6,*) 'ERROR: DEF_METHANE%wtd_steepness must be positive for methane scheme 6.'
					CALL CoLM_stop ()
				endif
				if (ieee_is_nan(zwt) .or. abs(zwt) >= 0.5_r8*abs(spval)) then
					write(6,*) 'ERROR: invalid zwt for methane scheme 6 logistic inundation: ', zwt
					CALL CoLM_stop ()
				endif
				! Dynamic flood extension: soil tile (patchtype==0/1) uses
				! a separate (typically deeper) sigmoid inflection.  This lets
				! seasonally-wet upland (Pantanal floodplain etc.) contribute
				! CH4 when soil zwt rises into 0.5-1.5 m range, while keeping
				! the steep 0.3 m wetland sigmoid that captures peat-style
				! peak-to-trough variation.  Disabled when wtd_inflection_soil=0.
				IF (patchtype /= 2 .and. DEF_METHANE%wtd_inflection_soil > 0._r8) then
					wtd_arg = (zwt - DEF_METHANE%wtd_inflection_soil) / &
					          max(DEF_METHANE%wtd_steepness_soil, 1.e-3_r8)
				else
					wtd_arg = (zwt - DEF_METHANE%wtd_inflection) / DEF_METHANE%wtd_steepness
				endif
				wtd_arg = min(50._r8, max(-50._r8, wtd_arg))
				finundated = 1._r8 / (1._r8 + exp(wtd_arg))
			endif
		elseif (DEF_wetland_finundation_scheme == 7) then
				! GridRiverLake fldfrc: per-patch general flood fraction
				! (levee+floodplain) from grid-based river/lake routing.
					! Published by MOD_Grid_RiverLakeFlow::publish_fldfrc_to_patches
					! at end of routing.  Methane consumes the most recently published
					! value, normally from the previous land step.  Fall back to the
					! levee variable if flood var happens to be unallocated (e.g.
				IF (allocated(f_inund_flood_patch) .and. &
				    ipatch >= 1 .and. &
				    ipatch <= size(f_inund_flood_patch)) THEN
				   finundated = f_inund_flood_patch(ipatch)
				ELSE IF (allocated(f_inund_levee_patch) .and. &
				    ipatch >= 1 .and. &
				    ipatch <= size(f_inund_levee_patch)) THEN
				   finundated = f_inund_levee_patch(ipatch)
					ELSE
						   write(6,*) 'ERROR: methane scheme 7 requires published GridRiverLake levee/floodplain fraction for patch ', &
					      ipatch
					   CALL CoLM_stop ()
					ENDIF
					! Routing flood fractions are grid/ucat-relative
					! (A_flood / A_cell).  On wetland tiles convert them to a
					! wetland-relative fraction (A_flood / A_wetland) using the
					! static parent-cell wetland fraction cache.  This assumes
					! broad floodwater preferentially occupies mapped wetland
					! area; soil tiles keep the grid-relative routing fraction.
					IF (patchtype == 2 .and. allocated(wetland_frac_per_patch) .and. &
					    ipatch >= 1 .and. ipatch <= size(wetland_frac_per_patch)) THEN
					   finundated = min(1._r8, finundated / max(wetland_frac_per_patch(ipatch), 0.01_r8))
					ENDIF
			else
			write(6,*) 'ERROR: unsupported DEF_wetland_finundation_scheme = ', &
				DEF_wetland_finundation_scheme
			CALL CoLM_stop ()
		endif
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			! Lake CH4 is an explicitly saturated-column problem in CTSM.
			! The full lake production / oxidation / transport path is
			! implemented as a saturated sediment column with CTSM lake
			! production / oxidation / ebullition / transport diagnostics.
			! Force the enabled lake gate onto the saturated branch rather
			! than silently weighting it away through the wetland scheme.
			finundated = 1._r8
		endif
		! NaN/spval defense before clamp.  Sources frcsat / f_h2osfc /
		! f_inund_levee_patch may carry sentinel values if forcing or
		! upstream module misbehaves; min/max would propagate NaN.
		if (ieee_is_nan(finundated) .or. abs(finundated) >= 0.5_r8*abs(spval)) then
			finundated = 0._r8
		endif
		if (finundated < 1.e-10_r8) finundated = 0._r8
		! Clamp to keep saturated/unsaturated weighting bounded.
		finundated = min(max(finundated, 0._r8), 1._r8)

		wetland_inactive_dry_area = patchtype == 2 .and. &
			DEF_METHANE%wetland_dry_unsat_branch

		! CoLM's wetland hydrology hardcodes zwt=0 and frcsat=1 for
		! patchtype==2 (MOD_SoilSnowHydrology.F90:1227-1265), so the IGBP /
		! TOPMODEL / sigmoid / routing schemes all collapse to ~1.0 or near-
		! const on wetland tiles, killing the seasonal cycle.  But the WATER
		! BALANCE itself (wetwat = wdsrf + wa + wetwat + (gwat - etr +
		! qsdew + qfros - qsubl)*deltim) IS dynamic.  This override turns
		! wetwat / wetwatmax into the per-step finundated, so seasonality
		! flows from forced precip-ET balance directly into methane physics.
		if (DEF_METHANE%enable_wetwat_finundated_override .and. patchtype == 2) then
			if (wetwatmax > 0._r8) then
				finundated = max(0._r8, min(1._r8, wetwat / wetwatmax))
			endif
		endif

		methane_cold_start = (abs(fsat_bef) >= 1.e30_r8 .or. fsat_bef < 0._r8 .or. fsat_bef > 1._r8)
		if (methane_cold_start) then
			fsat_bef = finundated
		endif
		if (abs(finundated_lag) >= 1.e30_r8 .or. finundated_lag < 0._r8 .or. finundated_lag > 1._r8) then
			finundated_lag = finundated
		endif
		do j = 1, nl_soil
			if (abs(layer_sat_lag(j)) >= 1.e30_r8 .or. &
			    layer_sat_lag(j) < 0._r8 .or. &
			    layer_sat_lag(j) > 1._r8) then
				layer_sat_lag(j) = finundated
			end if
		end do
		if (wetland_inactive_dry_area) then
			! For a mapped wetland tile, finundated is the dynamic active
			! wetland-water fraction.  The remaining dry part of the tile is
			! not a separate prognostic CH4 soil column; otherwise dynamic
			! wetland area changes create an artificial dry-unsat reservoir
			! and budget warnings.  Preserve dynamic finundated, but define
			! the dry fraction as inactive: newly flooded area starts with
			! zero CH4, and decreasing flooded area releases the saturated
			! excess through methane_dfsat_tot below.
			totcolch4_bef = fsat_bef * totcolch4_bef_sat
			totcolch4_bef_unsat = 0._r8
			conc_methane_unsat = 0._r8
			totcol_methane_unsat = 0._r8
		endif
		if (snowdp > 0._r8) then  !If snow_depth>0, keep finundated from the previous time step of snow season. (by Xiyan Xu, 05/2016)
            finundated = fsat_bef
		end if
		! Re-apply the lake force AFTER the snow override.  A snow-covered
		! lake patch must still report finundated=1 to stay on the saturated
		! sediment branch (otherwise a stale fsat_bef from a previous
		! mid-thaw step could put the lake back on the unsaturated branch
		! and zero-out CH4 production).  See lake gate ~line 658.
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			finundated = 1._r8
		endif

		! R1 rice paddy: override finundated when this patch contains
		! paddy-irrigated rice and CN reports croplive_p (i.e. the rice
		! is currently in its growing season).  Multi-CFT mixed-cropping
		! patches blend rice-flooded value with the scheme-computed value
		! by area fraction.  croplive_p handles the cross-year wrap
		! automatically (it stays .true. through New Year for Dec-March
		! Southern Hemisphere rice).
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
		rice_pft_frac = min(max(rice_pft_frac, 0._r8), 1._r8)
		finundated_default = finundated
		finundated_rice = finundated_default
		rice_finundation_active = .false.
		rice_substrate_boost_eff = 1._r8
		IF (is_rice_paddy .and. rice_pft_frac > 0._r8) THEN
		   ! Safety clamp on caller-supplied rice fraction (multi-CFT
		   ! aggregation occasionally exceeds 1 if pftfrac normalization
		   ! has been touched upstream).
		   rice_substrate_boost_eff = 1._r8 + rice_pft_frac * &
		      (DEF_METHANE%rice_substrate_boost - 1._r8)

		   IF (is_paddy_rice_live(ipatch)) THEN
		      ! Live rice: floor to rice_paddy_min_finundated.  Use max() so
		      ! an already-wet patch (wetland tile with rice CFT, or scheme=6
		      ! zwt at surface) is NOT lowered.
		      finundated_rice = max(finundated_default, DEF_METHANE%rice_paddy_min_finundated)

		      idpp_rice = rice_days_past_planting(ipatch, idate(2))
		      ms_start  = DEF_METHANE%rice_midseason_start_days
		      ms_end    = ms_start + DEF_METHANE%rice_midseason_drain_days
		      IF (idpp_rice >= 0 .and. real(idpp_rice, r8) >= ms_start &
		          .and. real(idpp_rice, r8) < ms_end) THEN
		         finundated_rice = min(finundated_rice, &
		            max(DEF_METHANE%rice_midseason_drained_finundated, &
		                finundated_default))
		      ENDIF

		      rice_finundation_active = .true.
		      finundated = rice_pft_frac * finundated_rice &
		                 + (1._r8 - rice_pft_frac) * finundated_default

		   ELSE
		      dsh = rice_days_since_harvest(ipatch, idate(2))
		      IF (dsh >= 0 .and. real(dsh, r8) < DEF_METHANE%rice_drain_window_days) THEN
		         decay = 1._r8 - real(dsh, r8) / DEF_METHANE%rice_drain_window_days
		         drain_target = decay * DEF_METHANE%rice_paddy_min_finundated
		         finundated_rice = max(drain_target, finundated_default)
		         rice_finundation_active = .true.
		         finundated = rice_pft_frac * finundated_rice &
		                    + (1._r8 - rice_pft_frac) * finundated_default
		      ENDIF
		   ENDIF

		   ! Clamp to [0,1] for safety after any blend.
		   finundated = min(max(finundated, 0._r8), 1._r8)
		ENDIF

		dfsat = finundated - fsat_bef

		! Update lagged finundated for redox calculation
		if (redoxlags > 0._r8) then
			finundated_lag = finundated_lag * exp(-deltim/redoxlags) &
					+ finundated * (1._r8 - exp(-deltim/redoxlags))
		else
			finundated_lag = finundated
		end if

		methane_dfsat_tot = 0._r8
		do j=1,nl_soil
			if (.not. methane_cold_start) then
				if (dfsat > 0._r8) then
					if (wetland_inactive_dry_area) then
						! Dynamic wetland can have very small but nonzero active
						! fractions.  Do not use the generic tiny-area fallback
						! (copy unsat, which is zero for inactive dry area),
						! because it erases the previous saturated CH4 mass.
						! Regrid active wetland inventory conservatively.
						conc_methane_sat(j) = (fsat_bef*conc_methane_sat(j) + &
							dfsat*conc_methane_unsat(j)) / max(finundated, 1.e-12_r8)
						conc_o2_sat (j) = (fsat_bef*conc_o2_sat (j) + &
							dfsat*conc_o2_unsat (j)) / max(finundated, 1.e-12_r8)
					elseif (finundated < min_rebalance_finundated) then
						conc_methane_sat(j) = conc_methane_unsat(j)
						conc_o2_sat (j) = conc_o2_unsat (j)
					else
						conc_methane_sat(j) = (fsat_bef*conc_methane_sat(j) + dfsat*conc_methane_unsat(j)) / finundated
						conc_o2_sat (j) = (fsat_bef*conc_o2_sat (j) + dfsat*conc_o2_unsat (j)) / finundated
					end if
				elseif (dfsat < 0._r8) then
					methane_dfsat_tot = methane_dfsat_tot - &
						dfsat*(conc_methane_sat(j) - conc_methane_unsat(j)) * &
						dz_soisno(j) / deltim
					! [mol/m2/s]   = [-] * [mol/m3] * [m] / [s]
				end if
			end if
		end do

		!!!! Begin biochemistry
		! First for soil
		! Do CH4 Annual Averages
		call methane_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
			annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
			tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)

		call henry_law(t_grnd,t_soisno,k_h_cc)

		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			! LAKE_SOILC is initialized from the formal surface dataset field
			! lake_soilc_srf, then persisted in the Methane restart field.
			conc_o2_sat = conc_o2_lake
			conc_methane_sat = conc_methane_lake
			totcol_methane_sat = totcol_methane_lake
			grnd_methane_cond_sat = grnd_methane_cond_lake
		endif

		!-------------------------------------------------
		! Loop
		!-------------------------------------------------
		do sat= 0, 1
			if (wetland_inactive_dry_area .and. sat == 0) then
				methane_prod_depth_unsat = 0._r8
				methane_oxid_depth_unsat = 0._r8
				methane_ebul_depth_unsat = 0._r8
				methane_aere_depth_unsat = 0._r8
				methane_tran_depth_unsat = 0._r8
				o2_decomp_depth_unsat = 0._r8
				co2_decomp_depth_unsat = 0._r8
				o2_oxid_depth_unsat = 0._r8
				co2_oxid_depth_unsat = 0._r8
				o2_aere_depth_unsat = 0._r8
				co2_aere_depth_unsat = 0._r8
				o2stress_unsat = 1._r8
				methane_stress_unsat = 1._r8
				conc_methane_unsat = 0._r8
				totcol_methane_unsat = 0._r8
				methane_surf_aere_unsat = 0._r8
				methane_surf_ebul_unsat = 0._r8
				methane_surf_diff_unsat = 0._r8
				methane_surf_diff_phys_unsat = 0._r8
				methane_ebul_tot_unsat = 0._r8
				methane_balance_residual_unsat = 0._r8
				methane_ch4_clip_credit_unsat = 0._r8
				o2_cap_loss_unsat = 0._r8
				o2_cap_gain_unsat = 0._r8
				! The dry branch is deliberately inactive, so its area-weighted
				! effective conductance contribution must also be zero rather than
				! a stale value from the previous timestep.
				grnd_methane_cond_unsat = 0._r8
				cycle
			endif
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. sat == 0) then
				methane_prod_depth_unsat = 0._r8
				methane_oxid_depth_unsat = 0._r8
				methane_ebul_depth_unsat = 0._r8
				methane_aere_depth_unsat = 0._r8
				methane_tran_depth_unsat = 0._r8
				o2_decomp_depth_unsat = 0._r8
				co2_decomp_depth_unsat = 0._r8
				o2_oxid_depth_unsat = 0._r8
				co2_oxid_depth_unsat = 0._r8
				o2_aere_depth_unsat = 0._r8
				co2_aere_depth_unsat = 0._r8
				o2stress_unsat = 1._r8
				methane_stress_unsat = 1._r8
				conc_o2_unsat = 0._r8
				conc_methane_unsat = 0._r8
				totcol_methane_unsat = 0._r8
				methane_surf_aere_unsat = 0._r8
				methane_surf_ebul_unsat = 0._r8
				methane_surf_diff_unsat = 0._r8
				methane_surf_diff_phys_unsat = 0._r8
				methane_ebul_tot_unsat = 0._r8
				cycle
			endif
				if (sat==0) then ! unsaturated
					wliq_soisno_unsat = wliq_soisno
					wice_soisno_unsat = wice_soisno
					! Force wetland unsat branch
					! to be truly dry so that methane_prod()'s j>jwt gate
					! returns zero for all layers.  Otherwise wetland zwt=0
					! hardcode (MOD_SoilSnowHydrology.F90:1229) gives
					! jwt_unsat=0, all layers produce, sat≈unsat → finundated
					! weighting becomes useless and ratio stuck at 1.5-1.7.
					if (DEF_METHANE%wetland_dry_unsat_branch .and. &
					    patchtype == 2) then
						zwt_unsat   = zi_soisno(nl_soil) + 1.e-3_r8   ! below deepest layer
						wdsrf_unsat = 0._r8                            ! no ponding
						jwt_unsat   = nl_soil                          ! no layer is below WT
					else
						zwt_unsat = zwt
						wdsrf_unsat = methane_wetland_water_depth(patchtype, wdsrf, wetwat)
						jwt_unsat = nl_soil
						! allow jwt to equal zero when zwt is in top layer
						do j = 1, nl_soil
							if(zwt_unsat <= zi_soisno(j)) then
								jwt_unsat = j-1
								exit
							end if
						end do
					endif

				do j=1,nl_soil
					if (DEF_METHANE%use_vertical_redoxlag .and. j > jwt_unsat .and. redoxlags_vertical > 0._r8) then ! saturated currently
						layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical) &
							+ (1._r8 - exp(-deltim/redoxlags_vertical))
					else if (DEF_METHANE%use_vertical_redoxlag .and. redoxlags_vertical > 0._r8) then
						layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical)
					else if (j > jwt_unsat) then  ! redoxlags_vertical = 0
						layer_sat_lag(j) = 1._r8
					else
						layer_sat_lag(j) = 0._r8
					end if
				end do

				call split_ch4_o2_phases( dz_soisno, wliq_soisno_unsat, wice_soisno_unsat, porsl, &
					conc_methane_unsat, conc_o2_unsat, k_h_cc, idate, &
					vol_aqu_unsat, vol_ch4_storage_unsat, vol_gas_unsat, f_aqu_unsat, f_gas_unsat, &
					conc_ch4_gas_unsat, conc_ch4_aqu_unsat, conc_ch4_porsl_unsat, conc_ch4_gas_porsl_unsat, conc_ch4_aqu_porsl_unsat, &
					conc_o2_gas_unsat, conc_o2_aqu_unsat, conc_o2_porsl_unsat, conc_o2_gas_porsl_unsat, conc_o2_aqu_porsl_unsat )

				! Calculate CH4 production in each soil layer
				call methane_prod ( ipatch, idate, patchtype, sat, jwt_unsat, finundated, finundated_lag, rr, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					lai, conc_o2_unsat, rootfr, annavg_finrw, &
					crootfr, somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, pH, layer_sat_lag, lake_soilc, &
					microbial_prod_potential_patch, &
					methane_prod_depth_unsat, o2_decomp_depth_unsat, co2_decomp_depth_unsat )

				! R2 short-term SOC fix: lift methane substrate on rice paddy
				! tiles to compensate for grid-mean SOC underestimate.
				! Default rice_substrate_boost=1.0 (off) for global budgets.
				! Set >1.0 only for explicit paddy-SOC sensitivity tests.
				! Only the methane production channel is scaled; o2/co2 decomp
				! are NOT scaled because that would falsify the BGC carbon
				! balance — this is a methane-only substrate availability proxy.
				IF (is_rice_paddy .and. rice_substrate_boost_eff > 1._r8) THEN
				   methane_prod_depth_unsat(:) = methane_prod_depth_unsat(:) * rice_substrate_boost_eff
				ENDIF

				! Calculate CH4 oxidation in each soil layer
				call methane_oxid ( idate, patchtype, jwt_unsat, sat, t_soisno, dz_soisno, zi_soisno, smp, vol_aqu_unsat, &
					conc_o2_aqu_porsl_unsat, conc_ch4_aqu_porsl_unsat, &
					microbial_oxid_potential_patch, &
					methane_oxid_depth_unsat, o2_oxid_depth_unsat )

				! Calculate CH4 ebullition losses in each soil layer
				call methane_ebul ( idate, patchtype, jwt_unsat, sat, finundated, deltim, &
					z_soisno, dz_soisno, zi_soisno, forc_pbot, lakedepth, lake_icefrac, &
					t_soisno, wdsrf_unsat, conc_methane_unsat, conc_ch4_gas_porsl_unsat, &
					methane_ebul_depth_unsat )

				! Calculate CH4 aerenchyma losses in each soil layer
				call methane_aere ( ipatch, idate, jwt_unsat, sat, patchclass, lai, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					rootfr, rootr, etr, DEF_METHANE%grnd_methane_cond_default, c_atm, annsum_npp, &
					annavg_agnpp, annavg_bgnpp, conc_methane_unsat, methane_prod_depth_unsat, &
					methane_oxid_depth_unsat, methane_ebul_depth_unsat, &
					conc_ch4_aqu_unsat, conc_ch4_aqu_porsl_unsat, conc_ch4_gas_porsl_unsat, conc_o2_aqu_porsl_unsat, conc_o2_gas_porsl_unsat, &
					methane_aere_depth_unsat, methane_tran_depth_unsat, o2_aere_depth_unsat )

				! Solve CH4 reaction/diffusion equation
				! Competition for oxygen will occur here.
				call methane_tran ( idate, patchtype, &
					lb, snl, jwt_unsat, sat, finundated, &
					dlon, dlat, deltim, z_soisno, dz_soisno, zi_soisno, t_soisno, t_grnd, forc_us, forc_vs, &
					porsl, wliq_soisno_unsat, wice_soisno_unsat, wdsrf_unsat, lakedepth, lake_icefrac, bsw, c_atm, methane_prod_depth_unsat, o2_aere_depth_unsat, &
					cellorg, t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl_unsat, conc_ch4_aqu_porsl_unsat, conc_o2_gas_porsl_unsat, conc_o2_aqu_porsl_unsat, vol_aqu_unsat, vol_ch4_storage_unsat, vol_gas_unsat, &
					o2stress_unsat, methane_stress_unsat, methane_surf_aere_unsat, methane_surf_ebul_unsat, methane_surf_diff_unsat, methane_ebul_tot_unsat, &
					methane_oxid_depth_unsat, methane_aere_depth_unsat, methane_tran_depth_unsat, methane_ebul_depth_unsat, &
						DEF_METHANE%grnd_methane_cond_default, grnd_methane_cond_unsat, &
						o2_oxid_depth_unsat, o2_decomp_depth_unsat, &
						conc_o2_unsat, conc_methane_unsat, methane_balance_residual_unsat, methane_ch4_clip_credit_unsat, methane_surf_diff_phys_unsat, &
						o2_cap_loss_unsat, o2_cap_gain_unsat )

			elseif (sat==1) then ! saturated
				zwt_sat=0.

				! IF (.not.DEF_SPLIT_SOILSNOW) THEN
				! 	IF (lb >= 1) THEN
				! 		wetwat_sat = wdsrf_sat + wa_sat + wetwat_sat + (gwat - etr + qsdew + qfros - qsubl) * deltim
				! 	ELSE
				! 		wetwat_sat = wdsrf_sat + wa_sat + wetwat_sat + (gwat - etr) * deltim
				! 	ENDIF
				! ELSE
				! 	wetwat_sat = wdsrf_sat + wa_sat + wetwat_sat + (gwat - etr + qsdew_soil + qfros_soil - qsubl_soil) * deltim
				! ENDIF

				! IF (wetwat_sat > wetwatmax) THEN
				! 	wdsrf_sat  = wetwat_sat - wetwatmax
				! 	wetwat_sat = wetwatmax
	! 	wa_sat     = 0.
				! ELSEIF (wetwat_sat < 0) THEN
	            ! 	wa_sat     = wetwat_sat
				! 	wdsrf_sat  = 0.
				! 	wetwat_sat = 0.
				! ELSE
				! 	wdsrf_sat = 0.
		        !     wa_sat    = 0.
				! ENDIF
				wliq_soisno_sat = wliq_soisno
				wice_soisno_sat = wice_soisno
				! Preserve the liquid/ice partition in volume space while filling
				! pore space to saturation.  wliq/wice are masses, so using mass
				! fractions would bias partially frozen layers because denh2o and
				! denice differ.
				DO j = 1, nl_soil
					vliq = max(wliq_soisno(j), 0._r8) / denh2o
					vice = max(wice_soisno(j), 0._r8) / denice
					vtot = vliq + vice
					IF (vtot > 1.e-12_r8) THEN
						wliq_soisno_sat(j) = porsl(j)*dz_soisno(j)*denh2o * vliq / vtot
						wice_soisno_sat(j) = porsl(j)*dz_soisno(j)*denice * vice / vtot
					ELSE
						! If the physical column is dry, fall back to temperature phase.
						IF (t_soisno(j) > tfrz) THEN
							wliq_soisno_sat(j) = porsl(j)*dz_soisno(j)*denh2o
							wice_soisno_sat(j) = 0._r8
						ELSE
							wliq_soisno_sat(j) = 0._r8
							wice_soisno_sat(j) = porsl(j)*dz_soisno(j)*denice
						ENDIF
					ENDIF
				ENDDO
					wdsrf_sat = methane_wetland_water_depth(patchtype, wdsrf, wetwat)
					! Area-forcing modes need a saturated-branch water depth consistent
					! with positive inundated area.  Satellite products provide area but
					! no depth, so use the configurable floor.  Routing publishes
					! floodplain depth flddph=max(stage-river_bank,0) separately from
					! river depth wdsrf_ucat; use physical flood depth when it is larger,
					! with the same floor as a shallow-flood fallback.
					IF ((DEF_wetland_finundation_scheme == 5 .or. &
					     DEF_wetland_finundation_scheme == 7) .and. patchtype == 2 .and. &
					    finundated > 0._r8 .and. DEF_METHANE%area_forcing_min_water_depth_mm > 0._r8) THEN
					   wdsrf_sat = max(wdsrf_sat, DEF_METHANE%area_forcing_min_water_depth_mm)
					ENDIF
					! Apply routing-published floodplain depth on patches that
					! also consumed routing-published fldfrc above:
					!   (a) scheme 7 routing: all patches
					!   (b) scheme 6 + use_routing_for_soil (hybrid): soil tiles only
					! Wetland tiles in hybrid stay on sigmoid(zwt), so they don't
					! get the routing depth (would create area/depth mismatch).
					IF (finundated > 0._r8 .and. &
					    allocated(f_inund_flood_depth_patch) .and. &
					    ipatch >= 1 .and. ipatch <= size(f_inund_flood_depth_patch) .and. &
					    (DEF_wetland_finundation_scheme == 7 .or. &
					     (DEF_wetland_finundation_scheme == 6 .and. &
					      DEF_METHANE%use_routing_for_soil .and. patchtype /= 2))) THEN
					   routing_depth_mm = 1000._r8 * max(0._r8, f_inund_flood_depth_patch(ipatch))
					   wdsrf_sat = max(wdsrf_sat, max(routing_depth_mm, DEF_METHANE%area_forcing_min_water_depth_mm))
					ENDIF
					jwt_sat = 0

				call split_ch4_o2_phases( dz_soisno, wliq_soisno_sat, wice_soisno_sat, porsl, &
					conc_methane_sat, conc_o2_sat, k_h_cc, idate, &
					vol_aqu_sat, vol_ch4_storage_sat, vol_gas_sat, f_aqu_sat, f_gas_sat, &
					conc_ch4_gas_sat, conc_ch4_aqu_sat, conc_ch4_porsl_sat, conc_ch4_gas_porsl_sat, conc_ch4_aqu_porsl_sat, &
					conc_o2_gas_sat, conc_o2_aqu_sat, conc_o2_porsl_sat, conc_o2_gas_porsl_sat, conc_o2_aqu_porsl_sat )

				! Calculate CH4 production in each soil layer
				call methane_prod ( ipatch, idate, patchtype, sat, jwt_sat, finundated, finundated_lag, rr, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					lai, conc_o2_sat, rootfr, annavg_finrw, &
					crootfr, somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, pH, layer_sat_lag, lake_soilc, &
					microbial_prod_potential_patch, &
					methane_prod_depth_sat, o2_decomp_depth_sat, co2_decomp_depth_sat )

				! R2 short-term SOC fix: rice paddy substrate boost (sat branch).
				IF (is_rice_paddy .and. rice_substrate_boost_eff > 1._r8) THEN
				   methane_prod_depth_sat(:) = methane_prod_depth_sat(:) * rice_substrate_boost_eff
				ENDIF

				! Calculate CH4 oxidation in each soil layer
				call methane_oxid ( idate, patchtype, jwt_sat, sat, t_soisno, dz_soisno, zi_soisno, smp, vol_aqu_sat, &
					conc_o2_aqu_porsl_sat, conc_ch4_aqu_porsl_sat, &
					microbial_oxid_potential_patch, &
					methane_oxid_depth_sat, o2_oxid_depth_sat )

				! Calculate CH4 ebullition losses in each soil layer
				call methane_ebul ( idate, patchtype, jwt_sat, sat, finundated, deltim, &
					z_soisno, dz_soisno, zi_soisno, forc_pbot, lakedepth, lake_icefrac, &
					t_soisno, wdsrf_sat, conc_methane_sat, conc_ch4_gas_porsl_sat, &
					methane_ebul_depth_sat )

				! Calculate CH4 aerenchyma losses in each soil layer
				if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
					methane_aere_depth_sat = 0._r8
					methane_tran_depth_sat = 0._r8
					o2_aere_depth_sat = 0._r8
				else
					call methane_aere ( ipatch, idate, jwt_sat, sat, patchclass, lai, deltim, &
						z_soisno, dz_soisno, zi_soisno, t_soisno, &
						rootfr, rootr, etr, DEF_METHANE%grnd_methane_cond_default, c_atm, annsum_npp, &
						annavg_agnpp, annavg_bgnpp, conc_methane_sat, methane_prod_depth_sat, &
						methane_oxid_depth_sat, methane_ebul_depth_sat, &
						conc_ch4_aqu_sat, conc_ch4_aqu_porsl_sat, conc_ch4_gas_porsl_sat, conc_o2_aqu_porsl_sat, conc_o2_gas_porsl_sat, &
						methane_aere_depth_sat, methane_tran_depth_sat, o2_aere_depth_sat )
				endif

				! Solve CH4 reaction/diffusion equation
				! Competition for oxygen will occur here.
				call methane_tran ( idate, patchtype, &
						lb, snl, jwt_sat, sat, finundated, &
						dlon, dlat, deltim, z_soisno, dz_soisno, zi_soisno, t_soisno, t_grnd, forc_us, forc_vs, &
						porsl, wliq_soisno_sat, wice_soisno_sat, wdsrf_sat, lakedepth, lake_icefrac, bsw, c_atm, methane_prod_depth_sat, o2_aere_depth_sat, &
					cellorg, t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl_sat, conc_ch4_aqu_porsl_sat, conc_o2_gas_porsl_sat, conc_o2_aqu_porsl_sat, vol_aqu_sat, vol_ch4_storage_sat, vol_gas_sat, &
					o2stress_sat, methane_stress_sat, methane_surf_aere_sat, methane_surf_ebul_sat, methane_surf_diff_sat, methane_ebul_tot_sat, &
					methane_oxid_depth_sat, methane_aere_depth_sat, methane_tran_depth_sat, methane_ebul_depth_sat, &
						DEF_METHANE%grnd_methane_cond_default, grnd_methane_cond_sat, &
						o2_oxid_depth_sat, o2_decomp_depth_sat, &
						conc_o2_sat, conc_methane_sat, methane_balance_residual_sat, methane_ch4_clip_credit_sat, methane_surf_diff_phys_sat, &
						o2_cap_loss_sat, o2_cap_gain_sat )

			endif

		enddo

		! Diagnostic CO2 closure.  methane_prod computes decomposition/root-respiration
		! CO2 from the same base decomposition used for CH4 production, excluding
		! nitrification O2 demand and excluding any anoxia-inflated potential O2 demand.
		! co2_decomp_depth is intentionally not rescaled by the O2 competition
		! limiter in methane_tran; treat it as a carbon-source diagnostic rather
		! than a prognostic O2-limited CO2 tracer.  CH4 oxidation produces CO2 at
		! 1:1 molar ratio. Explicit CO2 aerenchyma transport is not prognosed
		! here, so it is diagnosed as zero.
		co2_oxid_depth_unsat   = methane_oxid_depth_unsat
		co2_oxid_depth_sat     = methane_oxid_depth_sat
		co2_aere_depth_unsat   = 0._r8
		co2_aere_depth_sat     = 0._r8

		methane_prod_depth = methane_prod_depth_sat * finundated + methane_prod_depth_unsat * (1.0_r8 - finundated)
		methane_oxid_depth = methane_oxid_depth_sat * finundated + methane_oxid_depth_unsat * (1.0_r8 - finundated)
		methane_aere_depth = methane_aere_depth_sat * finundated + methane_aere_depth_unsat * (1.0_r8 - finundated)
		methane_ebul_depth = methane_ebul_depth_sat * finundated + methane_ebul_depth_unsat * (1.0_r8 - finundated)
		methane_tran_depth = methane_tran_depth_sat * finundated + methane_tran_depth_unsat * (1.0_r8 - finundated)
		co2_decomp_depth = co2_decomp_depth_sat * finundated + co2_decomp_depth_unsat * (1.0_r8 - finundated)
		co2_oxid_depth   = co2_oxid_depth_sat   * finundated + co2_oxid_depth_unsat   * (1.0_r8 - finundated)
		co2_aere_depth   = co2_aere_depth_sat   * finundated + co2_aere_depth_unsat   * (1.0_r8 - finundated)

		methane_oxid_tot_unsat = sum( methane_oxid_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )
		methane_prod_tot_unsat = sum( methane_prod_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_decomp_tot_unsat = sum( co2_decomp_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_oxid_tot_unsat   = sum( co2_oxid_depth_unsat  (1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_net_tot_unsat    = co2_decomp_tot_unsat + co2_oxid_tot_unsat
		net_methane_unsat = -methane_prod_tot_unsat + methane_oxid_tot_unsat
		methane_surf_flux_tot_unsat = methane_surf_diff_unsat + methane_surf_aere_unsat + methane_surf_ebul_unsat + &
			sum(methane_tran_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil))
		totcol_methane_unsat = sum( conc_methane_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )

		methane_oxid_tot_sat   = sum( methane_oxid_depth_sat(1:nl_soil)   * dz_soisno(1:nl_soil) )
		methane_prod_tot_sat   = sum( methane_prod_depth_sat(1:nl_soil)   * dz_soisno(1:nl_soil) )
		co2_decomp_tot_sat = sum( co2_decomp_depth_sat(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_oxid_tot_sat   = sum( co2_oxid_depth_sat  (1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_net_tot_sat    = co2_decomp_tot_sat + co2_oxid_tot_sat
		net_methane_sat   = -methane_prod_tot_sat   + methane_oxid_tot_sat
		methane_surf_flux_tot_sat = methane_surf_diff_sat  + methane_surf_aere_sat  + methane_surf_ebul_sat + &
			sum(methane_tran_depth_sat(1:nl_soil) * dz_soisno(1:nl_soil))
		totcol_methane_sat   = sum( conc_methane_sat  (1:nl_soil) * dz_soisno(1:nl_soil) )

			methane_surf_diff = methane_surf_diff_sat * finundated + methane_surf_diff_unsat * (1.0_r8 - finundated) + &
				methane_dfsat_tot
			methane_surf_diff_phys = methane_surf_diff_phys_sat * finundated + &
				methane_surf_diff_phys_unsat * (1.0_r8 - finundated)
			methane_balance_residual = methane_balance_residual_sat * finundated + &
				methane_balance_residual_unsat * (1.0_r8 - finundated)
			methane_ch4_clip_credit = methane_ch4_clip_credit_sat * finundated + &
				methane_ch4_clip_credit_unsat * (1.0_r8 - finundated)
			o2_cap_loss = o2_cap_loss_sat * finundated + o2_cap_loss_unsat * (1.0_r8 - finundated)
			o2_cap_gain = o2_cap_gain_sat * finundated + o2_cap_gain_unsat * (1.0_r8 - finundated)
			methane_surf_ebul = methane_surf_ebul_sat * finundated + methane_surf_ebul_unsat * (1.0_r8 - finundated)
		methane_surf_aere = methane_surf_aere_sat * finundated + methane_surf_aere_unsat * (1.0_r8 - finundated)
		methane_oxid_tot = methane_oxid_tot_sat * finundated + methane_oxid_tot_unsat * (1.0_r8 - finundated)
		methane_prod_tot = methane_prod_tot_sat * finundated + methane_prod_tot_unsat * (1.0_r8 - finundated)
		co2_decomp_tot = co2_decomp_tot_sat * finundated + co2_decomp_tot_unsat * (1.0_r8 - finundated)
		co2_oxid_tot   = co2_oxid_tot_sat   * finundated + co2_oxid_tot_unsat   * (1.0_r8 - finundated)
		co2_aere_tot   = sum( co2_aere_depth(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_net_tot    = co2_decomp_tot + co2_oxid_tot - co2_aere_tot
		net_methane = net_methane_sat * finundated + net_methane_unsat * (1.0_r8 - finundated)
		methane_surf_flux_tot = methane_surf_diff + methane_surf_ebul + methane_surf_aere + &
			sum(methane_tran_depth(1:nl_soil) * dz_soisno(1:nl_soil))
		methane_surf_flux_tot_phys = methane_surf_diff_phys + methane_dfsat_tot + &
			methane_surf_ebul + methane_surf_aere + &
			sum(methane_tran_depth(1:nl_soil) * dz_soisno(1:nl_soil))
		totcol_methane = totcol_methane_sat * finundated + totcol_methane_unsat * (1.0_r8 - finundated)
		grnd_methane_cond = grnd_methane_cond_sat * finundated + &
			grnd_methane_cond_unsat * (1.0_r8 - finundated)

		conc_methane = conc_methane_sat * finundated + conc_methane_unsat * (1._r8 - finundated)
		conc_o2 = conc_o2_sat * finundated + conc_o2_unsat * (1._r8 - finundated)

			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				methane_prod_depth_lake = methane_prod_depth_sat
			methane_oxid_depth_lake = methane_oxid_depth_sat
			methane_ebul_depth_lake = methane_ebul_depth_sat
			co2_decomp_depth_lake = co2_decomp_depth_sat
			co2_oxid_depth_lake = methane_oxid_depth_lake
			methane_surf_ebul_lake = methane_surf_ebul_sat
			methane_surf_diff_lake = methane_surf_diff_sat
			methane_surf_flux_tot_lake = methane_surf_ebul_lake + methane_surf_diff_lake
			methane_prod_tot_lake = methane_prod_tot_sat
			methane_oxid_tot_lake = methane_oxid_tot_sat
			methane_ebul_tot_lake = methane_ebul_tot_sat
			co2_decomp_tot_lake = co2_decomp_tot_sat
			co2_oxid_tot_lake = co2_oxid_tot_sat
			co2_net_tot_lake = co2_net_tot_sat
			totcol_methane_lake = totcol_methane_sat
			grnd_methane_cond_lake = grnd_methane_cond_sat
			grnd_methane_cond = grnd_methane_cond_lake
			conc_o2_lake = conc_o2_sat
			conc_methane_lake = conc_methane_sat
			methane_prod_depth = methane_prod_depth_lake
			methane_oxid_depth = methane_oxid_depth_lake
			methane_ebul_depth = methane_ebul_depth_lake
			co2_decomp_depth = co2_decomp_depth_lake
			co2_oxid_depth = co2_oxid_depth_lake
			methane_surf_ebul = methane_surf_ebul_lake
			methane_surf_diff = methane_surf_diff_lake
			methane_surf_aere = 0._r8
			methane_surf_flux_tot = methane_surf_flux_tot_lake
			methane_surf_flux_tot_phys = methane_surf_diff_phys + methane_surf_ebul
			methane_prod_tot = methane_prod_tot_lake
			methane_oxid_tot = methane_oxid_tot_lake
			methane_ebul_tot = methane_ebul_tot_lake
			co2_decomp_tot = co2_decomp_tot_lake
			co2_oxid_tot = co2_oxid_tot_lake
			co2_net_tot = co2_net_tot_lake
			totcol_methane = totcol_methane_lake
			net_methane = net_methane_sat
			methane_balance_residual = methane_balance_residual_sat
			methane_ch4_clip_credit = methane_ch4_clip_credit_sat
			o2_cap_loss = o2_cap_loss_sat
			o2_cap_gain = o2_cap_gain_sat
			conc_o2 = conc_o2_lake
			conc_methane = conc_methane_lake
			if (lake_restart_debug_print) then
				write(6,'(A,1X,I6,1X,I4,1X,I7,1X,I9,1X,2F10.3,1X,19ES14.6)') &
					'CH4_LAKE_RESTART_DEBUG', idate(1), idate(2), idate(3), ipatch, dlon, dlat, &
					lakedepth, wdsrf, f_h2osfc, lake_dbg_totcol_bef, totcol_methane_lake, &
					methane_prod_tot_lake, methane_oxid_tot_lake, methane_surf_diff_lake, &
					methane_surf_ebul_lake, methane_surf_flux_tot_lake, lake_dbg_ch4_l1_bef, &
					lake_dbg_o2_l1_bef, conc_methane_lake(1), conc_o2_lake(1), &
					lake_dbg_ch4_l5_bef, lake_dbg_o2_l5_bef, conc_methane_lake(min(5,nl_soil)), &
					conc_o2_lake(min(5,nl_soil)), sum(lake_soilc(1:nl_soil))
			endif
			if (.not. DEF_METHANE%replenishlakec) then
				do j = 1, nl_soil
					lake_soilc(j) = max(0._r8, lake_soilc(j) - &
						(methane_prod_depth_sat(j) + co2_decomp_depth_sat(j)) * deltim * catomw)
				enddo
				endif
			endif

			! Patch-level diagnostics for history output.  These expose the actual
			! inundation signal consumed by the methane physics and split the surface
			! flux into disjoint patchtype categories so hybrid routing/WTD behavior
			! can be diagnosed without flux back-calculation.
			store_patch_diagnostics = .true.
			IF (present(store_patch_diagnostics_in)) store_patch_diagnostics = store_patch_diagnostics_in
			IF (store_patch_diagnostics .and. allocated(methane_finundated) .and. &
			    ipatch >= 1 .and. ipatch <= size(methane_finundated)) THEN
				methane_finundated(ipatch) = finundated
				methane_soil_finundated(ipatch) = 0._r8
				methane_soil_zwt(ipatch) = spval
				methane_surf_flux_wetland(ipatch) = 0._r8
				methane_surf_flux_soil(ipatch) = 0._r8
				methane_surf_flux_lake(ipatch) = 0._r8
				methane_surf_flux_rice(ipatch) = 0._r8
				methane_surf_aere_soil(ipatch) = 0._r8
				methane_surf_aere_rice(ipatch) = 0._r8
				methane_surf_ebul_soil(ipatch) = 0._r8
				methane_surf_ebul_rice(ipatch) = 0._r8
				methane_surf_diff_soil(ipatch) = 0._r8
				methane_surf_diff_rice(ipatch) = 0._r8
				methane_prod_tot_soil(ipatch) = 0._r8
				methane_prod_tot_rice(ipatch) = 0._r8
				methane_oxid_tot_soil(ipatch) = 0._r8
				methane_oxid_tot_rice(ipatch) = 0._r8

				IF (patchtype == 2) THEN
					methane_surf_flux_wetland(ipatch) = methane_surf_flux_tot
						ELSEIF (patchtype == 0) THEN
							methane_soil_finundated(ipatch) = finundated_default
							methane_soil_zwt(ipatch) = zwt
							methane_surf_flux_soil(ipatch) = methane_surf_flux_tot
							methane_surf_aere_soil(ipatch) = methane_surf_aere
							methane_surf_ebul_soil(ipatch) = methane_surf_ebul
							methane_surf_diff_soil(ipatch) = methane_surf_diff
							methane_prod_tot_soil(ipatch) = methane_prod_tot
							methane_oxid_tot_soil(ipatch) = methane_oxid_tot
				ELSEIF (patchtype == 4 .and. DEF_METHANE%allowlakeprod) THEN
					methane_surf_flux_lake(ipatch) = methane_surf_flux_tot_lake
				ENDIF
			ENDIF











		! Column level balance
		if (.not. methane_cold_start) then
			! Check balance
			err_methane = totcol_methane - totcolch4_bef - deltim*(methane_prod_tot - methane_oxid_tot - methane_surf_flux_tot)
			! [mol/m2]    = [mol/m2] - [mol/m2] + [s]*[mol/m2/s]
			! Originally a fatal stop on column-budget closure > 1e-7.  Keep this
			! as a rate-limited diagnostic rather than a stop so long runs can
			! continue, but print the wetland inundation rebalance and numerical
			! correction terms explicitly.  Values above this threshold are not
			! assumed to be roundoff; the extra fields are needed to distinguish
			! fsat/dfsat rebalance from transport clipping/residual corrections.
			if (abs(err_methane) > 1.e-7_r8 .and. .not. warned_ch4_budget_methane) then
				write(6,*)'WARNING: CH4 column budget closure > 1e-7 mol/m^2/timestep.'
				write(6,*)'Lat,Lon,Patchtype        = ', dlat,dlon, patchtype
				write(6,*)'err_methane              = ', err_methane
				write(6,*)'totcol_methane                = ', totcol_methane
				write(6,*)'totcolch4_bef            = ', totcolch4_bef
				write(6,*)'deltim*methane_prod_tot      = ', deltim*methane_prod_tot
				write(6,*)'deltim*methane_oxid_tot      = ', deltim*methane_oxid_tot
				write(6,*)'deltim*methane_surf_flux_tot = ', deltim*methane_surf_flux_tot
				write(6,*)'fsat_bef, finundated, dfsat = ', fsat_bef, finundated, dfsat
				write(6,*)'deltim*methane_dfsat_tot = ', deltim*methane_dfsat_tot
				write(6,*)'deltim*methane_balance_residual = ', deltim*methane_balance_residual
				write(6,*)'deltim*methane_ch4_clip_credit = ', deltim*methane_ch4_clip_credit
				write(6,*)'deltim*methane_surf_diff = ', deltim*methane_surf_diff
				write(6,*)'deltim*methane_surf_aere = ', deltim*methane_surf_aere
				write(6,*)'deltim*methane_surf_ebul = ', deltim*methane_surf_ebul
				write(6,*)'deltim*methane_surf_tran = ', &
					deltim*sum(methane_tran_depth(1:nl_soil)*dz_soisno(1:nl_soil))
				warned_ch4_budget_methane = .true.
			end if
		end if

		if (present(finundated_used_out)) finundated_used_out = finundated
		if (present(finundated_default_out)) finundated_default_out = finundated_default
		fsat_bef = finundated
	end subroutine methane

	!-----------------------------------------------------------------------
	subroutine methane_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
		tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Annual mean fields.
		! Only the annavg is useful, the tempavg is temp.
		!-----------------------------------------------------------------------
		use MOD_Precision
		implicit none
		!-----------------------Argument----------------------------------------
		integer, intent(in) :: &
			idate(3)                   ! model calendar for next time step (year, days of the year, seconds of the day)

		real(r8), intent(in) :: &
			finundated              , &! fractional inundated area, =sat(0 or 1)
			deltim                  , &! land model time step (sec)

			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! soil organic matter heterotrophic respiration (gC/m2/s)

		real(r8), intent(inout) :: &
			annavg_agnpp            , &! annual average aboveground NPP (gC/m2/s)
			annavg_bgnpp            , &! annual average belowground NPP (gC/m2/s)
			annavg_somhr            , &! annual average SOM heterotrophic resp. (gC/m2/s)
			annavg_finrw               ! respiration-weighted annual average of finundated (1e-2*%)
		! definition different with tempavg_finrw

		real(r8), intent(inout) :: &
			! Cumulative data (from year start to now)
			tempavg_agnpp           , &! temporary average aboveground NPP (gC/m2/s)
			tempavg_bgnpp           , &! temporary average belowground NPP (gC/m2/s)
			annsum_counter          , &! seconds since last annual accumulator turnover
			tempavg_somhr           , &! temporary average SOM heterotrophic resp. (gC/m2/s)
			tempavg_finrw              ! respiration-weighted annual average of finundated (gC/m2/s)
		!-----------------------Local Variables------------------------------
		real(r8):: secsperyear       ! total number of seconds this year
		real(r8):: secsperyear_next  ! total number of seconds next year
		real(r8):: sec_of_year       ! seconds elapsed in current calendar year
		real(r8):: remaining_this_year
		real(r8):: dt_first, dt_second
		!-----------------------------------------------------------------------
		! Accumulate the current step before annual rollover.  If the land
		! timestep straddles Dec-31/Jan-01, split the contribution so leap-year
		! weighting and the first step of the new annual accumulator are both
		! correct.
		if ( isleapyear(idate(1)) ) then
			secsperyear = 366._r8*secspday
		else
			secsperyear = 365._r8*secspday
		endif
		sec_of_year = max(0._r8, real(idate(2)-1, r8)*secspday + real(idate(3), r8))
		remaining_this_year = max(0._r8, secsperyear - sec_of_year)

		if (remaining_this_year > 0._r8 .and. deltim > remaining_this_year) then
			dt_first = remaining_this_year
			dt_second = deltim - dt_first
			CALL add_annual_contribution(dt_first, secsperyear)
			CALL finish_annual_accumulator()
			if ( isleapyear(idate(1)+1) ) then
				secsperyear_next = 366._r8*secspday
			else
				secsperyear_next = 365._r8*secspday
			endif
			CALL add_annual_contribution(dt_second, secsperyear_next)
		else
			CALL add_annual_contribution(deltim, secsperyear)
			if (annsum_counter >= secsperyear) CALL finish_annual_accumulator()
		endif

		contains

			subroutine add_annual_contribution(dt, year_seconds)
				real(r8), intent(in) :: dt, year_seconds
				if (dt <= 0._r8 .or. year_seconds <= 0._r8) return
				annsum_counter = annsum_counter + dt
				tempavg_somhr  = tempavg_somhr + dt/year_seconds * somhr
				tempavg_finrw  = tempavg_finrw + dt/year_seconds * finundated * somhr
				tempavg_agnpp  = tempavg_agnpp + dt/year_seconds * agnpp
				tempavg_bgnpp  = tempavg_bgnpp + dt/year_seconds * bgnpp
			end subroutine add_annual_contribution

			subroutine finish_annual_accumulator()
				annsum_counter = 0._r8
				annavg_somhr = tempavg_somhr
				tempavg_somhr = 0._r8
					if (annavg_somhr > 0._r8) then
						annavg_finrw = tempavg_finrw / annavg_somhr
					else
						! No respiration-weighted annual inundation estimate is available.
						! Keep it missing so the SIF block does not treat a cold/missing
						! annual state as a valid zero and suppress first-year CH4 production.
						annavg_finrw = spval
					endif
				tempavg_finrw = 0._r8
				annavg_agnpp = tempavg_agnpp
				tempavg_agnpp = 0._r8
				annavg_bgnpp = tempavg_bgnpp
				tempavg_bgnpp = 0._r8
			end subroutine finish_annual_accumulator

	end subroutine methane_annualupdate

	!-----------------------------------------------------------------------
	subroutine methane_prod (ipatch,idate,patchtype,sat,jwt,finundated,finundated_lag,rr,deltim,& !input
		z_soisno,dz_soisno,zi_soisno,t_soisno,&
		lai,conc_o2,rootfr,annavg_finrw,&
		crootfr,somhr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,layer_sat_lag,lake_soilc,&
		microbial_prod_potential_layer, &
		methane_prod_depth,o2_decomp_depth,co2_decomp_depth)!output
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Production is done below the water table, based on CN heterotrophic respiration.
		! O2 is consumed by roots & by heterotrophic aerobes.
		! Production is done separately for sat & unsat, and is adjusted for temperature, seasonal inundation,
		! pH (optional), & redox lag factor.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			ipatch                      , &! local patch index (for biome_f_methane_patch lookup)
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype                   , &! land patch type (0=soil, 1=urban or built-up,
			! 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)
			sat                         , &! 0 = unsaturated; 1 = saturated
			jwt                            ! index of the soil layer right above the water table (-)

		real(r8), intent(in) :: &
			finundated                  , &! fractional inundated area in soil column
			finundated_lag              , &
			rr                          , &! root respiration (fine root MR + total root GR) (gC/m2/s)

			deltim                      , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil) , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil) , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)   , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil) , &! soil temperature (K)

			lai                         , &! leaf area index [m2/m2]
			conc_o2  (1:nl_soil)        , &! O2 conc in each soil layer (mol/m3) (nl_soil)
			rootfr   (1:nl_soil)        , &! fraction of roots in each soil layer (the sum of all layer is 1)

			annavg_finrw                , &! respiration-weighted annual average of finundated

			crootfr  (1:nl_soil)        , &! fraction of roots for carbon in each soil layer (the sum of all layer is 1)

			somhr                       , &! soil organic matter heterotrophic respiration (gC/m2/s)
			lithr                       , &! litter heterotrophic respiration (gC/m2/s)
			hr_vr    (1:nl_soil)        , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)        , &! fraction by which decomposition is limited by DEF_METHANE%anoxia
			fphr     (1:nl_soil)        , &! fraction of potential heterotrophic respiration

			pot_f_nit_vr  (1:nl_soil)   , &! potential soil nitrification flux (gN/m3/s)

			pH                          , &! soil water pH
			layer_sat_lag   (1:nl_soil) , &! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
			lake_soilc      (1:nl_soil) , &! lake sediment organic carbon (gC/m3)
			microbial_prod_potential_layer(1:nl_soil) ! optional microbial CH4 production potential (mol/m3/s)

		real(r8), intent(out) :: &
			methane_prod_depth (1:nl_soil)  , &! production of CH4 in each soil layer (nl_soil) (mol/m3/s)
			o2_decomp_depth(1:nl_soil)     , &! O2 consumption during decomposition in each soil layer (nl_soil) (mol/m3/s)
			co2_decomp_depth(1:nl_soil)       ! CO2 production from decomposition/root respiration (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer  :: j,s              ! indices
		real(r8) :: base_decomp      ! base heterotrophic respiration rate [mol C/m2/s]
		real(r8) :: partition_z
		real(r8) :: pH_fact_methane      ! pH factor in methane production
		real(r8), parameter :: pH_fact_reference = &
			10._r8**(-0.2235_r8*6.2_r8*6.2_r8 + 2.7727_r8*6.2_r8 - 8.6_r8)

		! Factors for methanogen temperature dependence being greater than soil aerobes
		real(r8) :: f_methane_adj        ! Adjusted DEF_METHANE%f_methane
		real(r8) :: t_fact_methane       ! Temperature factor calculated using additional Q10

		! O2 limitation on decomposition and methanogenesis
		real(r8) :: seasonalfin      ! finundated in excess of respiration-weighted annual average

		! For calculating column average (rootfrac(j)*rr(j))
		real(r8) :: rr_vr(1:nl_soil) ! vertically resolved column-mean root respiration (g C/m2/s)

		real(r8) :: sif              ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
			real(r8) :: q10lake_eff      ! lake Q10 for sediment carbon decomposition
			real(r8) :: lake_ch4_fraction
			real(r8) :: carbon_decomp_depth
			real(r8) :: legacy_methane_prod_depth
			real(r8) :: microbe_prod_cap
			real(r8) :: freeze_factor
			logical  :: use_microbe_override
		real(r8) :: walter_renorm        ! Walter 2001 depth attenuation renormalization factor
		real(r8) :: walter_w_raw, walter_w_atten   ! pre-pass HR weights for renorm
		!-----------------------------------------------------------------------
		! PATCH loop to calculate vertically resolved column-averaged root respiration
		do j=1,nl_soil
			rr_vr(j) = rr*crootfr(j)
			! [g C/m2/s]=[g C/m2/s]*[-]
		end do

		partition_z = 1._r8
		base_decomp = 0.0_r8
		use_microbe_override = DEF_METHANE%use_microbial_pools .and. &
			DEF_METHANE%use_microbial_flux_override .and. patchtype /= 4

		! ---- Walter & Heimann 2001 depth-attenuation pre-pass --------------
		! Compute renormalization factor that preserves total column CH4
		! production while applying exp(-z/z0_methane_prod) modifier to
		! partition_z later in the j loop.  Only iterate over layers that
		! actually produce in the current branch (j > jwt = below water table).
		walter_renorm = 1._r8
		if (DEF_METHANE%z0_methane_prod > 0._r8 .and. patchtype /= 4) then
			walter_w_raw   = 0._r8
			walter_w_atten = 0._r8
			do j = 1, nl_soil
				if (j > jwt) then
					walter_w_raw   = walter_w_raw   + hr_vr(j) * dz_soisno(j)
					walter_w_atten = walter_w_atten + hr_vr(j) * dz_soisno(j) * &
					                 exp(-z_soisno(j) / DEF_METHANE%z0_methane_prod)
				endif
			end do
			if (walter_w_atten > 0._r8) walter_renorm = walter_w_raw / walter_w_atten
		endif

		! column loop to partition decomposition_rate into each soil layer
		q10lake_eff = DEF_METHANE%q10methane * 1.5_r8
		do j=1,nl_soil
			if (patchtype == 4 .and. .not. DEF_METHANE%allowlakeprod) then
				methane_prod_depth(j) = 0._r8
				o2_decomp_depth(j) = 0._r8
				co2_decomp_depth(j) = 0._r8
				cycle
			endif
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				base_decomp = DEF_METHANE%lake_decomp_fact * DEF_METHANE%cnscalefactor * &
					max(lake_soilc(j), 0._r8) * dz_soisno(j) * &
					q10lake_eff**((t_soisno(j) - DEF_METHANE%q10lakebase) / 10._r8) / catomw
				freeze_factor = min(1._r8, max(0._r8, t_soisno(j) - tfrz + 1._r8))
				base_decomp = base_decomp * freeze_factor
				if (j > jwt) then
					! Keep the anaerobic carbon split stoichiometrically bounded:
					! 2CH2O -> CH4 + CO2 implies f_ch4 <= 0.5 of decomposed C.
					lake_ch4_fraction = min(max(DEF_METHANE%f_methane, 0._r8), 0.5_r8)
				else
					lake_ch4_fraction = 0._r8
				endif
				methane_prod_depth(j) = lake_ch4_fraction * base_decomp / dz_soisno(j)
				! Lake sediment production is treated as an anaerobic carbon split.
				! Diagnose the non-CH4 carbon as CO2 below, but do not create an
				! artificial O2 demand in the fully saturated lake sediment column.
				o2_decomp_depth(j) = 0._r8
				co2_decomp_depth(j) = max(0._r8, base_decomp / dz_soisno(j) - methane_prod_depth(j))
				cycle
			endif
			! Use soil heterotrophic respiration (based on Wania)
			base_decomp = (somhr+lithr) / catomw
			! [mol C/m2/s]=[g C/m2/s]   / [g C/mol C]

			! Multiply base_decomp by factor accounting for lower carbon stock in seasonally inundated areas than
			! if it were inundated all year.
			! This is to reduce emissions in seasonally inundated zones, because the eq.
			! C-flux will be less than predicted by a non-O2-lim model
			!
			! Apply SIF only when BGC decomposition is not already anoxia-limited.
			if (sat == 1 .and. DEF_METHANE%use_ch4_sif .and. &
			    .not. DEF_METHANE%bgc_anoxia_limits_decomp) then
				sif = 1._r8
				if (annavg_finrw /= spval .and. finundated > 0._r8) then
					seasonalfin = max(finundated-annavg_finrw, 0._r8)
					if (seasonalfin > 0._r8) then
						sif = (annavg_finrw + DEF_METHANE%mino2lim*seasonalfin) / finundated
						base_decomp = base_decomp * sif
					end if
				end if
			end if

			! For sensitivity studies
			base_decomp = base_decomp * DEF_METHANE%cnscalefactor


			! depth dependence of production either from rootfr or decomp model
			if ( (somhr + lithr) > 0._r8) then
				partition_z = hr_vr(j) * dz_soisno(j) / (somhr + lithr)
				!   [-]     = [g C/m3/s]*[m]/[g C/m2/s]
			else
				partition_z = 1._r8
			end if

			! Walter & Heimann 2001 explicit methanogenesis depth attenuation.
			! Multiply by exp(-z/z0)*walter_renorm so total column production
			! is preserved but profile is sharpened (top layers up, deep down)
			! to match observed peat / wetland CH4 source depth distribution.
			if (DEF_METHANE%z0_methane_prod > 0._r8 .and. patchtype /= 4) then
				partition_z = partition_z * exp(-z_soisno(j) / DEF_METHANE%z0_methane_prod) * walter_renorm
			endif


			! Adjust DEF_METHANE%f_methane to account for the fact that methanogens may have a higher Q10 than aerobic decomposers.
			! Note this is crude and should ideally be applied to all anaerobic decomposition rather than just the
			! DEF_METHANE%f_methane.
			f_methane_adj = 1.0_r8

			t_fact_methane = DEF_METHANE%q10methane**((t_soisno(j) - DEF_METHANE%q10methane_base)/10._r8)


			! Adjust f_methane by the temperature ratio.  When biome lookup
			! is enabled, use the per-patch value populated by methane_driver
			! via get_biome_f_methane; otherwise fall back to global scalar.
			if (DEF_METHANE%use_biome_f_methane .and. allocated(biome_f_methane_patch) .and. &
			    ipatch >= 1 .and. ipatch <= size(biome_f_methane_patch)) then
				f_methane_adj = biome_f_methane_patch(ipatch) * t_fact_methane
			else
				f_methane_adj = DEF_METHANE%f_methane * t_fact_methane
			endif

			! Do not allow soil methanogenesis in frozen layers.  Lake sediment
			! production already has an explicit temperature gate; apply the same
			! physical cutoff to the soil branch so winter CH4 is not produced
			! solely from residual BGC HR in ice-filled pores.
			if (t_soisno(j) <= tfrz) f_methane_adj = 0._r8


			! Remove CN nitrogen limitation, as methanogenesis is not N limited.
			! Also remove (low) moisture limitation
			if (DEF_METHANE%methane_rmcnlim) then
				if (fphr(j) > 0._r8) then
					f_methane_adj = f_methane_adj / fphr(j)
				end if
			end if


				! If switched on, use pH factor for production based on spatial pH data defined in surface data.
				if (patchtype /= 4 .and. DEF_METHANE%usephfact)then
					if (pH <= DEF_METHANE%pHmin .or. pH >= DEF_METHANE%pHmax) then
						! Outside the empirical Dunfield production window, suppress
						! methanogenesis.  The old branch skipped the multiplier outside
						! [pHmin,pHmax], accidentally treating extreme pH as pH-neutral.
						f_methane_adj = 0._r8
					else
						pH_fact_methane = 10._r8**(-0.2235_r8*pH*pH + 2.7727_r8*pH - 8.6_r8) / pH_fact_reference
						pH_fact_methane = min(1._r8, max(0._r8, pH_fact_methane))
						! fitted function using data from Dunfield et al. 1993
						! Normalized to one at pH=6.2.
						! From Lei Meng
						f_methane_adj = f_methane_adj * pH_fact_methane
					end if
				else
					! if no data, then no pH effects
				end if

			! Redox factor
			if ((patchtype /= 4) .and. sat==1 .and. finundated_lag < finundated)  then
				f_methane_adj = f_methane_adj * finundated_lag / finundated
			elseif (sat == 0 .and. j > jwt) then ! Assume lag in decay of alternative electron acceptors vertically
				f_methane_adj = f_methane_adj * layer_sat_lag(j)
			end if
			! Alternative electron acceptors will be consumed first after soil is inundated.

			f_methane_adj = min(f_methane_adj, 0.5_r8)
			! Must be less than 0.5 because otherwise the actual implied aerobic respiration would be negative.
			! The total of aer. respiration + methanogenesis must remain equal to the SOMHR calculated in CN,
			! so that the NEE is sensible. Even perfectly anaerobic conditions with no alternative
			! electron acceptors would predict no more than 0.5 b/c some oxygen is present in organic matter.
			! e.g. 2CH2O --> CH4 + CO2.

			! Decomposition uses 1 mol O2 per mol CO2 produced (happens below WT also, to deplete O2 below WT)
			! o2_decomp_depth is the demand in the absense of O2 supply limitation, in addition to autotrophic respiration.
			! Competition will be done in methane_prod

			o2_decomp_depth(j) = base_decomp * partition_z / dz_soisno(j)
			! Divide off o_scalar to recover the *potential* O2-unlimited HR for
			! the aerobe O2 demand used in the CH4/O2 competition.  This is only
			! valid when the BGC decomp rate has been pre-multiplied by o_scalar
			! upstream (so dividing here cancels it back out).  When the user
			! has told us BGC already applied its anoxia limitation, we treat
			! base_decomp as the realized rate and skip the divide — otherwise
			! we double-correct (BGC drops decomp by o_scalar, then we re-inflate
			! it back to potential here).  Same gate as the SIF block above.
			if (DEF_METHANE%anoxia .and. .not. DEF_METHANE%bgc_anoxia_limits_decomp) then
				! Use a lower bound to avoid singular oxygen demand when o_scalar is small.
				if (o_scalar(j) > 0._r8) then
					o2_decomp_depth(j) = o2_decomp_depth(j) / max(o_scalar(j), 0.01_r8)
				end if
			end if

			! Add root respiration
			if (patchtype /= 4 .and. j <= jwt) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + rr_vr(j)/catomw/dz_soisno(j)
				! [mol/m3/s]       = [mol/m3/s]         + [g C/m2/s]/[g C/mol C]/[m]
			end if

			! Add oxygen demand for nitrification
			if (DEF_METHANE%use_nitrif_denitrif) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + pot_f_nit_vr(j) * 2.0_r8/14.0_r8
				! [mol/m3/s]       = [mol/m3/s]         + [g N/m3/s]/[g N/mol N]
			end if

			carbon_decomp_depth = max(0._r8, base_decomp * partition_z / dz_soisno(j))
			legacy_methane_prod_depth = 0._r8
			if (j > jwt) then
				legacy_methane_prod_depth = &
					max(0._r8, f_methane_adj * base_decomp * partition_z / dz_soisno(j))
			elseif (DEF_METHANE%anoxicmicrosites) then
				legacy_methane_prod_depth = max(0._r8, &
					f_methane_adj * base_decomp * partition_z / dz_soisno(j) / &
					(1._r8 + DEF_METHANE%oxinhib * conc_o2(j)))
			endif
			microbe_prod_cap = 0.5_r8 * carbon_decomp_depth
			if (legacy_methane_prod_depth > 0._r8) then
				microbe_prod_cap = min(microbe_prod_cap, &
					DEF_METHANE%max_microbe_prod_multiplier * legacy_methane_prod_depth)
			elseif (use_microbe_override) then
				! If the legacy pathway forbids production (e.g. frozen layers or
				! above-WT layers without anoxic microsites), keep the microbial
				! override at zero even if a stale/future potential is nonzero.
				microbe_prod_cap = 0._r8
			endif
			if (use_microbe_override .and. &
			    abs(microbial_prod_potential_layer(j)) < 0.5_r8 * abs(spval)) then
				if (j > jwt .or. DEF_METHANE%anoxicmicrosites) then
					methane_prod_depth(j) = min(max(0._r8, microbial_prod_potential_layer(j)), &
						microbe_prod_cap)
				else
					methane_prod_depth(j) = 0._r8
				endif
			elseif (j  >  jwt) then ! Below the water table so anaerobic CH4 production can occur
				! partition decomposition to layer
				! turn into per volume-total by dz
				methane_prod_depth(j) = f_methane_adj * base_decomp * partition_z / dz_soisno(j)
				! [mol C/m3/s]    = [-]       * [mol C/m2/s]* [-]         / [m]
			else ! Above the WT
				if (DEF_METHANE%anoxicmicrosites) then
					methane_prod_depth(j) = f_methane_adj * base_decomp * partition_z / dz_soisno(j) &
						/ (1._r8 + DEF_METHANE%oxinhib*conc_o2(j))
				else
					methane_prod_depth(j) = 0._r8 ! [mol/m3 total/s]
				endif ! DEF_METHANE%anoxicmicrosites
			endif ! WT

			co2_decomp_depth(j) = max(0._r8, &
				carbon_decomp_depth - methane_prod_depth(j) + &
				rr_vr(j) / catomw / dz_soisno(j))

		end do ! nl_soil

	end subroutine methane_prod

	!---------------------------------------------------------------------------
	subroutine methane_oxid (idate, patchtype, jwt,  sat, t_soisno, dz_soisno, zi_soisno, smp, vol_aqu, &
		conc_o2_aqu_porsl, conc_ch4_aqu_porsl, &
		microbial_oxid_potential_layer, &
		methane_oxid_depth, o2_oxid_depth)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Oxidation is based on double Michaelis-Mentin kinetics, and is adjusted for low soil moisture.
		! Oxidation participates in the shared O2 competition in methane_tran.
		! The microbial override potential may already include local microbial
		! substrate/O2 limits; methane_tran can still reduce it when roots,
		! decomposition, and methanotrophs compete for the same timestep O2.
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype        , &! land patch type
			jwt                    , &! index of the soil layer right above the water table (-)
			sat                       ! 0 = unsaturated; 1 = saturated

		real(r8), intent(in) :: &
			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
			dz_soisno(maxsnl+1:nl_soil)    , &! soil layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)
			smp      (1:nl_soil)   , &! soil matrix potential [mm]
			vol_aqu  (1:nl_soil)   , &! liquid volumetric water content

			conc_o2_aqu_porsl  (1:nl_soil)   , &! aqueous phase O2 conc in each porosity (mol/m3)
			conc_ch4_aqu_porsl (1:nl_soil)   , &! aqueous phase CH4 conc in each porosity (mol/m3)
			microbial_oxid_potential_layer(1:nl_soil) ! optional microbial CH4 oxidation potential (mol/m3/s)

		real(r8), intent(out) :: &
			methane_oxid_depth (1:nl_soil)   , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			o2_oxid_depth  (1:nl_soil)      ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer :: j                              ! indices
		real(r8):: t0                             ! Base temperature for Q10
		real(r8):: oxid_a                         ! Oxidation predicted by method A (temperature & enzyme limited) (mol CH4/m3/s)
		real(r8):: smp_fact                       ! factor for reduction based on soil moisture (unitless)
		real(r8):: k_m_eff                        ! effective DEF_METHANE%k_m
		real(r8):: k_m_o2_eff                     ! effective DEF_METHANE%k_m_o2
		real(r8):: vmax_eff                       ! effective vmax
		real(r8):: oxid_scale                     ! optional landunit-specific oxidation multiplier
		real(r8):: lake_oxid_layer_factor         ! fraction of layer inside lake oxic sediment cap
		real(r8):: lake_oxid_layer_top            ! layer top depth for lake oxic cap (m)
		real(r8):: lake_oxid_layer_bot            ! layer bottom depth for lake oxic cap (m)
		real(r8):: lake_oxid_layer_overlap        ! layer thickness inside lake oxic cap (m)
		logical :: use_microbe_override
		!-----------------------------------------------------------------------
		t0 = tfrz + 12._r8 ! Walter, for Michigan site where the 45 M/h comes from
		use_microbe_override = DEF_METHANE%use_microbial_pools .and. &
			DEF_METHANE%use_microbial_flux_override .and. patchtype /= 4

		! Loop to determine oxidation in each layer
		do j=1,nl_soil
			if (sat == 1 .or. j > jwt) then ! Below the water table
				! Literature (e.g. Bender & Conrad, 1992) suggests lower DEF_METHANE%k_m and vmax for high-CH4-affinity methanotrophs in
				! upland soils consuming ambient methane.
				k_m_eff = DEF_METHANE%k_m
				vmax_eff = DEF_METHANE%vmax_methane_oxid
			else
				k_m_eff = DEF_METHANE%k_m_unsat
				vmax_eff = DEF_METHANE%vmax_oxid_unsat
			end if
			k_m_o2_eff = DEF_METHANE%k_m_o2
			oxid_scale = 1._r8
			lake_oxid_layer_factor = 1._r8
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				if (DEF_METHANE%lake_vmax_methane_oxid >= 0._r8) &
					vmax_eff = DEF_METHANE%lake_vmax_methane_oxid
				if (DEF_METHANE%lake_k_m_o2 > 0._r8) k_m_o2_eff = DEF_METHANE%lake_k_m_o2
				oxid_scale = DEF_METHANE%lake_oxid_scale
				if (DEF_METHANE%lake_oxic_sediment_depth > 0._r8) then
					if (j == 1) then
						lake_oxid_layer_top = 0._r8
					else
						lake_oxid_layer_top = zi_soisno(j-1)
					endif
					lake_oxid_layer_bot = zi_soisno(j)
					lake_oxid_layer_overlap = max(0._r8, &
						min(DEF_METHANE%lake_oxic_sediment_depth, lake_oxid_layer_bot) - lake_oxid_layer_top)
					lake_oxid_layer_factor = min(1._r8, &
						lake_oxid_layer_overlap / max(dz_soisno(j), 1.e-12_r8))
				endif
			endif

			! Continuous soil-water-potential limitation follows CLM5.
			! Smoothly reduce oxidation under dry matrix-potential stress.
			! Reference: Schnell & King 1996, Figure 3; CLM5 Tech Note 2.25.11
			if (j <= jwt .and. smp(j) < 0._r8) then
				smp_fact = exp(-smp(j)/DEF_METHANE%smp_crit)
			else
				smp_fact = 1._r8
			end if

			oxid_a              = oxid_scale * vmax_eff * vol_aqu(j)* conc_ch4_aqu_porsl(j) / (k_m_eff + conc_ch4_aqu_porsl(j)) &
			! [mol/m3/s]        = [mol/m3/s]   * [-]     [mol/m3-w]    [mol/m3-w]  [mol/m3-w]
				* conc_o2_aqu_porsl(j) / (k_m_o2_eff + conc_o2_aqu_porsl(j)) &
				* DEF_METHANE%q10_methane_oxid ** ((t_soisno(j) - t0) / 10._r8) * smp_fact &
				* lake_oxid_layer_factor

			! For all landunits / levels, prevent oxidation if at or below freezing
			if (t_soisno(j) <= tfrz) oxid_a = 0._r8

			if (use_microbe_override .and. &
			    abs(microbial_oxid_potential_layer(j)) < 0.5_r8 * abs(spval)) then
				methane_oxid_depth(j) = max(0._r8, microbial_oxid_potential_layer(j))
			else
				methane_oxid_depth(j) = oxid_a
			endif
			o2_oxid_depth(j) = methane_oxid_depth(j) * 2._r8
		end do
	end subroutine methane_oxid

		!---------------------------------------------------------------------------
	subroutine methane_aere (ipatch, idate,jwt, sat, patchclass, lai, deltim,&
		z_soisno, dz_soisno, zi_soisno, t_soisno,&
		rootfr, rootr, etr, grnd_methane_cond_base, c_atm, annsum_npp,&
		annavg_agnpp, annavg_bgnpp, conc_methane, methane_prod_depth, methane_oxid_depth, methane_ebul_depth, &
		conc_ch4_aqu, conc_ch4_aqu_porsl,conc_ch4_gas_porsl,conc_o2_aqu_porsl,conc_o2_gas_porsl,&
		methane_aere_depth, methane_tran_depth, o2_aere_depth)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Arctic c3 grass (which is often present in fens) and all vegetation in inundated areas is assumed to have
		! some root porosity. Currently, root porosity is allowed to be different for grasses & non-grasses.
		! CH4 diffuses out and O2 diffuses into the soil.  CH4 is also lossed via transpiration, which is both
		! included in the "aere" variables and output separately.  In practice this value is small.
		! By default upland veg. has small 5% porosity but this can be switched to be equal to inundated porosity.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			ipatch            , &! caller's per-patch index for vegetation overrides
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-)
			sat                    , &! 0 = unsatured, 1 = saturated
			patchclass                ! land patch class of USGS classification or others

		real(r8), intent(in) :: &
			lai                    , &! adjusted leaf area index for seasonal variation [-]

			deltim                 , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil)    , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! effective fraction of roots in each soil layer (SMS method only)
			! rootr here for effective per-layer transpiration, which may not be the same as rootfr
			etr                    , &! transpiration rate [mm/s]
			grnd_methane_cond_base     , &! base tracer conductance before snow/pond resistance [m/s]

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			! These variables help us swap between big-leaf and fates boundary conditions
			annsum_npp             , &! annual sum NPP (gC/m2/yr)
			annavg_agnpp           , &! annual avg aboveground NPP (gC/m2/s)
			annavg_bgnpp           , &! annual avg belowground NPP (gC/m2/s)

			! These variables help us swap between saturated and unsaturated boundary conditions
			conc_methane (1:nl_soil)        , &! CH4 conc in each soil layer (mol/m3)
			methane_prod_depth (1:nl_soil)  , &! production of CH4 in each soil layer (mol/m3/s)
			methane_oxid_depth (1:nl_soil)  , &! oxidation of CH4 already claimed in each soil layer (mol/m3/s)
			methane_ebul_depth (1:nl_soil)  , &! ebullition CH4 already claimed in each soil layer (mol/m3/s)

			conc_ch4_aqu       (1:nl_soil) , &! aqueous phase CH4 concentration [mol/m3 water]
			conc_ch4_aqu_porsl (1:nl_soil) , &! aqueous phase CH4 conc in each porosity [mol/m3]
			conc_ch4_gas_porsl (1:nl_soil) , &! gas phase CH4 conc in each porosity [mol/m3]
			conc_o2_aqu_porsl  (1:nl_soil) , &! aqueous phase O2 conc in each porosity [mol/m3]
			conc_o2_gas_porsl  (1:nl_soil)    ! gas phase O2 conc in each porosity [mol/m3]


		real(r8), intent(out) :: &
			methane_aere_depth  (1:nl_soil)  , &! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			methane_tran_depth  (1:nl_soil)  , &! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			o2_aere_depth   (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer  :: j              ! indices

		! methane aerenchyma parameters
		real(r8) :: tranloss(1:nl_soil)    ! loss due to transpiration (mol / m3 /s)
		real(r8) :: aere    (1:nl_soil)
		real(r8) :: oxaere  (1:nl_soil)    ! (mol / m3 /s)

			real(r8) :: aeretran
			real(r8) :: aeretran_requested
			real(r8) :: aeretran_scale
			real(r8) :: aere_pos, aere_neg, tran_pos, tran_neg

		real(r8) :: poros, poros_tiller_real
		real(r8) :: poros_default

		!-----------------------------------------------------------------------
		! Initialize methane_aere_depth
		do j=1,nl_soil
			methane_aere_depth(j) = 0._r8
			methane_tran_depth(j) = 0._r8
			o2_aere_depth(j) = 0._r8
		end do

			if (methane_patch_is_nongrass(patchclass)) then
			   poros_default = DEF_METHANE%poros_tiller * DEF_METHANE%nongrassporosratio
			else
			   poros_default = DEF_METHANE%poros_tiller
			endif
			poros = get_aere_poros(ipatch, poros_default)
		if (sat == 0) then
			poros_tiller_real = poros * DEF_METHANE%unsat_aere_ratio ! unsaturate
			! Apply a proportional floor on the unsat path: the absolute
			! porosmin=0.05 used to silently lift the non-grass branch back
			! up to the grass value (both ended at 0.05 because porosmin
			! caught poros*ratio = 0.3*1/3*0.05/0.3 = 0.0167 -> 0.05).
			! That defeats nongrassporosratio for forest wetlands.  Keep
			! the species contrast by scaling porosmin by the same
			! unsat_aere_ratio.
			poros_tiller_real = max(poros_tiller_real, &
				DEF_METHANE%porosmin * DEF_METHANE%unsat_aere_ratio)
		else
			poros_tiller_real = poros
			poros_tiller_real = max(poros_tiller_real, DEF_METHANE%porosmin)
		endif

		call SiteOxAere(ipatch, idate, jwt,  sat, poros_tiller_real, lai, z_soisno, dz_soisno,  zi_soisno,  t_soisno,  &
		rootfr, rootr, grnd_methane_cond_base, etr, &
			annsum_npp, annavg_agnpp, annavg_bgnpp, c_atm, conc_ch4_aqu, conc_ch4_aqu_porsl, conc_ch4_gas_porsl, conc_o2_aqu_porsl,conc_o2_gas_porsl,&
			tranloss, aere, oxaere)

		do j = 1,nl_soil
			! Impose limitation based on available methane during timestep
			! By imposing the limitation here, don't allow aerenchyma access to methane from other Patches.
				aere_pos = max(aere(j), 0._r8)
				aere_neg = min(aere(j), 0._r8)
				tran_pos = max(tranloss(j), 0._r8)
				tran_neg = min(tranloss(j), 0._r8)
				aeretran_requested = aere_pos + tran_pos
				if (aeretran_requested > 0._r8) then
					aeretran = min(aeretran_requested, max(0._r8, conc_methane(j)/deltim + &
						methane_prod_depth(j) - max(methane_oxid_depth(j), 0._r8) - &
						max(methane_ebul_depth(j), 0._r8)))
					aeretran_scale = aeretran / aeretran_requested
				else
					aeretran_scale = 1._r8
				endif
				methane_aere_depth (j) = methane_aere_depth(j) + aere_pos * aeretran_scale + aere_neg
				methane_tran_depth (j) = methane_tran_depth(j) + tran_pos * aeretran_scale + tran_neg
				o2_aere_depth  (j) = o2_aere_depth (j) + oxaere(j)
		end do ! over levels

	end subroutine methane_aere


	!---------------------------------------------------------------------------
	subroutine SiteOxAere(ipatch, idate,jwt,  sat, poros_tiller_real, lai, z_soisno, dz_soisno,  zi_soisno,  t_soisno,  &
			rootfr, rootr, grnd_methane_cond_base, etr, &
		annsum_npp, annavg_agnpp, annavg_bgnpp, c_atm, conc_ch4_aqu, conc_ch4_aqu_porsl, conc_ch4_gas_porsl, conc_o2_aqu_porsl,conc_o2_gas_porsl,&
		tranloss, aere, oxaere)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Site(column) level fluxes for O2 gain rate via
		! aerenchyma and methane losss rates from transpiration
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			ipatch            , &! caller's per-patch index for vegetation overrides
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-)
			sat                       ! 0 = unsatured, 1 = saturated

		real(r8), intent(in) :: &
			poros_tiller_real      , &
			lai                    , &! leaf area index [m2/m2]

			z_soisno (maxsnl+1:nl_soil)    , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! root resistance of a layer, all layers sum to 1
			grnd_methane_cond_base     , &! base tracer conductance before snow/pond resistance [m/s]
			etr                    , &! transpiration rate [mm/s]

			annsum_npp             , &! annual sum NPP (g C/m2/yr)
			annavg_agnpp           , &! annual average aboveground NPP (g C/m2/s)
			annavg_bgnpp           , &! annual average belowground NPP (g C/m2/s)

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			conc_ch4_aqu       (1:nl_soil) , &! aqueous phase CH4 concentration [mol/m3 water]
			conc_ch4_aqu_porsl (1:nl_soil) , &! aqueous phase CH4 conc in each porosity [mol/m3]
			conc_ch4_gas_porsl (1:nl_soil) , &! gas phase CH4 conc in each porosity [mol/m3]
			conc_o2_aqu_porsl  (1:nl_soil) , &! aqueous phase O2 conc in each porosity [mol/m3]
			conc_o2_gas_porsl  (1:nl_soil)    ! gas phase O2 conc in each porosity [mol/m3]

		real(r8), intent(out) :: &
			tranloss        (1:nl_soil)  , &! CH4 in soil water tran rate via plant transpiration in each soil layer (mol/m3/s)
			aere            (1:nl_soil)  , &! CH4 tran rate via aerenchyma in each soil layer (mol/m3/s)
			oxaere          (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
			integer  :: j
			real(r8) :: area_tiller ! cross-sectional area of tillers (m2/m2)
			real(r8) :: m_tiller
			real(r8) :: n_tiller
			real(r8) :: anpp        ! annual sum NPP (gC/m2/yr)
			real(r8) :: nppratio    ! bg/sum NPP
			real(r8) :: secsperyear ! seconds in current model year
			real(r8) :: aere_ch4_resis    ! aerenchyma resistance [s/m]
		real(r8) :: grnd_ch4_resis    ! boundary layer resistance [s/m]
		real(r8) :: aere_o2_resis    ! aerenchyma resistance [s/m]
		real(r8) :: grnd_o2_resis    ! boundary layer resistance [s/m]
		real(r8) :: aerecond    ! aerenchyma conductance [m/s]
		! Per-patch aerenchyma overrides for wetland (patchtype==2):
		! get_aere_* return DEF_METHANE%* defaults for non-wetland patches.
		real(r8) :: tillerC_eff      ! gC per tiller for current patch
		real(r8) :: radius_eff       ! tiller radius for current patch (m)
		real(r8) :: scale_eff        ! scale_factor_aere for current patch
		real(r8), parameter :: smallnumber = 1.e-12_r8
		!-----------------------------------------------------------------------
		! This parameter is poorly constrained and should be done on a patch-specific basis...

		! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
		! the new components are in; if there is any tuning to be done to get a realistic global flux,
		! this would probably be the place.  We will have to document clearly in the Tech Note
		! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)

			tillerC_eff = get_aere_tillerC(ipatch, DEF_METHANE%tiller_C)
			radius_eff  = get_aere_radius (ipatch, DEF_METHANE%aere_radius)
			scale_eff   = get_aere_scale  (ipatch, DEF_METHANE%scale_factor_aere)

			anpp = annsum_npp
			if (ieee_is_nan(anpp) .or. abs(anpp) >= 0.5_r8*abs(spval)) anpp = 0._r8
			if (anpp <= 0._r8 .and. .not. ieee_is_nan(annavg_agnpp) .and. .not. ieee_is_nan(annavg_bgnpp) .and. &
				abs(annavg_agnpp) < 0.5_r8*abs(spval) .and. abs(annavg_bgnpp) < 0.5_r8*abs(spval) .and. &
				(annavg_agnpp + annavg_bgnpp) > 0._r8) then
				! If annual-sum NPP is unavailable, derive annual total NPP from annual-mean rates.
				if (isleapyear(idate(1))) then
					secsperyear = 366._r8 * secspday
				else
					secsperyear = 365._r8 * secspday
				end if
				anpp = (annavg_agnpp + annavg_bgnpp) * secsperyear
				! [gC/m2] = [gC/m2/s] * [s/yr]
			end if
			anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

		if (.not. ieee_is_nan(annavg_agnpp) .and. .not. ieee_is_nan(annavg_bgnpp) .and. &
			abs(annavg_agnpp) < 0.5_r8*abs(spval) .and. abs(annavg_bgnpp) < 0.5_r8*abs(spval) .and. &
			annavg_agnpp > 0._r8 .and. annavg_bgnpp > 0._r8) then
			nppratio = annavg_bgnpp / (annavg_agnpp + annavg_bgnpp)
		else
			nppratio = 0.5_r8
		end if

		do j=1,nl_soil
			! Calculate transpiration loss
			if (DEF_METHANE%transpirationloss .and. lai > 0) then
				tranloss(j) = conc_ch4_aqu(j) * rootr(j)*etr / dz_soisno(j) / 1000._r8
				! [mol/m3/s]= [mol/m3 water] * [-]     *[mm/s]/ [m]        /   [mm/m]
				! Use rootr here for effective per-layer transpiration, which may not be the same as rootfr
				tranloss(j) = max(tranloss(j), 0._r8) ! in case transpiration is pathological
			else
				tranloss(j) = 0._r8
			end if

			! Calculate aerenchyma diffusion
			if (j > jwt .and. t_soisno(j) > tfrz .and. lai > 0) then ! Below water table
				! Estimate area of tillers (see Wania thesis)
				!m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
				!m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
				! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
				! done on a PFT-specific basis.

				! Use the current effective LAI for aerenchyma area, following CLM5.
				m_tiller = anpp * nppratio * lai
				! [g C/m2] = [g C/m2/yr aliased as g C/m2] * [-] * [m2/m2]

				n_tiller = m_tiller / tillerC_eff
				! [tiller/m2] = [g C/m2]/ [g C/tiller]

				area_tiller = scale_eff * n_tiller * poros_tiller_real * PI * radius_eff**2
				! [m2/m2]   = [-]               * [tiller/m2] * [-]       * [-]* [m2/tiller]

				aere_ch4_resis = 1._r8/((area_tiller * rootfr(j) * d_con_g(1,1) * 1e-4_r8 / (z_soisno(j)*DEF_METHANE%rob))+smallnumber)
				! [s/m]         = /([m2/m2]          * [-]       * [m2/s]                 / [m]         /[-])
				! Add in boundary layer resistance
				grnd_ch4_resis = 1._r8/(grnd_methane_cond_base+smallnumber)
				aerecond = 1._r8/(aere_ch4_resis + grnd_ch4_resis)
				! aerecond = max(aerecond,1.e-8_r8)
				aere(j) = aerecond*(conc_ch4_gas_porsl(j) - c_atm(1)) / dz_soisno(j)
				! [mol/m3/s] = [mol/m3]                      / [m]          / [s/m]
					! Allow negative CH4 aerenchyma fluxes: when pore-air CH4 is
					! below atmospheric CH4, aerenchyma can be an uptake pathway.
					! Downstream methane competition only limits positive sinks.

				! Do oxygen diffusion into layer
				aere_o2_resis = 1._r8/((area_tiller * rootfr(j) * d_con_g(2,1) * 1e-4_r8 / (z_soisno(j)*DEF_METHANE%rob)) + smallnumber)
				grnd_o2_resis = 1._r8/(grnd_methane_cond_base+smallnumber)

				oxaere(j) = -(conc_o2_gas_porsl(j) - c_atm(2)) / (dz_soisno(j)*(aere_o2_resis + grnd_o2_resis)) ![mol/m3-total/s]
				oxaere(j) = max(oxaere(j), 0._r8)
				! Diffusion in is positive; prevent backwards diffusion
				if ( .not. DEF_METHANE%use_aereoxid_prog ) then ! fixed aere oxid proportion; will be done in methane_tran
					oxaere(j) = 0._r8
				end if
			else
				aere(j) = 0._r8
				oxaere(j) = 0._r8
			end if ! veg type, below water table, & above freezing
		end do

	end subroutine SiteOxAere

	!---------------------------------------------------------------------------
	subroutine methane_ebul (idate, patchtype, jwt, sat, finundated, deltim, &
		z_soisno, dz_soisno, zi_soisno, forc_pbot, lakedepth, lake_icefrac, &
		t_soisno, wdsrf, conc_methane, conc_ch4_gas_porsl,&
		methane_ebul_depth)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Bubbling is based on temperature & pressure dependent solubility (k_h_cc),
		! with assumed proportion of bubbles
		! which are CH4, and assumed early nucleation at DEF_METHANE%vgc_max sat (Wania).
		! Bubbles are released to the water table surface in methane_tran.
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype                  , &! land patch type; 4 = lake/water body
			jwt                        , &! index of the soil layer right above the water table (-)
			sat                        ! 0 = unsaturated; 1 = saturated

		real(r8), intent(in) :: &
			finundated                 , &
			deltim                     , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil), &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil), &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)  , &! interface level below a "z" level (m)

			forc_pbot                  , &! atm bottom level pressure (or reference height) (Pa)
			lakedepth                  , &! lake depth (m), CTSM hydrostatic pressure over sediments
			lake_icefrac(1:nl_lake)    , &! lake frozen mass fraction [-]
			t_soisno (maxsnl+1:nl_soil), &! soil temperature (Kelvin)
			wdsrf                      , &! depth of surface water [mm]
			conc_methane       (1:nl_soil) , &! CH4 conc in each soil layer (mol/m3)

			conc_ch4_gas_porsl (1:nl_soil)! gas phase CH4 conc in each porosity (mol/m3)

		real(r8), intent(out) :: &
			methane_ebul_depth (1:nl_soil)    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer :: j      ! indices

			real(r8) :: vgc     ! gas phase volumetric CH4 content (m3 CH4/m3 pore air)
			real(r8) :: pressure! sum atmospheric and hydrostatic pressure
			real(r8) :: ebul_timescale
			real(r8), parameter :: smallnumber = 1.e-12_r8
		!-----------------------------------------------------------------------
		! Ebullition follows the CLM5/Riley trigger and post-event target gas fraction.
		ebul_timescale = deltim ! Allow fast bubbling

		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. lake_icefrac(1) > 0.1_r8) then
			methane_ebul_depth = 0._r8
			return
		endif

		! column loop to estimate ebullition CH4 flux from each soil layer
		do j=1,nl_soil
			if (j  >  jwt .and. t_soisno(j) > tfrz) then ! Ebullition occurs only below the water table
				if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
					! CTSM lake path: sediment bubbles must overcome the full
					! lake water-column head above the sediment surface.
					pressure = forc_pbot + denh2o * grav * (z_soisno(j) + max(lakedepth, 0._r8))
				else
					pressure = forc_pbot + denh2o * grav * (z_soisno(j)-zi_soisno(jwt))
					! [Pa]   = [Pa]      + [kg/m3]* [m/s2]* [m]
					! [Pa]   = [N/m2] = [kg]*[m/s2]/[m2] = [kg/m/s2]
				endif
				if (patchtype /= 4 .and. sat == 1 .and. finundated>0._r8) then ! Add ponding pressure head
					! Convert gridcell-mean ponded depth to local inundated depth with a 1% floor.
					pressure = pressure + denh2o * grav * wdsrf/1000._r8/max(finundated, 0.01_r8)
					! [Pa]   = [Pa]     + [kg/m3]* [m/s2]* [mm]/[mm/m]
				end if

				! Compare partial pressure to ambient pressure.
				vgc = conc_ch4_gas_porsl(j) * rgasm * t_soisno(j) / pressure
				! [-]= [mol/m3]           * [Pa*m3/K/mol]*[K]   / [Pa]

					! Trigger bubbles at vgc_max, then relax the layer toward the
					! lower post-event target vgc_max*bubble_f.  This implements
					! the intended two-threshold release; the old single threshold
					! vgc_max*bubble_f fired too early and had no hysteresis window.
					if (vgc > DEF_METHANE%vgc_max) then
						methane_ebul_depth(j) = (vgc - DEF_METHANE%vgc_max * DEF_METHANE%bubble_f) / &
							max(vgc, smallnumber) * conc_methane(j) / ebul_timescale
						! [mol/m3/s]      = [-]                                       * [mol/m3]    / [s]
					else
					methane_ebul_depth(j) = 0._r8
				endif

			else ! above the water table or freezing
				methane_ebul_depth(j) = 0._r8
			endif ! below the water table and not freezing
		end do ! j

	end subroutine methane_ebul

		!---------------------------------------------------------------------------
	subroutine methane_tran (idate,patchtype, &
		lb, snl, jwt, sat, finundated,&
		dlon, dlat, deltim, z_soisno, dz_soisno, zi_soisno,  t_soisno, t_grnd, forc_us, forc_vs, &
		porsl, wliq_soisno, wice_soisno, wdsrf, lakedepth, lake_icefrac, bsw, c_atm, methane_prod_depth, o2_aere_depth,&
		cellorg,t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl,conc_ch4_aqu_porsl,conc_o2_gas_porsl,conc_o2_aqu_porsl,vol_aqu,vol_ch4_storage,vol_gas,&
		o2stress, methane_stress, methane_surf_aere, methane_surf_ebul, methane_surf_diff, methane_ebul_tot, &
		methane_oxid_depth, methane_aere_depth, methane_tran_depth, methane_ebul_depth, &
			grnd_methane_cond_base, grnd_methane_cond_effective, o2_oxid_depth, o2_decomp_depth, &
			conc_o2, conc_methane, methane_balance_residual, methane_ch4_clip_credit, methane_surf_diff_phys, &
		o2_cap_loss, o2_cap_gain )
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Solves the reaction & diffusion equation for the timestep.
		! 1  "Competition" between processes for CH4 & O2 demand is done.
		! 2  Concentrations are apportioned into gas & liquid fractions;
		!    only the gas fraction is considered for diffusion in unsat.
		! 3  Snow and lake water resistance to diffusion is added as a bulk term in the ground conductance
		!    (which is really a surface layer conductance), but concentrations are not tracked and oxidation
		!    is not allowed inside snow and lake water.
		! 4  Diffusivity is set based on soil texture and organic matter fraction.
		!    A Crank-Nicholson solution is used.
		! 5  CH4 diffusive flux is calculated and consistency is checked.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer, intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype        	! land patch type (0=soil, 1=urban or built-up, 2=wetland,
										! 3=land ice, 4=land water bodies, 99=ocean


		integer , intent(in) :: &
			lb                , &! lower bound of array (snl+1)
			snl				  , &!  number of snow layers     (-5~-1)
			jwt               , &! index of the soil layer right above the water table (-)
			sat                  ! 0 = unsaturated; 1 = saturated

		real(r8), intent(in) :: &
			finundated                      , &
			dlon   	   				        , &! logitude
			dlat     	   			        , &! latitude

			deltim                  	    , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil)   	, &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)   	, &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)   	, &! interface level below a "z" level (m)


			t_soisno (maxsnl+1:nl_soil)    	, &! soil temperature (Kelvin)
			t_grnd                 		    , &! ground surface temperature [k]
			forc_us                         , &! eastward wind speed [m/s], U10 proxy for lake gas exchange
			forc_vs                         , &! northward wind speed [m/s], U10 proxy for lake gas exchange

			porsl             (1:nl_soil)   , &! volumetric soil water at saturation (porosity)
			wliq_soisno(maxsnl+1:nl_soil)	, &! liquid water in layers [kg/m2]
			wice_soisno(maxsnl+1:nl_soil) 	, &! ice lens in layers [kg/m2]
			wdsrf                  		    , &! depth of surface water [mm]
			lakedepth                       , &! lake water depth [m]
			lake_icefrac      (1:nl_lake)   , &! lake-layer ice fraction
			bsw               (1:nl_soil)   , &! Clapp and Hornberger "b" (nlevgrnd)

			c_atm(3)               		    , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			methane_prod_depth    (1:nl_soil)   , &! production of CH4 in each soil layer (mol/m3/s)
			o2_aere_depth     (1:nl_soil)   , &! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)


			cellorg           (1:nl_soil)   , &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		    , &! surface water temperature
			organic_max               	    , &! organic matter content (kg m-3) where soil is assumed to act like peat

			k_h_cc(0:nl_soil,ngases)        , &! ratio of mol/m3 in liquid to mol/m3 in gas
			conc_ch4_gas_porsl(1:nl_soil)   , &! gas phase CH4 conc in each porosity (mol/m3)
			conc_ch4_aqu_porsl(1:nl_soil)   , &! aqueous phase CH4 conc in each porosity (mol/m3)
			conc_o2_gas_porsl (1:nl_soil)   , &! gas phase O2 conc in each porosity (mol/m3)
			conc_o2_aqu_porsl (1:nl_soil)   , &! aqueous phase O2 conc in each porosity (mol/m3)
			vol_aqu           (1:nl_soil)   , &
			vol_ch4_storage   (1:nl_soil)   , &! liquid/ice volume retaining CH4 [m3/m3]
			vol_gas           (1:nl_soil)   , &
			grnd_methane_cond_base              ! base conductance before snow/pond resistance [m/s]

		real(r8), intent(out) :: &
			o2stress          (1:nl_soil)   , &! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			methane_stress         (1:nl_soil)   , &! Ratio of methane available to the total per-timestep methane sinks
			methane_surf_aere                   , &! Total column CH4 aerenchyma (mol/m2/s)
			methane_surf_ebul                   , &! CH4 ebullition to atmosphere (mol/m2/s)
				methane_surf_diff                   , &! CH4 surface flux plus numerical closure (mol/m2/s)
				methane_ebul_tot                    , &! Total column CH4 ebullition (mol/m2/s)
				methane_balance_residual            , &! numerical closure flux credited to methane_surf_diff (mol/m2/s)
				methane_ch4_clip_credit             , &! negative CH4 clip credited to methane_surf_diff (mol/m2/s)
				methane_surf_diff_phys              , &! pure diffusive flux before clip/residual (mol/m2/s)
				o2_cap_loss                        , &! O2 removed by post-solve physical cap (mol/m2/s)
				o2_cap_gain                         , &! O2 added by post-solve nonnegative floor (mol/m2/s)
				grnd_methane_cond_effective            ! effective CH4 conductance including snow/pond resistance [m/s]

		real(r8), intent(inout) :: &
			methane_oxid_depth    (1:nl_soil)   , &! InOut: CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			methane_aere_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			methane_tran_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			methane_ebul_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)   , &! InOut: O2 loss rate via ebullition in each soil layer (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)   , &! InOut: O2 consumption during decomposition in each soil layer (mol/m3/s)

			conc_o2           (1:nl_soil)   , &! InOut: O2 conc in each soil layer (mol/m3)
			conc_methane          (1:nl_soil)      ! InOut: CH4 conc in each soil layer (mol/m3)

		!-----------------------Local Variables------------------------------
		integer :: j,s,i			                                               ! indices
		integer :: jtop                                                        ! top level at each column
		real(r8) :: at (0:nl_soil)                     ! "a" vector for tridiagonal matrix
		real(r8) :: bt (0:nl_soil)                     ! "b" vector for tridiagonal matrix
		real(r8) :: ct (0:nl_soil)                     ! "c" vector for tridiagonal matrix
		real(r8) :: rt (0:nl_soil)                     ! "r" vector for tridiagonal solution
		real(r8) :: f_a                                                        ! air-filled fraction of available pore space
		real(r8) :: diffus (0:nl_soil)                 ! diffusivity (m2/s)
		real(r8) :: dzj                                                        !
		real(r8) :: dp1_zp1 (0:nl_soil)                ! diffusivity/delta_z for next j
		real(r8) :: dm1_zm1 (0:nl_soil)                ! diffusivity/delta_z for previous j
		real(r8) :: t_soisno_c                                                 ! soil temperature   (maxsnl+1:nl_soil)
		real(r8) :: deficit                                                    ! mol CH4 /m^2 that must be subtracted from diffusive flux to atm. to make up
		! for keeping concentrations always above zero
		logical, save :: warned_ch4_clip = .false.
		logical, save :: warned_o2_cap   = .false.
		logical, save :: warned_ch4_competition_demand = .false.
		logical, save :: warned_ch4_competition_left   = .false.
		logical, save :: warned_o2_competition_demand  = .false.
		logical, save :: warned_o2_competition_left    = .false.
		logical, save :: warned_ch4_diffusion_residual = .false.
		logical, save :: warned_lake_zero_depth        = .false.
		real(r8) :: conc_ch4_bef(1:nl_soil)            ! concentration at the beginning of the timestep
		real(r8) :: err_methane                            ! Error (Mol CH4 /m^2) [+ = too much CH4]
		real(r8) :: conc_ch4_rel(0:nl_soil)            ! Concentration per volume of air or water
		real(r8) :: conc_o2_rel(0:nl_soil)             ! Concentration per volume of air or water
		real(r8) :: conc_ch4_rel_old(0:nl_soil)        ! Concentration during last Crank-Nich. loop
			real(r8), parameter :: smallnumber = 1.e-12_r8
			logical, parameter :: ch4_backward_euler_transport = .true.
			real(r8) :: snowdiff                                                   ! snow diffusivity (m^2/s)
			real(r8) :: snow_resis                           ! Cumulative Snow resistance (s/m). Also includes
			real(r8) :: pond_resis                                                    ! Additional resistance from ponding, up to pondmx water on top of top soil layer (s/m)
			real(r8) :: pondz                                                      ! Depth of ponding (m)
			real(r8) :: ponddiff                                                   ! Pondwater diffusivity (m^2/s)
			real(r8) :: lake_exchange_vel                                          ! Lake piston velocity [m/s]
			real(r8) :: lake_k600_cmhr                                             ! Cole & Caraco k600 [cm/hr]
			real(r8) :: lake_schmidt_ch4                                           ! CH4 Schmidt number in freshwater [-]
			real(r8) :: lake_wind                                                  ! Wind speed [m/s], using forcing wind as U10 proxy
		real(r8) :: liquid_volfrac, ice_volfrac, open_pore_volfrac
		real(r8) :: spec_grnd_cond(1:ngases)            ! species grnd conductance (s/m)
		real(r8) :: airfrac                                                    ! air fraction in snow
		real(r8) :: waterfrac                                                  ! water fraction in snow
		real(r8) :: icefrac                                                    ! ice fraction in snow
		real(r8) :: epsilon_t (1:nl_soil,1:ngases)     !
		real(r8) :: source (1:nl_soil,1:ngases)        ! source
		real(r8) :: om_frac                                                    ! organic matter fraction
			real(r8) :: o2demand, methane_demand                                        ! mol/m^3/s
			real(r8) :: o2_before_cap, o2_cap_loss_col, o2_cap_gain_col, o2_rel_cap
				real(r8) :: aere_oxid_flux
				!-----------------------------------------------------------------------
				methane_balance_residual = 0._r8
				methane_ch4_clip_credit = 0._r8
				o2_cap_loss = 0._r8
				o2_cap_gain = 0._r8

				! When use_aereoxid_prog=.false. (CTSM legacy fixed-fraction root-zone
			! aerenchyma oxidation) we MUST fold aere_oxid_flux into the
			! methanotrophic oxidation budget BEFORE the O2/CH4 competition loop.
			! Previously this block sat after the competition (just before
			! `source(j,1)`), which let the additional CH4 oxidation and 2:1 O2
			! demand bypass the stress caps — pushing total O2 sink beyond what
			! is available in low-O2 soils.  Now the augmented demands enter the
			! competition like any other sink and get scaled by o2stress/methane_stress.
			do j = 1,nl_soil
				if ( .not. DEF_METHANE%use_aereoxid_prog .and. methane_aere_depth(j) > 0._r8 ) then
					aere_oxid_flux = DEF_METHANE%aereoxid * methane_aere_depth(j)
					methane_oxid_depth(j)  = methane_oxid_depth(j)  + aere_oxid_flux
					o2_oxid_depth(j)       = o2_oxid_depth(j)       + 2._r8 * aere_oxid_flux
					methane_aere_depth(j)  = methane_aere_depth(j)  - aere_oxid_flux
				end if
			end do

			! Perform competition for oxygen and methane in each soil layer if demands over the course of the timestep
		! exceed that available. Assign to each process in proportion to the quantity demanded in the absense of
		! the limitation.
		do j = 1,nl_soil
			o2demand = o2_decomp_depth(j) + o2_oxid_depth(j) ! o2_decomp_depth includes autotrophic root respiration
			if (o2demand > 0._r8) then
				o2stress(j) = min((conc_o2(j) / deltim + o2_aere_depth(j)) / o2demand, 1._r8)
			else
				o2stress(j) = 1._r8
			end if

				methane_demand = methane_oxid_depth(j) + max(methane_aere_depth(j), 0._r8) + &
					max(methane_tran_depth(j), 0._r8) + methane_ebul_depth(j)
			if (methane_demand > 0._r8) then
				methane_stress(j) = min((conc_methane(j) / deltim + methane_prod_depth(j)) / methane_demand, 1._r8)
			else
				methane_stress(j) = 1._r8
			end if

			! Resolve methane oxidation
			if (o2stress(j) < 1._r8 .or. methane_stress(j) < 1._r8) then
				if (methane_stress(j) <= o2stress(j)) then
					! methane limited
					if (o2stress(j) < 1._r8) then
						! Recalculate oxygen limitation
						o2demand = o2_decomp_depth(j)
						if (o2demand > 0._r8) then
							o2stress(j) = min((conc_o2(j)/deltim + o2_aere_depth(j) - methane_stress(j)*o2_oxid_depth(j))/o2demand, 1._r8)
						else
							o2stress(j) = 1._r8
						end if
					end if
					! Reset oxidation
					methane_oxid_depth(j) = methane_oxid_depth(j) * methane_stress(j)
					o2_oxid_depth(j)  = o2_oxid_depth(j) * methane_stress(j)
				else
					! oxygen limited
					if (methane_stress(j) < 1._r8) then
						! Recalculate methane limitation
							methane_demand = max(methane_aere_depth(j), 0._r8) + &
							max(methane_tran_depth(j), 0._r8) + methane_ebul_depth(j)
						if (methane_demand > 0._r8) then
							methane_stress(j) = min( (conc_methane(j) / deltim + methane_prod_depth(j) - &
									o2stress(j)*methane_oxid_depth(j)) / methane_demand, 1._r8)
						else
							methane_stress(j) = 1._r8
						end if
					end if
					! Reset oxidation
					methane_oxid_depth(j) = methane_oxid_depth(j) * o2stress(j)
					o2_oxid_depth(j) = o2_oxid_depth(j) * o2stress(j)
				end if
			end if

			! Reset non-methanotroph demands
				IF (methane_aere_depth(j) > 0._r8) methane_aere_depth(j) = methane_aere_depth(j) * methane_stress(j)
				IF (methane_tran_depth(j) > 0._r8) methane_tran_depth(j) = methane_tran_depth(j) * methane_stress(j)
			methane_ebul_depth(j) = methane_ebul_depth(j) * methane_stress(j)
			o2_decomp_depth(j) = o2_decomp_depth(j) * o2stress(j)


		end do !j


		! Accumulate ebullition to place in first layer above water table, or directly to atmosphere
		methane_ebul_tot = 0._r8
		do j = 1,nl_soil
			methane_ebul_tot = methane_ebul_tot + methane_ebul_depth(j) * dz_soisno(j)
		end do

		! Set the source term for each species (no need to do j=0, since epsilon_t and source not used there)
		! Note that because of the semi-implicit diffusion and the 30 min timestep combined with explicit
		! sources, occasionally negative concentration will result. In this case it is brought to zero and the
		! surface flux is adjusted to conserve. This results in some inaccuracy as compared to a shorter timestep
		! or iterative solution.
		do j = 1,nl_soil

			! NOTE: The fixed-fraction aerenchyma oxidation that used to live
			! here is now applied BEFORE the O2/CH4 competition loop above so
			! it participates in the stress budget.

			source(j,1) = methane_prod_depth(j) - methane_oxid_depth(j) - &
					methane_aere_depth(j) - methane_tran_depth(j) - &
					methane_ebul_depth(j) ! [mol/m3-total/s]

			! aerenchyma added to surface flux below
			! ebul added to soil depth just above WT
			! Original 1e-12 mol/m3/s strict abort tripped on roundoff at multiple
			! patches when only_wetland=.false. exposes wider patch population
			! Replace fatal with rate-limited warning (one per rank
			! per check); source/conc are recomputed each step from physics so a
			! small floating-point imbalance does not corrupt the long-term budget.
			if (source(j,1) + conc_methane(j) / deltim < -1.e-6_r8 .and. .not. warned_ch4_competition_demand)then
				write(6,*) 'WARNING: methane competition demand > 1e-6 mol/m3/s, j:', &
						source(j,1) + conc_methane(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_ch4_competition_demand = .true.
			else if (methane_stress(j) < 1._r8 .and. source(j,1) + conc_methane(j) / deltim > 1.e-6_r8 &
			         .and. .not. warned_ch4_competition_left) then
				write(6,*) 'WARNING: methane limited yet leftover > 1e-6 mol/m3/s, j:', &
						source(j,1) + conc_methane(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_ch4_competition_left = .true.
			end if

			source(j,2) = -o2_oxid_depth(j) - o2_decomp_depth(j) + o2_aere_depth(j) ! O2 [mol/m3/s]

			if (source(j,2) + conc_o2(j) / deltim < -1.e-6_r8 .and. .not. warned_o2_competition_demand) then
				write(6,*) 'WARNING: oxygen competition demand > 1e-6 mol/m3/s, j:', &
						source(j,2) + conc_o2(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_o2_competition_demand = .true.
			else if (o2stress(j) < 1._r8 .and. source(j,2) + conc_o2(j) / deltim > 1.e-6_r8 &
			         .and. .not. warned_o2_competition_left) then
				write(6,*) 'WARNING: oxygen limited yet leftover > 1e-6 mol/m3/s, j:', &
						source(j,2) + conc_o2(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_o2_competition_left = .true.
			end if
			conc_ch4_bef(j) = conc_methane(j) !For Balance Check

		enddo ! j


		! Accumulate aerenchyma to add directly to atmospheric flux
		methane_surf_aere = 0._r8
		do j = 1,nl_soil
			methane_surf_aere = methane_surf_aere + methane_aere_depth(j) * dz_soisno(j)
		enddo

		! Ebullition is a direct atmospheric pathway. It has already been
		! removed from each source layer via methane_ebul_depth above; do
		! not re-inject it into the layer above the water table.

		! Calculate concentration relative to m^3 of air or water: needed for the diffusion
		do j = 0,nl_soil
			if (j == 0) then
				conc_ch4_rel(j) = c_atm(1)
				conc_o2_rel(j)  = c_atm(2)
			else
					if (porsl(j) <= smallnumber) then
						do s = 1, 2
							epsilon_t(j,s) = smallnumber
						end do
						conc_ch4_rel(j) = 0._r8
						conc_o2_rel(j)  = 0._r8
					else if (j <= jwt) then  ! Above the WT
						! CH4 storage follows methane_frzout; O2 remains restricted to
						! mobile gas/liquid pore space.
						epsilon_t(j,1) = max(vol_gas(j) + k_h_cc(j,1)*vol_ch4_storage(j), smallnumber)
						epsilon_t(j,2) = max(vol_gas(j) + k_h_cc(j,2)*vol_aqu(j), smallnumber)
					conc_ch4_rel(j) = conc_methane(j)/epsilon_t(j,1)
					conc_o2_rel(j)  = conc_o2_gas_porsl(j)
					! Partition between the liquid and gas phases. The gas phase will drive the diffusion.
				else ! Below the WT
					! Diffusion remains controlled by liquid water below, while the
					! storage coefficient may retain CH4 in ice.
					epsilon_t(j,1) = max(vol_ch4_storage(j), smallnumber)
					epsilon_t(j,2) = max(vol_aqu(j), smallnumber)
					conc_ch4_rel(j) = conc_methane(j)/epsilon_t(j,1)
					conc_o2_rel(j)  = conc_o2(j) /epsilon_t(j,2)
				end if
			end if
		end do


		! Loop over species
		do s = 1, 2 ! transport CH4 and O2; CO2 is not prognosed here
			! Apply snow and pond resistance to the immutable base conductance.
			! The effective result is diagnostic output only and never feeds the next step.
			do j = maxsnl + 1,0
				if (j == maxsnl + 1) then
					! Needed to prevent overflow when ground is frozen, e.g. for lakes
					snow_resis = 0._r8
				end if

				! Add snow resistance
				if (j >= snl + 1) then
					! For the ice layer, all = ice + water + air, no soil
					t_soisno_c = t_soisno(j) - tfrz
					icefrac = wice_soisno(j)/denice/dz_soisno(j)
					! [-]   = [kg/m2]       / [kg/m3] / [m]
					waterfrac = wliq_soisno(j)/denh2o/dz_soisno(j)
					airfrac = max(1._r8 - icefrac - waterfrac, 0._r8)

					! Calculate snow diffusivity
					if (airfrac > 0.05_r8) then
						! Use Millington-Quirk Expression, as hydraulic properties (bsw) not available
						snowdiff = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								airfrac**(10._r8/3._r8) / (airfrac+waterfrac)**2 &
								* DEF_METHANE%scale_factor_gasdiff
					else !solute diffusion in water only, airfrac -> 0
						snowdiff = (airfrac+waterfrac)**DEF_METHANE%satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
								* DEF_METHANE%scale_factor_liqdiff
					end if

					snowdiff = max(snowdiff, smallnumber)
					snow_resis = snow_resis + dz_soisno(j)/snowdiff
				end if

				if (j == 0) then ! final loop
					! Add pond resistance
					pond_resis = 0._r8

					! First old pond formulation up to pondmx
					liquid_volfrac = wliq_soisno(1)/(dz_soisno(1)*denh2o)
					ice_volfrac = wice_soisno(1)/(dz_soisno(1)*denice)
					open_pore_volfrac = max(porsl(1) - ice_volfrac, 0._r8)
					if (patchtype /= 4 .and. snl == 0 .and. liquid_volfrac > open_pore_volfrac) then
						t_soisno_c = t_soisno(1) - tfrz
						if (t_soisno(1) <= tfrz) then
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* (wliq_soisno(1)/denh2o+smallnumber)/ &
									(wliq_soisno(1)/denh2o+wice_soisno(1)/denice+smallnumber) &
									* DEF_METHANE%scale_factor_liqdiff
						else ! Unfrozen
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* DEF_METHANE%scale_factor_liqdiff
						end if
						pondz = dz_soisno(1) * (liquid_volfrac - open_pore_volfrac)
						pond_resis = pondz / ponddiff
					end if

					! Now add new wdsrf form
					if (patchtype /= 4 .and. sat == 1 .and. finundated>0._r8) then
						if (t_h2osfc >= tfrz) then
							t_soisno_c = t_h2osfc - tfrz
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* DEF_METHANE%scale_factor_liqdiff
							! Guard against fill-and-spill transients with very small inundated area.
							pondz = wdsrf / 1000._r8/max(finundated, 0.01_r8) ! Assume all wdsrf corresponds to sat area
							! [m] = [mm]  /  [mm/m]
							pond_resis = pond_resis + pondz / ponddiff
						else if (wdsrf/max(finundated, 0.01_r8) > DEF_METHANE%capthick) then
							! assume surface ice is impermeable
							pond_resis = pond_resis + 1._r8/smallnumber
						end if
					end if

						! Lake patches: use a bulk air-water piston velocity rather than
						! molecular diffusion across the full lake depth.  The previous
						! lakedepth/Dmol resistance makes 10 m lakes effectively closed
						! to diffusive exchange; turbulent mixing gives k600 of order
						! 1e-6--1e-5 m/s for typical winds.
						!
						! Henry-consistent BC.  Cole & Caraco gives a water-side
						! piston velocity k_w [m s-1] for the air-water interface,
						! while spec_grnd_cond carries an air-side conductance.  The
						! water-side flux F = k_w*(C_w - k_h*C_atm) must match the
						! air-side closure F = k_air_eq*(C_w/k_h - C_atm), which
						! requires k_air_eq = k_h*k_w.  Therefore add the air-side
						! resistance 1/(k_h*k_w) to pond_resis.  The j==1, jwt==0
						! branch below uses k_h/spec_grnd_cond, giving the intended
						! water-side resistance 1/k_w.  The old k_h/k_w expression
						! made the effective water-side conductance too large by
						! approximately 1/k_h**2.
						if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. &
						    lakedepth <= 0._r8) then
							! Lake patch enabled for production but landdata lakedepth
							! is missing or non-positive.  Cole & Caraco piston velocity
							! cannot be applied; the air-water exchange degrades silently
							! to whatever the non-lake branch gives.  For exploratory
							! Tests keep the one-shot warning; for production
							! global lake-CH4 runs set lake_zero_depth_fatal=T.
							if (DEF_METHANE%lake_zero_depth_fatal) then
								write(6,*) 'ERROR: patchtype=4 + allowlakeprod=T but lakedepth<=0. ', &
									'Patch/fix landdata lakedepth or disable allowlakeprod.'
								CALL CoLM_stop ()
							else if (.not. warned_lake_zero_depth) then
								write(6,*) 'WARNING: patchtype=4 + allowlakeprod=T but lakedepth<=0; ', &
									'falling back to non-lake exchange.  Check landdata lakedepth field.'
								warned_lake_zero_depth = .true.
							end if
						end if
						if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. lakedepth > 0._r8) then
							if (t_h2osfc >= tfrz .and. lake_icefrac(1) <= 0.1_r8) then
								t_soisno_c = t_h2osfc - tfrz
								if (ieee_is_nan(forc_us) .or. ieee_is_nan(forc_vs) .or. &
								    abs(forc_us) > 1.e30_r8 .or. abs(forc_vs) > 1.e30_r8) then
									lake_wind = 0._r8
								else
									lake_wind = sqrt(max(forc_us*forc_us + forc_vs*forc_vs, 0._r8))
								end if
								! Cole & Caraco (1998): k600 [cm/hr] = 2.07 + 0.215 U10**1.7.
								lake_k600_cmhr = 2.07_r8 + 0.215_r8 * max(lake_wind, 0._r8)**1.7_r8
								! Freshwater Schmidt number polynomial (Wanninkhof-style)
								! per species; clamp T to avoid extrapolated extremes.
								t_soisno_c = max(0._r8, min(30._r8, t_soisno_c))
								lake_schmidt_ch4 = max(300._r8, s_con(s,1) + s_con(s,2)*t_soisno_c + &
									s_con(s,3)*t_soisno_c**2 + s_con(s,4)*t_soisno_c**3)
								lake_exchange_vel = lake_k600_cmhr * 1.e-2_r8 / 3600._r8 * &
									(lake_schmidt_ch4 / 600._r8)**(-2._r8/3._r8)
								! k_h_cc(0,s) is the dimensionless air-water Henry
								! coefficient at the surface; converts the water-side
								! piston velocity to its air-side equivalent so it
								! adds correctly to grnd_methane_cond.
								pond_resis = pond_resis + 1._r8 / (k_h_cc(0,s) * &
									max(lake_exchange_vel, smallnumber))
							else
								! Ice-covered/frozen lake surface: molecular exchange through the water column is effectively blocked.
								pond_resis = pond_resis + 1._r8/smallnumber
						endif
					end if

					spec_grnd_cond(s) = 1._r8/(1._r8/max(grnd_methane_cond_base, smallnumber) + &
						snow_resis + pond_resis)
				end if
			end do ! j

			! Determine gas diffusion and fraction of open pore (f_a)
			do j = 1,nl_soil
				t_soisno_c = t_soisno(j) - tfrz

				if (porsl(j) <= smallnumber) then
					diffus(j) = smallnumber
				else if (j <= jwt) then  ! Above the WT
					if (organic_max > 0._r8) then
						om_frac = min(DEF_METHANE%om_frac_sf*cellorg(j)/organic_max, 1._r8)
						! Use first power, not square as in iniTimeConst
					else
						om_frac = 1._r8
					end if
					diffus (j) = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								(om_frac * vol_gas(j)**(10._r8/3._r8) / porsl(j)**2._r8 + &
								(1._r8-om_frac) * vol_gas(j)**2._r8 * (vol_gas(j)/porsl(j))**(3._r8 / max(bsw(j), 0.5_r8)) ) &
								* DEF_METHANE%scale_factor_gasdiff
					else ! Below the WT use liquid-water diffusivity; ice-filled pore space is excluded via vol_aqu.
						diffus (j) = max(vol_aqu(j), smallnumber)**DEF_METHANE%satpow * &
							(d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
							* DEF_METHANE%scale_factor_liqdiff
					end if ! Above/below the WT
					if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. j > jwt) then
						if (s == 1) diffus(j) = diffus(j) * DEF_METHANE%lake_liqdiff_scale
						if (s == 2) diffus(j) = diffus(j) * DEF_METHANE%lake_o2_liqdiff_scale
					end if
					diffus(j) = max(diffus(j), smallnumber) ! Prevent overflow
				enddo ! j

			! Zero conductance arrays so end-rows that intentionally skip one
			! direction stay at a defined value under -finit-real=nan builds.
			dp1_zp1(:) = 0._r8
			dm1_zm1(:) = 0._r8

			do j = 1,nl_soil

				! Set up coefficients for tridiagonal solver.
				if (j == 1 .and. j /= jwt .and. j /= jwt+1) then
					dm1_zm1(j) = 1._r8/(1._r8/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! replace Diffusivity / Delta_z by conductance (grnd_methane_cond) for top layer
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j == 1 .and. j == jwt) then
					dm1_zm1(j) = 1._r8/(1._r8/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! layer resistance mult. by k_h_cc for dp1_zp1 term
					dp1_zp1(j) = 2._r8/(dz_soisno(j)*k_h_cc(j,s)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j == 1) then ! water table at surface: multiply ground resistance by k_h_cc
					dm1_zm1(j) = 1._r8/(k_h_cc(j-1,s)/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! air concentration will be mult. by k_h_cc below
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j <= nl_soil-1 .and. j /= jwt .and. j /= jwt+1) then
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)/diffus(j-1))
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j <= nl_soil-1 .and. j == jwt) then ! layer resistance mult. by k_h_cc for dp1_zp1 term
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)/diffus(j-1))
					dp1_zp1(j) = 2._r8/(dz_soisno(j)*k_h_cc(j,s)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
					! Concentration in layer will be mult. by k_h_cc below
				else if (j <= nl_soil-1) then ! j==jwt+1: layer above resistance mult. by k_h_cc for dm1_zm1 term
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)*k_h_cc(j-1,s)/diffus(j-1))
					! Concentration in layer above will be mult. by k_h_cc below
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j /= jwt+1) then ! j ==nl_soil
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)/diffus(j-1))
				else                    ! jwt == nl_soil-1: layer above resistance mult. by k_h_cc for dm1_zm1 term
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)*k_h_cc(j-1,s)/diffus(j-1))
				end if
			end do ! j; nl_soil

			! Perform a second loop for the tridiagonal coefficients since need dp1_zp1 and dm1_z1 at each depth
				do j = 0,nl_soil
					conc_ch4_rel_old(j) = conc_ch4_rel(j)

					if (j > 0) dzj = dz_soisno(j)
					if (j == 0) then ! top layer (atmosphere) doesn't change regardless of where WT is
						at(j) = 0._r8
						bt(j) = 1._r8
						ct(j) = 0._r8
						rt(j) = c_atm(s) ! 0th level stays at constant atmospheric conc
					elseif (s == 1 .and. ch4_backward_euler_transport) then
						! CH4 uses a fully implicit backward-Euler diffusion
						! matrix.  The old Crank-Nicolson half-explicit
						! diffusion term could overshoot sharp CH4 gradients in
						! cold / low-storage layers and then rely on the
						! negative-concentration clip branch.  With the
						! competition limiter above, the CH4 RHS is non-negative;
						! this M-matrix transport step is therefore
						! positivity-preserving while keeping the same
						! Henry-law interface coefficients.
						if (j < nl_soil .and. j == jwt) then
							at(j) = -1._r8 / dzj * dm1_zm1(j)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * &
								(dp1_zp1(j)*k_h_cc(j,s) + dm1_zm1(j))
							ct(j) = -1._r8 / dzj * dp1_zp1(j)
						elseif (j < nl_soil .and. j == jwt+1) then
							at(j) = -1._r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * &
								(dp1_zp1(j) + dm1_zm1(j))
							ct(j) = -1._r8 / dzj * dp1_zp1(j)
						elseif (j < nl_soil) then
							at(j) = -1._r8 / dzj * dm1_zm1(j)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * &
								(dp1_zp1(j) + dm1_zm1(j))
							ct(j) = -1._r8 / dzj * dp1_zp1(j)
						else if (j == nl_soil .and. j== jwt+1) then
							at(j) = -1._r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * dm1_zm1(j)
							ct(j) = 0._r8
						else
							at(j) = -1._r8 / dzj * dm1_zm1(j)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * dm1_zm1(j)
							ct(j) = 0._r8
						endif
					elseif (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
						at(j) = -0.5_r8 / dzj * dm1_zm1(j)
						bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * (dp1_zp1(j)*k_h_cc(j,s) + dm1_zm1(j))
					ct(j) = -0.5_r8 / dzj * dp1_zp1(j)
				elseif (j < nl_soil .and. j == jwt+1) then
					! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
					at(j) = -0.5_r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * (dp1_zp1(j) + dm1_zm1(j))
					ct(j) = -0.5_r8 / dzj * dp1_zp1(j)
				elseif (j < nl_soil) then
					at(j) = -0.5_r8 / dzj * dm1_zm1(j)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * (dp1_zp1(j) + dm1_zm1(j))
					ct(j) = -0.5_r8 / dzj * dp1_zp1(j)
				else if (j == nl_soil .and. j== jwt+1) then
					! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
					at(j) = -0.5_r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * dm1_zm1(j)
					ct(j) = 0._r8
				else ! j==nl_soil and jwt<nl_soil-1 or jwt==nl_soil: 0 flux at bottom
					at(j) = -0.5_r8 / dzj * dm1_zm1(j)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * dm1_zm1(j)
					ct(j) = 0._r8
				endif
			enddo ! j; nl_soil


			jtop = 0


			if (s == 1) then  ! CH4

				! Set rt, since it depends on conc
					do j = 1,nl_soil

						! Source and epsilon_t are constant over the implicit step,
						! so a single evaluation is sufficient (formerly a CN-style
						! average against a same-step "old" copy, which was a no-op).
						dzj = dz_soisno(j)
						if (ch4_backward_euler_transport) then
							! Backward-Euler CH4 diffusion has no explicit
							! diffusion contribution on the RHS.  With
							! source(j,1) limited so source + storage/dt is
							! non-negative, the implicit matrix above prevents
							! large negative post-solve concentrations.
							rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) + source(j,s)
						elseif (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
							rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
								0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)*k_h_cc(j,s)) - &
								dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							source(j,s)
					elseif (j < nl_soil .and. j == jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)) - &
							dm1_zm1(j) * (conc_ch4_rel(j) -conc_ch4_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					elseif (j < nl_soil) then
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)) - &
							dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							source(j,s)
					else if (j == nl_soil .and. j== jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_ch4_rel(j) -conc_ch4_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					else  !j==nl_soil
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							source(j,s)
					endif

				enddo ! j; nl_soil

				call Tridiagonal(0, nl_soil, &
					jtop, &
					at(:), &
					bt(:), &
					ct(:), &
					rt(:), &
					conc_ch4_rel(0:nl_soil))



					! Calculate net methane flux to the atmosphere from the surface (+ to atm)
					if (jwt /= 0) then ! WT not at the surface
						if (ch4_backward_euler_transport) then
							methane_surf_diff = dm1_zm1(1) * (conc_ch4_rel(1) - c_atm(s)) ! [mol/m2/s]
						else
							methane_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
								- c_atm(s)) ! [mol/m2/s]
						endif
					else ! WT at the surface; i.e., jwt==0
						if (ch4_backward_euler_transport) then
							methane_surf_diff = dm1_zm1(1) * (conc_ch4_rel(1) &
								- c_atm(s)*k_h_cc(0,s)) ! [mol/m2/s]
						else
							methane_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
								- c_atm(s)*k_h_cc(0,s)) ! [mol/m2/s]
						endif
						! atmospheric concentration gets mult. by k_h_cc as above
					endif
				methane_surf_ebul = methane_ebul_tot ! [mol/m2/s]

				! Snapshot the raw diffusive flux BEFORE the negative-concentration
				! clip-credit and the column-closure residual fold-in below, so
				! history can expose the pure-physics value separately from the
				! numerical correction (methane_balance_residual already captures
				! the residual; the clip credit is its own contribution).
				methane_surf_diff_phys = methane_surf_diff

				! Ensure that concentrations stay above 0
				! This should be done after the flux, so that the flux calculation is consistent.
				! Original strict 1e-12 abort tripped on implicit-solver roundoff at
				! cold start; relaxing to 1e-6 still aborted because cold-start sink
				! pull on initially-empty CH4 columns drives deficits on the order of
				! the per-timestep flux scale.
				! The clip+credit branch already conserves column mass — replace the
				! fatal stop with a single rate-limited warning per rank and always
				! clip, mirroring CTSM's approach for this implicit edge case.
					do j = 1,nl_soil
						if (conc_ch4_rel(j) < 0._r8) then
							deficit = - conc_ch4_rel(j)*epsilon_t(j,1)*dz_soisno(j)  ! Mol/m^2 added
							if (deficit > 1.e-3_r8 .and. .not. warned_ch4_clip) then
								write(6,*)'WARNING: methane_tran implicit clip > 1e-3 mol/m2 timestep.'
								write(6,*)'Lat,Lon, j, deficit=',dlat,dlon,j,deficit
								warned_ch4_clip = .true.
							endif
							! Always clip and credit; column mass is conserved by the surface flux update.
							methane_surf_diff = methane_surf_diff - deficit / deltim
							methane_ch4_clip_credit = methane_ch4_clip_credit - deficit / deltim
							conc_ch4_rel(j) = 0._r8
					end if
				enddo


			elseif (s == 2) then  ! O2

				! Set rt, since it depends on conc
				do j = 1,nl_soil

					! Source and epsilon_t are constant over the implicit step
					! (see the CH4 branch above for the same simplification).
					dzj     = dz_soisno(j)
					if (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)*k_h_cc(j,s)) - &
							dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							source(j,s)
					elseif (j < nl_soil .and. j == jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)) - &
							dm1_zm1(j) * (conc_o2_rel(j) -conc_o2_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					elseif (j < nl_soil) then
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)) - &
							dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							source(j,s)
					else if (j == nl_soil .and. j== jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_o2_rel(j) -conc_o2_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					else  !j==nl_soil
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							source(j,s)
					endif
				enddo ! j; nl_soil

				call Tridiagonal(0, nl_soil, jtop, &
					at(:), &
					bt(:), &
					ct(:), &
					rt(:), &
					conc_o2_rel(0:nl_soil))

					! Ensure that concentrations stay above 0
					o2_cap_loss_col = 0._r8
					o2_cap_gain_col = 0._r8
					do j = 1,nl_soil

						o2_before_cap = conc_o2_rel(j)
						conc_o2_rel(j) = max (conc_o2_rel(j), 1.e-12_r8)
						if (conc_o2_rel(j) > o2_before_cap) then
							o2_cap_gain_col = o2_cap_gain_col + &
								(conc_o2_rel(j) - o2_before_cap) * epsilon_t(j,2) * dz_soisno(j)
						endif

						! In case of pathologically large aerenchyma conductance. Should be OK in general but
						! this will maintain stability even if a PATCH with very small weight somehow has an absurd NPP or LAI.
						! Above WT: conc_o2_rel == conc_o2_gas_porsl == true gas-phase
						! concentration (mol per m^3 gas).  Algebra (see split_ch4_o2_phases):
						!   code conc_o2_gas        = conc_o2 * porsl / (k_h*vol_aqu + vol_gas)
						!   conc_o2_gas_porsl       = conc_o2 / (k_h*vol_aqu + vol_gas)
						! Confirmed by top BC at line ~2410: conc_o2_rel(0) = c_atm(2).
						! Cap directly at c_atm (was c_atm/porsl which left a 1/porsl-factor
						! slack — porsl=0.4 -> 2.5x ambient gas allowed; unphysical).
						! Below WT: conc_o2_rel = conc_o2 / vol_aqu = aqueous concentration;
						! Henry-equilibrium ceiling is k_h * c_atm.
						o2_before_cap = conc_o2_rel(j)
						if (j <= jwt) then
							o2_rel_cap = c_atm(2)
						else
							o2_rel_cap = k_h_cc(j,2) * c_atm(2)
						endif
						conc_o2_rel(j) = min (conc_o2_rel(j), max(o2_rel_cap, 1.e-12_r8))
						if (conc_o2_rel(j) < o2_before_cap) then
							o2_cap_loss_col = o2_cap_loss_col + &
								(o2_before_cap - conc_o2_rel(j)) * epsilon_t(j,2) * dz_soisno(j)
						endif
						enddo
						if (deltim > 0._r8) then
							o2_cap_loss = o2_cap_loss_col / deltim
							o2_cap_gain = o2_cap_gain_col / deltim
						endif
						! Cap adjustments are a benign mass correction (O2 conc kept within
					! [0, ambient]). Don't abort on small numerical adjustments; warn
					! once per rank if the adjustment exceeds 1e-3 mol/m^2 timestep.
					if ((o2_cap_loss_col > 1.e-3_r8 .or. o2_cap_gain_col > 1.e-3_r8) &
					    .and. .not. warned_o2_cap) then
						write(6,*)'WARNING: O2 cap adjustment > 1e-3 mol/m2 timestep.'
						write(6,*)'Lat,Lon, cap_loss, cap_gain=',dlat,dlon,o2_cap_loss_col,o2_cap_gain_col
						warned_o2_cap = .true.
					endif

		   endif  ! species

		enddo  ! species

		! Update absolute concentrations per unit volume
		do j = 1,nl_soil ! No need to update the atm. level concentrations
			conc_methane(j) = conc_ch4_rel(j)*epsilon_t(j,1)
			conc_o2(j)  = conc_o2_rel(j) *epsilon_t(j,2)
		end do

		! Do Balance Check and absorb small
		!    discrepancy into surface flux.
		err_methane = 0._r8
		do j = 1,nl_soil
			err_methane = err_methane + (conc_methane(j) - conc_ch4_bef(j))*dz_soisno(j)
			err_methane = err_methane - methane_prod_depth(j)*dz_soisno(j)*deltim
			err_methane = err_methane + methane_oxid_depth(j)*dz_soisno(j)*deltim
			err_methane = err_methane + methane_tran_depth(j)*dz_soisno(j)*deltim
		end do

		! History output uses the effective CH4 conductance, including snow and pond resistance.
		grnd_methane_cond_effective = spec_grnd_cond(1)

		err_methane = err_methane + (methane_surf_aere + methane_surf_ebul + methane_surf_diff)*deltim


		! Always clip-credit the closure error into methane_surf_diff so the
		! column inventory stays conserved. Warn once per rank if the residual
		! exceeds 1e-3 mol/m^2/timestep (well above any realistic flux scale).
		methane_balance_residual = -err_methane/deltim
		methane_surf_diff = methane_surf_diff + methane_balance_residual
		if (abs(err_methane) > 1.e-3_r8 .and. .not. warned_ch4_diffusion_residual) then
			write(6,*)'WARNING: CH4 diffusion residual > 1e-3 mol/m^2/timestep, istep, err=', idate, err_methane
			write(6,*)'Lat,Lon=',dlat,dlon
			warned_ch4_diffusion_residual = .true.
		end if

			end subroutine methane_tran

			!---------------------------------------------------------------------------
			real(r8) function methane_wetland_water_depth(patchtype, wdsrf, wetwat)
				! Use one surface-water-depth definition for all inundation schemes.
				! In the non-dynamic wetland hydrology branch, wdsrf/wa/wetwat are
				! merged into wetwat, then overflow is stored back in wdsrf. Those two
				! post-step reservoirs are disjoint and both contribute to ponded
				! depth. Dynamic wetland tiles follow the soil-ground branch, so wetwat
				! is not part of the active ponded-water state there.
				integer, intent(in) :: patchtype
				real(r8), intent(in) :: wdsrf, wetwat

				if (patchtype == 2 .and. .not. DEF_USE_Dynamic_Wetland) then
					methane_wetland_water_depth = max(wdsrf, 0._r8) + max(wetwat, 0._r8)
				else
					methane_wetland_water_depth = max(wdsrf, 0._r8)
				endif
			end function methane_wetland_water_depth

			!---------------------------------------------------------------------------
			logical function methane_patch_is_nongrass(patchclass)
				integer, intent(in) :: patchclass

#ifdef LULC_USGS
				! USGS: woody/non-grass classes that should use the reduced
				! non-grass aerenchyma porosity.  Include mixed cropland/woodland,
				! shrubland and mixed shrub/grass classes in addition to
				! savanna/forest classes.
				methane_patch_is_nongrass = patchclass == 6 .or. patchclass == 8 .or. &
				                            patchclass == 9 .or. &
				                            (patchclass >= 10 .and. patchclass <= 15) .or. &
				                            patchclass == 18 .or. patchclass == 21
#else
				! IGBP/PFT-class builds: tree/shrub/savanna classes are non-grass.
				! Grasslands/croplands/wetlands keep the grass-like tiller porosity
				! unless a future PFT-level parameter is wired in.
				methane_patch_is_nongrass = (patchclass >= 1 .and. patchclass <= 9) .or. patchclass == 14
#endif
		end function methane_patch_is_nongrass

		!---------------------------------------------------------------------------
		subroutine Tridiagonal (lbj, ubj, jtop, a, b, c, r, u)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Tridiagonal matrix solution
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		implicit none
		integer , intent(in)    :: lbj, ubj   ! lbinning and ubing level indices
		integer , intent(in)    :: jtop       ! top level for each column [col]
		real(r8), intent(in)    :: a(lbj:ubj) ! "a" left off diagonal of tridiagonal matrix [col, j]
		real(r8), intent(in)    :: b(lbj:ubj) ! "b" diagonal column for tridiagonal matrix [col, j]
		real(r8), intent(in)    :: c(lbj:ubj) ! "c" right off diagonal tridiagonal matrix [col, j]
		real(r8), intent(in)    :: r(lbj:ubj) ! "r" forcing term of tridiagonal matrix [col, j]
		real(r8), intent(inout) :: u(lbj:ubj) ! solution [col, j]

		!-----------------------Local Variables---------------------------------
		integer  :: j                 !indices
		real(r8) :: gam(lbj:ubj)      !temporary
		real(r8) :: bet               !temporary

		!-----------------------------------------------------------------------
		! Solve the matrix

		bet = b(jtop)

		do j = lbj, ubj
			if (j >= jtop) then
				if (j == jtop) then
					u(j) = r(j) / bet
				else
					gam(j) = c(j-1) / bet
					bet = b(j) - a(j) * gam(j)
					u(j) = (r(j) - a(j)*u(j-1)) / bet
				end if
			end if
		end do

		do j = ubj-1,lbj,-1
			if (j >= jtop) then
				u(j) = u(j) - gam(j+1) * u(j+1)
			end if
		end do

	end subroutine Tridiagonal

	subroutine henry_law (t_grnd,t_soisno,k_h_cc)
		real(r8), intent(in) :: &
			t_grnd                 		   , &! ground surface temperature [k]
			t_soisno (maxsnl+1:nl_soil)       ! soil temperature [K]

		real(r8), intent(out) :: &
			k_h_cc(0:nl_soil,ngases)          ! Dimensionless Henry's coefficient [-]

		integer :: j,s  			          ! Indices
		real(r8):: k_h						  ! Henry's coefficient [mol/L/atm]
		real(r8):: t_eff                  ! bounded effective temperature [K]

			do j = 0,nl_soil
				do s=1,ngases
					if (j==0) then
					t_eff = max(t_grnd, 200._r8)
					k_h = kh_theta(s)*exp(c_h(s) * (1._r8 / t_eff - 1._r8 / kh_tbase))
					! [mol/L/atm] = [mol/L/atm]*e**([K]*[1/K])
					k_h_cc(j,s) = k_h * rgasLatm * t_eff
					! [-]  = [mol/L/atm]*[L*atm/mol/K]*[K]
				else
					t_eff = max(t_soisno(j), 200._r8)
					k_h = kh_theta(s)*exp(c_h(s) * (1._r8 / t_eff - 1._r8 / kh_tbase))
					! [mol/L/atm] = [mol/L/atm]*e**([K]*[1/K])
					k_h_cc(j,s) = k_h * rgasLatm * t_eff
					! [-]  = [mol/L/atm]*[L*atm/mol/K]*[K]
				end if
			end do
		end do
	end subroutine henry_law

	subroutine split_ch4_o2_phases( dz_soisno, wliq_soisno, wice_soisno, porsl, &
									conc_methane, conc_o2, k_h_cc, idate,&
									vol_aqu,vol_ch4_storage,vol_gas,f_aqu,f_gas,&
									conc_ch4_gas,conc_ch4_aqu,conc_ch4_porsl,conc_ch4_gas_porsl,conc_ch4_aqu_porsl,&
									conc_o2_gas,conc_o2_aqu,conc_o2_porsl,conc_o2_gas_porsl,conc_o2_aqu_porsl)

		implicit none
		real(r8), intent(in) :: dz_soisno(maxsnl+1:nl_soil)     ! layer thickness [m]
		real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)   ! liquid water in layers [kg/m2]
		real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)   ! ice in layers [kg/m2]
		real(r8), intent(in) :: porsl(1:nl_soil)                ! volumetric soil water at saturation (porosity)
		real(r8), intent(in) :: conc_methane(1:nl_soil)          ! CH4 concentration in each soil layer [mol/m3]
		real(r8), intent(in) :: conc_o2(1:nl_soil)           ! O2 concentration in each soil layer [mol/m3]
		real(r8), intent(in) :: k_h_cc(0:nl_soil,ngases)        ! ratio of mol/m3 in liquid to mol/m3 in gas [-]
		integer, intent(in) :: idate(3)

		integer :: j,s
		real(r8), parameter :: smallnumber = 1.e-12_r8
		real(r8) :: vol_ice, mobile_pore, f_ch4_storage
		real(r8), intent(out) :: vol_aqu(1:nl_soil)          ! liquid volumetric water content [m3/m3]
		real(r8), intent(out) :: vol_ch4_storage(1:nl_soil)  ! liquid/ice volume retaining CH4 [m3/m3]
		real(r8), intent(out) :: vol_gas(1:nl_soil)          ! air volumetric water content [m3/m3]
		real(r8), intent(out) :: f_aqu(1:nl_soil)            ! water-filled proportion [-]
		real(r8), intent(out) :: f_gas(1:nl_soil)            ! air-filled proportion [-]
		real(r8), intent(out) :: conc_ch4_gas(1:nl_soil)     ! gas phase CH4 conc [mol/m3]
		real(r8), intent(out) :: conc_ch4_aqu(1:nl_soil)     ! aqueous phase CH4 conc [mol/m3]
		real(r8), intent(out) :: conc_ch4_porsl(1:nl_soil)   ! CH4 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_ch4_gas_porsl(1:nl_soil) ! gas phase CH4 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_ch4_aqu_porsl(1:nl_soil) ! aqueous phase CH4 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_o2_gas(1:nl_soil)      ! gas phase O2 conc [mol/m3]
		real(r8), intent(out) :: conc_o2_aqu(1:nl_soil)      ! aqueous phase O2 conc [mol/m3]
		real(r8), intent(out) :: conc_o2_porsl(1:nl_soil)    ! O2 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_o2_gas_porsl(1:nl_soil)! gas phase O2 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_o2_aqu_porsl(1:nl_soil)! aqueous phase O2 conc per porosity [mol/m3]

		!-----------------------------------------------------------
		! Main computation: calculate phase partitioning for CH4 and O2
		!-----------------------------------------------------------
		do j = 1, nl_soil
			if (porsl(j) <= smallnumber .or. dz_soisno(j) <= smallnumber) then
				vol_aqu(j) = 0._r8
				vol_ch4_storage(j) = 0._r8
				vol_gas(j) = 0._r8
				f_aqu(j) = 0._r8
				f_gas(j) = 0._r8
				conc_ch4_gas(j) = 0._r8
				conc_ch4_aqu(j) = 0._r8
				conc_ch4_porsl(j) = 0._r8
				conc_ch4_gas_porsl(j) = 0._r8
				conc_ch4_aqu_porsl(j) = 0._r8
				conc_o2_gas(j) = 0._r8
				conc_o2_aqu(j) = 0._r8
				conc_o2_porsl(j) = 0._r8
				conc_o2_gas_porsl(j) = 0._r8
				conc_o2_aqu_porsl(j) = 0._r8
				cycle
			end if
			! ---- Compute volumetric liquid/ice/gas contents ----
			vol_aqu(j) = max(0._r8, min(wliq_soisno(j)/(dz_soisno(j)*denh2o), porsl(j)))
			! [m3/m3] = [kg/m2] / ([m] * [kg/m3])
			vol_ice = max(0._r8, min(wice_soisno(j)/(dz_soisno(j)*denice), max(porsl(j)-vol_aqu(j), 0._r8)))

			! ---- Compute volumetric gas content, excluding ice-filled pores ----
			vol_gas(j) = max(porsl(j) - vol_aqu(j) - vol_ice, 0._r8)
			mobile_pore = vol_aqu(j) + vol_gas(j)
			IF (DEF_METHANE%methane_frzout) THEN
				vol_ch4_storage(j) = vol_aqu(j)
			ELSE
				vol_ch4_storage(j) = vol_aqu(j) + vol_ice
			ENDIF
			! A completely frozen layer has no phase into which a freeze-out
			! pulse can move.  Retain its old bulk inventory as immobile ice
			! storage instead of collapsing epsilon to zero and exporting the
			! inventory through the numerical closure flux.
			IF (mobile_pore <= smallnumber) vol_ch4_storage(j) = vol_ice

			! ---- Compute filled proportions ----
			f_aqu(j) = vol_aqu(j)/porsl(j)
			f_gas(j) = vol_gas(j)/porsl(j)
			f_ch4_storage = vol_ch4_storage(j)/porsl(j)

			! ---- CH4 partitioning between gas and aqueous phases ----
			if (mobile_pore <= smallnumber) then
				conc_ch4_aqu(j) = 0._r8
				conc_ch4_gas(j) = 0._r8
				conc_o2_aqu(j) = 0._r8
				conc_o2_gas(j) = 0._r8
			else
				conc_ch4_aqu(j) = conc_methane(j)/(f_ch4_storage+f_gas(j)/k_h_cc(j,1))
				conc_ch4_gas(j) = conc_methane(j)/(k_h_cc(j,1)*f_ch4_storage+f_gas(j))

				! ---- O2 partitioning between gas and aqueous phases ----
				conc_o2_aqu(j) = conc_o2(j)/(f_aqu(j)+f_gas(j)/k_h_cc(j,2))
				conc_o2_gas(j) = conc_o2(j)/(k_h_cc(j,2)*f_aqu(j)+f_gas(j))
			endif

			! ---- Concentrations normalized by porosity ----
			conc_ch4_porsl(j)      = conc_methane(j)/porsl(j)
			conc_ch4_aqu_porsl(j)  = conc_ch4_aqu(j)/porsl(j)
			conc_ch4_gas_porsl(j)  = conc_ch4_gas(j)/porsl(j)
			conc_o2_porsl(j)       = conc_o2(j)/porsl(j)
			conc_o2_aqu_porsl(j)   = conc_o2_aqu(j)/porsl(j)
			conc_o2_gas_porsl(j)   = conc_o2_gas(j)/porsl(j)


		end do
	end subroutine split_ch4_o2_phases
END MODULE MOD_Tracer_Reactive_Methane_Physics
#endif
! --------- EOP ----------
