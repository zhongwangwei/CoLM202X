#include <define.h>
#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Methane_Driver
!=======================================================================
! Driver wrapper for Methane reactive-tracer physics.
!=======================================================================

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: methane_driver

CONTAINS

	SUBROUTINE methane_driver (istep,i,idate,patchclass,patchtype,deltim,lb,snl,dlon,dlat,&!input
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai,rootr,fsatmax,fsatdcf,frcsat,f_h2osfc)

		use MOD_Precision
		use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
		use MOD_Tracer_Methane_Const
		use MOD_Namelist, only : DEF_USE_VariablySaturatedFlow
		use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad,z_soi,zi_soi,dz_soi
		use MOD_Tracer_Methane_Physics
		use MOD_SPMD_Task
			USE MOD_Tracer_Methane_BgcLink, only: tracer_ch4_bgc_patch_inputs, organic_max
		USE MOD_Tracer_Methane_Microbes, only: methane_microbes_step, &
		     microbial_prod_potential, microbial_oxid_potential

		! ----- Methane tracer state (module-level allocatables) -----
		USE MOD_Tracer_Methane_State, only: net_methane, methane_prod_depth, o2_decomp_depth, &
		     co2_decomp_depth, methane_oxid_depth, o2_oxid_depth, co2_oxid_depth, &
		     methane_aere_depth, methane_tran_depth, o2_aere_depth, co2_aere_depth, methane_ebul_depth, &
		     o2stress, methane_stress, &
		     methane_surf_flux_tot, methane_surf_aere, methane_surf_ebul, methane_surf_diff, &
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
		     methane_surf_ebul_lake, methane_surf_diff_lake, &
		     methane_prod_tot_lake, methane_oxid_tot_lake, methane_ebul_tot_lake, &
		     co2_decomp_tot_lake, co2_oxid_tot_lake, co2_net_tot_lake, &
		     totcol_methane_lake, grnd_methane_cond_lake, conc_o2_lake, conc_methane_lake, &
		     c_atm, forc_pmethanem, layer_sat_lag, lake_soilc, &
		     annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
		     tempavg_agnpp, tempavg_bgnpp, annsum_counter, &
		     tempavg_somhr, tempavg_finrw, &
		     fsat_bef, finundated_lag, methane_dfsat_tot

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

		real(r8):: &
				crootfr  (1:nl_soil)     , &! fraction of roots for carbon in each soil layer
				pH                       , &! soil water pH
				cellorg  (1:nl_soil)     , &! column 3D org (kg/m3 organic matter)
				t_h2osfc             	    ! surface water temperature

		real(r8) :: somhr_loc, lithr_loc, rr_loc, agnpp_loc, bgnpp_loc, annsum_npp_loc
		real(r8) :: hr_vr_loc(1:nl_soil), fphr_loc(1:nl_soil)
		real(r8) :: o_scalar_loc(1:nl_soil), pot_f_nit_vr_loc(1:nl_soil)
		logical  :: bgc_inputs_ready
		logical, save :: warned_missing_bgc_inputs = .false.

		t_h2osfc = max(t_grnd, tfrz)

			CALL tracer_ch4_bgc_patch_inputs (i, rootfr, crootfr, pH, cellorg, &
			     somhr_loc, lithr_loc, hr_vr_loc, rr_loc, agnpp_loc, bgnpp_loc, &
			     annsum_npp_loc, fphr_loc, o_scalar_loc, pot_f_nit_vr_loc, bgc_inputs_ready)
			IF (patchtype == 0 .and. .not. bgc_inputs_ready) THEN
				IF (p_is_master .and. .not. warned_missing_bgc_inputs) THEN
					write(6,*) 'WARNING: Soil methane running with sanitized BGC defaults before BGC/PFT inputs are initialized.'
					warned_missing_bgc_inputs = .true.
				ENDIF
			ENDIF
					! f_h2osfc is maintained by the water module before methane_driver.
				! Lake CH4 uses the CTSM-style sediment-carbon pathway, not the optional
				! soil microbial-pool override; keep lake microbial pools from evolving as
				! diagnostics-only state.
				IF (patchtype /= 4) THEN
					CALL methane_microbes_step (i, deltim, t_soisno(1:nl_soil), &
					     conc_o2(1:nl_soil,i), conc_methane(1:nl_soil,i), hr_vr_loc, cellorg)
				ENDIF

		CALL methane (istep,idate(1:3),patchclass,patchtype,lb,snl,dlon,dlat,deltim,&
		z_soisno(maxsnl+1:),dz_soisno(maxsnl+1:),zi_soisno(maxsnl:),t_soisno(maxsnl+1:),&
		t_grnd,wliq_soisno(maxsnl+1:),wice_soisno(maxsnl+1:),&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai,rootr,&
		annsum_npp_loc, rr_loc,&
		fsatmax,fsatdcf,frcsat,&
		agnpp_loc, bgnpp_loc, somhr_loc,&
		crootfr(1:nl_soil), lithr_loc, hr_vr_loc(1:nl_soil), o_scalar_loc(1:nl_soil), &
		fphr_loc(1:nl_soil), pot_f_nit_vr_loc(1:nl_soil), pH,&
			cellorg(1:nl_soil),t_h2osfc,organic_max,&
			microbial_prod_potential(1:nl_soil,i), microbial_oxid_potential(1:nl_soil,i), &
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
		methane_surf_flux_tot(i), methane_surf_aere(i), methane_surf_ebul(i), methane_surf_diff(i), &
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
		methane_surf_ebul_lake(i), methane_surf_diff_lake(i), &
		methane_prod_tot_lake(i), methane_oxid_tot_lake(i), methane_ebul_tot_lake(i), &
		co2_decomp_tot_lake(i), co2_oxid_tot_lake(i), co2_net_tot_lake(i), &
		totcol_methane_lake(i), grnd_methane_cond_lake(i), conc_o2_lake(1:nl_soil,i), conc_methane_lake(1:nl_soil,i), &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm(1:3,i), forc_pmethanem(i), layer_sat_lag(1:nl_soil,i), lake_soilc(1:nl_soil,i), &
		annavg_agnpp(i), annavg_bgnpp(i), annavg_somhr(i), annavg_finrw(i), &
		tempavg_agnpp(i), tempavg_bgnpp(i), annsum_counter(i), tempavg_somhr(i), &
		tempavg_finrw(i), fsat_bef(i), finundated_lag(i), methane_dfsat_tot(i), f_h2osfc)

	END SUBROUTINE methane_driver

END MODULE MOD_Tracer_Methane_Driver
#endif
