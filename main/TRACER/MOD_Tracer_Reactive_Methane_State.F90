#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_State
!=======================================================================
! Methane reactive-tracer state and restart fields.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_nan
   USE MOD_Vars_Global, only: nl_soil, spval, dz_soi

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Explicit exports: keep module data private by default and expose
   ! only the current cross-module methane interface/state fields.
   PUBLIC :: accumulate_methane_lake_substep_diagnostics
   PUBLIC :: allocate_methane_state
   PUBLIC :: annavg_agnpp
   PUBLIC :: annavg_bgnpp
   PUBLIC :: annavg_finrw
   PUBLIC :: annavg_somhr
   PUBLIC :: annsum_counter
   PUBLIC :: biome_f_methane_patch
   PUBLIC :: biome_redoxlag_patch
   PUBLIC :: c_atm
   PUBLIC :: co2_aere_depth
   PUBLIC :: co2_aere_depth_sat
   PUBLIC :: co2_aere_depth_unsat
   PUBLIC :: co2_aere_tot
   PUBLIC :: co2_decomp_depth
   PUBLIC :: co2_decomp_depth_lake
   PUBLIC :: co2_decomp_depth_sat
   PUBLIC :: co2_decomp_depth_unsat
   PUBLIC :: co2_decomp_tot
   PUBLIC :: co2_decomp_tot_lake
   PUBLIC :: co2_decomp_tot_sat
   PUBLIC :: co2_decomp_tot_unsat
   PUBLIC :: co2_net_tot
   PUBLIC :: co2_net_tot_lake
   PUBLIC :: co2_net_tot_sat
   PUBLIC :: co2_net_tot_unsat
   PUBLIC :: co2_oxid_depth
   PUBLIC :: co2_oxid_depth_lake
   PUBLIC :: co2_oxid_depth_sat
   PUBLIC :: co2_oxid_depth_unsat
   PUBLIC :: co2_oxid_tot
   PUBLIC :: co2_oxid_tot_lake
   PUBLIC :: co2_oxid_tot_sat
   PUBLIC :: co2_oxid_tot_unsat
   PUBLIC :: compute_f_h2osfc
   PUBLIC :: conc_methane
   PUBLIC :: conc_methane_lake
   PUBLIC :: conc_methane_sat
   PUBLIC :: conc_methane_unsat
   PUBLIC :: conc_o2
   PUBLIC :: conc_o2_lake
   PUBLIC :: conc_o2_sat
   PUBLIC :: conc_o2_unsat
   PUBLIC :: deallocate_methane_state
   PUBLIC :: f_h2osfc
   PUBLIC :: f_inund_flood_depth_patch
   PUBLIC :: f_inund_flood_patch
   PUBLIC :: f_inund_levee_patch
   PUBLIC :: finundated_lag
   PUBLIC :: forc_pmethanem
   PUBLIC :: fsat_bef
   PUBLIC :: grnd_methane_cond
   PUBLIC :: grnd_methane_cond_lake
   PUBLIC :: grnd_methane_cond_sat
   PUBLIC :: grnd_methane_cond_unsat
   PUBLIC :: init_methane_wetland_fraction_cache
   PUBLIC :: initialize_methane_lake_soilc_from_surface
   PUBLIC :: lake_soilc
   PUBLIC :: layer_sat_lag
   PUBLIC :: methane_aere_depth
   PUBLIC :: methane_aere_depth_sat
   PUBLIC :: methane_aere_depth_unsat
   PUBLIC :: methane_balance_residual
   PUBLIC :: methane_ch4_clip_credit
   PUBLIC :: restart_ch4_clip_credit_mass
   PUBLIC :: methane_dfsat_tot
   PUBLIC :: methane_ebul_depth
   PUBLIC :: methane_ebul_depth_lake
   PUBLIC :: methane_ebul_depth_sat
   PUBLIC :: methane_ebul_depth_unsat
   PUBLIC :: methane_ebul_tot
   PUBLIC :: methane_ebul_tot_lake
   PUBLIC :: methane_ebul_tot_sat
   PUBLIC :: methane_ebul_tot_unsat
   PUBLIC :: methane_finundated
   PUBLIC :: methane_oxid_depth
   PUBLIC :: methane_oxid_depth_lake
   PUBLIC :: methane_oxid_depth_sat
   PUBLIC :: methane_oxid_depth_unsat
   PUBLIC :: methane_oxid_tot
   PUBLIC :: methane_oxid_tot_lake
   PUBLIC :: methane_oxid_tot_sat
   PUBLIC :: methane_oxid_tot_unsat
   PUBLIC :: methane_prod_depth
   PUBLIC :: methane_prod_depth_lake
   PUBLIC :: methane_prod_depth_sat
   PUBLIC :: methane_prod_depth_unsat
   PUBLIC :: methane_prod_tot
   PUBLIC :: methane_prod_tot_lake
   PUBLIC :: methane_prod_tot_sat
   PUBLIC :: methane_prod_tot_unsat
   PUBLIC :: methane_soil_finundated
   PUBLIC :: methane_soil_zwt
   PUBLIC :: methane_stress
   PUBLIC :: methane_stress_sat
   PUBLIC :: methane_stress_unsat
   PUBLIC :: methane_surf_aere
   PUBLIC :: methane_surf_aere_sat
   PUBLIC :: methane_surf_aere_unsat
   PUBLIC :: methane_surf_diff
   PUBLIC :: methane_surf_diff_lake
   PUBLIC :: methane_surf_diff_phys
   PUBLIC :: methane_surf_diff_phys_sat
   PUBLIC :: methane_surf_diff_phys_unsat
   PUBLIC :: methane_surf_diff_sat
   PUBLIC :: methane_surf_diff_unsat
   PUBLIC :: methane_surf_ebul
   PUBLIC :: methane_surf_ebul_lake
   PUBLIC :: methane_surf_ebul_sat
   PUBLIC :: methane_surf_ebul_unsat
   PUBLIC :: methane_surf_flux_lake
   PUBLIC :: methane_surf_flux_rice
   PUBLIC :: methane_surf_flux_soil
   PUBLIC :: methane_surf_flux_tot
   PUBLIC :: methane_surf_flux_tot_lake
   PUBLIC :: methane_surf_flux_tot_phys
   PUBLIC :: methane_surf_flux_tot_sat
   PUBLIC :: methane_surf_flux_tot_unsat
   PUBLIC :: methane_surf_flux_wetland
   PUBLIC :: methane_tran_depth
   PUBLIC :: methane_tran_depth_sat
   PUBLIC :: methane_tran_depth_unsat
   PUBLIC :: net_methane
   PUBLIC :: net_methane_sat
   PUBLIC :: net_methane_unsat
   PUBLIC :: o2_aere_depth
   PUBLIC :: o2_aere_depth_sat
   PUBLIC :: o2_aere_depth_unsat
   PUBLIC :: o2_cap_gain
   PUBLIC :: o2_cap_loss
   PUBLIC :: o2_decomp_depth
   PUBLIC :: o2_decomp_depth_sat
   PUBLIC :: o2_decomp_depth_unsat
   PUBLIC :: o2_oxid_depth
   PUBLIC :: o2_oxid_depth_sat
   PUBLIC :: o2_oxid_depth_unsat
   PUBLIC :: o2stress
   PUBLIC :: o2stress_sat
   PUBLIC :: o2stress_unsat
   PUBLIC :: read_methane_restart
   PUBLIC :: remap_methane_lulcc_state
   PUBLIC :: publish_methane_levee_flood_patch
   PUBLIC :: publish_methane_flood_patch
   PUBLIC :: save_methane_lulcc_state
   PUBLIC :: tempavg_agnpp
   PUBLIC :: tempavg_bgnpp
   PUBLIC :: tempavg_finrw
   PUBLIC :: tempavg_somhr
   PUBLIC :: totcol_methane
   PUBLIC :: totcol_methane_lake
   PUBLIC :: totcol_methane_sat
   PUBLIC :: totcol_methane_unsat
   PUBLIC :: wetland_frac_per_patch
   PUBLIC :: write_methane_restart

   ! Public read-only data: external modules may inspect these arrays,
   ! but all writes stay inside this module through its APIs.
   PROTECTED :: f_inund_levee_patch, f_inund_flood_patch, f_inund_flood_depth_patch, &
      wetland_frac_per_patch, f_h2osfc

   ! -------------------- field declarations --------------------
   !!!! --------------------------------------------------------------------------------------------------------
   !!!!                                         sum data
   !!!! --------------------------------------------------------------------------------------------------------
	real(r8), allocatable :: net_methane           (:) ! average net methane correction to CO2 flux (mol/m2/s)
	real(r8), allocatable :: methane_prod_depth      (:,:) ! production of CH4 in each soil layer  (mol/m3/s)
	real(r8), allocatable :: o2_decomp_depth     (:,:) ! O2 consumption during decomposition in each soil layer (mol/m3/s)
	real(r8), allocatable :: co2_decomp_depth    (:,:) ! diagnostic CO2 from decomposition/methanogenesis before O2-stress scaling (mol/m3/s)
	real(r8), allocatable :: methane_oxid_depth      (:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
	real(r8), allocatable :: o2_oxid_depth       (:,:) ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)
	real(r8), allocatable :: co2_oxid_depth      (:,:) ! CO2 production from CH4 oxidation (mol/m3/s)
	real(r8), allocatable :: methane_aere_depth      (:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
	real(r8), allocatable :: methane_tran_depth      (:,:) ! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
	real(r8), allocatable :: o2_aere_depth       (:,:) ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)
	real(r8), allocatable :: co2_aere_depth      (:,:) ! CO2 aerenchyma diagnostic flux (mol/m3/s)
	real(r8), allocatable :: methane_ebul_depth      (:,:) ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)
	real(r8), allocatable :: o2stress            (:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
	real(r8), allocatable :: methane_stress           (:,:) ! Ratio of methane available to the total per-timestep methane sinks
	real(r8), allocatable :: methane_surf_flux_tot     (:) ! CH4 flux to atm incl. numerical corrections (mol/m2/s)
	real(r8), allocatable :: methane_surf_flux_tot_phys(:) ! CH4 physical flux before CH4 clip/residual corrections (mol/m2/s)
	real(r8), allocatable :: methane_surf_aere         (:) ! Output: Total column CH4 aerenchyma (mol/m2/s)
	real(r8), allocatable :: methane_surf_ebul         (:) ! Output: CH4 ebullition to atmosphere (mol/m2/s)
	real(r8), allocatable :: methane_surf_diff         (:) ! Output: CH4 diffusion flux plus numerical closure (mol/m2/s)
	real(r8), allocatable :: methane_surf_diff_phys    (:) ! Output: CH4 pure-physics diffusion before clip/residual (mol/m2/s)
	real(r8), allocatable :: methane_balance_residual (:) ! numerical CH4 closure flux credited to methane_surf_diff (mol/m2/s)
	real(r8), allocatable :: methane_ch4_clip_credit(:) ! negative CH4 clip credited to methane_surf_diff (mol/m2/s)
	real(r8), allocatable :: restart_ch4_clip_credit_mass(:) ! restart negative CH4 sanitation credit (mol/m2 impulse)
	real(r8), allocatable :: o2_cap_loss              (:) ! O2 removed by post-solve physical cap (mol/m2/s)
	real(r8), allocatable :: o2_cap_gain              (:) ! O2 added by post-solve nonnegative floor (mol/m2/s)
	real(r8), allocatable :: methane_ebul_tot        (:) ! Output: Total column CH4 ebullition (mol/m2/s)
	real(r8), allocatable :: methane_prod_tot        (:) ! Output: Total column CH4 production (mol/m2/s)
	real(r8), allocatable :: methane_oxid_tot        (:) ! Output: Total column CH4 oxidation (mol/m2/s)
	real(r8), allocatable :: co2_decomp_tot       (:) ! total diagnostic CO2 from decomposition/methanogenesis before O2-stress scaling (mol/m2/s)
	real(r8), allocatable :: co2_oxid_tot         (:) ! total CO2 from CH4 oxidation (mol/m2/s)
	real(r8), allocatable :: co2_aere_tot         (:) ! total CO2 aerenchyma diagnostic (mol/m2/s)
	real(r8), allocatable :: co2_net_tot          (:) ! net diagnosed CO2 source from CH4 module (mol/m2/s)
   real(r8), allocatable :: totcol_methane             (:) ! total methane in soil column, start of timestep (mol/m2)
	real(r8), allocatable :: grnd_methane_cond         (:) ! tracer conductance for boundary layer (m/s)
   real(r8), allocatable :: conc_o2             (:,:) ! O2 conc in each soil layer (mol/m3)
	real(r8), allocatable :: conc_methane            (:,:) ! CH4 conc in each soil layer (mol/m3)
   !!!! --------------------------------------------------------------------------------------------------------

   !!!! --------------------------------------------------------------------------------------------------------
   !!!!                                         sum data (unsaturated / saturated)
   !!!! --------------------------------------------------------------------------------------------------------
   real(r8), allocatable :: net_methane_unsat           (:)  ! average unsaturated net methane correction to CO2 flux (mol/m2/s)
   real(r8), allocatable :: net_methane_sat             (:)  ! average saturated   net methane correction to CO2 flux (mol/m2/s)

   real(r8), allocatable :: methane_prod_depth_unsat      (:,:)  ! production of CH4 in each soil layer (unsaturated)  (mol/m3/s)
   real(r8), allocatable :: methane_prod_depth_sat        (:,:)  ! production of CH4 in each soil layer (saturated)    (mol/m3/s)

   real(r8), allocatable :: o2_decomp_depth_unsat     (:,:)  ! O2 consumption during decomposition (unsaturated)    (mol/m3/s)
   real(r8), allocatable :: o2_decomp_depth_sat       (:,:)  ! O2 consumption during decomposition (saturated)      (mol/m3/s)
   real(r8), allocatable :: co2_decomp_depth_unsat    (:,:)  ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsat) (mol/m3/s)
   real(r8), allocatable :: co2_decomp_depth_sat      (:,:)  ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (sat)   (mol/m3/s)

   real(r8), allocatable :: methane_oxid_depth_unsat      (:,:)  ! CH4 oxidation rate (unsaturated)                     (mol/m3/s)
   real(r8), allocatable :: methane_oxid_depth_sat        (:,:)  ! CH4 oxidation rate (saturated)                       (mol/m3/s)

   real(r8), allocatable :: o2_oxid_depth_unsat       (:,:)  ! O2 oxidation rate (unsaturated)                      (mol/m3/s)
   real(r8), allocatable :: o2_oxid_depth_sat         (:,:)  ! O2 oxidation rate (saturated)                        (mol/m3/s)
   real(r8), allocatable :: co2_oxid_depth_unsat      (:,:)  ! CO2 from CH4 oxidation (unsaturated)                 (mol/m3/s)
   real(r8), allocatable :: co2_oxid_depth_sat        (:,:)  ! CO2 from CH4 oxidation (saturated)                   (mol/m3/s)

   real(r8), allocatable :: methane_aere_depth_unsat      (:,:)  ! CH4 loss via aerenchyma (unsaturated)                (mol/m3/s)
   real(r8), allocatable :: methane_aere_depth_sat        (:,:)  ! CH4 loss via aerenchyma (saturated)                  (mol/m3/s)

   real(r8), allocatable :: methane_tran_depth_unsat      (:,:)  ! CH4 loss via transpiration (unsaturated)             (mol/m3/s)
   real(r8), allocatable :: methane_tran_depth_sat        (:,:)  ! CH4 loss via transpiration (saturated)               (mol/m3/s)

   real(r8), allocatable :: o2_aere_depth_unsat       (:,:)  ! O2 gain via aerenchyma (unsaturated)                 (mol/m3/s)
   real(r8), allocatable :: o2_aere_depth_sat         (:,:)  ! O2 gain via aerenchyma (saturated)                   (mol/m3/s)
   real(r8), allocatable :: co2_aere_depth_unsat      (:,:)  ! CO2 aerenchyma diagnostic (unsaturated)              (mol/m3/s)
   real(r8), allocatable :: co2_aere_depth_sat        (:,:)  ! CO2 aerenchyma diagnostic (saturated)                (mol/m3/s)

   real(r8), allocatable :: methane_ebul_depth_unsat      (:,:)  ! CH4 loss via ebullition (unsaturated)                (mol/m3/s)
   real(r8), allocatable :: methane_ebul_depth_sat        (:,:)  ! CH4 loss via ebullition (saturated)                  (mol/m3/s)

   real(r8), allocatable :: o2stress_unsat            (:,:)  ! O2 availability/demand ratio (unsaturated)           (-)
   real(r8), allocatable :: o2stress_sat              (:,:)  ! O2 availability/demand ratio (saturated)             (-)

   real(r8), allocatable :: methane_stress_unsat           (:,:)  ! CH4 availability/sinks ratio (unsaturated)           (-)
   real(r8), allocatable :: methane_stress_sat             (:,:)  ! CH4 availability/sinks ratio (saturated)             (-)

   real(r8), allocatable :: methane_surf_flux_tot_unsat     (:)  ! CH4 flux to atmosphere (unsaturated)                 (mol/m2/s)
   real(r8), allocatable :: methane_surf_flux_tot_sat       (:)  ! CH4 flux to atmosphere (saturated)                   (mol/m2/s)

   real(r8), allocatable :: methane_surf_aere_unsat         (:)  ! Total column CH4 aerenchyma (unsaturated)            (mol/m2/s)
   real(r8), allocatable :: methane_surf_aere_sat           (:)  ! Total column CH4 aerenchyma (saturated)              (mol/m2/s)

   real(r8), allocatable :: methane_surf_ebul_unsat         (:)  ! CH4 ebullition to atmosphere (unsaturated)           (mol/m2/s)
   real(r8), allocatable :: methane_surf_ebul_sat           (:)  ! CH4 ebullition to atmosphere (saturated)             (mol/m2/s)

   real(r8), allocatable :: methane_surf_diff_unsat         (:)  ! CH4 surface diffusive flux (unsaturated)             (mol/m2/s)
   real(r8), allocatable :: methane_surf_diff_sat           (:)  ! CH4 surface diffusive flux (saturated)               (mol/m2/s)
   real(r8), allocatable :: methane_surf_diff_phys_unsat    (:)  ! CH4 pure-physics diff (unsaturated)                  (mol/m2/s)
   real(r8), allocatable :: methane_surf_diff_phys_sat      (:)  ! CH4 pure-physics diff (saturated)                    (mol/m2/s)

   real(r8), allocatable :: methane_ebul_tot_unsat          (:)  ! Total column CH4 ebullition (unsaturated)            (mol/m2/s)
   real(r8), allocatable :: methane_ebul_tot_sat            (:)  ! Total column CH4 ebullition (saturated)              (mol/m2/s)

   real(r8), allocatable :: methane_prod_tot_unsat          (:)  ! Total column CH4 production (unsaturated)            (mol/m2/s)
   real(r8), allocatable :: methane_prod_tot_sat            (:)  ! Total column CH4 production (saturated)              (mol/m2/s)

   real(r8), allocatable :: methane_oxid_tot_unsat          (:)  ! Total column CH4 oxidation (unsaturated)             (mol/m2/s)
   real(r8), allocatable :: methane_oxid_tot_sat            (:)  ! Total column CH4 oxidation (saturated)               (mol/m2/s)
   real(r8), allocatable :: co2_decomp_tot_unsat       (:)  ! total CO2 decomp/methanogenesis (unsaturated)        (mol/m2/s)
   real(r8), allocatable :: co2_decomp_tot_sat         (:)  ! total CO2 decomp/methanogenesis (saturated)          (mol/m2/s)
   real(r8), allocatable :: co2_oxid_tot_unsat         (:)  ! total CO2 from CH4 oxidation (unsaturated)           (mol/m2/s)
   real(r8), allocatable :: co2_oxid_tot_sat           (:)  ! total CO2 from CH4 oxidation (saturated)             (mol/m2/s)
   real(r8), allocatable :: co2_net_tot_unsat          (:)  ! net diagnosed CO2 source (unsaturated)               (mol/m2/s)
   real(r8), allocatable :: co2_net_tot_sat            (:)  ! net diagnosed CO2 source (saturated)                 (mol/m2/s)

   real(r8), allocatable :: totcol_methane_unsat             (:)  ! total methane in soil column, start (unsaturated)     (mol/m2)
   real(r8), allocatable :: totcol_methane_sat               (:)  ! total methane in soil column, start (saturated)       (mol/m2)

   real(r8), allocatable :: grnd_methane_cond_unsat         (:)  ! tracer conductance for boundary layer (unsaturated)  (m/s)
   real(r8), allocatable :: grnd_methane_cond_sat           (:)  ! tracer conductance for boundary layer (saturated)    (m/s)

   real(r8), allocatable :: conc_o2_unsat             (:,:)  ! O2 concentration in each soil layer (unsaturated)    (mol/m3)
   real(r8), allocatable :: conc_o2_sat               (:,:)  ! O2 concentration in each soil layer (saturated)      (mol/m3)

	   real(r8), allocatable :: conc_methane_unsat            (:,:)  ! CH4 concentration in each soil layer (unsaturated)   (mol/m3)
	   real(r8), allocatable :: conc_methane_sat              (:,:)  ! CH4 concentration in each soil layer (saturated)     (mol/m3)
	   !!!! --------------------------------------------------------------------------------------------------------

	   !!!! --------------------------------------------------------------------------------------------------------
	   !!!!                                         lake data (CTSM alignment)
	   !!!! --------------------------------------------------------------------------------------------------------
	   real(r8), allocatable :: methane_prod_depth_lake       (:,:)  ! lake CH4 production rate by layer (mol/m3/s)
	   real(r8), allocatable :: methane_oxid_depth_lake       (:,:)  ! lake CH4 oxidation rate by layer (mol/m3/s)
	   real(r8), allocatable :: methane_ebul_depth_lake       (:,:)  ! lake CH4 ebullition loss by layer (mol/m3/s)
	   real(r8), allocatable :: co2_decomp_depth_lake       (:,:)  ! lake CO2 from sediment decomposition/methanogenesis (mol/m3/s)
	   real(r8), allocatable :: co2_oxid_depth_lake         (:,:)  ! lake CO2 from CH4 oxidation (mol/m3/s)
	   real(r8), allocatable :: methane_surf_ebul_lake        (:)    ! lake surface ebullition flux (mol/m2/s)
	   real(r8), allocatable :: methane_surf_diff_lake        (:)    ! lake surface diffusive flux (mol/m2/s)
	   real(r8), allocatable :: methane_surf_flux_tot_lake    (:)    ! lake total surface CH4 flux = ebul + diff (mol/m2/s)
	   real(r8), allocatable :: methane_prod_tot_lake         (:)    ! lake total CH4 production (mol/m2/s)
	   real(r8), allocatable :: methane_oxid_tot_lake         (:)    ! lake total CH4 oxidation (mol/m2/s)
	   real(r8), allocatable :: methane_ebul_tot_lake         (:)    ! lake total ebullition (mol/m2/s)
	   real(r8), allocatable :: co2_decomp_tot_lake         (:)    ! lake total CO2 decomp/methanogenesis (mol/m2/s)
	   real(r8), allocatable :: co2_oxid_tot_lake           (:)    ! lake total CO2 oxidation product (mol/m2/s)
	   real(r8), allocatable :: co2_net_tot_lake            (:)    ! lake net diagnosed CO2 source (mol/m2/s)
	   real(r8), allocatable :: totcol_methane_lake           (:)    ! lake CH4 column stock (mol/m2)
	   real(r8), allocatable :: grnd_methane_cond_lake        (:)    ! lake-atmosphere CH4 conductance (m/s)
	   real(r8), allocatable :: conc_o2_lake                  (:,:)  ! lake O2 concentration by layer (mol/m3)
	   real(r8), allocatable :: conc_methane_lake             (:,:)  ! lake CH4 concentration by layer (mol/m3)
	   !!!! --------------------------------------------------------------------------------------------------------

	   real(r8), allocatable :: c_atm               (:,:) ! CH4, O2, CO2 atmospheric conc  (mol/m3)
	real(r8), allocatable :: forc_pmethanem            (:) ! CH4 concentration in atmos. (pascals)
	real(r8), allocatable :: layer_sat_lag       (:,:)
	real(r8), allocatable :: lake_soilc          (:,:) ! total soil organic matter found in level (gC / m3)
   real(r8), allocatable :: annavg_agnpp          (:) ! annual average above-ground NPP (gC/m2/s)
	real(r8), allocatable :: annavg_bgnpp          (:) ! annual average below-ground NPP (gC/m2/s)
	real(r8), allocatable :: annavg_somhr          (:) ! annual average SOM heterotrophic resp. (gC/m2/s)
	real(r8), allocatable :: annavg_finrw          (:) ! respiration-weighted annual average of finundated
   real(r8), allocatable :: tempavg_agnpp         (:) ! temporary average above-ground NPP (gC/m2/s)
	real(r8), allocatable :: tempavg_bgnpp         (:) ! temporary average below-ground NPP (gC/m2/s)
	real(r8), allocatable :: annsum_counter        (:) ! seconds since last annual accumulator turnover
	real(r8), allocatable :: tempavg_somhr         (:) ! temporary average SOM heterotrophic resp. (gC/m2/s)
	real(r8), allocatable :: tempavg_finrw         (:) ! respiration-weighted annual average of finundated

   real(r8), allocatable :: fsat_bef              (:) ! finundated from previous timestep
   real(r8), allocatable :: finundated_lag        (:) ! time-lagged fractional inundated area
   real(r8), allocatable :: methane_dfsat_tot         (:) ! CH4 flux to atm due to decreasing finundated [mol/m2/s]

   ! f_h2osfc: fractional area of surface water (0-1, dimensionless).
   ! Source-repo CLM5 microtopography-based prognostic h2osfc scheme.
   ! Maintained by compute_f_h2osfc (this module) before each methane_driver call.
	   real(r8), allocatable :: f_h2osfc              (:) ! fraction of surface water [-]
	   ! Diagnostics used to audit CH4 inundation choices and patchtype
	   ! contributions in history output.  These are not restart-critical state.
	   real(r8), allocatable :: methane_finundated        (:) ! actual CH4 finundated used by physics [-]
	   real(r8), allocatable :: methane_soil_finundated   (:) ! finundated on active soil/rice patches only [-]
	   real(r8), allocatable :: methane_soil_zwt          (:) ! zwt on active soil/rice patches only [m]
	   real(r8), allocatable :: methane_surf_flux_wetland (:) ! wetland contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_flux_soil    (:) ! non-rice soil contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_flux_lake    (:) ! lake contribution to CH4 surface flux [mol/m2/s]
	   real(r8), allocatable :: methane_surf_flux_rice    (:) ! rice-paddy contribution to CH4 surface flux [mol/m2/s]
	   ! Per-patch floodplain fraction (0-1) from GridRiverLakeFlow's levee
   ! diagnostic (levee_floodarea / topo_area), exposed to methane scheme 7.
   ! Default 0; populated by MOD_Grid_RiverLakeFlow via
   ! publish_levee_fldfrc_to_patches when GridRiverLakeFlow is compiled in and
   ! the catchment->patch mapping is available.
   real(r8), allocatable :: f_inund_levee_patch  (:) ! floodplain frac (-)

   ! General flood inundation fraction (levee+floodplain) for methane scheme 7.
   ! Published by MOD_Grid_RiverLakeFlow::publish_fldfrc_to_patches.
	   real(r8), allocatable :: f_inund_flood_patch  (:) ! flood frac (-)
	   real(r8), allocatable :: f_inund_flood_depth_patch (:) ! floodplain water depth [m]

	   ! Per-patch cache of the local wetland fraction in the parent
	   ! element/cell.  Used to convert grid/ucat-relative routing flood
	   ! fractions into wetland-relative finundated on patchtype==2:
	   !   fin_wet = min(1, flood_grid / max(wetland_frac_per_patch, 0.01)).
	   ! Built once after landpatch/patchtype are available, and rebuilt after
	   ! LULCC remaps.
	   real(r8), allocatable :: wetland_frac_per_patch (:) ! wetland area / active land area [-]

	   ! Biome-specific f_methane per-patch (set by Driver via get_biome_f_methane).
	   ! When DEF_METHANE%use_biome_f_methane=.true., methane_prod consumes this
	   ! array instead of the global DEF_METHANE%f_methane scalar.
	   real(r8), allocatable :: biome_f_methane_patch (:) ! [mol CH4 / mol CO2 anaerobic decomp]

	   ! Biome-specific redoxlag per-patch (set by Driver via get_biome_redoxlag).
	   ! When DEF_METHANE%use_biome_redoxlag=.true., methane uses this array
	   ! for finundated_lag time constant instead of global DEF_METHANE%redoxlag.
	   real(r8), allocatable :: biome_redoxlag_patch (:)  ! [days]

   ! LULCC remap snapshot. These arrays are private old-layout copies used
   ! to rebuild Methane state after landpatch/numpatch changes.
   logical :: methane_lulcc_snapshot_valid = .false.
   real(r8), allocatable :: lulcc_conc_o2_old(:,:), lulcc_conc_methane_old(:,:)
   real(r8), allocatable :: lulcc_co2_decomp_depth_old(:,:), lulcc_co2_oxid_depth_old(:,:), lulcc_co2_aere_depth_old(:,:)
   real(r8), allocatable :: lulcc_co2_decomp_depth_unsat_old(:,:), lulcc_co2_decomp_depth_sat_old(:,:)
   real(r8), allocatable :: lulcc_co2_oxid_depth_unsat_old(:,:), lulcc_co2_oxid_depth_sat_old(:,:)
   real(r8), allocatable :: lulcc_co2_aere_depth_unsat_old(:,:), lulcc_co2_aere_depth_sat_old(:,:)
   real(r8), allocatable :: lulcc_co2_decomp_depth_lake_old(:,:), lulcc_co2_oxid_depth_lake_old(:,:)
   real(r8), allocatable :: lulcc_conc_o2_unsat_old(:,:), lulcc_conc_o2_sat_old(:,:)
   real(r8), allocatable :: lulcc_conc_methane_unsat_old(:,:), lulcc_conc_methane_sat_old(:,:)
   real(r8), allocatable :: lulcc_conc_o2_lake_old(:,:), lulcc_conc_methane_lake_old(:,:)
   real(r8), allocatable :: lulcc_totcol_methane_old(:), lulcc_totcol_methane_unsat_old(:)
   real(r8), allocatable :: lulcc_co2_decomp_tot_old(:), lulcc_co2_oxid_tot_old(:), lulcc_co2_aere_tot_old(:), lulcc_co2_net_tot_old(:)
   real(r8), allocatable :: lulcc_co2_decomp_tot_unsat_old(:), lulcc_co2_decomp_tot_sat_old(:)
   real(r8), allocatable :: lulcc_co2_oxid_tot_unsat_old(:), lulcc_co2_oxid_tot_sat_old(:)
   real(r8), allocatable :: lulcc_co2_net_tot_unsat_old(:), lulcc_co2_net_tot_sat_old(:)
   real(r8), allocatable :: lulcc_co2_decomp_tot_lake_old(:), lulcc_co2_oxid_tot_lake_old(:), lulcc_co2_net_tot_lake_old(:)
   real(r8), allocatable :: lulcc_totcol_methane_sat_old(:), lulcc_totcol_methane_lake_old(:)
   real(r8), allocatable :: lulcc_grnd_methane_cond_old(:), lulcc_grnd_methane_cond_unsat_old(:)
   real(r8), allocatable :: lulcc_grnd_methane_cond_sat_old(:), lulcc_grnd_methane_cond_lake_old(:)
   real(r8), allocatable :: lulcc_layer_sat_lag_old(:,:), lulcc_lake_soilc_old(:,:)
   real(r8), allocatable :: lulcc_annavg_agnpp_old(:), lulcc_annavg_bgnpp_old(:)
   real(r8), allocatable :: lulcc_annavg_somhr_old(:), lulcc_annavg_finrw_old(:)
   real(r8), allocatable :: lulcc_tempavg_agnpp_old(:), lulcc_tempavg_bgnpp_old(:)
   real(r8), allocatable :: lulcc_annsum_counter_old(:), lulcc_tempavg_somhr_old(:)
   real(r8), allocatable :: lulcc_tempavg_finrw_old(:), lulcc_fsat_bef_old(:)
   real(r8), allocatable :: lulcc_finundated_lag_old(:), lulcc_methane_dfsat_tot_old(:)
   real(r8), allocatable :: lulcc_f_h2osfc_old(:), lulcc_forc_pmethanem_old(:)
	   real(r8), allocatable :: lulcc_c_atm_old(:,:)

	   ! Temporary lake-substep history buffers.  Lake methane runs on the
	   ! WATERBODY physics substep, while the normal history accumulator is
	   ! called once per land timestep.  These buffers keep time-weighted
	   ! diagnostic rates over the substeps and write the per-timestep mean
	   ! back to the module fields before accumulate_methane_fluxes samples
	   ! them.  Prognostic states such as concentrations and lake_soilc are
	   ! intentionally left at their final substep values.
	   integer, parameter :: methane_lake_substep_n2d = 44
	   integer, parameter :: methane_lake_substep_n1d = 49
	   real(r8), allocatable :: methane_lake_substep_acc2d(:,:)
	   real(r8), allocatable :: methane_lake_substep_acc1d(:)
	   integer :: methane_lake_substep_cached_ipatch = -1
	   integer :: methane_lake_substep_next_isub = 1

   ! -------------------- API --------------------

   INTERFACE credit_methane_restart_clip
      MODULE PROCEDURE credit_methane_restart_clip_1d
      MODULE PROCEDURE credit_methane_restart_clip_2d
   END INTERFACE credit_methane_restart_clip

CONTAINS

   SUBROUTINE credit_methane_restart_clip_1d (field, credit_mass)
      ! Accumulate positive CH4 mass introduced by restart sanitation.
      ! The in-step clip path reports this as a negative surface-flux credit;
      ! restart has no timestep length, so keep the boundary correction as a
      ! mol/m2 impulse budget item and fold its sign into methane_ch4_clip_credit
      ! for diagnostics that inspect the state immediately after restart.
      real(r8), intent(in)    :: field(:)
      real(r8), intent(inout) :: credit_mass(:)
      integer :: i, n

      n = min(size(field), size(credit_mass))
      DO i = 1, n
         IF (field(i) < 0._r8) credit_mass(i) = credit_mass(i) - field(i)
      ENDDO
   END SUBROUTINE credit_methane_restart_clip_1d

   SUBROUTINE credit_methane_restart_clip_2d (field, credit_mass)
      ! Convert negative concentration clips to a column impulse using the
      ! canonical soil-layer thickness. This mirrors the in-step CH4 clip-credit
      ! accounting while avoiding an unbudgeted restart mass insertion.
      real(r8), intent(in)    :: field(:,:)
      real(r8), intent(inout) :: credit_mass(:)
      integer :: i, j, nlev, npatch
      real(r8) :: dz

      nlev = min(size(field,1), nl_soil)
      npatch = min(size(field,2), size(credit_mass))
      DO i = 1, npatch
         DO j = 1, nlev
            dz = max(dz_soi(j), 0._r8)
            IF (field(j,i) < 0._r8) credit_mass(i) = credit_mass(i) - field(j,i) * dz
         ENDDO
      ENDDO
   END SUBROUTINE credit_methane_restart_clip_2d

   SUBROUTINE allocate_methane_state (numpatch)
      integer, intent(in) :: numpatch
      ! Keep zero-length arrays allocated on IO/control ranks.  Methane
      ! restart reads use collective vector I/O; if numpatch==0 ranks return
      ! here, workers enter ncio_read_vector while IO ranks skip it and MPI
      ! can deadlock during global startup.

      allocate (net_methane                 (numpatch)); net_methane            (:) = 0._r8
      allocate (methane_prod_depth      (nl_soil,numpatch)); methane_prod_depth       (:,:) = 0._r8
      allocate (o2_decomp_depth     (nl_soil,numpatch)); o2_decomp_depth      (:,:) = 0._r8
      allocate (co2_decomp_depth    (nl_soil,numpatch)); co2_decomp_depth     (:,:) = 0._r8
      allocate (methane_oxid_depth      (nl_soil,numpatch)); methane_oxid_depth       (:,:) = 0._r8
      allocate (o2_oxid_depth       (nl_soil,numpatch)); o2_oxid_depth        (:,:) = 0._r8
      allocate (co2_oxid_depth      (nl_soil,numpatch)); co2_oxid_depth       (:,:) = 0._r8
      allocate (methane_aere_depth      (nl_soil,numpatch)); methane_aere_depth       (:,:) = 0._r8
      allocate (methane_tran_depth      (nl_soil,numpatch)); methane_tran_depth       (:,:) = 0._r8
      allocate (o2_aere_depth       (nl_soil,numpatch)); o2_aere_depth        (:,:) = 0._r8
      allocate (co2_aere_depth      (nl_soil,numpatch)); co2_aere_depth       (:,:) = 0._r8
      allocate (methane_ebul_depth      (nl_soil,numpatch)); methane_ebul_depth       (:,:) = 0._r8
      allocate (o2stress            (nl_soil,numpatch)); o2stress             (:,:) = 1.0_r8
      allocate (methane_stress           (nl_soil,numpatch)); methane_stress            (:,:) = 1.0_r8
      allocate (methane_surf_flux_tot           (numpatch)); methane_surf_flux_tot      (:) = 0._r8
      allocate (methane_surf_flux_tot_phys      (numpatch)); methane_surf_flux_tot_phys (:) = 0._r8
      allocate (methane_surf_aere               (numpatch)); methane_surf_aere          (:) = 0._r8
      allocate (methane_surf_ebul               (numpatch)); methane_surf_ebul          (:) = 0._r8
      allocate (methane_surf_diff               (numpatch)); methane_surf_diff          (:) = 0._r8
      allocate (methane_surf_diff_phys          (numpatch)); methane_surf_diff_phys     (:) = 0._r8
      allocate (methane_balance_residual       (numpatch)); methane_balance_residual  (:) = 0._r8
      allocate (methane_ch4_clip_credit        (numpatch)); methane_ch4_clip_credit   (:) = 0._r8
      allocate (restart_ch4_clip_credit_mass (numpatch)); restart_ch4_clip_credit_mass(:) = 0._r8
      allocate (o2_cap_loss                    (numpatch)); o2_cap_loss               (:) = 0._r8
      allocate (o2_cap_gain                    (numpatch)); o2_cap_gain               (:) = 0._r8
      allocate (methane_ebul_tot                (numpatch)); methane_ebul_tot           (:) = 0._r8
      allocate (methane_prod_tot                (numpatch)); methane_prod_tot           (:) = 0._r8
      allocate (methane_oxid_tot                (numpatch)); methane_oxid_tot           (:) = 0._r8
      allocate (co2_decomp_tot                  (numpatch)); co2_decomp_tot             (:) = 0._r8
      allocate (co2_oxid_tot                    (numpatch)); co2_oxid_tot               (:) = 0._r8
      allocate (co2_aere_tot                    (numpatch)); co2_aere_tot               (:) = 0._r8
      allocate (co2_net_tot                     (numpatch)); co2_net_tot                (:) = 0._r8

      allocate (totcol_methane                   (numpatch)); totcol_methane              (:) = 0._r8
      allocate (grnd_methane_cond               (numpatch)); grnd_methane_cond          (:) = 1.e-6_r8
      ! Keep the combined cold-start O2 concentration consistent with the
      ! split unsat/sat pools.  Optional microbial pools may read conc_o2 on
      ! the first step before fsat_bef is valid.
      allocate (conc_o2             (nl_soil,numpatch)); conc_o2              (:,:) = 1.0_r8
      allocate (conc_methane            (nl_soil,numpatch)); conc_methane             (:,:) = 1.0e-6_r8
      !!!! --------------------------------------------------------------------------------------------------------

      !!!! --------------------------------------------------------------------------------------------------------
      !!!!                                         sum data (unsaturated / saturated)
      !!!! --------------------------------------------------------------------------------------------------------
      allocate (net_methane_unsat           (numpatch)); net_methane_unsat        (:)   = 0._r8
      allocate (net_methane_sat             (numpatch)); net_methane_sat          (:)   = 0._r8

      allocate (methane_prod_depth_unsat  (nl_soil,numpatch)); methane_prod_depth_unsat   (:,:) = 0._r8
      allocate (methane_prod_depth_sat    (nl_soil,numpatch)); methane_prod_depth_sat     (:,:) = 0._r8

      allocate (o2_decomp_depth_unsat (nl_soil,numpatch)); o2_decomp_depth_unsat  (:,:) = 0._r8
      allocate (o2_decomp_depth_sat   (nl_soil,numpatch)); o2_decomp_depth_sat    (:,:) = 0._r8
      allocate (co2_decomp_depth_unsat(nl_soil,numpatch)); co2_decomp_depth_unsat (:,:) = 0._r8
      allocate (co2_decomp_depth_sat  (nl_soil,numpatch)); co2_decomp_depth_sat   (:,:) = 0._r8

      allocate (methane_oxid_depth_unsat  (nl_soil,numpatch)); methane_oxid_depth_unsat   (:,:) = 0._r8
      allocate (methane_oxid_depth_sat    (nl_soil,numpatch)); methane_oxid_depth_sat     (:,:) = 0._r8

      allocate (o2_oxid_depth_unsat   (nl_soil,numpatch)); o2_oxid_depth_unsat    (:,:) = 0._r8
      allocate (o2_oxid_depth_sat     (nl_soil,numpatch)); o2_oxid_depth_sat      (:,:) = 0._r8
      allocate (co2_oxid_depth_unsat  (nl_soil,numpatch)); co2_oxid_depth_unsat   (:,:) = 0._r8
      allocate (co2_oxid_depth_sat    (nl_soil,numpatch)); co2_oxid_depth_sat     (:,:) = 0._r8

      allocate (methane_aere_depth_unsat  (nl_soil,numpatch)); methane_aere_depth_unsat   (:,:) = 0._r8
      allocate (methane_aere_depth_sat    (nl_soil,numpatch)); methane_aere_depth_sat     (:,:) = 0._r8

      allocate (methane_tran_depth_unsat  (nl_soil,numpatch)); methane_tran_depth_unsat   (:,:) = 0._r8
      allocate (methane_tran_depth_sat    (nl_soil,numpatch)); methane_tran_depth_sat     (:,:) = 0._r8

      allocate (o2_aere_depth_unsat   (nl_soil,numpatch)); o2_aere_depth_unsat    (:,:) = 0._r8
      allocate (o2_aere_depth_sat     (nl_soil,numpatch)); o2_aere_depth_sat      (:,:) = 0._r8
      allocate (co2_aere_depth_unsat  (nl_soil,numpatch)); co2_aere_depth_unsat   (:,:) = 0._r8
      allocate (co2_aere_depth_sat    (nl_soil,numpatch)); co2_aere_depth_sat     (:,:) = 0._r8

      allocate (methane_ebul_depth_unsat  (nl_soil,numpatch)); methane_ebul_depth_unsat   (:,:) = 0._r8
      allocate (methane_ebul_depth_sat    (nl_soil,numpatch)); methane_ebul_depth_sat     (:,:) = 0._r8

      allocate (o2stress_unsat        (nl_soil,numpatch)); o2stress_unsat         (:,:) = 1.0_r8
      allocate (o2stress_sat          (nl_soil,numpatch)); o2stress_sat           (:,:) = 1.0_r8

      allocate (methane_stress_unsat       (nl_soil,numpatch)); methane_stress_unsat        (:,:) = 1.0_r8
      allocate (methane_stress_sat         (nl_soil,numpatch)); methane_stress_sat          (:,:) = 1.0_r8

      allocate (methane_surf_flux_tot_unsat      (numpatch)); methane_surf_flux_tot_unsat  (:)   = 0._r8
      allocate (methane_surf_flux_tot_sat        (numpatch)); methane_surf_flux_tot_sat    (:)   = 0._r8

      allocate (methane_surf_aere_unsat          (numpatch)); methane_surf_aere_unsat      (:)   = 0._r8
      allocate (methane_surf_aere_sat            (numpatch)); methane_surf_aere_sat        (:)   = 0._r8

      allocate (methane_surf_ebul_unsat          (numpatch)); methane_surf_ebul_unsat      (:)   = 0._r8
      allocate (methane_surf_ebul_sat            (numpatch)); methane_surf_ebul_sat        (:)   = 0._r8

      allocate (methane_surf_diff_unsat          (numpatch)); methane_surf_diff_unsat      (:)   = 0._r8
      allocate (methane_surf_diff_sat            (numpatch)); methane_surf_diff_sat        (:)   = 0._r8
      allocate (methane_surf_diff_phys_unsat     (numpatch)); methane_surf_diff_phys_unsat (:)   = 0._r8
      allocate (methane_surf_diff_phys_sat       (numpatch)); methane_surf_diff_phys_sat   (:)   = 0._r8

      allocate (methane_ebul_tot_unsat           (numpatch)); methane_ebul_tot_unsat       (:)   = 0._r8
      allocate (methane_ebul_tot_sat             (numpatch)); methane_ebul_tot_sat         (:)   = 0._r8

      allocate (methane_prod_tot_unsat           (numpatch)); methane_prod_tot_unsat       (:)   = 0._r8
      allocate (methane_prod_tot_sat             (numpatch)); methane_prod_tot_sat         (:)   = 0._r8

      allocate (methane_oxid_tot_unsat           (numpatch)); methane_oxid_tot_unsat       (:)   = 0._r8
      allocate (methane_oxid_tot_sat             (numpatch)); methane_oxid_tot_sat         (:)   = 0._r8
      allocate (co2_decomp_tot_unsat             (numpatch)); co2_decomp_tot_unsat         (:)   = 0._r8
      allocate (co2_decomp_tot_sat               (numpatch)); co2_decomp_tot_sat           (:)   = 0._r8
      allocate (co2_oxid_tot_unsat               (numpatch)); co2_oxid_tot_unsat           (:)   = 0._r8
      allocate (co2_oxid_tot_sat                 (numpatch)); co2_oxid_tot_sat             (:)   = 0._r8
      allocate (co2_net_tot_unsat                (numpatch)); co2_net_tot_unsat            (:)   = 0._r8
      allocate (co2_net_tot_sat                  (numpatch)); co2_net_tot_sat              (:)   = 0._r8

      allocate (totcol_methane_unsat              (numpatch)); totcol_methane_unsat          (:)   = 0._r8
      allocate (totcol_methane_sat                (numpatch)); totcol_methane_sat            (:)   = 0._r8

      allocate (grnd_methane_cond_unsat          (numpatch)); grnd_methane_cond_unsat      (:)   = 1.e-6_r8
      allocate (grnd_methane_cond_sat            (numpatch)); grnd_methane_cond_sat        (:)   = 1.e-6_r8

      ! Physical cold-start concentrations for restart files without Methane state.
      allocate (conc_o2_unsat        (nl_soil,numpatch)); conc_o2_unsat           (:,:) = 1.0_r8
      allocate (conc_o2_sat          (nl_soil,numpatch)); conc_o2_sat             (:,:) = 1.0_r8

	      allocate (conc_methane_unsat       (nl_soil,numpatch)); conc_methane_unsat          (:,:) = 1.0e-6_r8
	      allocate (conc_methane_sat         (nl_soil,numpatch)); conc_methane_sat            (:,:) = 1.0e-6_r8
	      !!!! --------------------------------------------------------------------------------------------------------

	      !!!! --------------------------------------------------------------------------------------------------------
	      !!!!                                         lake data (CTSM alignment)
	      !!!! --------------------------------------------------------------------------------------------------------
	      allocate (methane_prod_depth_lake  (nl_soil,numpatch)); methane_prod_depth_lake     (:,:) = 0._r8
	      allocate (methane_oxid_depth_lake  (nl_soil,numpatch)); methane_oxid_depth_lake     (:,:) = 0._r8
	      allocate (methane_ebul_depth_lake  (nl_soil,numpatch)); methane_ebul_depth_lake     (:,:) = 0._r8
	      allocate (co2_decomp_depth_lake(nl_soil,numpatch)); co2_decomp_depth_lake     (:,:) = 0._r8
	      allocate (co2_oxid_depth_lake  (nl_soil,numpatch)); co2_oxid_depth_lake       (:,:) = 0._r8
	      allocate (methane_surf_ebul_lake          (numpatch)); methane_surf_ebul_lake       (:)   = 0._r8
	      allocate (methane_surf_diff_lake          (numpatch)); methane_surf_diff_lake       (:)   = 0._r8
	      allocate (methane_surf_flux_tot_lake      (numpatch)); methane_surf_flux_tot_lake   (:)   = 0._r8
	      allocate (methane_prod_tot_lake           (numpatch)); methane_prod_tot_lake        (:)   = 0._r8
	      allocate (methane_oxid_tot_lake           (numpatch)); methane_oxid_tot_lake        (:)   = 0._r8
	      allocate (methane_ebul_tot_lake           (numpatch)); methane_ebul_tot_lake        (:)   = 0._r8
	      allocate (co2_decomp_tot_lake           (numpatch)); co2_decomp_tot_lake        (:) = 0._r8
	      allocate (co2_oxid_tot_lake             (numpatch)); co2_oxid_tot_lake          (:) = 0._r8
	      allocate (co2_net_tot_lake              (numpatch)); co2_net_tot_lake           (:) = 0._r8
	      allocate (totcol_methane_lake             (numpatch)); totcol_methane_lake          (:)   = 0._r8
	      allocate (grnd_methane_cond_lake          (numpatch)); grnd_methane_cond_lake       (:)   = 1.e-6_r8
	      allocate (conc_o2_lake             (nl_soil,numpatch)); conc_o2_lake                (:,:) = 1.0_r8
	      allocate (conc_methane_lake        (nl_soil,numpatch)); conc_methane_lake           (:,:) = 0._r8
	      !!!! --------------------------------------------------------------------------------------------------------

	      allocate (c_atm                     (3,numpatch)); c_atm                (:,:) = 0._r8
      allocate (forc_pmethanem                  (numpatch)); forc_pmethanem             (:) = 0._r8
      allocate (layer_sat_lag       (nl_soil,numpatch)); layer_sat_lag        (:,:) = spval
      allocate (lake_soilc          (nl_soil,numpatch)); lake_soilc           (:,:) = 0._r8
      allocate (annavg_agnpp                (numpatch)); annavg_agnpp           (:) = 0._r8
      allocate (annavg_bgnpp                (numpatch)); annavg_bgnpp           (:) = 0._r8
      allocate (annavg_somhr                (numpatch)); annavg_somhr           (:) = 0._r8
      ! annavg_finrw is a prerequisite for the seasonal inundation factor (SIF).
      ! Use spval until a full annual accumulator has been completed; otherwise a
      ! cold start / missing CH4 restart is interpreted as annavg_finrw=0 and SIF
      ! suppresses first-year saturated CH4 production by mino2lim.
      allocate (annavg_finrw                (numpatch)); annavg_finrw           (:) = spval
      allocate (tempavg_agnpp               (numpatch)); tempavg_agnpp          (:) = 0._r8
      allocate (tempavg_bgnpp               (numpatch)); tempavg_bgnpp          (:) = 0._r8
      allocate (annsum_counter              (numpatch)); annsum_counter         (:) = 0._r8
      allocate (tempavg_somhr               (numpatch)); tempavg_somhr          (:) = 0._r8
      allocate (tempavg_finrw               (numpatch)); tempavg_finrw          (:) = 0._r8

      allocate (fsat_bef                    (numpatch)); fsat_bef               (:) = spval
      allocate (finundated_lag              (numpatch)); finundated_lag         (:) = spval
      allocate (methane_dfsat_tot               (numpatch)); methane_dfsat_tot          (:) = 0._r8
      allocate (f_h2osfc                        (numpatch)); f_h2osfc                   (:) = 0._r8
      allocate (methane_finundated              (numpatch)); methane_finundated         (:) = 0._r8
      allocate (methane_soil_finundated         (numpatch)); methane_soil_finundated    (:) = 0._r8
      allocate (methane_soil_zwt                (numpatch)); methane_soil_zwt           (:) = spval
      allocate (methane_surf_flux_wetland       (numpatch)); methane_surf_flux_wetland  (:) = 0._r8
      allocate (methane_surf_flux_soil          (numpatch)); methane_surf_flux_soil     (:) = 0._r8
      allocate (methane_surf_flux_lake          (numpatch)); methane_surf_flux_lake     (:) = 0._r8
      allocate (methane_surf_flux_rice          (numpatch)); methane_surf_flux_rice     (:) = 0._r8
	      allocate (f_inund_levee_patch             (numpatch)); f_inund_levee_patch        (:) = 0._r8
	      allocate (f_inund_flood_patch             (numpatch)); f_inund_flood_patch        (:) = 0._r8
	      allocate (f_inund_flood_depth_patch       (numpatch)); f_inund_flood_depth_patch  (:) = 0._r8
	      allocate (wetland_frac_per_patch           (numpatch)); wetland_frac_per_patch     (:) = 1._r8
	      allocate (biome_f_methane_patch            (numpatch)); biome_f_methane_patch      (:) = 0.20_r8  ! default = legacy DEF_METHANE%f_methane
	      allocate (biome_redoxlag_patch             (numpatch)); biome_redoxlag_patch       (:) = 30._r8   ! default = legacy DEF_METHANE%redoxlag (days)

	   END SUBROUTINE allocate_methane_state


	   SUBROUTINE init_methane_wetland_fraction_cache (numpatch)
	      ! Build a patch-local lookup for the wetland fraction of each parent
	      ! element/cell.  This is intentionally a static landdata cache:
	      ! routing/GIEMS provide an absolute/grid flood fraction, while methane
	      ! wetland tiles need a patch-relative inundation fraction.
	      USE MOD_LandPatch,           only: landpatch
	      USE MOD_Mesh,                only: numelm, mesh
	      USE MOD_Pixel,               only: pixel
	      USE MOD_Utils,               only: areaquad
	      USE MOD_SPMD_Task,           only: p_is_worker
	      USE MOD_Vars_TimeInvariants, only: patchtype
	      IMPLICIT NONE

	      integer, intent(in) :: numpatch
	      real(r8), allocatable :: elm_wet_area(:), elm_act_area(:)
	      integer :: ipatch, ie
	      real(r8) :: area

	      IF (.not. allocated(wetland_frac_per_patch)) THEN
	         allocate(wetland_frac_per_patch(numpatch))
	      ELSEIF (size(wetland_frac_per_patch) /= numpatch) THEN
	         deallocate(wetland_frac_per_patch)
	         allocate(wetland_frac_per_patch(numpatch))
	      ENDIF
	      wetland_frac_per_patch(:) = 1._r8

	      IF (.not. p_is_worker) RETURN
	      IF (numpatch <= 0 .or. numelm <= 0) RETURN
	      IF (.not. allocated(patchtype)) RETURN
	      IF (.not. allocated(landpatch%ielm)) RETURN

	      allocate(elm_wet_area(numelm), elm_act_area(numelm))
	      elm_wet_area(:) = 0._r8
	      elm_act_area(:) = 0._r8

	      DO ipatch = 1, min(numpatch, size(patchtype))
	         IF (ipatch > size(landpatch%ielm)) CYCLE
	         ie = landpatch%ielm(ipatch)
	         IF (ie < 1 .or. ie > numelm) CYCLE
	         area = methane_patch_area(ipatch)
	         IF (area <= 0._r8) CYCLE
	         IF (patchtype(ipatch) == 0 .or. patchtype(ipatch) == 2) THEN
	            elm_act_area(ie) = elm_act_area(ie) + area
	         ENDIF
	         IF (patchtype(ipatch) == 2) THEN
	            elm_wet_area(ie) = elm_wet_area(ie) + area
	         ENDIF
	      ENDDO

	      DO ipatch = 1, min(numpatch, size(patchtype))
	         IF (ipatch > size(landpatch%ielm)) CYCLE
	         ie = landpatch%ielm(ipatch)
	         IF (ie < 1 .or. ie > numelm) CYCLE
	         IF (elm_act_area(ie) > 0._r8) THEN
	            wetland_frac_per_patch(ipatch) = max(0._r8, min(1._r8, &
	               elm_wet_area(ie) / elm_act_area(ie)))
	         ELSE
	            wetland_frac_per_patch(ipatch) = 0._r8
	         ENDIF
	      ENDDO

	      deallocate(elm_wet_area, elm_act_area)

	   CONTAINS

	      real(r8) FUNCTION methane_patch_area (ip)
	         integer, intent(in) :: ip
	         integer :: ipxstt, ipxend, ipxl, ie_local

	         methane_patch_area = 0._r8
	         IF (ip < 1) RETURN
	         IF (ip > size(landpatch%ielm)) RETURN
	         IF (.not. allocated(landpatch%ipxstt)) RETURN
	         IF (.not. allocated(landpatch%ipxend)) RETURN
	         IF (ip > size(landpatch%ipxstt) .or. ip > size(landpatch%ipxend)) RETURN

	         ie_local = landpatch%ielm(ip)
	         IF (ie_local < 1 .or. ie_local > numelm) RETURN
	         ipxstt = landpatch%ipxstt(ip)
	         ipxend = landpatch%ipxend(ip)
	         IF (ipxstt == -1 .and. ipxend == -1) THEN
	            ipxstt = 1
	            ipxend = mesh(ie_local)%npxl
	         ENDIF
	         IF (ipxstt < 1 .or. ipxend < ipxstt) RETURN
	         IF (ipxend > mesh(ie_local)%npxl) RETURN

	         DO ipxl = ipxstt, ipxend
	            methane_patch_area = methane_patch_area + areaquad ( &
	               pixel%lat_s(mesh(ie_local)%ilat(ipxl)), &
	               pixel%lat_n(mesh(ie_local)%ilat(ipxl)), &
	               pixel%lon_w(mesh(ie_local)%ilon(ipxl)), &
	               pixel%lon_e(mesh(ie_local)%ilon(ipxl)) )
	         ENDDO
	         IF (landpatch%has_shared .and. allocated(landpatch%pctshared)) THEN
	            IF (ip <= size(landpatch%pctshared)) THEN
	               methane_patch_area = methane_patch_area * max(0._r8, landpatch%pctshared(ip))
	            ENDIF
	         ENDIF
	      END FUNCTION methane_patch_area

	   END SUBROUTINE init_methane_wetland_fraction_cache


	   SUBROUTINE deallocate_methane_state ()

      IF (allocated(net_methane)) deallocate (net_methane)
      IF (allocated(methane_prod_depth)) deallocate (methane_prod_depth)
      IF (allocated(o2_decomp_depth)) deallocate (o2_decomp_depth)
      IF (allocated(methane_oxid_depth)) deallocate (methane_oxid_depth)
      IF (allocated(o2_oxid_depth)) deallocate (o2_oxid_depth)
      IF (allocated(methane_aere_depth)) deallocate (methane_aere_depth)
      IF (allocated(methane_tran_depth)) deallocate (methane_tran_depth)
      IF (allocated(o2_aere_depth)) deallocate (o2_aere_depth)
      IF (allocated(methane_ebul_depth)) deallocate (methane_ebul_depth)
      IF (allocated(o2stress)) deallocate (o2stress)
      IF (allocated(methane_stress)) deallocate (methane_stress)
      IF (allocated(methane_surf_flux_tot)) deallocate (methane_surf_flux_tot)
      IF (allocated(methane_surf_flux_tot_phys)) deallocate (methane_surf_flux_tot_phys)
      IF (allocated(methane_surf_aere)) deallocate (methane_surf_aere)
      IF (allocated(methane_surf_ebul)) deallocate (methane_surf_ebul)
      IF (allocated(methane_surf_diff)) deallocate (methane_surf_diff)
      IF (allocated(methane_surf_diff_phys)) deallocate (methane_surf_diff_phys)
      IF (allocated(methane_balance_residual)) deallocate (methane_balance_residual)
      IF (allocated(methane_ch4_clip_credit)) deallocate (methane_ch4_clip_credit)
      IF (allocated(restart_ch4_clip_credit_mass)) deallocate (restart_ch4_clip_credit_mass)
      IF (allocated(o2_cap_loss)) deallocate (o2_cap_loss)
      IF (allocated(o2_cap_gain)) deallocate (o2_cap_gain)
      IF (allocated(methane_ebul_tot)) deallocate (methane_ebul_tot)
      IF (allocated(methane_prod_tot)) deallocate (methane_prod_tot)
      IF (allocated(methane_oxid_tot)) deallocate (methane_oxid_tot)
      IF (allocated(totcol_methane)) deallocate (totcol_methane)
      IF (allocated(grnd_methane_cond)) deallocate (grnd_methane_cond)
      IF (allocated(conc_o2)) deallocate (conc_o2)
      IF (allocated(conc_methane)) deallocate (conc_methane)
      !!!! --------------------------------------------------------------------------------------------------------

      !!!! --------------------------------------------------------------------------------------------------------
      !!!!                                         sum data (unsaturated / saturated)
      !!!! --------------------------------------------------------------------------------------------------------
      IF (allocated(net_methane_unsat)) deallocate (net_methane_unsat)
      IF (allocated(net_methane_sat)) deallocate (net_methane_sat)
      IF (allocated(methane_prod_depth_unsat)) deallocate (methane_prod_depth_unsat)
      IF (allocated(methane_prod_depth_sat)) deallocate (methane_prod_depth_sat)
      IF (allocated(o2_decomp_depth_unsat)) deallocate (o2_decomp_depth_unsat)
      IF (allocated(o2_decomp_depth_sat)) deallocate (o2_decomp_depth_sat)
      IF (allocated(methane_oxid_depth_unsat)) deallocate (methane_oxid_depth_unsat)
      IF (allocated(methane_oxid_depth_sat)) deallocate (methane_oxid_depth_sat)
      IF (allocated(o2_oxid_depth_unsat)) deallocate (o2_oxid_depth_unsat)
      IF (allocated(o2_oxid_depth_sat)) deallocate (o2_oxid_depth_sat)
      IF (allocated(methane_aere_depth_unsat)) deallocate (methane_aere_depth_unsat)
      IF (allocated(methane_aere_depth_sat)) deallocate (methane_aere_depth_sat)
      IF (allocated(methane_tran_depth_unsat)) deallocate (methane_tran_depth_unsat)
      IF (allocated(methane_tran_depth_sat)) deallocate (methane_tran_depth_sat)
      IF (allocated(o2_aere_depth_unsat)) deallocate (o2_aere_depth_unsat)
      IF (allocated(o2_aere_depth_sat)) deallocate (o2_aere_depth_sat)
      IF (allocated(methane_ebul_depth_unsat)) deallocate (methane_ebul_depth_unsat)
      IF (allocated(methane_ebul_depth_sat)) deallocate (methane_ebul_depth_sat)
      IF (allocated(o2stress_unsat)) deallocate (o2stress_unsat)
      IF (allocated(o2stress_sat)) deallocate (o2stress_sat)
      IF (allocated(methane_stress_unsat)) deallocate (methane_stress_unsat)
      IF (allocated(methane_stress_sat)) deallocate (methane_stress_sat)
      IF (allocated(methane_surf_flux_tot_unsat)) deallocate (methane_surf_flux_tot_unsat)
      IF (allocated(methane_surf_flux_tot_sat)) deallocate (methane_surf_flux_tot_sat)
      IF (allocated(methane_surf_aere_unsat)) deallocate (methane_surf_aere_unsat)
      IF (allocated(methane_surf_aere_sat)) deallocate (methane_surf_aere_sat)
      IF (allocated(methane_surf_ebul_unsat)) deallocate (methane_surf_ebul_unsat)
      IF (allocated(methane_surf_ebul_sat)) deallocate (methane_surf_ebul_sat)
      IF (allocated(methane_surf_diff_unsat)) deallocate (methane_surf_diff_unsat)
      IF (allocated(methane_surf_diff_sat)) deallocate (methane_surf_diff_sat)
      IF (allocated(methane_surf_diff_phys_unsat)) deallocate (methane_surf_diff_phys_unsat)
      IF (allocated(methane_surf_diff_phys_sat)) deallocate (methane_surf_diff_phys_sat)
      IF (allocated(methane_ebul_tot_unsat)) deallocate (methane_ebul_tot_unsat)
      IF (allocated(methane_ebul_tot_sat)) deallocate (methane_ebul_tot_sat)
      IF (allocated(methane_prod_tot_unsat)) deallocate (methane_prod_tot_unsat)
      IF (allocated(methane_prod_tot_sat)) deallocate (methane_prod_tot_sat)
      IF (allocated(methane_oxid_tot_unsat)) deallocate (methane_oxid_tot_unsat)
      IF (allocated(methane_oxid_tot_sat)) deallocate (methane_oxid_tot_sat)
      IF (allocated(totcol_methane_unsat)) deallocate (totcol_methane_unsat)
      IF (allocated(totcol_methane_sat)) deallocate (totcol_methane_sat)
      IF (allocated(grnd_methane_cond_unsat)) deallocate (grnd_methane_cond_unsat)
      IF (allocated(grnd_methane_cond_sat)) deallocate (grnd_methane_cond_sat)
      IF (allocated(conc_o2_unsat)) deallocate (conc_o2_unsat)
	      IF (allocated(conc_o2_sat)) deallocate (conc_o2_sat)
	      IF (allocated(conc_methane_unsat)) deallocate (conc_methane_unsat)
	      IF (allocated(conc_methane_sat)) deallocate (conc_methane_sat)
	      !!!! --------------------------------------------------------------------------------------------------------

	      !!!! --------------------------------------------------------------------------------------------------------
	      !!!!                                         lake data (CTSM alignment)
	      !!!! --------------------------------------------------------------------------------------------------------
	      IF (allocated(methane_prod_depth_lake)) deallocate (methane_prod_depth_lake)
	      IF (allocated(methane_oxid_depth_lake)) deallocate (methane_oxid_depth_lake)
	      IF (allocated(methane_ebul_depth_lake)) deallocate (methane_ebul_depth_lake)
	      IF (allocated(methane_surf_ebul_lake)) deallocate (methane_surf_ebul_lake)
	      IF (allocated(methane_surf_diff_lake)) deallocate (methane_surf_diff_lake)
	      IF (allocated(methane_surf_flux_tot_lake)) deallocate (methane_surf_flux_tot_lake)
	      IF (allocated(methane_prod_tot_lake)) deallocate (methane_prod_tot_lake)
	      IF (allocated(methane_oxid_tot_lake)) deallocate (methane_oxid_tot_lake)
	      IF (allocated(methane_ebul_tot_lake)) deallocate (methane_ebul_tot_lake)
	      IF (allocated(totcol_methane_lake)) deallocate (totcol_methane_lake)
	      IF (allocated(grnd_methane_cond_lake)) deallocate (grnd_methane_cond_lake)
	      IF (allocated(conc_o2_lake)) deallocate (conc_o2_lake)
	      IF (allocated(conc_methane_lake)) deallocate (conc_methane_lake)
	      !!!! --------------------------------------------------------------------------------------------------------

	      IF (allocated(c_atm)) deallocate (c_atm)
      IF (allocated(forc_pmethanem)) deallocate (forc_pmethanem)
      IF (allocated(layer_sat_lag)) deallocate (layer_sat_lag)
      IF (allocated(lake_soilc)) deallocate (lake_soilc)
      IF (allocated(annavg_agnpp)) deallocate (annavg_agnpp)
      IF (allocated(annavg_bgnpp)) deallocate (annavg_bgnpp)
      IF (allocated(annavg_somhr)) deallocate (annavg_somhr)
      IF (allocated(annavg_finrw)) deallocate (annavg_finrw)
      IF (allocated(tempavg_agnpp)) deallocate (tempavg_agnpp)
      IF (allocated(tempavg_bgnpp)) deallocate (tempavg_bgnpp)
      IF (allocated(annsum_counter)) deallocate (annsum_counter)
      IF (allocated(tempavg_somhr)) deallocate (tempavg_somhr)
      IF (allocated(tempavg_finrw)) deallocate (tempavg_finrw)
      IF (allocated(fsat_bef)) deallocate (fsat_bef)
      IF (allocated(finundated_lag)) deallocate (finundated_lag)
      IF (allocated(methane_dfsat_tot)) deallocate (methane_dfsat_tot)
	      IF (allocated(f_h2osfc)) deallocate (f_h2osfc)
	      IF (allocated(methane_finundated)) deallocate (methane_finundated)
	      IF (allocated(methane_soil_finundated)) deallocate (methane_soil_finundated)
	      IF (allocated(methane_soil_zwt)) deallocate (methane_soil_zwt)
	      IF (allocated(methane_surf_flux_wetland)) deallocate (methane_surf_flux_wetland)
	      IF (allocated(methane_surf_flux_soil)) deallocate (methane_surf_flux_soil)
	      IF (allocated(methane_surf_flux_lake)) deallocate (methane_surf_flux_lake)
	      IF (allocated(methane_surf_flux_rice)) deallocate (methane_surf_flux_rice)
	      IF (allocated(f_inund_levee_patch)) deallocate (f_inund_levee_patch)
	      IF (allocated(f_inund_flood_patch)) deallocate (f_inund_flood_patch)
	      IF (allocated(f_inund_flood_depth_patch)) deallocate (f_inund_flood_depth_patch)
	      IF (allocated(wetland_frac_per_patch)) deallocate (wetland_frac_per_patch)
	      IF (allocated(biome_f_methane_patch)) deallocate (biome_f_methane_patch)
	      IF (allocated(biome_redoxlag_patch)) deallocate (biome_redoxlag_patch)
      IF (allocated(co2_decomp_depth       )) deallocate (co2_decomp_depth       )
      IF (allocated(co2_oxid_depth         )) deallocate (co2_oxid_depth         )
      IF (allocated(co2_aere_depth         )) deallocate (co2_aere_depth         )
      IF (allocated(co2_decomp_tot         )) deallocate (co2_decomp_tot         )
      IF (allocated(co2_oxid_tot           )) deallocate (co2_oxid_tot           )
      IF (allocated(co2_aere_tot           )) deallocate (co2_aere_tot           )
      IF (allocated(co2_net_tot            )) deallocate (co2_net_tot            )
      IF (allocated(co2_decomp_depth_unsat )) deallocate (co2_decomp_depth_unsat )
      IF (allocated(co2_decomp_depth_sat   )) deallocate (co2_decomp_depth_sat   )
      IF (allocated(co2_oxid_depth_unsat   )) deallocate (co2_oxid_depth_unsat   )
      IF (allocated(co2_oxid_depth_sat     )) deallocate (co2_oxid_depth_sat     )
      IF (allocated(co2_aere_depth_unsat   )) deallocate (co2_aere_depth_unsat   )
      IF (allocated(co2_aere_depth_sat     )) deallocate (co2_aere_depth_sat     )
      IF (allocated(co2_decomp_tot_unsat   )) deallocate (co2_decomp_tot_unsat   )
      IF (allocated(co2_decomp_tot_sat     )) deallocate (co2_decomp_tot_sat     )
      IF (allocated(co2_oxid_tot_unsat     )) deallocate (co2_oxid_tot_unsat     )
      IF (allocated(co2_oxid_tot_sat       )) deallocate (co2_oxid_tot_sat       )
      IF (allocated(co2_net_tot_unsat      )) deallocate (co2_net_tot_unsat      )
      IF (allocated(co2_net_tot_sat        )) deallocate (co2_net_tot_sat        )
      IF (allocated(co2_decomp_depth_lake  )) deallocate (co2_decomp_depth_lake  )
      IF (allocated(co2_oxid_depth_lake    )) deallocate (co2_oxid_depth_lake    )
      IF (allocated(co2_decomp_tot_lake    )) deallocate (co2_decomp_tot_lake    )
	      IF (allocated(co2_oxid_tot_lake      )) deallocate (co2_oxid_tot_lake      )
	      IF (allocated(co2_net_tot_lake       )) deallocate (co2_net_tot_lake       )
	      IF (allocated(methane_lake_substep_acc2d)) deallocate (methane_lake_substep_acc2d)
	      IF (allocated(methane_lake_substep_acc1d)) deallocate (methane_lake_substep_acc1d)
	      methane_lake_substep_cached_ipatch = -1
	      methane_lake_substep_next_isub = 1

	   END SUBROUTINE deallocate_methane_state


	   SUBROUTINE accumulate_methane_lake_substep_diagnostics (ipatch, substep_dt, isub, nsub)
	      integer,  intent(in) :: ipatch
	      real(r8), intent(in) :: substep_dt
	      integer,  intent(in) :: isub
	      integer,  intent(in) :: nsub

	      real(r8) :: total_dt

	      IF (nsub <= 1) RETURN
	      IF (substep_dt <= 0._r8) RETURN
	      IF (ipatch <= 0) RETURN

	      IF (.not. allocated(methane_lake_substep_acc2d)) THEN
	         allocate (methane_lake_substep_acc2d(nl_soil, methane_lake_substep_n2d))
	      ENDIF
	      IF (.not. allocated(methane_lake_substep_acc1d)) THEN
	         allocate (methane_lake_substep_acc1d(methane_lake_substep_n1d))
	      ENDIF

	      IF (isub <= 1 .or. ipatch /= methane_lake_substep_cached_ipatch .or. &
	          isub /= methane_lake_substep_next_isub) THEN
	         methane_lake_substep_acc2d(:,:) = 0._r8
	         methane_lake_substep_acc1d(:)   = 0._r8
	      ENDIF
	      methane_lake_substep_cached_ipatch = ipatch

	      CALL add2d( 1, methane_prod_depth(:,ipatch))
	      CALL add2d( 2, o2_decomp_depth(:,ipatch))
	      CALL add2d( 3, co2_decomp_depth(:,ipatch))
	      CALL add2d( 4, methane_oxid_depth(:,ipatch))
	      CALL add2d( 5, o2_oxid_depth(:,ipatch))
	      CALL add2d( 6, co2_oxid_depth(:,ipatch))
	      CALL add2d( 7, methane_aere_depth(:,ipatch))
	      CALL add2d( 8, methane_tran_depth(:,ipatch))
	      CALL add2d( 9, o2_aere_depth(:,ipatch))
	      CALL add2d(10, co2_aere_depth(:,ipatch))
	      CALL add2d(11, methane_ebul_depth(:,ipatch))
	      CALL add2d(12, o2stress(:,ipatch))
	      CALL add2d(13, methane_stress(:,ipatch))
	      CALL add2d(14, methane_prod_depth_unsat(:,ipatch))
	      CALL add2d(15, o2_decomp_depth_unsat(:,ipatch))
	      CALL add2d(16, co2_decomp_depth_unsat(:,ipatch))
	      CALL add2d(17, methane_oxid_depth_unsat(:,ipatch))
	      CALL add2d(18, o2_oxid_depth_unsat(:,ipatch))
	      CALL add2d(19, co2_oxid_depth_unsat(:,ipatch))
	      CALL add2d(20, methane_aere_depth_unsat(:,ipatch))
	      CALL add2d(21, methane_tran_depth_unsat(:,ipatch))
	      CALL add2d(22, o2_aere_depth_unsat(:,ipatch))
	      CALL add2d(23, co2_aere_depth_unsat(:,ipatch))
	      CALL add2d(24, methane_ebul_depth_unsat(:,ipatch))
	      CALL add2d(25, o2stress_unsat(:,ipatch))
	      CALL add2d(26, methane_stress_unsat(:,ipatch))
	      CALL add2d(27, methane_prod_depth_sat(:,ipatch))
	      CALL add2d(28, o2_decomp_depth_sat(:,ipatch))
	      CALL add2d(29, co2_decomp_depth_sat(:,ipatch))
	      CALL add2d(30, methane_oxid_depth_sat(:,ipatch))
	      CALL add2d(31, o2_oxid_depth_sat(:,ipatch))
	      CALL add2d(32, co2_oxid_depth_sat(:,ipatch))
	      CALL add2d(33, methane_aere_depth_sat(:,ipatch))
	      CALL add2d(34, methane_tran_depth_sat(:,ipatch))
	      CALL add2d(35, o2_aere_depth_sat(:,ipatch))
	      CALL add2d(36, co2_aere_depth_sat(:,ipatch))
	      CALL add2d(37, methane_ebul_depth_sat(:,ipatch))
	      CALL add2d(38, o2stress_sat(:,ipatch))
	      CALL add2d(39, methane_stress_sat(:,ipatch))
	      CALL add2d(40, methane_prod_depth_lake(:,ipatch))
	      CALL add2d(41, methane_oxid_depth_lake(:,ipatch))
	      CALL add2d(42, methane_ebul_depth_lake(:,ipatch))
	      CALL add2d(43, co2_decomp_depth_lake(:,ipatch))
	      CALL add2d(44, co2_oxid_depth_lake(:,ipatch))

	      CALL add1d( 1, net_methane(ipatch))
	      CALL add1d( 2, methane_surf_flux_tot(ipatch))
	      CALL add1d(46, methane_surf_flux_tot_phys(ipatch))
	      CALL add1d( 3, methane_surf_aere(ipatch))
	      CALL add1d( 4, methane_surf_ebul(ipatch))
	      CALL add1d( 5, methane_surf_diff(ipatch))
	      CALL add1d(43, methane_balance_residual(ipatch))
	      CALL add1d(47, methane_ch4_clip_credit(ipatch))
	      CALL add1d(44, o2_cap_loss(ipatch))
	      CALL add1d(45, o2_cap_gain(ipatch))
	      CALL add1d( 6, methane_ebul_tot(ipatch))
	      CALL add1d( 7, methane_prod_tot(ipatch))
	      CALL add1d( 8, methane_oxid_tot(ipatch))
	      CALL add1d( 9, co2_decomp_tot(ipatch))
	      CALL add1d(10, co2_oxid_tot(ipatch))
	      CALL add1d(11, co2_aere_tot(ipatch))
	      CALL add1d(12, co2_net_tot(ipatch))
	      CALL add1d(13, net_methane_unsat(ipatch))
	      CALL add1d(14, net_methane_sat(ipatch))
	      CALL add1d(15, methane_surf_flux_tot_unsat(ipatch))
	      CALL add1d(16, methane_surf_flux_tot_sat(ipatch))
	      CALL add1d(17, methane_surf_aere_unsat(ipatch))
	      CALL add1d(18, methane_surf_aere_sat(ipatch))
	      CALL add1d(19, methane_surf_ebul_unsat(ipatch))
	      CALL add1d(20, methane_surf_ebul_sat(ipatch))
	      CALL add1d(21, methane_surf_diff_unsat(ipatch))
	      CALL add1d(22, methane_surf_diff_sat(ipatch))
	      CALL add1d(23, methane_ebul_tot_unsat(ipatch))
	      CALL add1d(24, methane_ebul_tot_sat(ipatch))
	      CALL add1d(25, methane_prod_tot_unsat(ipatch))
	      CALL add1d(26, methane_prod_tot_sat(ipatch))
	      CALL add1d(27, methane_oxid_tot_unsat(ipatch))
	      CALL add1d(28, methane_oxid_tot_sat(ipatch))
	      CALL add1d(29, co2_decomp_tot_unsat(ipatch))
	      CALL add1d(30, co2_decomp_tot_sat(ipatch))
	      CALL add1d(31, co2_oxid_tot_unsat(ipatch))
	      CALL add1d(32, co2_oxid_tot_sat(ipatch))
	      CALL add1d(33, co2_net_tot_unsat(ipatch))
	      CALL add1d(34, co2_net_tot_sat(ipatch))
	      CALL add1d(35, methane_surf_ebul_lake(ipatch))
	      CALL add1d(36, methane_surf_diff_lake(ipatch))
	      CALL add1d(48, methane_surf_flux_tot_lake(ipatch))
	      CALL add1d(49, methane_surf_flux_lake(ipatch))
	      CALL add1d(37, methane_prod_tot_lake(ipatch))
	      CALL add1d(38, methane_oxid_tot_lake(ipatch))
	      CALL add1d(39, methane_ebul_tot_lake(ipatch))
	      CALL add1d(40, co2_decomp_tot_lake(ipatch))
	      CALL add1d(41, co2_oxid_tot_lake(ipatch))
	      CALL add1d(42, co2_net_tot_lake(ipatch))

	      IF (isub == nsub) THEN
	         total_dt = substep_dt * real(nsub, r8)
	         IF (total_dt > 0._r8) THEN
	            CALL finish2d( 1, methane_prod_depth(:,ipatch))
	            CALL finish2d( 2, o2_decomp_depth(:,ipatch))
	            CALL finish2d( 3, co2_decomp_depth(:,ipatch))
	            CALL finish2d( 4, methane_oxid_depth(:,ipatch))
	            CALL finish2d( 5, o2_oxid_depth(:,ipatch))
	            CALL finish2d( 6, co2_oxid_depth(:,ipatch))
	            CALL finish2d( 7, methane_aere_depth(:,ipatch))
	            CALL finish2d( 8, methane_tran_depth(:,ipatch))
	            CALL finish2d( 9, o2_aere_depth(:,ipatch))
	            CALL finish2d(10, co2_aere_depth(:,ipatch))
	            CALL finish2d(11, methane_ebul_depth(:,ipatch))
	            CALL finish2d(12, o2stress(:,ipatch))
	            CALL finish2d(13, methane_stress(:,ipatch))
	            CALL finish2d(14, methane_prod_depth_unsat(:,ipatch))
	            CALL finish2d(15, o2_decomp_depth_unsat(:,ipatch))
	            CALL finish2d(16, co2_decomp_depth_unsat(:,ipatch))
	            CALL finish2d(17, methane_oxid_depth_unsat(:,ipatch))
	            CALL finish2d(18, o2_oxid_depth_unsat(:,ipatch))
	            CALL finish2d(19, co2_oxid_depth_unsat(:,ipatch))
	            CALL finish2d(20, methane_aere_depth_unsat(:,ipatch))
	            CALL finish2d(21, methane_tran_depth_unsat(:,ipatch))
	            CALL finish2d(22, o2_aere_depth_unsat(:,ipatch))
	            CALL finish2d(23, co2_aere_depth_unsat(:,ipatch))
	            CALL finish2d(24, methane_ebul_depth_unsat(:,ipatch))
	            CALL finish2d(25, o2stress_unsat(:,ipatch))
	            CALL finish2d(26, methane_stress_unsat(:,ipatch))
	            CALL finish2d(27, methane_prod_depth_sat(:,ipatch))
	            CALL finish2d(28, o2_decomp_depth_sat(:,ipatch))
	            CALL finish2d(29, co2_decomp_depth_sat(:,ipatch))
	            CALL finish2d(30, methane_oxid_depth_sat(:,ipatch))
	            CALL finish2d(31, o2_oxid_depth_sat(:,ipatch))
	            CALL finish2d(32, co2_oxid_depth_sat(:,ipatch))
	            CALL finish2d(33, methane_aere_depth_sat(:,ipatch))
	            CALL finish2d(34, methane_tran_depth_sat(:,ipatch))
	            CALL finish2d(35, o2_aere_depth_sat(:,ipatch))
	            CALL finish2d(36, co2_aere_depth_sat(:,ipatch))
	            CALL finish2d(37, methane_ebul_depth_sat(:,ipatch))
	            CALL finish2d(38, o2stress_sat(:,ipatch))
	            CALL finish2d(39, methane_stress_sat(:,ipatch))
	            CALL finish2d(40, methane_prod_depth_lake(:,ipatch))
	            CALL finish2d(41, methane_oxid_depth_lake(:,ipatch))
	            CALL finish2d(42, methane_ebul_depth_lake(:,ipatch))
	            CALL finish2d(43, co2_decomp_depth_lake(:,ipatch))
	            CALL finish2d(44, co2_oxid_depth_lake(:,ipatch))

	            CALL finish1d( 1, net_methane(ipatch))
	            CALL finish1d( 2, methane_surf_flux_tot(ipatch))
	            CALL finish1d(46, methane_surf_flux_tot_phys(ipatch))
	            CALL finish1d( 3, methane_surf_aere(ipatch))
	            CALL finish1d( 4, methane_surf_ebul(ipatch))
	            CALL finish1d( 5, methane_surf_diff(ipatch))
	            CALL finish1d(43, methane_balance_residual(ipatch))
	            CALL finish1d(47, methane_ch4_clip_credit(ipatch))
	            CALL finish1d(44, o2_cap_loss(ipatch))
	            CALL finish1d(45, o2_cap_gain(ipatch))
	            CALL finish1d( 6, methane_ebul_tot(ipatch))
	            CALL finish1d( 7, methane_prod_tot(ipatch))
	            CALL finish1d( 8, methane_oxid_tot(ipatch))
	            CALL finish1d( 9, co2_decomp_tot(ipatch))
	            CALL finish1d(10, co2_oxid_tot(ipatch))
	            CALL finish1d(11, co2_aere_tot(ipatch))
	            CALL finish1d(12, co2_net_tot(ipatch))
	            CALL finish1d(13, net_methane_unsat(ipatch))
	            CALL finish1d(14, net_methane_sat(ipatch))
	            CALL finish1d(15, methane_surf_flux_tot_unsat(ipatch))
	            CALL finish1d(16, methane_surf_flux_tot_sat(ipatch))
	            CALL finish1d(17, methane_surf_aere_unsat(ipatch))
	            CALL finish1d(18, methane_surf_aere_sat(ipatch))
	            CALL finish1d(19, methane_surf_ebul_unsat(ipatch))
	            CALL finish1d(20, methane_surf_ebul_sat(ipatch))
	            CALL finish1d(21, methane_surf_diff_unsat(ipatch))
	            CALL finish1d(22, methane_surf_diff_sat(ipatch))
	            CALL finish1d(23, methane_ebul_tot_unsat(ipatch))
	            CALL finish1d(24, methane_ebul_tot_sat(ipatch))
	            CALL finish1d(25, methane_prod_tot_unsat(ipatch))
	            CALL finish1d(26, methane_prod_tot_sat(ipatch))
	            CALL finish1d(27, methane_oxid_tot_unsat(ipatch))
	            CALL finish1d(28, methane_oxid_tot_sat(ipatch))
	            CALL finish1d(29, co2_decomp_tot_unsat(ipatch))
	            CALL finish1d(30, co2_decomp_tot_sat(ipatch))
	            CALL finish1d(31, co2_oxid_tot_unsat(ipatch))
	            CALL finish1d(32, co2_oxid_tot_sat(ipatch))
	            CALL finish1d(33, co2_net_tot_unsat(ipatch))
	            CALL finish1d(34, co2_net_tot_sat(ipatch))
	            CALL finish1d(35, methane_surf_ebul_lake(ipatch))
	            CALL finish1d(36, methane_surf_diff_lake(ipatch))
	            CALL finish1d(48, methane_surf_flux_tot_lake(ipatch))
	            CALL finish1d(49, methane_surf_flux_lake(ipatch))
	            CALL finish1d(37, methane_prod_tot_lake(ipatch))
	            CALL finish1d(38, methane_oxid_tot_lake(ipatch))
	            CALL finish1d(39, methane_ebul_tot_lake(ipatch))
	            CALL finish1d(40, co2_decomp_tot_lake(ipatch))
	            CALL finish1d(41, co2_oxid_tot_lake(ipatch))
	            CALL finish1d(42, co2_net_tot_lake(ipatch))
	         ENDIF
	         methane_lake_substep_cached_ipatch = -1
	         methane_lake_substep_next_isub = 1
	      ELSE
	         methane_lake_substep_next_isub = isub + 1
	      ENDIF

	   CONTAINS
	      SUBROUTINE add2d (icol, var)
	         integer,  intent(in) :: icol
	         real(r8), intent(in) :: var(1:nl_soil)
	         methane_lake_substep_acc2d(:,icol) = methane_lake_substep_acc2d(:,icol) + var(:) * substep_dt
	      END SUBROUTINE add2d

	      SUBROUTINE add1d (icol, var)
	         integer,  intent(in) :: icol
	         real(r8), intent(in) :: var
	         methane_lake_substep_acc1d(icol) = methane_lake_substep_acc1d(icol) + var * substep_dt
	      END SUBROUTINE add1d

	      SUBROUTINE finish2d (icol, var)
	         integer,  intent(in)    :: icol
	         real(r8), intent(inout) :: var(1:nl_soil)
	         var(:) = methane_lake_substep_acc2d(:,icol) / total_dt
	      END SUBROUTINE finish2d

	      SUBROUTINE finish1d (icol, var)
	         integer,  intent(in)    :: icol
	         real(r8), intent(inout) :: var
	         var = methane_lake_substep_acc1d(icol) / total_dt
	      END SUBROUTINE finish1d
	   END SUBROUTINE accumulate_methane_lake_substep_diagnostics


	   !-------------------------------------------------------------------
	   ! compute_f_h2osfc — CLM5 microtopography-based prognostic surface water fraction
   !
   ! Ported from the Methane source restart flow.
   !   main/MOD_Runoff.F90 (line 332-462) / MOD_Runoff_h2osfc.F90 (line 247)
   !   main/MOD_SoilSnowHydrology.F90 (line 1307-1322)
   ! Maintains f_h2osfc(i) as patch-level diagnostic state, called by CoLMDRIVER
   ! immediately before methane_driver. Implementation is minimal-invasion:
   ! does NOT modify CoLM202X WATER_2014/WATER_VSF/Runoff physics. f_h2osfc
   ! is recomputed each step from current wdsrf + slpratio + microtopography
   ! params (DEF_METHANE_hydrology%slopemax/slopebeta/pc).
   !
   ! Inputs:
   !   ipatch    — patch index
   !   slpratio  — slope ratio (from MOD_Vars_TimeInvariants)
   !   wdsrf     — surface water depth [mm] (from MOD_Vars_TimeVariables)
   ! Output:
   !   f_h2osfc(ipatch) is updated
   !-------------------------------------------------------------------
   SUBROUTINE compute_f_h2osfc (ipatch, slpratio_in, wdsrf_in)
      USE MOD_Vars_Global, only: PI
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE_hydrology
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: slpratio_in      ! slope ratio [-]
      real(r8), intent(in) :: wdsrf_in         ! surface water depth [mm]

      real(r8) :: micro_sigma, sigma_mm, d, slope_arg, slope_angle
      real(r8) :: fd, dfdd, d_lo, d_hi, d_mid, f_mid
      integer  :: p
      logical  :: converged
      real(r8), parameter :: pondmin = 1.e-8_r8
      real(r8), parameter :: fd_tol  = 1.e-10_r8

      IF (wdsrf_in <= pondmin .or. abs(slpratio_in) >= 1.e30_r8 .or. &
          DEF_METHANE_hydrology%slopemax <= 0._r8 .or. &
          DEF_METHANE_hydrology%slopebeta == 0._r8) THEN
         f_h2osfc(ipatch) = 0._r8
         RETURN
      END IF

      ! micro_sigma in [m], parametrised from a non-negative slope angle.
      ! CoLM landdata usually stores slpratio as tan(slope), while some raw
      ! topography products expose slope as radians, degrees, or percent.  Use
      ! the historical tan(slope) path for normal small values, but normalize
      ! obvious degree/percent encodings before applying the CLM expression.
      slope_arg = max(slpratio_in, 0._r8)
      IF (slope_arg <= 1._r8) THEN
         slope_angle = atan(slope_arg)
      ELSEIF (slope_arg <= 0.5_r8*PI) THEN
         slope_angle = slope_arg
      ELSEIF (slope_arg <= 90._r8) THEN
         slope_angle = slope_arg * PI / 180._r8
      ELSE
         slope_angle = atan(slope_arg / 100._r8)
      ENDIF
      slope_angle = max(0._r8, min(0.5_r8*PI, slope_angle))
      micro_sigma = (slope_angle + &
                     DEF_METHANE_hydrology%slopemax**(1._r8/DEF_METHANE_hydrology%slopebeta) &
                    )**DEF_METHANE_hydrology%slopebeta
      micro_sigma = max(0._r8, min(DEF_METHANE_hydrology%slopemax, micro_sigma))
      sigma_mm = 1.0e3_r8 * micro_sigma   ! m -> mm

	      IF (sigma_mm > 1.e-3_r8) THEN
	         IF (wdsrf_in >= 10._r8 * sigma_mm) THEN
	            f_h2osfc(ipatch) = 1._r8
	            RETURN
	         ENDIF
	         ! Newton iteration for the CLM fill-and-spill microtopography
         ! relation:
         !   W(d) = 0.5*d*(1+erf(d/(sigma*sqrt(2))))
         !        + sigma/sqrt(2*pi)*exp(-d**2/(2*sigma**2))
         ! where W is grid-cell mean surface-water depth [mm].  The old
         ! implementation solved an unrelated pc threshold and then replaced
         ! d by wdsrf, which made wdsrf=0 diagnose f_h2osfc=0.5.
		         d = min(max(0._r8, wdsrf_in), 10._r8 * sigma_mm)
		         converged = .false.
	         DO p = 1, 20
	            fd = 0.5_r8 * d * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8)))) &
	               + sigma_mm / sqrt(2.0_r8 * PI) &
	               * exp(-d**2 / (2.0_r8 * sigma_mm**2)) &
	               - wdsrf_in
	            dfdd = 0.5_r8 * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8))))
	            IF (abs(fd) < fd_tol) THEN
	               converged = .true.
	               EXIT
	            ENDIF
	            IF (dfdd < 1.e-12_r8) EXIT
	            d = d - fd / dfdd
	            d = max(-10._r8 * sigma_mm, min(10._r8 * sigma_mm, d))
	         END DO
	         IF (.not. converged) THEN
	            ! Monotone fallback for rare Newton failures.  This avoids
	            ! silently using the last Newton iterate when the derivative is
	            ! tiny near the dry tail.
	            d_lo = -10._r8 * sigma_mm
	            d_hi =  10._r8 * sigma_mm
	            DO p = 1, 60
	               d_mid = 0.5_r8 * (d_lo + d_hi)
	               f_mid = 0.5_r8 * d_mid * (1.0_r8 + erf(d_mid / (sigma_mm * sqrt(2.0_r8)))) &
	                  + sigma_mm / sqrt(2.0_r8 * PI) &
	                  * exp(-d_mid**2 / (2.0_r8 * sigma_mm**2)) &
	                  - wdsrf_in
	               IF (abs(f_mid) < fd_tol) EXIT
	               IF (f_mid > 0._r8) THEN
	                  d_hi = d_mid
	               ELSE
	                  d_lo = d_mid
	               ENDIF
	            END DO
	            d = d_mid
	         ENDIF
         f_h2osfc(ipatch) = 0.5_r8 * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8))))
         f_h2osfc(ipatch) = min(1._r8, max(0._r8, f_h2osfc(ipatch)))
      ELSE
         f_h2osfc(ipatch) = 0._r8
      END IF
   END SUBROUTINE compute_f_h2osfc


   SUBROUTINE write_methane_restart (file_restart, compress)
      USE MOD_LandPatch,     only: landpatch
      USE MOD_NetCDFVector,  only: ncio_write_vector
      character(len=*), intent(in) :: file_restart
      integer,          intent(in) :: compress

      IF (.not. allocated(conc_methane)) RETURN

	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2',          'soil', nl_soil, 'patch', landpatch, conc_o2,          compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_methane',     'soil', nl_soil, 'patch', landpatch, conc_methane,     compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_methane',   'patch', landpatch, totcol_methane,       compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond','patch', landpatch, grnd_methane_cond,    compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_unsat',    'soil', nl_soil, 'patch', landpatch, conc_o2_unsat,    compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_sat',      'soil', nl_soil, 'patch', landpatch, conc_o2_sat,      compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_unsat',   'soil', nl_soil, &
	                              'patch', landpatch, conc_methane_unsat, compress)
		      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_sat',     'soil', nl_soil, &
		                              'patch', landpatch, conc_methane_sat,   compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_methane_unsat','patch', landpatch, totcol_methane_unsat, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_methane_sat',  'patch', landpatch, totcol_methane_sat,   compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond_unsat','patch', landpatch, grnd_methane_cond_unsat, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond_sat',  'patch', landpatch, grnd_methane_cond_sat,   compress)
		      CALL ncio_write_vector (file_restart, 'ch4_conc_o2_lake',     'soil', nl_soil, 'patch', landpatch, conc_o2_lake,     compress)
	      CALL ncio_write_vector (file_restart, 'ch4_conc_ch4_lake',    'soil', nl_soil, 'patch', landpatch, conc_methane_lake, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_totcol_lake',      'patch', landpatch, totcol_methane_lake, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_grnd_methane_cond_lake','patch', landpatch, grnd_methane_cond_lake, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_layer_sat_lag',    'soil', nl_soil, 'patch', landpatch, layer_sat_lag,    compress)
	      CALL ncio_write_vector (file_restart, 'ch4_lake_soilc',       'soil', nl_soil, 'patch', landpatch, lake_soilc,       compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_agnpp',     'patch', landpatch, annavg_agnpp,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_bgnpp',     'patch', landpatch, annavg_bgnpp,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_somhr',     'patch', landpatch, annavg_somhr,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_annavg_finrw',     'patch', landpatch, annavg_finrw,     compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_agnpp',    'patch', landpatch, tempavg_agnpp,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_bgnpp',    'patch', landpatch, tempavg_bgnpp,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_annsum_counter',   'patch', landpatch, annsum_counter,   compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_somhr',    'patch', landpatch, tempavg_somhr,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_tempavg_finrw',    'patch', landpatch, tempavg_finrw,    compress)
      CALL ncio_write_vector (file_restart, 'ch4_fsat_bef',         'patch', landpatch, fsat_bef,         compress)
	      CALL ncio_write_vector (file_restart, 'ch4_finundated_lag',   'patch', landpatch, finundated_lag,   compress)
	      CALL ncio_write_vector (file_restart, 'ch4_methane_dfsat_tot','patch', landpatch, methane_dfsat_tot, compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_h2osfc',         'patch', landpatch, f_h2osfc,         compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_inund_levee_patch',       &
	                              'patch', landpatch, f_inund_levee_patch,       compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_inund_flood_patch',       &
	                              'patch', landpatch, f_inund_flood_patch,       compress)
	      CALL ncio_write_vector (file_restart, 'ch4_f_inund_flood_depth_patch', &
	                              'patch', landpatch, f_inund_flood_depth_patch, compress)
      IF (allocated(restart_ch4_clip_credit_mass)) CALL ncio_write_vector (file_restart, &
         'ch4_restart_ch4_clip_credit_mass', 'patch', landpatch, restart_ch4_clip_credit_mass, compress)
	   END SUBROUTINE write_methane_restart


   SUBROUTINE read_methane_restart (file_restart)
      USE MOD_LandPatch,     only: landpatch
	      USE MOD_NetCDFVector,  only: ncio_read_vector
	      character(len=*), intent(in) :: file_restart

	      IF (.not. allocated(conc_methane)) RETURN

	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2',          nl_soil, landpatch, conc_o2,          defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_methane',     nl_soil, landpatch, conc_methane,     defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_methane',   landpatch, totcol_methane,            defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond',landpatch, grnd_methane_cond,         defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2_unsat',    nl_soil, landpatch, conc_o2_unsat,    defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2_sat',      nl_soil, landpatch, conc_o2_sat,      defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_ch4_unsat',   nl_soil, landpatch, conc_methane_unsat, defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_ch4_sat',     nl_soil, landpatch, conc_methane_sat,   defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_methane_unsat', landpatch, totcol_methane_unsat,    defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_methane_sat',   landpatch, totcol_methane_sat,      defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond_unsat', landpatch, grnd_methane_cond_unsat, defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond_sat',   landpatch, grnd_methane_cond_sat,   defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_o2_lake',     nl_soil, landpatch, conc_o2_lake,       defval = 1._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_conc_ch4_lake',    nl_soil, landpatch, conc_methane_lake,  defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_totcol_lake',      landpatch, totcol_methane_lake,         defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_grnd_methane_cond_lake', landpatch, grnd_methane_cond_lake, defval = 1.e-6_r8)
	      CALL ncio_read_vector (file_restart, 'ch4_layer_sat_lag',    nl_soil, landpatch, layer_sat_lag,    defval = spval)
      CALL ncio_read_vector (file_restart, 'ch4_lake_soilc',       nl_soil, landpatch, lake_soilc,       defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_agnpp',     landpatch, annavg_agnpp,     defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_bgnpp',     landpatch, annavg_bgnpp,     defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_somhr',     landpatch, annavg_somhr,     defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annavg_finrw',     landpatch, annavg_finrw,     defval = spval)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_agnpp',    landpatch, tempavg_agnpp,    defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_bgnpp',    landpatch, tempavg_bgnpp,    defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_annsum_counter',   landpatch, annsum_counter,   defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_somhr',    landpatch, tempavg_somhr,    defval = 0._r8)
      CALL ncio_read_vector (file_restart, 'ch4_tempavg_finrw',    landpatch, tempavg_finrw,    defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_fsat_bef',         landpatch, fsat_bef,         defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_finundated_lag',   landpatch, finundated_lag,   defval = spval)
	      CALL ncio_read_vector (file_restart, 'ch4_methane_dfsat_tot',landpatch, methane_dfsat_tot, defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_h2osfc',         landpatch, f_h2osfc,         defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_inund_levee_patch',       landpatch, &
	                             f_inund_levee_patch,       defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_inund_flood_patch',       landpatch, &
	                             f_inund_flood_patch,       defval = 0._r8)
	      CALL ncio_read_vector (file_restart, 'ch4_f_inund_flood_depth_patch', landpatch, &
	                             f_inund_flood_depth_patch, defval = 0._r8)
      IF (allocated(restart_ch4_clip_credit_mass)) CALL ncio_read_vector (file_restart, &
         'ch4_restart_ch4_clip_credit_mass', landpatch, restart_ch4_clip_credit_mass, defval = 0._r8)

	      WHERE (invalid_restart_value(conc_o2))           conc_o2           = 1._r8
      WHERE (invalid_restart_value(conc_o2_unsat))     conc_o2_unsat     = 1._r8
      WHERE (invalid_restart_value(conc_o2_sat))       conc_o2_sat       = 1._r8
      WHERE (invalid_restart_value(conc_o2_lake))      conc_o2_lake      = 1._r8
      WHERE (invalid_restart_value(conc_methane))      conc_methane      = 1.e-6_r8
      WHERE (invalid_restart_value(conc_methane_unsat)) conc_methane_unsat = 1.e-6_r8
      WHERE (invalid_restart_value(conc_methane_sat))  conc_methane_sat  = 1.e-6_r8
      WHERE (invalid_restart_value(conc_methane_lake)) conc_methane_lake = 0._r8
      WHERE (invalid_restart_value(grnd_methane_cond) .or. grnd_methane_cond <= 0._r8) &
         grnd_methane_cond = 1.e-6_r8
      WHERE (invalid_restart_value(grnd_methane_cond_unsat) .or. grnd_methane_cond_unsat <= 0._r8) &
         grnd_methane_cond_unsat = 1.e-6_r8
      WHERE (invalid_restart_value(grnd_methane_cond_sat) .or. grnd_methane_cond_sat <= 0._r8) &
         grnd_methane_cond_sat = 1.e-6_r8
	      WHERE (invalid_restart_value(grnd_methane_cond_lake) .or. grnd_methane_cond_lake <= 0._r8) &
	         grnd_methane_cond_lake = 1.e-6_r8
	      WHERE (invalid_restart_value(f_inund_levee_patch) .or. f_inund_levee_patch < 0._r8) &
	         f_inund_levee_patch = 0._r8
	      WHERE (invalid_restart_value(f_inund_flood_patch) .or. f_inund_flood_patch < 0._r8) &
	         f_inund_flood_patch = 0._r8
	      WHERE (invalid_restart_value(f_inund_flood_depth_patch) .or. f_inund_flood_depth_patch < 0._r8) &
	         f_inund_flood_depth_patch = 0._r8
	      WHERE (f_inund_levee_patch > 1._r8) f_inund_levee_patch = 1._r8
	      WHERE (f_inund_flood_patch > 1._r8) f_inund_flood_patch = 1._r8

	      ! invalid_restart_value catches NaN/spval but not negatives.  A stale
      ! restart with sub-zero CH4 or O2 would feed phase-partition and
      ! Michaelis-Menten kinetics, yielding negative oxidation / production.
      ! Clip to zero defensively.  CH4 clips are credited before mutation so
      ! restart sanitation is visible in the same budget family as in-step
      ! negative-concentration clip credits.
      IF (allocated(restart_ch4_clip_credit_mass)) THEN
         CALL credit_methane_restart_clip(conc_methane,       restart_ch4_clip_credit_mass)
         CALL credit_methane_restart_clip(conc_methane_unsat, restart_ch4_clip_credit_mass)
         CALL credit_methane_restart_clip(conc_methane_sat,   restart_ch4_clip_credit_mass)
         CALL credit_methane_restart_clip(conc_methane_lake,  restart_ch4_clip_credit_mass)
      ENDIF
      WHERE (conc_o2           < 0._r8) conc_o2           = 0._r8
      WHERE (conc_o2_unsat     < 0._r8) conc_o2_unsat     = 0._r8
      WHERE (conc_o2_sat       < 0._r8) conc_o2_sat       = 0._r8
      WHERE (conc_o2_lake      < 0._r8) conc_o2_lake      = 0._r8
      WHERE (conc_methane      < 0._r8) conc_methane      = 0._r8
      WHERE (conc_methane_unsat < 0._r8) conc_methane_unsat = 0._r8
      WHERE (conc_methane_sat  < 0._r8) conc_methane_sat  = 0._r8
      WHERE (conc_methane_lake < 0._r8) conc_methane_lake = 0._r8

      ! Only force the cold-start sentinel when the *combined* column
      ! inventory (totcol_methane) is missing — that is the field the
      ! first-step balance check actually references.  A restart that
      ! lacks one of the sat/unsat/lake split inventories (e.g. an older
      ! file written before those were persisted) should NOT discard a
      ! valid fsat_bef; the per-variable defval at lines 1000-1003 will
      ! zero the missing split fields and physics will reconcile them on
      ! the next time step.
      WHERE (invalid_restart_value(totcol_methane))
         fsat_bef = spval
      END WHERE
      WHERE (invalid_restart_value(totcol_methane))       totcol_methane       = 0._r8
      WHERE (invalid_restart_value(totcol_methane_unsat)) totcol_methane_unsat = 0._r8
      WHERE (invalid_restart_value(totcol_methane_sat))   totcol_methane_sat   = 0._r8
      WHERE (invalid_restart_value(totcol_methane_lake))  totcol_methane_lake  = 0._r8
      IF (allocated(restart_ch4_clip_credit_mass)) THEN
         CALL credit_methane_restart_clip(totcol_methane,       restart_ch4_clip_credit_mass)
         CALL credit_methane_restart_clip(totcol_methane_unsat, restart_ch4_clip_credit_mass)
         CALL credit_methane_restart_clip(totcol_methane_sat,   restart_ch4_clip_credit_mass)
         CALL credit_methane_restart_clip(totcol_methane_lake,  restart_ch4_clip_credit_mass)
         IF (allocated(methane_ch4_clip_credit)) THEN
            methane_ch4_clip_credit(:) = methane_ch4_clip_credit(:) - restart_ch4_clip_credit_mass(:)
         ENDIF
      ENDIF
      WHERE (totcol_methane       < 0._r8) totcol_methane       = 0._r8
      WHERE (totcol_methane_unsat < 0._r8) totcol_methane_unsat = 0._r8
      WHERE (totcol_methane_sat   < 0._r8) totcol_methane_sat   = 0._r8
      WHERE (totcol_methane_lake  < 0._r8) totcol_methane_lake  = 0._r8

      ! Clean spval/NaN/negative lake_soilc.  Without this guard a stale
      ! restart sentinel (~1e36) would pass through sum(max(...,0)) as a
      ! huge positive value and the surface-data fallback below would
      ! skip the patch, leaving runaway lake CH4 production.
      WHERE (invalid_restart_value(lake_soilc) .or. lake_soilc < 0._r8) &
         lake_soilc = 0._r8
	   END SUBROUTINE read_methane_restart

   SUBROUTINE save_methane_lulcc_state ()
      IF (.not. allocated(conc_methane)) THEN
         methane_lulcc_snapshot_valid = .false.
         RETURN
      ENDIF

      CALL clear_methane_lulcc_snapshot ()

      allocate(lulcc_conc_o2_old(nl_soil,size(conc_o2,2))); lulcc_conc_o2_old = conc_o2
      allocate(lulcc_conc_methane_old(nl_soil,size(conc_methane,2))); lulcc_conc_methane_old = conc_methane
      allocate(lulcc_co2_decomp_depth_old(nl_soil,size(co2_decomp_depth,2))); lulcc_co2_decomp_depth_old = co2_decomp_depth
      allocate(lulcc_co2_oxid_depth_old(nl_soil,size(co2_oxid_depth,2))); lulcc_co2_oxid_depth_old = co2_oxid_depth
      allocate(lulcc_co2_aere_depth_old(nl_soil,size(co2_aere_depth,2))); lulcc_co2_aere_depth_old = co2_aere_depth
      allocate(lulcc_co2_decomp_depth_unsat_old(nl_soil,size(co2_decomp_depth_unsat,2))); lulcc_co2_decomp_depth_unsat_old = co2_decomp_depth_unsat
      allocate(lulcc_co2_decomp_depth_sat_old(nl_soil,size(co2_decomp_depth_sat,2))); lulcc_co2_decomp_depth_sat_old = co2_decomp_depth_sat
      allocate(lulcc_co2_oxid_depth_unsat_old(nl_soil,size(co2_oxid_depth_unsat,2))); lulcc_co2_oxid_depth_unsat_old = co2_oxid_depth_unsat
      allocate(lulcc_co2_oxid_depth_sat_old(nl_soil,size(co2_oxid_depth_sat,2))); lulcc_co2_oxid_depth_sat_old = co2_oxid_depth_sat
      allocate(lulcc_co2_aere_depth_unsat_old(nl_soil,size(co2_aere_depth_unsat,2))); lulcc_co2_aere_depth_unsat_old = co2_aere_depth_unsat
      allocate(lulcc_co2_aere_depth_sat_old(nl_soil,size(co2_aere_depth_sat,2))); lulcc_co2_aere_depth_sat_old = co2_aere_depth_sat
      allocate(lulcc_co2_decomp_depth_lake_old(nl_soil,size(co2_decomp_depth_lake,2))); lulcc_co2_decomp_depth_lake_old = co2_decomp_depth_lake
      allocate(lulcc_co2_oxid_depth_lake_old(nl_soil,size(co2_oxid_depth_lake,2))); lulcc_co2_oxid_depth_lake_old = co2_oxid_depth_lake
      allocate(lulcc_conc_o2_unsat_old(nl_soil,size(conc_o2_unsat,2))); lulcc_conc_o2_unsat_old = conc_o2_unsat
      allocate(lulcc_conc_o2_sat_old(nl_soil,size(conc_o2_sat,2))); lulcc_conc_o2_sat_old = conc_o2_sat
      allocate(lulcc_conc_methane_unsat_old(nl_soil,size(conc_methane_unsat,2))); lulcc_conc_methane_unsat_old = conc_methane_unsat
      allocate(lulcc_conc_methane_sat_old(nl_soil,size(conc_methane_sat,2))); lulcc_conc_methane_sat_old = conc_methane_sat
      allocate(lulcc_conc_o2_lake_old(nl_soil,size(conc_o2_lake,2))); lulcc_conc_o2_lake_old = conc_o2_lake
      allocate(lulcc_conc_methane_lake_old(nl_soil,size(conc_methane_lake,2))); lulcc_conc_methane_lake_old = conc_methane_lake
      allocate(lulcc_totcol_methane_old(size(totcol_methane))); lulcc_totcol_methane_old = totcol_methane
      allocate(lulcc_co2_decomp_tot_old(size(co2_decomp_tot))); lulcc_co2_decomp_tot_old = co2_decomp_tot
      allocate(lulcc_co2_oxid_tot_old(size(co2_oxid_tot))); lulcc_co2_oxid_tot_old = co2_oxid_tot
      allocate(lulcc_co2_aere_tot_old(size(co2_aere_tot))); lulcc_co2_aere_tot_old = co2_aere_tot
      allocate(lulcc_co2_net_tot_old(size(co2_net_tot))); lulcc_co2_net_tot_old = co2_net_tot
      allocate(lulcc_co2_decomp_tot_unsat_old(size(co2_decomp_tot_unsat))); lulcc_co2_decomp_tot_unsat_old = co2_decomp_tot_unsat
      allocate(lulcc_co2_decomp_tot_sat_old(size(co2_decomp_tot_sat))); lulcc_co2_decomp_tot_sat_old = co2_decomp_tot_sat
      allocate(lulcc_co2_oxid_tot_unsat_old(size(co2_oxid_tot_unsat))); lulcc_co2_oxid_tot_unsat_old = co2_oxid_tot_unsat
      allocate(lulcc_co2_oxid_tot_sat_old(size(co2_oxid_tot_sat))); lulcc_co2_oxid_tot_sat_old = co2_oxid_tot_sat
      allocate(lulcc_co2_net_tot_unsat_old(size(co2_net_tot_unsat))); lulcc_co2_net_tot_unsat_old = co2_net_tot_unsat
      allocate(lulcc_co2_net_tot_sat_old(size(co2_net_tot_sat))); lulcc_co2_net_tot_sat_old = co2_net_tot_sat
      allocate(lulcc_co2_decomp_tot_lake_old(size(co2_decomp_tot_lake))); lulcc_co2_decomp_tot_lake_old = co2_decomp_tot_lake
      allocate(lulcc_co2_oxid_tot_lake_old(size(co2_oxid_tot_lake))); lulcc_co2_oxid_tot_lake_old = co2_oxid_tot_lake
      allocate(lulcc_co2_net_tot_lake_old(size(co2_net_tot_lake))); lulcc_co2_net_tot_lake_old = co2_net_tot_lake
      allocate(lulcc_totcol_methane_unsat_old(size(totcol_methane_unsat))); lulcc_totcol_methane_unsat_old = totcol_methane_unsat
      allocate(lulcc_totcol_methane_sat_old(size(totcol_methane_sat))); lulcc_totcol_methane_sat_old = totcol_methane_sat
      allocate(lulcc_totcol_methane_lake_old(size(totcol_methane_lake))); lulcc_totcol_methane_lake_old = totcol_methane_lake
      allocate(lulcc_grnd_methane_cond_old(size(grnd_methane_cond))); lulcc_grnd_methane_cond_old = grnd_methane_cond
      allocate(lulcc_grnd_methane_cond_unsat_old(size(grnd_methane_cond_unsat))); lulcc_grnd_methane_cond_unsat_old = grnd_methane_cond_unsat
      allocate(lulcc_grnd_methane_cond_sat_old(size(grnd_methane_cond_sat))); lulcc_grnd_methane_cond_sat_old = grnd_methane_cond_sat
      allocate(lulcc_grnd_methane_cond_lake_old(size(grnd_methane_cond_lake))); lulcc_grnd_methane_cond_lake_old = grnd_methane_cond_lake
      allocate(lulcc_layer_sat_lag_old(nl_soil,size(layer_sat_lag,2))); lulcc_layer_sat_lag_old = layer_sat_lag
      allocate(lulcc_lake_soilc_old(nl_soil,size(lake_soilc,2))); lulcc_lake_soilc_old = lake_soilc
      allocate(lulcc_annavg_agnpp_old(size(annavg_agnpp))); lulcc_annavg_agnpp_old = annavg_agnpp
      allocate(lulcc_annavg_bgnpp_old(size(annavg_bgnpp))); lulcc_annavg_bgnpp_old = annavg_bgnpp
      allocate(lulcc_annavg_somhr_old(size(annavg_somhr))); lulcc_annavg_somhr_old = annavg_somhr
      allocate(lulcc_annavg_finrw_old(size(annavg_finrw))); lulcc_annavg_finrw_old = annavg_finrw
      allocate(lulcc_tempavg_agnpp_old(size(tempavg_agnpp))); lulcc_tempavg_agnpp_old = tempavg_agnpp
      allocate(lulcc_tempavg_bgnpp_old(size(tempavg_bgnpp))); lulcc_tempavg_bgnpp_old = tempavg_bgnpp
      allocate(lulcc_annsum_counter_old(size(annsum_counter))); lulcc_annsum_counter_old = annsum_counter
      allocate(lulcc_tempavg_somhr_old(size(tempavg_somhr))); lulcc_tempavg_somhr_old = tempavg_somhr
      allocate(lulcc_tempavg_finrw_old(size(tempavg_finrw))); lulcc_tempavg_finrw_old = tempavg_finrw
      allocate(lulcc_fsat_bef_old(size(fsat_bef))); lulcc_fsat_bef_old = fsat_bef
      allocate(lulcc_finundated_lag_old(size(finundated_lag))); lulcc_finundated_lag_old = finundated_lag
      allocate(lulcc_methane_dfsat_tot_old(size(methane_dfsat_tot))); lulcc_methane_dfsat_tot_old = methane_dfsat_tot
      allocate(lulcc_f_h2osfc_old(size(f_h2osfc))); lulcc_f_h2osfc_old = f_h2osfc
      allocate(lulcc_forc_pmethanem_old(size(forc_pmethanem))); lulcc_forc_pmethanem_old = forc_pmethanem
      allocate(lulcc_c_atm_old(3,size(c_atm,2))); lulcc_c_atm_old = c_atm

      methane_lulcc_snapshot_valid = .true.
   END SUBROUTINE save_methane_lulcc_state

   SUBROUTINE publish_methane_levee_flood_patch (fldfrc_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)

      integer :: ncopy

      IF (.not. allocated(f_inund_levee_patch)) RETURN

      ncopy = min(size(f_inund_levee_patch), size(fldfrc_patch))
      IF (ncopy > 0) f_inund_levee_patch(1:ncopy) = fldfrc_patch(1:ncopy)
      IF (size(f_inund_levee_patch) > ncopy) f_inund_levee_patch(ncopy+1:) = 0._r8

   END SUBROUTINE publish_methane_levee_flood_patch

   SUBROUTINE publish_methane_flood_patch (fldfrc_patch, flddph_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)
      real(r8), intent(in) :: flddph_patch(:)

      integer :: ncopy

      IF (.not. allocated(f_inund_flood_patch)) RETURN

      ncopy = min(size(f_inund_flood_patch), size(fldfrc_patch))
      IF (ncopy > 0) f_inund_flood_patch(1:ncopy) = fldfrc_patch(1:ncopy)
      IF (size(f_inund_flood_patch) > ncopy) f_inund_flood_patch(ncopy+1:) = 0._r8

      IF (allocated(f_inund_flood_depth_patch)) THEN
         ncopy = min(size(f_inund_flood_depth_patch), size(flddph_patch))
         IF (ncopy > 0) f_inund_flood_depth_patch(1:ncopy) = flddph_patch(1:ncopy)
         IF (size(f_inund_flood_depth_patch) > ncopy) f_inund_flood_depth_patch(ncopy+1:) = 0._r8
      ENDIF

   END SUBROUTINE publish_methane_flood_patch


	   SUBROUTINE remap_methane_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
	      lccpct_patches, old_patch_area, new_patch_area)
	      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
	      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
	      real(r8), intent(in), optional :: lccpct_patches(:,:)
	      real(r8), intent(in), optional :: old_patch_area(:)
	      real(r8), intent(in), optional :: new_patch_area(:)
      integer :: nnew

	      nnew = size(patchclass_new)
	      IF (allocated(conc_methane)) CALL deallocate_methane_state ()
	      CALL allocate_methane_state (nnew)
	      CALL init_methane_wetland_fraction_cache (nnew)

	      IF (.not. methane_lulcc_snapshot_valid) RETURN

      CALL remap2d(lulcc_conc_o2_old,              conc_o2)
      CALL remap2d(lulcc_conc_methane_old,         conc_methane)
      CALL remap2d(lulcc_co2_decomp_depth_old,     co2_decomp_depth)
      CALL remap2d(lulcc_co2_oxid_depth_old,       co2_oxid_depth)
      CALL remap2d(lulcc_co2_aere_depth_old,       co2_aere_depth)
      CALL remap2d(lulcc_co2_decomp_depth_unsat_old, co2_decomp_depth_unsat)
      CALL remap2d(lulcc_co2_decomp_depth_sat_old,   co2_decomp_depth_sat)
      CALL remap2d(lulcc_co2_oxid_depth_unsat_old, co2_oxid_depth_unsat)
      CALL remap2d(lulcc_co2_oxid_depth_sat_old,   co2_oxid_depth_sat)
      CALL remap2d(lulcc_co2_aere_depth_unsat_old, co2_aere_depth_unsat)
      CALL remap2d(lulcc_co2_aere_depth_sat_old,   co2_aere_depth_sat)
      CALL remap2d(lulcc_co2_decomp_depth_lake_old, co2_decomp_depth_lake)
      CALL remap2d(lulcc_co2_oxid_depth_lake_old,  co2_oxid_depth_lake)
      CALL remap2d(lulcc_conc_o2_unsat_old,        conc_o2_unsat)
      CALL remap2d(lulcc_conc_o2_sat_old,          conc_o2_sat)
      CALL remap2d(lulcc_conc_methane_unsat_old,   conc_methane_unsat)
      CALL remap2d(lulcc_conc_methane_sat_old,     conc_methane_sat)
      CALL remap2d(lulcc_conc_o2_lake_old,         conc_o2_lake)
      CALL remap2d(lulcc_conc_methane_lake_old,    conc_methane_lake)
      CALL remap2d(lulcc_layer_sat_lag_old,        layer_sat_lag)
      CALL remap2d(lulcc_lake_soilc_old,           lake_soilc)
      CALL remap2d(lulcc_c_atm_old,                c_atm)

      ! Concentrations, conductances, fractions, rates, and rolling diagnostics
      ! are remapped as intensive fields.  Column stocks and total per-area
      ! source/flux diagnostics are remapped as area-mass fields when both old
      ! and new patch areas are available, matching land-water tracer LULCC
      ! semantics for extensive per-area pools.
      CALL remap1d_mass(lulcc_totcol_methane_old,       totcol_methane)
      CALL remap1d_mass(lulcc_co2_decomp_tot_old,       co2_decomp_tot)
      CALL remap1d_mass(lulcc_co2_oxid_tot_old,         co2_oxid_tot)
      CALL remap1d_mass(lulcc_co2_aere_tot_old,         co2_aere_tot)
      CALL remap1d_mass(lulcc_co2_net_tot_old,          co2_net_tot)
      CALL remap1d_mass(lulcc_co2_decomp_tot_unsat_old, co2_decomp_tot_unsat)
      CALL remap1d_mass(lulcc_co2_decomp_tot_sat_old,   co2_decomp_tot_sat)
      CALL remap1d_mass(lulcc_co2_oxid_tot_unsat_old,   co2_oxid_tot_unsat)
      CALL remap1d_mass(lulcc_co2_oxid_tot_sat_old,     co2_oxid_tot_sat)
      CALL remap1d_mass(lulcc_co2_net_tot_unsat_old,    co2_net_tot_unsat)
      CALL remap1d_mass(lulcc_co2_net_tot_sat_old,      co2_net_tot_sat)
      CALL remap1d_mass(lulcc_co2_decomp_tot_lake_old,  co2_decomp_tot_lake)
      CALL remap1d_mass(lulcc_co2_oxid_tot_lake_old,    co2_oxid_tot_lake)
      CALL remap1d_mass(lulcc_co2_net_tot_lake_old,     co2_net_tot_lake)
      CALL remap1d_mass(lulcc_totcol_methane_unsat_old, totcol_methane_unsat)
      CALL remap1d_mass(lulcc_totcol_methane_sat_old,   totcol_methane_sat)
      CALL remap1d_mass(lulcc_totcol_methane_lake_old,  totcol_methane_lake)
      CALL remap1d(lulcc_grnd_methane_cond_old,    grnd_methane_cond)
      CALL remap1d(lulcc_grnd_methane_cond_unsat_old, grnd_methane_cond_unsat)
      CALL remap1d(lulcc_grnd_methane_cond_sat_old,   grnd_methane_cond_sat)
      CALL remap1d(lulcc_grnd_methane_cond_lake_old,  grnd_methane_cond_lake)
      CALL remap1d(lulcc_annavg_agnpp_old,         annavg_agnpp)
      CALL remap1d(lulcc_annavg_bgnpp_old,         annavg_bgnpp)
      CALL remap1d(lulcc_annavg_somhr_old,         annavg_somhr)
      CALL remap1d(lulcc_annavg_finrw_old,         annavg_finrw)
      CALL remap1d(lulcc_tempavg_agnpp_old,        tempavg_agnpp)
      CALL remap1d(lulcc_tempavg_bgnpp_old,        tempavg_bgnpp)
      CALL remap1d(lulcc_annsum_counter_old,       annsum_counter)
      CALL remap1d(lulcc_tempavg_somhr_old,        tempavg_somhr)
      CALL remap1d(lulcc_tempavg_finrw_old,        tempavg_finrw)
      CALL remap1d(lulcc_fsat_bef_old,             fsat_bef)
      CALL remap1d(lulcc_finundated_lag_old,       finundated_lag)
      CALL remap1d_mass(lulcc_methane_dfsat_tot_old,    methane_dfsat_tot)
      CALL remap1d(lulcc_f_h2osfc_old,             f_h2osfc)
      CALL remap1d(lulcc_forc_pmethanem_old,       forc_pmethanem)

      CALL clear_methane_lulcc_snapshot ()

   CONTAINS
      SUBROUTINE remap1d(old, new)
         real(r8), intent(in) :: old(:)
         real(r8), intent(inout) :: new(:)
         integer :: np, op, src
         real(r8) :: w, wsum, val
         DO np = 1, min(size(new), nnew)
            wsum = 0._r8
            val = 0._r8
            IF (present(lccpct_patches)) THEN
               DO op = 1, min(size(old), size(patchclass_old), size(eindex_old))
                  IF (eindex_old(op) /= eindex_new(np)) CYCLE
                  IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                      patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
	                  w = lulcc_source_weight(np, op)
	                  IF (w <= 0._r8) CYCLE
                  val = val + w * old(op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               new(np) = val / wsum
            ELSE
               src = fallback_source(np, size(old))
               IF (src > 0) new(np) = old(src)
            ENDIF
         ENDDO
      END SUBROUTINE remap1d

      SUBROUTINE remap1d_mass(old, new)
         real(r8), intent(in) :: old(:)
         real(r8), intent(inout) :: new(:)
         integer :: np, op, src
         real(r8) :: w, wsum, val, denom
         logical :: conserve_area_mass

         DO np = 1, min(size(new), nnew)
            wsum = 0._r8
            val = 0._r8
            IF (present(lccpct_patches)) THEN
               conserve_area_mass = area_mass_remap_available(np)
               DO op = 1, min(size(old), size(patchclass_old), size(eindex_old))
                  IF (eindex_old(op) /= eindex_new(np)) CYCLE
                  IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                      patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
                  IF (conserve_area_mass) THEN
                     w = lulcc_mass_transfer_area(np, op)
                  ELSE
                     w = lulcc_source_weight(np, op)
                  ENDIF
                  IF (w <= 0._r8) CYCLE
                  val = val + w * old(op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               denom = remap_denominator(np, wsum, conserve_area_mass)
               new(np) = val / denom
            ELSE
               src = fallback_source(np, size(old))
               IF (src > 0) new(np) = old(src)
            ENDIF
         ENDDO
      END SUBROUTINE remap1d_mass

      SUBROUTINE remap2d(old, new)
         real(r8), intent(in) :: old(:,:)
         real(r8), intent(inout) :: new(:,:)
         integer :: np, op, src
         real(r8) :: w, wsum
         real(r8) :: default_vals(size(new,1))
         DO np = 1, min(size(new,2), nnew)
            wsum = 0._r8
            default_vals(:) = new(:,np)
            IF (present(lccpct_patches)) THEN
               new(:,np) = 0._r8
               DO op = 1, min(size(old,2), size(patchclass_old), size(eindex_old))
                  IF (eindex_old(op) /= eindex_new(np)) CYCLE
                  IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
                      patchclass_old(op) > ubound(lccpct_patches,2)) CYCLE
	                  w = lulcc_source_weight(np, op)
	                  IF (w <= 0._r8) CYCLE
                  new(1:min(size(new,1),size(old,1)),np) = &
                     new(1:min(size(new,1),size(old,1)),np) + &
                     w * old(1:min(size(new,1),size(old,1)),op)
                  wsum = wsum + w
               ENDDO
            ENDIF
            IF (wsum > 0._r8) THEN
               new(:,np) = new(:,np) / wsum
            ELSE
               src = fallback_source(np, size(old,2))
               IF (src > 0) THEN
                  new(:,np) = 0._r8
                  new(1:min(size(new,1),size(old,1)),np) = old(1:min(size(new,1),size(old,1)),src)
               ELSE
                  new(:,np) = default_vals(:)
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE remap2d

	      INTEGER FUNCTION fallback_source(np, old_n) RESULT(src)
         integer, intent(in) :: np, old_n
         integer :: op
         src = 0
         DO op = 1, min(old_n, size(patchclass_old), size(eindex_old))
            IF (eindex_old(op) == eindex_new(np) .and. &
                patchclass_old(op) == patchclass_new(np)) THEN
               src = op
               RETURN
            ENDIF
         ENDDO
         ! Do not fall back across patch classes.  Methane lake/sediment
         ! inventories are class-specific and can be corrupted by copying from
         ! a same-eindex non-lake patch.
	      END FUNCTION fallback_source

	      REAL(r8) FUNCTION lulcc_source_weight(np, op) RESULT(w)
	         integer, intent(in) :: np, op
	         integer :: oq
	         real(r8) :: class_area

	         w = 0._r8
	         IF (.not. present(lccpct_patches)) RETURN
	         IF (np > size(lccpct_patches,1)) RETURN
	         IF (patchclass_old(op) < lbound(lccpct_patches,2) .or. &
	             patchclass_old(op) > ubound(lccpct_patches,2)) RETURN
	         w = max(0._r8, lccpct_patches(np, patchclass_old(op)))
	         IF (w <= 0._r8 .or. .not. present(old_patch_area)) RETURN
	         IF (op > size(old_patch_area)) RETURN

	         class_area = 0._r8
	         DO oq = 1, min(size(patchclass_old), size(eindex_old), size(old_patch_area))
	            IF (eindex_old(oq) == eindex_new(np) .and. &
	                patchclass_old(oq) == patchclass_old(op)) THEN
	               class_area = class_area + max(0._r8, old_patch_area(oq))
	            ENDIF
	         ENDDO
	         IF (class_area > 0._r8) THEN
	            w = w * max(0._r8, old_patch_area(op)) / class_area
	         ENDIF
	      END FUNCTION lulcc_source_weight

	      LOGICAL FUNCTION area_mass_remap_available(np) RESULT(ok)
	         integer, intent(in) :: np

	         ok = present(lccpct_patches) .and. present(old_patch_area) .and. present(new_patch_area)
	         IF (.not. ok) RETURN
	         ok = np <= size(new_patch_area)
	         IF (.not. ok) RETURN
	         ok = new_patch_area(np) > tiny(1._r8)
	      END FUNCTION area_mass_remap_available

	      REAL(r8) FUNCTION lulcc_mass_transfer_area(np, op) RESULT(w)
	         integer, intent(in) :: np, op
	         integer :: c, nq
	         real(r8) :: target_area, class_target_area

	         w = 0._r8
	         IF (.not. area_mass_remap_available(np)) RETURN
	         IF (op > size(old_patch_area)) RETURN
	         IF (op > size(patchclass_old) .or. op > size(eindex_old)) RETURN
	         c = patchclass_old(op)
	         IF (c < lbound(lccpct_patches,2) .or. c > ubound(lccpct_patches,2)) RETURN

	         target_area = lulcc_target_class_area(np, c)
	         IF (target_area <= tiny(1._r8)) RETURN

	         class_target_area = 0._r8
	         DO nq = 1, min(nnew, size(eindex_new), size(new_patch_area))
	            IF (eindex_new(nq) == eindex_new(np)) THEN
	               class_target_area = class_target_area + lulcc_target_class_area(nq, c)
	            ENDIF
	         ENDDO
	         IF (class_target_area <= tiny(1._r8)) RETURN

	         w = max(0._r8, old_patch_area(op)) * target_area / class_target_area
	      END FUNCTION lulcc_mass_transfer_area

	      REAL(r8) FUNCTION lulcc_target_class_area(np, c) RESULT(area)
	         integer, intent(in) :: np, c
	         integer :: cc
	         real(r8) :: class_sum

	         area = 0._r8
	         IF (.not. present(lccpct_patches)) RETURN
	         IF (.not. present(new_patch_area)) RETURN
	         IF (np > size(new_patch_area)) RETURN
	         IF (np > size(lccpct_patches,1)) RETURN
	         IF (c < lbound(lccpct_patches,2) .or. c > ubound(lccpct_patches,2)) RETURN

	         class_sum = 0._r8
	         DO cc = lbound(lccpct_patches,2), ubound(lccpct_patches,2)
	            class_sum = class_sum + max(0._r8, lccpct_patches(np, cc))
	         ENDDO
	         IF (class_sum <= tiny(1._r8)) RETURN

	         area = max(0._r8, new_patch_area(np)) * max(0._r8, lccpct_patches(np, c)) / class_sum
	      END FUNCTION lulcc_target_class_area

	      REAL(r8) FUNCTION remap_denominator(np, wsum, conserve_mass) RESULT(denom)
	         integer, intent(in) :: np
	         real(r8), intent(in) :: wsum
	         logical, intent(in) :: conserve_mass

	         denom = max(wsum, tiny(1._r8))
	         IF (conserve_mass .and. present(new_patch_area)) THEN
	            IF (np <= size(new_patch_area) .and. new_patch_area(np) > tiny(1._r8)) THEN
	               denom = new_patch_area(np)
	            ENDIF
	         ENDIF
	      END FUNCTION remap_denominator
	   END SUBROUTINE remap_methane_lulcc_state


   SUBROUTINE clear_methane_lulcc_snapshot ()
      IF (allocated(lulcc_conc_o2_old)) deallocate(lulcc_conc_o2_old)
      IF (allocated(lulcc_conc_methane_old)) deallocate(lulcc_conc_methane_old)
      IF (allocated(lulcc_co2_decomp_depth_old)) deallocate(lulcc_co2_decomp_depth_old)
      IF (allocated(lulcc_co2_oxid_depth_old)) deallocate(lulcc_co2_oxid_depth_old)
      IF (allocated(lulcc_co2_aere_depth_old)) deallocate(lulcc_co2_aere_depth_old)
      IF (allocated(lulcc_co2_decomp_depth_unsat_old)) deallocate(lulcc_co2_decomp_depth_unsat_old)
      IF (allocated(lulcc_co2_decomp_depth_sat_old)) deallocate(lulcc_co2_decomp_depth_sat_old)
      IF (allocated(lulcc_co2_oxid_depth_unsat_old)) deallocate(lulcc_co2_oxid_depth_unsat_old)
      IF (allocated(lulcc_co2_oxid_depth_sat_old)) deallocate(lulcc_co2_oxid_depth_sat_old)
      IF (allocated(lulcc_co2_aere_depth_unsat_old)) deallocate(lulcc_co2_aere_depth_unsat_old)
      IF (allocated(lulcc_co2_aere_depth_sat_old)) deallocate(lulcc_co2_aere_depth_sat_old)
      IF (allocated(lulcc_co2_decomp_depth_lake_old)) deallocate(lulcc_co2_decomp_depth_lake_old)
      IF (allocated(lulcc_co2_oxid_depth_lake_old)) deallocate(lulcc_co2_oxid_depth_lake_old)
      IF (allocated(lulcc_conc_o2_unsat_old)) deallocate(lulcc_conc_o2_unsat_old)
      IF (allocated(lulcc_conc_o2_sat_old)) deallocate(lulcc_conc_o2_sat_old)
      IF (allocated(lulcc_conc_methane_unsat_old)) deallocate(lulcc_conc_methane_unsat_old)
      IF (allocated(lulcc_conc_methane_sat_old)) deallocate(lulcc_conc_methane_sat_old)
      IF (allocated(lulcc_conc_o2_lake_old)) deallocate(lulcc_conc_o2_lake_old)
      IF (allocated(lulcc_conc_methane_lake_old)) deallocate(lulcc_conc_methane_lake_old)
      IF (allocated(lulcc_totcol_methane_old)) deallocate(lulcc_totcol_methane_old)
      IF (allocated(lulcc_co2_decomp_tot_old)) deallocate(lulcc_co2_decomp_tot_old)
      IF (allocated(lulcc_co2_oxid_tot_old)) deallocate(lulcc_co2_oxid_tot_old)
      IF (allocated(lulcc_co2_aere_tot_old)) deallocate(lulcc_co2_aere_tot_old)
      IF (allocated(lulcc_co2_net_tot_old)) deallocate(lulcc_co2_net_tot_old)
      IF (allocated(lulcc_co2_decomp_tot_unsat_old)) deallocate(lulcc_co2_decomp_tot_unsat_old)
      IF (allocated(lulcc_co2_decomp_tot_sat_old)) deallocate(lulcc_co2_decomp_tot_sat_old)
      IF (allocated(lulcc_co2_oxid_tot_unsat_old)) deallocate(lulcc_co2_oxid_tot_unsat_old)
      IF (allocated(lulcc_co2_oxid_tot_sat_old)) deallocate(lulcc_co2_oxid_tot_sat_old)
      IF (allocated(lulcc_co2_net_tot_unsat_old)) deallocate(lulcc_co2_net_tot_unsat_old)
      IF (allocated(lulcc_co2_net_tot_sat_old)) deallocate(lulcc_co2_net_tot_sat_old)
      IF (allocated(lulcc_co2_decomp_tot_lake_old)) deallocate(lulcc_co2_decomp_tot_lake_old)
      IF (allocated(lulcc_co2_oxid_tot_lake_old)) deallocate(lulcc_co2_oxid_tot_lake_old)
      IF (allocated(lulcc_co2_net_tot_lake_old)) deallocate(lulcc_co2_net_tot_lake_old)
      IF (allocated(lulcc_totcol_methane_unsat_old)) deallocate(lulcc_totcol_methane_unsat_old)
      IF (allocated(lulcc_totcol_methane_sat_old)) deallocate(lulcc_totcol_methane_sat_old)
      IF (allocated(lulcc_totcol_methane_lake_old)) deallocate(lulcc_totcol_methane_lake_old)
      IF (allocated(lulcc_grnd_methane_cond_old)) deallocate(lulcc_grnd_methane_cond_old)
      IF (allocated(lulcc_grnd_methane_cond_unsat_old)) deallocate(lulcc_grnd_methane_cond_unsat_old)
      IF (allocated(lulcc_grnd_methane_cond_sat_old)) deallocate(lulcc_grnd_methane_cond_sat_old)
      IF (allocated(lulcc_grnd_methane_cond_lake_old)) deallocate(lulcc_grnd_methane_cond_lake_old)
      IF (allocated(lulcc_layer_sat_lag_old)) deallocate(lulcc_layer_sat_lag_old)
      IF (allocated(lulcc_lake_soilc_old)) deallocate(lulcc_lake_soilc_old)
      IF (allocated(lulcc_annavg_agnpp_old)) deallocate(lulcc_annavg_agnpp_old)
      IF (allocated(lulcc_annavg_bgnpp_old)) deallocate(lulcc_annavg_bgnpp_old)
      IF (allocated(lulcc_annavg_somhr_old)) deallocate(lulcc_annavg_somhr_old)
      IF (allocated(lulcc_annavg_finrw_old)) deallocate(lulcc_annavg_finrw_old)
      IF (allocated(lulcc_tempavg_agnpp_old)) deallocate(lulcc_tempavg_agnpp_old)
      IF (allocated(lulcc_tempavg_bgnpp_old)) deallocate(lulcc_tempavg_bgnpp_old)
      IF (allocated(lulcc_annsum_counter_old)) deallocate(lulcc_annsum_counter_old)
      IF (allocated(lulcc_tempavg_somhr_old)) deallocate(lulcc_tempavg_somhr_old)
      IF (allocated(lulcc_tempavg_finrw_old)) deallocate(lulcc_tempavg_finrw_old)
      IF (allocated(lulcc_fsat_bef_old)) deallocate(lulcc_fsat_bef_old)
      IF (allocated(lulcc_finundated_lag_old)) deallocate(lulcc_finundated_lag_old)
      IF (allocated(lulcc_methane_dfsat_tot_old)) deallocate(lulcc_methane_dfsat_tot_old)
      IF (allocated(lulcc_f_h2osfc_old)) deallocate(lulcc_f_h2osfc_old)
      IF (allocated(lulcc_forc_pmethanem_old)) deallocate(lulcc_forc_pmethanem_old)
      IF (allocated(lulcc_c_atm_old)) deallocate(lulcc_c_atm_old)
      methane_lulcc_snapshot_valid = .false.
   END SUBROUTINE clear_methane_lulcc_snapshot


   SUBROUTINE initialize_methane_lake_soilc_from_surface (patchtype_in, lake_soilc_srf_in, allowlakeprod)
      integer,  intent(in) :: patchtype_in(:)
      real(r8), intent(in) :: lake_soilc_srf_in(:,:)
      logical,  intent(in) :: allowlakeprod

      integer :: ipatch, npatch
      integer :: lake_soilc_nlake, lake_soilc_nfilled, lake_soilc_missing
      logical, save :: lake_soilc_missing_warned = .false.
      real(r8), parameter :: smallnumber = 1.e-12_r8

      IF (.not. allocated(lake_soilc)) RETURN
      IF (.not. allowlakeprod) RETURN
      IF (size(lake_soilc_srf_in,1) < nl_soil) RETURN

      npatch = min(size(patchtype_in), size(lake_soilc,2), size(lake_soilc_srf_in,2))
      lake_soilc_nlake = 0
      lake_soilc_nfilled = 0
      lake_soilc_missing = 0
      DO ipatch = 1, npatch
         IF (patchtype_in(ipatch) /= 4) CYCLE
         lake_soilc_nlake = lake_soilc_nlake + 1
         IF (sum(max(lake_soilc(:,ipatch), 0._r8)) > smallnumber) CYCLE
         IF (sum(max(lake_soilc_srf_in(1:nl_soil,ipatch), 0._r8)) <= smallnumber) THEN
            lake_soilc_missing = lake_soilc_missing + 1
            CYCLE
         ENDIF
         lake_soilc(:,ipatch) = max(lake_soilc_srf_in(1:nl_soil,ipatch), 0._r8)
         lake_soilc_nfilled = lake_soilc_nfilled + 1
      END DO
      IF (.not. lake_soilc_missing_warned .and. lake_soilc_nlake > 0 .and. &
          lake_soilc_missing == lake_soilc_nlake .and. lake_soilc_nfilled == 0) THEN
         write(6,*) ' WARNING: lake CH4 production is enabled, but lake_soilc is zero/missing on this rank; ', &
            'lake CH4 production will remain zero for these lake patches.'
         lake_soilc_missing_warned = .true.
      ENDIF
   END SUBROUTINE initialize_methane_lake_soilc_from_surface


   ELEMENTAL LOGICAL FUNCTION invalid_restart_value (x)
      real(r8), intent(in) :: x

      invalid_restart_value = ieee_is_nan(x) .or. (abs(x) >= 0.5_r8 * abs(spval))
   END FUNCTION invalid_restart_value

END MODULE MOD_Tracer_Reactive_Methane_State
#endif
