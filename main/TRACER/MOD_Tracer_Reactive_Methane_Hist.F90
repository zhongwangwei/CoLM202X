#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_Hist

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: methane_reactive_history

CONTAINS

   SUBROUTINE methane_reactive_history (file_hist, itime_in_file, sumarea, filter, &
      nl_soil, forcing_has_missing_value, forcmask_pch)
   USE MOD_Tracer_Hist, only: HistForm, write_history_variable_2d, write_history_variable_3d
   USE MOD_HistGridded
   USE MOD_DataType
   USE MOD_LandPatch, only: numpatch
   USE MOD_SPMD_Task, only: p_is_worker
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Tracer_Methane_BgcLink, only: methane_patch_active_mask
            USE MOD_Tracer_Methane_Registry, only: igas_ch4
	            USE MOD_Tracer_Methane_Const,    only: mhist_on => methane_history_enabled, DEF_METHANE
            USE MOD_Tracer_Methane_AccFlux,  only: &
                 a_net_methane, a_methane_prod_depth, a_o2_decomp_depth, a_co2_decomp_depth, &
                 a_methane_oxid_depth, a_o2_oxid_depth, a_co2_oxid_depth, &
                 a_methane_aere_depth, a_methane_tran_depth, a_o2_aere_depth, a_co2_aere_depth, &
                 a_methane_ebul_depth, a_o2stress, a_methane_stress, &
                 a_methane_surf_flux_tot, a_methane_surf_flux_tot_phys, a_methane_surf_aere, a_methane_surf_ebul, &
                    a_methane_surf_diff, a_methane_surf_diff_phys, a_methane_balance_residual, a_methane_ch4_clip_credit, &
                    a_o2_cap_loss, a_o2_cap_gain, &
                    a_methane_ebul_tot, a_methane_prod_tot, a_methane_oxid_tot, &
                 a_co2_decomp_tot, a_co2_oxid_tot, a_co2_aere_tot, a_co2_net_tot, &
                 a_totcol_methane, a_grnd_methane_cond, a_conc_o2, a_conc_methane, &
                    a_net_methane_unsat, a_net_methane_sat, &
                    a_methane_prod_depth_unsat, a_methane_prod_depth_sat, &
                       a_o2_decomp_depth_unsat, a_o2_decomp_depth_sat, &
                       a_co2_decomp_depth_unsat, a_co2_decomp_depth_sat, &
                       a_methane_oxid_depth_unsat, a_methane_oxid_depth_sat, &
                       a_o2_oxid_depth_unsat, a_o2_oxid_depth_sat, &
                       a_co2_oxid_depth_unsat, a_co2_oxid_depth_sat, &
                       a_methane_aere_depth_unsat, a_methane_aere_depth_sat, &
                       a_methane_tran_depth_unsat, a_methane_tran_depth_sat, &
                       a_o2_aere_depth_unsat, a_o2_aere_depth_sat, &
                       a_co2_aere_depth_unsat, a_co2_aere_depth_sat, &
                       a_methane_ebul_depth_unsat, a_methane_ebul_depth_sat, &
                       a_o2stress_unsat, a_o2stress_sat, &
                       a_methane_stress_unsat, a_methane_stress_sat, &
                       a_methane_surf_flux_tot_unsat, a_methane_surf_flux_tot_sat, &
                       a_methane_surf_aere_unsat, a_methane_surf_aere_sat, &
                       a_methane_surf_ebul_unsat, a_methane_surf_ebul_sat, &
                       a_methane_surf_diff_unsat, a_methane_surf_diff_sat, &
                       a_methane_ebul_tot_unsat, a_methane_ebul_tot_sat, &
                       a_methane_prod_tot_unsat, a_methane_prod_tot_sat, &
                       a_methane_oxid_tot_unsat, a_methane_oxid_tot_sat, &
                       a_co2_decomp_tot_unsat, a_co2_decomp_tot_sat, &
                       a_co2_oxid_tot_unsat, a_co2_oxid_tot_sat, &
                       a_co2_net_tot_unsat, a_co2_net_tot_sat, &
                       a_totcol_methane_unsat, a_totcol_methane_sat, &
                       a_grnd_methane_cond_unsat, a_grnd_methane_cond_sat, &
                    a_conc_o2_unsat, a_conc_o2_sat, &
                    a_conc_methane_unsat, a_conc_methane_sat, &
                 a_methane_prod_depth_lake, a_methane_oxid_depth_lake, a_methane_ebul_depth_lake, &
                 a_co2_decomp_depth_lake, a_co2_oxid_depth_lake, &
                 a_methane_surf_ebul_lake, a_methane_surf_diff_lake, a_methane_surf_flux_tot_lake, &
                 a_methane_prod_tot_lake, a_methane_oxid_tot_lake, a_methane_ebul_tot_lake, &
                 a_co2_decomp_tot_lake, a_co2_oxid_tot_lake, a_co2_net_tot_lake, &
                 a_totcol_methane_lake, a_grnd_methane_cond_lake, &
                 a_conc_o2_lake, a_conc_methane_lake, &
                    a_forc_pmethanem, a_layer_sat_lag, &
                    a_annavg_agnpp, a_annavg_bgnpp, a_annavg_somhr, a_annavg_finrw, &
                    a_methane_dfsat_tot, a_f_h2osfc, &
                    a_methane_finundated, a_methane_soil_finundated, a_methane_soil_zwt, &
                    a_f_inund_flood_patch, a_f_inund_flood_depth_patch, a_wetland_frac_per_patch, &
                    a_methane_surf_flux_wetland, a_methane_surf_flux_soil, &
                    a_methane_surf_flux_lake, a_methane_surf_flux_rice, &
                    a_B_methanogen, a_B_methanotroph, &
                 a_B_methanogen_dormant, a_B_methanotroph_dormant, &
                 a_f_T_methanogen, a_f_S_methanogen, a_f_O2_methanogen, a_f_T_methanotroph, &
                 a_methanogen_growth_rate, a_methanotroph_growth_rate, &
                 a_microbial_prod_potential, a_microbial_oxid_potential, &
                 a_methane_acc_num, a_methane_acc_num_unsat, a_methane_acc_num_sat, &
                 a_methane_acc_num_lake, a_methane_acc_num_extra, a_methane_acc_num_microbe
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: itime_in_file
   type(block_data_real8_2d), intent(inout) :: sumarea
	   logical, intent(inout) :: filter(:)
	   integer, intent(in) :: nl_soil
	   logical, intent(in) :: forcing_has_missing_value
	   logical, intent(in) :: forcmask_pch(:)
	   real(r8), allocatable :: hist_ch4_active_without_lake(:)
	   real(r8), allocatable :: hist_ch4_global_with_lake(:)
	   logical, allocatable :: filter_active_without_lake(:)
	   logical, allocatable :: filter_all_land(:)

	            IF (igas_ch4 > 0) THEN
               ! Outer filter (patchtype==0) excludes wetland (patchtype==2) and
               ! lake (patchtype==4) — the only patches where the CH4 module
               ! actually runs.  Without overriding, every CH4 grid cell aggregates
               ! to zero because only soil patches are summed.  Re-derive a
               ! soil+wetland filter for the main/sat/unsat writes here, then
               ! swap to a lake-only filter before the lake section below.
                  ! Active CH4 mask mirrors the driver gate at
                  ! CoLMDRIVER.F90:381-384 — patchtype==2 and patchtype==0 are
                  ! active by default.  Urban (1) and ice (3) are NEVER active
                  ! even though "<=2" would include urban.
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  ! Use BgcLink helper so the history filter matches the
                  ! driver gate + rice paddy override exactly.  Same call
                  ! site pattern as MOD_Tracer_Methane_AccFlux.F90.
                  CALL methane_patch_active_mask (numpatch, DEF_METHANE%only_wetland, &
                                                  DEF_METHANE%enable_rice_paddy, &
                                                  patchtype, filter)
                  filter = filter .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
	               IF (HistForm == 'Gridded') THEN
	                  CALL mp2g_hist%get_sumarea (sumarea, filter)
	               ENDIF

	               ! Preserve un-normalized accumulators before write_history_variable_2d
	               ! divides them in place.  These two derived history fields use an
	               ! all-land denominator so they can be integrated directly with
	               ! landarea.  The legacy f_methane_surf_flux_tot below remains the
	               ! active CH4 patch mean and excludes lake patches.
	               allocate (hist_ch4_active_without_lake(numpatch))
	               allocate (hist_ch4_global_with_lake(numpatch))
	               allocate (filter_active_without_lake(numpatch))
	               allocate (filter_all_land(numpatch))
	               hist_ch4_active_without_lake(:) = 0._r8
	               hist_ch4_global_with_lake(:) = 0._r8
	               filter_active_without_lake(:) = .false.
	               filter_all_land(:) = .false.
	               IF ((p_is_worker) .and. (numpatch > 0)) THEN
	                  CALL methane_patch_active_mask (numpatch, DEF_METHANE%only_wetland, &
	                                                  DEF_METHANE%enable_rice_paddy, &
	                                                  patchtype, filter_active_without_lake)
	                  filter_active_without_lake = filter_active_without_lake .and. patchmask
	                  filter_all_land = (patchtype < 99) .and. patchmask
	                  IF (forcing_has_missing_value) THEN
	                     filter_active_without_lake = filter_active_without_lake .and. forcmask_pch
	                     filter_all_land = filter_all_land .and. forcmask_pch
	                  ENDIF
	                  hist_ch4_active_without_lake = a_methane_surf_flux_tot
	                  WHERE (filter_active_without_lake .or. &
	                         ((patchtype == 4) .and. DEF_METHANE%allowlakeprod .and. filter_all_land))
	                     hist_ch4_global_with_lake = a_methane_surf_flux_tot
	                  ELSEWHERE
	                     hist_ch4_global_with_lake = 0._r8
	                  END WHERE
	               ENDIF

	               CALL write_history_variable_2d (mhist_on('f_net_methane'), a_net_methane, file_hist, &
                  'f_net_methane', itime_in_file, sumarea, filter, &
                  'average net methane correction to CO2 flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth'), a_methane_prod_depth, file_hist, &
                  'f_methane_prod_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth'), a_o2_decomp_depth, file_hist, &
                  'f_o2_decomp_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 consumption during decomposition per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth'), a_co2_decomp_depth, file_hist, &
                  'f_co2_decomp_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth'), a_methane_oxid_depth, file_hist, &
                  'f_methane_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth'), a_o2_oxid_depth, file_hist, &
                  'f_o2_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 consumption via oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth'), a_co2_oxid_depth, file_hist, &
                  'f_co2_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CO2 production from CH4 oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_aere_depth'), a_methane_aere_depth, file_hist, &
                  'f_methane_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via aerenchyma per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_tran_depth'), a_methane_tran_depth, file_hist, &
                  'f_methane_tran_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via transpiration per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_aere_depth'), a_o2_aere_depth, file_hist, &
                  'f_o2_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 gain via aerenchyma per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_aere_depth'), a_co2_aere_depth, file_hist, &
                  'f_co2_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'diagnosed CO2 aerenchyma flux per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth'), a_methane_ebul_depth, file_hist, &
                  'f_methane_ebul_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via ebullition per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2stress'), a_o2stress, file_hist, &
                  'f_o2stress', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 availability/demand ratio', '-', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_stress'), a_methane_stress, file_hist, &
                  'f_methane_stress', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 availability/sinks ratio', '-', &
	                  acc_num=a_methane_acc_num)
	                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot'), a_methane_surf_flux_tot, file_hist, &
	                     'f_methane_surf_flux_tot', itime_in_file, sumarea, filter, &
	                     'legacy active CH4 total surface flux; lake excluded; active-area mean', 'mol/m2/s', &
	                     acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_phys'), a_methane_surf_flux_tot_phys, file_hist, &
                     'f_methane_surf_flux_tot_phys', itime_in_file, sumarea, filter, &
                     'CH4 physical total surface flux (before CH4 clip and closure residual)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere'), a_methane_surf_aere, file_hist, &
                  'f_methane_surf_aere', itime_in_file, sumarea, filter, &
                  'CH4 column aerenchyma flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul'), a_methane_surf_ebul, file_hist, &
                  'f_methane_surf_ebul', itime_in_file, sumarea, filter, &
                  'CH4 surface ebullition flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff'), a_methane_surf_diff, file_hist, &
                        'f_methane_surf_diff', itime_in_file, sumarea, filter, &
                        'CH4 diffusive surface flux plus inundated-fraction and closure adjustments', 'mol/m2/s', &
                        acc_num=a_methane_acc_num)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_phys'), a_methane_surf_diff_phys, file_hist, &
                        'f_methane_surf_diff_phys', itime_in_file, sumarea, filter, &
                        'CH4 pure-physics diffusive surface flux (before negative-conc clip and closure residual)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_balance_residual'), a_methane_balance_residual, file_hist, &
                     'f_methane_balance_residual', itime_in_file, sumarea, filter, &
                     'CH4 numerical closure flux credited to f_methane_surf_diff', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_ch4_clip_credit'), a_methane_ch4_clip_credit, file_hist, &
                     'f_methane_ch4_clip_credit', itime_in_file, sumarea, filter, &
                     'CH4 negative-concentration clip credit included in reported diffusive flux', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_o2_cap_loss'), a_o2_cap_loss, file_hist, &
                     'f_o2_cap_loss', itime_in_file, sumarea, filter, &
                     'O2 column sink from post-solve physical concentration cap', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_o2_cap_gain'), a_o2_cap_gain, file_hist, &
                     'f_o2_cap_gain', itime_in_file, sumarea, filter, &
                     'O2 column source from post-solve nonnegative concentration floor', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot'), a_methane_ebul_tot, file_hist, &
                  'f_methane_ebul_tot', itime_in_file, sumarea, filter, &
                  'CH4 column ebullition', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot'), a_methane_prod_tot, file_hist, &
                  'f_methane_prod_tot', itime_in_file, sumarea, filter, &
                  'CH4 column production', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot'), a_methane_oxid_tot, file_hist, &
                  'f_methane_oxid_tot', itime_in_file, sumarea, filter, &
                  'CH4 column oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot'), a_co2_decomp_tot, file_hist, &
                  'f_co2_decomp_tot', itime_in_file, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot'), a_co2_oxid_tot, file_hist, &
                  'f_co2_oxid_tot', itime_in_file, sumarea, filter, &
                  'CO2 column source from CH4 oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_aere_tot'), a_co2_aere_tot, file_hist, &
                  'f_co2_aere_tot', itime_in_file, sumarea, filter, &
                  'diagnosed CO2 column aerenchyma flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot'), a_co2_net_tot, file_hist, &
                  'f_co2_net_tot', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source from Methane module', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_totcol_methane'), a_totcol_methane, file_hist, &
                  'f_totcol_methane', itime_in_file, sumarea, filter, &
                  'total CH4 in soil column', 'mol/m2', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond'), a_grnd_methane_cond, file_hist, &
                  'f_grnd_methane_cond', itime_in_file, sumarea, filter, &
                  'boundary-layer CH4 conductance', 'm/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_conc_o2'), a_conc_o2, file_hist, &
                  'f_conc_o2', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 concentration per soil layer', 'mol/m3', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_conc_methane'), a_conc_methane, file_hist, &
                  'f_conc_methane', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per soil layer', 'mol/m3', &
                  acc_num=a_methane_acc_num)

                  ! optional microbial-pool diagnostics
                  IF (allocated(a_B_methanogen)) THEN
                     CALL write_history_variable_3d (mhist_on('f_methane_B_methanogen'), a_B_methanogen, file_hist, &
                        'f_methane_B_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'active methanogen biomass', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_B_methanotroph'), a_B_methanotroph, file_hist, &
                        'f_methane_B_methanotroph', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'active methanotroph biomass', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_B_methanogen_dormant'), a_B_methanogen_dormant, file_hist, &
                        'f_methane_B_methanogen_dormant', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'dormant methanogen biomass; active only when use_microbial_dormancy is true', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_B_methanotroph_dormant'), a_B_methanotroph_dormant, file_hist, &
                        'f_methane_B_methanotroph_dormant', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'dormant methanotroph biomass; active only when use_microbial_dormancy is true', 'gC/m3', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_T_methanogen'), a_f_T_methanogen, file_hist, &
                        'f_methane_microbe_f_T_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen microbial temperature factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_S_methanogen'), a_f_S_methanogen, file_hist, &
                        'f_methane_microbe_f_S_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen substrate factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_O2_methanogen'), a_f_O2_methanogen, file_hist, &
                        'f_methane_microbe_f_O2_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen oxygen inhibition factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_T_methanotroph'), a_f_T_methanotroph, file_hist, &
                        'f_methane_microbe_f_T_methanotroph', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanotroph microbial temperature factor', '-', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_methanogen_growth_rate'), a_methanogen_growth_rate, file_hist, &
                        'f_methane_methanogen_growth_rate', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanogen active-pool growth rate', 'day-1', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_methanotroph_growth_rate'), a_methanotroph_growth_rate, file_hist, &
                        'f_methane_methanotroph_growth_rate', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'methanotroph active-pool growth rate', 'day-1', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_microbial_prod_potential'), a_microbial_prod_potential, file_hist, &
                        'f_methane_microbial_prod_potential', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'carbon-capped microbial CH4 production potential', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_microbe)
                     CALL write_history_variable_3d (&
                        mhist_on('f_methane_microbial_oxid_potential'), a_microbial_oxid_potential, file_hist, &
                        'f_methane_microbial_oxid_potential', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'CH4/O2-capped microbial CH4 oxidation potential', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_microbe)
                  ENDIF

               ! unsat / sat
               CALL write_history_variable_2d (mhist_on('f_net_methane_unsat'), a_net_methane_unsat, file_hist, &
                  'f_net_methane_unsat', itime_in_file, sumarea, filter, &
                  'average net methane correction (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_net_methane_sat'), a_net_methane_sat, file_hist, &
                  'f_net_methane_sat', itime_in_file, sumarea, filter, &
                  'average net methane correction (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_unsat'), a_methane_prod_depth_unsat, file_hist, &
                  'f_methane_prod_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_sat'), a_methane_prod_depth_sat, file_hist, &
                  'f_methane_prod_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer (saturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_unsat'), a_methane_oxid_depth_unsat, file_hist, &
                  'f_methane_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_sat'), a_methane_oxid_depth_sat, file_hist, &
                     'f_methane_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 oxidation per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth_unsat'), a_o2_decomp_depth_unsat, file_hist, &
                     'f_o2_decomp_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'O2 consumption from decomposition/root respiration per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth_sat'), a_o2_decomp_depth_sat, file_hist, &
                        'f_o2_decomp_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 consumption from decomposition/root respiration per layer (saturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth_unsat'), a_o2_oxid_depth_unsat, file_hist, &
                        'f_o2_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 consumption from CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth_sat'), a_o2_oxid_depth_sat, file_hist, &
                        'f_o2_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 consumption from CH4 oxidation per layer (saturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_o2_aere_depth_unsat'), a_o2_aere_depth_unsat, file_hist, &
                        'f_o2_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 aerenchyma supply per layer (unsaturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2_aere_depth_sat'), a_o2_aere_depth_sat, file_hist, &
                        'f_o2_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 aerenchyma supply per layer (saturated)', 'mol/m3/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_o2stress_unsat'), a_o2stress_unsat, file_hist, &
                        'f_o2stress_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 availability stress factor (unsaturated)', '-', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_o2stress_sat'), a_o2stress_sat, file_hist, &
                        'f_o2stress_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'O2 availability stress factor (saturated)', '-', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_3d (mhist_on('f_methane_stress_unsat'), a_methane_stress_unsat, file_hist, &
                        'f_methane_stress_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'CH4 availability stress factor (unsaturated)', '-', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_3d (mhist_on('f_methane_stress_sat'), a_methane_stress_sat, file_hist, &
                        'f_methane_stress_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'CH4 availability stress factor (saturated)', '-', &
                        acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_methane_aere_depth_unsat'), a_methane_aere_depth_unsat, file_hist, &
                     'f_methane_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 aerenchyma loss per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_aere_depth_sat'), a_methane_aere_depth_sat, file_hist, &
                     'f_methane_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 aerenchyma loss per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_methane_tran_depth_unsat'), a_methane_tran_depth_unsat, file_hist, &
                     'f_methane_tran_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 transpiration loss per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_tran_depth_sat'), a_methane_tran_depth_sat, file_hist, &
                     'f_methane_tran_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 transpiration loss per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_unsat'), a_methane_ebul_depth_unsat, file_hist, &
                     'f_methane_ebul_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 ebullition loss per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_sat'), a_methane_ebul_depth_sat, file_hist, &
                     'f_methane_ebul_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CH4 ebullition loss per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_unsat'), a_co2_decomp_depth_unsat, file_hist, &
                     'f_co2_decomp_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                        'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_sat'), a_co2_decomp_depth_sat, file_hist, &
                  'f_co2_decomp_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (saturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_unsat'), a_co2_oxid_depth_unsat, file_hist, &
                  'f_co2_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CO2 from CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_sat'), a_co2_oxid_depth_sat, file_hist, &
                     'f_co2_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CO2 from CH4 oxidation per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_co2_aere_depth_unsat'), a_co2_aere_depth_unsat, file_hist, &
                     'f_co2_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CO2 aerenchyma diagnostic per layer (unsaturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_co2_aere_depth_sat'), a_co2_aere_depth_sat, file_hist, &
                     'f_co2_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'CO2 aerenchyma diagnostic per layer (saturated)', 'mol/m3/s', &
                     acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_unsat'), a_co2_decomp_tot_unsat, file_hist, &
                  'f_co2_decomp_tot_unsat', itime_in_file, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_sat'), a_co2_decomp_tot_sat, file_hist, &
                  'f_co2_decomp_tot_sat', itime_in_file, sumarea, filter, &
                     'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_unsat'), a_co2_oxid_tot_unsat, file_hist, &
                  'f_co2_oxid_tot_unsat', itime_in_file, sumarea, filter, &
                  'CO2 oxidation column source (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_sat'), a_co2_oxid_tot_sat, file_hist, &
                  'f_co2_oxid_tot_sat', itime_in_file, sumarea, filter, &
                  'CO2 oxidation column source (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_unsat'), a_co2_net_tot_unsat, file_hist, &
                  'f_co2_net_tot_unsat', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_sat'), a_co2_net_tot_sat, file_hist, &
                  'f_co2_net_tot_sat', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (&
                  mhist_on('f_methane_surf_flux_tot_unsat'), a_methane_surf_flux_tot_unsat, file_hist, &
                  'f_methane_surf_flux_tot_unsat', itime_in_file, sumarea, filter, &
                  'CH4 surface flux (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_sat'), a_methane_surf_flux_tot_sat, file_hist, &
                     'f_methane_surf_flux_tot_sat', itime_in_file, sumarea, filter, &
                     'CH4 surface flux (saturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_unsat'), a_methane_surf_aere_unsat, file_hist, &
                     'f_methane_surf_aere_unsat', itime_in_file, sumarea, filter, &
                     'CH4 column aerenchyma flux (unsaturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_sat'), a_methane_surf_aere_sat, file_hist, &
                     'f_methane_surf_aere_sat', itime_in_file, sumarea, filter, &
                     'CH4 column aerenchyma flux (saturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_unsat'), a_methane_surf_ebul_unsat, file_hist, &
                     'f_methane_surf_ebul_unsat', itime_in_file, sumarea, filter, &
                     'CH4 surface ebullition flux (unsaturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_sat'), a_methane_surf_ebul_sat, file_hist, &
                     'f_methane_surf_ebul_sat', itime_in_file, sumarea, filter, &
                     'CH4 surface ebullition flux (saturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_unsat'), a_methane_surf_diff_unsat, file_hist, &
                     'f_methane_surf_diff_unsat', itime_in_file, sumarea, filter, &
                     'CH4 diffusive surface flux (unsaturated)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_sat'), a_methane_surf_diff_sat, file_hist, &
                        'f_methane_surf_diff_sat', itime_in_file, sumarea, filter, &
                        'CH4 diffusive surface flux (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_unsat'), a_methane_ebul_tot_unsat, file_hist, &
                        'f_methane_ebul_tot_unsat', itime_in_file, sumarea, filter, &
                        'CH4 column ebullition (unsaturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_sat'), a_methane_ebul_tot_sat, file_hist, &
                        'f_methane_ebul_tot_sat', itime_in_file, sumarea, filter, &
                        'CH4 column ebullition (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_unsat'), a_methane_prod_tot_unsat, file_hist, &
                        'f_methane_prod_tot_unsat', itime_in_file, sumarea, filter, &
                        'CH4 column production (unsaturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_sat'), a_methane_prod_tot_sat, file_hist, &
                        'f_methane_prod_tot_sat', itime_in_file, sumarea, filter, &
                        'CH4 column production (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_unsat'), a_methane_oxid_tot_unsat, file_hist, &
                        'f_methane_oxid_tot_unsat', itime_in_file, sumarea, filter, &
                        'CH4 column oxidation (unsaturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_sat'), a_methane_oxid_tot_sat, file_hist, &
                        'f_methane_oxid_tot_sat', itime_in_file, sumarea, filter, &
                        'CH4 column oxidation (saturated)', 'mol/m2/s', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_totcol_methane_unsat'), a_totcol_methane_unsat, file_hist, &
                        'f_totcol_methane_unsat', itime_in_file, sumarea, filter, &
                        'total CH4 in unsaturated column', 'mol/m2', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_totcol_methane_sat'), a_totcol_methane_sat, file_hist, &
                        'f_totcol_methane_sat', itime_in_file, sumarea, filter, &
                        'total CH4 in saturated column', 'mol/m2', &
                        acc_num=a_methane_acc_num_sat)
                     CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_unsat'), a_grnd_methane_cond_unsat, file_hist, &
                        'f_grnd_methane_cond_unsat', itime_in_file, sumarea, filter, &
                        'CH4 boundary-layer conductance (unsaturated)', 'm/s', &
                        acc_num=a_methane_acc_num_unsat)
                     CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_sat'), a_grnd_methane_cond_sat, file_hist, &
                        'f_grnd_methane_cond_sat', itime_in_file, sumarea, filter, &
                        'CH4 boundary-layer conductance (saturated)', 'm/s', &
                        acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_conc_o2_unsat'), a_conc_o2_unsat, file_hist, &
                     'f_conc_o2_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'O2 concentration per layer (unsaturated)', 'mol/m3', &
                     acc_num=a_methane_acc_num_unsat)
                  CALL write_history_variable_3d (mhist_on('f_conc_o2_sat'), a_conc_o2_sat, file_hist, &
                     'f_conc_o2_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                     'O2 concentration per layer (saturated)', 'mol/m3', &
                     acc_num=a_methane_acc_num_sat)
                  CALL write_history_variable_3d (mhist_on('f_conc_methane_unsat'), a_conc_methane_unsat, file_hist, &
                     'f_conc_methane_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per layer (unsaturated)', 'mol/m3', &
                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_conc_methane_sat'), a_conc_methane_sat, file_hist, &
	                  'f_conc_methane_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 concentration per layer (saturated)', 'mol/m3', &
	                  acc_num=a_methane_acc_num_sat)

	               ! Land-area denominator diagnostics for global integrations.
	               ! active_total_without_lake = wetland + soil + rice, excluding lake.
	               ! global_total_with_lake    = wetland + soil + rice + lake.
	               IF (HistForm == 'Gridded') THEN
	                  CALL mp2g_hist%get_sumarea (sumarea, filter_all_land)
	               ENDIF
	               CALL write_history_variable_2d (&
	                  mhist_on('f_methane_surf_flux_active_total_without_lake'), &
	                  hist_ch4_active_without_lake, file_hist, &
	                  'f_methane_surf_flux_active_total_without_lake', itime_in_file, &
	                  sumarea, filter_active_without_lake, &
	                  'active CH4 total surface flux excluding lake; land-area mean contribution', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (&
	                  mhist_on('f_methane_surf_flux_global_total_with_lake'), &
	                  hist_ch4_global_with_lake, file_hist, &
	                  'f_methane_surf_flux_global_total_with_lake', itime_in_file, &
	                  sumarea, filter_all_land, &
	                  'global CH4 total surface flux including lake; land-area mean contribution', 'mol/m2/s')
	

	               ! lake diagnostics (CTSM lake CH4 path) — swap to lake-only filter
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  filter = (patchtype == 4) .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_lake'), a_methane_prod_depth_lake, file_hist, &
                  'f_methane_prod_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 production per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_lake'), a_methane_oxid_depth_lake, file_hist, &
                  'f_methane_oxid_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 oxidation per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_lake'), a_co2_decomp_depth_lake, file_hist, &
                  'f_co2_decomp_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CO2 decomp/methanogenesis per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_lake'), a_co2_oxid_depth_lake, file_hist, &
                  'f_co2_oxid_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CO2 from CH4 oxidation per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_lake'), a_methane_ebul_depth_lake, file_hist, &
                  'f_methane_ebul_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 ebullition loss per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_lake'), a_methane_surf_ebul_lake, file_hist, &
                  'f_methane_surf_ebul_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 surface ebullition flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_lake'), a_methane_surf_diff_lake, file_hist, &
                  'f_methane_surf_diff_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 surface diffusion flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_lake'), a_methane_surf_flux_tot_lake, file_hist, &
                     'f_methane_surf_flux_tot_lake', itime_in_file, sumarea, filter, &
                     'lake CH4 total surface flux (ebullition + diffusion)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num_lake)
	                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_lake'), a_methane_surf_flux_lake, file_hist, &
	                     'f_methane_surf_flux_lake', itime_in_file, sumarea, filter, &
	                     'lake CH4 total surface flux; lake-area intensive, multiply by area_lake for global total', 'mol/m2/s', &
	                     acc_num=a_methane_acc_num_lake)
                  CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_lake'), a_methane_prod_tot_lake, file_hist, &
                  'f_methane_prod_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column production', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_lake'), a_methane_oxid_tot_lake, file_hist, &
                  'f_methane_oxid_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_lake'), a_methane_ebul_tot_lake, file_hist, &
                  'f_methane_ebul_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column ebullition', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_lake'), a_co2_decomp_tot_lake, file_hist, &
                  'f_co2_decomp_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CO2 decomp/methanogenesis column source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_lake'), a_co2_oxid_tot_lake, file_hist, &
                  'f_co2_oxid_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CO2 oxidation column source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_lake'), a_co2_net_tot_lake, file_hist, &
                  'f_co2_net_tot_lake', itime_in_file, sumarea, filter, &
                  'lake net diagnosed CO2 source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_totcol_methane_lake'), a_totcol_methane_lake, file_hist, &
                  'f_totcol_methane_lake', itime_in_file, sumarea, filter, &
                  'total CH4 in lake sediment column', 'mol/m2', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_lake'), a_grnd_methane_cond_lake, file_hist, &
                  'f_grnd_methane_cond_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 boundary-layer conductance', 'm/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_conc_o2_lake'), a_conc_o2_lake, file_hist, &
                  'f_conc_o2_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake O2 concentration per sediment layer', 'mol/m3', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_conc_methane_lake'), a_conc_methane_lake, file_hist, &
                  'f_conc_methane_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 concentration per sediment layer', 'mol/m3', &
                  acc_num=a_methane_acc_num_lake)

               ! extras + microbe diagnostics — swap back to active CH4 mask
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  ! Use BgcLink helper so the history filter matches the
                  ! driver gate + rice paddy override exactly.  Same call
                  ! site pattern as MOD_Tracer_Methane_AccFlux.F90.
                  CALL methane_patch_active_mask (numpatch, DEF_METHANE%only_wetland, &
                                                  DEF_METHANE%enable_rice_paddy, &
                                                  patchtype, filter)
                  filter = filter .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_2d (mhist_on('f_forc_pmethanem'), a_forc_pmethanem, file_hist, &
                  'f_forc_pmethanem', itime_in_file, sumarea, filter, &
                  'atmospheric CH4 partial pressure forcing', 'Pa', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_3d (mhist_on('f_layer_sat_lag'), a_layer_sat_lag, file_hist, &
                  'f_layer_sat_lag', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'time-lagged saturation flag per layer', '-', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_annavg_agnpp'), a_annavg_agnpp, file_hist, &
                  'f_annavg_agnpp', itime_in_file, sumarea, filter, &
                  'annual avg above-ground NPP', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_bgnpp'), a_annavg_bgnpp, file_hist, &
                  'f_annavg_bgnpp', itime_in_file, sumarea, filter, &
                  'annual avg below-ground NPP', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_somhr'), a_annavg_somhr, file_hist, &
                  'f_annavg_somhr', itime_in_file, sumarea, filter, &
                  'annual avg SOM heterotrophic respiration', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_finrw'), a_annavg_finrw, file_hist, &
                  'f_annavg_finrw', itime_in_file, sumarea, filter, &
                  'respiration-weighted annual avg inundated fraction', '-', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_methane_dfsat_tot'), a_methane_dfsat_tot, file_hist, &
                  'f_methane_dfsat_tot', itime_in_file, sumarea, filter, &
                  'CH4 flux from decreasing inundated fraction', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_extra)
                  CALL write_history_variable_2d (mhist_on('f_f_h2osfc'), a_f_h2osfc, file_hist, &
                     'f_f_h2osfc', itime_in_file, sumarea, filter, &
                     'fraction of surface water (CLM5 microtopography)', '-', &
                     acc_num=a_methane_acc_num_extra)
                  CALL write_history_variable_2d (mhist_on('f_methane_finundated'), a_methane_finundated, file_hist, &
                     'f_methane_finundated', itime_in_file, sumarea, filter, &
                     'actual inundated fraction used by CH4 physics', '-', &
                     acc_num=a_methane_acc_num)
                  ! f_methane_soil_finundated moved to soil-only filter block below (avoid wetland zero dilution)
                  CALL write_history_variable_2d (mhist_on('f_inund_flood_patch'), a_f_inund_flood_patch, file_hist, &
                     'f_inund_flood_patch', itime_in_file, sumarea, filter, &
                     'routing-published flood inundation fraction on methane-active patches', '-', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_inund_flood_depth_patch'), a_f_inund_flood_depth_patch, file_hist, &
                     'f_inund_flood_depth_patch', itime_in_file, sumarea, filter, &
                     'routing-published floodplain water depth on methane-active patches', 'm', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_wetland_frac_patch'), a_wetland_frac_per_patch, file_hist, &
                     'f_wetland_frac_patch', itime_in_file, sumarea, filter, &
                     'static wetland fraction of parent cell for methane patch', '-', &
                     acc_num=a_methane_acc_num)
                  ! Sub-tile CH4 flux contributions: each variable is the active-CH4 area mean
                  ! of that tile's contribution; non-source patches contribute 0.
                  ! Sum (wetland+soil+rice+lake) equals f_methane_surf_flux_tot exactly.
                  ! NOT the intrinsic per-tile flux (divide by tile area fraction for that).
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_wetland'), a_methane_surf_flux_wetland, file_hist, &
                     'f_methane_surf_flux_wetland', itime_in_file, sumarea, filter, &
                     'wetland-patch contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_soil'), a_methane_surf_flux_soil, file_hist, &
                     'f_methane_surf_flux_soil', itime_in_file, sumarea, filter, &
                     'non-rice soil-patch contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_rice'), a_methane_surf_flux_rice, file_hist, &
                     'f_methane_surf_flux_rice', itime_in_file, sumarea, filter, &
                     'rice-paddy contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
                     acc_num=a_methane_acc_num)

                  IF ((p_is_worker) .and. (numpatch > 0)) THEN
                     filter = (patchtype == 0) .and. patchmask
                     IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
                  ENDIF
                  IF (HistForm == 'Gridded') THEN
                     CALL mp2g_hist%get_sumarea (sumarea, filter)
                  ENDIF
                  CALL write_history_variable_2d (mhist_on('f_methane_soil_zwt'), a_methane_soil_zwt, file_hist, &
                     'f_methane_soil_zwt', itime_in_file, sumarea, filter, &
                     'water-table depth on methane-active soil/rice patches', 'm', &
                     acc_num=a_methane_acc_num)
                  CALL write_history_variable_2d (mhist_on('f_methane_soil_finundated'), a_methane_soil_finundated, file_hist, &
                     'f_methane_soil_finundated', itime_in_file, sumarea, filter, &
                     'CH4 inundated fraction on methane-active soil/rice patches (soil-only mean, no wetland-zero dilution)', '-', &
                     acc_num=a_methane_acc_num)
               ENDIF

   END SUBROUTINE methane_reactive_history

END MODULE MOD_Tracer_Reactive_Methane_Hist
#endif
