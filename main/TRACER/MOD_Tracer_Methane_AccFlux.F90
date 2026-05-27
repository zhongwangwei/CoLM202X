#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Methane_AccFlux
!=======================================================================
! Accumulated Methane diagnostic buffers.
!=======================================================================

   USE MOD_Precision
	   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_nan
	   USE MOD_Vars_Global, only: nl_soil, spval
	   USE MOD_Tracer_Methane_Const, only: DEF_METHANE

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Explicit exports: keep module data private by default and expose
   ! only the current cross-module methane interface/state fields.
   PUBLIC :: a_B_methanogen
   PUBLIC :: a_B_methanogen_dormant
   PUBLIC :: a_B_methanotroph
   PUBLIC :: a_B_methanotroph_dormant
   PUBLIC :: a_annavg_agnpp
   PUBLIC :: a_annavg_bgnpp
   PUBLIC :: a_annavg_finrw
   PUBLIC :: a_annavg_somhr
   PUBLIC :: a_co2_aere_depth
   PUBLIC :: a_co2_aere_depth_sat
   PUBLIC :: a_co2_aere_depth_unsat
   PUBLIC :: a_co2_aere_tot
   PUBLIC :: a_co2_decomp_depth
   PUBLIC :: a_co2_decomp_depth_lake
   PUBLIC :: a_co2_decomp_depth_sat
   PUBLIC :: a_co2_decomp_depth_unsat
   PUBLIC :: a_co2_decomp_tot
   PUBLIC :: a_co2_decomp_tot_lake
   PUBLIC :: a_co2_decomp_tot_sat
   PUBLIC :: a_co2_decomp_tot_unsat
   PUBLIC :: a_co2_net_tot
   PUBLIC :: a_co2_net_tot_lake
   PUBLIC :: a_co2_net_tot_sat
   PUBLIC :: a_co2_net_tot_unsat
   PUBLIC :: a_co2_oxid_depth
   PUBLIC :: a_co2_oxid_depth_lake
   PUBLIC :: a_co2_oxid_depth_sat
   PUBLIC :: a_co2_oxid_depth_unsat
   PUBLIC :: a_co2_oxid_tot
   PUBLIC :: a_co2_oxid_tot_lake
   PUBLIC :: a_co2_oxid_tot_sat
   PUBLIC :: a_co2_oxid_tot_unsat
   PUBLIC :: a_conc_methane
   PUBLIC :: a_conc_methane_lake
   PUBLIC :: a_conc_methane_sat
   PUBLIC :: a_conc_methane_unsat
   PUBLIC :: a_conc_o2
   PUBLIC :: a_conc_o2_lake
   PUBLIC :: a_conc_o2_sat
   PUBLIC :: a_conc_o2_unsat
   PUBLIC :: a_f_O2_methanogen
   PUBLIC :: a_f_S_methanogen
   PUBLIC :: a_f_T_methanogen
   PUBLIC :: a_f_T_methanotroph
   PUBLIC :: a_f_h2osfc
   PUBLIC :: a_f_inund_flood_depth_patch
   PUBLIC :: a_f_inund_flood_patch
   PUBLIC :: a_forc_pmethanem
   PUBLIC :: a_grnd_methane_cond
   PUBLIC :: a_grnd_methane_cond_lake
   PUBLIC :: a_grnd_methane_cond_sat
   PUBLIC :: a_grnd_methane_cond_unsat
   PUBLIC :: a_layer_sat_lag
   PUBLIC :: a_methane_acc_num
   PUBLIC :: a_methane_acc_num_extra
   PUBLIC :: a_methane_acc_num_lake
   PUBLIC :: a_methane_acc_num_microbe
   PUBLIC :: a_methane_acc_num_sat
   PUBLIC :: a_methane_acc_num_unsat
   PUBLIC :: a_methane_aere_depth
   PUBLIC :: a_methane_aere_depth_sat
   PUBLIC :: a_methane_aere_depth_unsat
   PUBLIC :: a_methane_balance_residual
   PUBLIC :: a_methane_ch4_clip_credit
   PUBLIC :: a_methane_dfsat_tot
   PUBLIC :: a_methane_ebul_depth
   PUBLIC :: a_methane_ebul_depth_lake
   PUBLIC :: a_methane_ebul_depth_sat
   PUBLIC :: a_methane_ebul_depth_unsat
   PUBLIC :: a_methane_ebul_tot
   PUBLIC :: a_methane_ebul_tot_lake
   PUBLIC :: a_methane_ebul_tot_sat
   PUBLIC :: a_methane_ebul_tot_unsat
   PUBLIC :: a_methane_finundated
   PUBLIC :: a_methane_oxid_depth
   PUBLIC :: a_methane_oxid_depth_lake
   PUBLIC :: a_methane_oxid_depth_sat
   PUBLIC :: a_methane_oxid_depth_unsat
   PUBLIC :: a_methane_oxid_tot
   PUBLIC :: a_methane_oxid_tot_lake
   PUBLIC :: a_methane_oxid_tot_sat
   PUBLIC :: a_methane_oxid_tot_unsat
   PUBLIC :: a_methane_prod_depth
   PUBLIC :: a_methane_prod_depth_lake
   PUBLIC :: a_methane_prod_depth_sat
   PUBLIC :: a_methane_prod_depth_unsat
   PUBLIC :: a_methane_prod_tot
   PUBLIC :: a_methane_prod_tot_lake
   PUBLIC :: a_methane_prod_tot_sat
   PUBLIC :: a_methane_prod_tot_unsat
   PUBLIC :: a_methane_soil_finundated
   PUBLIC :: a_methane_soil_zwt
   PUBLIC :: a_methane_stress
   PUBLIC :: a_methane_stress_sat
   PUBLIC :: a_methane_stress_unsat
   PUBLIC :: a_methane_surf_aere
   PUBLIC :: a_methane_surf_aere_sat
   PUBLIC :: a_methane_surf_aere_unsat
   PUBLIC :: a_methane_surf_diff
   PUBLIC :: a_methane_surf_diff_lake
   PUBLIC :: a_methane_surf_diff_phys
   PUBLIC :: a_methane_surf_diff_sat
   PUBLIC :: a_methane_surf_diff_unsat
   PUBLIC :: a_methane_surf_ebul
   PUBLIC :: a_methane_surf_ebul_lake
   PUBLIC :: a_methane_surf_ebul_sat
   PUBLIC :: a_methane_surf_ebul_unsat
   PUBLIC :: a_methane_surf_flux_lake
   PUBLIC :: a_methane_surf_flux_rice
   PUBLIC :: a_methane_surf_flux_soil
   PUBLIC :: a_methane_surf_flux_tot
   PUBLIC :: a_methane_surf_flux_tot_lake
   PUBLIC :: a_methane_surf_flux_tot_phys
   PUBLIC :: a_methane_surf_flux_tot_sat
   PUBLIC :: a_methane_surf_flux_tot_unsat
   PUBLIC :: a_methane_surf_flux_wetland
   PUBLIC :: a_methane_tran_depth
   PUBLIC :: a_methane_tran_depth_sat
   PUBLIC :: a_methane_tran_depth_unsat
   PUBLIC :: a_methanogen_growth_rate
   PUBLIC :: a_methanotroph_growth_rate
   PUBLIC :: a_microbial_oxid_potential
   PUBLIC :: a_microbial_prod_potential
   PUBLIC :: a_net_methane
   PUBLIC :: a_net_methane_sat
   PUBLIC :: a_net_methane_unsat
   PUBLIC :: a_o2_aere_depth
   PUBLIC :: a_o2_aere_depth_sat
   PUBLIC :: a_o2_aere_depth_unsat
   PUBLIC :: a_o2_cap_gain
   PUBLIC :: a_o2_cap_loss
   PUBLIC :: a_o2_decomp_depth
   PUBLIC :: a_o2_decomp_depth_sat
   PUBLIC :: a_o2_decomp_depth_unsat
   PUBLIC :: a_o2_oxid_depth
   PUBLIC :: a_o2_oxid_depth_sat
   PUBLIC :: a_o2_oxid_depth_unsat
   PUBLIC :: a_o2stress
   PUBLIC :: a_o2stress_sat
   PUBLIC :: a_o2stress_unsat
   PUBLIC :: a_totcol_methane
   PUBLIC :: a_totcol_methane_lake
   PUBLIC :: a_totcol_methane_sat
   PUBLIC :: a_totcol_methane_unsat
   PUBLIC :: a_wetland_frac_per_patch
   PUBLIC :: accumulate_methane_fluxes
   PUBLIC :: allocate_methane_acc_fluxes
   PUBLIC :: deallocate_methane_acc_fluxes
   PUBLIC :: flush_methane_acc_fluxes
   PUBLIC :: read_methane_accflux_restart
   PUBLIC :: write_methane_accflux_restart

   !!!! ----------------------------------------------------------------
   !!!!                         sum data
   !!!! ----------------------------------------------------------------
   real(r8), allocatable :: a_net_methane           (:)
   real(r8), allocatable :: a_methane_prod_depth    (:,:)
   real(r8), allocatable :: a_o2_decomp_depth       (:,:)
   real(r8), allocatable :: a_co2_decomp_depth      (:,:)
   real(r8), allocatable :: a_methane_oxid_depth    (:,:)
   real(r8), allocatable :: a_o2_oxid_depth         (:,:)
   real(r8), allocatable :: a_co2_oxid_depth        (:,:)
   real(r8), allocatable :: a_methane_aere_depth    (:,:)
   real(r8), allocatable :: a_methane_tran_depth    (:,:)
   real(r8), allocatable :: a_o2_aere_depth         (:,:)
   real(r8), allocatable :: a_co2_aere_depth        (:,:)
   real(r8), allocatable :: a_methane_ebul_depth    (:,:)
   real(r8), allocatable :: a_o2stress              (:,:)
   real(r8), allocatable :: a_methane_stress        (:,:)
   real(r8), allocatable :: a_methane_surf_flux_tot (:)
   real(r8), allocatable :: a_methane_surf_flux_tot_phys (:)
   real(r8), allocatable :: a_methane_surf_aere     (:)
   real(r8), allocatable :: a_methane_surf_ebul     (:)
   real(r8), allocatable :: a_methane_surf_diff     (:)
   real(r8), allocatable :: a_methane_surf_diff_phys(:)
   real(r8), allocatable :: a_methane_balance_residual (:)
   real(r8), allocatable :: a_methane_ch4_clip_credit (:)
   real(r8), allocatable :: a_o2_cap_loss              (:)
   real(r8), allocatable :: a_o2_cap_gain              (:)
   real(r8), allocatable :: a_methane_ebul_tot      (:)
   real(r8), allocatable :: a_methane_prod_tot      (:)
   real(r8), allocatable :: a_methane_oxid_tot      (:)
   real(r8), allocatable :: a_co2_decomp_tot        (:)
   real(r8), allocatable :: a_co2_oxid_tot          (:)
   real(r8), allocatable :: a_co2_aere_tot          (:)
   real(r8), allocatable :: a_co2_net_tot           (:)
   real(r8), allocatable :: a_totcol_methane        (:)
   real(r8), allocatable :: a_grnd_methane_cond     (:)
   real(r8), allocatable :: a_conc_o2               (:,:)
   real(r8), allocatable :: a_conc_methane          (:,:)

   !!!! ----------------------------------------------------------------
   !!!!                 sum data (unsaturated / saturated)
   !!!! ----------------------------------------------------------------
   real(r8), allocatable :: a_net_methane_unsat           (:)
   real(r8), allocatable :: a_net_methane_sat             (:)
   real(r8), allocatable :: a_methane_prod_depth_unsat    (:,:)
   real(r8), allocatable :: a_methane_prod_depth_sat      (:,:)
   real(r8), allocatable :: a_o2_decomp_depth_unsat       (:,:)
   real(r8), allocatable :: a_o2_decomp_depth_sat         (:,:)
   real(r8), allocatable :: a_co2_decomp_depth_unsat      (:,:)
   real(r8), allocatable :: a_co2_decomp_depth_sat        (:,:)
   real(r8), allocatable :: a_methane_oxid_depth_unsat    (:,:)
   real(r8), allocatable :: a_methane_oxid_depth_sat      (:,:)
   real(r8), allocatable :: a_o2_oxid_depth_unsat         (:,:)
   real(r8), allocatable :: a_o2_oxid_depth_sat           (:,:)
   real(r8), allocatable :: a_co2_oxid_depth_unsat        (:,:)
   real(r8), allocatable :: a_co2_oxid_depth_sat          (:,:)
   real(r8), allocatable :: a_methane_aere_depth_unsat    (:,:)
   real(r8), allocatable :: a_methane_aere_depth_sat      (:,:)
   real(r8), allocatable :: a_methane_tran_depth_unsat    (:,:)
   real(r8), allocatable :: a_methane_tran_depth_sat      (:,:)
   real(r8), allocatable :: a_o2_aere_depth_unsat         (:,:)
   real(r8), allocatable :: a_o2_aere_depth_sat           (:,:)
   real(r8), allocatable :: a_co2_aere_depth_unsat        (:,:)
   real(r8), allocatable :: a_co2_aere_depth_sat          (:,:)
   real(r8), allocatable :: a_methane_ebul_depth_unsat    (:,:)
   real(r8), allocatable :: a_methane_ebul_depth_sat      (:,:)
   real(r8), allocatable :: a_o2stress_unsat              (:,:)
   real(r8), allocatable :: a_o2stress_sat                (:,:)
   real(r8), allocatable :: a_methane_stress_unsat        (:,:)
   real(r8), allocatable :: a_methane_stress_sat          (:,:)
   real(r8), allocatable :: a_methane_surf_flux_tot_unsat (:)
   real(r8), allocatable :: a_methane_surf_flux_tot_sat   (:)
   real(r8), allocatable :: a_methane_surf_aere_unsat     (:)
   real(r8), allocatable :: a_methane_surf_aere_sat       (:)
   real(r8), allocatable :: a_methane_surf_ebul_unsat     (:)
   real(r8), allocatable :: a_methane_surf_ebul_sat       (:)
   real(r8), allocatable :: a_methane_surf_diff_unsat     (:)
   real(r8), allocatable :: a_methane_surf_diff_sat       (:)
   real(r8), allocatable :: a_methane_ebul_tot_unsat      (:)
   real(r8), allocatable :: a_methane_ebul_tot_sat        (:)
   real(r8), allocatable :: a_methane_prod_tot_unsat      (:)
   real(r8), allocatable :: a_methane_prod_tot_sat        (:)
   real(r8), allocatable :: a_methane_oxid_tot_unsat      (:)
   real(r8), allocatable :: a_methane_oxid_tot_sat        (:)
   real(r8), allocatable :: a_co2_decomp_tot_unsat        (:)
   real(r8), allocatable :: a_co2_decomp_tot_sat          (:)
   real(r8), allocatable :: a_co2_oxid_tot_unsat          (:)
   real(r8), allocatable :: a_co2_oxid_tot_sat            (:)
   real(r8), allocatable :: a_co2_net_tot_unsat           (:)
   real(r8), allocatable :: a_co2_net_tot_sat             (:)
   real(r8), allocatable :: a_totcol_methane_unsat        (:)
   real(r8), allocatable :: a_totcol_methane_sat          (:)
   real(r8), allocatable :: a_grnd_methane_cond_unsat     (:)
   real(r8), allocatable :: a_grnd_methane_cond_sat       (:)
   real(r8), allocatable :: a_conc_o2_unsat               (:,:)
   real(r8), allocatable :: a_conc_o2_sat                 (:,:)
   real(r8), allocatable :: a_conc_methane_unsat          (:,:)
   real(r8), allocatable :: a_conc_methane_sat            (:,:)

   !!!! ----------------------------------------------------------------
   !!!!                         lake data
   !!!! ----------------------------------------------------------------
   real(r8), allocatable :: a_methane_prod_depth_lake (:,:)
   real(r8), allocatable :: a_methane_oxid_depth_lake (:,:)
   real(r8), allocatable :: a_methane_ebul_depth_lake (:,:)
   real(r8), allocatable :: a_co2_decomp_depth_lake  (:,:)
   real(r8), allocatable :: a_co2_oxid_depth_lake    (:,:)
   real(r8), allocatable :: a_methane_surf_ebul_lake  (:)
   real(r8), allocatable :: a_methane_surf_diff_lake  (:)
   real(r8), allocatable :: a_methane_surf_flux_tot_lake (:)
   real(r8), allocatable :: a_methane_prod_tot_lake   (:)
   real(r8), allocatable :: a_methane_oxid_tot_lake   (:)
   real(r8), allocatable :: a_methane_ebul_tot_lake   (:)
   real(r8), allocatable :: a_co2_decomp_tot_lake    (:)
   real(r8), allocatable :: a_co2_oxid_tot_lake      (:)
   real(r8), allocatable :: a_co2_net_tot_lake       (:)
   real(r8), allocatable :: a_totcol_methane_lake     (:)
   real(r8), allocatable :: a_grnd_methane_cond_lake  (:)
   real(r8), allocatable :: a_conc_o2_lake            (:,:)
   real(r8), allocatable :: a_conc_methane_lake       (:,:)

   !!!! ----------------------------------------------------------------
   !!!!                         extras
   !!!! ----------------------------------------------------------------
   real(r8), allocatable :: a_forc_pmethanem        (:)
   real(r8), allocatable :: a_layer_sat_lag         (:,:)
   real(r8), allocatable :: a_annavg_agnpp          (:)
   real(r8), allocatable :: a_annavg_bgnpp          (:)
   real(r8), allocatable :: a_annavg_somhr          (:)
	   real(r8), allocatable :: a_annavg_finrw          (:)
	   real(r8), allocatable :: a_methane_dfsat_tot     (:)
	   real(r8), allocatable :: a_f_h2osfc              (:)
	   real(r8), allocatable :: a_methane_finundated        (:)
	   real(r8), allocatable :: a_methane_soil_finundated   (:)
	   real(r8), allocatable :: a_methane_soil_zwt          (:)
	   real(r8), allocatable :: a_f_inund_flood_patch       (:)
	   real(r8), allocatable :: a_f_inund_flood_depth_patch (:)
	   real(r8), allocatable :: a_wetland_frac_per_patch    (:)
	   real(r8), allocatable :: a_methane_surf_flux_wetland (:)
	   real(r8), allocatable :: a_methane_surf_flux_soil    (:)
	   real(r8), allocatable :: a_methane_surf_flux_lake    (:)
	   real(r8), allocatable :: a_methane_surf_flux_rice    (:)

   !!!! ----------------------------------------------------------------
   !!!!                 optional microbial-pool diagnostics
   !!!! ----------------------------------------------------------------
   real(r8), allocatable :: a_B_methanogen              (:,:)
   real(r8), allocatable :: a_B_methanotroph            (:,:)
   real(r8), allocatable :: a_B_methanogen_dormant      (:,:)
   real(r8), allocatable :: a_B_methanotroph_dormant    (:,:)
   real(r8), allocatable :: a_f_T_methanogen            (:,:)
   real(r8), allocatable :: a_f_S_methanogen            (:,:)
   real(r8), allocatable :: a_f_O2_methanogen           (:,:)
   real(r8), allocatable :: a_f_T_methanotroph          (:,:)
   real(r8), allocatable :: a_methanogen_growth_rate    (:,:)
   real(r8), allocatable :: a_methanotroph_growth_rate  (:,:)
   real(r8), allocatable :: a_microbial_prod_potential  (:,:)
   real(r8), allocatable :: a_microbial_oxid_potential  (:,:)

   ! Per-patch sample counts used by history writes.  These prevent optional
   ! or mutually exclusive diagnostics (unsat/sat/lake/microbial) from being
   ! divided by the global accumulator count when they were spval for part of
   ! the history window.
   real(r8), allocatable :: a_methane_acc_num       (:)
   real(r8), allocatable :: a_methane_acc_num_unsat (:)
   real(r8), allocatable :: a_methane_acc_num_sat   (:)
   real(r8), allocatable :: a_methane_acc_num_lake  (:)
   real(r8), allocatable :: a_methane_acc_num_extra (:)
   real(r8), allocatable :: a_methane_acc_num_microbe (:)


CONTAINS

   !-------------------------------------------------------------------
   SUBROUTINE allocate_methane_acc_fluxes (numpatch)
      integer, intent(in) :: numpatch

      ! sum data
      allocate (a_net_methane           (numpatch))
      allocate (a_methane_prod_depth    (nl_soil,numpatch))
      allocate (a_o2_decomp_depth       (nl_soil,numpatch))
      allocate (a_co2_decomp_depth      (nl_soil,numpatch))
      allocate (a_methane_oxid_depth    (nl_soil,numpatch))
      allocate (a_o2_oxid_depth         (nl_soil,numpatch))
      allocate (a_co2_oxid_depth        (nl_soil,numpatch))
      allocate (a_methane_aere_depth    (nl_soil,numpatch))
      allocate (a_methane_tran_depth    (nl_soil,numpatch))
      allocate (a_o2_aere_depth         (nl_soil,numpatch))
      allocate (a_co2_aere_depth        (nl_soil,numpatch))
      allocate (a_methane_ebul_depth    (nl_soil,numpatch))
      allocate (a_o2stress              (nl_soil,numpatch))
      allocate (a_methane_stress        (nl_soil,numpatch))
      allocate (a_methane_surf_flux_tot (numpatch))
      allocate (a_methane_surf_flux_tot_phys (numpatch))
      allocate (a_methane_surf_aere     (numpatch))
      allocate (a_methane_surf_ebul     (numpatch))
      allocate (a_methane_surf_diff     (numpatch))
      allocate (a_methane_surf_diff_phys(numpatch))
      allocate (a_methane_balance_residual (numpatch))
      allocate (a_methane_ch4_clip_credit (numpatch))
      allocate (a_o2_cap_loss              (numpatch))
      allocate (a_o2_cap_gain              (numpatch))
      allocate (a_methane_ebul_tot      (numpatch))
      allocate (a_methane_prod_tot      (numpatch))
      allocate (a_methane_oxid_tot      (numpatch))
      allocate (a_co2_decomp_tot        (numpatch))
      allocate (a_co2_oxid_tot          (numpatch))
      allocate (a_co2_aere_tot          (numpatch))
      allocate (a_co2_net_tot           (numpatch))
      allocate (a_totcol_methane        (numpatch))
      allocate (a_grnd_methane_cond     (numpatch))
      allocate (a_conc_o2               (nl_soil,numpatch))
      allocate (a_conc_methane          (nl_soil,numpatch))

      ! unsat/sat data
      allocate (a_net_methane_unsat           (numpatch))
      allocate (a_net_methane_sat             (numpatch))
      allocate (a_methane_prod_depth_unsat    (nl_soil,numpatch))
      allocate (a_methane_prod_depth_sat      (nl_soil,numpatch))
      allocate (a_o2_decomp_depth_unsat       (nl_soil,numpatch))
      allocate (a_o2_decomp_depth_sat         (nl_soil,numpatch))
      allocate (a_co2_decomp_depth_unsat      (nl_soil,numpatch))
      allocate (a_co2_decomp_depth_sat        (nl_soil,numpatch))
      allocate (a_methane_oxid_depth_unsat    (nl_soil,numpatch))
      allocate (a_methane_oxid_depth_sat      (nl_soil,numpatch))
      allocate (a_o2_oxid_depth_unsat         (nl_soil,numpatch))
      allocate (a_o2_oxid_depth_sat           (nl_soil,numpatch))
      allocate (a_co2_oxid_depth_unsat        (nl_soil,numpatch))
      allocate (a_co2_oxid_depth_sat          (nl_soil,numpatch))
      allocate (a_methane_aere_depth_unsat    (nl_soil,numpatch))
      allocate (a_methane_aere_depth_sat      (nl_soil,numpatch))
      allocate (a_methane_tran_depth_unsat    (nl_soil,numpatch))
      allocate (a_methane_tran_depth_sat      (nl_soil,numpatch))
      allocate (a_o2_aere_depth_unsat         (nl_soil,numpatch))
      allocate (a_o2_aere_depth_sat           (nl_soil,numpatch))
      allocate (a_co2_aere_depth_unsat        (nl_soil,numpatch))
      allocate (a_co2_aere_depth_sat          (nl_soil,numpatch))
      allocate (a_methane_ebul_depth_unsat    (nl_soil,numpatch))
      allocate (a_methane_ebul_depth_sat      (nl_soil,numpatch))
      allocate (a_o2stress_unsat              (nl_soil,numpatch))
      allocate (a_o2stress_sat                (nl_soil,numpatch))
      allocate (a_methane_stress_unsat        (nl_soil,numpatch))
      allocate (a_methane_stress_sat          (nl_soil,numpatch))
      allocate (a_methane_surf_flux_tot_unsat (numpatch))
      allocate (a_methane_surf_flux_tot_sat   (numpatch))
      allocate (a_methane_surf_aere_unsat     (numpatch))
      allocate (a_methane_surf_aere_sat       (numpatch))
      allocate (a_methane_surf_ebul_unsat     (numpatch))
      allocate (a_methane_surf_ebul_sat       (numpatch))
      allocate (a_methane_surf_diff_unsat     (numpatch))
      allocate (a_methane_surf_diff_sat       (numpatch))
      allocate (a_methane_ebul_tot_unsat      (numpatch))
      allocate (a_methane_ebul_tot_sat        (numpatch))
      allocate (a_methane_prod_tot_unsat      (numpatch))
      allocate (a_methane_prod_tot_sat        (numpatch))
      allocate (a_methane_oxid_tot_unsat      (numpatch))
      allocate (a_methane_oxid_tot_sat        (numpatch))
      allocate (a_co2_decomp_tot_unsat        (numpatch))
      allocate (a_co2_decomp_tot_sat          (numpatch))
      allocate (a_co2_oxid_tot_unsat          (numpatch))
      allocate (a_co2_oxid_tot_sat            (numpatch))
      allocate (a_co2_net_tot_unsat           (numpatch))
      allocate (a_co2_net_tot_sat             (numpatch))
      allocate (a_totcol_methane_unsat        (numpatch))
      allocate (a_totcol_methane_sat          (numpatch))
      allocate (a_grnd_methane_cond_unsat     (numpatch))
      allocate (a_grnd_methane_cond_sat       (numpatch))
      allocate (a_conc_o2_unsat               (nl_soil,numpatch))
      allocate (a_conc_o2_sat                 (nl_soil,numpatch))
      allocate (a_conc_methane_unsat          (nl_soil,numpatch))
      allocate (a_conc_methane_sat            (nl_soil,numpatch))

      ! lake data
      allocate (a_methane_prod_depth_lake (nl_soil,numpatch))
      allocate (a_methane_oxid_depth_lake (nl_soil,numpatch))
      allocate (a_methane_ebul_depth_lake (nl_soil,numpatch))
      allocate (a_co2_decomp_depth_lake  (nl_soil,numpatch))
      allocate (a_co2_oxid_depth_lake    (nl_soil,numpatch))
      allocate (a_methane_surf_ebul_lake  (numpatch))
      allocate (a_methane_surf_diff_lake  (numpatch))
      allocate (a_methane_surf_flux_tot_lake (numpatch))
      allocate (a_methane_prod_tot_lake   (numpatch))
      allocate (a_methane_oxid_tot_lake   (numpatch))
      allocate (a_methane_ebul_tot_lake   (numpatch))
      allocate (a_co2_decomp_tot_lake    (numpatch))
      allocate (a_co2_oxid_tot_lake      (numpatch))
      allocate (a_co2_net_tot_lake       (numpatch))
      allocate (a_totcol_methane_lake     (numpatch))
      allocate (a_grnd_methane_cond_lake  (numpatch))
      allocate (a_conc_o2_lake            (nl_soil,numpatch))
      allocate (a_conc_methane_lake       (nl_soil,numpatch))

      ! extras
      allocate (a_forc_pmethanem    (numpatch))
      allocate (a_layer_sat_lag     (nl_soil,numpatch))
      allocate (a_annavg_agnpp      (numpatch))
      allocate (a_annavg_bgnpp      (numpatch))
      allocate (a_annavg_somhr      (numpatch))
	      allocate (a_annavg_finrw      (numpatch))
	      allocate (a_methane_dfsat_tot (numpatch))
	      allocate (a_f_h2osfc          (numpatch))
	      allocate (a_methane_finundated        (numpatch))
	      allocate (a_methane_soil_finundated   (numpatch))
	      allocate (a_methane_soil_zwt          (numpatch))
	      allocate (a_f_inund_flood_patch       (numpatch))
	      allocate (a_f_inund_flood_depth_patch (numpatch))
	      allocate (a_wetland_frac_per_patch    (numpatch))
	      allocate (a_methane_surf_flux_wetland (numpatch))
	      allocate (a_methane_surf_flux_soil    (numpatch))
	      allocate (a_methane_surf_flux_lake    (numpatch))
	      allocate (a_methane_surf_flux_rice    (numpatch))

	      ! optional microbial-pool diagnostics
	      IF (DEF_METHANE%use_microbial_pools .and. .not. allocated(a_B_methanogen)) THEN
	         allocate (a_B_methanogen             (nl_soil,numpatch))
	         allocate (a_B_methanotroph           (nl_soil,numpatch))
	         allocate (a_B_methanogen_dormant     (nl_soil,numpatch))
	         allocate (a_B_methanotroph_dormant   (nl_soil,numpatch))
	         allocate (a_f_T_methanogen           (nl_soil,numpatch))
	         allocate (a_f_S_methanogen           (nl_soil,numpatch))
	         allocate (a_f_O2_methanogen          (nl_soil,numpatch))
	         allocate (a_f_T_methanotroph         (nl_soil,numpatch))
	         allocate (a_methanogen_growth_rate   (nl_soil,numpatch))
	         allocate (a_methanotroph_growth_rate (nl_soil,numpatch))
	         allocate (a_microbial_prod_potential (nl_soil,numpatch))
	         allocate (a_microbial_oxid_potential (nl_soil,numpatch))
	         allocate (a_methane_acc_num_microbe  (numpatch))
	      ENDIF

	      allocate (a_methane_acc_num         (numpatch))
	      allocate (a_methane_acc_num_unsat   (numpatch))
	      allocate (a_methane_acc_num_sat     (numpatch))
	      allocate (a_methane_acc_num_lake    (numpatch))
	      allocate (a_methane_acc_num_extra   (numpatch))

      CALL flush_methane_acc_fluxes ()

   END SUBROUTINE allocate_methane_acc_fluxes

   !-------------------------------------------------------------------
   SUBROUTINE flush_methane_acc_fluxes ()
      ! sum data
      a_net_methane           (:)   = 0._r8
      a_methane_prod_depth    (:,:) = 0._r8
      a_o2_decomp_depth       (:,:) = 0._r8
      a_co2_decomp_depth      (:,:) = 0._r8
      a_methane_oxid_depth    (:,:) = 0._r8
      a_o2_oxid_depth         (:,:) = 0._r8
      a_co2_oxid_depth        (:,:) = 0._r8
      a_methane_aere_depth    (:,:) = 0._r8
      a_methane_tran_depth    (:,:) = 0._r8
      a_o2_aere_depth         (:,:) = 0._r8
      a_co2_aere_depth        (:,:) = 0._r8
      a_methane_ebul_depth    (:,:) = 0._r8
      a_o2stress              (:,:) = 0._r8
      a_methane_stress        (:,:) = 0._r8
      a_methane_surf_flux_tot (:)   = 0._r8
      a_methane_surf_flux_tot_phys (:) = 0._r8
      a_methane_surf_aere     (:)   = 0._r8
      a_methane_surf_ebul     (:)   = 0._r8
      a_methane_surf_diff     (:)   = 0._r8
      a_methane_surf_diff_phys(:)   = 0._r8
      a_methane_balance_residual (:) = 0._r8
      a_methane_ch4_clip_credit (:) = 0._r8
      a_o2_cap_loss              (:) = 0._r8
      a_o2_cap_gain              (:) = 0._r8
      a_methane_ebul_tot      (:)   = 0._r8
      a_methane_prod_tot      (:)   = 0._r8
      a_methane_oxid_tot      (:)   = 0._r8
      a_co2_decomp_tot        (:)   = 0._r8
      a_co2_oxid_tot          (:)   = 0._r8
      a_co2_aere_tot          (:)   = 0._r8
      a_co2_net_tot           (:)   = 0._r8
      a_totcol_methane        (:)   = 0._r8
      a_grnd_methane_cond     (:)   = 0._r8
      a_conc_o2               (:,:) = 0._r8
      a_conc_methane          (:,:) = 0._r8

      ! unsat/sat data
      a_net_methane_unsat        (:)   = 0._r8
      a_net_methane_sat          (:)   = 0._r8
      a_methane_prod_depth_unsat (:,:) = 0._r8
      a_methane_prod_depth_sat   (:,:) = 0._r8
      a_o2_decomp_depth_unsat    (:,:) = 0._r8
      a_o2_decomp_depth_sat      (:,:) = 0._r8
      a_co2_decomp_depth_unsat (:,:) = 0._r8
      a_co2_decomp_depth_sat   (:,:) = 0._r8
      a_methane_oxid_depth_unsat (:,:) = 0._r8
      a_methane_oxid_depth_sat   (:,:) = 0._r8
      a_o2_oxid_depth_unsat      (:,:) = 0._r8
      a_o2_oxid_depth_sat        (:,:) = 0._r8
      a_co2_oxid_depth_unsat   (:,:) = 0._r8
      a_co2_oxid_depth_sat     (:,:) = 0._r8
      a_methane_aere_depth_unsat (:,:) = 0._r8
      a_methane_aere_depth_sat   (:,:) = 0._r8
      a_methane_tran_depth_unsat (:,:) = 0._r8
      a_methane_tran_depth_sat   (:,:) = 0._r8
      a_o2_aere_depth_unsat      (:,:) = 0._r8
      a_o2_aere_depth_sat        (:,:) = 0._r8
      a_co2_aere_depth_unsat   (:,:) = 0._r8
      a_co2_aere_depth_sat     (:,:) = 0._r8
      a_methane_ebul_depth_unsat (:,:) = 0._r8
      a_methane_ebul_depth_sat   (:,:) = 0._r8
      a_o2stress_unsat           (:,:) = 0._r8
      a_o2stress_sat             (:,:) = 0._r8
      a_methane_stress_unsat     (:,:) = 0._r8
      a_methane_stress_sat       (:,:) = 0._r8
      a_methane_surf_flux_tot_unsat (:) = 0._r8
      a_methane_surf_flux_tot_sat   (:) = 0._r8
      a_methane_surf_aere_unsat     (:) = 0._r8
      a_methane_surf_aere_sat       (:) = 0._r8
      a_methane_surf_ebul_unsat     (:) = 0._r8
      a_methane_surf_ebul_sat       (:) = 0._r8
      a_methane_surf_diff_unsat     (:) = 0._r8
      a_methane_surf_diff_sat       (:) = 0._r8
      a_methane_ebul_tot_unsat      (:) = 0._r8
      a_methane_ebul_tot_sat        (:) = 0._r8
      a_methane_prod_tot_unsat      (:) = 0._r8
      a_methane_prod_tot_sat        (:) = 0._r8
      a_methane_oxid_tot_unsat      (:) = 0._r8
      a_methane_oxid_tot_sat        (:) = 0._r8
      a_co2_decomp_tot_unsat      (:) = 0._r8
      a_co2_decomp_tot_sat        (:) = 0._r8
      a_co2_oxid_tot_unsat        (:) = 0._r8
      a_co2_oxid_tot_sat          (:) = 0._r8
      a_co2_net_tot_unsat         (:) = 0._r8
      a_co2_net_tot_sat           (:) = 0._r8
      a_totcol_methane_unsat        (:) = 0._r8
      a_totcol_methane_sat          (:) = 0._r8
      a_grnd_methane_cond_unsat     (:) = 0._r8
      a_grnd_methane_cond_sat       (:) = 0._r8
      a_conc_o2_unsat               (:,:) = 0._r8
      a_conc_o2_sat                 (:,:) = 0._r8
      a_conc_methane_unsat          (:,:) = 0._r8
      a_conc_methane_sat            (:,:) = 0._r8

      ! lake data
      a_methane_prod_depth_lake (:,:) = 0._r8
      a_methane_oxid_depth_lake (:,:) = 0._r8
      a_methane_ebul_depth_lake (:,:) = 0._r8
      a_co2_decomp_depth_lake  (:,:) = 0._r8
      a_co2_oxid_depth_lake    (:,:) = 0._r8
      a_methane_surf_ebul_lake  (:)   = 0._r8
      a_methane_surf_diff_lake  (:)   = 0._r8
      a_methane_surf_flux_tot_lake (:) = 0._r8
      a_methane_prod_tot_lake   (:)   = 0._r8
      a_methane_oxid_tot_lake   (:)   = 0._r8
      a_methane_ebul_tot_lake   (:)   = 0._r8
      a_co2_decomp_tot_lake    (:)   = 0._r8
      a_co2_oxid_tot_lake      (:)   = 0._r8
      a_co2_net_tot_lake       (:)   = 0._r8
      a_totcol_methane_lake     (:)   = 0._r8
      a_grnd_methane_cond_lake  (:)   = 0._r8
      a_conc_o2_lake            (:,:) = 0._r8
      a_conc_methane_lake       (:,:) = 0._r8

      ! extras
      a_forc_pmethanem    (:)   = 0._r8
      a_layer_sat_lag     (:,:) = 0._r8
      a_annavg_agnpp      (:)   = 0._r8
      a_annavg_bgnpp      (:)   = 0._r8
      a_annavg_somhr      (:)   = 0._r8
      a_annavg_finrw      (:)   = 0._r8
	      a_methane_dfsat_tot (:)   = 0._r8
		      a_f_h2osfc          (:)   = 0._r8
	      a_methane_finundated        (:) = 0._r8
	      a_methane_soil_finundated   (:) = 0._r8
	      a_methane_soil_zwt          (:) = 0._r8
	      a_f_inund_flood_patch       (:) = 0._r8
	      a_f_inund_flood_depth_patch (:) = 0._r8
	      a_wetland_frac_per_patch    (:) = 0._r8
	      a_methane_surf_flux_wetland (:) = 0._r8
	      a_methane_surf_flux_soil    (:) = 0._r8
	      a_methane_surf_flux_lake    (:) = 0._r8
	      a_methane_surf_flux_rice    (:) = 0._r8

	      ! optional microbial-pool diagnostics
	      IF (allocated(a_B_methanogen))             a_B_methanogen             (:,:) = 0._r8
	      IF (allocated(a_B_methanotroph))           a_B_methanotroph           (:,:) = 0._r8
	      IF (allocated(a_B_methanogen_dormant))     a_B_methanogen_dormant     (:,:) = 0._r8
	      IF (allocated(a_B_methanotroph_dormant))   a_B_methanotroph_dormant   (:,:) = 0._r8
	      IF (allocated(a_f_T_methanogen))           a_f_T_methanogen           (:,:) = 0._r8
	      IF (allocated(a_f_S_methanogen))           a_f_S_methanogen           (:,:) = 0._r8
	      IF (allocated(a_f_O2_methanogen))          a_f_O2_methanogen          (:,:) = 0._r8
	      IF (allocated(a_f_T_methanotroph))         a_f_T_methanotroph         (:,:) = 0._r8
	      IF (allocated(a_methanogen_growth_rate))   a_methanogen_growth_rate   (:,:) = 0._r8
	      IF (allocated(a_methanotroph_growth_rate)) a_methanotroph_growth_rate (:,:) = 0._r8
	      IF (allocated(a_microbial_prod_potential)) a_microbial_prod_potential (:,:) = 0._r8
	      IF (allocated(a_microbial_oxid_potential)) a_microbial_oxid_potential (:,:) = 0._r8

	      a_methane_acc_num         (:) = 0._r8
      a_methane_acc_num_unsat   (:) = 0._r8
      a_methane_acc_num_sat     (:) = 0._r8
      a_methane_acc_num_lake    (:) = 0._r8
      a_methane_acc_num_extra   (:) = 0._r8
	      IF (allocated(a_methane_acc_num_microbe)) a_methane_acc_num_microbe (:) = 0._r8

   END SUBROUTINE flush_methane_acc_fluxes

   !-------------------------------------------------------------------
   SUBROUTINE accumulate_methane_fluxes ()
      ! Accumulate Methane state into a_* via acc1d/acc2d.
	      USE MOD_Tracer_Methane_State,    only: &
	           net_methane, methane_prod_depth, o2_decomp_depth, co2_decomp_depth, &
	           methane_oxid_depth, o2_oxid_depth, co2_oxid_depth, &
           methane_aere_depth, methane_tran_depth, o2_aere_depth, co2_aere_depth, &
           methane_ebul_depth, o2stress, methane_stress, &
           methane_surf_flux_tot, methane_surf_flux_tot_phys, methane_surf_aere, methane_surf_ebul, methane_surf_diff, methane_surf_diff_phys, &
           methane_balance_residual, methane_ch4_clip_credit, o2_cap_loss, o2_cap_gain, &
           methane_ebul_tot, methane_prod_tot, methane_oxid_tot, &
	           co2_decomp_tot, co2_oxid_tot, co2_aere_tot, co2_net_tot, &
           totcol_methane, grnd_methane_cond, &
           conc_o2, conc_methane, &
           ! unsat/sat
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
           o2stress_unsat, o2stress_sat, &
           methane_stress_unsat, methane_stress_sat, &
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
           conc_o2_unsat, conc_o2_sat, &
           conc_methane_unsat, conc_methane_sat, &
           ! lake data
           methane_prod_depth_lake, methane_oxid_depth_lake, methane_ebul_depth_lake, &
	           co2_decomp_depth_lake, co2_oxid_depth_lake, &
           methane_surf_ebul_lake, methane_surf_diff_lake, methane_surf_flux_tot_lake, &
           methane_prod_tot_lake, methane_oxid_tot_lake, methane_ebul_tot_lake, &
	           co2_decomp_tot_lake, co2_oxid_tot_lake, co2_net_tot_lake, &
           totcol_methane_lake, grnd_methane_cond_lake, &
           conc_o2_lake, conc_methane_lake, &
	           ! extras
	           forc_pmethanem, layer_sat_lag, &
	           annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
	           methane_dfsat_tot, f_h2osfc, &
	           methane_finundated, methane_soil_finundated, methane_soil_zwt, &
	           f_inund_flood_patch, f_inund_flood_depth_patch, wetland_frac_per_patch, &
	           methane_surf_flux_wetland, methane_surf_flux_soil, methane_surf_flux_lake, &
	           methane_surf_flux_rice
      USE MOD_Tracer_Methane_Const,    only: DEF_METHANE
      USE MOD_Tracer_Methane_Microbes, only: &
           B_methanogen, B_methanotroph, B_methanogen_dormant, B_methanotroph_dormant, &
           f_T_methanogen, f_S_methanogen, f_O2_methanogen, f_T_methanotroph, &
           methanogen_growth_rate, methanotroph_growth_rate, &
           microbial_prod_potential, microbial_oxid_potential
      USE MOD_SPMD_Task, only: p_is_worker
      USE MOD_Vars_TimeInvariants, only: patchtype
      USE MOD_Tracer_Methane_BgcLink, only: methane_patch_active_mask
      ! Use private acc1d/acc2d helpers (below). Cannot USE MOD_Vars_1DAccFluxes
      ! here because MOD_Vars_1DAccFluxes itself calls back into this module
      ! (BLOCK USE in accumulate_fluxes / FLUSH_acc_fluxes), which would
      ! create a circular module dependency at compile time.

      logical, save, allocatable :: methane_active_mask(:)

      ! Only worker ranks own local land-patch physics/invariants.  IO and
      ! aggregator ranks can still have methane accumulator buffers sized to
      ! the group's patches, but no local patchtype array.  They must not
      ! build patchtype-based activity masks.
      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(patchtype)) RETURN
      IF (size(patchtype) /= size(totcol_methane)) THEN
         write(6,*) 'accumulate_methane_fluxes: worker patchtype/state size mismatch=', &
            size(patchtype), size(totcol_methane)
         CALL abort
      ENDIF

	      ! Build the methane-active patch mask once via BgcLink helper so
	      ! AccFlux / history filters match the driver gate (CoLMDRIVER.F90
	      ! ~line 380) exactly, including the rice paddy override on patchtype==0
	      ! rice tiles.
	      ! Base wetland/soil mask: (patchtype == 2) .or. ((.not. DEF_METHANE%only_wetland) .and. patchtype == 0)
	      ! save attribute amortises allocation across CN timesteps.
	      IF (allocated(methane_active_mask)) THEN
	         IF (size(methane_active_mask) /= size(patchtype)) deallocate(methane_active_mask)
	      ENDIF
	      IF (.not. allocated(methane_active_mask)) allocate(methane_active_mask(size(patchtype)))
	      methane_active_mask = (patchtype == 2) .or. ((.not. DEF_METHANE%only_wetland) .and. patchtype == 0)
	      CALL methane_patch_active_mask (size(patchtype), DEF_METHANE%only_wetland, &
	                                       DEF_METHANE%enable_rice_paddy, patchtype, &
	                                       methane_active_mask)

      ! sum data
      CALL acc1d (net_methane           , a_net_methane           )
      CALL acc2d (methane_prod_depth    , a_methane_prod_depth    )
      CALL acc2d (o2_decomp_depth       , a_o2_decomp_depth       )
      CALL acc2d (co2_decomp_depth      , a_co2_decomp_depth      )
      CALL acc2d (methane_oxid_depth    , a_methane_oxid_depth    )
      CALL acc2d (o2_oxid_depth         , a_o2_oxid_depth         )
      CALL acc2d (co2_oxid_depth        , a_co2_oxid_depth        )
      CALL acc2d (methane_aere_depth    , a_methane_aere_depth    )
      CALL acc2d (methane_tran_depth    , a_methane_tran_depth    )
      CALL acc2d (o2_aere_depth         , a_o2_aere_depth         )
      CALL acc2d (co2_aere_depth        , a_co2_aere_depth        )
      CALL acc2d (methane_ebul_depth    , a_methane_ebul_depth    )
      CALL acc2d (o2stress              , a_o2stress              )
      CALL acc2d (methane_stress        , a_methane_stress        )
      CALL acc1d (methane_surf_flux_tot , a_methane_surf_flux_tot )
      CALL acc1d (methane_surf_flux_tot_phys, a_methane_surf_flux_tot_phys)
      CALL acc1d (methane_surf_aere     , a_methane_surf_aere     )
      CALL acc1d (methane_surf_ebul     , a_methane_surf_ebul     )
      CALL acc1d (methane_surf_diff     , a_methane_surf_diff     )
      CALL acc1d (methane_surf_diff_phys, a_methane_surf_diff_phys)
      CALL acc1d (methane_balance_residual, a_methane_balance_residual)
      CALL acc1d (methane_ch4_clip_credit, a_methane_ch4_clip_credit)
      CALL acc1d (o2_cap_loss, a_o2_cap_loss)
      CALL acc1d (o2_cap_gain, a_o2_cap_gain)
      CALL acc1d (methane_ebul_tot      , a_methane_ebul_tot      )
      CALL acc1d (methane_prod_tot      , a_methane_prod_tot      )
      CALL acc1d (methane_oxid_tot      , a_methane_oxid_tot      )
      CALL acc1d (co2_decomp_tot        , a_co2_decomp_tot        )
      CALL acc1d (co2_oxid_tot          , a_co2_oxid_tot          )
      CALL acc1d (co2_aere_tot          , a_co2_aere_tot          )
      CALL acc1d (co2_net_tot           , a_co2_net_tot           )
      CALL acc1d (totcol_methane        , a_totcol_methane        )
      CALL acc1d (grnd_methane_cond     , a_grnd_methane_cond     )
      CALL acc2d (conc_o2               , a_conc_o2               )
      CALL acc2d (conc_methane          , a_conc_methane          )
      ! Match the exact driver gate at CoLMDRIVER.F90:381-384:
      !   only_wetland=.true.  -> patchtype == 2 only
      !   only_wetland=.false. -> patchtype == 2 .or. patchtype == 0
      ! Urban (1) is NEVER an active CH4 patch; the previous "patchtype<=2"
      ! mask wrongly counted urban (and counted soil even under only_wetland),
      ! diluting per-grid-cell averages at history time.
      CALL acc_count1d_masked (totcol_methane, a_methane_acc_num, &
                               methane_active_mask)

      ! unsat/sat data
      CALL acc1d (net_methane_unsat        , a_net_methane_unsat        )
      CALL acc1d (net_methane_sat          , a_net_methane_sat          )
      CALL acc2d (methane_prod_depth_unsat , a_methane_prod_depth_unsat )
      CALL acc2d (methane_prod_depth_sat   , a_methane_prod_depth_sat   )
      CALL acc2d (o2_decomp_depth_unsat    , a_o2_decomp_depth_unsat    )
      CALL acc2d (o2_decomp_depth_sat      , a_o2_decomp_depth_sat      )
      CALL acc2d (co2_decomp_depth_unsat , a_co2_decomp_depth_unsat )
      CALL acc2d (co2_decomp_depth_sat   , a_co2_decomp_depth_sat   )
      CALL acc2d (methane_oxid_depth_unsat , a_methane_oxid_depth_unsat )
      CALL acc2d (methane_oxid_depth_sat   , a_methane_oxid_depth_sat   )
      CALL acc2d (o2_oxid_depth_unsat      , a_o2_oxid_depth_unsat      )
      CALL acc2d (o2_oxid_depth_sat        , a_o2_oxid_depth_sat        )
      CALL acc2d (co2_oxid_depth_unsat   , a_co2_oxid_depth_unsat   )
      CALL acc2d (co2_oxid_depth_sat     , a_co2_oxid_depth_sat     )
      CALL acc2d (methane_aere_depth_unsat , a_methane_aere_depth_unsat )
      CALL acc2d (methane_aere_depth_sat   , a_methane_aere_depth_sat   )
      CALL acc2d (methane_tran_depth_unsat , a_methane_tran_depth_unsat )
      CALL acc2d (methane_tran_depth_sat   , a_methane_tran_depth_sat   )
      CALL acc2d (o2_aere_depth_unsat      , a_o2_aere_depth_unsat      )
      CALL acc2d (o2_aere_depth_sat        , a_o2_aere_depth_sat        )
      CALL acc2d (co2_aere_depth_unsat   , a_co2_aere_depth_unsat   )
      CALL acc2d (co2_aere_depth_sat     , a_co2_aere_depth_sat     )
      CALL acc2d (methane_ebul_depth_unsat , a_methane_ebul_depth_unsat )
      CALL acc2d (methane_ebul_depth_sat   , a_methane_ebul_depth_sat   )
      CALL acc2d (o2stress_unsat           , a_o2stress_unsat           )
      CALL acc2d (o2stress_sat             , a_o2stress_sat             )
      CALL acc2d (methane_stress_unsat     , a_methane_stress_unsat     )
      CALL acc2d (methane_stress_sat       , a_methane_stress_sat       )
      CALL acc1d (methane_surf_flux_tot_unsat , a_methane_surf_flux_tot_unsat )
      CALL acc1d (methane_surf_flux_tot_sat   , a_methane_surf_flux_tot_sat   )
      CALL acc1d (methane_surf_aere_unsat     , a_methane_surf_aere_unsat     )
      CALL acc1d (methane_surf_aere_sat       , a_methane_surf_aere_sat       )
      CALL acc1d (methane_surf_ebul_unsat     , a_methane_surf_ebul_unsat     )
      CALL acc1d (methane_surf_ebul_sat       , a_methane_surf_ebul_sat       )
      CALL acc1d (methane_surf_diff_unsat     , a_methane_surf_diff_unsat     )
      CALL acc1d (methane_surf_diff_sat       , a_methane_surf_diff_sat       )
      CALL acc1d (methane_ebul_tot_unsat      , a_methane_ebul_tot_unsat      )
      CALL acc1d (methane_ebul_tot_sat        , a_methane_ebul_tot_sat        )
      CALL acc1d (methane_prod_tot_unsat      , a_methane_prod_tot_unsat      )
      CALL acc1d (methane_prod_tot_sat        , a_methane_prod_tot_sat        )
      CALL acc1d (methane_oxid_tot_unsat      , a_methane_oxid_tot_unsat      )
      CALL acc1d (methane_oxid_tot_sat        , a_methane_oxid_tot_sat        )
      CALL acc1d (co2_decomp_tot_unsat      , a_co2_decomp_tot_unsat      )
      CALL acc1d (co2_decomp_tot_sat        , a_co2_decomp_tot_sat        )
      CALL acc1d (co2_oxid_tot_unsat        , a_co2_oxid_tot_unsat        )
      CALL acc1d (co2_oxid_tot_sat          , a_co2_oxid_tot_sat          )
      CALL acc1d (co2_net_tot_unsat         , a_co2_net_tot_unsat         )
      CALL acc1d (co2_net_tot_sat           , a_co2_net_tot_sat           )
      CALL acc1d (totcol_methane_unsat        , a_totcol_methane_unsat        )
      CALL acc1d (totcol_methane_sat          , a_totcol_methane_sat          )
      CALL acc1d (grnd_methane_cond_unsat     , a_grnd_methane_cond_unsat     )
      CALL acc1d (grnd_methane_cond_sat       , a_grnd_methane_cond_sat       )
      CALL acc2d (conc_o2_unsat               , a_conc_o2_unsat               )
      CALL acc2d (conc_o2_sat                 , a_conc_o2_sat                 )
      CALL acc2d (conc_methane_unsat          , a_conc_methane_unsat          )
      CALL acc2d (conc_methane_sat            , a_conc_methane_sat            )
      CALL acc_count1d_masked (totcol_methane_unsat, a_methane_acc_num_unsat, &
                               methane_active_mask)
      CALL acc_count1d_masked (totcol_methane_sat,   a_methane_acc_num_sat,   &
                               methane_active_mask)

      ! lake data
      CALL acc2d (methane_prod_depth_lake , a_methane_prod_depth_lake )
      CALL acc2d (methane_oxid_depth_lake , a_methane_oxid_depth_lake )
      CALL acc2d (methane_ebul_depth_lake , a_methane_ebul_depth_lake )
      CALL acc2d (co2_decomp_depth_lake  , a_co2_decomp_depth_lake  )
      CALL acc2d (co2_oxid_depth_lake    , a_co2_oxid_depth_lake    )
      CALL acc1d (methane_surf_ebul_lake  , a_methane_surf_ebul_lake  )
      CALL acc1d (methane_surf_diff_lake  , a_methane_surf_diff_lake  )
      CALL acc1d (methane_surf_flux_tot_lake, a_methane_surf_flux_tot_lake)
      CALL acc1d (methane_prod_tot_lake   , a_methane_prod_tot_lake   )
      CALL acc1d (methane_oxid_tot_lake   , a_methane_oxid_tot_lake   )
      CALL acc1d (methane_ebul_tot_lake   , a_methane_ebul_tot_lake   )
      CALL acc1d (co2_decomp_tot_lake    , a_co2_decomp_tot_lake    )
      CALL acc1d (co2_oxid_tot_lake      , a_co2_oxid_tot_lake      )
      CALL acc1d (co2_net_tot_lake       , a_co2_net_tot_lake       )
      CALL acc1d (totcol_methane_lake     , a_totcol_methane_lake     )
      CALL acc1d (grnd_methane_cond_lake  , a_grnd_methane_cond_lake  )
      CALL acc2d (conc_o2_lake            , a_conc_o2_lake            )
      CALL acc2d (conc_methane_lake       , a_conc_methane_lake       )
      ! Only active lake patches contribute to the lake count.  Without this
      ! mask, every non-lake patch with totcol_methane_lake=0 (its initial
      ! value, never touched by methane_driver) would increment the denominator
      ! and dilute per-grid-cell averages once history aggregates by area.
      ! The history filter already excludes non-lake patches at write time,
      ! but a clean lake-only count avoids relying on the filter for math.
      CALL acc_count1d_masked (totcol_methane_lake, a_methane_acc_num_lake, &
                               patchtype == 4 .and. DEF_METHANE%allowlakeprod)

      ! extras
      CALL acc1d (forc_pmethanem    , a_forc_pmethanem    )
      CALL acc2d (layer_sat_lag     , a_layer_sat_lag     )
      CALL acc1d (annavg_agnpp      , a_annavg_agnpp      )
      CALL acc1d (annavg_bgnpp      , a_annavg_bgnpp      )
      CALL acc1d (annavg_somhr      , a_annavg_somhr      )
	      CALL acc1d (annavg_finrw      , a_annavg_finrw      )
	      CALL acc1d (methane_dfsat_tot , a_methane_dfsat_tot )
	      CALL acc1d (f_h2osfc          , a_f_h2osfc          )
	      CALL acc1d (methane_finundated        , a_methane_finundated        )
	      CALL acc1d (methane_soil_finundated   , a_methane_soil_finundated   )
	      CALL acc1d (methane_soil_zwt          , a_methane_soil_zwt          )
	      CALL acc1d (f_inund_flood_patch       , a_f_inund_flood_patch       )
	      CALL acc1d (f_inund_flood_depth_patch , a_f_inund_flood_depth_patch )
	      CALL acc1d (wetland_frac_per_patch    , a_wetland_frac_per_patch    )
	      CALL acc1d (methane_surf_flux_wetland , a_methane_surf_flux_wetland )
	      CALL acc1d (methane_surf_flux_soil    , a_methane_surf_flux_soil    )
	      CALL acc1d (methane_surf_flux_lake    , a_methane_surf_flux_lake    )
	      CALL acc1d (methane_surf_flux_rice    , a_methane_surf_flux_rice    )
	      CALL acc_count1d (forc_pmethanem, a_methane_acc_num_extra)

      ! optional microbial-pool diagnostics
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL acc2d (B_methanogen             , a_B_methanogen             )
         CALL acc2d (B_methanotroph           , a_B_methanotroph           )
         CALL acc2d (B_methanogen_dormant     , a_B_methanogen_dormant     )
         CALL acc2d (B_methanotroph_dormant   , a_B_methanotroph_dormant   )
         CALL acc2d (f_T_methanogen           , a_f_T_methanogen           )
         CALL acc2d (f_S_methanogen           , a_f_S_methanogen           )
         CALL acc2d (f_O2_methanogen          , a_f_O2_methanogen          )
         CALL acc2d (f_T_methanotroph         , a_f_T_methanotroph         )
         CALL acc2d (methanogen_growth_rate   , a_methanogen_growth_rate   )
         CALL acc2d (methanotroph_growth_rate , a_methanotroph_growth_rate )
         CALL acc2d (microbial_prod_potential , a_microbial_prod_potential )
         CALL acc2d (microbial_oxid_potential , a_microbial_oxid_potential )
         CALL acc_count2d (B_methanogen, a_methane_acc_num_microbe)
      ENDIF

   END SUBROUTINE accumulate_methane_fluxes

   !-------------------------------------------------------------------
   SUBROUTINE write_methane_accflux_restart (file_restart, compress)
      USE MOD_LandPatch,    only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector
      character(len=*), intent(in) :: file_restart
      integer,          intent(in) :: compress

      IF (.not. allocated(a_methane_acc_num)) RETURN

      IF (allocated(a_net_methane)) CALL ncio_write_vector (file_restart, 'ch4_a_net_methane', &
         'patch', landpatch, a_net_methane, compress)
      IF (allocated(a_methane_prod_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_depth', &
         'soil', nl_soil, 'patch', landpatch, a_methane_prod_depth, compress)
      IF (allocated(a_o2_decomp_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_decomp_depth', &
         'soil', nl_soil, 'patch', landpatch, a_o2_decomp_depth, compress)
      IF (allocated(a_co2_decomp_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_depth', &
         'soil', nl_soil, 'patch', landpatch, a_co2_decomp_depth, compress)
      IF (allocated(a_methane_oxid_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_depth', &
         'soil', nl_soil, 'patch', landpatch, a_methane_oxid_depth, compress)
      IF (allocated(a_o2_oxid_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_oxid_depth', &
         'soil', nl_soil, 'patch', landpatch, a_o2_oxid_depth, compress)
      IF (allocated(a_co2_oxid_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_depth', &
         'soil', nl_soil, 'patch', landpatch, a_co2_oxid_depth, compress)
      IF (allocated(a_methane_aere_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_aere_depth', &
         'soil', nl_soil, 'patch', landpatch, a_methane_aere_depth, compress)
      IF (allocated(a_methane_tran_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_tran_depth', &
         'soil', nl_soil, 'patch', landpatch, a_methane_tran_depth, compress)
      IF (allocated(a_o2_aere_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_aere_depth', &
         'soil', nl_soil, 'patch', landpatch, a_o2_aere_depth, compress)
      IF (allocated(a_co2_aere_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_aere_depth', &
         'soil', nl_soil, 'patch', landpatch, a_co2_aere_depth, compress)
      IF (allocated(a_methane_ebul_depth)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_depth', &
         'soil', nl_soil, 'patch', landpatch, a_methane_ebul_depth, compress)
      IF (allocated(a_o2stress)) CALL ncio_write_vector (file_restart, 'ch4_a_o2stress', &
         'soil', nl_soil, 'patch', landpatch, a_o2stress, compress)
      IF (allocated(a_methane_stress)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_stress', &
         'soil', nl_soil, 'patch', landpatch, a_methane_stress, compress)
      IF (allocated(a_methane_surf_flux_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_tot', &
         'patch', landpatch, a_methane_surf_flux_tot, compress)
      IF (allocated(a_methane_surf_flux_tot_phys)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_tot_phys', &
         'patch', landpatch, a_methane_surf_flux_tot_phys, compress)
      IF (allocated(a_methane_surf_aere)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_aere', &
         'patch', landpatch, a_methane_surf_aere, compress)
      IF (allocated(a_methane_surf_ebul)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_ebul', &
         'patch', landpatch, a_methane_surf_ebul, compress)
      IF (allocated(a_methane_surf_diff)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_diff', &
         'patch', landpatch, a_methane_surf_diff, compress)
      IF (allocated(a_methane_surf_diff_phys)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_diff_phys', &
         'patch', landpatch, a_methane_surf_diff_phys, compress)
      IF (allocated(a_methane_balance_residual)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_balance_residual', &
         'patch', landpatch, a_methane_balance_residual, compress)
      IF (allocated(a_methane_ch4_clip_credit)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ch4_clip_credit', &
         'patch', landpatch, a_methane_ch4_clip_credit, compress)
      IF (allocated(a_o2_cap_loss)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_cap_loss', &
         'patch', landpatch, a_o2_cap_loss, compress)
      IF (allocated(a_o2_cap_gain)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_cap_gain', &
         'patch', landpatch, a_o2_cap_gain, compress)
      IF (allocated(a_methane_ebul_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_tot', &
         'patch', landpatch, a_methane_ebul_tot, compress)
      IF (allocated(a_methane_prod_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_tot', &
         'patch', landpatch, a_methane_prod_tot, compress)
      IF (allocated(a_methane_oxid_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_tot', &
         'patch', landpatch, a_methane_oxid_tot, compress)
      IF (allocated(a_co2_decomp_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_tot', &
         'patch', landpatch, a_co2_decomp_tot, compress)
      IF (allocated(a_co2_oxid_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_tot', &
         'patch', landpatch, a_co2_oxid_tot, compress)
      IF (allocated(a_co2_aere_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_aere_tot', &
         'patch', landpatch, a_co2_aere_tot, compress)
      IF (allocated(a_co2_net_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_net_tot', &
         'patch', landpatch, a_co2_net_tot, compress)
      IF (allocated(a_totcol_methane)) CALL ncio_write_vector (file_restart, 'ch4_a_totcol_methane', &
         'patch', landpatch, a_totcol_methane, compress)
      IF (allocated(a_grnd_methane_cond)) CALL ncio_write_vector (file_restart, 'ch4_a_grnd_methane_cond', &
         'patch', landpatch, a_grnd_methane_cond, compress)
      IF (allocated(a_conc_o2)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_o2', &
         'soil', nl_soil, 'patch', landpatch, a_conc_o2, compress)
      IF (allocated(a_conc_methane)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_methane', &
         'soil', nl_soil, 'patch', landpatch, a_conc_methane, compress)
      IF (allocated(a_net_methane_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_net_methane_unsat', &
         'patch', landpatch, a_net_methane_unsat, compress)
      IF (allocated(a_net_methane_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_net_methane_sat', &
         'patch', landpatch, a_net_methane_sat, compress)
      IF (allocated(a_methane_prod_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_prod_depth_unsat, compress)
      IF (allocated(a_methane_prod_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_prod_depth_sat, compress)
      IF (allocated(a_o2_decomp_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_decomp_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_o2_decomp_depth_unsat, compress)
      IF (allocated(a_o2_decomp_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_decomp_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_o2_decomp_depth_sat, compress)
      IF (allocated(a_co2_decomp_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_co2_decomp_depth_unsat, compress)
      IF (allocated(a_co2_decomp_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_co2_decomp_depth_sat, compress)
      IF (allocated(a_methane_oxid_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_oxid_depth_unsat, compress)
      IF (allocated(a_methane_oxid_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_oxid_depth_sat, compress)
      IF (allocated(a_o2_oxid_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_oxid_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_o2_oxid_depth_unsat, compress)
      IF (allocated(a_o2_oxid_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_oxid_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_o2_oxid_depth_sat, compress)
      IF (allocated(a_co2_oxid_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_co2_oxid_depth_unsat, compress)
      IF (allocated(a_co2_oxid_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_co2_oxid_depth_sat, compress)
      IF (allocated(a_methane_aere_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_aere_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_aere_depth_unsat, compress)
      IF (allocated(a_methane_aere_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_aere_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_aere_depth_sat, compress)
      IF (allocated(a_methane_tran_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_tran_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_tran_depth_unsat, compress)
      IF (allocated(a_methane_tran_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_tran_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_tran_depth_sat, compress)
      IF (allocated(a_o2_aere_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_aere_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_o2_aere_depth_unsat, compress)
      IF (allocated(a_o2_aere_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2_aere_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_o2_aere_depth_sat, compress)
      IF (allocated(a_co2_aere_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_aere_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_co2_aere_depth_unsat, compress)
      IF (allocated(a_co2_aere_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_aere_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_co2_aere_depth_sat, compress)
      IF (allocated(a_methane_ebul_depth_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_depth_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_ebul_depth_unsat, compress)
      IF (allocated(a_methane_ebul_depth_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_depth_sat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_ebul_depth_sat, compress)
      IF (allocated(a_o2stress_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2stress_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_o2stress_unsat, compress)
      IF (allocated(a_o2stress_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_o2stress_sat', &
         'soil', nl_soil, 'patch', landpatch, a_o2stress_sat, compress)
      IF (allocated(a_methane_stress_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_stress_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_stress_unsat, compress)
      IF (allocated(a_methane_stress_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_stress_sat', &
         'soil', nl_soil, 'patch', landpatch, a_methane_stress_sat, compress)
      IF (allocated(a_methane_surf_flux_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_tot_unsat', &
         'patch', landpatch, a_methane_surf_flux_tot_unsat, compress)
      IF (allocated(a_methane_surf_flux_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_tot_sat', &
         'patch', landpatch, a_methane_surf_flux_tot_sat, compress)
      IF (allocated(a_methane_surf_aere_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_aere_unsat', &
         'patch', landpatch, a_methane_surf_aere_unsat, compress)
      IF (allocated(a_methane_surf_aere_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_aere_sat', &
         'patch', landpatch, a_methane_surf_aere_sat, compress)
      IF (allocated(a_methane_surf_ebul_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_ebul_unsat', &
         'patch', landpatch, a_methane_surf_ebul_unsat, compress)
      IF (allocated(a_methane_surf_ebul_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_ebul_sat', &
         'patch', landpatch, a_methane_surf_ebul_sat, compress)
      IF (allocated(a_methane_surf_diff_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_diff_unsat', &
         'patch', landpatch, a_methane_surf_diff_unsat, compress)
      IF (allocated(a_methane_surf_diff_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_diff_sat', &
         'patch', landpatch, a_methane_surf_diff_sat, compress)
      IF (allocated(a_methane_ebul_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_tot_unsat', &
         'patch', landpatch, a_methane_ebul_tot_unsat, compress)
      IF (allocated(a_methane_ebul_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_tot_sat', &
         'patch', landpatch, a_methane_ebul_tot_sat, compress)
      IF (allocated(a_methane_prod_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_tot_unsat', &
         'patch', landpatch, a_methane_prod_tot_unsat, compress)
      IF (allocated(a_methane_prod_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_tot_sat', &
         'patch', landpatch, a_methane_prod_tot_sat, compress)
      IF (allocated(a_methane_oxid_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_tot_unsat', &
         'patch', landpatch, a_methane_oxid_tot_unsat, compress)
      IF (allocated(a_methane_oxid_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_tot_sat', &
         'patch', landpatch, a_methane_oxid_tot_sat, compress)
      IF (allocated(a_co2_decomp_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_tot_unsat', &
         'patch', landpatch, a_co2_decomp_tot_unsat, compress)
      IF (allocated(a_co2_decomp_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_tot_sat', &
         'patch', landpatch, a_co2_decomp_tot_sat, compress)
      IF (allocated(a_co2_oxid_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_tot_unsat', &
         'patch', landpatch, a_co2_oxid_tot_unsat, compress)
      IF (allocated(a_co2_oxid_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_tot_sat', &
         'patch', landpatch, a_co2_oxid_tot_sat, compress)
      IF (allocated(a_co2_net_tot_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_net_tot_unsat', &
         'patch', landpatch, a_co2_net_tot_unsat, compress)
      IF (allocated(a_co2_net_tot_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_net_tot_sat', &
         'patch', landpatch, a_co2_net_tot_sat, compress)
      IF (allocated(a_totcol_methane_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_totcol_methane_unsat', &
         'patch', landpatch, a_totcol_methane_unsat, compress)
      IF (allocated(a_totcol_methane_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_totcol_methane_sat', &
         'patch', landpatch, a_totcol_methane_sat, compress)
      IF (allocated(a_grnd_methane_cond_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_grnd_methane_cond_unsat', &
         'patch', landpatch, a_grnd_methane_cond_unsat, compress)
      IF (allocated(a_grnd_methane_cond_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_grnd_methane_cond_sat', &
         'patch', landpatch, a_grnd_methane_cond_sat, compress)
      IF (allocated(a_conc_o2_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_o2_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_conc_o2_unsat, compress)
      IF (allocated(a_conc_o2_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_o2_sat', &
         'soil', nl_soil, 'patch', landpatch, a_conc_o2_sat, compress)
      IF (allocated(a_conc_methane_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_methane_unsat', &
         'soil', nl_soil, 'patch', landpatch, a_conc_methane_unsat, compress)
      IF (allocated(a_conc_methane_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_methane_sat', &
         'soil', nl_soil, 'patch', landpatch, a_conc_methane_sat, compress)
      IF (allocated(a_methane_prod_depth_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_depth_lake', &
         'soil', nl_soil, 'patch', landpatch, a_methane_prod_depth_lake, compress)
      IF (allocated(a_methane_oxid_depth_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_depth_lake', &
         'soil', nl_soil, 'patch', landpatch, a_methane_oxid_depth_lake, compress)
      IF (allocated(a_methane_ebul_depth_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_depth_lake', &
         'soil', nl_soil, 'patch', landpatch, a_methane_ebul_depth_lake, compress)
      IF (allocated(a_co2_decomp_depth_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_depth_lake', &
         'soil', nl_soil, 'patch', landpatch, a_co2_decomp_depth_lake, compress)
      IF (allocated(a_co2_oxid_depth_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_depth_lake', &
         'soil', nl_soil, 'patch', landpatch, a_co2_oxid_depth_lake, compress)
      IF (allocated(a_methane_surf_ebul_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_ebul_lake', &
         'patch', landpatch, a_methane_surf_ebul_lake, compress)
      IF (allocated(a_methane_surf_diff_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_diff_lake', &
         'patch', landpatch, a_methane_surf_diff_lake, compress)
      IF (allocated(a_methane_surf_flux_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_tot_lake', &
         'patch', landpatch, a_methane_surf_flux_tot_lake, compress)
      IF (allocated(a_methane_prod_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_prod_tot_lake', &
         'patch', landpatch, a_methane_prod_tot_lake, compress)
      IF (allocated(a_methane_oxid_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_oxid_tot_lake', &
         'patch', landpatch, a_methane_oxid_tot_lake, compress)
      IF (allocated(a_methane_ebul_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_ebul_tot_lake', &
         'patch', landpatch, a_methane_ebul_tot_lake, compress)
      IF (allocated(a_co2_decomp_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_decomp_tot_lake', &
         'patch', landpatch, a_co2_decomp_tot_lake, compress)
      IF (allocated(a_co2_oxid_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_oxid_tot_lake', &
         'patch', landpatch, a_co2_oxid_tot_lake, compress)
      IF (allocated(a_co2_net_tot_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_co2_net_tot_lake', &
         'patch', landpatch, a_co2_net_tot_lake, compress)
      IF (allocated(a_totcol_methane_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_totcol_methane_lake', &
         'patch', landpatch, a_totcol_methane_lake, compress)
      IF (allocated(a_grnd_methane_cond_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_grnd_methane_cond_lake', &
         'patch', landpatch, a_grnd_methane_cond_lake, compress)
      IF (allocated(a_conc_o2_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_o2_lake', &
         'soil', nl_soil, 'patch', landpatch, a_conc_o2_lake, compress)
      IF (allocated(a_conc_methane_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_conc_methane_lake', &
         'soil', nl_soil, 'patch', landpatch, a_conc_methane_lake, compress)
      IF (allocated(a_forc_pmethanem)) CALL ncio_write_vector (file_restart, 'ch4_a_forc_pmethanem', &
         'patch', landpatch, a_forc_pmethanem, compress)
      IF (allocated(a_layer_sat_lag)) CALL ncio_write_vector (file_restart, 'ch4_a_layer_sat_lag', &
         'soil', nl_soil, 'patch', landpatch, a_layer_sat_lag, compress)
      IF (allocated(a_annavg_agnpp)) CALL ncio_write_vector (file_restart, 'ch4_a_annavg_agnpp', &
         'patch', landpatch, a_annavg_agnpp, compress)
      IF (allocated(a_annavg_bgnpp)) CALL ncio_write_vector (file_restart, 'ch4_a_annavg_bgnpp', &
         'patch', landpatch, a_annavg_bgnpp, compress)
      IF (allocated(a_annavg_somhr)) CALL ncio_write_vector (file_restart, 'ch4_a_annavg_somhr', &
         'patch', landpatch, a_annavg_somhr, compress)
      IF (allocated(a_annavg_finrw)) CALL ncio_write_vector (file_restart, 'ch4_a_annavg_finrw', &
         'patch', landpatch, a_annavg_finrw, compress)
      IF (allocated(a_methane_dfsat_tot)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_dfsat_tot', &
         'patch', landpatch, a_methane_dfsat_tot, compress)
      IF (allocated(a_f_h2osfc)) CALL ncio_write_vector (file_restart, 'ch4_a_f_h2osfc', &
         'patch', landpatch, a_f_h2osfc, compress)
      IF (allocated(a_methane_finundated)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_finundated', &
         'patch', landpatch, a_methane_finundated, compress)
      IF (allocated(a_methane_soil_finundated)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_soil_finundated', &
         'patch', landpatch, a_methane_soil_finundated, compress)
      IF (allocated(a_methane_soil_zwt)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_soil_zwt', &
         'patch', landpatch, a_methane_soil_zwt, compress)
      IF (allocated(a_f_inund_flood_patch)) CALL ncio_write_vector (file_restart, 'ch4_a_f_inund_flood_patch', &
         'patch', landpatch, a_f_inund_flood_patch, compress)
      IF (allocated(a_f_inund_flood_depth_patch)) CALL ncio_write_vector (file_restart, 'ch4_a_f_inund_flood_depth_patch', &
         'patch', landpatch, a_f_inund_flood_depth_patch, compress)
      IF (allocated(a_wetland_frac_per_patch)) CALL ncio_write_vector (file_restart, 'ch4_a_wetland_frac_per_patch', &
         'patch', landpatch, a_wetland_frac_per_patch, compress)
      IF (allocated(a_methane_surf_flux_wetland)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_wetland', &
         'patch', landpatch, a_methane_surf_flux_wetland, compress)
      IF (allocated(a_methane_surf_flux_soil)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_soil', &
         'patch', landpatch, a_methane_surf_flux_soil, compress)
      IF (allocated(a_methane_surf_flux_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_lake', &
         'patch', landpatch, a_methane_surf_flux_lake, compress)
      IF (allocated(a_methane_surf_flux_rice)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_surf_flux_rice', &
         'patch', landpatch, a_methane_surf_flux_rice, compress)
      IF (allocated(a_B_methanogen)) CALL ncio_write_vector (file_restart, 'ch4_a_B_methanogen', &
         'soil', nl_soil, 'patch', landpatch, a_B_methanogen, compress)
      IF (allocated(a_B_methanotroph)) CALL ncio_write_vector (file_restart, 'ch4_a_B_methanotroph', &
         'soil', nl_soil, 'patch', landpatch, a_B_methanotroph, compress)
      IF (allocated(a_B_methanogen_dormant)) CALL ncio_write_vector (file_restart, 'ch4_a_B_methanogen_dormant', &
         'soil', nl_soil, 'patch', landpatch, a_B_methanogen_dormant, compress)
      IF (allocated(a_B_methanotroph_dormant)) CALL ncio_write_vector (file_restart, 'ch4_a_B_methanotroph_dormant', &
         'soil', nl_soil, 'patch', landpatch, a_B_methanotroph_dormant, compress)
      IF (allocated(a_f_T_methanogen)) CALL ncio_write_vector (file_restart, 'ch4_a_f_T_methanogen', &
         'soil', nl_soil, 'patch', landpatch, a_f_T_methanogen, compress)
      IF (allocated(a_f_S_methanogen)) CALL ncio_write_vector (file_restart, 'ch4_a_f_S_methanogen', &
         'soil', nl_soil, 'patch', landpatch, a_f_S_methanogen, compress)
      IF (allocated(a_f_O2_methanogen)) CALL ncio_write_vector (file_restart, 'ch4_a_f_O2_methanogen', &
         'soil', nl_soil, 'patch', landpatch, a_f_O2_methanogen, compress)
      IF (allocated(a_f_T_methanotroph)) CALL ncio_write_vector (file_restart, 'ch4_a_f_T_methanotroph', &
         'soil', nl_soil, 'patch', landpatch, a_f_T_methanotroph, compress)
      IF (allocated(a_methanogen_growth_rate)) CALL ncio_write_vector (file_restart, 'ch4_a_methanogen_growth_rate', &
         'soil', nl_soil, 'patch', landpatch, a_methanogen_growth_rate, compress)
      IF (allocated(a_methanotroph_growth_rate)) CALL ncio_write_vector (file_restart, 'ch4_a_methanotroph_growth_rate', &
         'soil', nl_soil, 'patch', landpatch, a_methanotroph_growth_rate, compress)
      IF (allocated(a_microbial_prod_potential)) CALL ncio_write_vector (file_restart, 'ch4_a_microbial_prod_potential', &
         'soil', nl_soil, 'patch', landpatch, a_microbial_prod_potential, compress)
      IF (allocated(a_microbial_oxid_potential)) CALL ncio_write_vector (file_restart, 'ch4_a_microbial_oxid_potential', &
         'soil', nl_soil, 'patch', landpatch, a_microbial_oxid_potential, compress)
      IF (allocated(a_methane_acc_num)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_acc_num', &
         'patch', landpatch, a_methane_acc_num, compress)
      IF (allocated(a_methane_acc_num_unsat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_acc_num_unsat', &
         'patch', landpatch, a_methane_acc_num_unsat, compress)
      IF (allocated(a_methane_acc_num_sat)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_acc_num_sat', &
         'patch', landpatch, a_methane_acc_num_sat, compress)
      IF (allocated(a_methane_acc_num_lake)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_acc_num_lake', &
         'patch', landpatch, a_methane_acc_num_lake, compress)
      IF (allocated(a_methane_acc_num_extra)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_acc_num_extra', &
         'patch', landpatch, a_methane_acc_num_extra, compress)
      IF (allocated(a_methane_acc_num_microbe)) CALL ncio_write_vector (file_restart, 'ch4_a_methane_acc_num_microbe', &
         'patch', landpatch, a_methane_acc_num_microbe, compress)
   END SUBROUTINE write_methane_accflux_restart

   !-------------------------------------------------------------------
   SUBROUTINE read_methane_accflux_restart (file_restart)
      USE MOD_LandPatch,    only: landpatch
      USE MOD_NetCDFVector, only: ncio_read_vector
      character(len=*), intent(in) :: file_restart

      IF (.not. allocated(a_methane_acc_num)) RETURN

      IF (allocated(a_net_methane)) CALL ncio_read_vector (file_restart, 'ch4_a_net_methane', &
         landpatch, a_net_methane, defval = 0._r8)
      IF (allocated(a_methane_prod_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_depth', &
         nl_soil, landpatch, a_methane_prod_depth, defval = 0._r8)
      IF (allocated(a_o2_decomp_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_decomp_depth', &
         nl_soil, landpatch, a_o2_decomp_depth, defval = 0._r8)
      IF (allocated(a_co2_decomp_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_depth', &
         nl_soil, landpatch, a_co2_decomp_depth, defval = 0._r8)
      IF (allocated(a_methane_oxid_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_depth', &
         nl_soil, landpatch, a_methane_oxid_depth, defval = 0._r8)
      IF (allocated(a_o2_oxid_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_oxid_depth', &
         nl_soil, landpatch, a_o2_oxid_depth, defval = 0._r8)
      IF (allocated(a_co2_oxid_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_depth', &
         nl_soil, landpatch, a_co2_oxid_depth, defval = 0._r8)
      IF (allocated(a_methane_aere_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_aere_depth', &
         nl_soil, landpatch, a_methane_aere_depth, defval = 0._r8)
      IF (allocated(a_methane_tran_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_tran_depth', &
         nl_soil, landpatch, a_methane_tran_depth, defval = 0._r8)
      IF (allocated(a_o2_aere_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_aere_depth', &
         nl_soil, landpatch, a_o2_aere_depth, defval = 0._r8)
      IF (allocated(a_co2_aere_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_aere_depth', &
         nl_soil, landpatch, a_co2_aere_depth, defval = 0._r8)
      IF (allocated(a_methane_ebul_depth)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_depth', &
         nl_soil, landpatch, a_methane_ebul_depth, defval = 0._r8)
      IF (allocated(a_o2stress)) CALL ncio_read_vector (file_restart, 'ch4_a_o2stress', &
         nl_soil, landpatch, a_o2stress, defval = 0._r8)
      IF (allocated(a_methane_stress)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_stress', &
         nl_soil, landpatch, a_methane_stress, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_tot', &
         landpatch, a_methane_surf_flux_tot, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_tot_phys)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_tot_phys', &
         landpatch, a_methane_surf_flux_tot_phys, defval = 0._r8)
      IF (allocated(a_methane_surf_aere)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_aere', &
         landpatch, a_methane_surf_aere, defval = 0._r8)
      IF (allocated(a_methane_surf_ebul)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_ebul', &
         landpatch, a_methane_surf_ebul, defval = 0._r8)
      IF (allocated(a_methane_surf_diff)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_diff', &
         landpatch, a_methane_surf_diff, defval = 0._r8)
      IF (allocated(a_methane_surf_diff_phys)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_diff_phys', &
         landpatch, a_methane_surf_diff_phys, defval = 0._r8)
      IF (allocated(a_methane_balance_residual)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_balance_residual', &
         landpatch, a_methane_balance_residual, defval = 0._r8)
      IF (allocated(a_methane_ch4_clip_credit)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ch4_clip_credit', &
         landpatch, a_methane_ch4_clip_credit, defval = 0._r8)
      IF (allocated(a_o2_cap_loss)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_cap_loss', &
         landpatch, a_o2_cap_loss, defval = 0._r8)
      IF (allocated(a_o2_cap_gain)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_cap_gain', &
         landpatch, a_o2_cap_gain, defval = 0._r8)
      IF (allocated(a_methane_ebul_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_tot', &
         landpatch, a_methane_ebul_tot, defval = 0._r8)
      IF (allocated(a_methane_prod_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_tot', &
         landpatch, a_methane_prod_tot, defval = 0._r8)
      IF (allocated(a_methane_oxid_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_tot', &
         landpatch, a_methane_oxid_tot, defval = 0._r8)
      IF (allocated(a_co2_decomp_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_tot', &
         landpatch, a_co2_decomp_tot, defval = 0._r8)
      IF (allocated(a_co2_oxid_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_tot', &
         landpatch, a_co2_oxid_tot, defval = 0._r8)
      IF (allocated(a_co2_aere_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_aere_tot', &
         landpatch, a_co2_aere_tot, defval = 0._r8)
      IF (allocated(a_co2_net_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_net_tot', &
         landpatch, a_co2_net_tot, defval = 0._r8)
      IF (allocated(a_totcol_methane)) CALL ncio_read_vector (file_restart, 'ch4_a_totcol_methane', &
         landpatch, a_totcol_methane, defval = 0._r8)
      IF (allocated(a_grnd_methane_cond)) CALL ncio_read_vector (file_restart, 'ch4_a_grnd_methane_cond', &
         landpatch, a_grnd_methane_cond, defval = 0._r8)
      IF (allocated(a_conc_o2)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_o2', &
         nl_soil, landpatch, a_conc_o2, defval = 0._r8)
      IF (allocated(a_conc_methane)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_methane', &
         nl_soil, landpatch, a_conc_methane, defval = 0._r8)
      IF (allocated(a_net_methane_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_net_methane_unsat', &
         landpatch, a_net_methane_unsat, defval = 0._r8)
      IF (allocated(a_net_methane_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_net_methane_sat', &
         landpatch, a_net_methane_sat, defval = 0._r8)
      IF (allocated(a_methane_prod_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_depth_unsat', &
         nl_soil, landpatch, a_methane_prod_depth_unsat, defval = 0._r8)
      IF (allocated(a_methane_prod_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_depth_sat', &
         nl_soil, landpatch, a_methane_prod_depth_sat, defval = 0._r8)
      IF (allocated(a_o2_decomp_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_decomp_depth_unsat', &
         nl_soil, landpatch, a_o2_decomp_depth_unsat, defval = 0._r8)
      IF (allocated(a_o2_decomp_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_decomp_depth_sat', &
         nl_soil, landpatch, a_o2_decomp_depth_sat, defval = 0._r8)
      IF (allocated(a_co2_decomp_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_depth_unsat', &
         nl_soil, landpatch, a_co2_decomp_depth_unsat, defval = 0._r8)
      IF (allocated(a_co2_decomp_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_depth_sat', &
         nl_soil, landpatch, a_co2_decomp_depth_sat, defval = 0._r8)
      IF (allocated(a_methane_oxid_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_depth_unsat', &
         nl_soil, landpatch, a_methane_oxid_depth_unsat, defval = 0._r8)
      IF (allocated(a_methane_oxid_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_depth_sat', &
         nl_soil, landpatch, a_methane_oxid_depth_sat, defval = 0._r8)
      IF (allocated(a_o2_oxid_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_oxid_depth_unsat', &
         nl_soil, landpatch, a_o2_oxid_depth_unsat, defval = 0._r8)
      IF (allocated(a_o2_oxid_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_oxid_depth_sat', &
         nl_soil, landpatch, a_o2_oxid_depth_sat, defval = 0._r8)
      IF (allocated(a_co2_oxid_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_depth_unsat', &
         nl_soil, landpatch, a_co2_oxid_depth_unsat, defval = 0._r8)
      IF (allocated(a_co2_oxid_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_depth_sat', &
         nl_soil, landpatch, a_co2_oxid_depth_sat, defval = 0._r8)
      IF (allocated(a_methane_aere_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_aere_depth_unsat', &
         nl_soil, landpatch, a_methane_aere_depth_unsat, defval = 0._r8)
      IF (allocated(a_methane_aere_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_aere_depth_sat', &
         nl_soil, landpatch, a_methane_aere_depth_sat, defval = 0._r8)
      IF (allocated(a_methane_tran_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_tran_depth_unsat', &
         nl_soil, landpatch, a_methane_tran_depth_unsat, defval = 0._r8)
      IF (allocated(a_methane_tran_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_tran_depth_sat', &
         nl_soil, landpatch, a_methane_tran_depth_sat, defval = 0._r8)
      IF (allocated(a_o2_aere_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_aere_depth_unsat', &
         nl_soil, landpatch, a_o2_aere_depth_unsat, defval = 0._r8)
      IF (allocated(a_o2_aere_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2_aere_depth_sat', &
         nl_soil, landpatch, a_o2_aere_depth_sat, defval = 0._r8)
      IF (allocated(a_co2_aere_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_aere_depth_unsat', &
         nl_soil, landpatch, a_co2_aere_depth_unsat, defval = 0._r8)
      IF (allocated(a_co2_aere_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_aere_depth_sat', &
         nl_soil, landpatch, a_co2_aere_depth_sat, defval = 0._r8)
      IF (allocated(a_methane_ebul_depth_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_depth_unsat', &
         nl_soil, landpatch, a_methane_ebul_depth_unsat, defval = 0._r8)
      IF (allocated(a_methane_ebul_depth_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_depth_sat', &
         nl_soil, landpatch, a_methane_ebul_depth_sat, defval = 0._r8)
      IF (allocated(a_o2stress_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2stress_unsat', &
         nl_soil, landpatch, a_o2stress_unsat, defval = 0._r8)
      IF (allocated(a_o2stress_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_o2stress_sat', &
         nl_soil, landpatch, a_o2stress_sat, defval = 0._r8)
      IF (allocated(a_methane_stress_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_stress_unsat', &
         nl_soil, landpatch, a_methane_stress_unsat, defval = 0._r8)
      IF (allocated(a_methane_stress_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_stress_sat', &
         nl_soil, landpatch, a_methane_stress_sat, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_tot_unsat', &
         landpatch, a_methane_surf_flux_tot_unsat, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_tot_sat', &
         landpatch, a_methane_surf_flux_tot_sat, defval = 0._r8)
      IF (allocated(a_methane_surf_aere_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_aere_unsat', &
         landpatch, a_methane_surf_aere_unsat, defval = 0._r8)
      IF (allocated(a_methane_surf_aere_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_aere_sat', &
         landpatch, a_methane_surf_aere_sat, defval = 0._r8)
      IF (allocated(a_methane_surf_ebul_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_ebul_unsat', &
         landpatch, a_methane_surf_ebul_unsat, defval = 0._r8)
      IF (allocated(a_methane_surf_ebul_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_ebul_sat', &
         landpatch, a_methane_surf_ebul_sat, defval = 0._r8)
      IF (allocated(a_methane_surf_diff_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_diff_unsat', &
         landpatch, a_methane_surf_diff_unsat, defval = 0._r8)
      IF (allocated(a_methane_surf_diff_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_diff_sat', &
         landpatch, a_methane_surf_diff_sat, defval = 0._r8)
      IF (allocated(a_methane_ebul_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_tot_unsat', &
         landpatch, a_methane_ebul_tot_unsat, defval = 0._r8)
      IF (allocated(a_methane_ebul_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_tot_sat', &
         landpatch, a_methane_ebul_tot_sat, defval = 0._r8)
      IF (allocated(a_methane_prod_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_tot_unsat', &
         landpatch, a_methane_prod_tot_unsat, defval = 0._r8)
      IF (allocated(a_methane_prod_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_tot_sat', &
         landpatch, a_methane_prod_tot_sat, defval = 0._r8)
      IF (allocated(a_methane_oxid_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_tot_unsat', &
         landpatch, a_methane_oxid_tot_unsat, defval = 0._r8)
      IF (allocated(a_methane_oxid_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_tot_sat', &
         landpatch, a_methane_oxid_tot_sat, defval = 0._r8)
      IF (allocated(a_co2_decomp_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_tot_unsat', &
         landpatch, a_co2_decomp_tot_unsat, defval = 0._r8)
      IF (allocated(a_co2_decomp_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_tot_sat', &
         landpatch, a_co2_decomp_tot_sat, defval = 0._r8)
      IF (allocated(a_co2_oxid_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_tot_unsat', &
         landpatch, a_co2_oxid_tot_unsat, defval = 0._r8)
      IF (allocated(a_co2_oxid_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_tot_sat', &
         landpatch, a_co2_oxid_tot_sat, defval = 0._r8)
      IF (allocated(a_co2_net_tot_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_net_tot_unsat', &
         landpatch, a_co2_net_tot_unsat, defval = 0._r8)
      IF (allocated(a_co2_net_tot_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_net_tot_sat', &
         landpatch, a_co2_net_tot_sat, defval = 0._r8)
      IF (allocated(a_totcol_methane_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_totcol_methane_unsat', &
         landpatch, a_totcol_methane_unsat, defval = 0._r8)
      IF (allocated(a_totcol_methane_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_totcol_methane_sat', &
         landpatch, a_totcol_methane_sat, defval = 0._r8)
      IF (allocated(a_grnd_methane_cond_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_grnd_methane_cond_unsat', &
         landpatch, a_grnd_methane_cond_unsat, defval = 0._r8)
      IF (allocated(a_grnd_methane_cond_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_grnd_methane_cond_sat', &
         landpatch, a_grnd_methane_cond_sat, defval = 0._r8)
      IF (allocated(a_conc_o2_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_o2_unsat', &
         nl_soil, landpatch, a_conc_o2_unsat, defval = 0._r8)
      IF (allocated(a_conc_o2_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_o2_sat', &
         nl_soil, landpatch, a_conc_o2_sat, defval = 0._r8)
      IF (allocated(a_conc_methane_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_methane_unsat', &
         nl_soil, landpatch, a_conc_methane_unsat, defval = 0._r8)
      IF (allocated(a_conc_methane_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_methane_sat', &
         nl_soil, landpatch, a_conc_methane_sat, defval = 0._r8)
      IF (allocated(a_methane_prod_depth_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_depth_lake', &
         nl_soil, landpatch, a_methane_prod_depth_lake, defval = 0._r8)
      IF (allocated(a_methane_oxid_depth_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_depth_lake', &
         nl_soil, landpatch, a_methane_oxid_depth_lake, defval = 0._r8)
      IF (allocated(a_methane_ebul_depth_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_depth_lake', &
         nl_soil, landpatch, a_methane_ebul_depth_lake, defval = 0._r8)
      IF (allocated(a_co2_decomp_depth_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_depth_lake', &
         nl_soil, landpatch, a_co2_decomp_depth_lake, defval = 0._r8)
      IF (allocated(a_co2_oxid_depth_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_depth_lake', &
         nl_soil, landpatch, a_co2_oxid_depth_lake, defval = 0._r8)
      IF (allocated(a_methane_surf_ebul_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_ebul_lake', &
         landpatch, a_methane_surf_ebul_lake, defval = 0._r8)
      IF (allocated(a_methane_surf_diff_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_diff_lake', &
         landpatch, a_methane_surf_diff_lake, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_tot_lake', &
         landpatch, a_methane_surf_flux_tot_lake, defval = 0._r8)
      IF (allocated(a_methane_prod_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_prod_tot_lake', &
         landpatch, a_methane_prod_tot_lake, defval = 0._r8)
      IF (allocated(a_methane_oxid_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_oxid_tot_lake', &
         landpatch, a_methane_oxid_tot_lake, defval = 0._r8)
      IF (allocated(a_methane_ebul_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_ebul_tot_lake', &
         landpatch, a_methane_ebul_tot_lake, defval = 0._r8)
      IF (allocated(a_co2_decomp_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_decomp_tot_lake', &
         landpatch, a_co2_decomp_tot_lake, defval = 0._r8)
      IF (allocated(a_co2_oxid_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_oxid_tot_lake', &
         landpatch, a_co2_oxid_tot_lake, defval = 0._r8)
      IF (allocated(a_co2_net_tot_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_co2_net_tot_lake', &
         landpatch, a_co2_net_tot_lake, defval = 0._r8)
      IF (allocated(a_totcol_methane_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_totcol_methane_lake', &
         landpatch, a_totcol_methane_lake, defval = 0._r8)
      IF (allocated(a_grnd_methane_cond_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_grnd_methane_cond_lake', &
         landpatch, a_grnd_methane_cond_lake, defval = 0._r8)
      IF (allocated(a_conc_o2_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_o2_lake', &
         nl_soil, landpatch, a_conc_o2_lake, defval = 0._r8)
      IF (allocated(a_conc_methane_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_conc_methane_lake', &
         nl_soil, landpatch, a_conc_methane_lake, defval = 0._r8)
      IF (allocated(a_forc_pmethanem)) CALL ncio_read_vector (file_restart, 'ch4_a_forc_pmethanem', &
         landpatch, a_forc_pmethanem, defval = 0._r8)
      IF (allocated(a_layer_sat_lag)) CALL ncio_read_vector (file_restart, 'ch4_a_layer_sat_lag', &
         nl_soil, landpatch, a_layer_sat_lag, defval = 0._r8)
      IF (allocated(a_annavg_agnpp)) CALL ncio_read_vector (file_restart, 'ch4_a_annavg_agnpp', &
         landpatch, a_annavg_agnpp, defval = 0._r8)
      IF (allocated(a_annavg_bgnpp)) CALL ncio_read_vector (file_restart, 'ch4_a_annavg_bgnpp', &
         landpatch, a_annavg_bgnpp, defval = 0._r8)
      IF (allocated(a_annavg_somhr)) CALL ncio_read_vector (file_restart, 'ch4_a_annavg_somhr', &
         landpatch, a_annavg_somhr, defval = 0._r8)
      IF (allocated(a_annavg_finrw)) CALL ncio_read_vector (file_restart, 'ch4_a_annavg_finrw', &
         landpatch, a_annavg_finrw, defval = 0._r8)
      IF (allocated(a_methane_dfsat_tot)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_dfsat_tot', &
         landpatch, a_methane_dfsat_tot, defval = 0._r8)
      IF (allocated(a_f_h2osfc)) CALL ncio_read_vector (file_restart, 'ch4_a_f_h2osfc', &
         landpatch, a_f_h2osfc, defval = 0._r8)
      IF (allocated(a_methane_finundated)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_finundated', &
         landpatch, a_methane_finundated, defval = 0._r8)
      IF (allocated(a_methane_soil_finundated)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_soil_finundated', &
         landpatch, a_methane_soil_finundated, defval = 0._r8)
      IF (allocated(a_methane_soil_zwt)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_soil_zwt', &
         landpatch, a_methane_soil_zwt, defval = 0._r8)
      IF (allocated(a_f_inund_flood_patch)) CALL ncio_read_vector (file_restart, 'ch4_a_f_inund_flood_patch', &
         landpatch, a_f_inund_flood_patch, defval = 0._r8)
      IF (allocated(a_f_inund_flood_depth_patch)) CALL ncio_read_vector (file_restart, 'ch4_a_f_inund_flood_depth_patch', &
         landpatch, a_f_inund_flood_depth_patch, defval = 0._r8)
      IF (allocated(a_wetland_frac_per_patch)) CALL ncio_read_vector (file_restart, 'ch4_a_wetland_frac_per_patch', &
         landpatch, a_wetland_frac_per_patch, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_wetland)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_wetland', &
         landpatch, a_methane_surf_flux_wetland, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_soil)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_soil', &
         landpatch, a_methane_surf_flux_soil, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_lake', &
         landpatch, a_methane_surf_flux_lake, defval = 0._r8)
      IF (allocated(a_methane_surf_flux_rice)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_surf_flux_rice', &
         landpatch, a_methane_surf_flux_rice, defval = 0._r8)
      IF (allocated(a_B_methanogen)) CALL ncio_read_vector (file_restart, 'ch4_a_B_methanogen', &
         nl_soil, landpatch, a_B_methanogen, defval = 0._r8)
      IF (allocated(a_B_methanotroph)) CALL ncio_read_vector (file_restart, 'ch4_a_B_methanotroph', &
         nl_soil, landpatch, a_B_methanotroph, defval = 0._r8)
      IF (allocated(a_B_methanogen_dormant)) CALL ncio_read_vector (file_restart, 'ch4_a_B_methanogen_dormant', &
         nl_soil, landpatch, a_B_methanogen_dormant, defval = 0._r8)
      IF (allocated(a_B_methanotroph_dormant)) CALL ncio_read_vector (file_restart, 'ch4_a_B_methanotroph_dormant', &
         nl_soil, landpatch, a_B_methanotroph_dormant, defval = 0._r8)
      IF (allocated(a_f_T_methanogen)) CALL ncio_read_vector (file_restart, 'ch4_a_f_T_methanogen', &
         nl_soil, landpatch, a_f_T_methanogen, defval = 0._r8)
      IF (allocated(a_f_S_methanogen)) CALL ncio_read_vector (file_restart, 'ch4_a_f_S_methanogen', &
         nl_soil, landpatch, a_f_S_methanogen, defval = 0._r8)
      IF (allocated(a_f_O2_methanogen)) CALL ncio_read_vector (file_restart, 'ch4_a_f_O2_methanogen', &
         nl_soil, landpatch, a_f_O2_methanogen, defval = 0._r8)
      IF (allocated(a_f_T_methanotroph)) CALL ncio_read_vector (file_restart, 'ch4_a_f_T_methanotroph', &
         nl_soil, landpatch, a_f_T_methanotroph, defval = 0._r8)
      IF (allocated(a_methanogen_growth_rate)) CALL ncio_read_vector (file_restart, 'ch4_a_methanogen_growth_rate', &
         nl_soil, landpatch, a_methanogen_growth_rate, defval = 0._r8)
      IF (allocated(a_methanotroph_growth_rate)) CALL ncio_read_vector (file_restart, 'ch4_a_methanotroph_growth_rate', &
         nl_soil, landpatch, a_methanotroph_growth_rate, defval = 0._r8)
      IF (allocated(a_microbial_prod_potential)) CALL ncio_read_vector (file_restart, 'ch4_a_microbial_prod_potential', &
         nl_soil, landpatch, a_microbial_prod_potential, defval = 0._r8)
      IF (allocated(a_microbial_oxid_potential)) CALL ncio_read_vector (file_restart, 'ch4_a_microbial_oxid_potential', &
         nl_soil, landpatch, a_microbial_oxid_potential, defval = 0._r8)
      IF (allocated(a_methane_acc_num)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_acc_num', &
         landpatch, a_methane_acc_num, defval = 0._r8)
      IF (allocated(a_methane_acc_num_unsat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_acc_num_unsat', &
         landpatch, a_methane_acc_num_unsat, defval = 0._r8)
      IF (allocated(a_methane_acc_num_sat)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_acc_num_sat', &
         landpatch, a_methane_acc_num_sat, defval = 0._r8)
      IF (allocated(a_methane_acc_num_lake)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_acc_num_lake', &
         landpatch, a_methane_acc_num_lake, defval = 0._r8)
      IF (allocated(a_methane_acc_num_extra)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_acc_num_extra', &
         landpatch, a_methane_acc_num_extra, defval = 0._r8)
      IF (allocated(a_methane_acc_num_microbe)) CALL ncio_read_vector (file_restart, 'ch4_a_methane_acc_num_microbe', &
         landpatch, a_methane_acc_num_microbe, defval = 0._r8)
   END SUBROUTINE read_methane_accflux_restart


   !-------------------------------------------------------------------
   SUBROUTINE deallocate_methane_acc_fluxes ()
      ! sum data
      IF (allocated(a_net_methane          )) deallocate (a_net_methane          )
      IF (allocated(a_methane_prod_depth   )) deallocate (a_methane_prod_depth   )
      IF (allocated(a_o2_decomp_depth      )) deallocate (a_o2_decomp_depth      )
      IF (allocated(a_co2_decomp_depth     )) deallocate (a_co2_decomp_depth     )
      IF (allocated(a_methane_oxid_depth   )) deallocate (a_methane_oxid_depth   )
      IF (allocated(a_o2_oxid_depth        )) deallocate (a_o2_oxid_depth        )
      IF (allocated(a_co2_oxid_depth       )) deallocate (a_co2_oxid_depth       )
      IF (allocated(a_methane_aere_depth   )) deallocate (a_methane_aere_depth   )
      IF (allocated(a_methane_tran_depth   )) deallocate (a_methane_tran_depth   )
      IF (allocated(a_o2_aere_depth        )) deallocate (a_o2_aere_depth        )
      IF (allocated(a_co2_aere_depth       )) deallocate (a_co2_aere_depth       )
      IF (allocated(a_methane_ebul_depth   )) deallocate (a_methane_ebul_depth   )
      IF (allocated(a_o2stress             )) deallocate (a_o2stress             )
      IF (allocated(a_methane_stress       )) deallocate (a_methane_stress       )
      IF (allocated(a_methane_surf_flux_tot)) deallocate (a_methane_surf_flux_tot)
      IF (allocated(a_methane_surf_flux_tot_phys)) deallocate (a_methane_surf_flux_tot_phys)
      IF (allocated(a_methane_surf_aere    )) deallocate (a_methane_surf_aere    )
      IF (allocated(a_methane_surf_ebul    )) deallocate (a_methane_surf_ebul    )
      IF (allocated(a_methane_surf_diff    )) deallocate (a_methane_surf_diff    )
      IF (allocated(a_methane_surf_diff_phys)) deallocate (a_methane_surf_diff_phys)
      IF (allocated(a_methane_balance_residual)) deallocate (a_methane_balance_residual)
      IF (allocated(a_methane_ch4_clip_credit)) deallocate (a_methane_ch4_clip_credit)
      IF (allocated(a_o2_cap_loss)) deallocate (a_o2_cap_loss)
      IF (allocated(a_o2_cap_gain)) deallocate (a_o2_cap_gain)
      IF (allocated(a_methane_ebul_tot     )) deallocate (a_methane_ebul_tot     )
      IF (allocated(a_methane_prod_tot     )) deallocate (a_methane_prod_tot     )
      IF (allocated(a_methane_oxid_tot     )) deallocate (a_methane_oxid_tot     )
      IF (allocated(a_co2_decomp_tot       )) deallocate (a_co2_decomp_tot       )
      IF (allocated(a_co2_oxid_tot         )) deallocate (a_co2_oxid_tot         )
      IF (allocated(a_co2_aere_tot         )) deallocate (a_co2_aere_tot         )
      IF (allocated(a_co2_net_tot          )) deallocate (a_co2_net_tot          )
      IF (allocated(a_totcol_methane       )) deallocate (a_totcol_methane       )
      IF (allocated(a_grnd_methane_cond    )) deallocate (a_grnd_methane_cond    )
      IF (allocated(a_conc_o2              )) deallocate (a_conc_o2              )
      IF (allocated(a_conc_methane         )) deallocate (a_conc_methane         )

      ! unsat/sat data
      IF (allocated(a_net_methane_unsat       )) deallocate (a_net_methane_unsat       )
      IF (allocated(a_net_methane_sat         )) deallocate (a_net_methane_sat         )
      IF (allocated(a_methane_prod_depth_unsat)) deallocate (a_methane_prod_depth_unsat)
      IF (allocated(a_methane_prod_depth_sat  )) deallocate (a_methane_prod_depth_sat  )
      IF (allocated(a_o2_decomp_depth_unsat   )) deallocate (a_o2_decomp_depth_unsat   )
      IF (allocated(a_o2_decomp_depth_sat     )) deallocate (a_o2_decomp_depth_sat     )
      IF (allocated(a_co2_decomp_depth_unsat)) deallocate (a_co2_decomp_depth_unsat)
      IF (allocated(a_co2_decomp_depth_sat  )) deallocate (a_co2_decomp_depth_sat  )
      IF (allocated(a_methane_oxid_depth_unsat)) deallocate (a_methane_oxid_depth_unsat)
      IF (allocated(a_methane_oxid_depth_sat  )) deallocate (a_methane_oxid_depth_sat  )
      IF (allocated(a_o2_oxid_depth_unsat     )) deallocate (a_o2_oxid_depth_unsat     )
      IF (allocated(a_o2_oxid_depth_sat       )) deallocate (a_o2_oxid_depth_sat       )
      IF (allocated(a_co2_oxid_depth_unsat  )) deallocate (a_co2_oxid_depth_unsat  )
      IF (allocated(a_co2_oxid_depth_sat    )) deallocate (a_co2_oxid_depth_sat    )
      IF (allocated(a_methane_aere_depth_unsat)) deallocate (a_methane_aere_depth_unsat)
      IF (allocated(a_methane_aere_depth_sat  )) deallocate (a_methane_aere_depth_sat  )
      IF (allocated(a_methane_tran_depth_unsat)) deallocate (a_methane_tran_depth_unsat)
      IF (allocated(a_methane_tran_depth_sat  )) deallocate (a_methane_tran_depth_sat  )
      IF (allocated(a_o2_aere_depth_unsat     )) deallocate (a_o2_aere_depth_unsat     )
      IF (allocated(a_o2_aere_depth_sat       )) deallocate (a_o2_aere_depth_sat       )
      IF (allocated(a_co2_aere_depth_unsat  )) deallocate (a_co2_aere_depth_unsat  )
      IF (allocated(a_co2_aere_depth_sat    )) deallocate (a_co2_aere_depth_sat    )
      IF (allocated(a_methane_ebul_depth_unsat)) deallocate (a_methane_ebul_depth_unsat)
      IF (allocated(a_methane_ebul_depth_sat  )) deallocate (a_methane_ebul_depth_sat  )
      IF (allocated(a_o2stress_unsat          )) deallocate (a_o2stress_unsat          )
      IF (allocated(a_o2stress_sat            )) deallocate (a_o2stress_sat            )
      IF (allocated(a_methane_stress_unsat    )) deallocate (a_methane_stress_unsat    )
      IF (allocated(a_methane_stress_sat      )) deallocate (a_methane_stress_sat      )
      IF (allocated(a_methane_surf_flux_tot_unsat)) deallocate (a_methane_surf_flux_tot_unsat)
      IF (allocated(a_methane_surf_flux_tot_sat  )) deallocate (a_methane_surf_flux_tot_sat  )
      IF (allocated(a_methane_surf_aere_unsat )) deallocate (a_methane_surf_aere_unsat )
      IF (allocated(a_methane_surf_aere_sat   )) deallocate (a_methane_surf_aere_sat   )
      IF (allocated(a_methane_surf_ebul_unsat )) deallocate (a_methane_surf_ebul_unsat )
      IF (allocated(a_methane_surf_ebul_sat   )) deallocate (a_methane_surf_ebul_sat   )
      IF (allocated(a_methane_surf_diff_unsat )) deallocate (a_methane_surf_diff_unsat )
      IF (allocated(a_methane_surf_diff_sat   )) deallocate (a_methane_surf_diff_sat   )
      IF (allocated(a_methane_ebul_tot_unsat  )) deallocate (a_methane_ebul_tot_unsat  )
      IF (allocated(a_methane_ebul_tot_sat    )) deallocate (a_methane_ebul_tot_sat    )
      IF (allocated(a_methane_prod_tot_unsat  )) deallocate (a_methane_prod_tot_unsat  )
      IF (allocated(a_methane_prod_tot_sat    )) deallocate (a_methane_prod_tot_sat    )
      IF (allocated(a_methane_oxid_tot_unsat  )) deallocate (a_methane_oxid_tot_unsat  )
      IF (allocated(a_methane_oxid_tot_sat    )) deallocate (a_methane_oxid_tot_sat    )
      IF (allocated(a_co2_decomp_tot_unsat  )) deallocate (a_co2_decomp_tot_unsat  )
      IF (allocated(a_co2_decomp_tot_sat    )) deallocate (a_co2_decomp_tot_sat    )
      IF (allocated(a_co2_oxid_tot_unsat    )) deallocate (a_co2_oxid_tot_unsat    )
      IF (allocated(a_co2_oxid_tot_sat      )) deallocate (a_co2_oxid_tot_sat      )
      IF (allocated(a_co2_net_tot_unsat     )) deallocate (a_co2_net_tot_unsat     )
      IF (allocated(a_co2_net_tot_sat       )) deallocate (a_co2_net_tot_sat       )
      IF (allocated(a_totcol_methane_unsat    )) deallocate (a_totcol_methane_unsat    )
      IF (allocated(a_totcol_methane_sat      )) deallocate (a_totcol_methane_sat      )
      IF (allocated(a_grnd_methane_cond_unsat )) deallocate (a_grnd_methane_cond_unsat )
      IF (allocated(a_grnd_methane_cond_sat   )) deallocate (a_grnd_methane_cond_sat   )
      IF (allocated(a_conc_o2_unsat           )) deallocate (a_conc_o2_unsat           )
      IF (allocated(a_conc_o2_sat             )) deallocate (a_conc_o2_sat             )
      IF (allocated(a_conc_methane_unsat      )) deallocate (a_conc_methane_unsat      )
      IF (allocated(a_conc_methane_sat        )) deallocate (a_conc_methane_sat        )

      ! lake data
      IF (allocated(a_methane_prod_depth_lake)) deallocate (a_methane_prod_depth_lake)
      IF (allocated(a_methane_oxid_depth_lake)) deallocate (a_methane_oxid_depth_lake)
      IF (allocated(a_methane_ebul_depth_lake)) deallocate (a_methane_ebul_depth_lake)
      IF (allocated(a_co2_decomp_depth_lake )) deallocate (a_co2_decomp_depth_lake )
      IF (allocated(a_co2_oxid_depth_lake   )) deallocate (a_co2_oxid_depth_lake   )
      IF (allocated(a_methane_surf_ebul_lake )) deallocate (a_methane_surf_ebul_lake )
      IF (allocated(a_methane_surf_diff_lake )) deallocate (a_methane_surf_diff_lake )
      IF (allocated(a_methane_surf_flux_tot_lake)) deallocate (a_methane_surf_flux_tot_lake)
      IF (allocated(a_methane_prod_tot_lake  )) deallocate (a_methane_prod_tot_lake  )
      IF (allocated(a_methane_oxid_tot_lake  )) deallocate (a_methane_oxid_tot_lake  )
      IF (allocated(a_methane_ebul_tot_lake  )) deallocate (a_methane_ebul_tot_lake  )
      IF (allocated(a_co2_decomp_tot_lake   )) deallocate (a_co2_decomp_tot_lake   )
      IF (allocated(a_co2_oxid_tot_lake     )) deallocate (a_co2_oxid_tot_lake     )
      IF (allocated(a_co2_net_tot_lake      )) deallocate (a_co2_net_tot_lake      )
      IF (allocated(a_totcol_methane_lake    )) deallocate (a_totcol_methane_lake    )
      IF (allocated(a_grnd_methane_cond_lake )) deallocate (a_grnd_methane_cond_lake )
      IF (allocated(a_conc_o2_lake           )) deallocate (a_conc_o2_lake           )
      IF (allocated(a_conc_methane_lake      )) deallocate (a_conc_methane_lake      )

      ! extras
      IF (allocated(a_forc_pmethanem   )) deallocate (a_forc_pmethanem   )
      IF (allocated(a_layer_sat_lag    )) deallocate (a_layer_sat_lag    )
      IF (allocated(a_annavg_agnpp     )) deallocate (a_annavg_agnpp     )
      IF (allocated(a_annavg_bgnpp     )) deallocate (a_annavg_bgnpp     )
      IF (allocated(a_annavg_somhr     )) deallocate (a_annavg_somhr     )
	      IF (allocated(a_annavg_finrw     )) deallocate (a_annavg_finrw     )
	      IF (allocated(a_methane_dfsat_tot)) deallocate (a_methane_dfsat_tot)
	      IF (allocated(a_f_h2osfc         )) deallocate (a_f_h2osfc         )
	      IF (allocated(a_methane_finundated       )) deallocate (a_methane_finundated       )
	      IF (allocated(a_methane_soil_finundated  )) deallocate (a_methane_soil_finundated  )
	      IF (allocated(a_methane_soil_zwt         )) deallocate (a_methane_soil_zwt         )
	      IF (allocated(a_f_inund_flood_patch      )) deallocate (a_f_inund_flood_patch      )
	      IF (allocated(a_f_inund_flood_depth_patch)) deallocate (a_f_inund_flood_depth_patch)
	      IF (allocated(a_wetland_frac_per_patch   )) deallocate (a_wetland_frac_per_patch   )
	      IF (allocated(a_methane_surf_flux_wetland)) deallocate (a_methane_surf_flux_wetland)
	      IF (allocated(a_methane_surf_flux_soil   )) deallocate (a_methane_surf_flux_soil   )
	      IF (allocated(a_methane_surf_flux_lake   )) deallocate (a_methane_surf_flux_lake   )
	      IF (allocated(a_methane_surf_flux_rice   )) deallocate (a_methane_surf_flux_rice   )

      ! optional microbial-pool diagnostics
      IF (allocated(a_B_methanogen            )) deallocate (a_B_methanogen            )
      IF (allocated(a_B_methanotroph          )) deallocate (a_B_methanotroph          )
      IF (allocated(a_B_methanogen_dormant    )) deallocate (a_B_methanogen_dormant    )
      IF (allocated(a_B_methanotroph_dormant  )) deallocate (a_B_methanotroph_dormant  )
      IF (allocated(a_f_T_methanogen          )) deallocate (a_f_T_methanogen          )
      IF (allocated(a_f_S_methanogen          )) deallocate (a_f_S_methanogen          )
      IF (allocated(a_f_O2_methanogen         )) deallocate (a_f_O2_methanogen         )
      IF (allocated(a_f_T_methanotroph        )) deallocate (a_f_T_methanotroph        )
      IF (allocated(a_methanogen_growth_rate  )) deallocate (a_methanogen_growth_rate  )
      IF (allocated(a_methanotroph_growth_rate)) deallocate (a_methanotroph_growth_rate)
      IF (allocated(a_microbial_prod_potential)) deallocate (a_microbial_prod_potential)
      IF (allocated(a_microbial_oxid_potential)) deallocate (a_microbial_oxid_potential)
      IF (allocated(a_methane_acc_num        )) deallocate (a_methane_acc_num        )
      IF (allocated(a_methane_acc_num_unsat  )) deallocate (a_methane_acc_num_unsat  )
      IF (allocated(a_methane_acc_num_sat    )) deallocate (a_methane_acc_num_sat    )
      IF (allocated(a_methane_acc_num_lake   )) deallocate (a_methane_acc_num_lake   )
      IF (allocated(a_methane_acc_num_extra  )) deallocate (a_methane_acc_num_extra  )
      IF (allocated(a_methane_acc_num_microbe)) deallocate (a_methane_acc_num_microbe)

   END SUBROUTINE deallocate_methane_acc_fluxes

   !-------------------------------------------------------------------
   ! Private acc helpers — mirror MOD_Vars_1DAccFluxes:acc1d/acc2d but kept
   ! local to break the USE cycle. spval guards skip uninitialised slots.
   !-------------------------------------------------------------------
   SUBROUTINE acc1d (var, s)
      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: s  (:)
      integer :: i
      DO i = lbound(var,1), ubound(var,1)
         IF (valid_acc_value(var(i))) THEN
            IF (.not. valid_acc_value(s(i))) s(i) = 0._r8
            s(i) = s(i) + var(i)
         ENDIF
      ENDDO
   END SUBROUTINE acc1d

   SUBROUTINE acc2d (var, s)
      real(r8), intent(in)    :: var(:,:)
      real(r8), intent(inout) :: s  (:,:)
      integer :: i1, i2
      DO i2 = lbound(var,2), ubound(var,2)
         DO i1 = lbound(var,1), ubound(var,1)
            IF (valid_acc_value(var(i1,i2))) THEN
               IF (.not. valid_acc_value(s(i1,i2))) s(i1,i2) = 0._r8
               s(i1,i2) = s(i1,i2) + var(i1,i2)
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE acc2d

   SUBROUTINE acc_count1d (var, n)
      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: n(:)
      integer :: i
      DO i = lbound(var,1), ubound(var,1)
         IF (valid_acc_value(var(i))) n(i) = n(i) + 1._r8
      ENDDO
   END SUBROUTINE acc_count1d

   SUBROUTINE acc_count1d_masked (var, n, mask)
      ! Like acc_count1d but the increment is also gated by an external mask
      ! (typically a patchtype/active flag).  Used for lake/wetland-only
      ! accumulators so non-active patches with a (legitimate but irrelevant)
      ! value of 0 don't inflate the denominator at history time.
      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: n(:)
      logical,  intent(in)    :: mask(:)
      integer :: k, ivar, inum, imask
      IF (size(var) /= size(n) .or. size(var) /= size(mask)) THEN
         write(6,*) 'acc_count1d_masked: size mismatch var/n/mask=', &
            size(var), size(n), size(mask)
         CALL abort
      ENDIF
      DO k = 1, size(var)
         ivar  = lbound(var,  1) + k - 1
         inum  = lbound(n,    1) + k - 1
         imask = lbound(mask, 1) + k - 1
         IF (mask(imask) .and. valid_acc_value(var(ivar))) n(inum) = n(inum) + 1._r8
      ENDDO
   END SUBROUTINE acc_count1d_masked

   SUBROUTINE acc_count2d (var, n)
      real(r8), intent(in)    :: var(:,:)
      real(r8), intent(inout) :: n(:)
      integer :: i1, i2
      logical :: any_valid
      DO i2 = lbound(var,2), ubound(var,2)
         any_valid = .false.
         DO i1 = lbound(var,1), ubound(var,1)
            any_valid = any_valid .or. valid_acc_value(var(i1,i2))
         ENDDO
         IF (any_valid) n(i2) = n(i2) + 1._r8
      ENDDO
   END SUBROUTINE acc_count2d

   ELEMENTAL LOGICAL FUNCTION valid_acc_value (x)
      real(r8), intent(in) :: x
      valid_acc_value = (.not. ieee_is_nan(x)) .and. abs(x) < 0.5_r8 * abs(spval)
   END FUNCTION valid_acc_value

END MODULE MOD_Tracer_Methane_AccFlux
#endif
