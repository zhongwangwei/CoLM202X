#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_BgcShim
!=======================================================================
! Reactive tracer wetland/BGC coupling shim.
!
! This module is the reactive-tracer boundary for invoking the BGC
! decomposition cascade needed by wetland CH4.  Reactive_Methane_Impl
! should orchestrate methane driver calls, not reach directly into
! MOD_BGC_* internals.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Vars_Global, only: nl_soil, z_soi, dz_soi, &
      ndecomp_pools, ndecomp_transitions
   USE MOD_BGC_Soil_BiogeochemDecompCascadeBGC, only: decomp_rate_constants_bgc
   USE MOD_BGC_Soil_BiogeochemPotential,        only: SoilBiogeochemPotential
   USE MOD_BGC_Soil_BiogeochemCompetition,      only: SoilBiogeochemCompetitionNoPlant
   USE MOD_BGC_Soil_BiogeochemDecomp,           only: SoilBiogeochemDecomp
   USE MOD_BGC_Vars_1DFluxes, only: decomp_cpools_sourcesink, decomp_npools_sourcesink, &
      decomp_hr_vr, decomp_ctransfer_vr, &
      decomp_ntransfer_vr, decomp_sminn_flux_vr, sminn_to_denit_decomp_vr, &
      pmnf_decomp, p_decomp_cpool_loss, net_nmin_vr, gross_nmin_vr, &
      net_nmin, gross_nmin, potential_immob_vr, phr_vr, pot_f_nit_vr, &
      decomp_hr, somc_fire, som_c_leached, som_n_leached, denit, f_n2o_nit, &
      smin_no3_leached, smin_no3_runoff, sminn_leached, sminn_to_plant
   USE MOD_BGC_Vars_TimeVariables, only: fpi_vr, o_scalar
   USE MOD_Vars_TimeVariables, only: t_soisno
   USE MOD_Vars_TimeInvariants, only: patchclass
   USE MOD_BGC_Vars_TimeInvariants, only: i_met_lit, i_cel_lit, i_lig_lit
   USE MOD_Vars_1DFluxes, only: assim
   USE MOD_Const_LC, only: rootfr
   USE MOD_Const_Physical, only: tfrz
   USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
   USE MOD_BGC_Soil_BiogeochemVerticalProfile, only: surfprof_exp
   USE MOD_Tracer_Reactive_Methane_WetlandVeg, only: wetveg_bg_frac

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: reactive_bgc_run_wetland_decomp


CONTAINS

   SUBROUTINE reactive_bgc_set_wetland_anoxia (ipatch)

! !DESCRIPTION:
!  Anoxia limitation on decomposition, which CoLM202X otherwise lacks.
!
!  The wetland hydrology holds every thawed layer at saturation, so those layers
!  are anaerobic and decompose at mino2lim of the aerobic rate -- the parameter
!  already carries exactly that definition. Without it a saturated tile has no
!  limiter left at all: t_scalar tracks temperature, depth_scalar is fixed, and
!  w_scalar is 1 precisely because the tile is waterlogged. Frozen layers keep 1
!  so the suppression is not counted twice against t_scalar and w_scalar --
!  unless frozen_anoxic_decomp (C-16) is on: just below freezing w_scalar is
!  still near 1, and 1 there made a freezing layer decompose five times faster
!  than a thawed one.
!
!  Set BEFORE decomposition runs -- decomp_rate_constants_bgc folds o_scalar into
!  decomp_k inside itself, and its patchtype 2 exemption is what preserves this.
!
!  CALLERS MUST RESTRICT THIS TO patchtype == 2. It lived inside the wetland
!  decomposition shim until 2026-08-02, and that shim is called for every patch,
!  so an upland soil was silently limited to 20% of its aerobic rate too. No
!  single-point run could show it -- all 44 towers are patchtype 2 -- but a
!  global run decomposes its entire land surface through this.

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      integer :: j

      IF (.not. allocated(o_scalar)) RETURN
      DO j = 1, nl_soil
         IF (t_soisno(j,ipatch) > tfrz .or. DEF_METHANE%frozen_anoxic_decomp) THEN
            o_scalar(j,ipatch) = max(DEF_METHANE%mino2lim, 1.e-6_r8)
         ELSE
            o_scalar(j,ipatch) = 1._r8
         ENDIF
      ENDDO

   END SUBROUTINE reactive_bgc_set_wetland_anoxia

   SUBROUTINE reactive_bgc_run_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      IF (.not. ieee_is_finite(deltim) .or. deltim <= 0._r8) THEN
         CALL CoLM_stop(' ***** ERROR: wetland CH4/BGC coupling requires a finite positive timestep')
      ENDIF

      ! Start from the same clean per-patch flux state as the full BGC driver.
      IF (allocated(decomp_cpools_sourcesink))   decomp_cpools_sourcesink  (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_npools_sourcesink))   decomp_npools_sourcesink  (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_hr_vr))              decomp_hr_vr             (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_ctransfer_vr))       decomp_ctransfer_vr      (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_ntransfer_vr))       decomp_ntransfer_vr      (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_sminn_flux_vr))      decomp_sminn_flux_vr     (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(sminn_to_denit_decomp_vr))  sminn_to_denit_decomp_vr (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(pmnf_decomp))               pmnf_decomp              (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(p_decomp_cpool_loss))       p_decomp_cpool_loss      (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(net_nmin_vr))               net_nmin_vr              (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(gross_nmin_vr))             gross_nmin_vr            (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(potential_immob_vr))        potential_immob_vr       (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(phr_vr))                    phr_vr                   (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(pot_f_nit_vr))              pot_f_nit_vr             (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(o_scalar))                  o_scalar                 (1:nl_soil,ipatch)   = 1._r8
      IF (allocated(fpi_vr))                    fpi_vr                   (1:nl_soil,ipatch)   = 1._r8
      IF (allocated(net_nmin))                  net_nmin                 (ipatch)             = 0._r8
      IF (allocated(gross_nmin))                gross_nmin               (ipatch)             = 0._r8
      IF (allocated(decomp_hr))                 decomp_hr                (ipatch)             = 0._r8
      IF (allocated(somc_fire))                 somc_fire                (ipatch)             = 0._r8
      IF (allocated(som_c_leached))             som_c_leached            (ipatch)             = 0._r8
      IF (allocated(som_n_leached))             som_n_leached            (ipatch)             = 0._r8
      IF (allocated(denit))                     denit                    (ipatch)             = 0._r8
      IF (allocated(f_n2o_nit))                 f_n2o_nit                (ipatch)             = 0._r8
      IF (allocated(smin_no3_leached))          smin_no3_leached         (ipatch)             = 0._r8
      IF (allocated(smin_no3_runoff))           smin_no3_runoff          (ipatch)             = 0._r8
      IF (allocated(sminn_leached))             sminn_leached            (ipatch)             = 0._r8
      IF (allocated(sminn_to_plant))            sminn_to_plant           (ipatch)             = 0._r8

      ! Plant carbon input of the tile (paper V2 C-12), on the source/sink
      ! that CDecompStateUpdate adds to the pools once they are not fixed.
      IF (DEF_METHANE%wetland_plant_input) CALL wetland_plant_litter_input (ipatch, deltim)

      ! Anoxia limit on the waterlogged tile, only when BGC owns it under the
      ! CH4/BGC contract (bgc_anoxia_limits_decomp); otherwise o_scalar stays 1.
      IF (DEF_METHANE%bgc_anoxia_limits_decomp) CALL reactive_bgc_set_wetland_anoxia (ipatch)

      CALL decomp_rate_constants_bgc (ipatch, nl_soil, z_soi)
      CALL SoilBiogeochemPotential   (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions)
      CALL SoilBiogeochemCompetitionNoPlant (ipatch, deltim, nl_soil, dz_soi)
      CALL SoilBiogeochemDecomp      (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions, dz_soi)

   END SUBROUTINE reactive_bgc_run_wetland_decomp

   SUBROUTINE wetland_plant_litter_input (ipatch, deltim)

      ! Canopy assimilation [mol CO2 m-2 s-1] times NPP/GPP gives the carbon
      ! the tile's plants return to the soil at steady state; it enters the
      ! litter pools along the land class's root profile with CoLM's grass
      ! litter split and a grass litter C:N (mean of leaf litter 50 and fine
      ! roots 42 in MOD_Const_PFT). With wetland_bg_frac < 1 (C-15) only that
      ! share follows the roots; the aboveground rest is laid on the surface
      ! along CoLM's leaf-litter profile, exp(-surfprof_exp z).
      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer  :: j
      real(r8) :: cin, s, bg, prof(1:nl_soil), sprof(1:nl_soil)

      IF (.not. allocated(decomp_cpools_sourcesink)) RETURN
      cin = assim(ipatch)
      IF (.not. ieee_is_finite(cin) .or. cin <= 0._r8 .or. cin > 1.e-2_r8) RETURN
      cin = cin * 12.011_r8 * DEF_METHANE%wetland_npp_frac * deltim       ! g C m-2 per step

      prof(:) = max(rootfr(1:nl_soil, patchclass(ipatch)), 0._r8)
      s = sum(prof)
      IF (s > 0._r8) THEN
         prof(:) = prof(:) / s
      ELSE
         prof(:) = 0._r8
         prof(1) = 1._r8
      ENDIF

      bg = wetveg_bg_frac(ipatch)      ! wetland_bg_frac, forest-weighted under C-13
      IF (bg < 1._r8) THEN
         sprof(:) = exp(-surfprof_exp * z_soi(1:nl_soil)) * dz_soi(1:nl_soil)
         sprof(:) = sprof(:) / sum(sprof)
         prof(:) = bg * prof(:) + (1._r8 - bg) * sprof(:)
      ENDIF

      DO j = 1, nl_soil
         decomp_cpools_sourcesink(j,i_met_lit,ipatch) = decomp_cpools_sourcesink(j,i_met_lit,ipatch) &
            + cin * prof(j) / dz_soi(j) * 0.25_r8
         decomp_cpools_sourcesink(j,i_cel_lit,ipatch) = decomp_cpools_sourcesink(j,i_cel_lit,ipatch) &
            + cin * prof(j) / dz_soi(j) * 0.50_r8
         decomp_cpools_sourcesink(j,i_lig_lit,ipatch) = decomp_cpools_sourcesink(j,i_lig_lit,ipatch) &
            + cin * prof(j) / dz_soi(j) * 0.25_r8
         IF (allocated(decomp_npools_sourcesink)) THEN
            decomp_npools_sourcesink(j,i_met_lit,ipatch) = decomp_npools_sourcesink(j,i_met_lit,ipatch) &
               + cin * prof(j) / dz_soi(j) * 0.25_r8 / DEF_METHANE%wetland_litter_cn
            decomp_npools_sourcesink(j,i_cel_lit,ipatch) = decomp_npools_sourcesink(j,i_cel_lit,ipatch) &
               + cin * prof(j) / dz_soi(j) * 0.50_r8 / DEF_METHANE%wetland_litter_cn
            decomp_npools_sourcesink(j,i_lig_lit,ipatch) = decomp_npools_sourcesink(j,i_lig_lit,ipatch) &
               + cin * prof(j) / dz_soi(j) * 0.25_r8 / DEF_METHANE%wetland_litter_cn
         ENDIF
      ENDDO

   END SUBROUTINE wetland_plant_litter_input

END MODULE MOD_Tracer_Reactive_BgcShim
#endif
