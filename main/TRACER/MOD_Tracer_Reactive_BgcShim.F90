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
   USE MOD_BGC_Vars_TimeVariables, only: fpi_vr, o_scalar, t_scalar, w_scalar, &
      depth_scalar
   USE MOD_Namelist, only: DEF_BGC_DEBUG_SCALARS
   USE MOD_Vars_TimeVariables, only: t_soisno, smp, wliq_soisno, wice_soisno
   USE MOD_Vars_TimeInvariants, only: porsl

   IMPLICIT NONE
   PRIVATE

   ! Call counter, used only to throttle the decomposition-scalar dump.
   integer :: dbg_calls = 0

   PUBLIC :: reactive_bgc_run_wetland_decomp

CONTAINS

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

      CALL decomp_rate_constants_bgc (ipatch, nl_soil, z_soi)
      CALL SoilBiogeochemPotential   (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions)
      CALL SoilBiogeochemCompetitionNoPlant (ipatch, deltim, nl_soil, dz_soi)
      CALL SoilBiogeochemDecomp      (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions, dz_soi)

      ! Decomposition-scalar dump.  decomp_k is the product of these five
      ! multipliers with a per-pool base rate, so printing the multipliers
      ! says which one differs between a site that emits and one that does
      ! not, without having to reason about the pool index.
      IF (DEF_BGC_DEBUG_SCALARS) CALL dump_decomp_scalars (ipatch, deltim)

   END SUBROUTINE reactive_bgc_run_wetland_decomp

   SUBROUTINE dump_decomp_scalars (ipatch, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer :: j, every

      ! roughly monthly, whatever the timestep
      every = max(1, nint(30._r8 * 86400._r8 / deltim))
      dbg_calls = dbg_calls + 1
      IF (mod(dbg_calls, every) /= 1) RETURN

      write(6,'(A,I8,A,I6)') 'BGCSCAL step=', dbg_calls, ' patch=', ipatch
      ! smp and the liquid/ice split are here because w_scalar turned out to be
      ! the multiplier that differs, and the floor it sits on is a test on smp:
      ! knowing w_scalar is 0.001 does not say whether the soil is genuinely at
      ! -200 m or whether smp itself is wrong.
      write(6,'(A)') 'BGCSCAL  lyr    t_soisno    t_scalar    w_scalar    o_scalar'&
                  // '  depth_scal      fpi_vr       hr_vr      smp_mm        wliq'&
                  // '        wice       porsl'
      DO j = 1, nl_soil
         write(6,'(A,I5,11E12.4)') 'BGCSCAL', j, &
            t_soisno(j,ipatch), t_scalar(j,ipatch), w_scalar(j,ipatch), &
            o_scalar(j,ipatch), depth_scalar(j,ipatch), fpi_vr(j,ipatch), &
            sum(decomp_hr_vr(j,:,ipatch)), smp(j,ipatch), &
            wliq_soisno(j,ipatch), wice_soisno(j,ipatch), porsl(j,ipatch)
      ENDDO

   END SUBROUTINE dump_decomp_scalars

END MODULE MOD_Tracer_Reactive_BgcShim
#endif
