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
   USE MOD_Namelist, only: DEF_USE_NITRIF, DEF_USE_SASU, DEF_USE_DiagMatrix
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Vars_TimeInvariants, only: BD_all
   USE MOD_Vars_Global, only: nl_soil, z_soi, dz_soi, spval, &
      ndecomp_pools, ndecomp_transitions
   USE MOD_BGC_Vars_TimeInvariants, only: donor_pool, receiver_pool
   USE MOD_BGC_Soil_BiogeochemDecompCascadeBGC, only: decomp_rate_constants_bgc
   USE MOD_BGC_Soil_BiogeochemPotential,        only: SoilBiogeochemPotential
   USE MOD_BGC_Soil_BiogeochemDecomp,           only: SoilBiogeochemDecomp
   USE MOD_BGC_CNSummary, only: CNDriverSummarizeNonvegetatedSoilStates
   USE MOD_BGC_Vars_1DFluxes, only: decomp_hr_vr, decomp_ctransfer_vr, &
      decomp_ntransfer_vr, decomp_sminn_flux_vr, sminn_to_denit_decomp_vr, &
      pmnf_decomp, p_decomp_cpool_loss, net_nmin_vr, gross_nmin_vr, &
      net_nmin, gross_nmin, potential_immob_vr, phr_vr, pot_f_nit_vr, &
      decomp_hr, somc_fire, som_c_leached, som_n_leached, denit, f_n2o_nit, &
      smin_no3_leached, smin_no3_runoff, sminn_leached, sminn_to_plant, &
      decomp_cpools_sourcesink, decomp_npools_sourcesink
   USE MOD_BGC_Vars_TimeVariables, only: fpi_vr, o_scalar, &
      decomp_cpools_vr, decomp_npools_vr, sminn_vr, smin_nh4_vr, smin_no3_vr, &
      totvegc, totvegn, ctrunc_veg, ntrunc_veg

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: reactive_bgc_run_wetland_decomp

CONTAINS

   SUBROUTINE reactive_bgc_run_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      CALL require_wetland_decomposition_state (ipatch, deltim)

      ! Start from the same clean per-patch flux state as the full BGC driver.
      decomp_cpools_sourcesink(:,:,ipatch) = 0._r8
      decomp_npools_sourcesink(:,:,ipatch) = 0._r8
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
      CALL SoilBiogeochemDecomp      (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions, dz_soi)
      CALL apply_wetland_decomposition_state (ipatch, deltim)

   END SUBROUTINE reactive_bgc_run_wetland_decomp

   SUBROUTINE require_wetland_decomposition_state (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      logical :: invalid_state
      real(r8), parameter :: aggregate_tol = 1.e-10_r8

      IF (deltim <= 0._r8) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition requires deltim > 0')
      ENDIF
      IF (DEF_USE_SASU .or. DEF_USE_DiagMatrix) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition-only debit does not support SASU/DiagMatrix')
      ENDIF
      IF (.not. allocated(decomp_cpools_vr) .or. .not. allocated(decomp_npools_vr) .or. &
          .not. allocated(decomp_cpools_sourcesink) .or. .not. allocated(decomp_npools_sourcesink) .or. &
          .not. allocated(sminn_vr) .or. .not. allocated(o_scalar) .or. .not. allocated(fpi_vr)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition state is not allocated')
      ENDIF
      IF (.not. allocated(BD_all) .or. size(BD_all,1) < nl_soil .or. &
          ipatch < 1 .or. ipatch > size(BD_all,2)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland soil bulk density is not allocated')
      ENDIF
      IF (.not. allocated(totvegc) .or. .not. allocated(totvegn) .or. &
          .not. allocated(ctrunc_veg) .or. .not. allocated(ntrunc_veg)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland vegetation aggregates are not allocated')
      ENDIF
      IF (.not. allocated(decomp_hr_vr) .or. .not. allocated(decomp_ctransfer_vr) .or. &
          .not. allocated(decomp_ntransfer_vr) .or. .not. allocated(decomp_sminn_flux_vr) .or. &
          .not. allocated(sminn_to_denit_decomp_vr) .or. .not. allocated(pmnf_decomp) .or. &
          .not. allocated(p_decomp_cpool_loss)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition flux state is not allocated')
      ENDIF
      IF (.not. allocated(net_nmin_vr) .or. .not. allocated(gross_nmin_vr) .or. &
          .not. allocated(potential_immob_vr) .or. .not. allocated(phr_vr) .or. &
          .not. allocated(net_nmin) .or. .not. allocated(gross_nmin) .or. &
          .not. allocated(decomp_hr)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition diagnostics are not allocated')
      ENDIF
      IF (ipatch < 1 .or. ipatch > size(decomp_cpools_vr,3) .or. &
          ipatch > size(decomp_npools_vr,3) .or. ipatch > size(sminn_vr,2)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland patch index is outside BGC state')
      ENDIF
      IF (any(BD_all(1:nl_soil,ipatch) /= BD_all(1:nl_soil,ipatch)) .or. &
          any(BD_all(1:nl_soil,ipatch) <= tiny(1._r8)) .or. &
          any(abs(BD_all(1:nl_soil,ipatch)) >= 0.5_r8*abs(spval))) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland soil bulk density is invalid')
      ENDIF
      IF (totvegc(ipatch) /= totvegc(ipatch) .or. totvegn(ipatch) /= totvegn(ipatch) .or. &
          ctrunc_veg(ipatch) /= ctrunc_veg(ipatch) .or. ntrunc_veg(ipatch) /= ntrunc_veg(ipatch) .or. &
          abs(totvegc(ipatch)) > aggregate_tol .or. abs(totvegn(ipatch)) > aggregate_tol .or. &
          abs(ctrunc_veg(ipatch)) > aggregate_tol .or. abs(ntrunc_veg(ipatch)) > aggregate_tol) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland patch unexpectedly contains vegetation C/N state')
      ENDIF
      IF (DEF_USE_NITRIF) THEN
         IF (.not. allocated(smin_nh4_vr) .or. .not. allocated(smin_no3_vr)) THEN
            CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: nitrification wetland N state is not allocated')
         ENDIF
      ENDIF

      invalid_state = any(decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch) /= &
                          decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch)) .or. &
                      any(decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch) /= &
                          decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch)) .or. &
                      any(abs(decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch)) >= 0.5_r8*abs(spval)) .or. &
                      any(abs(decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch)) >= 0.5_r8*abs(spval)) .or. &
                      any(decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch) < 0._r8) .or. &
                      any(decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch) < 0._r8) .or. &
                      any(sminn_vr(1:nl_soil,ipatch) /= sminn_vr(1:nl_soil,ipatch)) .or. &
                      any(abs(sminn_vr(1:nl_soil,ipatch)) >= 0.5_r8*abs(spval)) .or. &
                      any(sminn_vr(1:nl_soil,ipatch) < 0._r8)
      IF (DEF_USE_NITRIF) THEN
         invalid_state = invalid_state .or. &
            any(smin_nh4_vr(1:nl_soil,ipatch) /= smin_nh4_vr(1:nl_soil,ipatch)) .or. &
            any(smin_no3_vr(1:nl_soil,ipatch) /= smin_no3_vr(1:nl_soil,ipatch)) .or. &
            any(abs(smin_nh4_vr(1:nl_soil,ipatch)) >= 0.5_r8*abs(spval)) .or. &
            any(abs(smin_no3_vr(1:nl_soil,ipatch)) >= 0.5_r8*abs(spval)) .or. &
            any(smin_nh4_vr(1:nl_soil,ipatch) < 0._r8) .or. &
            any(smin_no3_vr(1:nl_soil,ipatch) < 0._r8)
      ENDIF
      IF (invalid_state) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition contains invalid C/N pool state')
      ENDIF

   END SUBROUTINE require_wetland_decomposition_state

   SUBROUTINE apply_wetland_decomposition_state (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim

      integer :: j, k, pool, donor, receiver
      real(r8) :: available, demand, supply, loss, scale
      real(r8) :: mineral_delta(1:nl_soil)
      real(r8) :: mineral_rate, take_nh4, remaining
      real(r8), parameter :: state_tol = 1.e-10_r8

      ! The full BGC driver builds these same source/sink identities in
      ! CStateUpdate1 and SoilBiogeochemNStateUpdate1, then applies them in
      ! SoilBiogeochemLittVertTransp.  Wetland patches have no PFT driver, so
      ! reproduce only the decomposition transaction here: no phenology,
      ! deposition, plant uptake, fire, or vertical mixing is repeated.

      ! Bound every transition by its beginning-of-step donor C.  Multiple
      ! pathways leaving one pool share one scale, preventing order-dependent
      ! overdraw while keeping respiration/transfer/N flux ratios together.
      DO j = 1, nl_soil
         DO pool = 1, ndecomp_pools
            loss = 0._r8
            DO k = 1, ndecomp_transitions
               IF (donor_pool(k) == pool) loss = loss + &
                  max(decomp_hr_vr(j,k,ipatch), 0._r8) + &
                  max(decomp_ctransfer_vr(j,k,ipatch), 0._r8)
            ENDDO
            IF (loss*deltim > decomp_cpools_vr(j,pool,ipatch) .and. loss > 0._r8) THEN
               scale = max(decomp_cpools_vr(j,pool,ipatch), 0._r8) / (loss*deltim)
               DO k = 1, ndecomp_transitions
                  IF (donor_pool(k) == pool) CALL scale_wetland_transition (j, k, ipatch, scale)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      ! Apply the same beginning-of-step bound to organic N donors.  Terminal
      ! transitions may return N to the mineral pool, so only positive donor
      ! losses participate in this limiter.
      DO j = 1, nl_soil
         DO pool = 1, ndecomp_pools
            loss = 0._r8
            DO k = 1, ndecomp_transitions
               IF (donor_pool(k) /= pool) CYCLE
               loss = loss + max(decomp_ntransfer_vr(j,k,ipatch) + &
                  merge(decomp_sminn_flux_vr(j,k,ipatch), 0._r8, receiver_pool(k) == 0), 0._r8)
            ENDDO
            IF (loss*deltim > decomp_npools_vr(j,pool,ipatch) .and. loss > 0._r8) THEN
               scale = max(decomp_npools_vr(j,pool,ipatch), 0._r8) / (loss*deltim)
               DO k = 1, ndecomp_transitions
                  IF (donor_pool(k) == pool) CALL scale_wetland_transition (j, k, ipatch, scale)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      ! Limit immobilization and explicit decomposition denitrification to
      ! mineral N available in this layer plus same-step mineralization.
      DO j = 1, nl_soil
         IF (DEF_USE_NITRIF) THEN
            available = max(smin_nh4_vr(j,ipatch), 0._r8) + max(smin_no3_vr(j,ipatch), 0._r8)
         ELSE
            available = max(sminn_vr(j,ipatch), 0._r8)
         ENDIF
         demand = 0._r8
         supply = 0._r8
         DO k = 1, ndecomp_transitions
            mineral_rate = wetland_mineral_n_rate (j, k, ipatch)
            demand = demand + max(-mineral_rate, 0._r8)
            supply = supply + max( mineral_rate, 0._r8)
         ENDDO
         IF (demand*deltim > available + supply*deltim .and. demand > 0._r8) THEN
            scale = min(1._r8, (available/deltim + supply) / demand)
            DO k = 1, ndecomp_transitions
               IF (wetland_mineral_n_rate(j,k,ipatch) < 0._r8) &
                  CALL scale_wetland_transition (j, k, ipatch, scale)
            ENDDO
         ENDIF
      ENDDO

      decomp_cpools_sourcesink(:,:,ipatch) = 0._r8
      decomp_npools_sourcesink(:,:,ipatch) = 0._r8
      mineral_delta(:) = 0._r8
      DO k = 1, ndecomp_transitions
         donor = donor_pool(k)
         receiver = receiver_pool(k)
         DO j = 1, nl_soil
            decomp_cpools_sourcesink(j,donor,ipatch) = &
               decomp_cpools_sourcesink(j,donor,ipatch) - &
               (decomp_hr_vr(j,k,ipatch) + decomp_ctransfer_vr(j,k,ipatch))*deltim
            decomp_npools_sourcesink(j,donor,ipatch) = &
               decomp_npools_sourcesink(j,donor,ipatch) - decomp_ntransfer_vr(j,k,ipatch)*deltim
            IF (receiver /= 0) THEN
               decomp_cpools_sourcesink(j,receiver,ipatch) = &
                  decomp_cpools_sourcesink(j,receiver,ipatch) + decomp_ctransfer_vr(j,k,ipatch)*deltim
               decomp_npools_sourcesink(j,receiver,ipatch) = &
                  decomp_npools_sourcesink(j,receiver,ipatch) + &
                  (decomp_ntransfer_vr(j,k,ipatch) + decomp_sminn_flux_vr(j,k,ipatch))*deltim
            ELSE
               decomp_npools_sourcesink(j,donor,ipatch) = &
                  decomp_npools_sourcesink(j,donor,ipatch) - decomp_sminn_flux_vr(j,k,ipatch)*deltim
            ENDIF
            mineral_delta(j) = mineral_delta(j) + wetland_mineral_n_rate(j,k,ipatch)*deltim
         ENDDO
      ENDDO

      ! Validate the entire prospective transaction before mutating any pool.
      IF (any(decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch) + &
              decomp_cpools_sourcesink(1:nl_soil,1:ndecomp_pools,ipatch) < -state_tol) .or. &
          any(decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch) + &
              decomp_npools_sourcesink(1:nl_soil,1:ndecomp_pools,ipatch) < -state_tol)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition would overdraw C/N donor state')
      ENDIF
      IF (DEF_USE_NITRIF) THEN
         IF (any(smin_nh4_vr(1:nl_soil,ipatch) + smin_no3_vr(1:nl_soil,ipatch) + &
                 mineral_delta(:) < -state_tol)) THEN
            CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition would overdraw mineral N')
         ENDIF
      ELSEIF (any(sminn_vr(1:nl_soil,ipatch) + mineral_delta(:) < -state_tol)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive_BgcShim: wetland decomposition would overdraw mineral N')
      ENDIF

      decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch) = max(0._r8, &
         decomp_cpools_vr(1:nl_soil,1:ndecomp_pools,ipatch) + &
         decomp_cpools_sourcesink(1:nl_soil,1:ndecomp_pools,ipatch))
      decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch) = max(0._r8, &
         decomp_npools_vr(1:nl_soil,1:ndecomp_pools,ipatch) + &
         decomp_npools_sourcesink(1:nl_soil,1:ndecomp_pools,ipatch))

      IF (DEF_USE_NITRIF) THEN
         DO j = 1, nl_soil
            IF (mineral_delta(j) >= 0._r8) THEN
               ! Full BGC places gross mineralization in NH4.  The wetland
               ! decomposition-only path has no plant/heterotroph competition
               ! allocation, so immobilization removes NH4 first and then NO3.
               ! This deterministic partition preserves total mineral N and
               ! non-negativity without inventing stale actual_immob fluxes.
               smin_nh4_vr(j,ipatch) = smin_nh4_vr(j,ipatch) + mineral_delta(j)
            ELSE
               remaining = -mineral_delta(j)
               take_nh4 = min(max(smin_nh4_vr(j,ipatch), 0._r8), remaining)
               smin_nh4_vr(j,ipatch) = max(0._r8, smin_nh4_vr(j,ipatch) - take_nh4)
               smin_no3_vr(j,ipatch) = max(0._r8, smin_no3_vr(j,ipatch) - (remaining-take_nh4))
            ENDIF
            sminn_vr(j,ipatch) = smin_nh4_vr(j,ipatch) + smin_no3_vr(j,ipatch)
         ENDDO
      ELSE
         sminn_vr(1:nl_soil,ipatch) = max(0._r8, sminn_vr(1:nl_soil,ipatch) + mineral_delta(:))
      ENDIF

      ! Source/sink work arrays have been consumed locally; leave no pending
      ! state for a later generic BGC transport call to apply a second time.
      decomp_cpools_sourcesink(:,:,ipatch) = 0._r8
      decomp_npools_sourcesink(:,:,ipatch) = 0._r8
      CALL CNDriverSummarizeNonvegetatedSoilStates (ipatch, nl_soil, dz_soi, ndecomp_pools)
      CALL refresh_wetland_decomposition_diagnostics (ipatch)

   END SUBROUTINE apply_wetland_decomposition_state

   REAL(r8) FUNCTION wetland_mineral_n_rate (j, k, ipatch)

      IMPLICIT NONE
      integer, intent(in) :: j, k, ipatch

      IF (receiver_pool(k) /= 0) THEN
         wetland_mineral_n_rate = -decomp_sminn_flux_vr(j,k,ipatch)
      ELSE
         wetland_mineral_n_rate = decomp_sminn_flux_vr(j,k,ipatch)
      ENDIF
      wetland_mineral_n_rate = wetland_mineral_n_rate - sminn_to_denit_decomp_vr(j,k,ipatch)

   END FUNCTION wetland_mineral_n_rate

   SUBROUTINE scale_wetland_transition (j, k, ipatch, scale)

      IMPLICIT NONE
      integer, intent(in) :: j, k, ipatch
      real(r8), intent(in) :: scale
      real(r8) :: bounded_scale

      bounded_scale = max(0._r8, min(1._r8, scale))
      p_decomp_cpool_loss(j,k,ipatch)      = p_decomp_cpool_loss(j,k,ipatch)      * bounded_scale
      pmnf_decomp(j,k,ipatch)              = pmnf_decomp(j,k,ipatch)              * bounded_scale
      decomp_hr_vr(j,k,ipatch)             = decomp_hr_vr(j,k,ipatch)             * bounded_scale
      decomp_ctransfer_vr(j,k,ipatch)      = decomp_ctransfer_vr(j,k,ipatch)      * bounded_scale
      decomp_ntransfer_vr(j,k,ipatch)      = decomp_ntransfer_vr(j,k,ipatch)      * bounded_scale
      decomp_sminn_flux_vr(j,k,ipatch)     = decomp_sminn_flux_vr(j,k,ipatch)     * bounded_scale
      sminn_to_denit_decomp_vr(j,k,ipatch) = sminn_to_denit_decomp_vr(j,k,ipatch) * bounded_scale

   END SUBROUTINE scale_wetland_transition

   SUBROUTINE refresh_wetland_decomposition_diagnostics (ipatch)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      integer :: j, k

      net_nmin(ipatch) = 0._r8
      gross_nmin(ipatch) = 0._r8
      decomp_hr(ipatch) = 0._r8
      DO j = 1, nl_soil
         net_nmin_vr(j,ipatch) = 0._r8
         gross_nmin_vr(j,ipatch) = 0._r8
         potential_immob_vr(j,ipatch) = 0._r8
         phr_vr(j,ipatch) = 0._r8
         DO k = 1, ndecomp_transitions
            net_nmin_vr(j,ipatch) = net_nmin_vr(j,ipatch) - pmnf_decomp(j,k,ipatch)
            gross_nmin_vr(j,ipatch) = gross_nmin_vr(j,ipatch) + max(-pmnf_decomp(j,k,ipatch), 0._r8)
            potential_immob_vr(j,ipatch) = potential_immob_vr(j,ipatch) + max(pmnf_decomp(j,k,ipatch), 0._r8)
            phr_vr(j,ipatch) = phr_vr(j,ipatch) + decomp_hr_vr(j,k,ipatch)
         ENDDO
         net_nmin(ipatch) = net_nmin(ipatch) + net_nmin_vr(j,ipatch)*dz_soi(j)
         gross_nmin(ipatch) = gross_nmin(ipatch) + gross_nmin_vr(j,ipatch)*dz_soi(j)
         decomp_hr(ipatch) = decomp_hr(ipatch) + phr_vr(j,ipatch)*dz_soi(j)
      ENDDO

   END SUBROUTINE refresh_wetland_decomposition_diagnostics

END MODULE MOD_Tracer_Reactive_BgcShim
#endif
