#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Methane_BgcLink
!=======================================================================
! BGC coupling bridge for the methane reactive tracer.
!
! This is the single place where Methane reaches into BGC internals.
! CoLM202X stores the required carbon fluxes in pool/PFT-level arrays, so
! this bridge reconstructs patch-level inputs for Methane without moving
! tracer-specific fields back into MOD_BGC_* modules.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_nan
   USE MOD_Vars_Global, only: nl_soil, dz_soi, spval
   USE MOD_Tracer_Methane_Const, only: DEF_METHANE
   USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
   USE MOD_Vars_PFTimeInvariants,  only: pftfrac
   USE MOD_BGC_Vars_TimeInvariants, only: organic_max, &
      i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3
   USE MOD_BGC_Vars_1DFluxes,       only: decomp_hr_vr, pot_f_nit_vr
   USE MOD_BGC_Vars_TimeVariables,  only: decomp_cpools_vr, o_scalar
   USE MOD_BGC_Vars_PFTimeVariables, only: annsum_npp_p, cinput_rootfr_p
   USE MOD_BGC_Vars_1DPFTFluxes,    only: froot_mr_p, &
      cpool_to_leafc_p, cpool_to_leafc_storage_p, &
      cpool_to_livestemc_p, cpool_to_livestemc_storage_p, &
      cpool_to_deadstemc_p, cpool_to_deadstemc_storage_p, &
      cpool_to_frootc_p, cpool_to_frootc_storage_p, &
      cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, &
      cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, &
      cpool_froot_gr_p, cpool_froot_storage_gr_p, &
      cpool_livecroot_gr_p, cpool_livecroot_storage_gr_p, &
      cpool_deadcroot_gr_p, cpool_deadcroot_storage_gr_p, &
      transfer_froot_gr_p, transfer_livecroot_gr_p, transfer_deadcroot_gr_p

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: tracer_ch4_bgc_patch_inputs
   PUBLIC :: organic_max

CONTAINS

   SUBROUTINE tracer_ch4_bgc_patch_inputs (ipatch, rootfr, crootfr, pH, cellorg, &
      somhr, lithr, hr_vr, rr, agnpp, bgnpp, annsum_npp, fphr, &
      o_scalar_out, pot_f_nit_vr_out, bgc_inputs_ready)
      integer,  intent(in)  :: ipatch
      real(r8), intent(in)  :: rootfr(1:nl_soil)
      real(r8), intent(out) :: crootfr(1:nl_soil)
      real(r8), intent(out) :: pH
      real(r8), intent(out) :: cellorg(1:nl_soil)
      real(r8), intent(out) :: somhr
      real(r8), intent(out) :: lithr
      real(r8), intent(out) :: hr_vr(1:nl_soil)
      real(r8), intent(out) :: rr
      real(r8), intent(out) :: agnpp
      real(r8), intent(out) :: bgnpp
      real(r8), intent(out) :: annsum_npp
      real(r8), intent(out) :: fphr(1:nl_soil)
      real(r8), intent(out) :: o_scalar_out(1:nl_soil)
      real(r8), intent(out) :: pot_f_nit_vr_out(1:nl_soil)
      logical,  intent(out), optional :: bgc_inputs_ready

      integer :: ps, pe, j, ipool, m
      integer :: decomp_pool_ids(7)
      real(r8) :: frac, profile_sum
      logical :: have_decomp_state, have_pft_state

      ps = -1
      pe = -1
      IF (allocated(patch_pft_s) .and. allocated(patch_pft_e)) THEN
         IF (ipatch >= 1 .and. ipatch <= size(patch_pft_s) .and. &
             ipatch <= size(patch_pft_e)) THEN
            ps = patch_pft_s(ipatch)
            pe = patch_pft_e(ipatch)
         ENDIF
      ENDIF

	      crootfr(:) = 0._r8
	      ! CoLM202X does not currently expose a spatial soil-pH field to this
	      ! bridge.  Only provide the Dunfield-fit optimum when pH scaling is
	      ! enabled; otherwise leave a sentinel so accidental downstream pH use
	      ! cannot apply a non-unit factor behind the switch.
	      if (DEF_METHANE%usephfact) then
	         pH = 6.2_r8
	      else
	         pH = spval
	      endif
      cellorg(:) = 0._r8
      somhr = 0._r8
      lithr = 0._r8
      hr_vr(:) = 0._r8
      rr = 0._r8
      agnpp = 0._r8
      bgnpp = 0._r8
      annsum_npp = 0._r8
      fphr(:) = 1.0_r8
      o_scalar_out(:) = 1.0_r8
      pot_f_nit_vr_out(:) = 0._r8

      decomp_pool_ids = (/ i_met_lit, i_cel_lit, i_lig_lit, i_cwd, &
         i_soil1, i_soil2, i_soil3 /)

      ! decomp_cpools_vr is [gC/m3].  Methane physics expects organic
      ! matter density [kg OM/m3], using 0.5 gC/gOM as in the source CH4
      ! conversion: gC/m3 * (1 gOM / 0.5 gC) * (1 kg / 1000 g).
      have_decomp_state = .false.
	      IF (allocated(decomp_cpools_vr) .and. allocated(decomp_hr_vr)) THEN
	         have_decomp_state = ipatch >= 1 .and. &
	                             ipatch <= size(decomp_cpools_vr,3) .and. &
	                             ipatch <= size(decomp_hr_vr,3)
	         IF (have_decomp_state) THEN
	            have_decomp_state = any(valid_bgc_value(decomp_cpools_vr(:,:,ipatch))) .and. &
	                                any(valid_bgc_value(decomp_hr_vr(:,:,ipatch)))
	         ENDIF
	      ENDIF
	      have_pft_state = pft_arrays_ready(ps, pe)
	      IF (present(bgc_inputs_ready)) bgc_inputs_ready = have_decomp_state .and. have_pft_state

	      IF (have_pft_state) THEN
	         DO m = ps, pe
	            frac = safe_nonnegative(pftfrac(m))
	            DO j = 1, nl_soil
	               IF (valid_bgc_value(cinput_rootfr_p(j,m))) THEN
	                  crootfr(j) = crootfr(j) + frac * max(0._r8, cinput_rootfr_p(j,m)) * dz_soi(j)
	               ENDIF
	            ENDDO
	         ENDDO
	      ENDIF
	      IF (sum(crootfr) <= 0._r8) THEN
	         DO j = 1, nl_soil
	            IF (valid_bgc_value(rootfr(j))) crootfr(j) = max(0._r8, rootfr(j))
	         ENDDO
	      ENDIF
	      profile_sum = sum(crootfr)
	      IF (profile_sum > 0._r8) THEN
	         crootfr(:) = crootfr(:) / profile_sum
	      ELSE
	         crootfr(1) = 1._r8
	      ENDIF

      IF (allocated(decomp_cpools_vr) .and. ipatch >= 1 .and. ipatch <= size(decomp_cpools_vr,3)) THEN
         DO j = 1, nl_soil
            DO m = 1, size(decomp_pool_ids)
               ipool = decomp_pool_ids(m)
               IF (ipool < 1 .or. ipool > size(decomp_cpools_vr,2)) CYCLE
               IF (valid_bgc_value(decomp_cpools_vr(j, ipool, ipatch))) THEN
                  cellorg(j) = cellorg(j) + max(0._r8, decomp_cpools_vr(j, ipool, ipatch))
               ENDIF
            END DO
         END DO
         cellorg(:) = cellorg(:) * 1.e-3_r8 / 0.5_r8
      ENDIF

      IF (allocated(decomp_hr_vr) .and. ipatch >= 1 .and. ipatch <= size(decomp_hr_vr,3)) THEN
         DO j = 1, nl_soil
            IF (valid_bgc_value(decomp_hr_vr(j, i_soil1, ipatch))) &
               somhr = somhr + max(0._r8, decomp_hr_vr(j, i_soil1, ipatch)) * dz_soi(j)
            IF (valid_bgc_value(decomp_hr_vr(j, i_soil2, ipatch))) &
               somhr = somhr + max(0._r8, decomp_hr_vr(j, i_soil2, ipatch)) * dz_soi(j)
            IF (valid_bgc_value(decomp_hr_vr(j, i_soil3, ipatch))) &
               somhr = somhr + max(0._r8, decomp_hr_vr(j, i_soil3, ipatch)) * dz_soi(j)

            IF (valid_bgc_value(decomp_hr_vr(j, i_met_lit, ipatch))) &
               lithr = lithr + max(0._r8, decomp_hr_vr(j, i_met_lit, ipatch)) * dz_soi(j)
            IF (valid_bgc_value(decomp_hr_vr(j, i_cel_lit, ipatch))) &
               lithr = lithr + max(0._r8, decomp_hr_vr(j, i_cel_lit, ipatch)) * dz_soi(j)
            IF (valid_bgc_value(decomp_hr_vr(j, i_lig_lit, ipatch))) &
               lithr = lithr + max(0._r8, decomp_hr_vr(j, i_lig_lit, ipatch)) * dz_soi(j)
            IF (valid_bgc_value(decomp_hr_vr(j, i_cwd,     ipatch))) &
               lithr = lithr + max(0._r8, decomp_hr_vr(j, i_cwd,     ipatch)) * dz_soi(j)

            DO ipool = 1, size(decomp_hr_vr,2)
               IF (valid_bgc_value(decomp_hr_vr(j, ipool, ipatch))) THEN
                  hr_vr(j) = hr_vr(j) + max(0._r8, decomp_hr_vr(j, ipool, ipatch))
               ENDIF
            END DO
         END DO
      ENDIF

      IF (allocated(o_scalar) .and. ipatch >= 1 .and. ipatch <= size(o_scalar,2)) THEN
         DO j = 1, nl_soil
            IF (valid_bgc_value(o_scalar(j, ipatch))) THEN
               o_scalar_out(j) = min(1._r8, max(0._r8, o_scalar(j, ipatch)))
            ENDIF
         END DO
      ENDIF

      IF (allocated(pot_f_nit_vr) .and. ipatch >= 1 .and. ipatch <= size(pot_f_nit_vr,2)) THEN
         DO j = 1, nl_soil
            IF (valid_bgc_value(pot_f_nit_vr(j, ipatch))) THEN
               pot_f_nit_vr_out(j) = max(0._r8, pot_f_nit_vr(j, ipatch))
            ENDIF
         END DO
      ENDIF

      IF (have_pft_state) THEN
         DO m = ps, pe
            frac = safe_nonnegative(pftfrac(m))

            agnpp = agnpp + frac * ( &
               safe_nonnegative(cpool_to_leafc_p(m)) &
             + safe_nonnegative(cpool_to_leafc_storage_p(m)) &
             + safe_nonnegative(cpool_to_livestemc_p(m)) &
             + safe_nonnegative(cpool_to_livestemc_storage_p(m)) &
             + safe_nonnegative(cpool_to_deadstemc_p(m)) &
             + safe_nonnegative(cpool_to_deadstemc_storage_p(m)) )

            bgnpp = bgnpp + frac * ( &
               safe_nonnegative(cpool_to_frootc_p(m)) &
             + safe_nonnegative(cpool_to_frootc_storage_p(m)) &
             + safe_nonnegative(cpool_to_livecrootc_p(m)) &
             + safe_nonnegative(cpool_to_livecrootc_storage_p(m)) &
             + safe_nonnegative(cpool_to_deadcrootc_p(m)) &
             + safe_nonnegative(cpool_to_deadcrootc_storage_p(m)) )

            rr = rr + frac * ( &
               safe_nonnegative(froot_mr_p(m)) &
             + safe_nonnegative(cpool_froot_gr_p(m)) &
             + safe_nonnegative(cpool_froot_storage_gr_p(m)) &
             + safe_nonnegative(cpool_livecroot_gr_p(m)) &
             + safe_nonnegative(cpool_livecroot_storage_gr_p(m)) &
             + safe_nonnegative(cpool_deadcroot_gr_p(m)) &
             + safe_nonnegative(cpool_deadcroot_storage_gr_p(m)) &
             + safe_nonnegative(transfer_froot_gr_p(m)) &
             + safe_nonnegative(transfer_livecroot_gr_p(m)) &
             + safe_nonnegative(transfer_deadcroot_gr_p(m)) )

            annsum_npp = annsum_npp + frac * safe_nonnegative(annsum_npp_p(m))
         END DO
      ENDIF

      ! Source CH4 treats fphr as a per-layer scaling factor. CoLM202X does
      ! not expose an equivalent field, so use full HR contribution while the
      ! rest of BGC coupling remains field-for-field reconstructed above.
   END SUBROUTINE tracer_ch4_bgc_patch_inputs


   REAL(r8) FUNCTION safe_nonnegative (x)
      real(r8), intent(in) :: x

      IF (valid_bgc_value(x)) THEN
         safe_nonnegative = max(0._r8, x)
      ELSE
         safe_nonnegative = 0._r8
      ENDIF
   END FUNCTION safe_nonnegative


   LOGICAL FUNCTION pft_arrays_ready (ps, pe)
      integer, intent(in) :: ps, pe

      pft_arrays_ready = .false.
      IF (ps <= 0 .or. pe < ps) RETURN
	      IF (.not. allocated(pftfrac)) RETURN
	      IF (.not. allocated(annsum_npp_p)) RETURN
	      IF (.not. allocated(cinput_rootfr_p)) RETURN
      IF (.not. allocated(froot_mr_p)) RETURN
      IF (.not. allocated(cpool_to_leafc_p)) RETURN
      IF (.not. allocated(cpool_to_leafc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_livestemc_p)) RETURN
      IF (.not. allocated(cpool_to_livestemc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_deadstemc_p)) RETURN
      IF (.not. allocated(cpool_to_deadstemc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_frootc_p)) RETURN
      IF (.not. allocated(cpool_to_frootc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_livecrootc_p)) RETURN
      IF (.not. allocated(cpool_to_livecrootc_storage_p)) RETURN
      IF (.not. allocated(cpool_to_deadcrootc_p)) RETURN
      IF (.not. allocated(cpool_to_deadcrootc_storage_p)) RETURN
      IF (.not. allocated(cpool_froot_gr_p)) RETURN
      IF (.not. allocated(cpool_froot_storage_gr_p)) RETURN
      IF (.not. allocated(cpool_livecroot_gr_p)) RETURN
      IF (.not. allocated(cpool_livecroot_storage_gr_p)) RETURN
      IF (.not. allocated(cpool_deadcroot_gr_p)) RETURN
      IF (.not. allocated(cpool_deadcroot_storage_gr_p)) RETURN
      IF (.not. allocated(transfer_froot_gr_p)) RETURN
      IF (.not. allocated(transfer_livecroot_gr_p)) RETURN
      IF (.not. allocated(transfer_deadcroot_gr_p)) RETURN

	      IF (pe > size(pftfrac)) RETURN
	      IF (pe > size(annsum_npp_p)) RETURN
	      IF (nl_soil > size(cinput_rootfr_p,1)) RETURN
	      IF (pe > size(cinput_rootfr_p,2)) RETURN
      IF (pe > size(froot_mr_p)) RETURN
      IF (pe > size(cpool_to_leafc_p)) RETURN
      IF (pe > size(cpool_to_leafc_storage_p)) RETURN
      IF (pe > size(cpool_to_livestemc_p)) RETURN
      IF (pe > size(cpool_to_livestemc_storage_p)) RETURN
      IF (pe > size(cpool_to_deadstemc_p)) RETURN
      IF (pe > size(cpool_to_deadstemc_storage_p)) RETURN
      IF (pe > size(cpool_to_frootc_p)) RETURN
      IF (pe > size(cpool_to_frootc_storage_p)) RETURN
      IF (pe > size(cpool_to_livecrootc_p)) RETURN
      IF (pe > size(cpool_to_livecrootc_storage_p)) RETURN
      IF (pe > size(cpool_to_deadcrootc_p)) RETURN
      IF (pe > size(cpool_to_deadcrootc_storage_p)) RETURN
      IF (pe > size(cpool_froot_gr_p)) RETURN
      IF (pe > size(cpool_froot_storage_gr_p)) RETURN
      IF (pe > size(cpool_livecroot_gr_p)) RETURN
      IF (pe > size(cpool_livecroot_storage_gr_p)) RETURN
      IF (pe > size(cpool_deadcroot_gr_p)) RETURN
      IF (pe > size(cpool_deadcroot_storage_gr_p)) RETURN
      IF (pe > size(transfer_froot_gr_p)) RETURN
      IF (pe > size(transfer_livecroot_gr_p)) RETURN
      IF (pe > size(transfer_deadcroot_gr_p)) RETURN

      pft_arrays_ready = .true.
   END FUNCTION pft_arrays_ready


	   ELEMENTAL LOGICAL FUNCTION valid_bgc_value (x)
	      real(r8), intent(in) :: x

	      valid_bgc_value = (.not. ieee_is_nan(x)) .and. &
         abs(x) < 0.5_r8 * abs(spval)
	   END FUNCTION valid_bgc_value

END MODULE MOD_Tracer_Methane_BgcLink
#endif
