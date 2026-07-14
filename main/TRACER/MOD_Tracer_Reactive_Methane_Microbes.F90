#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_Microbes
!=======================================================================
! Optional microbial-pool diagnostics for the Methane tracer.
!
! This module is deliberately off by default.  When enabled with
! DEF_METHANE%use_microbial_pools it advances simple active/dormant
! methanogen and methanotroph biomass pools and exposes guarded potential
! CH4 production / oxidation diagnostics.  The legacy CLM/CTSM-style
! methane fluxes remain unchanged unless DEF_METHANE%use_microbial_flux_override
! is explicitly enabled in the physics driver.
!
! Coupling note: methane_driver advances microbial pools before methane()
! physics and passes the resulting potentials into the same land step.  The
! potentials are therefore an explicit-Euler tendency based on start-of-step
! CH4/O2/carbon states, not an implicit solve against post-transport
! concentrations.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_nan
   USE MOD_Vars_Global,        only: nl_soil, spval
   USE MOD_Const_Physical,     only: tfrz
   USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE, secspday, catomw

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Prognostic biomass pools [gC biomass / m3 soil]
   real(r8), allocatable :: B_methanogen(:,:)
   real(r8), allocatable :: B_methanotroph(:,:)
   real(r8), allocatable :: B_methanogen_dormant(:,:)
   real(r8), allocatable :: B_methanotroph_dormant(:,:)

   ! Per-step diagnostics [-], [day-1], [mol m-3 s-1]
   real(r8), allocatable :: f_T_methanogen(:,:)
   real(r8), allocatable :: f_S_methanogen(:,:)
   real(r8), allocatable :: f_O2_methanogen(:,:)
   real(r8), allocatable :: f_T_methanotroph(:,:)
   real(r8), allocatable :: methanogen_growth_rate(:,:)
   real(r8), allocatable :: methanotroph_growth_rate(:,:)
   real(r8), allocatable :: microbial_prod_potential(:,:)
   real(r8), allocatable :: microbial_oxid_potential(:,:)

   logical :: methane_microbes_lulcc_snapshot_valid = .false.
   real(r8), allocatable :: lulcc_B_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_B_methanotroph_old(:,:)
   real(r8), allocatable :: lulcc_B_methanogen_dormant_old(:,:)
   real(r8), allocatable :: lulcc_B_methanotroph_dormant_old(:,:)
   real(r8), allocatable :: lulcc_f_T_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_f_S_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_f_O2_methanogen_old(:,:)
   real(r8), allocatable :: lulcc_f_T_methanotroph_old(:,:)
   real(r8), allocatable :: lulcc_methanogen_growth_rate_old(:,:)
   real(r8), allocatable :: lulcc_methanotroph_growth_rate_old(:,:)
   real(r8), allocatable :: lulcc_microbial_prod_potential_old(:,:)
   real(r8), allocatable :: lulcc_microbial_oxid_potential_old(:,:)

   PUBLIC :: allocate_methane_microbes_state
   PUBLIC :: deallocate_methane_microbes_state
   PUBLIC :: methane_microbes_step
   PUBLIC :: read_methane_microbes_restart
   PUBLIC :: write_methane_microbes_restart
   PUBLIC :: save_methane_microbes_lulcc_state
   PUBLIC :: remap_methane_microbes_lulcc_state
   PUBLIC :: B_methanogen, B_methanotroph
   PUBLIC :: B_methanogen_dormant, B_methanotroph_dormant
   PUBLIC :: f_T_methanogen, f_S_methanogen, f_O2_methanogen, f_T_methanotroph
   PUBLIC :: methanogen_growth_rate, methanotroph_growth_rate
   PUBLIC :: microbial_prod_potential, microbial_oxid_potential

CONTAINS

	   SUBROUTINE allocate_methane_microbes_state(numpatch)
	      integer, intent(in) :: numpatch

	      IF (.not. DEF_METHANE%use_microbial_pools) RETURN
	      IF (allocated(B_methanogen)) RETURN
	      ! Keep zero-length arrays allocated on ranks with numpatch==0 when
	      ! microbial pools are enabled; restart/vector I/O paths are collective.

	      allocate(B_methanogen(nl_soil,numpatch))
      allocate(B_methanotroph(nl_soil,numpatch))
      allocate(B_methanogen_dormant(nl_soil,numpatch))
      allocate(B_methanotroph_dormant(nl_soil,numpatch))
      allocate(f_T_methanogen(nl_soil,numpatch))
      allocate(f_S_methanogen(nl_soil,numpatch))
      allocate(f_O2_methanogen(nl_soil,numpatch))
      allocate(f_T_methanotroph(nl_soil,numpatch))
      allocate(methanogen_growth_rate(nl_soil,numpatch))
      allocate(methanotroph_growth_rate(nl_soil,numpatch))
      allocate(microbial_prod_potential(nl_soil,numpatch))
      allocate(microbial_oxid_potential(nl_soil,numpatch))

      B_methanogen(:,:) = DEF_METHANE%B_init_methanogen
      B_methanotroph(:,:) = DEF_METHANE%B_init_methanotroph
      B_methanogen_dormant(:,:) = 0._r8
      B_methanotroph_dormant(:,:) = 0._r8
      f_T_methanogen(:,:) = spval
      f_S_methanogen(:,:) = spval
      f_O2_methanogen(:,:) = spval
      f_T_methanotroph(:,:) = spval
      methanogen_growth_rate(:,:) = spval
      methanotroph_growth_rate(:,:) = spval
      microbial_prod_potential(:,:) = spval
      microbial_oxid_potential(:,:) = spval
   END SUBROUTINE allocate_methane_microbes_state


   SUBROUTINE deallocate_methane_microbes_state()
      IF (allocated(B_methanogen)) deallocate(B_methanogen)
      IF (allocated(B_methanotroph)) deallocate(B_methanotroph)
      IF (allocated(B_methanogen_dormant)) deallocate(B_methanogen_dormant)
      IF (allocated(B_methanotroph_dormant)) deallocate(B_methanotroph_dormant)
      IF (allocated(f_T_methanogen)) deallocate(f_T_methanogen)
      IF (allocated(f_S_methanogen)) deallocate(f_S_methanogen)
      IF (allocated(f_O2_methanogen)) deallocate(f_O2_methanogen)
      IF (allocated(f_T_methanotroph)) deallocate(f_T_methanotroph)
      IF (allocated(methanogen_growth_rate)) deallocate(methanogen_growth_rate)
      IF (allocated(methanotroph_growth_rate)) deallocate(methanotroph_growth_rate)
      IF (allocated(microbial_prod_potential)) deallocate(microbial_prod_potential)
      IF (allocated(microbial_oxid_potential)) deallocate(microbial_oxid_potential)
      CALL clear_methane_microbes_lulcc_snapshot()
   END SUBROUTINE deallocate_methane_microbes_state


   SUBROUTINE methane_microbes_step(ipatch, deltim, t_soisno, conc_o2, conc_ch4, hr_vr, cellorg)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: t_soisno(1:nl_soil)
      real(r8), intent(in) :: conc_o2(1:nl_soil)
      real(r8), intent(in) :: conc_ch4(1:nl_soil)
      real(r8), intent(in) :: hr_vr(1:nl_soil)
      real(r8), intent(in) :: cellorg(1:nl_soil)

      integer :: j
      real(r8) :: dt_day, tempfac, sub_pool, sub_rate, f_s, f_o2, ch4_mm, o2_mm
	      real(r8) :: mu_m, mu_o, loss_m, loss_o, to_dormant, from_dormant
	      real(r8) :: prod_pot, oxid_pot, carbon_cap, ch4_cap, o2_cap, freeze_loss
	      real(r8) :: growth_factor, dormant_loss
	      real(r8) :: B_cap_methanogen, B_cap_methanotroph, organic_c_layer
      real(r8), parameter :: small = 1.e-30_r8

      IF (.not. DEF_METHANE%use_microbial_pools) RETURN
      IF (.not. allocated(B_methanogen)) RETURN
      IF (ipatch < lbound(B_methanogen,2) .or. ipatch > ubound(B_methanogen,2)) RETURN
      IF (deltim <= 0._r8) RETURN

      dt_day = deltim / secspday

      DO j = 1, nl_soil
	         IF (t_soisno(j) <= tfrz) THEN
	            freeze_loss = min(1._r8, max(0._r8, DEF_METHANE%gamma_microbial_freeze * dt_day))
	            B_methanogen(j,ipatch) = max(DEF_METHANE%B_min_methanogen, &
	               B_methanogen(j,ipatch) * (1._r8 - freeze_loss))
	            B_methanotroph(j,ipatch) = max(DEF_METHANE%B_min_methanotroph, &
	               B_methanotroph(j,ipatch) * (1._r8 - freeze_loss))
	            B_methanogen_dormant(j,ipatch) = max(0._r8, &
	               B_methanogen_dormant(j,ipatch) * (1._r8 - freeze_loss))
	            B_methanotroph_dormant(j,ipatch) = max(0._r8, &
	               B_methanotroph_dormant(j,ipatch) * (1._r8 - freeze_loss))

            f_T_methanogen(j,ipatch) = 0._r8
            f_S_methanogen(j,ipatch) = 0._r8
            f_O2_methanogen(j,ipatch) = 0._r8
            f_T_methanotroph(j,ipatch) = 0._r8
            methanogen_growth_rate(j,ipatch) = 0._r8
            methanotroph_growth_rate(j,ipatch) = 0._r8
            microbial_prod_potential(j,ipatch) = 0._r8
            microbial_oxid_potential(j,ipatch) = 0._r8
            CYCLE
         ELSE
            tempfac = DEF_METHANE%q10_microbe_growth ** &
               ((t_soisno(j) - DEF_METHANE%T_ref_microbe) / 10._r8)
         ENDIF

         ! cellorg is [kg OM m-3]; convert to [mol C m-3] using the
         ! same 580 gC kgOM-1 convention as the wetland/lake surface-data
         ! fallbacks and BgcLink.
         sub_pool = max(cellorg(j), 0._r8) * 580._r8 / catomw
         sub_rate = max(hr_vr(j) / catomw, 0._r8)
         f_s = sub_pool / (DEF_METHANE%K_substrate_methanogen_pool + sub_pool + small)
         f_o2 = DEF_METHANE%K_inh_O2_methanogen / &
            (DEF_METHANE%K_inh_O2_methanogen + max(conc_o2(j), 0._r8) + small)

	         mu_m = DEF_METHANE%mu_max_methanogen * tempfac * f_s * f_o2
	         loss_m = DEF_METHANE%gamma_methanogen * tempfac
	         to_dormant = 0._r8
	         from_dormant = 0._r8
	         growth_factor = exp(max(-50._r8, min(50._r8, (mu_m - loss_m) * dt_day)))
	         B_methanogen(j,ipatch) = max(DEF_METHANE%B_min_methanogen, &
            B_methanogen(j,ipatch) * growth_factor)
	         IF (DEF_METHANE%use_microbial_dormancy) THEN
	            IF (f_s < DEF_METHANE%dormancy_threshold_methanogen_fS .or. &
	                f_o2 < DEF_METHANE%dormancy_threshold_methanogen_fO2) THEN
	               to_dormant = min(B_methanogen(j,ipatch), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_active * B_methanogen(j,ipatch) * dt_day))
	            ELSE
	               from_dormant = min(B_methanogen_dormant(j,ipatch), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_revive * &
	                  B_methanogen_dormant(j,ipatch) * dt_day))
		            ENDIF
		         ENDIF

		         B_methanogen(j,ipatch) = max(DEF_METHANE%B_min_methanogen, &
            B_methanogen(j,ipatch) - to_dormant + from_dormant)
	         dormant_loss = min(B_methanogen_dormant(j,ipatch), &
	            max(0._r8, DEF_METHANE%gamma_microbial_dormant * tempfac * &
	            B_methanogen_dormant(j,ipatch) * dt_day))
	         B_methanogen_dormant(j,ipatch) = B_methanogen_dormant(j,ipatch) + &
	            to_dormant - from_dormant - dormant_loss
	         B_methanogen_dormant(j,ipatch) = max(B_methanogen_dormant(j,ipatch), 0._r8)

         ch4_mm = max(conc_ch4(j), 0._r8) / (DEF_METHANE%k_m + max(conc_ch4(j), 0._r8) + small)
         o2_mm = max(conc_o2(j), 0._r8) / (DEF_METHANE%k_m_o2 + max(conc_o2(j), 0._r8) + small)
         mu_o = DEF_METHANE%mu_max_methanotroph * tempfac * ch4_mm * o2_mm
         loss_o = DEF_METHANE%gamma_methanotroph * tempfac

	         growth_factor = exp(max(-50._r8, min(50._r8, (mu_o - loss_o) * dt_day)))
	         B_methanotroph(j,ipatch) = max(DEF_METHANE%B_min_methanotroph, &
            B_methanotroph(j,ipatch) * growth_factor)
	         to_dormant = 0._r8
	         from_dormant = 0._r8
	         IF (DEF_METHANE%use_microbial_dormancy) THEN
	            IF (ch4_mm < DEF_METHANE%dormancy_threshold_methanotroph_fS .or. &
	                o2_mm < DEF_METHANE%dormancy_threshold_methanotroph_fO2) THEN
	               to_dormant = min(B_methanotroph(j,ipatch), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_active * B_methanotroph(j,ipatch) * dt_day))
	            ELSE
	               from_dormant = min(B_methanotroph_dormant(j,ipatch), &
	                  max(0._r8, DEF_METHANE%dormancy_rate_revive * &
	                  B_methanotroph_dormant(j,ipatch) * dt_day))
	            ENDIF
	         ENDIF
	         B_methanotroph(j,ipatch) = max(DEF_METHANE%B_min_methanotroph, &
            B_methanotroph(j,ipatch) - to_dormant + from_dormant)
	         dormant_loss = min(B_methanotroph_dormant(j,ipatch), &
	            max(0._r8, DEF_METHANE%gamma_microbial_dormant * tempfac * &
	            B_methanotroph_dormant(j,ipatch) * dt_day))
	         B_methanotroph_dormant(j,ipatch) = B_methanotroph_dormant(j,ipatch) + &
	            to_dormant - from_dormant - dormant_loss
	         B_methanotroph(j,ipatch) = max(B_methanotroph(j,ipatch), DEF_METHANE%B_min_methanotroph)
	         B_methanotroph_dormant(j,ipatch) = max(B_methanotroph_dormant(j,ipatch), 0._r8)

         ! Constrain microbial biomass by local organic carbon.  Without this
         ! bound, positive net growth can compound exponentially while the
         ! represented substrate pool is not depleted by microbial growth.
         organic_c_layer = max(cellorg(j), 0._r8) * 580._r8
         IF (DEF_METHANE%B_max_fraction_methanogen > 0._r8) THEN
            B_cap_methanogen = max(DEF_METHANE%B_min_methanogen, &
               DEF_METHANE%B_max_fraction_methanogen * organic_c_layer)
            B_methanogen(j,ipatch) = min(B_methanogen(j,ipatch), B_cap_methanogen)
         ENDIF
         IF (DEF_METHANE%B_max_fraction_methanotroph > 0._r8) THEN
            B_cap_methanotroph = max(DEF_METHANE%B_min_methanotroph, &
               DEF_METHANE%B_max_fraction_methanotroph * organic_c_layer)
            B_methanotroph(j,ipatch) = min(B_methanotroph(j,ipatch), B_cap_methanotroph)
         ENDIF

	         ! B_* pools are stored as [gC biomass m-3 soil].  Treat kappa_m_*
	         ! as first-order biomass-C turnover [day-1] and convert biomass C
	         ! to mol C with catomw before exposing molar CH4 potentials.  This
	         ! keeps the optional microbial override in the same [mol m-3 s-1]
	         ! units as the legacy production/oxidation physics and the caps
	         ! below.
	         prod_pot = DEF_METHANE%kappa_m_methanogen * B_methanogen(j,ipatch) / catomw * &
	            tempfac * f_s * f_o2 / secspday
         carbon_cap = sub_rate
         prod_pot = min(max(prod_pot, 0._r8), carbon_cap)

	         oxid_pot = DEF_METHANE%kappa_m_methanotroph * B_methanotroph(j,ipatch) / catomw * &
	            tempfac * ch4_mm * o2_mm / secspday
         ch4_cap = max(conc_ch4(j), 0._r8) / deltim
         o2_cap = max(conc_o2(j), 0._r8) / (2._r8 * deltim)
         oxid_pot = min(max(oxid_pot, 0._r8), ch4_cap, o2_cap)

         f_T_methanogen(j,ipatch) = tempfac
         f_S_methanogen(j,ipatch) = f_s
         f_O2_methanogen(j,ipatch) = f_o2
         f_T_methanotroph(j,ipatch) = tempfac
         methanogen_growth_rate(j,ipatch) = mu_m
         methanotroph_growth_rate(j,ipatch) = mu_o
         microbial_prod_potential(j,ipatch) = prod_pot
         microbial_oxid_potential(j,ipatch) = oxid_pot
      END DO
   END SUBROUTINE methane_microbes_step


   SUBROUTINE write_methane_microbes_restart(file_restart, compress)
      USE MOD_LandPatch,    only: landpatch
      USE MOD_NetCDFVector, only: ncio_write_vector
      character(len=*), intent(in) :: file_restart
      integer,          intent(in) :: compress

      IF (.not. allocated(B_methanogen)) RETURN
      IF (.not. DEF_METHANE%use_microbial_pools) RETURN

      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanogen_dormant', 'soil', nl_soil, &
         'patch', landpatch, B_methanogen_dormant, compress)
      CALL ncio_write_vector(file_restart, 'ch4_B_methanotroph_dormant', 'soil', nl_soil, &
         'patch', landpatch, B_methanotroph_dormant, compress)
   END SUBROUTINE write_methane_microbes_restart


   SUBROUTINE read_methane_microbes_restart(file_restart)
      USE MOD_LandPatch,    only: landpatch
      USE MOD_NetCDFVector, only: ncio_read_vector => ncio_read_vector_complete
      character(len=*), intent(in) :: file_restart

      IF (.not. allocated(B_methanogen)) RETURN
      IF (.not. DEF_METHANE%use_microbial_pools) RETURN

      CALL ncio_read_vector(file_restart, 'ch4_B_methanogen', nl_soil, landpatch, &
         B_methanogen, defval=DEF_METHANE%B_init_methanogen)
      CALL ncio_read_vector(file_restart, 'ch4_B_methanotroph', nl_soil, landpatch, &
         B_methanotroph, defval=DEF_METHANE%B_init_methanotroph)
      CALL ncio_read_vector(file_restart, 'ch4_B_methanogen_dormant', nl_soil, landpatch, &
         B_methanogen_dormant, defval=0._r8)
      CALL ncio_read_vector(file_restart, 'ch4_B_methanotroph_dormant', nl_soil, landpatch, &
         B_methanotroph_dormant, defval=0._r8)

      ! Clean spval/NaN/negative biomass from restart.  defval only fires
      ! on completely-missing fields, but a malformed restart could still
      ! contain sentinels which would propagate into methane prod/oxid.
      WHERE (ieee_is_nan(B_methanogen) .or. abs(B_methanogen) > 1.e30_r8 .or. B_methanogen < 0._r8) &
         B_methanogen = DEF_METHANE%B_init_methanogen
      WHERE (ieee_is_nan(B_methanotroph) .or. abs(B_methanotroph) > 1.e30_r8 .or. B_methanotroph < 0._r8) &
         B_methanotroph = DEF_METHANE%B_init_methanotroph
      WHERE (ieee_is_nan(B_methanogen_dormant) .or. abs(B_methanogen_dormant) > 1.e30_r8 .or. B_methanogen_dormant < 0._r8) &
         B_methanogen_dormant = 0._r8
      WHERE (ieee_is_nan(B_methanotroph_dormant) .or. abs(B_methanotroph_dormant) > 1.e30_r8 .or. B_methanotroph_dormant < 0._r8) &
         B_methanotroph_dormant = 0._r8
   END SUBROUTINE read_methane_microbes_restart


   SUBROUTINE save_methane_microbes_lulcc_state()
      IF (.not. allocated(B_methanogen)) THEN
         methane_microbes_lulcc_snapshot_valid = .false.
         RETURN
      ENDIF

      CALL clear_methane_microbes_lulcc_snapshot()
      allocate(lulcc_B_methanogen_old(nl_soil,size(B_methanogen,2)))
      lulcc_B_methanogen_old = B_methanogen
      allocate(lulcc_B_methanotroph_old(nl_soil,size(B_methanotroph,2)))
      lulcc_B_methanotroph_old = B_methanotroph
      allocate(lulcc_B_methanogen_dormant_old(nl_soil,size(B_methanogen_dormant,2)))
      lulcc_B_methanogen_dormant_old = B_methanogen_dormant
      allocate(lulcc_B_methanotroph_dormant_old(nl_soil,size(B_methanotroph_dormant,2)))
      lulcc_B_methanotroph_dormant_old = B_methanotroph_dormant
      allocate(lulcc_f_T_methanogen_old(nl_soil,size(f_T_methanogen,2)))
      lulcc_f_T_methanogen_old = f_T_methanogen
      allocate(lulcc_f_S_methanogen_old(nl_soil,size(f_S_methanogen,2)))
      lulcc_f_S_methanogen_old = f_S_methanogen
      allocate(lulcc_f_O2_methanogen_old(nl_soil,size(f_O2_methanogen,2)))
      lulcc_f_O2_methanogen_old = f_O2_methanogen
      allocate(lulcc_f_T_methanotroph_old(nl_soil,size(f_T_methanotroph,2)))
      lulcc_f_T_methanotroph_old = f_T_methanotroph
      allocate(lulcc_methanogen_growth_rate_old(nl_soil,size(methanogen_growth_rate,2)))
      lulcc_methanogen_growth_rate_old = methanogen_growth_rate
      allocate(lulcc_methanotroph_growth_rate_old(nl_soil,size(methanotroph_growth_rate,2)))
      lulcc_methanotroph_growth_rate_old = methanotroph_growth_rate
      allocate(lulcc_microbial_prod_potential_old(nl_soil,size(microbial_prod_potential,2)))
      lulcc_microbial_prod_potential_old = microbial_prod_potential
      allocate(lulcc_microbial_oxid_potential_old(nl_soil,size(microbial_oxid_potential,2)))
      lulcc_microbial_oxid_potential_old = microbial_oxid_potential
      methane_microbes_lulcc_snapshot_valid = .true.
   END SUBROUTINE save_methane_microbes_lulcc_state


	   SUBROUTINE remap_methane_microbes_lulcc_state(patchclass_new, eindex_new, &
	      patchclass_old, eindex_old, lccpct_patches, new_patch_area, old_patch_area)
	      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
	      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
	      real(r8), intent(in), optional :: lccpct_patches(:,:)
	      real(r8), intent(in), optional :: new_patch_area(:)
	      real(r8), intent(in), optional :: old_patch_area(:)
      integer :: nnew

	      nnew = size(patchclass_new)
	      CALL deallocate_methane_microbes_arrays_only()
	      CALL allocate_methane_microbes_state(nnew)
	      IF (.not. allocated(B_methanogen)) RETURN
	      IF (.not. methane_microbes_lulcc_snapshot_valid) RETURN

      CALL remap2d(lulcc_B_methanogen_old, B_methanogen)
      CALL remap2d(lulcc_B_methanotroph_old, B_methanotroph)
      CALL remap2d(lulcc_B_methanogen_dormant_old, B_methanogen_dormant)
      CALL remap2d(lulcc_B_methanotroph_dormant_old, B_methanotroph_dormant)
      CALL remap2d(lulcc_f_T_methanogen_old, f_T_methanogen)
      CALL remap2d(lulcc_f_S_methanogen_old, f_S_methanogen)
      CALL remap2d(lulcc_f_O2_methanogen_old, f_O2_methanogen)
      CALL remap2d(lulcc_f_T_methanotroph_old, f_T_methanotroph)
      CALL remap2d(lulcc_methanogen_growth_rate_old, methanogen_growth_rate)
      CALL remap2d(lulcc_methanotroph_growth_rate_old, methanotroph_growth_rate)
      CALL remap2d(lulcc_microbial_prod_potential_old, microbial_prod_potential)
      CALL remap2d(lulcc_microbial_oxid_potential_old, microbial_oxid_potential)
      CALL clear_methane_microbes_lulcc_snapshot()

   CONTAINS
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
               new(:,np) = default_vals(:)
               src = fallback_source(np, size(old,2))
               IF (src > 0) THEN
                  new(1:min(size(new,1),size(old,1)),np) = &
                     old(1:min(size(new,1),size(old,1)),src)
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE remap2d

      integer FUNCTION fallback_source(np, nold) RESULT(src)
         integer, intent(in) :: np, nold
         integer :: op
         src = 0
         DO op = 1, min(nold, size(patchclass_old), size(eindex_old))
            IF (eindex_old(op) == eindex_new(np) .and. &
                patchclass_old(op) == patchclass_new(np)) THEN
               src = op
               RETURN
            ENDIF
         ENDDO
         ! Microbial pools are patch-class specific; avoid copying a same-eindex
         ! pool from a different land-cover class when no class match exists.
	      END FUNCTION fallback_source

	      REAL(r8) FUNCTION lulcc_source_weight(np, op) RESULT(w)
	         integer, intent(in) :: np, op
	         integer :: oq
	         real(r8) :: class_area

	         w = 0._r8
	         IF (.not. present(lccpct_patches)) RETURN
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
	   END SUBROUTINE remap_methane_microbes_lulcc_state


   SUBROUTINE deallocate_methane_microbes_arrays_only()
      IF (allocated(B_methanogen)) deallocate(B_methanogen)
      IF (allocated(B_methanotroph)) deallocate(B_methanotroph)
      IF (allocated(B_methanogen_dormant)) deallocate(B_methanogen_dormant)
      IF (allocated(B_methanotroph_dormant)) deallocate(B_methanotroph_dormant)
      IF (allocated(f_T_methanogen)) deallocate(f_T_methanogen)
      IF (allocated(f_S_methanogen)) deallocate(f_S_methanogen)
      IF (allocated(f_O2_methanogen)) deallocate(f_O2_methanogen)
      IF (allocated(f_T_methanotroph)) deallocate(f_T_methanotroph)
      IF (allocated(methanogen_growth_rate)) deallocate(methanogen_growth_rate)
      IF (allocated(methanotroph_growth_rate)) deallocate(methanotroph_growth_rate)
      IF (allocated(microbial_prod_potential)) deallocate(microbial_prod_potential)
      IF (allocated(microbial_oxid_potential)) deallocate(microbial_oxid_potential)
   END SUBROUTINE deallocate_methane_microbes_arrays_only


   SUBROUTINE clear_methane_microbes_lulcc_snapshot()
      IF (allocated(lulcc_B_methanogen_old)) deallocate(lulcc_B_methanogen_old)
      IF (allocated(lulcc_B_methanotroph_old)) deallocate(lulcc_B_methanotroph_old)
      IF (allocated(lulcc_B_methanogen_dormant_old)) deallocate(lulcc_B_methanogen_dormant_old)
      IF (allocated(lulcc_B_methanotroph_dormant_old)) deallocate(lulcc_B_methanotroph_dormant_old)
      IF (allocated(lulcc_f_T_methanogen_old)) deallocate(lulcc_f_T_methanogen_old)
      IF (allocated(lulcc_f_S_methanogen_old)) deallocate(lulcc_f_S_methanogen_old)
      IF (allocated(lulcc_f_O2_methanogen_old)) deallocate(lulcc_f_O2_methanogen_old)
      IF (allocated(lulcc_f_T_methanotroph_old)) deallocate(lulcc_f_T_methanotroph_old)
      IF (allocated(lulcc_methanogen_growth_rate_old)) deallocate(lulcc_methanogen_growth_rate_old)
      IF (allocated(lulcc_methanotroph_growth_rate_old)) deallocate(lulcc_methanotroph_growth_rate_old)
      IF (allocated(lulcc_microbial_prod_potential_old)) deallocate(lulcc_microbial_prod_potential_old)
      IF (allocated(lulcc_microbial_oxid_potential_old)) deallocate(lulcc_microbial_oxid_potential_old)
      methane_microbes_lulcc_snapshot_valid = .false.
   END SUBROUTINE clear_methane_microbes_lulcc_snapshot

END MODULE MOD_Tracer_Reactive_Methane_Microbes
#endif
