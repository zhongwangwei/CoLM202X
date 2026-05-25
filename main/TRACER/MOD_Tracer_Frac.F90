#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Frac

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_TRACER_USE_FRACTIONATION, &
      DEF_TRACER_NSS_LEAF_WATER_PER_LAI, DEF_TRACER_NSS_LEAF_PATH_LENGTH, &
      DEF_TRACER_NSS_LEAF_RB
   USE MOD_Tracer_Defs, only: tracers, ntracers, tracer_is_isotope, &
      Rsmow_18O, Rsmow_D, trc_tiny

   IMPLICIT NONE
   SAVE

   integer, parameter :: iso_none = 0
   integer, parameter :: iso_o18  = 1
   integer, parameter :: iso_d    = 2

   real(r8), parameter :: tfrz = 273.15_r8
   real(r8), parameter :: eps_h2o_over_h218o_air = 1.03189_r8
   real(r8), parameter :: eps_h2o_over_hdo_air   = 1.01636_r8
   real(r8), parameter :: water_moles_per_mm = 1000._r8 / 18.01528_r8
   real(r8), parameter :: liquid_water_molar_density = 55.5e3_r8
   real(r8), parameter :: universal_gas_constant = 8.31446261815324_r8
   real(r8), parameter :: craig_gordon_relhum_max = 0.95_r8
   real(r8), parameter :: craig_gordon_relhum_fallback = 0.90_r8
   real(r8), parameter :: craig_gordon_min_net_ratio_frac = 0.75_r8

   PUBLIC :: tracer_fractionation_active
   PUBLIC :: tracer_alpha_liq_vap, tracer_alpha_ice_vap, tracer_alpha_ice_liq
   PUBLIC :: tracer_rayleigh_freezing_loss
   PUBLIC :: tracer_diffusivity_ratio_air
   PUBLIC :: tracer_alpha_kinetic_leaf, tracer_alpha_kinetic_soil
   PUBLIC :: tracer_craig_gordon_evap_ratio
   PUBLIC :: tracer_equilibrium_deposition_ratio
   PUBLIC :: tracer_transpiration_nss_ratio
   PUBLIC :: tracer_saturation_vapor_pressure, tracer_surface_relhum

CONTAINS

   logical FUNCTION tracer_fractionation_active (itrc)
      integer, intent(in) :: itrc

      tracer_fractionation_active = .false.
      IF (.not. DEF_TRACER_USE_FRACTIONATION) RETURN
      IF (.not. tracer_is_isotope(itrc)) RETURN
      tracer_fractionation_active = tracer_isotope_kind(itrc) /= iso_none
   END FUNCTION tracer_fractionation_active

   real(r8) FUNCTION tracer_alpha_liq_vap (itrc, temp_k)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tracer_alpha_liq_vap = 1._r8
      tk = max(temp_k, 150._r8)

      SELECT CASE (tracer_isotope_kind(itrc))
      CASE (iso_o18)
         ! Majoube liquid-vapor equilibrium, alpha = R_liquid / R_vapor.
         tracer_alpha_liq_vap = exp(1137._r8 / (tk * tk) - 0.4156_r8 / tk - 0.0020667_r8)
      CASE (iso_d)
         tracer_alpha_liq_vap = exp(24844._r8 / (tk * tk) - 76.248_r8 / tk + 0.052612_r8)
      END SELECT
   END FUNCTION tracer_alpha_liq_vap

   real(r8) FUNCTION tracer_alpha_ice_vap (itrc, temp_k)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tracer_alpha_ice_vap = 1._r8
      tk = max(temp_k, 150._r8)

      SELECT CASE (tracer_isotope_kind(itrc))
      CASE (iso_o18)
         ! Ice-vapor form used by IsoLESC microphysics; alpha = R_ice / R_vapor.
         tracer_alpha_ice_vap = exp(11.839_r8 / tk - 0.028224_r8)
      CASE (iso_d)
         tracer_alpha_ice_vap = exp(16289._r8 / (tk * tk) - 0.0945_r8)
      END SELECT
   END FUNCTION tracer_alpha_ice_vap

   real(r8) FUNCTION tracer_alpha_ice_liq (itrc, temp_k)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: temp_k

      ! alpha = R_ice / R_liquid, derived from the existing equilibrium
      ! vapor references so phase-change and vapor processes stay consistent.
      tracer_alpha_ice_liq = tracer_alpha_ice_vap(itrc, temp_k) / &
         max(tracer_alpha_liq_vap(itrc, temp_k), trc_tiny)
   END FUNCTION tracer_alpha_ice_liq

   real(r8) FUNCTION tracer_rayleigh_freezing_loss (itrc, pool_trc, pool_water, &
      freeze_water, temp_k)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: pool_trc
      real(r8), intent(in) :: pool_water
      real(r8), intent(in) :: freeze_water
      real(r8), intent(in) :: temp_k
      real(r8) :: alpha_il
      real(r8) :: liquid_fraction, remaining_water
      real(r8) :: source_ratio, remaining_ratio, remaining_trc

      tracer_rayleigh_freezing_loss = 0._r8
      IF (pool_water <= trc_tiny .or. freeze_water <= trc_tiny) RETURN
      IF (pool_trc <= trc_tiny) RETURN

      IF (freeze_water >= pool_water * (1._r8 - 1.e-12_r8)) THEN
         tracer_rayleigh_freezing_loss = max(pool_trc, 0._r8)
         RETURN
      ENDIF

      source_ratio = max(pool_trc, 0._r8) / pool_water
      IF (.not. tracer_fractionation_active(itrc)) THEN
         tracer_rayleigh_freezing_loss = min(freeze_water * source_ratio, max(pool_trc, 0._r8))
         RETURN
      ENDIF

      alpha_il = tracer_alpha_ice_liq(itrc, temp_k)
      IF (alpha_il <= trc_tiny .or. alpha_il /= alpha_il) THEN
         tracer_rayleigh_freezing_loss = min(freeze_water * source_ratio, max(pool_trc, 0._r8))
         RETURN
      ENDIF

      ! Equilibrium freezing of a finite mixed liquid pool:
      ! R_liq = R0 * f**(alpha_ice_liq - 1), with f = remaining liquid fraction.
      ! The tracer transferred to ice is the mass removed from the liquid pool.
      remaining_water = pool_water - freeze_water
      liquid_fraction = max(remaining_water / pool_water, trc_tiny)
      remaining_ratio = source_ratio * liquid_fraction ** (alpha_il - 1._r8)
      remaining_trc = remaining_water * remaining_ratio
      tracer_rayleigh_freezing_loss = min(max(pool_trc - remaining_trc, 0._r8), &
         max(pool_trc, 0._r8))
   END FUNCTION tracer_rayleigh_freezing_loss

   real(r8) FUNCTION tracer_diffusivity_ratio_air (itrc)
      integer, intent(in) :: itrc

      ! Return D(H2O) / D(isotopologue) in air. The O18/D values match the
      ! IsoLESC isotope_ET constants and are kept here instead of being
      ! duplicated in individual process modules.
      tracer_diffusivity_ratio_air = 1._r8
      SELECT CASE (tracer_isotope_kind(itrc))
      CASE (iso_o18)
         tracer_diffusivity_ratio_air = eps_h2o_over_h218o_air
      CASE (iso_d)
         tracer_diffusivity_ratio_air = eps_h2o_over_hdo_air
      END SELECT
   END FUNCTION tracer_diffusivity_ratio_air

   real(r8) FUNCTION tracer_alpha_kinetic_leaf (itrc, ra, rb, rc)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: ra, rb, rc
      real(r8) :: ra1, rb1, rc1, denom, diff_ratio

      ra1 = max(ra, 0._r8)
      rb1 = max(rb, 0._r8)
      rc1 = max(rc, 0._r8)
      denom = ra1 + rb1 + rc1
      tracer_alpha_kinetic_leaf = 1._r8
      IF (denom <= trc_tiny) RETURN

      diff_ratio = tracer_diffusivity_ratio_air(itrc)
      tracer_alpha_kinetic_leaf = (ra1 + rb1 * diff_ratio ** (2._r8 / 3._r8) + &
         rc1 * diff_ratio) / denom
   END FUNCTION tracer_alpha_kinetic_leaf

   real(r8) FUNCTION tracer_alpha_kinetic_soil (itrc, ra, rs)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: ra, rs
      real(r8) :: ra1, rs1, denom, diff_ratio

      ra1 = max(ra, 0._r8)
      rs1 = max(rs, 0._r8)
      denom = ra1 + rs1
      tracer_alpha_kinetic_soil = 1._r8
      IF (denom <= trc_tiny) RETURN

      diff_ratio = tracer_diffusivity_ratio_air(itrc)
      tracer_alpha_kinetic_soil = (ra1 + diff_ratio * rs1) / denom
   END FUNCTION tracer_alpha_kinetic_soil

   real(r8) FUNCTION tracer_craig_gordon_evap_ratio (itrc, source_ratio, vapor_ratio, &
      temp_k, relhum, alpha_k, from_ice)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: source_ratio
      real(r8), intent(in) :: vapor_ratio
      real(r8), intent(in) :: temp_k
      real(r8), intent(in) :: relhum
      real(r8), intent(in) :: alpha_k
      logical,  intent(in) :: from_ice

      real(r8) :: h, h_raw, one_minus_h, alpha_eq, ak
      real(r8) :: equilibrium_ratio, cg_numerator
      real(r8) :: cg_ratio

      tracer_craig_gordon_evap_ratio = source_ratio
      IF (.not. tracer_fractionation_active(itrc)) RETURN

      IF (from_ice) THEN
         alpha_eq = tracer_alpha_ice_vap(itrc, temp_k)
      ELSE
         alpha_eq = tracer_alpha_liq_vap(itrc, temp_k)
      ENDIF
      ak = max(alpha_k, trc_tiny)
      IF (alpha_eq <= trc_tiny .or. alpha_eq /= alpha_eq) RETURN

      equilibrium_ratio = source_ratio / alpha_eq
      tracer_craig_gordon_evap_ratio = equilibrium_ratio

      h_raw = min(max(relhum, 0._r8), 0.999999_r8)
      ! The one-way net-loss Craig-Gordon form becomes singular near
      ! saturation and can produce an almost tracer-free vapor when
      ! h*R_vapor balances R_source/alpha. In those cases the net flux
      ! isotope signature is poorly constrained, so use the equilibrium
      ! vapor in contact with the source pool instead.
      IF (h_raw >= craig_gordon_relhum_fallback) RETURN

      h = min(h_raw, craig_gordon_relhum_max)
      one_minus_h = 1._r8 - h
      IF (one_minus_h <= 1.e-8_r8) RETURN

      cg_numerator = equilibrium_ratio - h * vapor_ratio
      IF (cg_numerator <= trc_tiny) RETURN

      cg_ratio = cg_numerator / (ak * one_minus_h)
      IF (cg_ratio <= craig_gordon_min_net_ratio_frac * equilibrium_ratio) RETURN

      tracer_craig_gordon_evap_ratio = cg_ratio
      IF (tracer_craig_gordon_evap_ratio /= tracer_craig_gordon_evap_ratio) THEN
         tracer_craig_gordon_evap_ratio = source_ratio
      ENDIF
      tracer_craig_gordon_evap_ratio = max(tracer_craig_gordon_evap_ratio, 0._r8)
   END FUNCTION tracer_craig_gordon_evap_ratio

   real(r8) FUNCTION tracer_equilibrium_deposition_ratio (itrc, vapor_ratio, temp_k, from_ice)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: vapor_ratio
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: from_ice

      tracer_equilibrium_deposition_ratio = vapor_ratio
      IF (.not. tracer_fractionation_active(itrc)) RETURN

      IF (from_ice) THEN
         tracer_equilibrium_deposition_ratio = tracer_alpha_ice_vap(itrc, temp_k) * vapor_ratio
      ELSE
         tracer_equilibrium_deposition_ratio = tracer_alpha_liq_vap(itrc, temp_k) * vapor_ratio
      ENDIF
   END FUNCTION tracer_equilibrium_deposition_ratio

   SUBROUTINE tracer_transpiration_nss_ratio (itrc, source_ratio, vapor_ratio, &
      temp_k, relhum, psrf, transp_water, deltim, leaf_area, stomatal_resistance, &
      prev_delta_e, prev_delta_b, prev_peclet, prev_leaf_moles, &
      trans_ratio, new_delta_e, new_delta_b, new_peclet, new_leaf_moles)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: source_ratio
      real(r8), intent(in) :: vapor_ratio
      real(r8), intent(in) :: temp_k
      real(r8), intent(in) :: relhum
      real(r8), intent(in) :: psrf
      real(r8), intent(in) :: transp_water
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: leaf_area
      real(r8), intent(in) :: stomatal_resistance
      real(r8), intent(in) :: prev_delta_e
      real(r8), intent(in) :: prev_delta_b
      real(r8), intent(in) :: prev_peclet
      real(r8), intent(in) :: prev_leaf_moles
      real(r8), intent(out) :: trans_ratio
      real(r8), intent(out) :: new_delta_e
      real(r8), intent(out) :: new_delta_b
      real(r8), intent(out) :: new_peclet
      real(r8), intent(out) :: new_leaf_moles

      real(r8) :: h, one_minus_h, alpha_eq, alpha_k, eps_eq, eps_k
      real(r8) :: delta_x, delta_v, delta_es, delta_t
      real(r8) :: leaf_moles, transp_moles, transp_moles_s
      real(r8) :: prev_e, prev_p, prev_w, liquid_diff, peclet_number
      real(r8) :: leaf_moles_ratio
      real(r8) :: relax_b, denom, tk
      real(r8) :: gross_moles, conductance_gross_moles
      real(r8) :: vapor_molar_density_sat, total_resistance
	      real(r8) :: prev_leaf_water, new_leaf_water, prev_bulk_ratio, new_bulk_ratio
	      real(r8) :: storage_tracer_change, storage_tracer_change_used
	      real(r8) :: storage_scale, storage_bound, target_leaf_storage
	      real(r8), parameter :: max_storage_jump_for_fallback = 10._r8
	      real(r8), parameter :: max_storage_tendency_fraction = 0.95_r8

      trans_ratio = source_ratio
      new_delta_e = tracer_ratio_to_delta(itrc, source_ratio)
      new_delta_b = new_delta_e
      new_peclet = 1._r8
      new_leaf_moles = 0._r8

      IF (.not. tracer_fractionation_active(itrc)) RETURN
      IF (source_ratio <= trc_tiny .or. vapor_ratio <= trc_tiny) RETURN
      IF (transp_water <= trc_tiny .or. deltim <= trc_tiny) RETURN

      leaf_moles = max(leaf_area, 0._r8) * max(DEF_TRACER_NSS_LEAF_WATER_PER_LAI, 0._r8) &
         * water_moles_per_mm
      new_leaf_moles = leaf_moles
      IF (leaf_moles <= trc_tiny) RETURN

      tk = max(temp_k, 150._r8)
      h = min(max(relhum, 0._r8), 0.95_r8)
      one_minus_h = max(1._r8 - h, 1.e-6_r8)
      alpha_eq = tracer_alpha_liq_vap(itrc, tk)
      eps_eq = (1._r8 - 1._r8 / max(alpha_eq, trc_tiny)) * 1000._r8
      eps_k = tracer_leaf_kinetic_epsilon(itrc, 0._r8, DEF_TRACER_NSS_LEAF_RB, &
         stomatal_resistance)
      alpha_k = 1._r8 + eps_k / 1000._r8

      delta_x = tracer_ratio_to_delta(itrc, source_ratio)
      delta_v = tracer_ratio_to_delta(itrc, vapor_ratio)
      delta_es = delta_x + eps_k + eps_eq + h * (delta_v - eps_k - delta_x)

      transp_moles = transp_water * water_moles_per_mm
      transp_moles_s = transp_moles / deltim
      liquid_diff = tracer_leaf_liquid_diffusivity(itrc, tk)
      IF (liquid_diff > trc_tiny .and. DEF_TRACER_NSS_LEAF_PATH_LENGTH > 0._r8) THEN
         peclet_number = transp_moles_s * DEF_TRACER_NSS_LEAF_PATH_LENGTH / &
            (liquid_water_molar_density * liquid_diff)
         IF (peclet_number > 1.e-8_r8) THEN
            new_peclet = (1._r8 - exp(-peclet_number)) / peclet_number
         ELSE
            new_peclet = 1._r8 - 0.5_r8 * peclet_number
         ENDIF
      ENDIF
      new_peclet = min(max(new_peclet, 0._r8), 1._r8)

	      IF (prev_leaf_moles > trc_tiny) THEN
	         leaf_moles_ratio = prev_leaf_moles / max(leaf_moles, trc_tiny)
	         IF (leaf_moles_ratio > 0.25_r8 .and. leaf_moles_ratio < 4._r8) THEN
	            prev_w = prev_leaf_moles
	         ELSE
	            ! Ignore a storage-size jump from changed NSS parameters or legacy restart files.
	            prev_w = leaf_moles
	         ENDIF
	         prev_e = prev_delta_e
	         prev_p = min(max(prev_peclet, 0._r8), 1._r8)
      ELSE
         prev_w = leaf_moles
         prev_e = delta_x
         prev_p = new_peclet
      ENDIF

	      ! Farquhar-Cernusak NSS form: the turnover time is controlled by
	      ! the one-way gross flux through stomata/boundary (g_t * w_i), not
	      ! by net transpiration when canopy coupling or water stress makes
	      ! E/(1-h) inconsistent with the same humidity used here.
	      gross_moles = transp_moles / one_minus_h
	      conductance_gross_moles = 0._r8
	      total_resistance = max(stomatal_resistance, 0._r8) + max(DEF_TRACER_NSS_LEAF_RB, 0._r8)
	      IF (total_resistance > trc_tiny .and. psrf > trc_tiny) THEN
	         vapor_molar_density_sat = tracer_saturation_vapor_pressure(tk, .false.) / &
	            (universal_gas_constant * tk)
	         conductance_gross_moles = deltim * vapor_molar_density_sat / total_resistance
	      ENDIF
	      gross_moles = max(gross_moles, conductance_gross_moles)
	      relax_b = alpha_k * alpha_eq / max(gross_moles, trc_tiny)
      denom = 1._r8 + relax_b * leaf_moles * new_peclet
      IF (denom > trc_tiny) THEN
         new_delta_e = (delta_es + relax_b * (leaf_moles * new_peclet * delta_x + &
            prev_w * prev_p * (prev_e - delta_x))) / denom
      ELSE
         new_delta_e = delta_es
      ENDIF
      new_delta_b = delta_x + new_peclet * (new_delta_e - delta_x)

      ! Diagnose the transpiration flux from the leaf isotope storage
      ! budget, not directly from the Craig-Gordon end-member. Over one
      ! step:
      !   E * R_T = E * R_x - d(W_leaf * R_leaf_bulk)
      ! This keeps daily transpiration close to source water whenever
      ! leaf water storage is approximately cyclic; any departure is then
      ! the explicit net storage tendency rather than a hidden flux bias.
      IF (prev_leaf_moles > trc_tiny) THEN
         prev_leaf_water = prev_w / water_moles_per_mm
         prev_bulk_ratio = tracer_delta_to_ratio(itrc, prev_delta_b)
      ELSE
         prev_leaf_water = leaf_moles / water_moles_per_mm
         prev_bulk_ratio = source_ratio
      ENDIF
      new_leaf_water = leaf_moles / water_moles_per_mm
      new_bulk_ratio = tracer_delta_to_ratio(itrc, new_delta_b)
	      storage_tracer_change = new_leaf_water * new_bulk_ratio - prev_leaf_water * prev_bulk_ratio
	      storage_tracer_change_used = storage_tracer_change
	      storage_scale = max(transp_water * source_ratio, trc_tiny)
	      storage_bound = max_storage_tendency_fraction * storage_scale
	      IF (abs(storage_tracer_change) > max_storage_jump_for_fallback * storage_scale) THEN
	         ! A tiny transpiration step cannot physically determine a large
	         ! leaf-storage isotope tendency (usually LAI/restart/NSS memory
	         ! mismatch). Freeze the budget tendency for this step and adjust
	         ! the diagnostic bulk leaf delta to the same storage state.
	         storage_tracer_change_used = 0._r8
	      ELSE
	         storage_tracer_change_used = min(max(storage_tracer_change_used, -storage_bound), &
	            storage_bound)
	      ENDIF
	      IF (abs(storage_tracer_change_used - storage_tracer_change) > trc_tiny) THEN
	         target_leaf_storage = prev_leaf_water * prev_bulk_ratio + storage_tracer_change_used
	         IF (new_leaf_water > trc_tiny .and. target_leaf_storage > trc_tiny) THEN
	            new_bulk_ratio = target_leaf_storage / new_leaf_water
	            new_delta_b = tracer_ratio_to_delta(itrc, new_bulk_ratio)
	            IF (new_peclet > 1.e-6_r8) THEN
	               new_delta_e = delta_x + (new_delta_b - delta_x) / new_peclet
	            ELSE
	               new_delta_e = new_delta_b
	            ENDIF
	         ELSE
	            storage_tracer_change_used = 0._r8
	            new_delta_e = delta_x
	            new_delta_b = delta_x
	         ENDIF
	      ENDIF
	      trans_ratio = (transp_water * source_ratio - storage_tracer_change_used) / transp_water
	      IF (trans_ratio /= trans_ratio .or. trans_ratio <= 0._r8) THEN
	         trans_ratio = source_ratio
	         new_delta_e = prev_delta_e
	         new_delta_b = prev_delta_b
	      ENDIF
      IF (new_delta_b /= new_delta_b) new_delta_b = prev_delta_b
   END SUBROUTINE tracer_transpiration_nss_ratio

   real(r8) FUNCTION tracer_surface_relhum (qair, psrf, temp_k, over_ice)
      real(r8), intent(in) :: qair
      real(r8), intent(in) :: psrf
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: over_ice
      real(r8) :: eair, esat, qsafe, psafe

      qsafe = max(qair, 0._r8)
      psafe = max(psrf, 1._r8)
      eair = qsafe * psafe / (0.622_r8 + 0.378_r8 * qsafe)
      esat = tracer_saturation_vapor_pressure(temp_k, over_ice)
      IF (esat > trc_tiny) THEN
         tracer_surface_relhum = min(max(eair / esat, 0._r8), 0.999999_r8)
      ELSE
         tracer_surface_relhum = 0._r8
      ENDIF
   END FUNCTION tracer_surface_relhum

   real(r8) FUNCTION tracer_saturation_vapor_pressure (temp_k, over_ice)
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: over_ice
      real(r8) :: tc, tk

      tk = max(temp_k, 150._r8)
      tc = tk - tfrz
      IF (over_ice) THEN
         tracer_saturation_vapor_pressure = 611.2_r8 * exp(22.46_r8 * tc / (272.62_r8 + tc))
      ELSE
         tracer_saturation_vapor_pressure = 611.2_r8 * exp(17.67_r8 * tc / (243.5_r8 + tc))
      ENDIF
   END FUNCTION tracer_saturation_vapor_pressure

   real(r8) FUNCTION tracer_leaf_kinetic_epsilon (itrc, ra, rb, rc)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: ra, rb, rc
      real(r8) :: ra1, rb1, rc1, denom, eps_stomatal, eps_boundary

      ra1 = max(ra, 0._r8)
      rb1 = max(rb, 0._r8)
      rc1 = max(rc, 0._r8)
      denom = ra1 + rb1 + rc1
      tracer_leaf_kinetic_epsilon = 0._r8
      IF (denom <= trc_tiny) RETURN

      SELECT CASE (tracer_isotope_kind(itrc))
      CASE (iso_o18)
         eps_stomatal = 28._r8
         eps_boundary = 19._r8
      CASE (iso_d)
         eps_stomatal = 25._r8
         eps_boundary = 17._r8
      CASE DEFAULT
         eps_stomatal = 0._r8
         eps_boundary = 0._r8
      END SELECT
      tracer_leaf_kinetic_epsilon = (eps_boundary * rb1 + eps_stomatal * rc1) / denom
   END FUNCTION tracer_leaf_kinetic_epsilon

   real(r8) FUNCTION tracer_leaf_liquid_diffusivity (itrc, temp_k)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      SELECT CASE (tracer_isotope_kind(itrc))
      CASE (iso_d)
         tracer_leaf_liquid_diffusivity = 116.e-9_r8 * exp(-626._r8 / max(tk - 139._r8, 1._r8))
      CASE DEFAULT
         tracer_leaf_liquid_diffusivity = 119.e-9_r8 * exp(-637._r8 / max(tk - 137._r8, 1._r8))
      END SELECT
   END FUNCTION tracer_leaf_liquid_diffusivity

   real(r8) FUNCTION tracer_ratio_to_delta (itrc, ratio)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: ratio
      real(r8) :: ref_ratio

      tracer_ratio_to_delta = 0._r8
      ref_ratio = 0._r8
      IF (itrc >= 1 .and. itrc <= ntracers .and. allocated(tracers)) ref_ratio = tracers(itrc)%ref_ratio
      IF (ref_ratio > trc_tiny) tracer_ratio_to_delta = (ratio / ref_ratio - 1._r8) * 1000._r8
   END FUNCTION tracer_ratio_to_delta

   real(r8) FUNCTION tracer_delta_to_ratio (itrc, delta)
      integer,  intent(in) :: itrc
      real(r8), intent(in) :: delta
      real(r8) :: ref_ratio

      tracer_delta_to_ratio = 0._r8
      ref_ratio = 0._r8
      IF (itrc >= 1 .and. itrc <= ntracers .and. allocated(tracers)) ref_ratio = tracers(itrc)%ref_ratio
      IF (ref_ratio > trc_tiny) tracer_delta_to_ratio = ref_ratio * (1._r8 + delta / 1000._r8)
   END FUNCTION tracer_delta_to_ratio

   integer FUNCTION tracer_isotope_kind (itrc)
      integer, intent(in) :: itrc
      character(len=32) :: lname

      tracer_isotope_kind = iso_none
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      IF (.not. allocated(tracers)) RETURN

      lname = lower_string(trim(tracers(itrc)%name))
      IF (index(lname, '18o') > 0 .or. index(lname, 'o18') > 0) THEN
         tracer_isotope_kind = iso_o18
      ELSEIF (index(lname, 'hdo') > 0 .or. index(lname, '2h') > 0 .or. &
              index(lname, 'deuter') > 0) THEN
         tracer_isotope_kind = iso_d
      ELSEIF (tracers(itrc)%ref_ratio > 0._r8) THEN
         IF (abs(tracers(itrc)%ref_ratio - Rsmow_18O) / Rsmow_18O < 0.1_r8) THEN
            tracer_isotope_kind = iso_o18
         ELSEIF (abs(tracers(itrc)%ref_ratio - Rsmow_D) / Rsmow_D < 0.1_r8) THEN
            tracer_isotope_kind = iso_d
         ENDIF
      ENDIF
   END FUNCTION tracer_isotope_kind

   FUNCTION lower_string (raw) RESULT(out)
      character(len=*), intent(in) :: raw
      character(len=len(raw)) :: out
      integer :: i, code

      out = raw
      DO i = 1, len(raw)
         code = iachar(raw(i:i))
         IF (code >= iachar('A') .and. code <= iachar('Z')) THEN
            out(i:i) = achar(code + 32)
         ENDIF
      ENDDO
   END FUNCTION lower_string

END MODULE MOD_Tracer_Frac
#endif
