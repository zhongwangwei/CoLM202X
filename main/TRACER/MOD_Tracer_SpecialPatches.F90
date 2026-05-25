#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_SpecialPatches

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, trc_tiny, tracer_init_water_ratio, &
      tracer_can_use_fixed_signature
   USE MOD_Tracer_Forcing, only: tracer_forcing_precip_value, tracer_forcing_vapor_value
   USE MOD_Tracer_Frac, only: tracer_fractionation_active, tracer_surface_relhum, &
      tracer_diffusivity_ratio_air, tracer_craig_gordon_evap_ratio, &
      tracer_equilibrium_deposition_ratio
   USE MOD_Tracer_Conservation, only: tracer_save_storage, tracer_balance_check, &
      tracer_apply_reactive_processes
   USE MOD_Tracer_Hist, only: tracer_hist_accumulate
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, trc_scv, &
      trc_ldew_rain, trc_ldew_snow, trc_rnof_step, a_trc_precip, a_trc_evap, &
      a_trc_rsur, a_trc_rnof, trc_wetwat, trc_waterstorage, trc_storage_beg, &
      trc_runtime_forced, sync_tracer_patch_ratio

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: tracer_glacier_patch
   PUBLIC :: tracer_waterbody_patch

CONTAINS

   SUBROUTINE tracer_glacier_patch(ipatch, maxsnl, nl_soil, deltim, &
      prc_rain, prl_rain, prc_snow, prl_snow, rnof, qseva, qsubl, qsdew, qfros, &
      endwb, totwb, glacier_overflow_mass, errorw, wdsrf, scv, t_grnd, forc_q, &
      forc_psrf, wliq_soisno, wice_soisno)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch, maxsnl, nl_soil
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: prc_rain, prl_rain, prc_snow, prl_snow
      real(r8), intent(in) :: rnof, qseva, qsubl, qsdew, qfros
      real(r8), intent(in) :: endwb, totwb, glacier_overflow_mass, errorw
      real(r8), intent(in) :: wdsrf, scv, t_grnd, forc_q, forc_psrf
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)

      integer  :: itrc, j_trc, snl_trc
      real(r8) :: R_init, R_precip, R_vapor, R_out, R_final
      real(r8) :: R_dew, R_frost, R_evap_liq, R_evap_ice, R_runoff
      real(r8) :: precip_mass, rnof_mass
      real(r8) :: evap_mass, dep_mass
      real(r8) :: evap_liq_mass, evap_ice_mass, dep_liq_mass, dep_ice_mass
      real(r8) :: water_dS, water_end, water_beg, water_input
      real(r8) :: water_before_output, water_after_evap
      real(r8) :: trc_input, trc_evap, trc_rnof, trc_available, trc_final
      real(r8) :: relhum_liq, relhum_ice, alpha_k, xerr_tracer
      logical  :: mixed_signature, fixed_signature, frac_active

      DO itrc = 1, ntracers
         trc_ldew_rain(itrc, ipatch) = 0._r8
         trc_ldew_snow(itrc, ipatch) = 0._r8
         trc_wetwat   (itrc, ipatch) = 0._r8
         IF (allocated(trc_waterstorage)) trc_waterstorage(itrc, ipatch) = 0._r8
      ENDDO

      snl_trc = 0
      CALL tracer_save_storage(ipatch, snl_trc, nl_soil)
      precip_mass = (prc_rain + prl_rain + prc_snow + prl_snow) * deltim
      rnof_mass   = max(rnof, 0._r8) * deltim
      evap_liq_mass = max(qseva, 0._r8) * deltim
      evap_ice_mass = max(qsubl, 0._r8) * deltim
      dep_liq_mass  = max(qsdew, 0._r8) * deltim
      dep_ice_mass  = max(qfros, 0._r8) * deltim
      evap_mass = evap_liq_mass + evap_ice_mass
      dep_mass  = dep_liq_mass  + dep_ice_mass
      water_input = precip_mass + dep_mass
      water_dS = endwb - totwb - glacier_overflow_mass
      water_end = max(wdsrf, 0._r8) + max(scv, 0._r8)
      DO j_trc = 1, nl_soil
         water_end = water_end + max(wliq_soisno(j_trc), 0._r8) + max(wice_soisno(j_trc), 0._r8)
      ENDDO
      water_beg = water_end - water_dS

      DO itrc = 1, ntracers
         trc_rnof_step(itrc, ipatch) = 0._r8
         R_init = tracer_init_water_ratio(itrc)
         frac_active = tracer_fractionation_active(itrc)
         fixed_signature = tracer_can_use_fixed_signature(itrc) .and. .not. frac_active
         IF (allocated(trc_runtime_forced)) THEN
            fixed_signature = fixed_signature .and. .not. trc_runtime_forced(itrc)
         ENDIF
         mixed_signature = .not. fixed_signature

         IF (mixed_signature) THEN
            R_precip = tracer_forcing_precip_value(itrc, ipatch)
            R_vapor  = tracer_forcing_vapor_value (itrc, ipatch)
            R_dew = R_vapor
            R_frost = R_vapor
            IF (frac_active) THEN
               R_dew = tracer_equilibrium_deposition_ratio(itrc, R_vapor, t_grnd, .false.)
               R_frost = tracer_equilibrium_deposition_ratio(itrc, R_vapor, t_grnd, .true.)
            ENDIF
            trc_input = precip_mass * R_precip + dep_liq_mass * R_dew + dep_ice_mass * R_frost
            trc_available = max(trc_storage_beg(itrc, ipatch) + trc_input, 0._r8)
            water_before_output = water_beg + water_input
            IF (water_before_output > trc_tiny) THEN
               R_out = trc_available / water_before_output
            ELSE
               R_out = R_init
            ENDIF
            R_evap_liq = R_out
            R_evap_ice = R_out
            IF (frac_active) THEN
               alpha_k = tracer_diffusivity_ratio_air(itrc)
               relhum_liq = tracer_surface_relhum(forc_q, forc_psrf, t_grnd, .false.)
               relhum_ice = tracer_surface_relhum(forc_q, forc_psrf, t_grnd, .true.)
               R_evap_liq = tracer_craig_gordon_evap_ratio(itrc, R_out, R_vapor, &
                  t_grnd, relhum_liq, alpha_k, .false.)
               R_evap_ice = tracer_craig_gordon_evap_ratio(itrc, R_out, R_vapor, &
                  t_grnd, relhum_ice, alpha_k, .true.)
               R_evap_liq = min(R_evap_liq, max(R_out, 0._r8))
               R_evap_ice = min(R_evap_ice, max(R_out, 0._r8))
            ENDIF
            trc_evap = min(evap_liq_mass * R_evap_liq + evap_ice_mass * R_evap_ice, trc_available)
            water_after_evap = water_before_output - evap_mass
            IF (water_after_evap > trc_tiny) THEN
               R_runoff = max(trc_available - trc_evap, 0._r8) / water_after_evap
            ELSE
               R_runoff = R_out
            ENDIF
            trc_rnof = min(rnof_mass * R_runoff, max(trc_available - trc_evap, 0._r8))
            trc_final = max(trc_available - trc_evap - trc_rnof, 0._r8)
            IF (water_end > trc_tiny) THEN
               R_final = trc_final / water_end
            ELSE
               R_final = 0._r8
            ENDIF
         ELSE
            trc_input = water_input * R_init
            trc_evap  = evap_mass * R_init
            trc_rnof  = rnof_mass * R_init
            R_final   = R_init
         ENDIF

         IF (trc_input > 0._r8) a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_input
         IF (trc_evap > 0._r8) a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_evap
         IF (trc_rnof > 0._r8) THEN
            trc_rnof_step(itrc, ipatch) = trc_rnof
            a_trc_rsur(itrc, ipatch) = a_trc_rsur(itrc, ipatch) + trc_rnof
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_rnof
         ENDIF
         CALL sync_tracer_patch_ratio(itrc, ipatch, snl_trc, maxsnl, nl_soil, &
            wliq_soisno, wice_soisno, 0._r8, wdsrf, scv, R_final)
      ENDDO

      CALL tracer_apply_reactive_processes(ipatch, snl_trc, nl_soil, deltim)
      CALL tracer_balance_check(ipatch, snl_trc, nl_soil, deltim, xerr_tracer, &
         patchtype_in = 3, water_err_in = errorw, water_dS_in = water_dS, &
         water_input_in = precip_mass + dep_mass, water_output_in = evap_mass + rnof_mass, &
         water_evap_in = evap_mass, water_rnof_in = rnof_mass)
      CALL tracer_hist_accumulate(ipatch, snl_trc, maxsnl, nl_soil, 0._r8, 0._r8, &
         wliq_soisno(snl_trc+1:nl_soil), wice_soisno(snl_trc+1:nl_soil), 0._r8, wdsrf, 0._r8, scv)

   END SUBROUTINE tracer_glacier_patch

   SUBROUTINE tracer_waterbody_patch(ipatch, maxsnl, nl_soil, snl, deltim, &
      forc_rain, forc_snow, lake_deficit, rnof, qseva, qsubl, qsdew, qfros, &
      endwb, totwb, errorw, wa, wdsrf, scv, t_grnd, forc_q, forc_psrf, &
      wliq_soisno, wice_soisno, use_dynamic_lake)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch, maxsnl, nl_soil, snl
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: forc_rain, forc_snow, lake_deficit, rnof
      real(r8), intent(in) :: qseva, qsubl, qsdew, qfros
      real(r8), intent(in) :: endwb, totwb, errorw, wa, wdsrf, scv, t_grnd, forc_q, forc_psrf
      real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)
      logical,  intent(in) :: use_dynamic_lake

      integer  :: itrc, j_trc
      real(r8) :: R_init, R_precip, R_vapor, R_pool, R_out, R_final
      real(r8) :: R_dew, R_frost, R_evap_liq, R_evap_ice, R_runoff
      real(r8) :: atm_precip_mass, deficit_mass, precip_mass, rnof_mass
      real(r8) :: evap_mass, dep_mass
      real(r8) :: evap_liq_mass, evap_ice_mass, dep_liq_mass, dep_ice_mass
      real(r8) :: water_dS, water_end, water_beg, water_input
      real(r8) :: water_before_output, water_after_evap
      real(r8) :: trc_input, trc_evap, trc_rnof, trc_available, trc_final
      real(r8) :: relhum_liq, relhum_ice, alpha_k, xerr_tracer
      logical  :: mixed_signature, fixed_signature, frac_active

      DO itrc = 1, ntracers
         trc_ldew_rain(itrc, ipatch) = 0._r8
         trc_ldew_snow(itrc, ipatch) = 0._r8
         trc_wetwat   (itrc, ipatch) = 0._r8
         IF (allocated(trc_waterstorage)) trc_waterstorage(itrc, ipatch) = 0._r8
      ENDDO

      CALL tracer_save_storage(ipatch, maxsnl, nl_soil)
      atm_precip_mass = (forc_rain + forc_snow) * deltim
      IF (use_dynamic_lake) THEN
         precip_mass = atm_precip_mass
      ELSE
         precip_mass = (forc_rain + forc_snow + lake_deficit) * deltim
      ENDIF
      deficit_mass = precip_mass - atm_precip_mass
      rnof_mass   = max(rnof, 0._r8) * deltim
      evap_liq_mass = max(qseva, 0._r8) * deltim
      evap_ice_mass = max(qsubl, 0._r8) * deltim
      dep_liq_mass  = max(qsdew, 0._r8) * deltim
      dep_ice_mass  = max(qfros, 0._r8) * deltim
      evap_mass = evap_liq_mass + evap_ice_mass
      dep_mass  = dep_liq_mass  + dep_ice_mass
      water_input = precip_mass + dep_mass
      water_dS = endwb - totwb
      water_end = wa + max(wdsrf, 0._r8)
      DO j_trc = maxsnl + 1, nl_soil
         IF (j_trc >= snl + 1) THEN
            water_end = water_end + max(wliq_soisno(j_trc), 0._r8) + max(wice_soisno(j_trc), 0._r8)
         ENDIF
      ENDDO
      IF (snl >= 0) water_end = water_end + max(scv, 0._r8)
      water_beg = water_end - water_dS

      DO itrc = 1, ntracers
         trc_rnof_step(itrc, ipatch) = 0._r8
         R_init = tracer_init_water_ratio(itrc)
         frac_active = tracer_fractionation_active(itrc)
         fixed_signature = tracer_can_use_fixed_signature(itrc) .and. .not. frac_active
         IF (allocated(trc_runtime_forced)) THEN
            fixed_signature = fixed_signature .and. .not. trc_runtime_forced(itrc)
         ENDIF
         mixed_signature = .not. fixed_signature

         IF (mixed_signature) THEN
            R_precip = tracer_forcing_precip_value(itrc, ipatch)
            R_vapor  = tracer_forcing_vapor_value (itrc, ipatch)
            IF (water_beg > trc_tiny) THEN
               R_pool = max(trc_storage_beg(itrc, ipatch), 0._r8) / water_beg
            ELSE
               R_pool = R_precip
            ENDIF
            R_dew = R_vapor
            R_frost = R_vapor
            IF (frac_active) THEN
               R_dew = tracer_equilibrium_deposition_ratio(itrc, R_vapor, t_grnd, .false.)
               R_frost = tracer_equilibrium_deposition_ratio(itrc, R_vapor, t_grnd, .true.)
            ENDIF
            trc_input = atm_precip_mass * R_precip + dep_liq_mass * R_dew + &
               dep_ice_mass * R_frost + deficit_mass * R_pool
            trc_available = max(trc_storage_beg(itrc, ipatch) + trc_input, 0._r8)
            water_before_output = water_beg + water_input
            IF (water_before_output > trc_tiny) THEN
               R_out = trc_available / water_before_output
            ELSE
               R_out = R_init
            ENDIF
            R_evap_liq = R_out
            R_evap_ice = R_out
            IF (frac_active) THEN
               alpha_k = tracer_diffusivity_ratio_air(itrc)
               relhum_liq = tracer_surface_relhum(forc_q, forc_psrf, t_grnd, .false.)
               relhum_ice = tracer_surface_relhum(forc_q, forc_psrf, t_grnd, .true.)
               R_evap_liq = tracer_craig_gordon_evap_ratio(itrc, R_out, R_vapor, &
                  t_grnd, relhum_liq, alpha_k, .false.)
               R_evap_ice = tracer_craig_gordon_evap_ratio(itrc, R_out, R_vapor, &
                  t_grnd, relhum_ice, alpha_k, .true.)
               R_evap_liq = min(R_evap_liq, max(R_out, 0._r8))
               R_evap_ice = min(R_evap_ice, max(R_out, 0._r8))
            ENDIF
            trc_evap = min(evap_liq_mass * R_evap_liq + evap_ice_mass * R_evap_ice, trc_available)
            water_after_evap = water_before_output - evap_mass
            IF (water_after_evap > trc_tiny) THEN
               R_runoff = max(trc_available - trc_evap, 0._r8) / water_after_evap
            ELSE
               R_runoff = R_out
            ENDIF
            trc_rnof = min(rnof_mass * R_runoff, max(trc_available - trc_evap, 0._r8))
            trc_final = max(trc_available - trc_evap - trc_rnof, 0._r8)
            IF (water_end > trc_tiny) THEN
               R_final = trc_final / water_end
            ELSE
               R_final = 0._r8
            ENDIF
         ELSE
            trc_input = water_input * R_init
            trc_evap  = evap_mass * R_init
            trc_rnof  = rnof_mass * R_init
            R_final   = R_init
         ENDIF

         IF (trc_input > 0._r8) a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_input
         IF (trc_evap > 0._r8) a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_evap
         IF (trc_rnof > 0._r8) THEN
            trc_rnof_step(itrc, ipatch) = trc_rnof
            a_trc_rsur(itrc, ipatch) = a_trc_rsur(itrc, ipatch) + trc_rnof
            a_trc_rnof(itrc, ipatch) = a_trc_rnof(itrc, ipatch) + trc_rnof
         ENDIF
         CALL sync_tracer_patch_ratio(itrc, ipatch, snl, maxsnl, nl_soil, &
            wliq_soisno, wice_soisno, wa, wdsrf, scv, R_final)
      ENDDO

      CALL tracer_apply_reactive_processes(ipatch, maxsnl, nl_soil, deltim)
      CALL tracer_balance_check(ipatch, maxsnl, nl_soil, deltim, xerr_tracer, &
         patchtype_in = 4, water_err_in = errorw, water_dS_in = water_dS, &
         water_input_in = precip_mass + dep_mass, water_output_in = evap_mass + rnof_mass, &
         water_evap_in = evap_mass, water_rnof_in = rnof_mass)
      CALL tracer_hist_accumulate(ipatch, snl, maxsnl, nl_soil, 0._r8, 0._r8, &
         wliq_soisno(snl+1:nl_soil), wice_soisno(snl+1:nl_soil), wa, wdsrf, 0._r8, scv)

   END SUBROUTINE tracer_waterbody_patch

END MODULE MOD_Tracer_SpecialPatches
#endif
