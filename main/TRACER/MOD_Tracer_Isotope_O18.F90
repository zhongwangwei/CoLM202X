#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Isotope_O18

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: Rsmow_18O
   USE MOD_Tracer_Isotope_Registry, only: register_isotope_physics

   IMPLICIT NONE
   SAVE
   PRIVATE

   real(r8), parameter :: o18_diffusivity_ratio_air = 1.03189_r8

   PUBLIC :: register_o18_isotope_physics

CONTAINS

   SUBROUTINE register_o18_isotope_physics ()
      CALL register_isotope_physics(name='O18', name_patterns='18o,o18', &
         ref_ratio_hint=Rsmow_18O, legacy_forcing_kind=1, &
         default_soil_init_varname='soilwat_O18', &
         alpha_liq_vap_fn=o18_alpha_liq_vap, &
         alpha_ice_vap_fn=o18_alpha_ice_vap, &
         diffusivity_ratio_air_fn=o18_diffusivity_ratio, &
         leaf_kinetic_epsilon_fn=o18_leaf_kinetic_epsilon, &
         leaf_liquid_diffusivity_fn=o18_leaf_liquid_diffusivity)
   END SUBROUTINE register_o18_isotope_physics

   real(r8) FUNCTION o18_alpha_liq_vap (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      o18_alpha_liq_vap = exp(1137._r8 / (tk * tk) - 0.4156_r8 / tk - 0.0020667_r8)
   END FUNCTION o18_alpha_liq_vap

   real(r8) FUNCTION o18_alpha_ice_vap (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      o18_alpha_ice_vap = exp(11.839_r8 / tk - 0.028224_r8)
   END FUNCTION o18_alpha_ice_vap

   real(r8) FUNCTION o18_diffusivity_ratio ()
      o18_diffusivity_ratio = o18_diffusivity_ratio_air
   END FUNCTION o18_diffusivity_ratio

   real(r8) FUNCTION o18_leaf_kinetic_epsilon (ra, rb, rc)
      real(r8), intent(in) :: ra, rb, rc
      o18_leaf_kinetic_epsilon = weighted_leaf_epsilon(ra, rb, rc, 19._r8, 28._r8)
   END FUNCTION o18_leaf_kinetic_epsilon

   real(r8) FUNCTION o18_leaf_liquid_diffusivity (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      o18_leaf_liquid_diffusivity = 119.e-9_r8 * exp(-637._r8 / max(tk - 137._r8, 1._r8))
   END FUNCTION o18_leaf_liquid_diffusivity

   real(r8) FUNCTION weighted_leaf_epsilon (ra, rb, rc, eps_boundary, eps_stomatal)
      real(r8), intent(in) :: ra, rb, rc, eps_boundary, eps_stomatal
      real(r8) :: ra1, rb1, rc1, denom

      ra1 = max(ra, 0._r8)
      rb1 = max(rb, 0._r8)
      rc1 = max(rc, 0._r8)
      denom = ra1 + rb1 + rc1
      weighted_leaf_epsilon = 0._r8
      IF (denom <= 0._r8) RETURN
      weighted_leaf_epsilon = (eps_boundary * rb1 + eps_stomatal * rc1) / denom
   END FUNCTION weighted_leaf_epsilon

END MODULE MOD_Tracer_Isotope_O18
#endif
