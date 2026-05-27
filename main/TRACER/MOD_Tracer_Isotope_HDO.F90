#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Isotope_HDO

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: Rsmow_D
   USE MOD_Tracer_Isotope_Registry, only: register_isotope_physics

   IMPLICIT NONE
   SAVE
   PRIVATE

   real(r8), parameter :: hdo_diffusivity_ratio_air = 1.01636_r8

   PUBLIC :: register_hdo_isotope_physics

CONTAINS

   SUBROUTINE register_hdo_isotope_physics ()
      CALL register_isotope_physics(name='HDO', name_patterns='hdo,2h,deuter,=h2', &
         ref_ratio_hint=Rsmow_D, legacy_forcing_kind=2, &
         default_soil_init_varname='soilwat_H2', &
         alpha_liq_vap_fn=hdo_alpha_liq_vap, &
         alpha_ice_vap_fn=hdo_alpha_ice_vap, &
         diffusivity_ratio_air_fn=hdo_diffusivity_ratio, &
         leaf_kinetic_epsilon_fn=hdo_leaf_kinetic_epsilon, &
         leaf_liquid_diffusivity_fn=hdo_leaf_liquid_diffusivity)
   END SUBROUTINE register_hdo_isotope_physics

   real(r8) FUNCTION hdo_alpha_liq_vap (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      hdo_alpha_liq_vap = exp(24844._r8 / (tk * tk) - 76.248_r8 / tk + 0.052612_r8)
   END FUNCTION hdo_alpha_liq_vap

   real(r8) FUNCTION hdo_alpha_ice_vap (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      hdo_alpha_ice_vap = exp(16289._r8 / (tk * tk) - 0.0945_r8)
   END FUNCTION hdo_alpha_ice_vap

   real(r8) FUNCTION hdo_diffusivity_ratio ()
      hdo_diffusivity_ratio = hdo_diffusivity_ratio_air
   END FUNCTION hdo_diffusivity_ratio

   real(r8) FUNCTION hdo_leaf_kinetic_epsilon (ra, rb, rc)
      real(r8), intent(in) :: ra, rb, rc
      hdo_leaf_kinetic_epsilon = weighted_leaf_epsilon(ra, rb, rc, 17._r8, 25._r8)
   END FUNCTION hdo_leaf_kinetic_epsilon

   real(r8) FUNCTION hdo_leaf_liquid_diffusivity (temp_k)
      real(r8), intent(in) :: temp_k
      real(r8) :: tk

      tk = max(temp_k, 150._r8)
      hdo_leaf_liquid_diffusivity = 116.e-9_r8 * exp(-626._r8 / max(tk - 139._r8, 1._r8))
   END FUNCTION hdo_leaf_liquid_diffusivity

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

END MODULE MOD_Tracer_Isotope_HDO
#endif
