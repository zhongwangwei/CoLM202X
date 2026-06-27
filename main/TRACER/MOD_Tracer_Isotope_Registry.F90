#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Isotope_Registry

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: tracers, ntracers, tracer_is_isotope, trc_tiny, tracer_lower
   USE MOD_SPMD_Task, only: CoLM_stop

   IMPLICIT NONE
   SAVE
   PRIVATE

   integer, parameter :: isotope_name_len = 32
   integer, parameter :: isotope_patterns_len = 160
   integer, parameter :: isotope_varname_len = 64

   ABSTRACT INTERFACE
      real(r8) FUNCTION isotope_alpha_temp_if (temp_k)
         IMPORT :: r8
         real(r8), intent(in) :: temp_k
      END FUNCTION isotope_alpha_temp_if

      real(r8) FUNCTION isotope_scalar_if ()
         IMPORT :: r8
      END FUNCTION isotope_scalar_if

      real(r8) FUNCTION isotope_leaf_kinetic_epsilon_if (ra, rb, rc)
         IMPORT :: r8
         real(r8), intent(in) :: ra, rb, rc
      END FUNCTION isotope_leaf_kinetic_epsilon_if
   END INTERFACE

   TYPE :: isotope_physics_type
      character(len=isotope_name_len) :: name = ''
      character(len=isotope_patterns_len) :: name_patterns = ''
      real(r8) :: ref_ratio_hint = 0._r8
      real(r8) :: ref_ratio_tolerance = 0.1_r8
      integer :: legacy_forcing_kind = 0
      character(len=isotope_varname_len) :: default_soil_init_varname = 'null'
      procedure(isotope_alpha_temp_if), pointer, nopass :: alpha_liq_vap => null()
      procedure(isotope_alpha_temp_if), pointer, nopass :: alpha_ice_vap => null()
      procedure(isotope_scalar_if), pointer, nopass :: diffusivity_ratio_air => null()
      procedure(isotope_leaf_kinetic_epsilon_if), pointer, nopass :: leaf_kinetic_epsilon => null()
      procedure(isotope_alpha_temp_if), pointer, nopass :: leaf_liquid_diffusivity => null()
   END TYPE isotope_physics_type

   TYPE(isotope_physics_type), allocatable :: isotope_physics(:)
   integer :: n_isotope_physics = 0

   PUBLIC :: isotope_physics_type
   PUBLIC :: register_isotope_physics
   PUBLIC :: find_isotope_physics
   PUBLIC :: isotope_fractionation_registered
   PUBLIC :: isotope_alpha_liq_vap
   PUBLIC :: isotope_alpha_ice_vap
   PUBLIC :: isotope_diffusivity_ratio_air
   PUBLIC :: isotope_leaf_kinetic_epsilon
   PUBLIC :: isotope_leaf_liquid_diffusivity
   PUBLIC :: isotope_legacy_forcing_kind
   PUBLIC :: isotope_default_soil_init_varname
   PUBLIC :: isotope_weighted_leaf_epsilon

CONTAINS

   real(r8) FUNCTION isotope_weighted_leaf_epsilon (ra, rb, rc, eps_boundary, eps_stomatal)
      real(r8), intent(in) :: ra, rb, rc, eps_boundary, eps_stomatal
      real(r8) :: ra1, rb1, rc1, denom

      ra1 = max(ra, 0._r8)
      rb1 = max(rb, 0._r8)
      rc1 = max(rc, 0._r8)
      denom = ra1 + rb1 + rc1
      isotope_weighted_leaf_epsilon = 0._r8
      IF (denom <= 0._r8) RETURN
      isotope_weighted_leaf_epsilon = (eps_boundary * rb1 + eps_stomatal * rc1) / denom
   END FUNCTION isotope_weighted_leaf_epsilon

   SUBROUTINE register_isotope_physics (name, name_patterns, ref_ratio_hint, &
      legacy_forcing_kind, default_soil_init_varname, ref_ratio_tolerance, &
      alpha_liq_vap_fn, alpha_ice_vap_fn, diffusivity_ratio_air_fn, &
      leaf_kinetic_epsilon_fn, leaf_liquid_diffusivity_fn)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: name_patterns
      real(r8), intent(in), optional :: ref_ratio_hint
      integer, intent(in), optional :: legacy_forcing_kind
      character(len=*), intent(in), optional :: default_soil_init_varname
      real(r8), intent(in), optional :: ref_ratio_tolerance
      procedure(isotope_alpha_temp_if), optional :: alpha_liq_vap_fn
      procedure(isotope_alpha_temp_if), optional :: alpha_ice_vap_fn
      procedure(isotope_scalar_if), optional :: diffusivity_ratio_air_fn
      procedure(isotope_leaf_kinetic_epsilon_if), optional :: leaf_kinetic_epsilon_fn
      procedure(isotope_alpha_temp_if), optional :: leaf_liquid_diffusivity_fn
      integer :: idx

      IF (len_trim(name) <= 0) THEN
         CALL CoLM_stop ('MOD_Tracer_Isotope_Registry: cannot register isotope physics with empty name')
      ENDIF
      IF (find_registered_isotope_by_name(name) > 0) THEN
         CALL CoLM_stop ('MOD_Tracer_Isotope_Registry: duplicate isotope physics registration')
      ENDIF

      CALL ensure_isotope_physics_capacity(n_isotope_physics + 1)
      n_isotope_physics = n_isotope_physics + 1
      idx = n_isotope_physics

      isotope_physics(idx)%name = trim(name)
      isotope_physics(idx)%name_patterns = trim(name_patterns)
      IF (present(ref_ratio_hint)) isotope_physics(idx)%ref_ratio_hint = ref_ratio_hint
      IF (present(ref_ratio_tolerance)) isotope_physics(idx)%ref_ratio_tolerance = ref_ratio_tolerance
      IF (present(legacy_forcing_kind)) isotope_physics(idx)%legacy_forcing_kind = legacy_forcing_kind
      IF (present(default_soil_init_varname)) &
         isotope_physics(idx)%default_soil_init_varname = trim(default_soil_init_varname)
      IF (present(alpha_liq_vap_fn)) isotope_physics(idx)%alpha_liq_vap => alpha_liq_vap_fn
      IF (present(alpha_ice_vap_fn)) isotope_physics(idx)%alpha_ice_vap => alpha_ice_vap_fn
      IF (present(diffusivity_ratio_air_fn)) &
         isotope_physics(idx)%diffusivity_ratio_air => diffusivity_ratio_air_fn
      IF (present(leaf_kinetic_epsilon_fn)) &
         isotope_physics(idx)%leaf_kinetic_epsilon => leaf_kinetic_epsilon_fn
      IF (present(leaf_liquid_diffusivity_fn)) &
         isotope_physics(idx)%leaf_liquid_diffusivity => leaf_liquid_diffusivity_fn
   END SUBROUTINE register_isotope_physics

   SUBROUTINE ensure_isotope_physics_capacity (need)
      integer, intent(in) :: need
      TYPE(isotope_physics_type), allocatable :: grown(:)
      integer :: old_size, new_size

      IF (.not. allocated(isotope_physics)) THEN
         new_size = max(need, 2)
         allocate(isotope_physics(new_size))
         RETURN
      ENDIF
      old_size = size(isotope_physics)
      IF (need <= old_size) RETURN

      new_size = max(need, old_size * 2)
      allocate(grown(new_size))
      IF (n_isotope_physics > 0) grown(1:n_isotope_physics) = isotope_physics(1:n_isotope_physics)
      CALL move_alloc(grown, isotope_physics)
   END SUBROUTINE ensure_isotope_physics_capacity

   integer FUNCTION find_registered_isotope_by_name (name)
      character(len=*), intent(in) :: name
      integer :: i

      find_registered_isotope_by_name = 0
      IF (.not. allocated(isotope_physics)) RETURN
      DO i = 1, n_isotope_physics
         IF (trim(tracer_lower(isotope_physics(i)%name)) == trim(tracer_lower(name))) THEN
            find_registered_isotope_by_name = i
            RETURN
         ENDIF
      ENDDO
   END FUNCTION find_registered_isotope_by_name

   integer FUNCTION find_isotope_physics (itrc)
      integer, intent(in) :: itrc
      character(len=32) :: lname
      integer :: i

      find_isotope_physics = 0
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      IF (.not. allocated(tracers)) RETURN
      IF (.not. tracer_is_isotope(itrc)) RETURN
      IF (.not. allocated(isotope_physics)) RETURN

      lname = tracer_lower(trim(tracers(itrc)%name))
      DO i = 1, n_isotope_physics
         IF (isotope_name_matches(lname, isotope_physics(i)%name_patterns)) THEN
            find_isotope_physics = i
            RETURN
         ENDIF
      ENDDO

      IF (tracers(itrc)%ref_ratio > 0._r8) THEN
         DO i = 1, n_isotope_physics
            IF (isotope_physics(i)%ref_ratio_hint > trc_tiny .and. &
                abs(tracers(itrc)%ref_ratio - isotope_physics(i)%ref_ratio_hint) / &
                isotope_physics(i)%ref_ratio_hint < isotope_physics(i)%ref_ratio_tolerance) THEN
               find_isotope_physics = i
               RETURN
            ENDIF
         ENDDO
      ENDIF
   END FUNCTION find_isotope_physics

   logical FUNCTION isotope_fractionation_registered (itrc)
      integer, intent(in) :: itrc

      isotope_fractionation_registered = find_isotope_physics(itrc) > 0
   END FUNCTION isotope_fractionation_registered

   real(r8) FUNCTION isotope_alpha_liq_vap (itrc, temp_k)
      integer, intent(in) :: itrc
      real(r8), intent(in) :: temp_k
      integer :: idx

      isotope_alpha_liq_vap = 1._r8
      idx = find_isotope_physics(itrc)
      IF (idx > 0 .and. associated(isotope_physics(idx)%alpha_liq_vap)) &
         isotope_alpha_liq_vap = isotope_physics(idx)%alpha_liq_vap(temp_k)
   END FUNCTION isotope_alpha_liq_vap

   real(r8) FUNCTION isotope_alpha_ice_vap (itrc, temp_k)
      integer, intent(in) :: itrc
      real(r8), intent(in) :: temp_k
      integer :: idx

      isotope_alpha_ice_vap = 1._r8
      idx = find_isotope_physics(itrc)
      IF (idx > 0 .and. associated(isotope_physics(idx)%alpha_ice_vap)) &
         isotope_alpha_ice_vap = isotope_physics(idx)%alpha_ice_vap(temp_k)
   END FUNCTION isotope_alpha_ice_vap

   real(r8) FUNCTION isotope_diffusivity_ratio_air (itrc)
      integer, intent(in) :: itrc
      integer :: idx

      isotope_diffusivity_ratio_air = 1._r8
      idx = find_isotope_physics(itrc)
      IF (idx > 0 .and. associated(isotope_physics(idx)%diffusivity_ratio_air)) &
         isotope_diffusivity_ratio_air = isotope_physics(idx)%diffusivity_ratio_air()
   END FUNCTION isotope_diffusivity_ratio_air

   real(r8) FUNCTION isotope_leaf_kinetic_epsilon (itrc, ra, rb, rc)
      integer, intent(in) :: itrc
      real(r8), intent(in) :: ra, rb, rc
      integer :: idx

      isotope_leaf_kinetic_epsilon = 0._r8
      idx = find_isotope_physics(itrc)
      IF (idx > 0 .and. associated(isotope_physics(idx)%leaf_kinetic_epsilon)) &
         isotope_leaf_kinetic_epsilon = isotope_physics(idx)%leaf_kinetic_epsilon(ra, rb, rc)
   END FUNCTION isotope_leaf_kinetic_epsilon

   real(r8) FUNCTION isotope_leaf_liquid_diffusivity (itrc, temp_k)
      integer, intent(in) :: itrc
      real(r8), intent(in) :: temp_k
      integer :: idx

      isotope_leaf_liquid_diffusivity = 0._r8
      idx = find_isotope_physics(itrc)
      IF (idx > 0 .and. associated(isotope_physics(idx)%leaf_liquid_diffusivity)) &
         isotope_leaf_liquid_diffusivity = isotope_physics(idx)%leaf_liquid_diffusivity(temp_k)
   END FUNCTION isotope_leaf_liquid_diffusivity

   integer FUNCTION isotope_legacy_forcing_kind (itrc)
      integer, intent(in) :: itrc
      integer :: idx

      isotope_legacy_forcing_kind = 0
      idx = find_isotope_physics(itrc)
      IF (idx > 0) isotope_legacy_forcing_kind = isotope_physics(idx)%legacy_forcing_kind
   END FUNCTION isotope_legacy_forcing_kind

   SUBROUTINE isotope_default_soil_init_varname (itrc, varname)
      integer, intent(in) :: itrc
      character(len=*), intent(out) :: varname
      integer :: idx

      varname = 'null'
      idx = find_isotope_physics(itrc)
      IF (idx > 0) varname = trim(isotope_physics(idx)%default_soil_init_varname)
   END SUBROUTINE isotope_default_soil_init_varname

   logical FUNCTION isotope_name_matches (lname, patterns)
      character(len=*), intent(in) :: lname
      character(len=*), intent(in) :: patterns
      character(len=isotope_patterns_len) :: rest
      character(len=32) :: token
      integer :: comma

      isotope_name_matches = .false.
      rest = tracer_lower(trim(patterns))
      DO WHILE (len_trim(rest) > 0)
         comma = index(rest, ',')
         IF (comma > 0) THEN
            token = adjustl(rest(1:comma-1))
            rest = adjustl(rest(comma+1:))
         ELSE
            token = adjustl(rest)
            rest = ''
         ENDIF
         IF (len_trim(token) > 0) THEN
            IF (token(1:1) == '=') THEN
               IF (trim(lname) == trim(token(2:))) THEN
                  isotope_name_matches = .true.
                  RETURN
               ENDIF
            ELSEIF (index(lname, trim(token)) > 0) THEN
               isotope_name_matches = .true.
               RETURN
            ENDIF
         ENDIF
      ENDDO
   END FUNCTION isotope_name_matches


END MODULE MOD_Tracer_Isotope_Registry
#endif
