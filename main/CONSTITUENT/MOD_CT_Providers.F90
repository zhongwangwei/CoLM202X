#include <define.h>

#ifdef ConstituentTransport
MODULE MOD_CT_Providers
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  The lists of compiled constituent species code:
!
!  - isotope fractionation physics (ensure_isotope_physics_registered);
!  - species that own state or phase callbacks (ct_register_providers,
!    below the module).
!
!  Parameter-only isotopes and solutes are not listed here. Add a line only
!  when the species code exists.
!-----------------------------------------------------------------------

   USE MOD_CT_Isotope_O18, only: register_o18_isotope_physics
   USE MOD_CT_Isotope_HDO, only: register_hdo_isotope_physics
   IMPLICIT NONE
   SAVE
   PRIVATE

   logical :: isotope_physics_registered = .false.

   PUBLIC :: ensure_isotope_physics_registered
   PUBLIC :: register_all_isotope_physics

CONTAINS

   SUBROUTINE ensure_isotope_physics_registered ()
      IF (isotope_physics_registered) RETURN
      CALL register_all_isotope_physics ()
      isotope_physics_registered = .true.
   END SUBROUTINE ensure_isotope_physics_registered

   SUBROUTINE register_all_isotope_physics ()
      CALL register_o18_isotope_physics ()
      CALL register_hdo_isotope_physics ()
   END SUBROUTINE register_all_isotope_physics

END MODULE MOD_CT_Providers


SUBROUTINE ct_register_providers ()
!-----------------------------------------------------------------------
!  Species with lifecycle callbacks. An external subroutine rather than a
!  module procedure because the providers USE MOD_CT_Lifecycle, which calls
!  this; no species is attached yet.
!-----------------------------------------------------------------------
   IMPLICIT NONE
END SUBROUTINE ct_register_providers
#endif
