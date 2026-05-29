#include <define.h>

#ifdef TRACER
SUBROUTINE register_all_particle_callbacks ()
!=======================================================================
! Central particle-species registration list.
!
! Generic dispatch in MOD_Tracer_Particle does not name individual species.
! Adding a new particle species should add its registration here (and the
! Makefile/object list), not in MOD_Tracer_Particle.F90.
!=======================================================================

#define TRACER_PARTICLE_SPECIES(module_name, register_fn) USE module_name, only: register_fn
#include "tracer_particle_species.inc"
#undef TRACER_PARTICLE_SPECIES

   IMPLICIT NONE

#define TRACER_PARTICLE_SPECIES(module_name, register_fn) CALL register_fn ()
#include "tracer_particle_species.inc"
#undef TRACER_PARTICLE_SPECIES

END SUBROUTINE register_all_particle_callbacks
#endif
