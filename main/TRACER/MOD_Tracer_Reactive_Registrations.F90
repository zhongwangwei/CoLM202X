#include <define.h>

#ifdef TRACER
SUBROUTINE register_all_reactive_callbacks ()
!=======================================================================
! Central reactive-species registration list.
!
! Generic dispatch in MOD_Tracer_Reactive does not name individual species.
! Adding a new reactive species should add its registration here (and the
! Makefile/object list), not in MOD_Tracer_Reactive.F90.
!=======================================================================

#ifdef BGC
#define TRACER_REACTIVE_SPECIES(module_name, register_fn) USE module_name, only: register_fn
#include "tracer_reactive_species.inc"
#undef TRACER_REACTIVE_SPECIES
#endif

   IMPLICIT NONE

#ifdef BGC
#define TRACER_REACTIVE_SPECIES(module_name, register_fn) CALL register_fn ()
#include "tracer_reactive_species.inc"
#undef TRACER_REACTIVE_SPECIES
#endif

END SUBROUTINE register_all_reactive_callbacks
#endif
