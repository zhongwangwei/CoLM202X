#include <define.h>

#ifdef TRACER
! Generic particle registration stub for mkinidata.x.
! mkinidata can link generic river/TRACER dispatch objects through OBJS_BASIC,
! but it does not run particle tracer lifecycles. Keep species-private sediment
! implementations out of the mkinidata target and satisfy the generic dispatch
! layer with this empty registrar.
SUBROUTINE register_all_particle_callbacks ()
   IMPLICIT NONE
END SUBROUTINE register_all_particle_callbacks
#endif
