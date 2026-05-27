#include <define.h>

#ifdef TRACER
! Generic reactive registration stub for mkinidata.x.
! mkinidata links generic TRACER restart/time-variable objects but does not
! execute land reactive tracer lifecycles.  Keep this species-agnostic so
! adding a new reactive species does not require another mkinidata stub file.
SUBROUTINE register_all_reactive_callbacks ()
   IMPLICIT NONE
END SUBROUTINE register_all_reactive_callbacks
#endif
