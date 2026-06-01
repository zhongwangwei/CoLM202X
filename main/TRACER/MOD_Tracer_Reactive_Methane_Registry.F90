#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_Registry
!=======================================================================
! Methane species lookup against the TRACER registry.
!
! Methane is a reactive tracer in the TRACER framework. Its activation
! is controlled at runtime by whether a tracer with name "CH4" (or
! "METHANE", case-insensitive) is registered as category="reactive" in
! the namelist DEF_TRACER_NAMES / DEF_TRACER_TYPES. No #define switch
! is needed at the macro level: when CH4 is not in the namelist, the
! igas_ch4 index stays at the sentinel value (-1), and all gated
! methane logic (driver / accumulators / history) becomes a no-op.
!
! Additional oxidant species O2 and CO2 follow the same lookup so
! per-process oxidation/respiration channels can reference their
! tracer indices when present; both are optional and default to -1.
!=======================================================================

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, tracer_is_reactive, tracer_upper

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Sentinel value: tracer absent (or non-reactive) in the registry.
   integer, parameter :: METHANE_GAS_ABSENT = -1

   integer :: igas_ch4 = METHANE_GAS_ABSENT
   integer :: igas_o2  = METHANE_GAS_ABSENT
   integer :: igas_co2 = METHANE_GAS_ABSENT

   PUBLIC :: methane_registry_init
   PUBLIC :: methane_registry_refresh
   PUBLIC :: methane_is_active
   PUBLIC :: igas_ch4, igas_o2, igas_co2, METHANE_GAS_ABSENT

CONTAINS

   !-------------------------------------------------------------------
   ! methane_registry_init — populate igas_* from the tracer registry and
   ! report the resulting CH4 activation state.
   ! Must be called AFTER tracer_defs_init() because it reads
   ! tracers(:)%name and %category.
   !-------------------------------------------------------------------
   SUBROUTINE methane_registry_init ()
      IMPLICIT NONE

      CALL update_methane_registry (verbose=.true.)

   END SUBROUTINE methane_registry_init

   !-------------------------------------------------------------------
   ! methane_registry_refresh — refresh cached tracer indices without
   ! producing init-time status messages.  This is the callback used by
   ! the generic reactive dispatcher before active-species gating.
   !-------------------------------------------------------------------
   SUBROUTINE methane_registry_refresh ()
      IMPLICIT NONE

      CALL update_methane_registry (verbose=.false.)

   END SUBROUTINE methane_registry_refresh

   !-------------------------------------------------------------------
   SUBROUTINE update_methane_registry (verbose)
      USE MOD_SPMD_Task, only: p_is_master
      IMPLICIT NONE
      logical, intent(in) :: verbose
      integer :: i
      character(len=:), allocatable :: nm

      igas_ch4 = METHANE_GAS_ABSENT
      igas_o2  = METHANE_GAS_ABSENT
      igas_co2 = METHANE_GAS_ABSENT

      IF (.not. allocated(tracers)) RETURN
      IF (ntracers <= 0) RETURN

      DO i = 1, ntracers
         nm = tracer_upper(tracers(i)%name)
         IF (trim(nm) == 'CH4' .or. trim(nm) == 'METHANE') THEN
            IF (tracer_is_reactive(i)) THEN
               igas_ch4 = i
            ELSEIF (verbose .and. p_is_master) THEN
               write(*,'(A,A,A)') ' WARNING methane_registry_init: tracer "', &
                  trim(tracers(i)%name), '" found but not category="reactive"; methane disabled.'
            ENDIF
         ELSEIF (trim(nm) == 'O2') THEN
            IF (tracer_is_reactive(i)) THEN
               igas_o2 = i
            ELSEIF (verbose .and. p_is_master) THEN
               write(*,'(A,A,A)') ' WARNING methane_registry_init: tracer "', &
                  trim(tracers(i)%name), '" found but not category="reactive"; O2 channel disabled.'
            ENDIF
         ELSEIF (trim(nm) == 'CO2') THEN
            IF (tracer_is_reactive(i)) THEN
               igas_co2 = i
            ELSEIF (verbose .and. p_is_master) THEN
               write(*,'(A,A,A)') ' WARNING methane_registry_init: tracer "', &
                  trim(tracers(i)%name), '" found but not category="reactive"; CO2 channel disabled.'
            ENDIF
         ENDIF
      ENDDO

      IF (verbose .and. p_is_master) THEN
         IF (igas_ch4 > 0) THEN
            write(*,'(A,I0,A)') ' methane_registry_init: CH4 active at tracer index ', igas_ch4, '.'
         ELSE
            write(*,'(A)') ' methane_registry_init: CH4 tracer not registered; methane driver disabled.'
         ENDIF
      ENDIF

   END SUBROUTINE update_methane_registry

   !-------------------------------------------------------------------
   logical FUNCTION methane_is_active ()
      methane_is_active = (igas_ch4 > 0)
   END FUNCTION methane_is_active


END MODULE MOD_Tracer_Reactive_Methane_Registry
#endif
