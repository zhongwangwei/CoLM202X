#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_Preprocessing

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Tracer_Defs, only: tracer_defs_init, tracer_param_file_for_index
   USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE, read_methane_namelist
   USE MOD_Tracer_Reactive_Methane_Registry, only: igas_ch4, methane_is_active, &
      methane_registry_refresh

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: methane_preprocessing_requirements

CONTAINS

   SUBROUTINE methane_preprocessing_requirements (requires_lake_soilc, requires_spatial_ph)
      IMPLICIT NONE
      logical, intent(out) :: requires_lake_soilc, requires_spatial_ph

      character(len=512) :: file_param
      logical :: found

      requires_lake_soilc = .false.
      requires_spatial_ph = .false.

      CALL tracer_defs_init ()
      CALL methane_registry_refresh ()
      IF (.not. methane_is_active()) RETURN

      ! Use the runtime parser so keyed aliases, positional files, separators,
      ! and null handling cannot drift between mksrfdata and colm.x.
      CALL tracer_param_file_for_index (igas_ch4, 'METHANE', file_param, found)
      IF (.not. found) THEN
         CALL CoLM_stop (' ***** ERROR: CH4 requires DEF_TRACER_PARAM_FILES to include a CH4 parameter file')
      ENDIF
      CALL read_methane_namelist (file_param)

      requires_lake_soilc = DEF_METHANE%allowlakeprod
      requires_spatial_ph = DEF_METHANE%use_spatial_ph

   END SUBROUTINE methane_preprocessing_requirements

END MODULE MOD_Tracer_Reactive_Methane_Preprocessing
#endif
