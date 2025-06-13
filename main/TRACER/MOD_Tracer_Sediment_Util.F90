#include <define.h>

MODULE MOD_Tracer_Sediment_Util

USE MOD_Precision, only: r8
USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, DEF_Tracer_Number, tracer_info_type
! We will get CaMa's nsed and sDiam (numeric particle diameters) from yos_cmf_sed
USE yos_cmf_sed, only: nsed_cama_actual => nsed, sDiam_cama_values => sDiam 
! Renaming for clarity to avoid conflict if this module also had nsed/sDiam
USE MOD_SPMD_Task, only: p_is_master ! For conditional WRITE statements

IMPLICIT NONE
SAVE

PUBLIC :: colm_to_cama_sed_map
PUBLIC :: Initialize_Tracer_To_CaMa_Mapping

! Stores the CaMa-Flood sediment class index (1-based) for each CoLM tracer.
! Value is -1 if the CoLM tracer is not mapped (e.g., not 'Suspended' or no match found).
integer, allocatable, public :: colm_to_cama_sed_map(:)

! Characteristic diameters for CoLM's 'sand', 'silt', 'clay' tracers (in meters)
! These are defined here as parameters for now. Could be made configurable.
real(r8), parameter :: DIAM_SAND_COLM = 0.50e-3  ! 0.50 mm
real(r8), parameter :: DIAM_SILT_COLM = 0.02e-3  ! 0.02 mm
real(r8), parameter :: DIAM_CLAY_COLM = 0.001e-3 ! 0.001 mm

CONTAINS

SUBROUTINE Initialize_Tracer_To_CaMa_Mapping()
  integer :: i_colm_tracer, j_cama_class
  real(r8) :: colm_d, cama_d
  real(r8) :: min_diff, current_diff
  integer :: best_match_idx

  IF (.NOT. allocated(DEF_Tracers) .OR. DEF_Tracer_Number == 0) THEN
     IF (p_is_master) WRITE(*,*) "Initialize_Tracer_To_CaMa_Mapping: CoLM Tracers not defined. Skipping mapping."
     RETURN
  ENDIF

  IF (.NOT. allocated(sDiam_cama_values) .OR. nsed_cama_actual == 0) THEN
     IF (p_is_master) THEN
        WRITE(*,*) "Initialize_Tracer_To_CaMa_Mapping: CaMa sediment diameters (sDiam_cama_values) not available or nsed_cama_actual is 0. Skipping mapping."
        WRITE(*,*) "Ensure CaMa-Flood sediment module (yos_cmf_sed) is initialized first and data is available."
     ENDIF
     ! Ensure map is allocated even if we can't populate it, to avoid issues if accessed later
     IF (.NOT. allocated(colm_to_cama_sed_map)) THEN
        ALLOCATE(colm_to_cama_sed_map(DEF_Tracer_Number))
     ENDIF
     colm_to_cama_sed_map = -1 
     RETURN
  ENDIF

  IF (allocated(colm_to_cama_sed_map)) DEALLOCATE(colm_to_cama_sed_map)
  ALLOCATE(colm_to_cama_sed_map(DEF_Tracer_Number))
  colm_to_cama_sed_map = -1 ! Default to not mapped

  IF (p_is_master) THEN
     WRITE(*,*) "Initialize_Tracer_To_CaMa_Mapping: Attempting to map CoLM tracers to CaMa sediment classes..."
     WRITE(*,*) "CaMa nsed_cama_actual: ", nsed_cama_actual
     DO j_cama_class = 1, nsed_cama_actual
         WRITE(*,*) "CaMa Class ", j_cama_class, " Diameter (sDiam_cama_values): ", sDiam_cama_values(j_cama_class), " m"
     ENDDO
  ENDIF

  DO i_colm_tracer = 1, DEF_Tracer_Number
     IF (TRIM(ADJUSTL(DEF_Tracers(i_colm_tracer)%type)) == 'Suspended' .OR. &
         TRIM(ADJUSTL(DEF_Tracers(i_colm_tracer)%type)) == 'suspended') THEN ! Case-insensitive check
        IF (TRIM(ADJUSTL(DEF_Tracers(i_colm_tracer)%name)) == 'sand') THEN
           colm_d = DIAM_SAND_COLM
        ELSE IF (TRIM(ADJUSTL(DEF_Tracers(i_colm_tracer)%name)) == 'silt') THEN
           colm_d = DIAM_SILT_COLM
        ELSE IF (TRIM(ADJUSTL(DEF_Tracers(i_colm_tracer)%name)) == 'clay') THEN
           colm_d = DIAM_CLAY_COLM
        ELSE
           IF (p_is_master) WRITE(*,*) "Warning: Suspended tracer '", TRIM(DEF_Tracers(i_colm_tracer)%name), "' has no defined diameter. Cannot map to CaMa class."
           CYCLE ! Next CoLM tracer
        ENDIF

        min_diff = 1.0e30_r8
        best_match_idx = -1
        DO j_cama_class = 1, nsed_cama_actual
           cama_d = sDiam_cama_values(j_cama_class) ! CaMa sDiam should be in meters
           current_diff = ABS(LOG(colm_d) - LOG(cama_d)) ! Compare logs of diameters for relative closeness
           IF (current_diff < min_diff) THEN
              min_diff = current_diff
              best_match_idx = j_cama_class
           ENDIF
        END DO

        IF (best_match_idx /= -1) THEN
           ! Heuristic: if diameters are within a factor of 5 (log diff approx 1.6)
           IF (min_diff < LOG(5.0_r8)) THEN
              colm_to_cama_sed_map(i_colm_tracer) = best_match_idx
              IF (p_is_master) WRITE(*,*) "Mapped CoLM tracer '", TRIM(DEF_Tracers(i_colm_tracer)%name), "' (D=", colm_d, "m)", &
                         " to CaMa class ", best_match_idx, " (D=", sDiam_cama_values(best_match_idx), "m)"
           ELSE
              IF (p_is_master) WRITE(*,*) "Warning: No close diameter match for CoLM tracer '", TRIM(DEF_Tracers(i_colm_tracer)%name), &
                         "' (D=", colm_d, "m). Closest CaMa class ", best_match_idx, " (D=", sDiam_cama_values(best_match_idx), "m) is too different."
           ENDIF
        ELSE
           IF (p_is_master) WRITE(*,*) "Warning: No CaMa sediment class found to map CoLM tracer '", TRIM(DEF_Tracers(i_colm_tracer)%name), "'"
        ENDIF
     ENDIF 
  END DO

END SUBROUTINE Initialize_Tracer_To_CaMa_Mapping

END MODULE MOD_Tracer_Sediment_Util
