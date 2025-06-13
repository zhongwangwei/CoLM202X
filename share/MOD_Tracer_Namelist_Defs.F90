#include <define.h>

MODULE MOD_Tracer_Namelist_Defs

USE MOD_Precision, only: r8
USE MOD_DataType, only: tracer_info_type

IMPLICIT NONE
SAVE

PUBLIC :: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL
PUBLIC :: allocate_tracer_defs, initialize_tracer_forcing_nl_defaults, parse_tracer_names
PUBLIC :: MAX_TRACER_FORCING_VARS

! Max number of variables per tracer forcing (e.g., value for precip, value for humidity)
integer, parameter :: MAX_TRACER_FORCING_VARS = 10

! Array to hold all configured tracers' metadata
type(tracer_info_type), allocatable :: DEF_Tracers(:)

! Derived type for a single tracer's forcing namelist configuration
type :: nl_tracer_forcing_type
  character(len=64)  :: tracer_name        = 'null'
  integer            :: NVAR               = 0
  integer            :: dtime(MAX_TRACER_FORCING_VARS)    = 0
  integer            :: offset(MAX_TRACER_FORCING_VARS)   = 0
  character(len=256) :: fprefix(MAX_TRACER_FORCING_VARS)  = 'null'
  character(len=256) :: vname(MAX_TRACER_FORCING_VARS)    = 'null'
  character(len=256) :: timelog(MAX_TRACER_FORCING_VARS)  = 'instant'
  character(len=256) :: tintalgo(MAX_TRACER_FORCING_VARS) = 'linear'
END type nl_tracer_forcing_type

! Array to hold the forcing namelist configurations for tracers that need them
type(nl_tracer_forcing_type), allocatable :: DEF_Tracer_Forcings_NL(:)

CONTAINS

SUBROUTINE allocate_tracer_defs(num_tracers, num_forced_tracers)
  integer, intent(in) :: num_tracers
  integer, intent(in) :: num_forced_tracers

  IF (num_tracers > 0) THEN
    IF (allocated(DEF_Tracers)) DEALLOCATE(DEF_Tracers)
    ALLOCATE(DEF_Tracers(num_tracers))
  ELSE IF (allocated(DEF_Tracers)) THEN ! num_tracers is 0 or less
    DEALLOCATE(DEF_Tracers)
  ENDIF

  IF (num_forced_tracers > 0) THEN
    IF (allocated(DEF_Tracer_Forcings_NL)) DEALLOCATE(DEF_Tracer_Forcings_NL)
    ALLOCATE(DEF_Tracer_Forcings_NL(num_forced_tracers))
  ELSE IF (allocated(DEF_Tracer_Forcings_NL)) THEN ! num_forced_tracers is 0 or less
    DEALLOCATE(DEF_Tracer_Forcings_NL)
  ENDIF
END SUBROUTINE allocate_tracer_defs

SUBROUTINE initialize_tracer_forcing_nl_defaults(tracer_forcing_nl_item)
  type(nl_tracer_forcing_type), intent(out) :: tracer_forcing_nl_item
  
  tracer_forcing_nl_item%tracer_name = 'null'
  tracer_forcing_nl_item%NVAR = 0 
  tracer_forcing_nl_item%dtime(:) = 0
  tracer_forcing_nl_item%offset(:) = 0
  tracer_forcing_nl_item%fprefix(:) = 'null'
  tracer_forcing_nl_item%vname(:) = 'null'
  tracer_forcing_nl_item%timelog(:) = 'instant'
  tracer_forcing_nl_item%tintalgo(:) = 'linear'
END SUBROUTINE initialize_tracer_forcing_nl_defaults

SUBROUTINE parse_tracer_names(tracer_name_str, tracer_type_str, num_tracers_expected, tracers_array_out)
  character(len=*), intent(in) :: tracer_name_str
  character(len=*), intent(in) :: tracer_type_str ! New argument
  integer, intent(in) :: num_tracers_expected
  type(tracer_info_type), intent(out) :: tracers_array_out(:)

  character(len=256) :: temp_name_str, temp_type_str
  integer :: i, current_name_pos, next_name_comma, len_name
  integer :: current_type_pos, next_type_comma, len_type
  integer :: actual_tracers_found, actual_types_found
  
  IF(num_tracers_expected <= 0) RETURN
  IF(SIZE(tracers_array_out) < num_tracers_expected) THEN
    WRITE(*,*) 'Error in parse_tracer_names: output array too small.'
    RETURN
  ENDIF

  temp_name_str = TRIM(tracer_name_str)
  temp_type_str = TRIM(tracer_type_str)
  current_name_pos = 1
  current_type_pos = 1
  actual_tracers_found = 0
  actual_types_found = 0

  DO i = 1, num_tracers_expected
      ! Parse Name
      IF (current_name_pos > LEN_TRIM(temp_name_str)) THEN
          actual_tracers_found = i - 1
          ! If names run out, types should also ideally, or it's a mismatch
          IF (current_type_pos <= LEN_TRIM(temp_type_str)) THEN
               WRITE(*,*) 'Warning: More tracer types provided than tracer names.'
          ENDIF
          tracers_array_out(i)%name = 'NAME_MISSING' ! Indicate name was expected but not found
          tracers_array_out(i)%type = 'TYPE_NOT_PARSED'
          EXIT
      ENDIF
      next_name_comma = INDEX(temp_name_str(current_name_pos:), ',')
      IF (next_name_comma == 0) THEN 
          len_name = LEN_TRIM(temp_name_str(current_name_pos:))
          IF (len_name > 0) THEN
              tracers_array_out(i)%name = TRIM(ADJUSTL(temp_name_str(current_name_pos:current_name_pos+len_name-1)))
              actual_tracers_found = i
          ELSE
              tracers_array_out(i)%name = 'NAME_PARSE_ERROR' 
              actual_tracers_found = i ! Still counts as an attempt to parse
          ENDIF
          ! For name parsing, if it's the last one, we will process type then exit loop.
      ELSE
          len_name = next_name_comma - 1
          IF (len_name > 0) THEN
              tracers_array_out(i)%name = TRIM(ADJUSTL(temp_name_str(current_name_pos:current_name_pos+len_name-1)))
              actual_tracers_found = i
          ELSE
              tracers_array_out(i)%name = 'EMPTY_TRACER_NAME'
              actual_tracers_found = i ! Still counts as an attempt to parse
          ENDIF
          current_name_pos = current_name_pos + len_name + 1 
      ENDIF

      ! Parse Type
      IF (current_type_pos > LEN_TRIM(temp_type_str)) THEN
          tracers_array_out(i)%type = 'TYPE_MISSING' ! Default if types run out
          IF (actual_tracers_found == i .AND. tracers_array_out(i)%name /= 'NAME_PARSE_ERROR' .AND. tracers_array_out(i)%name /= 'EMPTY_TRACER_NAME') THEN 
              WRITE(*,*) 'Warning: Fewer tracer types provided than tracer names. Tracer: ', TRIM(tracers_array_out(i)%name)
          ENDIF
          IF (next_name_comma == 0) EXIT ! if names also ended
          CYCLE
      ENDIF
      next_type_comma = INDEX(temp_type_str(current_type_pos:), ',')
      IF (next_type_comma == 0) THEN
          len_type = LEN_TRIM(temp_type_str(current_type_pos:))
          IF (len_type > 0 .AND. len_type <= LEN(tracers_array_out(i)%type)) THEN
              tracers_array_out(i)%type = TRIM(ADJUSTL(temp_type_str(current_type_pos:current_type_pos+len_type-1)))
              actual_types_found = i
          ELSE IF (len_type > LEN(tracers_array_out(i)%type)) THEN
              tracers_array_out(i)%type = temp_type_str(current_type_pos:current_type_pos+LEN(tracers_array_out(i)%type)-1)
              WRITE(*,*) 'Warning: Tracer type string too long for: ', TRIM(tracers_array_out(i)%name), '. Truncated.'
              actual_types_found = i
          ELSE
              tracers_array_out(i)%type = 'TYPE_PARSE_ERROR'
              actual_types_found = i ! Still counts as an attempt to parse
          ENDIF
          IF (next_name_comma == 0) EXIT ! Both names and types are on their last item
      ELSE
          len_type = next_type_comma - 1
          IF (len_type > 0 .AND. len_type <= LEN(tracers_array_out(i)%type)) THEN
              tracers_array_out(i)%type = TRIM(ADJUSTL(temp_type_str(current_type_pos:current_type_pos+len_type-1)))
              actual_types_found = i
          ELSE IF (len_type > LEN(tracers_array_out(i)%type)) THEN
              tracers_array_out(i)%type = temp_type_str(current_type_pos:current_type_pos+LEN(tracers_array_out(i)%type)-1)
              WRITE(*,*) 'Warning: Tracer type string too long for: ', TRIM(tracers_array_out(i)%name), '. Truncated.'
              actual_types_found = i
          ELSE
              tracers_array_out(i)%type = 'EMPTY_TRACER_TYPE'
              actual_types_found = i ! Still counts as an attempt to parse
          ENDIF
          current_type_pos = current_type_pos + len_type + 1
      ENDIF
      
      IF (next_name_comma == 0) EXIT ! if names ended, exit after processing current type
  END DO
  
  IF (actual_tracers_found < num_tracers_expected .OR. actual_types_found < num_tracers_expected) THEN
      WRITE(*,*) 'Warning: Mismatch in parsed tracer names/types or items provided.'
      WRITE(*,*) 'Expected: ', num_tracers_expected, ' tracer items.'
      WRITE(*,*) 'Found: ', actual_tracers_found, ' names and ', actual_types_found, ' types from strings.'
      IF (actual_tracers_found < num_tracers_expected) THEN
          DO i = actual_tracers_found + 1, num_tracers_expected
             tracers_array_out(i)%name = 'UNDEFINED_NAME'
             tracers_array_out(i)%type = 'UNDEFINED_TYPE'
          ENDDO
      ELSE IF (actual_types_found < num_tracers_expected) THEN ! Names might be more, but types less
          DO i = actual_types_found + 1, num_tracers_expected
             IF (LEN_TRIM(tracers_array_out(i)%name) > 0 .AND. tracers_array_out(i)%name(1:4) /= 'NAME' &
                 .AND. tracers_array_out(i)%name(1:5) /= 'EMPTY' .AND. tracers_array_out(i)%name(1:9) /= 'UNDEFINED') THEN
                 tracers_array_out(i)%type = 'TYPE_MISSING'
             ELSE
                 tracers_array_out(i)%type = 'UNDEFINED_TYPE' ! If name was also problematic
             ENDIF
          ENDDO
      ENDIF
  ENDIF

END SUBROUTINE parse_tracer_names

END MODULE MOD_Tracer_Namelist_Defs
