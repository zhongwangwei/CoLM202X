#include <define.h>

MODULE MOD_Tracer_Namelist_Defs

USE MOD_Precision, only: r8

IMPLICIT NONE
SAVE

PUBLIC :: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL
PUBLIC :: allocate_tracer_defs, initialize_tracer_forcing_nl_defaults, parse_tracer_names
PUBLIC :: MAX_TRACER_FORCING_VARS
PUBLIC :: initialize_tracer_namelists, broadcast_tracer_namelists
PUBLIC :: DEF_Tracer_Number, DEF_USE_Tracer, DEF_Tracer_Init, DEF_Tracer_Init_file
PUBLIC :: debug_tracer_forcing_configurations

! Max number of variables per tracer forcing (e.g., value for precip, value for humidity)
integer, parameter :: MAX_TRACER_FORCING_VARS = 10
!------- Tracer information type -------
type :: tracer_info_type
    character(len=64) :: name         ! Name of the tracer (e.g., 'O18', 'sand')
    character(len=16) :: type         ! Type: 'dissolved' or 'suspended'
!    Add other relevant metadata if needed in the future
END type tracer_info_type

! Array to hold all configured tracers' metadata
type(tracer_info_type), allocatable :: DEF_Tracers(:)

! Derived type for a single tracer's forcing namelist configuration
type :: nl_tracer_forcing_type
  character(len=64)  :: tracer_name                       = 'null'
  character(len=256) :: dataset_name                      = 'null'
  character(len=16)  :: tracer_type                       = 'null'
  character(len=256) :: tracer_dir                        = 'null'
  integer            :: startyr                           = 0
  integer            :: startmo                           = 0
  integer            :: endyr                             = 0
  integer            :: endmo                             = 0
  integer            :: dtime(MAX_TRACER_FORCING_VARS)    = 0
  integer            :: offset(MAX_TRACER_FORCING_VARS)   = 0

  logical            :: leapyear                          = .false.
  logical            :: data2d                            = .false.
  logical            :: hightdim                          = .false.
  logical            :: dim2d                             = .false.
  character(len=256) :: latname                           = 'null'
  character(len=256) :: lonname                           = 'null'
  character(len=256) :: groupby                           = 'null'
  integer            :: NVAR                              = 0
  logical            :: has_missing_value                 = .true.
  character(len=256) :: missing_value_name                = 'missing_value'

  logical            :: regional                          = .false.
  real(r8)           :: regbnd(4)                         = (/-90.0, 90.0, -180.0, 180.0/)

  character(len=256) :: fprefix(MAX_TRACER_FORCING_VARS)  = 'null'
  character(len=256) :: vname(MAX_TRACER_FORCING_VARS)    = 'null'
  character(len=256) :: timelog(MAX_TRACER_FORCING_VARS)  = 'instant'
  character(len=256) :: tintalgo(MAX_TRACER_FORCING_VARS) = 'linear'
END type nl_tracer_forcing_type



! Array to hold the forcing namelist configurations for tracers that need them
  type(nl_tracer_forcing_type), allocatable :: DEF_Tracer_Forcings_NL(:)
  
  ! Variables that should be defined here but are currently in MOD_Namelist
  ! These are added for compatibility with other modules that use this module
  integer :: DEF_Tracer_Number = 0
  logical :: DEF_USE_Tracer = .false.
  logical :: DEF_Tracer_Init = .false.
  character(len=256) :: DEF_Tracer_Init_file = 'null'
  
  ! Module-level variables for reading tracer forcing namelist
  character(len=64)  :: nl_tracer_name        = 'null'
  character(len=256) :: nl_dataset_name       = 'null'
  character(len=16)  :: nl_tracer_type        = 'null'
  character(len=256) :: nl_tracer_dir         = 'null'
  integer            :: nl_startyr            = 0
  integer            :: nl_startmo            = 0
  integer            :: nl_endyr              = 0
  integer            :: nl_endmo              = 0
  logical            :: nl_leapyear          = .false.
  logical            :: nl_data2d            = .false.
  logical            :: nl_hightdim          = .false.
  logical            :: nl_dim2d             = .false.
  character(len=256) :: nl_latname           = 'null'
  character(len=256) :: nl_lonname           = 'null'
  character(len=256) :: nl_groupby           = 'null'
  integer            :: nl_NVAR               = 0
  logical            :: nl_has_missing_value  = .false.
  character(len=256) :: nl_missing_value_name = 'missing_value'
  integer            :: nl_dtime(MAX_TRACER_FORCING_VARS)    = 0
  integer            :: nl_offset(MAX_TRACER_FORCING_VARS)   = 0
  logical            :: nl_regional                         = .false.
  real(r8)           :: nl_regbnd(4)                        = (/-90.0, 90.0, -180.0, 180.0/)
  character(len=256) :: nl_fprefix(MAX_TRACER_FORCING_VARS)  = 'null'
  character(len=256) :: nl_vname(MAX_TRACER_FORCING_VARS)    = 'null'
  character(len=256) :: nl_timelog(MAX_TRACER_FORCING_VARS)  = 'instant'
  character(len=256) :: nl_tintalgo(MAX_TRACER_FORCING_VARS) = 'linear'
  
  ! Namelist declaration for reading tracer forcing configuration
  namelist /nl_colm_tracer_forcing/ nl_tracer_name, nl_dataset_name, nl_tracer_type, nl_tracer_dir,  &
   nl_startyr, nl_startmo, nl_endyr, nl_endmo, nl_leapyear, nl_data2d, nl_hightdim, nl_dim2d, nl_latname, nl_lonname, nl_groupby, &
   nl_NVAR, nl_has_missing_value, nl_missing_value_name, nl_dtime, nl_offset, nl_regional, nl_regbnd, nl_fprefix, nl_vname, nl_timelog, nl_tintalgo
  
  CONTAINS

SUBROUTINE allocate_tracer_defs(num_tracers, num_forced_tracers)
 ! num_tracers is the total number of tracers, num_forced_tracers is the number of tracers that need forcing
  integer, intent(in) :: num_tracers
  integer, intent(in) :: num_forced_tracers

  ! Always deallocate first to prevent memory leaks
  IF (allocated(DEF_Tracers)) DEALLOCATE(DEF_Tracers)
  IF (allocated(DEF_Tracer_Forcings_NL)) DEALLOCATE(DEF_Tracer_Forcings_NL)

  IF (num_tracers > 0) THEN
    ALLOCATE(DEF_Tracers(num_tracers))
  ENDIF

  IF (num_forced_tracers > 0) THEN
    ALLOCATE(DEF_Tracer_Forcings_NL(num_forced_tracers))
  ENDIF
END SUBROUTINE allocate_tracer_defs

SUBROUTINE initialize_tracer_forcing_nl_defaults(tracer_forcing_nl_item)
  type(nl_tracer_forcing_type), intent(out) :: tracer_forcing_nl_item
  
  tracer_forcing_nl_item%tracer_name = 'null'
  tracer_forcing_nl_item%dataset_name = 'null'
  tracer_forcing_nl_item%tracer_type = 'null'
  tracer_forcing_nl_item%tracer_dir = 'null'
  tracer_forcing_nl_item%startyr = 0
  tracer_forcing_nl_item%startmo = 0
  tracer_forcing_nl_item%endyr = 0
  tracer_forcing_nl_item%endmo = 0
  tracer_forcing_nl_item%leapyear = .false.
  tracer_forcing_nl_item%data2d = .false.
  tracer_forcing_nl_item%hightdim = .false.
  tracer_forcing_nl_item%dim2d = .false.
  tracer_forcing_nl_item%latname = 'null'
  tracer_forcing_nl_item%lonname = 'null'
  tracer_forcing_nl_item%groupby = 'null'
  tracer_forcing_nl_item%NVAR = 0 
  tracer_forcing_nl_item%has_missing_value = .true.
  tracer_forcing_nl_item%missing_value_name = 'missing_value'
  tracer_forcing_nl_item%dtime(:) = 0
  tracer_forcing_nl_item%offset(:) = 0
  tracer_forcing_nl_item%regional = .false.
  tracer_forcing_nl_item%regbnd(:) = (/-90.0, 90.0, -180.0, 180.0/)
  tracer_forcing_nl_item%fprefix(:) = 'null'
  tracer_forcing_nl_item%vname(:) = 'null'
  tracer_forcing_nl_item%timelog(:) = 'instant'
  tracer_forcing_nl_item%tintalgo(:) = 'linear'
END SUBROUTINE initialize_tracer_forcing_nl_defaults

SUBROUTINE parse_tracer_names(tracer_name_str, tracer_type_str, num_tracers_expected, tracers_array_out)
  ! num_tracers_expected is the number of tracers expected to be parsed
  ! tracers_array_out is the array to hold the parsed tracer names and types
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
          DO i = actual_tracers_found + 1, min(num_tracers_expected, SIZE(tracers_array_out))
             tracers_array_out(i)%name = 'UNDEFINED_NAME'
             tracers_array_out(i)%type = 'UNDEFINED_TYPE'
          ENDDO
      ELSE IF (actual_types_found < num_tracers_expected) THEN ! Names might be more, but types less
          DO i = actual_types_found + 1, min(num_tracers_expected, SIZE(tracers_array_out))
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

SUBROUTINE initialize_tracer_namelists(tracer_forcing_namelist,num_tracers)
  IMPLICIT NONE
  
  character(len=256), intent(in) :: tracer_forcing_namelist
  integer, intent(in) :: num_tracers
  integer :: i
  integer :: ierr_tracer_nl_open, ierr_tracer_nl_read
  integer :: unit_number

  IF (.NOT. allocated(DEF_Tracer_Forcings_NL) .OR. .NOT. allocated(DEF_Tracers)) RETURN

  ! Initialize all tracers with defaults first
  DO i = 1, SIZE(DEF_Tracer_Forcings_NL)
      CALL initialize_tracer_forcing_nl_defaults(DEF_Tracer_Forcings_NL(i))
      IF (i <= SIZE(DEF_Tracers)) THEN
         DEF_Tracer_Forcings_NL(i)%tracer_name = DEF_Tracers(i)%name
      ENDIF
  END DO

  WRITE(*,*) 'Reading tracer forcing configurations from: ', trim(tracer_forcing_namelist)
  
  ! Get a free unit number
  unit_number = 20

  ! Open the tracer forcing namelist file
  OPEN(unit_number, FILE=trim(tracer_forcing_namelist), STATUS='OLD', FORM='FORMATTED', IOSTAT=ierr_tracer_nl_open)
   
  IF (ierr_tracer_nl_open == 0) THEN
      ! Read tracer forcing configurations for each tracer
      DO i = 1, SIZE(DEF_Tracer_Forcings_NL)
          IF (i > SIZE(DEF_Tracers)) EXIT
          
          ! Initialize module-level variables with current defaults for this tracer
          nl_tracer_name = DEF_Tracer_Forcings_NL(i)%tracer_name
          nl_dataset_name = DEF_Tracer_Forcings_NL(i)%dataset_name
          nl_tracer_type = DEF_Tracer_Forcings_NL(i)%tracer_type
          nl_tracer_dir = DEF_Tracer_Forcings_NL(i)%tracer_dir
          nl_startyr = DEF_Tracer_Forcings_NL(i)%startyr
          nl_startmo = DEF_Tracer_Forcings_NL(i)%startmo
          nl_endyr = DEF_Tracer_Forcings_NL(i)%endyr
          nl_endmo = DEF_Tracer_Forcings_NL(i)%endmo
          nl_leapyear = DEF_Tracer_Forcings_NL(i)%leapyear
          nl_data2d = DEF_Tracer_Forcings_NL(i)%data2d
          nl_hightdim = DEF_Tracer_Forcings_NL(i)%hightdim
          nl_dim2d = DEF_Tracer_Forcings_NL(i)%dim2d
          nl_latname = DEF_Tracer_Forcings_NL(i)%latname
          nl_lonname = DEF_Tracer_Forcings_NL(i)%lonname
          nl_groupby = DEF_Tracer_Forcings_NL(i)%groupby
          nl_NVAR = DEF_Tracer_Forcings_NL(i)%NVAR
          nl_has_missing_value = DEF_Tracer_Forcings_NL(i)%has_missing_value
          nl_missing_value_name = DEF_Tracer_Forcings_NL(i)%missing_value_name
          nl_dtime = DEF_Tracer_Forcings_NL(i)%dtime
          nl_offset = DEF_Tracer_Forcings_NL(i)%offset
          nl_regional = DEF_Tracer_Forcings_NL(i)%regional
          nl_regbnd = DEF_Tracer_Forcings_NL(i)%regbnd
          nl_fprefix = DEF_Tracer_Forcings_NL(i)%fprefix
          nl_vname = DEF_Tracer_Forcings_NL(i)%vname
          nl_timelog = DEF_Tracer_Forcings_NL(i)%timelog
          nl_tintalgo = DEF_Tracer_Forcings_NL(i)%tintalgo
          
          ! Rewind to start of file for each tracer
          REWIND(unit_number, IOSTAT=ierr_tracer_nl_read)
          IF (ierr_tracer_nl_read /= 0) THEN
              WRITE(*,*) 'Error rewinding file for tracer ', i, ': ', ierr_tracer_nl_read
              CYCLE
          ENDIF
          
          ! Try to read a namelist that matches this tracer
          ! The namelist file should contain sections like:
          ! &nl_colm_tracer_forcing
          !   nl_tracer_name = 'O18'
          !   nl_NVAR = 2
          !   ...
          ! /
          DO WHILE (.TRUE.)
              READ(unit_number, NML=nl_colm_tracer_forcing, IOSTAT=ierr_tracer_nl_read)
              
              IF (ierr_tracer_nl_read /= 0) THEN
                  ! End of file or no more namelists found
                  EXIT
              ENDIF
              
              ! Check if this namelist is for the current tracer
              IF (TRIM(ADJUSTL(nl_tracer_name)) == TRIM(ADJUSTL(DEF_Tracers(i)%name)) .OR. &
                  TRIM(ADJUSTL(nl_tracer_name)) == 'null') THEN
                  
                  ! Copy the read values back to the derived type for this tracer
                  DEF_Tracer_Forcings_NL(i)%tracer_name = DEF_Tracers(i)%name  ! Use actual tracer name
                  DEF_Tracer_Forcings_NL(i)%dataset_name = nl_dataset_name
                  DEF_Tracer_Forcings_NL(i)%tracer_type = DEF_Tracers(i)%type
                  DEF_Tracer_Forcings_NL(i)%tracer_dir = nl_tracer_dir
                  DEF_Tracer_Forcings_NL(i)%startyr = nl_startyr
                  DEF_Tracer_Forcings_NL(i)%startmo = nl_startmo
                  DEF_Tracer_Forcings_NL(i)%endyr = nl_endyr
                  DEF_Tracer_Forcings_NL(i)%endmo = nl_endmo
                  DEF_Tracer_Forcings_NL(i)%leapyear = nl_leapyear
                  DEF_Tracer_Forcings_NL(i)%data2d = nl_data2d
                  DEF_Tracer_Forcings_NL(i)%hightdim = nl_hightdim
                  DEF_Tracer_Forcings_NL(i)%dim2d = nl_dim2d
                  DEF_Tracer_Forcings_NL(i)%latname = nl_latname
                  DEF_Tracer_Forcings_NL(i)%lonname = nl_lonname
                  DEF_Tracer_Forcings_NL(i)%groupby = nl_groupby
                  DEF_Tracer_Forcings_NL(i)%NVAR = nl_NVAR
                  DEF_Tracer_Forcings_NL(i)%has_missing_value = nl_has_missing_value
                  DEF_Tracer_Forcings_NL(i)%missing_value_name = nl_missing_value_name
                  DEF_Tracer_Forcings_NL(i)%dtime = nl_dtime
                  DEF_Tracer_Forcings_NL(i)%offset = nl_offset
                  DEF_Tracer_Forcings_NL(i)%regional = nl_regional
                  DEF_Tracer_Forcings_NL(i)%regbnd = nl_regbnd
                  DEF_Tracer_Forcings_NL(i)%fprefix = nl_fprefix
                  DEF_Tracer_Forcings_NL(i)%vname = nl_vname
                  DEF_Tracer_Forcings_NL(i)%timelog = nl_timelog
                  DEF_Tracer_Forcings_NL(i)%tintalgo = nl_tintalgo
                  
                  WRITE(*,*) 'Successfully read tracer forcing config for tracer ', i, ' (', &
                             TRIM(DEF_Tracers(i)%name), ') with ', nl_NVAR, ' variables'
                  WRITE(*,*) '  Directory: ', TRIM(nl_tracer_dir)
                  EXIT  ! Found configuration for this tracer, move to next
              ENDIF
          END DO
          
          ! If no specific configuration was found, use defaults
          IF (ierr_tracer_nl_read /= 0) THEN
              WRITE(*,*) 'No specific forcing config found for tracer ', i, ' (', &
                         TRIM(DEF_Tracers(i)%name), '), using defaults'
          ENDIF
      END DO
      
      CLOSE(unit_number, IOSTAT=ierr_tracer_nl_read)
      IF (ierr_tracer_nl_read /= 0) THEN
          WRITE(*,*) 'Warning: Error closing tracer forcing namelist file: ', ierr_tracer_nl_read
      ENDIF
  ELSE
      WRITE(*,*) 'Warning: Could not open tracer forcing namelist file: ', trim(tracer_forcing_namelist)
      WRITE(*,*) 'IOSTAT error code: ', ierr_tracer_nl_open
      WRITE(*,*) 'Using default forcing configurations for all tracers'
  ENDIF

END SUBROUTINE initialize_tracer_namelists

SUBROUTINE broadcast_tracer_namelists(num_tracers)
#ifdef USEMPI
  USE MOD_SPMD_Task, only: p_is_master, p_comm_glb, p_err, p_address_master, MPI_CHARACTER, MPI_INTEGER, MPI_LOGICAL, MPI_REAL8
#endif
  IMPLICIT NONE
  
  integer, intent(in) :: num_tracers
  
#ifdef USEMPI
  integer :: num_tracers_alloc_mpi, num_forced_tracers_alloc_mpi
  integer :: i_mpi, k_mpi
  integer :: master_allocated_status, worker_allocated_status

  ! Add barrier to ensure all processes reach this point together
  CALL mpi_barrier(p_comm_glb, p_err)
  IF (p_err /= 0) THEN
      WRITE(*,*) 'Error in mpi_barrier: ', p_err
      RETURN
  ENDIF
  
  ! First, broadcast whether arrays are allocated on master
  IF (p_is_master) THEN
      IF (allocated(DEF_Tracers) .AND. allocated(DEF_Tracer_Forcings_NL)) THEN
          master_allocated_status = 1
      ELSE
          master_allocated_status = 0
      ENDIF
  ENDIF
  
  CALL mpi_bcast(master_allocated_status, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
  IF (p_err /= 0) THEN
      WRITE(*,*) 'Error in mpi_bcast for master_allocated_status: ', p_err
      RETURN
  ENDIF
  
  ! Early return if master arrays are not allocated
  IF (master_allocated_status == 0) THEN
      IF (p_is_master) WRITE(*,*) 'broadcast_tracer_namelists: Master arrays not allocated, skipping broadcast'
      RETURN
  ENDIF

  ! Ensure all worker processes allocate their arrays consistently
  num_tracers_alloc_mpi = num_tracers
  num_forced_tracers_alloc_mpi = num_tracers ! Simplified
  
  IF (.NOT. p_is_master) THEN
      ! Deallocate any existing arrays first to prevent memory issues
      IF (allocated(DEF_Tracers)) deallocate(DEF_Tracers)
      IF (allocated(DEF_Tracer_Forcings_NL)) deallocate(DEF_Tracer_Forcings_NL)
      
      ! Allocate arrays with proper sizes
      CALL allocate_tracer_defs(num_tracers_alloc_mpi, num_forced_tracers_alloc_mpi)
      
      ! Verify allocation was successful
      IF (.NOT. (allocated(DEF_Tracers) .AND. allocated(DEF_Tracer_Forcings_NL))) THEN
          WRITE(*,*) 'broadcast_tracer_namelists: Worker allocation failed'
          RETURN
      ENDIF
      
      ! Verify array sizes are correct
      IF (SIZE(DEF_Tracers) /= num_tracers_alloc_mpi .OR. SIZE(DEF_Tracer_Forcings_NL) /= num_forced_tracers_alloc_mpi) THEN
          WRITE(*,*) 'broadcast_tracer_namelists: Worker array sizes incorrect'
          RETURN
      ENDIF
  ENDIF

  ! Add another barrier to ensure all allocations are complete
  CALL mpi_barrier(p_comm_glb, p_err)

  ! Broadcast basic tracer information
  IF (allocated(DEF_Tracers)) THEN
      DO i_mpi = 1, num_tracers
          IF (i_mpi > SIZE(DEF_Tracers)) EXIT
          CALL mpi_bcast(DEF_Tracers(i_mpi)%name, 64, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracers(i_mpi)%type, 16, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
      END DO
  ENDIF
  
  ! Broadcast forcing configurations with bounds checking
  IF (allocated(DEF_Tracer_Forcings_NL)) THEN
      DO i_mpi = 1, min(num_tracers, SIZE(DEF_Tracer_Forcings_NL))
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%tracer_name, 64, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%dataset_name, 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%tracer_type, 16, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%tracer_dir, 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%has_missing_value, 1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%missing_value_name, 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%NVAR,    1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%startyr, 1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%startmo, 1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%endyr,   1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%endmo,   1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%leapyear, 1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%data2d,  1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%hightdim, 1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%dim2d,  1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%latname, 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%lonname, 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%groupby, 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%dtime,   MAX_TRACER_FORCING_VARS, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%offset,  MAX_TRACER_FORCING_VARS, MPI_INTEGER,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%regional, 1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err)
          CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%regbnd, 4, MPI_REAL8,   p_address_master, p_comm_glb, p_err)
          DO k_mpi = 1, MAX_TRACER_FORCING_VARS
              CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%fprefix(k_mpi), 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
              CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%vname(k_mpi),   256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
              CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%timelog(k_mpi), 256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
              CALL mpi_bcast(DEF_Tracer_Forcings_NL(i_mpi)%tintalgo(k_mpi),256, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
          END DO
      END DO
  ENDIF
  
  ! Final barrier to ensure all broadcasts are complete
  CALL mpi_barrier(p_comm_glb, p_err)
#endif

END SUBROUTINE broadcast_tracer_namelists

SUBROUTINE debug_tracer_forcing_configurations()
  IMPLICIT NONE
  
  integer :: i, j
  
  IF (.NOT. allocated(DEF_Tracers) .OR. .NOT. allocated(DEF_Tracer_Forcings_NL)) THEN
      WRITE(*,*) '========== Tracer Debug Information =========='
      WRITE(*,*) 'Warning: Tracer arrays not allocated'
      WRITE(*,*) '=============================================='
      RETURN
  ENDIF
  
  WRITE(*,*) '========== Tracer Debug Information =========='
  WRITE(*,*) 'Number of tracers configured: ', SIZE(DEF_Tracers)
  WRITE(*,*) 'Number of forced tracers: ', SIZE(DEF_Tracer_Forcings_NL)
  WRITE(*,*) ''
  
  ! Display basic tracer information
  WRITE(*,*) '--- Basic Tracer Information ---'
  DO i = 1, SIZE(DEF_Tracers)
      WRITE(*,*) 'Tracer ', i, ':'
      WRITE(*,*) '  Name: ', TRIM(DEF_Tracers(i)%name)
      WRITE(*,*) '  Type: ', TRIM(DEF_Tracers(i)%type)
  END DO
  WRITE(*,*) ''
  
  ! Display detailed forcing configurations
  WRITE(*,*) '--- Tracer Forcing Configurations ---'
  DO i = 1, SIZE(DEF_Tracer_Forcings_NL)
      WRITE(*,*) 'Forcing Config ', i, ':'
      WRITE(*,*) '  Tracer Name: ', TRIM(DEF_Tracer_Forcings_NL(i)%tracer_name)
      WRITE(*,*) '  Dataset Name: ', TRIM(DEF_Tracer_Forcings_NL(i)%dataset_name)
      WRITE(*,*) '  Tracer Type: ', TRIM(DEF_Tracer_Forcings_NL(i)%tracer_type)
      WRITE(*,*) '  Start Year: ', DEF_Tracer_Forcings_NL(i)%startyr
      WRITE(*,*) '  Start Month: ', DEF_Tracer_Forcings_NL(i)%startmo
      WRITE(*,*) '  End Year: ', DEF_Tracer_Forcings_NL(i)%endyr
      WRITE(*,*) '  End Month: ', DEF_Tracer_Forcings_NL(i)%endmo
      WRITE(*,*) '  Leap Year: ', DEF_Tracer_Forcings_NL(i)%leapyear
      WRITE(*,*) '  Data 2D: ', DEF_Tracer_Forcings_NL(i)%data2d
      WRITE(*,*) '  Hight Dim: ', DEF_Tracer_Forcings_NL(i)%hightdim
      WRITE(*,*) '  Dim 2D: ', DEF_Tracer_Forcings_NL(i)%dim2d
      WRITE(*,*) '  Lat Name: ', TRIM(DEF_Tracer_Forcings_NL(i)%latname)
      WRITE(*,*) '  Lon Name: ', TRIM(DEF_Tracer_Forcings_NL(i)%lonname)
      WRITE(*,*) '  Group By: ', TRIM(DEF_Tracer_Forcings_NL(i)%groupby)
      WRITE(*,*) '  Directory: ', TRIM(DEF_Tracer_Forcings_NL(i)%tracer_dir)
      WRITE(*,*) '  Number of Variables (NVAR): ', DEF_Tracer_Forcings_NL(i)%NVAR
      WRITE(*,*) '  Has Missing Value: ', DEF_Tracer_Forcings_NL(i)%has_missing_value
      WRITE(*,*) '  Missing Value Name: ', TRIM(DEF_Tracer_Forcings_NL(i)%missing_value_name)
      
      IF (DEF_Tracer_Forcings_NL(i)%NVAR > 0) THEN
          WRITE(*,*) '  Time intervals (dtime): ', (DEF_Tracer_Forcings_NL(i)%dtime(j), j=1, DEF_Tracer_Forcings_NL(i)%NVAR)
          WRITE(*,*) '  Offsets: ', (DEF_Tracer_Forcings_NL(i)%offset(j), j=1, DEF_Tracer_Forcings_NL(i)%NVAR)
          
          DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
              WRITE(*,*) '    Variable ', j, ':'
              WRITE(*,*) '      File prefix: ', TRIM(DEF_Tracer_Forcings_NL(i)%fprefix(j))
              WRITE(*,*) '      Variable name: ', TRIM(DEF_Tracer_Forcings_NL(i)%vname(j))
              WRITE(*,*) '      Time log: ', TRIM(DEF_Tracer_Forcings_NL(i)%timelog(j))
              WRITE(*,*) '      Time interpolation: ', TRIM(DEF_Tracer_Forcings_NL(i)%tintalgo(j))
          END DO
      ELSE
          WRITE(*,*) '  No forcing variables configured for this tracer'
      ENDIF
      WRITE(*,*) ''
  END DO
  
  WRITE(*,*) '=============================================='
  
END SUBROUTINE debug_tracer_forcing_configurations

END MODULE MOD_Tracer_Namelist_Defs
