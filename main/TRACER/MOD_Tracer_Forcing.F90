#include <define.h>

MODULE MOD_Tracer_Forcing

USE MOD_Precision, only: r8
USE MOD_Namelist, only: DEF_forcing_namelist, DEF_dir_forcing, nl_forcing_type, NVAR => forcing_NVAR_proxy ! Proxy for NVAR from nl_forcing_type if needed
USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL, MAX_TRACER_FORCING_VARS
USE MOD_DataType, only: block_data_real8_2d, pointer_real8_1d, allocate_block_data, flush_block_data, block_data_copy, block_data_linear_interp
USE MOD_Grid, only: grid_type
USE MOD_SpatialMapping, only: spatial_mapping_type
USE MOD_TimeManager, only: timestamp, calendarday, adj2begin, ticktime, isleapyear, get_calday, julian2monthday
USE MOD_SPMD_Task, only: p_is_master, p_is_io, p_is_worker, p_comm_glb, p_err, p_address_master, &
                       gblock, numpatch, mpi_tag_mesg, p_root, p_iam_worker, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER, &
                       MPI_SUM, mpi_bcast, mpi_barrier, CoLM_stop
USE MOD_NetCDFSerial, only: ncio_get_attr, ncio_var_exist, ncio_read_serial, ncio_read_bcast_serial, ncio_read_site_time, ncio_read_block_time
USE MOD_NetCDFVector, only: landpatch_vector_type ! Assuming this type exists for patch data
USE MOD_NetCDFBlock, only: ncio_def_var_real8_2d_blk, ncio_put_var_real8_2d_blk
USE MOD_UserDefFun, only: isnan_ud

IMPLICIT NONE
SAVE

PUBLIC :: tracer_forcing_init
PUBLIC :: read_tracer_forcing
PUBLIC :: get_current_tracer_forcing_value ! To access data from other modules

! Data structures to hold tracer forcing data
type(grid_type), allocatable :: gforc_tracers(:,:) ! Grid for each tracer forcing variable [num_forced_tracers, MAX_TRACER_FORCING_VARS]
type(spatial_mapping_type), allocatable :: mg2p_forc_tracers(:,:) ! Mapping for each [num_forced_tracers, MAX_TRACER_FORCING_VARS]

! Forcing data at current time step on patches
type(pointer_real8_1d), allocatable :: current_tracer_patch_forcing(:,:) ![num_forced_tracers, MAX_TRACER_FORCING_VARS]

! Time stamps for Lower/Upper Bound data for each tracer forcing variable
type(timestamp), allocatable :: tstamp_LB_tracer(:,:)  ! [num_forced_tracers, MAX_TRACER_FORCING_VARS]
type(timestamp), allocatable :: tstamp_UB_tracer(:,:)  ! [num_forced_tracers, MAX_TRACER_FORCING_VARS]

! Data buffers for Lower/Upper Bound data
type(block_data_real8_2d), allocatable :: forcn_LB_tracer(:,:,:) ! [num_forced_tracers, MAX_TRACER_FORCING_VARS, (block_data)]
type(block_data_real8_2d), allocatable :: forcn_UB_tracer(:,:,:) ! [num_forced_tracers, MAX_TRACER_FORCING_VARS, (block_data)]
type(block_data_real8_2d), allocatable :: forcn_tracer_tmp(:,:) ! Temp buffer for reading [num_forced_tracers, MAX_TRACER_FORCING_VARS]

integer :: num_tracers_with_forcing = 0
real(r8) :: deltim_model_sec = 1800.0_r8 ! Model timestep, to be set in init

! Proxy for NVAR from nl_forcing_type, used for DEF_forcing%startyr etc.
! This is a simplification, assuming tracer forcing files might share some global properties
! with main forcing, or these properties need to be added to nl_tracer_forcing_type.
INTEGER, PARAMETER :: forcing_NVAR_proxy = 8 


CONTAINS

SUBROUTINE tracer_forcing_init (model_timestep_sec, landpatch_data)
  USE MOD_LandPatch, only: landpatch ! For spatial mapping
  real(r8), intent(in) :: model_timestep_sec
  type(landpatch_vector_type), intent(in) :: landpatch_data ! Or actual landpatch type from MOD_LandPatch

  integer :: i, j, k, ierr
  character(len=256) :: filename_template, first_file_to_check
  real(r8), allocatable :: lat_in_tr(:), lon_in_tr(:)
  real(r8) :: missing_value_tracer(MAX_TRACER_FORCING_VARS)
  logical, allocatable :: forcmask_pch_tracer(:,:)


  deltim_model_sec = model_timestep_sec
  num_tracers_with_forcing = 0
  IF (allocated(DEF_Tracer_Forcings_NL)) THEN
     num_tracers_with_forcing = SIZE(DEF_Tracer_Forcings_NL, 1)
  ENDIF

  IF (num_tracers_with_forcing == 0 .OR. .NOT. allocated(DEF_Tracers)) THEN
     IF (p_is_master) WRITE(*,*) "MOD_Tracer_Forcing: No tracers with forcing defined or DEF_Tracers not allocated."
     RETURN
  ENDIF

  IF (p_is_master) WRITE(*,*) "MOD_Tracer_Forcing: Initializing for ", num_tracers_with_forcing, " tracers."

  ALLOCATE(gforc_tracers(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
  ALLOCATE(mg2p_forc_tracers(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
  ALLOCATE(current_tracer_patch_forcing(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
  ALLOCATE(tstamp_LB_tracer(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
  ALLOCATE(tstamp_UB_tracer(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
  
  IF (p_is_io) THEN
     ALLOCATE(forcn_LB_tracer(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS, 1)) ! Third dim is dummy for type
     ALLOCATE(forcn_UB_tracer(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS, 1))
     ALLOCATE(forcn_tracer_tmp(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
  ENDIF
  
  IF (numpatch > 0) THEN
     ALLOCATE(forcmask_pch_tracer(num_tracers_with_forcing, MAX_TRACER_FORCING_VARS))
     forcmask_pch_tracer = .TRUE.
  ENDIF


  DO i = 1, num_tracers_with_forcing
     DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
        IF (TRIM(DEF_Tracer_Forcings_NL(i)%fprefix(j)) == 'null' .OR. TRIM(DEF_Tracer_Forcings_NL(i)%vname(j)) == 'null') CYCLE
        
        tstamp_LB_tracer(i,j) = timestamp(-1,-1,-1)
        tstamp_UB_tracer(i,j) = timestamp(-1,-1,-1)

        ! Determine grid from the first forcing file of the first variable of the first tracer
        ! This assumes all tracer forcings are on the same grid.
        IF (i == 1 .AND. j == 1) THEN
           ! Construct filename: Need year/month/day. Use main forcing start for now.
           ! This part needs robust determination of start year/month for tracer files.
           ! For now, using a placeholder or assuming it matches main forcing start year/month.
           ! Let's assume DEF_forcing (main namelist) is available or use hardcoded start.
           integer :: start_yr_main, start_mo_main
           start_yr_main = 2000 ! Placeholder, should come from main forcing namelist if possible
           start_mo_main = 1
           
           filename_template = TRIM(DEF_dir_forcing) // TRIM(DEF_Tracer_Forcings_NL(i)%fprefix(j))
           first_file_to_check = get_tracer_metfilename(filename_template, start_yr_main, start_mo_main, 1) ! Assuming monthly files for simplicity

           IF (p_is_master) THEN
             IF (ncio_var_exist(first_file_to_check, 'lat') .AND. ncio_var_exist(first_file_to_check, 'lon')) THEN
                 CALL ncio_read_bcast_serial (first_file_to_check, 'lat', lat_in_tr)
                 CALL ncio_read_bcast_serial (first_file_to_check, 'lon', lon_in_tr)
                 CALL gforc_tracers(i,j)%define_by_center(lat_in_tr, lon_in_tr)
                 DEALLOCATE(lat_in_tr, lon_in_tr)
             ELSEIF (ncio_var_exist(first_file_to_check, 'LATIXY') .AND. ncio_var_exist(first_file_to_check, 'LONGXY')) THEN ! CLM style
                 real(r8), allocatable :: latxy_tr(:,:), lonxy_tr(:,:)
                 CALL ncio_read_bcast_serial (first_file_to_check, 'LATIXY', latxy_tr)
                 CALL ncio_read_bcast_serial (first_file_to_check, 'LONGXY', lonxy_tr)
                 ALLOCATE(lat_in_tr(SIZE(latxy_tr,2)))
                 ALLOCATE(lon_in_tr(SIZE(lonxy_tr,1)))
                 lat_in_tr = latxy_tr(1,:)
                 lon_in_tr = lonxy_tr(:,1)
                 CALL gforc_tracers(i,j)%define_by_center(lat_in_tr, lon_in_tr)
                 DEALLOCATE(latxy_tr, lonxy_tr, lat_in_tr, lon_in_tr)
             ELSE
                 WRITE(*,*) "MOD_Tracer_Forcing: Lat/Lon not found in example file: ", TRIM(first_file_to_check)
                 CALL CoLM_stop()
             ENDIF
             CALL gforc_tracers(i,j)%set_rlon()
             CALL gforc_tracers(i,j)%set_rlat()
           ENDIF
#ifdef USEMPI
           CALL mpi_bcast_grid_type(gforc_tracers(i,j)) ! Custom broadcast for grid_type needed
#endif
        ELSE
           CALL gforc_tracers(i,j)%clone_grid(gforc_tracers(1,1))
        ENDIF

        IF (p_is_io) THEN
           CALL allocate_block_data(gforc_tracers(i,j), forcn_LB_tracer(i,j,1))
           CALL allocate_block_data(gforc_tracers(i,j), forcn_UB_tracer(i,j,1))
           CALL allocate_block_data(gforc_tracers(i,j), forcn_tracer_tmp(i,j))
        ENDIF
        
        IF (numpatch > 0) THEN
           ALLOCATE(current_tracer_patch_forcing(i,j)%val(numpatch))
           current_tracer_patch_forcing(i,j)%val = 0.0_r8
        ENDIF

        CALL mg2p_forc_tracers(i,j)%build_arealweighted(gforc_tracers(i,j), landpatch_data) 
        ! Assuming 'arealweight' for now, could be made configurable per tracer
        
        ! Simplified missing value handling - assume not present or use a common one
        ! CALL mg2p_forc_tracers(i,j)%set_missing_value (forcn_tracer_tmp(i,j), common_missing_value, forcmask_pch_tracer(i,j,:))

     END DO ! Loop NVAR for tracer
  END DO ! Loop tracers
  IF (p_is_master) WRITE(*,*) "MOD_Tracer_Forcing: Initialization complete."
END SUBROUTINE tracer_forcing_init

SUBROUTINE read_tracer_forcing (idate_current)
  type(timestamp), intent(in) :: idate_current
  integer :: i, j, k
  integer :: yr, mo, dy, sec_of_day, time_idx_in_file
  character(len=256) :: file_path_lb, file_path_ub
  real(r8) :: dtLB_ratio, dtUB_ratio
  integer :: dtLB_sec, dtUB_sec
  type(block_data_real8_2d) :: temp_metdata_block ! Temporary block for reading

  IF (num_tracers_with_forcing == 0) RETURN

  DO i = 1, num_tracers_with_forcing
     DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
         IF (TRIM(DEF_Tracer_Forcings_NL(i)%fprefix(j)) == 'null' .OR. TRIM(DEF_Tracer_Forcings_NL(i)%vname(j)) == 'null') CYCLE

         ! Check if current data is still valid
         IF ( .NOT. (tstamp_LB_tracer(i,j)%year == -1) .AND. .NOT. (tstamp_UB_tracer(i,j)%year == -1) .AND. &
              idate_current >= tstamp_LB_tracer(i,j) .AND. idate_current < tstamp_UB_tracer(i,j) ) THEN
             ! Data bounds still valid, just interpolate
         ELSE
             ! Need to read new data bounds
             IF (p_is_io) THEN
                 CALL tracer_metread_LBUB(idate_current, i, j, DEF_Tracer_Forcings_NL(i), &
                                        tstamp_LB_tracer(i,j), tstamp_UB_tracer(i,j), &
                                        forcn_LB_tracer(i,j,1), forcn_UB_tracer(i,j,1), &
                                        forcn_tracer_tmp(i,j) )
             ENDIF
#ifdef USEMPI
             CALL mpi_bcast_timestamp(tstamp_LB_tracer(i,j), p_is_io)
             CALL mpi_bcast_timestamp(tstamp_UB_tracer(i,j), p_is_io)
             IF (p_is_io) THEN
                CALL mpi_bcast_block_data(forcn_LB_tracer(i,j,1), .TRUE.)
                CALL mpi_bcast_block_data(forcn_UB_tracer(i,j,1), .TRUE.)
             ELSE
                CALL mpi_bcast_block_data(forcn_LB_tracer(i,j,1), .FALSE.)
                CALL mpi_bcast_block_data(forcn_UB_tracer(i,j,1), .FALSE.)
             ENDIF
#endif
         ENDIF

         ! Temporal Interpolation
         dtLB_sec = idate_current - tstamp_LB_tracer(i,j)
         dtUB_sec = tstamp_UB_tracer(i,j) - idate_current
         
         IF (p_is_io) THEN
             IF ((dtLB_sec + dtUB_sec) > 0) THEN
                 dtLB_ratio = REAL(dtLB_sec, r8) / REAL(dtLB_sec + dtUB_sec, r8)
                 dtUB_ratio = REAL(dtUB_sec, r8) / REAL(dtLB_sec + dtUB_sec, r8)
                 CALL block_data_linear_interp(forcn_LB_tracer(i,j,1), dtUB_ratio, &
                                               forcn_UB_tracer(i,j,1), dtLB_ratio, &
                                               forcn_tracer_tmp(i,j) ) 
             ELSE
                 CALL block_data_copy(forcn_LB_tracer(i,j,1), forcn_tracer_tmp(i,j))
             ENDIF
         ENDIF
#ifdef USEMPI
         IF (p_is_io) THEN
             CALL mpi_bcast_block_data(forcn_tracer_tmp(i,j), .TRUE.)
         ELSE
             CALL mpi_bcast_block_data(forcn_tracer_tmp(i,j), .FALSE.)
         ENDIF
#endif
         ! Spatial Interpolation (Mapping)
         IF (numpatch > 0 .AND. allocated(current_tracer_patch_forcing(i,j)%val)) THEN
              CALL mg2p_forc_tracers(i,j)%grid2pset(forcn_tracer_tmp(i,j), current_tracer_patch_forcing(i,j)%val)
              ! Add missing value handling if necessary using forcmask_pch_tracer
         ENDIF
     END DO ! NVAR loop
  END DO ! Tracer loop
END SUBROUTINE read_tracer_forcing

SUBROUTINE tracer_metread_LBUB(idate_model, tracer_idx, var_idx, tracer_nl_config, &
                             tstamp_LB_out, tstamp_UB_out, &
                             data_LB_out, data_UB_out, data_tmp_read_buf)
  type(timestamp), intent(in) :: idate_model
  integer, intent(in) :: tracer_idx, var_idx
  type(nl_tracer_forcing_type), intent(in) :: tracer_nl_config
  type(timestamp), intent(out) :: tstamp_LB_out, tstamp_UB_out
  type(block_data_real8_2d), intent(out) :: data_LB_out, data_UB_out
  type(block_data_real8_2d), intent(inout) :: data_tmp_read_buf ! Buffer for ncio_read_block_time

  integer :: yr_lb, mo_lb, dy_lb, time_idx_lb
  integer :: yr_ub, mo_ub, dy_ub, time_idx_ub
  character(len=256) :: fpath_lb, fpath_ub
  character(len=256) :: filename_template
  
  filename_template = TRIM(DEF_dir_forcing) // TRIM(tracer_nl_config%fprefix(var_idx))

  ! Set Lower Bound
  CALL get_tracer_forcing_timestamp_and_idx(idate_model, tracer_nl_config%dtime(var_idx), tracer_nl_config%offset(var_idx), &
                                         'LB', tstamp_LB_out, yr_lb, mo_lb, dy_lb, time_idx_lb)
  fpath_lb = get_tracer_metfilename(filename_template, yr_lb, mo_lb, dy_lb)
  CALL ncio_read_block_time(fpath_lb, TRIM(tracer_nl_config%vname(var_idx)), gforc_tracers(tracer_idx, var_idx), time_idx_lb, data_tmp_read_buf)
  CALL block_data_copy(data_tmp_read_buf, data_LB_out)

  ! Set Upper Bound
  CALL get_tracer_forcing_timestamp_and_idx(idate_model, tracer_nl_config%dtime(var_idx), tracer_nl_config%offset(var_idx), &
                                         'UB', tstamp_UB_out, yr_ub, mo_ub, dy_ub, time_idx_ub, tstamp_LB_out)
  fpath_ub = get_tracer_metfilename(filename_template, yr_ub, mo_ub, dy_ub)
  CALL ncio_read_block_time(fpath_ub, TRIM(tracer_nl_config%vname(var_idx)), gforc_tracers(tracer_idx, var_idx), time_idx_ub, data_tmp_read_buf)
  CALL block_data_copy(data_tmp_read_buf, data_UB_out)

END SUBROUTINE tracer_metread_LBUB

SUBROUTINE get_tracer_forcing_timestamp_and_idx(model_time, dtime_var, offset_var, bound_type, &
                                              tstamp_out, yr_out, mo_out, dy_out, time_idx_file_out, &
                                              tstamp_lb_ref)
  type(timestamp), intent(in) :: model_time
  integer, intent(in) :: dtime_var, offset_var ! dtime in seconds
  character(len=*), intent(in) :: bound_type ! "LB" or "UB"
  type(timestamp), intent(out) :: tstamp_out
  integer, intent(out) :: yr_out, mo_out, dy_out, time_idx_file_out
  type(timestamp), intent(in), optional :: tstamp_lb_ref ! For UB calculation

  integer :: model_total_seconds_of_year, target_total_seconds
  integer :: num_intervals_in_day
  integer :: day_of_year_model
  
  day_of_year_model = model_time%day
  model_total_seconds_of_year = (day_of_year_model - 1) * 86400 + model_time%sec

  IF (dtime_var == 0) THEN ! Daily data
     num_intervals_in_day = 1
  ELSE
     num_intervals_in_day = 86400 / dtime_var
  ENDIF

  IF (bound_type == 'LB') THEN
     target_total_seconds = INT((REAL(model_total_seconds_of_year - offset_var) / REAL(dtime_var))) * dtime_var + offset_var
     IF (target_total_seconds > model_total_seconds_of_year .AND. model_total_seconds_of_year >= offset_var) THEN
         target_total_seconds = target_total_seconds - dtime_var
     ENDIF
     tstamp_out%year = model_time%year
     tstamp_out%day  = target_total_seconds / 86400 + 1
     tstamp_out%sec  = MOD(target_total_seconds, 86400)
     CALL adj2begin(tstamp_out%year, tstamp_out%day, tstamp_out%sec) ! Normalize

     time_idx_file_out = ( (tstamp_out%day -1) * num_intervals_in_day ) + &
                           (tstamp_out%sec / dtime_var) + 1
     
  ELSE IF (bound_type == 'UB') THEN
     IF (.NOT. PRESENT(tstamp_lb_ref)) THEN
         WRITE(*,*) "Error: tstamp_lb_ref required for UB calculation in get_tracer_forcing_timestamp_and_idx"
         CALL CoLM_stop()
     ENDIF
     tstamp_out = tstamp_lb_ref + dtime_var 
     CALL adj2begin(tstamp_out%year, tstamp_out%day, tstamp_out%sec) ! Normalize

     time_idx_file_out = ( (tstamp_out%day -1) * num_intervals_in_day ) + &
                           (tstamp_out%sec / dtime_var) + 1
  ENDIF
  
  yr_out = tstamp_out%year
  CALL julian2monthday(yr_out, tstamp_out%day, mo_out, dy_out)
  ! Assuming file grouping is monthly for get_tracer_metfilename, so dy_out might not be used by it.
  ! If daily files, dy_out is important.
  
END SUBROUTINE get_tracer_forcing_timestamp_and_idx

FUNCTION get_tracer_metfilename(fprefix_template, year, month, day_of_month_or_year) RESULT(fname)
  character(len=*), intent(in) :: fprefix_template
  integer, intent(in) :: year, month, day_of_month_or_year
  character(len=256) :: fname
  character(len=4) :: cyear
  character(len=2) :: cmonth, cday

  WRITE(cyear,'(I4.4)') year
  WRITE(cmonth,'(I2.2)') month
  WRITE(cday,'(I2.2)') day_of_month_or_year ! Assuming day means day of month for daily files

  fname = fprefix_template ! This needs to be smarter, like in MOD_UserSpecifiedForcing
  ! Replace YYYY, MM, DD in template. For now, assume simple append or fixed name for example.
  ! This is a placeholder for robust file name generation.
  ! Example: if fprefix_template is 'TracerX_YYYYMM.nc', replace YYYY and MM.
  ! For now, assume fprefix_template IS the filename if not using date parts, or needs specific logic.
  ! Based on example: 'IsoGSM_prate.YYYY-MM.nc'
  integer :: idx_yyyy, idx_mm
  idx_yyyy = INDEX(fprefix_template, 'YYYY')
  idx_mm   = INDEX(fprefix_template, 'MM')

  IF (idx_yyyy > 0 .AND. idx_mm > 0) THEN
     fname = fprefix_template(1:idx_yyyy-1) // cyear // fprefix_template(idx_yyyy+4:idx_mm-1) // cmonth // fprefix_template(idx_mm+2:)
  ELSE IF (idx_yyyy > 0) THEN ! Only year
     fname = fprefix_template(1:idx_yyyy-1) // cyear // fprefix_template(idx_yyyy+4:)
  ELSE
     fname = TRIM(fprefix_template) ! Assumes template is full path or relative path that doesn't need date substitution
  ENDIF
  ! Add more logic if DD (day) is also part of the filename pattern.
END FUNCTION get_tracer_metfilename

FUNCTION get_current_tracer_forcing_value(tracer_name_in, forcing_var_idx, patch_idx) RESULT(val_out)
     character(len=*), intent(in) :: tracer_name_in
     integer, intent(in) :: forcing_var_idx ! 1-based index for the NVAR of that tracer
     integer, intent(in) :: patch_idx
     real(r8) :: val_out

     integer :: i
     val_out = 0.0_r8 ! Default / not found

     IF (num_tracers_with_forcing == 0 .OR. .NOT. allocated(current_tracer_patch_forcing)) RETURN

     DO i = 1, num_tracers_with_forcing
         IF (TRIM(ADJUSTL(DEF_Tracer_Forcings_NL(i)%tracer_name)) == TRIM(ADJUSTL(tracer_name_in))) THEN
             IF (forcing_var_idx > 0 .AND. forcing_var_idx <= MAX_TRACER_FORCING_VARS) THEN
                IF (allocated(current_tracer_patch_forcing(i, forcing_var_idx)%val)) THEN
                   IF (patch_idx > 0 .AND. patch_idx <= SIZE(current_tracer_patch_forcing(i, forcing_var_idx)%val)) THEN
                      val_out = current_tracer_patch_forcing(i, forcing_var_idx)%val(patch_idx)
                      RETURN
                   ENDIF
                ENDIF
             ENDIF
             EXIT
         ENDIF
     ENDDO
END FUNCTION get_current_tracer_forcing_value

SUBROUTINE mpi_bcast_grid_type(grid_var) ! Simplified broadcast for grid_type
     type(grid_type), intent(inout) :: grid_var
#ifdef USEMPI
     CALL mpi_bcast(grid_var%nlon, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
     CALL mpi_bcast(grid_var%nlat, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
     ! Add other members of grid_type that need broadcasting
     ! This is a placeholder; a full broadcast is more complex.
     IF (.NOT. p_is_master) THEN
         IF (grid_var%nlon > 0 .AND. grid_var%nlat > 0) THEN
             ! Allocate arrays based on nlon, nlat if they are allocatable in grid_type
         ENDIF
     ENDIF
     ! Then broadcast array contents if any
#endif
END SUBROUTINE mpi_bcast_grid_type

SUBROUTINE mpi_bcast_timestamp(ts_var, is_sender)
    type(timestamp), intent(inout) :: ts_var
    logical, intent(in) :: is_sender ! True if current rank is the sender (master/io)
#ifdef USEMPI
    CALL mpi_bcast(ts_var%year, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
    CALL mpi_bcast(ts_var%day,  1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
    CALL mpi_bcast(ts_var%sec,  1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
END SUBROUTINE mpi_bcast_timestamp

SUBROUTINE mpi_bcast_block_data(bd_var, is_sender)
    type(block_data_real8_2d), intent(inout) :: bd_var
    logical, intent(in) :: is_sender
#ifdef USEMPI
    integer :: iblk, jblk
    IF (.NOT. is_sender) THEN
        ! Non-senders need to allocate based on grid info received prior
        ! Assuming gforc_tracers(i,j) is already broadcast and bd_var is for that grid
        ! This allocation should happen after grid is known on non-sender
        ! CALL allocate_block_data(grid_ref, bd_var) ! Needs grid_ref
    ENDIF
    IF (allocated(bd_var%blk)) THEN
      DO jblk = 1, gblock%nyblk
         DO iblk = 1, gblock%nxblk
            IF (gblock%is_me(iblk,jblk)) THEN
               IF (allocated(bd_var%blk(iblk,jblk)%val)) THEN
                  CALL mpi_bcast(bd_var%blk(iblk,jblk)%val, &
                                 SIZE(bd_var%blk(iblk,jblk)%val), MPI_REAL8, &
                                 p_address_master, p_comm_glb, p_err)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
    ENDIF
#endif
END SUBROUTINE mpi_bcast_block_data


END MODULE MOD_Tracer_Forcing
