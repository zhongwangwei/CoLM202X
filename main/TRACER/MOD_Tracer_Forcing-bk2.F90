#include <define.h>

MODULE MOD_Tracer_Forcing

   !-----------------------------------------------------------------------
   ! !DESCRIPTION:
   !  read in the atmospheric forcing using user defined interpolation method or
   !  downscaling forcing
   !
   ! !REVISIONS:
   !  Yongjiu Dai and Hua Yuan, 04/2014: initial code from CoLM2014 (metdata.F90,
   !                                     GETMET.F90 and rd_forcing.F90
   !
   !  Shupeng Zhang, 05/2023: 1) porting codes to MPI parallel version
   !                          2) codes for dealing with missing forcing value
   !                          3) interface for downscaling
   !
   ! !TODO...(need complement)
   !-----------------------------------------------------------------------
   
   USE MOD_TimeManager, only: timestamp, calendarday, adj2begin, ticktime, isleapyear, get_calday, julian2monthday

   
   CONTAINS
   
   !-----------------------------------------------------------------------
      SUBROUTINE forcing_init (dir_forcing, deltatime, ststamp, lc_year, etstamp, lulcc_call)
   
      USE MOD_SPMD_Task
      USE MOD_Namelist
      USE MOD_Block
      USE MOD_DataType
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
   #ifdef CROP
      USE MOD_LandCrop
   #endif
      USE MOD_UserSpecifiedForcing
      USE MOD_NetCDFSerial
      USE MOD_NetCDFVector
      USE MOD_NetCDFBlock
      USE MOD_Vars_TimeInvariants
      USE MOD_Vars_1DForcing
      IMPLICIT NONE
   
      character(len=*), intent(in) :: dir_forcing
      real(r8),         intent(in) :: deltatime  ! model time step
      type(timestamp),  intent(in) :: ststamp
      integer,          intent(in) :: lc_year    ! which year of land cover data used
      type(timestamp),  intent(in), optional :: etstamp
      logical,          intent(in), optional :: lulcc_call ! whether it is a lulcc CALL
   
      ! Local variables
      integer            :: idate(3)
      type(timestamp)    :: tstamp
      character(len=256) :: filename, lndname, cyear
      integer            :: ivar, year, month, day, time_i
      real(r8)           :: missing_value
      integer            :: ielm, istt, iend
   
      integer :: iblkme, xblk, yblk, xloc, yloc
   
         CALL init_user_specified_forcing
   
         ! CO2 data initialization
         CALL init_monthly_co2_mlo
   
         ! get value of fmetdat and deltim
         deltim_int  = int(deltatime)
         ! deltim_real = deltatime
   
         ! set initial values
         IF (allocated(tstamp_LB)) deallocate(tstamp_LB)
         IF (allocated(tstamp_UB)) deallocate(tstamp_UB)
         allocate (tstamp_LB(NVAR))
         allocate (tstamp_UB(NVAR))
         tstamp_LB(:) = timestamp(-1, -1, -1)
         tstamp_UB(:) = timestamp(-1, -1, -1)
   
         idate = (/ststamp%year, ststamp%day, ststamp%sec/)
         CALL adj2begin (idate)
   
         CALL metread_latlon (dir_forcing, idate)
   
         IF (p_is_io) THEN
   
            IF (allocated(forcn   )) deallocate(forcn   )
            IF (allocated(forcn_LB)) deallocate(forcn_LB)
            IF (allocated(forcn_UB)) deallocate(forcn_UB)
            allocate (forcn    (NVAR))
            allocate (forcn_LB (NVAR))
            allocate (forcn_UB (NVAR))
   
            DO ivar = 1, NVAR
               CALL allocate_block_data (gforc, forcn   (ivar))
               CALL allocate_block_data (gforc, forcn_LB(ivar))
               CALL allocate_block_data (gforc, forcn_UB(ivar))
            ENDDO
   
            ! allocate memory for forcing data
            CALL allocate_block_data (gforc, metdata)  ! forcing data
            CALL allocate_block_data (gforc, avgcos )  ! time-average of cos(zenith)
   #if(defined URBAN_MODEL && defined SinglePoint)
            CALL allocate_block_data (gforc, rainf)
            CALL allocate_block_data (gforc, snowf)
   #endif
   
         ENDIF
   
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (forcmask_pch(numpatch));  forcmask_pch(:) = .true.
            ENDIF
         ENDIF
   
         IF (DEF_forcing%has_missing_value) THEN
   
            tstamp = idate
            CALL setstampLB(tstamp, 1, year, month, day, time_i)
            filename = trim(dir_forcing)//trim(metfilename(year, month, day, 1))
            tstamp_LB(1) = timestamp(-1, -1, -1)
   
            IF (p_is_master) THEN
               CALL ncio_get_attr (filename, vname(1), trim(DEF_forcing%missing_value_name), missing_value)
            ENDIF
   #ifdef USEMPI
            CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
   #endif
   
            CALL ncio_read_block_time (filename, vname(1), gforc, time_i, metdata)
   
         ENDIF
   
         IF (trim(DEF_Forcing_Interp_Method) == 'arealweight') THEN
            IF (present(lulcc_call)) CALL mg2p_forc%forc_free_mem
            CALL mg2p_forc%build_arealweighted (gforc, landpatch)
         ELSEIF (trim(DEF_Forcing_Interp_Method) == 'bilinear') THEN
            IF (present(lulcc_call)) CALL mg2p_forc%forc_free_mem
            CALL mg2p_forc%build_bilinear (gforc, landpatch)
         ENDIF
   
         IF (DEF_forcing%has_missing_value) THEN
            CALL mg2p_forc%set_missing_value (metdata, missing_value, forcmask_pch)
         ENDIF
   
         IF (DEF_USE_Forcing_Downscaling) THEN
   
            IF (p_is_worker .and. (numpatch > 0)) THEN
               forc_topo = topoelv
               WHERE(forc_topo == spval) forc_topo = 0.
            ENDIF
   
            IF (p_is_io) CALL allocate_block_data (gforc, topo_grid)
            CALL mg2p_forc%pset2grid (forc_topo, topo_grid)
   
            IF (p_is_io) CALL allocate_block_data (gforc, sumarea_grid)
            CALL mg2p_forc%get_sumarea (sumarea_grid)
   
            CALL block_data_division (topo_grid, sumarea_grid)
   
            IF (p_is_io) CALL allocate_block_data (gforc, maxelv_grid)
            CALL mg2p_forc%pset2grid_max (forc_topo, maxelv_grid)
   
   
            CALL mg2p_forc%allocate_part (forc_topo_grid  )
            CALL mg2p_forc%allocate_part (forc_maxelv_grid)
   
            CALL mg2p_forc%allocate_part (forc_t_grid     )
            CALL mg2p_forc%allocate_part (forc_th_grid    )
            CALL mg2p_forc%allocate_part (forc_q_grid     )
   #ifdef USE_ISOTOPE
            CALL mg2p_forc%allocate_part (forc_q_grid_O18 )
            CALL mg2p_forc%allocate_part (forc_q_grid_H2  )
   #endif
            CALL mg2p_forc%allocate_part (forc_pbot_grid  )
            CALL mg2p_forc%allocate_part (forc_rho_grid   )
            CALL mg2p_forc%allocate_part (forc_prc_grid   )
            CALL mg2p_forc%allocate_part (forc_prl_grid   )
   
   #ifdef USE_ISOTOPE
            CALL mg2p_forc%allocate_part (forc_prc_grid_O18 )
            CALL mg2p_forc%allocate_part (forc_prl_grid_O18 )
            CALL mg2p_forc%allocate_part (forc_prc_grid_H2  )
            CALL mg2p_forc%allocate_part (forc_prl_grid_H2  )
   #endif
            CALL mg2p_forc%allocate_part (forc_lwrad_grid )
            CALL mg2p_forc%allocate_part (forc_swrad_grid )
            CALL mg2p_forc%allocate_part (forc_hgt_grid   )
            CALL mg2p_forc%allocate_part (forc_us_grid    )
            CALL mg2p_forc%allocate_part (forc_vs_grid    )
   
            CALL mg2p_forc%allocate_part (forc_t_part     )
            CALL mg2p_forc%allocate_part (forc_th_part    )
            CALL mg2p_forc%allocate_part (forc_q_part     )
   #ifdef USE_ISOTOPE
            CALL mg2p_forc%allocate_part (forc_q_part_O18    )
            CALL mg2p_forc%allocate_part (forc_q_part_H2    )
   #endif
            CALL mg2p_forc%allocate_part (forc_pbot_part  )
            CALL mg2p_forc%allocate_part (forc_rhoair_part)
            CALL mg2p_forc%allocate_part (forc_prc_part   )
            CALL mg2p_forc%allocate_part (forc_prl_part   )
            CALL mg2p_forc%allocate_part (forc_frl_part   )
   #ifdef USE_ISOTOPE
            CALL mg2p_forc%allocate_part (forc_prc_part_O18  )
            CALL mg2p_forc%allocate_part (forc_prl_part_O18  )
            CALL mg2p_forc%allocate_part (forc_prc_part_H2  )
            CALL mg2p_forc%allocate_part (forc_prl_part_H2  )
   #endif
            CALL mg2p_forc%allocate_part (forc_swrad_part )
            CALL mg2p_forc%allocate_part (forc_us_part    )
            CALL mg2p_forc%allocate_part (forc_vs_part    )
   
            CALL mg2p_forc%grid2part (topo_grid,   forc_topo_grid  )
            CALL mg2p_forc%grid2part (maxelv_grid, forc_maxelv_grid)
   
            IF (p_is_worker .and. (numpatch > 0)) THEN
               allocate (glacierss(numpatch))
               glacierss(:) = patchtype(:) == 3
            ENDIF
   
         ENDIF
   
         forcing_read_ahead = .false.
         IF (trim(DEF_forcing%dataset) == 'POINT') THEN
            IF (USE_SITE_ForcingReadAhead .and. present(etstamp)) THEN
               forcing_read_ahead = .true.
               CALL metread_time (dir_forcing, ststamp, etstamp, deltatime)
            ELSE
               CALL metread_time (dir_forcing)
            ENDIF
            allocate (iforctime(NVAR))
         ENDIF
   
         IF (trim(DEF_forcing%dataset) == 'POINT') THEN
   
            filename = trim(dir_forcing)//trim(fprefix(1))
   
   
         ENDIF
   
      END SUBROUTINE forcing_init



      SUBROUTINE tracer_forcing_init (model_timestep_sec,ststamp, landpatch_data)
         USE MOD_LandPatch, only: landpatch ! For spatial mapping
         real(r8), intent(in) :: model_timestep_sec
         type(timestamp),  intent(in) :: ststamp

         type(landpatch_vector_type), intent(in) :: landpatch_data ! Or actual landpatch type from MOD_LandPatch
       
         integer :: i, j, k, ierr
         character(len=256) :: filename_template, first_file_to_check
         real(r8), allocatable :: lat_in_tr(:), lon_in_tr(:)
         real(r8) :: missing_value_tracer(MAX_TRACER_FORCING_VARS)
         logical, allocatable :: forcmask_pch_tracer(:,:)
       
       
         deltim_model_sec = int(model_timestep_sec)
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
                  
               idate = (/ststamp%year, ststamp%day, ststamp%sec/)
               CALL adj2begin (idate)
   
               CALL metread_latlon (DEF_Tracer_Forcings_NL(i)%tracer_dir, idate)
       
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
   END MODULE MOD_Tracer_Forcing










