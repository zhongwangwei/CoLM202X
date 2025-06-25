#include <define.h>

MODULE MOD_Tracer_Forcing



   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Grid
   USE MOD_SpatialMapping
   USE MOD_TimeManager
   USE MOD_SPMD_Task
   USE MOD_MonthlyinSituCO2MaunaLoa
   USE MOD_Vars_Global, only: pi
   USE MOD_OrbCoszen
   USE MOD_UserDefFun
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL, MAX_TRACER_FORCING_VARS

   IMPLICIT NONE

   !type (grid_type), PUBLIC :: gforc

   !type (spatial_mapping_type) :: mg2p_forc   ! area weighted mapping from forcing to model unit

   type(grid_type), allocatable :: gforc_tracers(:,:) ! Grid for each tracer forcing variable [num_forced_tracers, max_nvar]
   type(spatial_mapping_type), allocatable :: mg2p_forc_tracers(:,:) ! Mapping for each [num_forced_tracers, max_nvar]

   logical, allocatable :: forcmask_pch (:)
   ! local variables
   integer  :: deltim_int                ! model time step length
   ! real(r8) :: deltim_real             ! model time step length

   ! Variables for array allocation
   integer :: max_nvar                   ! maximum number of variables across all tracers

   !  for SinglePoint
   type(timestamp), allocatable :: forctime (:)
   integer,  allocatable :: iforctime(:)

   logical :: forcing_read_ahead
   real(r8), allocatable :: forc_disk(:,:)

   type(timestamp), allocatable :: tstamp_LB(:,:)  ! time stamp of low boundary data
   type(timestamp), allocatable :: tstamp_UB(:,:)  ! time stamp of up boundary data

   type(block_data_real8_2d),allocatable :: traceravgcos(:,:)   ! time-average of cos(zenith)
   type(block_data_real8_2d),allocatable :: tracerdata(:,:)  ! forcing data


   type(block_data_real8_2d), allocatable :: forcn_tracer(:,:)  ! forcing data
   type(block_data_real8_2d), allocatable :: forcn_LB_tracer(:,:)  ! forcing data at lower boundary
   type(block_data_real8_2d), allocatable :: forcn_UB_tracer(:,:)  ! forcing data at upper boundary
   
   ! Tracer patch data - stores tracer forcing data on patches
   ! Structure: [tracer_idx][var_idx] -> patch data array
   real(r8), allocatable, target :: tracer_patch_data(:,:,:)  ! [num_tracers, max_nvar, numpatch]
   
   PUBLIC :: tracer_forcing_init
   PUBLIC :: read_tracer_forcing
   PUBLIC :: tracer_forcing_final
   PUBLIC :: tracer_forcing_reset
   PUBLIC :: get_tracer_forcing_data

CONTAINS

   SUBROUTINE tracer_forcing_init (deltatime, ststamp, lc_year, etstamp, lulcc_call)

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
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing
   IMPLICIT NONE

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
   integer            :: i, j, k

   integer :: iblkme, xblk, yblk, xloc, yloc

   integer :: num_tracers_with_forcing
   
   ! Check allocation status of tracer arrays
   IF (p_is_master) THEN
      WRITE(*,*) "=== TRACER_FORCING_INIT START ==="
      WRITE(*,*) "  deltatime: ", deltatime
      WRITE(*,*) "  ststamp: ", ststamp%year, ststamp%day, ststamp%sec
      WRITE(*,*) "  lc_year: ", lc_year
      IF (present(etstamp)) WRITE(*,*) "  etstamp: ", etstamp%year, etstamp%day, etstamp%sec
      IF (present(lulcc_call)) WRITE(*,*) "  lulcc_call: ", lulcc_call
   ENDIF
   
   num_tracers_with_forcing = SIZE(DEF_Tracer_Forcings_NL, 1)
   
   IF (p_is_master) THEN
      WRITE(*,*) "  num_tracers_with_forcing: ", num_tracers_with_forcing
      WRITE(*,*) "  DEF_Tracer_Forcings_NL allocated: ", allocated(DEF_Tracer_Forcings_NL)
   ENDIF



   ! get value of deltim (moved outside loop)
   deltim_int  = int(deltatime)

   ! Calculate maximum number of variables across all tracers (moved outside loop)
   max_nvar = 0
   DO k = 1, num_tracers_with_forcing
      IF (DEF_Tracer_Forcings_NL(k)%NVAR > max_nvar) THEN
         max_nvar = DEF_Tracer_Forcings_NL(k)%NVAR
      ENDIF
   ENDDO
   
   ! Check if any tracer has variables defined
   IF (max_nvar == 0) THEN
      IF (p_is_master) THEN
         WRITE(*,*) "  WARNING: No tracer variables defined (max_nvar = 0)"
         DO k = 1, num_tracers_with_forcing
            WRITE(*,*) "    Tracer ", k, ": ", trim(DEF_Tracer_Forcings_NL(k)%tracer_name), &
                       " NVAR = ", DEF_Tracer_Forcings_NL(k)%NVAR
         ENDDO
      ENDIF
      RETURN  ! Exit early if no variables are defined
   ENDIF

   ! Initialize gforc_tracers and mg2p_forc_tracers
   allocate(gforc_tracers(num_tracers_with_forcing, max_nvar))
   allocate(mg2p_forc_tracers(num_tracers_with_forcing, max_nvar))
   

   IF (allocated(tstamp_LB)) deallocate(tstamp_LB)
   IF (allocated(tstamp_UB)) deallocate(tstamp_UB)

   allocate (tstamp_LB(num_tracers_with_forcing, max_nvar))
   allocate (tstamp_UB(num_tracers_with_forcing, max_nvar))


   tstamp_LB(:,:) = timestamp(-1, -1, -1)
   tstamp_UB(:,:) = timestamp(-1, -1, -1)

   ! Initialize global tracer patch data array (only once)
   ! It needs to be allocated on all processes to be passed as an argument.
   IF (allocated(tracer_patch_data)) deallocate(tracer_patch_data)
   IF (numpatch > 0) THEN
      allocate (tracer_patch_data(num_tracers_with_forcing, max_nvar, numpatch))
      tracer_patch_data(:,:,:) = 0.0_r8
   ELSE
      ! Allocate a zero-sized array for processes with no patches (e.g. master)
      allocate (tracer_patch_data(num_tracers_with_forcing, max_nvar, 0))
   ENDIF

   ! Allocate forcmask_pch once (moved outside loop)
   IF (p_is_worker) THEN
      IF (numpatch > 0) THEN
         allocate (forcmask_pch(numpatch))
         forcmask_pch(:) = .true.
      ENDIF
   ENDIF

   ! Set initial date (moved outside loop)
   idate = (/ststamp%year, ststamp%day, ststamp%sec/)
   CALL adj2begin (idate)

      IF (p_is_io) THEN
         IF (allocated(forcn_tracer)) deallocate(forcn_tracer)
         IF (allocated(forcn_LB_tracer)) deallocate(forcn_LB_tracer)
         IF (allocated(forcn_UB_tracer)) deallocate(forcn_UB_tracer)
         IF (allocated(tracerdata)) deallocate(tracerdata)
         IF (allocated(traceravgcos)) deallocate(traceravgcos)
         
         allocate (forcn_tracer(num_tracers_with_forcing, max_nvar))
         allocate (forcn_LB_tracer(num_tracers_with_forcing, max_nvar)) 
         allocate (forcn_UB_tracer(num_tracers_with_forcing, max_nvar))
         allocate (tracerdata(num_tracers_with_forcing, max_nvar))
         allocate (traceravgcos(num_tracers_with_forcing, max_nvar))

      ENDIF

   DO i = 1, num_tracers_with_forcing
      ! Initialize user specified tracer forcing for this tracer (moved inside loop)
  !    CALL init_user_specified_tracer_forcing(i)
      
      ! Initialize grid for each tracer variable
      DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
         CALL tracer_read_latlon (DEF_Tracer_Forcings_NL(i)%tracer_dir, idate, i, j)
      ENDDO

      ! Allocate block data for this tracer on I/O processes
      IF (p_is_io) THEN
         DO ivar = 1, DEF_Tracer_Forcings_NL(i)%NVAR
            CALL allocate_block_data (gforc_tracers(i,ivar), forcn_tracer(i,ivar))
            CALL allocate_block_data (gforc_tracers(i,ivar), forcn_LB_tracer(i,ivar))
            CALL allocate_block_data (gforc_tracers(i,ivar), forcn_UB_tracer(i,ivar))
            ! Allocate memory for forcing data (using first grid for this tracer)
            CALL allocate_block_data (gforc_tracers(i,ivar), tracerdata(i,ivar))  ! forcing data
            CALL allocate_block_data (gforc_tracers(i,ivar), traceravgcos(i,ivar))  ! time-average of cos(zenith)
         ENDDO
      ENDIF 

      ! Process missing value information if needed
      ! Temporarily disabled for debugging
      IF (.FALSE.) THEN
         IF (p_is_master) WRITE(*,*) "    Processing missing value information for tracer ", i

         tstamp = idate
         CALL setstampLB_tracer(tstamp, 1, i, year, month, day, time_i)
         filename = trim(DEF_Tracer_Forcings_NL(i)%tracer_dir)//trim(tracerfilename(year, month, day, i, 1))
         tstamp_LB(i,1) = timestamp(-1, -1, -1)  ! Fixed: use 2D array indexing

         IF (p_is_master) THEN
            CALL ncio_get_attr (filename, DEF_Tracer_Forcings_NL(i)%vname(1), trim(DEF_Tracer_Forcings_NL(i)%missing_value_name), missing_value)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif
         ! Check that arrays are properly allocated before use
         IF (.NOT. allocated(tracerdata) .OR. .NOT. allocated(gforc_tracers)) THEN
            IF (p_is_master) WRITE(*,*) "ERROR: Arrays not allocated for missing value processing"
            CALL CoLM_stop()
         ENDIF
         
         ! Check array bounds
         IF (i > SIZE(tracerdata, 1) .OR. SIZE(tracerdata, 2) < DEF_Tracer_Forcings_NL(i)%NVAR) THEN
            IF (p_is_master) THEN
               WRITE(*,*) "ERROR: Array bounds mismatch in missing value processing"
               WRITE(*,*) "  Tracer index: ", i, " max: ", SIZE(tracerdata, 1)
               WRITE(*,*) "  Variable count: ", DEF_Tracer_Forcings_NL(i)%NVAR, " max: ", SIZE(tracerdata, 2)
            ENDIF
            CALL CoLM_stop()
         ENDIF
         
         DO ivar = 1, DEF_Tracer_Forcings_NL(i)%NVAR
            CALL ncio_read_block_time (filename, DEF_Tracer_Forcings_NL(i)%vname(ivar), gforc_tracers(i,ivar), time_i, tracerdata(i,ivar))
         ENDDO
      ENDIF

      ! Initialize spatial mapping for each tracer variable
      ! Build on all processes (following the pattern from MOD_Forcing.F90)
      DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
         IF (trim(DEF_Forcing_Interp_Method) == 'arealweight') THEN
            IF (present(lulcc_call)) THEN
               CALL mg2p_forc_tracers(i,j)%forc_free_mem
            ENDIF
            CALL mg2p_forc_tracers(i,j)%build_arealweighted (gforc_tracers(i,j), landpatch)
         ELSEIF (trim(DEF_Forcing_Interp_Method) == 'bilinear') THEN
            IF (present(lulcc_call)) THEN
               CALL mg2p_forc_tracers(i,j)%forc_free_mem
            ENDIF
            CALL mg2p_forc_tracers(i,j)%build_bilinear (gforc_tracers(i,j), landpatch)
         ELSE
            IF (p_is_master) WRITE(*,*) "      WARNING: Unknown interpolation method: ", TRIM(DEF_Forcing_Interp_Method)
         ENDIF
      ENDDO

      ! Set missing value for spatial mapping if needed
      ! Temporarily disabled for debugging
      IF (.FALSE.) THEN
         DO ivar = 1, DEF_Tracer_Forcings_NL(i)%NVAR
            CALL mg2p_forc_tracers(i,ivar)%set_missing_value (tracerdata(i,ivar), missing_value, forcmask_pch)
         ENDDO
      ENDIF

      ! Handle POINT dataset specific setup
      IF (trim(DEF_Tracer_Forcings_NL(i)%dataset_name) == 'POINT') THEN
         filename = trim(DEF_Tracer_Forcings_NL(i)%tracer_dir)//trim(DEF_Tracer_Forcings_NL(i)%fprefix(1))
      ENDIF

      IF (p_is_master) WRITE(*,*) "    --- Completed processing tracer ", i, " ---"
   ENDDO



   END SUBROUTINE tracer_forcing_init



   ! ---- forcing finalize ----
   SUBROUTINE tracer_forcing_final ()

   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

      IF (allocated(forcmask_pch)) deallocate(forcmask_pch)
      IF (allocated(forctime    )) deallocate(forctime    )
      IF (allocated(iforctime   )) deallocate(iforctime   )
      IF (allocated(forc_disk   )) deallocate(forc_disk   )
      IF (allocated(tstamp_LB   )) deallocate(tstamp_LB   )
      IF (allocated(tstamp_UB   )) deallocate(tstamp_UB   )
      IF (allocated(tracer_patch_data)) deallocate(tracer_patch_data)
      IF (allocated(gforc_tracers)) deallocate(gforc_tracers)
      IF (allocated(mg2p_forc_tracers)) deallocate(mg2p_forc_tracers)
      
      ! These arrays are only allocated on I/O processes
      IF (p_is_io) THEN
         IF (allocated(tracerdata)) deallocate(tracerdata)
         IF (allocated(traceravgcos)) deallocate(traceravgcos)
         IF (allocated(forcn_tracer)) deallocate(forcn_tracer)
         IF (allocated(forcn_LB_tracer)) deallocate(forcn_LB_tracer)
         IF (allocated(forcn_UB_tracer)) deallocate(forcn_UB_tracer)
      ENDIF

      ! IF (DEF_USE_Forcing_Downscaling) THEN
      !    IF (p_is_worker) THEN
      !       IF (numpatch > 0) THEN
      !
      !          deallocate (forc_topo_grid  )
      !          deallocate (forc_maxelv_grid)
      !
      !          deallocate (forc_t_grid     )
      !          deallocate (forc_th_grid    )
      !          deallocate (forc_q_grid     )
      ! #ifdef USE_TRACER
      !          deallocate (forc_q_grid_O18    )
      !          deallocate (forc_q_grid_H2    )
      ! #endif
      !          deallocate (forc_pbot_grid  )
      !          deallocate (forc_rho_grid   )
      !          deallocate (forc_prc_grid   )
      !          deallocate (forc_prl_grid   )
      !
      ! #ifdef USE_TRACER
      !          deallocate (forc_prc_grid_O18  )
      !          deallocate (forc_prl_grid_O18  )
      !          deallocate (forc_prc_grid_H2  )
      !          deallocate (forc_prl_grid_H2  )
      ! #endif
      !
      !          deallocate (forc_lwrad_grid )
      !          deallocate (forc_swrad_grid )
      !          deallocate (forc_hgt_grid   )
      !
      !          deallocate (forc_t_part     )
      !          deallocate (forc_th_part    )
      !          deallocate (forc_q_part     )
      !
      ! #ifdef USE_TRACER
      !          deallocate (forc_q_part_O18    )
      !          deallocate (forc_q_part_H2    )
      ! #endif
      !
      !          deallocate (forc_pbot_part  )
      !          deallocate (forc_rhoair_part)
      !          deallocate (forc_prc_part   )
      !          deallocate (forc_prl_part   )
      !
      ! #ifdef USE_TRACER
      !          deallocate (forc_prc_part_O18  )
      !          deallocate (forc_prl_part_O18  )
      !          deallocate (forc_prc_part_H2  )
      !          deallocate (forc_prl_part_H2  )
      ! #endif
      !
      !          deallocate (forc_frl_part   )
      !          deallocate (forc_swrad_part )
      !
      !       ENDIF
      !    ENDIF
      ! ENDIF

   END SUBROUTINE tracer_forcing_final

   ! ------------
   SUBROUTINE tracer_forcing_reset ()

   IMPLICIT NONE

      tstamp_LB(:,:) = timestamp(-1, -1, -1)
      tstamp_UB(:,:) = timestamp(-1, -1, -1)

   END SUBROUTINE tracer_forcing_reset


!-----------------------------------------------------------------------
   SUBROUTINE read_tracer_forcing (idate)
   USE MOD_OrbCosazi
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Const_Physical, only: rgas, grav
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables, only: alb
   USE MOD_Vars_1DForcing
   USE MOD_Vars_2DForcing
   USE MOD_Block
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandPatch
   USE MOD_RangeCheck
   USE MOD_ForcingDownscaling, only: rair, cpair, downscale_forcings, downscale_wind
   USE MOD_NetCDFVector

   IMPLICIT NONE

      integer, intent(in) :: idate(3)



   ! local variables:
   integer  :: ivar, istt, iend, id(3)
   integer  :: iblkme, ib, jb, i, j, ilon, ilat, np, ipart, ne
   integer  :: tracer_idx, var_idx
   real(r8) :: calday                                             ! Julian cal day (1.xx to 365.xx)
   real(r8) :: sunang, cloud, difrat, vnrat
   real(r8) :: a, hsolar, ratio_rvrf
   type(block_data_real8_2d) :: forc_xy_solarin
   integer  :: ii
   character(10) :: cyear = "2005"
   character(256):: lndname

   type(timestamp) :: mtstamp
   integer  :: dtLB, dtUB
   real(r8) :: cosz, coszen(numpatch), cosa, cosazi(numpatch), balb
   integer  :: year, month, mday
   logical  :: has_u,has_v
   real solar, frl, prcp, tm, us, vs, pres, qm
   real(r8) :: pco2m
   real(r8), dimension(12, numpatch) :: spaceship !NOTE: 12 is the dimension size of spaceship
   integer target_server, ierr

   ! Local variables
   integer :: num_tracers_with_forcing
   real(r8), allocatable :: temp_tracer_data(:)
   character(len=50) :: tracer_label
   type(block_data_real8_2d) :: dummy_block_data  ! Declare dummy for grid2pset calls
   
   ! Initialize num_tracers_with_forcing properly
   IF (allocated(DEF_Tracer_Forcings_NL)) THEN
      num_tracers_with_forcing = SIZE(DEF_Tracer_Forcings_NL, 1)
   ELSE
      num_tracers_with_forcing = 0
   ENDIF
 
   IF (num_tracers_with_forcing == 0 .OR. .NOT. allocated(DEF_Tracers)) THEN
      IF (p_is_master) WRITE(*,*) "No tracer forcing data to process, returning..."
      RETURN
   ENDIF
   
   ! Check if DEF_Tracer_Forcings_NL is properly allocated
   IF (.NOT. allocated(DEF_Tracer_Forcings_NL)) THEN
      IF (p_is_master) WRITE(*,*) "ERROR: DEF_Tracer_Forcings_NL not allocated"
      RETURN
   ENDIF
   
   ! Verify max_nvar is consistent with current configuration
   ! (max_nvar should have been set in tracer_forcing_init)
   IF (max_nvar <= 0) THEN
      IF (p_is_master) THEN
         WRITE(*,*) "ERROR: max_nvar not properly initialized (max_nvar = ", max_nvar, ")"
         WRITE(*,*) "Tracer forcing configuration:"
         DO tracer_idx = 1, num_tracers_with_forcing
            WRITE(*,*) "  Tracer ", tracer_idx, ": ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name), &
                       " NVAR = ", DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
         ENDDO
      ENDIF
      RETURN
   ENDIF


   ! Loop over each tracer
   DO tracer_idx = 1, num_tracers_with_forcing

      !------------------------------------------------------------
      ! READ in THE TRACER FORCING
      ! read lower and upper boundary forcing data for this tracer
      IF (p_is_io) THEN
         CALL metreadLBUB_tracer(idate, TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir), tracer_idx)
         ! set model time stamp
         id(:) = idate(:)
         mtstamp = id
         
         ! loop for variables of this tracer
         DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
            IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) == 'NULL') CYCLE     ! no data, CYCLE
            IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx)) == 'NULL') CYCLE
            ! to make sure the forcing data calculated is in the range of time
            ! interval [LB, UB]. If not, read new data.
            
            ! First check if time boundaries are properly initialized
            IF (tstamp_LB(tracer_idx,var_idx)%year == -1 .or. tstamp_UB(tracer_idx,var_idx)%year == -1) THEN
               write(6, *) "WARNING: Tracer time boundaries not initialized, re-reading..."
               CALL metreadLBUB_tracer(idate, TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir), tracer_idx)
            ENDIF
            
            IF ( (mtstamp < tstamp_LB(tracer_idx,var_idx)) .or. (tstamp_UB(tracer_idx,var_idx) < mtstamp) ) THEN
               write(6, *) "========== TRACER FORCING TIME RANGE ERROR =========="
               write(6, *) "Current time stamp: ", mtstamp%year, mtstamp%day, mtstamp%sec
               write(6, *) "Tracer index: ", tracer_idx, " Variable index: ", var_idx
               write(6, *) "Lower bound: ", tstamp_LB(tracer_idx,var_idx)%year, tstamp_LB(tracer_idx,var_idx)%day, tstamp_LB(tracer_idx,var_idx)%sec
               write(6, *) "Upper bound: ", tstamp_UB(tracer_idx,var_idx)%year, tstamp_UB(tracer_idx,var_idx)%day, tstamp_UB(tracer_idx,var_idx)%sec
               write(6, *) "Variable name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx))
               write(6, *) "Tracer name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name)
               write(6, *) "Trying to reinitialize time boundaries..."
               
               ! Try to reinitialize boundaries
               tstamp_LB(tracer_idx,var_idx) = timestamp(-1, -1, -1)
               tstamp_UB(tracer_idx,var_idx) = timestamp(-1, -1, -1)
               
               ! Call the boundary reading function again
               CALL metreadLBUB_tracer(idate, TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir), tracer_idx)
               
               ! Check again
               IF ( (mtstamp < tstamp_LB(tracer_idx,var_idx)) .or. (tstamp_UB(tracer_idx,var_idx) < mtstamp) ) THEN
                  write(6, *) "After reinitializing boundaries:"
                  write(6, *) "Lower bound: ", tstamp_LB(tracer_idx,var_idx)%year, tstamp_LB(tracer_idx,var_idx)%day, tstamp_LB(tracer_idx,var_idx)%sec
                  write(6, *) "Upper bound: ", tstamp_UB(tracer_idx,var_idx)%year, tstamp_UB(tracer_idx,var_idx)%day, tstamp_UB(tracer_idx,var_idx)%sec
                  write(6, *) "the tracer forcing data required is still out of range! STOP!"
                  CALL CoLM_stop()
               ELSE
                  write(6, *) "Successfully fixed time range boundaries!"
               ENDIF
            ENDIF

            ! calculate distance to lower/upper boundary
            dtLB = mtstamp - tstamp_LB(tracer_idx,var_idx)
            dtUB = tstamp_UB(tracer_idx,var_idx) - mtstamp

            ! linear method, for most variables
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'linear') THEN
               IF ( (dtLB+dtUB) > 0 ) THEN
                  CALL block_data_linear_interp ( &
                     forcn_LB_tracer(tracer_idx,var_idx), real(dtUB,r8)/real(dtLB+dtUB,r8), &
                     forcn_UB_tracer(tracer_idx,var_idx), real(dtLB,r8)/real(dtLB+dtUB,r8), &
                     forcn_tracer(tracer_idx,var_idx))
               ELSE
                  CALL block_data_copy (forcn_LB_tracer(tracer_idx,var_idx), forcn_tracer(tracer_idx,var_idx))
               ENDIF
            ENDIF

            ! nearest method, for precipitation
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'nearest') THEN
               IF (dtLB <= dtUB) THEN
                  CALL block_data_copy (forcn_LB_tracer(tracer_idx,var_idx), forcn_tracer(tracer_idx,var_idx))
               ELSE
                  CALL block_data_copy (forcn_UB_tracer(tracer_idx,var_idx), forcn_tracer(tracer_idx,var_idx))
               ENDIF
            ENDIF

            ! uniform method
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'uniform') THEN
               IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%timelog(var_idx)) == 'forward') THEN
                  CALL block_data_copy (forcn_LB_tracer(tracer_idx,var_idx), forcn_tracer(tracer_idx,var_idx))
               ELSE
                  CALL block_data_copy (forcn_UB_tracer(tracer_idx,var_idx), forcn_tracer(tracer_idx,var_idx))
               ENDIF
            ENDIF

            ! coszen method, for SW radiation
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'coszen') THEN
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  DO j = 1, gforc_tracers(tracer_idx,var_idx)%ycnt(jb)
                     DO i = 1, gforc_tracers(tracer_idx,var_idx)%xcnt(ib)

                        ilat = gforc_tracers(tracer_idx,var_idx)%ydsp(jb) + j
                        ilon = gforc_tracers(tracer_idx,var_idx)%xdsp(ib) + i
                        IF (ilon > gforc_tracers(tracer_idx,var_idx)%nlon) ilon = ilon - gforc_tracers(tracer_idx,var_idx)%nlon

                        calday = calendarday(mtstamp)
                        cosz = orb_coszen(calday, gforc_tracers(tracer_idx,var_idx)%rlon(ilon), gforc_tracers(tracer_idx,var_idx)%rlat(ilat))
                        cosz = max(0.001, cosz)
                        ! 10/24/2024, yuan: deal with time log with backward or forward
                        IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%timelog(var_idx)) == 'forward') THEN
                           forcn_tracer(tracer_idx,var_idx)%blk(ib,jb)%val(i,j) = &
                              cosz / traceravgcos(tracer_idx,var_idx)%blk(ib,jb)%val(i,j) * forcn_LB_tracer(tracer_idx,var_idx)%blk(ib,jb)%val(i,j)
                        ELSE
                           forcn_tracer(tracer_idx,var_idx)%blk(ib,jb)%val(i,j) = &
                              cosz / traceravgcos(tracer_idx,var_idx)%blk(ib,jb)%val(i,j) * forcn_UB_tracer(tracer_idx,var_idx)%blk(ib,jb)%val(i,j)
                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ENDDO

      ENDIF

#ifdef USEMPI
      ! Ensure all processes complete I/O before proceeding to mapping
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! Map tracer forcing data to patches using the appropriate spatial mapping
      ! This must be called on ALL processes. Internally, grid2pset handles
      ! sending from I/O processes and receiving on worker processes.
      DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
         ! Check array bounds before accessing mg2p_forc_tracers
         IF (tracer_idx <= SIZE(mg2p_forc_tracers, 1) .AND. var_idx <= SIZE(mg2p_forc_tracers, 2)) THEN
            ! Skip NULL variables
            IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) == 'NULL') CYCLE
            
#ifdef USEMPI
            ! Ensure all processes are synchronized before spatial mapping
            CALL mpi_barrier (p_comm_glb, p_err)
#endif
            ! Call spatial mapping function - this handles MPI communication automatically
            ! Note: Only I/O processes have forcn_tracer data, but all processes must call grid2pset
            IF (p_is_io .AND. allocated(forcn_tracer)) THEN
               ! I/O processes send the data
               CALL mg2p_forc_tracers(tracer_idx,var_idx)%grid2pset (forcn_tracer(tracer_idx,var_idx), tracer_patch_data(tracer_idx,var_idx,:))
            ELSE
               ! Non-I/O processes receive the data (using dummy source)
               CALL mg2p_forc_tracers(tracer_idx,var_idx)%grid2pset (dummy_block_data, tracer_patch_data(tracer_idx,var_idx,:))
            ENDIF
            
#ifdef USEMPI
            ! Ensure all processes complete mapping before proceeding
            CALL mpi_barrier (p_comm_glb, p_err)
#endif
         ENDIF
      ENDDO

   ENDDO  ! End tracer loop

   IF (p_is_master) THEN
      WRITE(*,*) "=== READ_TRACER_FORCING COMPLETED SUCCESSFULLY ==="
   ENDIF
   

#ifdef RangeCheck
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) write(*,'(/, A25)') 'Checking tracer forcing ...'

      ! Check tracer patch data on worker processes, and receive messages on master
      IF (p_is_worker) THEN
          IF (allocated(tracer_patch_data) .AND. numpatch > 0) THEN
             DO tracer_idx = 1, num_tracers_with_forcing
                DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
                   IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) /= 'NULL') THEN
                      ! Allocate temporary array for range checking
                      IF (allocated(temp_tracer_data)) deallocate(temp_tracer_data)
                      allocate(temp_tracer_data(numpatch))
                      temp_tracer_data(:) = tracer_patch_data(tracer_idx,var_idx,:)
                      
                      ! Create a fixed-length descriptive label for the tracer variable.
                      tracer_label = repeat(' ', 25) ! Initialize with spaces
                      write(tracer_label, '(A,1X,A,1X,A)') 'Tracer', &
                         trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name), &
                         trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx))
                      
                      CALL check_vector_data(tracer_label, temp_tracer_data)
                      
                      deallocate(temp_tracer_data)
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
      ELSE IF (p_is_master) THEN
         ! Master process must call check_vector_data to receive MPI messages
         ! Allocate dummy array for master
         IF (.NOT. allocated(temp_tracer_data)) THEN
            allocate(temp_tracer_data(1))  ! Master only needs a dummy array
            temp_tracer_data(1) = 0.0_r8
         ENDIF
         
         DO tracer_idx = 1, num_tracers_with_forcing
            DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
               IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) /= 'NULL') THEN
                  ! Master calls with dummy data; it only performs MPI_Recv
                  CALL check_vector_data('dummy', temp_tracer_data)
               ENDIF
            ENDDO
         ENDDO
         
         ! Clean up dummy array
         IF (allocated(temp_tracer_data)) deallocate(temp_tracer_data)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif

   END SUBROUTINE read_tracer_forcing


   SUBROUTINE metreadLBUB_tracer (idate, dir_forcing, tracer_idx)

   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_Block
   USE MOD_NetCDFBlock
   USE MOD_RangeCheck
   IMPLICIT NONE

   integer, intent(in) :: idate(3)
   character(len=*), intent(in) :: dir_forcing
   integer, intent(in) :: tracer_idx

   ! Local variables
   integer         :: ivar, year, month, day, time_i
   integer         :: iblkme, ib, jb, i, j
   type(timestamp) :: mtstamp
   character(len=256) :: filename
   character(len=256) :: dataset_name_lower


      ! Arrays should already be allocated in tracer_forcing_init on I/O processes

      mtstamp = idate

      DO ivar = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR

         IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)) == 'NULL') THEN
            CYCLE     ! no data, CYCLE
         ENDIF



         ! Check array bounds before accessing
         IF (tracer_idx > SIZE(tstamp_LB, 1) .OR. ivar > SIZE(tstamp_LB, 2) .OR. &
             tracer_idx > SIZE(tstamp_UB, 1) .OR. ivar > SIZE(tstamp_UB, 2) .OR. &
             tracer_idx > SIZE(forcn_LB_tracer, 1) .OR. ivar > SIZE(forcn_LB_tracer, 2) .OR. &
             tracer_idx > SIZE(forcn_UB_tracer, 1) .OR. ivar > SIZE(forcn_UB_tracer, 2)) THEN
            CYCLE
         ENDIF


         ! lower and upper boundary data already exist, CYCLE
         IF ( .not.(tstamp_LB(tracer_idx,ivar)%year == -1) .and. .not.(tstamp_UB(tracer_idx,ivar)%year == -1) .and. &
            tstamp_LB(tracer_idx,ivar)%year<=mtstamp%year .and. mtstamp%year<tstamp_UB(tracer_idx,ivar)%year ) THEN
            CYCLE
         ENDIF

         ! set lower boundary time stamp and get data
         IF (tstamp_LB(tracer_idx,ivar)%year == -1) THEN
            IF (p_is_master) THEN
               write(6, *) "=== SETTING LB TIMESTAMP ==="
               write(6, *) "Tracer:", tracer_idx, "Variable:", ivar
               write(6, *) "Current mtstamp:", mtstamp%year, mtstamp%day, mtstamp%sec
            ENDIF
            CALL setstampLB_tracer(mtstamp, ivar, tracer_idx, year, month, day, time_i)
            IF (p_is_master) THEN
               write(6, *) "After setting LB:"
               write(6, *) "LB timestamp:", tstamp_LB(tracer_idx,ivar)%year, tstamp_LB(tracer_idx,ivar)%day, tstamp_LB(tracer_idx,ivar)%sec
               write(6, *) "File year/month/day/time_i:", year, month, day, time_i
            ENDIF

            filename = trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir)//trim(tracerfilename(year, month, day, tracer_idx, ivar))
            ! Convert dataset name to lowercase for case-insensitive comparison
            dataset_name_lower = adjustl(trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name))
            CALL to_lower(dataset_name_lower)
     
            
            IF (trim(dataset_name_lower) == 'point') THEN
               IF (forcing_read_ahead) THEN
                  tracerdata(tracer_idx,ivar)%blk(gblock%xblkme(1),gblock%yblkme(1))%val = forc_disk(time_i,ivar)
               ELSE
#ifndef URBAN_MODEL
                  CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata(tracer_idx,ivar))
#else
                  IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)) == 'Rainf') THEN
                     CALL ncio_read_site_time (filename, 'Rainf', time_i, rainf)
                     CALL ncio_read_site_time (filename, 'Snowf', time_i, snowf)

                     DO iblkme = 1, gblock%nblkme
                        ib = gblock%xblkme(iblkme)
                        jb = gblock%yblkme(iblkme)

                        tracerdata(tracer_idx,ivar)%blk(ib,jb)%val(1,1) = rainf%blk(ib,jb)%val(1,1) + snowf%blk(ib,jb)%val(1,1)
                     ENDDO
                  ELSE
                     CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata(tracer_idx,ivar))
                  ENDIF
#endif
               ENDIF
            ELSE

               ! Clear any previous allocation to ensure fresh read
               IF (allocated(tracerdata(tracer_idx,ivar)%blk)) THEN
                  DO ib = 1, SIZE(tracerdata(tracer_idx,ivar)%blk, 1)
                     DO jb = 1, SIZE(tracerdata(tracer_idx,ivar)%blk, 2)
                        IF (allocated(tracerdata(tracer_idx,ivar)%blk(ib,jb)%val)) THEN
                           deallocate(tracerdata(tracer_idx,ivar)%blk(ib,jb)%val)
                        ENDIF
                     ENDDO
                  ENDDO
                  deallocate(tracerdata(tracer_idx,ivar)%blk)
               ENDIF
               
               ! Reallocate block data structure
               CALL allocate_block_data (gforc_tracers(tracer_idx,ivar), tracerdata(tracer_idx,ivar))
               
               ! Perform the NetCDF read with explicit error checking
               CALL ncio_read_block_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), gforc_tracers(tracer_idx,ivar), time_i, tracerdata(tracer_idx,ivar))
               
               ! Explicit error checking after NetCDF read
               IF (.NOT. allocated(tracerdata(tracer_idx,ivar)%blk)) THEN
                  ! Re-allocate and fill with zeros as fallback
                  CALL allocate_block_data (gforc_tracers(tracer_idx,ivar), tracerdata(tracer_idx,ivar))
                  CALL flush_block_data (tracerdata(tracer_idx,ivar), 0.0_r8)
               ENDIF
            ENDIF

            CALL block_data_copy (tracerdata(tracer_idx,ivar), forcn_LB_tracer(tracer_idx,ivar))
            
         ENDIF

         ! set upper boundary time stamp and get data
         IF (tstamp_UB(tracer_idx,ivar)%year == -1) THEN
            IF (p_is_master) THEN
               write(6, *) "=== SETTING UB TIMESTAMP ==="
               write(6, *) "Tracer:", tracer_idx, "Variable:", ivar
            ENDIF
            CALL setstampUB_tracer(ivar, tracer_idx, year, month, day, time_i)
            IF (p_is_master) THEN
               write(6, *) "After setting UB:"
               write(6, *) "UB timestamp:", tstamp_UB(tracer_idx,ivar)%year, tstamp_UB(tracer_idx,ivar)%day, tstamp_UB(tracer_idx,ivar)%sec
               write(6, *) "File year/month/day/time_i:", year, month, day, time_i
            ENDIF

            IF (year <= DEF_Tracer_Forcings_NL(tracer_idx)%endyr) THEN
               ! read forcing data
               filename = trim(dir_forcing)//trim(tracerfilename(year, month, day,tracer_idx, ivar))


               
   
               ! Convert dataset name to lowercase for case-insensitive comparison
               dataset_name_lower = adjustl(trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name))
               CALL to_lower(dataset_name_lower)
               
               IF (trim(dataset_name_lower) == 'point') THEN

                  IF (forcing_read_ahead) THEN
                     tracerdata(tracer_idx,ivar)%blk(gblock%xblkme(1),gblock%yblkme(1))%val = forc_disk(time_i,ivar)
                  ELSE
#ifndef URBAN_MODEL
                     ! IF (p_is_master) WRITE(*,*) "    Calling ncio_read_site_time for UB variable: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar))
                     CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata(tracer_idx,ivar))
#else
                     IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)) == 'Rainf') THEN
                        CALL ncio_read_site_time (filename, 'Rainf', time_i, rainf)
                        CALL ncio_read_site_time (filename, 'Snowf', time_i, snowf)

                        DO iblkme = 1, gblock%nblkme
                           ib = gblock%xblkme(iblkme)
                           jb = gblock%yblkme(iblkme)

                           tracerdata(tracer_idx,ivar)%blk(ib,jb)%val(1,1) = rainf%blk(ib,jb)%val(1,1) + snowf%blk(ib,jb)%val(1,1)
                        ENDDO
                     ELSE
                        CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata(tracer_idx,ivar))
                     ENDIF
#endif
                  ENDIF
               ELSE

                  
                  ! Clear any previous allocation to ensure fresh read
                  IF (allocated(tracerdata(tracer_idx,ivar)%blk)) THEN
                     DO ib = 1, SIZE(tracerdata(tracer_idx,ivar)%blk, 1)
                        DO jb = 1, SIZE(tracerdata(tracer_idx,ivar)%blk, 2)
                           IF (allocated(tracerdata(tracer_idx,ivar)%blk(ib,jb)%val)) THEN
                              deallocate(tracerdata(tracer_idx,ivar)%blk(ib,jb)%val)
                           ENDIF
                        ENDDO
                     ENDDO
                     deallocate(tracerdata(tracer_idx,ivar)%blk)
                  ENDIF
                  
                  ! Reallocate block data structure
                  CALL allocate_block_data (gforc_tracers(tracer_idx,ivar), tracerdata(tracer_idx,ivar))
                  
                  ! Perform the UB NetCDF read with explicit error checking
                  CALL ncio_read_block_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), gforc_tracers(tracer_idx,ivar), time_i, tracerdata(tracer_idx,ivar))
                  
                  ! Explicit error checking after UB NetCDF read
                  IF (.NOT. allocated(tracerdata(tracer_idx,ivar)%blk)) THEN
                     ! Re-allocate and fill with zeros as fallback
                     CALL allocate_block_data (gforc_tracers(tracer_idx,ivar), tracerdata(tracer_idx,ivar))
                     CALL flush_block_data (tracerdata(tracer_idx,ivar), 0.0_r8)
                  ENDIF
               ENDIF

               CALL block_data_copy (tracerdata(tracer_idx,ivar), forcn_UB_tracer(tracer_idx,ivar))
               
            ELSE
               write(*,*) year, DEF_Tracer_Forcings_NL(tracer_idx)%endyr
               print *, 'NOTE: reaching the END of forcing data, always reuse the last time step data!'
            ENDIF

            
            IF (TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(ivar)) == 'coszen') THEN
               ! IF (p_is_master) WRITE(*,*) "    Calculating avgcos for tracer ", tracer_idx, " var ", ivar
               CALL calavgcos(idate, tracer_idx, ivar)
            ENDIF
         ENDIF

      ENDDO
   END SUBROUTINE metreadLBUB_tracer

!-----------------------------------------------------------------------
   SUBROUTINE tracer_read_latlon (dir_forcing, idate, ivar, jvar)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Namelist
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_forcing
   integer, intent(in) :: idate(3)
   integer, intent(in) :: ivar, jvar

   ! Local variables
   character(len=256) :: filename
   integer         :: year, month, day, time_i
   type(timestamp) :: mtstamp
   real(r8), allocatable :: latxy (:,:)    ! latitude values in 2d
   real(r8), allocatable :: lonxy (:,:)    ! longitude values in 2d
   real(r8), allocatable :: lon_in(:)
   real(r8), allocatable :: lat_in(:)
   character(len=256) :: dataset_name_lower

      ! Convert dataset name to lowercase for case-insensitive comparison
      dataset_name_lower = trim(DEF_Tracer_Forcings_NL(ivar)%dataset_name)
      CALL to_lower(dataset_name_lower)

      IF (trim(dataset_name_lower) == 'point' .or. trim(dataset_name_lower) == 'cpl7' ) THEN
         CALL gforc_tracers(ivar,jvar)%define_by_ndims (360, 180)
      ELSE
         ! Follow the exact same grid initialization logic as main forcing system
         mtstamp = idate
         CALL setstampLB_tracer(mtstamp, 1, ivar, year, month, day, time_i)
         filename = trim(DEF_Tracer_Forcings_NL(ivar)%tracer_dir)//trim(tracerfilename(year, month, day, ivar, 1))
         tstamp_LB(ivar,1) = timestamp(-1, -1, -1)

         ! Read actual grid from NetCDF file (same as main forcing system)
         if (DEF_Tracer_Forcings_NL(ivar)%dim2d) then
            CALL ncio_read_serial (filename, (DEF_Tracer_Forcings_NL(ivar)%latname), latxy)
            CALL ncio_read_serial (filename, (DEF_Tracer_Forcings_NL(ivar)%lonname), lonxy)
            allocate (lat_in (size(latxy,2)))
            allocate (lon_in (size(lonxy,1)))
            lat_in = latxy(1,:)
            lon_in = lonxy(:,1)
            deallocate (latxy)
            deallocate (lonxy)
         ELSE
            CALL ncio_read_serial (filename, (DEF_Tracer_Forcings_NL(ivar)%latname), lat_in)
            CALL ncio_read_serial (filename, (DEF_Tracer_Forcings_NL(ivar)%lonname), lon_in)
         ENDIF

         IF (.not. DEF_Tracer_Forcings_NL(ivar)%regional) THEN
            CALL gforc_tracers(ivar,jvar)%define_by_center (lat_in, lon_in)
         ELSE
            CALL gforc_tracers(ivar,jvar)%define_by_center (lat_in, lon_in, &
               south = DEF_forcing%regbnd(1), north = DEF_forcing%regbnd(2), &
               west  = DEF_forcing%regbnd(3), east  = DEF_forcing%regbnd(4))
         ENDIF

         IF (allocated(lat_in)) deallocate(lat_in)
         IF (allocated(lon_in)) deallocate(lon_in)

      ENDIF
      CALL gforc_tracers(ivar,jvar)%set_rlon ()
      CALL gforc_tracers(ivar,jvar)%set_rlat ()

   END SUBROUTINE tracer_read_latlon



!-----------------------------------------------------------------------
! !DESCRIPTION:
!    set the lower boundary time stamp and record information,
!    a KEY FUNCTION of this MODULE
!
! - for time stamp, set it regularly as the model time step.
! - for record information, account for:
!    o year alternation
!    o month alternation
!    o leap year
!    o required data just beyond the first record
!
! !REVISIONS:
!  04/2014, Hua Yuan: initial code
!
!-----------------------------------------------------------------------
   SUBROUTINE setstampLB_tracer(mtstamp, var_i, tracer_idx, year, month, mday, time_i)

   IMPLICIT NONE
   type(timestamp), intent(in)  :: mtstamp
   integer,         intent(in)  :: var_i
   integer,         intent(in)  :: tracer_idx
   integer,         intent(out) :: year
   integer,         intent(out) :: month
   integer,         intent(out) :: mday
   integer,         intent(out) :: time_i

   integer :: i, day, sec, ntime
   integer :: months(0:12)
   character(len=256) :: dataset_name_lower



      year = mtstamp%year
      day  = mtstamp%day
      sec  = mtstamp%sec

      ! Convert dataset name to lowercase for case-insensitive comparison
      dataset_name_lower = trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
      CALL to_lower(dataset_name_lower)

      IF (trim(dataset_name_lower) == 'point') THEN

         ! For POINT data, we need to handle this differently
         ! Since forctime is not initialized, we'll use a simplified approach
         tstamp_LB(tracer_idx,var_i)%year = year
         tstamp_LB(tracer_idx,var_i)%day  = day
         tstamp_LB(tracer_idx,var_i)%sec  = sec
         time_i = 1



         RETURN
      ENDIF

      tstamp_LB(tracer_idx,var_i)%year = year
      tstamp_LB(tracer_idx,var_i)%day  = day


      ! in the case of one year one file
      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'year' ) THEN

         ! calculate the initial second
         sec    = 86400*(day-1) + sec
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
         sec    = (time_i-1)*DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i) - 86400*(day-1)
         tstamp_LB(tracer_idx,var_i)%sec = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            tstamp_LB(tracer_idx,var_i)%sec = 86400 + sec
            tstamp_LB(tracer_idx,var_i)%day = day - 1
            IF (tstamp_LB(tracer_idx,var_i)%day == 0) THEN
               tstamp_LB(tracer_idx,var_i)%year = year - 1
               IF ( isleapyear(tstamp_LB(tracer_idx,var_i)%year) ) THEN
                  tstamp_LB(tracer_idx,var_i)%day = 366
               ELSE
                  tstamp_LB(tracer_idx,var_i)%day = 365
               ENDIF
            ENDIF
         ENDIF

         ! set record info (year, time_i)
         IF ( sec<0 .or. (sec==0 .and. DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i).NE.0) ) THEN

            ! IF the required data just behind the first record
            ! -> set to the first record
            IF ( year==DEF_Tracer_Forcings_NL(tracer_idx)%startyr .and. month==DEF_Tracer_Forcings_NL(tracer_idx)%startmo .and. day==1 ) THEN
               sec = DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)

               ! ELSE, set to one record backward
            ELSE
               sec = 86400 + sec
               day = day - 1
               IF (day == 0) THEN
                  year = year - 1
                  IF ( isleapyear(year) ) THEN
                     day = 366
                  ELSE
                     day = 365
                  ENDIF
               ENDIF
            ENDIF
         ENDIF ! ENDIF (sec <= 0)

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before after FEB 28th (Julian day 59).
         IF ( .not. DEF_Tracer_Forcings_NL(tracer_idx)%leapyear .and. isleapyear(year) .and. day>59 ) THEN
            day = day - 1
         ENDIF

         ! get record time index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
         ! Ensure month and mday are initialized for year grouping
         IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'year') THEN
             CALL julian2monthday(year, day, month, mday)
         ENDIF
      ENDIF

      ! in the case of one month one file
      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'month' ) THEN

         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)


         ! calculate initial second value
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
         sec    = (time_i-1)*DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i) - 86400*(mday-1)
         tstamp_LB(tracer_idx,var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            tstamp_LB(tracer_idx,var_i)%sec = 86400 + sec
            tstamp_LB(tracer_idx,var_i)%day = day - 1
            IF (tstamp_LB(tracer_idx,var_i)%day == 0) THEN
               tstamp_LB(tracer_idx,var_i)%year = year - 1
               IF ( isleapyear(tstamp_LB(tracer_idx,var_i)%year) ) THEN
                  tstamp_LB(tracer_idx,var_i)%day = 366
               ELSE
                  tstamp_LB(tracer_idx,var_i)%day = 365
               ENDIF
            ENDIF
         ENDIF

         ! set record info (year, month, time_i)
         IF ( sec<0 .or. (sec==0 .and. DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i).ne.0) ) THEN

            ! IF just behind the first record -> set to first record
            IF ( year==DEF_Tracer_Forcings_NL(tracer_idx)%startyr .and. month==DEF_Tracer_Forcings_NL(tracer_idx)%startmo .and. mday==1 ) THEN
               sec = DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)

               ! set to one record backward
            ELSE
               sec = 86400 + sec
               mday = mday - 1
               IF (mday == 0) THEN
                  month = month - 1
                  ! bug found by Zhu Siguang & Zhang Xiangxiang, 05/19/2014
                  ! move the below line in the 'ELSE' statement
                  !mday = months(month) - months(month-1)
                  IF (month == 0) THEN
                     month = 12
                     year = year - 1
                     mday = 31
                  ELSE
                     mday = months(month) - months(month-1)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before, i.e., FEB 28th.
         IF ( .not. DEF_Tracer_Forcings_NL(tracer_idx)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! get record time index
         sec = 86400*(mday-1) + sec
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
  
      ENDIF

      ! in the case of one day one file
      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'day' ) THEN

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)


         ! calculate initial second value
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
         sec    = (time_i-1)*DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)
         tstamp_LB(tracer_idx,var_i)%sec  = sec


         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            tstamp_LB(tracer_idx,var_i)%sec = 86400 + sec
            tstamp_LB(tracer_idx,var_i)%day = day - 1
            IF (tstamp_LB(tracer_idx,var_i)%day == 0) THEN
               tstamp_LB(tracer_idx,var_i)%year = year - 1
               IF ( isleapyear(tstamp_LB(tracer_idx,var_i)%year) ) THEN
                  tstamp_LB(tracer_idx,var_i)%day = 366
               ELSE
                  tstamp_LB(tracer_idx,var_i)%day = 365
               ENDIF
            ENDIF

            IF ( year==DEF_Tracer_Forcings_NL(tracer_idx)%startyr .and. month==DEF_Tracer_Forcings_NL(tracer_idx)%startmo .and. mday==1 ) THEN
               sec = DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)
            ! set to one record backward
            ELSE
               sec = 86400 + sec
               year = tstamp_LB(tracer_idx,var_i)%year
               CALL julian2monthday(tstamp_LB(tracer_idx,var_i)%year, tstamp_LB(tracer_idx,var_i)%day, month, mday)
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before, i.e., FEB 28th.
         IF ( .not. DEF_Tracer_Forcings_NL(tracer_idx)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! get record time index
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
      ENDIF

      IF (time_i <= 0) THEN
         IF (p_is_master) WRITE(*,*) "      ERROR: Invalid time_i = ", time_i
         write(6, *) "got the wrong time record of forcing! STOP!"; CALL CoLM_stop()
      ENDIF

      RETURN

   END SUBROUTINE setstampLB_tracer

!-----------------------------------------------------------------------
! !DESCRIPTION:
!    set the upper boundary time stamp and record information,
!    a KEY FUNCTION of this MODULE
!
! !REVISIONS:
!  04/2014, Hua Yuan: initial code
!
!-----------------------------------------------------------------------
   SUBROUTINE setstampUB_tracer(var_i, tracer_idx, year, month, mday, time_i)

   IMPLICIT NONE
   integer,         intent(in)  :: var_i
   integer,         intent(in)  :: tracer_idx
   integer,         intent(out) :: year
   integer,         intent(out) :: month
   integer,         intent(out) :: mday
   integer,         intent(out) :: time_i

   integer :: day, sec
   integer :: months(0:12)
   character(len=256) :: dataset_name_lower


      ! Convert dataset name to lowercase for case-insensitive comparison
      dataset_name_lower = adjustl(trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name))
      CALL to_lower(dataset_name_lower)

      IF (trim(dataset_name_lower) == 'point') THEN         
         IF ( tstamp_UB(tracer_idx,var_i)%year == -1 ) THEN
            tstamp_UB(tracer_idx,var_i) = tstamp_LB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
         ELSE
            tstamp_UB(tracer_idx,var_i) = tstamp_UB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
         ENDIF

         time_i = 1  ! Default time index for POINT data
         year = tstamp_UB(tracer_idx,var_i)%year
         RETURN
      ENDIF

      ! calculate the time stamp
      IF ( tstamp_UB(tracer_idx,var_i)%year == -1 ) THEN
         tstamp_UB(tracer_idx,var_i) = tstamp_LB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
      ELSE
         tstamp_LB(tracer_idx,var_i) = tstamp_UB(tracer_idx,var_i)
         tstamp_UB(tracer_idx,var_i) = tstamp_UB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
      ENDIF

      ! calculate initial year, day, and second values
      year = tstamp_UB(tracer_idx,var_i)%year
      day  = tstamp_UB(tracer_idx,var_i)%day
      sec  = tstamp_UB(tracer_idx,var_i)%sec

 
      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'year' ) THEN
         ! adjust year value
         IF ( sec==86400 .and. DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i).eq.0 ) THEN
            sec = 0
            day = day + 1
            IF( isleapyear(year) .and. day==367) THEN
               year = year + 1; day = 1
            ENDIF
            IF( .not. isleapyear(year) .and. day==366) THEN
               year = year + 1; day = 1
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before after FEB 28th (Julian day 59).
         IF ( .not. DEF_Tracer_Forcings_NL(tracer_idx)%leapyear .and. isleapyear(year) .and. day>59 ) THEN
            day = day - 1
         ENDIF

         ! set record index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
      ENDIF

      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'month' ) THEN

         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! record in the next day, adjust year, month and second values
         IF ( sec==86400 .and. DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i).eq.0 ) THEN
            sec  = 0
            mday = mday + 1
            IF ( mday > (months(month)-months(month-1)) ) THEN
               mday = 1
               ! bug found by Zhu Siguang, 05/25/2014
               ! move the below line in the 'ELSE' statement
               !month = month + 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! for day 29th Feb, USE the data 1 day before, i.e., 28th FEB.
         IF ( .not. DEF_Tracer_Forcings_NL(tracer_idx)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! set record index
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
      ENDIF

      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'day' ) THEN         
         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)
         !mday = day

         ! record in the next day, adjust year, month and second values
         IF ( sec==86400 .and. DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i).eq.0 ) THEN
            sec  = 0
            mday = mday + 1
            IF ( mday > (months(month)-months(month-1)) ) THEN
               mday = 1
               ! bug found by Zhu Siguang, 05/25/2014
               ! move the below line in the 'ELSE' statement
               !month = month + 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! for day 29th Feb, USE the data 1 day before, i.e., 28th FEB.
         IF ( .not. DEF_Tracer_Forcings_NL(tracer_idx)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! set record index
         time_i = floor( (sec-DEF_Tracer_Forcings_NL(tracer_idx)%offset(var_i)) *1. / DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i) ) + 1
      ENDIF

      IF (time_i < 0) THEN
         IF (p_is_master) WRITE(*,*) "      ERROR: Invalid time_i = ", time_i
         write(6, *) "got the wrong time record of forcing! STOP!"; CALL CoLM_stop()
      ENDIF


      RETURN

   END SUBROUTINE setstampUB_tracer

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  calculate time average coszen value between [LB, UB]
!
! !REVISIONS:
!  04/2014, Hua Yuan: this method is adapted from CLM
!
!-----------------------------------------------------------------------
   SUBROUTINE calavgcos(idate, tracer_idx, var_idx)

   USE MOD_Block
   USE MOD_DataType
   IMPLICIT NONE

   integer, intent(in) :: idate(3)
   integer, intent(in) :: tracer_idx, var_idx

   integer  :: ntime, iblkme, ib, jb, i, j, ilon, ilat
   real(r8) :: calday, cosz
   type(timestamp) :: tstamp

      tstamp = idate !tstamp_LB(tracer_idx,var_idx)
      ntime = 0
      DO WHILE (tstamp < tstamp_UB(tracer_idx,var_idx))
         ntime  = ntime + 1
         tstamp = tstamp + deltim_int
      ENDDO

      tstamp = idate !tstamp_LB(tracer_idx,var_idx)
      CALL flush_block_data (traceravgcos(tracer_idx,var_idx), 0._r8)

      DO WHILE (tstamp < tstamp_UB(tracer_idx,var_idx))

         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            DO j = 1, gforc_tracers(tracer_idx,var_idx)%ycnt(jb)
               DO i = 1, gforc_tracers(tracer_idx,var_idx)%xcnt(ib)

                  ilat = gforc_tracers(tracer_idx,var_idx)%ydsp(jb) + j
                  ilon = gforc_tracers(tracer_idx,var_idx)%xdsp(ib) + i
                  IF (ilon > gforc_tracers(tracer_idx,var_idx)%nlon) ilon = ilon - gforc_tracers(tracer_idx,var_idx)%nlon

                  calday = calendarday(tstamp)
                  cosz = orb_coszen(calday, gforc_tracers(tracer_idx,var_idx)%rlon(ilon), gforc_tracers(tracer_idx,var_idx)%rlat(ilat))
                  cosz = max(0.001, cosz)
                  traceravgcos(tracer_idx,var_idx)%blk(ib,jb)%val(i,j) = traceravgcos(tracer_idx,var_idx)%blk(ib,jb)%val(i,j) &
                     + cosz / real(ntime,r8) !  * deltim_real /real(tstamp_UB(tracer_idx,var_idx)-tstamp_LB(tracer_idx,var_idx))

               ENDDO
            ENDDO
         ENDDO

         tstamp = tstamp + deltim_int

      ENDDO

   END SUBROUTINE calavgcos

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Convert a string to lowercase
!
! !REVISIONS:
!  2024: Added for case-insensitive string comparisons
!
!-----------------------------------------------------------------------
   SUBROUTINE to_lower(str)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(INOUT) :: str
      INTEGER :: i, str_len
      
      str_len = LEN_TRIM(str)
      DO i = 1, str_len
         IF (str(i:i) >= 'A' .AND. str(i:i) <= 'Z') THEN
            str(i:i) = CHAR(ICHAR(str(i:i)) + 32)
         ENDIF
      ENDDO
      
      ! Clear any remaining characters to avoid buffer issues
      IF (str_len < LEN(str)) THEN
         str(str_len+1:) = ' '
      ENDIF
   END SUBROUTINE to_lower

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Get tracer forcing data for a specific tracer and variable
!
! !REVISIONS:
!  2024: Added for data access interface
!
!-----------------------------------------------------------------------
   FUNCTION get_tracer_forcing_data(tracer_idx, var_idx) RESULT(data_ptr)
   IMPLICIT NONE
   
   integer, intent(in) :: tracer_idx, var_idx
   real(r8), pointer :: data_ptr(:)
   
   nullify(data_ptr)
   
   IF (allocated(tracer_patch_data)) THEN
      IF (tracer_idx >= 1 .AND. tracer_idx <= SIZE(tracer_patch_data, 1) .AND. &
          var_idx >= 1 .AND. var_idx <= SIZE(tracer_patch_data, 2)) THEN
         data_ptr => tracer_patch_data(tracer_idx, var_idx, :)
      ELSE
         IF (p_is_master) THEN
            WRITE(*,*) "ERROR: get_tracer_forcing_data - Index out of bounds"
            WRITE(*,*) "  tracer_idx: ", tracer_idx, " max: ", SIZE(tracer_patch_data, 1)
            WRITE(*,*) "  var_idx: ", var_idx, " max: ", SIZE(tracer_patch_data, 2)
         ENDIF
      ENDIF
   ELSE
      IF (p_is_master) THEN
         WRITE(*,*) "ERROR: get_tracer_forcing_data - tracer_patch_data not allocated"
      ENDIF
   ENDIF
   
   END FUNCTION get_tracer_forcing_data


   FUNCTION tracerfilename(year, month, day, tracer_idx,var_i) RESULT(metfilename)

      USE MOD_Namelist
      USE MOD_SPMD_Task
      IMPLICIT NONE
   
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer, intent(in) :: tracer_idx
      integer, intent(in) :: var_i
      character(len=256)  :: metfilename
      character(len=256)  :: yearstr
      character(len=256)  :: monthstr
      character(len=256)  :: dataset_name_lower
   
         ! IF (p_is_master) THEN
         !    WRITE(*,*) "      === DEBUG: tracerfilename called ==="
         !    WRITE(*,*) "        Input: year=", year, " month=", month, " day=", day, " var_i=", var_i
         ! ENDIF
   
         write(yearstr, '(I4.4)') year
         write(monthstr, '(I2.2)') month
         
         dataset_name_lower = adjustl(trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name))
         CALL to_lower(dataset_name_lower)
   
         ! IF (p_is_master) THEN
         !    WRITE(*,*) "        Formatted strings: yearstr='", TRIM(yearstr), "' monthstr='", TRIM(monthstr), "'"
         !    WRITE(*,*) "        Dataset name (original): '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name), "'"
         !    WRITE(*,*) "        Dataset name (lower): '", dataset_name_lower, "'"
         ! ENDIF
   
         select CASE (dataset_name_lower)
         
         CASE ('isogsm')
            !DESCRIPTION
            !===========
               !--- Isotopes-incorporated Global Spectral Model (IsoGSM)
      
            !data source:
            !-------------------
               !---https://isotope.iis.u-tokyo.ac.jp/about-our-lab?lang=en
      
            !References:
            !-------------------
               !---Bong, H., Cauquoin, A., Okazaki, A., Chang, E.-C., Werner, M., Wei, Z., et al. (2024). 
               !   Process-based intercomparison of water isotope-enabled models and reanalysis nudging effects. 
               !   Journal of Geophysical Research: Atmospheres, 129, e2023JD038719. 
               !   https://doi.org/10.1029/2023JD038719
      
            !REVISION HISTORY
            !----------------
               !---2025.03.23   Zhongwang Wei @ SYSU: add the isotope forcing data
      
               ! Check if var_i is within bounds for fprefix array


               ! Construct filename based on the specific prefix for the variable, year, and .nc suffix
               metfilename = '/'//trim(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(var_i))//'_'//trim(yearstr)//'.nc'

   
         
         CASE ('POINT')
            metfilename = '/'//trim(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(1))
            ! IF (p_is_master) THEN
            !    WRITE(*,*) "        CASE POINT selected"
            !    WRITE(*,*) "        fprefix_tracer_forcing(1) = '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(1)), "'"
            !    WRITE(*,*) "        Generated filename: '", TRIM(metfilename), "'"
            ! ENDIF
            
         CASE DEFAULT
            ! IF (p_is_master) THEN
            !    WRITE(*,*) "        WARNING: Unknown dataset name '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name), "'"
            !    WRITE(*,*) "        Using default POINT format"
            ! ENDIF
            metfilename = '/'//trim(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(1))
         END select
         
         ! IF (p_is_master) THEN
         !    WRITE(*,*) "      === DEBUG: tracerfilename returning '", TRIM(metfilename), "' ==="
         ! ENDIF
         
         ! IF (DEF_USE_CBL_HEIGHT) THEN
         !    select CASE (var_i)
         !    CASE (9)
         !       metfilename = '/'//trim(fprefix_tracer_forcing(9))//'_'//trim(yearstr)//'_'//trim(monthstr)//&
         !          '_boundary_layer_height.nc4'
         !    END select
         ! ENDIF
      END FUNCTION tracerfilename




END MODULE MOD_Tracer_Forcing
