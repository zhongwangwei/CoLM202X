#include <define.h>

MODULE MOD_Tracer_Forcing

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  read in the atmospheric forcing using user defined interpolation method or
!  downscaling forcing
!

!-----------------------------------------------------------------------

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

   type (grid_type), PUBLIC :: gforc

   type (spatial_mapping_type) :: mg2p_forc   ! area weighted mapping from forcing to model unit

   type(grid_type), allocatable :: gforc_tracers(:,:) ! Grid for each tracer forcing variable [num_forced_tracers, max_nvar]
   type(spatial_mapping_type), allocatable :: mg2p_forc_tracers(:,:) ! Mapping for each [num_forced_tracers, max_nvar]

   logical, allocatable :: forcmask_pch (:)


   type(pointer_real8_1d), allocatable :: forc_prc_grid_O18   (:)
   type(pointer_real8_1d), allocatable :: forc_prl_grid_O18   (:)
   type(pointer_real8_1d), allocatable :: forc_prc_part_O18   (:)
   type(pointer_real8_1d), allocatable :: forc_prl_part_O18   (:)
   type(pointer_real8_1d), allocatable :: forc_q_part_O18     (:)
   type(pointer_real8_1d), allocatable :: forc_q_grid_O18     (:)

   type(pointer_real8_1d), allocatable :: forc_prc_grid_H2    (:)
   type(pointer_real8_1d), allocatable :: forc_prl_grid_H2    (:)
   type(pointer_real8_1d), allocatable :: forc_prc_part_H2    (:)
   type(pointer_real8_1d), allocatable :: forc_prl_part_H2    (:)
   type(pointer_real8_1d), allocatable :: forc_q_part_H2      (:)
   type(pointer_real8_1d), allocatable :: forc_q_grid_H2      (:)



   logical, allocatable :: glacierss (:)

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

   type(block_data_real8_2d) :: traceravgcos   ! time-average of cos(zenith)
   type(block_data_real8_2d) :: tracerdata  ! forcing data
#ifdef URBAN_MODEL
   type(block_data_real8_2d) :: rainf
   type(block_data_real8_2d) :: snowf
#endif

   type(block_data_real8_2d), allocatable :: forcn    (:)  ! forcing data
   type(block_data_real8_2d), allocatable :: forcn_LB (:)  ! forcing data at lower boundary
   type(block_data_real8_2d), allocatable :: forcn_UB (:)  ! forcing data at upper boundary
   
   ! Tracer patch data - stores tracer forcing data on patches
   ! Structure: [tracer_idx][var_idx] -> patch data array
   real(r8), allocatable, target :: tracer_patch_data(:,:,:)  ! [num_tracers, max_nvar, numpatch]
   
   PUBLIC :: tracer_forcing_init
   PUBLIC :: read_tracer_forcing
   PUBLIC :: tracer_forcing_final
   PUBLIC :: tracer_forcing_reset
   PUBLIC :: get_tracer_forcing_data

CONTAINS

!-----------------------------------------------------------------------
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
   num_tracers_with_forcing = SIZE(DEF_Tracer_Forcings_NL, 1)



   ! get value of deltim (moved outside loop)
   deltim_int  = int(deltatime)

   ! Calculate maximum number of variables across all tracers (moved outside loop)
   max_nvar = 0
   DO k = 1, num_tracers_with_forcing
      IF (DEF_Tracer_Forcings_NL(k)%NVAR > max_nvar) THEN
         max_nvar = DEF_Tracer_Forcings_NL(k)%NVAR
      ENDIF
   ENDDO

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
   IF (p_is_worker) THEN
      IF (numpatch > 0) THEN
         allocate (tracer_patch_data(num_tracers_with_forcing, max_nvar, numpatch))
         tracer_patch_data(:,:,:) = 0.0_r8
      ENDIF
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

         ! Allocate block data arrays once (moved outside loop)
      ! If no dedicated I/O processes exist, master process handles I/O
      IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN
         IF (allocated(forcn)) deallocate(forcn)
         IF (allocated(forcn_LB)) deallocate(forcn_LB)
         IF (allocated(forcn_UB)) deallocate(forcn_UB)
         allocate (forcn    (max_nvar))
         allocate (forcn_LB (max_nvar)) 
         allocate (forcn_UB (max_nvar))
      ENDIF

   DO i = 1, num_tracers_with_forcing
      ! Initialize user specified tracer forcing for this tracer (moved inside loop)
  !    CALL init_user_specified_tracer_forcing(i)
      
      ! Initialize grid for each tracer variable
      DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
         CALL tracer_read_latlon (DEF_Tracer_Forcings_NL(i)%tracer_dir, idate, i, j)
      ENDDO

               ! Allocate block data for this tracer
         ! If no dedicated I/O processes exist, master process handles I/O
         IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN
            DO ivar = 1, DEF_Tracer_Forcings_NL(i)%NVAR
               CALL allocate_block_data (gforc_tracers(i,ivar), forcn   (ivar))
               CALL allocate_block_data (gforc_tracers(i,ivar), forcn_LB(ivar))
               CALL allocate_block_data (gforc_tracers(i,ivar), forcn_UB(ivar))
            ENDDO

            ! Allocate memory for forcing data (using first grid for this tracer)
            CALL allocate_block_data (gforc_tracers(i,1), tracerdata)  ! forcing data
            CALL allocate_block_data (gforc_tracers(i,1), traceravgcos )  ! time-average of cos(zenith)
#if (defined URBAN_MODEL && defined SinglePoint)
            CALL allocate_block_data (gforc_tracers(i,1), rainf)
            CALL allocate_block_data (gforc_tracers(i,1), snowf)
#endif
         ENDIF

      ! Process missing value information if needed
      IF (DEF_Tracer_Forcings_NL(i)%has_missing_value) THEN
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
         CALL ncio_read_block_time (filename, DEF_Tracer_Forcings_NL(i)%vname(1), gforc_tracers(i,1), time_i, tracerdata)
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
      IF (DEF_Tracer_Forcings_NL(i)%has_missing_value) THEN
         CALL mg2p_forc_tracers(i,1)%set_missing_value (tracerdata, missing_value, forcmask_pch)
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
      IF (allocated(glacierss   )) deallocate(glacierss   )
      IF (allocated(forctime    )) deallocate(forctime    )
      IF (allocated(iforctime   )) deallocate(iforctime   )
      IF (allocated(forc_disk   )) deallocate(forc_disk   )
      IF (allocated(tstamp_LB   )) deallocate(tstamp_LB   )
      IF (allocated(tstamp_UB   )) deallocate(tstamp_UB   )
      IF (allocated(tracer_patch_data)) deallocate(tracer_patch_data)
      IF (allocated(gforc_tracers)) deallocate(gforc_tracers)
      IF (allocated(mg2p_forc_tracers)) deallocate(mg2p_forc_tracers)
      IF (allocated(forcn)) deallocate(forcn)
      IF (allocated(forcn_LB)) deallocate(forcn_LB)
      IF (allocated(forcn_UB)) deallocate(forcn_UB)

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

   ! Local variables for range checking
   real(r8), allocatable :: temp_tracer_data(:)
   character(len=64) :: var_name
   character(len=8) :: tracer_str, var_str
   integer :: num_tracers_with_forcing
   
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
   
   IF (p_is_master) THEN
      WRITE(*,*) "=== READ_TRACER_FORCING START ==="
      WRITE(*,*) "  Processing date: ", idate(1), idate(2), idate(3)
      WRITE(*,*) "  Number of tracers with forcing: ", num_tracers_with_forcing
   ENDIF
   
   ! Debug: Show process information
   IF (p_is_worker .AND. p_iam_glb <= 3) THEN  ! Only first few workers
      WRITE(*,*) "  WORKER ", p_iam_glb, ": Starting tracer forcing processing with ", numpatch, " patches"
   ENDIF

   ! Check if tracer_patch_data is properly allocated
   IF (p_is_worker .AND. numpatch > 0) THEN
      IF (.NOT. allocated(tracer_patch_data)) THEN
         IF (p_is_master) WRITE(*,*) "ERROR: tracer_patch_data not allocated!"
         CALL CoLM_stop()
      ENDIF
      
      ! Verify the dimensions
      IF (SIZE(tracer_patch_data, 1) < num_tracers_with_forcing .OR. &
          SIZE(tracer_patch_data, 3) /= numpatch) THEN
         IF (p_is_master) THEN
            WRITE(*,*) "ERROR: tracer_patch_data dimension mismatch!"
            WRITE(*,*) "  Expected: ", num_tracers_with_forcing, " x ", max_nvar, " x ", numpatch
            WRITE(*,*) "  Actual: ", SIZE(tracer_patch_data, 1), " x ", SIZE(tracer_patch_data, 2), " x ", SIZE(tracer_patch_data, 3)
         ENDIF
         CALL CoLM_stop()
      ENDIF
   ENDIF
   ! Loop over each tracer
   DO tracer_idx = 1, num_tracers_with_forcing
      
      ! Debug: Show process type and tracer info (visible from all processes)
      IF (p_is_master) THEN
         WRITE(*,*) "=== PROCESSING TRACER ", tracer_idx, " ==="
         WRITE(*,*) "  Tracer name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name)
         WRITE(*,*) "  Dataset: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
         WRITE(*,*) "  Directory: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir)
         WRITE(*,*) "  Number of variables: ", DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
         WRITE(*,*) "  p_is_io: ", p_is_io
         WRITE(*,*) "  p_is_worker: ", p_is_worker
         WRITE(*,*) "  p_is_master: ", p_is_master
      ENDIF
      
      !------------------------------------------------------------
      ! READ in THE TRACER FORCING
      ! read lower and upper boundary forcing data for this tracer
      ! If no dedicated I/O processes exist, master process handles I/O
      IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN
         IF (p_is_master) WRITE(*,*) "  >>> I/O process reading tracer forcing data"
         CALL metreadLBUB_tracer(idate, TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir), tracer_idx)
         
         ! Debug: Show timestamps after reading
         IF (p_is_master) THEN
            WRITE(*,*) "  After metreadLBUB_tracer for tracer ", tracer_idx, ":"
            DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
               IF (tracer_idx <= SIZE(tstamp_LB, 1) .AND. var_idx <= SIZE(tstamp_LB, 2) .AND. &
                   tracer_idx <= SIZE(tstamp_UB, 1) .AND. var_idx <= SIZE(tstamp_UB, 2)) THEN
                  WRITE(*,*) "    Variable ", var_idx, " LB: ", tstamp_LB(tracer_idx,var_idx)%year, tstamp_LB(tracer_idx,var_idx)%day, tstamp_LB(tracer_idx,var_idx)%sec
                  WRITE(*,*) "    Variable ", var_idx, " UB: ", tstamp_UB(tracer_idx,var_idx)%year, tstamp_UB(tracer_idx,var_idx)%day, tstamp_UB(tracer_idx,var_idx)%sec
               ENDIF
            ENDDO
         ENDIF
   
         ! set model time stamp
         id(:) = idate(:)
         mtstamp = id
         
         ! Debug: Show model timestamp
         IF (p_is_master) THEN
            WRITE(*,*) "  Model timestamp: ", mtstamp%year, mtstamp%day, mtstamp%sec
         ENDIF
         
         ! loop for variables of this tracer
         DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
            IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) == 'NULL') CYCLE     ! no data, CYCLE
            IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx)) == 'NULL') CYCLE

            ! Debug information for each variable
            IF (p_is_master) THEN
               WRITE(*,*) "  Processing tracer ", tracer_idx, " variable ", var_idx, ":"
               WRITE(*,*) "    Variable name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx))
               WRITE(*,*) "    Interpolation: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx))
               WRITE(*,*) "    Time log: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%timelog(var_idx))
            ENDIF

            ! Check array bounds before accessing
            IF (tracer_idx > SIZE(tstamp_LB, 1) .OR. var_idx > SIZE(tstamp_LB, 2) .OR. &
                tracer_idx > SIZE(tstamp_UB, 1) .OR. var_idx > SIZE(tstamp_UB, 2) .OR. &
                var_idx > SIZE(forcn_LB) .OR. var_idx > SIZE(forcn_UB) .OR. &
                var_idx > SIZE(forcn)) THEN
               IF (p_is_master) THEN
                  WRITE(*,*) "ERROR: Array index out of bounds in read_tracer_forcing"
                  WRITE(*,*) "  tracer_idx = ", tracer_idx, " var_idx = ", var_idx, " but array sizes are:"
                  WRITE(*,*) "  tstamp_LB size = ", SIZE(tstamp_LB, 1), "x", SIZE(tstamp_LB, 2)
                  WRITE(*,*) "  tstamp_UB size = ", SIZE(tstamp_UB, 1), "x", SIZE(tstamp_UB, 2)
                  WRITE(*,*) "  forcn_LB size = ", SIZE(forcn_LB)
                  WRITE(*,*) "  forcn_UB size = ", SIZE(forcn_UB)
                  WRITE(*,*) "  forcn size = ", SIZE(forcn)
               ENDIF
               CYCLE
            ENDIF

            ! to make sure the forcing data calculated is in the range of time
            ! interval [LB, UB]
            IF ( (mtstamp < tstamp_LB(tracer_idx,var_idx)) .or. (tstamp_UB(tracer_idx,var_idx) < mtstamp) ) THEN
               IF (p_is_master) THEN
                  WRITE(*,*) "=== TRACER FORCING DEBUG: OUT OF RANGE ERROR ==="
                  WRITE(*,*) "  Tracer index: ", tracer_idx
                  WRITE(*,*) "  Variable index: ", var_idx
                  WRITE(*,*) "  Tracer name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name)
                  WRITE(*,*) "  Variable name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx))
                  WRITE(*,*) "  Model timestamp: ", mtstamp%year, mtstamp%day, mtstamp%sec
                  WRITE(*,*) "  Lower boundary timestamp: ", tstamp_LB(tracer_idx,var_idx)%year, tstamp_LB(tracer_idx,var_idx)%day, tstamp_LB(tracer_idx,var_idx)%sec
                  WRITE(*,*) "  Upper boundary timestamp: ", tstamp_UB(tracer_idx,var_idx)%year, tstamp_UB(tracer_idx,var_idx)%day, tstamp_UB(tracer_idx,var_idx)%sec
                  WRITE(*,*) "  Condition check:"
                  WRITE(*,*) "    mtstamp < tstamp_LB: ", (mtstamp < tstamp_LB(tracer_idx,var_idx))
                  WRITE(*,*) "    tstamp_UB < mtstamp: ", (tstamp_UB(tracer_idx,var_idx) < mtstamp)
                  WRITE(*,*) "  Tracer forcing configuration:"
                  WRITE(*,*) "    Dataset: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
                  WRITE(*,*) "    Directory: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir)
                  WRITE(*,*) "    Start year/month: ", DEF_Tracer_Forcings_NL(tracer_idx)%startyr, DEF_Tracer_Forcings_NL(tracer_idx)%startmo
                  WRITE(*,*) "    End year/month: ", DEF_Tracer_Forcings_NL(tracer_idx)%endyr, DEF_Tracer_Forcings_NL(tracer_idx)%endmo
                  WRITE(*,*) "    Group by: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby)
                  WRITE(*,*) "================================================"
               ENDIF
               write(6, *) "the data required is out of range! STOP!"; CALL CoLM_stop()
            ENDIF

            ! calculate distance to lower/upper boundary
            dtLB = mtstamp - tstamp_LB(tracer_idx,var_idx)
            dtUB = tstamp_UB(tracer_idx,var_idx) - mtstamp

            ! linear method, for most variables
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'linear') THEN
               IF ( (dtLB+dtUB) > 0 ) THEN
                  CALL block_data_linear_interp ( &
                     forcn_LB(var_idx), real(dtUB,r8)/real(dtLB+dtUB,r8), &
                     forcn_UB(var_idx), real(dtLB,r8)/real(dtLB+dtUB,r8), &
                     forcn(var_idx))
               ELSE
                  CALL block_data_copy (forcn_LB(var_idx), forcn(var_idx))
               ENDIF
            ENDIF

            ! nearest method, for precipitation
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'nearest') THEN
               IF (dtLB <= dtUB) THEN
                  CALL block_data_copy (forcn_LB(var_idx), forcn(var_idx))
               ELSE
                  CALL block_data_copy (forcn_UB(var_idx), forcn(var_idx))
               ENDIF
            ENDIF

            ! uniform method
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'uniform') THEN
               IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%timelog(var_idx)) == 'forward') THEN
                  CALL block_data_copy (forcn_LB(var_idx), forcn(var_idx))
               ELSE
                  CALL block_data_copy (forcn_UB(var_idx), forcn(var_idx))
               ENDIF
            ENDIF

            ! coszen method, for SW radiation
            IF (DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(var_idx) == 'coszen') THEN
               ! Check if gforc_tracers is properly allocated and accessible
               IF (tracer_idx <= SIZE(gforc_tracers, 1) .AND. var_idx <= SIZE(gforc_tracers, 2)) THEN
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
                              forcn(var_idx)%blk(ib,jb)%val(i,j) = &
                                 cosz / traceravgcos%blk(ib,jb)%val(i,j) * forcn_LB(var_idx)%blk(ib,jb)%val(i,j)
                           ELSE
                              forcn(var_idx)%blk(ib,jb)%val(i,j) = &
                                 cosz / traceravgcos%blk(ib,jb)%val(i,j) * forcn_UB(var_idx)%blk(ib,jb)%val(i,j)
                           ENDIF

                        ENDDO
                     ENDDO
                  ENDDO
               ELSE
                  IF (p_is_master) THEN
                     WRITE(*,*) "ERROR: gforc_tracers array bounds exceeded in coszen interpolation"
                     WRITE(*,*) "  tracer_idx = ", tracer_idx, " var_idx = ", var_idx
                     WRITE(*,*) "  gforc_tracers dimensions = ", SIZE(gforc_tracers, 1), "x", SIZE(gforc_tracers, 2)
                  ENDIF
               ENDIF
            ENDIF

         ENDDO

         ! Debug: Check interpolated forcing data values
         IF (p_is_master) THEN
            WRITE(*,*) "  After interpolation for tracer ", tracer_idx, ":"
            DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
               IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) /= 'NULL' .AND. var_idx <= SIZE(forcn)) THEN
                  IF (allocated(forcn(var_idx)%blk)) THEN
                     ! Check first block for sample values
                     IF (SIZE(forcn(var_idx)%blk) > 0) THEN
                        DO ib = 1, MIN(2, SIZE(forcn(var_idx)%blk, 1))
                           DO jb = 1, MIN(2, SIZE(forcn(var_idx)%blk, 2))
                              IF (allocated(forcn(var_idx)%blk(ib,jb)%val)) THEN
                                 WRITE(*,*) "    Variable ", var_idx, " (", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)), ") block(", ib, ",", jb, ") sample values:", &
                                           forcn(var_idx)%blk(ib,jb)%val(1:MIN(3,SIZE(forcn(var_idx)%blk(ib,jb)%val,1)), &
                                                                        1:MIN(3,SIZE(forcn(var_idx)%blk(ib,jb)%val,2)))
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSE
                     WRITE(*,*) "    Variable ", var_idx, " (", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)), ") forcn%blk not allocated"
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

                  ! preprocess for forcing data
         CALL tracerpreprocess (gforc_tracers(tracer_idx,1), forcn)

      ENDIF

      ! Map tracer forcing data to patches using the appropriate spatial mapping
      ! (This must be done on I/O processes since they have the forcn arrays)
      ! If no dedicated I/O processes exist, master process handles I/O
      IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN
         IF (p_is_master) WRITE(*,*) "  >>> I/O process mapping tracer forcing data to patches"
         DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
            ! Check array bounds before accessing mg2p_forc_tracers
            IF (tracer_idx <= SIZE(mg2p_forc_tracers, 1) .AND. var_idx <= SIZE(mg2p_forc_tracers, 2)) THEN
               ! Skip NULL variables
               IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) == 'NULL') CYCLE
               
               ! Debug: Check conditions for mapping
               IF (p_is_master .AND. var_idx == 1) THEN
                  WRITE(*,*) "    Debug mapping conditions for tracer ", tracer_idx, ":"
                  WRITE(*,*) "      forcn allocated: ", allocated(forcn)
                  IF (allocated(forcn)) WRITE(*,*) "      forcn size: ", SIZE(forcn)
                  WRITE(*,*) "      var_idx: ", var_idx
                  WRITE(*,*) "      p_is_worker: ", p_is_worker
                  WRITE(*,*) "      numpatch: ", numpatch
                  IF (allocated(forcn) .AND. var_idx <= SIZE(forcn)) THEN
                     WRITE(*,*) "      forcn(", var_idx, ")%blk allocated: ", allocated(forcn(var_idx)%blk)
                  ENDIF
               ENDIF
               
               ! Only perform mapping if we have valid data and proper array dimensions
               IF (allocated(forcn) .AND. var_idx <= SIZE(forcn) .AND. &
                   allocated(forcn(var_idx)%blk) .AND. &
                   tracer_idx <= SIZE(tracer_patch_data, 1) .AND. &
                   var_idx <= SIZE(tracer_patch_data, 2)) THEN
                  
                  ! Debug: Report mapping attempt
                  IF (p_is_master) THEN
                     WRITE(*,*) "    Calling grid2pset for tracer ", tracer_idx, " variable ", var_idx, " (", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)), ")"
                  ENDIF
                  
                  ! Call spatial mapping function - this handles MPI communication automatically
                  CALL mg2p_forc_tracers(tracer_idx,var_idx)%grid2pset (forcn(var_idx), tracer_patch_data(tracer_idx,var_idx,:))
                  
                  ! Debug: Report successful mapping on worker processes
                  IF (p_is_worker .AND. p_iam_glb <= 3 .AND. numpatch > 0) THEN  ! Only first few workers to avoid spam
                     WRITE(*,*) "    WORKER ", p_iam_glb, ": Successfully mapped tracer ", tracer_idx, " variable ", var_idx, " - first 3 values: ", &
                                tracer_patch_data(tracer_idx,var_idx,1:min(3,numpatch))
                     WRITE(*,*) "    WORKER ", p_iam_glb, ": Min/Max values: ", &
                                MINVAL(tracer_patch_data(tracer_idx,var_idx,:)), "/", MAXVAL(tracer_patch_data(tracer_idx,var_idx,:))
                  ENDIF
               ELSE
                  IF (p_is_master) THEN
                     WRITE(*,*) "    Skipping mapping for tracer ", tracer_idx, " variable ", var_idx, " - invalid conditions"
                     WRITE(*,*) "      allocated(forcn): ", allocated(forcn)
                     IF (allocated(forcn)) WRITE(*,*) "      SIZE(forcn): ", SIZE(forcn)
                     WRITE(*,*) "      var_idx: ", var_idx
                     IF (allocated(forcn) .AND. var_idx <= SIZE(forcn)) THEN
                        WRITE(*,*) "      forcn(", var_idx, ")%blk allocated: ", allocated(forcn(var_idx)%blk)
                     ENDIF
                     WRITE(*,*) "      tracer_patch_data dimensions: ", SIZE(tracer_patch_data, 1), "x", SIZE(tracer_patch_data, 2), "x", SIZE(tracer_patch_data, 3)
                  ENDIF
               ENDIF
            ELSE
               IF (p_is_master) THEN
                  WRITE(*,*) "ERROR: mg2p_forc_tracers array bounds exceeded"
                  WRITE(*,*) "  tracer_idx = ", tracer_idx, " var_idx = ", var_idx
                  WRITE(*,*) "  mg2p_forc_tracers dimensions = ", SIZE(mg2p_forc_tracers, 1), "x", SIZE(mg2p_forc_tracers, 2)
               ENDIF
            ENDIF
         ENDDO
      ELSE
         IF (p_is_master) WRITE(*,*) "  >>> NOT an I/O process, skipping file operations for tracer ", tracer_idx
      ENDIF

   ENDDO

   IF (p_is_master) THEN
      WRITE(*,*) "=== READ_TRACER_FORCING COMPLETED SUCCESSFULLY ==="
   ENDIF
   
   ! Debug: Show worker completion status
   IF (p_is_worker .AND. p_iam_glb <= 3) THEN  ! Only first few workers
      WRITE(*,*) "  WORKER ", p_iam_glb, ": Completed tracer forcing processing for ", numpatch, " patches"
      IF (numpatch > 0 .AND. allocated(tracer_patch_data)) THEN
         WRITE(*,*) "    Sample tracer data for patch 1: tracer 1 var 1 = ", tracer_patch_data(1,1,1)
      ENDIF
   ENDIF

#ifdef RangeCheck
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) write(*,'(/, A25)') 'Checking tracer forcing ...'

      ! Check tracer forcing data for each tracer and variable (only on worker processes with patches)
      IF (p_is_worker .AND. numpatch > 0) THEN
         DO tracer_idx = 1, num_tracers_with_forcing
            IF (allocated(DEF_Tracer_Forcings_NL) .AND. tracer_idx <= SIZE(DEF_Tracer_Forcings_NL, 1)) THEN
               DO var_idx = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR
                  IF (var_idx <= max_nvar .AND. allocated(tracer_patch_data) .AND. &
                      var_idx <= SIZE(DEF_Tracer_Forcings_NL(tracer_idx)%vname) .AND. &
                      trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) /= 'NULL') THEN
                  ! Create a temporary array for the current tracer variable
                  IF (allocated(temp_tracer_data)) deallocate(temp_tracer_data)
                  allocate(temp_tracer_data(numpatch))
                  temp_tracer_data(:) = tracer_patch_data(tracer_idx, var_idx, :)
                  
                  ! Generate variable name for checking - use meaningful names when available
                  IF (var_idx <= SIZE(DEF_Tracer_Forcings_NL(tracer_idx)%vname) .AND. &
                      trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx)) /= 'NULL') THEN
                     var_name = 'Tracer ' // trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name) // &
                                ' ' // trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx))
                  ELSE
                     write(tracer_str, '(I0)') tracer_idx
                     write(var_str, '(I0)') var_idx
                     var_name = 'Tracer_' // trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name) // '_var_' // trim(var_str)
                  ENDIF
                  
                  ! Check the tracer forcing data
                  CALL check_vector_data (trim(var_name), temp_tracer_data)
                  
                  deallocate(temp_tracer_data)
                                 ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif

   END SUBROUTINE read_tracer_forcing

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  read lower and upper boundary forcing data for a specific tracer, a major interface of this
!  MODULE
!
! !REVISIONS:
!  04/2014, Hua Yuan: initial code
!  2024: Modified for multiple tracers
!
!-----------------------------------------------------------------------
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


      mtstamp = idate

      DO ivar = 1, DEF_Tracer_Forcings_NL(tracer_idx)%NVAR

         IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)) == 'NULL') THEN
            CYCLE     ! no data, CYCLE
         ENDIF



         ! Check array bounds before accessing
         IF (tracer_idx > SIZE(tstamp_LB, 1) .OR. ivar > SIZE(tstamp_LB, 2) .OR. &
             tracer_idx > SIZE(tstamp_UB, 1) .OR. ivar > SIZE(tstamp_UB, 2) .OR. &
             ivar > SIZE(forcn_LB) .OR. ivar > SIZE(forcn_UB)) THEN
            CYCLE
         ENDIF


         ! lower and upper boundary data already exist, CYCLE
         IF ( .not.(tstamp_LB(tracer_idx,ivar)%year == -1) .and. .not.(tstamp_UB(tracer_idx,ivar)%year == -1) .and. &
            tstamp_LB(tracer_idx,ivar)%year<=mtstamp%year .and. mtstamp%year<tstamp_UB(tracer_idx,ivar)%year ) THEN
            CYCLE
         ENDIF

         ! set lower boundary time stamp and get data
         IF (tstamp_LB(tracer_idx,ivar)%year == -1) THEN
            CALL setstampLB_tracer(mtstamp, ivar, tracer_idx, year, month, day, time_i)

            filename = trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir)//trim(tracerfilename(year, month, day, tracer_idx, ivar))
            
            IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN ! Debug on I/O ranks
                WRITE(*,*) "    DEBUG IO: Attempting to read LB for tracer_idx=", tracer_idx, " var_idx=", ivar
                WRITE(*,*) "      Tracer Name: ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name)
                WRITE(*,*) "      Directory:   ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir)
                WRITE(*,*) "      File Func Input (y,m,d,ti): ", year, month, day, time_i
                WRITE(*,*) "      Generated Filename: ", TRIM(filename)
                WRITE(*,*) "      Namelist Var Name: ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar))
                WRITE(*,*) "      Time Index (from setstamp): ", time_i
                WRITE(*,*) "      Dataset type: ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
            ENDIF
            
            ! Debug: Show file being read
            IF (p_is_master) THEN
               WRITE(*,*) "    Reading LB data from file: ", trim(filename)
               WRITE(*,*) "    Variable: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)), " time_i: ", time_i
            ENDIF
  
            
            ! Convert dataset name to lowercase for case-insensitive comparison
            dataset_name_lower = trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
            CALL to_lower(dataset_name_lower)
     
            
            IF (trim(dataset_name_lower) == 'point') THEN
               IF (forcing_read_ahead) THEN
                  tracerdata%blk(gblock%xblkme(1),gblock%yblkme(1))%val = forc_disk(time_i,ivar)
               ELSE
#ifndef URBAN_MODEL
                  CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata)
#else
                  IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)) == 'Rainf') THEN
                     CALL ncio_read_site_time (filename, 'Rainf', time_i, rainf)
                     CALL ncio_read_site_time (filename, 'Snowf', time_i, snowf)

                     DO iblkme = 1, gblock%nblkme
                        ib = gblock%xblkme(iblkme)
                        jb = gblock%yblkme(iblkme)

                        tracerdata%blk(ib,jb)%val(1,1) = rainf%blk(ib,jb)%val(1,1) + snowf%blk(ib,jb)%val(1,1)
                     ENDDO
                  ELSE
                     CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata)
                  ENDIF
#endif
               ENDIF
            ELSE
               CALL ncio_read_block_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), gforc_tracers(tracer_idx,ivar), time_i, tracerdata)
            ENDIF

            CALL block_data_copy (tracerdata, forcn_LB(ivar))
            
            ! Debug: Check values after reading LB data
            IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN ! Only on I/O ranks
               IF (allocated(tracerdata%blk)) THEN
                  IF (SIZE(tracerdata%blk) > 0) THEN
                     DO ib = 1, MIN(1, SIZE(tracerdata%blk, 1))
                        DO jb = 1, MIN(1, SIZE(tracerdata%blk, 2))
                           IF (allocated(tracerdata%blk(ib,jb)%val)) THEN
                              WRITE(*,*) "    DEBUG IO: After ncio_read (LB) for tracer ", tracer_idx, " var ", ivar, " (", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)), ")"
                              WRITE(*,*) "      tracerdata sample: ", tracerdata%blk(ib,jb)%val(1:MIN(3,SIZE(tracerdata%blk(ib,jb)%val,1)),1:MIN(3,SIZE(tracerdata%blk(ib,jb)%val,2)))
                              WRITE(*,*) "      tracerdata min/max: ", MINVAL(tracerdata%blk(ib,jb)%val), MAXVAL(tracerdata%blk(ib,jb)%val)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE
                  WRITE(*,*) "    DEBUG IO: tracerdata%blk not allocated after ncio_read (LB) for tracer ", tracer_idx, " var ", ivar
               ENDIF
            ENDIF

         ENDIF

         ! set upper boundary time stamp and get data
         IF (tstamp_UB(tracer_idx,ivar)%year == -1) THEN
            CALL setstampUB_tracer(ivar, tracer_idx, year, month, day, time_i)

            IF (year <= DEF_Tracer_Forcings_NL(tracer_idx)%endyr) THEN
               ! read forcing data
               filename = trim(dir_forcing)//trim(tracerfilename(year, month, day,tracer_idx, ivar))

               IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN ! Debug on I/O ranks
                  WRITE(*,*) "    DEBUG IO: Attempting to read UB for tracer_idx=", tracer_idx, " var_idx=", ivar
                  WRITE(*,*) "      Tracer Name: ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name)
                  WRITE(*,*) "      Directory:   ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_dir)
                  WRITE(*,*) "      File Func Input (y,m,d,ti): ", year, month, day, time_i
                  WRITE(*,*) "      Generated Filename: ", TRIM(filename)
                  WRITE(*,*) "      Namelist Var Name: ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar))
                  WRITE(*,*) "      Time Index (from setstamp): ", time_i
                  WRITE(*,*) "      Dataset type: ", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
               ENDIF
               
               ! Debug: Show UB file being read
               IF (p_is_master) THEN
                  WRITE(*,*) "    Reading UB data from file: ", trim(filename)
               ENDIF
               
   
               ! Convert dataset name to lowercase for case-insensitive comparison
               dataset_name_lower = trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
               CALL to_lower(dataset_name_lower)
               
               IF (trim(dataset_name_lower) == 'point') THEN

                  IF (forcing_read_ahead) THEN
                     tracerdata%blk(gblock%xblkme(1),gblock%yblkme(1))%val = forc_disk(time_i,ivar)
                  ELSE
#ifndef URBAN_MODEL
                     IF (p_is_master) WRITE(*,*) "    Calling ncio_read_site_time for UB variable: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar))
                     CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata)
#else
                     IF (trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)) == 'Rainf') THEN
                        CALL ncio_read_site_time (filename, 'Rainf', time_i, rainf)
                        CALL ncio_read_site_time (filename, 'Snowf', time_i, snowf)

                        DO iblkme = 1, gblock%nblkme
                           ib = gblock%xblkme(iblkme)
                           jb = gblock%yblkme(iblkme)

                           tracerdata%blk(ib,jb)%val(1,1) = rainf%blk(ib,jb)%val(1,1) + snowf%blk(ib,jb)%val(1,1)
                        ENDDO
                     ELSE
                        CALL ncio_read_site_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), time_i, tracerdata)
                     ENDIF
#endif
                  ENDIF
               ELSE
                  CALL ncio_read_block_time (filename, DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar), gforc_tracers(tracer_idx,ivar), time_i, tracerdata)
               ENDIF

               CALL block_data_copy (tracerdata, forcn_UB(ivar))
               
               ! Debug: Check values after reading UB data
               IF (p_is_io .OR. (p_is_master .AND. p_np_io == 0)) THEN ! Only on I/O ranks
                  IF (allocated(tracerdata%blk)) THEN
                     IF (SIZE(tracerdata%blk) > 0) THEN
                        DO ib = 1, MIN(1, SIZE(tracerdata%blk, 1))
                           DO jb = 1, MIN(1, SIZE(tracerdata%blk, 2))
                              IF (allocated(tracerdata%blk(ib,jb)%val)) THEN
                                 WRITE(*,*) "    DEBUG IO: After ncio_read (UB) for tracer ", tracer_idx, " var ", ivar, " (", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%vname(ivar)), ")"
                                 WRITE(*,*) "      tracerdata sample: ", tracerdata%blk(ib,jb)%val(1:MIN(3,SIZE(tracerdata%blk(ib,jb)%val,1)),1:MIN(3,SIZE(tracerdata%blk(ib,jb)%val,2)))
                                 WRITE(*,*) "      tracerdata min/max: ", MINVAL(tracerdata%blk(ib,jb)%val), MAXVAL(tracerdata%blk(ib,jb)%val)
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSE
                     WRITE(*,*) "    DEBUG IO: tracerdata%blk not allocated after ncio_read (UB) for tracer ", tracer_idx, " var ", ivar
                  ENDIF
               ENDIF
            ELSE
               write(*,*) year, DEF_Tracer_Forcings_NL(tracer_idx)%endyr
               print *, 'NOTE: reaching the END of forcing data, always reuse the last time step data!'
            ENDIF

            
            IF (TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%tintalgo(ivar)) == 'coszen') THEN
               IF (p_is_master) WRITE(*,*) "    Calculating avgcos for tracer ", tracer_idx, " var ", ivar
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

         mtstamp = idate

         CALL setstampLB_tracer(mtstamp, 1, ivar, year, month, day, time_i)
         filename = trim(DEF_Tracer_Forcings_NL(ivar)%tracer_dir)//trim(tracerfilename(year, month, day, ivar, 1))
         tstamp_LB(ivar,1) = timestamp(-1, -1, -1)  ! Fixed: use 2D array indexing

         IF (DEF_Tracer_Forcings_NL(ivar)%dim2d) THEN
            CALL ncio_read_bcast_serial (filename, DEF_Tracer_Forcings_NL(ivar)%latname, latxy)
            CALL ncio_read_bcast_serial (filename, DEF_Tracer_Forcings_NL(ivar)%lonname, lonxy)

            allocate (lat_in (size(latxy,2)))
            allocate (lon_in (size(lonxy,1)))
            lat_in = latxy(1,:)
            lon_in = lonxy(:,1)

            deallocate (latxy)
            deallocate (lonxy)
         ELSE
            CALL ncio_read_bcast_serial (filename, DEF_Tracer_Forcings_NL(ivar)%latname, lat_in)
            CALL ncio_read_bcast_serial (filename, DEF_Tracer_Forcings_NL(ivar)%lonname, lon_in)
         ENDIF

         IF (.not. DEF_Tracer_Forcings_NL(ivar)%regional) THEN
            CALL gforc_tracers(ivar,jvar)%define_by_center (lat_in, lon_in)
         ELSE
            CALL gforc_tracers(ivar,jvar)%define_by_center (lat_in, lon_in, &
               south = DEF_Tracer_Forcings_NL(ivar)%regbnd(1), north = DEF_Tracer_Forcings_NL(ivar)%regbnd(2), &
               west  = DEF_Tracer_Forcings_NL(ivar)%regbnd(3), east  = DEF_Tracer_Forcings_NL(ivar)%regbnd(4))
         ENDIF

         deallocate (lat_in)
         deallocate (lon_in)
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

   ! Debug information at entry
   IF (p_is_master) THEN
      WRITE(*,*) "    === SETSTAMPLB_TRACER DEBUG: ENTRY ==="
      WRITE(*,*) "      Input mtstamp: ", mtstamp%year, mtstamp%day, mtstamp%sec
      WRITE(*,*) "      Variable index: ", var_i
      WRITE(*,*) "      Tracer index: ", tracer_idx
      WRITE(*,*) "      Variable name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_i))
      WRITE(*,*) "      Dataset name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
      WRITE(*,*) "      Group by: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby)
      WRITE(*,*) "      Dtime: ", DEF_Tracer_Forcings_NL(tracer_idx)%dtime
      WRITE(*,*) "      ==================================="
   ENDIF

      year = mtstamp%year
      day  = mtstamp%day
      sec  = mtstamp%sec

      ! Convert dataset name to lowercase for case-insensitive comparison
      dataset_name_lower = trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
      CALL to_lower(dataset_name_lower)

      IF (trim(dataset_name_lower) == 'point') THEN
         IF (p_is_master) WRITE(*,*) "      Processing POINT dataset"

         ! For POINT data, we need to handle this differently
         ! Since forctime is not initialized, we'll use a simplified approach
         tstamp_LB(tracer_idx,var_i)%year = year
         tstamp_LB(tracer_idx,var_i)%day  = day
         tstamp_LB(tracer_idx,var_i)%sec  = sec
         time_i = 1

         IF (p_is_master) THEN
            WRITE(*,*) "      POINT dataset - simplified approach:"
            WRITE(*,*) "        tstamp_LB: ", tstamp_LB(tracer_idx,var_i)%year, tstamp_LB(tracer_idx,var_i)%day, tstamp_LB(tracer_idx,var_i)%sec
            WRITE(*,*) "        time_i: ", time_i
         ENDIF

         RETURN
      ENDIF

      tstamp_LB(tracer_idx,var_i)%year = year
      tstamp_LB(tracer_idx,var_i)%day  = day

      IF (p_is_master) WRITE(*,*) "      Initial tstamp_LB set to: ", tstamp_LB(tracer_idx,var_i)%year, tstamp_LB(tracer_idx,var_i)%day

      ! in the case of one year one file
      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'year' ) THEN
         IF (p_is_master) WRITE(*,*) "      Group by YEAR - using year-based file grouping"

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
      ENDIF

      ! in the case of one month one file
      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'month' ) THEN
         IF (p_is_master) WRITE(*,*) "      Group by MONTH - using month-based file grouping"

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
         IF (p_is_master) WRITE(*,*) "      Group by DAY - using day-based file grouping"

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

      IF (p_is_master) THEN
         WRITE(*,*) "      Final tstamp_LB: ", tstamp_LB(tracer_idx,var_i)%year, tstamp_LB(tracer_idx,var_i)%day, tstamp_LB(tracer_idx,var_i)%sec
         WRITE(*,*) "      Final time_i: ", time_i
         WRITE(*,*) "      Output values: year=", year, " month=", month, " mday=", mday
         WRITE(*,*) "    === SETSTAMPLB_TRACER DEBUG: EXIT ==="
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

   ! Debug information at entry
   IF (p_is_master) THEN
      WRITE(*,*) "    === SETSTAMPUB_TRACER DEBUG: ENTRY ==="
      WRITE(*,*) "      Variable index: ", var_i
      WRITE(*,*) "      Tracer index: ", tracer_idx
      WRITE(*,*) "      Variable name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_i))
      WRITE(*,*) "      Dataset name: ", trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
      WRITE(*,*) "      Current tstamp_LB: ", tstamp_LB(tracer_idx,var_i)%year, tstamp_LB(tracer_idx,var_i)%day, tstamp_LB(tracer_idx,var_i)%sec
      WRITE(*,*) "      Current tstamp_UB: ", tstamp_UB(tracer_idx,var_i)%year, tstamp_UB(tracer_idx,var_i)%day, tstamp_UB(tracer_idx,var_i)%sec
      WRITE(*,*) "      Dtime: ", DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
      WRITE(*,*) "      ==================================="
   ENDIF

      ! Convert dataset name to lowercase for case-insensitive comparison
      dataset_name_lower = trim(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
      CALL to_lower(dataset_name_lower)

      IF (trim(dataset_name_lower) == 'point') THEN
         IF (p_is_master) WRITE(*,*) "      Processing POINT dataset for UB"
         
         IF ( tstamp_UB(tracer_idx,var_i)%year == -1 ) THEN
            ! For POINT data, we need to handle this differently
            ! Since forctime is not initialized, we'll use a simplified approach
            tstamp_UB(tracer_idx,var_i) = tstamp_LB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
            IF (p_is_master) WRITE(*,*) "      Setting initial tstamp_UB for POINT data"
         ELSE
            ! For POINT data, we need to handle this differently
            ! For now, just increment the time
            tstamp_UB(tracer_idx,var_i) = tstamp_UB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
            IF (p_is_master) WRITE(*,*) "      Incrementing existing tstamp_UB for POINT data"
         ENDIF

         time_i = 1  ! Default time index for POINT data
         year = tstamp_UB(tracer_idx,var_i)%year
         
         IF (p_is_master) THEN
            WRITE(*,*) "      POINT dataset UB result:"
            WRITE(*,*) "        tstamp_UB: ", tstamp_UB(tracer_idx,var_i)%year, tstamp_UB(tracer_idx,var_i)%day, tstamp_UB(tracer_idx,var_i)%sec
            WRITE(*,*) "        time_i: ", time_i
            WRITE(*,*) "        year: ", year
         ENDIF
         
         RETURN
      ENDIF

      ! calculate the time stamp
      IF ( tstamp_UB(tracer_idx,var_i)%year == -1 ) THEN
         IF (p_is_master) WRITE(*,*) "      Setting initial tstamp_UB"
         tstamp_UB(tracer_idx,var_i) = tstamp_LB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
      ELSE
         IF (p_is_master) WRITE(*,*) "      Updating existing tstamp_UB"
         tstamp_LB(tracer_idx,var_i) = tstamp_UB(tracer_idx,var_i)
         tstamp_UB(tracer_idx,var_i) = tstamp_UB(tracer_idx,var_i) + DEF_Tracer_Forcings_NL(tracer_idx)%dtime(var_i)
      ENDIF

      ! calculate initial year, day, and second values
      year = tstamp_UB(tracer_idx,var_i)%year
      day  = tstamp_UB(tracer_idx,var_i)%day
      sec  = tstamp_UB(tracer_idx,var_i)%sec

      IF (p_is_master) THEN
         WRITE(*,*) "      Calculated from tstamp_UB:"
         WRITE(*,*) "        year: ", year, " day: ", day, " sec: ", sec
      ENDIF

      IF ( trim(DEF_Tracer_Forcings_NL(tracer_idx)%groupby) == 'year' ) THEN
         IF (p_is_master) WRITE(*,*) "      Group by YEAR - processing year-based grouping"

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
         IF (p_is_master) WRITE(*,*) "      Group by MONTH - processing month-based grouping"

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
         IF (p_is_master) WRITE(*,*) "      Group by DAY - processing day-based grouping"
         
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

      IF (p_is_master) THEN
         WRITE(*,*) "      Final results:"
         WRITE(*,*) "        tstamp_UB: ", tstamp_UB(tracer_idx,var_i)%year, tstamp_UB(tracer_idx,var_i)%day, tstamp_UB(tracer_idx,var_i)%sec
         WRITE(*,*) "        time_i: ", time_i
         WRITE(*,*) "        year: ", year, " month: ", month, " mday: ", mday
         WRITE(*,*) "    === SETSTAMPUB_TRACER DEBUG: EXIT ==="
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
      CALL flush_block_data (traceravgcos, 0._r8)

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
                  traceravgcos%blk(ib,jb)%val(i,j) = traceravgcos%blk(ib,jb)%val(i,j) &
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
      INTEGER :: i
      
      DO i = 1, LEN_TRIM(str)
         IF (str(i:i) >= 'A' .AND. str(i:i) <= 'Z') THEN
            str(i:i) = CHAR(ICHAR(str(i:i)) + 32)
         ENDIF
      ENDDO
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
   
         IF (p_is_master) THEN
            WRITE(*,*) "      === DEBUG: tracerfilename called ==="
            WRITE(*,*) "        Input: year=", year, " month=", month, " day=", day, " var_i=", var_i
         ENDIF
   
         write(yearstr, '(I4.4)') year
         write(monthstr, '(I2.2)') month
         
         dataset_name_lower = TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name)
         CALL to_lower(dataset_name_lower)
   
         IF (p_is_master) THEN
            WRITE(*,*) "        Formatted strings: yearstr='", TRIM(yearstr), "' monthstr='", TRIM(monthstr), "'"
            WRITE(*,*) "        Dataset name (original): '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name), "'"
            WRITE(*,*) "        Dataset name (lower): '", dataset_name_lower, "'"
         ENDIF
   
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
      
               metfilename = '/'//trim(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(var_i))//'_'//trim(yearstr)//'.nc'
               IF (p_is_master) THEN
                  WRITE(*,*) "        CASE IsoGSM selected"
                  WRITE(*,*) "        fprefix_tracer_forcing(", var_i, ") = '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(var_i)), "'"
                  WRITE(*,*) "        Generated filename: '", TRIM(metfilename), "'"
               ENDIF
   
         
         CASE ('POINT')
            metfilename = '/'//trim(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(1))
            IF (p_is_master) THEN
               WRITE(*,*) "        CASE POINT selected"
               WRITE(*,*) "        fprefix_tracer_forcing(1) = '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(1)), "'"
               WRITE(*,*) "        Generated filename: '", TRIM(metfilename), "'"
            ENDIF
            
         CASE DEFAULT
            IF (p_is_master) THEN
               WRITE(*,*) "        WARNING: Unknown dataset name '", TRIM(DEF_Tracer_Forcings_NL(tracer_idx)%dataset_name), "'"
               WRITE(*,*) "        Using default POINT format"
            ENDIF
            metfilename = '/'//trim(DEF_Tracer_Forcings_NL(tracer_idx)%fprefix(1))
         END select
         
         IF (p_is_master) THEN
            WRITE(*,*) "      === DEBUG: tracerfilename returning '", TRIM(metfilename), "' ==="
         ENDIF
         
         ! IF (DEF_USE_CBL_HEIGHT) THEN
         !    select CASE (var_i)
         !    CASE (9)
         !       metfilename = '/'//trim(fprefix_tracer_forcing(9))//'_'//trim(yearstr)//'_'//trim(monthstr)//&
         !          '_boundary_layer_height.nc4'
         !    END select
         ! ENDIF
      END FUNCTION tracerfilename


      SUBROUTINE tracerpreprocess(grid, forcn)

         USE MOD_Const_Physical
         USE MOD_Namelist
         USE MOD_SPMD_Task
         USE MOD_Block
         USE MOD_Grid
         USE MOD_DataType
         USE MOD_Qsadv
         IMPLICIT NONE
         type(grid_type), intent(in) :: grid
         type(block_data_real8_2d), intent(inout) :: forcn(:)
      
         integer  :: iblkme, ib, jb, i, j
         real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT, e, ea
      
            !----------------------------------------------------------------------------
            ! use polynomials to calculate saturation vapor pressure and derivative with
            ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
            ! required to convert relative humidity to specific humidity
            !----------------------------------------------------------------------------
            IF (trim(DEF_Tracer_Forcings_NL(1)%dataset_name) == 'POINT') THEN
#ifdef SinglePoint
               CALL qsadv(forcn(1)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1), &
                          forcn(3)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1), &
                          es,esdT,qsat_tmp,dqsat_tmpdT)
               IF (qsat_tmp < forcn(2)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1)) THEN
                  forcn(2)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1) = qsat_tmp
               ENDIF
#endif
            ELSE
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)
      
                  DO j = 1, grid%ycnt(jb)
                     DO i = 1, grid%xcnt(ib)
      
                        select CASE (trim(DEF_Tracer_Forcings_NL(1)%dataset_name))
      
      
                        CASE ('IsoGSM') ! IsoGSM forcing
                         !  CALL qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                         !     es,esdT,qsat_tmp,dqsat_tmpdT)
                         !    IF (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) THEN
                         !       forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                         !    ENDIF
                           print *, 'IsoGSM forcing is not implemented yet'
      
                        END select
      
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
      
         END SUBROUTINE tracerpreprocess

END MODULE MOD_Tracer_Forcing
