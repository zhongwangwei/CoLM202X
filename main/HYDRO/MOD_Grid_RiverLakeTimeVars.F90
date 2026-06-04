#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeTimeVars
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in gridded hydrological processes.
!
! Created by Shupeng Zhang, Oct 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
#ifdef TRACER
   USE MOD_Tracer_Particle, only: tracer_particle_write_restart
#endif
   USE MOD_Grid_RiverLakeBifurcation, only: write_bifurcation_restart
   IMPLICIT NONE

   real(r8), parameter :: RIVERLAKE_DRY_DEPTH = 1.e-5_r8

   ! -- state variables --
   real(r8), allocatable :: wdsrf_ucat (:) ! river or lake water depth [m]
   real(r8), allocatable :: veloc_riv  (:) ! river velocity            [m/s]
   real(r8), allocatable :: momen_riv  (:) ! unit river momentum       [m^2/s]
   real(r8), allocatable :: volresv    (:) ! reservoir water volume    [m^3]
   ! Prognostic routing volume [m^3]. For ordinary cells this is the total
   ! movable river/floodplain volume; for levee cells it is the visible
   ! river-side volume and levsto carries the protected-side volume.  This is
   ! owned by TimeVars because it is a hydrologic state/restart variable, not a
   ! levee-parameter variable; the levee module may read/write it through
   ! explicit arguments only.
   real(r8), allocatable :: volwater_ucat (:)
   logical               :: volwater_ucat_valid = .false.

   ! -- routing accumulator state (must persist across restart so that a
   !    restart written mid routing-period does not drop the land-runoff
   !    increments queued for the next flush). Owned here (rather than in
   !    MOD_Grid_RiverLakeFlow) so WRITE/READ_GridRiverLakeTimeVars can
   !    serialise them without a circular USE; Flow imports them via the
   !    existing USE MOD_Grid_RiverLakeTimeVars at the top of that module.
   real(r8), save       :: acctime_rnof = 0._r8     ! accumulated land time since last routing flush [s]
   real(r8), allocatable :: acc_rnof_uc (:)         ! accumulated runoff depth per ucatch [m]

   ! -- restart file path (saved for deferred particle-tracer restart read) --
   character(len=512) :: gridriver_restart_file = ''

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_GridRiverLakeTimeVars
   PUBLIC :: deallocate_GridRiverLakeTimeVars

   PUBLIC :: read_GridRiverLakeTimeVars
   PUBLIC :: write_GridRiverLakeTimeVars

CONTAINS

   SUBROUTINE allocate_GridRiverLakeTimeVars

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

   integer :: ncell_state

      ncell_state = 0
      IF (p_is_worker) ncell_state = numucat

      IF (p_is_worker) THEN

         IF (numucat > 0)  allocate (wdsrf_ucat (numucat))
         IF (numucat > 0)  allocate (veloc_riv  (numucat))
         IF (numucat > 0)  allocate (momen_riv  (numucat))
         IF (numresv > 0)  allocate (volresv    (numresv))

         ! Allocate on ALL workers (zero-length if numucat=0) so assumed-
         ! shape MPI wrappers downstream never see an unallocated array.
         allocate (acc_rnof_uc (numucat))
         acc_rnof_uc = 0._r8

      ENDIF

      ! Allocated on every rank because read_levee_restart is collective and
      ! receives this TimeVars-owned state via an assumed-shape argument.
      allocate (volwater_ucat(ncell_state))
      volwater_ucat = 0._r8

      acctime_rnof = 0._r8
      volwater_ucat_valid = .false.

   END SUBROUTINE allocate_GridRiverLakeTimeVars


   SUBROUTINE READ_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist,              only: DEF_Reservoir_Method
   USE MOD_Vector_ReadWrite
   USE MOD_NetCDFSerial,          only: ncio_var_exist, ncio_read_bcast_serial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, resv_data_address, totalnumresv
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical :: has_var

      gridriver_restart_file = trim(file_restart)

      CALL vector_read_and_scatter (file_restart, wdsrf_ucat, numucat, 'wdsrf_ucat', ucat_data_address)
      CALL vector_read_and_scatter (file_restart, veloc_riv,  numucat, 'veloc_riv',  ucat_data_address)

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN
            CALL vector_read_and_scatter (file_restart, volresv, numresv, 'volresv', resv_data_address)
         ENDIF
      ENDIF

      ! Routing accumulator recovery. Absent variables (old-format restart
      ! or first cold start) keep the zero defaults from allocate_*, so the
      ! old behaviour is preserved when the file predates this persistence.
      has_var = .false.
      IF (p_is_master) has_var = ncio_var_exist(file_restart, 'acctime_rnof', readflag = .false.)
#ifdef USEMPI
      CALL mpi_bcast (has_var, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      IF (has_var) THEN
         CALL ncio_read_bcast_serial (file_restart, 'acctime_rnof', acctime_rnof)
      ENDIF

      has_var = .false.
      IF (p_is_master) has_var = ncio_var_exist(file_restart, 'acc_rnof_uc', readflag = .false.)
#ifdef USEMPI
      CALL mpi_bcast (has_var, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      IF (has_var) THEN
         CALL vector_read_and_scatter (file_restart, acc_rnof_uc, numucat, 'acc_rnof_uc', ucat_data_address)
      ENDIF

      has_var = .false.
      IF (p_is_master) has_var = ncio_var_exist(file_restart, 'volwater_ucat', readflag = .false.)
#ifdef USEMPI
      CALL mpi_bcast (has_var, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      IF (has_var) THEN
         CALL vector_read_and_scatter (file_restart, volwater_ucat, numucat, 'volwater_ucat', ucat_data_address)
         volwater_ucat_valid = .true.
      ENDIF

      ! Note: levee restart (levsto) is read separately in grid_riverlake_flow_init
      ! after levee_init() allocates the levsto array. Same deferred-restart pattern.

      ! Note: bifurcation restart (pth_veloc, pth_momen) is read separately
      ! in grid_riverlake_flow_init after bifurcation_init() allocates the
      ! arrays. Same deferred-restart pattern as levee.

      ! Note: tracer restart is read separately in grid_riverlake_flow_init
      ! after river_lake_tracer_init() allocates the arrays. Same pattern as others.

      ! Note: particle-tracer restart is read through MOD_Tracer_Particle,
      ! called from grid_riverlake_flow_init after particle species are initialized.

   END SUBROUTINE READ_GridRiverLakeTimeVars


   SUBROUTINE WRITE_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist,              only: DEF_Reservoir_Method, DEF_USE_LEVEE, DEF_USE_BIFURCATION
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_data_address
   USE MOD_Grid_RiverLakeLevee,   only: levsto
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'ucatch', totalnumucat)
      ENDIF

      CALL vector_gather_and_write (&
         wdsrf_ucat, numucat, totalnumucat, ucat_data_address, file_restart, 'wdsrf_ucat', 'ucatch')

      CALL vector_gather_and_write (&
         veloc_riv, numucat, totalnumucat, ucat_data_address, file_restart, 'veloc_riv', 'ucatch')

      ! Persist the routing accumulator so a restart written mid routing
      ! period does not lose the queued land-runoff increments. Paired
      ! with the acctime_rnof / acc_rnof_uc reads in READ_*.
      IF (p_is_master) CALL ncio_write_serial (file_restart, 'acctime_rnof', acctime_rnof)
      CALL vector_gather_and_write (&
         acc_rnof_uc, numucat, totalnumucat, ucat_data_address, file_restart, 'acc_rnof_uc', 'ucatch')

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (p_is_master) CALL ncio_define_dimension(file_restart, 'reservoir', totalnumresv)

            CALL vector_gather_and_write (&
               volresv, numresv, totalnumresv, resv_data_address, file_restart, 'volresv', 'reservoir')
         ENDIF
      ENDIF

      CALL vector_gather_and_write (&
         volwater_ucat, numucat, totalnumucat, ucat_data_address, file_restart, 'volwater_ucat', 'ucatch')

      IF (DEF_USE_LEVEE) THEN
         CALL vector_gather_and_write (&
            levsto, numucat, totalnumucat, ucat_data_address, file_restart, 'levsto', 'ucatch')
      ENDIF

      IF (DEF_USE_BIFURCATION) THEN
         CALL write_bifurcation_restart(file_restart)
      ENDIF

#ifdef TRACER
      CALL tracer_particle_write_restart(file_restart)
#endif

   END SUBROUTINE WRITE_GridRiverLakeTimeVars

   SUBROUTINE deallocate_GridRiverLakeTimeVars

   IMPLICIT NONE

      IF (allocated (wdsrf_ucat )) deallocate (wdsrf_ucat )
      IF (allocated (veloc_riv  )) deallocate (veloc_riv  )
      IF (allocated (momen_riv  )) deallocate (momen_riv  )
      IF (allocated (volresv    )) deallocate (volresv    )
      IF (allocated (volwater_ucat)) deallocate (volwater_ucat)
      IF (allocated (acc_rnof_uc)) deallocate (acc_rnof_uc)
      volwater_ucat_valid = .false.

   END SUBROUTINE deallocate_GridRiverLakeTimeVars

END MODULE MOD_Grid_RiverLakeTimeVars
#endif
