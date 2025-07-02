#include <define.h>

!-----------------------------------------------------------------------
! Created by Yongjiu Dai, 03/2014
!-----------------------------------------------------------------------

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
MODULE MOD_Tracer_Vars_PFTimeInvariants
!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Define PFT time invariables
!
!  Added by Hua Yuan, 08/2019
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
   SAVE

  

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_PFTimeInvariants
   PUBLIC :: Tracer_READ_PFTimeInvariants
   PUBLIC :: Tracer_WRITE_PFTimeInvariants
   PUBLIC :: Tracer_deallocate_PFTimeInvariants
#ifdef RangeCheck
   PUBLIC :: Tracer_check_PFTimeInvariants
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_PFTimeInvariants
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM PFT 1d [numpft] variables
   ! -------------------------------------------------------------------

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_LandPFT,   only: numpft
   USE MOD_Precision
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN

         ENDIF


      ENDIF

   END SUBROUTINE Tracer_allocate_PFTimeInvariants

   SUBROUTINE Tracer_READ_PFTimeInvariants (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandPatch
   USE MOD_LandPFT
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

     ! CALL ncio_read_vector (file_restart, 'pftclass', landpft, pftclass) !
     ! CALL ncio_read_vector (file_restart, 'pftfrac ', landpft, pftfrac ) !
     ! CALL ncio_read_vector (file_restart, 'htop_p  ', landpft, htop_p  ) !
    !  CALL ncio_read_vector (file_restart, 'hbot_p  ', landpft, hbot_p  ) !

   END SUBROUTINE Tracer_READ_PFTimeInvariants

   SUBROUTINE Tracer_WRITE_PFTimeInvariants (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandPFT
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_Vars_Global
   IMPLICIT NONE

   ! Local variables
   character(len=*), intent(in) :: file_restart
   integer :: compress

      compress = DEF_REST_CompressLevel

     ! CALL ncio_create_file_vector (file_restart, landpft)
      !CALL ncio_define_dimension_vector (file_restart, landpft, 'pft')

     !  CALL ncio_write_vector (file_restart, 'pftclass', 'pft', landpft, pftclass, compress) !
     !  CALL ncio_write_vector (file_restart, 'pftfrac ', 'pft', landpft, pftfrac , compress) !
     !  CALL ncio_write_vector (file_restart, 'htop_p  ', 'pft', landpft, htop_p  , compress) !
     !  CALL ncio_write_vector (file_restart, 'hbot_p  ', 'pft', landpft, hbot_p  , compress) !



   END SUBROUTINE Tracer_WRITE_PFTimeInvariants

   SUBROUTINE Tracer_deallocate_PFTimeInvariants
   ! -------------------------------------------------------------------
   ! Deallocates memory for CoLM PFT 1d [numpft] variables
   ! -------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_LandPFT

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN
     !       deallocate (pftclass)
         ENDIF
      ENDIF

   END SUBROUTINE Tracer_deallocate_PFTimeInvariants

#ifdef RangeCheck
   SUBROUTINE Tracer_check_PFTimeInvariants ()

   USE MOD_RangeCheck
   IMPLICIT NONE

     ! CALL check_vector_data ('pftfrac', pftfrac) !
     ! CALL check_vector_data ('htop_p ', htop_p ) !
     ! CALL check_vector_data ('hbot_p ', hbot_p ) !

   END SUBROUTINE Tracer_check_PFTimeInvariants
#endif

END MODULE MOD_Tracer_Vars_PFTimeInvariants
#endif

MODULE MOD_Tracer_Vars_TimeInvariants
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

   USE MOD_Precision

   IMPLICIT NONE
   SAVE

! -----------------------------------------------------------------
! surface classification and soil information

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_TimeInvariants
   PUBLIC :: Tracer_deallocate_TimeInvariants
   PUBLIC :: Tracer_READ_TimeInvariants
   PUBLIC :: Tracer_WRITE_TimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_TimeInvariants ()
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM 1d [numpatch] variables
   ! -------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN

   
         ENDIF
      ENDIF

   END SUBROUTINE Tracer_allocate_TimeInvariants

   !---------------------------------------
   SUBROUTINE Tracer_READ_TimeInvariants (lc_year, casename, dir_restart)

   !====================================================================
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !====================================================================

   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_NetCDFVector
   USE MOD_NetCDFSerial
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_LandPatch
   USE MOD_Vars_Global
   USE MOD_Const_LC, only: patchtypes

   IMPLICIT NONE

   integer         , intent(in) :: lc_year
   character(len=*), intent(in) :: casename
   character(len=*), intent(in) :: dir_restart

   ! Local variables
   character(len=256) :: file_restart, cyear, lndname

      write(cyear,'(i4.4)') lc_year
      file_restart = trim(dir_restart) // '/const/' // trim(casename) //'_tracer_restart_const' // '_lc' // trim(cyear) // '.nc'


#ifdef RangeCheck
      CALL Tracer_check_TimeInvariants ()
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,'(A29)') 'Loading Time Invariants done.'
      ENDIF

   END SUBROUTINE Tracer_READ_TimeInvariants

   !---------------------------------------
   SUBROUTINE Tracer_WRITE_TimeInvariants (lc_year, casename, dir_restart)

   !====================================================================
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !====================================================================

   USE MOD_Namelist, only: DEF_REST_CompressLevel, DEF_USE_BEDROCK
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_LandPatch
   USE MOD_Vars_Global

   IMPLICIT NONE

   integer         , intent(in) :: lc_year
   character(len=*), intent(in) :: casename
   character(len=*), intent(in) :: dir_restart

   ! Local Variables
   character(len=256) :: file_restart, cyear
   integer :: compress

      compress = DEF_REST_CompressLevel

      write(cyear,'(i4.4)') lc_year

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(dir_restart)//'/const')
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      file_restart = trim(dir_restart) // '/const/' // trim(casename) //'_tracer_restart_const' //'_lc'// trim(cyear) // '.nc'


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_master) then

      end if
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif



   END SUBROUTINE Tracer_WRITE_TimeInvariants

   SUBROUTINE Tracer_deallocate_TimeInvariants ()

   USE MOD_Namelist, only: DEF_USE_Forcing_Downscaling
   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE

      ! --------------------------------------------------
      ! Deallocates memory for CoLM 1d [numpatch] variables
      ! --------------------------------------------------

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN

    

         ENDIF
      ENDIF


   END SUBROUTINE Tracer_deallocate_TimeInvariants

#ifdef RangeCheck
   SUBROUTINE Tracer_check_TimeInvariants ()

   USE MOD_SPMD_Task
   USE MOD_RangeCheck
   USE MOD_Namelist, only: DEF_USE_BEDROCK, DEF_USE_Forcing_Downscaling

   IMPLICIT NONE

      real(r8), allocatable :: tmpcheck(:,:)

      IF (p_is_master) THEN
         write(*,'(/,A29)') 'Checking Time Invariants ...'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN

      ENDIF


   END SUBROUTINE Tracer_check_TimeInvariants
#endif

END MODULE MOD_Tracer_Vars_TimeInvariants
! ---------- EOP ------------
