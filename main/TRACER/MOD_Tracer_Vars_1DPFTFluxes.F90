#include <define.h>

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

MODULE MOD_Tracer_Vars_1DPFTFluxes
!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Define PFT flux variables
!
!  Created by Hua Yuan, 08/2019
!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE
!----------------------------------------------------------------------
!  Fluxes
!-----------------------------------------------------------------------

   real(r8), allocatable :: qintr_p_O18(:) !interception (mm h2o/s)
   real(r8), allocatable :: qintr_rain_p_O18(:) !rainfall interception (mm h2o/s)
   real(r8), allocatable :: qintr_snow_p_O18(:) !snowfall interception (mm h2o/s)
   real(r8), allocatable :: qintr_p_H2(:) !interception (mm h2o/s)
   real(r8), allocatable :: qintr_rain_p_H2(:) !rainfall interception (mm h2o/s)
   real(r8), allocatable :: qintr_snow_p_H2(:) !snowfall interception (mm h2o/s)
! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_1D_PFTFluxes
   PUBLIC :: Tracer_deallocate_1D_PFTFluxes
   PUBLIC :: Tracer_set_1D_PFTFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_1D_PFTFluxes
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM PFT 1d [numpft] variables
   ! -------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_LandPFT
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN
            allocate (qintr_p_O18      (numpft)) ; qintr_p_O18      (:) = spval !interception (mm h2o/s)
            allocate (qintr_rain_p_O18 (numpft)) ; qintr_rain_p_O18 (:) = spval !rainfall interception (mm h2o/s)
            allocate (qintr_snow_p_O18 (numpft)) ; qintr_snow_p_O18 (:) = spval !snowfall interception (mm h2o/s)
            allocate (qintr_p_H2      (numpft)) ; qintr_p_H2      (:) = spval !interception (mm h2o/s)
            allocate (qintr_rain_p_H2 (numpft)) ; qintr_rain_p_H2 (:) = spval !rainfall interception (mm h2o/s)
            allocate (qintr_snow_p_H2 (numpft)) ; qintr_snow_p_H2 (:) = spval !snowfall interception (mm h2o/s)
         ENDIF
      ENDIF



   END SUBROUTINE Tracer_allocate_1D_PFTFluxes

   SUBROUTINE Tracer_deallocate_1D_PFTFluxes
   ! -------------------------------------------------------------------
   ! deallocates memory for CoLM PFT 1d [numpft] variables
   ! -------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_LandPFT

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN


            deallocate (qintr_p_O18      )
            deallocate (qintr_rain_p_O18 )
            deallocate (qintr_snow_p_O18 )
            deallocate (qintr_p_H2      )
            deallocate (qintr_rain_p_H2 )
            deallocate (qintr_snow_p_H2 )

         ENDIF
      ENDIF



   END SUBROUTINE Tracer_deallocate_1D_PFTFluxes

   SUBROUTINE Tracer_set_1D_PFTFluxes(Values, Nan)
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM PFT 1d [numpft] variables
   ! -------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_LandPFT
   IMPLICIT NONE

   real(r8),intent(in) :: Values
   real(r8),intent(in) :: Nan

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN

            qintr_p_O18      (:) = Values  !interception (mm h2o/s)
            qintr_rain_p_O18 (:) = Values  !rainfall interception (mm h2o/s)
            qintr_snow_p_O18 (:) = Values  !snowfall interception (mm h2o/s)
            qintr_p_H2      (:) = Values  !interception (mm h2o/s)
            qintr_rain_p_H2 (:) = Values  !rainfall interception (mm h2o/s)
            qintr_snow_p_H2 (:) = Values  !snowfall interception (mm h2o/s)

         ENDIF
      ENDIF


   END SUBROUTINE Tracer_set_1D_PFTFluxes

END MODULE MOD_Tracer_Vars_1DPFTFluxes

#endif
! ---------- EOP ------------
