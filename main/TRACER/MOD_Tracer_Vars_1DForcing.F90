#include <define.h>

MODULE MOD_Tracer_Vars_1DForcing
!-----------------------------------------------------------------------
!  Meteorological Forcing
!
!  Created by Yongjiu Dai, 03/2014
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
   real(r8), allocatable, target :: forc_q_O18(:)    ! atmospheric isotope vapor data [kg/m/s]
   real(r8), allocatable, target :: forc_q_H2(:)     ! atmospheric isotope vapor data [kg/m/s]
   real(r8), allocatable, target :: forc_prc_O18(:)  ! convective precipitation [kg/m/s]
   real(r8), allocatable, target :: forc_prl_O18(:)  ! large scale precipitation [kg/m/s]
   real(r8), allocatable, target :: forc_prc_H2(:)   ! convective precipitation [kg/m/s]
   real(r8), allocatable, target :: forc_prl_H2(:)   ! large scale precipitation [kg/m/s]
   real(r8), allocatable, target :: forc_rain_O18(:) ! rain [kg/m/s]
   real(r8), allocatable, target :: forc_snow_O18(:) ! snow [kg/m/s]
   real(r8), allocatable, target :: forc_rain_H2(:)  ! rain [kg/m/s]
   real(r8), allocatable, target :: forc_snow_H2(:)  ! snow [kg/m/s]
   real(r8), allocatable, target :: tracer_forc(:,:,:) ! tracer forcing [kg/m/s]
! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_1D_Forcing_O18
   PUBLIC :: Tracer_allocate_1D_Forcing_H2
   PUBLIC :: Tracer_deallocate_1D_Forcing_O18
   PUBLIC :: Tracer_deallocate_1D_Forcing_H2

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------



   SUBROUTINE Tracer_allocate_1D_Forcing_O18
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM 1d [numpatch] variables
   ! -------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_Mesh
   USE MOD_LandPatch
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            ! O18
            allocate (forc_q_O18(numpatch) ) ! atmospheric isotope vapor data [kg/m/s]
            allocate (forc_prc_O18(numpatch) ) ! convective precipitation [kg/m/s]
            allocate (forc_prl_O18(numpatch) ) ! large scale precipitation [kg/m/s]
            allocate (forc_rain_O18(numpatch) ) ! rain [mm/s]
            allocate (forc_snow_O18(numpatch) ) ! snow [mm/s]
         ENDIF
      ENDIF
   END SUBROUTINE Tracer_allocate_1D_Forcing_O18


   SUBROUTINE Tracer_allocate_1D_Forcing_H2
      ! -------------------------------------------------------------------
      ! Allocates memory for CoLM 1d [numpatch] variables
      ! -------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE MOD_LandPatch
      IMPLICIT NONE
   
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               ! H2
               allocate (forc_q_H2(numpatch) ) ! atmospheric isotope vapor data [kg/m/s]
               allocate (forc_prc_H2(numpatch) ) ! convective precipitation [kg/m/s]
               allocate (forc_prl_H2(numpatch) ) ! large scale precipitation [kg/m/s]
               allocate (forc_rain_H2(numpatch) ) ! rain [mm/s]
               allocate (forc_snow_H2(numpatch) ) ! snow [mm/s]
            ENDIF
   
         ENDIF
   
   END SUBROUTINE Tracer_allocate_1D_Forcing_H2



   SUBROUTINE Tracer_deallocate_1D_Forcing_O18 ()

   USE MOD_SPMD_Task
   USE MOD_Mesh
   USE MOD_LandPatch
   IMPLICIT NONE

   IF (p_is_worker) THEN
      IF (numpatch > 0) THEN
         deallocate ( forc_q_O18 ) ! atmospheric isotope vapor data [kg/m/s]
         deallocate ( forc_prc_O18 ) ! convective precipitation [kg/m/s]
         deallocate ( forc_prl_O18 ) ! large scale precipitation [kg/m/s]
         deallocate ( forc_rain_O18   ) ! rain [mm/s]
         deallocate ( forc_snow_O18   ) ! snow [mm/s]
      ENDIF
   ENDIF

   END SUBROUTINE Tracer_deallocate_1D_Forcing_O18

   SUBROUTINE Tracer_deallocate_1D_Forcing_H2 ()
      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE MOD_LandPatch
      IMPLICIT NONE
   
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               deallocate ( forc_q_H2 ) ! atmospheric isotope vapor data [kg/m/s]
               deallocate ( forc_prc_H2 ) ! convective precipitation [kg/m/s]
               deallocate ( forc_prl_H2 ) ! large scale precipitation [kg/m/s]
               deallocate ( forc_rain_H2   ) ! rain [mm/s]
               deallocate ( forc_snow_H2   ) ! snow [mm/s]
            ENDIF
         ENDIF
   END SUBROUTINE Tracer_deallocate_1D_Forcing_H2


   END MODULE MOD_Tracer_Vars_1DForcing
! ---------- EOP ------------
