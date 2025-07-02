#include <define.h>

MODULE MOD_Tracer_Vars_2DForcing
!-----------------------------------------------------------------------
!  Meteorogical Forcing
!
!  Created by Yongjiu Dai, 03/2014
!-----------------------------------------------------------------------

   USE MOD_DataType
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
   type(block_data_real8_2d) :: forc_xy_q_O18 ! atmospheric isotope precipitation data [kg/m/s]
   type(block_data_real8_2d) :: forc_xy_prc_O18    ! convective precipitation [mm/s]
   type(block_data_real8_2d) :: forc_xy_prl_O18    ! large scale precipitation [mm/s]

   type(block_data_real8_2d) :: forc_xy_q_H2 ! atmospheric isotope precipitation data [kg/m/s]
   type(block_data_real8_2d) :: forc_xy_prc_H2    ! convective precipitation [mm/s]
   type(block_data_real8_2d) :: forc_xy_prl_H2    ! large scale precipitation [mm/s]
   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_2D_Forcing

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_2D_Forcing (grid)
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM 2d [lon_points,lat_points] variables
   ! -------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   type(grid_type), intent(in) :: grid
      IF (p_is_io) THEN
         CALL allocate_block_data (grid, forc_xy_prc_O18 ) ! convective precipitation [mm/s]
         CALL allocate_block_data (grid, forc_xy_prl_O18 ) ! large scale precipitation [mm/s]
         CALL allocate_block_data (grid, forc_xy_q_O18 ) ! atmospheric isotope precipitation data [kg/m/s]

         CALL allocate_block_data (grid, forc_xy_prc_H2 ) ! convective precipitation [mm/s]
         CALL allocate_block_data (grid, forc_xy_prl_H2 ) ! large scale precipitation [mm/s]
         CALL allocate_block_data (grid, forc_xy_q_H2 ) ! atmospheric isotope precipitation data [kg/m/s]
      ENDIF

   END SUBROUTINE Tracer_allocate_2D_Forcing

END MODULE MOD_Tracer_Vars_2DForcing
! ---------- EOP ------------
