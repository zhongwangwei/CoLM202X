#include <define.h>

MODULE MOD_Tracer_Vars_Global
#ifdef USE_TRACER
!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!  Define some global variables for isotope calculations
!
!  Zhongwang Wei, 03/2025: initial version
!
!-----------------------------------------------------------------------
! !USES:
   USE MOD_Precision
   USE MOD_Namelist
   IMPLICIT NONE
   SAVE
   ! Isotope ratio constants
   real(r8), parameter :: Rsmow_O = 2005.2e-6_r8    ! Standard Mean Ocean Water ratio for O18/O16
   real(r8), parameter :: Rsmow_D = 155.76e-6_r8    ! Standard Mean Ocean Water ratio for D/H

   ! Diffusion coefficients
   real(r8), public :: dif18o  ! Diffusion coefficient for O18
   real(r8), public :: difhdo  ! Diffusion coefficient for HDO
   real(r8), public :: nn      ! Power law coefficient

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_Init_GlobalVars

CONTAINS

   SUBROUTINE Tracer_Init_GlobalVars
      IMPLICIT NONE
      ! Initialize diffusion coefficients
      dif18o = 1.02849_r8
      difhdo = 1.02512_r8
      nn     = 0.58_r8
   END SUBROUTINE Tracer_Init_GlobalVars

#endif
END MODULE MOD_Tracer_Vars_Global
! ---------- EOP ------------
