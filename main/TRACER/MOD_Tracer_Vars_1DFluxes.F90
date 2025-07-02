#include <define.h>

MODULE MOD_Tracer_Vars_1DFluxes
!-----------------------------------------------------------------------
!  Created by Yongjiu Dai, 03/2014
!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
!  Fluxes
!-----------------------------------------------------------------------
   real(r8), allocatable :: fevpa_O18  (:) !evapotranspiration from canopy to atmosphere [mm/s]
   real(r8), allocatable :: fevpa_H2  (:) !evapotranspiration from canopy to atmosphere [mm/s]
   real(r8), allocatable :: fevpl_O18  (:) !evapotranspiration from leaves to atmosphere [mm/s]
   real(r8), allocatable :: fevpl_H2  (:) !evapotranspiration from leaves to atmosphere [mm/s]
   real(r8), allocatable :: etr_O18  (:) !transpiration rate [mm/s]
   real(r8), allocatable :: etr_H2  (:) !transpiration rate [mm/s]
   real(r8), allocatable :: fevpg_O18  (:) !evaporation heat flux from ground [mm/s]
   real(r8), allocatable :: fevpg_H2  (:) !evaporation heat flux from ground [mm/s]
   real(r8), allocatable :: rsur_O18  (:) !surface runoff (mm h2o/s)
   real(r8), allocatable :: rsur_se_O18  (:) !saturation excess surface runoff (mm h2o/s)
   real(r8), allocatable :: rsur_ie_O18  (:) !infiltration excess surface runoff (mm h2o/s)
   real(r8), allocatable :: rsub_O18  (:) !subsurface runoff (mm h2o/s)
   real(r8), allocatable :: rnof_O18  (:) !total runoff (mm h2o/s)
   real(r8), allocatable :: qintr_O18  (:) !interception (mm h2o/s)
   real(r8), allocatable :: rsur_H2  (:) !surface runoff (mm h2o/s)
   real(r8), allocatable :: rsur_se_H2  (:) !saturation excess surface runoff (mm h2o/s)
   real(r8), allocatable :: rsur_ie_H2  (:) !infiltration excess surface runoff (mm h2o/s)
   real(r8), allocatable :: rsub_H2  (:) !subsurface runoff (mm h2o/s)
   real(r8), allocatable :: rnof_H2  (:) !total runoff (mm h2o/s)
   real(r8), allocatable :: qintr_H2  (:) !interception (mm h2o/s)
   real(r8), allocatable :: qinfl_O18  (:) !infiltration (mm h2o/s)
   real(r8), allocatable :: qinfl_H2  (:) !infiltration (mm h2o/s)
   real(r8), allocatable :: qdrip_O18  (:) !throughfall (mm h2o/s)
   real(r8), allocatable :: qdrip_H2  (:) !throughfall (mm h2o/s)
   real(r8), allocatable :: qcharge_O18(:) !groundwater recharge [mm/s]
   real(r8), allocatable :: qcharge_H2(:) !groundwater recharge [mm/s]

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_1D_Fluxes
   PUBLIC :: Tracer_deallocate_1D_Fluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------
CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_1D_Fluxes
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM 1d [numpatch] variables
   ! -------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN

            allocate ( fevpa_O18  (numpatch) )  ; fevpa_O18  (:) = spval ! evapotranspiration from canopy to atmosphere [mm/s]
            allocate ( fevpa_H2  (numpatch) )  ; fevpa_H2  (:) = spval ! evapotranspiration from canopy to atmosphere [mm/s]
            allocate ( fevpl_O18  (numpatch) )  ; fevpl_O18  (:) = spval ! evapotranspiration from leaves to atmosphere [mm/s]
            allocate ( fevpl_H2  (numpatch) )  ; fevpl_H2  (:) = spval ! evapotranspiration from leaves to atmosphere [mm/s]
            allocate ( etr_O18  (numpatch) )  ; etr_O18  (:) = spval ! transpiration rate [mm/s]
            allocate ( etr_H2  (numpatch) )  ; etr_H2  (:) = spval ! transpiration rate [mm/s]
            allocate ( fevpg_O18  (numpatch) )  ; fevpg_O18  (:) = spval ! evaporation heat flux from ground [mm/s]
            allocate ( fevpg_H2  (numpatch) )  ; fevpg_H2  (:) = spval ! evaporation heat flux from ground [mm/s]

            allocate ( rsur_O18  (numpatch) )  ; rsur_O18  (:) = spval ! surface runoff (mm h2o/s)
            allocate ( rsur_se_O18  (numpatch) )  ; rsur_se_O18  (:) = spval ! saturation excess surface runoff (mm h2o/s)
            allocate ( rsur_ie_O18  (numpatch) )  ; rsur_ie_O18  (:) = spval ! infiltration excess surface runoff (mm h2o/s)
            allocate ( rsub_O18  (numpatch) )  ; rsub_O18  (:) = spval ! subsurface runoff (mm h2o/s)
            allocate ( rnof_O18  (numpatch) )  ; rnof_O18  (:) = spval ! total runoff (mm h2o/s)
            allocate ( qintr_O18  (numpatch) )  ; qintr_O18  (:) = spval ! interception (mm h2o/s)
            allocate ( rsur_H2  (numpatch) )  ; rsur_H2  (:) = spval ! surface runoff (mm h2o/s)
            allocate ( rsur_se_H2  (numpatch) )  ; rsur_se_H2  (:) = spval ! saturation excess surface runoff (mm h2o/s)
            allocate ( rsur_ie_H2  (numpatch) )  ; rsur_ie_H2  (:) = spval ! infiltration excess surface runoff (mm h2o/s)
            allocate ( rsub_H2  (numpatch) )  ; rsub_H2  (:) = spval ! subsurface runoff (mm h2o/s)
            allocate ( rnof_H2  (numpatch) )  ; rnof_H2  (:) = spval ! total runoff (mm h2o/s)
            allocate ( qintr_H2  (numpatch) )  ; qintr_H2  (:) = spval ! interception (mm h2o/s)

            allocate ( qinfl_O18  (numpatch) )  ; qinfl_O18  (:) = spval ! infiltration (mm h2o/s)
            allocate ( qinfl_H2  (numpatch) )  ; qinfl_H2  (:) = spval ! infiltration (mm h2o/s)

            allocate ( qdrip_O18  (numpatch) )  ; qdrip_O18  (:) = spval ! throughfall (mm h2o/s)
            allocate ( qdrip_H2  (numpatch) )  ; qdrip_H2  (:) = spval ! throughfall (mm h2o/s)

            allocate ( qcharge_O18(numpatch) )  ; qcharge_O18(:) = spval ! groundwater recharge [mm/s]
            allocate ( qcharge_H2(numpatch) )  ; qcharge_H2(:) = spval ! groundwater recharge [mm/s]

         ENDIF
      ENDIF



   END SUBROUTINE Tracer_allocate_1D_Fluxes

   SUBROUTINE Tracer_deallocate_1D_Fluxes ()
   ! --------------------------------------------------------------------
   ! deallocates memory for CoLM 1d [numpatch] variables
   ! --------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_LandPatch

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            deallocate ( fevpa_O18   )  ! evapotranspiration from canopy to atmosphere [mm/s]
            deallocate ( fevpa_H2   )  ! evapotranspiration from canopy to atmosphere [mm/s]
            deallocate ( fevpl_O18   )  ! evapotranspiration from leaves to atmosphere [mm/s]
            deallocate ( fevpl_H2   )  ! evapotranspiration from leaves to atmosphere [mm/s]
            deallocate ( etr_O18   )  ! transpiration rate [mm/s]
            deallocate ( etr_H2   )  ! transpiration rate [mm/s]
            deallocate ( fevpg_O18   )  ! evaporation heat flux from ground [mm/s]
            deallocate ( fevpg_H2   )  ! evaporation heat flux from ground [mm/s]
            deallocate ( rsur_O18   )  ! surface runoff (mm h2o/s)
            deallocate ( rsur_se_O18   )  ! saturation excess surface runoff (mm h2o/s)
            deallocate ( rsur_ie_O18   )  ! infiltration excess surface runoff (mm h2o/s)
            deallocate ( rsub_O18   )  ! subsurface runoff (mm h2o/s)
            deallocate ( rnof_O18   )  ! total runoff (mm h2o/s)
            deallocate ( qintr_O18   )  ! interception (mm h2o/s)
            deallocate ( rsur_H2   )  ! surface runoff (mm h2o/s)
            deallocate ( rsur_se_H2   )  ! saturation excess surface runoff (mm h2o/s)
            deallocate ( rsur_ie_H2   )  ! infiltration excess surface runoff (mm h2o/s)
            deallocate ( rsub_H2   )  ! subsurface runoff (mm h2o/s)
            deallocate ( rnof_H2   )  ! total runoff (mm h2o/s)
            deallocate ( qintr_H2   )  ! interception (mm h2o/s)
            deallocate ( qinfl_O18   )  ! infiltration (mm h2o/s)
            deallocate ( qinfl_H2   )  ! infiltration (mm h2o/s)

            deallocate ( qdrip_O18   )  ! interception (mm h2o/s)
            deallocate ( qdrip_H2   )  ! interception (mm h2o/s)

            deallocate ( qcharge_O18   )  ! groundwater recharge [mm/s]
            deallocate ( qcharge_H2   )  ! groundwater recharge [mm/s]
         ENDIF
      ENDIF


   END SUBROUTINE Tracer_deallocate_1D_Fluxes

END MODULE MOD_Tracer_Vars_1DFluxes
! ---------- EOP ------------
