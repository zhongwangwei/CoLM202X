#include <define.h>

#ifdef LULCC
MODULE MOD_Lulcc_Driver

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccDriver

 CONTAINS


   SUBROUTINE LulccDriver (casename, dir_landdata, dir_restart, &
                           jdate, greenwich)

!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!  the main subroutine for Land use and land cover change simulation
!
!  Created by Hua Yuan, 04/08/2022
!
! !REVISIONS:
!  07/2023, Wenzong Dong: porting to MPI version.
!  08/2023, Wanyi Lin: add interface for Mass&Energy conserved scheme.
!
!-----------------------------------------------------------------------
!
!                  ***** For Development *****
!
!  Extra processes when adding a new variable and #define LULCC:
!
!  1. Save a copy of new variable (if called "var", save it to "var_")
!  with 2 steps:
!
!  1.1 main/LULCC/MOD_Lulcc_Vars_TimeVariables.F90's subroutine
!  "allocate_LulccTimeVariables":
!           allocate (var_(dimension))
!
!  1.2 main/LULCC/MOD_Lulcc_Vars_TimeVariables.F90's subroutine
!  "SAVE_LulccTimeVariables":
!           var_ = var
!
!  2. Reassignment for the next year
!
!  2.1 if used Same Type Assignment (SAT) scheme for variable recovery
!  main/LULCC/MOD_Lulcc_Vars_TimeVariables.F90's subroutine
!  "REST_LulccTimeVariables"
!           var(np) = var_(np_)
!
!  2.2 if using Mass and Energy conservation (MEC) scheme for variable
!  recovery
!
!  2.2.1 main/LULCC/MOD_Lulcc_Vars_TimeVariables.F90's subroutine
!  "REST_LulccTimeVariables":
!           var(np) = var_(np_)
!
!  2.2.2 [No need for PFT/PC scheme] Mass and Energy conserve
!  adjustment, add after line 519 of MOD_Lulcc_MassEnergyConserve.F90.
!
!  o if variable should be mass conserved:
!       var(:,np) = var(:,np) + &
!       var(:,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
!
!  o if variable should be energy conserved, take soil temperature
!  "t_soisno" as an example: [May neeed extra calculation]
!       t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) + &
!                t_soisno_(1:nl_soil,frnp_(k))*cvsoil_(1:nl_soil,k)* &
!                lccpct_np(patchclass_(frnp_(k)))/wgt(1:nl_soil)
!  where cvsoil_ is the heat capacity, wgt is the sum of
!  cvsoil_(1:nl_soil,k)*lccpct_np(patchclass_(frnp_(k))), which need to
!  be calculated in advance.
!
!  3. Deallocate the copy of new variable in
!  MOD_Lulcc_Vars_TimeVariables.F90's subroutine
!  "deallocate_LulccTimeVariables":
!           deallocate (var_)
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_Lulcc_Initialize
   USE MOD_Vars_TimeVariables
   USE MOD_Lulcc_TransferTraceReadin
   USE MOD_Lulcc_MassEnergyConserve
   USE MOD_Namelist
   USE MOD_LandPatch,          only: landpatch
   USE MOD_Vars_TimeInvariants, only: patchclass, patchtype, lake_soilc_srf
#if (defined TRACER) && (defined BGC)
   USE MOD_Tracer_Methane_Const, only: DEF_METHANE
   USE MOD_Tracer_Methane_Registry, only: igas_ch4
   USE MOD_Tracer_Methane_State, only: save_methane_lulcc_state, remap_methane_lulcc_state, &
      initialize_methane_lake_soilc_from_surface
   USE MOD_Tracer_Methane_AccFlux, only: deallocate_methane_acc_fluxes, allocate_methane_acc_fluxes
   USE MOD_Tracer_Methane_Microbes, only: save_methane_microbes_lulcc_state, remap_methane_microbes_lulcc_state
#endif

   IMPLICIT NONE

   character(len=256), intent(in) :: casename      !casename name
   character(len=256), intent(in) :: dir_landdata  !surface data directory
   character(len=256), intent(in) :: dir_restart   !case restart data directory

   logical, intent(in)    :: greenwich   !true: greenwich time, false: local time
   integer, intent(inout) :: jdate(3)    !year, julian day, seconds of the starting time
#if (defined TRACER) && (defined BGC)
   integer :: op_methane
   real(r8), allocatable :: methane_old_patch_area(:)
#endif
!-----------------------------------------------------------------------

      ! allocate Lulcc memory
      CALL allocate_LulccTimeInvariants
      CALL allocate_LulccTimeVariables

      ! SAVE variables
      CALL SAVE_LulccTimeInvariants
      CALL SAVE_LulccTimeVariables
#if (defined TRACER) && (defined BGC)
      IF (igas_ch4 > 0) THEN
         CALL save_methane_lulcc_state ()
         CALL save_methane_microbes_lulcc_state ()
      ENDIF
#endif

      ! =============================================================
      ! cold start for Lulcc
      ! =============================================================

      IF (p_is_master) THEN
         print *, ">>> LULCC: initializing..."
      ENDIF

      CALL LulccInitialize (casename, dir_landdata, dir_restart, &
                            jdate, greenwich)


      ! =============================================================
      ! 1. Same Type Assignment (SAT) scheme for variable recovery
      ! =============================================================

      IF (DEF_LULCC_SCHEME == 1) THEN
         IF (p_is_master) THEN
            print *, ">>> LULCC: Same Type Assignment (SAT) scheme for variable recovery..."
         ENDIF
         CALL REST_LulccTimeVariables
#if (defined TRACER) && (defined BGC)
         IF (igas_ch4 > 0) THEN
            CALL deallocate_methane_acc_fluxes ()
            CALL remap_methane_lulcc_state (patchclass, landpatch%eindex, patchclass_, landpatch_%eindex)
            CALL remap_methane_microbes_lulcc_state (patchclass, landpatch%eindex, patchclass_, landpatch_%eindex)
            CALL initialize_methane_lake_soilc_from_surface (patchtype, lake_soilc_srf, DEF_METHANE%allowlakeprod)
            CALL allocate_methane_acc_fluxes (size(patchclass))
         ENDIF
#endif
      ENDIF


      ! =============================================================
      ! 2. Mass and Energy conservation (MEC) scheme for variable recovery
      ! =============================================================

      IF (DEF_LULCC_SCHEME == 2) THEN
         IF (p_is_master) THEN
            print *, ">>> LULCC: Mass&Energy conserve (MEC) for variable recovery..."
         ENDIF
         CALL allocate_LulccTransferTrace()
         CALL REST_LulccTimeVariables
         CALL LulccTransferTraceReadin(jdate(1))
         CALL LulccMassEnergyConserve()
#if (defined TRACER) && (defined BGC)
         IF (igas_ch4 > 0) THEN
            CALL deallocate_methane_acc_fluxes ()
            allocate(methane_old_patch_area(size(patchclass_)))
	            methane_old_patch_area(:) = 1._r8
	            DO op_methane = 1, size(methane_old_patch_area)
	               IF (allocated(landpatch_%ipxstt) .and. allocated(landpatch_%ipxend) .and. &
	                   op_methane <= size(landpatch_%ipxstt) .and. op_methane <= size(landpatch_%ipxend)) THEN
	                  IF (landpatch_%ipxstt(op_methane) > 0 .and. &
	                      landpatch_%ipxend(op_methane) >= landpatch_%ipxstt(op_methane)) THEN
	                     methane_old_patch_area(op_methane) = real(landpatch_%ipxend(op_methane) - &
	                        landpatch_%ipxstt(op_methane) + 1, r8)
	                  ENDIF
	               ENDIF
	               IF (allocated(landpatch_%pctshared) .and. op_methane <= size(landpatch_%pctshared)) THEN
	                  methane_old_patch_area(op_methane) = methane_old_patch_area(op_methane) * &
	                     max(0._r8, landpatch_%pctshared(op_methane))
	               ENDIF
            ENDDO
            CALL remap_methane_lulcc_state (patchclass, landpatch%eindex, patchclass_, &
               landpatch_%eindex, lccpct_patches, methane_old_patch_area)
            CALL remap_methane_microbes_lulcc_state (patchclass, landpatch%eindex, patchclass_, &
               landpatch_%eindex, lccpct_patches, methane_old_patch_area)
            IF (allocated(methane_old_patch_area)) deallocate(methane_old_patch_area)
            CALL initialize_methane_lake_soilc_from_surface (patchtype, lake_soilc_srf, DEF_METHANE%allowlakeprod)
            CALL allocate_methane_acc_fluxes (size(patchclass))
         ENDIF
#endif
      ENDIF


      ! deallocate Lulcc memory
      CALL deallocate_LulccTimeInvariants()
      CALL deallocate_LulccTimeVariables()
      IF (DEF_LULCC_SCHEME == 2) THEN
         CALL deallocate_LulccTransferTrace()
      ENDIF

   END SUBROUTINE LulccDriver

END MODULE MOD_Lulcc_Driver
#endif
! ---------- EOP ------------
