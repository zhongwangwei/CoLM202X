#include <define.h>

MODULE MOD_Tracer_Vars_3DAccFluxes

   USE MOD_Precision
   USE MOD_Tracer_Vars_3DForcing, only: num_tracers_3d, max_vars_3d
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE
   SAVE

   ! Number of accumulation counters for 3D tracer data
   real(r8) :: tracer_3d_nac ! number of accumulation

   ! 3D accumulated tracer forcing array matching tracer_forc structure
   real(r8), allocatable :: a_tracer_forc_3d(:,:,:)  ! [tracer_idx, var_idx, patch_idx]

   ! Module control variables
   logical :: tracer_3d_acc_initialized = .false.

   PUBLIC :: allocate_tracer_3d_acc_fluxes
   PUBLIC :: deallocate_tracer_3d_acc_fluxes
   PUBLIC :: FLUSH_tracer_3d_acc_fluxes
   PUBLIC :: accumulate_tracer_3d_fluxes

CONTAINS

   SUBROUTINE allocate_tracer_3d_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   integer :: ierr

      IF (p_is_master) THEN
         WRITE(*,*) "Allocating 3D tracer accumulation arrays:"
         WRITE(*,*) "  num_tracers_3d:", num_tracers_3d
         WRITE(*,*) "  max_vars_3d:", max_vars_3d
         WRITE(*,*) "  numpatch:", numpatch
      ENDIF

      ! Allocate accumulation arrays only on worker processes
      IF (p_is_worker .AND. numpatch > 0) THEN
         
         IF (num_tracers_3d > 0 .AND. max_vars_3d > 0) THEN
            ! Allocate main 3D accumulation array with same structure as tracer_forc
            allocate (a_tracer_forc_3d(num_tracers_3d, max_vars_3d, numpatch), stat=ierr)
            IF (ierr /= 0) THEN
               WRITE(*,*) "ERROR: Failed to allocate a_tracer_forc_3d array"
               RETURN
            ENDIF
         ENDIF
      ENDIF

      ! Set initialization flag to true after successful allocation
      tracer_3d_acc_initialized = .true.
      
      IF (p_is_master) THEN
         WRITE(*,*) "3D tracer accumulation system initialized successfully"
      ENDIF

   END SUBROUTINE allocate_tracer_3d_acc_fluxes


   SUBROUTINE deallocate_tracer_3d_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE

      IF (p_is_worker .AND. numpatch > 0) THEN
         IF (allocated(a_tracer_forc_3d)) deallocate(a_tracer_forc_3d)
      ENDIF

      tracer_3d_acc_initialized = .false.

      IF (p_is_master) THEN
         WRITE(*,*) "3D tracer accumulation arrays deallocated"
      ENDIF

   END SUBROUTINE deallocate_tracer_3d_acc_fluxes


   SUBROUTINE FLUSH_tracer_3d_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE
   
   IF (p_is_worker) THEN
      tracer_3d_nac = 0.
      IF (numpatch > 0 ) THEN
         ! Initialize accumulation arrays with zero for proper accumulation
         a_tracer_forc_3d(:,:,:) = spval
      ENDIF
   ENDIF
   
   END SUBROUTINE FLUSH_tracer_3d_acc_fluxes


   SUBROUTINE accumulate_tracer_3d_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_Tracer_Vars_3DForcing, only: tracer_forc
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   integer :: itrace, ivar, ipatch

      IF (.NOT. tracer_3d_acc_initialized) THEN
         IF (p_is_master) WRITE(*,*) "WARNING: tracer_3d_acc_initialized is false, skipping 3D accumulation"
         RETURN
      ENDIF
      
      IF (p_is_worker .AND. numpatch > 0) THEN
        tracer_3d_nac = tracer_3d_nac + 1.
        call tracer_acc3d (tracer_forc, a_tracer_forc_3d)
      ENDIF
      
   END SUBROUTINE accumulate_tracer_3d_fluxes



   SUBROUTINE tracer_acc3d (var, s)

    USE MOD_Precision
    USE MOD_Vars_Global, only: spval
 
    IMPLICIT NONE
 
    real(r8), intent(in)    :: var(:,:,:)
    real(r8), intent(inout) :: s  (:,:,:)
    ! Local variables
    integer :: i1, i2, i3
 
       DO i3 = lbound(var,3), ubound(var,3)
          DO i2 = lbound(var,2), ubound(var,2)
             DO i1 = lbound(var,1), ubound(var,1)
                IF (var(i1,i2,i3) /= spval) THEN
                   IF (s(i1,i2,i3) /= spval) THEN
                      s(i1,i2,i3) = s(i1,i2,i3) + var(i1,i2,i3)
                   ELSE
                      s(i1,i2,i3) = var(i1,i2,i3)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
 
    END SUBROUTINE tracer_acc3d

END MODULE MOD_Tracer_Vars_3DAccFluxes 