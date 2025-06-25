#include <define.h>

MODULE MOD_Tracer_Vars_1DAccFluxes

   USE MOD_Precision
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, DEF_Tracer_Forcings_NL, DEF_Tracer_Number, MAX_TRACER_FORCING_VARS

   IMPLICIT NONE

   ! Number of accumulation counters
   real(r8) :: tracer_nac ! number of accumulation

   ! Accumulated tracer forcing variables for each tracer and variable
   ! Using dynamic allocation based on actual configured tracers
   real(r8), allocatable :: a_tracer_forc(:,:,:)  ! [tracer_idx, var_idx, patch_idx]
   
   ! More specific tracer accumulation variables for commonly used tracers
   ! These provide easier access for specific tracers like O18, H2, etc.
   real(r8), allocatable :: a_tracer_prate(:,:)   ! precipitation tracer [tracer_idx, patch_idx]
   real(r8), allocatable :: a_tracer_spfh(:,:)    ! humidity tracer [tracer_idx, patch_idx]

   ! Module control variables
   integer :: num_tracers_acc     ! number of tracers with accumulation
   integer :: max_vars_acc        ! maximum variables per tracer
   logical :: tracer_acc_initialized = .false.

   PUBLIC :: allocate_tracer_acc_fluxes
   PUBLIC :: deallocate_tracer_acc_fluxes
   PUBLIC :: FLUSH_tracer_acc_fluxes
   PUBLIC :: accumulate_tracer_fluxes

CONTAINS

   SUBROUTINE allocate_tracer_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE

   integer :: i, j, ierr

      ! Count the number of tracers with forcing
      num_tracers_acc = 0
      max_vars_acc = 0
      
      DO i = 1, SIZE(DEF_Tracer_Forcings_NL, 1)
         IF (trim(DEF_Tracer_Forcings_NL(i)%tracer_name) /= '') THEN
            num_tracers_acc = num_tracers_acc + 1
            IF (DEF_Tracer_Forcings_NL(i)%NVAR > max_vars_acc) THEN
               max_vars_acc = DEF_Tracer_Forcings_NL(i)%NVAR
            ENDIF
         ENDIF
      ENDDO

      IF (p_is_master) THEN
         WRITE(*,*) "Allocating tracer accumulation arrays:"
         WRITE(*,*) "  Number of tracers: ", num_tracers_acc
         WRITE(*,*) "  Maximum variables per tracer: ", max_vars_acc
         WRITE(*,*) "  Number of patches: ", numpatch
      ENDIF

      ! Allocate accumulation arrays only on worker processes
      IF (p_is_worker .AND. numpatch > 0) THEN
         
         ! Allocate main accumulation array
         allocate (a_tracer_forc(num_tracers_acc, max_vars_acc, numpatch), stat=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) "ERROR: Failed to allocate a_tracer_forc array"
            CALL CoLM_stop()
         ENDIF

         ! Allocate specific tracer arrays for common access patterns
         allocate (a_tracer_prate(num_tracers_acc, numpatch), stat=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) "ERROR: Failed to allocate a_tracer_prate array"
            CALL CoLM_stop()
         ENDIF

         allocate (a_tracer_spfh(num_tracers_acc, numpatch), stat=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) "ERROR: Failed to allocate a_tracer_spfh array"
            CALL CoLM_stop()
         ENDIF

      ENDIF

      tracer_acc_initialized = .true.
      
      IF (p_is_master) THEN
         WRITE(*,*) "Tracer accumulation arrays allocated successfully"
      ENDIF

   END SUBROUTINE allocate_tracer_acc_fluxes


   SUBROUTINE deallocate_tracer_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE

      IF (p_is_worker .AND. numpatch > 0) THEN
         IF (allocated(a_tracer_forc)) deallocate(a_tracer_forc)
         IF (allocated(a_tracer_prate)) deallocate(a_tracer_prate)
         IF (allocated(a_tracer_spfh)) deallocate(a_tracer_spfh)
      ENDIF

      tracer_acc_initialized = .false.

      IF (p_is_master) THEN
         WRITE(*,*) "Tracer accumulation arrays deallocated"
      ENDIF

   END SUBROUTINE deallocate_tracer_acc_fluxes


   SUBROUTINE FLUSH_tracer_acc_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

      tracer_nac = 0.

      IF (p_is_worker .AND. numpatch > 0 .AND. tracer_acc_initialized) THEN
         
         IF (allocated(a_tracer_forc)) THEN
            a_tracer_forc(:,:,:) = 0.
         ENDIF

         IF (allocated(a_tracer_prate)) THEN
            a_tracer_prate(:,:) = 0.
         ENDIF

         IF (allocated(a_tracer_spfh)) THEN
            a_tracer_spfh(:,:) = 0.
         ENDIF

      ENDIF

   END SUBROUTINE FLUSH_tracer_acc_fluxes


   SUBROUTINE accumulate_tracer_fluxes ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_Tracer_Forcing, only: get_tracer_forcing_data
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   integer :: i, j, k
   real(r8), pointer :: tracer_data_ptr(:)

      IF (.NOT. tracer_acc_initialized) RETURN
      IF (.NOT. p_is_worker .OR. numpatch == 0) RETURN

      tracer_nac = tracer_nac + 1.

      ! Loop through all tracers and accumulate their forcing data
      DO i = 1, num_tracers_acc
         DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
            IF (j <= max_vars_acc .AND. trim(DEF_Tracer_Forcings_NL(i)%vname(j)) /= 'NULL') THEN
               
               tracer_data_ptr => get_tracer_forcing_data(i, j)
               
               IF (associated(tracer_data_ptr)) THEN
                  IF (SIZE(tracer_data_ptr) == numpatch) THEN
                     DO k = 1, numpatch
                        IF (tracer_data_ptr(k) /= spval) THEN
                           a_tracer_forc(i, j, k) = a_tracer_forc(i, j, k) + tracer_data_ptr(k)
                           
                           ! Also accumulate to specific arrays for common variables
                           SELECT CASE (trim(DEF_Tracer_Forcings_NL(i)%vname(j)))
                           CASE ('prate1sfc', 'prate2sfc')
                              a_tracer_prate(i, k) = a_tracer_prate(i, k) + tracer_data_ptr(k)
                           CASE ('spfh12m', 'spfh22m')
                              a_tracer_spfh(i, k) = a_tracer_spfh(i, k) + tracer_data_ptr(k)
                           END SELECT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE accumulate_tracer_fluxes


END MODULE MOD_Tracer_Vars_1DAccFluxes
! ---------- EOP ------------
