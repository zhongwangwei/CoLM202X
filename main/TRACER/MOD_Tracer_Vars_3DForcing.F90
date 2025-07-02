#include <define.h>

MODULE MOD_Tracer_Vars_3DForcing
!-----------------------------------------------------------------------
!  3D Tracer Forcing Variables
!
!  Created by Yongjiu Dai, 03/2014
!  Modified for 3D tracer integration, 12/2024
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL, MAX_TRACER_FORCING_VARS
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
   ! 3D tracer forcing array: (num_tracers, max_vars, numpatch)
   real(r8), allocatable, target :: tracer_forc(:,:,:) ! tracer forcing data
   
   ! Global parameters for tracer dimensions
   integer, public :: num_tracers_3d = 0      ! number of active tracers
   integer, public :: max_vars_3d = 0         ! maximum variables per tracer
   
   ! Tracer names and variable names for reference
   character(len=64), allocatable, public :: tracer_names_3d(:)
   character(len=64), allocatable, public :: tracer_var_names_3d(:,:)

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_3D_Forcing_Tracer
   PUBLIC :: Tracer_deallocate_3D_Forcing_Tracer

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
   SUBROUTINE Tracer_allocate_3D_Forcing_Tracer
      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE MOD_LandPatch
      IMPLICIT NONE
      
      integer :: i, j, itrace, ivar
      
      IF (p_is_master) THEN
         WRITE(*,*) "=== Allocating 3D Tracer Forcing Arrays ==="
         
         ! Count active tracers and maximum variables on master process only
         num_tracers_3d = 0
         max_vars_3d = 0
         
         DO i = 1, SIZE(DEF_Tracer_Forcings_NL, 1)
            IF (trim(DEF_Tracer_Forcings_NL(i)%tracer_name) /= '') THEN
               num_tracers_3d = num_tracers_3d + 1
               IF (DEF_Tracer_Forcings_NL(i)%NVAR > max_vars_3d) THEN
                  max_vars_3d = DEF_Tracer_Forcings_NL(i)%NVAR
               ENDIF
            ENDIF
         ENDDO
         
         WRITE(*,*) "Number of tracers:", num_tracers_3d
         WRITE(*,*) "Maximum variables per tracer:", max_vars_3d
      ENDIF
      
#ifdef USEMPI
      ! Broadcast the computed values to all processes to ensure consistency
      CALL mpi_bcast(num_tracers_3d, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(max_vars_3d, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      
      IF (p_is_worker) THEN
         IF (numpatch > 0 .and. num_tracers_3d > 0 .and. max_vars_3d > 0) THEN
            ! Allocate main 3D array: (num_tracers, max_vars, numpatch)
            allocate(tracer_forc(num_tracers_3d, max_vars_3d, numpatch))
            tracer_forc = spval
         ENDIF
      ENDIF
      
      ! Print allocation info from master only
      IF (p_is_master .and. num_tracers_3d > 0 .and. max_vars_3d > 0) THEN
         WRITE(*,*) "3D tracer array allocated: (", num_tracers_3d, ",", max_vars_3d, ",", numpatch, ")"
      ENDIF
      
      ! Allocate and fill tracer names arrays on all processes
      IF (num_tracers_3d > 0) THEN
         allocate(tracer_names_3d(num_tracers_3d))
         allocate(tracer_var_names_3d(num_tracers_3d, max_vars_3d))
         
         tracer_names_3d = ''
         tracer_var_names_3d = ''
         
         IF (p_is_master) THEN
            ! Fill tracer names and variable names on master process only
            itrace = 0
            DO i = 1, SIZE(DEF_Tracer_Forcings_NL, 1)
               IF (trim(DEF_Tracer_Forcings_NL(i)%tracer_name) /= '') THEN
                  itrace = itrace + 1
                  tracer_names_3d(itrace) = DEF_Tracer_Forcings_NL(i)%tracer_name
                  
                  DO ivar = 1, DEF_Tracer_Forcings_NL(i)%NVAR
                     tracer_var_names_3d(itrace, ivar) = DEF_Tracer_Forcings_NL(i)%vname(ivar)
                  ENDDO
                  
                  WRITE(*,*) "Tracer", itrace, ":", trim(tracer_names_3d(itrace)), &
                            " with", DEF_Tracer_Forcings_NL(i)%NVAR, "variables"
               ENDIF
            ENDDO
         ENDIF
          
#ifdef USEMPI
         ! Broadcast tracer names to all processes
         DO itrace = 1, num_tracers_3d
            CALL mpi_bcast(tracer_names_3d(itrace), 64, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
            DO ivar = 1, max_vars_3d
               CALL mpi_bcast(tracer_var_names_3d(itrace, ivar), 64, MPI_CHARACTER, p_address_master, p_comm_glb, p_err)
            ENDDO
         ENDDO
#endif
      ENDIF
      
      IF (p_is_master) THEN
         WRITE(*,*) "3D Tracer forcing allocation completed!"
      ENDIF
      
   END SUBROUTINE Tracer_allocate_3D_Forcing_Tracer

   !-----------------------------------------------------------------------
   SUBROUTINE Tracer_deallocate_3D_Forcing_Tracer
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      IF (p_is_master) THEN
         WRITE(*,*) "Deallocating 3D Tracer Forcing Arrays..."
      ENDIF
      
      IF (allocated(tracer_forc)) deallocate(tracer_forc)
      IF (allocated(tracer_names_3d)) deallocate(tracer_names_3d)
      IF (allocated(tracer_var_names_3d)) deallocate(tracer_var_names_3d)
      
      num_tracers_3d = 0
      max_vars_3d = 0
      
      IF (p_is_master) THEN
         WRITE(*,*) "3D Tracer forcing deallocation completed!"
      ENDIF
      
   END SUBROUTINE Tracer_deallocate_3D_Forcing_Tracer
END MODULE MOD_Tracer_Vars_3DForcing
! ---------- EOP ------------
