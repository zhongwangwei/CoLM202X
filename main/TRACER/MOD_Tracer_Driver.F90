#include <define.h>

MODULE MOD_Tracer_Driver

USE MOD_Precision, only: r8
USE MOD_Tracer_Namelist_Defs, only: DEF_Tracer_Number, DEF_USE_Tracer
USE MOD_Tracer_State, only: allocate_tracer_state, deallocate_tracer_state
USE MOD_Tracer_Forcing, only: tracer_forcing_init, read_tracer_forcing 
! Add tracer_forcing_final if it gets created in MOD_Tracer_Forcing
USE MOD_Tracer_Initialize, only: tracer_initialize_concentrations

USE MOD_SPMD_Task, only: numpatch, p_is_master ! For numpatch in allocate_tracer_state
USE MOD_Vars_Soil, only: nl_soil      ! For nl_soil in allocate_tracer_state
USE MOD_TimeManager, only: timestamp   ! For current_model_time argument type
USE MOD_LandPatch, only: landpatch_vector_type ! For landpatch_data argument type in tracer_forcing_init
USE MOD_SPMD_Task, only: CoLM_stop ! For error handling

IMPLICIT NONE
SAVE

PUBLIC :: Tracer_Initialize_Master
PUBLIC :: Tracer_Advance_Timestep
PUBLIC :: Tracer_Finalize

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE Tracer_Initialize_Master(model_timestep_sec, landpatch_data_for_forcing)
!-----------------------------------------------------------------------
! Initializes all tracer-related components.
! Assumes main model namelists (including tracer flags) have been read by MOD_Namelist.
!-----------------------------------------------------------------------
  real(r8), intent(in) :: model_timestep_sec
  type(landpatch_vector_type), intent(in) :: landpatch_data_for_forcing ! Or actual type from MOD_LandPatch

  IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: Starting tracer initialization..."

  IF (DEF_USE_Tracer .AND. DEF_Tracer_Number > 0) THEN
     ! 1. Allocate tracer state variables
     !    numpatch and nl_soil are assumed to be available from their respective modules.
     IF (numpatch > 0 .AND. nl_soil > 0) THEN
        CALL allocate_tracer_state(DEF_Tracer_Number, numpatch, nl_soil)
        IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: Tracer state allocated for ", DEF_Tracer_Number, " tracers."
     ELSE
        IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: numpatch or nl_soil is zero, cannot allocate tracer state."
        CALL CoLM_stop("Tracer_Initialize_Master: numpatch or nl_soil is zero.")
        RETURN
     ENDIF

     ! 2. Initialize tracer forcing system
     CALL tracer_forcing_init(model_timestep_sec, landpatch_data_for_forcing)
     IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: Tracer forcing system initialized."

     ! 3. Initialize tracer concentrations (from file or defaults)
     CALL tracer_initialize_concentrations()
     IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: Tracer concentrations initialized."

     IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: Tracer system initialization complete."
  ELSE
     IF (p_is_master) WRITE(*,*) "Tracer_Initialize_Master: Tracers are not used (DEF_USE_Tracer is false or DEF_Tracer_Number is 0)."
  ENDIF

END SUBROUTINE Tracer_Initialize_Master

!-----------------------------------------------------------------------
SUBROUTINE Tracer_Advance_Timestep(current_model_time, dt_seconds)
!-----------------------------------------------------------------------
! Called each model timestep to advance tracer processes.
!-----------------------------------------------------------------------
  type(timestamp), intent(in) :: current_model_time
  real(r8), intent(in) :: dt_seconds ! Current model timestep length

  IF (DEF_USE_Tracer .AND. DEF_Tracer_Number > 0) THEN
     ! 1. Read and process tracer forcing data for the current timestep
     CALL read_tracer_forcing(current_model_time)

     ! 2. Tracer transport logic (e.g., in soil) is embedded within the
     !    respective physics modules (e.g., MOD_Hydro_SoilWater) and is
     !    called as part of their timestep advancement.
     !    No explicit call needed here for that part based on current design.

     ! Other per-timestep tracer processes could be called here if they exist
     ! (e.g., chemical reactions, decay - not in current scope).
  ENDIF

END SUBROUTINE Tracer_Advance_Timestep

!-----------------------------------------------------------------------
SUBROUTINE Tracer_Finalize()
!-----------------------------------------------------------------------
! Finalizes tracer components at the end of the simulation.
!-----------------------------------------------------------------------
  IF (p_is_master) WRITE(*,*) "Tracer_Finalize: Finalizing tracer components..."

  IF (DEF_USE_Tracer .AND. DEF_Tracer_Number > 0) THEN
     ! 1. Deallocate tracer state variables
     CALL deallocate_tracer_state()
     IF (p_is_master) WRITE(*,*) "Tracer_Finalize: Tracer state deallocated."

     ! 2. Call finalize for tracer forcing (if it exists)
     ! CALL tracer_forcing_finalize() ! Add if MOD_Tracer_Forcing gets a finalize routine

     IF (p_is_master) WRITE(*,*) "Tracer_Finalize: Tracer system finalization complete."
  ELSE
      IF (p_is_master) WRITE(*,*) "Tracer_Finalize: Tracers were not used."
  ENDIF

END SUBROUTINE Tracer_Finalize

END MODULE MOD_Tracer_Driver
