! Reactive tracer species manifest.
!
! Add one line per implemented reactive species:
!   TRACER_REACTIVE_SPECIES(Module_Name, register_callbacks_subroutine)
!
! MOD_Tracer_Reactive_Registrations.F90 expands this single manifest into
! both USE statements and registration CALLs. Do not add placeholder species;
! add a line only when the species implementation exists.
TRACER_REACTIVE_SPECIES(MOD_Tracer_Reactive_Methane, ch4_register_reactive_callbacks)
