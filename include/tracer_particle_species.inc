! Particle tracer species manifest.
!
! Add one line per implemented particle species:
!   TRACER_PARTICLE_SPECIES(Module_Name, register_callbacks_subroutine)
!
! MOD_Tracer_Particle_Registrations.F90 expands this single manifest into
! both USE statements and registration CALLs. Do not add placeholder particle
! species; add a line only when the implementation exists.
TRACER_PARTICLE_SPECIES(MOD_Tracer_Particle_Sediment, register_sediment_particle_callbacks)
