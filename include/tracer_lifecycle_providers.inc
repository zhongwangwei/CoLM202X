! Implemented TRACER lifecycle providers.
!
! Parameter-only isotope/solute tracers do not belong here. Add one entry
! only when a compiled species owns state or phase callbacks.
#ifdef BGC
TRACER_LIFECYCLE_PROVIDER(MOD_Tracer_Reactive_Methane, ch4_register_tracer_provider)
#endif
TRACER_LIFECYCLE_PROVIDER(MOD_Tracer_Particle_Sediment, register_sediment_tracer_provider)
