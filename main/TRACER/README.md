# TRACER reactive architecture

Reactive tracers are kept flat at the species layer and generic above it:

- `MOD_Tracer_Reactive` is the generic dispatcher. It must not name concrete
  species such as CH4, N2O, O3, or isotope tracers.
- Implemented reactive species register callbacks through
  `register_reactive_callbacks`.
- Every species registration must provide `has_fn`; the generic dispatcher uses
  it to gate lifecycle callbacks, including `init`.
- Species with cached registry/index state should also provide `refresh_fn`;
  the generic dispatcher calls refresh before `has_fn` only when callback state
  is marked dirty (initial registration and after LULCC remap), not on every
  per-patch lifecycle dispatch.
- The callback registry grows dynamically; do not introduce a fixed maximum
  number of reactive species.
- `include/tracer_reactive_species.inc` is the single implemented-species
  manifest. Add a line there only when the species implementation exists.
- Do not add placeholder species modules or objects.

For a future reactive tracer, add its species-owned modules and one manifest
line, and add the species facade object to the Makefile species object
manifest `TRACER_REACTIVE_REGISTERED_SPECIES_OBJS`:

```fortran
TRACER_REACTIVE_SPECIES(MOD_Tracer_Reactive_<Species>, <species>_register_reactive_callbacks)
```

The future species module should own all species-specific wiring and expose one
registration routine. It should call `register_reactive_callbacks` with its
init/restart/step/history/LULCC/HYDRO callbacks. The generic dispatcher should
not change for a new species.

Species-private state modules should use default `PRIVATE`. If raw arrays must
remain externally visible for legacy internal code, they must be explicitly
`PUBLIC` and named with a species prefix such as `methane_*`. New cross-layer
data exchange should prefer small species-owned API routines over direct writes
to raw state arrays.

Current CH4 legacy raw-array boundaries are intentionally narrow:

- Methane_State raw state must not be imported outside the CH4 internal
  boundary. External coupling such as LULCC/HYDRO must go through
  `MOD_Tracer_Reactive` lifecycle/publish hooks and species-owned API routines.
- Methane_AccFlux raw accumulators must not be imported outside the CH4
  facade/history boundary. History output is the only current external reader
  of those accumulator arrays.

## Isotope fractionation physics

Conservative/isotope tracer transport stays generic over `itrc = 1, ntracers`.
Species-specific isotope fractionation physics is registered separately:

- `MOD_Tracer_Frac` owns the generic fractionation algorithms and must not
  classify O18/HDO locally.
- `MOD_Tracer_Isotope_Registry` stores implemented isotope physics callbacks
  and metadata such as name patterns, reference-ratio hints, legacy forcing
  shortcut kind, and default soil-init variable name.
- `include/tracer_isotope_species.inc` is the single implemented-isotope
  manifest. It currently registers only the existing O18 and HDO physics.
- `MOD_Tracer_Forcing` and `MOD_Tracer_SoilInit` get default isotope
  auto-detection from the isotope registry instead of keeping duplicate
  O18/HDO name-matching chains.
- Do not add placeholder isotope species. Add a new isotope only when its
  physics module exists, then add one manifest line and one Makefile object
  entry under `TRACER_ISOTOPE_REGISTERED_SPECIES_OBJS`.
