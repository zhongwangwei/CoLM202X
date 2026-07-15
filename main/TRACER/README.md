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

## Methane model contract

The implemented CH4 model is an offline land diagnostic. It predicts soil,
wetland, rice-paddy, and reduced-order lake-sediment CH4 fluxes, but it does
not feed those fluxes back into an atmospheric CH4 state. Wetland production
uses the BGC decomposition cascade and debits its donor C/N pools at source;
the emitted CH4 flux does not apply a second carbon-budget debit.
`DEF_METHANE%methane_offline` must therefore remain true; requesting online
coupling is rejected instead of silently running the offline equations.

Lake CH4 is a reduced-order sediment-column-to-atmosphere calculation. The
soil-like lake sediment column represents production, oxidation, diffusion,
and ebullition below the interface. There is no resolved lake-water CH4/O2
column, storage, overturning, ice-cover transport, or air-water piston model.
`allowlakeprod` enables this reduced-order pathway only; it must not be read as
enabling a complete aquatic methane model.

The physical surface-flux diagnostics exclude numerical concentration clipping
and column-balance corrections. Corrected diagnostics include those terms so
that restart-to-restart inventory budgets close. Both forms, together with the
explicit correction terms, must remain observable; numerical closure must not
be mistaken for a physical emission process.

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

## Particle tracers

Particle tracers use the same flat TRACER-owned pattern. `MOD_Tracer_Particle`
is the generic HYDRO-facing dispatcher; it does not name sediment or any other
concrete particle species. Implemented species register callbacks through
`register_particle_callbacks`, and `include/tracer_particle_species.inc` is the
single implemented-particle manifest.

Suspended sediment is currently the only particle species. Its implementation,
restart, history, and accumulator flushing live in
`MOD_Tracer_Particle_Sediment`. HYDRO modules should call only the generic
`tracer_particle_*` APIs so species-private sediment state does not leak back
into `main/HYDRO`.

Sediment is enabled exactly like any other tracer species: by listing it in
the generic tracer registry as a particle. There is no top-level sediment
``USE`` switch.

```fortran
DEF_TRACER_NAMES = '...,SEDIMENT'
DEF_TRACER_TYPES = '...,particle'
DEF_TRACER_PARAM_FILES = 'SEDIMENT:./standard_sediment_parameter.nml'
```

The optional sediment parameter file is read by `MOD_Tracer_Particle_Sediment`
from `&nl_colm_sediment_parameter`; generic tracer metadata in the same file is
read by `MOD_Tracer_Defs`.

For a future particle tracer, add its species-owned module, one manifest line,
and one Makefile object entry under `TRACER_PARTICLE_REGISTERED_SPECIES_OBJS`.
Do not add placeholder particle species.
