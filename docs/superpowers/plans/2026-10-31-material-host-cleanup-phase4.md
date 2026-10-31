# Comprehensive material-host extraction: phase 4

User directs one comprehensive migration pass before tests; do not repeat long
integration tests per module. Prior verified working tree is frozen at
/tmp/colm-host-phase4-8xvcb9cc with SHA256 source-manifest.json; the completed
phase-3 model job codex-00988b4f47864faba508af0815e1aa3e is the numerical baseline.
Preserve all uncommitted work from earlier rounds. No production edits precede
this freeze. Existing coverage and frozen sources lock behavior for comparison.

## Batched ownership

1. Thermal interfaces and bookkeeping: Thermal, LeafTemperature/PC,
   GroundTemperature, PhaseChange, GroundFluxes, MAIN and lightweight land work.
2. Runtime lifecycle/history/variable bookkeeping: CoLM, Hist, Vars time/flux
   modules; material-side history/runtime helper with no dependency cycles.
3. Water/snow/interception host intrusions: SoilSnowHydrology,
   SnowLayersCombineDivide, LeafInterception, Hydro_SoilWater; material helpers.
4. Routing/driver/LULCC intrusions: CoLMDRIVER, HYDRO grid routing modules,
   LULCC driver and material routing helpers.
5. BGC, preprocessing and initialization: BGC, mksrfdata, mkinidata,
   share/Namelist; move material-only implementations, retain host-only features.

Parent also moves material config/build dependency rules into MATERIAL/Makefile.material.
Configuration entrypoint guards stay in define.h: they validate the selected
compiler options and are necessary at that boundary.

Parent owns Makefile/build integration, final whole-repo residual inventory and
validation. Each owner writes a bounded proposal before edits, migrates all safe
blocks in their slice, records callers/config dependencies and what must stay.
Do not move host water/energy/BGC algorithms or duplicate complete host solvers.
No result, event ordering, arithmetic order, guard, optional-presence, lower-bound,
error handling, configuration, history/restart schema or MPI ownership changes.
Prefer existing material modules; new lightweight modules only break real cycles
or own cohesive material-only storage. No new runtime dependencies.

## One verification phase after all migration

- Independent read-only review and residual audit (including untouched files).
- Compile and run focused baseline comparisons adapted to final interfaces;
  inspect all tests rather than weakening assertions to accommodate changes.
- Compare effective MT-off host code to frozen sources wherever migration should
  leave it identical; verify CPP variants affected by rewritten interfaces.
- Compile MT-on/off, relevant PFT/PC/CROP/snow/routing modes; diagnose dependencies.
- One final T7920 continuous96 and split48+48 fixture against phase-3 outputs;
  strict finite/file/variable/time coverage and all eight output comparisons.
- Long builds and model integration only in private run-7920 snapshots. Keep all
  existing shared server sources, old results and external inputs read-only.
