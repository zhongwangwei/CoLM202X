# Material host cleanup phase 2 validation

## Scope and baseline

This pass preserves the uncommitted phase-1 work. Its baseline is the complete
phase-1 working tree, whose coupled/restart outputs were already validated,
not Git HEAD. The five pre-edit files and SHA256 manifest are retained under
`.omx/artifacts/material-host-cleanup-phase2-20261031/baseline/`.

- CoLMMAIN: 2,052 -> 1,960 lines (92 fewer in this pass; 240 fewer across both).
- Post-interception release/gross-interception calculations and storage/precip
  dispatch moved to `MOD_MT_LandPhase::mt_land_interception`.
- Work storage moved to lightweight `MOD_MT_LandWork`, with bound methods for
  snow capture, pre-thermal preparation, post-thermal water capture and release.
- Recharge/infiltration conversion moved to the land facade; flood evaporation
  ledger access moved to `MOD_MT_RiverLakeCoupling::mt_rlc_land_flood_evap`.
- The exact decay -> balance -> history -> release sequence moved to
  `MOD_MT_Hist::mt_land_finish_step`. The separate Work module avoids the real
  Hist -> 1DAccFluxes -> LandPhase dependency cycle.
- Unsupported CaMa generic-water feedback is rejected at the original event
  through a land-facade helper, with the original error message preserved.

Water/thermal algorithms, phase sampling points, state fields, configuration,
restart schema and snow-layer topology are unchanged.

## Local and independent checks

| Check | Result |
| --- | --- |
| Existing focused material physics | PASS before/after: 9 + 15 + 24 numerical checks; 5 wiring checks |
| Phase bookkeeping/lifecycle | PASS: 384 CROP/flood/snow/patch cases; baseline callback values/order bitwise identical |
| CaMa capability guard | PASS: 5 inactive/provider/unsupported paths, exact error retained |
| Land boundary wiring | PASS: 64 cases; real Work type, all arguments and optional presence identical |
| Flood evaporation | PASS: 7 cases and 7 baseline bitwise comparisons; guards, units, clear-before-call, patch isolation |
| Workspace compilation | PASS with GridRiverLakeFlow on/off and runtime bounds checks; MT-off empty module compiles |
| MT-off CoLMMAIN preprocessing | PASS: effective source identical to phase-1 baseline |
| Naming inventory / Python syntax / whitespace | PASS |
| Independent production review | APPROVE; no findings; new dependency graph has no cycle |

The callback tests use real storage and production helper bodies with recording
service stubs. They verify bookkeeping/forwarding/order, not the full physical
model. The coupled validation below is separate evidence.

## Remote builds and integration

The 12.3 MB code-only source snapshot includes the current local edits and a
per-file SHA256 manifest (also retained in the local artifact directory).
All data inputs were preflighted on T7920 and bound read-only; no datasets were
uploaded. Source-content checks identify changed/new files; Make recompiles
them and their dependencies in private copies of the previous build. Both main
program objects are rebuilt with the selected flags.

| Run | Job | Result |
| --- | --- | --- |
| Complete MT-off compile/link | `codex-29b6d1aefad248499796634ea28e3e33` | SUCCEEDED, exit 0 |
| MT-on compile + continuous/split integration | `codex-46c411fb98924fdb8d9e031d1155c31c` | SUCCEEDED, exit 0; 96 + 48+48 steps returned success |
| Independent numerical verification | `codex-5205d65ae5064e78972b68cd1fd274cc` | SUCCEEDED, exit 0; all 8 comparisons have zero differences |

The integration uses the same NEC IsoGSM inputs and PFT+BGC+CROP+MaterialTransport
configuration as phase 1, with isotope fractionation, CL, sediment, levees,
bifurcation, reservoir, flood feedback and dynamic wetlands, 16 MPI processes.
Dates are 2005-06-01 to 2005-06-03, dt=1,800 seconds: 96 continuous steps and
48 + restart + 48 steps. Compiler flags match phase 1 and retain bounds checks;
global FPE traps remain disabled because the same HDF5 initialization probe
raises invalid (documented in the phase-1 report).

Reference output is the completed phase-1 result from
`codex-206b3f72fae9489691bc150a1e1be58b`, including its successful MPI-segment
markers. Strict verification requires complete file/variable/time coverage,
finite non-missing floating values, raw-bit equality, and identical semantic
attributes. It includes all five final restart files, including the history
accumulator sidecar. Every check passed, including finite values and successful
MPI-segment exits.

| Comparison | History variables | Restart variables | Differences |
| --- | ---: | ---: | ---: |
| Phase 1 vs phase 2, continuous | 505 | 1425 (all dates/constant state) | 0 |
| Phase 1 vs phase 2, split | 505 | 1425 (all dates/constant state) | 0 |
| Phase 1 continuous vs split | 505 | 481 (final state) | 0 |
| Phase 2 continuous vs split | 505 | 481 (final state) | 0 |

All eight comparisons have exact nonempty file, variable, dimension and time
coverage, matching masking/scaling semantics and raw data bits. Global creation
times and output-path attributes are excluded. These 5,832 variable comparisons
include repeated variables across comparisons; they are not independent
scientific calibration tests.

Fetched summaries, JSON counts, full logs, configs, binary hashes and terminal
verification-job metadata are retained under
`.omx/artifacts/material-host-cleanup-phase2-20261031/validation/`.
The MT-off build evidence is under `mt-off/` alongside it. Large NetCDF outputs
remain in their private T7920 job directories.

## Limits

The scientific claim is equivalence for this fixed two-day case and 16-process
layout. Long integrations, other MPI decompositions, and extreme snow topology
transitions are not covered. Thermal interfaces and snow remapping were not
refactored in this pass.
