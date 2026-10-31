# Host optional scaffolding cleanup: phase 5

Baseline: complete verified phase-4 working tree, frozen before this pass at
`/tmp/colm-host-phase5-_s7lztm7`, with durable source and SHA256 manifest under
`.omx/artifacts/material-host-cleanup-phase5-20261031/baseline`.
Git HEAD was not used because earlier refactors remain uncommitted.

## Changes and scope

- MAIN computes irrigation once at the common join after either water routine.
  Sampling remains after mutable irrigation output arguments; addition order stays
  drip + flood + paddy. The common WATER_VSF argument tail is shared.
- LandWork no longer duplicates thermal diagnostic initialization. Its unique
  production preparation call is followed unconditionally by THERMAL, which still
  initializes its output. Water snapshots/reset and standalone PhaseChange reset
  are retained.
- Eight host files share common argument tails outside material guards: snow layer
  topology, GroundFluxes, LeafTemperature/PC, PhaseChange, GroundTemperature,
  Hydro_SoilWater and SoilSnowHydrology. Twenty duplicate tails removed without
  changing effective interfaces, argument order, declarations or algorithms.
- HistAccRestart calls the material sample-count handler directly in the generated
  field list. Its one-use forwarding procedure is deleted and the generator emits
  the exact previous call and callbacks in the same position.
- No Makefile changes in this pass, no dependencies or schema changes.

Eleven production Fortran files changed, net **34 fewer physical lines** (33 in
host files, one in MATERIAL). MAIN: 1888 → 1883 lines. Using the preceding audit's
physical nonblank/noncomment material-code/interface counting rules, host residual
code: **456 → 429**, material-related CPP: **679 → 654**. These counts include
necessary optional integration boundaries, not 429 removable algorithm lines.
Inventory and exact per-file counts are retained in the phase-5 artifact directory.

## Verification

- Baseline lifecycle check passed before editing. Adapted lifecycle check uses
  THERMAL's real initializer after preparation: 384 cases and 5 CaMa guards pass;
  all six trace hashes match the pre-edit run.
- New focused comparison: 8 configurations (CROP/Grid × O0/O2), 128 water branches,
  96 thermal lifecycles, 72 history calls; **32,288 bytes identical**. Mutations that
  prematurely sample irrigation or swap write/read callbacks are rejected.
- Interface-only changes: 24 macro configurations × 8 files = **192 exact token
  and statement comparisons**, including material on/off, Campbell/VGM,
  IGBP/PFT/PC and BGC/CROP. All equal. A further 36 MT-off comparisons for MAIN,
  HistAccRestart and LandWork also match exactly.
- Actual thermal phase regression: 1728 calls per configuration, Campbell/VGM ×
  O0/O2; each 490,824-byte trace matches the pre-phase-4 reference.
- Actual standard/SNICAR snow remapping: 18,720 cases across layer/material bounds,
  16,074,240 bytes match this pass's frozen baseline. Snow adapters: 960 calls,
  1,107,200 bytes match pre-phase-4 direct host calls.
- Soil-column smoke: six columns × three steps, finite states, balance and root
  pathways pass. History-count compatibility: 16 reference cases pass.
- Real NetCDF/MPI: 18 cases at one/four ranks, including legacy names, malformed
  fields and expected rejection, pass. Identifier inventory and exceptions pass.
- Generated history registry: all 525 fields match; --write is byte-idempotent.
- Independent review of 14 changed files: no issues. Python syntax and
  `git diff --check` pass.

## Build evidence

Clean T7920 builds of colm.x, mkinidata.x and mksrfdata.x run from the private
`/tmp/colm-host-phase5-i66fxlzt` source snapshot, verified byte-identical to final
local production files. Modes and terminal results are recorded in `jobs.json`;
full commands, generated headers, build logs and binary hashes are retained in
phase-5 artifacts. All three jobs **SUCCEEDED, exit 0**, with all three programs linked:

| Mode | Job ID |
| --- | --- |
| GRID+PC+BGC+MT on | codex-e040b3dd867c47618bda06bf56aa308d |
| GRID+PC+BGC+MT off | codex-8b4f1bc3489e46ec9f24f4badd36e06e |
| GRID+PFT+BGC+CROP+MT on | codex-f51c204b784542a9a05799d88c8e82a7 |

## Limits

This pass did not repeat the long continuous/split full-model integration; phase 4's
successful 5,832-variable comparison predates these changes. Current bitwise proof
covers focused boundaries, actual phase/snow routines and unchanged CPP interfaces;
recording stubs replace the large water solvers and history count handler in the
new focused test. Real history I/O and complete program compilation are separate
checks. No new full-model numerical equivalence claim is made for phase 5.
