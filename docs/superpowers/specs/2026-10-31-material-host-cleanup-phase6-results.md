# Phase 6: consolidate material stages in existing modules

Baseline is the complete verified phase-5 working tree, frozen before edits at
`/tmp/colm-host-phase5-t4jdp4xy`; durable baseline and SHA256 manifest are retained
under `.omx/artifacts/material-host-cleanup-phase6-20261031/baseline`.
All previous uncommitted refactors are preserved; this report describes only phase 6.

## Production changes

- MAIN makes one water-stage call. Existing LandPhase water entry now owns the
  original irrigation sum, water-flux correction, flood material extraction,
  wetland/soil transport dispatch and flood publication, in that exact order.
  Its automatic flood species array is local to MATERIAL. MAIN no longer imports
  nmts or owns that material-only array. Flood water credit keeps its allocatable
  dummy contract and all host water/forcing arguments remain read-only.
- Existing LandPhase thermal entry now performs snowmelt carry and saves the
  postthermal water snapshot after evaporation transport. MAIN's lb is snl+1
  through this stage; array lower bounds and snapshot data are preserved.
- SnowHost adds one normalization entry using the existing combine/divide routines.
  It checks snl after combine and calls divide only if snow remains, preserving
  combine lb:1 versus divide lb:0 and interface-level legacy entries.
- Makefile.material adds the direct LandPhase -> RiverLakeCoupling dependency.
  No new files implementing algorithms, packages, configuration or restart/history
  schemas are introduced. Host solver algorithms are unchanged.

MAIN: **1883 -> 1859 physical lines**. LandPhase adds 25 lines, SnowHost adds 22,
reflecting moved code and the unified entry's interface declarations. Total production
Fortran changes are +23 physical lines; the goal here is reducing host exposure.
By the previous residual audit's counting rules, material-related host code/interface
lines decrease **429 -> 415**, and related CPP lines **654 -> 648**.
Exact diffs, per-file counts and final production hashes are in phase-6 artifacts.

## Validation

- New thermal/water-stage proof: 8 configurations x 256 stage pairs,
  **3,427,840 bytes identical** to frozen phase 5. Actual old/new MAIN calls and
  LandPhase entry bodies execute with real canopy partition, snowmelt carry and
  work snapshots. Recording evaporation/soil/flood services check all argument bits,
  optional presence and ordering. Covers material counts 0/3, 0–3 snow layers,
  soil/wetland, dynamic wetland, Grid/CROP, O0/O2, stale snapshots, cancellation
  in the irrigation sum, and allocated/unallocated flood credit.

- Existing bookkeeping/lifecycle regression: 384 cases plus five CaMa guards pass,
  with the same six trace hashes as the pre-edit run.
- Actual standard/SNICAR snow adapters: original 960 cases / 1,107,200 bytes and
  new normalized-stage 576 cases / 664,320 bytes match their frozen references.
  Includes zero/one/three materials, two bounds, empty/collapsed snow and patch isolation.
- Routing/driver/history/publication reference comparison: 104 cases / 6,720 bytes
  match pre-phase-4 code; publication is now extracted from its new LandPhase owner.
- Actual snowmelt carry: six scenarios pass. Flood input: nine scenarios plus
  rejection without remaining water credit pass. Flood evaporation: seven pass.
- Material-off preprocessing: 12 macro configurations x three changed Fortran files
  = 36 token/statement comparisons; all exactly equal to frozen phase 5.
- Legacy physics-boundary argument mapping: 64 cases pass after adapting stage
  preparation. Phase-5 optional-scaffolding proof: 8 configurations, 128 water
  branches, 96 thermal lifecycles, 72 history calls, 32,288 bytes still identical.
- Generated history registry still matches all 525 fields.
- Separate production review: four changed production/build files, zero findings.
  Reviewer independently repeated MT-off and snow-stage comparisons successfully.

## Build and verification evidence

Builds run on T7920 private source `/tmp/colm-host-phase5-sykm6gn3` snapshot,
verified byte-identical to final local production files. Each mode compiles and
links colm.x, mkinidata.x and mksrfdata.x. Job terminal states, headers, build logs
and binary hashes are retained under the phase-6 artifact directory. All jobs
**SUCCEEDED, exit 0**:

| Mode | Job ID | Linked programs |
| --- | --- | --- |
| GRID+PC+BGC+MT on | codex-d5f46c36448c45c9930dd31558c46b52 | 3 |
| GRID+PC+BGC+MT off | codex-78790649b215400fabaefa6111057261 | 3 |
| GRID+PFT+BGC+CROP+MT on | codex-5276f362aff349a3a9685f9480f879d0 | 3 |

## Limits

The new phase orchestration uses changed internal LandPhase signatures; all
repository callers are updated. External callers of these material facade routines
must add the new phase inputs. Existing SnowHost combine/divide entry signatures
remain available. No full continuous/split model simulation was rerun in this pass;
targeted bit comparisons and full compilation do not imply a new complete-model
numerical equivalence result. Phase 4's complete-model verification remains prior
historical evidence only.
