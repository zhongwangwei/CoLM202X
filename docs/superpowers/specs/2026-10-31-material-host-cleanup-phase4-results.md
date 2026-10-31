# Comprehensive material-host migration: phase 4

The user requested a complete migration batch before unified testing. All
production edits were completed against a frozen verified phase-3 working tree;
no per-module long model runs were performed. Earlier uncommitted changes remain.

## Scope and footprint

28 host Fortran files changed; the host footprint shrinks by **1300 lines**.
Two complete material-only preprocessing source files (768 lines) are included in
that reduction; their implementations are relocated, not deleted.

| Host file | Before | After | Change |
|---|---:|---:|---:|
| `main/BGC/MOD_BGC_CNBalanceCheck.F90` | 285 | 287 | +2 |
| `main/BGC/MOD_BGC_CNNStateUpdate3.F90` | 288 | 283 | -5 |
| `main/BGC/MOD_BGC_CNSummary.F90` | 1137 | 1140 | +3 |
| `main/BGC/MOD_BGC_CNZeroFluxes.F90` | 853 | 850 | -3 |
| `main/BGC/MOD_BGC_Soil_BiogeochemNLeaching.F90` | 197 | 178 | -19 |
| `main/BGC/MOD_BGC_Vars_1DFluxes.F90` | 1153 | 1151 | -2 |
| `main/CoLM.F90` | 831 | 787 | -44 |
| `main/CoLMDRIVER.F90` | 412 | 391 | -21 |
| `main/CoLMMAIN.F90` | 1960 | 1888 | -72 |
| `main/HYDRO/MOD_Grid_RiverLakeHist.F90` | 747 | 743 | -4 |
| `main/HYDRO/MOD_Hydro_SoilWater.F90` | 3701 | 3680 | -21 |
| `main/MOD_GroundFluxes.F90` | 278 | 281 | +3 |
| `main/MOD_GroundTemperature.F90` | 466 | 468 | +2 |
| `main/MOD_Hist.F90` | 5747 | 5670 | -77 |
| `main/MOD_HistAccRestart.F90` | 1641 | 1499 | -142 |
| `main/MOD_LeafInterception.F90` | 781 | 763 | -18 |
| `main/MOD_LeafTemperature.F90` | 1651 | 1650 | -1 |
| `main/MOD_LeafTemperaturePC.F90` | 2214 | 2212 | -2 |
| `main/MOD_PhaseChange.F90` | 861 | 857 | -4 |
| `main/MOD_SoilSnowHydrology.F90` | 2611 | 2601 | -10 |
| `main/MOD_Thermal.F90` | 1572 | 1542 | -30 |
| `main/MOD_Vars_1DAccFluxes.F90` | 3548 | 3537 | -11 |
| `main/MOD_Vars_1DFluxes.F90` | 314 | 317 | +3 |
| `main/MOD_Vars_TimeVariables.F90` | 1719 | 1688 | -31 |
| `mksrfdata/Aggregation_LakeSoilC.F90` | 202 | 0 | -202 |
| `mksrfdata/Aggregation_MethanePH.F90` | 566 | 0 | -566 |
| `mksrfdata/MKSRFDATA.F90` | 606 | 595 | -11 |
| `mksrfdata/MOD_RegionClip.F90` | 1082 | 1065 | -17 |

Nine cohesive material modules now own thermal/canopy diagnostics, snow adapters,
NH4 state/operations, methane aggregation/clipping, runtime orchestration, history
state and material restart metadata. Existing material modules own the remaining
moved blocks. Material config/dependency rules live in MATERIAL/Makefile.material.
No external package dependency, user namelist or restart/history schema was added.

## Boundary audit

Independent review inspected 46 changed source/build files and all32 remaining
host/config entries. No large material-only implementation remains outside its
subsystem. Remaining host content is justified:

- Exact acquisition/event hooks in MAIN, thermal, interception, water and snow.
- Optional diagnostic forwarding at physical phase boundaries; snow layer hooks
  remain at each original merge, split, shift and collapse.
- Runtime, routing, LULCC and BGC lifecycle/callback locations, including exact
  NH4 allocation/zeroing order. Generic host water/BGC algorithms stay in place.
- Generic history/restart engines and clipping callbacks. Compiler compatibility
  guards remain at the configuration entrypoint. share has only a comment;
  mkinidata has no residual material implementation.

The two complete methane aggregation procedure bodies are byte-for-byte equal
to the frozen versions (198 and562 procedure lines). Their surrounding guards
and entry dispatch remain equivalent.

## Focused validation

All checks below pass against frozen production source or independent assertions:

- 4 thermal configurations (Campbell/VGM × O0/O2), 1728 real phase calls each,
  490824 compared bytes per configuration, plus diagnostic reset/PFT SUM checks.
- 960 snow adapter calls,1107200 bytes exact;18720 snow topology calls,
  16074240 bytes exact, full optional/bounds/0-species/repeated-call coverage.
- 290 real canopy calls and9 guards,92944 bytes exact;684 soil-column real values.
- 108 NH4 cases,17280 bytes exact;32 clipping cases,19448 callback bytes exact.
- 104 routing/driver/publication cases,6720 bytes exact;16 runtime/historycount
  cases and lifetime/error/empty-rank assertions.
- 384 phase-bookkeeping cases+5 CaMa guards;64 land boundary cases;18 real
  NetCDF/MPI history compatibility cases at1/4 ranks;525 generated history fields.
- Existing material physics48 numerical checks+5 wiring checks, snowmelt6,
  flood-input9+error guard, flood-evaporation7, canopy isotope mixing,
  naming inventory and narrow upstream BGC exception checks.

## Unified remote validation

Source snapshot: /tmp/colm-host-phase4-ygewzvwc; SHA256 manifest and frozen baseline
are in .omx/artifacts/material-host-cleanup-phase4-20261031/.

- PFT+BGC+CROP+MT build,96-step continuous and48+48 restart:
  codex-2e7181c3a4674af0854b2c217a1a736c.
- PC+BGC+MT build: codex-c00cca0def5b490c88463c5346bc000d.
- IGBP,noBGC,LULCC+MT build: codex-40033e33f8784a90af9c921928d09493.
- PC+BGC,MT-off build: codex-5023378cc04340c5beddb7e7d1de153b.
- SinglePoint,IGBP,noBGC+MT build: codex-4406d54df4c24dedb0ee947250226488.
- UNSTRUCTURED,PC+BGC+MT build: codex-8171dbb996f84e3db70e8f57135f8927.

The two added spatial modes use /tmp/colm-host-phase4-3edkebsn, whose production
hashes are identical; only the matrix runner adds those mode choices.

All six build configurations succeeded for all three executables. Five build-only
jobs have terminal succeeded/exit0 evidence, fetched headers/build logs and three
binary hashes each; matrix-results.json records them. The final PFT model job
passed its build subprocess and all96/48/48 model segments (terminal exit0).
The final PFT job completed 96 continuous steps and 48+48 restart steps normally.
Strict verification job `codex-89e0eccacd534da7af379093ad774757` succeeded with exit 0:
**all eight groups, 5,832 variable comparisons, zero differences**. Full file sets,
variable/dimension/encoding semantics, two expected daily records, material and
sediment fields, and finite valid floating-point values all pass.

| Comparison | Variables | Differences |
|---|---:|---:|
| revisions-continuous-history | 505 | 0 |
| revisions-continuous-restart | 1425 | 0 |
| revisions-split-history | 505 | 0 |
| revisions-split-restart | 1425 | 0 |
| baseline-split-history | 505 | 0 |
| baseline-split-restart | 481 | 0 |
| current-split-history | 505 | 0 |
| current-split-restart | 481 | 0 |

Complete evidence is under
`.omx/artifacts/material-host-cleanup-phase4-20261031/validation/`: comparison JSON,
eight detailed logs, model/build logs, input namelists, actual compiler header,
three binary hashes and source-change records. Separate matrix directories retain
all five build-only terminal-success reports. Production source hashes still match
the final model snapshot.

An early whole-output check of the completed continuous run also passed while
restart segments ran (job `codex-307a06fa1a0a42898712514a641f29aa`); its evidence is
under `continuous-check/`. No numerical mismatch was found in any check.
The numerical baseline is completed phase-3 job
`codex-00988b4f47864faba508af0815e1aa3e`.

The first MT-on build attempt exposed an omitted explicit `nmts` import in MAIN.
It was fixed without changing arithmetic. Initial failed jobs e46b9409/f25d4096/
7f28b69d never reached a model run. The replacement PFT job reuses the completed
private object cache from its failed compile and rebuilds the changed source;
PC/LULCC replacements are clean builds. The MT-off job used the preceding source,
but the only follow-up import is inside MT guards, leaving effective MT-off code
identical. No completed model run was repeated for this fix.

MT-off effective-code checks covered26 changed existing host files ×3 CPP variants:
75/78 were identical; only HistAccRestart differs because its always-available
legacy-name function moved behind an import alias. Its body matches after renaming,
and the real NetCDF/MPI cases exercise the migrated function.

## Limits

The actual coupled fixture spans two June days, one NEC cell and16 MPI ranks,
with isotope fractionation,CL andsediment enabled. Dedicated fixtures cover
cold/snow/NH4 cases absent from that run. Bulk surface-data generation and longer
simulations/other MPI layouts are not executed. Internal optional diagnostic APIs
were intentionally grouped; repository callers are converted, but external callers
using old material-only keyword parameters must adopt the new diagnostic object.
