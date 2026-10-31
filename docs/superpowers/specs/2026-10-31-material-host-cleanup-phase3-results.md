# Material host cleanup: phase 3 results

## Changes

This round targets upstream host modules outside CoLMMAIN. It preserves the
uncommitted phase-1 and phase-2 changes. Baseline sources and SHA256 hashes are
under `.omx/artifacts/material-host-cleanup-phase3-20261031/baseline/`.

- Soil-water exchange diagnostic helpers move into existing HydroDiagnostics.
  Their bodies and solver observation positions remain unchanged.
- Snow material remapping moves into a Precision-only material module; the four
  public standard/SNICAR interfaces and host water algorithms remain unchanged.
- New object and direct host dependency are wired in Makefile.

## Verification

Targeted checks completed:

- Snow host: 1901 → 1550 lines, **351 fewer host lines**. The shared helper is
  212 lines, so snow host plus helper also shrink by 139 lines overall.
- Soil-water host: **118 fewer lines**, retaining the two unchanged helper bodies
  in HydroDiagnostics. Combined host reduction outside CoLMMAIN: **469 lines**.
- 18,720 actual snow calls and 684 soil-column real outputs are bitwise equal
  to frozen pre-phase-3 sources. Snow comparison covers 16,074,240 output bytes,
  including water and eight aerosol arrays, plus finite/budget assertions.
- Both MT-off host modules are identical after preprocessing, ignoring blank
  and comment-only lines. Existing diagnostics, snowmelt, flood input/evaporation,
  384 phase-bookkeeping cases, five CaMa guards, 48 material physics cases and
  naming inventory checks pass. `git diff --check` passes.

Independent final review: APPROVE, zero remaining findings.
MT-on compile/link and all three model segments completed successfully (runner
terminal state succeeded, exit 0). MT-off compile/link completed successfully
(terminal state succeeded, exit 0). Strict comparison succeeded (terminal state succeeded, exit 0): all eight
comparison groups have zero differences, covering **5,832 variable comparisons**.
This includes revision-to-revision continuous/split history and all restart
files, and continuous-versus-split history/final restart in each revision.
All valid floating-point outputs are finite; file sets, variables, dimensions,
time records, material/sediment coverage and encoding semantics pass.

Fetched evidence: `.omx/artifacts/material-host-cleanup-phase3-20261031/validation/`
contains `results/comparison.json`, eight detailed comparison logs, build/model
logs and configs; MT-off evidence is in the adjacent `mt-off/` directory.
Current production source hashes still match the verified immutable snapshot. Numerical baseline is completed phase-2 coupled job
`codex-46c411fb98924fdb8d9e031d1155c31c`, with MT-off build cache
`codex-29b6d1aefad248499796634ea28e3e33`.

Remote jobs:

- MT-on build and continuous/split model: `codex-00988b4f47864faba508af0815e1aa3e`.
- MT-off compile/link: `codex-aa2a33eecd3a43dcaabf4df2c7f1a545`.
- Strict numerical verification: `codex-47b47a0890e94411af96885a5760b3fb`.
- Immutable source snapshot: `/tmp/colm-host-phase3-sxwnh10z`, with SHA256 manifest.

The initial MT-on job `codex-180ca117373e4ad7bee89860a92dfb87` was cancelled
before verification because review requested private visibility for helper imports.
The replacement above builds the final source. MT-off job predates that
visibility-only change, whose entire block is removed with MT disabled; effective
MT-off source equality was checked after the fix.

## Scope limits

The coupled fixture covers two June days, one NEC cell, PFT+BGC+CROP and 16 MPI
ranks. Targeted tests exercise snow cases absent from that warm coupled fixture.
Longer simulations and other MPI layouts remain outside this round.
