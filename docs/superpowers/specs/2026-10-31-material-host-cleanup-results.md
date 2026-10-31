# Material host cleanup validation (2026-10-31)

Baseline: `e766527a` (the working tree was clean before this task).

## Changes

- `CoLMMAIN`: 2,200 -> 2,052 lines. Material snapshots/scratch are held in a
  patch-call-local `mt_land_work`; water diagnostics in `mt_water_diagnostics`.
- `WATER_VSF`: 18 material diagnostic arguments -> one diagnostic object.
  `soil_water_vertical_movement`: seven diagnostic arguments -> the same object.
- The soil solver resets only its own root/runoff fields, retaining snow,
  top-boundary, evaporation and water-balance diagnostics sampled by WATER_VSF.
- Soil/static-wetland material dispatch lives in `mt_land_water_step` in the
  existing LandPhase module. The low-level material routines are unchanged.
- Thin-snow carry bookkeeping and flood-infiltration material calculation moved
  to the existing Snow and RiverLakeCoupling modules, at the original call times.
- No new dependencies, namelist settings, history variables or restart fields.
  Snow-layer topology remapping is outside this pass.

## Local evidence

| Check | Result |
| --- | --- |
| Existing runoff handoff regression | PASS: 9 numerical + 5 wiring checks, before and after |
| Existing isotope physics regression | PASS: 15 numerical checks, before and after |
| Existing solute ice recovery regression | PASS: 24 numerical checks, before and after |
| Real soil solver, 6 columns x 3 steps | PASS: 684 host/diagnostic real outputs bitwise equal to baseline |
| Diagnostic lifetime and reset | PASS: full/partial reset, resize, negative/empty snow bounds, object isolation |
| Land boundary argument forwarding | PASS: 64 combinations, all argument values and optional presence bitwise equal to baseline |
| Thin-snow helper | PASS: 6 cases, mass conservation, species/patch isolation, zero materials |
| Flood-input helper | PASS: 9 independent expected-value cases, 9 bitwise baseline comparisons and both rejection paths |
| MT-off preprocessed host source | PASS: CoLMMAIN, SoilSnowHydrology and Hydro_SoilWater unchanged after removing comments/blank lines |
| Python syntax and git diff whitespace | PASS |
| Independent production review | No blocking findings; identified Git-dependent test was corrected |

Commands and coverage details are in `tests/material_host_cleanup/README.md`.
Baseline source files were exported before production edits; optional baseline
arguments make the numerical comparison reproducible without embedding old
production implementations in tests.

The initial naming inventory failed on upstream `conc_trcr_c` / `conc_trcr_n`.
The authorized follow-up now exempts only those two identifiers in their original
BGC file; that production file remains identical to upstream. The inventory and
its new exception regression pass (same names elsewhere and additional legacy
names are still rejected). The history-I/O regression also passes all 18 cases
across one and four MPI processes.

## T7920 full builds

Private snapshots of the edited local source; 16 CPUs per job, 1,800 s timeout.
Both use `tests/run_grid_transport_build.sh`, compile and link `colm.x` and
`mkinidata.x` with bounds and floating-point runtime checks, GRID + IGBP-PC + BGC.

| Configuration | Job | Result |
| --- | --- | --- |
| MTON | `codex-e3d571ce51e94843943b66b04de85182` | SUCCEEDED, exit 0 |
| MTOFF | `codex-c329d424e21d41ca953aae194e1124ae` | SUCCEEDED, exit 0 |

Fetched build logs, generated headers, binary SHA256 sums and terminal job
metadata are under `.omx/artifacts/material-host-cleanup-20261031/mt-on/` and
`mt-off/`. Both report `PASS complete GRID+BGC ... compile and link`.

## Limits

The focused solver test executes real hydraulic/root/aquifer physics. The land
boundary test records arguments through the actual new wrapper and MAIN call;
its material consumers are stubs, so it proves forwarding, not coupled physics.
The coupled follow-up below adds full CoLMMAIN integration and continuous versus
restarted runs. Coverage remains one two-day fixture and one 16-process layout;
no long scientific integration, alternative MPI decomposition, extreme-snow
campaign or performance benchmark was run. Snow-layer remapping is unchanged.

## Coupled follow-up

The input preflight confirmed existing, readable T7920 CASE, RAW, RT and ISO
directories. No external datasets were uploaded. A 12.6 MB code-only snapshot
was assembled from the current local code, retaining a SHA256 source manifest.
Baseline mode restores the changed production files from `e766527a` with their
own SHA256 manifest; both variants use identical compiler flags and parameters.

- GRID + PFT + BGC + CROP + MaterialTransport, 16 MPI processes, bounds checks.
- NEC 125--127.5 E, 47.5--50 N, 0.25 degrees, IsoGSM forcing.
- Water isotopes with fractionation, CL, sediment, levees, bifurcation, reservoir,
  flood feedback and dynamic wetlands.
- 2005-06-01 to 2005-06-03: continuous 96 steps, and 48 steps + restart + 48 steps.
- Strict comparison covers nonempty complete file sets, groups, dimensions,
  variables, raw floating-point bits, masking/scaling and time semantics.
- Extra fixture checks require both expected daily records, each day's valid
  isotope/CL/sediment output, and the five expected final restart files (including
  the history-accumulator sidecar). All non-missing floating-point values in all
  history/restart files must be finite. Runner markers prove each MPI command
  returned success, beyond merely printing a model completion message.

| Check | Job | Result |
| --- | --- | --- |
| Strict comparator, 23 microchecks | `codex-101e35532eb34c25a34aaf748eccaebb` | SUCCEEDED, exit 0 |
| Historical fixture validation | `codex-e099e8d131204aa29451608954206c0a` | SUCCEEDED, exit 0 |
| Initial baseline build + init | `codex-a93ccd4959e54e9e9c176634efdfeb2c` | Build passed; HDF5 invalid trap during init |
| Initial edited build + init | `codex-62355f49cfbf403da8c5f51bd7c7fdee` | Build passed; same HDF5 invalid trap |
| Baseline coupled continuous + split retry | `codex-353570a8be144589b3c893416d005aa3` | All model segments passed; old harness comparator launch failed |
| Edited coupled continuous + split retry | `codex-206b3f72fae9489691bc150a1e1be58b` | All model segments passed; old harness comparator launch failed |
| Independent full output verification | `codex-582152fefae644728e1b49f8d43297f5` | SUCCEEDED, exit 0; all 8 comparisons have zero differences |

The first attempts failed before model physics at `nf90_open(pixel.nc)`.
A minimal same-library probe (`codex-46f398d76f9a4e8a96a345612d00e606`, succeeded)
located SIGFPE in HDF5's `H5T__init_native_float_types`. Opening the same input
without global invalid traps succeeds and reads edges=47.5; merely relinking
without traps does not remove GFortran's startup trap setup. Retries compare
all compiled source files against the read-only previous build, copy its objects
privately, then recompile both program entry points and relink without global
FPE traps, retaining bounds checks. No production logic was changed for this.
Probe logs are preserved in the local artifact directory.

The frozen retry harness invoked its comparator with system Python, which lacks
NumPy/netCDF4. Both model commands had already returned success for all 96 and
48+48 steps; the wrapper jobs therefore ended failed only at comparator launch.
The delivered runner now leaves comparisons to the independent verifier using
the scientific Python. The verifier independently checks all successful-exit
markers, complete output coverage, finite values and raw bits. Its first attempt
(`codex-0ea5ae581e244a388823395c6449fdcb`) rejected the new history sidecar absent
from the older historical fixture's four-file list; the explicit expected set
was corrected to include the fifth file, and all five are compared.

| Comparison | History variables | Restart variables | Differences |
| --- | ---: | ---: | ---: |
| Baseline vs edited, continuous | 505 | 1425 (all dates/constant state) | 0 |
| Baseline vs edited, split | 505 | 1425 (all dates/constant state) | 0 |
| Baseline continuous vs split | 505 | 481 (final state) | 0 |
| Edited continuous vs split | 505 | 481 (final state) | 0 |

All eight directory comparisons pass with identical file/variable/dimension
coverage and masking/scaling/time semantics. Global creation-time/path attributes
are deliberately excluded. The 5,832 variable comparisons are not independent
scientific calibration tests; they establish numerical equivalence for this case.

Fetched summaries, JSON comparison counts, configs, model logs, binary hashes
and successful verification-job metadata are under
`.omx/artifacts/material-host-cleanup-20261031/coupled-validation/`. Large NetCDF
outputs remain in the two private T7920 retry result directories recorded above.
