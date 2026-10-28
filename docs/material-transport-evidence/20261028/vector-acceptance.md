# Vector / real restart acceptance handoff

## Verdict — PASS (bounded numerical/runtime scope)
Read-only independent verification against source `6d99543861c38c708f6c833e6349bcd3d7d3b879` and actual completed job artifacts. No models rerun; no production or tracked documentation changes. This does not close the whole constituent-transport project.

| Claim | Authoritative evidence |
|---|---|
| Five-day Vector continuation and MPI invariance | `codex-6d036be5acbd4ef5bfec97ab3c0605d3`, succeeded/0. Full5 vs split5, full9, split59: each comparison covers **2 history files / 193 variables** and **115 restart files / 4,690 variables**, exact raw variable bits and complete file/variable sets. |
| Correct full windows and prefixes | 2005-01-01..05, 240 steps per full case; split segments 96+144. All eight segment logs end normally with zero SIG/BAL/hard/host water/energy failure messages. Both split prefixes retain their first two records exactly; four cases have five increasing daily records, finite unmasked material values and **60,551 nonzero CL values each**. |
| Real half-day count persistence | Jan1 noon split, NP5→9: 24+24 steps; six actual NetCDF blocks / **2,952 patches**, `hacc_nac=24`, valid-count values within 0..24, version1. All **193 variables** at the completed full-day endpoint equal continuous output. The extra partial endpoint is explicitly accounted: 55225440 vs full endpoint55226160 minutes since1900. |
| Actual dynamic state after half-day split | Supplemental read-only job **`codex-ca35eb4112334ef991fe7ad03f6b8d00`**, succeeded/0: explicit **2005-002-00000**, **18 files / 750 variables** completely raw-bit identical. See correction below. |
| Actual final-code legacy compatibility | `codex-005e25fe7c6948e69428ca773dcc762b`, succeeded/0. Normal/legacy each run 72 steps (Jan1 noon..Jan2 end), no conservation/host failures. **61 non-history-accumulator restart files / 3,136 variables**, plus **2 history files / 192 variables** excluding intended `f_frcsat` difference, raw-bit identical. Legacy first-day `f_frcsat` is all missing; second day data and mask recover exactly. |
| Actual corrupt sidecar rejection | Same IO job: orphan count, missing count and mixed-NaN marker each exit1 before any model step, with the respective required fatal text; no END. Script checks all MPI model processes have exited after every case, not merely any nonzero return. Shared seed content manifest unchanged. |

## Corrected verification defect — preserve this history
Original `p9t2/validate.py` selected `max(all restart directory names)` and therefore selected **`const`**, not a dated checkpoint. Its original **13 files / 364 variables** result must **not** be cited as dynamic-state proof. `p9t2-validate-original.py` preserves that exact validator and matches the executed script hash. The local `.tmp/ct/p9t2/validate.py` now filters date-form directory names and asserts `2005-002-00000`; `py_compile` passes. The original remote job is untouched. New supplemental analysis explicitly compares the correct dynamic checkpoint and passes (**18/750**, above), without rerunning physics.

## Source traceability
- Normal Vector runtime HACC SHA: `f4b1b773f71a4ad1d99792a794124731b06117ec7f3296a81cb6c2c3f6077d30`.
- Final/real-IO HACC SHA: `7a3dbb0e7413c58a8c1441af2696c677d5e33ba56214e3f6896591e8017ef4d2`; rebuilt binary SHA `6bdcf1527d329789a32219fb71ef796780a46557c3cb88b472f62039c6f91777`.
- Other five relevant production/generator SHA values match final source byte-for-byte. The only HACC difference is the additional count-without-version metadata probe/rejection, recorded in `final-invalid-guard.diff`; normal arithmetic is unchanged.
- Executed remote drivers/validators match local pre-correction scripts. Actual statuses, JSON results, log date/error counts and source hashes are preserved in `vector-status.json`, `vector-evidence.json`, `io-status.json`, `io-evidence.json`, `checkpoint-status.json`, `dynamic-checkpoint.json`, `actual-fatal-excerpts.txt`. `collect.py` reproduces the small evidence collection; no bulk NetCDF/build download was used.

## Limits / do not overclaim
UNSTRUCTURED/PC/VG, materials on, static isotopes (fractionation off), nonzero CL, BGC and lateral river routing off. No isotope forcing, LULCC, scientific/observational accuracy, or broad configuration matrix is proved here. Legacy format is an actual NetCDF/MPI test made by removing the two new field names from a valid seed, not an archived old forcing-schema restart. All raw-bit claims exclude attributes/file encodings. IO physical comparison intentionally excludes history-accumulator sidecars and one diagnostic; no claim those intentionally different windows are identical. Existing SinglePoint same-month buffered overwrite and empty-worker mapping boundaries are unchanged. Full latest-code multi-day MPI matrix was not rerun solely for the final invalid-input guard; final-code actual MPI IO paths were exercised.

## Repository archive names

The original private filenames above are retained as provenance. In this Git
archive, `vector-evidence.json` is `vector-original-evidence.json`,
`io-evidence.json` is `vector-io-results.json`, and `dynamic-checkpoint.json`
is `vector-dynamic-checkpoint.json`. Status files have the corresponding
`vector-job-`, `vector-io-job-`, and `vector-checkpoint-job-status.json` names;
fatal excerpts are in `vector-actual-fatal-excerpts.txt`.
Unarchived scripts remain in the stated private directories or the respective
runner's frozen `work/`; no bulk result files are included here.
