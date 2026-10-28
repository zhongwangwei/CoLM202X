# p9q read-only validation handoff

No new model simulation. Original model job `codex-1f593daee0854b2587cee991b7996e45` remains failed and untouched.

- First analysis: `codex-5d025b6634c442778093eca0089788fd`, terminal failed. Checker 17 positive/negative fixtures passed; raw-bits comparator fixture passed. Analysis stopped on an over-strict assumption that a split run's saved final nml explicitly retained DEF_LULCC_SCHEME. It does not, while the actual transition log records MEC. This is not a model/physics failure.
- Second analysis: `codex-bc9d9644688a4ea0b68817319dcc242c`. Validates any explicit saved scheme and requires actual MEC/SAT transition log. Records omitted saved scheme as null, never inventing an nml setting. Read terminal status/report before claiming full verification.
- Fetch report using `run-7920 fetch JOBID --out NEW_DIRECTORY`; outputs are under results/report.json, metadata separate.

Checker frozen: `tests/material_transport/check_forcing_cache.py`, SHA256 `6851188840ce13ee878ea291cc2c65dec297a64763322a59b033236fb07d623e`; `py_compile` and `git diff --check` passed. README command section changed. No production edits or commits by this child.

Scope: nine existing p9q runs, segment completion/expected timestep sequence, NBAD=0, real LULCC scheme logs, all later schema-2 caches, four complete restart+history variable-bit comparisons, unmasked finite material/isotope outputs. Does not validate attributes or observational/scientific accuracy. No latest frcsat/RegionClip combined-source LULCC test. Original source provenance is in original_source_versions.txt; this analysis does not independently rehash retained build sources against that manifest.

Cold exception proof: initial land file contains complete=1/schema=5/transport_count=0/aquifer=-1 metadata, not a material pool or forcing cache. MOD_MT_Rest's zero-transport write path commits this empty transaction then returns before forcing writing. Explicit first initial-date alone permits it (or truly host-only file); later missing cache, any real initial material pools, invalid schema/count/commit marker are errors.

Final: second analysis terminal succeeded exit 0; report PASS (9 runs, 8 restart/history comparisons). Fetched at `fetched/results/report.json`; metadata in `fetched/metadata`. No third job needed.
