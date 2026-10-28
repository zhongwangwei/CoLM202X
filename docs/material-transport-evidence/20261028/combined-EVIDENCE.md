# combined10 — latest combined-source acceptance (frcsat + forced isotopes + LULCC + CL + sediment/river)

Status record: 2026-10-29 (Asia/Shanghai). All heavy work ran on T7920 through `run-7920`;
this directory holds only the small JSON/text evidence pulled back to the Mac.

## Verdict

**The required full combination does NOT pass on the latest production source.**
The new frcsat protocol itself passes its isolation and protocol checks, but the requested
combined configuration (IsoGSM isotope forcing + fractionation + LULCC IGBP + CL +
sediment + river routing with levee/bifurcation/reservoir/flood-feedback) aborts on the
LULC_IGBP/NEC build inside the first day:

```
TIMESTEP = 4 | DATE = 2005-12-25-05400
TRC_BAL residual note (excluded from tol check): n_entries=2 worst_abs= 0.90882E-04
TRC_BAL residual hard failure: n_entries=1 worst_abs= 0.90882E-04 tol= 0.50000E-04 @ipatch=213 itrc=1 ptype=0 owner=0
 TRC_BAL hard failure: numerical residual exceeds independent threshold
Abort(1): application called MPI_Abort, exit 1
```

This is the explicitly booked `trc_numerical_residual_step` (`numerical_source_sink`) from the
near-dry/flood material paths in `MOD_MT_SoilWater.F90`, checked against the independent hard
limit `trc_resid_abs_tol = 5.0e-5` in `MOD_MT_Conservation.F90`. It is **not** the ordinary
aggregate `TRC_BAL`/`TRC_SIG` check. No threshold was relaxed and no delta was clipped for this
report.

Single-factor ablation (all on the same frozen source and the same build inputs) localizes the
trigger to one option: **`DEF_GridRiverLake_FloodFeedback = .true.` together with runtime isotope
forcing (IsoGSM `prate1sfc/prate2sfc`, `spfh12m/spfh22m`)**. Removing flood feedback makes the
rest of the combined configuration run clean for 5 days; enabling only flood feedback on the
p9q-like forced-isotope configuration reproduces the same step-4 abort.

Trigger classification (job `codex-69298c1322fd4c749b69e9999c7b74d8`): the same 12-step trigger
case run on the two **pre-built, no-LULCC** p9r binaries — old = `28478b18`
(`results/bld/old/run/colm.x`, sha256 `eb30ab5b…e90e`) and new = `c315c849`
(`results/bld/new/run/colm.x`, sha256 `c41e8853…f3ea`) — aborts in both cases with the **identical**
`worst_abs = 9.0765E-05 > 5.0E-05 @ipatch=213 itrc=1 ptype=0` residual at step 4. Therefore this is
**not a regression introduced by 6d995438 or the recent fixes**; it is a pre-existing interaction
defect of the branch base (present in `28478b18`), only *exposed* here because previous acceptance
work never combined flood feedback with runtime isotope forcing (p9q had no flood feedback; p9r's
flood-feedback runs had no runtime isotope forcing).

## Source freeze (what exactly was accepted)

| item | value |
|---|---|
| repo / branch | `/Users/zhongwangwei/Desktop/Github/CoLM-constituent`, `feature/constituent-transport` |
| production source | `6d99543861c38c708f6c833e6349bcd3d7d3b879`, snapshot = `git archive 6d995438` into `new/` |
| HEAD at freeze time | `10e9fbd204a99ff81bd60a98198b364d489e10d9` |
| HEAD later (docs/tests only) | `74cb9175a239cb3fdb0270c0d1e6eb7e441333b0`; `git diff --name-only 6d995438..HEAD -- . ':(exclude)docs' ':(exclude)tests'` = empty |
| snapshot integrity | 620 regular files; in-job `sha256sum -c source.sha256` rc=0 (`logs/source-check.log` on 7920) |
| tools | `compare.py`, `iso_check.py` = frozen p9q copies (byte-identical to HEAD tests); `check_forcing_cache.py` = delivery HEAD test tool required by the task |
| build | `make all FF=mpif90`, conda GNU/MPI/NetCDF, `create_defineh.bash GRID <LULC> URBANOFF vanGenu CaMaOFF BGCOFF CROPOFF MTON LULCCON` |

## Jobs and attempts (preserve, do not misread)

| # | job id | purpose | terminal | outcome |
|---|---|---|---|---|
| 1 | `codex-2e853460a0984d3c93ce2e6cdaa12b5d` | combined10 main acceptance, attempt 1 | failed/1 | Stopped in the build gate before any model ran: the script check required a literal "#define USEMPI" but create_defineh writes "#define  USEMPI"; fixed the checker (no model was run, no physics result exists for this job). |
| 2 | `codex-1c0b12b7efa44efc9eaf2540062163ae` | frcsat isolation attempt 1 | failed/1 | Both builds OK; runs failed because the helper used `local c=$OUT/runs/$tag bld=$2`, so $tag/$bld expanded empty before assignment and cases landed in the wrong directory; validator then found no case dirs. Fixed. |
| 3 | `codex-6d761e65ddd54ccf987df1d0bc68b2f0` | frcsat isolation attempt 2 | failed/1 | Both builds OK; runs aborted with "No &nl_colm_isotope_parameter in .../iso_off.nml" because the no-forcing isotope namelist file was empty; replaced with an explicit empty namelist group. |
| 4 | `codex-6284ea5dd86e486ab6fa97288b96dff1` | option-subset localization attempt 1 | failed/1 | run_variant was not exported to the xargs subshells ("run_variant: command not found"); no model ran. Fixed with export -f. |
| 5 | `codex-17e7132400c44adeb691e66e3dc2c586` | option-subset localization round 1 (9 variants) | succeeded/0 | PASS as a diagnostic job. Isolated the trigger to DEF_GridRiverLake_FloodFeedback=.true.; V4-fb and every variant containing flood feedback aborted at step 4. |
| 6 | `codex-2434412927a94e38a44ad21db5643a38` | combined10 main acceptance, attempt 2 | EXPECTED_FAILED_1 | All six IGBP/NEC main cases abort at TIMESTEP=4 with TRC_BAL residual hard failure 9.0882e-5 > 5.0e-5; f-* PFT and i-AMZ cases run on. Runner state will be failed/1 because the validator gate fails; this is the acceptance verdict, not an infrastructure error. |
| 7 | `codex-05c172686e434105884e3a0c7297f6fb` | frcsat fix isolation: c315c849 vs 6d995438 on an identical no-LULCC/no-forcing MT case | succeeded/0 | PASS: 120 steps each, history identical raw bits, all 13 restart files and 767 common variables identical; only the intended two new sidecar fields at the mid-window checkpoint. |
| 8 | `codex-0c007acf2ceb4863aebe3d3dbf812cd2` | ablation round 2: fractionation / LULCC scheme / long no-feedback control | succeeded/0 | W2 (fractionation off) and W3 (SAT) still abort at step 4 -> fractionation and scheme are not the trigger. W4 (full config without flood feedback, 5 days) passes with zero diagnostics. W1 was a namelist configuration error (IsoGSM isotope-forcing groups left in a JRA-forced run), redone correctly in job 9. |
| 9 | `codex-3c6c916b96a142f6a17f88cbfe25c374` | ablation round 3: JRA/static-isotope full config and LULCC-off build | succeeded/0 | X1 (JRA dynamics + static isotopes + full river options incl. flood feedback + CL + sediment, 3 days) passes; X2 (flood feedback + IsoGSM isotope forcing on the LULCC-off build) still aborts at step 4 -> isotope forcing is required for the failure and the LULCC compile switch is irrelevant. |
| 10 | `codex-69298c1322fd4c749b69e9999c7b74d8` | trigger regression classification on pre-built p9r binaries (old=28478b18, new=c315c849) | succeeded/0 | Both binaries abort identically at step 4 with the same residual -> the interaction is pre-existing at the branch base 28478b18, not introduced by the latest commits. |
| 11 | `codex-81a0f423c967451abfcb5adb5ef28881` | read-only supplementary audit (corrected iso_check path, field-level forcing-cache audit, frcsat boundary counts) | failed/1 | Two explained audit flags (cold-start exemption; isotope-only positivity assumption of the delivery cache tool) plus the expected no-history result for the aborted i-NEC-MEC case; supplies the corrected subset facts used in the gate discussion. |

Every job except the main one was a bounded probe/ablation using existing builds; none of them
re-ran the failed main configuration in the hope of a different result.

## Main acceptance job `codex-2434412927a94e38a44ad21db5643a38`

### Gate results as reported by the main job

Job: `codex-2434412927a94e38a44ad21db5643a38` — terminal state `failed/1` (validation gate verdict, not infrastructure).
Full machine report: `evidence/main/report.json`; log: `evidence/main/validation.log`.

| gate | result | key evidence from report.json |
|---|---|---|
| G1 两种构建（LULC_IGBP / LULC_IGBP_PFT，MT+LULCC+GridRiverLakeFlow，BGC/CaMa/CROP off） | **PASS** | no errors recorded |
| G2 组合配置成立 + 连续/拆分 restart/history 全文件全变量原始位一致 | **FAIL** | full history comparison i-NEC-MEC vs i-NEC-MEC-s0101 failed: FILESET a=0 b=0 unmatched=0; full history comparison i-NEC-MEC vs i-NEC-MEC-s0103 failed: FILESET a=0 b=0 unmatched=0; noon-split restart alignment: expected exactly the mid-window checkpoint 2006-001-43200 extra, got [] |
| G3 全部 END / 步数时间戳精确 / NBAD=0 / 无 SIG/BAL/LULCC_BAL/residual hard/STOP/abort / 无未掩码 NaN/Inf | **FAIL** | i-NEC-MEC: log.colm: no "CoLM Execution Completed"; i-NEC-MEC: log.colm: diagnostic "residual hard failure" x1; i-NEC-MEC: log.colm: diagnostic "MPI_Abort" x1 |
| G4 每个 restart 陆地块 schema-2 forcing cache 有效（check_forcing_cache.py） | **FAIL** | i-NEC-MEC: forcing cache check failed: AssertionError('no initialized forced checkpoints'); i-NEC-MEC-s0101: forcing cache check failed: AssertionError('no initialized forced checkpoints'); i-NEC-MEC-s0103: forcing cache check failed: AssertionError('no initialized forced checkpoints') |
| G5 frcsat 有效样本协议（history 掩码/计数一致、mid-window 计数版本保存/边界重置/恢复、部分窗口重算） | **FAIL** | i-NEC-MEC: history f_frcsat has no valid sample at all (cannot claim the protocol produced real data); i-NEC-MEC-s0101: history f_frcsat has no valid sample at all (cannot claim the protocol produced real data); i-NEC-MEC-s0101: 2006-001-00000: history sidecar missing at the split checkpoint |
| G6 真实 MEC/SAT 过渡 + isotope/CL 账本诊断为零 + tracer 注册 | **FAIL** | i-NEC-MEC: LULCC transition banner count 0 != 1; i-NEC-MEC: transition log does not contain "(MEC)"; i-NEC-MEC: no scheme-specific transition line found |
| G7 同位素 δ 宽范围诊断（不作观测准确度判据） | **FAIL** | i-NEC-MEC: iso_check rc=1; f-NEC-MEC: iso_check rc=1; i-AMZ-MEC-cl2e5: iso_check rc=1 |

Key counts (from the same report):

| case | steps | segments | final timestamp | diagnostics | history files | restart files |
|---|---|---|---|---|---|---|
| i-NEC-MEC | 4 | [4] | 2005-12-25-05400 | 2 | 0 | 4 |
| i-NEC-MEC-s0101 | 4 | [4] | 2005-12-25-05400 | 2 | 0 | 4 |
| i-NEC-MEC-s0103 | 4 | [4] | 2005-12-25-05400 | 2 | 0 | 4 |
| i-NEC-SAT | 4 | [4] | 2005-12-25-05400 | 2 | 0 | 4 |
| f-NEC-MEC | 576 | [576] | 2006-01-05-84600 | 0 | 6 | 57 |
| f-NEC-MEC-s0101 | 576 | [336, 240] | 2006-01-05-84600 | 0 | 6 | 57 |
| f-NEC-MEC-s0103 | 576 | [432, 144] | 2006-01-05-84600 | 0 | 6 | 57 |
| f-NEC-SAT | 576 | [576] | 2006-01-05-84600 | 0 | 6 | 57 |
| i-NEC-MEC-noon | 4 | [4] | 2005-12-25-05400 | 2 | 0 | 4 |

### Subset behaviour (not a pass of the whole combination)

- `i-NEC-MEC`: FAIL — log.colm: no "CoLM Execution Completed"
- `i-NEC-MEC-s0101`: FAIL — expected 2 colm segment logs, found 1
- `i-NEC-MEC-s0103`: FAIL — expected 2 colm segment logs, found 1
- `i-NEC-SAT`: FAIL — log.colm: no "CoLM Execution Completed"
- `f-NEC-MEC`: model-completion checks PASS locally (steps=576, final=2006-01-05-84600, diagnostics=0).
- `f-NEC-MEC-s0101`: model-completion checks PASS locally (steps=576, final=2006-01-05-84600, diagnostics=0).
- `f-NEC-MEC-s0103`: model-completion checks PASS locally (steps=576, final=2006-01-05-84600, diagnostics=0).
- `f-NEC-SAT`: model-completion checks PASS locally (steps=576, final=2006-01-05-84600, diagnostics=0).
- `i-NEC-MEC-noon`: FAIL — expected 2 colm segment logs, found 1

The f-* (LULC_IGBP_PFT/NEC) and i-AMZ (LULC_IGBP/AMZ) cases do not show the step-4
residual in the portion of the 12-day window they reached; that is subset evidence only.
The main acceptance verdict remains FAIL because the required IGBP/NEC combined cases abort.

### Validator reporting corrections (read-only, no physics rerun)

The main job's `report.json` is the executed validator's output and is preserved as-is. Two
*reporting* defects were found afterwards and corrected in `combined10/validate.py`; neither
changes the FAIL verdict (the IGBP/NEC abort is in the model logs):

1. `iso_check.py` was invoked on `<run>/` instead of the model case directory `<run>/case/`, so the
   G7 rows for `f-NEC-MEC` and `i-AMZ-MEC-cl2e5` show `rc=1` spuriously. Re-run correctly inside
   `codex-81a0f423c967451abfcb5adb5ef28881`: `f-NEC-MEC`, `f-NEC-SAT` and `i-AMZ-MEC-cl2e5` all
   return rc=0 with 50 isotope variables, ~18,219,880 / 18,219,870 / 18,217,609 valid values,
   0 unmasked non-finite values and 19–21 d-excess variable pairs (broad diagnostics only).
   `i-NEC-MEC` still fails because the aborted run produced no isotope history at all.
2. The delivery tool `tests/material_transport/check_forcing_cache.py` assumes an isotope-only run:
   it requires `mt_forcing_precip_last/vapor_last > 0` for *every* column of
   `mt_forcing_species`. A combined run registers CL and SEDIMENT as tracers, so
   `mt_forcing_species = ntracers = 4`, while the CL column is legitimately 0 (standard chloride
   namelist: init/precip/vapour default 0) or 3e-6 (nonzero-CL extra) and SEDIMENT is 0. The
   tool therefore rejects otherwise valid schema-2 caches. The same read-only job audited the
   fields directly on the completed subset: `f-NEC-MEC`, `-s0101`, `-s0103`, `-SAT` each have 12
   initialized checkpoints with schema=2, forcing_variables=6, identity shape (7, 1544),
   reliability flags in {0,1}, both isotope columns strictly positive and all non-isotope columns
   exactly 0; `i-AMZ-MEC-cl2e5` (nonzero CL) matches with CL precip column = 3e-6 as configured.
   The only audit flag left is the declared cold-start checkpoint `2005-359-00000` (no cache by
   design, allowed exception). Consequence: gate G4 cannot be certified by the required delivery
   tool for a combined isotope+CL+SEDIMENT cache without extending the tool; the caches
   themselves are structurally consistent on the completed subset, and the IGBP/NEC cases have no
   initialized checkpoints at all because they abort at step 4.

### Positive subset facts (still not a pass of the full combination)

From the same preserved report/logs, the four LULC_IGBP_PFT/NEC cases are internally clean over
the full 12-day window: 576 steps each (split segments 336+240 / 432+144), final timestamp
`2006-01-05-84600`, model END, zero SIG/BAL/LULCC_BAL/residual-hard/STOP/abort messages, zero
unmasked NaN/Inf, real LULCC MEC and SAT transitions at the 2005→2006 boundary, registered
transport species `H2_18O, HDO, CL`, `f_frcsat` valid for all 1,200 element-days per case with no
masked value, 12 window-boundary sidecars per case with `hacc_nac=0` (window reset before the
restart), and **continuous runs bit-identical to both split runs** (restart and history, all
files/variables: 4 strict comparisons in the report). The `i-AMZ-MEC-cl2e5` extra also reaches
576 steps / END with the same clean diagnostics (its only G4 issue is the tool limitation above).
Subset results do not turn the overall verdict into a pass.

## Trigger localization (single-factor ablations)

Round 1 (`codex-17e7132400c44adeb691e66e3dc2c586`, 12 steps each, IGBP/NEC, IsoGSM forcing,
fractionation on, LULCC MEC build; `V*` rows add exactly the listed options on top of the
p9q-like forced-isotope setup):

| variant | added options | status |
|---|---|---|
| V0-control | — (p9q-like) | PASS (steps=12) |
| V1-cl | CL | PASS (steps=12) |
| V2-sed | sediment | PASS (steps=12) |
| V3-rivaropts | reservoir+levee+bifurcation | PASS (steps=12) |
| V4-fb | flood feedback only | DIAG_RESIDUAL (steps=4 ; resid=TRC_BAL residual hard failure: n_entries=       1 worst_abs= 0.90765E-04 tol= 0.50000E-04 @ipatch=     213 itrc=   1 ptype=  0 owner=   0) |
| V5-riveropts-fb | river options + flood feedback | DIAG_RESIDUAL (steps=4 ; resid=TRC_BAL residual hard failure: n_entries=       1 worst_abs= 0.90882E-04 tol= 0.50000E-04 @ipatch=     213 itrc=   1 ptype=  0 owner=   0) |
| V6-full-no-cl | river options + sediment + flood feedback | DIAG_RESIDUAL (steps=4 ; resid=TRC_BAL residual hard failure: n_entries=       1 worst_abs= 0.90882E-04 tol= 0.50000E-04 @ipatch=     213 itrc=   1 ptype=  0 owner=   0) |
| V7-full | river options + sediment + flood feedback + CL | DIAG_RESIDUAL (steps=4 ; resid=TRC_BAL residual hard failure: n_entries=       1 worst_abs= 0.90882E-04 tol= 0.50000E-04 @ipatch=     213 itrc=   1 ptype=  0 owner=   0) |
| V8-sed-cl | sediment + CL | PASS (steps=12) |

Round 2 (`codex-0c007acf2ceb4863aebe3d3dbf812cd2`):

| variant | configuration | status |
|---|---|---|
| W2-fb-nofrac | flood feedback only, **fractionation off**, IsoGSM forcing | DIAG_RESIDUAL 9.09e-5 > 5.0e-5 at step 4 -> fractionation not required |
| W3-fb-sat | flood feedback only, fractionation on, **SAT scheme** | DIAG_RESIDUAL 9.0765e-5 at step 4 -> LULCC scheme not required |
| W4-nofb-full-frac | full river options + sediment + CL **without flood feedback**, 5 days (240 steps) | **PASS**, zero diagnostics, ends 2005-12-29-84600 |
| W1-jra-full | JRA + flood feedback, but isotope-forcing groups accidentally left in the isotope namelist | CONFIG ERROR ("tracer_forcing_add_var: raw/total forcing temporal config mismatch"), redone as X1 in round 3 |

Round 3 (`codex-3c6c916b96a142f6a17f88cbfe25c374`):

| variant | configuration | status |
|---|---|---|
| X1-jra-static-full | JRA dynamics forcing, **no isotope forcing**, static isotopes, full river options incl. flood feedback + CL + sediment, 3 days | **PASS** 144 steps, no diagnostics |
| X2-fb-lulccoff | flood feedback + IsoGSM isotope forcing, **LULCC-off build** | DIAG_RESIDUAL 9.0882e-5 at step 4 -> LULCC compile switch not required |

Net minimal trigger: `DEF_GridRiverLake_FloodFeedback=.true.` plus runtime isotope forcing
(IsoGSM material forcing fields), on the LULC_IGBP/NEC build. It is independent of CL,
sediment, levee/bifurcation/reservoir, fractionation setting, LULCC scheme and the LULCC
compile switch; it is build/region specific (PFT/NEC and AMZ/IGBP do not show it so far).

Interpretation (from the source comments in `MOD_MT_SoilWater.F90` / `MOD_MT_Conservation.F90`):
with flood feedback the code deliberately disables the R_init-host-water correction for a
flooded soil column and books any unresolved dry-surface isotope pool as an explicit numerical
residual. With runtime-forced isotopes that booked residual reaches 9.09e-5 for patch 213,
above the code's own independent hard limit 5.0e-5, so the run stops. This is a real numerical
path in the latest source, not a harness artifact; it needs a root-cause fix or an explicit
user-level decision about the residual policy (relaxing `DEF_MT_RESID_ABORT_NBAD` would only
hide it and is NOT counted as a pass here).

## Trigger regression classification (pre-existing vs latest-commit regression)

Job `codex-69298c1322fd4c749b69e9999c7b74d8` reused the two already-built binaries from p9r
(`codex-4ba194bf50074f68aaa3d76178c49830/results/bld/{old,new}`; both are
`MT+GridRiverLakeFlow`, `LULC_IGBP`, `LULCC` off) and ran the **minimal trigger case** for 12 steps:
IsoGSM isotope forcing (fractionation on) + `DEF_GridRiverLake_FloodFeedback=.true.`, no CL, no
sediment, no levee/bifurcation/reservoir, one 30-minute step stream, 5 ranks, NEC.

| binary | source | colm.x sha256 | step-4 result |
|---|---|---|---|
| old | `28478b18` (branch base) | `eb30ab5ba3f0b81088181c1d7e4612915c668166ecf5d4ffe02fa3d2ed80e90e` | `TRC_BAL residual hard failure: worst_abs= 0.90765E-04 tol= 0.50000E-04 @ipatch=213 itrc=1 ptype=0` |
| new | `c315c849` (pre-frcsat-fix) | `c41e88539990eb90bcdd094f25043ed09736a767f467275e1d58ea8b48fcf3ea` | same line, same numbers, same patch/tracer |

Both aborted at `TIMESTEP = 4 / DATE = 2005-12-25-05400`; the log excerpts are
`evidence/iso5/OLD-28478b18.log.colm` and `evidence/iso5/NEW-c315c849.log.colm`.

**Conclusion: pre-existing at the branch base, not a regression from the latest commits.** No bisect
is needed for this finding (and none was attempted). The defect is a never-covered combination:
flood-feedback material coupling together with runtime-forced isotopes. It is outside the scope of
this acceptance task to fix; the follow-up decision is either a root-cause fix of the booked
`trc_numerical_residual_step` path or an explicit, documented residual-policy decision.

## frcsat fix isolation: same case on c315c849 vs 6d995438

Job `codex-05c172686e434105884e3a0c7297f6fb` (independent builds of both commits, same case:
NEC / LULC_IGBP / MT on / no LULCC / no isotope forcing / JRA forcing / CL / river routing;
120 steps ending 2005-12-27 12:00 so a restart is written in the middle of a daily history window):

| check | result |
|---|---|
| both runs | 120 steps, final `2005-12-27-41400`, `CoLM Execution Completed`, zero diagnostics |
| history | **identical raw bits for every file/variable including `f_frcsat`** |
| restart | 13 files, **767 common variables compared, all bit-identical**; the fixed source's only difference is the intended extra mid-window sidecar fields |
| new sidecar fields | `hacc_nac=24`, `nac_frcsat=24` on all 587 patches, `hacc_frcsat_count_version=1` |

This closes the "frcsat fix does not change the numbers" question for a real NetCDF/MPI case
and is the only part of this work that can be cited as a positive frcsat number-invariance result.

## Limits / not covered

- The full combined configuration is a **hard FAIL** (see above); do not cite the frcsat
  isolation or the PFT/AMZ subset as acceptance of the full combination.
- Reproduction coverage of the trigger is only for LULC_IGBP/NEC and, through the pre-built p9r
  binaries, for the branch base 28478b18 and c315c849 (both fail). PFT/NEC and IGBP/AMZ did not
  show the abort in the main job's 12-day subset runs; no other region, LULC_IGBP_PC build,
  other scheme, or other forcing dataset was tested for this trigger.
- The root cause inside the booked residual path was not fixed and no bisect beyond the existing
  p9r old/new binaries was performed (the defect predates both).
- PFT/NEC and IGBP/AMZ reached the full window cleanly in the main job, but only NEC/IGBP (and the
  p9r old/new binaries) were tested for the step-4 trigger; other regions, LULC_IGBP_PC, other
  forcing datasets and any additional scheme combinations were not re-run.
- Gate G4 cannot be certified with the delivered `check_forcing_cache.py` for a combined
  isotope+CL+SEDIMENT restart because that tool encodes the isotope-only positivity invariant;
  the field-level audit in this directory is supplementary evidence, not a replacement tool.
- The dry-lake invalid-sample branch of the frcsat protocol was not observed in a completed
  run: the extra `DEF_USE_Dynamic_Lake=.true.` case aborted at the same flood-feedback step 4,
  and the standard combined runs keep `DEF_USE_Dynamic_Lake=.false.` (so valid counts equal
  step counts by construction). The protocol-sidecar save/restore semantics are demonstrated
  by the mid-window sidecar and the old/new isolation; mask semantics themselves rely on the
  frozen-source `check_frcsat_history.py` self check (PASS, 14 scenarios x -O0/-O2, 10
  corruption classes) and are not claimed from the failed combined run.
- No hosted-CI compiler matrix and no methane work is part of this directory.
- `compare.py` compares variable bits/dimensions/shapes/dtypes; file encoding and attributes
  are outside the contract.
- Delta ranges are diagnostics only; no observational accuracy is claimed anywhere.

## Evidence map

Local (this directory):

- `run.sh`, `lib.sh`, `validate.py`, `new/`, `source.sha256`, `source_versions.txt`, `regions.txt` — the executed acceptance harness and frozen source.
- `iso/` — fcrsat isolation harness (`iso/run.sh`, `iso/validate.py`, `iso/old` = c315c849 archive).
- `iso2/`, `iso3/`, `iso4/`, `iso5/` — ablation/classification harnesses (locate the flood-feedback trigger; iso5 runs the minimal trigger on the pre-built p9r old/new binaries).
- `iso6/` — read-only supplementary audit harness (corrected iso_check path, field-level forcing-cache audit, frcsat boundary counts).
- `results.json` — machine summary of the verdict, jobs, gates, localization and limits.
- `evidence/main/` — main job status/report/validation log, build define.h/logs, failed-case log excerpts.
- `evidence/iso*/` — small reports of the isolation and ablation jobs; `evidence/attempts.json`.

Remote (kept on 7920, read-only for later audits):

- `codex-2434412927a94e38a44ad21db5643a38/results` — main job full outputs (failed verdict).
- `codex-17e7132400c44adeb691e66e3dc2c586/results`, `codex-0c007acf2ceb4863aebe3d3dbf812cd2/results`, `codex-3c6c916b96a142f6a17f88cbfe25c374/results` — ablations.
- `codex-05c172686e434105884e3a0c7297f6fb/results` — frcsat isolation.
- `codex-69298c1322fd4c749b69e9999c7b74d8/results` — trigger regression classification on the pre-built p9r binaries.
- `codex-81a0f423c967451abfcb5adb5ef28881/results` — read-only supplementary audit.
- Builds reused read-only by the probes: `codex-2434412927.../results/bld/i-on` (LULCC on), `codex-05c17268.../results/bld/new6` (LULCC off).
