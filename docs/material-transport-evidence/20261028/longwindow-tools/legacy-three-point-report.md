# AMZ IsoGSM evaporation trace: evidence, not a parameter fix

## Result

Job `codex-4f7bfe583e234dfe8ac35115d15fee56` **succeeded**. Frozen numerical source `3a1c391a`, private logging instrumentation only; GRID/IGBP/MT/vanGenu, PHS default, NP8, fractionation on, IsoGSM 2004-01-01..08. Actual final step 384 at 2004-01-08-84600 and model END; TRC_SIG/TRC_BAL/residual hard messages all zero.

Instrumented vs original P9L history: **all 3 files and all 225 variables bit-identical**, including the three reported extremes. This directly establishes transparency for the history being diagnosed (not a claim about every restart field or attribute).

All three extremes reconstruct from actual gross-water/signed-mass bookkeeping and the finite-pool Craig–Gordon formula. The Jan 8 HDO FillValue is caused by **negative signed isotope transfer / positive gross water**, not a missing water denominator, small-water mask, or integration explosion. This is an allowed atmospheric heavy-isotope uptake in the implemented CG exchange formulation; it is not a positive outgoing-water isotope ratio. Numerical consistency does not establish independent physical calibration of the humidity/kinetic assumptions.

## Actual target mapping and daily ledgers

Indices supplied originally were zero-based. Exact history intersection areas, filters, patch centroids, and both eligible/excluded contributions are in `replay-results/patches.json`. Identifiers below are process-local at NP8, not global patch IDs.

| Target/date | Grid center lat,lon | Rank/patch/class | Contributing patches / eligible | Actual patch gross water mm | Signed O18 mass | Signed HDO mass |
|---|---|---|---|---|---|---|
| soil Jan 2 | -1.125,-59.125 | 3/17/2 | 1/1 | 0.105683392176150637 | 1.87166866787154967e-4 | 4.03545517257660261e-6 |
| wetland Jan 6 | -2.375,-58.375 | 6/22/11 | 4/1 | 0.169442996809659546 | 3.07849921055591835e-4 | 9.07468674438138834e-6 |
| wetland Jan 8 | -2.375,-58.125 | 6/28/11 | 6/1 | 0.100588829590334772 | 1.67730165096191541e-4 | -2.95883268333629879e-6 |

Each eligible patch strictly exceeds the real 0.1 mm daily threshold. Jan 2 and Jan 8 are very close to it. Other mapped patches have zero/insufficient wetland accumulations and are explicitly retained in the evidence, not silently dropped.

The grid-level logged numerator/denominator are **intersection-area-weighted sums**, not patch mm. Their common area scale cancels in the isotope ratio. No `f_fevpg`/`f_fevpa` proxy was used. Reconstructing all filtered patch contributions reproduces each valid grid delta within 6.51e-14 permil.

| Target | O18 delta permil | HDO delta / mask | Actual finite-pool events per isotope | Negative HDO events |
|---|---|---|---|---|
| soil Jan 2 | -116.788877341969 | -754.851159068499 | 6 | 1 |
| wetland Jan 6 | -93.9383960657421 | -656.163418857909 | 33 | 13 |
| wetland Jan 8 | -168.420597675629 | FillValue: negative ratio | 32 | 27 |

Jan 8 signed HDO diagnostic before masking corresponds algebraically to -1188.849010603698 permil, i.e. a negative isotope ratio. `mass_to_delta` correctly masks it rather than publishing a physically interpretable water composition. The O18 net transfer stays positive. No constituent has been clipped or changed.

## High-precision finite-pool replay

`replay.py` uses Python standard-library Decimal at 50 digits and Fortran ES26.17E3 round-trippable inputs. Each event records pool water/mass, water loss, pool/vapor ratios, CG humidity before/after cap, equilibrium/kinetic alpha, slope/offset, remaining water/ratio, signed isotope loss and remaining mass. Replay evaluates:

```
W1 = W0 - loss
x = ln(W1/W0); y = (slope-1)*x
R1 = R0*exp(y) + offset*x*(exp(y)-1)/y
Mloss = M0 - R1*W1
```

It preserves the model's residual-ratio bounds, with the analytic y=0 limit. No production solver changes or parameter retuning.

- 3,600 captured finite-pool events over the three spatial windows / eight days; 188 signed negative losses.
- Maximum model vs Decimal50 mass-loss error: **8.3211e-15** over all events.
- Selected-day maxima: soil O18/HDO 4.01e-18 / 3.46e-19; Jan6 wetland 7.11e-16 / 7.83e-17; Jan8 wetland 7.37e-16 / 7.78e-17.
- Every selected-day BOOK water and signed isotope mass is exactly accounted for by the logged finite-pool events (6, 33, 32 per isotope); no unrecorded deficit pathway is required for these targets.
- Reaccumulated BOOK vs actual patch daily arrays differ at most 1.265e-17 mm water and 4.61e-22 isotope mass (normal summation rounding).
- Decimal50 day-aggregate deltas differ by at most about 6.68e-9 permil from the reported binary64 diagnostics, far too small to explain the extrema; Jan8 HDO remains negative.
- `finite_pool_events.json` contains every step's exact inputs and replay residual; `finite_daily_replay.json` contains independently recomputed day totals.

Example Jan8 HDO uptake, rank6 patch28 day8 sec3600:
- pool water 1208.74199951231753; loss 0.00180643164672671954 mm;
- pool ratio 1.45541817173543684e-4, vapor ratio 1.39964866197246976e-4;
- RH entering CG 0.992055348878099097, capped to 0.989999999999999991;
- alpha_eq 1.08024069894586283, alpha_k 1.00530802637928263;
- pool mass before 0.175922507103005332, after 0.175923196021014944;
- signed outgoing mass -6.88918009611461102e-7; Decimal50 -6.8891800968597952761e-7 (7.46e-17 difference).

## Reproduce and scope

```
python3 .tmp/ct/iso_trace9/replay.py \
  .tmp/ct/iso_trace9/evidence \
  .tmp/ct/iso_trace9/replay-results
```

Exit 0 asserts complete target mapping, all six target/isotope ledgers, all finite-loss contributions, small replay residuals, correct reconstructed valid deltas, and exactly one negative-ratio mask. Fresh run passed.

Artifacts: `instrument.py`, `run.sh`, `replay.py`, `evidence/status.json`, `evidence/model-verification.txt`, three `trace_rank_*.txt`, and `replay-results/*.json`. Full instrumented build, source and native NetCDF remain under the job's results directory on 7920l.

Limits: RH logged as entering the CG function (already subject to the upstream surface-relative-humidity safeguard) and after the CG max-RH cap; raw uncapped meteorological eair/esat is not separately logged. This replay diagnoses actual CG inputs, not the upstream saturation formula. Neither broad physical plausibility nor independent validation of the CG humidity closure is claimed. No source outside this private task directory changed; no commit, push, merge, or edits to prior jobs.

### Separate long-window evidence is not resolved by this trace

The history-drift reviewer reports additional p9k AMZ Jan–Mar wetland extrema (O18 -368.91..214.58; HDO -997.51..1708.50) in `.tmp/ct/p9o/p9k-strict/amz-wet-isotope.txt`. This three-cell/eight-day replay does **not** validate or explain every extreme in that longer window. In particular it must not be generalized as “all wetland outliers are caused by near-zero host ET”; our demonstrated mechanism uses the actual gross denominator and signed uptake for these exact three targets only. Further long-window targets require their own gross/signed/input ledgers.
