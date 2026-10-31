# Material host cleanup: phase 2

Baseline is the complete, uncommitted phase-1 working tree, already verified by
the two-day coupled/restart comparison. Preserve that work; do not reset to HEAD.

## Scope

1. Move post-interception material bookkeeping into the existing land facade.
2. Move scratch allocation/reset/snapshot/release into a lightweight work module.
   It cannot depend on the facade or history: MT_Hist uses 1DAccFluxes, which
   already uses LandPhase, so importing MT_Hist into LandPhase is a real cycle.
3. Move recharge/infiltration conversion into the land facade; move flood
   evaporation ledger access into the existing RiverLakeCoupling module.
4. Move decay -> balance -> history -> scratch release into the existing
   material history module, importing only the lightweight work type.

Do not change water/thermal algorithms, event timing, material update order,
configuration, restart fields or snow-layer topology. Keep snapshots on both
sides of thermal and interception. Do not introduce shared mutable scratch.

## Verification and review

- Existing focused physics checks pass before editing (completed).
- Independently review this plan and retain separate implementation/review roles.
- Add numerical/callback-order/bounds regression coverage for the moved phase
  bookkeeping and workspace lifecycle; preserve existing boundary checks.
- Confirm MT-off CoLMMAIN is unchanged after preprocessing.
- Compile MT-on and MT-off with dependency ordering checked.
- Run the same declared T7920 NEC two-day continuous and split case for the new
  source; compare strictly against the completed phase-1 results and internally
  across continuous/restarted execution. Retain source hashes, exit markers,
  complete time/file/variable coverage and finite-value checks.
- Keep heavy builds/integration on run-7920 private snapshots; do not modify
  shared server code, old results, or input datasets.
