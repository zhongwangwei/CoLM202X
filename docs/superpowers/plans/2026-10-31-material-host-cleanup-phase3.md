# Material host cleanup: phase 3

Baseline: the uncommitted phase-2 tree, with successful MT-on/off builds and
strict two-day continuous/restart comparison. Preserve all earlier changes.

## Scope and invariants

- Extract the two material-only soil-water exchange diagnostic helpers into
  existing MOD_MT_HydroDiagnostics; retain observations at the same solver sites.
- Extract snow material scratch and remapping operations from all four standard
  and SNICAR combine/divide routines into a Precision-only material helper.
  Keep water layer topology, thresholds, sampling sites, arithmetic order,
  lower bounds, optional argument behavior and solid/ice/liquid phases unchanged.
- Parent owns build wiring; independent owners handle snow, hydro and regression.
  A separate reviewer checks the final change and numerical coverage.
- No new dependencies, global scratch, configuration or restart changes.

## Verification

Before production edits, run actual solver regression and freeze current sources.
Test all snow routine variants, layer disappearance/merge/shift/split thresholds,
material gradients, and absent/partial/full optional arrays against that baseline.
Check soil-column outputs bitwise against the phase-2 interface. Check effective
MT-off snow/soil sources unchanged. Run existing focused material checks.
Build MT-on/off on T7920 private snapshots. Run the established 96-step continuous
and 48+48 restart case, strictly compare all outputs to completed phase-2 job
codex-46c411fb98924fdb8d9e031d1155c31c and across restart execution. Require
terminal success, finite values and complete time/file/variable coverage.
