# Material host interface cleanup

## Scope and invariants

Reduce material-specific declarations and parameter forwarding in CoLMMAIN,
WATER_VSF and soil_water_vertical_movement. Keep water calculations, diagnostic
sampling points, material update order, restart/history fields and configuration
unchanged. Snow-layer remapping and unrelated host fixes are outside this pass.

## Sequence

1. Review this plan independently and establish existing focused numerical
   regression results before editing production code.
2. Group water-path diagnostics in a local derived type shared by the water
   routines and the material consumer. Keep allocation/reset explicit and do
   not introduce shared mutable state. Update all callers and Makefile ordering.
3. Group CoLMMAIN's material-only snapshots and scratch values in a local work
   type under MATERIAL. Move self-contained material bookkeeping out of the
   host where that reduces its responsibilities without changing event order.
4. Run focused numerical checks, a diagnostic reset/bounds check, and independent
   review. Compile full MT-on/off configurations on T7920 using private snapshots.

## Validation

- Existing material physics regression checks pass before and after edits.
- Compare pre/post water-solver host outputs and path diagnostics numerically
  using identical synthetic column inputs, including impermeable layers and
  root/aquifer transfers, where a focused harness can cover the real routines.
- Check reset on repeated calls, changed column sizes, snow lower bounds and
  independent diagnostic objects with compiler runtime bounds checks.
- Confirm no diagnostic references remain as long argument lists at the two
  host interfaces, and compile MT off as well as MT on.
- Report any missing coupled land-model/restart coverage explicitly; existing
  historical test reports alone are not evidence for this refactor.

## Risks

Lower bounds for snow arrays and lifetime of patch-local allocatables must be
preserved. Reset must not discard thermal/snow diagnostics before their consumer.
Passing a diagnostic object must not make it an input to host water calculations.

## Follow-up authorized after the initial cleanup

- Fix the naming inventory's false positive with an exception restricted to
  the original BGC file and its two upstream concentration identifiers. Prove
  the exception does not permit those names elsewhere or other legacy names.
- Reuse declared, read-only T7920 NEC IsoGSM inputs for same-compiler baseline
  and edited PFT+BGC+CROP builds. Run continuous two-day and one-day/restart/
  one-day cases for each revision. Check completion and step counts.
- Require exact nonempty file/variable/time coverage, raw values and semantic
  NetCDF attributes; fail on missing files or differences. Compare both revisions'
  histories/restart outputs and each revision's continuous/split final state.
