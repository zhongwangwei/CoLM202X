# Signed soil-tracer debt: isolated experiment

**Off by default; not a production fix.** `TRACER_SIGNED_SOIL_DEBT_EXPERIMENT`
enables a direct `tracer_soil_water` kernel experiment only. No hydrology,
non-TRACER source, public interface, namelist or restart schema is changed.

## What is being tested

For closed, liquid-only advection of inert, generic-water-owned isotopes,
retain signed tracer inventory when the carrier water is slightly negative,
and use equal-and-opposite inter-layer transfers. The water deficit is bounded
by `0.001 / (2 * nl_soil)` mm, matching the current solver tolerance.

**This trades exact conservation for bounded concentrations.** The debt is
water the solver over-withdrew; it never existed. The old rule kept the error
in the concentration, as an unbounded, unrecorded extrapolation. This rule
moves it to the mass ledger instead: repaying debt water of ratio `c_debt`
with inflow of ratio `c_in` books

    dM = min(inflow, debt) * (c_debt - c_in)

in `trc_numerical_residual_step`. This is a real numerical source or sink
even in exact arithmetic, not rounding. The ledger is a diagnostic, and the
mass is never repaid.

Three rules keep every mixture a nonnegative combination of real water:

- **Debt is repaid first, at its own composition.** Inflow into a layer in
  debt first cancels that debt; only the rest mixes with the layer or passes
  on. Inflow with a different composition leaves the mismatch `dM` above.
  It is explicit and at most the repaid water times the concentration range.
  Without this rule the debt acts as a negative mixing weight and
  extrapolates the concentration.
- **Bounds are checked on mass, not on mass / water.** Each layer must hold
  between `lowest ratio * water` and `highest ratio * water`. The rounding
  allowance scales with the tracer that passed through the layer. The rounding
  is then booked so the next snapshot is sign-compatible and in range.
- **An unresolved (empty) layer borrows its receiver's ratio.** Its outflow
  does not change the receiver's concentration, and a refill from the same
  receiver repays the debt exactly. This is a numerical closure assumption,
  not the unique physical answer: an empty layer has no composition of its own.

Species 1 uses a unit ratio (exact arithmetic); species 2 uses an O18-like
ratio, so rounding near a tiny carrier is exercised. Uniform cases drain, hold
and rewet, including repayment to exactly zero and to 1e-10 mm. Nonuniform
cases check throughflow mixing in both directions, repayment by richer and
leaner water, partial repayment, a debt layer in throughflow, convergent
repayment, and an empty layer that drains and is refilled. A deterministic
400-step random walk with a nonuniform profile checks, after every step,
conservation including the booked residual and the initial ratio range.
Every call, in every mode, must book at most the water it repaid that call
times the concentration range, plus rounding. A call without repayment may
book rounding only, and bookings of opposite sign cannot cancel across calls.
These are kernel tests, not evidence of a measured error or correction in a
real-site simulation.

## Deliberate limits

- No blanket removal of existing nonnegative clamps or validation.
- Snow, ice, evaporation, roots, boundary exchange, diffusion, fractionation,
  reactions, provider-owned species and solutes are not supported while debt
  is present. Production keeps nonvolatile solutes on the clamped, unmixed
  sweep, and signed transport would change that, so solutes stop explicitly.
- Production-coupling optional arguments are rejected during debt, even if
  their values are zero. The real full-model caller supplies these arguments;
  enabling the macro does **not** make production runs support signed debt.
- Incoming signed inventory is rejected by `tracer_save_storage` before a full
  timestep, and by restart validation before checkpoint acceptance. History and
  restart persistence of signed debt are not implemented.
- Exact signed conservation with no residual is given up when repayment
  water differs in composition (see `dM` above).
- Ratios large enough that a product with any layer water or face volume
  could overflow are rejected up front.

Full-model support needs separate lifecycle and process design; do not enable
this macro for scientific production runs.

## Reproduce the targeted checks

Start with a successful GNU Fortran TRACER build of the same revision and
configuration. This compiles three override objects into a private directory;
it does not modify the existing build. `BUILD` must contain `include/define.h`,
`.bld/*.mod` and `lib/libcolm.a`. Compiler, MPI and NetCDF must match that build.
Run from the repository root:

```bash
export BUILD=/absolute/path/to/verified-tracer-build
export OUT=$(mktemp -d)
mkdir -p "$OUT/modules"
export FC='mpif90 -fopenmp -ffpe-trap=invalid,zero,overflow'
export NETCDF_LIBS="$(nf-config --flibs)"
export EXTRA_LIBS='-llapack -lblas'
for name in MOD_Tracer_Conservation MOD_Tracer_Rest MOD_Tracer_SoilWater; do
  $FC -c -cpp -DTRACER_SIGNED_SOIL_DEBT_EXPERIMENT \
    -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 \
    -ffree-line-length-none -fcheck=all \
    -I"$BUILD/include" -I"$OUT/modules" -I"$BUILD/.bld" \
    -J"$OUT/modules" "main/TRACER/$name.F90" -o "$OUT/modules/$name.o"
done
ar rcs "$OUT/modules/patched.a" "$OUT/modules/"*.o
bash tests/check_tracer_signed_soil_debt.sh \
  "$BUILD" "$OUT/results" "$OUT/modules/patched.a"
```

The 31 modes include 19 accepted cases and 12 explicit safety rejections.
A generic crash is not accepted as a safety rejection. The script checks the
specific diagnostic. The full-timestep guard is exercised; restart rejection
is compile-reviewed, not exercised through a checkpoint round trip.

As a negative control, the `2d73bd4a` kernel fails 13 modes, for three
reasons. `repayexact`, `repaynear` and `drainlarge` are uniform and should
pass under the old rule too; they fail on rounding over a tiny carrier, a
defect of the old rule. `enriched`, `different`, `partialsettle`,
`debtthrough`, `convergent`, `emptybackflow` and `randomwalk` reflect the rule
change: the old rule stops, or it conserves exactly at an extrapolated
concentration (`convergent` gives 1.25 instead of 4/3). `solute`,
`solubility` and `hugeratio` reflect the narrower scope and the new guard.

With the macro **off**, the preprocessed source of all three modules (blank
lines excluded) is identical to `ed94372f`, so the default build is unchanged.
Non-TRACER source files are unchanged.
