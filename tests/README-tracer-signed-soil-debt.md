# Signed soil-tracer debt: isolated experiment

**Off by default; not a production fix.** `TRACER_SIGNED_SOIL_DEBT_EXPERIMENT`
enables a direct `tracer_soil_water` kernel experiment only. No hydrology,
non-TRACER source, public interface, namelist or restart schema is changed.

## What is being tested

For closed, liquid-only advection of inert, generic-water-owned isotopes and
passive solutes, retain signed tracer inventory when the carrier water is
slightly negative. Recover concentration from signed inventory / signed water,
and use equal-and-opposite inter-layer transfers. The water deficit is bounded
by `0.001 / (2 * nl_soil)` mm, matching the current solver tolerance.

The synthetic uniform-concentration case drains, holds and rewets repeatedly.
Both species retain concentration 1 and conserve signed total inventory without
adding numerical source/sink entries. Nonuniform three-layer cases separately
check incoming-water mixing in both directions. These are kernel tests, not
evidence of a measured error or correction in a real-site simulation.

## Deliberate limits

- No blanket removal of existing nonnegative clamps or validation.
- Snow, ice, evaporation, roots, boundary exchange, diffusion, fractionation,
  reactions, provider-owned species and finite-solubility species are not
  supported while debt is present.
- Production-coupling optional arguments are rejected during debt, even if
  their values are zero. The real full-model caller supplies these arguments;
  enabling the macro does **not** make production runs support signed debt.
- Incoming signed inventory is rejected by `tracer_save_storage` before a full
  timestep, and by restart validation before checkpoint acceptance. History and
  restart persistence of signed debt are not implemented.
- Repayment with a different composition can leave incompatible signs or cause
  a concentration overshoot. Such cases stop explicitly rather than being
  claimed as a general solution.

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

The 21 modes include nine accepted cases and twelve explicit safety rejections.
A generic crash is not accepted as a safety rejection. The script checks the
specific diagnostic. The full-timestep guard is exercised; restart rejection
is compile-reviewed, not exercised through a checkpoint round trip.

Validation also compared the three modules with the macro **off**: preprocessed
source excluding blank lines matched the baseline, and a 24-row actual-module
runtime probe matched byte-for-byte. Non-TRACER source files were unchanged.
