# Material transport checks

These focused checks do not replace the MPI and data-driven regression matrix.
Current integration status and outstanding work are in
[the status record](../../docs/material-transport-status.md).

## Run the self checks

The forcing check needs Python 3 and `gfortran`. It extracts the production
Fortran routines and compiles them with small host stubs; no model data or MPI
is needed. Fractionation coefficients are stubbed, so this checks control flow
and cache remapping, not fractionation physics.

```sh
python3 tests/material_transport/check_forcing_edges.py
```

The NetCDF tools also need NumPy and netCDF4, already present in the T7920 test
environment. Their self checks create only temporary, tiny fixtures.

```sh
python3 tests/material_transport/compare.py --self-test
python3 tests/material_transport/iso_check.py --self-test
```

## Check CI configurations

This standard-library check needs a C preprocessor (`cc`). It generates headers
in temporary directories, checks effective macros and unique case names, and
never changes the working `include/define.h`.

```sh
python3 tests/material_transport/check_defineh.py
```

The generator accepts an optional ninth argument, `LULCCON` or `LULCCOFF`
(default off); the optional eighth material-transport argument is unchanged.
BGC is supported with both PFT and PC, as in the production header.

## Check simulation outputs

```sh
python3 tests/material_transport/compare.py RESTART_A RESTART_B
python3 tests/material_transport/compare.py HISTORY_A HISTORY_B
python3 tests/material_transport/iso_check.py CASE_DIRECTORY
```

The comparator requires identical NetCDF file and variable sets, dimensions,
shapes, types and raw variable bits. It returns nonzero on a mismatch. It does
not compare file encodings or attributes. Compare equivalent output layouts;
intentional schema additions or segmented river-history filenames need separate
accounting, not an ignored missing-file warning.

The isotope report distinguishes declared fill values from unmasked NaN/Inf,
and computes d-excess using jointly valid positions in the same file. Unmasked
non-finite values and missing isotope histories fail the check. Its broad delta
ranges are diagnostics, not observation-based acceptance criteria.

## Additional regression checks

```sh
python3 tests/material_transport/check_singlebucket_precip.py
python3 tests/material_transport/check_site_history.py CASE/history --days 365
```

The canopy test compiles extracted production routines and needs `gfortran`.
The site-history validator checks complete daily timelines, masks and exact
numerical equality against an optional reference; it is not a raw-bit comparator.
Use `compare.py` separately when both NetCDF layouts should match bit for bit.
The old host same-month site-writeback restart limitation is recorded in the
status document and is not hidden by these checks.
