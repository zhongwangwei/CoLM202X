# Three-day integration validation (2026-09-23)

This is a bounded functional and conservation check of the development branch,
not a calibrated or fully validated model release. Runs used real IsoGSM forcing,
real canopy-structure products, a 1,800 s land step, MPI=3, and OMP=1. No
canopy morphology, snowfall, carbon input, or groundwater mixing parameter was
fabricated to make a case pass.

| Build / case | Three-day result | Evidence checked |
| --- | --- | --- |
| Plain IGBP, scheme 1/8 × `DEF_VEG_SNOW` on/off | 4/4 runs completed 144/144 steps; four required broadleaf patches had valid structure | All 257 initial files matched by hash across controls; scheme-8 initial capacity was 1.36–1.43× scheme 1; after including land `wdsrf`, each coupled three-day water residual was below 0.002 m³; 80/80 NetCDF files passed `ncks --chk_nan`. |
| PFT-IGBP, scheme 8 | Full BIF/levee/reservoir/flood-feedback configuration, all-off control, and `DEF_SPLIT_SOILSNOW`-on control each completed 144/144 steps; 32/32 required tree PFT records valid | Full-on continuous versus 1+2-day restart and history states were bitwise equal; finite-solubility, dry-pool, nonnegative-inventory and independent river-water checks passed. |
| PC, scheme 8 | Ordinary and extended land runs, plus extended `GridRiverLakeFlow`, each completed 144/144 steps; 32/32 required broadleaf PFT records valid | CH4 local/global balance residuals were at most about 7e-27/2e-25; the route-enabled days-2–3 coupled water residual including land `wdsrf` was 0.006805 m³; H2_18O, HDO and FINITE had positive river mass. A same-initial scheme-1 control gave 21.66% less three-day interception. |
| USGS, scheme 8 | Warm-domain ordinary and extended runs completed; a real-snow Adirondack scheme-8 run and same-initial scheme-1 control each completed 144/144 steps | Nine required forest patches had valid structure. Snow canopy and ground stores were nonzero and nonnegative; `ldew = ldew_rain + ldew_snow` held within 1.11e-16 mm; independent river-water closure had maximum relative residual 5.58e-15; 1+2-day restart ended with all 308 numeric variables identical to the continuous run. |

All completed cases exited normally. Scanned output had no NaN/Inf or negative
checked water/tracer inventory. Scheme 8 was not inferred from the namelist
alone: the input structure gate and same-initial scheme-1 comparisons showed
its capacity/process branch was active. The warm USGS comparison had 11.85%
more gross interception with scheme 8 than with scheme 1.

The repository-wide Python suite was rerun before publication:
`1067 passed, 27 skipped, 35 subtests passed` (`python -m pytest -q`). Skips
remain environment/fixture gates; this count is not substituted for the
real-data, three-day runs above.

## Accounting boundaries and remaining tests

- Land-patch `wdsrf` and routing `volwater_ucat` are different stores. Omitting
  land `wdsrf` initially produced a false 15,134 m³ residual in the plain-IGBP
  audit; including it reduced all four residuals below 0.002 m³. A cold
  `mkinidata` routing restart encodes initial stage with a zero volume
  placeholder; the first routing call materializes the effective volume from
  the stage curve. Raw day-zero `volwater_ucat=0` is not a physical inventory.
- Adirondack genuinely exercised snow interception, but selected histories did
  not expose phase-resolved snow unloading, throughfall and sublimation water
  fluxes. A term-by-term canopy-snow mass identity remains unverified.
- Enabling BIF, levees and reservoirs is not proof that overtopping or every
  dam-construction/operation threshold occurred during these three days.
  Complete whole-domain tracer and sediment budgets cannot be reconstructed
  from the selected history outputs. CL had no source and stayed at zero.
- Plain IGBP and the completed PC/PFT domains exercised broadleaf formulae.
  A real PNW PC surface product had 54/54 valid required tree PFT structures,
  including needleleaf, but its isotope cold start correctly stopped because
  real `wa` reached -7,459.743 mm while the **uncalibrated QA** effective
  aquifer-mixing water was 1,000 mm. No needleleaf three-day claim is made.
- Full soil/wetland/rice CH4 was not run in PFT/USGS cases lacking suitable
  nonzero BGC carbon/NPP inputs. The completed PC CH4 restart differed only
  in six inactive component-lag fields; aggregate state and subsequent daily
  CH4 histories matched. The 1,000 mm aquifer-mixing setting used elsewhere
  is a QA value, not a physical calibration.

Detailed logs, NetCDF outputs, input/binary hashes, and machine-readable
audits were retained locally under `/tmp/colm_*_20260923*`; they are not part
of this Git branch because forcing and restart products are large and
machine-specific. The report records their conclusions and limitations; it
does not substitute for longer 30/90-day, active overtopping/dam-switch, or
calibration runs.
