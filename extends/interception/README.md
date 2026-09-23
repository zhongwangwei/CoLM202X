# Extended canopy interception schemes

This directory keeps the extended canopy interception parameterizations
(CLM4/CLM5/Noah-MP/MATSIRO/VIC/JULES/CoLM2024) outside the `main/` physics
path.

The current repository default defines `extend_interception` in
`include/define.h`.  The Makefile therefore replaces the four same-name
`main/` modules with:

- `MOD_LeafInterception_Extended.F90`
- `MOD_LeafTemperature_Extended.F90`
- `MOD_LeafTemperaturePC_Extended.F90`
- `MOD_Thermal_CanopyPhase_Extended.F90`

Set `DEF_Interception_scheme = 6` in the runtime namelist to select VIC.  VIC
uses its native vegetation-tile convention (`F=1`); `fsno` is passed separately
for snow-process activation.  The repository default is now
`DEF_VEG_SNOW = .true.`; set it explicitly to `.false.` only for experiments
that need the old no-vegetation-snow behavior.  Change the switch back to
`#undef extend_interception` and run a clean rebuild to restore the `main/`
CoLM2014 path.

`DEF_Interception_scheme = 8` (CoLM2024) uses a morphology- and wind-dependent
liquid canopy-water capacity.  In `LULC_USGS`, the unambiguous classes are
mapped as shrubland (8), broadleaf forest (11/13), needleleaf forest (12/14),
and mixed forest (15).  Mosaic, savanna, wooded-wetland/tundra, and other
classes retain the explicit LAI/SAI capacity because USGS does not provide the
woody fraction or leaf-type composition needed to choose a morphology formula.
USGS crown inputs are aggregated from the same `canopy_data` product used by
IGBP.  Initialization must select scheme 1 for the whole run, with a clear
warning, when the complete input group is unavailable; it must never mix
invented per-patch structure into a scheme-8 run.

USGS broadleaf and mixed-forest capacity currently uses the classification
mean `htop0` shared by every interception scheme.  It is not an observed crown
height.  Keeping the common height avoids changing aerodynamic roughness and
displacement only in the scheme-8 experiment; measured-height support requires
a separate all-scheme design and validation.
