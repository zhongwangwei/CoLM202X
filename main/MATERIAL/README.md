# Material Transport

Material transport carries dissolved and suspended materials with the water
that CoLM already moves: water isotopes, solutes, river sediment and methane.
It is compiled in with `#define MaterialTransport` in `include/define.h`
(it needs `vanGenuchten_Mualem_SOIL_MODEL`).

Compiled in but with nothing switched on, the model gives bit-identical
results to a build without it. With materials on, host water, energy and
carbon are not changed by the materials; the only host difference is the
nitrogen-deposition / fire start-year fix in `MOD_MT_HostFixes`
(`mt_bgc_data_init_year`), which acts whenever any material is on.

## Switches

The main namelist file takes one extra group; leaving it out keeps every
default and carries no material.

```fortran
&nl_colm_material
   DEF_USE_ISOTOPE         = .true.
   DEF_ISOTOPE_SPECIES     = 'H2_18O,HDO'
   DEF_ISOTOPE_PARAM_FILE  = 'run/standard_isotope_parameter.nml'
   DEF_SOLUTE_NAMES        = 'CL'
   DEF_SOLUTE_PARAM_FILES  = 'CL:run/standard_chloride_parameter.nml'
   DEF_SEDIMENT_PARAM_FILE = 'run/standard_sediment_parameter.nml'
   DEF_USE_METHANE         = .true.
   DEF_METHANE_PARAM_FILE  = 'run/standard_methane_parameter.nml'
/
```

- Isotopes and solutes: land water everywhere, plus the grid river-lake
  network when `GridRiverLakeFlow` is compiled in. Isotopes need a positive
  `DEF_MT_AQUIFER_MIXING_WATER_MM` (the effective mixed aquifer water depth,
  mm) because vanGenuchten turns on variably saturated flow; there is no
  default, set a measured or calibrated value. With fractionation on,
  liquid evaporation follows Craig-Gordon and transpiration the
  non-steady-state leaf model; sublimation of snow and ice does not
  fractionate (the vapour leaves at the ratio of the ice skin,
  `DEF_ISO_SUBL_SKIN_MM`), and the leaf water keeps its ratio when LAI
  changes its size. Settings are in `run/standard_isotope_parameter.nml`.
- Sediment: `DEF_USE_SEDIMENT = .true.` in `&nl_colm`, needs
  `GridRiverLakeFlow`; parameters in `&nl_colm_sediment_parameter`.
- Methane: needs `BGC` with `LULC_IGBP_PFT` or `LULC_IGBP_PC`; all settings
  in `&nl_colm_methane_parameter` of `DEF_METHANE_PARAM_FILE`.

Framework settings (`DEF_MT_*`: balance abort thresholds, aquifer mixing
depth, soil initial state file, forcing thresholds) are in the same group;
see `MOD_MT_Namelist.F90`.

- `DEF_MT_FORCING_MIN_PRECIP` (default 5e-6 kg m-2 s-1): total
  precipitation rate below which a heavy/total precipitation forcing ratio
  is not used, because IsoGSM ratios are noise at smaller rates. For an
  isotope with both a precipitation and a vapour forcing, rain then takes
  the ratio in equilibrium with the forcing vapour at the air temperature
  and snow the ratio of the last reliable record; other materials keep the
  last valid ratio.
- `DEF_MT_FORCING_MEMORY_PRECIP` (default 2e-5 kg m-2 s-1): only records
  above this rate are remembered as that snow ratio. Before the first
  one (cold start, or a new patch after land cover change) the snow ratio
  is the ratio of ice deposited from the forcing vapour. Land cover changes
  mix only reliable source ratios, weighted by their contributing area;
  a patch with no reliable source uses the deposition-equilibrium estimate.

Every parameter file entry is optional, the values in
`run/standard_*_parameter.nml` are the code defaults unless stated there.
A full case is in
`run/examples/Amazon_Grid_25km_IGBP_VG_RiverLake_Material.nml`.

## Layout

| Directory | Contents |
| --- | --- |
| this directory | switches (`MOD_MT_Namelist`), material table (`MOD_MT_Defs`), provider lifecycle (`MOD_MT_Lifecycle`, `MT_RegisterProviders`), land driver (`MOD_MT_LandPhase`), state, restart, history, LULCC (`MOD_MT_Lulcc`), river-lake glue (`MOD_MT_RiverLakeCoupling`), host fixes |
| `GENERIC_WATER/` | code shared by every material that moves with water: precipitation, snow, soil water, evaporation, special patches, forcing, conservation, river-lake transport |
| `FAMILY_ISOTOPE/` | fractionation and the two isotope species |
| `FAMILY_SOLUTE/` | solute-only rules (dissolution limit, first-order decay) |
| `FAMILY_PARTICLE/` | river sediment (`MOD_MT_Sediment`) |
| `FAMILY_GAS/` | methane (`MOD_MT_Methane_*`) |

Each material has a family (ISOTOPE, SOLUTE, PARTICLE, GAS), a state owner
(the generic water code, or a provider such as methane or sediment) and a
reaction mode (none, first order, provider). Providers attach their
restart, history and LULCC callbacks through `MOD_MT_Lifecycle`, so the host
calls one set of `mt_*` routines whatever is switched on.

## Methane and the host BGC

Wetland decomposition for methane runs on methane-owned C/N shadow pools.
Each step, on a wetland patch, the host BGC column is saved, the shadow pools
are swapped in, the upstream decomposition cascade runs, methane keeps what
it needs, and the host column is restored. Host `decomp_hr`, `er` and `ar`
are not changed; methane writes its own CO2-only respiration
`f_co2_hr_methane`.

At a cold start, a wetland layer without carbon gets its shadow pools from
the organic matter of the soil data (`OM_density`, 580 gC/kg, C:N 15). The
CoLM202X `TRACER` code did this only with `DEF_USE_CN_INIT`, so its wetlands
started without carbon and emitted almost no methane by default; here it is
done at every cold start.

The surface data for lake methane production (`lake_soilc_patches`) and
spatial soil pH (`methane_ph_patches`) are made by mksrfdata
(`Aggregation_LakeSoilC`, `Aggregation_MethanePH`) only when the methane
settings ask for them.

## Output

- History: material fields go to `<case>_hist_material_<date>.nc` and follow
  the host history switches (`DEF_hist_vars`) of the water flux or state
  they belong to. Methane has its own `history_mode` ('core' by default).
- Restart: material state is written next to the host restart; a restart
  without material state starts the materials cold.
- Balance: per-step land-water balance checks; a cell over tolerance is
  reported, and the run stops when more than `DEF_MT_BALANCE_ABORT_NBAD`
  cells fail.

## Land cover change

With `LULCC`, `LulccDriver` calls `mt_lulcc_save` before the new year's
land surface is loaded and `mt_lulcc_remap` after the host variables are
recovered. With the mass and energy conservation scheme
(`DEF_LULCC_SCHEME = 2`) the transfer trace is checked against the old and
new patch areas, the material pools are remapped by area, and the remapped
land-water material mass must match the mass before the change.

The materials then follow the host, so that every material stays attached
to its water:

- The host does not carry `wetwat` across the change (it is reset to the
  cold-start value). The material in the wetland water is reset with it at
  the initial ratio.
- A new patch with no source patch (with `DEF_LULCC_SCHEME = 1`, no old
  patch of the same class in the element) keeps the host cold-start state,
  and its materials are cold-started from that water.
- Elsewhere the old host water is remapped with the same weights as the
  material pools and compared with the new host water, compartment by
  compartment (canopy, snow column, each soil layer and phase, surface
  water). Where they differ (MEC re-bins snow layers to the dominant source
  patch and drops the pre-layer snow of the others; with PFT or PC the
  patch canopy water is rebuilt from the new PFT fractions), the material
  is set on the host water at the ratio it had; a snow column that differs
  anywhere gets one column ratio. Where they agree nothing changes.

All three are counted as a mass entering or leaving at the change, so the
mass check still closes. With isotopes and `DEF_MT_AQUIFER_MIXING_WATER_MM`,
the aquifer reference carrier is removed from a patch that became a lake,
glacier or urban patch, and created on a new soil or wetland patch that
received none, at the patch's own aquifer ratio (`FAMILY_ISOTOPE/MOD_MT_Isotope_Lulcc`).
