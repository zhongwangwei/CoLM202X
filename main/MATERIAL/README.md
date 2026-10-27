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
  network when `GridRiverLakeFlow` is compiled in.
- Sediment: `DEF_USE_SEDIMENT = .true.` in `&nl_colm`, needs
  `GridRiverLakeFlow`; parameters in `&nl_colm_sediment_parameter`.
- Methane: needs `BGC` with `LULC_IGBP_PFT` or `LULC_IGBP_PC`; all settings
  in `&nl_colm_methane_parameter` of `DEF_METHANE_PARAM_FILE`.

Framework settings (`DEF_MT_*`: balance abort thresholds, aquifer mixing
depth, soil initial state file) are in the same group; see
`MOD_MT_Namelist.F90`. Every parameter file entry is optional, the values in
`run/standard_*_parameter.nml` are the code defaults unless stated there.

## Layout

| Directory | Contents |
| --- | --- |
| `MATERIAL/` | switches (`MOD_MT_Namelist`), material table (`MOD_MT_Defs`), provider lifecycle (`MOD_MT_Lifecycle`, `MT_RegisterProviders`), land driver (`MOD_MT_LandPhase`), state, restart, history, LULCC (`MOD_MT_Lulcc`), river-lake glue (`MOD_MT_RiverLakeCoupling`), host fixes |
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
