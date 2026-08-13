# CONFIG_CARD — frozen methane baseline

    anchor: methane-lean @ 7e05e0c4 (2026-08-13)
    rule:   FROZEN. Values never change. New code parameters are appended
            at neutral values (switch=off, scale=1, increment=0), marked +.
            Removed parameters keep their row, annotated "removed in <commit>".
    config: effective single-point batch config (code defaults overridden by
            scripts/sites/ch4_parameter_baseline.nml in the Methane repo).
    marks:  * = explicit override of the code default (appendix A lists both).

## 1 Macros

| macro | state |
| --- | --- |
| SinglePoint, LULC_IGBP_PFT, vanGenuchten_Mualem_SOIL_MODEL | on |
| BGC, CROP, TRACER, LULC_IGBP_WFT | on |
| CatchLateralFlow, GridRiverLakeFlow (forced off by SinglePoint) | off |
| URBAN_MODEL, LULCC, RangeCheck, CoLMDEBUG, USEMPI (forced off) | off |

## 2 Classes

Wetland patches (patchtype 2) carry one WFT tile. Methane parameters are
class-blind at baseline: f_methane/redoxlag resolve by site latitude
(biome tree), substrate and inundation are identical across classes.

| class | idx | base PFT | template | substrate | inundation |
| --- | --- | --- | --- | --- | --- |
| permafrost peatland | 79 | 12 | moss | frozen SOM + root resp | saturated |
| wet tundra | 80 | 12 | graminoid | frozen SOM + root resp | saturated |
| bog | 81 | 12 | moss | frozen SOM + root resp | saturated |
| fen | 82 | 13 | graminoid | frozen SOM + root resp | saturated |
| marsh | 83 | 13 | graminoid | frozen SOM + root resp | saturated |
| salt marsh | 84 | 13 | graminoid | frozen SOM + root resp | saturated |
| tropical swamp | 85 | 4 | woody | frozen SOM + root resp | saturated |
| temperate swamp | 86 | 7 | woody | frozen SOM + root resp | saturated |

Non-wetland routing: rice -> CROP path (enable_rice_paddy); lake ->
sediment module (allowlakeprod); upland soil -> soil path (only_wetland=F).
Class source: SITE_wetland_class (sites) / global WFT map (gridded).

## 3 Vegetation templates (WFT columns in MOD_Const_PFT)

Overrides on the base-PFT copy; empty = pure copy. Source: scripts/wft/
wft_pft_templates.yaml (generated columns, checked by --check).

| param | moss (79,81) | graminoid (80,82-84) | woody (85,86) |
| --- | --- | --- | --- |
| slatop | 0.00781 | base | base |
| leafcn | 35.56 | base | base |
| lflitcn | 66.0 | base | base |
| leaf_long | 0.9744 | base | base |
| froot_leaf | 0.3944 | base | base |
| htop0_p / hbot0_p | 0.05 / 0.0 | base | base |
| sai0_p | 0.0 | base | base |
| vmax25_p | 4.3 | base | base |
| roota / rootb | base | base | 11.0 / 5.0 |
| woody flag | 0 | 0 | 1 |

## 4 Production

P = f_methane x (SOMHR + root resp) x Q10(T) x pH factor x O2 inhibition
x redox lag. Wetland SOM pools are frozen at their initial stock
(wetland_fixed_substrate): decomposition is bookkept, never debited.

| param | value | note |
| --- | --- | --- |
| f_methane | 0.2 | fallback fraction |
| use_biome_f_methane | T * | biome tree on |
| f_methane_tropical_peat | 0.12 | tree leaf |
| f_methane_tropical_floodplain | 0.07 | tree leaf |
| f_methane_floodplain | 0.10 | tree leaf |
| f_methane_temperate_marsh | 0.15 | tree leaf |
| f_methane_boreal_fen | 0.12 | tree leaf |
| f_methane_boreal_bog | 0.08 | tree leaf |
| f_methane_rice_paddy | 0.10 | tree leaf |
| f_methane_upland_soil | 0.05 | tree leaf |
| q10methane | 2.0 | extra Q10 vs aerobes |
| q10methane_base | 295.0 | K, Q10 reference |
| z0_methane_prod | 0.30 * | m, depth e-folding |
| cnscalefactor | 1.0 | substrate scaling |
| mino2lim | 0.2 | anaerobic decomp floor |
| anoxia | T | anoxia framework |
| anoxicmicrosites | F | unsat microsites |
| oxinhib | 400.0 | O2 inhibition scale |
| usephfact | T | pH factor |
| pHmin / pHmax | 2.2 / 9.0 | pH factor range |
| use_spatial_ph | F | no spatial pH map |
| methane_rmcnlim | F | keep CN limit |
| wetland_fixed_substrate | T * | SOM frozen at init |
| use_ch4_sif | F * | SIF substrate path off |
| bgc_anoxia_limits_decomp | T * | anoxia feeds back on BGC |
| use_nitrif_denitrif | T | O2 competition |
| enable_rice_paddy | T * | rice channel on |
| rice_substrate_boost | 1.0 | neutral (retired 68a507f8) |
| rice_drain_window_days | 30 | drain window, d |
| methane_frzout | F | no freeze-out flux |

## 5 Redox lag

| param | value | note |
| --- | --- | --- |
| use_vertical_redoxlag | T | per-layer lag |
| redoxlag | 30 | d, base |
| redoxlag_vertical | 30 | d, per-layer |
| use_biome_redoxlag | T * | biome tree on |
| redoxlag_tropical_peat | 12 | tree leaf |
| redoxlag_tropical_floodplain | 7 | tree leaf |
| redoxlag_temperate_marsh | 20 | tree leaf |
| redoxlag_boreal_fen | 30 | tree leaf |
| redoxlag_boreal_bog | 45 | tree leaf |
| redoxlag_rice_paddy | 25 | tree leaf |
| redoxlag_upland_soil | 30 | tree leaf |

## 6 Oxidation

| param | value | note |
| --- | --- | --- |
| vmax_methane_oxid | 1.25e-5 | mol/m3/s, saturated |
| vmax_oxid_unsat | 1.25e-6 | mol/m3/s, unsaturated |
| k_m | 5e-3 | CH4 half-sat, saturated |
| k_m_unsat | 5e-4 | CH4 half-sat, unsat |
| k_m_o2 | 2e-2 | O2 half-sat |
| q10_methane_oxid | 1.9 | — |
| smp_crit | -2.4e5 | mm, drought cutoff |

## 7 Transport

| param | value | note |
| --- | --- | --- |
| aere_radius | 2.9e-3 | m, aerenchyma radius |
| rob | 3.0 | root obliquity |
| poros_tiller | 0.3 | tiller porosity |
| nongrassporosratio | 1/3 | non-grass porosity ratio |
| unsat_aere_ratio | 1/6 | unsat aere efficiency |
| porosmin | 0.05 | porosity floor |
| tiller_C | 0.22 | tiller carbon |
| scale_factor_aere | 1.0 | aere flux scaling |
| aereoxid | 0.0 | fixed aere oxidation (unused) |
| use_aereoxid_prog | T | prognostic aere oxidation |
| transpirationloss | T | transpiration export |
| vgc_max | 0.15 | ebullition VGC threshold |
| bubble_f | 0.57 | ebullition partition |
| scale_factor_gasdiff | 1.0 | gas diffusion scaling |
| scale_factor_liqdiff | 1.0 | liquid diffusion scaling |
| satpow | 2.0 | diffusivity sat. power |
| capthick | 100 | mm, capillary fringe |
| grnd_methane_cond_default | 1e-6 | surface conductance |
| om_frac_sf | 1.0 | OM fraction scaling |

## 8 Hydrology-inundation

Saturated mode keeps the wetland column at full saturation; the mapping
parameters below act on soil patches / other modes only.

| param | value | note |
| --- | --- | --- |
| inundation_mode | saturated * | wetland column saturated |
| wtd_inflection | 0.30 | m, WTD S-curve midpoint |
| wtd_steepness | 0.05 | m, S-curve width |
| wtd_inflection_soil | 0.0 | soil-patch variant |
| wtd_steepness_soil | 0.3 | soil-patch variant |
| slopemax | 0.4 | slope cap on finundated |
| slopebeta | -3 | slope decay exponent |
| use_routing_for_soil | F | no routed depth |
| hybrid_soil_threshold | 0.05 * | hybrid-mode threshold |
| enable_wetwat_finundated_override | F | wetwat override off |
| wetland_dry_unsat_branch | F | dry-wetland branch off |
| only_wetland | F | soil patches also active |

## 9 Lake

| param | value | note |
| --- | --- | --- |
| allowlakeprod | T * | lake production on |
| lake_decomp_fact | 9e-11 | /s, sediment decomp |
| q10lakebase | 298 | K, reference |
| lake_oxid_scale | 1.0 | oxidation scaling |
| lake_k_m_o2 | -1 | inherit soil value |
| lake_vmax_methane_oxid | -1 | inherit soil value |
| lake_oxic_sediment_depth | -1 | inherit soil value |
| lake_liqdiff_scale | 1.0 | CH4 liq diffusion |
| lake_o2_liqdiff_scale | 1.0 | O2 liq diffusion |
| replenishlakec | F | no C replenishment |
| lake_zero_depth_fatal | T * | zero depth aborts |
| lake_restart_debug | F | probe, unused |
| lake_restart_debug_year / lake_restart_debug_start_doy | 1996 / 1 | probe, unused |
| lake_restart_debug_end_doy / lake_restart_debug_sec | 3 / 1800 | probe, unused |

## 10 Microbial pools (all off at baseline)

| param | value | param | value |
| --- | --- | --- | --- |
| use_microbial_pools | F | use_microbial_flux_override | F |
| use_microbial_dormancy | F | max_microbe_prod_multiplier | 3.0 |
| B_init_methanogen | 1.0 | B_init_methanotroph | 1.0 |
| B_min_methanogen | 1e-3 | B_min_methanotroph | 1e-3 |
| B_max_fraction_methanogen | -1 | B_max_fraction_methanotroph | -1 |
| mu_max_methanogen | 0.2 | mu_max_methanotroph | 0.5 |
| gamma_methanogen | 0.05 | gamma_methanotroph | 0.10 |
| gamma_microbial_dormant | 0.005 | gamma_microbial_freeze | 0.001 |
| k_substrate_methanogen_pool | 0.04 | k_inh_o2_methanogen | 1e-3 |
| kappa_m_methanogen | 1e-2 | kappa_m_methanotroph | 1e-2 |
| q10_microbe_growth | 3.0 | t_ref_microbe | 298.15 |
| dormancy_rate_active | 0.1 | dormancy_rate_revive | 1.0 |
| dormancy_threshold_methanogen_fs | 0.1 | dormancy_threshold_methanogen_fo2 | 0.1 |
| dormancy_threshold_methanotroph_fs | 0.1 | dormancy_threshold_methanotroph_fo2 | 0.1 |

## 11 Boundary & output

| param | value | note |
| --- | --- | --- |
| methane_offline | T | fixed atmosphere |
| atm_methane | 1.7e-6 | mol/mol |
| use_transient_atm_methane | F | no yearly series |
| atm_methane_file | null | unused |
| atm_methane_file_units | auto | unused |
| write_ch4_history | T | — |
| ch4_history_vars | all * | full variable set |
| numerical_correction_fatal_threshold | 1e-3 * | conservation abort |

## 12 Legacy

| param | value | note |
| --- | --- | --- |
| vdcf | 2.0 | backward compat only |
| pc | 0.4 | backward compat only |

## 13 Tracer registration & forcing

| item | value |
| --- | --- |
| DEF_TRACER unit_kind / mol_weight | species_owned / 16.04 |
| DEF_TRACER ref_ratio / init_delta / reactive_decay_rate | 1.0 / 0.0 / 0.0 |
| forcing 1 | GIEMS_MC_v1.1 wetland_frac, monthly, nearest, direct |
| forcing 2 | atm_ch4_global ch4_molmol, monthly, linear, direct |

## A Explicit overrides (the * rows)

| param | baseline | code default |
| --- | --- | --- |
| inundation_mode | saturated | hybrid |
| wetland_fixed_substrate | T | F |
| use_biome_f_methane | T | F |
| use_biome_redoxlag | T | F |
| allowlakeprod | T | F |
| enable_rice_paddy | T | F |
| hybrid_soil_threshold | 0.05 | 0.0 |
| bgc_anoxia_limits_decomp | T | F |
| use_ch4_sif | F | T |
| z0_methane_prod | 0.30 | 0.0 |
| numerical_correction_fatal_threshold | 1e-3 | -1 (off) |
| lake_zero_depth_fatal | T | F |
| ch4_history_vars | all | core |
| enable_rice_paddy | T | F |
| hybrid_soil_threshold | 0.05 | 0.0 |
