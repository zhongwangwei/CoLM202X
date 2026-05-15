#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

fail() { echo "ERROR: $*" >&2; exit 1; }
need_file() { [[ -f "$1" ]] || fail "missing required file: $1"; }
need_rg() { local pat="$1"; shift; rg -n -- "$pat" "$@" >/dev/null || fail "missing pattern '$pat' in $*"; }

core_files=(
  main/TRACER/MOD_Tracer_Methane_Const.F90
  main/TRACER/MOD_Tracer_Methane_Physics.F90
  main/TRACER/MOD_Tracer_Methane_Driver.F90
  main/TRACER/MOD_Tracer_Methane_State.F90
  main/TRACER/MOD_Tracer_Methane_Microbes.F90
  main/TRACER/MOD_Tracer_Methane_AccFlux.F90
  main/TRACER/MOD_Tracer_Methane_Registry.F90
  main/TRACER/MOD_Tracer_Methane_BgcLink.F90
  run/standard_ch4_parameter.nml
  run/examples/Global_Grid_2x2_PFT_VG_BGC_CH4.nml
  methane-dssat-csm-reference-20260515.md
)
for f in "${core_files[@]}"; do need_file "$f"; done

# Build graph: all runtime pieces must be linked into colm.x.
for obj in MOD_Tracer_Methane_Const.o MOD_Tracer_Methane_Registry.o \
           MOD_Tracer_Methane_State.o MOD_Tracer_Methane_Microbes.o \
           MOD_Tracer_Methane_AccFlux.o \
           MOD_Tracer_Methane_Physics.o MOD_Tracer_Methane_BgcLink.o \
           MOD_Tracer_Methane_Driver.o; do
  need_rg "${obj}" Makefile
done

# Startup/finalize/driver/history/accumulator hooks.
need_rg '^MODULE MOD_Tracer_Methane_Driver' main/TRACER/MOD_Tracer_Methane_Driver.F90
need_rg 'CALL methane_registry_init\(' main/CoLM.F90
need_rg 'CALL allocate_methane_state' main/CoLM.F90
need_rg 'CALL allocate_methane_microbes_state' main/CoLM.F90
need_rg 'CALL read_methane_restart' main/CoLM.F90
need_rg 'CALL read_methane_microbes_restart' main/CoLM.F90
need_rg 'CALL read_methane_namelist' main/CoLM.F90
need_rg 'open\(newunit=unit_nml' main/TRACER/MOD_Tracer_Methane_Const.F90
need_rg 'DEF_USE_METHANE_para' share/MOD_Namelist.F90 main/CoLM.F90
need_rg 'DEF_file_METHANE_para' share/MOD_Namelist.F90 main/CoLM.F90
need_rg 'DEF_METHANE_only_wetland' share/MOD_Namelist.F90 main/CoLMDRIVER.F90
need_rg 'DEF_wetland_finundation_scheme' share/MOD_Namelist.F90 main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'mpi_bcast .*DEF_USE_METHANE_para' share/MOD_Namelist.F90
need_rg 'mpi_bcast .*DEF_file_METHANE_para' share/MOD_Namelist.F90
need_rg 'CALL deallocate_methane_state' main/CoLM.F90
need_rg 'CALL deallocate_methane_microbes_state' main/CoLM.F90
need_rg 'CALL methane_driver' main/CoLMDRIVER.F90
need_rg 'IF \(igas_ch4 > 0\)' main/CoLMDRIVER.F90 main/CoLM.F90 main/MOD_Hist.F90 main/MOD_Vars_1DAccFluxes.F90
need_rg 'CALL flush_methane_acc_fluxes' main/MOD_Vars_1DAccFluxes.F90
need_rg 'CALL accumulate_methane_fluxes' main/MOD_Vars_1DAccFluxes.F90
need_rg 'CALL write_methane_restart' main/MOD_Vars_TimeVariables.F90
need_rg 'CALL write_methane_microbes_restart' main/MOD_Vars_TimeVariables.F90
need_rg 'f_methane_surf_flux_tot' main/MOD_Hist.F90
need_rg 'DEF_TRACER_NAMES[[:space:]]*=[[:space:]]*"CH4,O2,CO2"' run/examples/Global_Grid_2x2_PFT_VG_BGC_CH4.nml
need_rg 'DEF_file_METHANE_para' run/examples/Global_Grid_2x2_PFT_VG_BGC_CH4.nml

# Methane defensive-fix regressions from review.
need_rg 'forc_pmethanem <= 0\._r8' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'forc_pmethanem[[:space:]]*=[[:space:]]*DEF_METHANE%atm_methane\*forc_pbot' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'real\(r8\), intent\(out\) :: crootfr' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
need_rg 'pH[[:space:]]*=[[:space:]]*6\.2_r8' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
if rg -n 'usephfact requires valid soil pH|pH[[:space:]]*=[[:space:]]*spval' main/TRACER/MOD_Tracer_Methane_BgcLink.F90 main/TRACER/MOD_Tracer_Methane_Driver.F90; then
  fail 'usephfact must use the BgcLink neutral default until a spatial soil-pH field is wired'
fi
need_rg 'o_scalar_out\(j\)[[:space:]]*=[[:space:]]*min\(1\._r8, max\(0\._r8, o_scalar\(j, ipatch\)\)\)' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
need_rg 'do s=1,ngases' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'c_atm\(3\)[[:space:]]*=[[:space:]]*max\(forc_pco2m, 0\._r8\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'porsl\(j\) <= smallnumber \.or\. dz_soisno\(j\) <= smallnumber' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'slope_angle[[:space:]]*=[[:space:]]*atan\(slope_arg\)' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'slope_angle[[:space:]]*=[[:space:]]*slope_arg \* PI / 180\._r8' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'slope_angle[[:space:]]*=[[:space:]]*atan\(slope_arg / 100\._r8\)' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'conc_ch4_porsl\(j\)[[:space:]]*=[[:space:]]*0\._r8' main/TRACER/MOD_Tracer_Methane_Physics.F90

# Methane lake/CTSM-alignment regressions.
need_rg 'lakedepth,lake_icefrac,wdsrf' main/TRACER/MOD_Tracer_Methane_Driver.F90 main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'lake_soilc_srf' main/MOD_Vars_TimeInvariants.F90 mkinidata/MOD_Initialize.F90 main/CoLM.F90
need_rg 'lake_soilc_patches' mkinidata/MOD_Initialize.F90 mksrfdata/Aggregation_LakeSoilC.F90 mksrfdata/MKSRFDATA.F90
need_rg 'q10lake_eff[[:space:]]*=[[:space:]]*DEF_METHANE%q10methane[[:space:]]*\*[[:space:]]*1\.5_r8' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'lake_ch4_fraction[[:space:]]*=[[:space:]]*min\(max\(DEF_METHANE%f_methane' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'methane_prod_depth\(j\)[[:space:]]*=[[:space:]]*lake_ch4_fraction \* base_decomp / dz_soisno\(j\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'co2_decomp_depth\(j\)[[:space:]]*=[[:space:]]*max\(0\._r8, base_decomp / dz_soisno\(j\) - methane_prod_depth\(j\)\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'lake_soilc\(j\)[[:space:]]*=[[:space:]]*max\(0\._r8, lake_soilc\(j\) - &' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'methane_prod_depth_sat\(j\) \+ co2_decomp_depth_sat\(j\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
if rg -n 'DEF_METHANE%q10lake([^a-zA-Z0-9_]|$)' main/TRACER/MOD_Tracer_Methane_Physics.F90; then
  fail 'physics must derive effective lake Q10 from q10methane, not the legacy q10lake namelist field'
fi
need_rg 'patchtype == 4 \.and\. \.not\. DEF_METHANE%allowlakeprod' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'patchtype == 4 \.and\. DEF_METHANE%allowlakeprod' main/TRACER/MOD_Tracer_Methane_Physics.F90
if rg -n 'lake_soilc[[:space:]]*=[[:space:]]*580\._r8 \* organic_max \* max\(vf_om' main/TRACER/MOD_Tracer_Methane_Physics.F90; then
  fail 'runtime lake_soilc fallback must use formal lake_soilc_srf surface/restart field, not organic_max/vf_om'
fi
need_rg 'SUBROUTINE initialize_methane_lake_soilc_from_surface' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'CALL initialize_methane_lake_soilc_from_surface' main/CoLM.F90 main/LULCC/MOD_Lulcc_Driver.F90
need_rg 'lake_soilc_srf_in' main/TRACER/MOD_Tracer_Methane_State.F90
if rg -n 'organic_max_in|vf_om_in|580\._r8 \* organic_max' main/TRACER/MOD_Tracer_Methane_State.F90 main/TRACER/MOD_Tracer_Methane_Physics.F90; then
  fail 'lake_soilc initialization must not depend on legacy organic_max/vf_om fallback'
fi
need_rg 'CALL read_methane_namelist' main/CoLM.F90
need_rg 'CALL save_methane_lulcc_state' main/LULCC/MOD_Lulcc_Driver.F90
need_rg 'CALL save_methane_microbes_lulcc_state' main/LULCC/MOD_Lulcc_Driver.F90
need_rg 'CALL remap_methane_lulcc_state' main/LULCC/MOD_Lulcc_Driver.F90
need_rg 'CALL remap_methane_microbes_lulcc_state' main/LULCC/MOD_Lulcc_Driver.F90
need_rg 'lccpct_patches' main/LULCC/MOD_Lulcc_Driver.F90
need_rg 'SUBROUTINE remap_methane_lulcc_state' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'CALL allocate_methane_state[[:space:]]*\(nnew\)' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'lake_icefrac\(1\)[[:space:]]*>[[:space:]]*0\.1_r8' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'vgc[[:space:]]*>[[:space:]]*DEF_METHANE%vgc_max[[:space:]]*\*[[:space:]]*DEF_METHANE%bubble_f' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'methane_ebul_depth\(j\)[[:space:]]*=[[:space:]]*\(vgc - DEF_METHANE%vgc_max[[:space:]]*\*[[:space:]]*DEF_METHANE%bubble_f\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'Legacy compatibility field; effective lake Q10 is q10methane \* 1\.5' run/standard_ch4_parameter.nml main/TRACER/MOD_Tracer_Methane_Const.F90

# CO2 diagnostic closure must be carried through physics, accumulation, and history.
need_rg 'carbon_decomp_depth - methane_prod_depth\(j\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
if rg -n 'co2_decomp_depth_unsat[[:space:]]*=[[:space:]]*o2_decomp_depth_unsat[[:space:]]*\+[[:space:]]*methane_prod_depth_unsat' main/TRACER/MOD_Tracer_Methane_Physics.F90; then
  fail 'CO2 decomposition diagnostic must not use o2_decomp + methane_prod; that includes nitrification/anoxia terms and has the wrong CH4 sign'
fi
need_rg 'lulcc_co2_decomp_depth_old' main/TRACER/MOD_Tracer_Methane_State.F90
need_rg 'co2_oxid_depth_unsat[[:space:]]*=[[:space:]]*methane_oxid_depth_unsat' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'co2_net_tot[[:space:]]*=[[:space:]]*co2_decomp_tot[[:space:]]*\+[[:space:]]*co2_oxid_tot[[:space:]]*-[[:space:]]*co2_aere_tot' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'a_co2_decomp_tot' main/TRACER/MOD_Tracer_Methane_AccFlux.F90 main/MOD_Hist.F90
need_rg 'f_co2_net_tot' main/MOD_Hist.F90
need_rg 'f_co2_net_tot_lake' main/MOD_Hist.F90
need_rg 'Aggregation_LakeSoilC' Makefile mksrfdata/MKSRFDATA.F90
need_rg 'lake_soilc.nc' mksrfdata/Aggregation_LakeSoilC.F90
need_rg 'vf_om_s_l' mksrfdata/Aggregation_LakeSoilC.F90
need_rg 'carbon_per_kg_om \* organic_max_default \* vf_om_s_patches' mksrfdata/Aggregation_LakeSoilC.F90
need_rg 'CALL Aggregation_SoilParameters[[:space:]]*\(grid_soil' mksrfdata/MKSRFDATA.F90
need_rg 'CALL Aggregation_LakeSoilC[[:space:]]*\(grid_soil' mksrfdata/MKSRFDATA.F90


# Optional microbial-pool dynamics are present but default-off; when both
# microbial switches are enabled, capped microbial potentials replace legacy
# production/oxidation rates inside the Methane physics.
need_rg 'MODULE MOD_Tracer_Methane_Microbes' main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'use_microbial_pools[[:space:]]*=[[:space:]]*\.false\.' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml
need_rg 'use_microbial_flux_override[[:space:]]*=[[:space:]]*\.false\.' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml
need_rg 'CALL methane_microbes_step' main/TRACER/MOD_Tracer_Methane_Driver.F90
need_rg 'cellorg' main/TRACER/MOD_Tracer_Methane_Driver.F90 main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'ch4_B_methanogen' main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'f_methane_B_methanogen' main/MOD_Hist.F90
need_rg 'a_microbial_prod_potential' main/TRACER/MOD_Tracer_Methane_AccFlux.F90 main/MOD_Hist.F90
need_rg 'microbial_prod_potential_patch' main/TRACER/MOD_Tracer_Methane_Driver.F90 main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'microbial_oxid_potential_patch' main/TRACER/MOD_Tracer_Methane_Driver.F90 main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'methane_prod_depth\(j\)[[:space:]]*=[[:space:]]*min\(max\(0\._r8, microbial_prod_potential_layer\(j\)\), carbon_decomp_depth\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'methane_oxid_depth\(j\)[[:space:]]*=[[:space:]]*max\(0\._r8, microbial_oxid_potential_layer\(j\)\)' main/TRACER/MOD_Tracer_Methane_Physics.F90
need_rg 'IF[[:space:]]*\(patchtype[[:space:]]*/=[[:space:]]*4\)[[:space:]]*THEN' main/TRACER/MOD_Tracer_Methane_Driver.F90
need_rg 'carbon_cap' main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'o2_cap' main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'K_substrate_methanogen_pool' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'gamma_microbial_freeze' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'dormancy_threshold_methanogen_fS' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'dormancy_threshold_methanotroph_fS' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'dormancy_threshold_methanotroph_fO2' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml main/TRACER/MOD_Tracer_Methane_Microbes.F90
need_rg 'DEPRECATED: legacy rate-based substrate half-saturation' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml
need_rg 'DEPRECATED: legacy shared' main/TRACER/MOD_Tracer_Methane_Const.F90 run/standard_ch4_parameter.nml

# External rice/paddy reference is documented separately from the CTSM
# wetland/lake alignment path and must remain default-off conceptually.
need_rg 'DSSAT-CSM rice Methane reference' methane-dssat-csm-reference-20260515.md
need_rg 'BSD-3-Clause' methane-dssat-csm-reference-20260515.md
need_rg 'rice/paddy extension' ctsm-ch4-alignment-plan-20260515.md methane-dssat-csm-reference-20260515.md

# BGC isolation: Methane driver must go through BgcLink, not direct BGC USE.
if rg -n '^\s*USE\s+MOD_BGC_' main/TRACER/MOD_Tracer_Methane_Driver.F90; then
  fail 'MOD_Tracer_Methane_Driver.F90 must access BGC only through MOD_Tracer_Methane_BgcLink'
fi
need_rg '^\s*USE\s+MOD_BGC_' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
need_rg 'tracer_ch4_bgc_patch_inputs' main/TRACER/MOD_Tracer_Methane_Driver.F90 main/TRACER/MOD_Tracer_Methane_BgcLink.F90
need_rg 'valid_bgc_value' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
need_rg 'pft_arrays_ready' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
need_rg 'safe_nonnegative' main/TRACER/MOD_Tracer_Methane_BgcLink.F90
if rg -n 'PUBLIC :: o_scalar|PUBLIC :: pot_f_nit_vr' main/TRACER/MOD_Tracer_Methane_BgcLink.F90; then
  fail 'BGC bridge must not re-export raw BGC internals'
fi

# BGC driver remains in native soil scope; wetland Methane is protected by
# sanitized bridge inputs rather than executing full BGC on uninitialized
# wetland pools/PFT mappings.
if rg -n 'CALL bgc_driver' main/CoLMDRIVER.F90 -B 8 | rg 'patchtype\(i\).*2'; then
  fail 'bgc_driver must not be widened directly to wetland patchtype==2'
fi

# Methane is not a compile macro guard. Mentions in comments are allowed.
if rg -n '^\s*#\s*(ifn?def|if).*\bCH4\b' include main share mkinidata mksrfdata postprocess Makefile; then
  fail 'found active CH4 preprocessor guard; use TRACER+BGC compile gate and igas_ch4 runtime gate'
fi

# Optional source comparison when the legacy source tree is present.
source_root="${1:-/Users/zhongwangwei/Desktop/CH4}"
if [[ -d "$source_root" ]]; then
  need_file "$source_root/main/MOD_ch4.F90"
  need_file "$source_root/main/MOD_ch4_driver.F90"
  need_file "$source_root/main/MOD_Const_ch4.F90"
fi

echo 'Methane integration checks passed.'
