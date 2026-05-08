#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

check_main() {
  file=$1

  test "$(rg -c "wet_area\\s*=\\s*1\\._r8" "$file")" -ge 2
  rg -U "wet_area\\s*=\\s*1\\._r8\\n\\s*ec_pot_VIC" "$file" >/dev/null
  rg -U "wet_area\\s*=\\s*1\\._r8\\n\\s*wet_area_cfw = wet_area" "$file" >/dev/null
  test "$(rg -c "MaxInt_VIC\\s*=\\s*max\\(0\\.1_r8 \\* .*DEF_VIC_WDMAX_SCALE, 1\\.e-10_r8\\)" "$file")" -eq 2
  rg -n "satcap_rain_eff\\s*=\\s*sigf_safe \\* 0\\.1_r8 \\* max\\(lai_perveg, 0\\._r8\\) \\* DEF_VIC_WDMAX_SCALE" "$file" >/dev/null
  rg -n "satcap_snow_eff\\s*=\\s*sigf_safe \\* 0\\.5_r8 \\* Lr \\* max\\(lai_perveg, 0\\._r8\\)" "$file" >/dev/null
}

check_pc() {
  file=$1

  test "$(rg -c "wet_area\\s*=\\s*1\\._r8" "$file")" -ge 2
  rg -U "wet_area\\s*=\\s*1\\._r8\\n\\s*ec_pot_VIC" "$file" >/dev/null
  rg -U "wet_area\\s*=\\s*1\\._r8\\n\\s*wet_area_cfw = wet_area" "$file" >/dev/null
  rg -n "wet_area_cfw\\s*=\\s*1\\._r8" "$file" >/dev/null
  test "$(rg -c "MaxInt_VIC\\s*=\\s*max\\(0\\.1_r8 \\* .*DEF_VIC_WDMAX_SCALE, 1\\.e-10_r8\\)" "$file")" -eq 2
  rg -n "satcap_rain_eff\\s*=\\s*sigf_safe \\* 0\\.1_r8 \\* max\\(lai_perveg, 0\\._r8\\) \\* DEF_VIC_WDMAX_SCALE" "$file" >/dev/null
  rg -n "satcap_snow_eff\\s*=\\s*sigf_safe \\* 0\\.5_r8 \\* Lr \\* max\\(lai_perveg, 0\\._r8\\)" "$file" >/dev/null
}

check_phs() {
  file=$1

  rg -n "FUNCTION wet_coupling_area_twoleaf" "$file" >/dev/null
  rg -n "FUNCTION wet_coupling_conductance_twoleaf" "$file" >/dev/null
  rg -U "CASE \\(6\\)\\n\\s*wet_area\\s*=\\s*1\\._r8" "$file" >/dev/null
  rg -n "wet_cond_cfw\\s*=\\s*wet_coupling_conductance_twoleaf\\(laisun, laisha, sai, raw, gb_mol, cf\\)" "$file" >/dev/null
}

check_namelist() {
  file=$1

  rg -n "real\\(r8\\) :: DEF_VIC_WDMAX_SCALE = 1\\.0_r8" "$file" >/dev/null
  rg -n "DEF_VIC_WDMAX_SCALE," "$file" >/dev/null
  rg -n "mpi_bcast \\(DEF_VIC_WDMAX_SCALE\\s*,1\\s*,mpi_real8" "$file" >/dev/null
}

check_interception() {
  file=$1

  rg -n "MaxInt\\s*=\\s*0\\.1 \\* lai_perveg \\* DEF_VIC_WDMAX_SCALE" "$file" >/dev/null
}

check_main "$repo_dir/main/MOD_LeafTemperature.F90"
check_pc "$repo_dir/main/MOD_LeafTemperaturePC.F90"
check_phs "$repo_dir/main/MOD_PlantHydraulic.F90"
check_namelist "$repo_dir/share/MOD_Namelist.F90"
check_interception "$repo_dir/main/MOD_LeafInterception.F90"
