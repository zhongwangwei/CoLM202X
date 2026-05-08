#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

check_thermal_main() {
  file=$1

  rg -U -P "ELSEIF \\(DEF_Interception_scheme == 7\\) THEN(?:(?:\\n|.){0,240})wet_area = 1\\._r8" "$file" >/dev/null
  rg -n "epot_JULES\\s*=\\s*rhoair / max\\(raw, 1\\.e-10_r8\\) \\* \\(qsatl - qaf\\)" "$file" >/dev/null
  rg -n "wet_cond\\s*=\\s*1\\._r8 / max\\(raw, 1\\.e-10_r8\\)" "$file" >/dev/null
}

check_thermal_pc() {
  file=$1

  rg -U -P "ELSEIF \\(DEF_Interception_scheme == 7\\) THEN(?:(?:\\n|.){0,240})wet_area = 1\\._r8" "$file" >/dev/null
  rg -n "epot_JULES\\s*=\\s*rhoair / max\\(raw, 1\\.e-10_r8\\) \\* \\(qsatl\\(i\\) - qaf\\(clev\\)\\)" "$file" >/dev/null
  rg -n "wet_cond\\s*=\\s*1\\._r8 / max\\(raw, 1\\.e-10_r8\\)" "$file" >/dev/null
}

check_phs() {
  file=$1

  rg -U "CASE \\(7\\)\\n\\s*wet_area\\s*=\\s*1\\._r8" "$file" >/dev/null
  rg -U "CASE \\(7\\)\\n\\s*wet_cond\\s*=\\s*1\\._r8 / max\\(raw, 1\\.e-10_r8\\)" "$file" >/dev/null
}

check_interception() {
  file=$1

  rg -n "r_rain_ls\\s*=\\s*MAX\\(0\\.0_r8, prl_rain \\+ qflx_irrig_sprinkler\\)" "$file" >/dev/null
  rg -n "r_rain_con\\s*=\\s*MAX\\(0\\.0_r8, prc_rain\\)" "$file" >/dev/null
  rg -n "ldew_snow_pre\\s*=\\s*ldew_snow" "$file" >/dev/null
  rg -n "unload_snow\\s*=\\s*snowunloadfact \\* melt_rate \\* deltim" "$file" >/dev/null
  rg -n "intercept_snow\\s*=\\s*snowinterceptfact \\* \\(can_cpy_snow - ldew_snow_pre\\)" "$file" >/dev/null
  rg -n "ldew_snow\\s*=\\s*ldew_snow_pre \\+ intercept_snow - unload_snow" "$file" >/dev/null
}

check_thermal_main "$repo_dir/main/MOD_LeafTemperature.F90"
check_thermal_pc "$repo_dir/main/MOD_LeafTemperaturePC.F90"
check_phs "$repo_dir/main/MOD_PlantHydraulic.F90"
check_interception "$repo_dir/main/MOD_LeafInterception.F90"
