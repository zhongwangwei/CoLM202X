#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

check_main() {
  file=$1

  rg -U "wet_area\\s*=\\s*min\\(1\\._r8, max\\(lai, 0\\._r8\\)\\)\\n\\s*wet_area_cfw = wet_area" "$file" >/dev/null
  rg -n "dewmx_MATSIRO\\s*=\\s*0\\.2_r8 \\* DEF_MATSIRO_CWCAP_SCALE" "$file" >/dev/null
  rg -n "satcap_snow_eff\\s*=\\s*dewmx_MATSIRO \\* max\\(lai, 0\\._r8\\)" "$file" >/dev/null
}

check_pc() {
  file=$1

  rg -U "wet_area\\s*=\\s*min\\(1\\._r8, max\\(lai\\(i\\), 0\\._r8\\)\\)\\n\\s*wet_area_cfw = wet_area" "$file" >/dev/null
  rg -n "wet_area_cfw\\s*=\\s*min\\(1\\._r8, max\\(lai\\(i\\), 0\\._r8\\)\\)" "$file" >/dev/null
  rg -n "dewmx_MATSIRO\\s*=\\s*0\\.2_r8 \\* DEF_MATSIRO_CWCAP_SCALE" "$file" >/dev/null
  rg -n "satcap_snow_eff\\s*=\\s*dewmx_MATSIRO \\* max\\(lai, 0\\._r8\\)" "$file" >/dev/null
}

check_namelist() {
  file=$1

  rg -n "real\\(r8\\)\\s*::\\s*DEF_MATSIRO_CWCAP_SCALE\\s*=" "$file" >/dev/null
  rg -n "DEF_MATSIRO_CWCAP_SCALE" "$file" >/dev/null
  rg -n "CALL mpi_bcast \\(DEF_MATSIRO_CWCAP_SCALE" "$file" >/dev/null
}

check_interception() {
  file=$1

  rg -n "USE MOD_Namelist, only: .*DEF_MATSIRO_CWCAP_SCALE" "$file" >/dev/null
  rg -n "dewmx_MATSIRO\\s*=\\s*0\\.2_r8 \\* DEF_MATSIRO_CWCAP_SCALE" "$file" >/dev/null
}

check_main "$repo_dir/main/MOD_LeafTemperature.F90"
check_pc "$repo_dir/main/MOD_LeafTemperaturePC.F90"
check_namelist "$repo_dir/share/MOD_Namelist.F90"
check_interception "$repo_dir/main/MOD_LeafInterception.F90"
