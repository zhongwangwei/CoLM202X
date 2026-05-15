#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

check_main() {
  file=$1

  rg -n "ELSEIF \\(DEF_Interception_scheme == 4\\) THEN" "$file" >/dev/null
  rg -n "wet_area_cfw" "$file" >/dev/null
  rg -n "wet_cond_cfw" "$file" >/dev/null
  rg -n "wet_area\\s*=\\s*min\\(6\\._r8, max\\(lai \\+ sai, 0\\._r8\\)\\)" "$file" >/dev/null
  rg -n "wet_area_cfw\\s*=\\s*wet_area" "$file" >/dev/null
  rg -n "wet_cond_cfw\\s*=\\s*wet_area_cfw / rb" "$file" >/dev/null
  rg -n "cfw\\s*=\\s*\\(1\\.-delta\\*\\(1\\.-fwet\\)\\)\\*wet_cond_cfw" "$file" >/dev/null
}

check_pc() {
  file=$1

  rg -n "ELSEIF \\(DEF_Interception_scheme == 4\\) THEN" "$file" >/dev/null
  rg -n "wet_area_cfw" "$file" >/dev/null
  rg -n "wet_cond_cfw" "$file" >/dev/null
  rg -n "wet_area\\s*=\\s*min\\(6\\._r8, max\\(lsai\\(i\\), 0\\._r8\\)\\)" "$file" >/dev/null
  rg -n "wet_area_cfw\\s*=\\s*wet_area" "$file" >/dev/null
  rg -n "wet_cond_cfw\\s*=\\s*wet_area_cfw / rb\\(i\\)" "$file" >/dev/null
  rg -n "cfw\\(i\\)\\s*=\\s*\\(1\\.-delta\\(i\\)\\*\\(1\\.-fwet\\(i\\)\\)\\)\\*wet_cond_cfw" "$file" >/dev/null
}

check_main "$repo_dir/main/MOD_LeafTemperature.F90"
check_pc "$repo_dir/main/MOD_LeafTemperaturePC.F90"
