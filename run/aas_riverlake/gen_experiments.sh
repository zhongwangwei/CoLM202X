#!/bin/bash
# ==============================================================================
#  从模板 aas_E7_ALL.nml 生成全部 8 组实验的 namelist
# ------------------------------------------------------------------------------
#  只替换 5 处：DEF_CASE_NAME、DEF_dir_output、DEF_dir_restart、
#  三个实验开关、以及 E0_CTL 的历史变量清单。
#  其余字段 8 组逐字节相同 —— 这是"差异只来自三个开关"这一论断的技术保证。
#
#  用法：  ./gen_experiments.sh
#  重跑安全：每次都从模板重新生成，不会累积修改。
# ==============================================================================
set -euo pipefail

cd "$(dirname "$0")"
TEMPLATE=aas_E7_ALL.nml

[ -f "$TEMPLATE" ] || { echo "错误：找不到模板 $TEMPLATE" >&2; exit 1; }

# 实验名:Reservoir_Method:USE_LEVEE:USE_BIFURCATION
EXPS="
E0_CTL:0:.false.:.false.
E1_D:1:.false.:.false.
E2_L:0:.true.:.false.
E3_B:0:.false.:.true.
E4_DL:1:.true.:.false.
E5_DB:1:.false.:.true.
E6_LB:0:.true.:.true.
E7_ALL:1:.true.:.true.
"

for line in $EXPS; do
    name=$(echo "$line" | cut -d: -f1)
    dam=$( echo "$line" | cut -d: -f2)
    lev=$( echo "$line" | cut -d: -f3)
    bif=$( echo "$line" | cut -d: -f4)

    out="aas_${name}.nml"
    [ "$name" = "E7_ALL" ] && continue     # 模板本身即 E7_ALL，不覆盖

    sed -e "s|DEF_CASE_NAME = 'E7_ALL'|DEF_CASE_NAME = '${name}'|" \
        -e "s|/output/E7_ALL/|/output/${name}/|" \
        -e "s|/restart/E7_ALL/|/restart/${name}/|" \
        -e "s|^   DEF_Reservoir_Method = 1|   DEF_Reservoir_Method = ${dam}|" \
        -e "s|^   DEF_USE_LEVEE        = .true.|   DEF_USE_LEVEE        = ${lev}|" \
        -e "s|^   DEF_USE_BIFURCATION  = .true.|   DEF_USE_BIFURCATION  = ${bif}|" \
        "$TEMPLATE" > "$out"

    # E0_CTL 额外输出陆面场（8 组陆面解相同，只需一份）
    # 只改 DEF_HIST_vars_namelist 那一行，不要动注释里出现的同名文字
    if [ "$name" = "E0_CTL" ]; then
        sed -i.bak '/^   DEF_HIST_vars_namelist/s|hist_vars_routing\.nml|hist_vars_land_full.nml|' "$out"
        rm -f "${out}.bak"
    fi

    echo "生成 $out  (dam=${dam} levee=${lev} bif=${bif})"
done

# md5 工具在 Linux 是 md5sum，macOS 是 md5
if command -v md5sum >/dev/null 2>&1; then MD5=md5sum; else MD5=md5; fi

# 实际生效的配置行（行首三空格 + 变量名），排除注释
CFG='^   (DEF_CASE_NAME|DEF_Reservoir_Method|DEF_USE_LEVEE|DEF_USE_BIFURCATION)'
VAR='^   (DEF_CASE_NAME|DEF_dir_output|DEF_dir_restart|DEF_Reservoir_Method|DEF_USE_LEVEE|DEF_USE_BIFURCATION|DEF_HIST_vars_namelist)'

ALL="aas_E0_CTL.nml aas_E1_D.nml aas_E2_L.nml aas_E3_B.nml \
     aas_E4_DL.nml aas_E5_DB.nml aas_E6_LB.nml aas_E7_ALL.nml"

echo ""
echo "=== 实验矩阵核对 ==="
printf "%-14s %-22s %-22s %s\n" 实验 Reservoir_Method USE_LEVEE USE_BIFURCATION
for f in $ALL; do
    printf "%-14s %-22s %-22s %s\n" \
      "$(grep -E "$CFG" "$f" | grep DEF_CASE_NAME       | sed "s/.*= *//; s/'//g")" \
      "$(grep -E "$CFG" "$f" | grep Reservoir_Method    | sed 's/.*= *//')" \
      "$(grep -E "$CFG" "$f" | grep DEF_USE_LEVEE       | sed 's/.*= *//')" \
      "$(grep -E "$CFG" "$f" | grep DEF_USE_BIFURCATION | sed 's/.*= *//')"
done

echo ""
echo "=== 一致性核对：除 case 名/路径/三开关/变量清单外，8 组必须逐字节相同 ==="
n=$(for f in $ALL; do grep -vE "$VAR" "$f" | $MD5 | awk '{print $1}'; done | sort -u | wc -l | tr -d ' ')
if [ "$n" = "1" ]; then
    echo "通过：8 份 namelist 的其余部分完全一致"
else
    echo "！！失败：其余部分存在 $n 种差异，实验将无法归因到三个开关。请检查。"
    exit 1
fi
