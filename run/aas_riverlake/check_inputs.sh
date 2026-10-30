#!/bin/bash
# ==============================================================================
#  前置检查：确认汇流参数文件含有 levee / bifurcation / dam 所需的全部字段
# ------------------------------------------------------------------------------
#  为什么必须先跑这个：
#    缺字段不会在 namelist 解析阶段报错，而是在 river_lake_network_init 读参数
#    时才失败 —— 那时 480 核作业已经排队、启动、读完地表数据，浪费的是小时级
#    的机时和排队时间。本脚本几秒钟就能给出答案。
#
#    仓库里 run/Amazon.nml:39 用的是 grid_routing_data_15min_sed_levee.nc（含 levee），
#    而 run/Amazon_iso_0p5_2020.nml:55 用的是 grid_routing_data_15min.nc ——
#    文件名不同说明不是同一份数据，不能假定后者也含 levee/bifurcation 字段。
#
#  用法：  ./check_inputs.sh [参数文件路径]
#          不给参数时使用下面的默认路径。
# ==============================================================================
set -uo pipefail

PARA=${1:-/share/home/dq013/zhwei/colm/data/CoLMruntime/unitcatchment/grid_routing_data_15min.nc}

command -v ncdump >/dev/null 2>&1 || { echo "错误：找不到 ncdump" >&2; exit 1; }
[ -f "$PARA" ] || { echo "错误：参数文件不存在 -> $PARA" >&2; exit 1; }

echo "=============================================================="
echo " 检查参数文件：$PARA"
echo "=============================================================="

HDR=$(ncdump -h "$PARA")

# 变量名:用途:对应开关
CHECKS="
seq_x:河网拓扑:总是
seq_y:河网拓扑:总是
seq_next:河网拓扑:总是
inpmat_x:陆面->汇流映射:总是
inpmat_y:陆面->汇流映射:总是
inpmat_area:陆面->汇流映射:总是
topo_rivelv:河道地形:总是
topo_rivhgt:河道地形:总是
topo_rivlen:河道地形:总是
topo_rivman:河道地形:总是
topo_rivwth:河道地形:总是
topo_rivstomax:河道地形:总是
topo_area:洪泛地形:总是
topo_fldhgt:洪泛地形:总是
levee_frc:堤防:DEF_USE_LEVEE
levee_hgt:堤防:DEF_USE_LEVEE
bifurcation_upst:分汊:DEF_USE_BIFURCATION
bifurcation_down:分汊:DEF_USE_BIFURCATION
bifurcation_distance:分汊:DEF_USE_BIFURCATION
bifurcation_elevation:分汊:DEF_USE_BIFURCATION
bifurcation_width:分汊:DEF_USE_BIFURCATION
bifurcation_manning:分汊:DEF_USE_BIFURCATION
dam_seq:水库:DEF_Reservoir_Method
dam_year:水库:DEF_Reservoir_Method
dam_GRAND_ID:水库:DEF_Reservoir_Method
dam_TotalVol_mcm:水库:DEF_Reservoir_Method
dam_ConVol_mcm:水库:DEF_Reservoir_Method
dam_Qn:水库:DEF_Reservoir_Method
dam_Qf:水库:DEF_Reservoir_Method
"

MISS_CORE=0; MISS_LEV=0; MISS_BIF=0; MISS_DAM=0
printf "%-24s %-16s %-24s %s\n" 变量 用途 需要它的开关 状态
printf -- "--------------------------------------------------------------------------------\n"

for line in $CHECKS; do
    v=$(  echo "$line" | cut -d: -f1)
    use=$(echo "$line" | cut -d: -f2)
    sw=$( echo "$line" | cut -d: -f3)

    if echo "$HDR" | grep -qE "[[:space:]]${v}\(|[[:space:]]${v}[[:space:]]*;"; then
        st="OK"
    else
        st="缺失"
        case "$sw" in
            总是)                  MISS_CORE=$((MISS_CORE+1)) ;;
            DEF_USE_LEVEE)         MISS_LEV=$((MISS_LEV+1))   ;;
            DEF_USE_BIFURCATION)   MISS_BIF=$((MISS_BIF+1))   ;;
            DEF_Reservoir_Method)  MISS_DAM=$((MISS_DAM+1))   ;;
        esac
    fi
    printf "%-24s %-16s %-24s %s\n" "$v" "$use" "$sw" "$st"
done

echo ""
echo "=============================================================="
RC=0
if [ $MISS_CORE -gt 0 ]; then
    echo "！！致命：缺 $MISS_CORE 个核心字段，任何实验都跑不了。"; RC=1
else
    echo "核心河网/地形字段：完整（E0_CTL 可跑）"
fi

[ $MISS_DAM -gt 0 ] && { echo "！！缺 $MISS_DAM 个水库字段 -> E1/E4/E5/E7 无法运行"; RC=1; } \
                    || echo "水库字段：完整（E1_D E4_DL E5_DB E7_ALL 可跑）"
[ $MISS_LEV -gt 0 ] && { echo "！！缺 $MISS_LEV 个堤防字段 -> E2/E4/E6/E7 无法运行"; RC=1; } \
                    || echo "堤防字段：完整（E2_L E4_DL E6_LB E7_ALL 可跑）"
[ $MISS_BIF -gt 0 ] && { echo "！！缺 $MISS_BIF 个分汊字段 -> E3/E5/E6/E7 无法运行"; RC=1; } \
                    || echo "分汊字段：完整（E3_B E5_DB E6_LB E7_ALL 可跑）"
echo "=============================================================="

if [ $RC -ne 0 ]; then
    echo ""
    echo "处理建议："
    echo "  1) 先看看同目录下是否有更完整的版本，例如："
    echo "       ls -l $(dirname "$PARA")/"
    echo "     仓库 run/Amazon.nml 用的是 grid_routing_data_15min_sed_levee.nc"
    echo "  2) 若确实没有，需从 MERIT-Hydro / CaMa-Flood v4 的全球 15 min 数据集重新制备。"
    echo "  3) 在补齐之前，只能跑核心字段支持的那几组实验。"
fi

exit $RC
