#!/bin/bash
# ==============================================================================
#  把共享 spinup (1996-1999) 的末期 restart 分发给 8 组实验
# ------------------------------------------------------------------------------
#  原理：汇流对陆面是单向的（push_ucat2grid 只用于历史输出，不回写 patch），
#        所以 bif/levee/dam 三个开关不影响陆面解，一次 spinup 的陆面状态对
#        8 组完全等价。这一步省下 4x7 = 28 个模式年。
#
#  restart 命名规则（main/MOD_Vars_TimeVariables.F90:1119-1129）：
#        {dir_restart}/{cdate}/{casename}_restart[_子类]_{cdate}_lc{cyear}.nc
#        cdate 格式为 YYYY-DDD-SSSSS（年-儒略日-当日秒数），
#        例如 1999 年末为 1999-365-86400。
#        子类包括 _pft_、_gridriver_ 等，脚本用通配符统一处理。
#
#  用法：  ./distribute_restart.sh            # 预览，不改动任何文件
#          ./distribute_restart.sh --run      # 实际执行
# ==============================================================================
set -euo pipefail

CASE_ROOT=/share/home/dq013/zhwei/colm/cases/aas_riverlake
SPINUP_CASE=SPINUP
EXPS="E0_CTL E1_D E2_L E3_B E4_DL E5_DB E6_LB E7_ALL"

DRYRUN=1
[ "${1:-}" = "--run" ] && DRYRUN=0

SRC_ROOT="${CASE_ROOT}/restart/${SPINUP_CASE}"
[ -d "$SRC_ROOT" ] || { echo "错误：找不到 spinup restart 目录 $SRC_ROOT" >&2; exit 1; }

# 自动取最新的日期子目录（不硬编码，避免闰年/秒数写错）
CDATE=$(ls -1 "$SRC_ROOT" | grep -E '^[0-9]{4}-[0-9]{3}-[0-9]{5}$' | sort | tail -1)
[ -n "$CDATE" ] || { echo "错误：$SRC_ROOT 下没有形如 YYYY-DDD-SSSSS 的子目录" >&2; exit 1; }

SRC="${SRC_ROOT}/${CDATE}"
NFILE=$(ls -1 "$SRC"/${SPINUP_CASE}_restart*.nc 2>/dev/null | wc -l | tr -d ' ')

echo "=============================================================="
echo " 源目录 : $SRC"
echo " 日期戳 : $CDATE"
echo " 文件数 : $NFILE"
echo "=============================================================="
[ "$NFILE" -gt 0 ] || { echo "错误：源目录下没有 ${SPINUP_CASE}_restart*.nc" >&2; exit 1; }

echo ""
echo "源文件清单："
ls -1 "$SRC"/${SPINUP_CASE}_restart*.nc | sed 's|.*/|  |'
echo ""

# 确认这确实是 1999 年末的 restart，而不是中途某年
if [ "${CDATE%%-*}" != "1999" ]; then
    echo "！！警告：日期戳年份是 ${CDATE%%-*}，而共享 spinup 应当结束于 1999。"
    echo "   若 spinup 尚未跑完或配置有变，请先确认后再执行。"
    [ $DRYRUN -eq 0 ] && { echo "已中止。"; exit 1; }
fi

for exp in $EXPS; do
    DST="${CASE_ROOT}/restart/${exp}/${CDATE}"
    echo "--> $exp"
    if [ $DRYRUN -eq 1 ]; then
        echo "    mkdir -p $DST"
        for f in "$SRC"/${SPINUP_CASE}_restart*.nc; do
            base=$(basename "$f")
            echo "    cp $base  ->  ${base/${SPINUP_CASE}_/${exp}_}"
        done
    else
        mkdir -p "$DST"
        for f in "$SRC"/${SPINUP_CASE}_restart*.nc; do
            base=$(basename "$f")
            cp -f "$f" "${DST}/${base/${SPINUP_CASE}_/${exp}_}"
        done
        echo "    完成，$(ls -1 "$DST" | wc -l | tr -d ' ') 个文件"
    fi
done

echo ""
if [ $DRYRUN -eq 1 ]; then
    echo "以上为预览。确认无误后执行：  ./distribute_restart.sh --run"
else
    echo "分发完成。8 组实验现在共享同一份 1999-12-31 陆面状态。"
    echo ""
    echo "下一步：提交 8 组作业"
    echo "  for e in $EXPS; do sbatch submit_\${e}.slurm; done"
fi
