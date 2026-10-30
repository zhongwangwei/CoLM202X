#!/bin/bash
# ==============================================================================
#  生成 LSF 作业脚本：1 份 spinup + 8 组实验
# ------------------------------------------------------------------------------
#  运行环境栈直接复用 run/Amazon_iso_0p5_2020_run.lsf —— 那是本机（/share/home/dq013）
#  上已经跑通的配置，不要替换成 miniconda 里的 HDF5（其 HDF5 链接的是另一套 MPI，
#  会出现 MPI_Info_f2c undefined）。
#
#  用法：  ./gen_jobs.sh              # 用默认核数
#          NCORE=480 ./gen_jobs.sh    # 自定义核数
# ==============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# ---- 可调参数 ----
NCORE=${NCORE:-480}            # 总 MPI rank 数（全球 0.5 度陆面 + 0.25 度汇流）
PTILE=${PTILE:-48}             # 每节点核数
QUEUE=${QUEUE:-normal}
MEMREQ=${MEMREQ:-200G}
COLM_ROOT=/share/home/dq013/zhwei/colm/script/CoLM202X
NMLDIR=${COLM_ROOT}/run/aas_riverlake

emit () {   # $1=作业名  $2=nml 文件名  $3=是否含 mksrfdata/mkinidata (yes/no)
    local job="$1" nml="$2" prep="$3"
    local out="submit_${job}.lsf"

    cat > "$out" <<EOF
#!/bin/bash
#BSUB -J ${job}
#BSUB -q ${QUEUE}
#BSUB -o ${job}-o.out
#BSUB -e ${job}-e.out
#BSUB -n ${NCORE}
#BSUB -R rusage[mem=${MEMREQ}]
#BSUB -R span[ptile=${PTILE}]

# ------------------------------------------------------------------------------
#  AAS 河湖流实验 —— ${job}
#  自动生成于 gen_jobs.sh，请勿手工编辑；改参数后重跑 gen_jobs.sh。
# ------------------------------------------------------------------------------
export OMP_NUM_THREADS=1
ulimit -s unlimited

# 复用 Amazon_iso_0p5_2020_run.lsf 里已验证的运行时栈。
# 不要加入 miniconda/lib：其 HDF5 引用了另一套 MPI，会导致 MPI_Info_f2c undefined。
MV2=/share/home/dq013/software/mvapich2-2.3.7-1-intel
export PATH=\$MV2/bin:\$PATH
export LD_LIBRARY_PATH=/share/home/dq013/software/netcdf_parallel_N/lib:/share/home/dq013/software/hdf5-1.14.3-Intel_MPI_N/lib:/share/home/dq013/software/mvapich2-2.3.7-1-intel/lib:/share/home/dq013/software/zlib-1.2.11/lib:/share/home/dq013/software/szip-2.1.1/lib:/share/home/dq013/software/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin:/usr/lib64:\${LD_LIBRARY_PATH:-}

NML="${NMLDIR}/${nml}"
cd ${COLM_ROOT}/run

echo "job start: \$(date)  host: \$(hostname)"
echo "nml: \${NML}"
EOF

    if [ "$prep" = "yes" ]; then
        cat >> "$out" <<EOF

echo ">>> mksrfdata.x  \$(date)"
\$MV2/bin/mpirun -np ${NCORE} ./mksrfdata.x "\${NML}" &> ${job}_mksrf.log
RC=\$?; echo "mksrfdata rc=\${RC} \$(date)"
[ \${RC} -ne 0 ] && { echo "mksrfdata FAILED"; tail -120 ${job}_mksrf.log; exit \${RC}; }

echo ">>> mkinidata.x  \$(date)"
\$MV2/bin/mpirun -np ${NCORE} ./mkinidata.x "\${NML}" &> ${job}_mkini.log
RC=\$?; echo "mkinidata rc=\${RC} \$(date)"
[ \${RC} -ne 0 ] && { echo "mkinidata FAILED"; tail -120 ${job}_mkini.log; exit \${RC}; }
EOF
    fi

    cat >> "$out" <<EOF

# rank 0 单独启动作 master/IO，其余为 worker（与已验证的生产作业写法一致）
echo ">>> colm.x (${NCORE} ranks)  \$(date)"
\$MV2/bin/mpirun -np 1 ./colm.x "\${NML}" : -np \$((${NCORE}-1)) ./colm.x "\${NML}" &> ${job}_colm.log
RC=\$?; echo "colm rc=\${RC} \$(date)"
[ \${RC} -ne 0 ] && { echo "colm FAILED"; tail -120 ${job}_colm.log; exit \${RC}; }

echo "${job} COMPLETE: \$(date)"
exit 0
EOF

    chmod +x "$out"
    echo "生成 $out"
}

# 阶段 0：地表数据 + 初始场 + 共享 spinup（含 mksrfdata/mkinidata）
emit SPINUP aas_spinup.nml yes

# 阶段 1+2：8 组实验（复用同一套地表数据，只跑 colm.x）
for e in E0_CTL E1_D E2_L E3_B E4_DL E5_DB E6_LB E7_ALL; do
    emit "$e" "aas_${e}.nml" no
done

cat <<'TXT'

==============================================================================
 提交顺序（必须按序，不能并行）
------------------------------------------------------------------------------
 1) bsub < submit_SPINUP.lsf
      地表数据 + 初始场 + 1996-1999 共享 spinup。8 组复用其产物。

 2) ./distribute_restart.sh          # 先预览
    ./distribute_restart.sh --run    # 确认后执行
      把 1999 年末 restart 分发给 8 组。

 3) 先单独提交 E7_ALL 做单年标定（把 nml 的 end_year 临时改为 2000）：
      bsub < submit_E7_ALL.lsf
      ALL 是最慢的一组（BIF 强制全局 dt 同步，MOD_Grid_RiverLakeFlow.F90:974），
      用它标定得到的是墙钟时间上界。同时采集 Section 4.5 需要的性能数据。

 4) 标定通过后恢复 end_year=2005，提交全部 8 组：
      for e in E0_CTL E1_D E2_L E3_B E4_DL E5_DB E6_LB E7_ALL; do
          bsub < submit_${e}.lsf
      done
==============================================================================
TXT
