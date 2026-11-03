# MPICH调试构建耗时修复（2026-11-03）

基于个人仓库83da4d3e；PR516的base仍为ac659a3d09cf378b3e90f7b34880ceb78c22a40d。本次尚未推送，不合并上游PR。

## 原因与最小变更

上游运行41558745865有22项成功，Ubuntu22.04/MPICH/GNU12超过6小时被取消，已完成87条配置、正在最后一条。原workflow每条均make clean && make all，共88次串行全量构建。

该作业安装的Ubuntu包mpich4.0-3中，mpifort.mpich默认在用户参数前加入-O2、-flto=auto、-ffat-lto-objects；链接还加入-flto。模型GNU配置原来有-g、-C、边界/浮点异常检查，却没有显式优化/LTO设置，因此MPICH继承了额外优化，而其他MPI组合不一定如此。同轮首19个已完成配置的耗时中位数为MPICH322.9秒、OpenMPI72.3秒；这些跨runner观察不是严格对照。

本次仅修改include/Makeoptions.github两处默认FOPTS，显式加入-O0 -fno-lto，并扩展tests/check_makeoptions.py：保留GNU9/14选项分支、所有既有检查和显式FOPTS覆盖。CaMa自己的FFLAGS=-O3保持原样，其他用户Makeoptions不改。官方workflow、88条配置、模型对象、源码和依赖自检均不改，不增加分片或make并行度。

这是统一默认GNU调试构建的优化行为，不是科学代码修复。隐式-O2/LTO不再是此调试配置的验证覆盖；自定义发布/性能构建仍可显式给出完整FOPTS。默认调试程序的运行速度可能低于此前隐式优化的MPICH版本，不能把编译提速说成模拟运行提速或科学结果逐位一致。

## 对照验证

全部重计算使用run-7920私有快照；没有修改服务器系统或上传科研数据。GNU12.5、Ubuntu26.04.1、系统NetCDF-C4.9.3及已验证的匹配Fortran库。没有声称这是完整GitHub镜像复刻。

1. 大模块配对：任务codex-3839fef8cbc1497db53388cc500b27b2，终态succeeded。使用CI同版本原始mpifort.mpich脚本，实际编译器选择GNU12；使用已有匹配模块接口，编译而不链接MPI运行库。校验两份源文件SHA256与当前仓库一致，两个次序相反的重复中只有候选的-O0 -fno-lto不同。

| 源文件 | wrapper默认平均秒数 | 显式调试平均秒数 |
|---|---:|---:|
| main/LULCC/MOD_Lulcc_Vars_TimeVariables.F90 | 5.686 | 0.866 |
| mkinidata/MOD_Initialize.F90 | 4.301 | 1.119 |

2. 整配置配对：任务codex-ee8342c1dee243a6984d4025ad0d9b1a，终态succeeded。使用上述MPICH4.0脚本的相同编译/链接参数，测试环境仅将其MPI include/lib目录及库名称映射到已有Conda MPICH4.3.2，避免把MPI运行库版本差异带入A/B；这不是Ubuntu22.04原安装环境的完整复刻。两边同一GNU12、源码、MPI库、CaMa/NetCDF路径，先准备兼容库，再分别make clean。原始第一配置GRID/PFT/URBANON/Campbell/CaMaON/BGC/CROP均完整make all，包含三个主程序、后处理及静态库；基线显式覆盖旧FOPTS，候选使用仓库新默认值。

| 版本 | 完整make all | 边界/浮点异常运行检查 | 代表对象LTO段 |
|---|---:|---|---|
| wrapper原默认 | 422.024秒 | 正常触发 | 有 |
| -O0 -fno-lto | 132.122秒 | 正常触发 | 无 |

编译墙钟时间减少68.7%，此单次完整配置约3.19倍加速；不是88配置或GitHub全部runner的耗时承诺。微基准使用交错次序补充，不能完全排除整配置先后运行的文件缓存或外部负载影响。

测试分别编译并运行故意越界和除零程序，必须非零退出且给出相应诊断；不是仅检查选项字符串。两边模型ldd均确认系统libnetcdf.so.22。未运行真实数据科学模拟，没有输出逐位一致的结论。

3. 本地tests/check_makeoptions.py通过，新增断言先在旧默认下因缺少-O0失败，修改后通过；保留全部既有检查、用户FOPTS=-O3 -g覆盖。OpenMP目标、七项依赖助手测试、make准备流程检查通过。git diff --check以及与上游base整个.github/workflows字节一致检查通过。

完整日志、展开后的wrapper命令、对象ELF段和timings.json已取回/tmp/colm-opt-checks/module-bench及/tmp/colm-opt-checks/full-bench。准备好后须用户确认推送至zhongwangwei/CoLM202X:master，才能验证全部hosted矩阵；尚不能宣称六小时时限问题已在23个组合全部解决。
