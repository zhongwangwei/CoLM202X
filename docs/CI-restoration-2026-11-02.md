# 官方 CI 恢复与验证记录（2026-11-02）

## 来源与交付边界

实时查询 PR https://github.com/CoLM-SYSU/CoLM/pull/516 （GitHub connector 与 gh API）：base 为 master，SHA `ac659a3d09cf378b3e90f7b34880ceb78c22a40d`；head 为 zhongwangwei/CoLM202X:master，SHA `c545a96c0aa94f5323861bd3ad7929b7dcfd27f2`。本地从该 head 开始，初始工作区干净。未推送、未合并。拟推送目标为个人仓库 master，须用户确认。

## 完整目录比较、用途与引用

`git diff --name-status ac659a3d -- .github/workflows` 最初仅列出以下四项；上游目录本身只有前三个文件。

| 文件 | 原差异与用途 | 恢复结果、引用及覆盖影响 |
|---|---|---|
| build_CoLM_gnu.yml | GNU/MPI 构建入口；个人版新增私有 NetCDF、缓存、自检、MPI I/O 探针、四分片、显式编译器传递、历史字段检查 | 恢复上游全部内容，保留 ubuntu-latest/ubuntu-22.04。引用 actions/checkout@v4、mpi4py/setup-mpi@v1、fortran-lang/setup-fortran@v1、apt 系统包及下列两个文件；不再执行上述个人版检查 |
| TestCaseLists | 原88配置，个人版新增5个 MT/LULCC 配置 | 恢复88配置；上游 workflow 读取，不再由官方矩阵覆盖新增5配置 |
| create_defineh.bash | 七参数宏生成器；个人版新增 MT/LULCC 参数并允许 PC+BGC | 恢复原脚本；测试调用改为 tests 专用辅助器，防止附加参数被静默忽略 |
| hist_acc_restart_fields.py | 新增历史累积字段检查/生成器 | 原样移动到 tests/history_restart/；不再自动由官方 CI 执行；MOD_HistAccRestart 两处说明同步修改 |

完整目录按上游 blob 逐字节核验；`git diff --exit-code ac659a3d -- .github/workflows` 通过。上游原有空白被保留，因此普通 `git diff --check` 会报告三处恢复引入的上游尾空白；目录外 diff-check 通过，不为消除尾空白改变官方文件。

目录外引用检索使用 `rg --hidden`，包含 material_naming/run_build.sh、material_host_cleanup/run_coupled.py，以及本地被忽略的 tests/run_grid_transport_build.sh；三者均转用 tests/create_material_defineh.py。新增辅助器保留七参数生成，再明确补充 MT/LULCC、PC+BGC，并用 C 预处理器核验有效宏；无效配置或上游 guard 改动立即失败。run/scripts/SummaryTest.bash 中 TestCaseLists 是传入路径变量，并非硬编码 CI 文件。历史字段脚本仅移动，不删除功能。

## 构建修复与历史澄清

- include/Makeoptions 是指向 Makeoptions.github 的受版本控制符号链接，实际修改文件为 include/Makeoptions.github。
- 原文件 CaMa 默认 FC/FCMP 为 /usr/bin/gfortran，来自2024年上游；此前个人 workflow 已通过 make 命令行覆盖，因此不能说此前个人 CI 中 CaMa 一直使用系统编译器。恢复原 workflow 的无参数 make all 后才重新暴露默认值问题。
- 构建配置现在保留环境/命令行 FC，默认采用 PATH 的 gfortran，CaMa FCMP 跟随 FC；OpenMPI、MPICH、IntelMPI 后端默认跟随 FC，已有显式 backend 覆盖保留。模型 FF 与其他用户构建配置不整体覆盖。
- 版本判断改为读取实际 FF 的 Fortran 主版本，并按数字比较，避免读取另一套 gcc；版本查询显式传入同一 MPI backend 环境。保留命令行 GCC_VERSION 覆盖。
- NetCDF 默认 /usr/include、/usr/lib/x86_64-linux-gnu 与链接库保持原样；没有编辑 mpi.mod，没有引入私有 NetCDF 路径。
- c545a96c 的三个甲烷文件修复完整保留，经 MOD_SPMD_Task 的 mpif.h 接口使用 MPI。其他 USE MPI 仅见 CaMa 的条件 MPI 模块与独立 ICI 模块；不能把未激活的路径当作本次失败原因。
- 除历史脚本路径注释外，模型源代码未修改，模型对象/依赖未回退。

## 旧运行诊断（不是本次恢复后的验证）

上游 https://github.com/CoLM-SYSU/CoLM/actions/runs/41314111396 有55成功、37失败。36个 latest 作业失败于私有 NetCDF-Fortran 构建检查；代表 job124870426870 的 nf_test4/ftst_vars 与 f03tst_vars 失败，但没有展开 test-suite.log，依赖根因尚未证明。实际 latest 镜像为 Ubuntu24.04.5 / ubuntu-24.04 / 20261025.353.1。

Ubuntu22.04/mpich/GNU13/shard3 的 job124870427174 约98分钟失去通信，check-run annotation 为 “The hosted runner lost communication with the server.”，不是已证明的编译错误或六小时超时。

个人 https://github.com/zhongwangwei/CoLM202X/actions/runs/41314106116 同head的旧92项全部成功；代表作业使用私有 NetCDF cache，不能证明全新依赖构建已修复，也不能证明恢复后的官方 CI 已通过。最长分片约112.47分钟，四分片去除后存在接近/超过六小时的风险；需要新 hosted 全程实测，不能按分片简单相加宣称必然超时。

## 本次验证

- 上游三文件字节一致、YAML解析、矩阵展开通过：23个组合。latest GNU12/13/14 × mpich/openmpi/intelmpi；22.04 GNU9–13 ×三MPI，排除GNU9+mpich。
- 88条配置逐条执行原 create_defineh、cpp、make -n all：全部通过；这是配置/依赖解析，不是88次编译。
- `bash -n .github/workflows/create_defineh.bash tests/material_naming/run_build.sh tests/run_grid_transport_build.sh` 通过。
- `python3 tests/create_material_defineh.py --self-test`：9种有效宏配置及无效参数拒绝通过。
- `python3 tests/check_makeoptions.py`：模拟GNU9/14、环境/命令行FC、MPI backend 与 CaMa 传递通过；模拟器不是实际GNU9/14编译。
- `python3 tests/history_restart/hist_acc_restart_fields.py --check`：525字段匹配。
- material_naming/check_inventory.py、check_inventory_exceptions.py 通过；check_history_io.py 的单rank/四rank共18项真实I/O检查通过。
- material_host_cleanup/run_coupled.py --self-test 通过，仅检查辅助函数，不是耦合科学验证。
- Python语法编译及目录外 git diff-check 通过。

远程重计算经 run-7920，以本地代码创建私有快照，无外部数据上传：

1. `codex-55cff4b0ae0444029afbd33c9a7b56b0`：succeeded，`bash tests/material_naming/run_build.sh pc`，MT+PC+BGC，colm/mkinidata/mksrfdata 三程序编译链接通过。该早期快照的宏检查使用 cpp；Linux执行成功，最终辅助器改用跨平台 cc 预处理，并在后续任务重验。构建选项使用测试的显式debug flags，不代表官方默认 flags。
2. `codex-8b3f76045ccb40be972f0c93826b39f8`：最终构建配置快照。GNU15.2 / Conda MPICH / Conda NetCDF，Ubuntu26.04.1；运行宏自检、host phase bookkeeping 和原七参数 GRID/PFT/vanGenu/CaMaON/BGCON/CROPON 的 make all。模型与后处理程序编译链接成功；make all 在静态库打包失败，原因是快照带入本地旧 lib/*.o 悬空生成链接，尚不能归为代码缺陷。host phase bookkeeping 的6组各64个案例及5个CaMa guard案例通过。私有快照清理生成链接后复验，见补充。

上述远程环境不等同 Ubuntu24.04 hosted runner 或 Ubuntu22.04，且没有系统 NetCDF-Fortran 的验证证据。恢复后的23个 hosted组合均尚未运行，GNU9–14各MPI真实编译仍待确认后推送个人仓库触发；没有把已有缓存依赖测试、模型编译和科学正确性混为一谈。

没有运行完整带真实数据的 continuous/restart 耦合科学验证；本次不宣称科学正确性得到新增证明。推送需用户确认，不合并官方PR。

最终补充：`codex-95b948f874a943acae5725362aebc768` 终态 succeeded，复用上述已编译对象，在私有快照清除旧 lib/*.o 生成链接及 macOS 旧 lib/libcolm.a 后，原配置默认 flags 的 `make all` 完整通过（含静态库），不修改 Makefile。未将此增量收尾称为全新编译；新编译证据来自8b3f任务。结果及二进制摘要已取回 `/tmp/colm-ci-check/cama-final/`；初次编译证据在 `/tmp/colm-ci-check/cama-build/`，MT编译在 `/tmp/colm-ci-check/pc-build/`。

复验中间记录：`codex-ded9041bd6434f3583fc33a5e9fba169` 在命令启动前因同步期间报告变动而被快照一致性校验拒绝；`codex-fa7ded63c6c44e28a2ff6f1864259bbe` 清理悬空链接后仍因带入的macOS旧归档导致GNU ranlib报malformed archive；最终任务同时清理两类生成产物后成功。上述清理仅限私有快照，不改变用户本地旧产物或官方CI。
