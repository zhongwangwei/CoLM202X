# 完整 hosted CI 编译器/MPI 矩阵（个人 fork workflow_dispatch）

状态记录：2026-10-29 03:29 CST（最后核对）。**本文件记录的是“已触发 + 终态前”的真实状态，不是最终通过结论。**

## 触发方式与身份

| 项 | 值 |
| --- | --- |
| 仓库 | 个人 fork `zhongwangwei/CoLM202X`（public；admin 权限） |
| 分支 / ref | `feature/constituent-transport` |
| head SHA | `10e9fbd204a99ff81bd60a98198b364d489e10d9`（= 交付 HEAD；生产目录等于 `6d995438`） |
| 事件 | `workflow_dispatch`（fork 上手动触发，**不是**向 SYSU 开 PR、不是合并 master/colm30、未改触发过滤） |
| run | `40725517808` — https://github.com/zhongwangwei/CoLM202X/actions/runs/40725517808 |
| 触发时间 | 2026-10-28T16:23:40Z（2026-10-29 00:23:40 CST） |
| 实际矩阵 | **23 个 job**（YAML：2 OS × 3 MPI × 6 GCC = 36，exclude 去重后 23；任务简报的“25”与 YAML 不符，以实际为准） |
| 分支上的 workflow blob | `2390fad0…`（含 `workflow_dispatch`、`check_defineh.py` 步骤、`TestCaseLists` 93 行、逐例 `make clean && make all` 循环），与本地 HEAD 一致 |
| 旧上游 run 对照 | 上游 PR #514（已 CLOSED）的 run `40682560947` 用的是 `1f4cc9da`（早于 CI 用例补充与全部最新修复），**不作为**本轮 CI 验收；未触碰该 run |

## CI 内容核对（来自首个完成 job 的完整日志，逐项核对）

- `check_defineh.py` 输出 `PASS: optional arguments and effective macros for 93 unique CI cases`（出现 1 次）。
- 历史累加量字段检查：`main/MOD_HistAccRestart.F90: 525 fields match main/MOD_Vars_1DAccFluxes.F90`。
- 93 个唯一用例各 1 次 `Create test case … Complete!` / `Start Compilation`，逐例 `make clean && make all`，无缺失、无 error marker。
- `define.h` 中 `#define MaterialTransport` 出现 5 次（5 个 MTON 例）、`#define LULCC` 2 次（IGBP、PFT 的 MTON+LULCCON 例），与 `TestCaseLists` 一致。
- 本地 `check_defineh.py` 的覆盖声明与 `TestCaseLists` 实际内容一致（93 例、唯一名、参数个数 7/8/9、PC+BGC+MT、IGBP/PFT+LULCC+MT、BGC+LULCC 被排除、非法 LULCC 参数负例）。独立脚本 `verify_matrix.py` 输出见 `.tmp/ct/ci10/coverage-check.txt`。

## 当前状态（2026-10-29 03:29 CST 快照）

| 指标 | 值 |
| --- | --- |
| 已完成 | **16 / 23**，全部 **success**，0 失败、0 取消 |
| 仍需等待 | **7 个 mpich job** 全部 in progress（ubuntu-latest mpich gcc12/13/14；ubuntu-22.04 mpich gcc10/11/12/13） |
| 已耗时 | 约 3.0 小时（run 于 16:23:40Z 启动；7 个 mpich 已跑约 2.98 h） |
| 历史参照 | 上游同类 run（`40682560947`）的 mpich job 单例耗时 **3h06m–5h48m**，明显慢于 openmpi/intelmpi（约 1.5 h）；本 run 的 mpich 预计 19:30Z–22:15Z 陆续完成 |
| 单 job 上限 | GitHub Actions 默认 **6 h** → 最迟 22:23Z（**CST 06:23**）；若个别 job 超时被取消，按“未完成/缺口”保留，不算通过 |
| 结论口径 | **只有 23/23 全部 success 才算完整矩阵通过**；当前只能报告“16 个非 mpich job 全通过、mpich 待定”，不能写成完整通过 |

## 证据

- GitHub Actions run：https://github.com/zhongwangwei/CoLM202X/actions/runs/40725517808
- 本地私有（不进 Git、含 API 快照/日志标记/独立覆盖核对）：`/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/ci10/`
  - `run.json`、`jobs.json`、`coverage-check.txt`、`verify_matrix.py`、`monitor.log`、`source-snapshots/`
- 本目录归档：`ci-coverage-check.txt`（独立矩阵展开与 check_defineh 覆盖核对输出）

## 未完成 / 终止条件

- 终态若为 23/23 success：可在状态文档记录“完整 hosted 编译器/MPI 矩阵通过（commit、run URL、23 job）”。
- 终态若有失败/取消/超时：保留原始状态，逐 job 记录结论与关键错误，不做“部分成功=完整通过”的替换。
- 本轮 CI 不覆盖：真实模型运行（CI 只做构建与配置检查）、甲烷、数据依赖测试。
