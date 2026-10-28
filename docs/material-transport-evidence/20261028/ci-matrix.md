# 完整 hosted CI 编译器/MPI 矩阵（个人 fork workflow_dispatch）

状态记录：2026-10-29 07:10 CST。本文件记录两次 run：旧 run（`10e9fbd2`，**20/23 通过、3 个 6h 超时取消，
完整矩阵未通过**）与新 run（`2bc7510a`，CI 分片修复后触发，**进行中，92 个 job**）。

## 触发方式与身份

| 项 | 值 |
| --- | --- |
| 仓库 | 个人 fork `zhongwangwei/CoLM202X`（public；admin 权限） |
| 分支 / ref | `feature/constituent-transport` |
| 触发方式 | fork `workflow_dispatch`（**不是**向上游开 PR、不是合并 master/colm30、未改触发过滤） |
| 证据目录 | `docs/material-transport-evidence/20261028/ci-matrix/`（REST/web 快照、逐 job 汇总、取消 job 摘录） |

## Run 1：`40725517808`（head `10e9fbd2`）— **未通过**

| 项 | 值 |
| --- | --- |
| run | https://github.com/zhongwangwei/CoLM202X/actions/runs/40725517808 |
| head | `10e9fbd204a99ff81bd60a98198b364d489e10d9`（event=workflow_dispatch，attempt 1） |
| 官方 REST 终态 | `status=completed`、**`conclusion=cancelled`**、`run_started_at=2026-10-28T16:23:40Z`、`updated_at=2026-10-28T22:23:59Z` |
| job 统计 | **23 job = 20 success + 3 cancelled + 0 failure**（REST 快照 `jobs-final.json`、`run.json`） |
| 耗时 | 中位 99.0 min；最长 360.2 min（= 6h 硬上限） |
| 3 个 cancelled | `ubuntu-latest × mpich × gcc{12,13,14}`，均在 22:23:56–57Z 触发单 job 6h 上限，93 例逐例构建步骤未完成 |
| 结论口径 | **完整矩阵未通过**：20/23 成功、3 个超时未完成，必须原样记为未完成；不能写成“完整通过”。 |

CI 内容与覆盖核对（已完成）：

- 实际矩阵 = **23 job 而不是 25**：YAML 2 OS × 3 MPI × 6 GCC = 36，exclude 去重 13 后 = 23；YAML 穷举名单与 API 名单完全一致；
  历史上游 run `40682560947`（head `1f4cc9da`，已 CLOSED 的 PR #514 触发）同为 23。
- 23/23 job 的 `Check the history-accumulator restart field list` 与 `Check effective test case macros` 步骤均 success；
  `check_defineh.py` 输出 `PASS: optional arguments and effective macros for 93 unique CI cases`；
  2 个 job 有日志行级 marker 复核（`log-markers/`，含完整日志缓存）。
- 本地独立核对（`ci-coverage-check.txt`）：93 唯一用例、MTON=5、LULCCON=2（IGBP/PFT）、BGC+LULCC=0、
  PC+BGC+MT 与 IGBP/PFT+LULCC+MT 断言成立。

**根因**：不是宏/语法错误，而是 `ubuntu-latest`(24.04)+mpich 组合下 93 例串行 `make clean && make all` 的累计耗时超过
GitHub 单 job 6h 硬上限（ubuntu-22.04 mpich 组 3h37m–5h52m 都能通过；openmpi/intelmpi 组约 0.9–1.8h）。

## CI 结构修复（提交 `2bc7510a`）与 Run 2：`40784782071`（head `2bc7510a`）— 进行中

- 修复：`.github/workflows/build_CoLM_gnu.yml` 增加 `matrix.shard: [0,1,2,3]`，把 93 例构建**分片**到 4 个 job
  （`index % 4 == shard`），每个 os/mpi/toolchain 组合仍完整覆盖 93 例；`check_defineh.py` 仍在每个 job 全量执行。
  **不减少覆盖、不改触发过滤**，只把每个 job 的工作量压到 6h 以内。矩阵变为 23×4 = **92 job**。
- Run 2：`40784782071`，head `2bc7510a2f9bddcf294b4a8999f98540086d4316`，event=workflow_dispatch，
  2026-10-28T23:05:05Z 创建；首次核对为 92 job（40 in_progress / 52 queued）。
- **Run 2 终态尚未取得**：在 92/92 job 终态之前，仍不能声称完整矩阵通过。终态、逐 job 结论与失败/取消摘录
  写入 `ci-matrix/`（由 `ci-matrix` 子任务持续跟进；`run-fixed` 证据目录）。

## 残余缺口（如实保留）

1. Run 1 的 3 个 cancelled job 无完整日志（cancelled 步骤未产出）；已保存步骤级记录与官方终态，缺口为“这 3 个
   组合在该 run 中未完成构建”，不得用其它 20 个 job 推断其通过。
2. Run 2 需要 92 job 全部终态；若仍有 job 因 6h/基础设施失败，按“未通过/未完成”记录并考虑进一步分片或拆分 job。
3. CI 只验证构建与配置（93 例宏/编译），不运行模型、不覆盖真实数据物理验收。
4. 新 run 的 dispatch 与本机网络：本机代理故障期间使用本机 SSH SOCKS 中继（`ssh -f -N -D 11080 7920l`）恢复
   GitHub 访问；未修改系统代理、未把 token 传到任何远端。

## 证据

- `ci-matrix/run.json`、`ci-matrix/jobs-final.json`：Run 1 官方 REST 终态快照（23 job：20 success + 3 cancelled）。
- `ci-matrix/job-summary.tsv`、`job-stats.txt`：逐 job 结论与耗时。
- `ci-matrix/EVIDENCE-run-40725517808.md`：子任务终态报告（含 web 页面证据、逐 job 表、缺口与自动补齐机制）。
- `ci-matrix/web-terminal.json`、`job-summary-web.tsv`、`job-stats-web.txt`：web 页面终态推断与对照。
- `ci-matrix/failure-excerpts/`：cancelled/失败 job 摘录（含 `PENDING-cancelled-jobs.md`）。
- `ci-coverage-check.txt`：独立矩阵展开与 `check_defineh.py` 覆盖核对。
