# 分片 CI 矩阵 — 新 run（40784782071，head 2bc7510a）终态证据（run-fixed）

状态：**run 终态 = completed / success；92/92 job 全部 success（0 failure、0 cancelled、0 timeout）；92/92 job 日志行级 marker 核对通过。分片矩阵通过。**
旧 run（40725517808，head 10e9fbd2）证据见上级目录 `../EVIDENCE.md`：该 run 为 **20/23、3 个 6h 超时 cancelled，未通过**（历史记录保留，不被本 run 改写）。两套证据分开记录。
观测时点：2026-10-29 CST（UTC+8）；表格内为 UTC。执行者：teammate `ci-matrix`（Team task-3）。

## 1. 触发记录

| 项目 | 值 |
| --- | --- |
| run id | **40784782071** |
| run URL | https://github.com/zhongwangwei/CoLM202X/actions/runs/40784782071 |
| event / attempt | `workflow_dispatch` / 1 |
| created_at / updated_at | 2026-10-28T23:05:05Z / 2026-10-29T01:06:43Z |
| **run 结论** | **completed / success**（官方 REST） |
| head_branch / head_sha | `feature/constituent-transport` / `2bc7510a2f9bddcf294b4a8999f98540086d4316` |
| workflow | id 78251663 `Build_CoLM202X_gnu`，path `.github/workflows/build_CoLM_gnu.yml` |
| job 总数 | **92**（= 23 matrix cell × 4 shard） |
| 实际墙钟耗时 | **2h01m38s**（7298s） |

由 Lead 触发；本会话用 REST API 观测与取证（SOCKS 中继），未 push/改仓库/动历史 run。

## 2. 源码身份（dispatch ref 上 CI 文件，已冻结到 `source-snapshots/`）

| 文件 | git blob SHA（commit 2bc7510a） | 与旧 run（10e9fbd2）相比 |
| --- | --- | --- |
| `.github/workflows/build_CoLM_gnu.yml` | `079eefc054aa0c750185b25c185509c691a73fed`（5050 B） | **已变更：新增 `matrix.shard: [0,1,2,3]` + 分片构建循环** |
| `.github/workflows/TestCaseLists` | `612b55722dacb63b58f80b27dac495f8f239604e`（93 行） | 未变 |
| `tests/material_transport/check_defineh.py` | `dc486706d83bf8663954b3604dd8369d0e03d17c` | 未变 |
| `.github/workflows/hist_acc_restart_fields.py` | `0b25fd85afcac7088693b1adaf57fc48fed56988` | 未变 |
| `.github/workflows/create_defineh.bash` | `d70d48cbee0d67dd50b49f542e0923f4d8bc7e68` | 未变 |

blob SHA 由 `verify_sharded_matrix.py` 以 git 对象格式独立复算（见 `coverage-check.txt`），与 contents API 返回一致。

## 3. 分片矩阵展开与覆盖（独立核对）

- 维度：`os`(2) × `mpi`(3) × `gcc`(6) × `shard`(4) = **144** 笛卡尔组合；5 条 `exclude` 在 4 个 shard 上同样生效，去重后排除 52 个 → **92 job**。
- YAML 独立展开的 job 名（含第 5 维 shard）与 Actions API 返回的 92 个 job 名集合**完全相等**。
- TestCaseLists 分片规则（workflow 内嵌 python，`index % 4 == shard`）：**shard0=24 例、shard1/2/3=23 例**；4 个 shard 并集 = 完整 93 例、无重复、无遗漏。
- 每 shard 期望宏命中（`per-shard-expectations.json`，用于日志 marker 校验）：
  - MTON：shard0=2、shard1=1、shard2=1、shard3=1（合计 5）；
  - LULCCON：shard0=1、shard1=0、shard2=0、shard3=1（合计 2）。
- workflow 每个 job 均执行：`Check the history-accumulator restart field list`、`Check effective test case macros`（check_defineh.py 对**全量 93 例清单**校验，不是分片子集）、`Build CoLM202X`（`set -euo pipefail` + shard 循环 + 每例 `make clean && make all`）。
- 分片目的（workflow 注释）：93 例编译拆到 4 个 job，使每个 job 低于 GitHub 6 小时硬上限，**不减少覆盖**。

## 4. CI 检查内容与覆盖声明（与旧 run 同一套，未变）

`check_defineh.py`、`TestCaseLists`、`hist_acc_restart_fields.py`、`create_defineh.bash` blob 均未变，因此旧 run 的本地独立核对继续有效：

- 本地 `python3 tests/material_transport/check_defineh.py` → `PASS: optional arguments and effective macros for 93 unique CI cases`；
- 本地 `python3 .github/workflows/hist_acc_restart_fields.py --check` → `main/MOD_HistAccRestart.F90: 525 fields match main/MOD_Vars_1DAccFluxes.F90`；
- 93 唯一用例、参数 7/8/9 个 = 88/3/2、MTON=5、LULCCON=2（IGBP/PFT）、BGC+LULCC=0；PC+BGC+MT 与 IGBP/PFT+LULCC+MT 断言成立，非法第 9 参数被拒。

## 5. 运行结果（官方 REST 终态）

| 项目 | 值 |
| --- | --- |
| status / conclusion | **completed / success** |
| job 总数 | **92**（全部 completed） |
| **success / failure / cancelled / timed_out** | **92 / 0 / 0 / 0** |
| job 耗时 | min **11.2 min**；median **24.4 min**；max **101.9 min**（无 job 接近 6h 上限） |
| 23 个 matrix cell | **每个 cell 的 4 个 shard 全部 success（23/23 cell 完整、4/4 shard 各成功）** |
| 关键步骤 | 92/92 job 的 `Check the history-accumulator restart field list`、`Check effective test case macros`、`Build CoLM202X` 均 success |

- 每 cell 的 4 个 shard 合计构建该 cell 的完整 93 例，且每个 job 内 check_defineh 对全量 93 例清单校验；23 cell × 93 例 = 全部 2,139 个「cell-shard 用例构建」均成功。
- 按 gcc 版本分布的 job 数：gcc9=8、gcc10=12、gcc11=12、gcc12=24、gcc13=24、gcc14=12（与 exclude 后矩阵一致）。
- MPI 生效路径：60 个 job 用 `/usr/bin/mpif90`（mpich/openmpi），32 个 job 用 `/opt/intel/oneapi/mpi/2021.19/bin/mpif90`（intelmpi，8 cell×4），均按矩阵生效。

## 6. 逐 job 日志行级核对（92/92 通过）

来源：REST job logs（92 个，缓存于 `cache-logs/*.log.gz`），解析结果 `log-markers.tsv` + `log-markers/<job_id>.txt`。

- **92/92 job** 的 `Check effective test case macros` PASS 行存在（93 unique CI cases）；
- **92/92 job** 的 hist 字段检查行存在（525 fields match）；
- 每 job 按 shard 期望：`Create test case … Complete!` / `Start Compilation …` 数量 = **24（shard0）/23（shard1,2,3）**，唯一名同数，`make clean && make all` 全数完成；
- 每 job 的 `#define MaterialTransport` 数量 = shard 期望（2/1/1/1），`#define LULCC` 数量 = shard 期望（1/0/0/1）；
- 每 job 出现 `Shard X/4: N test cases are compiled successfully!` 与 `All 93 cases are covered across the four shards of this matrix cell.`；
- error marker = 0（92/92）；无失败/取消 job，`failure-excerpts/` 为空（仅保留旧 run 的历史摘录在上级目录）。
- 抽样示例（`cold=ok` 全部为 True）：shard0 job cases=24 mt=2 lulcc=1；shard1/2 job cases=23 mt=1 lulcc=0；shard3 job cases=23 mt=1 lulcc=1。
- 汇总：`jobs_with_full_expected_markers=92/92 (processed=92 terminal=True) missing_logs=0`（`collection-report.txt`、`collected.ok`）。

## 7. 结论与口径

1. **分片 CI 矩阵通过：92/92 job success，run conclusion = success，墙钟 2h01m38s，无 job 接近 6h 上限（最大 101.9 min）。**
2. **覆盖等价**：23 个 os/mpi/gcc cell（与旧矩阵一致）每个都由 4 个 shard job 覆盖完整 93 例；`check_defineh.py` 仍对全量 93 例清单在每个 job 执行并通过；hist 字段检查在每个 job 执行并通过。
3. 与旧 run 的关系：旧 run 40725517808（commit 10e9fbd2）**20/23 success、3 个 ubuntu-latest×mpich×gcc{12,13,14} 因单 job 6h 上限 cancelled，仍未通过**（历史记录保留，不因本 run 改写）；本 run 在 **commit 2bc7510a** 上通过分片解决了该耗时瓶颈，3 个此前未完成的 cell 均由 4 个 shard 成功覆盖。
4. 范围说明：本 run 验证的是「分片后的编译矩阵」（每个 cell 的 93 例分散到 shard），不是单 job 内一次性跑完 93 例；这是 workflow 的显式设计（注释与日志均说明 all 93 cases covered across four shards），不减少用例/宏覆盖。CI 仅为编译覆盖，不替代模型运行/科学验收。
5. 无残余失败项；本 run 证据完整。唯一需注意：run 属 commit `2bc7510a`（仓库在该提交上新增分片），旧 run 结果仍按旧 commit 记录。

## 8. 证据文件索引

| 文件 | 内容 |
| --- | --- |
| `run.json` / `jobs.json` | 官方终态 run / jobs 快照（completed/success；92 jobs） |
| `jobs-final.json` | collector 写入的官方终态 jobs 快照（含每 job steps） |
| `job-summary.tsv` / `job-stats.txt` | 逐 job 官方结论/起止/耗时与统计（92/92 success；min 11.2m / median 24.4m / max 101.9m） |
| `cells-summary.json` | 23 个 cell × 4 shard 的结论与耗时汇总（0 缺失、0 失败） |
| `log-markers.tsv` / `log-markers/<job_id>.txt` | 92 job 日志行级 marker 计数 + shard/期望值列 + 摘录 |
| `cache-logs/*.log.gz` | 92 个 job 完整日志缓存（可离线复核） |
| `collection-report.txt` / `collected.ok` | 采集报告：92/92 全 marker 通过、missing_logs=0 |
| `coverage-check.txt` / `verify_sharded_matrix.py` / `per-shard-expectations.json` | 92 job 展开、分片覆盖 24/23/23/23、per-shard 宏期望、blob 身份 |
| `source-snapshots/` | commit 2bc7510a 的 5 个 CI 文件冻结快照 |
| `monitor-curl.log` / `snapshots/` | 轮询时间线与各时点快照（23:05Z→01:06Z） |
| `collect-partial.stdout` / `collection-report.partial.txt` | 终态前对已完成 shard job 的抽样核对（验证分片 marker 逻辑） |
| `../EVIDENCE.md` | 旧 run（40725517808，20/23 未通过）终态证据 |
