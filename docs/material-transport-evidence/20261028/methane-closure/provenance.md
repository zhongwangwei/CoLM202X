# 证据与来源说明（provenance）

任务：task-2「甲烷候选输入闭合、参数冻结与水稻/两年全球运行设计（只读）」。
本任务**只做只读调查**：未启动模型、未提交/取消任何 7920 任务（除本任务自己的两个只读 survey 作业，其一因超时被自己取消，另一成功）、未上传/下载任何大数据、未修改仓库、未 git commit/push、未改动任何他人 7920 作业。

## 1. 代码基线与仓库状态

- 分析基线：`/Users/zhongwangwei/Desktop/Github/CoLM-constituent`，分支 `feature/constituent-transport`，HEAD 在任务开始时为 **`10e9fbd204a99ff81bd60a98198b364d489e10d9`**。
- 分析期间其他 agent 在该仓库新增了两个**纯文档/证据**提交：`2d5c15b7`（p9r isotope/CL 回归记录）、`74cb9175`（p9r 源码版本核对）。`git diff --name-only 10e9fbd2..74cb9175` 不含任何 `.F90/.h/.nml` 文件 → 本任务所有源码结论仍对 `10e9fbd2` 成立。
- `/Users/zhongwangwei/Desktop/Github/CoLM202X` 工作树有他人既定未提交修改（`main/*.F90` 等），与本任务无关；本任务只在其 `.tmp/ct/ch4_closure10/` 下写文件。
- 远程代码：BD `script/CoLM202X` 当前 HEAD `080a09fb`（2026-09-17），**不是** 2026-05-27 旧运行的代码；旧运行使用的工作树 define.h 未提交、不可恢复，见 GAP_ANALYSIS §2.4。

## 2. 本地证据身份核对（下载/复制前先比哈希）

| 证据 | 本地路径 | sha256 | 与 7920 原件的核对 |
| --- | --- | --- | --- |
| 甲烷输入审计全文 | `evidence/audit.json` | `1d6fb3cbdbacc3d68e7af0e1d1a69d945db003a0b237bfed3b995cd4d98bdca1` | 与 `…/codex-ad6b48e0cb0044a9bbf0ce2903195234/results/audit.json` 相同 |
| 边界请求台账 | `evidence/local-boundary-requests.csv` | `d9e3d523b1c1c3fee34659a9569825104a161d3c475b41fb5b46ae52bf785b7b` | 与远端 `boundary-requests.csv` 相同 |
| 审计脚本 | `audit.py` / `probe_boundaries.py` | `689f3fbe…` / `53ddd41f…` | 与远端结果目录相同 |
| 审计源码快照 | `audit-source.sha256` | — | 5 个 F90 与 CoLM-constituent HEAD `10e9fbd2` 逐字节相同（本轮再次核对） |

`audit-source.sha256` 中的 `MOD_Forcing.F90/MOD_Precision.F90/MOD_TimeManager.F90/MOD_UserDefFun.F90/MOD_UserSpecifiedForcing.F90` 与 HEAD 对比：5/5 相同 → 15 文件请求集在 HEAD 上仍然有效。

## 3. 远程只读调查（7920）与作业记录

- 主 survey 作业：**`codex-fa5a1a057cc848e3a2468dd9797a73c0`**（`run-7920 submit --agent codex --source-local <scan> --cpus 8 --memory-gb 16 --timeout 1500`，**succeeded/0**）。输入以 `--input-remote` 只读声明：候选 staging、`CoLMruntime`、`CoLMrawdata`。产出：`scan.json`、`probe_windows.json`、`probe_window.log`、`source.sha256`（本地副本在 `evidence/`，`scan/results/metadata/job-result.json` 保留调度元数据）。
- 第一次 survey 作业 `codex-7ebe46ca0d4e48ee9ffc1312113f0462` 因全盘 find 过慢且脚本有 pathlib bug，由本任务自行 `run-7920 cancel`（终态 cancelled/-15），**不属于证据**；其脚本已修正后重跑。
- 只读 SSH 读取（无写操作）：7920 上的 staging 目录、`CoLMruntime`、`CoLMrawdata`、`data-stage` 列表、候选 history 头；以及经 7920l→`bd`（BD 集群，只读 `ls`/`cat`/`tail`/`grep`）确认 1994/1997 CRUJRA 文件与旧运行日志/参数。
- 未使用 `run-7920 fetch` 拉取任何模型结果；只取回本任务的 200KB 级 survey JSON。
- `scan/probe_window.py` 在 7920 上用系统 `gfortran 15.2` 编译「只提取生产 `setstampLB/setstampUB` 例程 + 桩」的探针；源码 4 个 F90 的 sha256 与审计快照相同（见 `scan/sources/`、`evidence/source.sha256`）。探针 W0 重放出审计的 15 条外部请求（23,407 行、35,088 步），证明探针与审计一致。

## 4. 数据来源（未复制、未上传）

- 候选包来源与字节数见审计 manifest 与 GAP_ANALYSIS §2.1；盘上复核 40,450 文件 / 26,618,457,915 B 完全一致。
- 15 个缺失文件的权威来源在 BD：`bd:/share/home/dq013/zhwei/colm/data/CoLM_Forcing/crujra_v2.5/<var>/crujra.v2.5.5d.<var>.{1994,1997}.365d.noc.nc`，只读 `ls -l` 实测存在且尺寸与 1995 同变量文件一致。**本轮没有传输或复制这些文件。**
- 旧运行日志：BD `…/cases/globe/g_global_ch4_flat_2yr_20260527_200428/` 下 `logs/lsf.168594.out`（11 KB，完整读取）、`mkinidata_global_current.log`（47 KB，只读 head/tail）、`colm_global_current.log`（1.28 GB，仅只读 `head/tail/grep -m`，**未复制**）。基准数字见 `evidence/bd_old_run_benchmark.json`。
- 历史参数/构建线索：BD git `show 39a8b691:include/define.h`、`git show 39a8b691:run/standard_ch4_parameter.nml` 的只读 grep/快照（`evidence/bd_historical_source.txt`、`bd_define_h_sep2026_snapshot.h`、`bd_standard_ch4_parameter_historical.nml`）。
- 7920 上 `training2026/CoLM_Forcing/crujra_v2.5` 是**悬空 symlink**（指向 BD 未挂载路径），不代表 7920 有 CRUJRA 数据。

## 5. 本地复现脚本与工具

- 环境：macOS 27.0.1、Python 3.12.12、`netCDF4`/`pandas`（本机分析）；7920 python = conda `/home/zhwei/miniconda3/bin/python3`（netCDF4 1.7.2 / numpy 2.3.5），gfortran 15.2.0。
- 脚本：
  - `scan/remote_scan.py`：7920 只读清点/搜索（sections A–F，逐节增量写 scan.json）
  - `scan/probe_window.py`：生产时间索引例程的窗口重放
  - `evidence/gen_closure.py`：从 audit.json/compat/scan/probe 生成 `closure.json` 与 `boundary_request_table.md`
- 命令入口：`~/.local/bin/run-7920`（README 规则：>2 min 计算走 7920；本任务无长跑）。

## 6. 未覆盖 / 不可得

- 7920 上没有旧运行日志或其它 CoLM 全球甲烷基准；旧运行在 BD，未经用户授权不传输其大文件。
- 旧运行的实际编译宏（define.h）未提交，无法从源码管理恢复；只能从 landdata 内容与日志推断部分宏（CROP 参与 landdata 生成；甲烷 active；hybrid；spatial pH 未生效；lake prod off；rice off）。
- 运行时依赖（soilstate/cnsteadystate/routing/crop）与 BD 历史算例的**内容身份**无哈希证据，标注 UNPROVEN。
- 候选输入包的唯一性、版本正确性与科学 ground truth 均未建立；旧运行完成 ≠ 输入闭合 ≠ 科学正确。
