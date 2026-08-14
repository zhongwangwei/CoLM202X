# River history 分片输出 — 第 1 步基线报告与 go/no-go

对应 `.omx/plans/river-history-sharded-output.md` 第 1 步。
测量机器：macOS arm64 笔记本（Apple Silicon），gfortran 16.1 + Open MPI 5.0.9，
本地 SSD。**所有 rank 均为超订阅（oversubscribed）**。

## 交付物

| 文件 | 作用 |
|---|---|
| `tests/river_hist_baseline_harness.F90` | 驱动**真实** `one` 写出路径（`MOD_Vector_ReadWrite`）的 MPI harness，产出参考文件并分阶段计时 |
| `tests/run_river_hist_baseline.sh` | 对着 `lib/libcolm.a` + `.bld` 编译并在两个 rank 规模 × 3 次运行取中位数 |
| `tests/river_hist_schema_lock.py` | 变量名/维序/属性/missing value/时间轴 + **按全局 ID 重算**的数值校验 |
| `tests/test_river_hist_schema_lock.py` | 用合成损坏样本证明能抓到变量缺失、维度转置、ID 错位、时间轴与属性差异 |
| `tests/river_hist_memory_baseline.py` | 按计划 1.5 公式估算 root 峰值内存 |
| `tests/test_river_hist_baseline_step1.py` | 固定公式、校验 harness/驱动接线 |
| `tests/artifacts/river_hist_ref_r{4,8}.nc` | 参考输出（**未提交**：可由驱动脚本重新生成，不把二进制放进仓库） |

## 1.1 / 1.2 参考输出与 schema 锁定 — 部分通过

**锁住的是写出机制，不是生产变量集合。** harness 写的是合成变量 `f_ucat_synthNN`，
不是 `hist_grid_riverlake_out` 实际产生的 river / Levee / reservoir / sediment /
tracer 全集——那需要 landdata 和 forcing，本机没有。因此下面验证的是"gather → 重映射
→ 串行写"这条路径本身与 ID 归位，生产 schema 的锁定要靠
`tests/run_river_hist_mode_parity.sh` 在目标机器上完成。

harness 走的是真实 writer（`vector_gather_to_master` /
`vector_gather_map2grid_and_write` / `vector_gather_matrix_to_master`）与真实
SPMD 布局（`spmd_init` + `divide_processes_into_groups`），并**始终保留一个零长度
worker**。数值编码为 `global_id + 1e-3*ivar + 1e-6*itime`，因此 ID 错位表现为数值
不符而不是"看着合理"的场。

结果：

- 4 ranks 与 8 ranks 的输出 **schema 完全一致**（维度、维序、属性、missing value、
  时间轴）；
- 两个规模下 16 个 unitcat 场 + 2 条 BIF 记录的**每个全局 ID 都落在正确格点**；
- 校验器自测 7 项全过，能识别计划点名的三类错误。

**发现的一个事实**（不是缺陷，但聚合器必须遵守）：`f_bifflw_lev` 在磁盘上是
`(time, bifurcation_pathway, bifurcation_level)`。CoLM 的 `ncio_write_serial_time`
按 Fortran 快维在前接收维名，落盘为 C 序，所以 Fortran 的 `(level, pathway)` 数组
变成 `(pathway, level)`。第 4 步重建该变量时必须按维名定位轴，不能假设顺序。

## 1.3 / 1.4 分阶段计时与 rank 标度

`large` profile：`totalnumucat=200000`、grid `1440x720`、16 个 unitcat 变量、
`totalnpthout=20000`、`npthlev=3`、2 个 history 时间点。3 次取中位数。

| ranks | workers | gather_only | full_write | netcdf(推算) | bif_matrix |
|---|---|---|---|---|---|
| 4  | 2 | 0.0213 s | 0.694 s | 0.674 s | 0.00040 s |
| 8  | 5 | 0.112–0.641 s | 0.725–0.988 s | 0.34–0.61 s | 0.0087–0.052 s |
| 12 | 8 | 1.746 s | 1.996 s | 0.250 s | 0.105 s |

标度倍数——**这些数字不可用作门槛或立项依据**：

ranks=8 有两次独立测量，中位数相差 **5.7×**（0.112 vs 0.641 s）。取哪一个，结论
完全不同：

| 取 8-rank 值 | 4→8 | 8→12 |
|---|---|---|
| 0.641 s | 30.1× | 2.7× |
| 0.112 s | 5.3× | 15.6× |

本报告早先版本写的"4→8 增 30×、8→12 增 16×"是**把两次不同的 8-rank 测量拼接**
得到的（30× 用 0.641，16× 用 0.112），两者不可能同时成立。该表述已从这里和计划
正文中撤除。

**能成立的只有定性结论**：master gather 随 rank 数增长且快于线性；4 ranks 时它只占
`hist_grid_riverlake_out` 的 3%，12 ranks 时已占多数。这与计划的机制诊断一致（每个
变量一次串行 `mpi_recv` 循环 + 一次全局 `mpi_barrier`，本例每次 history 写出 21 次）。

**到此为止。** 笔记本超订阅 + 5.7× 的运行间波动，使任何具体倍数都不可外推。真实标度
必须在目标机器的生产 rank 数上重测——这也是 go/no-go 仍未闭合的原因之一。

## 1.5 root 内存基线 — **未达到 go/no-go 门槛**

按计划公式，`workers=512`、每节点 256 GB：

| 分辨率 | totalnumucat | grid | peak root 下界 | 占单节点 |
|---|---|---|---|---|
| 15 arcmin | 86,400 | 1440×720 | 9 MB | 0.00% |
| 6 arcmin | 540,000 | 3600×1800 | 56 MB | 0.02% |
| 3 arcmin | 2,160,000 | 7200×3600 | **225 MB** | **0.08%** |

（估算值。实测 RSS/HWM 仍待在目标机器上采集；结论"内存不是瓶颈"在三个数量级的
裕度下不会被实测推翻，但严格说这一栏是下界而非观测。）

`vector_gather_map2grid_and_write` 在每个变量末尾 `deallocate` 掉 `wdata` 和
`wdata2d`，因此峰值是**单变量**的，不随变量数累积。即使全球 3 arcmin，root 侧也
只有约 225 MB——**远低于"节点可用内存 25%"的门槛，不构成 OOM 上限**。

> 这一条**推翻了我在计划评审里的建议**。我当时主张把内存能力提到首要收益，理由是
> "加机器也解决不了的硬天花板"。按公式实算，这个天花板在 3 arcmin 全球尺度上仍只有
> 225 MB，不成立。计划目标一节的收益优先级应据此调整：**主要收益是全局停顿，不是内存**。

## 1.6 go/no-go

三条判据：

| 判据 | 状态 |
|---|---|
| 河网 history 输出 ≥ 端到端墙钟 5% | **未测**——需要目标机器上的真实 case |
| 单次 history 全局停顿已影响作业稳定性/超时 | **未测**——需要你的作业记录 |
| root 内存 ≥ 节点可用 25%，或已 OOM | **不满足**（3 arcmin 全球仅 0.08%） |

**判定：暂不进入第 2 步。** 按计划"若三项都不满足，停止实施并保留测量报告"。

但 gather 的超线性标度是**支持继续**的实质证据，只是缺少能触发判据的那一个数字。
需要在目标机器 + 真实 case 上补测的三项：

1. `hist_grid_riverlake_out` 占端到端墙钟的比例（用现有 `MPI_Wtime` 计时点，
   注意 `-fdefault-real-8` 陷阱，见下）；
2. 同一 case 在两个生产 rank 规模下的 gather 耗时，确认笔记本上的超线性在真机复现；
3. routing compute（含 BIF global-dt lockstep）的耗时，用于区分"开 BIF 后慢"到底
   是算的慢还是写的慢——计划第 39-42 行已要求这一点。

若第 1 项达到 5%，或第 2 项显示 gather 在生产规模下已是主导项，则 go；否则本计划应
让位给 BIF 求解器优化。

## 附带发现：`MPI_WTIME` 在本仓库的编译选项下返回 0

`share/MOD_SPMD_Task.F90:34` 用 `include 'mpif.h'`，其中 `MPI_WTIME` 声明为
`DOUBLE PRECISION`。构建使用 `-fdefault-real-8` 而**没有** `-fdefault-double-8`
（`include/Makeoptions.Mac-arm`、`Makeoptions.github` 均如此），gfortran 会把该声明
提升为 `real(16)`，于是按 16 字节读取库返回的 8 字节值，结果恒为 0。

harness 已改用 `system_clock` 规避。**任何在模型内部加 `MPI_Wtime` 计时点的人都会
踩到这个坑**——第 1.3 步要在真机加计时时必须注意。当前模型代码没有使用 `MPI_WTIME`，
所以这不是现存缺陷，但会直接影响接下来的测量工作。
