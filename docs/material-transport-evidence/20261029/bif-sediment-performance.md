# N01 BIF 泥沙子步化性能记录（2026-10-29）

- 分支：`feature/constituent-transport`
- 基线：`87c0435b`（优化前）；优化提交见文末 HEAD
- 环境：本机 macOS（gfortran 16 + OpenMPI + NetCDF-Fortran）用于编译/回归；性能实跑在 7920。

## 1. 已获得的真实运行测量（优化前）

7920 short smoke（job `codex-014c8aeb…`）：
- 配置：GRID + USEMPI(5 ranks) + MaterialTransport，JRA 强迫，2005-01-01…01-03（3 天，dt=1800 s），
  开启 isotope + CL + sediment + GridRiverLakeFlow + BIF + levee + flood feedback；
- 结果：**job 在 3600 s timeout 时只输出 2 条 TIMESTEP 记录**（`state=timed_out`），
  日志尾部仍停留在首个/第二个时间步的材料诊断。构建与 mkinidata 均 PASS。
- 结论：BIF 路径的子步化在该配置下不可接受地慢（或停滞），必须优化后才能用于验收。

## 2. 已实施的低风险优化

1. **超小载水单元不再驱动全局子步数**（本次关键修复）
   - 位置：`grid_sediment_calc` 与 `calc_sediment_advection` 的 CFL 循环；
   - 规则：`rivsto_donor > SED_MIN_SUBSTEP_CARRIER (=0.1 m³)` 才参与全局 `nsub` 取最大；
   - 依据：近干涸/微小载水单元的周转比可趋于无穷，会把所有单元的 BIF 子步数拉到上限；
     这类单元的泥沙库存近零，其出流已被 advection core 的 donor limiter 限制，不构成
     “已解析载水单元”的恒浓度约束；
   - 不变量：对 `V > 0.1 m³` 的单元仍严格执行 gross-outflow ≤ donor-water 的子步；
     仍无静默截断（N09 的 >240 显式报错保留）。
2. **工作区复用**：核心子步的 `w_don/w_end` 改为持久分配（仅维度变化时重分配）；
   provider 的 `riv_don_step/riv_step` 仍是每次 `grid_sediment_calc` 分配一次、子步内复用。
3. **低开销探针**（无扫描、无通信）：模块级整数计数器
   `sed_adv_substep_calls/sum/max`、`sed_bif_substep_calls/sum/max`，通过
   `mt_sediment_substep_stats()` 读取，供后续性能对照。

## 3. 性能对照表（按任务要求，未实测项标 NOT_MEASURED）

| 指标 | 优化前 | 优化后 | 改善 |
|---|---:|---:|---:|
| BIF 总耗时 | NOT_MEASURED（仅有 timeout 证据） | NOT_MEASURED | — |
| prepare_bif_sediment 耗时 | NOT_MEASURED | NOT_MEASURED | — |
| 平流耗时 | NOT_MEASURED | NOT_MEASURED | — |
| 平均子步数 | NOT_MEASURED | NOT_MEASURED | — |
| MPI 等待 | NOT_MEASURED | NOT_MEASURED | — |
| 总运行时间 | 3 天 case >3600 s 超时（2 步） | NOT_MEASURED | — |
| 最大质量预算误差 | NOT_MEASURED | 本地回归 0 failure | — |
| 最大浓度相对误差 | NOT_MEASURED | 本地回归 0 failure | — |

说明：本轮没有可用的优化后完整模型运行数据，绝不把理论估算当性能结果。

## 4. 正确性回归（本轮后）

- `python3 tests/material_transport/check_audit_p1_fixes.py` → 0 failure；
- `check_sediment_joint_substep.py`、`check_terminal_sed_history.py`、`check_push_mpi.py`（真实 2-rank MPI）继续通过；
- GRID 与 BGC+CROP 两配置 `make all` → MAKE_RC=0；
- N10 水量随子步推进已验证编译与结构/不变量检查（完整模型恒浓度验证仍待 7920）。

## 5. 未完成/下一步

1. 用新 HEAD 重跑同一 7920 smoke，记录：实际 wall time、TIMESTEP 数、`mt_sediment_substep_stats` 的
   平均/最大子步数、`t_adv` 计时与 MPI 等待；这是验证“性能确有下降”的必需证据。
2. 若仍慢：对 `prepare_bif_sediment` 做静态拓扑/几何缓存与“每子步合并通信”，并评估按单元自适应
   （非全局最大）的子步方案；快速模式需与完整算法并行对照并报告误差。
3. 完整 BIF/levee 多分支数值验收（本任务第六节矩阵）。

## 6. 决策更新：BIF/levee 子步已按项目决定回退为单步

- 背景：provider 全局子步（commit `bed6193b`）被单个近干涸单元拖垮；7920 实测 3 天 case
  在 3600 s 内只推进 2 个 TIMESTEP。
- 决定（2026-10-29）：**BIF/levee 路径退回单步执行**，只保留普通（非 BIF/levee）advection 的
  材料子步。回退后：
  - `grid_sediment_calc` 对 BIF/levee 算子每个 morphology interval 只执行一次（原 dt_morph、
    原 reservation/credits/levee 顺序）；
  - `calc_sediment_advection` 普通路径子步、N09 普通路径 >240 显式失败、N10 普通路径水量插值保留；
  - `sed_bif_substep_*` 计数器保留（恒为 1 次/区间），用于确认没有意外的重复调用。
- 正确性影响：BIF/levee 高周转下的恒浓度不变量**不再强制执行**（已知限制，回归测试
  `check_sediment_joint_substep.py` case4 记录该单步结果，不再声称 BIF 路径满足 N01 恒浓度）。
- 性能影响：预期 BIF 路径回到优化前的一次算子开销（未实测完整模型，标记 NOT_MEASURED）。
- 回退后验证：GRID 与 BGC+CROP 两配置 `make all` 均 MAKE_RC=0；本地回归 0 failure。

## 7. 退回单步后的真实速度测试（2026-10-29）

同源同配置（GRID+USEMPI 5 ranks、JRA 强迫、MT+BIF+levee+sediment+isotope+CL）在 7920 上：

| 项目 | 优化前（provider 全局子步） | 退回单步后 |
| --- | --- | --- |
| job | codex-014c8aeb…（3 天 case，timeout 3600 s） | codex-fef2905b…（2 天 case，timeout 1800 s） |
| 结果 | timed_out：**2 个 TIMESTEP / 3600 s** | **succeeded，exit 0** |
| 步数 | 2 | **95**（2005-01-01 00:00 → 01-02 23:30，dt=1800 s） |
| wall time | ≥3600 s（未完成） | **976.5 s**（started_at 1793252759.05 → finished_at 1793253735.55，含构建+mkinidata） |
| 每步成本 | ≥1800 s/step（未完成，下限） | ≤10.3 s/step（含构建/初始化摊薄，实际模型步更快） |
| 错误 | — | 无 ERROR STOP / MPI_Abort / residual hard failure；`CoLM Execution Completed`=1 |
| 加速比 | 1× | **约 175×**（按保守的每步 wall 时间；旧数据未跑完，只能取下限） |

- 结论：退回单步后 BIF/levee 配置回到可用速度，2 天 case 约 **16.3 分钟**跑完（含构建与初始化）；
  按线性外推 3 天 case 约 **25 分钟**（外推值，非实测）。
- 注意：两次 case 长度不同（3 天 vs 2 天），但都使用相同域、进程数、物理开关与 dt，每步成本可比。
- 仍未实测：`prepare_bif_sediment`/advection 分项耗时、MPI wait、内存峰值、不同河网分区与粒级矩阵；
  这些在单步回退后属于“原单步算法成本”，不再有全局子步放大。
- 数值验收边界：BIF/levee 恒浓度不再强制（已知限制，见 `check_sediment_joint_substep.py` case4）；
  本次成功运行证明“无崩溃/无硬失败”，不是完整 BIF 物理验收。
