# p9r: correct-baseline unforced isotope/CL regression (accepted)

接续观察时点 2026-10-29 00:26:59 CST。本文只记录该作业的实际终态与结果，不推广到未测组合。

## 作业终态

| 项 | 值 |
| --- | --- |
| Job | `codex-4ba194bf50074f68aaa3d76178c49830` |
| 终态 | `succeeded`, exit 0（2026-10-28 21:39:34 → 2026-10-29 00:26:59 CST，2h47m） |
| 配额 | 24 CPU / 48 GB / timeout 14400 s |
| 源码 | old = `28478b18e6bb87db7bff1af7eb539812563110a3`；new = `c315c84973d0672d120b2b6926343a1b42b6ce21`（git archive，无未提交修改） |
| 配置 | GRID / LULC_IGBP / MT / GridRiverLakeFlow / BGCOFF / CROPOFF；无 LULCC，无同位素强迫；生产默认 PHS；9 MPI ranks；最多两例并发 |
| 输入 | RAW=`/media/zhwei/data02/zhwei/CoLMrawdata`, RT=`.../CoLMruntime`, JRA=`.../JRA3Q`, UCAT=`.../grid_routing_data_15min.nc`, GRIDJOB=`claude-2026f9ac...`（landdata 只读声明） |

## 结果（`report.json`, `validation.log`）

- **10/10 模型例 PASS**：旧/新源码各重跑 5 种无强迫配置——isotope NEC、isotope AMZ、CL NEC（2005-12 整月，各 1,488 步/2005-12-31-84600），isotope+CL+sediment+routing NEC/AMZ（2005-01-01..03，各 144 步/2005-01-03-84600）。
- 每例均验证：描述符身份（iso→{H2_18O,HDO}；CL→{CL}；all→三者）、诊断失败 0、无 SIG/BAL/residual hard/ERROR STOP/MPI_Abort，forcing cache 按无强迫配置预期为 inactive。
- **10/10 严格比较 PASS**（完整文件集、变量集、形状、dtype、原始变量字节）：
  - restart：各 7 文件 / 289（iso）、283（CL）、438（all）变量；**每例恰好 1 处标量 `mt_forcing_cache_schema` 1→2 迁移**，且两侧 `mt_forcing_cache_count=0`、无 `mt_forcing_precip_reliable`；其余全部原始位相同。
  - history：各 3 文件 / 225（iso）、191（CL）、276（all）变量；0 处迁移，全部原始位相同。
- `failures: []`。验证器自测先 PASS：接受"仅 schema 标量迁移"，拒绝 signed-zero、active-cache、额外变量。
- 两次独立构建 PASS（old/new，相同 Conda MPI/GNU/NetCDF 与 Fortran 选项）；构建日志与源清单随任务 `work/` 保存。

## 证据

- `p9r-report.json`（完整 report.json，含比较计数）
- `p9r-job-status.json`（runner 终态与输入声明，已去掉 cpu_affinity 长列表）
- `p9r-task-log.txt`（runner task.log：构建/逐例/比较记录）
- `p9r-source-versions.txt`（不可变源码提交）
- 远程完整结果：`/media/zhwei/data01/compute_runs/codex/codex-4ba194bf50074f68aaa3d76178c49830/results/`；构建/驱动冻结于同目录 `work/`。

## 边界（不得外推）

- 只证明 `c315c849` 相对 `28478b18` 在五种**无强迫**配置下回归一致；**不包含**随后 `6d995438` 的 frcsat 诊断字段增补，也不等于最新组合源码（frcsat+强迫+LULCC）验收——后者由独立作业覆盖。
- CL 使用标准零浓度参数，只验证 CL 描述符真实注册与输运数值一致；非零 CL 有单独历史证据，不由本作业替代。
- 不覆盖强迫同位素、LULCC、BGC、甲烷、非 9-rank 全部矩阵或观测精度；不比较 NetCDF 属性/编码。
