# （已撤销的第一次尝试）无载体干表层同位素质量的“隔离库+复水再溶解”方案

> **状态：本方案已撤销（提交 `198da79` 撤销 `9f7c9a52`）。** 水同位素质量就是水本身，不存在“干残留相、
> 复水再溶解”的物理机制；该机制只对非挥发溶质成立。撤销后重新诊断并实现“示踪随水”的最终修复，
> 见 [late-carrier-fix.md](late-carrier-fix.md)（提交 `3b0af4f6`）。本文保留为失败/被撤销方案的完整记录，
> 其组件级验证与回归差异数字仍可作为“该方案在数值上可行但物理上不可接受”的反证。

状态记录：2026-10-29 DeepSeek 接续。第一次尝试提交 `9f7c9a52`（分支 `feature/constituent-transport`，
基于 `6d995438`/HEAD `10e9fbd2` 的生产代码）。本文给出根因诊断、修复内容、组件级验证证据、
边界与未完成项。**该方案已撤销；最终修复与验收见 late-carrier-fix.md。**

## 缺陷（来自组合验收 FAIL）

`DEF_GridRiverLake_FloodFeedback=.true.` + 运行期 IsoGSM 同位素强迫的组合在 LULC_IGBP/NEC 上
step 4 / 2005-12-25-05400 中止：

```
TRC_BAL residual hard failure: n_entries=1 worst_abs=9.0882E-05 tol=5.0000E-05 @ipatch=213 itrc=1 ptype=0
 TRC_BAL hard failure: numerical residual exceeds independent threshold
```

## 根因定位（两个 7920 私有诊断作业，只加打印、不改数值）

作业 `codex-3f19455717724b49af2ec0a85bad17ac`（逐点 `trc_numerical_residual_step` 写入打印）与
`codex-5119aba74f9a4ee2a3bd926451d695e2` / `codex-7d77108954fc4101b70480bf7979ae5c`（物理量打印）：

- 记账位置唯一：`main/MATERIAL/GENERIC_WATER/MOD_MT_SoilWater.F90` 的**表层混合池“无载体塌陷”分支**
  （原 1322 行，`trc_surface_collapse_residual`）。
- 触发时实际量（patch 213, itrc=1/H2_18O）：`wdsrf=10.0 mm`、`rsur=1.7476e-3 mm/s`、
  `qinfl=-7.3146e-3 mm/s`、`flood_water=5.7582e-07 mm`、`late_surface_water=0`、
  `surface_base_balance=-2.0647e-02` → `surface_base_water=0`（可分辨载体水=0）；
  同时该 patch 表层同位素库存 `trc_wdsrf=9.0765e-05 kg/m²`、`trc_soil_upflow=2.5973e-02 kg/m²`。
- 即：表层**可分辨载体水被洪水/入渗记账清零**，但对应的 **9.0765e-05 kg/m² 同位素质量仍在表层池**；
  旧代码把全部剩余 tracer 池记成 `trc_numerical_residual_step`（数值源汇），该值超过代码自身
  独立硬限 `trc_resid_abs_tol=5.0e-05`（`MOD_MT_Conservation.F90`），于是 fatal。
- 该路径对**非挥发溶质**早已有正确处置：存入既有 `trc_surface_residue`（不可移动隔离库），
  复水时溶解回混合池。**同位素此前被排除在外**，是本次缺陷的根因。

## 修复（`9f7c9a52`）

把“无载体质量隔离”规则从非挥发溶质**推广到水同位素**，复用既有状态量，不放松任何阈值、不裁 δ：

| 位置 | 修改 |
| --- | --- |
| 表层塌陷（soil surf） | 同位素剩余池 → `trc_surface_residue`，不再记数值源汇；复水时溶解 |
| 洪水 late-pool 残余 | 同位素小残余 → `trc_surface_residue`（超限 fatal 分支保留） |
| dew/frost late surface 无载体 | 由 fatal 改为隔离入库（该场景质量合法、无载体） |
| 湿地池塌陷 | 同位素 → `trc_surface_residue` |
| 4 处复水/混合池入口 | 溶解条件扩展到同位素（表层混合池、洪水 late pool、late surface pool、湿地混合池） |
| `MOD_MT_Hist.F90` | 同位素也输出 `f_mt_surface_residue_<species>` / `f_mt_subsurface_residue_<species>`（可审计；dry-layer / solid inventory 仍仅溶质） |
| `MOD_MT_Vars.F90` | 更新 `trc_surface_residue` 语义注释 |

`trc_surface_residue` 本来就是守恒/历史/restart/LULCC remap/LULCC 审计都覆盖的 storage 组件，
因此该修复**不引入新状态字段、不改变 restart schema**；把质量从“数值源汇”搬到“物理隔离库存”，
使 `TRC_BAL` 严格闭合。

## 组件级验证（作业 `codex-6696730ffec34c518f011b395657e983`）

- 构建：`LULC_IGBP`（无 LULCC）与 `LULC_IGBP+LULCCON` 两套 MT+GridRiverLakeFlow 配置 `make all` 均 PASS；
  define.h 实测 `MaterialTransport`、`GridRiverLakeFlow`，LULCC 仅第二套定义；二进制 SHA256 见
  `residual-fix/binary-*.sha256`。
- 触发案例（IsoGSM 2005-2006 同位素强迫 + `DEF_GridRiverLake_FloodFeedback=.true.`，NEC，
  2005-12-25 00:00–23:00，47 步）：两套构建均 **`CoLM Execution Completed.`**，
  **零** `TRC_SIG step report` / `TRC_BAL step report` / `residual hard failure` / `residual note` /
  `ERROR STOP` / `MPI_Abort`（修复前同一案例在 step 4 abort）。
- 同位素隔离库真实出现并被审计：restart 中 `mt_surface_residue` 形状 `(587 patches, 2 tracers)` 有限非负；
  history 出现 `f_mt_surface_residue_H2_18O/HDO` 与 subsurface 变量、无缺测；日末数值很小，
  说明 step-4 隔离的质量在后续复水/洪水 late-pool 中已重新溶解并继续输运（守恒闭合）。
- 仓库自检 9 项全部 rc=0：`compare --self-test`、`iso_check --self-test`、`check_forcing_cache --self-test`、
  `check_forcing_edges`、`check_singlebucket_precip`、`check_orphan_leaf_inventory`、
  `check_frcsat_history`、`check_regionclip_paths`、`check_defineh`（93 唯一用例）。
- 证据：`residual-fix/validation.json`、`residual-fix/*-log-markers.txt`、`residual-fix/define-*.h`、
  `residual-fix/environment.txt`；诊断作业日志在 7920 对应 job 的 `results/logs/`。

## 修复回归对照（已完成：舍入级惰性，差异已逐项定性）

作业 `codex-d64e0d0b17264f8db7d81038056a6030`：用同一套 `MT+GridRiverLakeFlow`（无 LULCC）配置分别构建
fixed（`9f7c9a52`）与 base（`6d995438`），各跑两个 3 天案例：

- A `all`：JRA 动力强迫 + 静态同位素 + CL + 泥沙 + levee/bifurcation/reservoir + flood feedback；
- B `iso`：JRA 动力强迫 + 静态同位素（无 flood feedback、无 CL/泥沙/河网选项）。

结果（作业终态 failed/1 只是严格比较器/复核脚本对“预期差异”返回非零；4 个模型 run 全部 143 步 END、无诊断）：

| 比较 | 结果 |
| --- | --- |
| restart 宿主/网格/水量状态 | **逐位相同**（`case_restart_gridriver_*`、`case_restart_hist_*` 0 差异；`mt_wliq/wice/wdsrf` 等仅 4 个 tracer 状态差） |
| restart tracer 状态差异 | `all`：`mt_wliq_soisno` ≤8.7e-19、`mt_wice_soisno` ≤3.5e-18、`mt_wdsrf` ≤2.1e-20、`mt_surface_residue` ≤1.8e-18（kg/m²）；`iso`：相应 ≤5.0e-17 / 3.5e-18 / 1.1e-16 / 1.3e-16 |
| history 变量集 | fixed 仅**新增** 4 个诊断变量：`f_mt_surface_residue_{H2_18O,HDO}`、`f_mt_subsurface_residue_{H2_18O,HDO}`；base 无多余变量 |
| history 公共变量差异 | `all` 8 个、`iso` 12 个，全部是 tracer 浓度/δ/evap 质量诊断；max abs 4.4e-13 / 1.1e-11（源于极小的 tracer 差与除法放大），**宿主水/能量/网格状态无差异** |

**结论**：该修复对不触发“无载体干池”路径的配置是**舍入级惰性**（差异 ≤1.3e-16 kg/m² tracer 状态、
≤1.1e-11 诊断），只新增 4 个可审计的同位素隔离库 history 变量；它只改变“本来会被记成数值源汇并可能 fatal”
的隔离质量去向（这正是修复目的）。严格原始位比较器的非零返回被如实保留，不能写成“逐位一致通过”。

## 当时计划的完整验收（**本方案已撤销，相关作业均已取消/过时，记录仅作历史**）

| 作业 | 目的 | 状态 |
| --- | --- | --- |
| `codex-d48a500ecf994f2fb89a6e2856a1cc6a` | combined10 全量验收重跑（两套 LULCC 构建 + 10 案例 + 全部 G1–G7 门禁，源码 = `9f7c9a52`） | 运行中（IGBP/NEC 原 abort 案例已越过 step 4，接近 12 天窗口尾部） |
| `codex-d64e0d0b17264f8db7d81038056a6030` | 修复 vs `6d995438` 回归对照 | **已完成**（上表；终态 failed/1 是严格比较器对预期差异返回） |
| 待办 | 修复提交后的完整 hosted CI 矩阵；p9r 类验收在修复源码上的进一步结论 | 待办 |

## 边界与未覆盖

- 组件级验证只证明：原 abort 消失、模型 END、无 TRC 诊断、隔离库可运行、仓库自检不破。
  完整组合验收（G1–G7）与回归对照结论以两个进行中作业的终态为准。
- 修复只影响“无载体干池”这一条路径。**在会触发该路径的配置里，同位素状态/历史会与修复前不同
  （这是修复目的）**；不触发该路径的配置应逐位不变，由回归对照作业验证。
- 其它数值残余路径（土壤液相/冰相边界 1527/1747/2007/2062/2090 行、雪层 180/197 行、蒸发 375/397/426/448 行）
  **未修改**，仍按原独立硬限记账；若进一步组合覆盖触发，会再次显式暴露而非隐藏。
- 未修改阈值、未裁 δ、未改水/能量物理、未改 restart schema；isolated residue 不进入水收支。
