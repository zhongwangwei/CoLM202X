# 最新组合源码整体验收（frcsat + 强迫同位素 + LULCC + CL + 泥沙/河网）：FAIL，触发项已定位

状态记录：2026-10-29 DeepSeek 接续。生产源码 `6d99543861c38c708f6c833e6349bcd3d7d3b879`（对交付 HEAD `10e9fbd2`
的生产目录逐一相同，之后只有 docs/tests 提交）。本文记录严格验收的**实际结果**，不把子集通过写成整体通过。

## 验收结论

**要求的完整组合不通过。** 主验收作业 `codex-2434412927a94e38a44ad21db5643a38` 终态 **failed/1**
（2026-10-29 02:39 CST；runner 失败来自验收门禁 FAIL，不是基础设施故障）。LULC_IGBP/NEC 的六个组合案例
都在同一位置中止：

```
TIMESTEP = 4 | DATE = 2005-12-25-05400
TRC_BAL residual note (excluded from tol check): n_entries=       2 worst_abs= 0.90882E-04
TRC_BAL residual hard failure: n_entries=       1 worst_abs= 0.90882E-04 tol= 0.50000E-04 @ipatch=     213 itrc=   1 ptype=  0 owner=   0
 TRC_BAL hard failure: numerical residual exceeds independent threshold
Abort(1) on node 1 ... application called MPI_Abort
```

这是 `MOD_MT_Conservation.F90` 中被显式记账的 `numerical_source_sink`（来自 `MOD_MT_SoilWater.F90` 近干/洪水
耦合路径）超过代码自身独立硬限 `trc_resid_abs_tol = 5.0e-5`，**不是**普通聚合 `TRC_BAL`/`TRC_SIG` 检查。
没有放宽阈值、没有裁剪 δ。

| 门禁 | 结果 |
| --- | --- |
| G1 构建（LULC_IGBP / LULC_IGBP_PFT，MT+LULCC+GridRiverLakeFlow） | **PASS**（两套 `make all` 成功，define.h 宏核对一致） |
| G2 配置与完整位比较 | FAIL（IGBP/NEC 中止无法比较）；**子集事实**：PFT/NEC 的 4 组拆分对连续比较（restart+history 全文件全变量原始位）通过 |
| G3 模型完成 | FAIL（IGBP/NEC 中止）；**子集事实**：PFT 4 例 + AMZ 非零 CL 例均 576 步、`2006-01-05-84600` END、零诊断、无未掩码非有限值 |
| G4 forcing cache | FAIL（IGBP/NEC 无后续 checkpoint；交付工具对 isotope+CL+SEDIMENT 组合 cache 存在工具缺口，见下）；**子集事实**：PFT 每例 12 个初始化 checkpoint 字段审计合格（schema=2、forcing_variables=6、identity(7,1544)、reliability∈{0,1}、同位素列>0、非同位素列=0），AMZ 非零 CL 列=3e-6 符合配置 |
| G5 frcsat 协议 | FAIL（IGBP/NEC 无 history）；**子集事实**：PFT 每例 1,200 个有效要素-日 `f_frcsat`、0 掩码、12 个窗口边界 sidecar（`hacc_nac=0` 已重置） |
| G6 LULCC 注册/守恒 | FAIL（IGBP/NEC 未到过渡）；**子集事实**：PFT 4 例为真实 MEC/SAT transition，注册 `H2_18O,HDO,CL`，守恒诊断 0 |
| G7 同位素诊断 | FAIL（报告值，主作业有路径 bug）；**只读补核更正**：f-NEC-MEC、f-NEC-SAT、i-AMZ 的 iso_check rc=0（各 50 个同位素变量、约 1,822 万有效值、0 未掩码非有限、19–21 组 d-excess）；i-NEC-MEC 因中止无 history |
| 总体 | **FAIL**（`validation` report：status=FAIL，failures=89） |

### 校验器报告瑕疵（只读补核 `codex-81a0f423c967451abfcb5adb5ef28881`，不改变 FAIL）

1. G7 的 `iso_check.py` 被主作业在 `<run>/` 而非 `<run>/case/` 上调用，导致 f-NEC-MEC 与
   i-AMZ 出现假 rc=1；只读补核按正确路径重跑后上述子集 rc=0（数字见上表）。i-NEC-MEC 仍因中止无 history。
2. 交付工具 `tests/material_transport/check_forcing_cache.py` 假定纯同位素运行，要求
   `mt_forcing_species` 每一列的 `mt_forcing_precip_last/vapor_last > 0`；组合运行 `ntracers=4`
   （isotope+CL+SEDIMENT），CL/SEDIMENT 列合法为 0 或 3e-6，工具因此拒绝结构有效的 schema-2 cache。
   只读补核直接审计字段确认子集 cache 结构一致；**G4 在组合配置下无法仅凭该交付工具认证**，
   需要扩展工具（本轮未改工具）。IGBP/NEC 中止案例没有任何 initialized checkpoint。

## 触发项定位（单因素消融，全部复用冻结源码与既有构建，只改 namelist）

| 轮次 | 变体 | 结果 |
| --- | --- | --- |
| 1 (`codex-17e71324…`) | V0 控制 / V1 CL / V2 sediment / V3 reservoir+levee+bif / V8 sediment+CL | PASS（12 步） |
| 1 | V4 仅 `DEF_GridRiverLake_FloodFeedback=.true.`，及 V5/V6/V7 含 flood feedback | **DIAG_RESIDUAL**（step 4，9.0765–9.0882e-5 > 5e-5） |
| 2 (`codex-0c007acf…`) | W2 分馏 off、W3 SAT scheme | 仍 DIAG_RESIDUAL → 分馏/LULCC scheme 无关 |
| 2 | W4 完整配置去掉 flood feedback，5 天 240 步 | **PASS**，零诊断 |
| 3 (`codex-3c6c916b…`) | X1 JRA 动力强迫 + 静态同位素 + 全部 river 选项（含 flood feedback）+ CL + sediment | **PASS**（144 步） |
| 3 | X2 flood feedback + IsoGSM 运行期同位素强迫，**LULCC-off 构建** | 仍 DIAG_RESIDUAL → LULCC 编译开关无关 |

**最小触发组合 = `DEF_GridRiverLake_FloodFeedback=.true.` + IsoGSM 运行期同位素强迫**；
与 CL、泥沙、levee/bifurcation/reservoir、分馏、LULCC scheme 及 LULCC 编译开关无关。当前只在
LULC_IGBP/NEC 复现；PFT/NEC 与 IGBP/AMZ 子集未出现。

## 既有缺陷分类（不是最新提交的回归）

`codex-69298c1322fd4c749b69e9999c7b74d8` 用 p9r 已构建的两套无 LULCC 二进制跑同一 12 步触发案例：

| 二进制 | 源码 | 结果 |
| --- | --- | --- |
| old | `28478b18`（分支基线） | `worst_abs=0.90765E-04 > 5.0E-04` @ipatch=213 itrc=1 ptype=0，step 4（同一位置） |
| new | `c315c849` | 同一行、同一数值、同一 patch/tracer |

两套都失败 → **该交互缺陷在分支基线 `28478b18` 已存在，不是最新提交引入的回归**；旧验收从未把
flood feedback 与运行期同位素强迫同时打开（p9q 无 flood feedback；p9r 的 flood-feedback 案例无运行期同位素强迫），
因此这是本轮新测试**暴露的既有缺陷**。根因修复未做；需要根因修复或明确的残余策略决策（放宽
`DEF_MT_RESID_ABORT_NBAD` 只掩盖问题，不算通过）。

## frcsat 修复隔离（正面的数值不变证据）

`codex-05c172686e434105884e3a0c7297f6fb`（succeeded/0）：同一 NEC / IGDB / MT / 无 LULCC / 无同位素强迫 /
JRA / CL / 河网案例分别用 `c315c849` 与 `6d995438` 独立构建跑 120 步（结束于 2005-12-27 12:00，制造窗口中途 restart）：

- history **所有文件/变量（含 `f_frcsat`）原始位相同**；
- restart 13 文件、**767 个公共变量全部原始位相同**；唯一差异是修复新增的两个 mid-window 侧车字段
  （`hacc_nac=24`、`nac_frcsat=24` 全部 587 patches、`hacc_frcsat_count_version=1`）。
- 结论：frcsat 修复在真实 NetCDF/MPI 案例上**不改数值**；这是本轮唯一可正面引用的 frcsat 数值不变证据。

## 子集结果（不得当作整体通过）

四个 LULC_IGBP_PFT/NEC 案例在完整 12 天窗口内部干净：各 576 步（拆分 336+240 / 432+144）、终点
`2006-01-05-84600`、模型 END、零 SIG/BAL/LULCC_BAL/residual hard/STOP/abort、零未掩码 NaN/Inf；
2005→2006 边界为真实 MEC/SAT transition；注册物种 `H2_18O,HDO,CL`；每例 1,200 个有效要素-日
`f_frcsat` 且 0 掩码，12 个窗口边界 sidecar 且 `hacc_nac=0`（窗口在 restart 前已重置）；
**连续运行与两种拆分运行 restart/history 全文件全变量原始位一致（报告内 4 组严格比较）**。
IGBP/AMZ 非零 CL 案例 `i-AMZ-MEC-cl2e5` 同样 576 步 END、诊断干净（其 G4 问题仅为上述工具缺口）。
这些子集只说明触发不是全局的，不证明完整组合通过。

## 失败/尝试保留记录

主验收作业与所有探针作业编号、失败原因见 `combined-EVIDENCE.md`（attempts 表）与 `combined-attempts.json`。
早期失败均为脚本/检查器错误或配置错误（无模型物理结果），原始记录保留，不冒充成功。

## 证据

- `combined-report.json`：主作业完整验收报告（status=FAIL、逐门禁、逐案例、89 条 failures、limits）。
- `combined-results.json`：执行子任务的机器可读结论（verdict、门禁、触发定位、回归分类、未覆盖项）。
- `combined-job-status.json`：runner 终态（failed/1）与输入/资源声明。
- `combined-EVIDENCE.md`（由执行子任务维护）：源码冻结、作业/尝试表、消融与分类细节。
- `combined-attempts.json`：作业尝试记录。
- 远程：`codex-2434412927a94e38a44ad21db5643a38/results`（主作业）、`codex-17e7132400c44adeb691e66e3dc2c586`、
  `codex-0c007acf2ceb4863aebe3d3dbf812cd2`、`codex-3c6c916b96a142f6a17f88cbfe25c374`（消融）、
  `codex-05c172686e434105884e3a0c7297f6fb`（frcsat 隔离）、`codex-69298c1322fd4c749b69e9999c7b74d8`（分类）。

## 未覆盖 / 未做

- 根因未修复，残余策略未决策。
- 触发复现只覆盖 LULC_IGBP/NEC 与 28478b18/c315c849 两套二进制；未测 LULC_IGBP_PC 构建、其它区域/scheme/强迫产品。
- 干湖无效样本分支未在**完成的**组合案例中出现（标准案例 `DEF_USE_Dynamic_Lake=.false.`，valid count 等于步数）；
  协议侧车保存/恢复语义由 mid-window sidecar 与新旧隔离证明，mask 语义依赖冻结源码的 `check_frcsat_history.py` 自检，
  不声称由失败的主验收证明。
- 不比较 NetCDF 属性/文件编码；δ 范围只作诊断，不验证观测准确度。
