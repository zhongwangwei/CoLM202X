# 组合验收缺陷的根因修复（最终方案）：无载体表层同位素质量随其 resolved 表层水输运

状态记录：2026-10-29。提交 `3b0af4f6`（分支 `feature/constituent-transport`，基于 `6d995438`/HEAD `10e9fbd2`
的生产代码）。本文记录**最终采用**的根因修复。此前尝试 `9f7c9a52`（把无载体同位素质量存入
`trc_surface_residue`、复水再溶解）已被**撤销**（`198da79`）：水同位素质量就是水本身，不存在“干残留相、
复水再溶解”的物理机制；该方案只对非挥发溶质成立。撤销后重新诊断，找到并实现了“示踪随水”的修复。

## 缺陷（来自组合验收 FAIL）

`DEF_GridRiverLake_FloodFeedback=.true.` + 运行期 IsoGSM 同位素强迫的组合在 LULC_IGBP/NEC 上
step 4 / 2005-12-25-05400 中止：

```
TRC_BAL residual hard failure: n_entries=1 worst_abs=9.0882E-05 tol=5.0E-05 @ipatch=213 itrc=1 ptype=0
```

## 诊断（逐点透明打印作业 `codex-01984a85605e40fb90f24a7bd1d18ffe`，只加打印、不改数值）

对 patch 213（H2_18O）逐 step 打印进入状态、宿主水通量、池组成与终态（原始行在
`residual-fix/v2/origin-dbg-lines.txt`）：

| 时刻 | 打印量 | 值 |
| --- | --- | --- |
| step 3 末 | `wdsrf`, `trc_wdsrf` | 0.457898 mm, 9.0765e-5 kg/m²（浓度 1.98e-4） |
| step 4 初 | `wdsrf_bef`, `trc_wdsrf` | 0.457898 mm, 9.0765e-5 kg/m² |
| step 4 通量 | `qinfl`, `rsur`, `flood_water`, `top_boundary_out_water` | −7.3146e-3 mm/s, 1.7476e-3 mm/s, 5.758e-7 mm, 13.1664 mm |
| step 4 终 | `wdsrf` | **10.0 mm（全部来自土壤 exfiltration，`top_exfil_water`=13.0999 mm）** |
| step 4 重建 | `surface_base_water = wdsrf+rsur·dt+qinfl·dt−flood` | −0.0206 → 0 |
| step 4 池 | `trc_pool_total`, `trc_soil_upflow` | **9.0765e-5**, 0.025973 kg/m² |

**根因**：普通混合池是按“洪水前的普通水”重建的分配工具，但在该步它重建为 0 水（旧水被蒸发/入渗消耗，
最终 10 mm 塘水是 exfiltration+flood 的 late 池）。同位素 tracer 本身仍对应**同一个物理表层水体**
（late 池的 10 mm 水），而旧代码在普通池塌陷时把 tracer 全部记成数值源汇并清零，
`late_tracer = trc_pool_total + trc_soil_upflow` 因此丢掉这 9.0765e-5，累计超过 5e-5 硬限而 fatal。
late 池混合后的浓度是合理的（(9.0765e-5+0.025973)/10 ≈ 2.6e-3），并不会产生非物理值。

## 修复（`3b0af4f6`）

在普通混合池“水为 0 但 tracer 非 0”的塌陷分支：

- **非挥发溶质**：保持原行为（存入 `trc_surface_residue`，复水再溶解）——这是溶质的真实不溶相。
- **水同位素**：若存在 **resolved late 载体水**（`flood_water > 0` 或 `late_surface_water > trc_tiny`），
  **不清零 tracer 池、也不记数值源汇**；只把普通池 ratio 置 0（普通 runoff/infiltration 不能在没有载体水时
  输出 tracer），随后既有 late-pool 混合把该 tracer 与该步真实的表层水（exfiltration/flood/dew）混合，
  `trc_wdsrf` 得到单塘混合浓度。示踪质量随它所属的水走，不新增任何“干残留”状态。
- **真正无载体**（既无 flood water 也无 late surface water）：保留原有显式数值 residual 记账，
  使真实不一致继续在 balance 报告里暴露（fail-closed 只留在这里）。

没有放松任何阈值、没有裁 δ、没有新增 restart/history 字段、没有改水/能量物理；`trc_surface_residue`
对同位素完全不再使用（验证作业里同位素 residue 全为 0）。

## 已验证（组件级，作业 `codex-fc0e28c1093d4a968a5f36df4cdc6849`）

- 构建 PASS（带 MT+GridRiverLakeFlow 的 LULC_IGBP）；二进制 SHA256：`colm.x=589b0fc3…`、
  `mkinidata.x=4a36bd83…`（`residual-fix/v2/binary.sha256`）。
- 触发案例（IsoGSM 2005-2006 同位素强迫 + `DEF_GridRiverLake_FloodFeedback=.true.`，NEC，
  2005-12-25 00:00–05:00，12 步）：**`CoLM Execution Completed.`，步数 12，零
  `TRC_SIG/TRC_BAL/residual hard/residual note/STOP/abort`**；restart 中同位素
  `mt_surface_residue`/`mt_subsurface_residue` **全为 0**（确认未用干残留）；仓库 9 项自检全部 PASS
  （`residual-fix/v2/validation.json`、`v2/log-markers.txt`）。
- 原组合验收的失败点在修复后不再出现（12 步含原 step 4）。

## 进行中的最终验收（终态前不写通过）

| 作业 | 目的 | 状态 |
| --- | --- | --- |
| `codex-64f6cca13a0d487ba438bab454d6d413` | combined10 全量验收重跑（两套 LULCC 构建 + 10 案例 + G1–G7 全门禁，源码 `3b0af4f6`，48 CPU/8h） | 运行中 |
| `codex-8ff87f9ecf494c0a9754cc9333553f11` | 修复 vs `6d995438` 回归对照（相同配置构建/运行；A: JRA+同位素+CL+泥沙+河网+flood feedback 3 天；B: JRA+同位素 3 天） | 运行中 |

回归对照的预期：宿主水/能量/网格状态逐位相同；差异只可能出现在塌陷路径触发时的 tracer 池
（新方案把质量留在池中而非记 residual），需按实测值逐项核对后写结论。

## 边界与未覆盖

- 本修复只处理“普通表层池水重建为 0、但存在 resolved late 载体水”的同位素情形；真正无载体的塌陷仍记数值
  residual（可能仍 fatal，这是有意的可见性保留）。
- 其它数值残余路径（土壤液相/冰相边界、雪层、蒸发边界、含水层 orphan）未修改。
- 组件级验证不等于组合验收：以 `codex-64f6cca1…`（G1–G7）与回归对照终态为准。
- 未覆盖 LULC_IGBP_PC、其它区域/强迫产品；观测准确度与全球代表性未验证。
