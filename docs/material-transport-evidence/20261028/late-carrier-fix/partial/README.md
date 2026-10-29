# 取消前 partial 证据（**不是验收通过**）

2026-10-29 10:03 CST 按用户要求“先停掉手头全部测试（已发现明显 bugs）”，以下作业被取消；
本目录只归档**取消前已完成的部分**的可复核结果，**不作为组合验收通过、也不作为回归对照通过**。

## 回归对照 partial（作业 `codex-8ff87f9ecf494c0a9754cc9333553f11`，cancelled/-15）

配置：fixed = `3b0af4f6` vs base = `6d995438`，同一 `MT+GridRiverLakeFlow`（无 LULCC）构建；
A = JRA 动力强迫 + 静态同位素 + CL + 泥沙 + levee/bifurcation/reservoir + flood feedback；B = JRA + 静态同位素。

| 比较 | 结果（取消前） |
| --- | --- |
| B(iso) history（3 天完整，143 步） | **逐位相同**：3 文件 / 225 变量 |
| B(iso) 最终 restart（2005-003-84600） | **逐位相同**：3 文件 / 311 变量 |
| A(all) restart（day-1，2005-002-00000，48 步后） | 3 个 tracer 变量有差异：`mt_wliq_soisno` ≤8.7e-19（1/26415 点）、`mt_wice_soisno` ≤3.5e-18（3/26415）、`mt_wdsrf` ≤7.6e-21（8/1761）；`case_restart_gridriver_*` 与 `case_restart_hist_*` **0 差异**（宿主/网格状态相同） |

解释：B 配置不触发“无载体塌陷+late 载体”路径 → 最终修复数值上完全惰性；A 配置触发该路径时只有 tracer
库存在舍入级差异（≤1e-18 kg/m²），宿主水/能量/网格无差异。这支持“修复不破坏已验收数值”，但 A 案例只跑完
1 天且对照作业被取消，**不构成完整回归通过**。

## 组合验收 partial（作业 `codex-64f6cca13a0d487ba438bab454d6d413`，cancelled/-15）

取消时 10 个案例跑到 **step 64–70**（2005-12-26），均**无** `TRC_BAL` / `TRC_SIG` / `residual hard` /
`residual note` / `MPI_Abort` / `ERROR STOP` / `CoLM_stop` / `unresolved dry isotope` / `no resolved store`；
原缺陷在 step 4 触发，修复后已越过 step 4 到约 step 70 且零诊断。日志仅有预期 warning（sediment grain
diameter、isotope fractionation experimental）与 flood material land/vapor/donor/residual ~1e-14–1e-16 舍入。
**G1–G7 校验未运行（作业取消），因此组合验收仍未通过、结论未定。**

## 可复核命令

```sh
# 在 7920 上（results 目录按 job id 定位）
python3 /tmp/compare.py <job>/results/runs/fixed-iso/case/history <job>/results/runs/base-iso/case/history
python3 /tmp/compare.py <job>/results/runs/fixed-iso/case/restart/2005-003-84600 \
                       <job>/results/runs/base-iso/case/restart/2005-003-84600
python3 /tmp/vardiff.py <job>/results/runs/fixed-all/case/restart/2005-002-00000 \
                        <job>/results/runs/base-all/case/restart/2005-002-00000
```
