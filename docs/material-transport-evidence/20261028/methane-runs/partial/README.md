# 30 天水稻运行取消前的 partial 证据（**不是验收通过**）

2026-10-29 09:29 CST 启动的 30 天水稻验收作业 `codex-923e86aad16940228b7def12bae54d4d`
（dynamic_wtd + enable_rice_paddy + irrigation，delivery-default tolerance，1995-06-01 → 1995-07-01，
1440 步，64 CPU/128GB，复用 v2 build + irrigation-on restart）在 **10:03 CST 按用户要求取消**
（用户“已发现明显 bugs”，先停全部测试）。以下为取消前已完成部分的只读提取，**不构成 30 天验收通过**。

| 项 | 取消前实测 |
| --- | --- |
| 进度 | **336 / 1440 步（1995-06-07 23:00，7 天）** |
| 致命/诊断 | `TRC_BAL` / `TRC_SIG` / `residual hard` / `ERROR STOP` / `MPI_Abort` / `Water balance violation` / `inundation exceeds` 计数 **0** |
| 水量平衡 | `Total water balance` 最近三条 2.3e-2 / 4.6e-3 / 1.6e-2 m³（全球存量 ~1.9e12 m³） |
| host-water | worst disagreement 0.8843 → 0.9466 → **0.9615** 孔隙体积（capped，不 abort；report 行如实记录） |
| 警告 | 110 条，全部为 `WARNING: methane host-water disagreement ...; tolerated`（无其它 warning 类型） |
| 输出 | `case_hist_1995.nc`（部分 6 月，10:03 仍在写）、`case_hist_material_1995.nc`；restart 至 1995-152 |
| 结论 | 水稻 30 天**未完成**；取消前没有失败迹象，但**不得写成通过**。重启后需从头/按 restart 续跑完整 30 天。 |

## 可复核命令

```sh
J=codex-923e86aad16940228b7def12bae54d4d
C=/media/zhwei/data01/compute_runs/codex/$J/work/pilot/rice30d/case
grep -a "TIMESTEP" $C/log.colm | tail -2
grep -a -c -E "TRC_BAL|TRC_SIG|residual hard|ERROR STOP|MPI_Abort|Water balance violation|inundation exceeds" $C/log.colm
grep -a "Total water balance" $C/log.colm | tail -3
grep -a "methane host-water: worst" $C/log.colm | tail -3
```
