# 甲烷真实构建、pilot 与两年全球准备（methane10 证据归档）

状态记录：2026-10-29。本目录归档甲烷子任务的**小型证据**：真实 CROP+BGC+MT 构建、mkinidata、
交付默认 hybrid 口径被 fail-closed 守卫阻断的根因、替代 wetwat 口径 1/3/30 天运行、两年全球 W1 配置与预算。
**完整原始报告在 7920 作业结果目录与本机私有 `.tmp/ct/methane10/`；这里只保留摘要性证据。**

## 结论速览（截至本归档）

| 子项 | 判定 | 摘要 |
| --- | --- | --- |
| 真实构建 MT+BGC+CROP+LULC_IGBP_PFT+VG+GridRiverLakeFlow | **READY** | `make all` 成功；define.h 与三二进制 SHA256 固定；CROP+BGC+MT 无 CI 覆盖，这是首份真实编译证据 |
| mkinidata 冷启动 1995-06-01 | **READY** | 4m12s / 15,714 core-s / restart 642MB / `CoLM Initialization Execution Completed` |
| 交付默认 **hybrid** 全球口径 | **NOT_READY（根因已定位并修复）** | 旧代码 fail-closed 守卫 `methane inundation exceeds host soil water` 全球冷启动 step 2/3 中止；根因：finundated 是独立水位估计，与土壤孔隙水不是同一量；修复 `3778404f` 用 host-supported saturation 封顶分配替代 abort（分配本身守恒、不造水） |
| 替代 **wetwat** 口径 | **READY（有口径变更声明）** | 1 天/3 天/30 天真实运行全部 END；30 天 wall 53m12s、cgroup 峰值 37.1GB、水量平衡 8.7e-3 m³（全球存量 1.9e12 m³）、CH4 balance residual 4.5e-21、通量有限非零 |
| 水稻甲烷 | **未完成（等待 guard 修复后重跑）** | rice 可启动（irrigation-on restart），但被同一守卫在 step 36 以同一 patch/数值中止；水稻不是触发因素；wetwat 口径下稻田 finundated≡0、rice 通量=0，不能替代水稻验收 |
| 两年全球 W1 | **配置+预算 READY，运行待 guard 修复验证** | W1 1995-01-01 03:00 → 1996-12-31 18:00，35,070 步；拆 Y1 17,514 步 / Y2 17,556 步（YEARLY restart），每作业 9.1–11.5h @64 CPU（实测 1.84 s/step、峰值内存 33–37GB），远低于 24h 上限 |

## 关键实测数字

- 构建（作业 `codex-9db71e90…`）：`define.h=017da942…`；`colm.x=3f83278c…`、`mkinidata.x=e47a908b…`、`mksrfdata.x=3141da87…`。
- mkinidata：4m12s、15,714 core-s、maxRSS 703MB/rank、cgroup 28.1GB、restart 642MB；irrigation-on 版 5m16s、restart 661MB。
- hybrid 守卫诊断（1 天 48 步）：默认 tol 0.05→step2 中止；0.20→step3（0.258）；0.95→step36（失配增长到 0.96≈1 孔隙体积）；
  hybrid + flood feedback 仍 step3；dynamic_wtd 0.95→step36（0.9466，同一 patch/layer）；水稻同 patch/数值 → 结构性不一致，
  调大容差不可救（上限 <1）。
- wetwat 30 天（作业 `codex-c5057cf1…`）：1,440 步、53m12s、`CoLM Execution Completed`、worst host-water 0.083 稳定、
  水量平衡 8.7e-3 m³、CH4 balance residual |max| 4.5e-21、CH4 30 天面积分 tot=3.755 mol/s（wetland 3.764、soil −0.0089、
  lake −5.7e-5、rice 0；冷启动量级 sanity，不是科学结论）。
- 步速：wet30d 1.84 s/step @64 rank（≈118 core-s/step）；3 天短跑 2.3 s/step（启动开销占比高）。
- W1 预算：Y1 17,514 步 ≈9.1h（保守 11.3h）；Y2 17,556 步同量级；两年合计 ≈18–23h，需两个年度作业串行。

## 参数口径（相对交付默认的差异，全部记录）

`inundation_mode`：全球运行目标为交付默认 hybrid（修复后）；诊断替代 wetwat（可运行口径，丢 routing floodplain/稻田淹水物理）。
`enable_rice_paddy`：全球 false / 水稻 true；`only_wetland=false`；`use_spatial_ph=true`；`allowlakeprod=true`；
`ch4_history_vars` 自定义 28 项；`host_water_tolerance` 记录值 0.20（诊断另用 0.95）；`Dynamic_Wetland`：hybrid=true / wetwat=false；
`irrigation`：全球 false / 水稻 true（mkinidata 同步开）；pilot WRST=MONTHLY、生产 W1 用 YEARLY。

## 边界

- 候选输入仍为 `UNVERIFIED_CANDIDATE`；运行时依赖与 BD 历史算例内容身份未证实；1995/96 使用 2005 landdata 是历史惯例，需在正式报告中声明。
- 旧 BD 甲烷结果（TRACER 谱系、spH off、lake off、rice off）**只作参照，不是正确基准**；新实现 MATERIAL 谱系不可逐位/逐量比较。
- wetwat 口径的科学代价必须随结果一起声明；不得把 wetwat 的通过写成 hybrid 或水稻通过。
- 本归档不是甲烷最终验收。2026-10-29 已完成的实测：hybrid 1/3 天 PASS、水稻 1/3 天 PASS（非零、sign 记录）；
  **30 天水稻与两年完整窗口运行于 10:03 按用户要求暂停/取消**（用户已发现明显 bugs，等待清单）。
- 用户同日批准：从 `bd:` 补齐 15 个缺失边界文件（≈21.15 GiB）以使用**完整窗口 1995-01-01 00:00 →
  1996-12-31 21:00（35,088 步）**、两年全球用**hybrid 交付默认口径**、只用个人 fork 留档；
  **补数据与 Y1/Y2 均尚未启动**（暂停中）。

## 文件

| 文件 | 内容 |
| --- | --- |
| `EVIDENCE.md` | 甲烷子任务完整报告（构建、mkinidata、守卫诊断、wetwat 运行、W1 预算、未解决项） |
| `results.json` | 机器可读结果汇总（作业终态、测量值、参数与判定） |
| `PLAN_global2yr_7920.md` | 两年 W1 运行配置与预算（含 run_year.sh 设计） |
| `run_year.sh` / `mkparam.py` | 年度运行 runbook 与参数生成脚本（7920 私有作业使用） |
| `job-attempts.txt` | 所有提交/失败/重试记录 |
| `source_versions.txt` | 该阶段源码版本说明 |
