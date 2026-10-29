# 甲烷 host-water fail-closed 守卫的根因修复（提交 `3778404f`）

状态记录：2026-10-29 DeepSeek 接续。本文记录交付默认 **hybrid** 全球甲烷口径被 fail-closed 守卫阻断的根因、
修复内容、设计依据与验证状态。**修复后 hybrid/水稻/两年运行的实测验收仍在进行；在那些结果终态前，
不能声称甲烷验收通过。**

## 缺陷

真实构建（MT+BGC+CROP+LULC_IGBP_PFT+VG+GridRiverLakeFlow）与 mkinidata 均成功，但交付默认 hybrid 全球
2° 冷启动 3 天运行在 step 2 中止：

```
ERROR: methane inundation exceeds host soil water: <patch> <layer> finundated vtot pore_volume
       deficit/pore = 0.0509 ... tolerance = 0.05
CALL CoLM_stop('methane inundation exceeds host soil water')
```

1 天诊断（48 步）显示是**结构性不一致、不是数值噪声**：默认 tol 0.05 → step2；tol 0.20 → step3（0.258）；
tol 0.95 → step36 且失配增长到 0.96（≈1 孔隙体积）；hybrid + flood feedback 仍 step3（0.2577）；
dynamic_wtd 0.95 → step36（0.9466，同一 patch/layer）；水稻与不含水稻在同一 patch/同一数值中止。
只有 wetwat（`Dynamic_Wetland=false`，finundated≡0）能通过——但 wetwat 丢掉了 routing floodplain/稻田淹水物理，
稻田 CH4 恒为 0，不能替代默认口径或水稻验收。

## 根因（源码级，`MOD_MT_Methane_Physics.F90` 的水淹分配块）

- 守卫比较 `finundated × 孔隙体积` 与 `wliq+wice` 土壤水，**只计土壤孔隙水**。
- 但 `finundated` 在 hybrid/scheme 6-7 中是**独立由水位（zwt）S 曲线/河网洪泛比例计算的面积分数**；
  与其对应的水可能存在于塘水（wdsrf）、湿地/含水层（wetwat/wa）或 routing 存储里，不在土壤孔隙水统计内。
- 因此 `deficit/pore` 可以合法地接近 1；把两个独立估计当作必须一致并 fail-closed，会让交付默认 hybrid 在全球
  冷启动下不可用。**调大容差不能解决**（失配上限 <1，且失配随步数增长）。
- 关键：同一代码块里已有 `host_water_scale = min(pore/vtot, 1/finundated)` 的**封顶分配**，原注释已证明
  它使面积加权 host 水精确守恒（saturated+unsaturated 分配之和等于 host 总量，不造水、不丢水）；
  守卫只是“有效性策略”，不是守恒所必需。

## 修复（`3778404f`）

把该 fail-closed abort 改为**记录 + host-supported 封顶**：

- 超过 `DEF_METHANE%host_water_tolerance` 时不再 `CoLM_stop`，而是记录 `host_water_resid` 到
  `host_water_max_resid`，首次出现时打印一次 WARNING（含 patch/layer、deficit/pore、tolerance，并注明
  “capped to host-supported saturation (host water conserved)”）。
- 随后的封顶分配保持原逻辑：`host_water_scale = min(pore/vtot, 1/finundated)`，面积加权 host 水精确守恒；
  不修改 `wliq_soisno/wice_soisno` 的 host 状态。
- `methane_host_water_report` 仍按原机制在所有 worker 上 reduce 并打印最坏失配，因此该结构不一致**不会被隐藏**。
- 未改任何物理参数、未放宽 `host_water_tolerance` 的语义（它仍是宽容差，不再决定 fatal）、未改 CH4 产生/氧化/输运。

**性质声明**：这是策略修复（把“独立估计不一致 → 中止”改为“按 host 可支撑饱和度封顶 + 记录”），
不是“调参通过”；正确性依赖封顶分配的面积加权守恒（代码内已有注释证明，需实测水量平衡验证，见下）。

## 验证状态（部分通过；30 天水稻与两年运行已按用户要求暂停）

- [x] 真实构建 `CROP+BGC+MT+LULC_IGBP_PFT+VG+GridRiverLakeFlow`：PASS（首份该组合真实编译证据）。
- [x] mkinidata：PASS（4m12s、15,714 core-s、restart 642MB、END banner）。
- [x] 修复后 hybrid 1 天 PASS（48 步，guard abort 0、water-balance violation 0、host-water worst 0.0383→0.9615
  由报告行记录、CH4 residual |max| 1.8e-21、通量有限非零 `codex-6b46c0ce…`）。
- [x] 修复后 hybrid 3 天 PASS（144 步，同 gates；CH4 residual 2.5e-21；`codex-6b46c0ce…`）。
- [x] 修复后水稻（dynamic_wtd+rice+irrigation）1 天/3 天 PASS：rice 通量非零（积分 −4.63e-4 / −4.45e-4 mol/s，
  冷启动首日为小净汇、sign 已记录；prod_rice/oxid_rice 非零；`codex-683d918f…`）。
- [ ] 30 天水稻：作业 `codex-923e86aa…` 于 2026-10-29 09:29 启动，**10:03 按用户要求取消** → 未完成。
  取消前 partial：336/1440 步（7 天），零 fatal/诊断，water balance 1.6e-2–2.3e-2 m³，host-water worst 0.9615
  （capped、report 行记录），110 条 warning 全为 tolerated host-water；**不构成 30 天验收通过**。见
  [methane-runs/partial/README.md](methane-runs/partial/README.md)。
- [ ] 两年完整窗口（用户已批准补 15 文件、hybrid 默认、完整 35,088 步）：**补数据与 Y1/Y2 均未启动**
  （10:03 暂停；W1 少 9h 方案不再采用）。
- [x] 用实测水量平衡残差证明“不造水”：hybrid/水稻短跑 routing 水量平衡 |残差| ≤0.068 m³（存量 1.9e12 m³）、
  violation 0；但**未逐 patch 检查 0 水饱和列**（列为自审候选问题 #2）。

以上由 `methane-pilot` 在 7920 执行（task-5），证据写 `.tmp/ct/methane10/`；结果出来后在
`methane-runs/README.md` 与状态文档更新。**在终态通过前，不把 hybrid、水稻或两年全球写成已完成。**

## 备选与边界

- 若修复后仍有其他甲烷守卫失败（如 rice component fraction、BGC 合约等）：保留失败、按根因范围报告，
  不得绕过或改判成功；不得把 wetwat 替代口径写成 hybrid/水稻通过。
- wetwat 口径的 1/3/30 天实测（已 PASS，最坏 host-water 0.083、CH4 balance residual 4.5e-21）保留为
  “可运行替代口径”的证据与科学代价声明，但默认口径以本修复后的 hybrid 结果为准。
- 两年 W1 窗口较原始 1995-1996 少 9h（零新增边界文件）；实验定义变更需在正式结果中声明。
- 候选输入仍 `UNVERIFIED_CANDIDATE`；旧 BD 甲烷结果只作参照。
