# 甲烷 guard 修复后的短跑实测证据（hybrid / 水稻，1 天与 3 天）

来源作业：`codex-6b46c0ceb804474d8abea8ec7d661a1c`（v2 build + hybrid 1/3 天 + rice 1/3 天）与
`codex-683d918f970e4730b6e8f1347463d7da`（v2r：修好 phase-JSON 命名后重跑 rice 验收）。
源码快照：`3778404f`（甲烷 guard 修复后）。原始 JSON 名称：`h1d_full.json`/`h3d_full.json`/
`r1d_full.json`/`r3d_full.json`、`report.json`、`rice_gate.json`、`a_fix_present.json`。

## 机器门禁（全部 PASS，rc=0）

| 运行 | 配置 | 步数 | END | guard abort | water-balance violation | 最大水量平衡残差 | host-water worst |
| --- | --- | --- | --- | --- | --- | --- | --- |
| h1d | hybrid（交付默认，rice off，host_water_tolerance 默认 0.05） | 48 | ✓ | **0** | **0** | 0.045 m³ | 0.9615（capped，report 行记录） |
| h3d | hybrid | 144 | ✓ | **0** | **0** | 0.069 m³ | 0.9615 |
| r1d | dynamic_wtd + enable_rice_paddy + irrigation | 48 | ✓ | **0** | **0** | 0.044 m³ | 0.9615 |
| r3d | 同上 | 144 | ✓ | **0** | **0** | 0.068 m³ | 0.9615 |

## 通量/守恒（作业报告实测）

| 运行 | CH4 balance residual | 通量积分 | 水稻通道 |
| --- | --- | --- | --- |
| h1d | \|max\| 1.8e-21 | tot 3.450 mol/s（wetland 3.459、soil −0.0085、lake −5.1e-5、rice 0） | 未启用 |
| h3d | \|max\| 2.5e-21 | tot 3.111 mol/s | 未启用 |
| r1d | \|max\| 1.8e-21 | rice 积分 −4.63e-4 mol/s | prod_rice 1.889e-5、oxid_rice 9.96e-4（**非零**；冷启动首日为小净汇，sign 记录） |
| r3d | \|max\| 2.5e-21 | rice 积分 −4.45e-4 mol/s | prod_rice 5.34e-5、oxid_rice 6.79e-4、`f_methane_surf_flux_rice` n_nonzero=2526 |

## 边界

- 冷启动 1995-06-01、全球 2°；通量量级只作“有限/非零/守恒”的 sanity，不是科学结论或观测验证。
- 1/3 天短跑不等于 30 天/两年验收：30 天水稻取消前只到 336/1440 步（见 `../partial/README.md`）；
  两年完整窗口（用户已批准补数据、hybrid 默认、35,088 步）尚未启动。
- `host_water_worst=0.9615` 是封顶分配后仍记录的最坏结构不一致（不再 fatal）；不是“通过阈值”的证明。
