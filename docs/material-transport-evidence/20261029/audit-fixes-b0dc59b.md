# 审查基线 b0dc59b7 的逐项修复记录（2026-10-29）

- 仓库/分支：`zhongwangwei/CoLM202X` `feature/constituent-transport`
- 审查基线：`b0dc59b70f49a61a5fc837ca5c22648909e6305f`
- 本轮工作起点 HEAD：`2a9573b9`（本地；按要求未拉远端）
- 附件 `CoLM-constituent-audit-b0dc59b.zip` 不可用 → 按任务描述自行构造最小复现。
- **F04、K01 按用户指示本轮不处理**（“不是我们管的”）；不修改 `create_defineh.bash`（已回退该脚本修改，
  CATCHMENT 限制改放在受控 `include/define.h` + `MOD_MT_Defs.F90`）。

## 状态总表

| 项 | 最新代码中是否仍存在 | 本轮状态 |
| --- | --- | --- |
| F01 湖泊多子步只交付最后子步示踪径流 | 是（未改） | **未修复**（下一步） |
| F02 非 BIF 内部单元排空 dry cleanup 删除真实质量 | 是（未改） | **未修复**（下一步） |
| F03 甲烷名义饱和列含气相但输送用纯液相容量 | 是（未改） | **未修复**（下一步） |
| F05 NSS 叶水历史项错用当前源水组成 | 否（已修） | **已修复并本地验证** |
| F06 水汽扩散误用于非挥发溶质 | 否（已修） | **已修复并本地验证** |
| F07 冻土冰体积按液水密度换算 | 否（已修） | **已修复并本地验证** |
| F08 甲烷 6 个汇总诊断漏更新 | 是（未改） | **未修复**（下一步） |
| F04 CATCHMENT 侧向水同物质输送 | 是 | **按指示不处理** |
| K01 SinglePoint history 续跑覆盖 | 是 | **按指示不处理** |

## F05：NSS 叶水历史项（已修复）

- 位置：`main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90`，`tracer_transpiration_nss_ratio`。
- 旧式：`prev_w * prev_p * (prev_e - delta_x)`；新式：`prev_w * (prev_delta_b - delta_x)`；
  有记忆用存储的 `prev_delta_b`，无记忆初始化 `prev_bulk_delta = delta_x`（历史项 0）。
- 验证（`tests/material_transport/check_mt_frac_physics.py`，生产例程逐字提取 + gfortran stub 编译）：
  - 首步对 `prev_delta_e` 完全不变（旧代码 3 项 FAIL 中第 2 项即此项）；
  - 有记忆时改变 `prev_delta_b` 必须改变结果（通过）；
  - 无记忆时 `prev_delta_b` 被忽略（通过）；
  - 恒定强迫 400 步收敛，稳态蒸腾比 reldiff < 5%。
- 兼容性：无记忆路径与旧代码一致（历史项本来就为 0）；有记忆首步会改变，属修复目标。
- 局限：未在 7920 做完整模型集成（测试暂停）；未复现附件中 −22.38003316‰ 的具体数值
  （附件不可用，改用不变量验证）。

## F06：水汽扩散物种限制（已修复）

- 位置：`MOD_MT_Frac.F90` 新增 `tracer_isotope_soil_vapor_equivalent_diffusivity` /
  `tracer_isotope_snow_vapor_equivalent_diffusivity`（家族检查 `tracer_is_isotope`）；
  `MOD_MT_SoilWater.F90` **三个入口**（土壤、普通雪、湿地雪）改调新 helper 并在入口条件加
  `.and. tracer_is_isotope(itrc)`。
- 验证：新测试中同位素 soil/firn 扩散 > 0；非挥发溶质两类 helper **精确 0**；饱和冻土两类均 0。
  旧代码在饱和冻土上仍给 8.6e-14（F07 的旧行为）。
- 兼容性：溶质只失去非物理的水汽路径；液相/水动力路径未改；无 restart/schema 变化。

## F07：冻土冰体积密度（已修复）

- 位置：`MOD_MT_Frac.F90` `tracer_soil_vapor_equivalent_diffusivity`：
  `theta_ice = ice_mm/(rho_ice*dz)`（`rho_ice=917 kg/m³`，与宿主 `denice` 一致），
  并对 `air_porosity < 1e-8` 视为饱和（数值零），避免舍入残余孔隙产生虚假扩散。
- 验证：液水 10 kg/m² + 冰 32.095 kg/m²、dz=0.1 m、porsl=0.45 时空气孔隙精确零、D_vap 精确 0；
  未冻、零液水边界保持正通路。旧代码同例给出 8.6e-14（`/1000` 换算）。
- 物理影响：冻结土壤的水汽等效扩散在接近饱和时显著减小（预期改变冷季同位素剖面）；
  不涉及质量源汇，守恒形式不变。

## 本地验证命令与结果

```sh
PATH=/opt/homebrew/bin:$PATH python3 tests/material_transport/check_mt_frac_physics.py
# PASS: 15 checks; FAIL: 0 checks   (gfortran 16.1, -O0 -fcheck=bounds)
```

```sh
# 旧代码复现（git show 2a9573b9:main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90 提取旧例程）
# PASS: 6 checks; FAIL: 3 checks  (exit 1) —— 失败项：
#   F07 saturated/pore-filling: D_vap != 0（8.6e-14）
#   F05 first step is independent of prev_delta_e（旧代码依赖 prev_delta_e）
```

```sh
# CATCHMENT 限制（本轮从 create_defineh.bash 移到受控 define.h + MOD_MT_Defs）
# 临时 define.h: #define MaterialTransport + #define CATCHMENT
gfortran -cpp -E -I <tmp>/cat main/MATERIAL/MOD_MT_Defs.F90   # rc=1，报出限制消息
# 仅 #define MaterialTransport
gfortran -cpp -E -I <tmp>/grid main/MATERIAL/MOD_MT_Defs.F90  # rc=0
```

## 尚未完成（下一轮）

- F01、F02、F03、F08：代码未改；已定位入口（F01：`MOD_MT_SpecialPatches` 子步清零/赋值、
  `CoLMDRIVER` 湖泊子步；F02：`MOD_MT_RiverLake`/`MOD_MT_RiverLakeCoupling` 的 `trc_dry_drain`；
  F03：`MOD_MT_Methane_Physics` `jwt_sat=0`/`Crel=Ctot/theta_l`；F08：6 个汇总诊断的 sat/unsat 聚合）。
- 全量构建/CI/模型集成验证：需测试恢复授权后执行（本轮只做本地纯函数/片段测试）。
- 7920 侧尚未运行任何新作业。
