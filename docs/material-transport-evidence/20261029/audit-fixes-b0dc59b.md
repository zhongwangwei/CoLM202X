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
| F01 湖泊多子步只交付最后子步示踪径流 | 否（已修） | **已修复并本地验证** |
| F02 非 BIF 内部单元排空 dry cleanup 删除真实质量 | 是（未改） | **未修复**（根因/下一步见下） |
| F03 甲烷名义饱和列含气相但输送用纯液相容量 | 是（未改） | **未修复**（根因/下一步见下） |
| F05 NSS 叶水历史项错用当前源水组成 | 否（已修） | **已修复并本地验证** |
| F06 水汽扩散误用于非挥发溶质 | 否（已修） | **已修复并本地验证** |
| F07 冻土冰体积按液水密度换算 | 否（已修） | **已修复并本地验证** |
| F08 甲烷 6 个汇总诊断漏更新 | 是（未改） | **未修复**（根因/下一步见下） |
| F04 CATCHMENT 侧向水同物质输送 | 是 | **按指示不处理** |
| K01 SinglePoint history 续跑覆盖 | 是 | **按指示不处理** |

## F01：湖泊/水体子步示踪径流交接（已修复）

- 位置：`main/MATERIAL/GENERIC_WATER/MOD_MT_SpecialPatches.F90`（`mt_glacier_patch`/`mt_waterbody_patch`）、
  `MOD_MT_Conservation.F90`、`MOD_MT_Vars.F90`、`main/CoLM.F90`。
- 根因：每个 1800 s 子步的 `mt_save_storage`/`mt_waterbody_patch` 清零 `trc_rnof_step`，末尾又覆盖为
  最后子步的 `trc_rnof`；`MOD_MT_Conservation` 的守恒快照也每步清零。河网 hand-off 只读一次，
  于是只收到最后子步质量。
- 修复：
  - `MOD_MT_Vars` 新增 `mt_add_rnof_step`（子步累加，忽略非 generic-water 与越界/负值）和
    `mt_reset_rnof_step`（只清 generic-water）；`trc_rnof_step` 语义改为“本主步累计 hand-off 量”，
    瞬时状态、不进 restart；
  - 特殊斑块子步改为 `CALL mt_add_rnof_step`，删除两处子步清零；
  - `MOD_MT_Conservation` 删除每次守恒快照的清零；
  - `CoLM.F90` 在 `CoLMDRIVER` 之前**每主步一次** `CALL mt_reset_rnof_step()`（在 land 物理与
    river coupling 读取之前）。
- 验证（`tests/material_transport/check_rnof_handoff.py`，生产例程逐字提取 + gfortran）：
  1/2/3 子步恒流、变流量、最后子步为零：`river_received == land_accumulated`；
  provider-owned 材料不受影响；负值/越界被忽略；reset 语义正确；另有 5 项源码结构断言
  （无子步清零、两处 `mt_add_rnof_step`、CoLM.F90 主步 reset 顺序等）。本地 **14/14 PASS**。
- 兼容性：宿主水量、能量、BGC 未改；`trc_rnof_step` 非 restart 状态；provider-owned（甲烷）路径未涉及。
- 局限：未跑 7920 全模型集成（测试暂停）；变流量下的水量积分仍沿用宿主原有 `rnof(last)*deltim`，
  本轮只让**示踪质量**按子步精确累加，不引入隐式水量变化。

## 未修复：F02（已定位，未改动代码）

- 位置：`main/MATERIAL/GENERIC_WATER/MOD_MT_RiverLake.F90` 的 `mt_riverlake_diag_accumulate`
  dry-cell 分支（`volwater <= trc_v_dry_off`）：把 `trc_mass + trc_inp_buf` 正值记入
  `trc_dry_drain` 并清零；`MOD_MT_RiverLakeCoupling.F90` 再把 `trc_dry_drain` 汇总进 `trc_mass_dis`。
- 复现条件（按任务描述）：非 BIF 串联单元 `[100,1,100] m3`、比值 `[.002,.0018,.002]`、出流 `[1,2,0] m3/s`、
  dt=1 s → 中间单元排空，`.0002 R*m3` 被删除、下游少收同样质量。
- 为什么本轮不改：正确修复要求把内部排空单元的剩余真实质量沿**同一水通量**送到真实下游接收单元，
  需要与 `mt_riverlake_route_calc`/donor limiter 的拓扑（非 BIF 与 BIF/levee 分支、双向流）统一处理并保证
  正性/守恒/不重复；在缺少 7920 集成验证（测试暂停）时盲改风险高（可能把质量送错单元或破坏 BIF 路径）。
- 下一步（需恢复测试授权后）：先加“直接调用生产 transport 例程”的三单元复现单元测试（当前没有文件级
  Fortran 源码可独立编译 RiverLake 模块），再改 limiter 使同一步入流可被继续输送（pass-through），
  或把 dry cleanup 剩余质量按 outflux 边界分给下游；验证恒比值保持、异质混合正确、剔除数值 sink 后收支闭合，
  并记录迭代次数/步数/耗时。F02 目前记为**未修复**。

## 未修复：F03 与 F08（甲烷，已定位，未改动代码）

- F03：`MOD_MT_Methane_Physics.F90` 中 `jwt_sat=0` 固定、`split_ch4_o2_phases` 用两相容量
  `Ctot=Catm*(theta_g+KH*theta_l)` 初始化/相分配，而 below-WT 输送使用 `Crel=Ctot/theta_l`；
  CH4 和 O2 都要统一为含气相的等效容量（若以液相浓度为未知量则为 `theta_l+theta_g/KH`），并同步矩阵、
  回写与边界。复现：两个 0.1 m 层、298.15 K、无反应/植物/冰雪，`theta_l=.2, theta_g=.3, jwt=0`
  基线虚假扩散 `~1.9198e-13 mol/m2/s`；`theta_l=.5, theta_g=0` 为对照。本轮未改。
- F08：6 个汇总诊断 `o2_decomp_depth`/`o2_oxid_depth`/`o2_aere_depth`/`o2stress`/`methane_stress`/
  `methane_ebul_tot` 在 sat/unsat 聚合处未更新（总鼓泡只有湖分支另行更新）；需要按变量定义
  （面积平均 vs 总供给/需求比）在统一聚合点补齐，并检查 soil/rice 再聚合、湖泊路径和输出掩码。
  本轮未改。
- 未改原因：两者都需要甲烷模块的完整依赖（BGC 宿主）与集成/最小耦合测试；当前测试暂停且本地无
  netCDF/MPI 全量构建环境，无法给出满足“质量残差小不能替代相平衡/深度积分验证”的验收证据。

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

- F02、F03、F08：代码未改；已定位入口（F02：`MOD_MT_RiverLake`/`MOD_MT_RiverLakeCoupling` 的 `trc_dry_drain`；
  F03：`MOD_MT_Methane_Physics` `jwt_sat=0`/`Crel=Ctot/theta_l`；F08：6 个汇总诊断的 sat/unsat 聚合）。
- 全量构建/CI/模型集成验证：需测试恢复授权后执行（本轮只做本地纯函数/片段测试）。
- 7920 侧尚未运行任何新作业。
