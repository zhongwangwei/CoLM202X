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
| F02 非 BIF 内部单元排空 dry cleanup 删除真实质量 | 否（已修） | **已修复并本地验证**（结构与公式不变量；集成待授权） |
| F03 甲烷名义饱和列含气相但输送用纯液相容量 | 否（已修） | **已修复并本地验证**（容量公式/结构；集成待授权） |
| F05 NSS 叶水历史项错用当前源水组成 | 否（已修） | **已修复并本地验证** |
| F06 水汽扩散误用于非挥发溶质 | 否（已修） | **已修复并本地验证** |
| F07 冻土冰体积按液水密度换算 | 否（已修） | **已修复并本地验证** |
| F08 甲烷 6 个汇总诊断漏更新 | 否（已修） | **已修复并本地验证**（结构与加权断言；集成待授权） |
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

## F02：内部排空单元真实质量下传（已修复）

- 位置：`main/MATERIAL/GENERIC_WATER/MOD_MT_RiverLake.F90` `mt_riverlake_diag_accumulate`。
- 修复：dry-cell cleanup 区分“内部下游存在（`ucat_next(i)>0`）”与“终端/河口”：
  - 内部单元排空后剩余 `trc_mass + max(trc_inp_buf,0)` 不再记入 `trc_dry_drain`，而是累计到新的
    `trc_dry_fwd(itrc,i)`；本子步清理循环结束后用 `worker_push_data(push_ups2ucat, ..., mode='sum')`
    沿主河道下游边转发，并加到接收单元 `trc_mass`（若接收单元也干，则后续子步按同一规则继续转发，
    最终到湿单元或真实河口）；转发是内部转移，不记 `a_trc_out`/`trc_dry_drain`；
  - 终端单元（无下游 ucat）保留原有“计入 discharge 的 budgeted sink”语义。
  这样既不重设比值、也不删除真实库存，且不改变宿主水量步长。BIF/levee 路径不改动：转发只走主河道边，
  protected/levee 池仍按其原 donor-limiter 路径处理。
- 验证（`tests/material_transport/check_audit_p1_fixes.py`，8 项 F02 检查）：
  内部 drained mass 被缓冲而非记 sink；终端仍记 sink；转发使用下游 push；接收单元库存增加；
  buffer 分配/清零与 tracer 状态一致；3 单元链条例子（中间单元残留 .0002）转发后**总量守恒**且到达下游。
  另有 6 个改动文件 `gfortran -cpp -E` 预处理 rc=0。
- 兼容性：无 restart 新字段（`trc_dry_fwd` 为瞬时工作状态）；不改变水量/步长；终端 discharge 预算语义保持。
- 局限：未跑 7920 集成/守恒对比（测试暂停）；BIF/levee 路径的 drained-mass 转发未扩展（主河道修复优先）。

## F03/F08：甲烷二相容量与 6 个汇总诊断（已修复）

- F03：`MOD_MT_Methane_Physics.F90` `methane_tran` below-WT 分支在保留气相时改用等效容量
  `epsilon_t(j,1)=vol_ch4_storage+vol_gas/k_h_cc(j,1)`、`epsilon_t(j,2)=vol_aqu+vol_gas/k_h_cc(j,2)`，
  相对浓度用该容量（CH4/O2 同步）；无气相时退化为原液相容量。公式不变量：Henry 平衡
  `Ctot=Cgas*(θg+KH*θl)` 下液相浓度 `=KH*Cgas` 在每层一致 → 零净扩散；旧的纯液相容量会产生伪梯度。
- F08：在同一聚合点补齐 6 个诊断：
  `o2_decomp_depth`、`o2_oxid_depth`、`o2_aere_depth` 按 sat/unsat 面积加权；
  `o2stress`、`methane_stress` 按定义取面积平均（dimensionless 可用/需求比）；
  `methane_ebul_tot = ebul_sat*f + ebul_unsat*(1-f)`，湖泊分支随后仍覆盖为 lake 值。
- 验证（`check_audit_p1_fixes.py`）：F03 3 个结构断言 + 2 个数值不变量；F08 6 个聚合模式 + 湖泊覆盖顺序，
  全部 PASS；6 个改动文件预处理 rc=0。
- 兼容性：无 restart/history schema 变化；湿层（无气相）数值行为与旧代码一致；干/饱和含气列会改变
  扩散强度（这是修复目标，物理原因：容量包含气相存储）。
- 局限：未跑 7920 甲烷集成（相平衡三组对照、sat/unsat 深度积分）——需测试恢复授权后执行；本轮以
  生产源码结构与公式不变量验证，不宣称集成通过。

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
