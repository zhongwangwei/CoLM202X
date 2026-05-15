# CoLM CH4 对齐 CTSM 计划

日期：2026-05-15
目标仓库：`/Users/zhongwangwei/Desktop/Github/CoLM202X`
基准：CTSM `master` = `e1f563b8345400d33421263f05b237f13e6b0e3d`；可选扩展：CTSM `wetlands` = `10690e577ed0f2a93d9a02c85f2f146babe892ae` / PR #3893。

## 1. 目标

把 CoLM CH4 从“土壤/湿地 CH4 主链基本可运行”推进到“与 CTSM master 的 CH4 物理过程、状态变量、诊断和守恒检查基本对齐”。

分两层：

1. **CTSM master 对齐**：lake CH4、CO2 diagnostics、balance check、DynamicColumnAdjustments/面积变化、参数/restart/history 整理。
2. **CTSM wetlands / PR #3893 可选对齐**：surface-layer transport、WT 到地表时 unsat oxidation、snow+h2osfc 面积协调、surface water drainage 新方案。

非目标：CTSM 对齐主线不做 acetoclastic/hydrogenotrophic 拆分、完整 NetCDF paramfile 重构，也不把 rice CH4 混入自然湿地/lake 主链。rice/paddy CH4 作为独立扩展参考 DSSAT-CSM，见 `methane-dssat-csm-reference-20260515.md`。

## 2. 当前证据

- CoLM driver 只有 `igas_ch4 > 0` 且 `patchtype==2/0` 时运行 methane，lake `patchtype==4` 当前不进入 CH4 driver：`main/CoLMDRIVER.F90:244-316`。
- CoLM 有 `lake_soilc` state/restart 骨架，但缺少 lake production/oxidation/diffusion/ebullition 变量族：`main/MOD_Tracer_Methane_State.F90:128,265,470,501`。
- CoLM 主流程已有 production / oxidation / aerenchyma / ebullition / transport / annual update：`main/MOD_Tracer_Methane_Physics.F90`。
- CoLM inline column balance 存在，但不是 CTSM 式独立 init/gridcell balance：`main/MOD_Tracer_Methane_Physics.F90:797-803`。
- CoLM 有 `use_ch4_sif` 与 `bgc_anoxia_limits_decomp` 双开关：`main/MOD_Tracer_Methane_Const.F90:219-220`，使用点：`main/MOD_Tracer_Methane_Physics.F90:1005-1010`。
- PR #3893 是 `wetlands` 分支内容，不是当前 CTSM master 内容；应作为扩展目标。

## 3. 决策原则

1. 先对齐物理闭环，再对齐工程外观。
2. 小步可验证，每阶段必须能跑单点/小域 smoke test。
3. 保持 CoLM 架构，不机械复制 CTSM 类型系统。
4. 防止 BGC anoxia 与 CH4 SIF 双重压低 SOMHR。
5. master 与 wetlands 分离：先 master，再 PR #3893。

## 4. 工作包

### WP0：建立 baseline 与测试壳

目的：改动前冻结当前 upland/wetland CH4 行为。

工作：
- 建立 CH4 单点运行说明。
- 记录 wetland `patchtype==2` 与 upland `patchtype==0` baseline 输出。
- 检查 flux、prod、oxid、aere、ebul、diff、mass balance。

验收：
- CH4 打开时能跑过单点/小域。
- CH4 关闭时模型路径不变。
- 有 baseline 输出可回归比较。

### WP1：lake CH4 driver gate

目的：让 lake patch 在明确开关下进入 methane。

工作：
- 修改 `main/CoLMDRIVER.F90`：`DEF_METHANE%allowlakeprod .and. patchtype==4` 时允许运行 methane。
- 默认 `allowlakeprod=.false.` 保持现有 lake 零贡献。
- 确认 `main/MOD_Tracer_Methane_Const.F90` namelist 可控。

验收：
- `allowlakeprod=.false.`：lake 不产生 CH4，现有结果不变。
- `allowlakeprod=.true.`：lake 进入 CH4 driver，不 NaN、不崩溃。

### WP2：lake CH4 物理闭环

目的：补齐 CTSM master 的 lake production / oxidation / diffusion / ebullition。

工作：
- `main/MOD_Tracer_Methane_Physics.F90`：新增 lake-specific 分支，避免 lake 走陆地 root/aerenchyma/nitrification 路径。
- `main/MOD_Tracer_Methane_State.F90`：新增变量：
  - `methane_prod_depth_lake`
  - `methane_oxid_depth_lake`
  - `methane_ebul_depth_lake`
  - `methane_surf_ebul_lake`
  - `methane_surf_diff_lake`
  - `conc_methane_lake`
  - `conc_o2_lake`
- 接入 `lake_soilc`、`lake_decomp_fact`、`q10lakebase`。
- 补 allocation/deallocation/restart/history。
- 找 CoLM lake temperature / lake resistance 可用接口；如果没有，先实现最小 conductance 近似并标注。

验收：
- `lake_soilc>0` 且 `allowlakeprod=.true.` 时 lake production 非零。
- `lake_soilc==0` 或 `allowlakeprod=.false.` 时 lake CH4 为零。
- lake CH4 column balance 闭合。
- lake 输出变量单位、维度、restart/history 正确。

### WP3：CO2 diagnostics

目的：补齐 CTSM 的 CO2 诊断闭环。

工作：
- `main/MOD_Tracer_Methane_State.F90` 新增：
  - `co2_decomp_depth_sat/unsat`
  - `co2_oxid_depth_sat/unsat`
  - `co2_aere_depth_sat/unsat`
- `main/MOD_Tracer_Methane_Physics.F90` 在 production/oxidation/aerenchyma 中记录 CO2 诊断。
- 第一阶段只做 diagnostics，不直接改 BGC carbon pool / NEE。

验收：
- CO2 diagnostics 与 CH4 production/oxidation 量纲自洽。
- 默认 NEE 不改变，除非显式打开 coupled correction。
- 输出可用于 CH4-C / CO2-C budget 检查。

### WP4：独立 balance check

目的：把 inline balance 升级为可复用、可定位的检查框架。

工作：
- 抽出 `methane_column_balance_check`。
- 增加 lake balance check。
- 增加 initialization balance check。
- 必要时新增 `main/MOD_Methane_Balance.F90`。

验收：
- upland/wetland/lake 均有 residual 输出。
- residual 阈值和单位明确。
- 报错包含 patch id、patchtype、sat/lake、layer/time。

### WP5：DynamicColumnAdjustments / 面积变化

目的：LULCC 或 patch/column 面积变化时 CH4 库存不凭空丢失。

工作：
- 找 CoLM 现有土地利用/patch 面积变化入口。
- 对 CH4/O2 concentration、column stock、diagnostic flux、lake_soilc 分别定义面积变化规则。

验收：
- 面积变化前后 CH4 stock 守恒。
- 新增 patch 初始化合理，无 NaN。
- 消失 patch 库存处理有明确诊断。

### WP6：参数与 namelist 整理

目的：清理 dead parameter，明确 active/lake-only 参数。

工作：
- `main/MOD_Tracer_Methane_Const.F90` 标注 active / lake-only / diagnostic / deprecated。
- 决定 `q10lake` 是独立参数还是 `q10methane*1.5`。
- 确保 `allowlakeprod`、`replenishlakec`、`lake_decomp_fact`、`q10lakebase` 可控。
- 更新 `run/standard_ch4_parameter.nml`。

验收：
- 参数有默认值、单位、来源说明。
- lake disabled 时 lake 参数不影响 soil/wetland CH4。

### WP7：可选对齐 CTSM wetlands / PR #3893

前置：WP0-WP6 通过后再做。

工作：
- surface-layer methane transport。
- WT 到表面时 unsat oxidation 薄层。
- snow+h2osfc fractional area 协调。
- surface water drainage matric-potential-gradient 方案。

验收：
- 新功能默认关闭时结果与 master-aligned 版本一致。
- 打开后，在 WT near surface、融雪、surface water case 中 flux 变化物理可解释。
- 无负面积、负孔隙、负浓度、NaN。


## 4b. 外部参考：DSSAT-CSM rice Methane

DSSAT-CSM 的 CH4 模块是 rice/paddy 管理场景参考，不是 CTSM wetland/lake 对齐对象。已整理到 `methane-dssat-csm-reference-20260515.md`。可借鉴点包括：

- decomposition carbon substrate cap；
- flood/drainage and stored-CH4 release；
- alternate-electron-acceptor redox buffer；
- crop-root-density-driven plant transport；
- DSSAT/Arah-style output decomposition（production/consumption/plant/ebullition/diffusion/leaching/storage）。

集成原则：默认关闭、仅 rice/paddy 专用、复用 CoLM 现有 transport/oxidation/aerenchyma 框架，不直接移植 DSSAT solver。

## 5. 测试矩阵

| Case | patchtype | CH4 | lake | 目的 |
|---|---:|---:|---:|---|
| upland off | 0 | off | no | 确认无 CH4 时不影响原模型 |
| upland on | 0 | on | no | 检查 unsat/sat 混合与 oxidation |
| wetland on | 2 | on | no | 检查主 CH4 flux |
| lake disabled | 4 | on | off | 确认默认 lake 为 0 |
| lake enabled no carbon | 4 | on | on, `lake_soilc=0` | 确认无碳不生产 |
| lake enabled carbon | 4 | on | on, `lake_soilc>0` | 确认 lake CH4 非零且守恒 |

边界条件：WT at surface、`finundated=0/0.01/0.5/1`、snow depth > 0、surface water > 0、frozen top layer、very low O2、high CH4 ebullition。

## 6. 里程碑

| 里程碑 | 内容 | 可交付物 |
|---|---|---|
| M0 | baseline + test harness | 当前 CH4 单点说明、baseline 输出 |
| M1 | lake driver gate | patchtype 4 可控进入/跳过 CH4 |
| M2 | lake CH4 core | lake production/oxidation/transport/ebullition + diagnostics |
| M3 | CO2 diagnostics | CO2 变量族 + budget 输出 |
| M4 | balance framework | 独立 balance check，覆盖 soil/lake |
| M5 | dynamic area handling | 面积变化下 CH4 state 规则 |
| M6 | parameter cleanup | namelist 与参数说明整理 |
| M7 | optional PR #3893 | surface-layer 等新物理可选开关 |

## 7. 推荐执行顺序

1. WP0 baseline。
2. WP1 lake driver gate。
3. WP2 lake core。
4. WP4 balance check，尽早插入守恒约束。
5. WP3 CO2 diagnostics。
6. WP6 参数整理。
7. WP5 DynamicColumnAdjustments。
8. WP7 可选 wetlands / PR #3893。

## 8. 完成标准

- `igas_ch4 <= 0` 时 CH4 完全 no-op。
- lake disabled 时 upland/wetland CH4 不破坏。
- `allowlakeprod=.true.` 时 lake patch 可以产生、氧化、释放 CH4。
- lake CH4 有独立 diagnostics 与 restart/history 支持。
- soil/wetland/lake CH4 每步 mass balance residual 在阈值内。
- CO2 diagnostics 与 CH4 production/oxidation 自洽。
- 新增参数有默认值、单位、namelist 或明确硬编码理由。
- 至少通过单点 upland/wetland/lake smoke tests。

## 9. 分工建议

- Lane A：lake driver + lake physics，负责 `CoLMDRIVER.F90`、`MOD_Tracer_Methane_Physics.F90` lake 分支。
- Lane B：state/restart/history/diagnostics，负责 `MOD_Tracer_Methane_State.F90` 和输出注册。
- Lane C：balance/test harness，负责守恒检查与单点测试。
- Lane D：parameter/namelist/docs，负责 `MOD_Tracer_Methane_Const.F90`、`run/standard_ch4_parameter.nml`、运行说明。

## 10. 停止条件

- 若 CoLM 无可靠 lake conductance 来源，先停在 lake production + diagnostic-only，不宣称 lake 完全对齐。
- 若 baseline wetland/upland 回归显著漂移，暂停扩展，先修主链。
- 若 PR #3893 hydrology 改动影响面超过 CH4 模块，单独立项。
