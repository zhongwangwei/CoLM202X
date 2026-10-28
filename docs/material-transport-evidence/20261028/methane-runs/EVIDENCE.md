# 甲烷 CROP+BGC+MT 真实构建与全球 2° pilot：证据报告（methane10）

任务：task-4（Team Lead: lead；执行: methane-pilot）。只覆盖 7920l 私有作业；不改仓库、不提交、不取消/重提他人任务。
所有重计算经 `~/.local/bin/run-7920 submit --agent codex`（显式只读 `--input-remote`）；本机只做只读分析与文本证据取回（不 fetch NetCDF）。
源码基线：`CoLM-constituent` HEAD **`e33d64f81508014b175592df4359474ede4712ab`**（生产树 = 交付 HEAD `10e9fbd2`；`git diff 10e9fbd2 e33d64f8 -- . ':(exclude)docs' ':(exclude)tests'` 为空）。

## 0. 结论速览（2026-10-29）

| 子项 | 判定 | 一句话 |
| --- | --- | --- |
| A 真实构建 `CROP+BGC+MaterialTransport+LULC_IGBP_PFT+VG+GridRiverLakeFlow` | **READY** | 首次真实 `make all` 成功；define.h 与三二进制 sha256 固定 |
| mkinidata（候选 landdata+运行时依赖，冷启动 1995-06-01） | **READY** | 4m12s / 15,714 core-s / 642MB restart，`CoLM Initialization Execution Completed` |
| B 3 天真实短跑（交付默认 **hybrid** 口径） | **NOT_READY** | step 2 被交付版 fail-closed 守卫 `methane inundation exceeds host soil water` 中止；hybrid/dynamic_wtd 都不可运行（详见 §5） |
| B′ 3 天短跑（替代 **wetwat** 口径） | **READY（有口径变更声明）** | 144 步全部完成、`CoLM Execution Completed`、worst host-water 失配平台 0.083<0.20、CH4 通量有限非零、平衡残差 ~1e-21 |
| C 水稻 3 天 / 30 天 | **NOT_READY** | 水稻路径可执行（irrigation-on restart + 36 步），但被同一 host-water 守卫在 step 36 中止；wetwat 口径下稻田 finundated≡0（无水稻 CH4 物理） |
| D 两年全球 W1（1995-01-01 03:00 → 1996-12-31 18:00，35,070 步） | **配置+预算 READY / 口径 NOT_READY** | 拆两个年度作业（各 ≈9.1–11.5 h @64 CPU，均 <24 h；runbook `run_year.sh`）；hybrid 口径需先修守卫或改口径；wetwat 口径可跑 |
| 与旧 BD 甲烷结果可比性 | 仅参照 | 旧作业=TRACER 谱系 + spH off/lake off/rice off；本实现=MATERIAL 谱系 + 新默认参数，不可逐位/逐量比较 |

## 1. 合规、源码与构建身份

- 只在 `/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/methane10/` 写文件；仓库只读（`git archive e33d64f8` 冻结到 `new/`，686 文件，job 内 `sha256sum -c source.sha256` rc=0）。
- **注意（冻结后的仓库推进）**：本任务冻结源之后，仓库又多了两个提交（`9f7c9a52` isotope residual 修复、`d6347f6c` 文档，HEAD 变为 `d6347f6c`）。
  `git diff --name-only e33d64f8 HEAD` 的生产文件变化仅 `main/MATERIAL/GENERIC_WATER/MOD_MT_SoilWater.F90`、`main/MATERIAL/MOD_MT_Hist.F90`、`main/MATERIAL/MOD_MT_Vars.F90`（isotope dry-pool 修复），
  **`main/MATERIAL/FAMILY_GAS/`（甲烷实现）逐字节未变**；diff 中对 "methane" 的命中为 0 → 本报告的甲烷构建/运行结论对当前 HEAD 同样适用，但构建哈希对应 e33d64f8。
- 构建宏（实际 define.h）：`MaterialTransport+BGC+CROP+LULC_IGBP_PFT+vanGenuchten_Mualem_SOIL_MODEL+GridRiverLakeFlow+USEMPI`；
  `LULCC/CaMa_Flood/Campbell_SOIL_MODEL/URBAN_MODEL` 未定义（**不可 LULCC**、Campbell 与 MT 互斥）。
- 构建命令：`.github/workflows/create_defineh.bash GRID LULC_IGBP_PFT URBANOFF vanGenu CaMaOFF BGCON CROPON MTON LULCCOFF` +
  `make all FF=mpif90 NETCDF_INC=$(nf-config --includedir) NETCDF_LIB=$(nc-config --libdir) FOPTS='-fdefault-real-8 -fdefault-double-8 -ffree-form -O2 -cpp -ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium'`；
  工具链 conda gfortran 15.2.0 / MPICH hydra 4.3.2 / netCDF4-python 1.7.2。构建 ~3 min（-j8 并行 + 增量串行补齐）。
- 二进制 sha256（作业 A 实测，后续所有 job 复用同一份只读 build）：
  - `colm.x     3f83278c468e21bb64595ab91b0e892b1b6e7b4ebbc3de51067be675e122f90f`
  - `mkinidata.x e47a908b537a02deb9bf34e1cf3220aed3075e70fec770349d97f75eef81ce7e`
  - `mksrfdata.x 3141da8773db5afa635720388d71fdad74356b201192fcbc60bd09f668650dba`
  - `define.h    017da94293fffa8ba6e07d288fb6684ecce77a9df17e9164bd30b5cb87f7728a`
- CI 覆盖：`.github/workflows/TestCaseLists` 中 CROPON 用例全为 MTOFF、MTON 用例全为 CROPOFF → **CROP+BGC+MT 组合没有 CI 覆盖**；本构建是该组合的第一份真实编译证据。

## 2. 输入与运行环境（候选包 UNVERIFIED_CANDIDATE）

- 候选根 `/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`：
  landdata 40,429 文件 2.36 GB，forcing 16 文件 24.2 GB（CRUJRA v2.5 0.5°, 1995/1996 各 8 变量 365d 文件）；传输一致，科学身份未验证。
- landdata 实测：全球 2°、仅 2005 年；`landpatch` 144 块 / **37,347 patches**、`landpft` 137 块（含 CFT `nrice=61`/`nirrig_rice=62`）；
  `soil/2005/methane_ph_patches_*`（pH 3.0–10.0）、`lake_soilc_patches_*`、`topography/2005/fsatmax|fsatdcf_patches_*`。
- 运行时 `/media/zhwei/data02/zhwei/CoLMruntime/`（soilstate 34.2MB、cnsteadystate 46.2MB、routing 240MB、crop 4 文件、ndep/nitrif 可读；与 BD 历史算例同内容**未证实**）。
- 窗口：pilot 1995-06-01 00:00 → 06-04 00:00（144 步）/ 06-01 → 07-01（1,440 步），全部落在 1995 forcing 文件内部，**不需要 15 个缺失边界文件**。

## 3. 作业台账（全部 agent=codex，私有作业）

| # | job id | 用途 | 终态 | 结论 |
| --- | --- | --- | --- | --- |
| A0 | `codex-af7f004937a344e1b9c9ac1fd269d333` | 提交 | failed | 提交时在 source 目录内写日志导致 manifest 变化，未启动（已修） |
| A | `codex-9db71e90e32846b59b3cb0d7ecd29137` | 构建+mkinidata+3天 hybrid | failed/1 | **构建 PASS、mkinidata PASS**；3 天 run step 2 被 host-water 守卫中止 |
| A2 | `codex-a214907adda84f1ab08c5d849e5c7167` | 参数准备 | failed | mkparam 包装丢第 6 参数（host tolerance），未跑模型（已修） |
| A3 | `codex-bdde92b67e8c4fdea9caf52e7eca937f` | 3 天 hybrid tol 0.20 | failed | step 3 中止，deficit/pore 0.2578；worst 0.038→0.142→0.258 |
| D1 | `codex-90b01488287a414e96cc983d0002b1c7` | 1 天模式诊断（3 变体） | failed/1（脚本门禁） | hybrid 0.95/36 步中止；dynamic_wtd 0.20/1 步中止；**wetwat 0.20 全 48 步 PASS** |
| D2 | `codex-da8dce57739443e68401a0932d5bc6f8` | 1 天诊断 2 | succeeded（脚本） | hybrid+floodfeedback step 3 中止；dynamic_wtd 0.95 step 36 中止；rice 因 restart 缺 `irrig_rate` 启动失败 |
| P3 | `codex-c5057cf13eaf43c89ca89e6153eae4dc` | 最终 pilot（4 阶段） | **succeeded** | mkinidata(irrigation) PASS；rice 1 天 36 步被守卫中止（同 patch/同数值）；**wet3d 3 天 PASS**；**wet30d 30 天 PASS**（§4.2-4.3） |

## 4. 实测证据与数字

### 4.1 mkinidata（作业 A，64 rank，1995-06-01 冷启动）

- wall **4m12s**；User 15,551.99 s + Sys 162.28 s（≈15,714 core-s）；CPU 6228%；maxRSS **703 MB/rank**；cgroup 峰值 **28.1 GB**；
- 日志 `CoLM Initialization Execution Completed`；restart 树 **642 MB / 1,430 文件**（含 const）；
- 读入 40,429 landdata 文件；`Check vector data` 显示 `plant date value for rice: 4–365`、`plantdate_pfts` 部分缺测 fallback、`fert nitro` fallback 正常。
- 另：作业 P3 的 **M1 mkinidata with `DEF_USE_IRRIGATION=.true.`** PASS，wall 5m16s，restart **661 MB**（含 `irrig_rate` 变量）——水稻 irrigation 物理的真实运行前置。

### 4.2 3 天短跑（wetwat 替代口径，P3-W3，1995-06-01 → 06-04）

- **144/144 步完成，rc=0，`CoLM Execution Completed`，wall 855 s（含 init 与最终 restart 写盘）**；
- 日志：`CH4 inundation mode: wetwat -> scheme=1 wetwat_override=T dry_unsat=T`；`PH-LOAD patch pH: active=T, patches=37347, range=[3.00,10.00]`；
- **host-water 最坏失配 8.30e-2（< 0.20 容差）且全程平台化**（1 天时同样 0.083 → 稳定，不是增长型不一致）；
- CH4 变量（history `case_hist_material_1995.nc`）：关键通量/产生/氧化/平衡残差均有限；3 天窗口的 global 面积分见 `report.json`（冷启动 3 天，仅量级 sanity，不外推科学结论）；
- 内存 cgroup 峰值 ≈33 GB；单 rank maxRSS ≈1.7 GB；输出：host hist ≈24 MB/slice、material hist ≈0.6 MB、unitcat ≈18 MB/slice、restart ≈1.0–1.6 GB（1,430 文件）；
- 步速（日志采样）：48 步/≈110.5 s → **0.434 步/s = 2.30 s/step @64 rank**（当时服务器有其他 agent 作业、load≈100；core 成本 ≈147 core-s/step）。

### 4.3 30 天全球 2° run（wetwat 替代口径，1995-06-01 → 1995-07-01，1,440 步）

- **PASS**：1,440/1,440 步，rc=0，`CoLM Execution Completed`，最终时间戳 `1995-06-30-84600`；
- wall **53m12s（3,191.6 s）**，User 179,540 s + Sys 343.7 s（CPU 5,640%）；内存 cgroup 峰值 **37.1 GB**（单 rank maxRSS 1.68 GB）；
- host-water 最坏失配**全程稳定在 8.30e-2**（< 0.20 容差，无中止）；河网水量平衡 8.7e-03 m³（全球存量 1.9e12 m³）；
- CH4 通量（30 天平均、面积分，冷启动仅量级 sanity）：`f_methane_surf_flux_tot` **3.755 mol/s**（≈1.9 Tg/yr 当量，**不可外推**）、
  wetland 3.764、soil −0.0089、lake −5.7e-5、**rice 0**；prod 4.206、oxid 0.0405；`f_methane_balance_residual` 30 天最坏 |值| **4.52e-21**（数值闭合）；
- `f_methane_finundated` 面积平均 0.0322（最大 1.0）；输出：host hist 587.6 MB + material 25.1 MB + unitcat 426.4 MB（DAILY 诊断输出）、restart 1.63 GB/1,430 文件；
- 步速：运行中 log 采样 **1,440 步 / 2,650 s = 0.543 步/s = 1.84 s/step @64 rank**（逐区间中位 1.90 s/step）；单步 core 成本 **≈118 core-s/step**。
- 结论：**全球 2° 冷启动在 wetwat 口径下 30 天可稳定运行**（这是本任务唯一可运行的甲烷口径）；1 天/3 天/30 天三档均通过。

### 4.4 水稻 pilot 的实测边界

- `enable_rice_paddy=.true.` 需要 CROP 编译（本构建满足）；paddy 识别接受 `nrice/nirrig_rice` CFT + flood/paddy 灌溉方式。
- irrigation-on restart 已真实生成（M1）；`dynrice1d` 模型**成功启动并运行 36 步**，说明水稻路径、CROP+irrigation 数据链真实可用；
  但随后被同一 host-water 守卫中止（patch 172 layer 2，finundated 0.9939、vtot 0.000476 m（0.476 mm）、pore 14.68 mm、deficit/pore 0.9615 > tol 0.95）——
  与**不含水稻**的 `dynamic_wtd` 1 天测试完全同一位置/同一数值 → **水稻本身不是触发因素**。
- wetwat 口径（Dynamic_Wetland=false）下土壤 patch finundated≡0 → 稻田无淹水驱动，`f_methane_surf_flux_rice` 积分=0（1 天实测）
  → **wetwat 口径不能做水稻**。水稻需要 scheme 6（hybrid/dynamic_wtd），而 scheme 6 在全球冷启动被守卫阻断。

## 5. 关键缺陷：交付版 host-water fail-closed 守卫 vs 全球 2° 冷启动（根因分析）

源码：`main/MATERIAL/FAMILY_GAS/MOD_MT_Methane_Physics.F90:1096-1124`（由 P5 提交 `4d033986` 2026-10-28 引入；该提交自述只在 7920 的
NEC/AMZ 1 个月、`wetwat` 与 `routing levee+FB` 配置上验证过，**从未验证全球 hybrid/dynamic_wtd 冷启动**）。

判据：对每个甲烷活跃 patch 的**每一层**要求 `finundated × 孔隙体积 ≤ 宿主土壤水体积(wliq+wice)`，容差 `DEF_METHANE%host_water_tolerance`（默认 0.05）；
超出即 `CoLM_stop`。守卫只计土壤孔隙水，**不计 wetwat/wdsrf 表面水、也不计含水层 wa**；而 `Dynamic_Wetland=true` 的宿主水位 zwt 来自含水层+土壤水，
与土壤孔隙水本来就可能不一致（湿地/淹水 patch 表层可未饱和）。

实测矩阵（同一次冷启动 1995-06-01，1 天=48 步，64 rank）：

| 变体 | 模式/附加 | tolerance | 结果 | 观测 |
| --- | --- | --- | --- | --- |
| A-b3 | hybrid（交付默认） | 0.05（交付默认） | FAIL **step 2** | deficit/pore 0.0509/0.1417/0.0505/0.0572；tolerated 警告最大 2.9e-2 |
| A3 | hybrid | 0.20 | FAIL **step 3** | 0.2578 > 0.20；worst 0.038→0.142→0.258（增长型） |
| hyb_t095 | hybrid | 0.95 | FAIL **step 36** | worst→**0.9615**；失配持续增长到 ≈1 个孔隙体积（routing 洪水不在宿主土壤水，flood feedback 关闭） |
| hybff_t020 | hybrid + `GridRiverLake_FloodFeedback=.true.` | 0.20 | FAIL **step 3** | 0.2577；flood feedback 早期步内未消除不一致（失败 patch 为湿地类） |
| dyn_t020 | dynamic_wtd | 0.20 | FAIL **step 1** | 0.3037（patch 663 layer 1） |
| dyn_t095 | dynamic_wtd | 0.95 | FAIL **step 36** | worst→0.9466；**同一 patch 172 layer 2、同一数值 0.9615** |
| wet_t020 | **wetwat（Dynamic_Wetland=false）** | 0.20 | **PASS 48 步** | worst→0.083 平台；河网水量平衡 2.4e-2 m³（全球存量 1.9e12 m³） |

结论：
1. **所有 scheme-6 模式（hybrid/dynamic_wtd）在候选全球 2° 冷启动下都无法通过交付默认守卫；把 tolerance 调到 0.95 也只是把中止从 step 2 推到 step 36，
   失配是增长型（→≈1 pore volume），tolerance 合法上限 <1 → 调参数不是可行修复。**
2. 根因是交付版新增守卫与全球 2° 冷启动初态/宿主流域状态的**结构性不一致**（守卫不计表面水/含水层；Dynamic_Wetland 的 zwt 与土壤孔隙水不一致），
   不是输入缺项、不是编译错误、也不是水稻/irrigation 引起。
3. 可用替代：**wetwat 口径（Dynamic_Wetland=false）** 通过 1 天/3 天（30 天见 §4.3），代价是丢掉 routing floodplain 与稻田淹水物理。
4. 若两年全球必须保留交付默认 hybrid 物理，需要（超出本任务写权限，需 Lead/用户决策）：
   a) 修源码（守卫计入 wetwat/wdsrf/wa，或按宿主可支撑量对 finundated 封顶而不是中止）；或
   b) 加 host spinup（先用 host-only 运行把 Dynamic_Wetland 水位与土壤水积分到平衡，再冷启动/续跑甲烷）；或
   c) 明确接受 `host_water_tolerance≈0.95` 并声明守卫近似失效（不推荐，hybrid 仍会到 step 36 中止）。

## 6. 参数口径冻结（与交付默认的逐项差异）

基准：`new/run/standard_methane_parameter.nml`（HEAD `e33d64f8` 的 `run/standard_methane_parameter.nml`，其余字节不变）。允许差异只有：

| 项 | 交付默认 | pilot 用值 | 用途/理由 |
| --- | --- | --- | --- |
| `DEF_METHANE%enable_rice_paddy` | `.true.` | baseline/wet = `.false.`；rice 诊断 = `.true.` | 全球两年不带水稻；水稻 pilot 单独打开 |
| `DEF_METHANE%only_wetland` | `.false.` | `.false.`（rice_wetonly 变体 `.true.`） | 与历史/生产口径一致（湿地+土壤） |
| `DEF_METHANE%ch4_history_vars` | `'core'` | 自定义 28 项列表 | 保留 flux 分项 + 守恒诊断（balance residual/clip credit/finundated/zwt） |
| `DEF_METHANE%host_water_tolerance` | `0.05` | `0.20`（wet3d/wet30d）；`0.95` 仅诊断 | 交付默认在冷启动 step 2 中止；0.20 是 wetwat 实测可运行且失配平台 0.083；0.95 只用于证明 scheme-6 不可救 |
| `DEF_METHANE%inundation_mode` | `'hybrid'` | **`'wetwat'`（可运行口径）**；诊断另测 hybrid/dynamic_wtd | 交付默认 hybrid 被守卫阻断（§5）；wetwat 是唯一实测通过的模式 |
| 主 namelist `DEF_USE_Dynamic_Wetland` | 默认 `.false.` | wetwat 口径 `.false.`（与模式要求一致）；诊断 hybrid/dynamic_wtd 用 `.true.` | wetwat 模式要求 `.false.` |
| 主 namelist `DEF_USE_IRRIGATION` | 默认 `.false.` | 全球 wetwat 口径 `.false.`；水稻诊断 `.true.`（mkinidata+colm 同步） | 水稻需要灌溉物理；irrigation-on restart 单独生成 |
| `DEF_WRST_FREQ` | 历史 `YEARLY` | pilot `MONTHLY`（避免每日 restart 大 I/O；末步 restart 总会写） | 性能优化，不影响物理；两年 W1 配置仍用 `YEARLY` |
| 其余物理参数 | HEAD 默认 | **未改** | q10、f_methane biome、redoxlag、z0_methane_prod、anoxia、aere、bubble、lake…… 全部 delivery 默认 |

构建宏冻结：`MaterialTransport+BGC+LULC_IGBP_PFT+vanGenuchten+GridRiverLakeFlow`，CROP 打开（水稻需要，全球也兼容），**LULCC 关闭**。

## 7. 两年全球 W1 配置与预算（不启动；详见 `PLAN_global2yr_7920.md`）

- 窗口 W1：1995-01-01 **03:00** → 1996-12-31 **18:00**，**35,070 步**，候选包外文件请求 **0**（生产时间索引例程实测；比原窗口少 9h 模拟时间，需 Lead 确认实验定义）。
- 拆两个年度作业（单作业上限 24h）：Y1 1995-01-01 03:00→1996-01-01 00:00 = **17,514 步**（含 mkinidata 一次）；Y2 1996-01-01 00:00→1996-12-31 18:00 = **17,556 步**（复用 Y1 的 YEARLY restart）。
- 实测预算（64 CPU，基于 30 天 run 实测）：**≈1.84 s/step（≈118 core-s/step）**；含固定开销上界 124.9 core-s/step（1.95 s/step）；保守取 2.3 s/step。
  年度作业 Y1/Y2 ≈ **9.1–11.5 h**（含负载/保守区间），两年 ≈ **18–23 h 墙钟**，拆两个 <24 h 作业；建议 64 CPU / 128 GB / timeout 72,000–86,400 s。
- 内存：cgroup 峰值 **≈33 GB**（64 rank）→ 建议 128 GB 上限即可；output：history ≈0.3–0.6 GB/年 + material，restart YEARLY ≈1.6 GB/checkpoint。
- **口径阻塞**：W1 若按交付默认 hybrid 跑会立刻撞 §5 守卫；需先决定 hybrid 修复/口径方案，或明确用 wetwat 口径（可跑，但物理不同）。

## 8. 未解决 / 需要 Lead/用户决策

1. **交付默认 hybrid 全球运行被守卫阻断（本任务最高优先级发现）**：需选择 (a) 修源码（守卫计入 `wetwat/wdsrf/wa` 或对 finundated 按宿主支撑封顶）；(b) host spinup 方案；(c) wetwat 口径变更；或 (d) 其它残余策略。未修前两年全球 hybrid NOT_READY。
2. **wetwat 口径的科学代价**：丢失 routing floodplain 与稻田淹水物理；与旧 TRACER 结果/交付默认参数不可逐位/逐量比较。
3. **水稻路径**：在守卫修复前不可完成（3 天/30 天均无法过 step 36）；irrigation-on restart 资产已具备，修复后可立即复用重试。
4. **早期窗口证据缺口**：候选包精确 1995-1996 原窗口缺 15 个 CRUJRA 边界文件（21.15 GiB，来源已实测存在于 BD，待授权）；W1 是无传输替代（少 9h）。
5. 运行时依赖身份未证实；2005 landdata 用于 1995/1996 是历史既有偏差；冷启动无 spinup，CH4 池未平衡 → 本轮所有 CH4 数字都不是科学结论。
6. 7920 上同时有其他 agent 的作业（pilot 期间 load≈70–130，未取消/未重提）；年度作业排期前需查负载。

## 9. 证据索引（本目录）

| 内容 | 路径 |
| --- | --- |
| 冻结源清单 / 来源 | `source.sha256`、`source_versions.txt`、`new/` |
| 作业台账（含所有失败原样记录） | `job-attempts.txt` |
| 构建证据 | `evidence/jobA-failed-03/results/build/`（define.h、binaries.sha256、make 日志、defines 日志） |
| 3 天 hybrid 失败证据 | `evidence/jobA-failed-03/worklogs/key-evidence.txt`、`evidence/jobA3-failed-06/worklogs/key-evidence.txt` |
| 1 天模式诊断（D1/D2） | `evidence/job-diag1/`、`evidence/job-diag2/`（phase JSON + worklogs/key-evidence.txt） |
| 最终 pilot P3（mkini irrigation / rice attempt / wet3d / wet30d） | `evidence/job-pilots3/`（fetch 后） |
| 步速采样 | `evidence/steprate/samples.tsv` |
| 机器可读汇总 | `results.json`（构建哈希、逐作业/逐阶段数字、步速、两年预算、verdicts；由 `make_results.py` 生成） |
| 复现脚本 | `run.sh`（构建+mkinidata+3天）、`pilots3.sh`（最终 pilot）、`diag.sh`/`diag2.sh`（模式诊断）、`run30.sh`（30 天）、`run_year.sh`（W1 Y1/Y2 年度作业 runbook 脚本）、`lib.sh`、`mkparam.py`、`memwatch.py`、`validate.py` |
| 两年全球设计/runbook | `PLAN_global2yr_7920.md` |
