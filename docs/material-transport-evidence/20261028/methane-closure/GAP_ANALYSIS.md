# 甲烷候选输入闭合、参数冻结与水稻/两年全球运行设计：缺口分析（GAP_ANALYSIS）

任务：task-2（只读）。源码基线：`CoLM-constituent` HEAD **`10e9fbd2`**（分支 `feature/constituent-transport`；分析期间其他 agent 又提交了两个纯文档提交 `2d5c15b7`、`74cb9175`，`git diff --name-only 10e9fbd2..74cb9175` 无任何 `.F90/.h/.nml` 变化）。
候选输入根：`/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`。
审计权威：job `codex-ad6b48e0cb0044a9bbf0ce2903195234`（succeeded/0，`results/audit.json` 与本地副本 sha256 相同：`1d6fb3cb…`）。
**本轮未启动任何模型、未传输/复制任何数据、未修改仓库、未动别人 7920 任务。**

## 0. 结论速览

| 子项 | 判定 | 一句话 |
| --- | --- | --- |
| A 精确 1995–1996 输入闭合 | **NOT_READY** | 15 个边界文件不在候选包；但**已实测确认存在的权威来源在 BD**（原传输的同一路径），待用户授权传输 |
| A′ 无传输替代（窗口改写） | **可用（需批准）** | 生产时间索引例程实测：起点 03:00 + 终点 18:00 时**外部文件需求为 0**；只损失 9 小时模拟时间 |
| B 候选包身份 | TRANSFER_VERIFIED / UNVERIFIED_CANDIDATE | 40,450 文件 / 26,618,457,915 B 在盘上逐组复核一致；2 namelist + 1 LSF 确认 |
| B′ namelist 兼容 | **旧 case namelist 不可直接用** | `&nl_colm` 有 15 个 TRACER 时代键在 HEAD 已删除（会 namelist 读失败）；forcing namelist 26 键全兼容 |
| B″ 甲烷参数冻结 | **未冻结（且已发现旧运行实际值与当前默认不同）** | 旧运行 log 实测 spatial pH=off、lake prod=off、rice=off；HEAD 标准参数文件是 pH=on、lake=on、rice=on |
| C 运行时依赖 | 可读，**身份未证实** | 8 类路径可读；CROP 运行文件 4/5（第 5 个仅 allocation=3 才需要）；与 BD 历史算例同内容无哈希证据 |
| D 水稻甲烷设计 | **设计就绪、可运行性未验证** | 候选 landdata 含 crop CFT/pH/lake soilc，够用；但 **CROP+BGC+MT 无任何 CI 组合**，需先真实编译 |
| E 两年全球估时 | **有权威基准** | 旧全球 2° 甲烷两年作业真实日志：720 ranks、mkinidata 4m23s、mode 3h05m、35,088 步、成功完成；折算 7920 @64 CPU ≈ **17 h/年、需拆两个年度作业** |

**总体：READY 条件未满足（精确窗口 NOT_READY）；但两件事已可立即推进**：(1) 用户授权后用现成 `bd:` 通道补 15 个文件即可闭合；(2) 或批准窗口改写方案 W1 立即零传输进入短跑。水稻 pilot 与两年全球的设计已给出，均需先做真实编译/短跑标定。

## 1. A：15 个缺失边界文件（精确到变量/年份/记录）

生产 reader 的时间请求逻辑（`setstampLB/setstampUB`）在审计快照与 HEAD `10e9fbd2` 的 5 个源文件**逐字节相同**，因此 15 文件缺口在 HEAD 上仍然成立。完整逐行请求表见 `evidence/boundary_request_table.md`；原始 23,407 条请求台账 `evidence/local-boundary-requests.csv`。

| 组 | 变量 | 文件年 | 记录（1-based） | 类型 | 触发模型时刻 |
| --- | --- | --- | --- | --- | --- |
| 7 | `dlwrf,pre,pres,spfh,tmp,ugrd,vgrd`（无 dswrf） | 1994 | 1460（年末最后一条 6h 记录 = 1994-12-31 18:00） | LB | 1995-01-01 00:00 第一步 |
| 8 | 全部 8 变量（含 dswrf） | 1997 | 1（1997-01-01 00:00） | UB | 1996-12-31 最后一步（6h 变量 21:00 / dswrf 18:00） |

- `dswrf` 1994 不需要：其 `offset=0`+forward/coszen，第一步 LB 直接取 1995 第 1 条。
- 每文件只读 1 条记录，但 reader 按绝对 record 下标寻址 → 安全单位是整年文件；按 1995 同变量尺寸外推，15 文件共 **22,706,115,044 B ≈ 21.15 GiB**（7 个 1994 = 9.87 GiB；8 个 1997 = 11.28 GiB）。
- **权威来源已实测存在**（只读 `ls`，经 7920l→`bd` 到原传输源机）：`bd:/share/home/dq013/zhwei/colm/data/CoLM_Forcing/crujra_v2.5/<var>/crujra.v2.5.5d.<var>.{1994,1997}.365d.noc.nc` 全部存在，尺寸与 1995 同变量文件**逐字节一致**（例如 dlwrf 1,513,741,012；tmp 1,513,741,024）；该目录覆盖 1849–2024。候选暂存的原始传输就是 `rsync -rlt … bd:/share/home/dq013/…`（见 staging 的 `*.rsync-status.json` argv），所以补传通道现成，**只差用户授权**。
- 7920 本机**没有**任何 CRUJRA 实际数据：只找到悬空符号链接 `/media/zhwei/data02/zhwei/training2026/CoLM_Forcing/crujra_v2.5 -> /share/home/dq013/...`（目标未挂载）；`/media/zhwei/data01`、`/home/zhwei`、`/media/zhwei/data02`（裁剪 CoLMrawdata 后）对 `*crujra*` 的 bounded find 除该链接和代码/文档名外无数据命中。JRA-3Q 在 7920 有 1947-09..2024-09 共 8 变量（每月一文件，~926 文件/变量），但**是不同强迫产品**，不能当成 CRUJRA 的边界补齐（除非用户决定换产品，属另一实验）。

### A′ 窗口改写实测（生产例程，非猜测）

用与审计相同的生产 `setstampLB/setstampUB` 例程（`probe_window.py`，gfortran 15.2 在 7920 上编译，源码哈希与审计一致）对四种窗口重放请求：

| 窗口 | start | end | 步数 | 候选包外请求 |
| --- | --- | --- | --- | --- |
| W0（原窗口，重放校验） | 1995-01-01 00:00 | 1997-01-01 00:00 | 35,088 | 15（7×1994 + 8×1997）与审计完全一致 ✅ |
| **W1（双端平移）** | 1995-01-01 **03:00** | 1996-12-31 **18:00** | 35,070 | **0** |
| W2（仅移起点） | 1995-01-01 03:00 | 1997-01-01 00:00 | 35,082 | 8（仅 1997） |
| W3（仅移终点） | 1995-01-01 00:00 | 1996-12-31 18:00 | 35,076 | 7（仅 1994） |

W1 损失 9 h 模拟时间（起点 3 h + 终点 6 h），换来零新增数据；W2 只缺 1997 的 8 个文件（11.28 GiB）。三种方案均需用户确认是否改变实验定义。

## 2. B：候选包来源、身份与参数/宏

### 2.1 传输一致性（在盘复核，非仅 manifest）

| 组 | manifest files/bytes | 7920 盘上 files/bytes | 一致 |
| --- | --- | --- | --- |
| forcing | 16 / 24,219,856,144 | 16 / 24,219,856,144 | ✅ |
| landdata | 40,429 / 2,362,078,324 | 40,429 / 2,362,078,324 | ✅ |
| baseline | 5 / 36,523,447 | 5 / 36,523,447 | ✅ |
| 合计 | **40,450 / 26,618,457,915** | **40,450 / 26,618,457,915** | ✅ |

- baseline 5 文件 = `config/` 3 文件（**2 份 namelist + 1 份 LSF**，确认不是 3 份 namelist）+ `history/` 2 个 nc。
- 来源（manifest 原值）：forcing 与 landdata 来自 BD `dq013` 的 `g_global_ch4_flat_2yr_20260527_200428` / `g_global_ch4_test`；传输方式 `rsync -rlt`，三组 `checksum dry-run clean`。
- baseline history 是 **TRACER 谱系** `hist_tracer`（64 变量、12 月），不是当前 MATERIAL 谱系的 `hist_material`。
- 身份仍是 **UNVERIFIED_CANDIDATE**：不证明唯一、版本正确或科学 ground truth；旧运行完成 ≠ 输入闭合 ≠ 正确。

### 2.2 候选 landdata 的实测内容（决定可用构建）

- 全球 2°（mesh：18×9 blocks，180×90 格；patch/pft 按 2005 年）。只有 **2005** 一年。
- patch/PFT/CFT：`landpatch` settyp 1..17；`landpft` settyp 0..**78**，其中含 **61/62**（`nrice`/`nirrig_rice`）；`pctpft/2005/pct_crops_*.nc` 144 块、`pct_pfts` 144 块。→ 该 landdata 由 **定义 CROP 的构建**生成（`N_PFT=15, N_CFT=64` 分支），且含稻 CFT。PFT vs PC 两种 subgrid 宏从文件无法唯一区分（要求构建侧冻结并用 mkinidata 验证）。
- 甲烷预处理向量齐全且 variable 名/维度与 HEAD reader 一致：
  - `soil/2005/methane_ph_patches_*`（patch；pH 5.0–8.2，fallback 6.2 在范围内）→ 当前 fail-closed reader 可读；
  - `soil/2005/lake_soilc_patches_*`（soil×patch，10×patch；mean 202.8、max 7697.5 gC/m³）→ 当前 reader 可读；
  - `topography/2005/fsatmax_patches_*`、`fsatdcf_patches_*` → Dynamic_Wetland 所需。
- 其他族：`LAI/SAI_patches01..12`+`LAI/SAI_pfts01..12`、`htop_patches/htop_pfts`、`lakedepth`、`soil`（215 个基名：vanGenuchten 参数、OM_density、csol、土壤质地/反射率等）、`topography`、`mesh`、`landelm/landpatch/landpft`。
- **缺少**：`cstructure`、`dbedrock`、`urban`、`riverlake(unitcatchment_regional)`、`lulcc` 族 → 冻结配置必须保持 `DEF_RS_CROWN_STRUCTURE=.false.`（默认）、`DEF_USE_BEDROCK=.false.`（默认）、不编译 URBAN_MODEL、`DEF_UnitCatchment_regional=.false.`（默认）、不编译/不启用 LULCC。
- 结论：候选 landdata 与 HEAD 的 mkinidata 输入布局（按块文件 `<base>_<block>.nc`，`ncio_read_vector`/`pixelset_load_from_file`）**形式匹配**；但这是"可读+布局匹配"，尚未做真实 mkinidata 运行验证。

### 2.3 namelist 兼容性（HEAD `10e9fbd2`，机器可读结果 `evidence/namelist_compat.json`）

- `CRUJRA_local_abs.nml`（`&nl_colm_forcing`）：26 个键 + 全部 `DEF_forcing%` derived-type 字段**全部存在** → 兼容，只需改 `DEF_dir_forcing` 绝对路径（BD→候选根）。
- 旧 case `&nl_colm`：65 键中 **15 键在 HEAD 已不存在**，Fortran namelist 读取会 `CoLM_Stop`（未知变量名），旧文件**不能直接运行**：
  `DEF_USE_METHANE_para, DEF_file_METHANE_para, DEF_METHANE_enable_rice_paddy, DEF_METHANE_only_wetland, DEF_METHANE_inundation_mode, DEF_USE_SpatialpH, DEF_TRACER_NUM, DEF_TRACER_NAMES, DEF_TRACER_TYPES, DEF_TRACER_MRAT, DEF_TRACER_REF_RATIO, DEF_TRACER_INIT_DELTA, DEF_TRACER_REACTIVE_PARAM_FILES, DEF_TRACER_REACTIVE_OPTIONS, DEF_TRACER_REACTIVE_DECAY_RATE`。
- 映射（改写要点）：甲烷开关迁到 `&nl_colm_material` 的 `DEF_USE_METHANE=.true.` + `DEF_METHANE_PARAM_FILE='<path>'`；CH4 私有开关迁到参数文件 `&nl_colm_methane_parameter` 的 `DEF_METHANE%…`；`DEF_USE_SpatialpH`→`DEF_METHANE%use_spatial_ph`；TRACER 注册整体删除（HEAD 由 `DEF_USE_METHANE` 驱动）。

### 2.4 构建宏要求（甲烷；HEAD 源码证据）

`MaterialTransport`（需 `vanGenuchten_Mualem_SOIL_MODEL`，否则 `define.h` 编译错误）+ `BGC`（`MOD_MT_Namelist.check_settings`：否则 `CoLM_Stop`）+ `LULC_IGBP_PFT` **或** `LULC_IGBP_PC`（同上）+ `GridRiverLakeFlow`（标准 `inundation_mode='hybrid'` 运行时硬性要求，`MOD_MT_Methane_Const::configure_methane_inundation_mode`）+ `CROP`（仅水稻模式需要；`MOD_MT_Methane.F90` 编译期检查）。**不能加 LULCC**（上游 BGC+LULCC 不支持）。
历史运行的宏无法从源码管理恢复：BD `script/CoLM202X` 在运行日附近的提交（`39a8b691` 2026-05-25、`1c7c365d` 2026-05-27 21:55）里 `include/define.h` 是模板（BGC/CROP/TRACER 均 `#undef`、`LULC_IGBP_PC`），说明实际构建用了本地未提交的 define.h；该工作树当前 HEAD 是 `080a09fb`（2026-09-17，晚于运行、且已改为 TRACER 谱系），**不是运行代码**。→ 宏必须重新冻结；候选 landdata 的 CFT 证据只能证明"生成 landdata 的构建带 CROP"，不能证明运行构建的完整宏集。

### 2.5 甲烷参数冻结：旧运行实际值 vs HEAD 默认（**关键差异**）

旧运行 log（BD，只读 grep）实测：`methane_registry_init: CH4 active at tracer index 1`；`CH4 inundation mode: hybrid -> scheme=6 wetwat_override=F dry_unsat=T`；`PH-LOAD patch pH: active=F, nonfallback=0/0, range=[6.20,6.20]`（**spatial pH 未生效，用常数 6.2**）；运行前最近提交 `39a8b691:run/standard_ch4_parameter.nml` 中 `DEF_METHANE%allowlakeprod = .false.`（**湖泊甲烷关闭**）；主 namelist 中 `enable_rice_paddy=.false.`、`only_wetland=.false.`。旧运行**没有** CROPOFF 的 "Fertilization is on when CROP is off" 警告 → 推断其构建 CROP=ON（未证明）。

HEAD `run/standard_methane_parameter.nml` 默认是：`use_spatial_ph=.true.`、`allowlakeprod=.true.`、`enable_rice_paddy=.true.`。→ **用 HEAD 默认参数文件不能复现旧结果**；冻结清单必须显式二选一（复现口径 vs 新默认），并记录在案。

### 2.6 候选 LSF 的资源/环境（仅作历史记录，不是 7920 基准）

`BSUB -n 720, mem=150G/task, span[ptile=48]`；BD 集群 mvapich2-2.3.7-intel + netcdf_parallel；`mkinidata -np 48` → `colm.x -np 720`；无 `--network`/容器，路径为 BD `/share/home/dq013/...`。

## 3. C：运行时依赖可读性（7920）

全部可读（`evidence/scan.json` → `sections.C_runtime`）：

| 路径 | 状态 | 备注 |
| --- | --- | --- |
| `CoLMruntime/soilstate.nc` | 34,249,742 B，可读 | dims month=12, lat=90, lon=180, layer=10 |
| `CoLMruntime/cnsteadystate.nc` | 46,243,066 B，可读 | dims lon=360, lat=180, soil=10 |
| `CoLMruntime/unitcatchment/grid_routing_data_15min.nc` | 240,217,739 B，可读 | 1440×720 河网；bifurcation/levee/reservoir 参数同文件 |
| `CoLMruntime/{nitrif,ndep,crop,fire,snicar}` | 目录可读 | 见下 |
| `CoLMruntime/crop/` | 4/5 文件 | `plantdt-colm-64cfts-rice2_fillcoast.nc`、`fertnitro_fillcoast.nc`、`fertilizer_2015soc.nc`、`surfdata_irrigation_method_96x144.nc` 均在；`surfdata_irrigation_allocation.nc` 缺失，但仅 `DEF_IRRIGATION_ALLOCATION==3` 需要（默认 1） |
| `CoLMrawdata/CFT/…`、`landtypes/`、`plant_15s/`(30,387)、`lai_15s_8day/`(21)、`soil/`(26)、`urban*`、`Forest_Height.nc`、`bedrock.nc`、`topography.nc`(41 GB)、`TWI.nc`(58 GB) 等 | 可读 | 仅重生成 landdata 时才用 |
| `CoLMrawdata/global_CFT_surface_data.nc` | **悬空符号链接**（指向 BD，未挂载） | 如要重生成 landdata 需用 `CFT/global_0.25x0.25.MOD2000_V4.5_CFT_lf-20260120.nc` 或修复链接 |
| 原始 pH `PHH2O1.nc` | 不在 CoLMrawdata | 但已在独立 staging：`…/data-stage/mt-ph-20261028-2053-1793191930/soil/PHH2O1.nc`，2,903,280,596 B，可读（PHH2O2 当前代码不需要） |

**身份未证实**：这些 7920 文件与 BD 历史算例使用的同名文件是否内容相同，没有任何哈希证据（BD 侧无备份哈希、本轮未传数据）→ 标注 UNPROVEN；不能声称"同输入复现"。3 个输入文件尺寸/时间维度（12×90×180×10 等）与 global 2° 配置量级一致。

另：7920 上搜索到的 "methane" 相关目录（`data-stage/methane-runtime-smoke-20261028`）是 **OpenBench/OBRef raw-engine 的甲烷 smoke 测试**，与 CoLM 甲烷运行无关；`training2026/ex04_bgc` 是 BD 训练算例的副本（见 §5），不是 7920 上的 CoLM 甲烷基准。

## 4. D：水稻甲烷最小真实运行（详见 `PLAN_rice.md`）

- 构建：`MaterialTransport+BGC+CROP+(LULC_IGBP_PFT|PC)+vanGenuchten+GridRiverLakeFlow`；**CI 无此组合**（`TestCaseLists` 89–93 的 MTON 全为 CROPOFF；CROPON 用例全 MTOFF）→ 必须先真实 `make all`。
- landdata：候选 2005 全球 2° 已含 `pct_crops`（CFT 分数、含 61/62 稻 CFT）、`methane_ph_patches`、`lake_soilc_patches`；mkinidata 可直接读该布局。稻 CFT 常数 `nrice=61`、`nirrig_rice=62`、`irrig_method_flood=3`、`paddy=4`（`MOD_Vars_Global`）。
- 最小真实运行窗口**不需要任何缺失边界文件**：取完全落在文件内部的窗口（如 1995-06-01 00:00 → 1995-06-30 24:00，水稻季）；起点/终点均在 1995 文件内部。
- 参数：`enable_rice_paddy=.true.`；`only_wetland` 建议 pilot 用 `.true.`（湿地+稻，排除普通 upland 土壤）；`inundation_mode='hybrid'`（需 `DEF_USE_Dynamic_Wetland=.true.`）；`use_spatial_ph` 与 `allowlakeprod` 按冻结决定（候选数据都支持）；`DEF_USE_IRRIGATION=.true.` 才有真实稻田淹水物理（4 个 crop 文件均可用；默认 allocation=1 不需要第 5 个）。灌溉关也能走代码路径（flood/paddy 均被 CH4 识别），但只能算 code-path pilot。
- 区域：全球 2° 网格本身的 patch/CFT 已含稻田分布（活跃阈值 `PADDY_RICE_FRAC_MIN=0.01`）；东亚/东南亚/南亚大格通常满足。区域裁剪（RegionClip）是可选 Phase-2，需要额外预处理。
- 估时：30 天 = 1,440 步；按 §5 的旧运行甲烷成本（≈227 core-s/步）→ 64 CPU 约 **1.4 h**；90 天约 4–5 h；加 mkinidata/构建，30 天 pilot 整体 2–6 h @64 CPU。**先跑 3–5 天标定**。

## 5. E：两年全球——已找到旧运行的权威日志（7920 本身无基准）

7920 上**没有** CoLM 全球甲烷运行日志/结果（E 段搜索只命中 OBRef methane smoke 与 BD 训练材料副本）。但**旧运行在 BD 的完整日志仍在**，本轮只读读取了 LSF 输出与日志关键行（未复制大数据）：

- job `168594`（`ch4_flat_2yr`，BD `myCluster`）：**成功完成**。720 tasks = 48×**15 nodes**；`rusage[mem=150G]`/task；提交 2026-05-27 20:04:29，结束 23:14:14。
- **mkinidata**（48 ranks）：20:04:31 → 20:08:59 = **4m28s**；mkinidata 日志 "Overall system time used: 4 minutes 23 seconds"；"CoLM Initialization Execution Completed"。
- **colm**（720 ranks）：20:08:59 → 23:14:04 = **3h05m05s**；日志末端 **"CoLM Execution Completed."**（最后 rank "Time elapsed: 3h04m49s"）；`accumnstep` 最大 **35,088**（与 W0 探测步数一致；1995–1996 两年、30 min 步长）。
- LSF 记账：**CPU time 7,952,842 s ≈ 2,209 core-hour**；Max Memory 1,288,799 MB、Average 741,290 MB（720 tasks 聚合）。
- 输出量：`history/` 899 MB（含 2 年 hist/hist_tracer/hist_unitcat），`restart/` 2.3 GB（1995/1996/1997-001 三个 checkpoint + const）。
- 旧运行自带分析（只读）：1996 年 CH4 总量 **119.813 Tg CH4/yr**（production 145.848、oxidation 6.871、ebullition 96.217、diffusion 28.646、aerenchyma 3.004；面积/时间平均 finundated 0.02924）；`analysis/ch4_1996_monthly_summary.csv` 有逐月分解。**这是候选对照数据，不是 ground truth**；且与当前 MATERIAL 谱系变量集不同。

### 7920 两年估时（用真实基准重算）

- 单位步成本（旧运行）：7,952,842 / 35,088 ≈ **227 core-s/步**（约 0.316 s/步 @720 ranks）。两年 = **≈2,210 core-hour**。
- 7920 单作业最多 64 CPU（10 槽、96 逻辑核、750 GB 上限）：**64 CPU 时 ≈ 34.5 h 计算**（理想）；96 ranks 时 ≈ 23 h。→ **单个 64-CPU 作业无法在 24 h 超时内完成两年**，必须**按年拆两个作业**（各 ≈17.3 h + 开销，每年 < 24 h ✅），用 YEARLY restart 续跑。
- 交叉校验：候选 history 的 `create_time` 差 1h28m45s/年；实际年均 ≈1h32m32s（colm 墙钟 3h05m05s/2）→ 与 1h28m45s 相差约 4%，口径一致，说明该推断口径成立。
- mkinidata：BD 4m28s@48 ranks（读 40k landdata 文件）→ 7920 估 5–20 min。
- 内存是 7920 上最大不确定项：BD 聚合 max 1.29 TB / 720 tasks ≈ 1.8 GB/task；若 7920 64 tasks 保持 ~1.8 GB/task → ~115 GB 可行；但若内存随每 rank patch 数增长（64 ranks 时 826 patches/rank vs BD 73），可能逼近/超过 750 GB 服务器上限 → **必须用 3–5 天短跑实测内存与步速**，必要时提高 ranks（>64）或先用 128 ranks 测。申请建议：短跑标定 64 CPU / 256 GB / 7200 s；年度作业 64 CPU / 256–512 GB / 72,000 s。
- 构建 `make all`：10–30 min（4 个源文件组合，具体待实操）。
- **完整两年 campaign 估：建 0.5 h + mkinidata 0.3 h + 年1 ≈17 h + 年2 ≈17 h + 分析 ≈ 35–40 h 墙钟（约 2,200 core-hour），跨 3 个 7920 作业**；若 7920 单核性能低于 BD，按 15–30 h/年区间准备。**没有 7920 甲烷基准，这是外推，不是实测。**

## 6. READY / NOT_READY 判定与逐项缺口

**NOT_READY（相对"精确 1995–1996 两年全球甲烷对照"）**，缺：

1. 15 个 CRUJRA 边界文件（变量/年/记录见 §1 表）——来源已实测存在，待授权传输（≈21.15 GiB 到新私有 data-stage；通道 `bd:` 现成）。
2. 冻结的甲烷参数文件（旧实际值 ≠ HEAD 默认：spatial pH off/on、lake prod off/on、rice off/on）。
3. 冻结的构建宏矩阵：`LULC_IGBP_PFT` vs `PC`、`CROP` 是否开（全球两年可不开，水稻必须开）、以及历史宏无法从源码恢复 → 需重建并验证。
4. 运行时依赖的身份（soilstate/cnsteadystate/routing/crop 与 BD 历史算例相同）——无哈希，未证实。
5. 7920 上的短跑标定（步速/内存/输出完整性）；当前所有耗时都是外推。
6. 拟合/对照定义：旧结果是 TRACER 谱系、spatial pH off、lake off、rice off；与 HEAD MATERIAL 谱系不可逐位比较 → "平场对照"需先定义口径。

**READY 的子项（证据充分，可立即执行）**：
- 候选传输一致性（40,450/26,618,457,915 B）✅
- 15 文件请求集与 HEAD 一致性 ✅（源文件逐字节相同 + W0 重放一致）
- 无传输窗口改写 W1/W2/W3 的请求性质量化 ✅
- 候选 landdata/参数/依赖可读性 ✅
- 水稻与两年全球的设计草案 ✅（PLAN_rice.md / PLAN_global2yr.md）
- 旧运行完成性/资源/估值基准 ✅（LSF+日志只读）

## 7. 需要用户决策（按优先级）

1. **输入闭合路线**：A 授权从 `bd:` 补 15 文件（21.15 GiB，精确原窗口）；B 只补 1997 的 8 文件（11.28 GiB，起点 03:00）；C 零传输窗口 W1（1995-01-01 03:00→1996-12-31 18:00，少 9 h）。
2. **参数冻结口径**：复现旧运行（spatial pH=off、lake prod=off、rice=off、hybrid、only_wetland=false）还是用 HEAD 新默认（pH=on、lake=on、rice=on）？两者结果不可互比，选择必须写死在冻结清单。
3. **构建宏矩阵与是否补 CI 用例**（CROP+BGC+MTON 目前无 CI 覆盖）。
4. **两年全球拆两个年度作业**是否接受（24 h 超时约束；每年 ~17 h @64 CPU）。
5. 用于 1995/1996 的 2005 年 landdata 是否接受（历史 case 即如此，属既有偏差）。
6. 运行时依赖身份（未证实）是否接受，或要求先从 BD 取哈希/重传（属新增数据操作）。

## 8. 证据文件（本目录 + `evidence/`）

- `evidence/audit.json`（审计全文，sha256 `1d6fb3cb…` 与 7920 一致）、`evidence/local-boundary-requests.csv`（23,407 请求）、`evidence/boundary_request_table.md`
- `evidence/probe_windows.json`、`evidence/probe_window.log`、`evidence/scan.json`（7920 只读调查原始结果）、`evidence/scan/results/metadata/job-result.json`（survey job `codex-fa5a1a057cc848e3a2468dd9797a73c0`）
- `evidence/g_global_ch4_flat_2yr_20260527_200428.nml`、`evidence/CRUJRA_local_abs.nml`、`evidence/run_g_global_ch4_flat_2yr_20260527_200428.lsf`
- `evidence/namelist_compat.json`
- `evidence/bd_historical_source.txt`（BD git HEAD/define.h/历史 `standard_ch4_parameter.nml` 只读快照）、`evidence/bd_define_h_sep2026_snapshot.h`、`evidence/bd_standard_ch4_parameter_historical.nml`、`evidence/bd_old_run_benchmark.json`
- `PLAN_rice.md`、`PLAN_global2yr.md`、`closure.json`、`provenance.md`
- 复现脚本：`scan/remote_scan.py`、`scan/probe_window.py`、`scan/sources/*`、`evidence/gen_closure.py`
