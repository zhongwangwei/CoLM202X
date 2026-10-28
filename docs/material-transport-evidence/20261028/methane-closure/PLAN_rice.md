# 水稻甲烷最小真实运行设计草案（PLAN_rice，未运行）

状态：**DESIGN_ONLY_NOT_RUN**。本文件不含任何模型结果；所有"可用"均指输入已可读，不代表科学验收。
源码基线：`CoLM-constituent` HEAD `10e9fbd204a99ff81bd60a98198b364d489e10d9`（分支 `feature/constituent-transport`）。
候选输入根：`/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`。

## 1. 构建宏（必须同时满足）

| 宏 | 原因 / 证据 |
| --- | --- |
| `MaterialTransport` | methane provider 属 MATERIAL 框架；`main/MATERIAL/README.md` |
| `BGC` | `main/MATERIAL/MOD_MT_Namelist.F90: check_settings`：`DEF_USE_METHANE` 无 BGC 直接 `CoLM_Stop` |
| `LULC_IGBP_PFT` 或 `LULC_IGBP_PC` | 同上；甲烷只在这些 subgrid 上构建/测试 |
| `CROP` | `FAMILY_GAS/MOD_MT_Methane.F90`：`enable_rice_paddy=.true.` 而无 CROP 编译时 `CoLM_stop`；define.h 中 CROP 依赖 BGC |
| `vanGenuchten_Mualem_SOIL_MODEL` | `include/define.h`：MaterialTransport 与 Campbell 互斥 |
| `GridRiverLakeFlow` | 标准参数 `inundation_mode='hybrid'` 需要（`MOD_MT_Methane_Const.F90::configure_methane_inundation_mode`，`#ifndef GridRiverLakeFlow` 直接报错）；给出 `DEF_UnitCatchment_file` |
| 不要 `LULCC` | 上游 BGC+LULCC 不支持；甲烷依赖 BGC，因此甲烷+LULCC 不是可运行组合 |

**CI 缺口（已实测）**：`.github/workflows/TestCaseLists` 第 89–93 行 MTON 用例全部为 `CROPOFF`；`CROPON` 用例（1–32 行）全部 MT/OFF。**没有任何 CI 组合是 CROP+BGC+MT**，水稻甲烷构建没有任何 CI 覆盖，因此短跑前必须先在 7920 上做一次真实 `make all`（且不能把"其它 CROP 构建通过"当成此组合通过的证据）。

## 2. 运行配置（基于候选 &nl_colm，按 HEAD 改写）

候选包唯一 case namelist 是 TRACER 时代的 `g_global_ch4_flat_2yr_20260527_200428.nml`，其中有 15 个键 HEAD 已不存在（见 `evidence/namelist_compat.json`）；不能直接使用，须改写为：

```fortran
&nl_colm
   ! 域与时间：水稻 pilot 建议全球 2 度域、窗口完全落在 1995 文件内部
   DEF_domain%edges=-90. edgen=90. edgew=-180. edgee=180.
   DEF_nx_blocks = 18,  DEF_ny_blocks = 9,  DEF_PIO_groupsize = 12   ! 7920 上按进程数重设
   DEF_GRIDBASED_lon_res = 2, DEF_GRIDBASED_lat_res = 2
   DEF_simulation_time%start_year=1995 start_month=6 start_day=1 start_sec=0
   DEF_simulation_time%end_year=1995 end_month=6 end_day=30 end_sec=86400
   DEF_simulation_time%timestep=1800.
   DEF_dir_rawdata = '<7920 路径>/CoLMrawdata/'      ! mkinidata 若重生成 landdata 才用到
   DEF_dir_runtime = '/media/zhwei/data02/zhwei/CoLMruntime/'
   DEF_dir_output  = '$OBREF_JOB_RESULTS/cases/'
   DEF_LC_YEAR = 2005                                 ! 候选 landdata 只有 2005
   DEF_LAI_MONTHLY=.true.  DEF_LAI_CHANGE_YEARLY=.false.
   DEF_USE_SoilInit=.true.  DEF_file_SoilInit='/media/zhwei/data02/zhwei/CoLMruntime/soilstate.nc'
   DEF_USE_CN_INIT=.true.   DEF_file_cn_init='/media/zhwei/data02/zhwei/CoLMruntime/cnsteadystate.nc'
   DEF_USE_NITRIF=.true.
   DEF_USE_VariablySaturatedFlow=.true.
   DEF_USE_Dynamic_Wetland=.true.                     ! hybrid 模式硬性要求
   DEF_UnitCatchment_file='/media/zhwei/data02/zhwei/CoLMruntime/unitcatchment/grid_routing_data_15min.nc'
   DEF_Runoff_SCHEME=0
   DEF_USE_IRRIGATION=.true.                          ! 真实稻田水管理；见 §4
   DEF_forcing_namelist='<case>/config/CRUJRA_local_abs_7920.nml'
   DEF_WRST_FREQ='DAILY'   DEF_HIST_FREQ='DAILY'  DEF_HIST_groupby='YEAR'
   DEF_HIST_mode='one'  DEF_REST_CompressLevel=1  DEF_HIST_CompressLevel=1
/
&nl_colm_material
   DEF_USE_METHANE = .true.
   DEF_METHANE_PARAM_FILE = '<case>/config/standard_methane_parameter_rice.nml'
/
```

参数文件（复制 `run/standard_methane_parameter.nml` 后只改水稻相关项）：

```fortran
&nl_colm_methane_parameter
   DEF_METHANE%inundation_mode   = 'hybrid'   ! 需要 Dynamic_Wetland=.true. + GridRiverLakeFlow
   DEF_METHANE%enable_rice_paddy = .true.     ! 需要 CROP 编译 + pftclass 61/62 数据
   DEF_METHANE%only_wetland      = .true.     ! 建议 pilot: 湿地 + 稻田，排除一般 upland 土壤 CH4
   DEF_METHANE%use_spatial_ph    = .true.     ! 候选 landdata 有 methane_ph_patches
   DEF_METHANE%allowlakeprod     = .true.     ! 候选 landdata 有 lake_soilc_patches；湖泊 CH4 可保留
   DEF_METHANE%write_ch4_history = .true.
   DEF_METHANE%ch4_history_vars  = 'core'
/
```

候选 forcing namelist `CRUJRA_local_abs.nml` 与 HEAD 的 `&nl_colm_forcing` **完全兼容**（26 个键 + 全部 derived-type 字段逐项核对），只需把 `DEF_dir_forcing` 改到候选根：`.../forcing/`。

## 3. 数据可用性（实测，见 `evidence/scan.json` 与 GAP_ANALYSIS）

- landdata（候选根 `landdata/`，2005 年，全球 2°）：含 `pctpft/2005/pct_crops_*.nc`（144 块，crop CFT 分数）、`soil/2005/methane_ph_patches_*.nc`（144 块）、`soil/2005/lake_soilc_patches_*.nc`（144 块），并有 `landpatch/landpft/landelm/mesh/LAI/htop/lakedepth/soil/topography` 等 mkinidata 输入族。**说明历史 case 是带 CROP + methane 预处理构建的**（推断：构建宏未随包记录）。
- CROP 运行时文件（`CoLMruntime/crop/`）：`plantdt-colm-64cfts-rice2_fillcoast.nc`、`fertnitro_fillcoast.nc`、`fertilizer_2015soc.nc`、`surfdata_irrigation_method_96x144.nc` 四件可读；`surfdata_irrigation_allocation.nc` 缺失，但只在 `DEF_IRRIGATION_ALLOCATION == 3` 时需要，默认 1 → 不影响。
- 水稻 CFT 常数：`nrice=61`、`nirrig_rice=62`；`irrig_method_flood=3`、`paddy=4`（`main/MOD_Vars_Global.F90`）。当 `DEF_USE_IRRIGATION=.true.` 时 `MOD_Irrigation` 会把 ivt=62 的 flood 改判为 paddy（`main/MOD_Irrigation.F90:295-299`）；CH4 侧 `ch4_rice_pft_is_paddy` 同时接受 flood/paddy，因此灌溉开关不是识别稻田的必要条件，但**没有灌溉物理就没有稻田淹水位**，`finundated` 只能来自 hybrid 的 WTD/routing，属于"能跑通"而非"物理真实"。两种口径要分开报告。
- forcing：1995、1996 两个完整年可读。**水稻 pilot 只要窗口完全落在文件内部就与 15 个缺失边界文件无关**：例如 1995-06-01 00:00 起步，LB=1995-05-31 18:00、UB 均在 1995 文件内（同理由 `probe_windows.py` 的生产时间索引例程确认边界请求只发生在窗口端点）。因此**水稻 pilot 不需要任何新数据传输**。
- 初始场：`soilstate.nc`、`cnsteadystate.nc`、`grid_routing_data_15min.nc` 可读；与历史 BD 算例内容是否相同**未证实**（无 BD 侧哈希）。

## 4. 可用区域与"最小真实运行"

- 全球 2° 网格的 landdata 已含全部 land patch/PFT/CFT；稻田分布由 2005 年 CFT 分数给出。甲烷激活阈值 `PADDY_RICE_FRAC_MIN`（0.01）在 `MOD_MT_Methane_BgcLink`，即 patch 内水稻 CFT 面积分数 >1% 才进入水稻甲烷路径。东亚/东南亚/南亚 2° 大格通常满足；具体活跃 patch 数应由模型日志/`ch4_*` history 实测统计，本设计不预先声称。
- **推荐 Phase-1（最小真实运行）**：全球 2° 域，1995-06-01 00:00 → 1995-06-30 24:00（30 天，水稻生长季，完全内部窗口）；CROP+BGC+MT 构建；`enable_rice_paddy=.true.`、`only_wetland=.true.`、hybrid、spatial pH、灌溉开。目标：验证 (a) `PH-LOAD` 与稻田激活日志、(b) `hist_material` 中 `f_*ch4*` 有限且非零、(c) 质量/残差诊断为零、(d) 与关闭 `enable_rice_paddy` 的对照差异确实来自稻田 patch。**这不是科学验收**：初始 CH4 池是冷启动，30 天不足以做年总量。
- **Phase-2（可选）**：把域裁剪到东亚/东南亚（mksrfdata `RegionClip` 从同一 2005 landdata 生成区域 landdata + 从全球 CRUJRA 抽区域 forcing 文件），可显著降低单步成本；需要新的区域 landdata/forcing 预处理，属于新增工作。
- 区域三个月（如 NEC/AMZ）对照不能替代水稻全球/区域真实运行；本设计不把旧 NEC/AMZ 结果当基准。

## 5. 资源估计（用旧运行真实日志重算，仍需短跑标定）

旧全球 2° 甲烷运行在 BD 的 LSF 记账（job 168594，成功完成）：720 ranks、两年 35,088 步用 **CPU time 7,952,842 s ≈ 2,209 core-hour**，即 ≈**227 core-s/步**；mkinidata 48 ranks 用 4m28s（见 PLAN_global2yr §3 与 `evidence/bd_old_run_benchmark.json`）。

- 30 天 = 1,440 步 → 约 **91 core-hour**；7920 @64 CPU 理想 ≈ **1.4 h**，按 2–3 倍安全系数估 **2–5 h**。90 天 ≈ 273 core-hour → **4–8 h @64 CPU**。
- mkinidata 5–20 min + 构建 10–30 min（CROP+BGC+MT 未在 CI 覆盖）。
- 内存：旧运行 720 tasks 聚合 max 1.29 TB / avg 741 GB ≈ 每 task 1–1.8 GB；水稻 pilot 的 patch 数相同（全球 2°），64 tasks 估 100–200 GB，建议先申请 256 GB。
- **必须先跑 3–5 天标定**（240–480 步，64 CPU / 256 GB / 7200 s）测步速、峰值内存、log 完整性，再决定 30 天；不要直接上两年或长窗口。
- 注意：旧运行的实际甲烷参数是 spatial pH=off（常数 6.2）、lake prod=off、rice=off（主 namelist 与运行日志实测）；水稻 pilot 会打开 rice 且大概率打开 spatial pH，**成本与结果都不与旧运行相同**，不能拿旧 119.8 Tg/yr 当对照。

## 6. 运行前必须回答的问题（用户决策）

1. 冻结的构建宏矩阵（是否就是 §1 的 6 个宏；是否允许新增 CI 用例覆盖 CROP+BGC+MT）。
2. 是否接受"灌溉开/关"两种口径；若只做 code-path pilot 可关，但要写明不声称稻田水位物理。
3. 水稻 pilot 的参照/验收指标（田间通量？GRWP？还是仅模型内部守恒+敏感性对照）。
4. 是否允许 2005 年 landdata 用于 1995 forcing（历史 case 就是这样，但这是已知科学偏差，不是本轮新引入）。
5. 若要求区域裁剪，需要另行授权 mksrfdata 预处理与区域 forcing 生成。
