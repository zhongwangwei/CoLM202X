# 两年全球甲烷运行设计草案（PLAN_global2yr，未运行）

状态：**DESIGN_ONLY_NOT_RUN**；输入闭合 **NOT_READY**（精确 1995–1996 原始窗口），存在两个**不传新数据**的窗口改写方案（见 §2）。
源码基线：`CoLM-constituent` HEAD `10e9fbd204a99ff81bd60a98198b364d489e10d9`。
候选输入根：`/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`（传输一致 **TRANSFER_VERIFIED**，科学身份 **UNVERIFIED_CANDIDATE**）。

## 1. 输入闭合现状（精确到变量/年份/记录）

审计任务 `codex-ad6b48e0cb0044a9bbf0ce2903195234`（succeeded/0）用**真实生产时间索引例程**（`setstampLB/setstampUB`，与 HEAD `10e9fbd2` 逐字节相同）得到 23,407 条请求；其中 15 条落在候选包之外：

| 组 | 变量 | 文件年 | 记录（1-based） | 请求类型 | 模型时间 |
| --- | --- | --- | --- | --- | --- |
| 7 个 | `dlwrf,pre,pres,spfh,tmp,ugrd,vgrd`（除 dswrf） | 1994 | 1460（该年最后一条 6h 记录） | LB | 1995-01-01 00:00 第一步 |
| 8 个 | 全部 8 变量（含 dswrf） | 1997 | 1（该年第一条记录） | UB | 1996-12-31 最后一步（6h 变量 21:00，dswrf 18:00） |

- 1994 年 `dswrf` **不需要**：dswrf `offset=0` 且 forward/coszen，第一步 LB 直接用 1995 第一条。
- 完整表格与逐条证据：`evidence/boundary_request_table.md`、`local-boundary-requests.csv`（67 万字节请求台账）。
- 读取只会各取 1 条记录，但生产 reader 按绝对 record 下标寻址，**安全做法是补齐整年文件**：15 文件 × ~1.51 GB ≈ **22,706,115,044 B ≈ 21.15 GiB**（按 1995 同变量文件尺寸外推）。
- 若文件缺失，第一步/最后一步在打开对应年份文件时即失败（生产 reader 不做截断；`endyr` 只是警告，不能省文件）。
- **窗口改写实测**（`evidence/probe_windows.json`，同一生产例程，改动仅模型起止时刻）：

| 窗口 | 起 | 止 | 步数 | 需要候选包外文件 |
| --- | --- | --- | --- | --- |
| W0 原窗口 | 1995-01-01 00:00 | 1997-01-01 00:00 | 35,088 | 15（7×1994 + 8×1997） |
| W1 双端平移 | 1995-01-01 03:00 | 1996-12-31 18:00 | 35,070 | **0** ✅ |
| W2 只移起点 | 1995-01-01 03:00 | 1997-01-01 00:00 | 35,082 | 8（仅 1997 首记录，~12,109,928,072 B ≈ 11.28 GiB） |
| W3 只移终点 | 1995-01-01 00:00 | 1996-12-31 18:00 | 35,076 | 7（仅 1994 末记录，~10,596,186,972 B ≈ 9.87 GiB） |

结论：**精确原始窗口 NOT_READY；W1 无需任何新数据即可跑（少模拟 9 小时）；W2/W3 只需一半边界文件。** 这三者都必须由用户确认是否可接受，才可从"输入未闭合"转为"可进入短跑"。

## 2. 运行配置（候选 &nl_colm 改写，global 2°）

构建宏：`MaterialTransport + BGC + LULC_IGBP_PFT`（或 PC）+ `vanGenuchten_Mualem_SOIL_MODEL` + `GridRiverLakeFlow`；**不加 LULCC**（BGC+LULCC 上游不支持）。两年全球甲烷**不需要 CROP**（水稻模式关闭时）；若同一构建也跑水稻，见 `PLAN_rice.md`。

关键 namelist（HEAD 键名；候选旧 TRACER 键 15 个已失效，映射见 `evidence/namelist_compat.json`）：

```fortran
&nl_colm
   DEF_domain [-90,90,-180,180]; nx_blocks=18, ny_blocks=9; PIO_groupsize 按 7920 进程数重设
   DEF_GRIDBASED_lon_res=2, lat_res=2
   start 1995-01-01 03:00（W1/W2）或 00:00（W3/补齐数据后）; end 1996-12-31 18:00（W1/W3）或 24:00
   timestep=1800
   DEF_LC_YEAR=2005                       ! 候选 landdata 只有 2005 一年
   DEF_LAI_MONTHLY=.true., DEF_LAI_CHANGE_YEARLY=.false.
   DEF_USE_SoilInit + soilstate.nc; DEF_USE_CN_INIT + cnsteadystate.nc; DEF_USE_NITRIF=.true.
   DEF_USE_VariablySaturatedFlow=.true.
   DEF_USE_Dynamic_Wetland=.true.         ! hybrid 甲烷模式硬性要求
   DEF_UnitCatchment_file=.../grid_routing_data_15min.nc; DEF_Runoff_SCHEME=0
   DEF_forcing_namelist=<case>/config/CRUJRA_local_abs_7920.nml   ! 旧文件兼容，仅需改目录
   DEF_WRST_FREQ='YEARLY'; DEF_HIST_FREQ='MONTHLY'; DEF_HIST_groupby='YEAR'
/
&nl_colm_material
   DEF_USE_METHANE=.true.; DEF_METHANE_PARAM_FILE='.../standard_methane_parameter.nml'
/
```

甲烷参数（`run/standard_methane_parameter.nml`，hybrid 默认）：`inundation_mode='hybrid'`、`enable_rice_paddy=.false.`（两年全球对照）、`only_wetland=.false.`（湿地+土壤；若要"自然湿地仅"口径设 `.true.`，这是科学定义决策）、`use_spatial_ph=.true.`、`allowlakeprod=.true.`、`write_ch4_history='core'`。**旧 `standard_ch4_parameter.nml` 未随候选包传输**，历史参数值不可重现；当前默认值必须显式冻结。

`mkinidata` 需要用候选 landdata 生成初始 restart；实践中把候选 landdata 以只读 symlink 挂到 case 目录（`<dir_output>/<case>/landdata`），避免复制 2.36 GB。mkinidata 读 40,429 个 landdata 文件 + soilstate/cnsteadystate。

## 3. 资源与耗时估计（已有旧运行真实日志基准）

7920 上**没有** CoLM 全球甲烷运行日志/基准；但旧运行在 BD 上的 **LSF 记账和模型日志仍然存在**，本轮只读读取了关键行（未复制大文件）：

| 项 | 实测值 | 证据 |
| --- | --- | --- |
| 作业 | BD job 168594 `ch4_flat_2yr`，**成功完成**（"Successfully completed."） | `logs/lsf.168594.out` |
| 并行 | 720 tasks（15 节点 × 48），`rusage[mem=150G]`/task | 同上 |
| mkinidata | 48 ranks，20:04:31→20:08:59 = **4m28s**（日志 4m23s） | `mkinidata_global_current.log` |
| 模型 | 720 ranks，20:08:59→23:14:04 = **3h05m05s**；最后 rank "3h04m49s"；**"CoLM Execution Completed."**；`accumnstep=35,088` | `colm_global_current.log` 末尾 |
| CPU / 内存 | CPU time **7,952,842 s ≈ 2,209 core-hour**；Max Memory 1,288,799 MB，Average 741,290 MB（720 tasks 聚合） | LSF resource summary |
| 输出体积 | history 899 MB + restart 2.3 GB（3 个年度 checkpoint + const） | `du -sh` 只读 |

- 单位成本：7,952,842 / 35,088 ≈ **227 core-s/步**（≈0.316 s/步@720 ranks）；两年 ≈ **2,210 core-hour**。
- 交叉校验：候选 history `create_time` 差 1h28m45s/年，与实际年均 1h31m09s 一致（差 2.6%）。
- **7920 外推**：单作业 64 CPU 上限下，两年 ≈ **34.5 h 计算**（96 ranks ≈23 h）→ **单个 64-CPU 作业无法在 24 h 超时内完成两年，必须按年拆两个作业**（每年 ≈17.3 h + 开销 < 24 h ✅，用 YEARLY restart 续跑）。若 7920 单核性能低于 BD，按 **15–30 h/年** 区间准备。
- mkinidata：7920 估 5–20 min（BD 48 ranks 4m28s 读 40k landdata 文件）。
- 构建 `make all`：10–30 min（CROP+BGC+MT 组合未在 CI 覆盖，需真实编译）。
- **内存是最大不确定项**：BD 聚合 max 1.29 TB / 720 tasks ≈ 1.8 GB/task；若 64 tasks 保持量级 → ~115 GB 可行；若随每 rank patch 数增长（64 ranks 时 826 patches/rank vs BD 73）则可能逼近 750 GB 服务器上限。**必须先用 3–5 天短跑实测步速+内存**，必要时提高 rank 数（>64）。
- 完整两年 campaign 估：构建 0.5 h + mkinidata 0.3 h + 年度作业 ×2 ≈17 h + 分析 ≈ **35–40 h 墙钟、约 2,200 core-hour，跨 3 个 7920 作业**。这是**外推**，不是 7920 实测。

短跑标定建议：64 CPU / 256 GB / 7200 s；年度作业 64 CPU / 256–512 GB / 72,000 s；报告必须记录实测 step/s、峰值内存、输出完整性。

## 4. 科学/验收边界（不得夸大）

- 旧运行"完成"≠ 输入闭合 ≠ 科学正确。旧 history 是 **TRACER 代码谱系**输出（`hist_tracer`，64 变量），当前是 MATERIAL 谱系（`hist_material`）；两者不是逐位可比的"平场对照"，只是同参数包不同代码代的定性比较。
- 初始条件是 `soilstate.nc` + `cnsteadystate.nc` 气候态，**无模式 spinup**（候选 namelist `spinup_repeat=0`）；CH4/BGC 冷启动，头几年池子未平衡，不能用年总量做科学结论。
- 候选 landdata 是 **2005 年** land cover/LAI（DEF_LC_YEAR 默认 2005，landdata 只有 2005），用于 1995–1996 forcing 是历史 case 的既有偏差，不是本轮引入，但必须在报告中声明。
- 运行时依赖 `soilstate.nc`、`cnsteadystate.nc`、`grid_routing_data_15min.nc`、`nitrif/ndep/fire/snicar` 在 7920 可读，但**与 BD 历史算例内容是否相同未证实**（无 BD 侧哈希）；不能声称"同一输入复现"。
- 验收门禁：所有进程 END；NBAD=0；SIG/BAL/LULCC/residual hard=0；无未掩码 NaN/Inf；history/restart 文件集与日期集完整；甲烷闭合残差量级（此前空间 pH 两天测试最大 ~2e-22 mol/m²/s，仅为该测试证据，不外推全球两年）；**不把 finite 当 scientific pass**。

## 5. 需要用户决策

1. 输入闭合路线：**A** 补齐 15 文件（~21.15 GiB，需给出 7920 已有副本或授权从 BD/官方传输）；**B** 只补 8 个 1997 文件（~11.28 GiB，起点 03:00）；**C** 零传输、W1 窗口（1995-01-01 03:00 → 1996-12-31 18:00，少 9 h）。
2. 是否接受"两年全球"拆成两个年度作业（旧运行实测两年 ≈2,210 core-hour；7920 64 CPU 需 ≈34.5 h，超 24 h 单作业上限 → 必须拆年）。
3. 冻结的甲烷参数文件版本与 `only_wetland` 口径（湿地+土壤 vs 仅湿地）——影响与旧结果的可比性。
4. 是否要求同构建同时支持水稻（CROP 开）与两年全球；若共用构建需先补 CI 覆盖。
5. 2005 landdata 用于 1995–1996 是否接受（历史 case 已如此）。
