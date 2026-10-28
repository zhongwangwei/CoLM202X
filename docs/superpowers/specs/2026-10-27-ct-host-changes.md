# 组分输运移植：宿主修改登记表

> **历史记录。** 本文保留原阶段的设计、决定和验收记录；“待做”“未修”“进行中”以及旧命名不代表当前状态。最新进度以[物质输送 Draft PR 状态](../../material-transport-status.md)为准。旧比较器存在缺文件、缺变量和失败退出码处理缺陷，本文中的“逐位一致”记录须重新核验，不能直接作为本次合并依据。历史命令未作为本次可复现入口重新执行。

分支 `feature/constituent-transport`（基于 `feature/riverlake-flow` 915bdf45）。规则见计划第 4 节：
- 只有三种情况才改上游代码：(a) 组分要的宿主量没有别的途径拿到；(b) 不改就编译不过或运行出错；(c) 用户已拍板的替换。
- 每处都写理由和证据，没有证据的改动不进提交。

新文件（不算改上游代码，列在这里方便审阅）：
- `main/MOD_HistAccRestart.F90`（history 累加量 restart）；
- `.github/workflows/hist_acc_restart_fields.py`（字段清单检查）。

## P1

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `main/CoLM.F90` USE 段 | `USE MOD_HistAccRestart, only: read_history_acc_restart, save_history_acc_restart` | c | 侧车 restart 是用户定的移植内容（计划第 8 节第 3 条，默认启用） | — |
| `main/CoLM.F90` 河湖流初始化之后 | `CALL read_history_acc_restart (jdate)` | c | 河流累加量在 `grid_riverlake_flow_init` 里才分配，所以读要放在它之后；陆面和河流一起读，只占一个调用点（fork 用了 2 个） | 不动模型状态：只写 `a_*`、`nac` 和河流累加量；没有侧车时直接返回（P1 job 里 d-base/d-full 跟基线比，restart 逐位一致） |
| `main/CoLM.F90` `WRITE_TimeVariables` 之后 | `CALL save_history_acc_restart (jdate)` | c | 侧车要跟模型 restart 同时写；写侧车、标记完成在一个调用里做完（fork 用了 3 个调用，还要往模型 restart 里写标记） | 只写侧车文件，模型 restart 文件不变（P1 job 里 sha256 跟基线比） |
| `main/MOD_Hist.F90` `hist_out` USE 段 | `USE MOD_HistAccRestart, only: hist_acc_final_snapshot` | c | 同下 | — |
| `main/MOD_Hist.F90` `accumulate_fluxes` 之后 | `CALL hist_acc_final_snapshot (idate, deltim, itstamp, etstamp)` | c | 最后一步如果窗口没走完，`hist_out` 会先写出部分窗口、清零累加量，然后才轮到 CoLM.F90 存 restart。所以必须在清零前、在 `hist_out` 里面先存一份原始累加量。fork 改了 `hist_out` 的接口（加 3 个可选参数）并有 2 个调用点；这里不改接口，只有 1 个调用点 | 只在最后一步、窗口没走完时写文件，其余时候直接返回；不改任何模型或 history 变量 |
| `Makefile` `OBJS_MAIN` 和依赖行 | 加 `MOD_HistAccRestart.o`，以及两行依赖 | b | 新文件要编译；依赖行保证 `MOD_Hist.o` 在它之后编译 | — |
| `include/define.h` 末尾 | `#undef ConstituentTransport` + van Genuchten 检查 | c | 组分输运的编译开关（计划 3.2） | 宏未定义，不影响任何现有代码 |
| `.github/workflows/create_defineh.bash` | 可选的第 8 个参数 `CTON/CTOFF`（缺省 CTOFF）；`echo` 多打印 `$8`；生成的 define.h 末尾加宏 | c | CI 要能编译组分配置 | 原来 7 个参数的调用生成的 define.h 只多了 `#undef ConstituentTransport` |
| `.github/workflows/build_CoLM_gnu.yml` | awk 多取 `$9`；加一步 `hist_acc_restart_fields.py --check` | c | 把第 8 列传给脚本；上游以后新增 `a_*` 累加量时，提示同步侧车字段清单 | 原有配置行没有第 9 列，awk 输出不变 |
| `.github/workflows/TestCaseLists` | 加 2 个 CTON 配置（PFT、PFT+BGC） | c | 编译检查组分配置。名字把 `_CTON_` 放在 `_CROPOFF` 前面，避免 `grep $CaseName` 同时匹配到原有配置。BGC 配置用 PFT，因为 CI 脚本在非 PFT 时会关掉 BGC | 已用脚本核对，每个配置名在清单里只匹配一行 |

## P2（进行中）

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `Makefile` VPATH | 加 `main/CONSTITUENT` | b | 新目录要编译 | — |
| `Makefile` `OBJS_BASIC` / `OBJS_MAIN` | 加 21 个 `MOD_CT_*.o`，以及按 `USE` 生成的依赖行 | b | 组分文件要编译；各文件都包在 `#ifdef ConstituentTransport` 里，宏关时是空对象 | job `claude-383676353e1b44c89e0d5f3f21549114`：宏关（igbp）、宏开（igbpct、pftbgcct）三个配置都编译通过 |
| `Makefile` | `MOD_UserSpecifiedForcing.o` 从 `OBJS_MAIN` 挪到 `OBJS_BASIC`（`MOD_Qsadv.o` 之后） | b | `MOD_CT_Forcing` 要用 `metfilename`，而它会被 restart 和初始化钩子带进 `OBJS_BASIC`，也就链进 `mkinidata.x`；不挪就链接不上。它只依赖共享模块和 `MOD_Qsadv` | 同上 job 编译通过；宏关结果是否一致，留待 P2 验收 |
| `main/MOD_HistAccRestart.F90` `hacc_fields` | `#ifdef ConstituentTransport`：组分累加量 `a_trc_*`/`a_water_*` 的 44 个字段；只在组分 restart 读进来（`ct_land_restart_loaded`）时才读 | c | 取代 fork 的 tracer 侧车（`ncio_read_vector` 带了 fork 才有的 `known_present` 参数，没法在上游编译）；组分冷启动时，旧累加量属于别的状态 | 新文件；宏关不编译这段 |
| `main/HYDRO/MOD_Hydro_SoilWater.F90` `soil_water_vertical_movement` | `#ifdef`：7 个诊断输出（各层根系实取水、含水层/积水的根系份额、各来源的基流）；单次 aquifer exchange 前后快照，新私有子程序 `ct_exchange_layer_water`、`ct_split_aquifer_exchange` 按比例拆分 | a | 组分要知道每条路径的水量；宿主语句一行未改，诊断量不回写宿主 | 待 P2 验收（开组分 vs 不开逐位一致）。fork 的两次 exchange、Richards 投影、`qinfl` 扣 pond 都是改宿主物理，按方案 1 丢弃 |
| `main/MOD_SoilSnowHydrology.F90` `WATER_VSF`、`snowwater`、`snowwater_snicar` | `#ifdef`：诊断输出（根系、基流来源、不透水层蒸发分配、雪层出流、`qgtop`、透水标志、`wblc` 冰汇）；`WATER_2014` 不动 | a | 同上 | 同上。fork 的露水封顶、wdsrf 每步并入 gwat、霜冰溢出、湿地 smp/hk 重算、`1.e-3_r8` 容差都丢弃 |
| `main/MOD_LeafTemperature.F90`、`MOD_LeafTemperaturePC.F90` | `#ifdef`：可选输出冠层融雪/冻结质量（`qmelt*deltim`、`qfrz*deltim`）和 `raw` | a | 组分的冠层相变迁移 | fork 的 gssun 重赋值（缺 1e6）、`etr<0` 并入 `evplwet`、fwet 改写、`rstfacsun` intent 都丢弃 |
| `main/MOD_Thermal.F90`、`MOD_GroundTemperature.F90`、`MOD_PhaseChange.F90` | `#ifdef`：上述诊断量逐级传出；`meltf`/`meltf_snicar` 输出各层融冻质量 | a | 同上 | fork 的洪水蒸发重构、PFT 级 `ldew_p` 重分配丢弃 |
| `main/MOD_SnowLayersCombineDivide.F90` | `#ifdef`：四个合并/分裂子程序加可选的组分层数组，组分跟着水一起合并、分裂 | a | 雪层重排时组分要跟着走 | 与 fork 对比只差空行和换行（宿主语句未动） |
| `main/MOD_LeafInterception.F90` | 不改 | — | 截留的 xsc/gross 由 CoLMMAIN 用输出反推：`xsc = pg − (输入 − qintr)`、`gross = max(qintr,0)`；本版只有方案 1、8，fork 的冠层相变输出恒为 0 | — |
| `main/CoLMMAIN.F90` | `#ifdef`：组分调用（截留后 `ct_precip`、`ct_newsnow`、`ct_evapo`、`ct_soil_water`/`ct_wetland`、守恒检查、冰川/水体）；`DEF_VEG_SNOW` 关时冠层组分池由 `ct_canopy_lump/split` 只在组分侧重分配；`ct_soil_water` 收到的 `qinfl` 扣掉 aquifer exchange 拿走的积水；洪水反馈 + 组分暂时报错（P3 再接） | a | 组分记账入口；宿主变量不改 | fork 在 CoLMMAIN 改写宿主 `ldew_rain/snow`、`errorw`、`relocate_soil_frost_ice` 都丢弃 |
| `main/CoLMMAIN.F90` 调 `WATER_VSF` | 上游 bug 修复（不分宏）：`qflx_irrig_flood=qflx_irrig_drip` 改成 `=qflx_irrig_flood`，加注释 | b | 漫灌从未进入地面、滴灌算了两次，而 `CalIrrigationApplicationFluxes` 已从 `waterstorage` 扣掉漫灌量；`WATER_2014` 那一路本来就传对了 | 2026-10-27 用户确认修。只影响 CROP + 灌溉 + VSF 且有滴灌/漫灌的算例 |

| `main/CoLM.F90` BGC 初始化（ndep 年度/月度、`init_fire_data`） | 新局部变量 `bgc_data_year`，替换三处 `sdate(1)`；`#ifdef ConstituentTransport` 时取 `ct_bgc_data_init_year(sdate(1), s_year)`（`MOD_CT_HostFixes`） | b | `sdate` 经过 adj2end，1 月 1 日 00:00 起算时会变成上一年，整年都用上一年的 ndep；而在年内重启又会读到当年的，所以分段运行和连续运行不一致。做法与 fork 相同，fork 在 TRACER 下用 `s_year`。只在开组分时生效 | 组分关闭时 `bgc_data_year = sdate(1)`，与上游逐位一致；job `claude-43a513843b89434b8a5917f1e663b8fc` 验证 on0 = base、开组分拆分 = 连续。注意：开组分的 BGC 算例从此与 BGC 基线在 ndep 上有差，这是有意的 |
| `main/MOD_Hist.F90` 两处 `hist_write_time` 之前 | `#ifdef`：`CALL ct_hist_continue_file(file, idate, HistForm, file_last)` | b | 上游 `file_last` 启动时是 `'null'`，在组中间重启会重新创建本组 history 文件，重启前写的记录就丢了（拆分运行的 `case_hist_2005-01.nc` 只剩 2 个时刻）。fork 也没修这个。新子程序的判断是：开了组分、本次运行第一次写、起始时刻在本组开头之后，并且文件已存在；满足时把 `file_last` 设成这个文件，宿主就会接着追加（`ncio_write_time` 按时刻找记录，找不到就追加）。文件是否存在：网格 history 在 `DEF_HIST_mode='block'` 下由各 IO 进程检查自己的全部分块文件（`get_filename_block`），`MPI_LAND` 归约；其余情况（`one` 模式、Vector/Single）由 master 检查单个文件再广播。其他情况直接返回（bea42d4e 加了 block 模式） | 组分关闭时第一行就返回；job `claude-43a51384…` 验证 `one` 模式，job `claude-859bd66f…` 验证 `block` 模式：拆分运行（第 2 天拆开，后半段 5 或 9 进程）restart 与陆面/组分 history 都和连续运行逐位一致，只有 unitcat 分片多一个（分片按运行段命名，属设计如此）。尚未测：Vector 形式；CaMa history 只在缺失时创建，不需要处理 |
| `main/MOD_Hist.F90` 组分 history 段 | `#ifdef` 段内再包 `IF (ct_any_enabled())`：组分 history 的时间轴、窗口变量和 `ct_hist_out` 只在开了组分时执行；窗口变量 `history_window_seconds/_end_minutes` 只在 `HistForm=='Gridded'`、master、`DEF_HIST_mode='one'`、不开 `DEF_HIST_WriteBack` 时写进单个文件（bea42d4e） | b | 不开组分时原来也会建一个空的 `_hist_tracer_` 文件；`block` 模式下没有单个文件，master 写窗口变量时报 "No such file or directory"，每次运行都死在第一次写 history（与上游无关，是移植代码的 bug） | 同上 job；上游 `block` 模式对照（job `claude-9d98a82c…`）本来就能跑 |
| `Makefile` | 加 `MOD_CT_HostFixes.o`（在 `OBJS_MAIN`，`MOD_CT_Hist.o` 之后），依赖 `MOD_Hist.o: MOD_HistAccRestart.o MOD_CT_HostFixes.o`、`MOD_CT_HostFixes.o: MOD_CT_Namelist.o MOD_HistGridded.o` | b | 新文件；用到 `ghist`（`MOD_HistGridded`），所以从 `OBJS_BASIC` 挪到 `OBJS_MAIN` | — |
| `run/forcing/IsoGSM.nml` | 删掉 20 行 `DEF_forcing%precipitation_O18_*`/`precipitation_H2_*`/`water_vapor_O18_*`/`water_vapor_H2_*`，在 `/` 后面加注释 | b | 上游的 `DEF_forcing` 里没有这些成员，读这个文件会因 namelist 错误停下，IsoGSM 强迫根本用不了。fork 0f2020ae 删的就是这些行；同位素强迫改到同位素参数文件里，按物种写在 `&nl_colm_ct_forcing` 组中 | 同上 job：IsoGSM 基线（上游代码）和分馏组都直接读仓库里的这个文件 |

## P3（河湖流组分，溶解态：同位素、溶质）

新文件：`main/CONSTITUENT/MOD_CT_RiverLake.F90`（fork `MOD_Tracer_RiverLake` 改名；分汊改成全局路径顺序，见下）、`main/CONSTITUENT/MOD_CT_RiverLakeCoupling.F90`（宿主调用的胶水层，`ct_rlc_*`，加上洪水反馈的组分数组）。

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `main/HYDRO/MOD_Grid_RiverLakeFlow.F90` | `#ifdef`：17 处单行调用 `ct_rlc_*`（初始化/restart、径流输入、汇入缓冲、3 处 levee 重分配、子步输运、内流洼地溢出、history 累加、CoLMDEBUG 平衡、调用末尾清零；洪水反馈的分配、withhold、publish 三处、debit 两处；final） | a | 组分要跟着每一步水量走；逻辑都在胶水层，宿主只给参数 | job `claude-b743a51b…`：NEC levee+bif+dams 3 天，组分关、开组分无组分、开组分（同位素+CL）三种宿主 restart/history 都与上游逐位一致 |
| 同上 | fork 1330 速度截断后重算 `momen_riv` | — | 丢弃：改宿主物理 | — |
| 同上 | fork 的泥沙钩子（`route_*`、`particle_*`、降水推送）与甲烷的淹没面积发布（`publish_*_fldfrc_to_patches`） | — | 留给 P4（泥沙）、P5（甲烷） | — |
| `MOD_CT_RiverLake` 分汊 | 本地下游单元不再直接加，所有路径都经 `push_bif_influx` 按全局路径顺序求和（与水量 915bdf45 相同） | — | fork 的写法让结果随进程数变化 | 同上 job：开组分 5/9/17 进程全部变量（含 `trc_*`）逐位一致 |
| `main/HYDRO/MOD_Grid_RiverLakeHist.F90` | `#ifdef`：写 history 时调 `ct_riverlake_write_history`；`flush_acc_fluxes_riverlake` 里调 `ct_riverlake_flush_acc` | a | 河流组分 history（`f_trc_*`） | 同上 |
| `main/MOD_Vars_TimeVariables.F90` | `#ifdef`：`WRITE_GridRiverLakeTimeVars` 与 `commit_GridRiverLakeRestart` 之间调 `ct_riverlake_write_restart` | a | 河流组分 restart（`trc_*`）写进同一个河流 restart 文件，在提交之前 | 同上 job：第 2 天拆分重启全部变量逐位一致；CoLMDEBUG 河流组分平衡最大相对残差 3e-14 |
| `main/CoLMMAIN.F90` | `#ifdef`：P2 的"洪水反馈 + 组分报错"只留给 CaMa；网格河湖流下接上陆面侧：洪水蒸发的同位素交换（`ct_flood_evap_loss`）、入渗带进土壤的组分（传给 `ct_soil_water` 的 `flood_tracer_input`/`flood_infil_water`，记到 `flood_tracer_land_patch`）、守恒检查加 `flood_heterogeneous_in` | a | 洪水反馈的组分账：陆面拿走多少，下个路由调用就从河流组分池扣多少（`ct_rlc_flood_debit` 核账，对不上就停） | job `claude-dbc48696…`：AMZ 开洪水反馈宿主与上游逐位一致，账本闭合（残差 1e-13）；NEC 开洪水反馈原先停在陆面检查，查明是下面两个组分模块问题，修后 job `claude-b006eb89…`：NEC/AMZ 开洪水反馈、径流方案 0 不开洪水反馈，宿主组分开关逐位一致；NEC 洪水反馈 5 vs 9 进程、第 2 天拆分重启（后段 9 进程）全部变量逐位一致；账本闭合（残差 1e-14） |
| `Makefile` | `OBJS_BASIC` 加 `MOD_CT_RiverLake.o`，`OBJS_MAIN` 在 `MOD_Grid_RiverLakeFlow.o` 前加 `MOD_CT_RiverLakeCoupling.o`，加依赖行 | b | 新文件要编译 | job `claude-0ffba1b6…`：组分开（有/无 CoLMDEBUG）、组分关都编译通过 |

组分模块修复（不改宿主，`main/CONSTITUENT/MOD_CT_SoilWater.F90`，fork 同样写法）：
- 含水层外排上限（`aquifer_export_tracer`，用于根系取水、基流、qcharge 外排三处）：只有含水层实际存有的水按示踪剂质量设上限，超出部分是宿主造成的含水层亏欠，按同一比值记账。原写法在没有混合参考水量的 patch（如 patchtype 1）上，`wa=0` 加 6e-14 mm 根系回补就把 1.9e-3 mm TOPMODEL 基流的示踪剂夹成 6e-14（径流方案 0 才有含水层基流，方案 3 碰不到）。
- 洪水晚到池的积水溢流：有洪水入渗时，超出普通地表池的 rsur 只能是求解后 `wdsrf > pondmx` 的溢流，改从晚到池按其比值取并记入 rsur/rnof。原写法从普通池取，`max(.,0)` 截断凭空造出示踪剂（NEC patch：上涌 13.09 mm、普通池 0.812 mm、rsur 3.97 mm，误差 3.378R = 6.77e-3）；晚到池原本也没有给溢流记账。
- 证据：job `claude-2403abfc…`（只含第一处）与 `claude-b006eb89…`（两处）：径流方案 3 与 AMZ 开洪水反馈修前修后组分全部变量逐位一致；径流方案 0 5 vs 9 进程逐位一致；TRC_SIG/TRC_BAL 全部通过；河流组分平衡最大相对残差 3.3e-14。
- `run/standard_chloride_parameter.nml`：CL 溶解上限 1.0e-2 改为 0.218 kg Cl/kg 水（25 °C 岩盐饱和），用户确认。

## P4（河湖流组分，泥沙）

新文件：`main/CONSTITUENT/FAMILY_PARTICLE/MOD_CT_Sediment.F90`（CoLM202X 泥沙模块改名，默认参数改成 CaMa-Flood v4，见下）、`run/standard_sediment_parameter.nml`（全部参数的 CaMa 默认值，注释里有旧 `DEF_SED_*` 名字对照表）。同位素四个文件挪进 `main/CONSTITUENT/FAMILY_ISOTOPE/`（只是搬家）。

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `main/HYDRO/MOD_Grid_RiverLakeSediment.F90` | 删除 | c | 用户定：用 CoLM202X 泥沙替换上游泥沙（计划 §8） | job `claude-ad216af1…`：NEC/AMZ（levee+bif+dams+洪水反馈）、NEC 无 levee/bif，开泥沙、泥沙+同位素+CL 时宿主 restart/history 与组分关逐位一致；溶解态组分与 P3 提交逐位一致，开不开泥沙不变；泥沙 5/9/17 进程、第 2 天拆分重启（后段 9 进程）全部变量逐位一致；`standard_sediment_parameter.nml` 与默认值逐位一致；河流组分平衡最大相对残差 3.3e-14；AMZ 72 次路由调用里 61 次有泥沙输入（NEC 1 月降水 < 10 mm/天 门槛，无产沙） |
| `include/define.h` | 删 `GridRiverLakeSediment` 块 | c | 同上；泥沙改由 `ConstituentTransport` + `GridRiverLakeFlow` 编译 | 同上 |
| `share/MOD_Namelist.F90` | 删 16 个 `DEF_SED_*`（声明、namelist 列表、bcast 共 48 行），保留 `DEF_USE_SEDIMENT`；没编组分或河湖流时 `DEF_USE_SEDIMENT = .true.` 打印说明并 `CoLM_Stop` | c | 参数改由 `&nl_colm_sediment_parameter`（`DEF_SEDIMENT_PARAM_FILE`）读；没编进来时原样忽略开关会让用户以为跑了泥沙（用户要求“加警告并退出”） | 同上 job `off-sedstop`：组分关编译 + `DEF_USE_SEDIMENT=.true.` 在 mkinidata 停下并打印说明；`DEF_USE_SEDIMENT=.false.` 的算例不受影响 |
| `main/HYDRO/MOD_Grid_RiverLakeFlow.F90` | `ct_rlc_diag_accumulate` 多传 `hflux_fc`、`levee_floodarea`；删掉上游泥沙留下的 `sed_clock_*`/`sed_elapsed` 局部变量 | a | 泥沙要用洪泛区通量和堤防返回的淹没面积；局部变量随上游泥沙一起失效 | 同上 |
| `main/HYDRO/MOD_Grid_RiverLakeHist.F90` | `#ifdef`：`ct_riverlake_write_history` 后调 `ct_lifecycle_route_write_history`，`ct_riverlake_flush_acc` 后调 `ct_lifecycle_route_flush_history` | a | 泥沙 history（沿用上游 `sedcon/sedinp/sedlayer/sedout` 开关名） | 同上 |
| `main/MOD_Vars_TimeVariables.F90` | `#ifdef`：`ct_riverlake_write_restart` 后、提交前调 `ct_lifecycle_route_write_restart` | a | 泥沙 restart 写进同一个河流 restart 文件 | 同上 |
| `main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90` | 只改注释（延迟读取的列表加上组分） | — | 上游注释提到泥沙模块 | — |
| `MOD_CT_RiverLakeCoupling`（胶水层，不算宿主） | 初始化/restart、降水推送（patch→ucat，按有效面积）、levee 重分配、分汊累加、子步前后储量与保护区面积、`ct_lifecycle_route_calc`（每次路由调用末尾）、final | — | 宿主只给参数 | 同上 |
| `Makefile` | VPATH 加 `FAMILY_ISOTOPE`、`FAMILY_PARTICLE`；`MOD_CT_RiverLakeCoupling.o` 依赖加 `MOD_CT_Lifecycle.o` | b | 新目录和新依赖 | 同上 job 三个编译通过 |
| `run/Canopy_interception_0p1deg_era5land_Amazon_sed.nml` | 删 `DEF_SED_*`，同值写进 `&nl_colm_sediment_parameter`，`&nl_colm_constituent` 指向本文件 | c | 上游示例参数原样保留，换成新的读法 | 同上 job `sedamz-nec` 读到的参数 |

泥沙默认参数（用户 2026-10-27 定：“改成 cama 的”）：标量全部与 CaMa-Flood v4 `cmf_ctrl_sed_mod` 相同（`lyrdph 5e-5`、`totlyrnum 5`、`psedD 2.65`、`lambda 0.4`、`visKin 1e-6`、`vonKar 0.4`、`IGNORE_DPH 0.05`、`sedDT 3600`、`pyld 0.01/2/2`、`dsylunit 1e-6`）；粒径 CaMa 代码里没有默认值，取 `etc/sediment/run_sediment.sh` 的 `0.001,0.000025,0.0000001`（原上游默认粗 20~200000 倍）。fork 原来硬性要求 `lyrdph ≥ 最大粒径`（CaMa 没有这个检查，活动层只是体积），改成只打印一次警告。仍与 CaMa 不同、留待以后（2026-10-28 对照 CaMa-Flood_v4 源码核过）：没有 `revEgia=.false.` 分支（CaMa 默认 `.true.`，与我们一样）；地貌步长每次调用内取 min(剩余, `max_timestep_s`)，CaMa 是固定 `sedDT=3600` 且须为 DT 整数倍；浅水单元（河道水深 < 5 cm）CaMa 在 `calc_entrainment` 里整段跳过、连当步坡面输入也不加（产沙丢失），我们把输入直接放进床层并记为沉降（守恒）；fork 自加的堤防保护区泥沙、分汊转移和转移限幅 CaMa 没有。原来记的四条不是差异：1% 浓度上限两边都对每个水深 ≥ 5 cm 的单元无条件做（我们在 `calc_sediment_exchange` 交换后做一次、加坡面输入后再做一次，CaMa 在交换+输入后做一次，总量相同）；地貌步长默认同为 3600 s（汇流每 `DEF_GRIDBASED_ROUTING_MAX_DT`=3600 s 调一次）；产沙的 10 mm/d 阈值和幂律两边都作用在当时段的降水率上（CaMa 按 `IFRQ_INP` 的输入时段，我们按模型步），降水同为 1 小时时算法相同；CaMa 示例脚本 `IFRQ_INP=24` 喂日降水，那只是配置不同（公式按 mm/h 雨强定义，`pyld=0.01` 是在什么降水时间分辨率下率定的，Hatono & Yoshimura 2020 原文：式 1 产沙 = βθ^c r^d，θ 度、r mm/h，β=0.1、c=1、d=2，用 Kim et al. 2009 的 1° 日降水、步长约 15 分钟、对亚马孙 Obidos 悬沙人工试凑率定，文中没有 10 mm/d 阈值；CaMa v4 代码默认 0.01/2/2 加阈值，出处未知。系数是日降水下率定的，按小时降水用会偏多；2026-10-28 用户定：先按小时用，不改）；`crocdph` 在 CaMa v4 里只读入和打印、计算不用；`MAX_SED_CONC=0.01` 就是 CaMa 的 1% 上限（`cmf_calc_sedflw_mod` 第 307 行）。`sed_frc` 三类顺序已确认是砂、粉、黏（与粒径 1 mm、25 µm、0.1 µm 和 CaMa-sed `mkSedfrc.py` 一致）：第 1 类在撒哈拉/阿拉伯最高（0.54/0.52），第 2 类在黄土高原最高（0.44），第 3 类在亚马孙/刚果偏高（0.30/0.29），全球均值 0.43/0.31/0.23。

## 河湖流上游修复（2026-10-27 用户要求直接修在本分支）

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `main/HYDRO/MOD_Grid_RiverLakeFlow.F90` 分汊块 | `sync_global_routing_dt` 从限幅之后挪到算 `ordinary_scale` 之前；开动量 dt 约束时约束之后再同步一次 | b | 限幅用的是本 worker 的 dt，同步只会把 dt 变短，已按长 dt 压低的通量不会恢复，结果依赖分区 | job `claude-519daf20…`：修后 NEC/AMZ 5/9/17 进程 restart/history 逐位一致；本次修前 5 vs 9 也一致（没碰上）；修前后差：`volwater_ucat` 最大 1.5e7 m³（量级 2e10）、`f_discharge` 最大 146 m³/s（量级 1e6）；耗时不变（1728 s vs 1772 s） |
| 同上 history 累加 | 有堤防单元一律用堤防返回的淹没面积（0 也用），不再回退洪泛曲线 | b | 岸顶到堤顶之间 `fldfrc = levee_frc`，`levee_frc = 0` 合法，回退会把受保护的洪泛区记成淹没 | 同上 job；只改 history `f_floodarea`（最大差 6.4e8 m²，量级 1.6e9） |
| `postprocess/RiverHistConcatenate.F90` | 复制 `time` 属性；按全局记录合并 `history_window_seconds/_end_minutes`；重建 `(unitcat_local)` 静态场（`mask_complete_upstream`） | b | 原来只收 `(unitcat_local,time)` 变量，这些都丢了 | 同上 job：block + 合并（1、2 进程）与 one 模式比，变量清单、维度、数值全同，属性只差 `f_bifflw_lev` 的 `missing_value`（one 模式本来就不写，数据里也没有缺测值，未改） |
| 同上 | 所有进程一起 `read_namelist`，之后非 master 再退出 | b | 原来非 master 在 barrier、master 在 bcast，同一通信域集合通信不匹配 | 同上 job：2 进程合并 rc=0，结果与 1 进程相同 |

## 改名与重组（2026-10-28，提交 0f98622c）

Constituent Transport 改名为 Material Transport（用户定），同时按“跟水走 / 按族”重组目录。纯机械改名，结果不变。

- 目录：`main/CONSTITUENT` → `main/MATERIAL`；跟水走的 11 个文件进 `GENERIC_WATER/`；各族进 `FAMILY_ISOTOPE/`、`FAMILY_SOLUTE/`（新 `MOD_MT_Solute`）、`FAMILY_PARTICLE/`（泥沙）；根目录留 Defs、Hist、HostFixes、LandPhase、Lifecycle、Namelist、Providers、Rest、RiverLakeCoupling、Vars。
- 名字：`MOD_CT_*` → `MOD_MT_*`，`ct_*` → `mt_*`（含 `ct_rlc_*` → `mt_rlc_*`），`DEF_CT_*` → `DEF_MT_*`，宏 `ConstituentTransport` → `MaterialTransport`，namelist 组 `&nl_colm_constituent` → `&nl_colm_material`、`&nl_colm_ct_forcing` → `&nl_colm_mt_forcing`，CI 参数 `CTON/CTOFF` → `MTON/MTOFF`，history 文件 `*_hist_tracer_*` → `*_hist_material_*`，NetCDF `trc_*` → `mt_*`。上游自己的 "constituent phases"/"soil constituents" 不动。
- 宿主文件只是跟着改调用名（CoLM.F90、CoLMDRIVER、CoLMMAIN、MOD_Hist、MOD_Namelist、HYDRO Flow/Hist/TimeVars 等），上面各节登记的改动内容不变。
- 证据：job `claude-0b4c90ed…`：按名字对照，与 P4 提交 ea874545 全部 restart/history 变量逐位一致（组分关的宿主、NEC 洪水反馈泥沙+同位素+CL、径流方案 3、AMZ），拆分重启逐位一致，无 TRC_SIG/TRC_BAL。（第一次 job `claude-1ea9e5a7` 编译失败：重组脚本把注释里提到的 `MOD_MT_Solute` 当成依赖，造成循环；已改为只认 `USE` 语句。）

## 网格河湖流：出口面宽、水库上游（2026-10-28）

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `main/HYDRO/MOD_Grid_RiverLakeNetwork.F90` 出口面宽 | `outletwth` 由本单元与下游单元河宽的平均改为本单元河宽（提交 090c5df8） | b | 用户定（比较了 min 与本单元宽度后选后者）。窄支流汇入宽河/湖时平均宽度是本单元的几百倍，有效长度 `topo_rivare/outletwth` 远小于 `rivlen`，CFL 不成立，出流隔步翻号：1e-13 扰动一天内放大到 O(1)；分汊时不同进程数结果不同，原因就在这里。CaMa-Flood `cmf_calc_outflw_mod` 用上游单元自身的 `D2RIVWTH` | job `claude-7d69d572…`、`claude-805ad7b5…`（0.25° 全球河网）：1e-13 扰动下 ≥10 m³/s 且差 >1% 的单元 73–89 → 0；30 天子步数少 27%（dt 一直在 60 s 上限）；5/9 进程逐位一致。改变结果：两侧河宽不同的面占 69%（0.25°）/50%（0.1°）；冷启动 30 天 10 大河口平均流量低 0–13%，河道储量最多高 1.6%；min(w, w_next) 稳定性和耗时相同，但结果改变约为两倍 |
| 新文件 `main/HYDRO/MOD_Grid_ReservoirUpstream.F90` + `MOD_Grid_RiverLakeFlow.F90` 钩子（提交 620e521d） | 下游是已建水库的单元，出流改为本单元河道+洪泛区的运动波（CaMa `cmf_ctrl_damout_mod` UPDATE_INFLOW，I1DAM=10；PMINSLP 1e-5、PMANFLD 0.10、洪泛坡度上限 0.005；坡度用岸顶高程差/`topo_distance`）；这些单元速度取运动波速度，不受动量 dt 约束 | b | 用户报告（0.1° 下密西西比开水库：普通河段流量下降、储量升两个量级、水深几十米、淹没比例近 100%），用户定“照 CaMa 改”。水库全部 GRanD 库容放在坝所在一个单元，水位按该单元曲线算，远高于真实湖面（0.1° 正常库容下坝单元水深中位数 16 m、P90 64 m），经 HLL 面把上游单元顶住、倒灌。冷启动 `volresv = volresv_normal`，第一天就出现 | job `claude-69b7d73b…`（静态统计：0.1° 有 1666 座水库的上游单元岸顶低于正常水位，共 72150 个单元）、`claude-c75c092f…`：不开水库逐位一致；5/9 进程逐位一致；拆分重启 restart 逐位一致；开关组分宿主一致。第 6 天紧邻上游单元：水深 P90 7.1 → 0.8 m、最大 109 → 16 m、>10 m 的单元 260 → 1、储量 3.1e11 → 1.3e10 m³；水量平衡不变（约 0.04 m³）。未照搬：CaMa 的出流 ≤ 储量/dt 上限（这里由储量 dt 约束或分汊限幅保证非负）、CaMa 在坝上游单元也关分汊路径（这里仍只在水库单元关） |
| `main/MOD_Hist.F90` history 续写 | 不改：仍只在开物质输送时生效（P2 的 `mt_hist_continue_file`） | — | 2026-10-28 曾改成对所有运行生效（新模块 `MOD_HistContinue`），用户决定不修，已撤回；补丁留在 `.tmp/ct/histfix/always_on.patch`。不开物质输送时拆分运行的 history 仍与连续运行不同（上游原有问题），restart 逐位一致 | — |

综合验证（HEAD 620e521d，job `claude-9344c40e…`，NEC/AMZ levee+bif+dams+洪水反馈 3 天）：不开物质输送拆分重启 restart 逐位一致（history 不同，即上游组中重启覆盖问题，按用户决定不修）；开物质输送（泥沙+同位素+CL）宿主 restart/history 与不开逐位一致；5/9/17 进程、NEC/AMZ 第 2 天拆分重启（后段 9 进程）全部变量逐位一致；无 TRC_SIG/TRC_BAL；不编物质输送开泥沙在 mkinidata 停下并打印说明。两个提交已推送 origin。

## P5（甲烷，FAMILY_GAS）

湿地分解用甲烷自己的影子库（用户定“影子库，宿主不动”）：每个湿地 patch 每步先存下宿主 BGC 这一列，换进影子库，跑上游 `decomp_rate_constants_bgc`/`SoilBiogeochemPotential`/无植物竞争/`SoilBiogeochemDecomp`，取出甲烷要读的量后把宿主列原样恢复；库的更新在甲烷收尾时按 fork 的运算顺序做在影子库上。`decomp_hr`/`er`/`ar` 不改，甲烷另外输出只含 CO2 的呼吸 `f_co2_hr_methane`。fork 改宿主 BGC 的地方（湿地列分解、NStateUpdate1 重构、`decomp_hr` 扣甲烷碳）全部按方案 1 丢弃。

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `mksrfdata/MKSRFDATA.F90` | `#if MaterialTransport && BGC`：`Aggregation_SoilParameters` 之后按甲烷设置调 `Aggregation_LakeSoilC`（湖底碳）、`Aggregation_MethanePH`（土壤 pH） | c | 甲烷湖泊产甲烷、空间 pH 要的地表数据；只在甲烷设置要求时生成，别的地表数据不动 | 编译通过；运行时未跑（7920 的 rawdata 是旧目录结构，新 mksrfdata 跑不起来，按用户意见用旧 landdata，湖底碳文件用脚本按 OM 密度生成） |
| `mksrfdata/MOD_RegionClip.F90` | `#if`：裁剪 `soil/<年>/lake_soilc_patches.nc`、`methane_ph_patches.nc`（存在才裁） | c | 区域裁剪时带上上面两个文件 | 编译通过；未跑区域裁剪 |
| 新文件 `mksrfdata/Aggregation_LakeSoilC.F90`、`Aggregation_MethanePH.F90` | 从 fork 搬来，宏改名；湖底碳读 OM 密度用 `MOD_MT_NcIO` | c | 同上 | 同上 |
| `Makefile` | VPATH 加 `FAMILY_GAS`；新组 `OBJS_MT_CONFIG`（`MOD_MT_Namelist`、`Defs`、`NcIO`、甲烷 Registry/Const/Preprocessing/PHMapping），按共享模块的规则编译，链进三个程序（mksrfdata 也要读甲烷设置）；`MOD_MT_Namelist.o`、`MOD_MT_Defs.o` 从 `OBJS_BASIC` 移过去；`OBJS_MKSRFDATA` 加两个 Aggregation；`OBJS_MAIN` 加 13 个甲烷对象和 `MT_RegisterProviders.o`，依赖行按 `USE` 写 | b | 甲烷要用 `OBJS_MAIN` 里的 BGC 土壤模块，所以放在 `OBJS_MAIN`；`mt_register_providers` 从 `MOD_MT_Providers.F90` 拆到 `MT_RegisterProviders.F90`，由 `MOD_MT_LandPhase` 当参数传给 `tracer_lifecycle_init`，这样链进 mkinidata 的 `MOD_MT_Lifecycle` 不再引用它 | job `claude-697e7b6b…`：PFT+BGC 物质输送关/开、IGBP 开物质输送三个配置编译通过 |
| `MOD_MT_RiverLakeCoupling`（胶水层，不算宿主） | 最后一个子步的淹没面积/水深/堤防淹没面积存下来，路由调用末尾映射到 patch 发布给甲烷（fork 的 `publish_*_fldfrc_to_patches`，宿主 Flow 文件不改）；面积与河流 history 一致（有堤防单元用堤防结果，0 也用） | — | 甲烷淹没方案 7 | job `claude-b8b84148…`：NEC/AMZ 堤防+FB 路由下甲烷方案 7 的宿主与只开 CL 逐位一致；甲烷湿地排放非零 |

测试（旧 PFT landdata，lc 2005，1 个月，NEC 东北/AMZ 亚马孙两区）：

- job `claude-779f4200…`：物质输送编进去但不开组分 == 不编物质输送，逐位一致；甲烷 np 5/9/17 逐位一致；NEC/AMZ 第 2 天断点续算逐位一致；路由（方案 7，堤防+FB）np 5 vs 9 逐位一致。
- job `claude-b8b84148…`：甲烷+CL vs 只开 CL，宿主 restart/history 共同变量（含 `mt_CL`）在 NEC/AMZ × wetwat/路由堤防+FB 四组逐位一致；续算逐位一致；`f_co2_hr_methane` 非零；甲烷平衡残差 ≤1e-22。
- job `claude-1248e644…`、`claude-4b7b0404…`：PFT+BGC 开/关物质输送、IGBP/USGS 开物质输送都能编译。注意 define.h 只允许 BGC 配 PFT/PC（CI 生成脚本只允许 PFT），所以 IGBP/USGS 那几组里 BGC 实际是关的；PC+BGC 放到 P6 的 job 里单独编。非 PFT/PC 开 `DEF_USE_METHANE` 的运行期报错只是兜底，正常情况下没编 BGC 那条检查先报错。
- 宿主对比只能拿“甲烷+另一组分”对“只开另一组分”：P2 的 ndep 起始年修复（`mt_bgc_data_init_year`）在任何组分打开时生效，所以开甲烷 vs 不编物质输送从第 2 天起 ndep 就不同，这是 P2 登记过的有意改动，不是甲烷问题。
- 没测：mksrfdata 两个 Aggregation 的运行、空间 pH（7920 没有 `soil/PHH2O1.nc`）、水稻（要 CROP）、与 fork 两年全球甲烷结果对比。

## P6（LULCC）

fork 把物质输送的 LULCC 处理直接写在 `MOD_Lulcc_Driver` 里（`#ifdef TRACER`，约 200 行：旧/新 patch 面积、转移矩阵核对、各池重映射）。这里全部搬进新文件 `main/MATERIAL/MOD_MT_Lulcc.F90`，上游 driver 只加两个调用。fork 对 `MOD_Lulcc_Initialize` 的改动（旧/新进程 patch 有无一致的检查）在 `mt_lulcc_remap` 里已经覆盖，不搬；fork 的 `rewind` 参数和 `zwt_init` 改动属于 fork 宿主，不搬。

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
| `main/LULCC/MOD_Lulcc_Driver.F90` `SAVE_LulccTimeVariables` 之后 | `#ifdef MaterialTransport`：`CALL mt_lulcc_save ()` | c | 在 `LulccInitialize` 换掉 landpatch/mesh/pixel 之前存下物质池和旧 patch 面积；只读宿主变量 | job `claude-6247c29c…`+`claude-bc43eaf3…`：IGBP/PFT/PC × NEC/AMZ × SAT/MEC，物质输送编进不开 = 不编、只开 CL 宿主 = 不编，restart/history 逐位一致（跨 2005→2006 土地覆盖变化） |
| `main/LULCC/MOD_Lulcc_Driver.F90` 释放 LULCC 内存之前 | `#ifdef MaterialTransport`：`CALL mt_lulcc_remap (jdate(1), dir_landdata)` | c | 宿主变量恢复（SAT/MEC）之后、`landpatch_`/`patchclass_`/转移矩阵释放之前，把物质池映射到新 patch；只写物质状态 | 同上；CL 在 MEC 下的面积核对和陆地水物质守恒断言通过；5/9 进程、01-01（正好跨 LULCC）和 01-03 续算逐位一致 |
| `Makefile` | `OBJS_MAIN` 在 `MOD_Lulcc_Driver.o` 前加 `MOD_MT_Lulcc.o`，加依赖行 | b | 新模块要 USE LULCC 模块 | 同上；PC+BGC+物质输送+LULCC 也编过（job `claude-a6e83f6d…`，目标文件里有 methane 和 MT LULCC 符号） |

测试数据：7920 的 rawdata 是旧结构，2005/2006 的 landdata 和转移矩阵用合并前 fork 的 mksrfdata 做（旧脚本 create_defineh 第 8 个参数要 `TRACEROFF`）。NEC、AMZ 在 IGBP 下各有 8 个单元的 patch 类型组成变化。

两个原有问题（已报告，等用户定，没改）：
- 同位素 12 月冷启动在第 39 步 `TRC_SIG` 失败（NEC patch 38，1.7e-12 > 容差 1.0e-12，相对约 3e-8）；有无 LULCC、新旧 landdata 都一样，1 月起不出现，只开 CL 不出现（job `claude-224c4d82…`）。
- 同位素 + LULCC 必停：物质输送要求 vanGenuchten，上游据此强制打开 VSF，同位素就要求土壤/湿地 patch 带含水层参考水量；`MOD_MT_Vars` 的重映射在 SAT（新土壤类没有参考）和 MEC（参考流到特殊 patch）下都有意 `CoLM_stop`。所以 P6 的物质验证只用了 CL。

上游 `LULCC` 不支持 BGC（`MOD_Namelist` 里 `LULC_USGS || BGC` 直接报错），所以甲烷（要 BGC）在上游碰不到 LULCC，甲烷的 LULCC 回调目前走不到。

## 计划中、尚未实施

| 文件:位置 | 改动 | 类 | 理由 | 证据 |
| --- | --- | --- | --- | --- |
