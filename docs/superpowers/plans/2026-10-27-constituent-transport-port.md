# 组分输运（Constituent Transport）移植计划

> **历史记录。** 本文保留原阶段的设计、决定和验收记录；“待做”“未修”“进行中”以及旧命名不代表当前状态。最新进度以[物质输送 Draft PR 状态](../../material-transport-status.md)为准。旧比较器存在缺文件、缺变量和失败退出码处理缺陷，本文中的“逐位一致”记录须重新核验，不能直接作为本次合并依据。历史命令未作为本次可复现入口重新执行。

日期：2026-10-27
状态：用户 2026-10-27 已审定（决定见第 8 节）
源：CoLM202X `master`（fcb2d1eb + 工作区）的 `main/TRACER/`
目标：CoLM-SYSU/CoLM，基于 `feature/riverlake-flow`（PR #512）新开分支 `feature/constituent-transport`

---

## 1. 已定决策（2026-10-27）

| 项 | 决定 |
| --- | --- |
| 对外名称 | Constituent Transport；宏 `ConstituentTransport`，目录 `main/CONSTITUENT/`，模块前缀 `MOD_CT_*`；内部变量 `trc`/`tracers` 可以保留 |
| 兼容性 | 不兼容 CoLM202X 旧结果：history/restart 变量名一起改 |
| `.inc` | 一个都不留（上游只 include `define.h`） |
| Makefile | 不再用 cpp 探测 TRACER；所有文件都编译，文件内容用 `#ifdef` 包住，跟上游 BGC/HYDRO 一样 |
| namelist | 不再出现 `tracer`：`DEF_USE_METHANE` / `DEF_USE_SEDIMENT`（沿用上游）/ `DEF_USE_ISOTOPE` + `DEF_ISOTOPE_SPECIES` / `DEF_SOLUTE_NAMES` + `DEF_SOLUTE_PARAM_FILES`；框架级设置用 `DEF_CT_*`；删掉 `DEF_TRACER_NUM/NAMES/TYPES` |
| 参数 | 参数文件可选，默认值写在代码里；同位素物理参数全部移到 `&nl_colm_isotope_parameter` |
| 溶质 | 不拆 family，反应性只看 `reaction_mode`（NONE / FIRST_ORDER / PROVIDER） |
| 泥沙 | 废弃上游 `MOD_Grid_RiverLakeSediment`（宏 `GridRiverLakeSediment`、主 namelist 的 `DEF_SED_*`），改用 CoLM202X 的 `MOD_Tracer_Particle_Sediment`；在本 PR 里删上游模块，PR 说明写清替代关系；开关保留 `DEF_USE_SEDIMENT`，参数移进泥沙参数组 |
| land-history 侧车 restart | 移植为新模块 `MOD_HistAccRestart.F90`，`USE MOD_Vars_1DAccFluxes`；清单手写，CI 脚本比对声明；上游文件只加两个调用 |
| 宿主改动 | 遵守 [改上游代码要写注释]：每处写清修改理由和“关掉组分时为何不影响结果”，并用前后对比验证逐位一致；别动别人的代码，规模大就新建文件或子程序 |

## 2. 现状盘点（CoLM202X master）

### 2.1 组分代码本体

`main/TRACER/` 共 42 个 F90，约 43,500 行：

| 组 | 文件 | 行数 |
| --- | --- | --- |
| 框架 | Defs, Vars, Lifecycle(+Registrations/Stubs), LandPhase, Conservation, Rest, Hist, SpecialPatches, SoilInit, EvapLimit | ≈ 8,700 |
| 陆面水过程 | Forcing, ForcingInput, Precip, Evapo, SoilWater, Snow | ≈ 6,600 |
| 同位素 | Frac, Isotope_Registry/Registrations/O18/HDO | ≈ 1,450 |
| 河湖 | RiverLake | 3,109 |
| 泥沙 | Particle_Sediment | 4,396 |
| 甲烷 | Reactive_Methane_* (17 个) + BgcShim | ≈ 19,000 |

另有 4 个 `.inc`：

| 文件 | 行数 | 作用 | 去向 |
| --- | --- | --- | --- |
| `include/land_history_restart.inc` | 7,116 | 522 个 `a_*` 的侧车 restart（清单 + 读写） | 新模块 `MOD_HistAccRestart.F90` |
| `include/tracer_land_history_restart.inc` | 366 | 组分 `a_trc_*`/`a_water_*` 清单 | 并入 `MOD_CT_Vars`（普通 contains 子程序） |
| `include/tracer_lifecycle_providers.inc` | 14 | 甲烷/泥沙 provider 注册宏 | `MOD_CT_Providers.F90` 普通子程序 |
| `include/tracer_isotope_species.inc` | 10 | O18/HDO 物理注册宏 | 同上，`MOD_CT_Providers.F90` |

### 2.2 宿主钩子

宿主侧 34 个文件约 430 处 `#if…TRACER`，13 个文件 `USE MOD_Tracer_*`。主要分布：

| 文件 | `#if` 数 | 性质 |
| --- | --- | --- |
| MOD_SnowLayersCombineDivide | 84 | 雪层合并/分裂时的组分层间搬运（内联，跟水量逐层同步） |
| HYDRO/MOD_Grid_RiverLakeFlow | 54 | 河湖组分输运；PR #512 里已经全部剥掉，需要按新结构重新挂 |
| MOD_SoilSnowHydrology | 41 | 下渗、径流、地下水的组分通量 |
| CoLMMAIN | 31 | 各过程调用的 optional 参数传递 |
| MOD_Thermal / LeafTemperature(PC) | 18 / 11 / 11 | 蒸散通量拆分 |
| CoLM.F90 | 18 | 初始化、restart、history 侧车 |
| MOD_LeafInterception / HYDRO/MOD_Hydro_SoilWater | 17 / 17 | 截留、VSF 土壤水 |
| MOD_Hist / CoLMDRIVER | 15 / 12 | history 输出、驱动 |
| LULCC (Driver/Initialize) | 13 | 土地利用变化时重映射组分库 |
| BGC (CNCStateUpdate1, NStateUpdate1, Competition, CNSummary) | 16 | 甲烷与 BGC 的耦合 |
| mksrfdata (MKSRFDATA, RegionClip, Aggregation_MethanePH, Aggregation_LakeSoilC) | 7 | 甲烷用的 pH 和湖泊土壤碳 |
| 其余（PhaseChange, NewSnow, GroundFluxes, GroundTemperature, Vars_*, Namelist, Initialize, extends/interception…） | 约 60 | 零散 |

内联块（以雪层为代表）没法改成一行调用：组分搬运跟水量的合并/分裂逐步交织，拆出去就得把宿主算法复制一份。这类保留内联，但必须写注释。

### 2.3 上游现状（feature/riverlake-flow）

- 没有甲烷，没有同位素，没有 `Aggregation_MethanePH` / `Aggregation_LakeSoilC`。
- 有上游泥沙 `MOD_Grid_RiverLakeSediment.F90`（1,974 行）。`define.h` 里 `GridRiverLakeSediment` 已被 `#error` 挡住（跟新河湖流不兼容）。`MOD_Namelist` 里有 `DEF_USE_SEDIMENT` 和 `DEF_SED_*`（约 50 处引用），`DEF_hist_vars` 有 `sedcon/sedinp/sedlayer/sedout`。
- CI（`create_defineh.bash`）只有 7 个参数，没有组分配置。

## 3. 目标形态

### 3.1 文件映射

规则：`main/TRACER/MOD_Tracer_<X>.F90` → `main/CONSTITUENT/MOD_CT_<X>.F90`，再按下表去掉多余的层级词：

| 旧 | 新 |
| --- | --- |
| `MOD_Tracer_Reactive_Methane_*` / `MOD_Tracer_Reactive_Methane` | `MOD_CT_Methane_*` / `MOD_CT_Methane` |
| `MOD_Tracer_Reactive_BgcShim` | `MOD_CT_Methane_BgcShim` |
| `MOD_Tracer_Particle_Sediment` | `MOD_CT_Sediment` |
| `MOD_Tracer_Isotope_{Registry,O18,HDO}` | `MOD_CT_Isotope_{Registry,O18,HDO}` |
| `MOD_Tracer_Isotope_Registrations` + `MOD_Tracer_Lifecycle_Registrations` + 两个 `.inc` | `MOD_CT_Providers`（`ct_register_isotopes()`、`ct_register_providers()`，`#ifdef BGC` / `#ifdef GridRiverLakeFlow` 写在子程序里） |
| `MOD_Tracer_Lifecycle_Registrations_Stubs` | 去掉（见 3.3：所有可执行文件都链接同一组对象之后就不需要桩了；如果 mkinidata 的链接顺序不允许，就在 `MOD_CT_Providers` 里加一个 `#ifdef` 分支，不另建文件） |
| `include/land_history_restart.inc` | `main/MOD_HistAccRestart.F90` |
| `README.md` | `main/CONSTITUENT/README.md`（改成新名字和新 namelist） |

模块内的 public 名字：对外接口（宿主会调的）改成 `ct_*`；只在组分内部用的 `tracer_*`、`trc_*` 不动，免得改动量失控。

### 3.2 宏

- `define.h`：`#undef ConstituentTransport`，加原来的两条检查（必须 vanGenuchten；甲烷要 BGC + PFT/PC）。注意：甲烷改成 namelist 开关之后，“BGC + PFT/PC” 这条检查要从编译期挪到 `DEF_USE_METHANE=.true.` 时的运行期检查。
- 删 `GridRiverLakeSediment` 宏及其 `#error`。
- 所有组分文件全文用 `#ifdef ConstituentTransport … #endif` 包住；甲烷文件再套 `#ifdef BGC`；泥沙、河湖组分再套 `#ifdef GridRiverLakeFlow`。

### 3.3 Makefile

- 删 `TRACER_ENABLED` / `METHANE_ENABLED` 的 cpp 探测和所有 `ifeq` 分支。
- `VPATH` 加 `main/CONSTITUENT`；组分对象按依赖顺序直接写进 `OBJS_BASIC` / `OBJS_MAIN`（跟 BGC 一样），在 7920 上用串行 `make all` 验证（`-j` 在 Linux 上本来就缺模块依赖）。
- 甲烷的 mksrfdata 依赖（Const、Registry、Preprocessing、PHMapping）放进 `OBJS_MKSRFDATA`。

### 3.4 namelist

主 namelist 文件里单独一个 `&nl_colm_constituent` 组，由 `MOD_CT_Namelist` 自己读（master 读、广播；`MOD_Namelist.F90` 不改）。整个组可以不写，不写就什么组分都不开。`DEF_USE_SEDIMENT` 留在上游原来的 `&nl_colm` 里（2026-10-27 用户同意）：

```fortran
&nl_colm_constituent
   DEF_USE_ISOTOPE           = .false.
   DEF_ISOTOPE_SPECIES       = 'H2_18O,HDO'     ! 只能取已实现的物种
   DEF_ISOTOPE_PARAM_FILE    = 'null'           ! 可选，&nl_colm_isotope_parameter + 各物种 &nl_colm_ct_forcing
   DEF_SOLUTE_NAMES          = ''               ! 例 'CL,NO3'
   DEF_SOLUTE_PARAM_FILES    = 'null'           ! 例 'CL:run/standard_chloride_parameter.nml'，可选
   DEF_CT_BALANCE_ABORT_NBAD = 0                ! 框架级（原 DEF_TRACER_*_ABORT_NBAD 等）
   DEF_CT_RESID_ABORT_NBAD   = 0
   DEF_CT_LULCC_ABORT_NBAD   = 0
   DEF_CT_AQUIFER_MIXING_WATER_MM = -1.         ! VSF 时同位素必须 > 0；Vars/Rest 对所有陆面水组分都用，所以放框架级
   DEF_CT_USE_SOIL_INIT      = .false.
   DEF_CT_SOIL_INIT_FILE     = 'null'
   DEF_CT_SOIL_INIT_VARS     = ''
/
```

`DEF_USE_METHANE`（P5，第二个 PR）以后也加进这个组。

参数组：

| 组 | 内容 | 原位置 |
| --- | --- | --- |
| `&nl_colm_isotope_parameter` | `DEF_ISO_*`：分馏开关、动力学方案、冰过饱和斜率、CG 湿度上限、开放水体/土壤动力学、升华表皮、冠层/融雪平衡、NSS 叶片参数；各物种的 `DEF_ISO_REF_RATIO_<物种>`、`DEF_ISO_INIT_DELTA_<物种>` | 主 namelist 的 15 个 `DEF_TRACER_*` + `standard_O18/HDO_parameter.nml` |
| `&nl_colm_solute_parameter` | `DEF_SOLUTE%unit_kind/mol_weight/init_conc/precip_default_conc/vapor_default_conc/max_dissolved_conc/decay_rate/charge` | `&nl_colm_tracer_parameter` |
| `&nl_colm_methane_parameter` | 不变（已经是这个名字） | `standard_ch4_parameter.nml`；`DEF_file_GIEMS`、`DEF_wetland_finundation_scheme` 一并挪进来 |
| `&nl_colm_sediment_parameter` | CoLM202X 现有的泥沙参数，加上游 `DEF_SED_*` 中仍有意义的 | 同名组 + 上游主 namelist |
| `&nl_colm_ct_forcing` | 每个物种自己的驱动输入（原 `&nl_colm_tracer_forcing`），多了 `forcing_species`：同位素文件是共用的，每组必须写；溶质自己的文件可以不写。一个物种只能有一组 | 物种参数文件 |

每个参数都有代码内默认值，文件缺省时直接用默认值跑（O18/HDO 的默认值取当前 `standard_*_parameter.nml`）。

`DEF_TRACER_TYPES` 的 `conservative`、`reactive` 两个旧写法随 `DEF_TRACER_TYPES` 一起删掉；descriptor 的 family 直接由开关决定（同位素开关 → ISOTOPE，溶质列表 → SOLUTE，甲烷 → GAS，泥沙 → PARTICLE）。

### 3.5 history / restart 命名

- history：`f_trc_<qty>_<species>` → `f_ct_<qty>_<species>`（机械替换，物种后缀不变）。甲烷自己的 `f_ch4_*` 和泥沙自己的字段保持现名。
- `DEF_hist_vars`：上游的 `sedcon/sedinp/sedlayer/sedout` 开关名保留不动，由新泥沙模块输出同名字段。只有含义确实对不上的字段才改，并在登记表里写明差别。
- restart：`trc_*` 变量 → `ct_*`；restart schema 版本号从 1 重新开始，读到旧 schema 直接报错并提示冷启动组分。
- 侧车 restart 文件名中的 `tracer` 字样同步改掉。

## 4. 宿主钩子规则

总原则（用户 2026-10-27 重申）：**除非必须，不改别人的代码；必须改的，要说明理由并给出证据。**

“必须”只包括三种情况：
- (a) 组分要读写的宿主量没有别的途径拿到（例如雪层合并时的层间搬运，只能在宿主循环里做）；
- (b) 不改就编译不过或运行出错；
- (c) 用户已经拍板的替换（上游泥沙模块）。

其余一律放进 `MOD_CT_*` 新文件或新子程序。宿主只加 optional 参数和一行调用。

每处宿主改动都登记进“宿主修改登记表”（`docs/superpowers/specs/…-ct-host-changes.md`，PR 说明附同一张表）：

| 文件:行 | 改动 | 属于 a/b/c | 理由 | 证据 |
| --- | --- | --- | --- | --- |

“证据”栏要写能复现的东西：
- 为什么非改不可：例如“该变量是 SUBROUTINE 局部量，循环结束即丢失，见 X.F90:123”；
- 为什么不影响结果：cpp 输出 diff 为空，或某次 job ID 的 restart/history 哈希逐位一致。

没有证据的改动不进提交。每个阶段结束时扫一遍全部被改的上游文件，把不在表里的 diff 退回去。

1. 一律 `#ifdef ConstituentTransport` 包住。宏没定义时，宿主文件预处理后的结果要跟 `feature/riverlake-flow` 完全相同（用 `cpp` 输出 diff 检查，比编译后对比更直接）。
2. 优先用“optional 参数 + 调一个 `MOD_CT_*` 子程序”。只有像雪层合并这种跟宿主算法逐步交织的地方才保留内联。
3. 每个钩子块上方写注释：做什么、为什么必须在这里、宏打开但所有组分开关为假时为何不改变宿主结果（例如“只读宿主变量”“只在 `present()` 时执行”）。
4. 宏打开、所有开关为假，结果也必须跟基线逐位一致。因为 Makefile 以后全编译，这个组合会是很多人的默认配置。
5. 不顺手改上游写法。需要宿主额外导出变量时，先找现成的 public 变量，找不到再加，并注明原因。
6. 上游泥沙相关的宿主调用（Flow、Hist、TimeVars、Namelist）随模块一起删，每处在 PR 说明里列出。

### 4.x fork 在 TRACER 块里改写的宿主物理（2026-10-27 用户定：方案 1）

证据：job `claude-4163d8a26f5d4b71a804507f7cac144c`，fork HEAD 编进 TRACER、不开组分，宿主 history 有 68（AMZ）/ 90（NEC）个变量和不编 TRACER 时不同；静态筛查约 70 个 `#ifdef TRACER … #else` 块改写宿主。两种写法的区域水量平衡误差都在 1e-2 m³ 以内，fork 的改写不是守恒修复。

规则：
- 宿主计算一行不改，`#else` 一侧（上游写法）原样保留；
- 组分要的分路径水量，在宿主算完后用前后差分 + 比例分配得到（仿 fork 不开 PHS 时的 `et_fraction`/`rsub_fraction`）；
- 宿主截断/清零处，组分按同比例一起截断，比值不变；
- fork 只为组分方便改写的宿主块（露水封顶、`wdsrf` 每步并入 `gwat`、含水层交换拆两次、不透水面分配、`relocate_soil_frost_ice`、`gssun` 重新赋值等）一律不搬。
- 验收加一条：开组分时宿主变量和不开组分时逐位一致。

## 5. 分阶段步骤

每个阶段一个或几个提交；每个阶段结束都要在 7920 上编译 + 回归（命令见第 6 节）通过后才进下一阶段。

### P0 准备（不改代码）
- 从 `feature/riverlake-flow`（915bdf45）拉 `feature/constituent-transport`，worktree `../CoLM-constituent`。
- 写钩子清单脚本：从 CoLM202X master 抽出宿主侧每个 `#if…TRACER` 块（文件、行号、内容），生成 `docs/superpowers/specs/…-ct-hook-inventory.md`。之后每个阶段在清单上逐条打勾，防止漏挂。
- 跑出基线：`feature/riverlake-flow` 在 AMZ 0.25° 上冷启动 2 天 + 拆分重启，保存 restart/history 的哈希。

### P1 骨架 + 侧车 restart
- 建 `main/CONSTITUENT/`；`define.h` 加宏；Makefile 只做加法，不改现有对象表，只追加组分对象和依赖行。CI `create_defineh.bash` 加第 8 个可选参数 `CTON/CTOFF`，缺省就是 OFF，原有 7 参数调用行为不变；工作流加 2 个配置（GRID_PFT + CT，GRID_PC + BGC + CT）。
- `MOD_HistAccRestart.F90`：把 `land_history_restart.inc` 改写成独立模块（522 个 `a_*`，手写清单），上游 `MOD_Hist`、`CoLM.F90` 各加一个调用。附 `.github/workflows/check_hist_acc_restart.py`：比对清单与 `MOD_Vars_1DAccFluxes` 的声明，有遗漏就让 CI 失败。
- fork 的侧车还改了 `MOD_NetCDFVector`（+295 行）、`MOD_Grid_RiverLakeHistState`（+335 行）和 `MOD_Hist`（`history_saved_raw`）。按用户 2026-10-27 的决定，这些功能都放进新文件，`MOD_NetCDFVector` 和 `HistState` 尽量不改。
- 侧车默认启用（第 8 节第 3 条）。
- 默认启用会多写文件，要验证两件事：
  - 写侧车不改变 land/river restart 本身的内容（哈希比对）；
  - 用上游旧 restart（没有侧车）启动时，按“新开 history 窗口”处理，不报错。
- 验收：宏关、宏开两种都跟基线逐位一致；侧车开启时拆分重启的 history 跟连续运行逐位一致。

### P2 框架 + 同位素 + 溶质（陆面）
- 移植 Defs、Vars（并入原 tracer 侧车清单）、Lifecycle、Providers、LandPhase、Frac、Isotope_*、EvapLimit、Forcing、ForcingInput、Precip、Evapo、SoilWater、Snow、SoilInit、Rest、Hist、Conservation、SpecialPatches。
- 新 namelist（`DEF_USE_ISOTOPE`、`DEF_SOLUTE_*`、`DEF_CT_*`、`&nl_colm_isotope_parameter`、`&nl_colm_solute_parameter`、`&nl_colm_ct_forcing`），删 `DEF_TRACER_*`。
- 挂陆面钩子：CoLMMAIN、CoLMDRIVER、CoLM.F90、SoilSnowHydrology、Hydro_SoilWater、SnowLayersCombineDivide、LeafInterception、Thermal、LeafTemperature(PC)、PhaseChange、NewSnow、GroundFluxes、GroundTemperature、Vars_*、Initialize、Hist、extends/interception。
- 注意：上游宿主跟 fork 宿主不一样（截留、火、上游写法），钩子按上游当前代码重挂，不照搬 fork 的上下文。
- **宿主 bug 修复：年度强迫数据滞后一年**（P0 基线发现，用户 2026-10-27 决定修）。
  - 问题：
    - 上游 `CoLM.F90:373` 用 `init_ndep_data_annually (sdate(1))` 初始化，1 月 1 日 00:00 冷启动时 `sdate` 已经 `adj2end` 成上一年；
    - `CoLM.F90:471` 跨年时用 `update_ndep_data_annually (idate(1))`，`idate` 还停在旧年份。
    - 两处叠加，运行时一直落后一年。火模块的 `init_fire_data (sdate(1))`（:383）和 `update_hdm_data (idate(1))`（:485）是同样写法。
  - 证据：基线 job `claude-75bf183e6d0e4caf928f3800e7c59fb5`。
    - 连续运行的 ndep 等于文件里 2004 年的值，从 01-02 拆分重启后等于 2005 年的值；
    - BGC 拆分重启因此不逐位一致（bgc restart 15 个变量、pft restart 28 个变量不同）。
  - 修法：年份一律取 `jdate(1)`。
    - 涉及 4 处：ndep 年度初始化、ndep 月度初始化的年份、ndep 年度跨年更新、火模块的初始化和 hdm 跨年更新。
    - **只在“编译了 `ConstituentTransport` 且至少开了一个组分”时生效**，否则走上游原来的调用，保证第 4 节第 4 条（宏开全关逐位一致）。
    - 用 `#ifdef ConstituentTransport` + 运行期判断 `ct_any_enabled()`，登记表归 (b) 类，证据栏引上面的 job。
  - fork `CoLM.F90:415–437` 的 TRACER 块只修了初始化、没修跨年，清单标为 `drop`，不照搬。
- **宿主 bug 修复：组中间重启会覆盖 history 文件**（P1 开头发现，用户 2026-10-27 决定修，同样只在开组分时生效）。
  - 问题：`MOD_Hist.F90:45` 的 `file_last='null'` 每次启动都会重置。重启后第一次写 history 时，`hist_*_write_time`（`MOD_HistGridded.F90:398` 等）判断 `filename /= filelast`，就重建本组文件，重启前写进去的时刻全部丢失。上游和 fork 都这样。
  - 证据：基线 job `claude-75bf183e6d0e4caf928f3800e7c59fb5`，拆分运行的 `case_hist_2005-01.nc` 只剩 01-02、01-03，连续运行有 3 天。
  - 修法（优先只动 `MOD_Hist`）：组分打开时，重启后第一次写之前，如果本组文件已存在，而且时间轴最后一个时刻早于当前时刻，就把 `file_last` 设成该文件名，接着追加。追加时 `itime` 怎么取，P2 先核对 `hist_*_write_time` 再定。如果必须改 `MOD_HistGridded/Vector/Single`，在登记表里写理由。
  - 验收：组分开时，拆分运行的 history 跟连续运行逐位一致；组分全关时跟基线一致。
  - **同时处理多余的部分窗口记录**（P1 验证发现，用户 2026-10-27 同意）。
    - 现象：运行结束时窗口没走完，`hist_out` 照常写一条部分窗口记录。有了侧车以后，重启后会把整个窗口续完，前一段那条部分记录就多余了。
    - 证据：job `claude-d87d778d6e53427ea9ff45a7065585b4`。拆分运行的 `case_hist_unitcat_2005-01.nc` 多一条 `history_window_seconds=86400` 的记录；最后一条跟连续运行逐位一致。河流 history 本来就追加所以已经能看到；陆面 history 修好追加之后也会出现。
    - 修法：侧车记下“前一段末尾写过部分记录”，重启后第一次写时覆盖这条记录，不再追加。陆面和河流都改。
- 验收标准更正：模型 restart 文件不能要求 sha256 相同，因为 NetCDF 的全局属性 `create_time` 每次写都不一样（P0 基线自己的连续和拆分运行之间也是这样）。改为“除 `create_time` 外内容相同”。
  - 验收：
    - 组分开 + BGC：拆分重启逐位一致，ndep 等于当年文件值；
    - 组分全关：跟基线逐位一致；
    - 火模块开：另跑一次跨年（12-31 → 01-02）检查 hdm 年份。
- 验收：开关全关逐位一致；O18/HDO + 1 个溶质 AMZ 2 天，陆面守恒残差 < 1e-10；拆分重启逐位一致；5/9/17 进程陆面逐位一致；跟 CoLM202X 同配置结果比较，差异只能来自宿主物理差别（列出量级）。

### P3 河湖组分
- 移植 `MOD_CT_RiverLake`，按 PR #512 的新结构（`MOD_Grid_RiverLakePush` 的 CSR、`MOD_Grid_RiverLakeVectorIO`、分汊全局通道顺序求和）重挂 Flow / Hist / TimeVars 钩子。
- 验收：开关全关时河流 restart 逐位一致；组分开启时河网守恒 < 1e-10；5/9/17 进程逐位一致（分汊组分通量也要用全局通道顺序求和，否则会重现 915bdf45 修过的进程数依赖）。

#### P3 细化（2026-10-27，开工前核对 fork `MOD_Grid_RiverLakeFlow` 的 54 个 TRACER 块）

范围：只做溶解组分（同位素、溶质）在河湖里的输运及其与陆面的交换。fork Flow 里另外两类块不在 P3：
- route 钩子（`route_forcing_put` / `route_diag_accumulate` / `route_calc` / `route_sediment_*`、`particle_*` 数组、降水上推到 ucat）：只有泥沙 provider 注册，归 P4；
- 洪泛面积 publisher（`publish_levee_fldfrc_to_patches`、`publish_fldfrc_to_patches`、`total_floodarea/flooddepth`）：只有甲烷 provider 注册，归 P5。

| fork 块 | 内容 | P3 处理 |
| --- | --- | --- |
| 172–212、133–138 的 restart 部分 | `river_lake_tracer_init`，读组分 restart，缺的按水量初值 | 挂到上游 `grid_riverlake_flow_init` 读完河流 restart 之后，一个 `ct_riverlake_init(...)` 调用 |
| 445–479 | 陆面 `trc_rnof_step` 上推到 ucat，`tracer_input_from_runoff` | 挂在上游径流上推之后；上推照抄宿主同一套 remap/push |
| 627–654 | 输入缓冲并入 `trc_inp_buf` | 同位置 |
| 703–708、1296–1300、2181–2186 | 堤防水量重分配时组分跟着走（`levee_tracer_repartition`） | 挂在上游三处 `levee_repartition_storage` 之后，只调组分子程序 |
| 1215–1227 | `tracer_substep`（主河道 + 分汊） | 挂在上游子步内水量更新之前同位置；**分汊组分通量改成和水量一样：下游一律经 `push_bif_influx` 按全局通道顺序求和**，不再本地直接加减 |
| 1272–1287 | 内陆洼地溢出时组分按比例移走 | 挂在上游同一分支，只动组分 |
| 1330–1333 | 速度截断后重算 `momen_riv` | **不搬**：改宿主物理（方案 1），登记为 drop |
| 1337–1343 | 组分诊断累加 | 同位置 |
| 565–573、1351–1356、1499–1563 | CoLMDEBUG 组分质量平衡 | 同位置，仍在 `CoLMDEBUG` 下 |
| 1577–1605（除 route_calc） | 限制器统计、清零累加量 | 同位置 |
| 73–79、219–259、1815–2186 的组分部分 | 洪泛反馈组分账（河→陆 credit、陆面蒸发/下渗回账、河侧扣减与总账检查） | 挂到上游 `init_flood_feedback` / `publish_flood_feedback` / 回账子程序里；陆面一侧 CoLMMAIN 1179–1247、1346–1356、1433 一起挂（现在开组分 + `DEF_GridRiverLake_FloodFeedback` 会因缺陆面组分停机）；LULCC 重建时组分数组跟着重建 |
| 1860–1862、2170–2176 | 有溶解上限的组分在格点内平衡 | 同位置 |
| 2207–2227 | final | 同位置 |
| Hist 577、680 | 组分河流 history、累加量清零 | 上游分片输出同位置 |
| TimeVars 683；`MOD_Vars_TimeVariables` 1341 | 写组分河流 restart | 放进 `WRITE_GridRiverLakeTimeVars` 之后、`commit_GridRiverLakeRestart` 之前（同 fork） |

验收（7920，NEC 0.25°，levee + bif + dams，3 天）：组分关 = 基线逐位一致；开组分时宿主河流/陆面 restart、history 与关组分逐位一致；河网组分守恒 < 1e-10；拆分重启逐位一致；5/9/17 进程组分逐位一致；开洪泛反馈时陆河组分总账闭合。

### P4 泥沙替换
- 移植 `MOD_CT_Sediment`，适配新河湖流结构；`DEF_USE_SEDIMENT` + `&nl_colm_sediment_parameter`。
- 删除上游 `MOD_Grid_RiverLakeSediment.F90`、`GridRiverLakeSediment` 宏、主 namelist 的 `DEF_SED_*` 及其广播、相关 history 开关的旧实现。这是 4.(c) 类改动，证据写进登记表：上游模块在新河湖流下已被 `define.h` 的 `#error` 挡住，现在本来就用不了；逐项列出新模块覆盖的功能和参数。
- 写“上游泥沙 → 新泥沙”对照表（参数、输出字段、物理方案差别），放进 README 和 PR 说明。
- 验收：AMZ 泥沙 1 天冒烟（对照 `Amazon_sediment_0p5_2020_smoke_1day`），泥沙质量守恒、拆分重启、多进程一致。

#### P4 细化（2026-10-27，开工前核对 fork `MOD_Tracer_Particle_Sediment` 与 fork `MOD_Grid_RiverLakeFlow` 的泥沙钩子）

模块：`MOD_Tracer_Particle_Sediment`（4,396 行）→ `MOD_CT_Sediment`，`#ifdef ConstituentTransport` 再套 `#ifdef GridRiverLakeFlow`。依赖只有河网（`MOD_Grid_RiverLakeNetwork` 30 处）、`MOD_WorkerPushData`、levee、分汊、`MOD_Grid_RiverLakeHistRoute`、`MOD_Tracer_Defs/Lifecycle`（改成 `MOD_CT_Defs/Lifecycle`，P2 已移植，`route_*` 回调接口都在）。泥沙是 provider 自管状态的 particle 组分，河湖溶解态代码已经用 `ct_uses_land_water_transport` 跳过它。

宿主钩子（fork Flow 里 13 处，全部收进胶水层 `MOD_CT_RiverLakeCoupling`，Flow 只加单行调用）：
- 初始化 + 读 restart（`route_init` / `route_read_restart`），final；
- 径流输入后把降水推到单元集水区（`prcp_*` 数组的 patch→inpm→ucat 推送，约 45 行，搬进胶水层）→ `route_forcing_put`；
- 3 处 levee 重分配（与溶解态同位置）→ `route_sediment_levee`；
- 分汊块里 → `route_sediment_bif`；
- 子步开始记水量（`particle_water_storage_start`/`particle_protected_start`），子步末尾组 `particle_*` 数组 → `route_diag_accumulate`；
- 一次路由调用结束 → `route_calc`；
- Hist：写 history（`route_history`）、`flush_acc_fluxes_riverlake` 里 `route_flush_history`；TimeVars：写 restart。
- `particle_floodarea` 用 history 累加里算出的淹没面积（fork 的 `total_floodarea`，甲烷 P5 也要），取上游修过 levee 回退之后的值。

输入：`sed_frc(nseqmax, sed_n)`、`sed_slope(nseqmax, slope_layers)` 从河网文件读；7920 的 `grid_routing_data_15min.nc` 和 `_sed.nc` 都有。

删除上游泥沙（4.(c)，用户已定）：`main/HYDRO/MOD_Grid_RiverLakeSediment.F90`（1,974 行）；`define.h` 的 `GridRiverLakeSediment` 与 `#error`；`MOD_Namelist` 的 17 个 `DEF_SED_*`（声明、namelist 列表、bcast）；`run/Canopy_interception_0p1deg_era5land_Amazon_sed.nml` 里的 `DEF_SED_*`。`DEF_USE_SEDIMENT` 留在 `&nl_colm`；开了但没编进 `ConstituentTransport`+`GridRiverLakeFlow` 就报错停。`DEF_hist_vars` 的 7 个开关（`sedcon/sedout/sedinp/sedlayer/bedout/netflw/shearvel`）上游已有，名字和 fork 输出一致，不改。

参数对照（上游 `DEF_SED_*` → `&nl_colm_sediment_parameter` 的 `DEF_SEDIMENT%*`）：

| 上游 | 上游默认 | 新名 | fork 标准文件 | 备注 |
| --- | --- | --- | --- | --- |
| `DEF_SED_DIAMETER` | "0.0002,0.002,0.02" m | `grain_diameter` + `nsed` | 2e-5, 2e-4, 2e-3 m | 上游比 fork 粗 10 倍 |
| `DEF_SED_DENSITY` | 2.65 g/cm³ | `grain_density` | 2650 kg/m³ | 单位不同 |
| `DEF_SED_WATER_DENSITY` | 1.0 g/cm³ | `water_density` | 1000 kg/m³ | 单位不同 |
| `DEF_SED_LAMBDA` | 0.4 | `porosity` | 0.4 | |
| `DEF_SED_LYRDPH` | 5e-5 m | `active_layer_depth` | 5e-2 m | 差 1000 倍 |
| `DEF_SED_TOTLYRNUM` | 5 | `ndeposit_layers` | 5 | |
| `DEF_SED_IGNORE_DPH` | 0.05 m | `ignore_depth_m` | 0.01 m | |
| `DEF_SED_VISKIN` | 1e-6 | `viscosity` | 1e-6 | |
| `DEF_SED_VONKAR` | 0.4 | `von_karman` | 0.4 | |
| `DEF_SED_PSET` | 1.0 | `settling_multiplier` | 1.0 | |
| `DEF_SED_PYLD/PYLDC/PYLDPC` | 0.01/2/2 | `yield_coefficient/slope_exponent/precipitation_exponent` | 0.01/2/2 | |
| `DEF_SED_DSYLUNIT` | 1e-6 | `unit_conversion` | 1e-6 | |
| `DEF_SED_CFL_ADV` | 0.5 | `cfl_adv` | 0.5 | |
| `DEF_SED_DT_MAX` | 3600 s | `max_timestep_s` | 3600 s | |
| — | — | `bed_depth` | （fork 新增） | |

fork 的参数没有代码默认值（-1，必须写在参数文件里）。

用户定（2026-10-27）：代码默认值用上游 `DEF_SED_*` 的默认值（换成新单位：密度 2650/1000 kg/m³），参数文件可选；上游示例 `run/Canopy_interception_0p1deg_era5land_Amazon_sed.nml` 改写成 `&nl_colm_sediment_parameter`，保留原取值并换算单位。

验收：AMZ 泥沙 3 天（对照 fork `Amazon_sediment_0p5_2020_smoke_1day`）+ NEC levee+bif+dams 泥沙 3 天：宿主开泥沙与组分全关逐位一致（泥沙不反馈河床）；泥沙质量守恒；拆分重启、5/9/17 进程全部变量逐位一致；泥沙 + 同位素 + CL 同开时溶解态结果与只开溶解态逐位一致。

### P5 甲烷
- 移植 Methane_*（17 个）+ BgcShim；BGC 钩子（CNCStateUpdate1、NStateUpdate1、Competition、CNSummary）；mksrfdata 的 MethanePH、LakeSoilC、MKSRFDATA、RegionClip 钩子。
- `DEF_USE_METHANE` 运行期检查 BGC + PFT/PC；`DEF_file_GIEMS`、`DEF_wetland_finundation_scheme` 挪进甲烷参数组。
- 验收：
  - BGC 开、甲烷关时跟基线逐位一致（这是 BGC 用户的默认配置，最重要）；
  - 甲烷开启时拆分重启逐位一致（依赖 P2 的年度强迫修复）；
  - 甲烷开启 2 年全球平场（对照 `g_global_ch4_flat_2yr`）的守恒检查与 fork 量级一致。

#### P5 细化设计（2026-10-28）

盘点结论（fork `main/TRACER/MOD_Tracer_Reactive_Methane*`，17 个文件 + BgcShim，共 18904 行）：
- 甲烷对宿主只有读，唯一的写是 `Methane_BgcLink` 的 `tracer_ch4_bgc_finalize_step`（改 `decomp_hr`、`er`，湿地还把 `ar` 置 0）和 BgcShim（湿地 patch 上用宿主 BGC 数组跑整套分解：清零并写约 30 个 BGC 通量数组，`CDecompStateUpdate(...,.true.)` 直接改 `decomp_cpools_vr`/`decomp_npools_vr`）。
- 宿主钩子已在 P2 挂好：`CoLMDRIVER` 的 `mt_lake_step`、`mt_wetland_decomp`、`mt_soil_step`、`mt_report`；restart、history、LULCC 都走 `MOD_MT_Lifecycle` 的 provider 回调（`tracer_lifecycle_hooks_type` 的字段与 fork 相同）。
- fork 改了 4 个宿主 BGC 文件（抽出 `CDecompStateUpdate`、`SoilBiogeochemNDecompStateUpdate`、`SoilBiogeochemCompetitionNoPlant`、`CNDriverSummarizeNonvegetatedSoilStates`），只为甲烷所用。

用户定（2026-10-28）：影子库，宿主不动。
- 湿地分解：甲烷自己保存湿地 patch 的 C/N 库（`decomp_cpools_vr`、`decomp_npools_vr`、`sminn_vr` 等，进甲烷 restart）。每步在湿地 patch i 上：先存下宿主数组第 i 列，换入影子库，调用上游原有的 `decomp_rate_constants_bgc`、`SoilBiogeochemPotential`、`SoilBiogeochemDecomp`，以及搬到甲烷侧的 NoPlant 竞争和库更新，再把结果存回影子库，最后把宿主第 i 列原样恢复。宿主 BGC 文件一行不改。
- `decomp_hr`、`er`、`ar` 不改；甲烷另出只含 CO2 的呼吸（history `f_mt_ch4_co2_hr` 一类）。
- 验收：BGC 开、甲烷关 = 基线；甲烷开时宿主全部 restart/history 与甲烷关逐位一致；甲烷开拆分重启、5/9/17 进程逐位一致；碳守恒（甲烷侧）；与 fork 2 年全球平场的 CH4 通量量级对照（湿地分解从“宿主库”换成“影子库”，初值相同，预期数值一致或很接近，差异要能解释）。

步骤：
1. 文件搬家改名：`main/MATERIAL/FAMILY_GAS/MOD_MT_Methane_*.F90`（BgcShim 改成 `MOD_MT_Methane_WetlandBGC`，实现上面的换入换出）；宏 `(TRACER && BGC)` → `(MaterialTransport && BGC)`；`MOD_Tracer_Defs/Hist/Lifecycle` → `MOD_MT_*`（函数名不变）。
2. 开关：`&nl_colm_material` 加 `DEF_USE_METHANE`、`DEF_METHANE_PARAM_FILE`；`DEF_file_GIEMS`、`DEF_wetland_finundation_scheme`（fork 在主 namelist）放进甲烷参数组；`MOD_MT_Defs` 加 CH4 条目（FAMILY_GAS / STATE_OWNER_PROVIDER / REACTION_PROVIDER），没编 BGC 或不是 PFT/PC 时报错退出；`mt_register_providers` 加甲烷。
3. 洪水：网格河湖流把淹没面积/深度发布给甲烷（fork 的 `publish_levee_fldfrc_to_patches`、`publish_fldfrc_to_patches`），走 `MOD_MT_RiverLakeCoupling` 胶水层，宿主只加一处调用。
4. mksrfdata：`Aggregation_MethanePH`、`Aggregation_LakeSoilC` 新文件；`MKSRFDATA`、`MOD_RegionClip` 加调用（只在开甲烷时生成）。
5. 示例 namelist 与参数文件；Makefile；CI 加一个 BGC+MTON 配置的甲烷编译。

#### P5 完成情况（2026-10-28）
- 1–4 步完成；第 5 步 CI 已有 `GRID_PFT_..._BGCON_MTON` 配置，现在会编进甲烷。`DEF_file_GIEMS` 等仍在甲烷参数组里（沿用 fork 的 `nl_colm_methane_parameter`）。
- 验收结果：见 `docs/superpowers/specs/2026-10-27-ct-host-changes.md` 的 P5 节（job `claude-779f4200…`、`claude-b8b84148…`）。宿主对比用“甲烷+CL 对只开 CL”，因为 P2 的 ndep 起始年修复在任何组分打开时都会生效。
- 还没做：与 fork 两年全球平场对照（数据在 HPC）；mksrfdata 两个 Aggregation 的运行测试（7920 rawdata 是旧结构）；空间 pH、水稻没跑；湿地影子库的 LULCC 重映射只按同 eindex+同 patchclass 拷贝（面积加权留给 P6）；`run/standard_methane_parameter.nml` 默认值沿用 fork（hybrid + 水稻，要 CROP 和 `DEF_USE_Dynamic_Wetland`）。

### P6 LULCC
- Lulcc_Driver / Lulcc_Initialize 钩子：组分库随土地利用变化重映射（含 2026-05-31 审计修过的 B2 问题）。
- 验收：LULCC 开 + 组分开，跨年守恒；LULCC 开 + 组分关逐位一致。

#### P6 完成情况（2026-10-28，提交 44edfec1）
- `MOD_MT_Lulcc.F90` + `LulccDriver` 两个调用；验收见宿主改动登记表 P6 节。
- 上游 LULCC 拒绝 BGC，所以甲烷不涉及 LULCC。
- 同位素 + LULCC 在 fork 代码里有意停下（含水层参考），待用户定；12 月冷启动同位素 TRC_SIG 失败，待用户定。

### P7 收尾
- README、示例 namelist（`run/` 下 iso / solute / ch4 / sediment 各一个，路径改相对路径）、`DEF_hist_vars` 默认清单补齐新字段。
- 全仓 `grep -i tracer` 只允许出现在内部变量名和注释的历史说明里。
- 最后一轮：CI 全配置 + 7920 完整回归矩阵。

#### P7 进展（2026-10-28，提交 9dccacb5）
- 字符串里的 tracer → material（报错、单位、长名），数值逐位一致；`main/MATERIAL/README.md`。
- `DEF_hist_vars` 不用补：物质字段跟宿主对应通量/状态的开关走。
- 未做：CI 生成脚本只允许 PFT+BGC（上游 define.h 允许 PC），CI 没有 PC+BGC+MT 和 LULCC+MT 配置。
- 示例 `run/examples/Amazon_Grid_25km_IGBP_VG_RiverLake_Material.nml`（同位素+CL+泥沙，IsoGSM 强迫）和 README 补充已写，未提交、未实跑（7920 上没有 IsoGSM 强迫）。

#### 原有问题：12 月冷启动 TRC_SIG（2026-10-28 定位，未修，等用户定）
- 现象：NEC 2005-12-01 冷启动，第 39 步起 patch 38 签名检查超差（误差 1.7e-12~1.2e-11，容差约 1e-12）；TRC_BAL 质量守恒始终为 0，只是成分漂移。
- 关植物水力（`DEF_USE_PLANTHYDRAULICS=.false.`）后完全消失（claude-8b0eb3f0）。
- 插桩（claude-2cecb56f，`.tmp/ct/p6e`）：第 3 层土壤液态水降到约 0（冻结），PHS 仍每步往这一层回流 -1.8e-5 mm、并从中吸水（如第 27 步 +4.5e-4）。这一层的物质与水对不上：有时液水为正而示踪物为 0（dev=-1），有时液水≈0 而示踪物还在（dev 达 1e14~1e26）。第 38 步后这部分错配随层间通量进入第 1、2 层（偏差 -3.3e-8、-4.7e-5），再由土壤蒸发（kind 2）带出，签名检查超差。
- 代码：`MOD_MT_SoilWater` 的根系回流块和残差回补块与 fork `MOD_Tracer_SoilWater` 逐字相同，是原有问题，不是移植引入的。

## 6. 验证矩阵

全部在 7920 上跑（`run-7920 --agent claude --source-local`，7920 预装编译环境，串行 make，每个配置独立副本）。

| 配置 | 宏 | 开关 | 比对 |
| --- | --- | --- | --- |
| A 基线 | 无 | — | `feature/riverlake-flow` 哈希 |
| B 宏关 | 无 | — | = A（cpp 输出 + 结果） |
| C 宏开全关 | `ConstituentTransport` | 全 false | = A |
| D 同位素 + 溶质 | 同上 | iso, solute | 守恒、拆分、5/9/17 进程 |
| E 泥沙 | + GridRiverLakeFlow | sediment | 同上 |
| F 甲烷 | + BGC | methane | 守恒、拆分；BGC 开甲烷关 = BGC 基线 |
| G 全开 + LULCC | 全部 | 全 true | 守恒、拆分 |
| SinglePoint | 有/无 | iso | 编译 + 1 年（不含河湖、泥沙） |

已知坑：AMZ 区域网验证不了梯级水库；np2 起不来，进程数对比用 5/9/17；SinglePoint 的 `river_hist_concatenate.x` 链接失败是旧问题。

## 7. 风险

1. **规模**：约 4.4 万行组分代码 + 约 430 处宿主钩子，所以拆成两个 PR（第 8 节第 5 条）。
2. **宿主漂移**：fork 的钩子挂在 fork 宿主上，上游宿主写法不同（截留方案、火、BGC），P2 是最容易出错的阶段，要靠钩子清单逐条核对。
3. **河湖结构变了**：PR #512 重构了推送和分汊求和，组分输运必须跟着用同一套顺序，否则进程数依赖会回来。
4. **PR #512 还没合并**：本分支叠在它上面，#512 再改就要 rebase；P3/P4 尽量等 #512 稳定后再做。
5. **甲烷变成运行期开关**：以前 `TRACER && BGC` 就编译甲烷，现在 BGC 构建都会编译甲烷代码，必须保证甲烷关时 BGC 结果逐位一致（配置 F）。
6. **侧车清单手写**：522 个字段，靠 CI 脚本兜底，上游以后加 `a_*` 时脚本要给出清楚的报错信息。

## 8. 用户决定（2026-10-27）

1. **命名**：确认。宏 `ConstituentTransport`，目录 `main/CONSTITUENT/`，模块前缀 `MOD_CT_*`，history 前缀 `f_ct_`，restart 前缀 `ct_`。
2. **同位素参数文件**：所有同位素共用一个文件 `DEF_ISOTOPE_PARAM_FILE`。文件里一个 `&nl_colm_isotope_parameter`，包含共享物理参数和按物种的数组；驱动场放 `&nl_colm_ct_forcing`。
3. **侧车 restart**：`MOD_HistAccRestart` 默认启用，跟组分开不开无关。它只保证拆分重启后的 history 正确，不改物理结果。代价是上游所有用户的 restart 目录里会多出侧车文件，这一点写进 PR 说明和登记表。
4. **旧类型名**：`conservative`、`reactive` 直接删，不留兼容。
5. **PR 拆分**：同一个分支分阶段提交，对上游开两个 PR：
   - ① 框架 + 同位素 + 溶质 + 河湖 + 泥沙替换 + LULCC（P1–P4、P6）；
   - ② 甲烷（P5），叠在 ① 上。

   PR 都由用户自测后自己开。
6. **子目录按 family 分**（2026-10-27）：物种自己的代码按 `family_id` 放子目录——`main/CONSTITUENT/FAMILY_ISOTOPE/`（O18、HDO、Registry、Frac）、`FAMILY_PARTICLE/`（泥沙）、`FAMILY_GAS/`（P5 的甲烷）；目录名保留 `FAMILY_` 前缀，与 `family_id` 常量同名；随水输运的公共代码（同位素和溶质共用）留在 `main/CONSTITUENT/` 根目录。溶质目前没有自己的代码，等有了再建 `FAMILY_SOLUTE/`。
7. **随水输运公共代码叫 `GENERIC_WATER`**（2026-10-27）：与 `STATE_OWNER_GENERIC_WATER` 同名，放 `main/CONSTITUENT/GENERIC_WATER/`（Precip、Forcing、ForcingInput、Snow、Evapo、EvapLimit、SoilWater、SoilInit、SpecialPatches、Conservation、RiverLake）；框架留根目录（Defs、Namelist、Lifecycle、Providers、LandPhase、Vars、Hist、Rest、HostFixes、RiverLakeCoupling）。`ct_uses_land_water_transport` 统一改名为 `ct_uses_generic_water`。Defs 里的溶质专用代码抽到 `FAMILY_SOLUTE/MOD_CT_Solute.F90`，同位素专用函数抽到 `FAMILY_ISOTOPE/`。顺序：先完成并提交 P4，再做一个“只移动、不改结果”的提交，7920 上验证与 P4 提交逐位一致。
8. **改名为“物质输送” MaterialTransport**（2026-10-27，替代第 1 条的命名）：宏 `MaterialTransport`，目录 `main/MATERIAL/`，模块 `MOD_MT_*`、子程序前缀 `mt_*`，namelist 组 `&nl_colm_material`、`DEF_MT_*`，history 前缀 `f_mt_`，restart 前缀 `mt_`；CI 参数 `MTON/MTOFF`。与第 7 条的纯移动一起做成机械改名提交，验证与 P4 提交逐位一致（变量名不同，按改名映射比对）。
9. **河湖流宿主修复，单独提交在本分支**（2026-10-28）：
   - 出口面宽 `outletwth` 取本单元河宽，与 CaMa-Flood 相同（比较过 min(w, w_next)，结果改变约为本单元宽度的两倍）。
   - 已建水库正上游单元的出流改为 CaMa 运动波（I1DAM=10），新模块 `MOD_Grid_ReservoirUpstream`。
   - 组中间重启的 history 续写修复：曾改成对所有运行生效，用户决定不修，保持只在开物质输送时生效。

**进度（2026-10-28）**：P1–P4 已推送；第 7、8 条的改名重组已提交（0f98622c，与 P4 逐位一致）；第 9 条两项已提交（090c5df8、620e521d）。待做：P5（甲烷进 `FAMILY_GAS/`）、P6、P7，以及 `main/MATERIAL/README.md`（族说明）。`sed_frc` 三类顺序 2026-10-28 已确认为砂、粉、黏（见宿主改动登记表 P4 节）；泥沙与 CaMa 的其余差异以后再处理。
