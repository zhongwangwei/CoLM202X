# 两年全球甲烷作业（W1 零传输窗口）配置、预算与 runbook —— 基于 7920 实测

状态：**配置/预算 READY；口径 NOT_READY（交付默认 hybrid 被 host-water 守卫阻断，见 EVIDENCE.md §5）**。
本文件不启动两年全球；启动由 Lead 决定。可运行替代口径 = `wetwat`；若决定修源码或接受其它口径，按 §5 调整 `CASE_MODE`。

源码：`CoLM-constituent` HEAD `e33d64f8`（生产树 = 交付 HEAD `10e9fbd2`）。
构建：`MaterialTransport+BGC+LULC_IGBP_PFT+vanGenuchten_Mualem_SOIL_MODEL+GridRiverLakeFlow`，LULCC off；sha256 见 EVIDENCE.md §1。
候选输入：`/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`（40,450 文件 / 26,618,457,915 B；UNVERIFIED_CANDIDATE）。
运行时依赖：`/media/zhwei/data02/zhwei/CoLMruntime/`（可读；与 BD 历史算例内容身份未证实）。

## 1. 窗口 W1（零新增边界文件；生产时间索引例程实测）

| 项 | 值 |
| --- | --- |
| 起点 | 1995-01-01 **03:00** |
| 终点 | 1996-12-31 **18:00** |
| 步数 | **35,070**（30 min 步长；含 1996 闰日） |
| 候选包外文件需求 | **0**（1994 末记录 / 1997 首记录都不被请求） |
| 相对原 W0（1995-01-01 00:00 → 1997-01-01 00:00，35,088 步） | 少 18 步 = 9 h 模拟时间 |

W1 是**改变实验定义**的替代方案（少 9 h），是否接受由 Lead/用户决定；精确原窗口需要补 15 个 CRUJRA 边界文件
（21.15 GiB，来源已实测存在于 BD，待用户授权）。

## 2. 运行配置（pilot 已验证的配置形状）

- 域：全球 2°（`DEF_nx_blocks=18, DEF_ny_blocks=9`，`DEF_GRIDBASED_lon_res=lat_res=2`），`DEF_LC_YEAR=2005`，
  `DEF_LAI_MONTHLY=.true.`、`DEF_LAI_CHANGE_YEARLY=.false.`。
- 模型开关：`DEF_USE_SoilInit` + soilstate.nc、`DEF_USE_CN_INIT` + cnsteadystate.nc、`DEF_USE_NITRIF=.true.`、
  `DEF_USE_VariablySaturatedFlow=.true.`、`DEF_UnitCatchment_file=…/grid_routing_data_15min.nc`、`DEF_Runoff_SCHEME=0`、
  `DEF_TOPMOD_method=0`、`DEF_USE_IRRIGATION=.false.`（全球两年不带水稻）。
- 甲烷模式：**`wetwat` 口径（可运行）**：`DEF_METHANE%inundation_mode='wetwat'`、`enable_rice_paddy=.false.`、
  `only_wetland=.false.`、`use_spatial_ph=.true.`、`allowlakeprod=.true.`、其余 HEAD 默认；
  `DEF_METHANE%host_water_tolerance=0.20`；主 namelist `DEF_USE_Dynamic_Wetland=.false.`（wetwat 模式要求）。
  - 若改用交付默认 `hybrid`：主 namelist `DEF_USE_Dynamic_Wetland=.true.`，但**当前源码会在 step 2-36 被守卫中止**（必须先在 §5 决策/修复后重跑）。
- restart：`DEF_WRST_FREQ='YEARLY'`；history：`DEF_HIST_FREQ='MONTHLY'`、`DEF_HIST_groupby='YEAR'`、压缩级别 1。
- **mkinidata 归属**：只在 Y1 作业里跑一次（起点 1995-01-01 **03:00** 的冷启动 restart）；Y2 从 Y1 的 YEARLY restart 续跑，
  两个作业使用同一构建二进制（sha256 固定），Y2 把 Y1 restart 目录作为只读输入声明。

## 3. 年度拆分（单作业 24 h 上限）

| 作业 | 窗口 | 步数 | 说明 |
| --- | --- | --- | --- |
| Y1 | 1995-01-01 03:00 → 1996-01-01 00:00 | **17,514** | mkinidata + colm；年界写 YEARLY restart（`y1_..._restart_1995-...`→ 供 Y2） |
| Y2 | 1996-01-01 00:00 → 1996-12-31 18:00 | **17,556** | colm only；复用 Y1 restart；末步写 restart/年度 history |
| 合计 | W1 | **35,070** | 与生产时间索引 probe 的 W1 步数一致 |

## 4. 实测资源与预算（2026-10-29 pilot，7920l，64 rank）

### 4.1 实测固定开销与步速

| 项 | 实测 | 来源 |
| --- | --- | --- |
| 构建（make all，含并行+增量串行补齐） | ≈3 min | 作业 A |
| mkinidata（40,429 landdata 文件，64 rank） | wall **4m12s**，15,714 core-s，restart 642 MB | 作业 A |
| mkinidata（irrigation on） | wall 5m16s，restart 661 MB | 作业 P3 |
| colm 初始化（含 BGC/CROP/landdata 读入） | ≈4–6.5 min | 作业 A/P3 |
| 步速（全球 2°，64 rank） | **1.84 s/step**（wet30d 运行中 log 采样：1,440 步/2,650 s，0.543 步/s；逐区间中位 1.90 s/step，区间 1.4–2.3；wet3d 为 2.30 s/step。当时服务器有其他 agent 作业、load≈100；独占时会更快） | `evidence/steprate/samples.tsv`、`results.json` |
| 单步 core 成本 | **≈118 core-s/step**（1.84 s×64；/usr/bin/time 直接除得 124.9 含固定开销，为保守上界） | `evidence/job-pilots3/results/phases/wet30d_full.json` |
| 30 天 run 实测总账 | 1,440 步，wall **53m12s**（3,191.6 s），User 179,540 s + Sys 344 s，CPU 5,640%；`CoLM Execution Completed`；内存 cgroup 峰值 **37.1 GB** | wet30d |
| 30 天输出体积（DAILY 诊断输出） | host hist 587.6 MB + material hist 25.1 MB + unitcat 426.4 MB = **1,039 MB**；restart 1.63 GB（1,430 文件） | wet30d；生产用 MONTHLY hist 时 history 约 1/30 |
| cgroup 峰值内存（64 rank） | **≈33 GB**（750 GB 的 4.4%）；单 rank maxRSS ≈1.7 GB | memwatch + cgroup |
| 最终 restart 写盘 | 1,430 文件 / 1.0–1.6 GB / ≈5–10 min（I/O 为主） | 作业 A/P3 |
| 输出体积 | host history ≈24 MB/日 slice、unitcat ≈18 MB/日、material ≈0.6 MB；restart ≈1.0–1.6 GB/checkpoint | 作业 P3 |

### 4.2 年度作业预算（W1，64 CPU）

- 年作业 wall ≈ init 0.15 h + steps×1.84 s/3600 + restart 写盘 0.15 h：
  - 按 30 天实测 1.84 s/step：Y1 17,514 步 ≈ **9.1 h**、Y2 17,556 步 ≈ **9.1 h**（合计 ≈18.2 h）；
  - 含固定开销上界 1.95 s/step：Y1/Y2 各 ≈ **9.6 h**；
  - 保守口径（负载重时 2.3 s/step）：Y1/Y2 各 ≈ **11.3 h**（合计 22.6 h）。
- 结论：**年度作业 9–11.5 h，均 < 24 h 上限；必须且只需拆两个年度作业。**
- 建议资源：**64 CPU / 128 GB 内存 / `--timeout 72000`（20 h，≈1.8–2.2× 预算）**；若负载高或保守取值，用 `--timeout 86400`。
- 两年 campaign 实测外推：Y1 + Y2 ≈ **18–23 h 计算墙钟**（另加 build 0.05 h、mkinidata 0.1 h），跨 2 个作业、每个 <24 h。
- mkinidata 只在 Y1（4–5 min）；Y2 不跑 mkinidata。
- 两年 campaign 预计墙钟：build 0.05 h + Y1 8–11 h + Y2 8–11 h + 分析 ≈ **17–23 h，跨 2 个 7920 作业**（外加一次 0.1 h 的 build 作业，或把 build 放进 Y1）。
- 风险：a) **口径阻塞**（§5，最高优先级）；b) 服务器同时负载（pilot 期间其他 agent 作业使步速慢约 2×，排期前先查 `uptime`/槽位）；
  c) 候选输入科学身份未证实、2005 landdata 用于 1995/1996、无 spinup → CH4 年总量不能当科学结论；d) restart 写盘慢，年度作业 timeout 已留余量。

## 5. 口径决策（Lead/用户；详见 EVIDENCE.md §5）

1. **交付默认 hybrid 在候选全球 2° 冷启动不可运行**：交付版 host-water guard（`MOD_MT_Methane_Physics.F90:1096-1124`）step 2 即中止；
   调大 tolerance 只推迟到 step 36（失配增长到 ≈1 孔隙体积，上限 <1）→ **调参数不是修复**。
2. 可选路线：
   a) **修源码**（推荐）：守卫计入 `wetwat/wdsrf/wa`，或把 `finundated` 按宿主可支撑量封顶而非中止；修后 hybrid 可用、水稻也可用。需 Lead/用户授权改仓库。
   b) **host spinup**：先跑 host-only（甲烷关）若干时段，让 `Dynamic_Wetland` 水位与土壤水积分平衡，再开甲烷；需要额外 run，且不保证消除结构不一致。
   c) **改用 `wetwat` 口径**（本文件 §2/§4 预算即按此给出）：立即可跑、物理不同（无 routing floodplain、无稻田淹水），与旧结果/交付默认不可比。
3. 水稻：任何方案在 scheme 6 guard 修复前都无法完成；irrigation-on restart 资产已生成（`evidence/` 与 §4.1），修复后可立即复用。

## 6. Runbook（Lead 决策后可直接用本目录脚本；不自动启动）

本目录 `run_year.sh` 已实现 Y1/Y2 两阶段（mkinidata 归属、YEARLY restart、验证）。示例提交（Y1，wetwat 口径）：

```sh
~/.local/bin/run-7920 submit --agent codex \
  --source-local /Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/methane10 \
  --cpus 64 --memory-gb 128 --timeout 72000 \
  --input-remote LANDDATA=/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/landdata \
  --input-remote FORCING=/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/forcing \
  --input-remote RT=/media/zhwei/data02/zhwei/CoLMruntime \
  --input-remote RAWROOT=/media/zhwei/data02/zhwei/CoLMrawdata \
  --input-remote BUILD=<Y1 build dir (run_year.sh 首次需自行 build，或用 EVIDENCE §1 的构建流程)> \
  --env CASE_TAG=y1 --env CASE_MODE=wetwat --env CASE_TOL=0.20 -- bash run_year.sh
```

- Y2：同一命令，`--env CASE_TAG=y2`，并加 `--input-remote MKINI_RESTART=<Y1 work>/pilot/y1/case/restart`。
- 首次运行需要先构建：本目录 `run.sh` 的 `build_coLM` 已冻结构建命令（也可直接复用本任务作业 A 的 build，sha256 见 EVIDENCE §1；
  注意 7920 任务 work 目录不会自动清理，路径 `/media/zhwei/data01/compute_runs/codex/codex-9db71e90e32846b59b3cb0d7ecd29137/work/build/mt-bgc-crop-pft` 仍可只读声明使用）。
- 验证：`run_year.sh` 末尾自动跑 `validate.py`（步数/最终时间戳/END/CH4 变量有限非零/平衡残差/文件清单），小 JSON 写 `$OBREF_JOB_RESULTS`。
- 年度作业前请确认：(i) 口径决策已定；(ii) 服务器负载/槽位；(iii) 结果证据只回取小 JSON（history/restart 不 fetch）。
