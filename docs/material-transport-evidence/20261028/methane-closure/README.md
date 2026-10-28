# 甲烷候选输入闭合与运行设计（只读分析，未运行模型）

2026-10-29 DeepSeek 接续调查。**本轮未启动任何甲烷模型、未传输数据、未修改仓库。**
候选输入仍是 **incomplete candidate package**，不得把旧甲烷结果当正确基准。完整原报告在
7920 的 `codex-ad6b48e0cb0044a9bbf0ce2903195234/results/`（审计）、
`codex-fa5a1a057cc848e3a2468dd9797a73c0/results/`（只读调查）与候选暂存根
`/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`。

- `GAP_ANALYSIS.md`：总报告（判定、15 个缺失边界文件清单、窗口改写实测、候选包身份、参数/宏、水稻与两年设计、用户决策点）。
- `PLAN_global2yr.md`：两年全球运行设计草案（**未运行**），含旧全球作业真实基准与资源外推。
- `PLAN_rice.md`：水稻甲烷最小真实运行设计草案（**未运行**）；CROP+BGC+MT 无任何 CI 组合，须先真实编译。
- `provenance.md`：候选包来源/身份与证据范围。
- `closure.json`：机器可读缺口清单（含逐文件请求索引）。
- `MANIFEST.sha256`：本地分析工件哈希。
- `boundary_request_table.md` / `namelist_compat.json` / `probe_windows.json` /
  `bd_old_run_benchmark.json` / `CRUJRA_local_abs.nml` / `run_g_global_ch4_flat_2yr_20260527_200428.lsf` /
  `bd_standard_ch4_parameter_historical.nml`：关键证据与历史配置只读快照。

## 关键判定

| 子项 | 判定 |
| --- | --- |
| 精确 1995-1996 输入闭合 | **NOT_READY**：15 个边界文件（1994 年 7 个末记录、1997 年 8 个首记录）不在候选包；7920 本机无 CRUJRA 数据。权威源已实测存在于 `bd:` 原路径，补传 21.15 GiB 需用户授权。 |
| 无传输替代 | 生产时间索引例程实测：窗口改为 1995-01-01 03:00 → 1996-12-31 18:00 时外部文件需求为 **0**（少 9 h 模拟时间）；只补 1997（11.28 GiB）或只补 1994（9.87 GiB）也可行。均需用户确认改变实验定义。 |
| 候选包身份 | 传输一致 40,450 文件 / 26,618,457,915 B；仍为 **UNVERIFIED_CANDIDATE**（不证明唯一/版本/科学正确）。 |
| 参数冻结 | **未冻结**：旧运行实测 spatial pH=off（常数 6.2）、`allowlakeprod=.false.`、rice=off；HEAD 标准参数默认 pH=on、lake=on、rice=on。旧 case namelist 含 15 个已删除 TRACER 键，不能直接跑。 |
| 运行时依赖 | 8 类路径可读；与 BD 历史算例**内容身份未证实**（无哈希）→ UNPROVEN。 |
| 水稻甲烷 | 设计就绪、可运行性未验证；候选 landdata 含 crop CFT/pH/lake_soilc。 |
| 两年全球估时 | 旧全球 2° 作业真实日志：720 ranks、mkinidata 4m28s、model 3h05m、35,088 步、成功；7920 @64CPU 约 34.5 h/两年 → **必须拆年度作业**（每年约 17.3 h）。内存需短跑标定。 |
