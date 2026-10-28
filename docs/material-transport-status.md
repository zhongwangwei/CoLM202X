# 物质输送 Draft PR 状态

状态记录：2026-10-28 启动、2026-10-29 整理交接；**2026-10-29 起 DeepSeek 接续验收**，后台状态以注明的观测为准。本次累计 Draft PR 用于保存实现、历史决定和未完成清单，**不是可合并版本**。

用户已要求整理提交并交接 DeepSeek；接手入口为[交接文档](material-transport-handoff-deepseek-20261028.md)，小型验收证据见[证据目录](material-transport-evidence/20261028/README.md)。已有后台计算保留，不因交接重复提交。

本次仅在个人 fork `zhongwangwei/CoLM202X` 留档，来源分支为 `feature/constituent-transport`，草稿基线为 `checkpoint/material-transport-base-20261028`（冻结于 `f9ea408e`）。不合并至 `CoLM-SYSU/CoLM`，也不改动个人 fork 的 `master`。前置河湖流 [PR #512](https://github.com/CoLM-SYSU/CoLM/pull/512) 已合并。整理起点为 `28478b18`，当时相对该基线包含 19 个分支提交；下述最新修复尚待完整集成验收。

原计划将框架与甲烷拆为两个 PR。本次先累计留档，不改变原始决定；正式合并前是否拆分仍待确定。

`checkpoint/material-transport-base-20261028` 是 PR #29 当前仍引用的冻结对比基线，不是另一开发分支；个人 master 与它不同，因此目前保留。PR 仍为未合并 Draft，既没有向 SYSU 交付，也没有合并到个人 master。

## 已实现范围

- P1–P4：框架、history 累加量 restart、同位素、溶质、河湖输运及泥沙替换。
- P5：甲烷 provider、湿地 C/N 影子库、洪水耦合和地表数据预处理模块。
- P6：土地利用变化的物质池重映射、守恒核对及含水层参考载体处理。
- P7 主体：命名整理、README、示例和矩阵中发现的问题修复；CI 与最终验收尚未完成。
- 最近提交还包括区域河网裁剪位置调整和有限水库同位素蒸发的精确积分。

“已实现”不等于所有配置已验收。历史测试记录见下方存档，最新证明以逐项完成的验收为准。

## 最新降水强迫修复

仅对同时配置降水、水汽强迫的同位素启用相态回退：近干记录的雨取水汽液相平衡值，雪取最近可靠记录；只有降水率 **大于** `DEF_MT_FORCING_MEMORY_PRECIP`（默认 `2e-5 kg m-2 s-1`）才更新雪的可靠记忆。没有可靠记录时，雪取冰沉积平衡值。

缓存 schema 升为 2，增加 `mt_forcing_precip_reliable`。本轮接续审查确认并修正以下问题；代码与定向自检完成，冻结 p9q 的跨年模型、逐位对照及严格缓存补核已通过，最新组合源码仍需按交接清单核对。

| 问题 | 修正 |
| --- | --- |
| 无效或缺测的湿降水也进入弱雨平衡回退 | 只有 `DECODE_NEAR_DRY` 使用新雨回退；其他无效记录保留缓存回退 |
| 新阈值关联约束误拒绝未启用相态回退的旧配置 | `MEMORY >= MIN` 仅在启用相态回退时检查；保留参数非负校验 |
| LULCC 将多数可靠来源与不可靠估计混合后整体标为可靠 | 按用户决定，仅按可靠来源面积混合；无可靠来源则使用冰沉积平衡，并保存相态配置跨越 forcing teardown |
| 比较器只比较共同文件与变量，差异仍可成功退出 | 新比较器核对完整文件集、变量集、形状、类型和原始变量字节；失败返回非零，不比较属性 |
| 测试驱动可能吞掉编译、模型或比较失败 | 新 p9m 增加失败传播和逐案例完成、守恒、LULCC、比较门禁 |
| 同位素检查漏报未掩码的 NaN/Inf | 保留缺测掩码，计数真实非有限值并使检查失败 |
| d-excess 独立筛选两个物种后错配位置 | 在同一文件、同一形状上使用联合有效掩码后计算 |

### 已完成的轻量验证

- 真实生产例程提取后的 Fortran 自检：雨雪回退、阈值兼容、LULCC save/remap/restore 生命周期，以及 SAT/MEC 的可靠来源混合，均通过。
- 整理后的三组定向自检及 history 累加量清单检查在 7920 通过（`codex-404145c3b2594619aceefe9373d2616b`，终态 `succeeded`）；524 个字段与声明匹配。严格比较器和联合有效位置的 d-excess 自测均通过。
- p9m 的 IGBP 与 PFT 两种 LULCC 构建通过；语法检查和 `git diff --check` 通过。
- Fortran 定向自检的分馏系数使用桩；它不能替代完整 MPI、真实强迫或物理合理性验收。测试工具整理在 [tests/material_transport](../tests/material_transport/README.md)，不随仓库保存大数据或运行产物。

## 验收结果与仍在运行的任务

以下为已完成模拟的验收状态，不是自动更新的仪表盘。**模型计算结束、runner 终态成功、验收通过是三个不同条件。**

| 任务 | 目的 | 本次状态 |
| --- | --- | --- |
| p9k | 精确蒸发积分的长窗口合理性：NEC 冬季、夏季和 AMZ 湿季的 δ、d-excess 与守恒 | 模拟及严格复核已结束；70,227 张二维平面证明区域外无有效值，无未掩码 NaN/Inf；长期通量极值未全面解释，不能记为科学总通过 |
| p9l | 弱降水相态回退、schema-1 兼容、拆分续跑，以及无同位素强迫配置的回归 | 模拟及严格复核已结束；Jan4 split 完整位一致，legacy schema 兼容通过；P9B 混入更早 history 掩码变化，不是本次修复的单因素基线 |
| p9m / p9q | IsoGSM 2005-12-25 至 2006-01-05 跨年 LULCC；IGBP/PFT、SAT/MEC、断点与进程数对比 | 首轮 p9m 的 MEC 守恒失败已修。p9q 九模型及所有分段完成、诊断为零；四组完整原始位比较通过。原 runner 因冷启动缓存误判仍为 failed/1；新只读补核 succeeded/0、report=PASS，关闭该冻结配置的检查缺口，不改写旧终态，也不代表后来的 frcsat/LULCC 最新源码组合已测 |

任务号用于找回既有结果，不能据此代替读取结果：

- p9k：`claude-8d1bd4b1d9554b86a9c77e8d2b266458`
- p9l：`claude-102542ded83c43a78b73e8f320c67271`
- p9m：`codex-f92acf9040bf46b0bdf82d0782ba8c4f`
- p9q：`codex-1f593daee0854b2587cee991b7996e45`，`6666cdbb` 加叶片库存修复；不包含随后独立整理的 RegionClip 或进行中的 Vector 输出改动。
- p9q 只读补核：首个 `codex-5d025b6634c442778093eca0089788fd` 因错误要求 split 最终 namelist 保留显式 scheme 而 failed/1；实际 transition 日志为 MEC。修订审计 `codex-bc9d9644688a4ea0b68817319dcc242c` 为 succeeded/0，report=PASS；九例 108 个后续缓存日期与八次完整比较通过，没有重跑模型。缓存 checker 的 17 项正负 fixture 通过；仅显式初始日期允许 host-only 或 schema5/complete1/transport_count0 的空事务，后续必须有 schema2 cache。
- p9r：`codex-4ba194bf50074f68aaa3d76178c49830`，正确基线 `28478b18` 对冻结 `c315c849` 的五种无强迫配置，每种均重跑旧/新源码，保持原完整月或三日窗口。**2026-10-29 00:26:59 CST 终态 succeeded/0**：10/10 模型例 END、步数/最终时间戳精确、描述符身份正确、诊断失败 0；10/10 严格比较通过（restart 各 7 文件/289·283·438 变量，history 各 3 文件/225·191·276 变量原始位一致），唯一差异为设计允许的标量 `mt_forcing_cache_schema` 1→2（两侧 cache count=0、无可靠性向量）；`failures: []`，验证器自测通过。证据见[验收记录](material-transport-evidence/20261028/p9r-acceptance.md)。**边界：只覆盖 `c315c849` 相对 `28478b18` 的无强迫回归，不含 `6d995438` frcsat 诊断字段，也不等于最新组合源码验收。**
- p9s：`codex-23354a2d716f45059da97e74b357b906`，冻结 `c315c849` 重跑三个完整长窗口，补齐正确区域配置。**2026-10-29 01:28 CST 终态 succeeded/0**：NEC 冬季 2004-01-01..02-29（60 日/2,880 步）、NEC 夏季 05-01..07-31（92 日/4,416 步）、AMZ 湿季 01-01..03-31（91 日/4,368 步）三窗口全部 END、最终时间戳精确、本格 10×10、未掩码非有限值 0、诊断失败 0，`failures: []`。`report.json` 已定位真实尾部（AMZ 湿地蒸发 O18 −368.91/+214.60‰、HDO −997.50/+1708.59‰；NEC 夏季土壤蒸发 HDO −815.04‰、O18 −125.13‰；NEC 冬季土壤蒸发 O18 +119.10‰ 等，均只作宽范围诊断触发，不是物理通过线）。尾部解释由透明探针 `codex-452fcc287377488daef69b1491144c1b` 从对应月 restart 重跑目标月、逐位核对 history 透明性并做 Decimal50 gross-water/signed-mass 账本回放；在其结论前不宣称长期极值已解释。详见[p9s 报告证据](material-transport-evidence/20261028/p9s-report.json)。
- Vector p9t2：`codex-6d036be5acbd4ef5bfec97ab3c0605d3`，终态 succeeded/0。三组五日对照各覆盖 history 2 文件/193 变量、restart 115 文件/4690 变量，全量原始位一致；真实半日 NP5→9 的计数与完整日输出恢复通过。不含最后一个纯损坏输入护栏，正常数值路径未改变。
- Vector 动态 checkpoint 补核：`codex-ca35eb4112334ef991fe7ad03f6b8d00`，succeeded/0，明确 `2005-002-00000` 的 18 文件/750 变量全量原始位一致。原 p9t2 validator 误用目录名最大值选到 `const`，其 13 文件/364 变量不能证明动态状态；旧记录保留，由此补核关闭缺口。
- 实际 NetCDF/MPI sidecar：`codex-005e25fe7c6948e69428ca773dcc762b`，succeeded/0，最终 `6d995438` 代码的新旧计数兼容、三类指定损坏拒绝和 MPI 全退出通过。旧窗口的首日诊断缺测、翌日恢复；物理 restart 61 文件/3136 变量、history 排除预期 `f_frcsat` 差异后 2 文件/192 变量全量原始位一致。
- 甲烷候选输入只读审计：`codex-ad6b48e0cb0044a9bbf0ce2903195234`，succeeded/0，但结论为 `NOT_READY / UNVERIFIED_CANDIDATE`，不是模型验收；缺 15 个边界年份文件。

先前 mask-only Vector 快照 `codex-6629489828ac4f269bc375b628c9d31a` 因缺少有效样本归一化被主动取消，不作为通过证据。

p9k 的源快照遗漏 `regions.txt`，生成的区域边界字段为空，history 使用 `720 × 1440` 全球尺寸；后处理因此变慢，其配置和输出仍需复核，不能据此宣称完成了全球模拟或预期区域验收。p9k/p9l 的旧检查器存在上表中的漏洞，旧 `IDENTICAL` 和 runner 成功状态不能单独作为验收证据。

IsoGSM 2005–2006 的 12 个文件（约 25.4 GB）已传至 7920，供 p9m 使用；数据未加入仓库。

## 未完成清单

- [x] 关闭冻结 p9q 的冷启动缓存检查误判，并严格复核九模型跨年、断点、进程数、守恒与缓存。
- [x] 读取 p9r 完整结果并完成同位素无强迫回归验收（终态 succeeded/0，10/10 例 + 10/10 严格比较，`failures: []`，唯一样本差异为允许的 schema 标量迁移；见[验收记录](material-transport-evidence/20261028/p9r-acceptance.md)）。
- [ ] 读取 p9s 完整结果并完成同位素长期尾部账本诊断；最新组合源码的 LULCC 覆盖由独立组合验收作业确认。
- [x] 补 CI 的 `PC+BGC+MT` 与 IGBP/PFT 的 `LULCC+MT`，并修正精确案例名匹配和全参数传递。93 个配置的有效宏检查通过；四项定向 GNU 构建通过，完整 GitHub 编译器/MPI 矩阵未运行。
- [x] 补最新修复的非零 CL 回归：以 `28478b18` 对当前数值版本，NEC 两天 CL-only 与 isotope+CL 的四组实际运行、26 文件和 1,296 变量完整位比较通过；仅明确的缓存 schema `1→2` 例外。此项不替代强迫/LULCC 验证。
- [ ] 完成甲烷与 fork 的两年全球平场对照。候选 1995–1996 输入不闭合：按真实生产时间索引需要 1994 年除 dswrf 外的 7 个末记录文件及 1997 年全部 8 个首记录文件；旧参数与构建宏也未冻结。必须先补齐来源与输入闭合性，再短跑估时；没有启动新的全球甲烷模型。NEC/AMZ 三个月区域对照不能替代这一项。
- [x] 完成湖底碳、土壤 pH 两个 Aggregation 的真实预处理与区域裁剪验收。RegionClip 年份路径、元数据及向量子集修复通过真实数据检查；TOPMODEL 常数模式与合法无 PFT 块两处复审发现已修正，独立复审通过。未扩展至所有可选地表模式或空 MPI worker。
- [x] 空间 pH：真实生成、加载、强制裁剪及两天连续/拆分运行通过。27 个相关斑块全部非默认，19 restart 文件/1356 变量与 6 history 文件/706 变量逐位一致；甲烷闭合残差最大约 `2e-22 mol/m²/s`。这不是观测或长期平衡验证。
- [ ] 水稻甲烷真实运行验收。
- [x] SinglePoint 默认缓存一年验收：修复 host/material 共用单流缓存越界，全年 365 日/流；与直接写出 history 24 文件/2340 变量、restart 27 文件/1352 变量原始位一致。
- [x] Vector 完整宿主 history 与 NP5/NP9/5→9 续跑验收；真实半日计数保存和最终动态状态补核通过。范围为 UNSTRUCTURED/PC/VG、静态同位素与非零 CL、无 BGC/侧向河网，不替代强迫同位素、LULCC 或全部配置验收。
- [ ] 完成最新代码的回归汇总与最终审查，再决定正式 PR 的拆分和合并。

甲烷、CL、泥沙、河湖输运和非同位素 LULCC 的主体已有实现与历史测试记录，不应列为“尚未开发”。当前上游 LULCC 禁止 BGC，而甲烷需要 BGC，所以甲烷加 LULCC 属于不支持组合，不是本轮漏测。泥沙与 CaMa 的其余差异按历史决定延期，不在本次补修中扩大范围。

## 保留的兼容边界

- schema-1 缓存按既定策略全部视为可靠。这是迁移约定，不能证明旧缓存满足新记忆阈值；旧文件无法恢复这一来源信息。
- restart 的 forcing identity 不包含两个降水阈值。改阈值后续跑会沿用历史可靠标记，不会按新阈值重新判定；本轮未改变这一行为。
- 旧部分 history 窗口缺少新增的 `frcsat` 有效样本计数时，保留旧累计量及全部物理状态并警告；只将该窗口的 `f_frcsat` 输出为缺测，正常清零后恢复。不能从旧累计和反推干湖无效步数，亦不应为此拒绝整条旧物理续跑。

## 历史记录

以下文档保留原阶段的决定和证据；其中过时的“待做”状态以本页为准：

- [移植计划](superpowers/plans/2026-10-27-constituent-transport-port.md)
- [宿主修改登记表](superpowers/specs/2026-10-27-ct-host-changes.md)
- [宿主钩子清单](superpowers/specs/2026-10-27-ct-hook-inventory.md)
- [物质输送使用说明](../main/MATERIAL/README.md)

## 接续验收记录

- **p9r 同位素/CL 无强迫回归**（DeepSeek 接续，`codex-4ba194bf50074f68aaa3d76178c49830`）：2026-10-29 00:26:59 CST 终态 succeeded/0。旧 `28478b18` 对冻结 `c315c849` 的五种无强迫配置各重跑，10/10 例 END 且步数/最终时间戳精确、描述符身份正确、诊断失败 0；10/10 完整文件集/变量集/原始位比较通过，唯一差异为设计允许的标量 `mt_forcing_cache_schema 1→2`（两侧 count=0、无可靠性向量）；`failures: []`，验证器自测通过。不覆盖 `6d995438` frcsat 诊断字段或最新组合源码。详见[验收记录](material-transport-evidence/20261028/p9r-acceptance.md)。
- 定向编译 `codex-d479cc7d5bb04bf190119a0427536a24` 终态成功：PC+BGC+MT、IGBP+LULCC+MT、PFT+LULCC+MT 的 `make all`，以及 SinglePoint 的 `colm.x/mkinidata.x`。构建输入为 `3a1c391a` 加 CI 改动；不包含后续运行暴露问题的修复。首轮环境链接失败已保留，重试显式统一 Conda GNU/NetCDF 和编译选项，不改生产构建文件。
- 非零 CL 验收 `codex-577f4c4aa667427f8b04ceeca34e5946` 终态成功。每例 96 步，签名、质量与硬残差诊断为零；实际 CL 初值 `2e-5`、降水浓度 `3e-6`，每例有 11,742 个正 CL 值。
- p9l 严格复核确认 Jan4 拆分的 history/restart 完整位一致。与 P9B 的四组 history 不同，冻结源码证明 P9B 混入先前 `9d351369` 的 history patch 掩码阈值变化，不能当作本次降水修复的单因素基线；保留原始失败记录，不回滚正确的输出掩码。尚未逐点消融所有末位差异。
- p9k/p9l 模型均已完成，所有模型日志的守恒诊断为零；区域裁剪与科学极值的进一步复核仍在进行。无 NaN/Inf 不等于所有 δ 极值已解释。
- 首轮 SinglePoint 默认缓存真实运行在 Jan17 越界，根因是 host/material 共用单流缓存；后已修复并通过默认缓存全年/直接写出对照，不回写旧失败为成功。
- 首轮 Vector isotope 首日签名断言关闭 PHS 后仍复现，后定位并修复单桶冠层 donor；继而修复干湖输出有效性，最新限定范围验证见上表。
- 首轮真实 `mksrfdata` 湖底碳已满足 `580×OM`，随后旧 RegionClip 路径不兼容失败；该问题后已修复并完成预处理/裁剪复核，不能把原失败任务当成成功任务。

本次 Draft PR 只在个人 fork 的冻结留档基线上记录。该基线不匹配现有 GitHub Actions 的 master/colm30 触发过滤，因此上面的定向远程结果不代表托管 CI 矩阵已执行。

## 本次修复及输入边界

- SinglePoint 缓存隔离与 PC/PFT 单桶混合雨雪释放已完成独立双路复核；后者 17 个提取生产例程场景通过，7 个不受影响场景六项输出与旧版逐位一致。修 donor 后的第一轮 Vector 仅剩宿主 `f_frcsat` 差异；后续有效性与计数修复的完整结果见上表。
- LULCC 继承到湖/冰川的 signed 叶片 NSS 异常已经包含在步初库存；混入接收水池后清除叶片库存和失效记忆，避免下一步重复计入。LULCC 前后审计也补计该库存。24 个提取生产代码用例及真实旧 restart 两步重放 `codex-2e5c81aaa21e422a9e1f4c44378dbb25` 通过、NBAD=0 且无 SIG/BAL，独立代码审查通过。非零异常无安全载体、负接收库存或 fixed-signature 不兼容时明确停止，不凭空增水或丢弃异常；真实冰川及完整跨年组合不由两步湖泊重放证明。
- Vector `f_frcsat` 差异精确对应 131 个干湖元素、393 个日值；干湖分支未赋值，连续运行沿用旧值而 restart 初始化为缺测。修复只记录步初有效性并修改材料启用时的诊断累加/均值，不改原始 `frcsat` 或水文/BGC/甲烷消费者。有效样本计数随窗口重置及 sidecar 保存；真实 marker reader 新增逐项有限性检查，防止 MIN/MAX 掩盖混合 NaN。定向测试、独立审查、525 字段清单及真实多进程/半日续跑已通过，详见[独立复核](material-transport-evidence/20261028/vector-acceptance.md)。
- SinglePoint 原有 host 同月缓冲续跑覆盖问题未在窄修复中改变：材料记录正确，宿主 January 时间轴可能被覆盖。月底续跑数值对照通过，不能宣称任意同月 host 缓冲续跑已修复。`check_site_history.py` 检查时间/掩码/数值精确相等，原始位结论来自另行严格 `compare.py`，不混用两种证据。
- 三个短期蒸发通量极值已用真实 gross-water/signed-mass 账本和 Decimal50 回放：3600 个有限池事件最大质量差约 `8.32e-15`，三个目标日 δ 差最多 `6.68e-9‰`；负 HDO 比值来自净同位素吸收而被遮罩。结论只覆盖这些样本，不推广为长期所有极值的物理校准。
- RegionClip 修复按一次调用的真实土地覆盖年份选择一致的现代或旧目录布局，保留 packed/padded mesh、像元映射、`pctshared` 及变量属性；必需的湖碳/pH 缺失时明确停止。`codex-87e0723ea9a94135a18147b0d9820771` 验证 TOPMOD method 0、合法无 PFT 块及对应负例，混合块 541 个变量精确对应源子集；MT+BGC 与无 MT+LULCC 构建、四进程初始化另有通过证据。八进程空 worker 映射崩溃未修复；其他未实现的地表模式仍明确拒绝，旧目录支持不等于旧变量 schema 转换。
- 经用户授权，pH 单文件及全球候选包已传入 7920 新私有 staging，并完成内容一致性校验。全球包实际净文件大小为 26,618,457,915 字节、40,450 文件；目录项口径与文件净大小不同。**只证明传输一致，仍为 `UNVERIFIED_CANDIDATE`，不证明输入唯一、完整、科学正确或旧参数可重现。**
- 不能声称“这些数据只在 BD”。准确边界是在已查位置找到一组候选输入/旧结果；旧运行完成不等于可信基准，接近旧结果不是甲烷通过标准。长跑前仍须核验版本、时间/单位、网格、边界年份依赖、参数和守恒。
- 候选包实际配置是两份 namelist 与一份 LSF 脚本，不是三份 namelist。元数据、完整生产边界索引探针和依赖可读性已检查；这仍不是全科学场质量检查，完整[小型审计摘要](material-transport-evidence/20261028/methane-input-audit-summary.json)保留 `NOT_READY` 结论。
