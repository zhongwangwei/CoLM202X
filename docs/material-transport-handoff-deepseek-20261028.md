# 物质输送：DeepSeek 接续交接（2026-10-28）

交付整理跨至 2026-10-29；下文每项后台状态使用注明的观测时点，不把交付日期误当作业完成时间。

## 1. 交接结论与边界

本次目标是整理并保存已有实现、证据和未完成项，不是宣布全部验收完成。
用户要求停止扩展实现和新增模型计算，把剩余工作交给 DeepSeek；健康后台作业继续，不因会话/agent 重启而重提。
用户已明确选择：**GitHub 交接提交完成后暂停 Codex 的持续执行，保留后台计算**。
暂停后本会话不再主动收尾；DeepSeek 应主动查询台账中的任务，不依赖旧 shell/agent 自动通知。
本文源码基线为 `6d99543861c38c708f6c833e6349bcd3d7d3b879`；后续检查器和交接文档提交不等于重新验证过全部模型组合。

| 项目 | 权威位置 / 规则 |
| --- | --- |
| 本地工作树 | `/Users/zhongwangwei/Desktop/Github/CoLM-constituent`；目录名不是分支名 |
| 正确分支 | **`feature/constituent-transport`**，不是 `CoLM-constituent` |
| 个人远端 | `origin=https://github.com/zhongwangwei/CoLM202X.git` |
| 上游远端 | `upstream=https://github.com/CoLM-SYSU/CoLM.git`；另有 `colm202x=CoLM-SYSU/CoLM202X` |
| GitHub 留档 | [个人 Draft PR #29](https://github.com/zhongwangwei/CoLM202X/pull/29)，未合并 |
| 比较基线 | `checkpoint/material-transport-base-20261028`，冻结于 `f9ea408e` |
| 禁止动作 | 不推送/合并到 SYSU，不改变任何 `master`，不为整理记录重启健康作业 |

**checkpoint 目前有实际用途，因此暂留。** 本次查询 PR #29 仍为 OPEN/Draft、`mergedAt=null`，以它为 base；个人 `master=fcb2d1e` 与它不同。
它是 PR 的冻结比较基线，不是“已经合并”的证明。以后清理前必须重新核对 PR base、提交可达性与留档需求。

阅读入口：[当前状态](material-transport-status.md)、[测试说明](../tests/material_transport/README.md)、[材料模块说明](../main/MATERIAL/README.md)、[小型验收证据](material-transport-evidence/20261028/README.md)。
[原移植计划](superpowers/plans/2026-10-27-constituent-transport-port.md)及[宿主修改登记](superpowers/specs/2026-10-27-ct-host-changes.md)保留历史决定；旧目录名、旧“进行中”和旧宽松比较器结论不能覆盖本次逐项证据。

## 2. 已完成的代码与 family 组织

| 提交 | 生产位置 | 决策与已验证范围 |
| --- | --- | --- |
| `6666cdbb` | `main/MOD_Hist.F90`、`main/MATERIAL/MOD_MT_Hist.F90`、`GENERIC_WATER/MOD_MT_Precip.F90` | 隔离 host/material SinglePoint 单流缓存；PC/PFT 单桶冠层只有名义相态 donor 不足时才合并旧库存，避免混合雨雪释放伪失败 |
| `8a8c444e` | `mksrfdata/MOD_RegionClip.F90`、`MKSRFDATA.F90` | 使用一致土地覆盖年份/目录布局，保留 packed/padded mesh、属性与 `pctshared`；必需湖碳/pH 缺失明确失败，支持合法零 PFT 和 TOPMOD 常数模式 |
| `c315c849` | `GENERIC_WATER/MOD_MT_SpecialPatches.F90`、`main/MATERIAL/MOD_MT_Vars.F90` | 继承的 signed 叶片 NSS 库存进入特殊斑块水池后只清除一次，LULCC 两侧审计纳入该库存；不造水、不裁异常，拒绝无安全载体等不兼容情形 |
| `6d995438` | `main/CoLMMAIN.F90`、`MOD_Vars_1DFluxes.F90`、`MOD_Vars_1DAccFluxes.F90`、`MOD_Hist.F90`、`MOD_HistAccRestart.F90` 及字段生成器 | 干湖 `frcsat` 按有效步累加/归一化；计数与版本随 sidecar 保存，旧部分窗口未知计数只影响该诊断窗口；损坏 marker/count 明确拒绝 |

上表 `GENERIC_WATER/` 均位于 `main/MATERIAL/` 下。修复不重定义原始 `frcsat`，不改变其水文/BGC/甲烷消费者。
对应定向测试已纳入 `tests/material_transport/`，不是只有私有脚本才可检查。

| 目录 | 内容与边界 |
| --- | --- |
| `main/MATERIAL/FAMILY_ISOTOPE/` | 稳定水同位素、分馏和物种注册 |
| `main/MATERIAL/FAMILY_SOLUTE/` | 溶质，包括 CL；CL 的 `charge=-1` 是元数据，不代表实现电中性或电化学求解器 |
| `main/MATERIAL/FAMILY_PARTICLE/` | 泥沙/颗粒物；不是尚未开发 |
| `main/MATERIAL/FAMILY_GAS/` | 甲烷 provider 与相关过程 |
| `main/MATERIAL/GENERIC_WATER/` 与根目录 | 共享水载体输运、生命周期、守恒、history/restart；不是每个 family 各复制一套 |

降水修复已实现：仅“同位素 + 降水强迫 + 水汽强迫”齐全时启用相态回退。
近干雨用水汽液相平衡，雪用最近可靠降水记忆；仅降水率 **大于** `DEF_MT_FORCING_MEMORY_PRECIP=2e-5 kg m-2 s-1` 才更新记忆。
冷启动或 LULCC 无可靠来源时用冰沉积平衡；混合继承时**只混合可靠来源**，不用旧“可靠面积 ≥50% 即混入全部来源”的规则。
schema 2 新增 `mt_forcing_precip_reliable`；无效/缺测强迫与近干记录分开处理。

## 3. 任务台账：先看终态，再读验收内容

时间口径：2026-10-28 CST（UTC+8）。p9r/p9s 最新观测为 23:58:18（15:58:18 UTC）；其他行是已取到的终态，不是实时仪表盘。
**模型 END、runner 成功、配置验收通过、科学合理性通过是不同条件。** `failed/cancelled` 原记录永远保留。
下表结果均在第 6 节远程结果根内；换机器先恢复这些脚本/小报告，不依赖 Mac 私有目录。

| 任务 / 完整 ID | 源码与目的 | 已知状态 / 读取与下一步 |
| --- | --- | --- |
| p9q `codex-1f593daee0854b2587cee991b7996e45` | `6666cdbb` 加叶片修复；跨年 IsoGSM+LULCC 九例，不含后续 RegionClip/frcsat | 原 runner **failed/1** 保留；后续只读补核 **PASS**，九模型 END、SIG/BAL/LULCC/hard=0，四组 rawbits 全同，见下一段 |
| p9r `codex-4ba194bf50074f68aaa3d76178c49830` | 旧 `28478b18` 对新 `c315c849`，五种无强迫配置共十模型 | 最近 **running**；读 `report.json`，确认完整原窗口和全部变量比较，不用 P9B 代替基线 |
| p9s `codex-23354a2d716f45059da97e74b357b906` | `c315c849`，分馏开启、正确区域，三个完整长窗口 | 最近 **running**；读 **`report.json`、`audit.log`**，不是 `audit.json`；尚无 END，长期尾部未解释 |
| Vector p9t2 `codex-6d036be5acbd4ef5bfec97ab3c0605d3` | `c315c849` 加 frcsat 修复；最后损坏输入 guard 前的 HACC | **succeeded/0**；五天 full5/split5/full9/split59 与真实半日续跑通过，须连同下一行补核解释 |
| Vector 补核 `codex-ca35eb4112334ef991fe7ad03f6b8d00` | 只读原输出，不重跑模型 | **succeeded/0**；明确 `2005-002-00000`，18 动态 restart 文件/750 变量原始位相同 |
| 真实 sidecar IO `codex-005e25fe7c6948e69428ca773dcc762b` | 最终 `6d995438`，真实 NetCDF/MPI 旧格式和损坏输入 | **succeeded/0**；物理状态/正常诊断一致，legacy 一窗口缺测后恢复；三个坏格式指定 fatal 且全 MPI 退出 |
| 甲烷输入审计 `codex-ad6b48e0cb0044a9bbf0ce2903195234` | 元数据/时间请求审计，**不是甲烷模型计算** | **succeeded/0**，内容判定 **NOT_READY / UNVERIFIED_CANDIDATE**；读 `audit.json`、`boundary-requests.csv`，见第 5 节 |

### p9q：缺失并不自动等于失败，也不自动等于合法

四组比较为 IGBP MEC Jan1 split、Jan3 split、NP5/9，以及 PFT Jan1 split。
每个 IGBP 比较覆盖 restart 42 文件/1934 变量、history 6/450；PFT 覆盖 restart 57/2241、history 6/450，均完整 rawbits 相同。
原失败为初始 `2005-359-00000` 找不到 forcing cache：这是 2005-12-25 的 `mkinidata` 冷启动文件，早于 material/forcing 初始化；下一日开始已有 cache。
仅在显式指定的初始日期允许真 host-only，或 material schema=5、complete=1、transport_count=0 的已提交空事务（aquifer 元数据为 -1）。
后者由 `MOD_MT_Rest` 在 `ntransport<=0` 写元数据后返回，forcing writer 尚未执行；不能只因出现 `mt_` 字段就断言必须有 cache。
material restart schema 5 与 forcing cache schema 2 是不同协议；有实际 transport 的初始化文件及全部后续日期仍须完整有效 cache，不能全局跳过缺失。

新检查器 17 项正负 fixture 通过，已保存为 `tests/material_transport/check_forcing_cache.py`。
第一只读复核 **`codex-5d025b6634c442778093eca0089788fd`** 为 failed/1：split 的最终 `case.nml` 没有显式保留 `DEF_LULCC_SCHEME`，审计错误要求该文本必须存在。真实 transition 日志明确使用 MEC，不是新模型失败。
修订审计改为校验实际 transition 日志，若 namelist 显式给值则同时核对，缺项记录 null；第二只读复核 **`codex-bc9d9644688a4ea0b68817319dcc242c`** 已 **succeeded/0，report.json=PASS**。九例各 576 步（split 336+240 或 432+144）、全部实际 LULCC/分段完成/NBAD0、108 个后续完整 schema2 缓存日期及有效输出有限性通过；四组共八次完整 restart/history 比较原始位一致。两次都没有重跑模型，原失败不改写。
小型[通过报告](material-transport-evidence/20261028/crossyear-cache-results.json)和[审计交接](material-transport-evidence/20261028/crossyear-cache-acceptance.md)已进 Git；完整原报告在新分析 job 的 `results/report.json`。
本次补核不包含对冻结构建与源 manifest 的再次逐文件 SHA 审计，也不证明后来的 frcsat/LULCC 最新源码组合。
新工具入口见测试说明中的 `check_forcing_cache.py`；需独立核对模型完成和预期日期集合，缓存检查不替代这两个门禁。
原 p9m `codex-f92acf9040bf46b0bdf82d0782ba8c4f` 是真实叶片库存重复计入失败；不得与本次初始检查器误判混为一谈。
修复先通过 24 个提取生产例程用例及旧 restart 两步重放 `codex-2e5c81aaa21e422a9e1f4c44378dbb25`，后者 succeeded、NBAD=0。

### p9r / p9s：尚需收结果，不缩短窗口制造通过

p9r 五例：NEC/AMZ isotope、NEC CL 的 2005-12 月完整窗口，以及 NEC/AMZ isotope+CL+sediment+routing 的三日窗口；旧/新都完整运行。
比较必须覆盖文件集、变量集、shape、type 和原始位，仅明确记录 schema 1→2 例外；P9B 混入先前 history 掩码阈值变化，不是单因素基线。
23:58:18 快照：十例中六例已有 END（五个 old + new isoNEC）；new isoAMZ 到 Dec27/step1258，new CLNEC 到 Dec9/step416，两例 new allJan 尚待驱动启动。
p9r 配额 24 CPU / 48 GB / 14400 秒；保留原窗口，不因观察等待时间长而重提。

p9s：NEC 冬季 2004-01-01..02-29（60 日）、NEC 夏季 05-01..07-31（92 日）、AMZ 01-01..03-31（91 日）。
23:58:18 快照：NEC 冬季 step2218/Feb16，夏季 step2208/Jun15，AMZ step2238/Feb16；三例仍无 END。此前 SIG/BAL/hard=0，不代表后半程通过。
原作业配额 16 CPU / 32 GB / 21600 秒。正确 10×10 区域及 NBAD=0 已设置。
最终从 `report.json` 选真实极值/日期/网格，再用对应月 restart 做最小透明探针，不重跑整个 90 日来找晚期样本。
私有 `p9s/trace/instrument.py` 已记录 gross/signed 账本、真实日期及 raw `eair/esat`；**targets、runner、replay 尚未适配长期目标/月份**，不是即开即用工具。
探针先证明同段 history 全变量位一致，再核对 gross-water/signed-isotope 账本和 Decimal50 回放；不能拿净 ET 代替真实分母，不能裁 δ 或调参制造通过。
旧 `codex-4f7bfe583e234dfe8ac35115d15fee56` 仅解释三个短期点：3600 事件最大质量差约 `8.32e-15`，日 δ 差最大 `6.68e-9‰`。
其中负 HDO 是有符号净吸收后被遮罩，不是正向出水浓度；这不证明所有长期尾部、湿度闭合或观测精度。

### Vector / IO：限定范围已通过，保留检查器修正记录

Vector 五日各比较完整 2 history 文件/193 变量、115 restart 文件/4690 变量，原始位一致；时间轴/前缀保留、有限材料值和非零 CL 均检查。
真实半日 Jan1 noon NP5→9 保存 `hacc_nac=24`、有效计数和版本，完成整日的全部 193 history 变量与连续运行相同；额外半日端点明确单列。
原 validator 对 restart 目录取 `max()`，错误选中 `const`；其 13 文件/364 变量不能证明动态终态。
新只读补核明确日期后 18/750 全同，关闭这一检查缺口；原脚本/报告保留，旧远程输出不改写。
正常 Vector HACC SHA 为 `f4b1b773...`；最终 IO HACC SHA 为 `7a3dbb0e...`，差别是“count 存在但 version 缺失”的 guard，正常算术不变，其他五个相关源文件匹配最终版。

最终 IO：61 非累计侧车 restart 文件/3136 变量和 2 history 文件/192 变量物理比较一致；有意排除累计侧车及诊断 `f_frcsat` 的窗口差异。
去掉两个新字段的真实 legacy 输入：首日该诊断缺测，次日数据与 mask 恢复位一致，物理状态不变。
真实 orphan count、missing count、mixed-NaN marker 各在模型步前报对应 fatal，检查全部 MPI 进程退出，不以任意非零退出替代验收。
定向生产片段覆盖 MT-off/on 各 14 场景、`-O0/-O2` 十类损坏输入；实际字段生成器检查 525 项匹配。
范围仅 UNSTRUCTURED/PC/VG、静态同位素、非零 CL、无 BGC/横向河道；不证明强迫/LULCC/科学精度或完整配置矩阵。
原 mask-only `codex-6629489828ac4f269bc375b628c9d31a` 为 **cancelled/-15**，不是通过。

## 4. 已关闭的其他验收与历史遗留

| 项目 | 证据 / 范围 |
| --- | --- |
| 非零 CL | `codex-577f4c4aa667427f8b04ceeca34e5946`；旧 `28478b18` 对当时 `3a1c391a`，四模型各两日，26 文件/1296 变量完整位比较通过 |
| SinglePoint 年缓存 | `codex-630bcd15e7c841309ecdc5ebfa8652bc`；365 日缓存对直接写出，24 history/2340 变量、27 restart/1352 变量原始位相同 |
| 湖碳/pH 生成与裁剪 | `codex-ba20ec3d796349c0b7a84fd2a5bd3f68`；真实数据生成/裁剪，不只是桩测试 |
| RegionClip 独立复核 | `codex-87e0723ea9a94135a18147b0d9820771`；541 变量精确子集、TOPMOD0/合法零 PFT 及负例，MT+BGC/无 MT+LULCC 构建、4-rank 初始化通过 |
| 空间 pH 续跑 | `codex-b0d9511d68ed435ca9f6c669fda78ece`；27 非默认 patch，两日连续/拆分 19 restart/1356 变量、6 history/706 变量原始位相同，甲烷闭合残差约 `2e-22 mol/m²/s` |
| CI 定向覆盖 | 93 个宏配置、四项 GNU 构建通过；完整 hosted compiler/MPI 矩阵未跑，冻结 PR base 不匹配原 master/colm30 触发过滤 |

旧 p9k `claude-8d1bd4b1d9554b86a9c77e8d2b266458` 漏 `regions.txt`，得到全球尺寸输出，不等于完成全球模拟。
严格裁剪复核 `codex-6a7a6fbef16a4927b8a2a9d39cda8b53` 检查 70227 平面、无区域外有效值和未掩码 NaN/Inf；长期极值仍待解释。
旧 p9l `claude-102542ded83c43a78b73e8f320c67271` 的 Jan4 连续/拆分完整位比较及 schema1 导入通过；P9B 回归不能沿用为正确基线。
早期 5/9/17-rank 矩阵有历史记录，但旧比较器存在漏报，不能自动宣称最新源码完整矩阵通过。

## 5. 甲烷候选输入：传输完成不等于数据闭合

候选暂存根：`/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/`。
传输核对 40,450 文件 / 26,618,457,915 字节，只证明指定传输完整；不证明唯一、年份边界齐全、版本正确或科学 ground truth。
审计覆盖 16 个 1995–1996 CRUJRA 文件、两份旧 history 的元数据/坐标及生产时间请求。
生产 reader 还请求 **15 个缺失的边界文件**：1994 年除 dswrf 外七变量的末记录，1997 年八变量的首记录；`endyr` 只警告，不截断读取。
旧包是 **2 份 namelist + 1 份 LSF 脚本**，不是 3 份 namelist；旧甲烷参数和构建宏仍未冻结。
`/media/zhwei/data02/zhwei/CoLMruntime/` 下 soilstate、cnsteadystate、grid_routing、nitrif、ndep、crop、fire、snicar 路径已证实可读，尚未证实与历史算例内容相同。
先寻找 7920 已有匹配边界数据、核清来源/坐标/日历/缺测和宏参数，再决定 pilot 与两年运行；新增数据传输需按指定范围授权，不自行补传整套目录。
已有 pH 单文件在 `mt-ph-20261028-2053-1793191930/soil/PHH2O1.nc`（同一 data-stage 根），2,903,280,596 字节；SHA256 `bac6936e6b3322b3fe2eff5bdae9d8a0f8cc3c5ef9181c707c57b85396b7e965`。
**本轮未启动全球两年或水稻甲烷模型。** 元数据审计成功不能计为甲烷模型验收成功。

## 6. 接手操作与证据恢复

先核对本机工作树，不把另一个 `CoLM202X` 工作树误当当前代码：

```sh
cd /Users/zhongwangwei/Desktop/Github/CoLM-constituent
git status --short
git branch --show-current
git rev-parse HEAD
git remote -v
```

以下 CLI 参数已核对当前帮助；用台账中的完整 ID 替换 `JOB_ID`。`wait` 可能长时间等待，不启动新模型。

```sh
~/.local/bin/run-7920 status JOB_ID
~/.local/bin/run-7920 logs JOB_ID
~/.local/bin/run-7920 wait JOB_ID
~/.local/bin/run-7920 fetch JOB_ID --out NEW_EMPTY_LOCAL_DIRECTORY
```

**不要直接 fetch 大型模型任务。** 它复制整个 results；先通过 SSH 读取/定向复制小 JSON、文本和 manifest，避免把全部 NetCDF 拉回 Mac。
Codex 任务根为 `/media/zhwei/data01/compute_runs/codex/<JOB_ID>/`；`results/` 是输出，`work/` 是提交时冻结的代码/脚本。
Claude 任务归 `claude/`；更早位置可能不同，以 runner 返回的路径为准。
本机私有工件在 `/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/`，**不进 Git**：

| 子目录 | 用途 |
| --- | --- |
| `p9q_verify/` | 初始 cache 判定和只读补核，最终以新增报告为准 |
| `p9r/`、`p9s/` | 原 driver、validator/audit、源码 manifest |
| `p9s_tail_followup/HANDOFF.md` | 长窗口截至 23:46:53 的状态、探针剩余适配 |
| `iso_trace9/REPORT.md` | 三个短期点账本，不外推为长期全部通过 |
| `p9t2/`、`frcsat_io9/`、`vector_acceptance9/HANDOFF.md` | Vector/真实 IO 与明确动态日期补核 |
| `ch4_input_audit9/evidence/audit.json` | 甲烷候选缺口和时间请求证据 |

DeepSeek 换机器后从对应 runner `work/` 取原执行脚本，从 `results/` 取小报告。
长期探针 `trace/instrument.py` 是任务提交后在 Mac `.tmp` 准备的文件，不在最初远程快照；现已归档到 `docs/material-transport-evidence/20261028/longwindow-tools/`，换机器从该归档恢复，不能假定旧 runner `work/` 已包含。
旧结果只读，补核写新私有目录并保留原报告/脚本哈希。

### 远程计算规约与编译环境

先读 `~/.local/share/run-7920/README.md`。所有预计超过 2 分钟/8GB 的扫描、完整构建和模型运行走 `run-7920`、SSH 别名 `7920l`，不能悄悄改为本地重计算。
提交明确选择 `--source-local` 或 `--source-remote` 权威来源，不覆盖服务器共享工作副本；使用已配置 `--agent codex`，当前 CLI 只接受 codex/claude，不能凭模型名发明新 agent 值。
明确 `--cpus`、`--memory-gb`、`--timeout`；现有限额 64 CPU/750GB/86400s、10 槽，大任务先看负载。
先盘点输入，用 `check-data` 预检；每项外部输入通过 `--input-remote NAME=PATH` 只读声明，脚本用 `$OBREF_INPUT_NAME`，输出用 `$OBREF_JOB_RESULTS`。
共享 data01/data02 与代码保持只读；新输入只进已授权的唯一私有 data-stage 目录，不占根盘、不传凭证、不把大数据塞进 source。
保存提交 ID；查询超时/SSH 断线不是任务失败，不能盲目重提；取消只用 `run-7920 cancel JOB_ID`，不按进程名杀、不动他人任务。

既有远程成功构建使用同一 Conda GNU/MPI/NetCDF；旧 Makeoptions 的系统路径曾造成链接失败。以下为已用构建形式，本次文档整理没有再次执行重型构建：

```sh
export PATH=/home/zhwei/miniconda3/bin:$PATH OMP_NUM_THREADS=1
make all FF=mpif90 \
  NETCDF_INC="$(nf-config --includedir)" \
  NETCDF_LIB="$(nc-config --libdir)" \
  FOPTS='-fdefault-real-8 -fdefault-double-8 -ffree-form -O2 -cpp -ffree-line-length-0 -fallow-argument-mismatch -mcmodel=medium'
```

此命令只能在私有远程 job 的正确配置目录中使用；SinglePoint 既有验收只建 `colm.x mkinidata.x`，不能把其他附属目标的历史链接问题隐藏成 `make all` 成功。
定向检查命令与依赖以[测试 README](../tests/material_transport/README.md)为准；提取生产例程的桩测试不替代真实 MPI/数据运行。

## 7. 剩余工作、兼容边界与完成判据

建议按顺序接续，不另开无关实现：

1. 阅读已通过的 p9q 新 cache 补核，收齐 p9r/p9s 原作业终态及报告；核对完整窗口、日期集、变量集、退出码和原始失败记录，不重复提交旧模型。
2. 对 p9s 真实尾部适配最小探针和 replay，先证明透明性再作科学解释；有限数值/守恒不是观测准确度验收。
3. 评估并安排最终 `6d995438` 及交付 HEAD 的整体组合回归，尤其新 frcsat 与强迫/LULCC 同开的完整源码；旧 p9q 不能独自证明该组合。
4. 完成 hosted CI 编译器/MPI 矩阵及最终审查，明确当前 PR base 的触发策略；不得靠改成 SYSU/master 合并来触发。
5. 同位素优先事项关闭后，闭合甲烷候选输入/宏/参数，短跑估时，再按证据设计两年全球对照与水稻甲烷真实运行；区域三个月不替代全球两年。
6. 更新状态和交接证据，再由用户决定正式 PR 拆分/合并。当前是个人 Draft 留档，不是正式上游合并授权。

以下是已知边界，不应被“全部完成”掩盖：

- schema1 缓存全部当可靠是迁移约定；不能恢复旧记录是否满足新阈值。
- forcing identity 不包含两个降水阈值；改阈值续跑仍沿用旧可靠标记，不重新判定。
- 旧部分 history 窗口无 frcsat count 时仅该窗口缺测，清零后恢复；不能从旧累计量推回有效步数。
- 原 host SinglePoint 同月缓存续跑可能覆盖 January 时间轴；本轮修复材料共享缓存，不等于修复该宿主限制。
- 8-rank 空 MPI worker 映射崩溃未修；4-rank 通过不能代替它。
- BGC+LULCC 是当前上游不支持组合；甲烷依赖 BGC，因此不是可运行但漏测的一格。
- 泥沙/CaMa 其余差异按原决定延期；CL、泥沙、甲烷主体都已实现，并非只剩同位素之外的“空文件夹”。
- 日诊断极值、湿度闭合、长期科学合理性与观测精度尚无全面结论；不裁 δ、不调阈值迎合 RMSE，不把 finite 当 scientific pass。

完成判据：每项承诺都能对应冻结源码/配置、完整输入身份、真实终态和严格结果报告；失败有保留原因，未支持/延期/未验证有独立列表。
在上述未完成项关闭前，允许交付代码与交接，**不允许声称整个物质输送项目、全部同位素或全部甲烷验收已经完成**。
