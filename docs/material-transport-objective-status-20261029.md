# 目标 A–D 现状覆盖矩阵（2026-10-29 10:14 CST，暂停等待 bug 清单）

> 说明：本表只陈述**现有证据可支持**的状态。“完成”指本目标要求已有终态验证证据；
> `partial` 指取消/中断前有可复核片段但**不构成通过**；`未启动` 指尚未运行。所有作业在 2026-10-29 10:03
> 按用户要求停止，等待其提供的明显 bug 清单。

## A. 组合验收缺陷（flood feedback + 运行期同位素强迫）

| 子项 | 状态 | 证据/说明 |
| --- | --- | --- |
| 根因定位（不靠放宽阈值） | **完成** | 逐点打印作业 `codex-01984a85…`：step 4 普通池重建为 0 水而同位素质量 9.0765e-05 kg/m² 仍在；旧代码清零并记数值源汇 → 超过 5e-5 硬限。见 `late-carrier-fix.md`、`residual-fix/v2/origin-dbg-lines.txt` |
| 根因修复（最终方案） | **完成（代码）** | `3b0af4f6`：存在 resolved late 载体水时保留 tracer 池随水混合输运；真正无载体塌陷仍记显式 residual；不放松阈值、不裁 δ、不新增 restart 字段。撤销了不符合水同位素机理的 `9f7c9a52`（`198da79`） |
| 组件级验证（触发案例 END/零诊断/自检） | **完成** | `codex-fc0e28c1…`：12 步 END、零 TRC 诊断、同位素 residue 全 0、9 项自检 PASS |
| combined10 全部门禁 PASS | **未达成（作业被取消）** | `codex-64f6cca1…` 取消前 10 案例到 step 64–70、**零** TRC_BAL/TRC_SIG/residual/abort（原缺陷 step 4 已越过），但 G1–G7 未运行；partial 见 `late-carrier-fix/partial/`。**恢复时需从头重跑全部 G1–G7** |
| 不破坏已验收结果（回归） | **partial** | 回归对照 `codex-8ff87f9e…` 取消前：iso 配置 3 天 history 3 文件/225 变量、最终 restart 3 文件/311 变量**逐位相同**；all 配置仅 3 个 tracer 变量差 ≤3.5e-18 kg/m²（1–8/26415 点），gridriver/hist restart 0 差异。**完整对照被取消，未完成** |
| p9r 十例在最终修复上的逐位重跑 | **未做** | p9r 验收是在 `c315c849` 上完成的；最终 `3b0af4f6` 上只有 partial 回归对照。恢复时若要严格“数值不回归”结论需重跑 p9r 类十例 |

## B. Hosted CI 矩阵

| 子项 | 状态 | 证据/说明 |
| --- | --- | --- |
| 旧 run `40725517808` 终态记录 | **完成** | 官方 REST：completed/**cancelled**；23 job = 20 success + 3 cancelled（ubuntu-latest×mpich×gcc12/13/14 触发 6h 上限）+ 0 failure；证据 `ci-matrix/`、`ci-matrix/EVIDENCE-run-40725517808.md` |
| CI 结构修复 | **完成** | `2bc7510a`：93 例构建分片到 23×4=92 job，单 job <6h，覆盖/触发不变 |
| 修复后提交触发并记录完整矩阵 | **partial** | Run 2 `40784782071`（head `2bc7510a`）：**92/92 success、0 failure/cancelled/timeout，墙钟 2h01m38s**，日志 marker 92/92 核对通过（`ci-matrix/run-fixed/`）。但之后生产修复 `3b0af4f6` 改变了最终 head；Run 3（head `613e0a49`）取消前仅 9 success → **最终 head 的完整矩阵未记录，待重跑** |

## C. 甲烷（真实数据、CROP+BGC+MT）

| 子项 | 状态 | 证据/说明 |
| --- | --- | --- |
| 真实构建 CROP+BGC+MT | **完成** | 首份真实 `make all`（`MT+BGC+CROP+LULC_IGBP_PFT+VG+GridRiverLakeFlow`），define.h 与二进制 SHA256 固定；`methane-runs/EVIDENCE.md` |
| mkinidata 冷启动 | **完成** | 4m12s、restart 642MB、END banner；irrigation-on 版 5m16s |
| 交付默认 hybrid 全球口径 | **修复完成 + 1/3 天 PASS；30 天/两年未完成** | guard 修复 `3778404f`（封顶分配替代 abort）；hybrid 1/3 天 END、guard abort 0、water-balance violation 0、CH4 residual ~1e-21、通量有限非零（`methane-runs/evidence-final/`） |
| 水稻甲烷 pilot/验收 | **1/3 天 PASS；30 天未完成** | rice 1/3 天 END、通量非零（−4.63e-4 / −4.45e-4 mol/s，sign 记录）；30 天作业取消前 336/1440 步零 fatal（`methane-runs/partial/`） |
| 两年全球运行 | **未启动** | 用户 2026-10-29 已批准：补 15 个 CRUJRA 边界文件（≈21.15 GiB，`bd:` 源）、**完整窗口 1995-01-01 00:00 → 1996-12-31 21:00（35,088 步）**、hybrid 默认、Y1 17,520 步 → Y2 17,568 步串行。W1 少 9h 方案不再采用。恢复计划见 `material-transport-resume-plan-20261029.md` |
| 旧 BD 结果 | **仅参照** | 旧作业 TRACER 谱系/spH off/lake off/rice off，不作为正确基准 |

## D. 证据/文档/推送与合规

| 子项 | 状态 | 证据/说明 |
| --- | --- | --- |
| 文档与状态更新 | **持续** | `material-transport-status.md`、continuation、各 evidence 目录均记录取消/暂停真实状态 |
| 推送个人 fork | **完成（到当前 head）** | head `01229c7f`；本地=远端；本轮后续提交继续推送 |
| 不合并 SYSU/master | **遵守** | PR #29 open/draft、merged=false，base=checkpoint 分支；未向上游推送 |
| 不删 checkpoint | **遵守** | 分支 `checkpoint/material-transport-base-20261028` 仍存在，sha `f9ea408e…` |
| 未完成/未验证项保留 | **完成** | 见本表与 `material-transport-status.md` 未完成清单；取消作业只按 partial 记录 |

## 恢复前必须由用户提供

1. 明显 bug 的**位置/现象/期望行为**清单；
2. 确认重启测试后，按 `material-transport-resume-plan-20261029.md` 的顺序执行；
3. 若 bug 修复会改变生产代码，恢复时以修复后的新 head 重新触发 CI 与全部验收。
