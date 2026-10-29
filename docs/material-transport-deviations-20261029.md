# 对原始目标口径的偏离/待确认登记（2026-10-29 10:15 CST）

本文件列出本轮工作相对“原始目标 A–D 字面表述”的**有意偏离**与**待用户确认点**，避免把口径变化埋进结论里。
状态：全部测试自 2026-10-29 10:03 暂停（用户发现明显 bugs，清单待给）；未启动任何新作业。

## 偏离 1（需用户确认）：CI 矩阵从 23 job 变为 92 job（分片）

- 目标 B 原文：*“修复后提交再触发并记录完整 hosted CI 矩阵终态（23 job）”*。
- 实际：原 23-job 矩阵在 `ubuntu-latest × mpich × gcc{12,13,14}` 三个 job 上超过 GitHub **单 job 6h 硬上限**
  （旧 run `40725517808` 实测 6h00m14–15s 被 cancelled，20/23 success）。为完成“所有组合都真实构建”，
  提交 `2bc7510a` 把 93 例构建**分片**为 4 个 shard：矩阵变为 **23×4 = 92 job**；每个 cell 的 4 个 shard
  合起来仍覆盖全部 93 例，`check_defineh.py` 每个 job 全量执行；**不减少覆盖、不改触发过滤**。
- 结果：Run 2 `40784782071`（head `2bc7510a`）**92/92 success、0 failure/cancelled/timeout、墙钟 2h01m38s**；
  92/92 job 日志 marker 核对通过；分片用例数 24/23/23/23=93 无重漏。
- **待确认**：
  - 方案 A（当前）：接受 92-job 分片矩阵作为“完整矩阵通过”的等价证据（覆盖相同，job 结构不同）；
  - 方案 B：坚持以 23 job 为验收口径，则必须接受 3 个 mpich+ubuntu-latest 组合**无法在 6h 内完成**，
    矩阵只能记为“20/23 通过 + 3 个因平台硬上限未完成”。
  - 若用户选择 B，恢复时只需不改 workflow、重跑并记录 20/23；若选 A，最终 head 上需重跑 92-job 并记录 92/92。
- 另注意：最终 head 的矩阵（Run 3）在 10:03 被取消，**最终 head 目前没有完整矩阵终态记录**（无论 A/B）。

## 偏离 2（已有用户明确批准，替代目标原文）：甲烷输入窗口

- 目标 C 原文：*“输入用零传输 W1 窗口”*。
- 实际：用户 2026-10-29 明确批准 **“可以，补上”**：从 `bd:` 补 15 个 CRUJRA 边界文件（≈21.15 GiB），
  使用**完整窗口 1995-01-01 00:00 → 1996-12-31 21:00（35,088 步）**，两年全球用 **hybrid 默认口径**。
- 因此 W1（少 9h）方案作废；**数据补传与 Y1/Y2 尚未启动**（10:03 暂停）。恢复时按
  `material-transport-resume-plan-20261029.md` 执行，并在证据里记录补传来源/清单/校验。

## 偏离 3（未完成，不是偏离）：同位素验收状态

- combined10 G1–G7 在最终修复 `3b0af4f6` 上**未完成**（作业取消时只到 step 64–70，零 TRC 诊断，但 G1–G7 未跑）；
- p9r 十例逐位验收仍停留在 `c315c849`，**未在最终修复上重跑**；
- 修复回归对照只到 partial（iso 逐位相同、all 仅 ≤3.5e-18 tracer 差）。
- 这些都在 `material-transport-objective-status-20261029.md` 与状态表中标为未完成。

## 偏离 4（代码历史）：第一次修复方案被撤销

- `9f7c9a52` 曾把无载体同位素质量存入 `trc_surface_residue` 并复水再溶解；用户指出水同位素无此机理后，
  以 `198da79` 撤销，最终采用 `3b0af4f6`（随 resolved late 载体水输运）。撤销记录保留在
  `residual-fix.md`，作为失败方案证据。
- 生产代码最终只含两处修复（`MOD_MT_SoilWater.F90` 同位素 late-carrier、`MOD_MT_Methane_Physics.F90`
  甲烷 host-water 封顶）与一处 CI 结构改动（`build_CoLM_gnu.yml` 分片）；其余为 docs/evidence。

## 偏离 5（合规，无偏离）

- 只在个人 fork `feature/constituent-transport` 提交/推送；PR #29 open/draft、未合并；checkpoint 分支
  `checkpoint/material-transport-base-20261028`（`f9ea408e…`）保留未删；未向 SYSU/master 推送；
  未动历史 CI run（Run 3 为自己 dispatch 后按用户要求取消）。
