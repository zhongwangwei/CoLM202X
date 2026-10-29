# 暂停后的恢复运行计划（2026-10-29，待用户 bug 清单后执行）

状态：2026-10-29 10:03 CST 按用户要求 **停止全部在途测试**（用户已发现明显 bugs，清单待给）。
本文是“恢复时怎么跑”的**可执行计划**，不是已执行记录；在用户确认 bug 修复完成前**不启动任何作业**。
所有重计算走 `~/.local/bin/run-7920 submit --agent codex --source-local "$PWD" ...`，私有目录不在 Git。

## 0. 当前已取消/未启动的作业（恢复时需重跑）

| 作业 | 内容 | 取消时状态 | 恢复方式 |
| --- | --- | --- | --- |
| `codex-64f6cca13a0d487ba438bab454d6d413` | combined10 + G1–G7（源码 `3b0af4f6`） | cancelled/-15，10 案例到 step 64–70，**零 TRC 诊断**，G1–G7 未跑 | 修 bug 后新提交上**从头**跑全部门禁 |
| `codex-8ff87f9ecf494c0a9754cc9333553f11` | 回归对照 fixed vs `6d995438` | cancelled/-15；iso 案例 143/143 已 END 且逐位相同；all 案例 ~70/143 | 从头重跑（约 20–40 分钟） |
| `codex-923e86aad16940228b7def12bae54d4d` | 30 天水稻 | cancelled/-15；336/1440 步，零 fatal | 从头（或从 1995-152 restart 续跑后重新验收） |
| GitHub CI Run 3 `40801158901` | 最终 head 完整矩阵（分片 92 job） | completed/cancelled，取消前 9 success | 修 bug 后在新 head 重新 `workflow_dispatch` |
| 甲烷补数据 + 两年 Y1/Y2 | 用户已批准：补 15 文件（21.15 GiB）、完整窗口 35,088 步、hybrid 默认 | **未启动** | bug 清单处理完后启动 |

## 1. 恢复前置检查（每次启动作业前）

1. `git status` clean；`git rev-parse HEAD` = 远端 head；记录源码 commit。
2. 7920 负载：`ssh 7920l uptime`；确认没有他人/自己的在跑作业；高峰时降低并行或错峰。
3. 数据输入只读声明：`--input-remote RAW=... RT=... UCAT=... P6B=... ISOGSM2=...`（回归/组合验收）；
   甲烷另加 `LANDDATA/FORCING/BUILD`（复用已验收的 v2 build 与 irrigation-on restart）。
4. 作业超时按实测设置：combined 8h、regression 2h、rice30 6h、两年 Y1/Y2 各 24h；
   **Y2 必须等 Y1 成功终态**（不得并行）。

## 2. 恢复顺序（按依赖与价值）

1. **同位素修复回归（快）**：`residue_regress` 套件（fixed vs `6d995438`；A=JRA+iso+CL+sediment+routing+FB 3 天、
   B=JRA+iso 3 天），约 20–40 分钟；判定：宿主/网格/hist 状态逐位相同，tracer 差异逐项量化并归档。
2. **同位素组合验收（长）**：`combined10_fix` 套件（两套 LULCC 构建 + 10 案例 + G1–G7），48 CPU/8h；
   判定：G1–G7 全 PASS、验证报告 status=PASS、零 TRC 诊断、失败项为 0。
3. **CI 矩阵（GitHub）**：在修复后的最终 head 上 `workflow_dispatch`（92 job、约 2h）；判定 92/92 success。
4. **甲烷补数据**：rsync 15 文件（bd: → 候选 forcing 根），保存 argv/日志/尺寸对照/可读性/哈希抽查；
   生产时间索引确认精确窗口 35,088 步边界请求全部闭合。
5. **甲烷 30 天水稻**：复用 v2 build + irrigation-on restart，1440 步；判定 END、无 fatal、
   water-balance violation 0、CH4 residual 有限、rice 通道非零。
6. **甲烷两年 W1–full**：Y1（1995，17,520 步，含 mkinidata）→ Y2（1996，17,568 步，复用 Y1 YEARLY restart）。
   判定：两个作业 END、步数精确、水量平衡与 CH4 守恒/通量证据保存；旧 BD 只作参照。
7. 每项完成后：证据 JSON/日志摘录 → `docs/material-transport-evidence/20261028/` → 更新
   `material-transport-status.md` / continuation → 本地提交 → 推送到个人 fork 并核对远端 head。

## 3. 验收判定口径（不得放宽）

- combined10：G1–G7 按 `combined10_fix/validate.py` 原定义；**不得靠调大 `DEF_MT_RESID_ABORT_NBAD`、
  不得裁 δ、不得跳过案例**。
- 回归对照：逐位比较以 `tests/material_transport/compare.py` 为准；差异必须逐变量量化并解释；
  “新增诊断变量”与“数值差异”分开表述。
- CI：92/92 job success 才算完整矩阵通过；分片覆盖必须 93 例无重漏、`check_defineh.py` 每 job 全量。
- 甲烷：默认 hybrid、`host_water_tolerance` 保持 0.05；wetwat 只作对照；30 天/两年不得用短跑代替。

## 4. 若用户 bug 涉及代码修改

1. 先静态核对/最小复现（必要时 1-day 级最小作业，需用户同意“重启测试”后）；
2. 修改后先跑**对应最小验证**（触发案例/单配置），再跑上表完整套件；
3. 所有失败/取消/未完成原样保留在 docs 与证据目录，不写成通过。

## 5. 资源与时间预算（7920 实测）

| 项目 | 配置 | 实测/预估 |
| --- | --- | --- |
| combined10 全量 | 48 CPU / 96 GB / 8h 超时 | 约 3–5h（受负载影响；上次取消前 10 案例 ~2 步/分钟/案例） |
| 回归对照 | 24–32 CPU / 64 GB / 2h | 约 20–40 分钟 |
| CI 矩阵 | GitHub 92 job | 约 2h |
| 甲烷补数据 | 7920 登录节点或小作业 | 21.15 GiB @ ~82 MB/s ≈ 4–5 分钟 + 验证 |
| 30 天水稻 | 64 CPU / 128 GB / 6h | 1440 步；历史 30 天 wetwat ≈53 分钟，rice 更慢 |
| 两年 Y1 / Y2 | 64 CPU / 128 GB / 24h 各 | 17,520 / 17,568 步；实测 1.84 s/step → 约 9–11.5h 各 |
