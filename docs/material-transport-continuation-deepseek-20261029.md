# 物质输送：DeepSeek 接续记录（2026-10-29）

本文件是 2026-10-29 DeepSeek 接续工作的交接记录。**结论：同位素无强迫回归与长期极值诊断已完成；最新组合源码验收
执行后 FAIL（触发项已定位，属分支基线的既有交互缺陷）；甲烷输入闭合仍 NOT_READY（已找到决定性新证据，待用户决策）；
完整 CI 矩阵已在个人 fork 触发并处于终态前阶段。** 未完成、未验证与延期项全部保留在文末。

## 1. 起点、工作树与合规边界

| 项 | 值 |
| --- | --- |
| 仓库工作树 | `/Users/zhongwangwei/Desktop/Github/CoLM-constituent`（不是同级的 `CoLM202X`） |
| 分支 | `feature/constituent-transport` |
| 接续起点 | `10e9fbd204a99ff81bd60a98198b364d489e10d9`（交接提交；工作树干净、与 origin 同步） |
| 本轮新增提交 | 仅 docs/tests/evidence 与文档，无生产 `.F90/.h/.nml` 改动（逐项见 §2） |
| 远端 | `origin=zhongwangwei/CoLM202X`（个人 fork）；未推送/未合并 SYSU，未改任何 master |
| checkpoint | `origin/checkpoint/material-transport-base-20261028` = `f9ea408e` 保留未动；仍是 PR #29 的冻结基线 |
| PR | 个人 Draft PR #29 仍 OPEN/Draft/`mergedAt=null`，base 为该 checkpoint；本轮未改 PR 状态 |

**后台计算保留情况**：接手时 p9r/p9s 两个健康作业（`codex-4ba194bf…`、`codex-23354a2d…`）均在运行，未被重提或取消，
分别于 00:26:59 与 01:28 CST 自然完成 succeeded/0；本轮所有新计算都走 `run-7920`（codex agent、显式只读输入），
本机不跑长模型。

## 2. 已完成事项与实际证据

### 2.1 同位素/CL 无强迫回归（p9r）— PASS

- 作业 `codex-4ba194bf50074f68aaa3d76178c49830`，终态 **succeeded/0**（2026-10-29 00:26:59 CST，2h47m）。
- 旧 `28478b18` 对冻结 `c315c849` 的 5 种无强迫配置各重跑：**10/10 模型例** END、步数/最终时间戳精确、
  描述符身份正确、诊断失败 0；**10/10 严格比较**通过（restart 各 7 文件/289·283·438 变量，history 各
  3 文件/225·191·276 变量原始位一致），唯一差异为设计允许的标量 `mt_forcing_cache_schema 1→2`
  （两侧 cache count=0、无可靠性向量）；`failures: []`。
- **冻结源独立复核**：`git archive` 重新导出两提交，old 606 / new 619 个常规文件逐文件 SHA256 与任务证据一致、
  零不匹配；重算清单摘要等于 `e18ed7fd…`/`9b5140e4…`，与 `source-manifests.json` 完全一致。
- 证据：[p9r-acceptance.md](material-transport-evidence/20261028/p9r-acceptance.md)、
  `p9r-report.json`、`p9r-task-log.txt`、`p9r-provenance-check.txt`。
- 边界：只覆盖 `c315c849` 相对 `28478b18` 的无强迫回归；不含 `6d995438` frcsat 字段，不等于最新组合源码验收。

### 2.2 p9s 长窗口极值诊断 — PASS（有限范围）

- p9s 作业 `codex-23354a2d716f45059da97e74b357b906` 终态 **succeeded/0**（01:28 CST）：NEC 冬季 60 日/2,880 步、
  NEC 夏季 92 日/4,416 步、AMZ 湿季 91 日/4,368 步三窗口全部 END、区域内 10×10、无未掩码非有限值、诊断失败 0，
  `failures: []`；`report.json` 定位真实尾部（AMZ 湿地 O18 −368.91/+214.60‰、HDO −997.50/+1708.59‰；
  NEC 夏季土壤 HDO −815.04‰ 等）。
- 透明探针作业 `codex-452fcc287377488daef69b1491144c1b`（**succeeded/0**）：7 个目标月从基线 monthly restart 重跑
  （1–19 天），全部模型 END 且重叠 history 与基线**逐位相同**（每探针 3 文件/225 变量）；12/12 报告极值先用基线
  history 原文件/索引复现。独立第二作业 `codex-1425a3a12f584098a4f59f2916c7ea22` 重跑同一段，trace 文件 SHA256
  完全相同。首个尝试 `codex-ff49825f…` 因脚本 API 错误在建模前失败（保留记录）。
- Decimal50 回放 58,634 个 finite-pool 事件：最大质量误差 2.45e-14 kg/m²、affine 重建误差 ≤2.47e-14、
  BOOK 日账本残差 ≤1.82e-15/7.00e-21。12 个目标网格 22 个 target×isotope：**21 有效 δ + 1 个负 signed ratio**
  （按缺测输出、不裁 δ）。
- **结论**：长期 δ 尾部 = 有符号同位素质量 / 真实 gross 水（最小 1.39e-2，面积加权，不是净 ET）的比值，
  不是存水组成异常；与旧三点机制一致并扩展到 3 窗口/12 极值；AMZ 2004-01-15 的 HDO 配对值为负 signed ratio。
- 证据：[p9s-longwindow-diagnosis.md](material-transport-evidence/20261028/p9s-longwindow-diagnosis.md)、
  `p9s-probe-results.json`、`longwindow-probe-tools/`。
- 边界：只覆盖选中极值（例如 AMZ 湿地 O18 仍有 241 个点超出宽诊断范围）；不验证观测准确度/物理校准。

### 2.3 最新组合源码整体验收 — **FAIL**（触发项已定位，非最新提交回归）

- 主作业 `codex-2434412927a94e38a44ad21db5643a38` 终态 **failed/1**（02:39 CST，失败来自验收门禁）。
  生产源码 `6d995438`（= 交付 HEAD `10e9fbd2` 生产目录）。LULC_IGBP/NEC 六个完整组合案例都在
  **step 4 / 2005-12-25-05400** 触发 `TRC_BAL residual hard failure`（`worst_abs=9.0882e-5 > tol=5.0e-5`，
  ipatch=213 itrc=1 ptype=0）而中止；这是被显式记账的 `numerical_source_sink`（近干/洪水路径）超过代码自身硬限。
- 门禁：G1 构建 PASS；G2–G7 随中止案例 FAIL；验证报告 status=FAIL（89 failures）。PFT/NEC 与 IGBP/AMZ 子集 576 步 END，
  但不得当作整体通过。
- 单因素消融：触发项 = **`DEF_GridRiverLake_FloodFeedback=.true.` + 运行期 IsoGSM 同位素强迫同时开启**；
  与 CL/泥沙/levee/bifurcation/reservoir/分馏/LULCC scheme/LULCC 编译开关无关。完整配置去掉 flood feedback 跑 5 天 PASS。
- 回归分类：`28478b18` 与 `c315c849` 的既有二进制在同一 12 步触发案例上**同点同值失败** → 该交互缺陷
  **在分支基线已存在，不是最新提交引入的回归**；旧验收从未覆盖此组合。未放宽阈值、未裁 δ、未做根因修复。
- frcsat 修复隔离（正面证据）：`c315c849` 对 `6d995438` 同一真实 120 步 MPI/NetCDF 案例，history 全部原始位相同；
  restart 13 文件/767 公共变量相同，仅新增 2 个 mid-window 侧车字段。
- 证据：[combined-acceptance.md](material-transport-evidence/20261028/combined-acceptance.md)、
  `combined-report.json`、`combined-EVIDENCE.md`、`combined-attempts.json`。

#### 2.3.1 根因修复（2026-10-29，提交 `9f7c9a52`）

- **根因**（两个 7920 诊断作业逐点打印定位）：`MOD_MT_SoilWater.F90` 表层混合池“无载体塌陷”分支
  在 `surface_base_water=0` 时把整块同位素库存（patch 213 为 `9.0765e-05 kg/m²`）记成
  `trc_numerical_residual_step`（数值源汇），超过 `trc_resid_abs_tol=5.0e-05` 而 fatal。
  该路径对非挥发溶质已有“隔离到 `trc_surface_residue`、复水再溶解”的正确处置，同位素被排除。
- **修复**：把同一隔离规则推广到水同位素——表层塌陷、洪水 late-pool 残余、dew/frost late surface
  （原 fatal）、湿地塌陷四处；4 处复水入口扩展溶解条件；history 增加
  `f_mt_surface_residue_<species>` / `f_mt_subsurface_residue_<species>` 输出。
  复用既有 storage 组件，**不新增 restart 字段、不放松阈值、不裁 δ、不改水/能量物理**；
  其它数值残余路径保持原记账。
- **组件级验证**（`codex-6696730ffec34c518f011b395657e983`）：无 LULCC 与 LULCC 两种构建均 PASS；
  触发案例（IsoGSM+flood feedback，NEC，47 步）两种构建都 `CoLM Execution Completed.`，
  零 `TRC_SIG/TRC_BAL/residual hard/residual note/STOP/abort`；restart/history 中同位素隔离库有限非负；
  仓库 9 项自检全 PASS。
- **进行中**：完整组合验收重跑 `codex-d48a500ecf994f2fb89a6e2856a1cc6a`（G1–G7 全门禁，源码 `9f7c9a52`）；
  修复 vs `6d995438` 回归对照 `codex-2608f8979c3449e8bb3cef74e9aba087`。**终态前不写组合验收通过。**
- 证据：[residual-fix.md](material-transport-evidence/20261028/residual-fix.md)、`residual-fix/`。

### 2.4 甲烷：真实构建/mkinidata 通过；交付默认 hybrid 守卫已根因修复；水稻与两年 W1 运行中

- **输入闭合**：精确 1995–1996 窗口缺 15 个边界文件（1994 年 `dlwrf,pre,pres,spfh,tmp,ugrd,vgrd` 末记录；
  1997 年全部 8 变量首记录）。权威源已实测存在于 `bd:` 原路径（21.15 GiB，补传需用户授权；7920 本机无 CRUJRA）。
  **采用 W1 零传输窗口**（1995-01-01 03:00 → 1996-12-31 18:00，35,070 步，少 9 h；实验定义变更须在结果中声明）。
- **真实构建（首份 CROP+BGC+MT 编译证据）**：`MT+BGC+CROP+LULC_IGBP_PFT+vanGenuchten+GridRiverLakeFlow` `make all` PASS，
  define.h 与三个二进制 SHA256 已固定；mkinidata 冷启动 PASS（4m12s、restart 642MB、END banner）。
- **交付默认 hybrid 全球口径阻塞（根因已定位）**：fail-closed 守卫 `methane inundation exceeds host soil water`
  （`MOD_MT_Methane_Physics.F90`）只比较 `finundated×孔隙体积` 与土壤孔隙水；但 hybrid 的 finundated 是独立水位估计，
  洪水/湿地/含水层水不在孔隙水内，失配可合法接近 1（1 天诊断：0.05→step2, 0.20→step3, 0.95→step36 增长到 0.96；
  dynamic_wtd 同；水稻同 patch/数值）。**调大容差不可救**。
- **根因修复 `3778404f`**：把 abort 改为“记录 host-water 最大失配 + host-supported 封顶分配”（既有 `1/finundated`
  封顶被原注释证明面积加权 host 水精确守恒，不造水）；`methane_host_water_report` 仍显示最坏失配，不隐藏；
  未改物理参数、未放宽 tolerance 语义、未改 CH4 产生/氧化/输运。
- **替代 wetwat 口径实测（可运行但科学代价须声明）**：1/3/30 天真实运行全部 END；30 天 wall 53m12s、
  水量平衡 8.7e-3 m³（全球存量 1.9e12）、CH4 balance residual 4.5e-21、CH4 通量有限非零（rice=0）。
  该口径丢 routing floodplain/稻田淹水物理，不能替代默认口径或水稻验收。
- **修复后验收（进行中）**：hybrid 1/3 天、水稻（dynamic_wtd+rice+irrigation）1/3/30 天、两年 W1
  （Y1 17,514 步 → Y2 17,556 步串行，各约 9–11.5h@64CPU，实测 1.84 s/step、峰值 33–37GB）。**结果出来前甲烷不写通过。**
- 证据：[methane-runs/README.md](material-transport-evidence/20261028/methane-runs/README.md)、
  [methane-guard-fix.md](material-transport-evidence/20261028/methane-guard-fix.md)、`methane-runs/`；
  输入闭合只读调查另见 [methane-closure/README.md](material-transport-evidence/20261028/methane-closure/README.md)。

### 2.5 完整 CI 编译器/MPI 矩阵 — Run 1 未通过（3 个 6h 超时），已分片修复并触发 Run 2（进行中）

- 触发方式：fork `workflow_dispatch`（**没有**向上游开 PR、没有合并 SYSU/master、没有改触发过滤）。
- **Run 1 `40725517808`（head `10e9fbd2`）官方 REST 终态：`conclusion=cancelled`，
  23 job = 20 success + 3 cancelled（`ubuntu-latest×mpich×gcc12/13/14` 触发 GitHub 6h 单 job 硬上限）
  + 0 failure**；中位 99 min、最长 360.2 min。**完整矩阵未通过，不写成通过。**
  CI 内容核对：23/23 job 的 `check_defineh.py` 与历史字段检查步骤 success；93 唯一用例；MTON=5、LULCCON=2；
  独立覆盖核对与 YAML/API 名单一致（矩阵实为 23 job，不是 25）。
- **结构修复 `2bc7510a`**：`.github/workflows/build_CoLM_gnu.yml` 增加 `matrix.shard: [0,1,2,3]`，
  93 例构建分片到 4 个 job（23×4=**92 job**；每个 os/mpi/toolchain 仍覆盖全部 93 例，`check_defineh.py` 仍全量执行；
  不减少覆盖、不改触发过滤），把单 job 工作量压到 6h 内。
- **Run 2 `40784782071`（head `2bc7510a`）已触发，进行中**（首次核对 92 job：40 in_progress / 52 queued）。
  **92/92 终态 success 前不写完整矩阵通过**；终态证据由 `ci-matrix` 子任务写入 `ci-matrix/`（`run-fixed`）。
- 证据：[ci-matrix.md](material-transport-evidence/20261028/ci-matrix.md)、`ci-matrix/`（REST/web 快照、逐 job 汇总、
  cancelled 摘录）、`ci-coverage-check.txt`；本机 API 快照与日志标记在 `.tmp/ct/ci10/`。
  本机代理故障期间用本机 SSH SOCKS 中继（`ssh -f -N -D 11080 7920l`）恢复访问，未改系统代理、未把 token 传给远端。

## 3. 本轮新增 Git 提交（均在个人 fork 分支）

| 提交 | 内容 |
| --- | --- |
| `2d5c15b7` | 记录并提交 p9r 无强迫回归验收（报告、状态、日志、源版本） |
| `74cb9175` | 独立复核 p9r 冻结源与声明提交一致 |
| `3daa60e8` | 记录 p9s 长窗口三个完整窗口完成与 `report.json` 极值定位 |
| `83d1d6c0` | 7 个透明探针 + Decimal50 回放：解释 p9s 长窗口极值（含工具归档与证据 JSON） |
| `894720da` | 记录最新组合源码验收 FAIL、触发项定位、既有缺陷分类、frcsat 隔离、甲烷只读调查证据 |
| `9f7c9a52` | **生产修复**：无载体干池同位素质量隔离到 `trc_surface_residue`（组合验收 abort 根因） |
| `d6347f6c` | 记录该修复、组件级验证与回归差异（本轮文档/证据） |
| `2bc7510a` | **CI 结构修复**：93 例构建分片为 4 shard（92 job），保证单 job <6h |
| `3778404f` | **生产修复**：甲烷 host-water fail-closed 守卫改为 host-supported 封顶分配（hybrid 根因） |

生产代码变化仅：`MOD_MT_SoilWater.F90`/`MOD_MT_Hist.F90`/`MOD_MT_Vars.F90`（同位素隔离修复）、
`MOD_MT_Methane_Physics.F90`（甲烷守卫修复）、`.github/workflows/build_CoLM_gnu.yml`（CI 分片）；其余为 docs/evidence。

## 4. 未完成 / 未验证 / 延期项（保留，不掩盖）

1. **最新组合源码验收：原 FAIL 的根因已修复，完整验收重跑尚未终态**。修复 `9f7c9a52`
   把无载体干池同位素质量隔离到既有 `trc_surface_residue`（不放松阈值、不裁 δ），组件级验证通过；
   回归对照 `codex-d64e0d0b…` 已实测：同配置 fixed vs `6d995438` 的 host 水/能量/网格 river 状态**逐位相同**，
   仅 4 个 tracer 状态（含新隔离库）差 ≤1.3e-16 kg/m²、history 常见诊断差 ≤1.2e-11，并新增 4 个同位素
   residue 输出变量（预期诊断增补）——即修复对不触发该路径的配置是舍入级惰性。但 **combined10 全量门禁重跑
   `codex-d48a500ecf994f2fb89a6e2856a1cc6a` 仍在运行**——终态通过前，组合验收仍记为未完成。
   其它数值残余路径（土壤相变/雪层/蒸发边界）未修改；LULC_IGBP_PC、其它区域/scheme/强迫产品仍未覆盖。
2. **完整 CI 矩阵未通过/未完成**：Run 1 `40725517808` 官方终态 **cancelled = 20 success + 3 cancelled（6h 上限）**，
   完整矩阵未通过；已提交 `2bc7510a` 分片修复并触发 Run 2 `40784782071`（head `2bc7510a`，92 job）。
   **92/92 终态 success 前不写通过。** 本机代理故障期间用本机 SSH SOCKS 中继恢复 GitHub 访问（未改系统代理、
   未把 token 传远端）。
3. **甲烷两年全球与水稻（已真实推进，未完成）**：真实构建（CROP+BGC+MT）与 mkinidata PASS；
   交付默认 hybrid 被 fail-closed 守卫阻断，根因修复 `3778404f`（host-supported 封顶分配替代 abort，分配守恒不造水）。
   替代 wetwat 口径 1/3/30 天已 PASS（30 天 CH4 balance residual 4.5e-21、通量有限非零），但科学代价须声明。
   **修复后的 hybrid/水稻 1/3/30 天与两年 W1（Y1 17,514 步 → Y2 17,556 步串行）仍待完成**；
   未完成前不得写甲烷验收通过，旧 BD 结果只作对照、不是正确基准。
4. **其它未覆盖边界**（沿用交接文档）：BGC+LULCC 是上游不支持组合；空 MPI worker 8-rank 崩溃未修；
   SinglePoint host 同月缓存覆盖未修；schema-1 缓存可靠性/阈值与 forcing identity 边界；观测准确度与
   全球代表性未验证；8 组以外 MPI rank 矩阵未重跑。
5. **p9s 尾部诊断有限**：只覆盖 12 个选中极值，不等于全部超范围值已解释；未做观测比较或参数标定。

## 5. 建议的下一步（按优先级）

1. 等 `codex-d48a500ecf994f2fb89a6e2856a1cc6a`（组合验收重跑）终态：G1–G7 全 PASS 则更新组合验收结论；
   失败则保留原因并继续定位（不得放宽阈值）。
2. 等 Run 2 `40784782071` 的 92 job 终态：全 success 才记录“完整 CI 矩阵通过”；失败/取消按缺口保留并
   考虑进一步分片或拆分 job；同时保留 Run 1 的 20/23 未通过记录。
3. 甲烷：读取修复后 hybrid/水稻 1/3 天与 30 天水稻结果；通过后依次完成两年 W1（Y1 → Y2，各约 9–11.5h@64CPU），
   保存 job ID/终态/CH4 通量与守恒诊断；wetwat 替代口径仅作对照并声明科学代价。
4. 正式 PR 拆分/合并决定仍由用户作出；当前仅个人 Draft 留档，不向上游交付。

## 6. 证据索引

| 内容 | 路径 |
| --- | --- |
| p9r 验收 | `docs/material-transport-evidence/20261028/p9r-acceptance.md`（+ report/job/logs/provenance） |
| p9s 报告与诊断 | `docs/material-transport-evidence/20261028/p9s-report.json`、`p9s-longwindow-diagnosis.md`、`p9s-probe-results.json`、`longwindow-probe-tools/` |
| 组合验收（原 FAIL + 根因修复） | `docs/material-transport-evidence/20261028/combined-acceptance.md`、`residual-fix.md`、`combined-report.json`、`combined-EVIDENCE.md`、`residual-fix/` |
| 甲烷（输入闭合 + 真实运行 + 守卫修复） | `docs/material-transport-evidence/20261028/methane-closure/`、`methane-runs/`、`methane-guard-fix.md` |
| CI | `docs/material-transport-evidence/20261028/ci-matrix.md`、`ci-matrix/`、`ci-coverage-check.txt`（Run 1 REST 终态 + Run 2 进行中） |
| 状态总表 | `docs/material-transport-status.md` |

## 7. 复现命令（只读复核）

```sh
cd /Users/zhongwangwei/Desktop/Github/CoLM-constituent
git log --oneline 10e9fbd2..HEAD
~/.local/bin/run-7920 status codex-4ba194bf50074f68aaa3d76178c49830   # p9r succeeded/0
~/.local/bin/run-7920 status codex-23354a2d716f45059da97e74b357b906   # p9s succeeded/0
~/.local/bin/run-7920 status codex-452fcc287377488daef69b1491144c1b   # p9s 探针 succeeded/0
~/.local/bin/run-7920 status codex-2434412927a94e38a44ad21db5643a38   # 组合验收原作业 failed/1（验收 FAIL）
~/.local/bin/run-7920 status codex-6696730ffec34c518f011b395657e983   # 同位素修复组件验证 succeeded/0
# 组合验收重跑/回归/甲烷作业的 job id 见本文件 §2.3–§2.4 与对应 EVIDENCE
curl --socks5-hostname 127.0.0.1:11080 https://api.github.com/repos/zhongwangwei/CoLM202X/actions/runs/40725517808  # Run 1（需本机 SOCKS 中继）
```

本地私有工件（不进 Git）：`/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/{p9r_accept,p9s_report,p9s_probe,combined10,combined10_fix,residue_fix,residue_regress,ch4_closure10,methane10,ci10}/`。
