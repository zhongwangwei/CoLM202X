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

### 2.4 甲烷候选输入闭合与运行设计 — 只读调查，仍 NOT_READY

- **未运行任何甲烷模型、未传输数据。** 精确 1995–1996 窗口缺 15 个边界文件（1994 年 `dlwrf,pre,pres,spfh,tmp,
  ugrd,vgrd` 末记录；1997 年全部 8 变量首记录）。**权威源已实测存在于 `bd:` 原传输路径**
  （`/share/home/dq013/zhwei/colm/data/CoLM_Forcing/crujra_v2.5/`，尺寸与 1995 同变量一致，合计 21.15 GiB），
  补传通道现成、**待用户授权**；7920 本机无 CRUJRA 数据。
- **零传输替代**：生产时间索引例程实测窗口 1995-01-01 03:00 → 1996-12-31 18:00 时新增外部文件需求为 0
  （少 9 h 模拟时间，需用户确认实验定义）；只补 1997（11.28 GiB）或只补 1994（9.87 GiB）也可行。
- 候选包传输一致 40,450 文件/26,618,457,915 B，但仍是 **UNVERIFIED_CANDIDATE**；参数未冻结
  （旧运行实测 spatial pH=off、`allowlakeprod=.false.`、rice=off；HEAD 默认相反；旧 case namelist 15 个 TRACER 键已删除）。
- 两年全球有真实旧作业基准（720 ranks、35,088 步、mkinidata 4m28s、model 3h05m、成功；7920 @64CPU ≈34.5 h/两年，
  **必须拆年度作业**，每年约 17.3 h）；内存需短跑标定。水稻设计就绪但 CROP+BGC+MT 无任何 CI 组合，需先真实编译。
- 证据：[methane-closure/README.md](material-transport-evidence/20261028/methane-closure/README.md) 及同目录
  `GAP_ANALYSIS.md`、`PLAN_global2yr.md`、`PLAN_rice.md`、`closure.json`。

### 2.5 完整 CI 编译器/MPI 矩阵 — 已在个人 fork 触发（终态前状态，见 §4）

- 触发方式：fork `workflow_dispatch`（head `10e9fbd2`，与交付 HEAD 的生产目录一致），**没有**向上游开 PR、
  没有合并 SYSU/master、没有改 workflow 触发过滤——这正面解决“冻结 PR base 不匹配 master/colm30 触发”的限制。
- run `40725517808`：https://github.com/zhongwangwei/CoLM202X/actions/runs/40725517808
  实际矩阵展开为 **23 个 job**（YAML：2 OS × 3 MPI × 6 GCC = 36，exclude 去重后 23；不是 25）。
- 截至本文最后核对（2026-10-29 03:29 CST）：**已完成 16/23 全部 success、0 失败、7 个 mpich job 仍 in_progress**（已跑约 3.0 h）；CI 内容核对已确认
  `check_defineh.py` PASS（93 唯一用例）、525 字段历史累加量检查、93 例逐例 `make clean && make all`、
  5 个 MTON/2 个 LULCCON 宏出现次数正确。mpich 历史上单例 3h06–5h48，可能拖到 6h 单 job 上限（06:24 CST）。
- 终态判定口径：全部 job success 才算完整矩阵通过；若有 job 被 6h 超时取消或失败，按下述“未完成/缺口”保留。
- 证据：[ci-matrix.md](material-transport-evidence/20261028/ci-matrix.md)、`ci-coverage-check.txt`；本地 API 快照与日志标记在 `.tmp/ct/ci10/`。

## 3. 本轮新增 Git 提交（均在个人 fork 分支）

| 提交 | 内容 |
| --- | --- |
| `2d5c15b7` | 记录并提交 p9r 无强迫回归验收（报告、状态、日志、源版本） |
| `74cb9175` | 独立复核 p9r 冻结源与声明提交一致 |
| `3daa60e8` | 记录 p9s 长窗口三个完整窗口完成与 `report.json` 极值定位 |
| `83d1d6c0` | 7 个透明探针 + Decimal50 回放：解释 p9s 长窗口极值（含工具归档与证据 JSON） |
| `894720da` | 记录最新组合源码验收 FAIL、触发项定位、既有缺陷分类、frcsat 隔离、甲烷只读调查证据 |

无生产源码改动；只 `docs/` 与 `tests/`（本提交串）变化。

## 4. 未完成 / 未验证 / 延期项（保留，不掩盖）

1. **最新组合源码验收：原 FAIL 的根因已修复，完整验收重跑尚未终态**。修复 `9f7c9a52`
   把无载体干池同位素质量隔离到既有 `trc_surface_residue`（不放松阈值、不裁 δ），组件级验证通过；
   但 **combined10 全量门禁重跑 `codex-d48a500ecf994f2fb89a6e2856a1cc6a` 与修复-vs-`6d995438`
   回归对照 `codex-2608f8979c3449e8bb3cef74e9aba087` 仍在运行**——终态通过前，组合验收仍记为未完成。
   其它数值残余路径（土壤相变/雪层/蒸发边界）未修改；LULC_IGBP_PC、其它区域/scheme/强迫产品仍未覆盖。
2. **完整 CI 矩阵未到终态**：旧 run `40725517808` 在 2026-10-29 04:0x CST 仍是 16/23 success、0 失败、
   7 个 mpich in_progress；**本机 GitHub 网络在 04:08 后不可达（DNS 被解析到 198.18.0.33 代理假地址，
   git/curl 均 SSL_ERROR_SYSCALL）**，因此修复提交 `9f7c9a52` 的推送与新的完整 CI 触发被阻塞，
   恢复后必须补做：推送 → 触发新 run → 跟到 23/23 终态。旧 run 的 mpich 是否在 6 h 上限内成功同样以终态为准。
3. **甲烷两年全球与水稻**：已开始真实推进（用户已授权 7920l 作为测试平台）。已真实编译
   `MT+BGC+CROP+LULC_IGBP_PFT+vanGenuchten+GridRiverLakeFlow` 并通过 mkinidata（4m12s、峰值 28 GB）；
   3 天首次运行触发 fail-closed `methane inundation exceeds host soil water`（deficit/pore 最大 0.0509
   > 交付默认容差 0.05），已记录并按源码自身设计的合法容差（0.20）重跑 3 天（rice off / rice on /
   rice on+only_wetland），后续 30 天水稻 pilot 与两年全球 W1 运行仍待完成。**未完成前不得写甲烷验收通过，
   旧甲烷结果只作对照、不是正确基准。**
4. **其它未覆盖边界**（沿用交接文档）：BGC+LULCC 是上游不支持组合；空 MPI worker 8-rank 崩溃未修；
   SinglePoint host 同月缓存覆盖未修；schema-1 缓存可靠性/阈值与 forcing identity 边界；观测准确度与
   全球代表性未验证；8 组以外 MPI rank 矩阵未重跑。
5. **p9s 尾部诊断有限**：只覆盖 12 个选中极值，不等于全部超范围值已解释；未做观测比较或参数标定。

## 5. 建议的下一步（按优先级）

1. 等 `codex-d48a500ecf994f2fb89a6e2856a1cc6a`（组合验收重跑）与 `codex-2608f8979c3449e8bb3cef74e9aba087`
   （修复回归对照）终态：通过则更新组合验收结论；失败则保留原因并继续定位（不得放宽阈值）。
2. GitHub 网络恢复后：推送 `9f7c9a52` 及后续文档提交，触发**修复提交的**完整 CI 矩阵并跟到 23/23 终态；
   同时读取旧 run `40725517808` 的最终 mpich 结果。
3. 甲烷：读取 3 天三口径结果；通过后跑 30 天水稻 pilot；再按 W1 零传输窗口准备/提交两个年度全球作业
   （每作业约 17 h，注意 24 h 上限与内存标定）；参数口径（pH/lake/rice/only_wetland/irrigation）
   逐项记录，旧结果只作参照。
4. 正式 PR 拆分/合并决定仍由用户作出；当前仅个人 Draft 留档，不向上游交付。

## 6. 证据索引

| 内容 | 路径 |
| --- | --- |
| p9r 验收 | `docs/material-transport-evidence/20261028/p9r-acceptance.md`（+ report/job/logs/provenance） |
| p9s 报告与诊断 | `docs/material-transport-evidence/20261028/p9s-report.json`、`p9s-longwindow-diagnosis.md`、`p9s-probe-results.json`、`longwindow-probe-tools/` |
| 组合验收 | `docs/material-transport-evidence/20261028/combined-acceptance.md`、`combined-report.json`、`combined-EVIDENCE.md`、`combined-attempts.json` |
| 甲烷闭合 | `docs/material-transport-evidence/20261028/methane-closure/` |
| CI | run `40725517808`（fork；证据在 `.tmp/ct/ci10/`，原报告在 GitHub Actions） |
| 状态总表 | `docs/material-transport-status.md` |

## 7. 复现命令（只读复核）

```sh
cd /Users/zhongwangwei/Desktop/Github/CoLM-constituent
git log --oneline 10e9fbd2..HEAD
~/.local/bin/run-7920 status codex-4ba194bf50074f68aaa3d76178c49830   # p9r succeeded/0
~/.local/bin/run-7920 status codex-23354a2d716f45059da97e74b357b906   # p9s succeeded/0
~/.local/bin/run-7920 status codex-452fcc287377488daef69b1491144c1b   # p9s 探针 succeeded/0
~/.local/bin/run-7920 status codex-2434412927a94e38a44ad21db5643a38   # 组合验收 failed/1（验收 FAIL）
gh run view 40725517808 --repo zhongwangwei/CoLM202X                  # CI 矩阵
```

本地私有工件（不进 Git）：`/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/{p9r_accept,p9s_report,p9s_probe,combined10,ch4_closure10,ci10}/`。
