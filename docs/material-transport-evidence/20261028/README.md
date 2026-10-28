# 交接用小型验收证据

这里只保留状态、比较摘要与必要错误摘录，不包含模型数据、NetCDF、构建产物或凭证。
完整结果在 7920 的 `/media/zhwei/data01/compute_runs/codex/<job-id>/results/`，
冻结驱动与源快照在相应 `work/`。任务号和读法见[交接文档](../../material-transport-handoff-deepseek-20261028.md)。

- [Vector / 实际 restart 独立复核](vector-acceptance.md)：限定配置的完整变量原始位比较。
- [p9r 无强迫同位素/CL 回归验收](p9r-acceptance.md)：旧 `28478b18` 对新 `c315c849` 的十例完整窗口，
  10/10 模型例与 10/10 严格比较通过，唯一允许差异为标量 schema 1→2；原始记录为
  `p9r-report.json`、`p9r-task-log.txt`、`p9r-job-status.json`、`p9r-source-versions.txt`。
- `vector-original-evidence.json` 忠实保留原记录；其中 partial `final_checkpoint`
  选中了 `const`，不能用原 13 文件 / 364 变量证明动态状态。
  `vector-dynamic-checkpoint.json` 与对应 job 状态补证真实 `2005-002-00000`
  的 18 文件 / 750 变量。原模型未重跑，旧记录未覆盖。
- `vector-io-results.json` 与 fatal 摘录：真实新旧计数 sidecar 以及损坏输入拒绝。
- `crossyear-cache-results.json` 与 `crossyear-cache-acceptance.md`：冻结 p9q 九模型的
  只读严格补核通过；原 runner 的冷启动误判和首轮审计失败仍保留，不假称旧任务成功。
- [p9s 长窗口同位素极值诊断](p9s-longwindow-diagnosis.md)：7 个透明探针（每月 restart 重跑目标月、
  重叠 history 逐位一致）+ 12 条真实极值复现 + Decimal50 账本回放；21 个有效 δ 由 gross water /
  有符号质量重建，1 个负 signed ratio 按缺测输出。原始记录为 `p9s-probe-results.json`、
  `p9s-probe-job-status.json`，工具在 `longwindow-probe-tools/`。
- [最新组合源码整体验收](combined-acceptance.md)：frcsat + 强迫同位素 + LULCC + CL + 泥沙/河网的完整组合
  **FAIL**（LULC_IGBP/NEC 在 step 4 触发已记账数值残余硬限 9.0882e-5 > 5.0e-5）；触发项定位为
  flood feedback + 运行期同位素强迫，且 `28478b18` 基线同样复现（既有交互缺陷）。frcsat 修复的数值不变性
  另有真实 120 步隔离证据。原始记录为 `combined-report.json`、`combined-job-status.json`、
  `combined-attempts.json`、`combined-EVIDENCE.md`。
- [组合验收缺陷的根因修复](residual-fix.md)：把“无载体干池质量隔离”从非挥发溶质推广到水同位素
  （`trc_surface_residue`，不放松阈值、不裁 δ）；组件级验证：两种构建触发案例均 END、零 TRC 诊断、
  仓库 9 项自检 PASS。完整组合验收重跑与修复-vs-基线回归对照的作业编号见该文；**在它们终态前不写成组合通过**。
  原始记录为 `residual-fix/validation.json`、`residual-fix/*-log-markers.txt`、`residual-fix/define-*.h`。
- [完整 hosted CI 矩阵记录](ci-matrix.md)：个人 fork `workflow_dispatch` 触发（run `40725517808`，head `10e9fbd2`，
  23 job）。记录时 16/23 success、0 失败、7 个 mpich 待终态；**未到终态，不写成完整通过**。
  `ci-coverage-check.txt` 是独立矩阵展开与 `check_defineh.py` 覆盖核对输出。
- [甲烷候选输入闭合与运行设计](methane-closure/README.md)：只读调查，结论 `NOT_READY`；15 个缺失边界文件
  权威源已实测存在于 `bd:` 路径（补传需授权），窗口改写 W1 可零新增数据（少 9 h）。仍未运行甲烷模型。
- `methane-input-audit-summary.json`：只读候选包元数据审计，结论为
  `NOT_READY / UNVERIFIED_CANDIDATE`。传输一致不代表来源唯一、输入闭合或科学正确；
  没有启动甲烷模型。完整原报告仍在审计任务结果目录。

原始位结论不比较 NetCDF 属性或文件编码；数值一致不等于观测准确性或全部配置验收通过。
