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
- `methane-input-audit-summary.json`：只读候选包元数据审计，结论为
  `NOT_READY / UNVERIFIED_CANDIDATE`。传输一致不代表来源唯一、输入闭合或科学正确；
  没有启动甲烷模型。完整原报告仍在审计任务结果目录。

原始位结论不比较 NetCDF 属性或文件编码；数值一致不等于观测准确性或全部配置验收通过。
