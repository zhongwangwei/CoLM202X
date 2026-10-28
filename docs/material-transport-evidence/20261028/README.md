# 交接用小型验收证据

这里只保留状态、比较摘要与必要错误摘录，不包含模型数据、NetCDF、构建产物或凭证。
完整结果在 7920 的 `/media/zhwei/data01/compute_runs/codex/<job-id>/results/`，
冻结驱动与源快照在相应 `work/`。任务号和读法见[交接文档](../../material-transport-handoff-deepseek-20261028.md)。

- [Vector / 实际 restart 独立复核](vector-acceptance.md)：限定配置的完整变量原始位比较。
- `vector-original-evidence.json` 忠实保留原记录；其中 partial `final_checkpoint`
  选中了 `const`，不能用原 13 文件 / 364 变量证明动态状态。
  `vector-dynamic-checkpoint.json` 与对应 job 状态补证真实 `2005-002-00000`
  的 18 文件 / 750 变量。原模型未重跑，旧记录未覆盖。
- `vector-io-results.json` 与 fatal 摘录：真实新旧计数 sidecar 以及损坏输入拒绝。
- `crossyear-cache-results.json` 与 `crossyear-cache-acceptance.md`：冻结 p9q 九模型的
  只读严格补核通过；原 runner 的冷启动误判和首轮审计失败仍保留，不假称旧任务成功。
- `methane-input-audit-summary.json`：只读候选包元数据审计，结论为
  `NOT_READY / UNVERIFIED_CANDIDATE`。传输一致不代表来源唯一、输入闭合或科学正确；
  没有启动甲烷模型。完整原报告仍在审计任务结果目录。

原始位结论不比较 NetCDF 属性或文件编码；数值一致不等于观测准确性或全部配置验收通过。
