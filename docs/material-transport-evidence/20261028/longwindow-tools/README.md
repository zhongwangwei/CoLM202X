# 尚未完成适配的长窗口诊断工具

保留这些小脚本，避免 DeepSeek 换机器后丢失仅存在于 Mac `.tmp` 的准备工作。
它们不是生产代码，也不属于已经通过的全窗口测试。

- `instrument.py`：已参数化的 Fortran 稀疏记账补丁器；需要 `targets.json` 与冻结
  `new/` 源码。**先复制到新私有探针目录，再阅读其路径假设后执行**，不可直接修改
  当前工作树或正在运行的源快照。
- `legacy-three-point-replay.py`：旧 AMZ 三个短期点的 Decimal50 回放，目标、日期、
  六项账本断言仍硬编码；需要接手者按真实长期极值适配，不能直接宣称长期验证完成。
- `legacy-three-point-report.md`：已验证的旧三点范围与误差，不推广到全部长期极值。
- `longwindow-handoff.md`：当前 p9s 任务、准确结果路径、状态快照与接手步骤。
- `original-plan.md`：原假设与诊断计划，仍是待证假设，不是科学结论。

原短期 runner、日志和账本位于
`/media/zhwei/data01/compute_runs/codex/codex-4f7bfe583e234dfe8ac35115d15fee56/`。
新的大计算必须使用 `run-7920` 私有快照与显式只读输入，不在本机或旧结果目录运行。
