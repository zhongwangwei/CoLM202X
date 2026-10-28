# p9s long-window probe tools

这些脚本实现 2026-10-29 的 p9s 长窗口极值诊断（结果见[诊断记录](../p9s-longwindow-diagnosis.md)、
[p9s-probe-results.json](../p9s-probe-results.json)）。它们不是生产代码，也不属于已通过的全窗口测试。

- `instrument.py`：私有 Fortran 稀疏打点补丁器（与 `../longwindow-tools/instrument.py` 相同）。
  需要 `targets.json` 与冻结的 `new/` 源码；只写 trace 文本，不改数值路径。先复制到新私有探针
  目录再执行；会修改该目录下的 `new/`。
- `select_targets.py`：从 p9s `report.json` 选择真实 min/max 极值，去重相同 (cell,doy,host)，
  生成 `targets.json`、`probes.json`、`selection.json`（每条保留文件/记录索引/δ/lat/lon）。
- `run.sh`：探针作业。重建打点源码、`sha256`-free 构建、先在基线 history 上复现 12 条极值，
  再从基线 monthly restart 起跑目标月，逐位比较 probe 与基线重叠 history，并保存 trace 文件。
  需要只读输入 `P9S`（p9s results 根）、`ISOGSM`、`RAW`、`RT`、`JRA`、`UCAT`。
- `compare_subset.py`：probe 与基线重叠 history 的严格比较（文件集/变量集/形状/dtype/原始位/时间轴子集）。
- `verify_extremes.py`：用基线 history 的原文件/变量/索引复现审计报告的极值。
- `replay_longwindow.py`：Decimal50 回放（CG affine、finite-pool 残余水公式、BOOK 日账本、
  GRID 面积加权重建），目标从 `targets.json` 读取。
- `run_replay_all.py`：对每个已取回的 probe 目录执行回放并汇总。
- `gen_evidence.py`：从本地工件生成本目录证据 JSON 与诊断 markdown。

`targets.json` / `probes.json` / `selection.json` 是本次实际使用的 12 条极值记录与 7 个目标月定义，
保留以便审计。回放依赖 Python 标准库 Decimal（50 位）；NetCDF 比较依赖 netCDF4/NumPy。
