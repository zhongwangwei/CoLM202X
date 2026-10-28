# p9s 长窗口同位素诊断交接

## 当前状态（不是验收完成）

2026-10-28 **23:46:53 CST / 15:46:53 UTC** 实测：作业
`codex-23354a2d716f45059da97e74b357b906` 为 `running`。三个算例尚无 END；截至各自当前步，TRC_SIG、TRC_BAL 和 residual hard failure 计数均为 0。

| 算例 / 完整目标窗口 | 当前最后一步与日期 |
|---|---|
| NEC 冬季，2004-01-01..02-29（60 日） | 1992，2004-02-11-41400 |
| NEC 夏季，2004-05-01..07-31（92 日） | 1984，2004-06-11-27000 |
| AMZ 湿季，2004-01-01..03-31（91 日） | 2010，2004-02-11-73800 |

启动时间为 2026-10-28 13:42:14 UTC；配额 16 CPU / 32 GB，timeout 21600 秒。**保持现有健康作业，不因观察超时重新提交。** 本目录 `status-handoff.json`、`progress-handoff.txt` 为这次原始状态证据。

## 源码与结果

- 数值源码冻结于 `c315c84973d0672d120b2b6926343a1b42b6ce21`。正确 NEC/AMZ 10×10 区域、IsoGSM 2000–2004、分馏启用、NBAD=0；不是旧 p9k 缺失区域边界的全球尺寸输出。
- 原私有测试工件：`/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/p9s/`。`run.sh`、`audit.py`、`source_sha.json`、`source_commit.txt`、`PLAN.md` 已存在。
- 远程结果根：`/media/zhwei/data01/compute_runs/codex/codex-23354a2d716f45059da97e74b357b906/results`。
- 模型日志/数据：上述根下 `runs/{nec-win,nec-sum,amz-wet}/`；monthly restart 在各 `case/restart/`。
- 全窗口完成后脚本输出的是 **`report.json` 和 `audit.log`**，不是 `audit.json`。`audit.py` 检查完整日数/步数、全部材料浮点有效值、区域大小并列出 δ/d-excess 极值。范围阈值只是诊断触发，不是物理通过线。

## 下一步，不重跑整段窗口

```sh
~/.local/bin/run-7920 status codex-23354a2d716f45059da97e74b357b906
~/.local/bin/run-7920 logs codex-23354a2d716f45059da97e74b357b906
```

1. 等原作业实际终态；读取 `report.json` 的 failures 和实际完整范围，不能只看退出码。
2. 用完整窗口的真实 min/max/显著尾部选择日期、网格与变量，区分储量/叶片诊断和有符号通量比。不要沿用旧 p9k 极值，或把此前三个短期点推广为全部长期尾部。
3. 从相应月份前的已存 restart 做最小透明探针，避免为晚期目标重跑 90 日。重计算仍必须走 `run-7920 --agent codex`，输入显式声明只读；模型/配置保持相同。
4. 核对探针对应时段所有 history 变量位一致，再检查真实 gross-water / signed-isotope 日账本和有限池高精度回放；不能使用 `f_fevpa` 净 ET 代替真实分母，不能裁 δ 或调参数制造通过。

## 已准备工具及未完成适配

- `p9s/trace/instrument.py` 已参数化 `targets.json`（lat、lon、doy），记录 AFFINE/LIMIT/BOOK/PATCH/GRID、真实 history 日期及未截断 `eair/esat`。它会修改同目录 `new/`；**只在新私有探针目录复制后执行，不能修改旧源码/运行作业**。
- 可复用 `iso_trace9/{run.sh,replay.py,REPORT.md}` 的透明性与 Decimal50 方法。但旧 runner 固定 AMZ Jan1–8，旧 replay 固定三个点、日序和六个目标结果；**目前尚未适配长期新目标、跨月份段与 raw-RH 对应关系，不能直接当作已可运行/已完成的探针**。
- `p9s/PLAN.md` 的 H1 仍是假设：强尾部可能是高湿度下有符号净重同位素交换，而非储水组成异常。需要上述实际输入与账本证明；本文未得出长期极值根因，也未验证观测准确度/物理校准。

本接续子任务未启动新模型/探针、未改生产、未提交或推送。已按用户要求停止新实现，留给 DeepSeek 接续。
