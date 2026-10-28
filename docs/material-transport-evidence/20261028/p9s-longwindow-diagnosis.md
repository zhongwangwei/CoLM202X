# p9s 长窗口同位素极值诊断：透明探针 + Decimal50 账本回放

状态记录：2026-10-29 DeepSeek 接续。本文只覆盖选中的真实极值记录，不把有限/无 NaN 当科学通过。

## 结论

- **透明性成立**：7 个目标月探针全部模型 END（步数/最终时间戳精确），每个探针从基线月 restart 重跑，重叠段 history 与基线**逐位相同**（每探针 3 文件 / 225 变量）。12/12 报告极值先由基线 history 原文件/变量/索引复现（`verify_extremes.py`），再用于探针。
- **22 个 target×isotope 网格**：21 个有效 δ，1 个为负 signed ratio（模型输出 FillValue，不裁 δ、不调参）；所有有效 δ 由真实 gross water 与有符号同位素质量重建，patch 面积加权重建误差 < 1e-9‰。
- **Decimal50 回放**：58,634 个 finite-pool 事件（1407 个 signed 负损失）；最大质量回放误差 2.449e-14 kg/m²；affine 斜率/截距重建误差 ≤ 2.470e-14 / 4.317e-17；BOOK 日账本残差 ≤ 1.820e-15（水）/ 7.001e-21（同位素质量）。
- **物理解释**：长期窗口最大 δ 尾部（AMZ 湿地 O18 −368.91/+214.60‰、HDO −997.50/+1708.59‰；NEC 夏季土壤 HDO −815.04‰、O18 −125.13‰；NEC 冬季土壤 O18 +119.10‰ 等）是"有符号同位素质量 / 真实 gross 水"的比值，gross 水最小可至 1.394e-02（面积加权求和，不是净 ET）；不是存水组成异常。其中 AMZ 2004-01-15 的 HDO 配对值为负 signed ratio（Decimal50 回放 −2764.25‰），模型按缺测输出——这是有符号净吸收被遮罩，不是正向出水浓度。
- **与既有证据一致**：AMZ 2004-01-02 土壤 O18 −116.788877 / HDO −754.851159 与旧三点报告（`longwindow-tools/legacy-three-point-report.md`）一致；本次把同一机制从 3 个短期点/8 天扩展到 3 个窗口/12 个真实极值。
- **边界**：只覆盖选中的极值记录，不等于全部超范围值已解释（例如 AMZ 湿地 O18 有 241 个点超出宽诊断范围，本次只诊断其中 min/max 代表点）；不验证观测准确度、物理校准或全球代表性；网格为 10×10 区域。

## 透明性与独立复现

| 探针（目标月） | 模型 END | 步数/最终时间 | history 文件/变量 | 重叠记录 | trace rank 文件 | bad 诊断 |
| --- | --- | --- | --- | --- | --- | --- |
| amz-wet-2004-01 | PASS | 720 / 2004-01-15-84600 | 3/225 | 15 | 4 | none |
| amz-wet-2004-02 | PASS | 912 / 2004-02-19-84600 | 3/225 | 19 | 4 | none |
| amz-wet-2004-03 | PASS | 480 / 2004-03-10-84600 | 3/225 | 10 | 3 | none |
| nec-sum-2004-05 | PASS | 384 / 2004-05-08-84600 | 3/225 | 8 | 3 | none |
| nec-sum-2004-07 | PASS | 144 / 2004-07-03-84600 | 3/225 | 3 | 5 | none |
| nec-win-2004-01 | PASS | 864 / 2004-01-18-84600 | 3/225 | 18 | 5 | none |
| nec-win-2004-02 | PASS | 480 / 2004-02-10-84600 | 3/225 | 10 | 2 | none |

主探针作业 `codex-452fcc287377488daef69b1491144c1b`（succeeded/0）。另有独立作业 `codex-1425a3a12f584098a4f59f2916c7ea22` 重跑同一 `nec-win-2004-01` 段，5 个 trace 文件 SHA256 与主作业完全相同（独立提交的确定性交叉验证）。首个探针尝试 `codex-ff49825f37bf45f18201f49006db53c1` 因 Lead 的极值复现脚本 netCDF4 API 错误在建模前失败（无模型计算），修正后重提；该失败保留记录。

## 极值账本（模型 δ vs Decimal50 重建 vs 审计报告值）

| target | 报告值‰ | 变量 | 模型 δ‰ | Decimal50 δ‰ | mask | gross water | signed mass | patch 重建误差‰ |
| --- | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: |
| amz-wet-wet-HDO-max-2004-03-10 |  | wetland_evap_H2_18O | 214.597079 | 214.597079 | valid | 5.637e-02 | 1.373e-04 | 2.417e-13 |
| amz-wet-wet-HDO-max-2004-03-10 | 1708.586078 | wetland_evap_HDO | 1708.586078 | 1708.586078 | valid | 5.637e-02 | 2.378e-05 | 1.597e-13 |
| amz-wet-wet-HDO-min-2004-02-19 |  | wetland_evap_H2_18O | -133.582509 | -133.582509 | valid | 6.519e-02 | 1.133e-04 | 4.213e-14 |
| amz-wet-wet-HDO-min-2004-02-19 | -997.502658 | wetland_evap_HDO | -997.502658 | -997.502658 | valid | 6.519e-02 | 2.536e-08 | 2.498e-14 |
| amz-wet-soil-HDO-min-2004-01-02 |  | soilevap_H2_18O | -116.788877 | -116.788877 | valid | 8.166e+01 | 1.446e-01 | 6.504e-14 |
| amz-wet-soil-HDO-min-2004-01-02 | -754.851159 | soilevap_HDO | -754.851159 | -754.851159 | valid | 8.166e+01 | 3.118e-03 | 3.231e-14 |
| amz-wet-wet-H2_18O-min-2004-01-15 | -368.909653 | wetland_evap_H2_18O | -368.909653 | -368.909653 | valid | 6.200e-02 | 7.846e-05 | 7.035e-15 |
| amz-wet-wet-H2_18O-min-2004-01-15 |  | wetland_evap_HDO | FillValue | -2764.248717 | negative_signed_ratio | 6.200e-02 | -1.704e-05 |  |
| nec-sum-soil-HDO-min-2004-07-03 | -125.128595 | soilevap_H2_18O | -125.128595 | -125.128595 | valid | 1.394e-02 | 2.445e-05 | 3.325e-15 |
| nec-sum-soil-HDO-min-2004-07-03 | -815.040002 | soilevap_HDO | -815.040002 | -815.040002 | valid | 1.394e-02 | 4.015e-07 | 2.146e-14 |
| nec-sum-wet-HDO-max-2004-07-03 |  | wetland_evap_H2_18O | 73.218290 | 73.218290 | valid | 1.967e-01 | 4.233e-04 | 3.428e-14 |
| nec-sum-wet-HDO-max-2004-07-03 | 470.802937 | wetland_evap_HDO | 470.802937 | 470.802937 | valid | 1.967e-01 | 4.507e-05 | 4.104e-14 |
| nec-sum-wet-HDO-min-2004-05-08 |  | wetland_evap_H2_18O | -38.750575 | -38.750575 | valid | 1.956e+00 | 3.769e-03 | 7.355e-14 |
| nec-sum-wet-HDO-min-2004-05-08 | -279.245682 | wetland_evap_HDO | -279.245682 | -279.245682 | valid | 1.956e+00 | 2.195e-04 | 7.244e-14 |
| nec-win-soil-HDO-min-2004-02-10 |  | soilevap_H2_18O | -62.272686 | -62.272686 | valid | 6.089e-01 | 1.145e-03 | 9.324e-14 |
| nec-win-soil-HDO-min-2004-02-10 | -366.805730 | soilevap_HDO | -366.805730 | -366.805730 | valid | 6.089e-01 | 6.006e-05 | 4.641e-14 |
| nec-win-wet-HDO-min-2004-01-01 |  | wetland_evap_H2_18O | -25.341559 | -25.341559 | valid | 5.290e-01 | 1.034e-03 | 8.336e-14 |
| nec-win-wet-HDO-min-2004-01-01 | -154.992308 | wetland_evap_HDO | -154.992308 | -154.992308 | valid | 5.290e-01 | 6.962e-05 | 3.366e-14 |
| nec-win-soil-H2_18O-max-2004-01-18 | 119.097740 | soilevap_H2_18O | 119.097740 | 119.097740 | valid | 2.254e-01 | 5.058e-04 | 7.174e-14 |
| nec-win-soil-H2_18O-max-2004-01-18 |  | soilevap_HDO | -13.385477 | -13.385477 | valid | 2.254e-01 | 3.464e-05 | 3.757e-14 |
| nec-win-soil-HDO-max-2004-01-17 |  | soilevap_H2_18O | 29.091691 | 29.091691 | valid | 4.210e+01 | 8.687e-02 | 1.438e-13 |
| nec-win-soil-HDO-max-2004-01-17 | 84.463658 | soilevap_HDO | 84.463658 | 84.463658 | valid | 4.210e+01 | 7.111e-03 | 1.403e-13 |

## 每探针回放统计

| 探针 | finite-pool 事件 | signed 负损失 | 最大质量误差 kg/m² | 最大 affine 斜率误差 | BOOK 水/质量残差 | masks |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| amz-wet-2004-01 | 2354 | 304 | 1.935e-14 | 2.415e-14 | 1.265e-17 / 3.400e-22 | {'valid': 3, 'negative_signed_ratio': 1, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |
| amz-wet-2004-02 | 3472 | 169 | 2.255e-14 | 2.255e-14 | 9.202e-18 / 3.930e-22 | {'valid': 2, 'negative_signed_ratio': 0, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |
| amz-wet-2004-03 | 3134 | 33 | 2.449e-14 | 2.385e-14 | 9.677e-17 / 3.583e-22 | {'valid': 2, 'negative_signed_ratio': 0, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |
| nec-sum-2004-05 | 12550 | 75 | 3.162e-15 | 2.251e-14 | 1.820e-15 / 1.360e-21 | {'valid': 2, 'negative_signed_ratio': 0, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |
| nec-sum-2004-07 | 29430 | 731 | 3.082e-15 | 2.470e-14 | 6.050e-17 / 7.001e-21 | {'valid': 4, 'negative_signed_ratio': 0, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |
| nec-win-2004-01 | 7534 | 95 | 9.065e-15 | 1.442e-14 | 3.227e-17 / 1.000e-21 | {'valid': 6, 'negative_signed_ratio': 0, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |
| nec-win-2004-02 | 160 | 0 | 2.626e-15 | 8.079e-16 | 1.050e-17 / 4.900e-22 | {'valid': 2, 'negative_signed_ratio': 0, 'gross_water_nonpositive': 0, 'delta_sanity': 0, 'missing_patches': 0} |

## 方法与口径

1. 冻结源码 `c315c84973d0672d120b2b6926343a1b42b6ce21`（与 p9s 数值源码相同；`p9s/new` 与 `git archive c315c849` 619 文件逐字节一致），用私有 `instrument.py` 打点（只写 trace 文本，不改数值路径）。
2. 目标来自 p9s `report.json` 的真实 min/max 记录（文件/变量/时间索引/lat/lon），先用基线 history 复现后才用于探针。
3. 从 `case/restart` 的对应月 monthly checkpoint 起跑目标月（1–19 天）；`compare_subset.py` 比较 probe 与基线重叠 history 的完整文件集/变量集/形状/dtype/原始位与时间轴。
4. `replay_longwindow.py` 用 Decimal50 重算 Craig–Gordon affine 与 finite-pool 残余水公式、核对 BOOK 日账本与 GRID 面积加权重建；不做参数重标定。
5. 日账本区分 patch BOOK 累积与 finite-pool LIMIT 事件：部分湿润 patch 的蒸发并非全部走有限池分支，因此报告 LIMIT 覆盖率而不强求 100%；严格恒等式是 patch 累积 == BOOK 之和（残差见上表）。

## 复现

```sh
~/.local/bin/run-7920 status codex-452fcc287377488daef69b1491144c1b
# 本地工件（不进 Git）：/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/p9s_probe/
python3 run_replay_all.py fetched
```

## 证据

- `p9s-probe-results.json`（本目录）：作业终态、透明性计数、12 条极值复现、22 个 target×isotope 账本、回放汇总。
- `p9s-probe-job-status.json`：runner 终态与输入声明。
- `longwindow-probe-tools/`：探针与回放工具（instrument 打点、目标选择、subset 比较、极值复现、Decimal50 回放、汇总）。
- 远程：`codex-452fcc287377488daef69b1491144c1b/results/{report.json,logs/extreme-verify.json,trace/<probe>/}`（trace_rank_*.txt、compare.json、log.colm、case.nml 均保留）。
- 旧三点机制证据：`longwindow-tools/legacy-three-point-report.md`（本次独立交叉验证其 AMZ 土壤点数值）。
