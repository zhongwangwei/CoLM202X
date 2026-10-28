#!/usr/bin/env python
"""Generate the p9s probe evidence JSON and diagnosis markdown from local artifacts."""
import datetime
import json
from pathlib import Path

root = Path(__file__).parent
repo = Path('/Users/zhongwangwei/Desktop/Github/CoLM-constituent')
evid = repo / 'docs/material-transport-evidence/20261028'

job = json.loads((root / 'job-status.json').read_text())
job_id = job['job_id']
report = json.loads((root / 'job-report.json').read_text())
xver = json.loads((root / 'extreme-verify.json').read_text())
agg = json.loads((root / 'replay-results/summary.json').read_text())


def fmax(values):
    return float(max(float(v) for v in values))


def fmin(values):
    return float(min(float(v) for v in values))


probes = list(report['probes'])
mask = {}
for p in agg['probes'].values():
    for k, v in p['summary']['mask_reasons'].items():
        mask[k] = mask.get(k, 0) + v

out = {
    'job': {
        'id': job_id, 'state': job.get('state'), 'exit_code': job.get('exit_code'),
        'created_at': datetime.datetime.fromtimestamp(job['created_at']).isoformat() if job.get('created_at') else None,
        'finished_at': datetime.datetime.fromtimestamp(job['finished_at']).isoformat() if job.get('finished_at') else None,
        'resources': job.get('resources'), 'inputs': job.get('inputs_remote'),
    },
    'source': {'commit': 'c315c84973d0672d120b2b6926343a1b42b6ce21',
               'note': 'p9s/new verified equal to git archive c315c849 (619 files, diff empty)'},
    'secondary_job': {'id': 'codex-1425a3a12f584098a4f59f2916c7ea22',
                      'purpose': 'independent rerun of nec-win-2004-01; trace_rank_*.txt SHA256 identical to the main job'},
    'failed_attempt': {'id': 'codex-ff49825f37bf45f18201f49006db53c1',
                       'reason': 'Lead extreme-reproduction script used a netCDF4 Variable.astype call that does not exist; '
                                 'the job failed before any model run. Fixed and resubmitted; kept for the record.'},
    'transparency': {n: {'end': report['probes'][n]['end'], 'steps': report['probes'][n]['steps'],
                         'last': report['probes'][n]['last'],
                         'trace_rank_files': report['probes'][n]['trace_files'],
                         'history_files': report['probes'][n]['transparency']['files'],
                         'history_variables': report['probes'][n]['transparency']['variables'],
                         'records': report['probes'][n]['transparency']['records'],
                         'bad_diagnostics': report['probes'][n]['bad']}
                     for n in probes},
    'extremes_reproduced': {'checked': xver['checked'], 'failures': xver['failures'], 'detail': xver['detail']},
    'replay_totals': {
        'probes': len(agg['probes']),
        'unique_targets': len(agg['expected']),
        'target_isotope_grids': sum(len(v) for v in agg['targets'].values()),
        'finite_pool_events': sum(p['summary']['finite_pool_events'] for p in agg['probes'].values()),
        'negative_finite_pool_events': sum(p['summary']['negative_finite_pool_events'] for p in agg['probes'].values()),
        'max_finite_pool_mass_abs_error': str(fmax(p['summary']['max_finite_pool_mass_abs_error'] for p in agg['probes'].values())),
        'max_affine_slope_abs_error': str(fmax(p['summary']['max_affine_slope_abs_error'] for p in agg['probes'].values())),
        'max_affine_offset_abs_error': str(fmax(p['summary']['max_affine_offset_abs_error'] for p in agg['probes'].values())),
        'max_daily_water_book_residual': str(fmax(p['summary']['max_daily_water_book_residual'] for p in agg['probes'].values())),
        'max_daily_signed_mass_book_residual': str(fmax(p['summary']['max_daily_signed_mass_book_residual'] for p in agg['probes'].values())),
        'mask_reasons': mask,
    },
    'targets': [],
    'valid_target_gross_water_min': None,
    'notes': [
        'history files of every probe are bit-identical to the overlapping baseline records (225 variables per probe).',
        'every reported extreme was reproduced from the baseline history file/variable/index before probing.',
        'model delta, Decimal50 replay delta and the reported audit extreme agree for every valid grid.',
        'negative signed ratios are reported as filled/missing, never clipped.',
        'the daily ledger separates patch BOOK accumulation from finite-pool LIMIT events; '
        'LIMIT coverage of the daily book is reported rather than forced to 100% '
        '(a partly wet patch evaporates through the non-finite-pool branch too).',
    ],
}
valid_w = []
for name, exp in agg['expected'].items():
    for var, rec in sorted(agg['targets'][name].items()):
        out['targets'].append({
            'target': name, 'variable': var,
            'reported_expected': exp['expected'].get('H2_18O' if var.endswith('H2_18O') else 'HDO'),
            'model_delta': rec['model_delta'], 'decimal50_delta': rec['decimal_delta'],
            'mask_reason': rec['mask_reason'], 'gross_water': rec['gross_water'],
            'signed_mass': rec['signed_mass'],
            'patch_reconstruction_error': rec['patch_reconstruction_delta_error'],
            'patches': rec['patches'], 'eligible_patches': rec['eligible_patches'],
            'records': [{k: r[k] for k in ('name', 'species', 'expected', 'var')} for r in exp['records']],
        })
        if rec['mask_reason'] == 'valid':
            valid_w.append(float(rec['gross_water']))
out['valid_target_gross_water_min'] = str(fmin(valid_w))
(root / 'p9s-probe-results.json').write_text(json.dumps(out, indent=2) + '\n')
(evid / 'p9s-probe-results.json').write_text(json.dumps(out, indent=2) + '\n')

# ---------------- diagnosis markdown ----------------
t = out['replay_totals']
L = []
L.append('# p9s 长窗口同位素极值诊断：透明探针 + Decimal50 账本回放\n')
L.append('状态记录：2026-10-29 DeepSeek 接续。本文只覆盖选中的真实极值记录，不把有限/无 NaN 当科学通过。\n')
L.append('## 结论\n')
L.append('- **透明性成立**：7 个目标月探针全部模型 END（步数/最终时间戳精确），每个探针从基线月 restart 重跑，'
         '重叠段 history 与基线**逐位相同**（每探针 3 文件 / 225 变量）。12/12 报告极值先由基线 history 原文件/变量/索引复现'
         '（`verify_extremes.py`），再用于探针。')
L.append(f"- **22 个 target×isotope 网格**：21 个有效 δ，1 个为负 signed ratio（模型输出 FillValue，不裁 δ、不调参）；"
         f"所有有效 δ 由真实 gross water 与有符号同位素质量重建，patch 面积加权重建误差 < 1e-9‰。")
L.append(f"- **Decimal50 回放**：{t['finite_pool_events']:,} 个 finite-pool 事件（{t['negative_finite_pool_events']} 个 signed 负损失）；"
         f"最大质量回放误差 {float(t['max_finite_pool_mass_abs_error']):.3e} kg/m²；affine 斜率/截距重建误差 "
         f"≤ {float(t['max_affine_slope_abs_error']):.3e} / {float(t['max_affine_offset_abs_error']):.3e}；"
         f"BOOK 日账本残差 ≤ {float(t['max_daily_water_book_residual']):.3e}（水）/ "
         f"{float(t['max_daily_signed_mass_book_residual']):.3e}（同位素质量）。")
L.append('- **物理解释**：长期窗口最大 δ 尾部（AMZ 湿地 O18 −368.91/+214.60‰、HDO −997.50/+1708.59‰；'
         'NEC 夏季土壤 HDO −815.04‰、O18 −125.13‰；NEC 冬季土壤 O18 +119.10‰ 等）是"有符号同位素质量 / 真实 gross 水"的比值，'
         f"gross 水最小可至 {float(out['valid_target_gross_water_min']):.3e}（面积加权求和，不是净 ET）；不是存水组成异常。"
         '其中 AMZ 2004-01-15 的 HDO 配对值为负 signed ratio（Decimal50 回放 −2764.25‰），模型按缺测输出——'
         '这是有符号净吸收被遮罩，不是正向出水浓度。')
L.append('- **与既有证据一致**：AMZ 2004-01-02 土壤 O18 −116.788877 / HDO −754.851159 与旧三点报告'
         '（`longwindow-tools/legacy-three-point-report.md`）一致；本次把同一机制从 3 个短期点/8 天扩展到 3 个窗口/12 个真实极值。')
L.append('- **边界**：只覆盖选中的极值记录，不等于全部超范围值已解释（例如 AMZ 湿地 O18 有 241 个点超出宽诊断范围，'
         '本次只诊断其中 min/max 代表点）；不验证观测准确度、物理校准或全球代表性；网格为 10×10 区域。\n')

L.append('## 透明性与独立复现\n')
L.append('| 探针（目标月） | 模型 END | 步数/最终时间 | history 文件/变量 | 重叠记录 | trace rank 文件 | bad 诊断 |')
L.append('| --- | --- | --- | --- | --- | --- | --- |')
for p in probes:
    d = out['transparency'][p]
    rec = sum(d['records'].values()) // d['history_files'] if d['history_files'] else 0
    L.append(f"| {p} | {'PASS' if d['end'] else 'FAIL'} | {d['steps']} / {d['last'][1] if d['last'] else '-'} | "
             f"{d['history_files']}/{d['history_variables']} | {rec} | {d['trace_rank_files']} | "
             f"{', '.join(d['bad_diagnostics']) if d['bad_diagnostics'] else 'none'} |")
L.append('')
L.append('主探针作业 `codex-452fcc287377488daef69b1491144c1b`（succeeded/0）。另有独立作业 '
         '`codex-1425a3a12f584098a4f59f2916c7ea22` 重跑同一 `nec-win-2004-01` 段，5 个 trace 文件 SHA256 与主作业完全相同'
         '（独立提交的确定性交叉验证）。首个探针尝试 `codex-ff49825f37bf45f18201f49006db53c1` 因 Lead 的极值复现脚本 '
         'netCDF4 API 错误在建模前失败（无模型计算），修正后重提；该失败保留记录。\n')

L.append('## 极值账本（模型 δ vs Decimal50 重建 vs 审计报告值）\n')
L.append('| target | 报告值‰ | 变量 | 模型 δ‰ | Decimal50 δ‰ | mask | gross water | signed mass | patch 重建误差‰ |')
L.append('| --- | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: |')
for row in out['targets']:
    exp = row['reported_expected']
    exp_s = f"{float(exp):.6f}" if exp not in (None, '') else ''
    md = float(row['model_delta'])
    md_s = f"{md:.6f}" if abs(md) < 1e6 else 'FillValue'
    err = row['patch_reconstruction_error']
    err_s = f"{float(err):.3e}" if err not in (None, '', 'None') else ''
    L.append(f"| {row['target']} | {exp_s} | {row['variable'].split('f_mt_delta_')[-1]} | {md_s} | "
             f"{float(row['decimal50_delta']):.6f} | {row['mask_reason']} | {float(row['gross_water']):.3e} | "
             f"{float(row['signed_mass']):.3e} | {err_s} |")
L.append('')
L.append('## 每探针回放统计\n')
L.append('| 探针 | finite-pool 事件 | signed 负损失 | 最大质量误差 kg/m² | 最大 affine 斜率误差 | BOOK 水/质量残差 | masks |')
L.append('| --- | ---: | ---: | ---: | ---: | --- | --- |')
for name in probes:
    s = agg['probes'][name]['summary']
    L.append(f"| {name} | {s['finite_pool_events']} | {s['negative_finite_pool_events']} | "
             f"{float(s['max_finite_pool_mass_abs_error']):.3e} | {float(s['max_affine_slope_abs_error']):.3e} | "
             f"{float(s['max_daily_water_book_residual']):.3e} / {float(s['max_daily_signed_mass_book_residual']):.3e} | "
             f"{s['mask_reasons']} |")
L.append('')
L.append('## 方法与口径\n')
L.append('1. 冻结源码 `c315c84973d0672d120b2b6926343a1b42b6ce21`（与 p9s 数值源码相同；`p9s/new` 与 `git archive c315c849` '
         '619 文件逐字节一致），用私有 `instrument.py` 打点（只写 trace 文本，不改数值路径）。')
L.append('2. 目标来自 p9s `report.json` 的真实 min/max 记录（文件/变量/时间索引/lat/lon），先用基线 history 复现后才用于探针。')
L.append('3. 从 `case/restart` 的对应月 monthly checkpoint 起跑目标月（1–19 天）；`compare_subset.py` 比较 probe 与基线重叠 '
         'history 的完整文件集/变量集/形状/dtype/原始位与时间轴。')
L.append('4. `replay_longwindow.py` 用 Decimal50 重算 Craig–Gordon affine 与 finite-pool 残余水公式、核对 BOOK 日账本与 GRID '
         '面积加权重建；不做参数重标定。')
L.append('5. 日账本区分 patch BOOK 累积与 finite-pool LIMIT 事件：部分湿润 patch 的蒸发并非全部走有限池分支，因此报告 LIMIT '
         '覆盖率而不强求 100%；严格恒等式是 patch 累积 == BOOK 之和（残差见上表）。\n')
L.append('## 复现\n')
L.append('```sh')
L.append('~/.local/bin/run-7920 status codex-452fcc287377488daef69b1491144c1b')
L.append('# 本地工件（不进 Git）：/Users/zhongwangwei/Desktop/Github/CoLM202X/.tmp/ct/p9s_probe/')
L.append('python3 run_replay_all.py fetched')
L.append('```')
L.append('')
L.append('## 证据\n')
L.append('- `p9s-probe-results.json`（本目录）：作业终态、透明性计数、12 条极值复现、22 个 target×isotope 账本、回放汇总。')
L.append('- `p9s-probe-job-status.json`：runner 终态与输入声明。')
L.append('- `longwindow-probe-tools/`：探针与回放工具（instrument 打点、目标选择、subset 比较、极值复现、Decimal50 回放、汇总）。')
L.append('- 远程：`codex-452fcc287377488daef69b1491144c1b/results/{report.json,logs/extreme-verify.json,trace/<probe>/}`'
         '（trace_rank_*.txt、compare.json、log.colm、case.nml 均保留）。')
L.append('- 旧三点机制证据：`longwindow-tools/legacy-three-point-report.md`（本次独立交叉验证其 AMZ 土壤点数值）。')
(evid / 'p9s-longwindow-diagnosis.md').write_text('\n'.join(L) + '\n')
print('wrote evidence JSON and diagnosis markdown')
print('maxima:', {k: t[k] for k in t if k.startswith('max')})
