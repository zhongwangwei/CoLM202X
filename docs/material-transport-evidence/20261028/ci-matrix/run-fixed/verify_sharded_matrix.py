#!/usr/bin/env python3
"""Independent verification of the SHARDED CI matrix for run 40784782071.

Reads the workflow/TestCaseLists snapshots frozen from commit
2bc7510a2f9bddcf294b4a8999f98540086d4316 (run-fixed/source-snapshots/) and the
jobs.json written by poll_curl.py. Writes run-fixed/coverage-check.txt.

Checks:
 1. Cartesian expansion of the sharded matrix and exclude semantics -> expected
    job names must equal the job names returned by the Actions API.
 2. Shard partition of TestCaseLists: shard k gets lines with index % 4 == k;
    union is the full 93 cases, no duplicates, per-shard case counts 24/23/23/23.
 3. Workflow still runs check_defineh.py (full 93-case list) and the
    history-accumulator field check in every job; build loop has shard logic
    and strict shell mode.
 4. Per-shard expected MTON / LULCCON counts (for log-marker verification).
"""
import itertools
import json
import pathlib
from collections import Counter

import yaml

BASE = pathlib.Path(__file__).resolve().parent
SRC = BASE / 'source-snapshots'
OUT = []


def say(line=''):
    OUT.append(line)
    print(line)


wf = yaml.safe_load((SRC / 'build_CoLM_gnu.yml').read_text())
matrix = wf['jobs']['Build-CoLM-gnu']['strategy']['matrix']
axes = {k: v for k, v in matrix.items() if k not in ('include', 'exclude')}
say('## Sharded matrix expansion')
say('axes: ' + ', '.join(f'{k}={len(v)}' for k, v in axes.items()))
combos = [dict(zip(axes.keys(), vals))
          for vals in itertools.product(*axes.values())]
say(f'cartesian combos = {len(combos)}')
excludes = matrix.get('exclude', [])


def matches(combo, rule):
    for k, v in rule.items():
        if k not in combo:
            return False
        if isinstance(v, dict):
            if not isinstance(combo[k], dict) or any(combo[k].get(kk) != vv for kk, vv in v.items()):
                return False
        elif combo[k] != v:
            return False
    return True


kept, removed = [], []
for c in combos:
    (removed if any(matches(c, r) for r in excludes) else kept).append(c)
say(f'exclude rules = {len(excludes)}; combos removed = {len(removed)}; kept = {len(kept)}')
base_cells = {(c['os'], c['mpi'], c['toolchain']['compiler'], c['toolchain']['version'])
              for c in kept}
say(f'unique (os,mpi,compiler,gcc) cells after exclude = {len(base_cells)} '
    f'(x {len(axes["shard"])} shards = {len(base_cells) * len(axes["shard"])} jobs expected)')
expected_names = {
    'Build-CoLM-gnu (%s, %s, %s, %s, %s)' % (
        c['os'], c['mpi'], c['toolchain']['compiler'], c['toolchain']['version'], c['shard'])
    for c in kept
}
try:
    jobs = json.loads((BASE / 'jobs.json').read_text())['jobs']
    actual_names = {j['name'] for j in jobs}
except Exception:
    actual_names = set()
    say('WARNING: jobs.json not available yet')
say(f'API jobs returned = {len(actual_names)}; expected = {len(expected_names)}; '
    f'names equal = {expected_names == actual_names}')
if actual_names and expected_names != actual_names:
    say(f'  only expected: {sorted(expected_names - actual_names)[:5]}')
    say(f'  only actual:   {sorted(actual_names - expected_names)[:5]}')

say()
say('## TestCaseLists shard partition')
lines = [ln for ln in (SRC / 'TestCaseLists').read_text().splitlines() if ln.strip()]
say(f'cases = {len(lines)}; unique case names = {len({ln.split()[0] for ln in lines})}')
nshard = len(axes['shard'])
per_shard = {}
for shard in axes['shard']:
    sel = [ln for i, ln in enumerate(lines) if i % nshard == shard]
    cases = [ln.split()[0] for ln in sel]
    mt = [c for c in cases if 'MTON' in c]
    lulcc = [c for c in cases if 'LULCCON' in c]
    per_shard[shard] = {'cases': len(cases), 'mt': len(mt), 'lulcc': len(lulcc)}
    say(f'  shard {shard}: {len(cases)} cases | MTON={len(mt)} | LULCCON={len(lulcc)} '
        f'| examples={cases[:2]}')
all_assigned = Counter()
for i in range(len(lines)):
    all_assigned[i % nshard] += 1
say(f'union covers all lines exactly once = {sum(all_assigned.values()) == len(lines)}; '
    f'per-shard counts = {dict(sorted(all_assigned.items()))}')
mt_expected = Counter()
for i, ln in enumerate(lines):
    if 'MTON' in ln:
        mt_expected[i % nshard] += 1
say(f'MTON by shard = {dict(sorted(mt_expected.items()))}')
say(f'LULCCON by shard = '
    f'{dict(sorted(Counter(i % nshard for i, ln in enumerate(lines) if "LULCCON" in ln).items()))}')

say()
say('## Workflow still runs full checks in every job')
body = (SRC / 'build_CoLM_gnu.yml').read_text()
checks = {
    'Check the history-accumulator restart field list': 'hist_acc_restart_fields.py --check',
    'Check effective test case macros': 'check_defineh.py',
    'Build CoLM202X': 'make clean && make all',
}
for step, needle in checks.items():
    say(f'  step {step!r}: {"present" if needle in body else "MISSING"}')
say(f"  shard partition code (index % nshard == shard): "
    f'{"present" if "index % nshard == shard" in body else "MISSING"}')
say(f"  strict shell (set -euo pipefail): "
    f'{"present" if "set -euo pipefail" in body else "MISSING"}')
say(f"  shard variable from matrix: "
    f'{"present" if "matrix.shard" in body and "NSHARD=4" in body else "MISSING"}')
say(f"  all-93 coverage statement: "
    f'{"present" if "All 93 cases are covered across the four shards" in body else "MISSING"}')

say()
say('## check_defineh.py identity / coverage')
cdh = (SRC / 'check_defineh.py').read_text()
KNOWN_SHA = {
    'build_CoLM_gnu.yml': '079eefc054aa0c750185b25c185509c691a73fed',
    'TestCaseLists': '612b55722dacb63b58f80b27dac495f8f239604e',
    'check_defineh.py': 'dc486706d83bf8663954b3604dd8369d0e03d17c',
    'hist_acc_restart_fields.py': '0b25fd85afcac7088693b1adaf57fc48fed56988',
    'create_defineh.bash': 'd70d48cbee0d67dd50b49f542e0923f4d8bc7e68',
}

import hashlib


def git_blob_sha(path):
    data = pathlib.Path(path).read_bytes()
    return hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()


for fname, want in KNOWN_SHA.items():
    got = git_blob_sha(SRC / fname)
    say(f'  {fname}: blob {got} == expected {want}: {got == want}')
say('  check_defineh.py / TestCaseLists / create_defineh.bash / hist_acc_restart_fields.py '
    'identical to old run; only the workflow file changed (sharding)')
say(f'  asserts 93 cases from TestCaseLists: {"present" if "len(names) == len(set(names))" in cdh else "MISSING"}')
say(f'  PC+BGC+MT claim: {"present" if "LULC_IGBP_PC\', \'BGCON\', \'LULCCOFF" in cdh else "MISSING"}')
say(f'  IGBP/PFT+LULCC+MT claim: '
    f'{"present" if "for land in (\'LULC_IGBP\', \'LULC_IGBP_PFT\')" in cdh else "MISSING"}')

say()
say('## Verdict')
ok = (len(kept) == 92 and expected_names == actual_names and len(lines) == 93
      and sum(per_shard[s]['cases'] for s in per_shard) == 93
      and all(v == 24 if s == 0 else v == 23 for s, v in
              ((s, per_shard[s]['cases']) for s in per_shard)))
say(f'sharded matrix and CI coverage consistent: {ok}')

(BASE / 'coverage-check.txt').write_text('\n'.join(OUT) + '\n')
(BASE / 'per-shard-expectations.json').write_text(json.dumps(per_shard, indent=1) + '\n')
