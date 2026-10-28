#!/usr/bin/env python3
"""Freeze a methane parameter file from the delivery standard file.

Only the explicitly named keys are rewritten; everything else stays
byte-identical to run/standard_methane_parameter.nml from the frozen source,
so the parameter provenance stays auditable.

usage: mkparam.py SRC OUT RICE(true|false) ONLYWET(true|false) HISTVARS
"""
import sys

src, out, rice, onlywet, histvars = sys.argv[1:6]
host_tol = sys.argv[6] if len(sys.argv) > 6 and sys.argv[6] else None
mode_override = sys.argv[7] if len(sys.argv) > 7 and sys.argv[7] else None

with open(src) as fh:
    lines = fh.readlines()

want = {
    'DEF_METHANE%enable_rice_paddy': 'DEF_METHANE%enable_rice_paddy = .' + rice + '.',
    'DEF_METHANE%only_wetland': 'DEF_METHANE%only_wetland = .' + onlywet + '.',
    'DEF_METHANE%ch4_history_vars': "DEF_METHANE%ch4_history_vars = '" + histvars + "'",
}

seen = set()
new = []
mode_seen = False
for line in lines:
    stripped = line.strip()
    if mode_override is not None and stripped.startswith('DEF_METHANE%inundation_mode') and '=' in stripped:
        lead = line[: len(line) - len(line.lstrip())]
        new.append(lead + "DEF_METHANE%inundation_mode = '" + mode_override + "'" + '\n')
        mode_seen = True
        continue
    key = None
    for k in want:
        if stripped.startswith(k) and '=' in stripped:
            key = k
            break
    if key is not None:
        lead = line[: len(line) - len(line.lstrip())]
        new.append(lead + want[key] + '\n')
        seen.add(key)
    else:
        new.append(line)

if mode_override is not None and not mode_seen:
    sys.stderr.write('mkparam: inundation_mode line not found for override\n')
    sys.exit(1)
missing = set(want) - seen
if missing:
    sys.stderr.write('mkparam: missing keys %s\n' % sorted(missing))
    sys.exit(1)
if rice not in ('true', 'false') or onlywet not in ('true', 'false'):
    sys.stderr.write('mkparam: logical value must be true/false\n')
    sys.exit(1)

if host_tol is not None:
    text = ''.join(new)
    marker = '\n/\n'
    if 'DEF_METHANE%host_water_tolerance' in text:
        raise SystemExit('mkparam: host_water_tolerance already present in source file')
    idx = text.rindex(marker)
    add = ('    ! methane10 pilot override: documented host-water consistency tolerance\n'
           '    ! (delivery default is 0.05; global 2deg cold start exceeds it, see EVIDENCE.md)\n'
           '    DEF_METHANE%host_water_tolerance = ' + host_tol + '\n')
    text = text[:idx] + '\n' + add + text[idx + 1:]
    new = [text]

with open(out, 'w') as fh:
    fh.writelines(new)
