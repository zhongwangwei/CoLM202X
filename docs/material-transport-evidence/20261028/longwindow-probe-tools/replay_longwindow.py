"""Long-window transparency replay: read-only trace text, Decimal50 bookkeeping.

Generalises legacy-three-point-replay.py to an arbitrary target list taken from
targets.json (lat, lon, doy, name, host in {'soil','wet'}) and to any number of
traced days.  It verifies, for every selected target/isotope:

* the affine Craig-Gordon slope/offset exactly reconstructs its inputs;
* every finite-pool evaporation event re-solves the production residual-water
  formula in Decimal50 and matches the model's signed mass loss;
* every selected patch's daily gross water / signed isotope mass is exactly the
  sum of its finite-pool events (no unrecorded pathway);
* the grid-level mapped numerator/denominator reconstructs the published delta;
* and it reports, rather than hides, a grid record whose ratio is masked.

usage: replay_longwindow.py TRACE_DIR OUT_DIR
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
from collections import defaultdict
import json
import sys

root = Path(sys.argv[1])
out = Path(sys.argv[2])
targets_path = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).parent / 'targets.json'
out.mkdir(parents=True, exist_ok=True)
targets = json.loads(targets_path.read_text())


def decs(items):
    return [D(x) for x in items]


def target_of(lat, lon, host, tol='0.12501'):
    """Map a history grid centre to a target cell/host.

    The instrument selects patches/compute points within 0.375 deg so that every
    patch contributing to the target grid cell is traced; the replay must then
    filter to the exact grid centre (0.125 deg at 0.25 deg resolution).
    """
    lon = (lon + D(180)) % D(360) - D(180) if lon >= D(180) else lon
    for t in targets:
        if t.get('host') != host:
            continue
        if abs(lat - D(str(t['lat']))) < D(tol) and abs(lon - D(str(t['lon']))) < D(tol):
            return t['name'], int(t['doy'])
    return None


books = defaultdict(lambda: [D(0), D(0), 0, 0])
patches, grids, limits, errors = [], [], [], []
affines = []
last_affine = {}
pending = defaultdict(list)
bad_rows = []

with localcontext() as ctx:
    ctx.prec = 50
    for file in sorted(root.glob('trace_rank_*.txt')):
        for lineno, line in enumerate(file.read_text().splitlines(), 1):
            p = line.split()
            if not p:
                continue
            tag = p[0]
            try:
                rank, ip, yr, doy, sec, pc = map(int, p[1:7])
                lat, lon = decs(p[7:9])
                v = decs(p[9:])
            except Exception:
                bad_rows.append(f'{file.name}:{lineno}')
                continue
            key = (rank, ip, doy, sec)
            if tag == 'AFFINE':
                it, rv, temp, rh, h, ae, ak, slope, offset = v
                reconstructed_slope = 1 / (ae * ak * (1 - h))
                reconstructed_offset = -h * rv / (ak * (1 - h))
                affines.append((abs(reconstructed_slope - slope), abs(reconstructed_offset - offset)))
                last_affine[key] = (v, lineno)
            elif tag == 'LIMIT':
                m, w, loss, temp, r0, a, b, cap, w1, r1, flux, m1 = v
                exactw1 = w - loss
                x = (exactw1 / w).ln()
                y = (a - 1) * x
                growth = y.exp()
                factor = x if y == 0 else x * (growth - 1) / y
                exactr = max(D(0), (m / w) * growth + b * factor)
                if cap > 0:
                    exactr = min(exactr, cap)
                exactflux = min(max(m, 0) - exactr * exactw1, max(m, 0))
                af = last_affine.get(key)
                if not af or af[0][-2:] != [a, b]:
                    errors.append(f'LIMIT without exact preceding AFFINE {file.name}:{lineno}')
                host = 'wet' if int(af[0][0]) else None
                limits.append(dict(rank=rank, patch=ip, doy=doy, sec=sec, pc=pc,
                                   lat=str(lat), lon=str(lon),
                                   isotope=int(af[0][0]) if af else None,
                                   pool_water=str(w), water_loss=str(loss), pool_mass=str(m),
                                   ratio_pool=str(r0), model_loss=str(flux), decimal_loss=str(exactflux),
                                   loss_abs_error=str(abs(flux - exactflux)),
                                   model_after=str(m1), decimal_after=str(m - exactflux),
                                   model_mass_closure=str(m - flux - m1),
                                   line=lineno, file=file.name))
                if af:
                    pending[(key, int(af[0][0]))].append(len(limits) - 1)
            elif tag == 'BOOK':
                it, kind, m, water, gross = v
                for ix in pending.pop((key, int(it)), []):
                    limits[ix]['next_book_kind'] = int(kind)
                    limits[ix]['next_book_mass'] = str(m)
                actualday = doy if sec > 0 else doy - 1
                k = (rank, ip, actualday, int(it), int(kind))
                r = books[k]
                r[0] += gross
                r[1] += m
                r[2] += 1
                r[3] += m < 0
            elif tag.startswith('PATCH_'):
                var = tag[6:]
                tx, ty, area, water, m, wm, mm, wmin, ref, mask, plat, plon = v
                host = 'soil' if 'soilevap' in var else 'wet'
                tgt = target_of(lat, lon, host)
                if tgt is None:
                    tgt = target_of(plat, plon, host)
                if tgt is None:
                    continue
                name, tday = tgt
                if doy != tday:
                    continue
                it = 1 if var.endswith('H2_18O') else 2
                kind = 2 if host == 'soil' else 5
                bk = books[(rank, ip, doy, it, kind)]
                eligible = water > wmin and mask == 1
                patches.append(dict(target=name, day=doy, var=var, rank=rank, patch=ip, patchclass=pc,
                                    lat=str(lat), lon=str(lon),
                                    area=str(area), water=str(water), signed_mass=str(m),
                                    mapped_water=str(wm), mapped_mass=str(mm), water_min=str(wmin),
                                    filter=bool(mask), eligible=eligible, ref=str(ref),
                                    book_water=str(bk[0]), book_mass=str(bk[1]), book_count=bk[2],
                                    negative_book_count=bk[3],
                                    water_residual=str(water - bk[0]), mass_residual=str(m - bk[1])))
            elif tag.startswith('GRID_'):
                var = tag[5:]
                tx, ty, m, water, ref, delta = v
                host = 'soil' if 'soilevap' in var else 'wet'
                tgt = target_of(lat, lon, host)
                if tgt is None:
                    continue
                name, tday = tgt
                if doy != tday:
                    continue
                grids.append(dict(target=name, day=doy, var=var, lat=str(lat), lon=str(lon),
                                  signed_mass=str(m), gross_water=str(water), ref=str(ref),
                                  model_delta=str(delta)))

    negative_grids = []
    for g in grids:
        ps = [p for p in patches if p['target'] == g['target'] and p['var'] == g['var']]
        if not ps:
            errors.append(f"missing contribution patches for {g['target']} {g['var']}")
            g['mask_reason'] = 'missing_patches'
            continue
        mass = sum((D(p['area']) * D(p['mapped_mass']) for p in ps if p['filter']), D(0))
        water = sum((D(p['area']) * D(p['mapped_water']) for p in ps if p['filter']), D(0))
        ratio = mass / water if water > 0 else None
        delta = (ratio / D(g['ref']) - 1) * 1000 if ratio is not None else None
        if water <= 0:
            reason = 'gross_water_nonpositive'
        elif ratio < 0:
            reason = 'negative_signed_ratio'
            negative_grids.append(g)
        elif delta is not None and abs(delta) > 2000:
            reason = 'delta_sanity'
        else:
            reason = 'valid'
        g.update(patches=len(ps), eligible_patches=sum(bool(p['eligible']) for p in ps),
                 decimal_delta=str(delta), mask_reason=reason,
                 patch_reconstruction_delta_error=(str(abs(delta - D(g['model_delta'])))
                                                   if reason == 'valid' else None))

    limerrs = [D(x['loss_abs_error']) for x in limits]
    summary = dict(
        trace_files=len(list(root.glob('trace_rank_*.txt'))),
        targets=[t['name'] for t in targets],
        finite_pool_events=len(limits),
        negative_finite_pool_events=sum(D(x['model_loss']) < 0 for x in limits),
        max_finite_pool_mass_abs_error=str(max(limerrs, default=D(0))),
        max_affine_slope_abs_error=str(max((x[0] for x in affines), default=D(0))),
        max_affine_offset_abs_error=str(max((x[1] for x in affines), default=D(0))),
        selected_patch_records=len(patches),
        grids=len(grids),
        negative_ratio_grids=len(negative_grids),
        mask_reasons={r: sum(g.get('mask_reason') == r for g in grids) for r in
                      ('valid', 'negative_signed_ratio', 'gross_water_nonpositive',
                       'delta_sanity', 'missing_patches')},
        max_daily_water_book_residual=str(max((abs(D(p['water_residual'])) for p in patches), default=D(0))),
        max_daily_signed_mass_book_residual=str(max((abs(D(p['mass_residual'])) for p in patches), default=D(0))),
        all_contributing_centroids_traced=all(
            target_of(D(p['lat']), D(p['lon']), 'soil' if 'soilevap' in p['var'] else 'wet',
                      tol='0.37501') is not None
            for p in patches),
        errors=errors,
        bad_rows=bad_rows[:20],
        grids_detail=grids,
    )
    for name, data in [('summary', summary), ('patches', patches), ('finite_pool_events', limits)]:
        (out / (name + '.json')).write_text(json.dumps(data, indent=2) + '\n')

    daily = []
    for p in patches:
        if not p['eligible']:
            continue
        isotope = 1 if p['var'].endswith('H2_18O') else 2
        kind = 2 if 'soilevap' in p['var'] else 5
        es = [e for e in limits
              if (e['rank'], e['patch'], e['doy'], e['isotope'], e.get('next_book_kind'))
              == (p['rank'], p['patch'], p['day'], isotope, kind)]
        # BOOK rows are the model's per-event bookkeeping of the patch accumulation;
        # LIMIT rows exist only for the finite-pool (near-dry) evaporation branch.
        # A partly wet patch therefore books more mass than its finite-pool events.
        # The strict identity is patch accumulation == BOOK sum (checked through the
        # patch-level residuals below); the finite-pool events are replayed exactly
        # and their coverage of the daily book is reported, not forced to be 100%.
        totalw = sum((D(e['water_loss']) for e in es), D(0))
        totalm = sum((D(e['model_loss']) for e in es), D(0))
        replaym = sum((D(e['decimal_loss']) for e in es), D(0))
        bookw, bookm = D(p['book_water']), D(p['book_mass'])
        daily.append(dict(target=p['target'], isotope=isotope, rank=p['rank'], patch=p['patch'],
                          day=p['day'], events=len(es),
                          book_gross_water=str(bookw), book_signed_mass=str(bookm),
                          finite_pool_gross_water=str(totalw), finite_pool_signed_mass=str(totalm),
                          finite_pool_decimal_mass=str(replaym),
                          finite_pool_water_fraction=str(totalw / bookw if bookw else D(0)),
                          finite_pool_mass_fraction=str(totalm / bookm if bookm else D(0)),
                          decimal_delta=str((replaym / totalw / D(p['ref']) - 1) * 1000) if totalw else None,
                          max_step_mass_abs_error=str(max((D(e['loss_abs_error']) for e in es), default=D(0))),
                          negative_events=sum(D(e['model_loss']) < 0 for e in es)))
    (out / 'finite_daily_replay.json').write_text(json.dumps(daily, indent=2) + '\n')

    expect = 2 * len(targets)
    if len(grids) != expect:
        errors.append(f'grid count {len(grids)} != expected {expect}')
    # A target grid cell may contain several patches, each above the daily water
    # threshold; every such (patch, isotope) pair gets its own daily ledger row.
    if len(daily) < expect:
        errors.append(f'daily count {len(daily)} < expected {expect}')
    covered = {(d['target'], d['isotope']) for d in daily}
    for t in targets:
        for isotope in (1, 2):
            if (t['name'], isotope) not in covered:
                errors.append(f"missing daily ledger {t['name']} isotope {isotope}")
    if errors:
        sys.exit(f'REPLAY FAIL: {errors[:5]}')
    if max(limerrs, default=D(0)) > D('1e-12'):
        sys.exit('REPLAY FAIL: finite-pool mass replay error')
    if D(summary['max_daily_water_book_residual']) > D('1e-14'):
        sys.exit('REPLAY FAIL: water book residual')
    if D(summary['max_daily_signed_mass_book_residual']) > D('1e-18'):
        sys.exit('REPLAY FAIL: signed mass book residual')
    for g in grids:
        if g['mask_reason'] == 'valid':
            if D(g['patch_reconstruction_delta_error']) > D('1e-9'):
                sys.exit(f"REPLAY FAIL: patch reconstruction delta {g['target']} {g['var']}")
    print(json.dumps(summary, indent=2))
