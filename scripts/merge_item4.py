#!/usr/bin/env python3
"""Merge item-4 outputs: frame-level paired cluster bootstrap + seed-level t.

Verdict per cell:
  robust            both intervals exclude zero and point the same way
  frame_only        only the paired cluster bootstrap excludes zero
  seed_only         only the seed-level t interval excludes zero
  not_established   neither excludes zero
  boot_cross / ...  (kept explicit in the CSV)
"""
import csv
from collections import defaultdict
from pathlib import Path

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')


def verdict(boot_cross, seed_cross, agree):
    if not boot_cross and not seed_cross:
        return 'robust' if agree else 'direction_conflict'
    if not boot_cross and seed_cross:
        return 'frame_only'
    if boot_cross and not seed_cross:
        return 'seed_only'
    return 'not_established'


def same_sign(lo, hi):
    return (lo > 0 and hi > 0) or (lo < 0 and hi < 0)


def main():
    boot = {(r['group'], r['comparison'], r['condition']): r
            for r in csv.DictReader((OUT / 'paired_diff_ci.csv').open())}
    seed = {(r['group'], r['comparison'], r['condition']): r
            for r in csv.DictReader((OUT / 'seed_level_diff.csv').open())}
    assert boot.keys() == seed.keys(), (len(boot), len(seed))

    rows = []
    for k in boot:
        b, s = boot[k], seed[k]
        blo, bhi = float(b['ci_lo_pp']), float(b['ci_hi_pp'])
        slo, shi = float(s['seed_t_lo_pp']), float(s['seed_t_hi_pp'])
        pt = float(b['point_diff_pp'])
        # A zero-crossing interval can still be one-sided if both bounds share a sign;
        # here cross flags are computed as lo<0<hi (strict), so not-cross == same sign.
        bc, sc = int(b['crosses_zero']), int(s['seed_t_cross'])
        agree = same_sign(blo, bhi) and same_sign(slo, shi) and (blo > 0) == (slo > 0)
        rows.append(dict(
            group=k[0], comparison=k[1], condition=k[2],
            arm_a=b['arm_a'], arm_b=b['arm_b'], protocol=b['protocol'],
            point_diff_pp=b['point_diff_pp'],
            boot_lo_pp=b['ci_lo_pp'], boot_hi_pp=b['ci_hi_pp'], boot_cross=bc,
            seed_t_lo_pp=s['seed_t_lo_pp'], seed_t_hi_pp=s['seed_t_hi_pp'], seed_t_cross=sc,
            seed_diffs_pp=f"{s['d_s0_pp']} {s['d_s1_pp']} {s['d_s2_pp']}",
            seed_sd_pp=s['sd_ddof1_pp'],
            n_positive_seeds=s['n_positive_seeds'],
            n_negative_seeds=s['n_negative_seeds'],
            verdict=verdict(bc, sc, agree),
            both_agree_direction=int(agree),
        ))
    order = {('FLIR_P26', 'DAT_vs_FMAug'): 0,
             ('M3FD_P26', 'DAT_vs_FMAug'): 1,
             ('M3FD_P26', 'MidFusion_vs_CleanFT'): 2,
             ('M3FD_P26', 'MidFusionDAT_vs_DAT'): 3,
             ('M3FD_P26', 'DAT_vs_CleanFT'): 4,
             ('M3FD_P26', 'FMAug_vs_CleanFT'): 5,
             ('M3FD_P26', 'MidFusionDAT_vs_CleanFT'): 6,
             ('M3FD_RS13', 'MidFusion_vs_CGCMFclean'): 7,
             ('M3FD_RS13', 'MidFusionDAT_vs_CGCMFDAT'): 8}
    cond_order = {c: i for i, c in enumerate(
        ['p2_blur_31', 'p2_misalign_60', 'p2_misalign_100', 'p2_missing_rect90',
         'p2_vcombo_hotpx_misalign', 'p2_vcombo_noise_misalign', 'p2_vcombo_blur_quant',
         'combined', 'misalign_20', 'missing_0p', 'missing_50p', 'missing_75p',
         'noise_s30', 'noise_s50'])}
    rows.sort(key=lambda r: (order[(r['group'], r['comparison'])], cond_order[r['condition']]))

    dst = OUT / 'paired_diff_merged.csv'
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f'wrote {dst} ({len(rows)} rows)')

    print('\n== per comparison: verdict counts over its conditions ==')
    groups = defaultdict(list)
    for r in rows:
        groups[(r['group'], r['comparison'])].append(r)
    for k, rs in groups.items():
        c = defaultdict(int)
        for r in rs:
            c[r['verdict']] += 1
        print(f"{k[0]:10s} {k[1]:28s} n={len(rs):2d} "
              f"robust={c['robust']} frame_only={c['frame_only']} "
              f"seed_only={c['seed_only']} not_established={c['not_established']} "
              f"direction_conflict={c['direction_conflict']}")

    print('\n== cells where the two intervals disagree in direction ==')
    for r in rows:
        if r['both_agree_direction'] == 0 and r['boot_cross'] == 0 and r['seed_t_cross'] == 0:
            print('DISAGREE', r['group'], r['comparison'], r['condition'])


if __name__ == '__main__':
    main()
