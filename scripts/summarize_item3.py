#!/usr/bin/env python3
"""Compact summary of the 28-cell reweighted intervals."""
import csv
from collections import Counter
from pathlib import Path

import numpy as np

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')
SRC = OUT / 'reweight_ci_28.csv'


def classify(r, lo, hi):
    if float(r[hi]) < 0:
        return 'neg'
    if float(r[lo]) > 0:
        return 'pos'
    return 'cross'


def main():
    rows = list(csv.DictReader(SRC.open()))
    print(f'rows={len(rows)}')
    for name, lo, hi in [('old', 'old_ci_lo_pp', 'old_ci_hi_pp'),
                         ('pps_ht', 'pps_lo_pp', 'pps_hi_pp'),
                         ('cluster_mean', 'cmean_lo_pp', 'cmean_hi_pp')]:
        c = Counter(classify(r, lo, hi) for r in rows)
        rel = np.array([float(r[{'old': 'old_rel', 'pps_ht': 'pps_rel',
                                 'cluster_mean': 'cmean_rel'}[name]]) for r in rows])
        hw = np.array([(float(r[hi]) - float(r[lo])) / 2 for r in rows])
        print(f'\n{name}: neg={c["neg"]} pos={c["pos"]} cross={c["cross"]}')
        print(f'  rel: median={np.median(rel):.3f} mean={rel.mean():.3f} '
              f'min={rel.min():.3f} max={rel.max():.3f} | rel>0.5: {np.sum(rel>0.5)}/28 '
              f'rel>0.7: {np.sum(rel>0.7)}/28 rel>0.9: {np.sum(rel>0.9)}/28')
        print(f'  half-width: median={np.median(hw):.3f} mean={hw.mean():.3f} '
              f'min={hw.min():.3f} max={hw.max():.3f}')
        for ds in ('FLIR', 'M3FD'):
            cds = Counter(classify(r, lo, hi) for r in rows if r['dataset'] == ds)
            rd = np.array([float(r[{'old': 'old_rel', 'pps_ht': 'pps_rel',
                                    'cluster_mean': 'cmean_rel'}[name]])
                           for r in rows if r['dataset'] == ds])
            print(f'  {ds}: neg={cds["neg"]} pos={cds["pos"]} cross={cds["cross"]} '
                  f'rel median={np.median(rd):.3f}')

    print('\n== classification changes vs old ==')
    for r in rows:
        old = classify(r, 'old_ci_lo_pp', 'old_ci_hi_pp')
        pps = classify(r, 'pps_lo_pp', 'pps_hi_pp')
        cme = classify(r, 'cmean_lo_pp', 'cmean_hi_pp')
        if pps != old or cme != old:
            print(f"{r['dataset']:4s} {r['arm']:12s} {r['condition']:30s} "
                  f"old={old:5s} pps={pps:5s} cmean={cme:5s} "
                  f"point={float(r['point_pooled_pp']):+8.3f} "
                  f"old[{float(r['old_ci_lo_pp']):+8.3f},{float(r['old_ci_hi_pp']):+8.3f}] "
                  f"pps[{float(r['pps_lo_pp']):+8.3f},{float(r['pps_hi_pp']):+8.3f}]")

    print('\n== cells with the largest PPS rel (still upper-half) ==')
    for r in sorted(rows, key=lambda x: -float(x['pps_rel']))[:6]:
        print(f"{r['dataset']:4s} {r['arm']:12s} {r['condition']:30s} "
              f"rel={float(r['pps_rel']):.3f} point={float(r['point_pooled_pp']):+8.3f} "
              f"pps[{float(r['pps_lo_pp']):+8.3f},{float(r['pps_hi_pp']):+8.3f}]")

    print('\n== bootstrap mean/median vs point (centering check) ==')
    for name, bm, bmed in [('pps_ht', 'pps_boot_mean', 'pps_boot_median'),
                           ('cluster_mean', 'cmean_boot_mean', 'cmean_boot_median')]:
        dif = np.array([float(r[bm]) - float(r['point_pooled_pp']) for r in rows])
        print(f'{name}: mean(boot)-point  median={np.median(dif):+.3f} '
              f'min={dif.min():+.3f} max={dif.max():+.3f}')


if __name__ == '__main__':
    main()
