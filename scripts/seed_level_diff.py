#!/usr/bin/env python3
"""Item 4 companion: seed-level paired t intervals (n=3, df=2).

Reads the per-seed mAP50 values already extracted from the evaluation JSONs
(point_diff_check.csv, bit-exact against the JSONs) and computes, for every
comparison x condition:

    d_s = arm_A(s) - arm_B(s)            (pp)
    mean, SD(ddof=1), SEM, 95% t interval with t_{0.975,2} = 4.3026527

This is the interval that answers "the ordering rests on 3 seeds" directly.
The paired cluster bootstrap (paired_diff_ci.py) answers a different question:
test-frame/cluster uncertainty given the three trained model pairs.  The two
are complementary; neither replaces the other.
"""
import csv
from pathlib import Path

import numpy as np

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')
SRC = OUT / 'point_diff_check.csv'
T975_DF2 = 4.302652729911275


def main():
    rows = list(csv.DictReader(SRC.open()))
    out = []
    for r in rows:
        d = np.array([float(r[f'a_s{s}']) - float(r[f'b_s{s}']) for s in (0, 1, 2)]) * 100.0
        m = float(d.mean())
        sd = float(d.std(ddof=1))
        sem = sd / np.sqrt(3)
        hw = T975_DF2 * sem
        out.append(dict(
            group=r['group'], comparison=r['comparison'], condition=r['condition'],
            arm_a=r['a'], arm_b=r['b'],
            d_s0_pp=f'{d[0]:.6f}', d_s1_pp=f'{d[1]:.6f}', d_s2_pp=f'{d[2]:.6f}',
            mean_diff_pp=f'{m:.6f}', sd_ddof1_pp=f'{sd:.6f}',
            sem_pp=f'{sem:.6f}', t95_halfwidth_pp=f'{hw:.6f}',
            seed_t_lo_pp=f'{m - hw:.6f}', seed_t_hi_pp=f'{m + hw:.6f}',
            seed_t_cross=int((m - hw) < 0 < (m + hw)),
            n_positive_seeds=int((d > 0).sum()),
            n_negative_seeds=int((d < 0).sum()),
            all3_same_sign=int(((d > 0).all() or (d < 0).all())),
        ))
    cols = list(out[0].keys())
    dst = OUT / 'seed_level_diff.csv'
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    print(f'wrote {dst} ({len(out)} rows)')

    print('\n== per comparison: conditions with all 3 seeds same sign / with seed-level t interval excluding 0 ==')
    seen = []
    for r in out:
        k = (r['group'], r['comparison'])
        if k not in seen:
            seen.append(k)
    for k in seen:
        rs = [r for r in out if (r['group'], r['comparison']) == k]
        n_same = sum(r['all3_same_sign'] for r in rs)
        n_t = sum(1 for r in rs if r['seed_t_cross'] == 0)
        print(f"{k[0]:10s} {k[1]:28s} same-sign {n_same}/{len(rs)}  "
              f"seed-t excludes 0: {n_t}/{len(rs)}")


if __name__ == '__main__':
    main()
