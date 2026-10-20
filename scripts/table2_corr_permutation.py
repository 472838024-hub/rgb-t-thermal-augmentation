#!/usr/bin/env python3
"""Item 5: how surprising is the implied within-seed correlation r ~ -1?

With n=3 seeds, the pairing is one of only 6 permutations.  For each arm x
condition row of Table 2 we report the observed Pearson r between the fusion
mAP50 and the RGB reference mAP50, plus all 6 permutation values of r and the
exact permutation tail count |r_perm| >= |r_obs|.  This shows that |r| ~ 1 is
a small-sample artefact of three points, not evidence of a computation error.
"""
import csv
import itertools
import math
from pathlib import Path

import numpy as np

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')
SRC = OUT / 'table2_per_seed_verify.csv'


def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def main():
    rows = list(csv.DictReader(SRC.open()))
    out = []
    for r in rows:
        f = [float(r[f'map50_s{s}']) for s in (0, 1, 2)]
        rgb = [float(r[f'rgb_s{s}']) for s in (0, 1, 2)]
        g = [float(r[f'gap_s{s}']) for s in (0, 1, 2)]
        r_obs = pearson(f, rgb)
        r_gap = pearson(g, rgb)
        perms = [pearson(f, [rgb[i] for i in p]) for p in itertools.permutations(range(3))]
        tail = sum(1 for x in perms if abs(x) >= abs(r_obs) - 1e-12)
        out.append(dict(
            arm=r['arm'], condition=r['condition'],
            r_obs_fusion_rgb=f'{r_obs:.9f}',
            r_obs_gap_rgb=f'{r_gap:.9f}',
            r_perm_min=f'{min(perms):.9f}', r_perm_max=f'{max(perms):.9f}',
            n_perm_abs_ge_obs=f'{tail}/6',
            perm_p_two_sided=f'{tail/6:.4f}',
            n_perm_abs_ge_0p9=sum(1 for x in perms if abs(x) >= 0.9),
            n_perm_abs_ge_0p99=sum(1 for x in perms if abs(x) >= 0.99),
            perm_values=' '.join(f'{x:+.6f}' for x in perms),
        ))
    dst = OUT / 'table2_corr_permutation.csv'
    with open(dst, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    print(f'wrote {dst} ({len(out)} rows)')
    for r in out:
        print(f"{r['arm']:16s} {r['condition']:30s} r={r['r_obs_fusion_rgb']} "
              f"gap_r={r['r_obs_gap_rgb']} perm[{r['r_perm_min']},{r['r_perm_max']}] "
              f"|r|>=obs {r['n_perm_abs_ge_obs']} >=0.9:{r['n_perm_abs_ge_0p9']} "
              f">=0.99:{r['n_perm_abs_ge_0p99']}")


if __name__ == '__main__':
    main()
