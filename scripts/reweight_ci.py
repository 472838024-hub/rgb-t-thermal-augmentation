#!/usr/bin/env python3
"""Item 3: recompute the 28 Table-5 intervals with cluster weighting.

Primary scheme (pps_ht): draw each cluster with probability proportional to its
frame count, then weight every selected cluster by 1/(K*p_i) (Horvitz-Thompson).
The resample therefore targets the same frame-weighted estimand as the reported
pooled point estimate, while giving small clusters their correct influence.

Secondary scheme (cluster_mean): replace the pooled AP with the cluster-size
weighted mean of per-cluster gaps and bootstrap clusters uniformly.  This is a
different estimand (macro over clusters, size-weighted) and is reported only as
a sensitivity check.

Read-only over existing prediction caches.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ci_core as C

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')
OLD = Path('/root/autodl-tmp/bootstrap_ci_20261018/table5_seven_ci.csv')
ARMS = {
    'FLIR': [('cgcmf_dat', '/root/autodl-tmp/p26_flir_eval/preds', '{arm}_s{seed}'),
             ('cgcmf_fmaug', '/root/autodl-tmp/p26_flir_eval/preds', '{arm}_s{seed}')],
    'M3FD': [('cgcmf_dat', '/autodl-fs/data/cviu_closeout_p26_eval/m3fd/preds',
              '{arm}_ub_s{seed}_final'),
             ('cgcmf_fmaug', '/autodl-fs/data/cviu_closeout_p26_eval/m3fd/preds',
              '{arm}_ub_s{seed}_final')],
}
SEED = {'FLIR': 20261024, 'M3FD': 20261025}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', choices=['FLIR', 'M3FD', 'both'], default='both')
    ap.add_argument('--n-boot', type=int, default=1000)
    ap.add_argument('--out', default=str(OUT / 'reweight_ci_28.csv'))
    a = ap.parse_args()

    old = {(r['dataset'], r['arm'], r['condition']): r
           for r in csv.DictReader(OLD.open())}
    rows = []
    for ds in (['FLIR', 'M3FD'] if a.dataset == 'both' else [a.dataset]):
        ctx = C.make_ctx(ds)
        C.log(f"=== {ds}: {len(ctx['stems'])} frames, {ctx['ncl']} clusters, "
              f"sizes {ctx['sizes'].tolist()}")
        rgb = C.load_rgb(ctx)
        p = ctx['sizes'] / ctx['sizes'].sum()
        rng_b = np.random.default_rng(SEED[ds])
        rng_c = np.random.default_rng(SEED[ds] + 100)
        # One draw matrix per dataset/scheme, shared across arms and conditions.
        W_b = C.draw_pps_ht(ctx, a.n_boot, rng_b)
        C.log('  precomputing shared RGB reference APs for all resamples')
        rgb_ap = C.rgb_ap_matrix(rgb, W_b)
        cnt_c = rng_c.multinomial(ctx['ncl'], np.full(ctx['ncl'], 1.0 / ctx['ncl']),
                                  size=a.n_boot).astype(np.float64)

        for arm, cache, fmt in ARMS[ds]:
            for cond in C.P26_SEVEN:
                pre = C.load_pre(ctx, cache, fmt, arm, cond)
                pt = C.point_gap(pre, rgb)
                o = old[(ds, arm, cond)]
                old_pt = float(o['mean_gap_pp'])
                if abs(pt - old_pt) > 1e-6:
                    raise SystemExit(f'point mismatch {ds} {arm} {cond}: {pt} vs {old_pt}')

                gaps_b = C.gap_matrix(pre, rgb, ctx, W_b, rgb_ap=rgb_ap)
                lo_b, hi_b = C.percentile_ci(gaps_b)

                gk = C.per_cluster_gaps(pre, rgb, ctx)
                wpt = float(np.dot(p, gk))
                # Unbiased resample statistic for theta = sum_i p_i g_i:
                # E[c_i] = 1 under uniform cluster draws, so sum_i c_i p_i g_i
                # has expectation theta.  (Dividing by K would give theta/K.)
                lin = (cnt_c * p * gk).sum(axis=1)
                lo_c, hi_c = C.percentile_ci(lin)

                old_lo, old_hi = float(o['ci_lo_pp']), float(o['ci_hi_pp'])
                row = dict(
                    dataset=ds, arm=arm, condition=cond,
                    n_frames=len(ctx['stems']), n_clusters=ctx['ncl'],
                    cluster_sizes=' '.join(str(x) for x in ctx['sizes'].tolist()),
                    point_pooled_pp=f'{pt:.6f}',
                    old_ci_lo_pp=f'{old_lo:.6f}', old_ci_hi_pp=f'{old_hi:.6f}',
                    old_rel=f'{(old_pt-old_lo)/(old_hi-old_lo):.4f}',
                    old_cross=int(old_lo < 0 < old_hi),
                    pps_lo_pp=f'{lo_b:.6f}', pps_hi_pp=f'{hi_b:.6f}',
                    pps_rel=f'{(pt-lo_b)/(hi_b-lo_b):.4f}',
                    pps_cross=int(lo_b < 0 < hi_b),
                    pps_halfwidth=f'{(hi_b-lo_b)/2:.6f}',
                    pps_boot_mean=f'{float(np.mean(gaps_b)):.6f}',
                    pps_boot_median=f'{float(np.median(gaps_b)):.6f}',
                    wpoint_cluster_mean_pp=f'{wpt:.6f}',
                    cmean_lo_pp=f'{lo_c:.6f}', cmean_hi_pp=f'{hi_c:.6f}',
                    cmean_rel=f'{(wpt-lo_c)/(hi_c-lo_c):.4f}',
                    cmean_cross=int(lo_c < 0 < hi_c),
                    cmean_boot_mean=f'{float(np.mean(lin)):.6f}',
                    cmean_boot_median=f'{float(np.median(lin)):.6f}',
                    cmean_delta_point_vs_pooled=f'{wpt-pt:+.6f}',
                    n_resamples=a.n_boot,
                )
                rows.append(row)
                C.log(f"{ds} {arm} {cond:30s} point {pt:+9.4f} | old [{old_lo:+8.3f},{old_hi:+8.3f}] "
                      f"rel {(old_pt-old_lo)/(old_hi-old_lo):.2f} | pps [{lo_b:+8.3f},{hi_b:+8.3f}] "
                      f"rel {(pt-lo_b)/(hi_b-lo_b):.2f} cross={int(lo_b<0<hi_b)} | "
                      f"cmean point {wpt:+8.3f} [{lo_c:+8.3f},{hi_c:+8.3f}] cross={int(lo_c<0<hi_c)}")
                del pre

    cols = list(rows[0].keys())
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    C.log(f'wrote {a.out} ({len(rows)} rows)')

    for name, lo, hi, cross in [
            ('old', 'old_ci_lo_pp', 'old_ci_hi_pp', 'old_cross'),
            ('pps_ht', 'pps_lo_pp', 'pps_hi_pp', 'pps_cross'),
            ('cluster_mean', 'cmean_lo_pp', 'cmean_hi_pp', 'cmean_cross')]:
        nneg = sum(1 for r in rows if float(r[hi]) < 0)
        npos = sum(1 for r in rows if float(r[lo]) > 0)
        ncross = sum(1 for r in rows if r[cross] == 1)
        C.log(f'SUMMARY {name}: fully-negative {nneg}/28, fully-positive {npos}/28, '
              f'cross-zero {ncross}/28')


if __name__ == '__main__':
    main()
