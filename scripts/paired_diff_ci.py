#!/usr/bin/env python3
"""Item 4: paired cluster-bootstrap intervals for arm differences.

Delta = arm A - arm B, paired by seed and by bootstrap replicate.  Because both
arms share the same frame multiplicities and the same RGB reference, the
difference of gaps equals the difference of mAP50; the RGB reference cancels
exactly, so it is not evaluated here.

Primary scheme: PPS cluster draws + Horvitz-Thompson frame weights, identical to
the item-3 reweighting (one draw matrix per dataset, shared by every arm and
condition).  Read-only over existing prediction caches.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ci_core as C

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')
SEED = {'FLIR': 20261024, 'M3FD': 20261025}

# (label, comparison, arm_a, arm_b, conditions, cache spec per arm)
P26 = '/autodl-fs/data/cviu_closeout_p26_eval/m3fd/preds'
MF = '/root/autodl-tmp/p26_midfusion_20261018/out/preds'
FLIR = '/root/autodl-tmp/p26_flir_eval/preds'
RS13 = '/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/preds'

# arm -> (cache_root, run_fmt)
ARMS = {
    'FLIR_P26': {
        'cgcmf_dat': (FLIR, '{arm}_s{seed}'),
        'cgcmf_fmaug': (FLIR, '{arm}_s{seed}'),
    },
    'M3FD_P26': {
        'cgcmf_dat_ub': (P26, '{arm}_s{seed}_final'),
        'cgcmf_fmaug_ub': (P26, '{arm}_s{seed}_final'),
        'cgcmf_cleanft_ub': (P26, '{arm}_s{seed}_final'),
        'midfusion_v2': (MF, '{arm}_s{seed}_final'),
        'midfusion_dat_v2': (MF, '{arm}_s{seed}_final'),
    },
    'M3FD_RS13': {
        'midfusion': (RS13, '{arm}_s{seed}'),
        'midfusion_dat': (RS13, '{arm}_s{seed}'),
        'cgcmf_clean': (RS13, '{arm}_s{seed}'),
        'cgcmf_dat': (RS13, '{arm}_s{seed}'),
    },
}

COMPARISONS = [
    # group, comparison label, arm A, arm B, conditions, protocol
    ('FLIR_P26', 'DAT_vs_FMAug', 'cgcmf_dat', 'cgcmf_fmaug', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_P26', 'DAT_vs_FMAug', 'cgcmf_dat_ub', 'cgcmf_fmaug_ub', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_P26', 'MidFusion_vs_CleanFT', 'midfusion_v2', 'cgcmf_cleanft_ub', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_P26', 'MidFusionDAT_vs_DAT', 'midfusion_dat_v2', 'cgcmf_dat_ub', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_P26', 'DAT_vs_CleanFT', 'cgcmf_dat_ub', 'cgcmf_cleanft_ub', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_P26', 'FMAug_vs_CleanFT', 'cgcmf_fmaug_ub', 'cgcmf_cleanft_ub', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_P26', 'MidFusionDAT_vs_CleanFT', 'midfusion_dat_v2', 'cgcmf_cleanft_ub', C.P26_SEVEN, 'P26 last.pt'),
    ('M3FD_RS13', 'MidFusion_vs_CGCMFclean', 'midfusion', 'cgcmf_clean', C.RS13_SEVEN, 'RS13 best.pt'),
    ('M3FD_RS13', 'MidFusionDAT_vs_CGCMFDAT', 'midfusion_dat', 'cgcmf_dat', C.RS13_SEVEN, 'RS13 best.pt'),
]


def ap_boot(pre, weights):
    """AP50 matrix (n_boot, 3 seeds) for one arm/condition."""
    out = np.empty((weights.shape[0], len(C.SEEDS)), np.float64)
    for b in range(weights.shape[0]):
        w = weights[b]
        for si, s in enumerate(C.SEEDS):
            out[b, si] = C.B.ap50_weighted(pre[s], w)
    return out


def ap_point(pre):
    ones = np.ones(len(next(iter(pre.values()))['gidx']), np.float64)
    return np.array([C.B.ap50_weighted(pre[s], ones) for s in C.SEEDS])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', choices=['FLIR', 'M3FD', 'both'], default='both')
    ap.add_argument('--n-boot', type=int, default=1000)
    ap.add_argument('--out', default=str(OUT / 'paired_diff_ci.csv'))
    a = ap.parse_args()

    # Point differences already extracted from the evaluation JSONs; used as an
    # independent check that the cache-based AP chain reproduces them.
    point_json = {}
    with open(OUT / 'point_diff_check.csv') as f:
        for r in csv.DictReader(f):
            point_json[(r['group'], r['comparison'], r['condition'])] = float(r['point_diff_pp'])

    rows = []
    cols = ['group', 'dataset', 'protocol', 'comparison', 'arm_a', 'arm_b',
            'condition', 'n_frames', 'point_diff_pp', 'ci_lo_pp', 'ci_hi_pp',
            'excludes_zero', 'crosses_zero', 'n_resamples', 'method', 'seed_list']

    def dump():
        tmp = Path(str(a.out) + '.tmp')
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader(); w.writerows(rows)
        tmp.replace(a.out)

    if a.dataset == 'both':
        groups = ['FLIR_P26', 'M3FD_P26', 'M3FD_RS13']
    elif a.dataset == 'FLIR':
        groups = ['FLIR_P26']
    else:
        groups = ['M3FD_P26', 'M3FD_RS13']
    for group in groups:
        ds = 'FLIR' if group == 'FLIR_P26' else 'M3FD'
        ctx = C.make_ctx(ds)
        C.log(f"=== {group} ({ds}): {len(ctx['stems'])} frames, {ctx['ncl']} clusters")
        W = C.draw_pps_ht(ctx, a.n_boot, np.random.default_rng(SEED[ds]))
        # group arm registry
        arm_specs = ARMS[group]
        conds_all = []
        for g, _lab, _a, _b, cs, _p in COMPARISONS:
            if g == group:
                for c in cs:
                    if c not in conds_all:
                        conds_all.append(c)
        for cond in conds_all:
            # Load only the arms needed by comparisons on this condition.
            needed = set()
            pairs = []
            for g, lab, aa, bb, cs, proto in COMPARISONS:
                if g != group or cond not in cs:
                    continue
                needed.add(aa); needed.add(bb); pairs.append((lab, aa, bb, proto))
            # One arm at a time: keep only its small (n_boot, 3) AP matrix.
            apm, pt = {}, {}
            for arm in sorted(needed):
                cache, fmt = arm_specs[arm]
                pre = C.load_pre(ctx, cache, fmt, arm, cond)
                pt[arm] = ap_point(pre)
                apm[arm] = ap_boot(pre, W)
                del pre
            for lab, aa, bb, proto in pairs:
                d_boot = (apm[aa] - apm[bb]).mean(axis=1) * 100.0
                d_point = float((pt[aa] - pt[bb]).mean()) * 100.0
                want = point_json[(group, lab, cond)]
                if abs(d_point - want) > 1e-6:
                    raise SystemExit(f'point mismatch {group} {lab} {cond}: '
                                     f'caches {d_point:.9f} vs JSON {want:.9f}')
                lo, hi = C.percentile_ci(d_boot)
                rows.append(dict(
                    group=group, dataset=ds, protocol=proto, comparison=lab,
                    arm_a=aa, arm_b=bb, condition=cond,
                    n_frames=len(ctx['stems']),
                    point_diff_pp=f'{d_point:.6f}',
                    ci_lo_pp=f'{lo:.6f}', ci_hi_pp=f'{hi:.6f}',
                    excludes_zero=int(lo < 0 < hi) == 0,
                    crosses_zero=int(lo < 0 < hi),
                    n_resamples=a.n_boot, method='paired-cluster-bootstrap-pps-ht',
                    seed_list='0,1,2',
                ))
                C.log(f"{group:9s} {lab:28s} {cond:30s} d {d_point:+8.3f} "
                      f"[{lo:+8.3f},{hi:+8.3f}] cross={int(lo<0<hi)}")
            del apm, pt
            dump()
            C.log(f"  partial: {len(rows)} rows -> {a.out}")

    dump()
    C.log(f'wrote {a.out} ({len(rows)} rows)')


if __name__ == '__main__':
    main()
