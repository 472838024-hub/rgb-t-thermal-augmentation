#!/usr/bin/env python3
"""Shared plumbing for the 2026-10-20 cluster-weighted CI recomputation.

Read-only over existing prediction caches.  The AP50 chain, GT loading and
cluster definitions are imported verbatim from bootstrap_ci_20261017 (B6 port,
validated bit-exact against the recorded eval JSONs there).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, '/root/autodl-tmp/bootstrap_ci_20261017')
import bootstrap_ci as B  # noqa: E402

SEEDS = [0, 1, 2]

P26_SEVEN = ['p2_blur_31', 'p2_misalign_60', 'p2_misalign_100',
             'p2_missing_rect90', 'p2_vcombo_hotpx_misalign',
             'p2_vcombo_noise_misalign', 'p2_vcombo_blur_quant']
RS13_SEVEN = ['combined', 'misalign_20', 'missing_0p', 'missing_50p',
              'missing_75p', 'noise_s30', 'noise_s50']


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def _flir_stems(root):
    return sorted((p.stem for p in (root / 'labels/test').glob('*.txt')),
                  key=lambda s: int(s.split('_')[-1]))


def _m3fd_stems(root):
    return sorted(p.stem for p in (root / 'labels/test').glob('*.txt'))


def make_ctx(dataset):
    if dataset == 'FLIR':
        root = Path('/root/autodl-tmp/FLIR_aligned_m3fd_v2')
        ctx = dict(label='FLIR', root=root, W=640, H=512, NC=6,
                   rgb_root=Path('/root/pcb_project/experiments/cviu_2026/test_eval/flir/preds'),
                   cluster_fn=B.flir_clusters, stems_fn=_flir_stems)
    elif dataset == 'M3FD':
        root = Path('/root/autodl-tmp/M3FD/M3FD_split_v2')
        ctx = dict(label='M3FD', root=root, W=1024, H=768, NC=6,
                   rgb_root=Path('/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/preds'),
                   cluster_fn=B.m3fd_clusters, stems_fn=_m3fd_stems)
    else:
        raise ValueError(dataset)
    ctx['stems'] = ctx['stems_fn'](ctx['root'])
    ctx['cluster'], ctx['ncl'] = ctx['cluster_fn'](ctx['stems'])
    ctx['sizes'] = np.bincount(ctx['cluster'])
    ctx['gidx'] = {B.image_id_of(s): i for i, s in enumerate(ctx['stems'])}
    ctx['gt'], ctx['counts'] = B.load_gt(ctx['root'], ctx['stems'],
                                         ctx['W'], ctx['H'], ctx['NC'])
    return ctx


def load_rgb(ctx):
    rgb = {}
    for s in SEEDS:
        p = ctx['rgb_root'] / f'rgb_s{s}' / 'clean' / 'predictions.json'
        t0 = time.time()
        rgb[s] = B.precompute(p, ctx['stems'], ctx['gidx'], ctx['gt'],
                              ctx['counts'], ctx['NC'])
        log(f"  pre rgb s{s}: {rgb[s]['n_rec']} records ({time.time()-t0:.1f}s)")
    return rgb


def load_pre(ctx, cache_root, run_fmt, arm, cond):
    pre = {}
    for s in SEEDS:
        p = Path(cache_root) / run_fmt.format(arm=arm, seed=s) / cond / 'predictions.json'
        if not p.exists():
            raise FileNotFoundError(p)
        t0 = time.time()
        pre[s] = B.precompute(p, ctx['stems'], ctx['gidx'], ctx['gt'],
                              ctx['counts'], ctx['NC'])
        log(f"  pre {arm} s{s} {cond}: {pre[s]['n_rec']} records "
            f"({time.time()-t0:.1f}s)")
    return pre


def point_gap(pre, rgb):
    ones = np.ones(len(next(iter(rgb.values()))['gidx']), np.float64)
    return float(np.mean([B.ap50_weighted(pre[s], ones) -
                          B.ap50_weighted(rgb[s], ones) for s in SEEDS])) * 100.0


def per_cluster_gaps(pre, rgb, ctx):
    out = np.empty(ctx['ncl'], np.float64)
    for k in range(ctx['ncl']):
        w = np.zeros(len(ctx['stems']), np.float64)
        w[ctx['cluster'] == k] = 1.0
        out[k] = float(np.mean([B.ap50_weighted(pre[s], w) -
                                B.ap50_weighted(rgb[s], w) for s in SEEDS])) * 100.0
    return out


def draw_uniform(ctx, n_boot, rng):
    """Integer frame multiplicities, clusters drawn with equal probability."""
    K, cl = ctx['ncl'], ctx['cluster']
    cnt = rng.multinomial(K, np.full(K, 1.0 / K), size=n_boot)
    return cnt[:, cl].astype(np.float64)          # (n_boot, n_frames)


def draw_pps_ht(ctx, n_boot, rng):
    """PPS cluster draws + Horvitz-Thompson inverse-probability frame weights.

    Each of the K draws selects cluster i with probability p_i = n_i / n.
    A selected cluster's frames get multiplicity c_i / (K * p_i); in expectation
    every frame keeps weight 1, so the resample targets the same frame-weighted
    estimand as the reported pooled point estimate.
    """
    K, cl = ctx['ncl'], ctx['cluster']
    p = ctx['sizes'] / ctx['sizes'].sum()
    cnt = rng.multinomial(K, p, size=n_boot).astype(np.float64)
    wclus = cnt / (K * p)                          # (n_boot, K)
    return wclus[:, cl]


def draw_cluster_pairs(ctx, n_boot, rng, scheme):
    """Return (n_boot, n_frames) frame-weight matrix."""
    if scheme == 'uniform':
        return draw_uniform(ctx, n_boot, rng)
    if scheme == 'pps_ht':
        return draw_pps_ht(ctx, n_boot, rng)
    raise ValueError(scheme)


def rgb_ap_matrix(rgb, weights):
    """RGB reference AP50 for every replicate row and seed; shape (n_boot, 3)."""
    out = np.empty((weights.shape[0], len(SEEDS)), np.float64)
    for b in range(weights.shape[0]):
        w = weights[b]
        for si, s in enumerate(SEEDS):
            out[b, si] = B.ap50_weighted(rgb[s], w)
    return out


def gap_matrix(pre, rgb, ctx, weights, rgb_ap=None):
    """Gap (pp) for every replicate row of `weights`; shape (n_boot,).

    Pass `rgb_ap` from rgb_ap_matrix() to avoid recomputing the shared reference.
    """
    if rgb_ap is None:
        rgb_ap = rgb_ap_matrix(rgb, weights)
    out = np.empty(weights.shape[0], np.float64)
    for b in range(weights.shape[0]):
        w = weights[b]
        out[b] = np.mean([B.ap50_weighted(pre[s], w) for s in SEEDS] -
                         rgb_ap[b]) * 100.0
    return out


def percentile_ci(a, lo=2.5, hi=97.5):
    return tuple(float(x) for x in np.percentile(np.asarray(a), [lo, hi]))
