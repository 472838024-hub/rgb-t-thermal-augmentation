#!/usr/bin/env python3
"""Render the 28-cell item-3 result as compact Markdown tables."""
import csv
from pathlib import Path

OUT = Path('/root/autodl-tmp/ci_reweight_20261020')
ROWS = list(csv.DictReader((OUT / 'reweight_ci_28.csv').open()))
ORDER = ['p2_blur_31', 'p2_misalign_60', 'p2_misalign_100', 'p2_missing_rect90',
         'p2_vcombo_hotpx_misalign', 'p2_vcombo_noise_misalign',
         'p2_vcombo_blur_quant']


def main():
    lines = []
    for ds in ('FLIR', 'M3FD'):
        for arm in ('cgcmf_dat', 'cgcmf_fmaug'):
            rs = {r['condition']: r for r in ROWS
                  if r['dataset'] == ds and r['arm'] == arm}
            lines.append(f'### {ds} / {arm}')
            lines.append('')
            lines.append('| condition | point (pp) | old 95% CI | PPS+HT 95% CI | rel old→pps | cmean 95% CI (sensitivity) |')
            lines.append('|---|---:|---|---|---:|---|')
            for c in ORDER:
                r = rs[c]
                lines.append(
                    f"| `{c}` | {float(r['point_pooled_pp']):+.3f} | "
                    f"[{float(r['old_ci_lo_pp']):+.3f}, {float(r['old_ci_hi_pp']):+.3f}] | "
                    f"[{float(r['pps_lo_pp']):+.3f}, {float(r['pps_hi_pp']):+.3f}] | "
                    f"{float(r['old_rel']):.2f}→{float(r['pps_rel']):.2f} | "
                    f"[{float(r['cmean_lo_pp']):+.2f}, {float(r['cmean_hi_pp']):+.2f}] |")
            lines.append('')
    dst = OUT / 'item3_tables.md'
    dst.write_text('\n'.join(lines))
    print(f'wrote {dst} ({len(lines)} lines)')


if __name__ == '__main__':
    main()
