#!/usr/bin/env python3
"""Verify Table-2 mean±SD (ddof=1) from raw P26 Mid-Fusion eval JSONs; per-seed raw values."""
import csv, json, math, statistics as st
from pathlib import Path
OUT=Path('/root/autodl-tmp/p26_midfusion_20261018/out')
ARMS={'midfusion_v2':'Mid-Fusion','midfusion_dat_v2':'Mid-Fusion+DAT'}
CONDS=['p2_blur_31','p2_misalign_60','p2_misalign_100','p2_missing_rect90',
       'p2_vcombo_blur_quant','p2_vcombo_hotpx_misalign','p2_vcombo_noise_misalign']
rows=[]
for arm,label in ARMS.items():
    for cond in CONDS:
        fs=[];rs=[];gs=[]
        for s in (0,1,2):
            j=json.loads((OUT/f'm3fd_test_{arm}_s{s}_final_p26.json').read_text())
            f=j['conditions'][cond]['map50']; r=j['rgb_flat_reference']['map50']
            fs.append(f); rs.append(r); gs.append((f-r)*100)
        sf=st.stdev(fs); sr=st.stdev(rs); sg=st.stdev(gs)
        r_imp=(sf*sf+sr*sr-(sg/100)**2)/(2*sf*sr) if sf*sr else float('nan')
        rows.append(dict(arm=arm,label=label,condition=cond,
            map50_s0=fs[0],map50_s1=fs[1],map50_s2=fs[2],
            rgb_s0=rs[0],rgb_s1=rs[1],rgb_s2=rs[2],
            gap_s0=gs[0],gap_s1=gs[1],gap_s2=gs[2],
            mean=st.mean(gs),sd_ddof1=sg,sem=sg/math.sqrt(3),
            t95_halfwidth_df2=4.302652729911275*sg/math.sqrt(3),
            implied_corr_fusion_rgb=r_imp,
            implied_corr_gap_rgb=(r_imp*sf*sr-sr*sr)/((sg/100)*sr) if sg*sr else float('nan'),
            fusion_sd_pp=sf*100,rgb_sd_pp=sr*100))
        print(f'{label:14s} {cond:30s} fusion {fs[0]:.6f} {fs[1]:.6f} {fs[2]:.6f} | rgb {rs[0]:.6f} {rs[1]:.6f} {rs[2]:.6f} | gap {gs[0]:+.4f} {gs[1]:+.4f} {gs[2]:+.4f} -> {st.mean(gs):+.4f}±{sg:.4f}  r_implied={r_imp:+.4f}')
# compare to published summary
pub={ (r['arm'],r['condition']):(float(r['gap_pp_mean']),float(r['gap_pp_sd_ddof1'])) for r in csv.DictReader(open('/root/autodl-tmp/p26_midfusion_20261018/p26_midfusion_seven_summary.csv'))}
print('\n--- comparison with p26_midfusion_seven_summary.csv ---')
bad=0
for r in rows:
    pm,ps=pub[(r['arm'],r['condition'])]
    dm=abs(pm-r['mean']); ds=abs(ps-r['sd_ddof1'])
    flag='OK' if dm<1e-9 and ds<1e-9 else 'DIFF'
    if flag!='OK': bad+=1
    print(f"{r['label']:14s} {r['condition']:30s} mean diff {dm:.2e} sd diff {ds:.2e} {flag}")
print('mismatches:',bad)
Path('/root/autodl-tmp/ci_reweight_20261020/table2_per_seed_verify.json').write_text(json.dumps(rows,indent=1))
with open('/root/autodl-tmp/ci_reweight_20261020/table2_per_seed_verify.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('wrote CSV/JSON')
