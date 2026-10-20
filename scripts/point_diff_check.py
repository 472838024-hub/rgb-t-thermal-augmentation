#!/usr/bin/env python3
"""Fast point-level validation of item-4 arm differences against eval JSONs."""
import csv, json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ci_core as C

J = {
 'FLIR': {
   'cgcmf_dat': '/root/autodl-tmp/p26_flir_eval/flir_p26_test_{arm}_s{seed}.json',
   'cgcmf_fmaug': '/root/autodl-tmp/p26_flir_eval/flir_p26_test_{arm}_s{seed}.json',
 },
 'M3FD_P26': {
   'cgcmf_dat_ub': '/autodl-fs/data/cviu_closeout_p26_eval/m3fd/m3fd_test_{arm}_s{seed}_final_p26.json',
   'cgcmf_fmaug_ub': '/autodl-fs/data/cviu_closeout_p26_eval/m3fd/m3fd_test_{arm}_s{seed}_final_p26.json',
   'cgcmf_cleanft_ub': '/autodl-fs/data/cviu_closeout_p26_eval/m3fd/m3fd_test_{arm}_s{seed}_final_p26.json',
   'midfusion_v2': '/root/autodl-tmp/p26_midfusion_20261018/out/m3fd_test_{arm}_s{seed}_final_p26.json',
   'midfusion_dat_v2': '/root/autodl-tmp/p26_midfusion_20261018/out/m3fd_test_{arm}_s{seed}_final_p26.json',
 },
 'M3FD_RS13': {
   'midfusion': '/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/m3fd_test_{arm}_s{seed}.json',
   'midfusion_dat': '/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/m3fd_test_{arm}_s{seed}.json',
   'cgcmf_clean': '/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/m3fd_test_{arm}_s{seed}.json',
   'cgcmf_dat': '/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/m3fd_test_{arm}_s{seed}.json',
 },
}
COMPS = [
 ('FLIR_P26','DAT_vs_FMAug','cgcmf_dat','cgcmf_fmaug',C.P26_SEVEN),
 ('M3FD_P26','DAT_vs_FMAug','cgcmf_dat_ub','cgcmf_fmaug_ub',C.P26_SEVEN),
 ('M3FD_P26','MidFusion_vs_CleanFT','midfusion_v2','cgcmf_cleanft_ub',C.P26_SEVEN),
 ('M3FD_P26','MidFusionDAT_vs_DAT','midfusion_dat_v2','cgcmf_dat_ub',C.P26_SEVEN),
 ('M3FD_P26','DAT_vs_CleanFT','cgcmf_dat_ub','cgcmf_cleanft_ub',C.P26_SEVEN),
 ('M3FD_P26','FMAug_vs_CleanFT','cgcmf_fmaug_ub','cgcmf_cleanft_ub',C.P26_SEVEN),
 ('M3FD_P26','MidFusionDAT_vs_CleanFT','midfusion_dat_v2','cgcmf_cleanft_ub',C.P26_SEVEN),
 ('M3FD_RS13','MidFusion_vs_CGCMFclean','midfusion','cgcmf_clean',C.RS13_SEVEN),
 ('M3FD_RS13','MidFusionDAT_vs_CGCMFDAT','midfusion_dat','cgcmf_dat',C.RS13_SEVEN),
]
rows=[]
for g,lab,a,b,conds in COMPS:
    ds = 'FLIR' if g=='FLIR_P26' else ('M3FD_P26' if g=='M3FD_P26' else 'M3FD_RS13')
    for c in conds:
        va=[]; vb=[]
        for s in (0,1,2):
            da=json.load(open(J[ds][a].format(arm=a,seed=s)))
            db=json.load(open(J[ds][b].format(arm=b,seed=s)))
            va.append(da['conditions'][c]['map50']); vb.append(db['conditions'][c]['map50'])
        d=[(x-y)*100 for x,y in zip(va,vb)]
        rows.append(dict(group=g,comparison=lab,condition=c,a=a,b=b,
                         a_s0=va[0],a_s1=va[1],a_s2=va[2],b_s0=vb[0],b_s1=vb[1],b_s2=vb[2],
                         point_diff_pp=sum(d)/3))
for r in rows:
    print(f"{r['group']:9s} {r['comparison']:28s} {r['condition']:30s} {r['point_diff_pp']:+9.4f}")
with open('/root/autodl-tmp/ci_reweight_20261020/point_diff_check.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
