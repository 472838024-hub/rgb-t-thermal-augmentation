# RECEIPT — 事项 3/4/5 复核（2026-10-20）

- 计算日期：2026-10-20（服务器，纯 CPU）
- 性质：**只读统计**。没有训练、没有新评测、没有改任何评测 JSON / yaml / 权重 / 预测缓存；test 未触碰。
- 运行环境：cgroup 0.5 CPU / 2 GiB；numpy `default_rng`。

## 数据源

- FLIR P26：`/root/autodl-tmp/p26_flir_eval/preds/{cgcmf_dat,cgcmf_fmaug}_s{0,1,2}/<cond>/predictions.json`
- M3FD P26：`/autodl-fs/data/cviu_closeout_p26_eval/m3fd/preds/{cgcmf_dat_ub,cgcmf_fmaug_ub,cgcmf_cleanft_ub}_s{0,1,2}_final/<cond>/predictions.json`
- M3FD Mid-Fusion：`/root/autodl-tmp/p26_midfusion_20261018/out/preds/{midfusion_v2,midfusion_dat_v2}_s{0,1,2}_final/<cond>/predictions.json`
- M3FD RS13：`/root/pcb_project/experiments/cviu_2026/test_eval/m3fd/preds/{midfusion,midfusion_dat,cgcmf_clean,cgcmf_dat}_s{0,1,2}/<cond>/predictions.json`
- 点差核对：既有评测 JSON（`point_diff_check.csv`），缓存链路逐格断言一致（容差 1e-6）。
- 表 2：`/root/autodl-tmp/p26_midfusion_20261018/out/m3fd_test_{arm}_s{seed}_final_p26.json`
  与 RGB 参照 `/root/autodl-tmp/step2_delivery_20261008/04_eval_json/m3fd_test_rgb_s{seed}.json`
  （6/6 sha256 与评测 JSON 内记录一致）。

## 方法

- 事项 3：28 格（7 条件 × 2 数据集 × 2 配方）cluster bootstrap，1,000 次，2.5/97.5 分位。
  主口径 **PPS+HT**：簇抽中概率 `p_i=n_i/n`，帧权重 `c_i/(K p_i)`；
  敏感性口径 **cluster_mean**：目标量 `theta=sum_i p_i g_i`，统计量 `sum_i c_i p_i g_i`。
  AP50 链路复用 B6 Ultralytics numpy 移植（与 2026-10-17/18 交付同一实现）。
- 事项 4：帧级**配对** cluster bootstrap（A/B 同轮共用权重，RGB 参照精确抵消），
  加种子级 n=3、df=2、t=4.30265 的配对 t 区间；两套区间并列判定。
- 事项 5：14 个原始评测 JSON 重算 mean/SD；`r` 用 `r=(s_f^2+s_r^2-s_g^2)/(2 s_f s_r)`
  反解，并用 6 种种子的置换分布校准。

## 结果速览

- **事项 3**：PPS+HT 下 28 格 **27 全负、1 跨零**（M3FD DAT `p2_blur_31`
  [−10.613,+0.321]）；cluster_mean 同为 27/28，跨零格换成 `p2_missing_rect90`
  [−3.430,+3.002]；**26 格三口径一致全负**。正文“26 格全负”改为 27/28（主口径）。
- **事项 4**：63 格 **robust 30 / frame_only 24 / not_established 9**。
  各臂 > CleanFT 最稳（FM-Aug 6/7、MidFusion+DAT 6/7、DAT 5/7 robust）；
  第二载体 > CGCMFclean（RS13）仅 2/7 robust；DAT−FM-Aug 方向随条件翻转；
  MidFusion+DAT−DAT 与 MidFusion+DAT−CGCMFDAT 均 0/7 robust。
- **事项 5**：±2.69 精确复现（ddof=1；ddof=0 为 2.2004）。`p2_blur_31` 逐种子
  gap（Mid-Fusion）−4.6706/−7.5754/−2.1911；r=−0.9021 是 `gap=fusion−RGB`
  且 `s_g >> s_f,s_r` 的代数结果，置换 p=1/3。种子级 t 区间 [−11.51,+1.88] 跨零；
  Mid-Fusion 6/7 格在零下，Mid-Fusion+DAT 7/7 跨零。

## 文件

全部 24 个文件的 sha256 见 `MANIFEST.sha256`（生成时按路径排除清单自身）。
逐格数据：`reweight_ci_28.csv`、`paired_diff_ci.csv`、`seed_level_diff.csv`、
`paired_diff_merged.csv`、`point_diff_check.csv`、`table2_per_seed_verify.csv`、
`table2_corr_permutation.csv`；报告：`README_items_3_4_5_20261020.md`、
`item3_tables.md`、`table2_seed_intervals.md`。

## 复现（只读；输出只写工作目录）

```bash
cd /root/autodl-tmp/ci_reweight_20261020
python3 -u -I reweight_ci.py --dataset both --n-boot 1000 > reweight_ci.log 2>&1
python3 -I point_diff_check.py
python3 -u -I paired_diff_ci.py --n-boot 1000 > paired_diff_ci.log 2>&1
python3 -I seed_level_diff.py
python3 -I merge_item4.py
python3 -I verify_table2.py
python3 -I table2_corr_permutation.py
python3 -I summarize_item3.py
```

脚本在打包目录内，可直接 `cd` 到解包目录后按同样顺序运行（输入路径为绝对路径）。

## 已知限制

- FLIR 6 簇（609 帧占 60.1%）、M3FD 3 簇；区间是簇级不确定性，不是逐图独立样本精度。
- cluster_mean 是不同估计量（macro-over-clusters），其点估计不能当论文数字。
- 种子级 t 只有 df=2；3 个种子的符号检验最小双侧 p=0.25，“显著”须注明抽样层。
- 两套协议（P26 `last.pt` / RS13 `best.pt`）分别判定，未合并。
- 帧级与种子级回答不同问题；只有两层都排除零才写成不依赖抽样层的排序。
