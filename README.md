# Accuracy recovery and residual failures in thermal augmentation for RGB-T detection

Code and result files for the manuscript

> **Thermal Augmentation in RGB-T Detection: Accuracy Recovery and Residual Failures**

The package contains the statistical analysis used in the paper: cluster-weighted
confidence intervals for the 28 condition-level deficit cells, paired-difference
intervals for 63 comparisons, and the per-seed verification of the continuation
results. Training and evaluation code, configuration files and checkpoints are
released separately (see *What is not here*).

## What is here

| Folder | Contents |
|---|---|
| `data/` | Result tables as CSV: the 28-cell intervals (`reweight_ci_28.csv`), the 63 paired differences (`paired_diff_merged.csv`, `paired_diff_ci.csv`, `seed_level_diff.csv`), the per-seed continuation values (`table2_per_seed_verify.csv`) and the permutation check for the implied correlation (`table2_corr_permutation.csv`). |
| `scripts/` | The analysis scripts that produce every table above from the condition-level scores. |
| `docs/` | The analysis report (`README_items_3_4_5_20261020.md`), the delivery receipt and the sha256 manifest of the original package. |

## What is not here

- Training and evaluation code, model configuration files and checkpoints. These live in the
  project repository and are released with the article.
- The M3FD and FLIR Aligned images. They remain subject to the terms of their providers
  (Liu et al., 2022; Teledyne FLIR; Zhang et al., 2020).
- The corruption generators for the two condition grids.

## Reproducing the tables

```bash
pip install -r requirements.txt
python scripts/ci_core.py            # shared cluster-bootstrap utilities
python scripts/reweight_ci.py        # data/reweight_ci_28.csv
python scripts/paired_diff_ci.py     # data/paired_diff_ci.csv
python scripts/seed_level_diff.py    # data/seed_level_diff.csv
python scripts/merge_item4.py        # data/paired_diff_merged.csv
python scripts/verify_table2.py      # data/table2_per_seed_verify.csv
```

The scripts expect the frozen condition-level scores that accompany the article. They are
CPU-only and run in a few minutes; no GPU is required.

## Intervals used in the paper

- **Main intervals.** Cluster resampling with probabilities proportional to unit size and
  Horvitz--Thompson frame weights, shared across the three seeds and paired with the RGB
  reference; 1,000 replicates, 2.5 and 97.5 percentiles.
- **Sensitivity intervals.** A macro-over-clusters estimator, reported to show that a
  different estimand does not remove the small-cluster limitation (six FLIR sequence units,
  three M3FD scene clusters).
- **Paired differences.** Frame-level paired bootstrap and a seed-level paired t interval
  (df = 2). A comparison is called robust only when both exclude zero with the same sign.

## License

Code under the MIT license (see `LICENSE`). Result tables under CC BY 4.0.

## Contact

Dongyuan Ge, gordon399@gxust.edu.cn
College of Mechanical and Automotive Engineering, Guangxi University of Science and Technology
