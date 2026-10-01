# Computation chain for the paper

*A multi-scale framework for reliable surface–function correlation.*

Every number, table and figure of the paper derived from simulation can be
regenerated from this folder. Each file header carries its version history.

## Requirements

- MATLAB R2025b (earlier releases from R2022a are expected to work but were not tested).
- Statistics and Machine Learning Toolbox (`fitlme`, `anova1`, `kruskalwallis`).
- Parallel Computing Toolbox, optional. Without it, `parfor` runs serially and
  the results are identical, only slower.

## Quick start

```matlab
run_all             % everything cited in the paper
run_all('court')    % short studies only (a few minutes)
run_all('long')     % Monte Carlo studies only
run_all('figures')  % redraw the figures from existing ws_*.mat files
run_all('sec4')     % only what feeds Section 4; likewise 'sec3', 'sec5', 'sec6'
run_all('appB')     % only what feeds Appendix B; likewise 'appC', 'appE'
```

Each step saves its workspace in `ws_<script>.mat` and the console is logged in
`journal_<date>.txt`. A failed step does not stop the chain; the final summary
lists what did not complete.

## Reproducibility

All random draws go through `paper_seed(scenario, replicate, level)`, which
returns a fixed seed for a `RandStream('Threefry', ...)`. The same scenario,
level and replicate always give the same data, whatever the order of execution
or the number of parallel workers. Seed blocks are listed below and in the
header of `paper_seed.m`. The reference draw of the worked example uses the
fixed seed 40007.

## Core

| File | Role |
|---|---|
| `paper_params.m` | Parameters and default thresholds: G = 5, n = 5, R = 5, S = 10, K = 5; q = 0.05, η² = 0.06, ICC = 0.60, R² = 0.30, C = 0.90. |
| `gen_volet1.m` | Synthetic data generator, all scenarios (`nominal`, `null`, `twoscale`, `unreliable`, `nonnormal`, `heteroscedastic`, `nonlinear`, `autocorr`, `curvature`, `randslope`, `ushape`, `fragile`). |
| `gen_fragile.m` | Anomaly injected on a single specimen. |
| `paper_pipeline.m` | The three gates, Spearman's ρ, bootstrap, leave-one-specimen-out check. |
| `attenuation_correct.m` | Attenuation correction of the slope (Section 3.9.4). |
| `paper_seed.m` | Seed allocation. |
| `paper_figstyle.m` | Common figure format (15 cm, 9 pt). |
| `swtest.m` | Shapiro–Wilk test, called by `paper_pipeline`. |

## Paper item → script

| Paper item | Section | Script | Output | Seed block |
|---|---|---|---|---|
| Descriptor design, ground truth (tables) | 4 | `paper_params`, `gen_volet1` | — | — |
| Discrimination and correlation curves and matrices (figures) | 4 | `run_volet1` → `make_results_figures` | `Fig_curves_*`, `Fig_matrices_*` | 40007 |
| Per-descriptor results, retained pairs (table) | 4 | `run_volet1`, `make_bootstrap_table(out)` | console | 40007 |
| Decision matrix (figure) | 4 | `run_volet1` | `Fig_decision_matrix.png` | 40007 |
| Bootstrap against model-based standard error (figure) | 4 | `run_boot_vs_param` | `Fig_boot_vs_param.png` | 0 |
| Coverage and stability of the slope interval (text) | 4 | `run_bca_coverage`, `run_bca_stability` | `ws_*.mat` | 0 |
| Branch coverage under assumption violations (table) | 4 | `run_robustness_compare` | console | 0 |
| Branching against permutation (two tables) | 4 | `run_perm_compare` | `perm_compare_results.csv` | 0 |
| Stability of the permutation level at n = 5 (text) | 4 | `run_perm_seeds` | `ws_*.mat` | 120–139 |
| Fragile flag, calibration (figure, text) | 4 | `run_fragile`, `run_loso_sweep` | `Fig_loso_sweep.png` | 0, 80–99 |
| Correlated but unreliable descriptor (figure, text) | 4 | `run_unreliable`, `make_unreliable_figure` | `Fig_unreliable_curves.png` | 40007, 0 |
| Null and two-scale scenarios (figure) | 4 | `run_scenarios` | `Fig_decision_matrices_BC.png` | 0 |
| Attenuation correction, configurations (b) and (c) (text) | 3, 5 | `run_attenuation`, `run_attenuation_c` | `ws_*.mat` | 260–269, 440–459 |
| Bias-correction parameter when no descriptor is exact (text) | 4 | `run_latent` | `ws_*.mat` | 270–279 |
| Stress tests (table) | 5 | `run_stress_tests` | `ws_*.mat` | 1–19 |
| Sample size and permutation (figure) | 5 | `run_perm_size` | `Fig_perm_size.png` | 100–119 |
| Curvature against noise (figure) | 5 | `run_curvature_noise` | `Fig_curvature_noise.png` | 140–199 |
| Response with an interior optimum (table) | 5 | `run_ushape` | `ws_*.mat` | 400–439 |
| Design study, 81 configurations (two tables, three figures) | 5, C | `run_volet2_pipeline` → `make_sensitivity_figures` | `volet2_results_paper.csv`, `Fig_sensitivity_*.png` | 300–399 |
| Contribution of the upstream gates (text) | 6 | `run_value_added` | `value_added_results.csv` | 0 |
| Chance R² with no descriptor–response link (text) | 6 | `run_r2_null` | `r2_null.csv` | 510–519 |
| Composite indicator (figure) | 6 | `run_volet1` | `Fig_Q_indicator.png` | 40007 |
| Aggregation example (table) | 6 | illustrative values, no script | — | — |
| ICC against contrast and number of conditions (figure) | 6 | `run_contrast_sweep` | `Fig_contrast_sweep.png` | 20–59 |
| Simplified path, R² level on specimen means | B | `run_r2_means_calib` | `r2_means_calib.csv`, `Fig_r2_means_calib.png` | 460–504 |
| Random-slope selection (text, table, figure) | E | `run_randslope`, `run_randslope_R` → `plot_randslope_R` | `Fig_randslope_R.png` | 60–79, 230–249 |
| Flowcharts | 3, B | TikZ in the manuscript | — | — |

`make_all_figures` redraws every figure from the saved `ws_*.mat` files in a few
seconds, without recomputation.

## License

This code is distributed under the GNU General Public License v3.0 (see
`LICENSE`). Integration into a closed-source distributed product requires a
separate agreement with the author.
