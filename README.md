# NSEI Sablefish MSE Analyses

Management strategy evaluation for Northern Southeast Inside (NSEI) sablefish, using
closed-loop simulation to identify biomass reference points.
[SPoRC](https://github.com/chengmatt/SPoRC/tree/dev-popn-seasons) is both the operating
model and the estimation model.

## Installation

```r
# install.packages("devtools")
devtools::install_github("chengmatt/SPoRC@dev-popn-seasons")
```

## Design

Three operating models, each 150 replicates x 111 years (1975-2085). Years 1-51
(1975-2025) are the conditioning period, a single deterministic trajectory shared by every
replicate. The first closed-loop assessment runs in year 51 and sets the 2026 catch, so
years 52-111 are HCR-controlled. The OM stores biomass in kg and the figures report tonnes.

| OM | File | Recruitment |
|---|---|---|
| Baseline | `outputs/base_mse_om.RDS` | resampled from the 1975-2025 SCAA estimates |
| Regime | `outputs/bh_regime_mse_om.RDS` | Beverton-Holt, 5-year 5x R0 booms every 20 years |
| Crash | `outputs/bh_crash_mse_om.RDS` | Beverton-Holt, R0 x 0.1 over 2040-2050 |

Each OM is run under two information assumptions, sharing operating models, seeds and
recruitment deviations, so replicate *j* is the same underlying future in both.

| Mode | Script | HCR grid | Results |
|---|---|---|---|
| Perfect information | `R/run_mse_loop_noasmt.r` | 298 rules | `outputs/*_results_noasmt.RDS` |
| With assessment | `R/run_mse_loop_asmt.R` | 7 rules | `outputs/*_results_asmt.RDS` |

The grids overlap on 7 rules (constant F50, and F40/F50 x alpha 0, 0.25 and 0.5, all
asymmetric). Every assessment-vs-no-assessment comparison is restricted to those, and the
SPR x alpha screening profiles come from the no-assessment runs only.

## Scripts

Run in this order.

| Script | Writes |
|---|---|
| `R/condition_om.R` | `outputs/*_mse_om.RDS` |
| `R/run_mse_loop_noasmt.r` | `outputs/*_results_noasmt.RDS` |
| `R/run_mse_loop_asmt.R` | `outputs/*_results_asmt.RDS` |
| `R/extract_mse_summaries.R` | `outputs/mse_summaries.RDS` |
| `R/extract_grid_metrics.R` | `outputs/grid_metrics.csv` |
| `R/make_figures.R` | everything in `figs/`, plus `outputs/asmt_vs_noasmt_performance.csv` |

Re-run both extract scripts whenever the MSE runs change. `make_figures.R` reads only the
two extract outputs and `outputs/base_mse_om.RDS`, and its sections are independent, so a
single figure can be re-made on its own.

`R/functions/` holds the helpers, `utils.R` and `setup_em.R` for the run loops and
`mse_viz_utils.R` for the figures. `old/` holds work no longer in the pipeline.
