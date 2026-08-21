# NSEI Sablefish MSE Analyses

Management strategy evaluation (MSE) for Northern Southeast Inside (NSEI) sablefish to identify biomass reference points using closed-loop simulation.

## Overview

This analysis uses [SPoRC](https://github.com/chengmatt/SPoRC/tree/dev-popn-seasons) as the operating model (OM) and estimation model (EM) framework to evaluate harvest control rules and biomass reference points for NSEI sablefish through closed-loop feedback simulation.

## Installation

```r
# install.packages("devtools")
devtools::install_github("chengmatt/SPoRC@dev-popn-seasons")
```

## Workflow

1. **Condition OM** on historical data (years prior to `feedback_start_yr`)
2. **Run closed-loop feedback** with annual assessment fitting via the EM
3. **Derive reference points** and project catch advice each assessment year
4. **Evaluate performance** of candidate biomass reference points across simulations

## General Simulation design

Three operating models, each 150 replicates x 111 years, with `feedback_start_yr = 51`.
Years 1-50 are the conditioning period (1977-2026, a single deterministic trajectory
shared by every replicate); years 51-111 are the projection (2027-2087).

| OM | File | Recruitment |
|---|---|---|
| Baseline | `outputs/base_mse_om.RDS` | mean recruitment |
| Regime | `outputs/bh_regime_mse_om.RDS` | Beverton-Holt with regime shifts |
| Crash | `outputs/bh_crash_mse_om.RDS` | Beverton-Holt with a recruitment crash |

Each OM was run twice, under two information assumptions:

| Mode | Script | HCR grid | Result files |
|---|---|---|---|
| **Perfect information** (no assessment; HCR sees true stock state) | `R/run_mse_loop_noasmt.r` | full: SPR 0.30-0.70 x alpha 0-0.50 x 3 stability options + constant F50 = 298 | `outputs/*_results_noasmt.RDS` (~7 GB each) |
| **With assessment** (EM fit every year; HCR sees estimates) | `R/run_mse_loop_asmt.R` | reduced: SPR {0.40, 0.50} x alpha {0, 0.25, 0.50}, asymmetric only + constant F50 = 7 | `outputs/*_results_asmt.RDS` (~200 MB each) |

Both run sets share operating models, seeds and recruitment deviations, so replicate *j*
is the same underlying future in both. Conditioning-period SSB is bit-identical between
them, which is what makes the comparison attributable to assessment error alone.

**The two grids only overlap on 7 rules** (constant F50, and F40/F50 x alpha {0, 0.25, 0.5},
all asymmetric). Every assessment-vs-no-assessment comparison is restricted to those.
The SPR x alpha heatmaps exist for the no-assessment runs only -- the assessment grid is
too coarse for them.

## Scripts

Run the following scripts in this order:

| Script | Purpose |
|---|---|
| `R/condition_om.R` | Condition the OM on historical data, write the `*_mse_om.RDS` files |
| `R/run_mse_loop_noasmt.r` | Closed loop with perfect information, full HCR grid |
| `R/run_mse_loop_asmt.R` | Closed loop with an EM fit each year, reduced HCR grid |
| `R/extract_mse_summaries.R` | Read all six result files once (~2 min), write `outputs/mse_summaries.RDS` (~35 MB) |
| `R/functions/mse_viz_utils.R` | Shared metric definitions, year mapping, palettes, plotting helpers |
| `R/visualize_mse_results_noasmt.R` | No-assessment figures, including the SPR x alpha heatmaps (reads the raw 7 GB files directly) |
| `R/visualize_mse_results_asmt.R` | Assessment-run figures + EM diagnostics -> `figs/asmt_*.png` |
| `R/compare_asmt_noasmt.R` | Assessment vs perfect information -> `figs/cmp_*.png`, `outputs/asmt_vs_noasmt_performance.csv` |
| `R/visualize_hcr_shapes.R` | Illustrative HCR shapes |
| `R/visualize_mse_params_conditoning.R` | OM conditioning fits and demographics |
| `R/southeast_sabie_prices.R` | Price model used for economic metrics |

Re-run `R/extract_mse_summaries.R` whenever the MSE runs change.

## Figures

`figs/` prefixes:

- **no prefix** -- no-assessment runs (`grid_*`, `hcr_comparison_*`, `base_hcr_*`, `om_*`, `hcr_shape_*`)
- **`asmt_`** -- assessment runs
- **`cmp_`** -- assessment vs perfect information

Assessment-run figures (`asmt_*`):

| Figure | Shows |
|---|---|
| `recruitment_timeseries` | Age-2 recruitment by OM |
| `base_hcr_{ssb,cat,f,status}_timeseries` | Reference rule (constant F50), median + 90% interval, terminal-year density |
| `hcr_comparison_{ssb,cat,f,status}_timeseries[_zoom]` | All 7 rules overlaid |
| `em_ssb_relative_error` | Relative error of the terminal SSB estimate against the OM -- the bias the rule acts on |
| `em_perceived_vs_true_status` | SSB/Bx as the rule saw it vs as it truly was |
| `em_convergence` | Positive-definite Hessian rate and worst gradient by year |
| `retro_fan` | Retrospective fan (squid plot) for one replicate against OM truth, one panel per rule x OM |

Comparison figures (`cmp_*`):

| Figure | Shows |
|---|---|
| `base_hcr_{ssb,cat,f,status}_timeseries` | Reference rule, both modes overlaid with intervals |
| `all_hcr_{ssb,cat,f,status}_timeseries` | Same, one row per rule (7 x 3 panels) |
| `performance_metrics` | Median SSB, median catch, catch AAV, P(SSB < 10% historical max), P(closure), both modes |
| `performance_change` | The same metrics as change attributable to assessment error, coloured better/worse |
| `tradeoff` | % change in catch vs % change in AAV, shared axes, 50% simulation intervals from paired replicates |
| `tradeoff_medians` | Same axes, medians only, labelled, free scales |
