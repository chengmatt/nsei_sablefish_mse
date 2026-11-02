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

Three operating models, each 150 replicates x 111 years (1975-2085), with `feedback_start_yr = 51`.
Years 1-51 are the conditioning period (1975-2025: the ADF&G SCAA's 51 model years, a single
deterministic trajectory shared by every replicate). The first closed-loop assessment is run in
year 51 (2025) and sets the 2026 catch, so years 52-111 (2026-2085) are HCR-controlled.
Performance summaries use model years 51-111 (`fb:n_yrs`), i.e. they include the shared 2025 row.
The OM stores biomass in kg; figures and `outputs/asmt_vs_noasmt_performance.csv` report SSB and
catch in t.

| OM | File | Recruitment |
|---|---|---|
| Baseline | `outputs/base_mse_om.RDS` | resampled with replacement from the 1975-2025 SCAA estimates (no stock-recruit relationship) |
| Regime | `outputs/bh_regime_mse_om.RDS` | Beverton-Holt (h = 0.8, sigmaR = 1) with 5-year 5x R0 booms every 20 years (2041-45, 2061-65, 2081-85) |
| Crash | `outputs/bh_crash_mse_om.RDS` | Beverton-Holt (h = 0.8, sigmaR = 1) with R0 x 0.1 over 2040-2050 |

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
The SPR x alpha screening profiles exist for the no-assessment runs only -- the assessment
grid is too coarse for them.

## Scripts

Run in this order:

| Script | Purpose |
|---|---|
| `R/condition_om.R` | Condition the OM on historical data, write the `*_mse_om.RDS` files |
| `R/run_mse_loop_noasmt.r` | Closed loop with perfect information, full HCR grid |
| `R/run_mse_loop_asmt.R` | Closed loop with an EM fit each year, reduced HCR grid |
| `R/extract_mse_summaries.R` | Read all six raw result files once (~2 min), write `outputs/mse_summaries.RDS` |
| `R/extract_grid_metrics.R` | Per-rule metrics for all 298 rules under perfect information (reads the raw 7 GB files once, ~1 min) -> `outputs/grid_metrics.csv` |
| `R/make_figures.R` | **Every figure in `figs/`**, plus `outputs/asmt_vs_noasmt_performance.csv` |

Helpers: `R/functions/utils.R` and `R/functions/setup_em.R` support the run loops;
`R/functions/mse_viz_utils.R` holds the metric definitions, year mapping and palettes
that `make_figures.R` draws on.

`R/make_figures.R` reads only `outputs/mse_summaries.RDS`, `outputs/grid_metrics.csv`
and `outputs/base_mse_om.RDS`, so it does not touch the raw result files. Its sections are
independent and labelled by figure, so a single figure can be re-made without running the
whole file. The one ordering constraint is internal: the performance summary near the top
feeds Figure 4 and the comparison figures.

Re-run `R/extract_mse_summaries.R` and `R/extract_grid_metrics.R` whenever the MSE runs change.

**Known break:** the last section of `make_figures.R` (`om_conditioning_fit`, S1) is the only
one that runs the OM forward, and the committed `outputs/base_mse_om.RDS` predates the
installed SPoRC, so `run_annual_cycle()` fails in `generate_initial_age_structure()`. It is
last in the file for that reason -- the other 18 figures are unaffected. Re-run
`R/condition_om.R` to rebuild the OM before re-making S1.

`old/` holds work no longer in the pipeline: the sablefish price model
(`southeast_sabie_prices.R` and its outputs — no economic performance metric was used in
the end) and a scratch explorer for the SCAA output.

## Figures

All 19 figures in `figs/` are written by `R/make_figures.R`. Main-text figures:

| File | Figure | Shows |
|---|---|---|
| `fig1_recruitment_scenarios` | 1 | Age-2 recruitment by OM (constant F50, perfect information) |
| `hcr_shape_both` | 2 | Threshold rule shapes across SPR targets and closure thresholds |
| `fig3_screening_profiles` | 3 | Screening metrics vs SPR target, one line per alpha, asymmetric constraint |
| `fig4_candidate_metrics` | 4 | The seven candidate rules, all metrics, perfect information vs with assessment |
| `fig5_candidate_timeseries` | 5 | SSB and catch time series of the candidate rules with assessment |

Supplementary figures are embedded with captions in `notes/supplementary_material_1.docx`
(edited in Word; a condensed version is `notes/supplementary_material_1_condensed.docx`),
in this order:

| File | Figure | Shows |
|---|---|---|
| `om_conditioning_fit` | S1 | OM dynamics over 1975-2025 against the SCAA estimates |
| `om_demographics` | S2 | Weight-, maturity- and natural-mortality-at-age |
| `om_selectivity` | S3 | Fishery and survey selectivity and retention, terminal block |
| `asmt_em_ssb_relative_error` | S4 | Relative error of the terminal SSB estimate against the OM -- the bias the rule acts on |
| `asmt_em_perceived_vs_true_status` | S5 | SSB/Bx as the rule saw it vs as it truly was |
| `figS_screening_profiles_sym` | S6 | Screening profiles, symmetric stability constraint |
| `figS_screening_profiles_none` | S7 | Screening profiles, no stability constraint |
| `cmp_all_hcr_ssb_timeseries` | S8 | SSB, both modes overlaid, one row per rule (7 x 3 panels) |
| `cmp_all_hcr_cat_timeseries` | S9 | Catch, both modes overlaid |
| `cmp_all_hcr_status_timeseries` | S10 | SSB/Bx, both modes overlaid |
| `cmp_performance_change` | S11 | Median SSB, median catch, catch AAV, P(SSB < 10% historical max) and P(closure) as change attributable to assessment error, coloured better/worse |

One further figure is produced but is not in the manuscript or supplement:

| File | Shows |
|---|---|
| `cmp_all_hcr_f_timeseries` | Fishing mortality, both modes overlaid |
