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
