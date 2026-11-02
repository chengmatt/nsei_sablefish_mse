# Purpose: Performance metrics for every screening-grid scenario (no assessment)
# Creator: Matthew LH. Cheng
# Same metric definitions as performance_metrics() in R/functions/mse_viz_utils.R.

library(here)
library(dplyr)

om_files <- list(
  Baseline = list(results = "base_results_noasmt.RDS",   om = "base_mse_om.RDS"),
  Regime   = list(results = "regime_results_noasmt.RDS", om = "bh_regime_mse_om.RDS"),
  Crash    = list(results = "crash_results_noasmt.RDS",  om = "bh_crash_mse_om.RDS")
)

# Stack one quantity across replicates: rows = year, cols = replicate
replicate_matrix <- function(replicates, pull_value) {
  columns <- vector("list", length(replicates))
  for (i in seq_along(replicates)) {
    columns[[i]] <- as.vector(pull_value(replicates[[i]]))
  }
  matrix(unlist(columns), nrow = length(columns[[1]]), ncol = length(columns))
}

# Catch AAV for one replicate: year-on-year change over total catch
catch_aav_one_replicate <- function(catch) {
  if (sum(head(catch, -1)) == 0) return(NA)
  sum(abs(diff(catch))) / sum(head(catch, -1))
}

all_rows <- list()
for (om_scenario in names(om_files)) {
  cat("reading", om_files[[om_scenario]]$results, "\n")
  scenarios <- readRDS(here("outputs", om_files[[om_scenario]]$results))
  feedback_start_yr <- readRDS(here("outputs", om_files[[om_scenario]]$om))$feedback_start_yr

  for (scenario_entry in scenarios) {
    scenario <- scenario_entry$scenario
    ssb   <- replicate_matrix(scenario_entry$results, function(replicate) replicate$om$SSB)
    catch <- replicate_matrix(scenario_entry$results, function(replicate) replicate$om$Catch)

    # everything below is scored over the projection period only
    projection_rows  <- feedback_start_yr:nrow(ssb)
    projection_ssb   <- ssb[projection_rows, , drop = FALSE]
    projection_catch <- catch[projection_rows, , drop = FALSE]

    # depletion threshold: 10% of each replicate's conditioning-period max
    crash_threshold <- apply(ssb[seq_len(feedback_start_yr - 1), , drop = FALSE], 2, max) * 0.1

    all_rows[[length(all_rows) + 1]] <- data.frame(
      om_scenario = om_scenario,
      scenario_id = scenario$scenario_id,
      hcr_type    = scenario$hcr_type,
      spr_x       = scenario$spr_x,
      alpha       = if (is.null(scenario$alpha) || is.na(scenario$alpha)) NA else scenario$alpha,
      dyn_b0      = scenario$dyn_b0,
      stability   = scenario$stability,
      n_reps      = ncol(ssb),
      # median over all projection years and replicates (OM stores kg; report t)
      med_ssb      = median(projection_ssb)   / 1e3,
      med_catch    = median(projection_catch) / 1e3,
      # median of the annual medians, so each year weighs the same
      med_ssb_yr   = median(apply(projection_ssb,   1, median)) / 1e3,
      med_catch_yr = median(apply(projection_catch, 1, median)) / 1e3,
      catch_aav    = median(apply(projection_catch, 2, catch_aav_one_replicate), na.rm = TRUE),
      # fraction of replicate-years below the threshold / with no catch
      p_crash      = mean(sweep(projection_ssb, 2, crash_threshold, "<")),
      p_closure    = mean(projection_catch == 0),
      # fraction of replicates that close at least once
      p_rep_any_closure = mean(apply(projection_catch == 0, 2, any)),
      min_ssb_med  = median(apply(projection_ssb, 2, min)) / 1e3,
      stringsAsFactors = FALSE
    )
  }

  rm(scenarios); gc(verbose = FALSE)   # the raw files are large, so free each one
}

grid_metrics <- bind_rows(all_rows)
write.csv(grid_metrics, here("outputs", "grid_metrics.csv"), row.names = FALSE)
cat("wrote outputs/grid_metrics.csv,", nrow(grid_metrics), "rows\n")
