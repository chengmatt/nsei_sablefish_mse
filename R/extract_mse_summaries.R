# Purpose: Condense the raw MSE result files (assessment and no-assessment)
# Creator: Matthew LH. Cheng

library(here)

terminal_conditioning_year <- 2025

# HCRs
hcr_set <- list(
  "Constant F50" = list(spr_x = 0.5, alpha = NA,   stability = "asymmetric", hcr_type = "constant"),
  "F40, a=0"     = list(spr_x = 0.4, alpha = 0,    stability = "asymmetric", hcr_type = "threshold"),
  "F40, a=0.25"  = list(spr_x = 0.4, alpha = 0.25, stability = "asymmetric", hcr_type = "threshold"),
  "F40, a=0.5"   = list(spr_x = 0.4, alpha = 0.5,  stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0"     = list(spr_x = 0.5, alpha = 0,    stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.25"  = list(spr_x = 0.5, alpha = 0.25, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.5"   = list(spr_x = 0.5, alpha = 0.5,  stability = "asymmetric", hcr_type = "threshold")
)

# the grids are built with seq(), so alpha/spr_x carry floating-point noise --
# match on tolerance rather than identity
match_hcr <- function(scenario, hcr_definition) {
  scenario_type <- if (is.null(scenario$hcr_type)) "threshold" else scenario$hcr_type
  if (scenario_type != hcr_definition$hcr_type) return(FALSE)
  if (as.character(scenario$stability) != hcr_definition$stability) return(FALSE)
  if (!isTRUE(all.equal(scenario$spr_x, hcr_definition$spr_x))) return(FALSE)
  # constant F has no ramp, so alpha is not part of its identity
  if (hcr_definition$hcr_type == "constant") return(TRUE)
  isTRUE(all.equal(scenario$alpha, hcr_definition$alpha))
}

# Position of the first scenario in a results file that is this HCR, else NA
find_hcr <- function(all_scenarios, hcr_definition) {
  for (i in seq_along(all_scenarios)) {
    if (match_hcr(all_scenarios[[i]]$scenario, hcr_definition)) return(i)
  }
  NA
}

# Stack one quantity across replicates: rows = year, cols = replicate
replicate_matrix <- function(replicates, pull_value) {
  columns <- vector("list", length(replicates))
  for (i in seq_along(replicates)) {
    columns[[i]] <- as.vector(pull_value(replicates[[i]]))
  }
  matrix(unlist(columns), nrow = length(columns[[1]]), ncol = length(columns))
}

# SSB the assessment reported for its own terminal year, one value per
# assessment year. This is the quantity the HCR actually saw, so it carries the
# assessment-error signal. Years where the fit failed stay NA.
terminal_assessed_ssb <- function(replicate) {
  terminal <- rep(NA, length(replicate$ssb_em))
  for (i in seq_along(replicate$ssb_em)) {
    series <- replicate$ssb_em[[i]]
    if (!is.null(series)) terminal[i] <- series[length(series)]
  }
  terminal
}

# Compact per HCR record per entry
summarize_hcr <- function(scenario_entry, mode) {
  replicates <- scenario_entry$results
  record <- list(
    scenario_id = scenario_entry$scenario$scenario_id,
    SSB   = replicate_matrix(replicates, function(replicate) replicate$om$SSB),
    Catch = replicate_matrix(replicates, function(replicate) replicate$om$Catch),
    Fmort = replicate_matrix(replicates, function(replicate) replicate$om$Fmort),
    Rec   = replicate_matrix(replicates, function(replicate) replicate$om$Rec),
    bx    = replicate_matrix(replicates, function(replicate) replicate$bx)
  )
  if (mode == "asmt") {
    record$ssb_em_term <- replicate_matrix(replicates, terminal_assessed_ssb)
    # convergence diagnostics, kept for reporting rather than for a figure
    record$pdHess <- replicate_matrix(replicates, function(replicate) replicate$pdHess)
    record$grad   <- replicate_matrix(replicates, function(replicate) replicate$grad)
  }
  record
}

# helper to help extract file
extract_file <- function(path, om_scenario, mode, feedback_start_yr) {
  all_scenarios <- readRDS(path)
  hcrs <- vector("list", length(hcr_set))
  names(hcrs) <- names(hcr_set)
  for (hcr_name in names(hcr_set)) {
    match_index <- find_hcr(all_scenarios, hcr_set[[hcr_name]])
    if (is.na(match_index)) {
      cat("  HCR not found:", hcr_name, "\n")
      next
    }
    hcrs[[hcr_name]] <- summarize_hcr(all_scenarios[[match_index]], mode)
  }

  rm(all_scenarios); gc(verbose = FALSE)   # the raw files are large so remove garbo
  list(om = om_scenario, mode = mode, fb = feedback_start_yr, hcrs = hcrs)
}


# Run ----------------------------------------------------------------------

om_files <- list(
  Baseline = "base_mse_om.RDS",
  Regime   = "bh_regime_mse_om.RDS",
  Crash    = "bh_crash_mse_om.RDS"
)
results_files <- list(
  Baseline = c(asmt = "base_results_asmt.RDS",   noasmt = "base_results_noasmt.RDS"),
  Regime   = c(asmt = "regime_results_asmt.RDS", noasmt = "regime_results_noasmt.RDS"),
  Crash    = c(asmt = "crash_results_asmt.RDS",  noasmt = "crash_results_noasmt.RDS")
)

summaries <- list()
for (om_scenario in names(results_files)) {
  feedback_start_yr <- readRDS(here("outputs", om_files[[om_scenario]]))$feedback_start_yr
  for (mode in c("asmt", "noasmt")) {
    cache_key <- paste(om_scenario, mode, sep = "_")
    summaries[[cache_key]] <- extract_file(here("outputs", results_files[[om_scenario]][[mode]]),
                                           om_scenario, mode, feedback_start_yr)
  }
}

summaries$meta <- list(term_cond_yr = terminal_conditioning_year, hcr_set = hcr_set)
saveRDS(summaries, here("outputs", "mse_summaries.RDS"))
