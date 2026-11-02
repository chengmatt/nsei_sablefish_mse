# Shared helpers for make_figures.R
#
# outputs/mse_summaries.RDS holds one entry per OM scenario and mode, each with
# one record per harvest control rule. Inside a record every quantity is a
# matrix with rows = model year and cols = replicate. The conditioning period
# runs from 1975 to the terminal conditioning year; from then on the rule is in
# control, and the quantities the rule produced (bx, ssb_em_term) have rows for
# those projection years only.

library(here)
library(ggplot2)
library(dplyr)

mse_summaries <- readRDS(here("outputs", "mse_summaries.RDS"))

terminal_conditioning_year <- mse_summaries$meta$term_cond_yr  # also the first assessment year
hcr_names         <- names(mse_summaries$meta$hcr_set)         # reference rule (constant F50) first
om_scenario_names <- c("Baseline", "Regime", "Crash")
reference_hcr     <- hcr_names[1]

# mode = the information the harvest control rule had to work with
mode_labels <- c(noasmt = "Perfect information", asmt = "With assessment")

# metrics that only exist once the rule is in control
projection_only_metrics <- c("status", "status_perceived", "ssb_re")

# the reference rule is drawn in black with a grey interval; the other six rules
# take the colourblind-safe palette (black and orange dropped)
hcr_colours      <- c("black",  ggthemes::colorblind_pal()(8)[-c(1, 2)])
hcr_fill_colours <- c("grey50", ggthemes::colorblind_pal()(8)[-c(1, 2)])
mode_colours     <- setNames(c("#0072B2", "#D55E00"), unname(mode_labels))


# One run ------------------------------------------------------------------

# Everything one OM scenario, mode and rule needs, under plain names. The OM
# carries biomass and catch in kg and recruitment in numbers, while every
# figure reports tonnes and millions of fish, so convert here and nowhere else.
get_run <- function(om_scenario, mode, hcr) {
  scenario_entry <- mse_summaries[[paste(om_scenario, mode, sep = "_")]]
  record <- scenario_entry$hcrs[[hcr]]

  first_projection_row <- scenario_entry$fb
  n_years <- nrow(record$SSB)

  list(
    ssb               = record$SSB   / 1e3,   # t
    catch             = record$Catch / 1e3,   # t
    fishing_mortality = record$Fmort,
    recruits          = record$Rec   / 1e6,   # millions of age-2 fish

    # projection years only, in the same units as ssb
    target_ssb   = record$bx / 1e3,   # the Bx the rule aimed at that year
    assessed_ssb = if (is.null(record$ssb_em_term)) NULL else record$ssb_em_term / 1e3,

    conditioning_rows = seq_len(first_projection_row - 1),
    projection_rows   = first_projection_row:n_years,

    # row 1 is this calendar year, so the rows run 1975 ... 2085
    years = seq(terminal_conditioning_year - first_projection_row + 1,
                length.out = n_years)
  )
}

# Depletion threshold: 10% of each replicate's conditioning-period maximum SSB
crash_threshold <- function(run) {
  apply(run$ssb[run$conditioning_rows, , drop = FALSE], 2, max) * 0.1
}


# Metrics ------------------------------------------------------------------

# One metric as a replicate matrix, together with the calendar year of each row.
#   SSB, Catch, Fmort, Rec  straight out of the OM, whole series
#   status                  true SSB / the target the rule used that year
#   status_perceived        SSB the assessment reported / that same target
#   ssb_re                  relative error of the assessment's terminal SSB
get_metric <- function(run, metric) {
  true_ssb <- run$ssb[run$projection_rows, , drop = FALSE]

  values <- switch(metric,
    SSB   = run$ssb,
    Catch = run$catch,
    Fmort = run$fishing_mortality,
    Rec   = run$recruits,
    status           = true_ssb / run$target_ssb,
    status_perceived = run$assessed_ssb / run$target_ssb,
    ssb_re           = (run$assessed_ssb - true_ssb) / true_ssb
  )

  years <- if (metric %in% projection_only_metrics) {
    run$years[run$projection_rows]
  } else {
    run$years
  }
  list(values = values, years = years)
}

# Median and 90% interval across replicates, one row per year. A year whose
# replicates all failed to converge stays NA.
summarize_replicates <- function(values, years) {
  yearly <- data.frame(year = years, median = NA, lower = NA, upper = NA)
  for (i in seq_len(nrow(values))) {
    replicates <- values[i, ]
    replicates <- replicates[!is.na(replicates)]
    if (length(replicates) > 0) {
      yearly[i, c("median", "lower", "upper")] <-
        unname(quantile(replicates, c(0.5, 0.05, 0.95)))
    }
  }
  yearly
}

# Factor levels fix the facet, legend and panel order in every figure
set_plot_order <- function(summary_table, hcrs, modes) {
  summary_table$om_scenario <- factor(summary_table$om_scenario, om_scenario_names)
  summary_table$hcr         <- factor(summary_table$hcr, hcrs)
  summary_table$mode        <- factor(summary_table$mode, unname(mode_labels[modes]))
  summary_table
}


# Long frames --------------------------------------------------------------

# Median and interval by year for one metric: one row per OM scenario, mode,
# rule and year.
build_timeseries <- function(metric, hcrs = hcr_names, modes = names(mode_labels),
                             om_scenarios = om_scenario_names) {
  all_rows <- list()
  for (om_scenario in om_scenarios) {
    for (mode in modes) {
      for (hcr in hcrs) {
        run <- get_run(om_scenario, mode, hcr)
        metric_data <- get_metric(run, metric)

        all_rows[[paste(om_scenario, mode, hcr)]] <- data.frame(
          om_scenario = om_scenario,
          mode        = unname(mode_labels[mode]),
          hcr         = hcr,
          summarize_replicates(metric_data$values, metric_data$years)
        )
      }
    }
  }
  set_plot_order(bind_rows(all_rows), hcrs, modes)
}

# Catch AAV for one replicate: year-on-year change over total catch
catch_aav_one_replicate <- function(catch) {
  if (sum(head(catch, -1)) == 0) return(NA)
  sum(abs(diff(catch))) / sum(head(catch, -1))
}

# One row per OM scenario, mode and rule, all scored over the projection
# period. The column names are the ones grid_metrics.csv already uses.
#   med_ssb, med_catch  median over every projection year and replicate
#   catch_aav           the median replicate's catch variability
#   p_crash             share of replicate-years under the depletion threshold
#   p_closure           share of replicate-years with no catch at all
#   term_status         median SSB/Bx in the final year
performance_metrics <- function(hcrs = hcr_names, modes = names(mode_labels),
                                om_scenarios = om_scenario_names) {
  all_rows <- list()
  for (om_scenario in om_scenarios) {
    for (mode in modes) {
      for (hcr in hcrs) {
        run              <- get_run(om_scenario, mode, hcr)
        projection_ssb   <- run$ssb[run$projection_rows, , drop = FALSE]
        projection_catch <- run$catch[run$projection_rows, , drop = FALSE]
        status           <- get_metric(run, "status")$values

        all_rows[[paste(om_scenario, mode, hcr)]] <- data.frame(
          om_scenario = om_scenario,
          mode        = unname(mode_labels[mode]),
          hcr         = hcr,
          med_ssb     = median(projection_ssb),
          med_catch   = median(projection_catch),
          catch_aav   = median(apply(projection_catch, 2, catch_aav_one_replicate),
                               na.rm = TRUE),
          p_crash     = mean(sweep(projection_ssb, 2, crash_threshold(run), "<")),
          p_closure   = mean(projection_catch == 0),
          term_status = median(status[nrow(status), ], na.rm = TRUE)
        )
      }
    }
  }
  set_plot_order(bind_rows(all_rows), hcrs, modes)
}
