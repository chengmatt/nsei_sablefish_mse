# Purpose: Shared helpers for plotting MSE results out of outputs/mse_summaries.RDS
#          (built by R/extract_mse_summaries.R). Keeps the assessment-run and
#          assessment-vs-perfect-information scripts working off one definition
#          of the metrics, year mapping and palettes.
# Creator: Matthew LH. Cheng

library(here)
library(ggplot2)
library(dplyr)

cache <- readRDS(here("outputs", "mse_summaries.RDS"))

term_cond_yr <- cache$meta$term_cond_yr   # terminal conditioning (assessment) year
hcr_levels   <- names(cache$meta$hcr_set) # reference rule (constant F50) first
om_levels    <- c("Baseline", "Regime", "Crash")
ref_hcr      <- hcr_levels[1]

# mode = information available to the harvest control rule
mode_labels <- c(noasmt = "Perfect information", asmt = "With assessment")
mode_levels <- unname(mode_labels)

# palettes: the reference rule is drawn in black with a grey interval; the other
# six rules take the colourblind-safe palette (black and orange dropped)
hcr_pal  <- c("black",  ggthemes::colorblind_pal()(8)[-c(1, 2)])
hcr_fill <- c("grey50", ggthemes::colorblind_pal()(8)[-c(1, 2)])
mode_pal <- setNames(c("#0072B2", "#D55E00"), mode_levels)

rec_of  <- function(om, mode) cache[[paste(om, mode, sep = "_")]]$hcrs
fb_of   <- function(om, mode) cache[[paste(om, mode, sep = "_")]]$fb
# calendar years for the full (conditioning + projection) series
years_of <- function(om, mode) {
  fb <- fb_of(om, mode); n_yr <- nrow(rec_of(om, mode)[[1]]$SSB)
  (term_cond_yr - (fb - 1) + 1):(term_cond_yr - (fb - 1) + n_yr)
}
proj_start_yr <- function(om, mode) term_cond_yr + 1

# Replicate matrix (rows = year, cols = replicate) for a metric. Derived metrics
# only exist over the projection period and are padded with NA so every matrix
# shares the same year rows.
#   status           true SSB / the Bx the rule used that year
#   status_perceived SSB the assessment reported / that same Bx -- what the rule saw
#   ssb_re           relative error of the assessment's terminal SSB vs the OM
metric_mat <- function(om, mode, hcr, metric) {
  r  <- rec_of(om, mode)[[hcr]]
  fb <- fb_of(om, mode)
  n_yr <- nrow(r$SSB)
  if (metric %in% c("SSB", "Catch", "Fmort", "Rec")) return(r[[metric]])

  pad <- function(m) rbind(matrix(NA_real_, fb - 1, ncol(m)), m)
  ssb_proj <- r$SSB[fb:n_yr, , drop = FALSE]

  switch(metric,
    status           = pad(ssb_proj / r$bx),
    status_perceived = pad(r$ssb_em_term / r$bx),
    ssb_em_term      = pad(r$ssb_em_term),
    ssb_re           = pad((r$ssb_em_term - ssb_proj) / ssb_proj),
    stop("unknown metric: ", metric)
  )
}

# median + 90% interval by year
summarize_mat <- function(mat, yrs) {
  q <- function(p) apply(mat, 1, function(v) {
    v <- v[!is.na(v)]; if (!length(v)) NA_real_ else unname(quantile(v, p))
  })
  data.frame(year = yrs, med = q(0.5), lwr = q(0.05), upr = q(0.95))
}

# Long time-series frame across OM scenarios / modes / HCRs
build_ts <- function(metric, scale = 1, hcrs = hcr_levels,
                     modes = names(mode_labels), oms = om_levels) {
  df <- do.call(rbind, lapply(oms, function(om)
    do.call(rbind, lapply(modes, function(md)
      do.call(rbind, lapply(hcrs, function(h) {
        m <- metric_mat(om, md, h, metric) / scale
        cbind(data.frame(om_scenario = om, mode = unname(mode_labels[md]), hcr = h),
              summarize_mat(m, years_of(om, md)),
              fb_year = proj_start_yr(om, md))
      }))))))
  finalize(df, hcrs, modes)
}

# Terminal projection-year replicate values (for the side density panels)
build_term <- function(metric, scale = 1, hcrs = hcr_levels,
                       modes = names(mode_labels), oms = om_levels) {
  df <- do.call(rbind, lapply(oms, function(om)
    do.call(rbind, lapply(modes, function(md)
      do.call(rbind, lapply(hcrs, function(h) {
        m <- metric_mat(om, md, h, metric) / scale
        data.frame(om_scenario = om, mode = unname(mode_labels[md]), hcr = h,
                   value = m[nrow(m), ])
      }))))))
  finalize(df, hcrs, modes)
}

finalize <- function(df, hcrs, modes) {
  df$om_scenario <- factor(df$om_scenario, levels = om_levels)
  df$hcr         <- factor(df$hcr, levels = hcrs)
  df$mode        <- factor(df$mode, levels = unname(mode_labels[modes]))
  df
}

# Per-replicate performance metrics over the projection period ---------------

# catch AAV: sum|dC| / sum C within a replicate
catch_aav_reps <- function(catch_mat, proj_rows) {
  apply(catch_mat[proj_rows, , drop = FALSE], 2, function(cc) {
    if (sum(head(cc, -1)) == 0) return(NA_real_)
    sum(abs(diff(cc))) / sum(head(cc, -1))
  })
}

# One row per OM x mode x HCR. med_ssb / med_catch are medians over all
# projection years and replicates; p_crash is the fraction of projection
# replicate-years with SSB below 10% of that replicate's conditioning-period max.
perf_metrics <- function(hcrs = hcr_levels, modes = names(mode_labels), oms = om_levels) {
  df <- do.call(rbind, lapply(oms, function(om)
    do.call(rbind, lapply(modes, function(md)
      do.call(rbind, lapply(hcrs, function(h) {
        r    <- rec_of(om, md)[[h]]
        fb   <- fb_of(om, md)
        proj <- fb:nrow(r$SSB)
        thr  <- apply(r$SSB[seq_len(fb - 1), , drop = FALSE], 2, max) * 0.1
        data.frame(
          om_scenario = om, mode = unname(mode_labels[md]), hcr = h,
          med_ssb     = median(r$SSB[proj, ])   / 1e5,
          med_catch   = median(r$Catch[proj, ]) / 1e5,
          catch_aav   = median(catch_aav_reps(r$Catch, proj), na.rm = TRUE),
          p_crash     = mean(sweep(r$SSB[proj, , drop = FALSE], 2, thr, "<")),
          p_closure   = mean(r$Catch[proj, ] == 0),
          term_status = median(metric_mat(om, md, h, "status")[nrow(r$SSB), ], na.rm = TRUE)
        )
      }))))))
  finalize(df, hcrs, modes)
}

metric_labels <- c(
  med_ssb     = "Median SSB (100k t)",
  med_catch   = "Median catch (100k t)",
  catch_aav   = "Catch AAV",
  p_crash     = "P(SSB < 10% of historical max)",
  p_closure   = "P(fishery closure)",
  term_status = expression(Terminal~SSB/B[x])
)

# Time-series plot ----------------------------------------------------------

# colour_by: "hcr" or "mode" (or NULL for a single series). The 90% interval is
# drawn only for the first level of the colour variable so overlapping ribbons
# stay readable; the conditioning period is shared and drawn once in black.
plot_ts <- function(df, ylab, colour_by = NULL, term_df = NULL, zoom = FALSE,
                    facet = "wrap", pal = NULL, fill_pal = NULL, ribbon_all = FALSE) {

  facet_vars <- c("om_scenario", if (facet == "grid") "hcr")

  # the conditioning period is a single deterministic trajectory shared by every
  # replicate, HCR and mode -- keep one copy per panel and draw it in black
  hist_df <- df %>%
    filter(year <= fb_year, !is.na(med)) %>%
    distinct(across(all_of(c(facet_vars, "year"))), .keep_all = TRUE)
  proj_df <- df %>% filter(year >= fb_year, !is.na(med))
  if (!is.null(term_df)) term_df <- term_df %>% filter(!is.na(value))

  p <- ggplot(mapping = aes(year, med))

  if (!zoom && dplyr::n_distinct(hist_df$year) > 1) {
    p <- p +
      geom_vline(data = distinct(df, om_scenario, fb_year), aes(xintercept = fb_year),
                 linetype = 2, colour = "grey50") +
      geom_ribbon(data = hist_df, aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.20) +
      geom_line(data = hist_df, linewidth = 0.9)
  } else {
    p <- p + geom_vline(data = distinct(df, om_scenario, fb_year), aes(xintercept = fb_year),
                        linetype = 2, colour = "grey50")
  }

  if (!is.null(colour_by)) {
    # with many overlapping series only the first level gets an interval, so the
    # ribbons stay readable; ribbon_all is for two-way (mode) comparisons
    base_lvl <- levels(proj_df[[colour_by]])[1]
    ribbon_df <- if (ribbon_all) proj_df else proj_df[proj_df[[colour_by]] == base_lvl, ]
    p <- p +
      geom_ribbon(data = ribbon_df,
                  aes(ymin = lwr, ymax = upr, fill = .data[[colour_by]],
                      group = interaction(hcr, mode)),
                  colour = NA, alpha = if (ribbon_all) 0.15 else 0.20) +
      geom_line(data = proj_df,
                aes(colour = .data[[colour_by]], group = interaction(hcr, mode)),
                linewidth = 0.9)
    if (!is.null(pal)) p <- p + scale_colour_manual(values = pal) +
        scale_fill_manual(values = if (is.null(fill_pal)) pal else fill_pal)
  } else {
    p <- p +
      geom_ribbon(data = proj_df, aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.20) +
      geom_line(data = proj_df, linewidth = 0.9)
  }

  if (!is.null(term_df)) {
    dens_aes <- if (is.null(colour_by)) {
      aes(y = value, x = after_stat(scaled))
    } else {
      aes(y = value, x = after_stat(scaled),
          colour = .data[[colour_by]], fill = .data[[colour_by]])
    }
    p <- p +
      ggside::geom_ysidedensity(data = term_df, dens_aes, alpha = 0.20, linewidth = 0.6) +
      ggside::scale_ysidex_continuous(breaks = NULL)
  }

  p <- p +
    (if (facet == "grid") facet_grid(hcr ~ om_scenario, scales = "free_y")
     else facet_wrap(~ om_scenario, scales = "free_y")) +
    labs(x = "Year", y = ylab, colour = NULL, fill = NULL) +
    coord_cartesian(ylim = if (zoom) NULL else c(0, NA)) +
    theme_bw(base_size = 20) +
    theme(legend.position = if (is.null(colour_by)) "none" else "top")

  # ggside theme bits must come after theme_bw() (complete themes reset them)
  if (!is.null(term_df)) p <- p + theme(ggside.panel.scale.x = 0.2)

  p
}

# Paired per-replicate comparison --------------------------------------------

# Projection-period summary for every replicate individually.
rep_perf <- function(om, mode, hcr) {
  r    <- rec_of(om, mode)[[hcr]]
  fb   <- fb_of(om, mode)
  proj <- fb:nrow(r$SSB)
  data.frame(
    med_ssb   = apply(r$SSB[proj, , drop = FALSE],   2, median) / 1e5,
    med_catch = apply(r$Catch[proj, , drop = FALSE], 2, median) / 1e5,
    catch_aav = catch_aav_reps(r$Catch, proj)
  )
}

# The two run sets share operating model, seeds and recruitment deviations, so
# replicate j is the same underlying future in both. Differencing replicate by
# replicate removes that shared variability, leaving the effect of assessment
# error plus a simulation interval for it.
# Replicates whose perfect-information value is zero (fully closed fishery) give
# no defined percent change and are dropped; n_used reports how many remain.
# probs defaults to the interquartile range: the 5-95% range of a per-replicate
# ratio has a very heavy tail under the Crash OM (single replicates swing from
# -97% to +780% when a threshold rule opens in one run and closes in the other),
# which is real but swamps the median effect on any shared axis.
paired_change <- function(vars = c("med_catch", "catch_aav", "med_ssb"),
                          hcrs = hcr_levels, oms = om_levels,
                          probs = c(0.25, 0.75)) {
  df <- do.call(rbind, lapply(oms, function(om)
    do.call(rbind, lapply(hcrs, function(h) {
      a <- rep_perf(om, "asmt", h)
      n <- rep_perf(om, "noasmt", h)
      do.call(rbind, lapply(vars, function(v) {
        ch <- 100 * (a[[v]] - n[[v]]) / n[[v]]
        ch[!is.finite(ch)] <- NA_real_
        data.frame(om_scenario = om, hcr = h, var = v,
                   med = median(ch, na.rm = TRUE),
                   lwr = unname(quantile(ch, probs[1], na.rm = TRUE)),
                   upr = unname(quantile(ch, probs[2], na.rm = TRUE)),
                   n_used = sum(!is.na(ch)))
      }))
    }))))
  df$om_scenario <- factor(df$om_scenario, levels = om_levels)
  df$hcr         <- factor(df$hcr, levels = hcrs)
  df
}
