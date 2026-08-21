library(vegan)
library(ggplot2)
library(here)
library(dplyr)

# Set Up ------------------------------------------------------------------

# Read in pricing stuff
price_mod <- readRDS(here("outputs", 'price_model.RDS'))
grade_df <- read.csv(here("outputs", 'rel_grade_price.csv'))

# no-assessment (perfect information) runs -- these files are ~7 GB each
om_scenarios <- list(
  list(name = "Baseline",
       results = readRDS(here('outputs', 'base_results_noasmt.RDS')),
       sim_list = readRDS(here("outputs", 'base_mse_om.RDS'))),
  list(name = "Regime",
       results = readRDS(here('outputs', 'regime_results_noasmt.RDS')),
       sim_list = readRDS(here("outputs", 'bh_regime_mse_om.RDS'))),
  list(name = "Crash",
       results = readRDS(here('outputs', 'crash_results_noasmt.RDS')),
       sim_list = readRDS(here("outputs", 'bh_crash_mse_om.RDS')))
)

# Time series plots -------------------------------------------------------

term_cond_yr <- 2026  # terminal conditioning (assessment) year

# find a scenario id by HCR settings (matched in om_scenarios[[1]], ids are shared)
# hcr_type "constant" matches the constant-F rule (alpha is ignored there)
find_scenario_id <- function(spr_x, alpha, dyn_b0, stability, hcr_type = "threshold") {
  m <- sapply(om_scenarios[[1]]$results, function(s) {
    sc <- s$scenario
    type_ok  <- if (is.null(sc$hcr_type)) hcr_type == "threshold" else sc$hcr_type == hcr_type
    alpha_ok <- if (hcr_type == "constant") TRUE else isTRUE(all.equal(sc$alpha, alpha))
    type_ok && alpha_ok && isTRUE(all.equal(sc$spr_x, spr_x)) &&
      identical(sc$dyn_b0, dyn_b0) && sc$stability == stability
  })
  sapply(om_scenarios[[1]]$results, function(s) s$scenario$scenario_id)[which(m)[1]]
}

# extract a metric summary for one scenario id, across all OMs
extract_ts <- function(scenario_id, metric, hcr_label, scale = 1) {
  df <- do.call(rbind, lapply(om_scenarios, function(om) {
    all_results <- om$results
    fb <- om$sim_list$feedback_start_yr
    ids <- sapply(all_results, function(s) s$scenario$scenario_id)
    ref <- which(ids == scenario_id)[1]

    mat <- sapply(all_results[[ref]]$results, function(x) as.vector(x$om[[metric]])) / scale
    n_yr <- nrow(mat)
    cal_yr <- (term_cond_yr - (fb - 1) + 1):(term_cond_yr - (fb - 1) + n_yr)

    data.frame(
      om_scenario = om$name,
      hcr         = hcr_label,
      year        = cal_yr,
      med         = apply(mat, 1, median),
      lwr         = apply(mat, 1, quantile, 0.05),
      upr         = apply(mat, 1, quantile, 0.95),
      fb_year     = term_cond_yr - (fb - 1) + fb  # first feedback (projection) year
    )
  }))
  df$om_scenario <- factor(df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
  df
}

# terminal-year replicate values for one scenario id, across all OMs
extract_term <- function(scenario_id, metric, hcr_label, scale = 1) {
  df <- do.call(rbind, lapply(om_scenarios, function(om) { #
    all_results <- om$results
    ids <- sapply(all_results, function(s) s$scenario$scenario_id)
    ref <- which(ids == scenario_id)[1]

    mat <- sapply(all_results[[ref]]$results, function(x) as.vector(x$om[[metric]])) / scale
    if (metric == "Fmort") mat[mat == 1] <- 0 # closure sentinel (see extract_ts)
    data.frame(
      om_scenario = om$name,
      hcr         = hcr_label,
      value       = mat[nrow(mat), ]  # terminal projection year, all reps
    )
  }))
  df$om_scenario <- factor(df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
  df
}

# colourblind-safe palette (ggthemes); black dropped since history uses it
cb_pal <- ggthemes::colorblind_pal()(8)

# function for plotting
plot_ts <- function(df, ylab, zoom = FALSE, term_df = NULL) {
  multi <- length(unique(df$hcr)) > 1 # check to see if multiple hcrs

  # filter historical and proj period
  hist_df <- df %>% filter(year <= fb_year, hcr == df$hcr[1])
  proj_df <- df %>% filter(year >= fb_year)

  p <- ggplot(mapping = aes(year, med))

  if (!zoom) { # if not zooming in
    p <- p +
      geom_vline(data = df, aes(xintercept = fb_year),
                 linetype = 2, colour = "grey50") +
      geom_ribbon(data = hist_df, aes(ymin = lwr, ymax = upr),
                  colour = NA, alpha = 0.20) +
      geom_line(data = hist_df, linewidth = 0.9)
  }

  # if multiple hcrs: interval only for the base (first) HCR, medians for all
  if (multi) {
    base_hcr <- if (is.factor(proj_df$hcr)) levels(proj_df$hcr)[1] else unique(proj_df$hcr)[1]
    p <- p +
      geom_ribbon(data = proj_df %>% filter(hcr == base_hcr),
                  aes(ymin = lwr, ymax = upr),
                  colour = NA, alpha = 0.20) +
      geom_line(data = proj_df, aes(colour = hcr), linewidth = 0.9) +
      scale_colour_manual(values = cb_pal) +
      scale_fill_manual(values = cb_pal)
  } else {
    p <- p +
      geom_ribbon(data = proj_df, aes(ymin = lwr, ymax = upr),
                  colour = NA, alpha = 0.20) +
      geom_line(data = proj_df, linewidth = 0.9)
  }

  # add density on side
  if (!is.null(term_df)) {
    # normalize densities within each scenario
    dens_aes <- if (multi) {
      aes(y = value, x = after_stat(scaled), colour = hcr, fill = hcr)
    } else {
      aes(y = value, x = after_stat(scaled))
    }
    p <- p +
      ggside::geom_ysidedensity(data = term_df, dens_aes,
                                alpha = 0.20, linewidth = 0.6) +
      ggside::scale_ysidex_continuous(breaks = NULL)
  }

  p <- p +
    facet_wrap(~ om_scenario, scales = "free_y") +
    labs(x = "Year", y = ylab, colour = NULL, fill = NULL) +
    coord_cartesian(ylim = if (zoom) NULL else c(0, NA)) +
    theme_bw(base_size = 20) +
    theme(legend.position = if (multi) "top" else "none")

  # ggside theme bits must come after theme_bw() (complete themes reset them)
  if (!is.null(term_df)) p <- p + theme(ggside.panel.scale.x = 0.2)

  p
}

# HCRs to show: label = matching criteria
# first entry is the reference (current management): constant F50, asymmetric
hcr_set <- list(
  "Current (constant F50)" = list(spr_x = 0.5, alpha = NA,   dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "constant"),
  "F50, a=0"               = list(spr_x = 0.5, alpha = 0,    dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.05"            = list(spr_x = 0.5, alpha = 0.05, dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.1"             = list(spr_x = 0.5, alpha = 0.1,  dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.15"            = list(spr_x = 0.5, alpha = 0.15, dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.20"            = list(spr_x = 0.5, alpha = 0.20, dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.25"            = list(spr_x = 0.5, alpha = 0.25, dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.5"             = list(spr_x = 0.5, alpha = 0.5,  dyn_b0 = FALSE, stability = "asymmetric", hcr_type = "threshold")
)

hcr_ids <- sapply(hcr_set, function(h)
  find_scenario_id(h$spr_x, h$alpha, h$dyn_b0, h$stability, h$hcr_type))

# stack all HCRs for one metric; hcr factor keeps the list order
extract_hcrs <- function(ids, metric, scale = 1, fn = extract_ts) {
  df <- do.call(rbind, Map(function(id, label)
    fn(id, metric, label, scale), ids, names(ids)))
  df$hcr <- factor(df$hcr, levels = names(ids))
  df
}

# Recruitment Time-Series -------------------------------------------------

rec_df <- extract_ts(1, "Rec", "ref", scale = 1e6)  # first hcr for reference

p_rec <- ggplot(rec_df %>% filter(year >= 2028), aes(year, med)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ om_scenario, scales = "free") +
  labs(x = "Year", y = "Age-2 Recruitment (millions)") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

print(p_rec)
ggsave(here("figs", "recruitment_timeseries.png"), p_rec,
       width = 15, height = 5, dpi = 300)


# SSB Time-Series ---------------------------------------------------------

# SSB
ssb_base <- extract_ts(hcr_ids[[1]], "SSB", names(hcr_ids)[1], scale = 1e5)
ssb_all  <- extract_hcrs(hcr_ids, "SSB", scale = 1e5)

# terminal-year replicate distributions (right-hand density panels)
ssb_term_base <- extract_term(hcr_ids[[1]], "SSB", names(hcr_ids)[1], scale = 1e5)
ssb_term_all  <- extract_hcrs(hcr_ids, "SSB", scale = 1e5, fn = extract_term)

p_ssb_base <- plot_ts(ssb_base, "SSB", term_df = ssb_term_base)
p_ssb_all  <- plot_ts(ssb_all,  "SSB", term_df = ssb_term_all)
p_ssb_zoom <- plot_ts(ssb_all, "SSB", zoom = TRUE, term_df = ssb_term_all)

print(p_ssb_base)
print(p_ssb_all)
print(p_ssb_zoom)

ggsave(here("figs", "base_hcr_ssb_timeseries.png"), p_ssb_base,
       width = 15, height = 5, dpi = 300)
ggsave(here("figs", "hcr_comparison_ssb_timeseries.png"), p_ssb_all,
       width = 15, height = 8, dpi = 300)
ggsave(here("figs", "hcr_comparison_ssb_timeseries_zoom.png"), p_ssb_zoom,
       width = 15, height = 8, dpi = 300)


# Catch Time-Series -------------------------------------------------------

# Catch
cat_base <- extract_ts(hcr_ids[[1]], "Catch", names(hcr_ids)[1], scale = 1e5)
cat_all  <- extract_hcrs(hcr_ids, "Catch", scale = 1e5)

# terminal-year replicate distributions (right-hand density panels)
cat_term_base <- extract_term(hcr_ids[[1]], "Catch", names(hcr_ids)[1], scale = 1e5)
cat_term_all  <- extract_hcrs(hcr_ids, "Catch", scale = 1e5, fn = extract_term)

p_cat_base <- plot_ts(cat_base, "Catch", term_df = cat_term_base)
p_cat_all  <- plot_ts(cat_all,  "Catch", term_df = cat_term_all)
p_cat_zoom <- plot_ts(cat_all, "Catch", zoom = TRUE, term_df = cat_term_all)

print(p_cat_base)
print(p_cat_all)
print(p_cat_zoom)

ggsave(here("figs", "base_hcr_cat_timeseries.png"), p_cat_base,
       width = 15, height = 5, dpi = 300)
ggsave(here("figs", "hcr_comparison_cat_timeseries.png"), p_cat_all,
       width = 15, height = 8, dpi = 300)
ggsave(here("figs", "hcr_comparison_cat_timeseries_zoom.png"), p_cat_zoom,
       width = 15, height = 8, dpi = 300)

# F Time-Series -----------------------------------------------------------

# F
f_base <- extract_ts(hcr_ids[[1]], "Fmort", names(hcr_ids)[1])
f_all  <- extract_hcrs(hcr_ids, "Fmort")

# terminal-year replicate distributions (right-hand density panels)
f_term_base <- extract_term(hcr_ids[[1]], "Fmort", names(hcr_ids)[1])
f_term_all  <- extract_hcrs(hcr_ids, "Fmort", fn = extract_term)

p_f_base <- plot_ts(f_base, "Fishing mortality (F)", term_df = f_term_base)
p_f_all  <- plot_ts(f_all,  "Fishing mortality (F)", term_df = f_term_all)
p_f_zoom <- plot_ts(f_all, "Fishing mortality (F)", zoom = TRUE, term_df = f_term_all)

print(p_f_base)
print(p_f_all)
print(p_f_zoom)

ggsave(here("figs", "base_hcr_f_timeseries.png"), p_f_base,
       width = 15, height = 5, dpi = 300)
ggsave(here("figs", "hcr_comparison_f_timeseries.png"), p_f_all,
       width = 15, height = 8, dpi = 300)
ggsave(here("figs", "hcr_comparison_f_timeseries_zoom.png"), p_f_zoom,
       width = 15, height = 8, dpi = 300)

# Stock Status (SSB / Bx%) Time-Series ------------------------------------

# SSB relative to the rule's realized biomass reference point (Bx%). bx is
# stored per feedback year, so status exists only in the projection period
extract_status <- function(scenario_id, hcr_label, terminal = FALSE) {
  df <- do.call(rbind, lapply(om_scenarios, function(om) {
    all_results <- om$results
    fb <- om$sim_list$feedback_start_yr
    ids <- sapply(all_results, function(s) s$scenario$scenario_id)
    ref <- which(ids == scenario_id)[1]

    ssb <- sapply(all_results[[ref]]$results, function(x) as.vector(x$om$SSB))
    bx  <- sapply(all_results[[ref]]$results, function(x) x$bx)
    status <- ssb[fb:nrow(ssb), , drop = FALSE] / bx

    if (terminal) {
      data.frame(om_scenario = om$name, hcr = hcr_label,
                 value = status[nrow(status), ])
    } else {
      data.frame(om_scenario = om$name, hcr = hcr_label,
                 year = (term_cond_yr + 1):(term_cond_yr + nrow(status)),
                 med = apply(status, 1, median),
                 lwr = apply(status, 1, quantile, 0.05),
                 upr = apply(status, 1, quantile, 0.95),
                 fb_year = term_cond_yr + 1)
    }
  }))
  df$om_scenario <- factor(df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
  df
}

stack_status <- function(ids, terminal = FALSE) {
  df <- do.call(rbind, Map(function(id, label)
    extract_status(id, label, terminal), ids, names(ids)))
  df$hcr <- factor(df$hcr, levels = names(ids))
  df
}

status_base <- extract_status(hcr_ids[[1]], names(hcr_ids)[1])
status_all  <- stack_status(hcr_ids)
status_term_base <- extract_status(hcr_ids[[1]], names(hcr_ids)[1], terminal = TRUE)
status_term_all  <- stack_status(hcr_ids, terminal = TRUE)

# reference line at B = Bx (status = 1, where the threshold ramp engages)
p_status_base <- plot_ts(status_base, expression(SSB / B[x]), term_df = status_term_base) +
  geom_hline(yintercept = 1, linetype = 3)
p_status_all  <- plot_ts(status_all, expression(SSB / B[x]), term_df = status_term_all) +
  geom_hline(yintercept = 1, linetype = 3)

print(p_status_base)
print(p_status_all)

ggsave(here("figs", "base_hcr_status_timeseries.png"), p_status_base,
       width = 15, height = 5, dpi = 300)
ggsave(here("figs", "hcr_comparison_status_timeseries.png"), p_status_all,
       width = 15, height = 8, dpi = 300)

# HCR Grid Visualization ------------------------------------

# Grab SSB results
spr_x_checkpoints <- seq(0.3, 0.7, by = 0.05)
alpha_checkpoints <- seq(0, 0.5, 0.05)

range_df <- do.call(rbind, lapply(om_scenarios, function(om) {
  fb <- om$sim_list$feedback_start_yr
  do.call(rbind, lapply(om$results, function(s) {
    sc <- s$scenario
    if (!(sc$spr_x %in% spr_x_checkpoints) || !(sc$alpha %in% alpha_checkpoints)) return(NULL)
    ssb   <- sapply(s$results, function(x) as.vector(x$om$SSB))
    catch <- sapply(s$results, function(x) as.vector(x$om$Catch))
    n_yr  <- nrow(ssb); cal_yr <- (term_cond_yr - (fb - 1) + 1):(term_cond_yr - (fb - 1) + n_yr)
    base  <- data.frame(om_scenario = om$name, scenario_id = sc$scenario_id,
                        spr_x = sc$spr_x, alpha = sc$alpha,
                        dyn_b0 = sc$dyn_b0, stability = sc$stability)
    rbind(cbind(base, metric = "SSB",   year = cal_yr, med = apply(ssb,   1, median)),
          cbind(base, metric = "Catch", year = cal_yr, med = apply(catch, 1, median)))
  }))
}))

range_df$om_scenario <- factor(range_df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
proj_range <- range_df %>% filter(year > term_cond_yr)   # projection only

# Summarize mean ssb and catch
sum_proj_range <- proj_range %>%
  group_by(om_scenario, spr_x, alpha, stability, metric) %>%
  summarize(med = median(med) / 1e5)

# Catch AAV across the rule space: per-replicate AAV over the projection period
# (sum |C_y+1 - C_y| / sum C_y), summarized as the median across replicates.
# needs the replicate-level series, so computed separately from the yearly medians
catch_aav_fn <- function(catch_mat, proj_rows) {
  apply(catch_mat[proj_rows, , drop = FALSE], 2, function(cc) {
    if (sum(head(cc, -1)) == 0) return(NA_real_)
    sum(abs(diff(cc))) / sum(head(cc, -1))
  })
}

aav_df <- do.call(rbind, lapply(om_scenarios, function(om) {
  fb <- om$sim_list$feedback_start_yr
  do.call(rbind, lapply(om$results, function(s) {
    sc <- s$scenario
    if (!(sc$spr_x %in% spr_x_checkpoints) || !(sc$alpha %in% alpha_checkpoints)) return(NULL)
    catch <- sapply(s$results, function(x) as.vector(x$om$Catch))
    data.frame(om_scenario = om$name, spr_x = sc$spr_x, alpha = sc$alpha,
               stability = sc$stability, metric = "Catch AAV",
               med = median(catch_aav_fn(catch, fb:nrow(catch)), na.rm = TRUE))
  }))
}))
aav_df$om_scenario <- factor(aav_df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
sum_proj_range <- bind_rows(sum_proj_range, aav_df)

# Reference (current management): constant F50, asymmetric. Its realized
# median projection SSB / Catch per OM is drawn as a red contour on the grids
# (the rule itself has no alpha, so it has no location in threshold space)
ref_id <- find_scenario_id(0.5, NA, FALSE, "asymmetric", hcr_type = "constant")

ref_vals <- do.call(rbind, lapply(om_scenarios, function(om) {
  fb  <- om$sim_list$feedback_start_yr
  ids <- sapply(om$results, function(s) s$scenario$scenario_id)
  s   <- om$results[[which(ids == ref_id)[1]]]
  ssb   <- sapply(s$results, function(x) as.vector(x$om$SSB))
  catch <- sapply(s$results, function(x) as.vector(x$om$Catch))
  n_yr  <- nrow(ssb)
  cal_yr <- (term_cond_yr - (fb - 1) + 1):(term_cond_yr - (fb - 1) + n_yr)
  proj  <- cal_yr > term_cond_yr
  # same summary as sum_proj_range: median over projection years of yearly medians
  data.frame(om_scenario = om$name, metric = c("SSB", "Catch"),
             ref = c(median(apply(ssb,   1, median)[proj]),
                     median(apply(catch, 1, median)[proj])) / 1e5)
}))
# reference rule catch AAV (same summary as aav_df)
ref_aav <- do.call(rbind, lapply(om_scenarios, function(om) {
  fb  <- om$sim_list$feedback_start_yr
  ids <- sapply(om$results, function(s) s$scenario$scenario_id)
  s   <- om$results[[which(ids == ref_id)[1]]]
  catch <- sapply(s$results, function(x) as.vector(x$om$Catch))
  data.frame(om_scenario = om$name, metric = "Catch AAV",
             ref = median(catch_aav_fn(catch, fb:nrow(catch)), na.rm = TRUE))
}))
ref_vals <- rbind(ref_vals, ref_aav)
ref_vals$om_scenario <- factor(ref_vals$om_scenario, levels = c("Baseline", "Regime", "Crash"))

# plot grid: one subplot per OM so the fill scale is computed within each panel.
# the red contour marks the reference rule (constant F50) realized value within
# that panel; when its level sits off the panel's surface (no contour possible),
# it is reported as a red label instead
grid_plot <- function(metric, stability) {
  dat  <- sum_proj_range %>%
    filter(metric == .env$metric, stability == .env$stability)
  refs <- ref_vals %>% filter(metric == .env$metric)
  oms  <- levels(droplevels(dat$om_scenario))

  panels <- lapply(seq_along(oms), function(i) {
    d   <- dat %>% filter(om_scenario == oms[i])
    ref <- refs$ref[refs$om_scenario == oms[i]]

    p <- ggplot(d, aes(x = spr_x, y = alpha, z = med)) +
      geom_tile(aes(fill = med)) +
      geomtextpath::geom_textcontour(color = 'white', size = 5, linewidth = 1) +
      # Catch AAV: lower is better, so reverse the scale (yellow = low)
      scale_fill_viridis_c(direction = if (metric == "Catch AAV") -1 else 1) +
      scale_x_continuous(breaks = seq(0.3, 0.7, 0.1), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
      facet_wrap(~om_scenario) + # single facet, kept for the strip label
      guides(fill = guide_colourbar(barwidth = grid::unit(6, 'cm'))) +
      theme_bw(base_size = 20) +
      theme(legend.position = 'top') +
      labs(x = 'SPR rate (Fx%)', y = if (i == 1) 'Alpha' else NULL, fill = metric)

    if (length(ref)) {
      if (ref >= min(d$med) && ref <= max(d$med)) {
        p <- p + geomtextpath::geom_textcontour(breaks = ref, color = 'red',
                                                size = 5, linewidth = 1.2)
      } else {
        p <- p + annotate('label', x = 0.31, y = 0.49, hjust = 0, vjust = 1,
                          label = sprintf('Constant F50 = %.3g', ref),
                          color = 'red', size = 5, alpha = 0.7)
      }
    }
    p
  })

  patchwork::wrap_plots(panels, nrow = 1)
}

ggsave(here('figs', 'grid_ssb_asym.png'), grid_plot('SSB', 'asymmetric'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_cat_asym.png'), grid_plot('Catch', 'asymmetric'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_ssb_none.png'), grid_plot('SSB', 'none'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_cat_none.png'), grid_plot('Catch', 'none'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_ssb_sym.png'), grid_plot('SSB', 'symmetric'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_cat_sym.png'), grid_plot('Catch', 'symmetric'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_aav_asym.png'), grid_plot('Catch AAV', 'asymmetric'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_aav_none.png'), grid_plot('Catch AAV', 'none'), width = 19, height = 7, dpi = 300)
ggsave(here('figs', 'grid_aav_sym.png'), grid_plot('Catch AAV', 'symmetric'), width = 19, height = 7, dpi = 300)

# # Process Results for NMDS ---------------------------------------------------------
#
# results_list <- lapply(om_scenarios, function(om) {
#
#   all_results <- om$results
#   sim_list <- om$sim_list
#
#   perf <- do.call(rbind, lapply(1:length(all_results), function(i) {
#
#     sc <- all_results[[i]]$scenario
#     ssb_mat <- sapply(all_results[[i]]$results, function(x) x$om$SSB) # get ssb
#     catch_mat <- sapply(all_results[[i]]$results, function(x) x$om$Catch) # get catch
#     f_mat <- sapply(all_results[[i]]$results, function(x) x$om$Fmort) # get fmort
#     catch_mat[catch_mat == 0] <- 1e-3 # add 0 to guard against aav calcs
#     naa_array <- abind::abind(lapply(all_results[[i]]$results, function(x) x$om$NAA), along = 4) # get NAA
#     caa_array <- abind::abind(lapply(all_results[[i]]$results, function(x) x$om$CAA), along = 4) # get CAA
#     fb <- sim_list$feedback_start_yr # get feedback year
#     proj_rows <- fb:nrow(ssb_mat) # get projection rows
#     ssb_threshold <- max(all_results[[i]]$results[[1]]$om$SSB[seq_len(fb - 1)]) * 0.1 # SSB threshold = max * 0.1
#
#     data.frame(
#       # get hcr type
#       scenario_id = sc$scenario_id, spr_x = sc$spr_x, alpha = sc$alpha,
#       dyn_b0 = sc$dyn_b0, stability = sc$stability,
#
#       # median ssb
#       med_ssb = median(ssb_mat[proj_rows, ]),
#       # mean age
#       mean_age = mean(apply(naa_array[proj_rows,,,,drop = F], c(1, 4), function(mat) {
#         n_total <- sum(mat)
#         if (n_total == 0) return(NA)
#         sum(mat * 2:31) / n_total
#       })),
#       # age eveness
#       pielou_j = mean(apply(naa_array[proj_rows,,,,drop = F], c(1, 4), function(mat) {
#         n_total <- sum(mat)
#         if (n_total == 0) return(NA)
#         p <- mat / n_total
#         p <- p[p > 0]
#         H <- -sum(p * log(p))
#         H / log(length(p))
#       })),
#
#       # economic value
#       econ_value = mean(sapply(1:ncol(ssb_mat), function(j) {
#         proj_vals <- sapply(proj_rows, function(y) {
#           n_total <- sum(naa_array[y,,,j]) # get total N
#           rel_price_base <- as.numeric(predict(price_mod, newdata = data.frame(n = n_total), type = "response")) # predict price based on N
#           caa_female <- caa_array[y,,1,j] # get female catch-at-age
#           caa_male   <- caa_array[y,,2,j] # get male catch-at-age
#           # multiply, grade relative pricing, and predicted relative pricing
#           value_male   <- sum(caa_male * grade_df$rel_price_male * rel_price_base, na.rm = TRUE)
#           value_female <- sum(caa_female * grade_df$rel_price_female * rel_price_base, na.rm = TRUE)
#           value_male + value_female
#         })
#         median(proj_vals)
#       })),
#
#       # proportion of years below ssb threshold
#       p_crash = mean(ssb_mat[proj_rows, ] < ssb_threshold),
#
#       # median catch
#       med_catch = median(catch_mat[proj_rows, ]),
#       # catch aav
#       catch_aav = median(sapply(1:ncol(catch_mat), function(j) {
#         cc <- catch_mat[proj_rows, j]
#         mean(abs(diff(cc)) / head(cc, -1), na.rm = TRUE)
#       }))
#     )
#   }))
#
#   # Include p_crash for non-baseline scenarios
#   if (om$name == "Crash") {
#     use_metrics <- c("med_ssb", "med_catch", "catch_aav", "p_crash",
#                      "mean_age", "pielou_j", 'econ_value')
#   } else {
#     use_metrics <- c("med_ssb", "med_catch", "catch_aav", "mean_age", "pielou_j", 'econ_value')
#   }
#
#   # fit NMDS
#   set.seed(42)
#   perf_scaled <- scale(as.matrix(perf[, use_metrics])) # rescale everything
#   nmds <- metaMDS(dist(perf_scaled), k = 2, trymax = 200, autotransform = FALSE) # 2 axes
#   ef <- envfit(nmds, perf_scaled, permutations = 999) # get vectors
#
#   # set up plotting for vectors
#   vec_df <- as.data.frame(ef$vectors$arrows * sqrt(ef$vectors$r))
#   vec_df$label <- rownames(vec_df)
#   arrow_scale <- 0.6 * max(abs(nmds$points))
#   vec_df$NMDS1 <- vec_df$NMDS1 * arrow_scale
#   vec_df$NMDS2 <- vec_df$NMDS2 * arrow_scale
#   vec_df$om_scenario <- om$name
#   df <- data.frame(NMDS1 = nmds$points[,1], NMDS2 = nmds$points[,2], perf, om_scenario = om$name)
#
#   # figure out potential outliers
#   dists <- sqrt(df$NMDS1^2 + df$NMDS2^2)
#   df$label <- ifelse(dists > quantile(dists, 0.90), df$scenario_id, NA)
#
#   list(points = df, vectors = vec_df)
# })
#
# # combine results into a df
# nmds_df <- do.call(rbind, lapply(results_list, `[[`, "points"))
# vec_df  <- do.call(rbind, lapply(results_list, `[[`, "vectors"))
#
# # get scenarios
# nmds_df$om_scenario <- factor(nmds_df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
# vec_df$om_scenario  <- factor(vec_df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
#
# nmds_df$dyn_b0    <- factor(nmds_df$dyn_b0, levels = c(FALSE, TRUE),
#                             labels = c("Static B0", "Dynamic B0"))
# nmds_df$stability <- factor(nmds_df$stability,
#                             levels = c("none", "symmetric", "asymmetric"),
#                             labels = c("No Constraint", "Symmetric", "Asymmetric"))
#
#
# # Plot NMDS ---------------------------------------------------------------
#
# pos <- position_jitter(width = 0.15, height = 0.15, seed = 42)
#
# nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
#   geom_point(aes(colour = stability, shape = dyn_b0, size = spr_x),
#              alpha = 0.5, position = pos) +
#   geom_text(aes(label = label), size = 2.5, vjust = -1, colour = "grey30",
#             na.rm = TRUE, position = pos) +
#   geom_segment(data = vec_df,
#                aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
#                arrow = arrow(length = unit(0.2, "cm")),
#                colour = "grey40", linewidth = 0.5, inherit.aes = FALSE) +
#   ggrepel::geom_label_repel(data = vec_df,
#                             aes(x = NMDS1, y = NMDS2, label = label),
#                             size = 3, colour = "grey30", fontface = "italic",
#                             fill = "white", label.size = 0.2,
#                             nudge_x = vec_df$NMDS1 * 0.3,
#                             nudge_y = vec_df$NMDS2 * 0.3,
#                             segment.color = "grey70", segment.size = 0.3,
#                             box.padding = 0.4, point.padding = 0.2,
#                             min.segment.length = 0, inherit.aes = FALSE) +
#   facet_wrap(~ om_scenario) +
#   scale_shape_manual(values = c("Static B0" = 16, "Dynamic B0" = 17)) +
#   scale_size_continuous(range = c(1.5, 5), breaks = c(0.3, 0.5, 0.7)) +
#   scale_linetype_discrete() +
#   labs(colour = "Stability Type", shape = "Reference Point",
#        size = expression(SPR[x]), linetype = expression(alpha)) +
#   theme_bw(base_size = 14) +
#   coord_cartesian(ylim = c(-7.5, 7.5), xlim = c(-7.5, 7.5))
#
# print(nmds_plot)
#
# ggsave(here("figs", "nmds_plot.png"), nmds_plot,
#        width = 13, height = 7, dpi = 300)

