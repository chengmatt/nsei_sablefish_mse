# Purpose: Figures for the closed-loop MSE runs that include an estimation model
#          (assessment). Mirrors R/visualize_mse_results.R for the reduced HCR grid
#          the assessment runs used, and adds the assessment-error diagnostics that
#          only exist once an EM is in the loop.
#          Requires outputs/mse_summaries.RDS (see R/extract_mse_summaries.R).
# Creator: Matthew LH. Cheng

source(here::here("R", "functions", "mse_viz_utils.R"))

MODE  <- "asmt"
PREF  <- "asmt_"     # figure prefix
W1    <- 15; H1 <- 5 # single-row figures
W2    <- 15; H2 <- 8 # figures with a legend row

save_fig <- function(name, plot, width = W1, height = H1) {
  ggsave(here("figs", paste0(PREF, name, ".png")), plot,
         width = width, height = height, dpi = 300)
}

# Recruitment Time-Series -------------------------------------------------

rec_df <- build_ts("Rec", scale = 1e6, hcrs = ref_hcr, modes = MODE)

p_rec <- ggplot(rec_df %>% filter(year >= term_cond_yr + 2), aes(year, med)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ om_scenario, scales = "free") +
  labs(x = "Year", y = "Age-2 Recruitment (millions)") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")

save_fig("recruitment_timeseries", p_rec)

# SSB / Catch / F Time-Series ---------------------------------------------

# one entry per figure family: metric, axis label, scale
ts_specs <- list(
  ssb    = list(metric = "SSB",    lab = "SSB",                     scale = 1e5),
  cat    = list(metric = "Catch",  lab = "Catch",                   scale = 1e5),
  f      = list(metric = "Fmort",  lab = "Fishing mortality (F)",   scale = 1),
  status = list(metric = "status", lab = expression(SSB / B[x]),    scale = 1)
)

for (nm in names(ts_specs)) {
  sp <- ts_specs[[nm]]

  base_ts   <- build_ts(sp$metric,   sp$scale, hcrs = ref_hcr,    modes = MODE)
  base_term <- build_term(sp$metric, sp$scale, hcrs = ref_hcr,    modes = MODE)
  all_ts    <- build_ts(sp$metric,   sp$scale, hcrs = hcr_levels, modes = MODE)
  all_term  <- build_term(sp$metric, sp$scale, hcrs = hcr_levels, modes = MODE)

  p_base <- plot_ts(base_ts, sp$lab, term_df = base_term)
  p_all  <- plot_ts(all_ts,  sp$lab, colour_by = "hcr", pal = hcr_pal,
                    fill_pal = hcr_fill, term_df = all_term)

  # SSB / Bx: mark the point where the threshold ramp engages
  if (nm == "status") {
    p_base <- p_base + geom_hline(yintercept = 1, linetype = 3)
    p_all  <- p_all  + geom_hline(yintercept = 1, linetype = 3)
  }

  save_fig(paste0("base_hcr_", nm, "_timeseries"), p_base, W1, H1)
  save_fig(paste0("hcr_comparison_", nm, "_timeseries"), p_all, W2, H2)

  # zoomed view (projection only, free y) for the stock/fishery metrics
  if (nm %in% c("ssb", "cat", "f")) {
    p_zoom <- plot_ts(all_ts, sp$lab, colour_by = "hcr", pal = hcr_pal,
                      fill_pal = hcr_fill, term_df = all_term, zoom = TRUE)
    save_fig(paste0("hcr_comparison_", nm, "_timeseries_zoom"), p_zoom, W2, H2)
  }
}

# Assessment Error --------------------------------------------------------

# Relative error of the SSB the assessment reported for its own terminal year --
# the quantity the HCR acts on. Positive = the assessment saw more fish than the
# OM actually had.
re_df <- build_ts("ssb_re", hcrs = hcr_levels, modes = MODE)

p_re <- ggplot(re_df %>% filter(!is.na(med)), aes(year, med, colour = hcr)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_ribbon(data = re_df %>% filter(!is.na(med), hcr == ref_hcr),
              aes(ymin = lwr, ymax = upr, fill = hcr), colour = NA, alpha = 0.20) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ om_scenario) +
  scale_colour_manual(values = hcr_pal) +
  scale_fill_manual(values = hcr_fill) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = "none") +   # the interval is the reference rule's; one legend is enough
  labs(x = "Year", y = "Terminal SSB relative error", colour = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top")

save_fig("em_ssb_relative_error", p_re, W2, H2)

# What the rule saw vs what was there: perceived stock status (EM SSB / Bx)
# against true status (OM SSB / the same Bx), for the reference rule
status_cmp <- rbind(
  transform(build_ts("status",           hcrs = ref_hcr, modes = MODE), series = "True (OM)"),
  transform(build_ts("status_perceived", hcrs = ref_hcr, modes = MODE), series = "Perceived (EM)")
)
status_cmp$series <- factor(status_cmp$series, levels = c("True (OM)", "Perceived (EM)"))

p_perc <- ggplot(status_cmp %>% filter(!is.na(med)), aes(year, med, colour = series)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = series), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ om_scenario, scales = "free_y") +
  scale_colour_manual(values = unname(mode_pal)) +
  scale_fill_manual(values = unname(mode_pal)) +
  labs(x = "Year", y = expression(SSB / B[x]), colour = NULL, fill = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top")

save_fig("em_perceived_vs_true_status", p_perc, W2, H2)

# Convergence -------------------------------------------------------------

conv_df <- do.call(rbind, lapply(om_levels, function(om)
  do.call(rbind, lapply(hcr_levels, function(h) {
    r  <- rec_of(om, MODE)[[h]]
    fb <- fb_of(om, MODE)
    yrs <- (term_cond_yr + 1):(term_cond_yr + nrow(r$pdHess))
    data.frame(om_scenario = om, hcr = h, year = yrs,
               pd_rate  = rowMeans(r$pdHess, na.rm = TRUE),
               max_grad = apply(r$grad, 1, max, na.rm = TRUE))
  }))))
conv_df$om_scenario <- factor(conv_df$om_scenario, levels = om_levels)
conv_df$hcr         <- factor(conv_df$hcr, levels = hcr_levels)

p_pd <- ggplot(conv_df, aes(year, pd_rate, colour = hcr)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ om_scenario) +
  scale_colour_manual(values = hcr_pal) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(NA, 1)) +
  labs(x = NULL, y = "Converged (pd Hessian)", colour = NULL) +
  theme_bw(base_size = 16) + theme(legend.position = "top")

p_gr <- ggplot(conv_df, aes(year, max_grad, colour = hcr)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ om_scenario) +
  scale_colour_manual(values = hcr_pal) +
  scale_y_log10() +
  labs(x = "Year", y = "Max |gradient| (worst replicate)", colour = NULL) +
  theme_bw(base_size = 16) + theme(legend.position = "none")

save_fig("em_convergence", patchwork::wrap_plots(p_pd, p_gr, ncol = 1), W2, 11)

cat(sprintf("\nassessment convergence: %.3f%% of fits positive-definite; worst gradient %.3g\n",
            100 * mean(conv_df$pd_rate), max(conv_df$max_grad)))
cat("wrote figs/", PREF, "*.png\n", sep = "")

# Retrospective Bias ------------------------------------------------------

# Every closed-loop assessment is a peel of every later one, so the retrospective
# pattern is visible exactly as it would be in practice -- without reference to
# the operating model. This is the complement to the relative-error figure
# above, which compares the assessment against the truth it could not see.
# Retrospective fan for a single replicate, one panel per management procedure:
# each line is one assessment's view of history, drawn over the projection
# period, against the operating model
fan_df <- do.call(rbind, lapply(om_levels, function(om)
  do.call(rbind, lapply(hcr_levels, function(h) {
    ex <- rec_of(om, MODE)[[h]]$ssb_em_example
    fb <- fb_of(om, MODE)
    keep <- seq(1, length(ex), by = 6)          # every 6th assessment
    do.call(rbind, lapply(keep, function(i) {
      v <- ex[[i]]
      if (is.null(v)) return(NULL)
      t <- seq_along(v)
      data.frame(om_scenario = om, hcr = h,
                 asmt_year = term_cond_yr - (fb - 1) + (fb + i - 1),
                 year = term_cond_yr - (fb - 1) + t, ssb = v / 1e5)[t >= fb, ]
    }))
  }))))
fan_df$om_scenario <- factor(fan_df$om_scenario, levels = om_levels)
fan_df$hcr         <- factor(fan_df$hcr, levels = hcr_levels)

truth_df <- do.call(rbind, lapply(om_levels, function(om)
  do.call(rbind, lapply(hcr_levels, function(h) {
    r  <- rec_of(om, MODE)[[h]]
    fb <- fb_of(om, MODE)
    data.frame(om_scenario = om, hcr = h, year = years_of(om, MODE),
               ssb = r$SSB[, 1] / 1e5)[fb:nrow(r$SSB), ]
  }))))
truth_df$om_scenario <- factor(truth_df$om_scenario, levels = om_levels)
truth_df$hcr         <- factor(truth_df$hcr, levels = hcr_levels)

p_fan <- ggplot(fan_df, aes(year, ssb, group = asmt_year, colour = asmt_year)) +
  geom_line(linewidth = 1.1) +
  geom_point(data = fan_df %>% group_by(om_scenario, hcr, asmt_year) %>%
               slice_max(year, n = 1) %>% ungroup(), size = 3.2) +
  geom_line(data = truth_df, aes(year, ssb), inherit.aes = FALSE,
            colour = "black", linewidth = 1.5) +
  facet_grid(hcr ~ om_scenario, scales = "free_y") +
  scale_colour_viridis_c(option = "viridis") +
  labs(x = "Year", y = "SSB", colour = "Assessment year") +
  theme_bw(base_size = 20) +
  theme(legend.position = "top",
        legend.key.width = grid::unit(3, "cm"))

save_fig("retro_fan", p_fan, 15, 22)
