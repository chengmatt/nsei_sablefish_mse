# Purpose: Head-to-head comparison of the closed-loop MSE runs with an estimation
#          model ("With assessment") against the perfect-information runs
#          ("Perfect information"). Both run sets share the same operating models,
#          seeds and conditioning period, so differences are attributable to
#          assessment error alone.
#          Restricted to the 7 harvest control rules present in both run sets
#          (the assessment runs used a reduced grid).
#          Requires outputs/mse_summaries.RDS (see R/extract_mse_summaries.R).
# Creator: Matthew LH. Cheng

source(here::here("R", "functions", "mse_viz_utils.R"))
library(tidyr)

PREF <- "cmp_"
save_fig <- function(name, plot, width = 15, height = 8) {
  ggsave(here("figs", paste0(PREF, name, ".png")), plot,
         width = width, height = height, dpi = 300)
}

ts_specs <- list(
  ssb    = list(metric = "SSB",    lab = "SSB",                   scale = 1e5),
  cat    = list(metric = "Catch",  lab = "Catch",                 scale = 1e5),
  f      = list(metric = "Fmort",  lab = "Fishing mortality (F)", scale = 1),
  status = list(metric = "status", lab = expression(SSB / B[x]),  scale = 1)
)

# Time series: reference rule, both modes -----------------------------------

for (nm in names(ts_specs)) {
  sp <- ts_specs[[nm]]

  base_ts   <- build_ts(sp$metric,   sp$scale, hcrs = ref_hcr)
  base_term <- build_term(sp$metric, sp$scale, hcrs = ref_hcr)

  p <- plot_ts(base_ts, sp$lab, colour_by = "mode", pal = mode_pal,
               term_df = base_term, ribbon_all = TRUE)
  if (nm == "status") p <- p + geom_hline(yintercept = 1, linetype = 3)

  save_fig(paste0("base_hcr_", nm, "_timeseries"), p, 15, 8)

  # every rule, one row per rule
  all_ts <- build_ts(sp$metric, sp$scale)
  p_all <- plot_ts(all_ts, sp$lab, colour_by = "mode", pal = mode_pal,
                   facet = "grid", ribbon_all = TRUE)
  if (nm == "status") p_all <- p_all + geom_hline(yintercept = 1, linetype = 3)

  save_fig(paste0("all_hcr_", nm, "_timeseries"), p_all, 15, 20)
}

# Performance metrics --------------------------------------------------------

perf <- perf_metrics()

metric_order <- c("med_ssb", "med_catch", "catch_aav", "p_crash", "p_closure")
perf_long <- perf %>%
  pivot_longer(all_of(metric_order), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_order,
                         labels = c("Median SSB\n(100k t)", "Median catch\n(100k t)",
                                    "Catch AAV", "P(SSB < 10%\nhistorical max)",
                                    "P(closure)")))

p_perf <- ggplot(perf_long, aes(value, hcr, colour = mode, shape = mode)) +
  geom_line(aes(group = hcr), colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3.5) +
  facet_grid(om_scenario ~ metric, scales = "free_x") +
  scale_colour_manual(values = mode_pal) +
  scale_y_discrete(limits = rev(hcr_levels)) +
  labs(x = NULL, y = NULL, colour = NULL, shape = NULL) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top")

save_fig("performance_metrics", p_perf, 17, 10)

# Change attributable to assessment error ------------------------------------

# stock/fishery scale metrics are compared as a percent change; probabilities are
# compared as a difference in percentage points, since the perfect-information
# value is often zero
pct_metrics <- c("med_ssb", "med_catch", "catch_aav", "term_status")
delta <- perf %>%
  pivot_longer(all_of(c(metric_order, "term_status")),
               names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = mode, values_from = value) %>%
  rename(perfect = `Perfect information`, assessment = `With assessment`) %>%
  mutate(
    change = ifelse(metric %in% pct_metrics,
                    100 * (assessment - perfect) / perfect,
                    100 * (assessment - perfect)),
    unit = ifelse(metric %in% pct_metrics, "% change", "percentage points")
  )

# direction of preference, so a bar can be labelled better/worse rather than
# just positive/negative
better_when_higher <- c(med_ssb = TRUE, med_catch = TRUE, catch_aav = FALSE,
                        p_crash = FALSE, p_closure = FALSE)

delta_long <- delta %>%
  filter(metric %in% metric_order) %>%
  mutate(direction = ifelse(change == 0, "No change",
                     ifelse((change > 0) == better_when_higher[metric],
                            "Better with assessment", "Worse with assessment"))) %>%
  mutate(metric = factor(metric, levels = metric_order,
                         labels = c("Median SSB\n(% change)", "Median catch\n(% change)",
                                    "Catch AAV\n(% change)", "P(SSB < 10% hist. max)\n(pp change)",
                                    "P(closure)\n(pp change)")))

p_delta <- ggplot(delta_long, aes(change, hcr, fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_col(width = 0.65, alpha = 0.85) +
  facet_grid(om_scenario ~ metric, scales = "free_x") +
  scale_fill_manual(values = c("Better with assessment" = "#0072B2",
                               "Worse with assessment"  = "#D55E00",
                               "No change"              = "grey70")) +
  scale_y_discrete(limits = rev(hcr_levels)) +
  labs(x = "Assessment run relative to perfect information", y = NULL, fill = NULL) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top")

save_fig("performance_change", p_delta, 17, 10)

# Trade-off space ------------------------------------------------------------

# Percent change from adding an assessment, rule by rule, on axes that are shared
# across operating models so the panels are directly comparable. Changes are
# computed replicate by replicate against the matching perfect-information
# replicate, so the shaded box and crosshairs are simulation intervals on the
# change itself rather than on either run in isolation.
# Reading: right = more catch, up = more variable catch. The lower-right quadrant
# is the only place assessment error helps on both axes.
TRADE_PROBS <- c(0.25, 0.75)   # shaded interval; c(0.05, 0.95) for the full spread
chg <- paired_change(vars = c("med_catch", "catch_aav"), probs = TRADE_PROBS)

chg_wide <- chg %>%
  select(om_scenario, hcr, var, med, lwr, upr) %>%
  pivot_wider(names_from = var, values_from = c(med, lwr, upr), names_sep = "_")

p_trade <- ggplot(chg_wide, aes(med_med_catch, med_catch_aav, colour = hcr)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_rect(aes(xmin = lwr_med_catch, xmax = upr_med_catch,
                ymin = lwr_catch_aav, ymax = upr_catch_aav, fill = hcr),
            colour = NA, alpha = 0.10) +
  geom_linerange(aes(xmin = lwr_med_catch, xmax = upr_med_catch), linewidth = 0.7, alpha = 0.8) +
  geom_linerange(aes(ymin = lwr_catch_aav, ymax = upr_catch_aav), linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 3.5) +
  facet_wrap(~ om_scenario) +
  scale_colour_manual(values = hcr_pal) +
  scale_fill_manual(values = hcr_fill) +
  labs(x = "Change in median catch (%)", y = "Change in catch AAV (%)",
       colour = NULL, fill = NULL) +
  theme_bw(base_size = 18) +
  theme(legend.position = "top")

save_fig("tradeoff", p_trade, 16, 7)

# same view without the intervals, for when only the central tendency is wanted
p_trade_med <- ggplot(chg_wide, aes(med_med_catch, med_catch_aav, colour = hcr)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = hcr), size = 4.2, show.legend = FALSE,
                           min.segment.length = 0, box.padding = 0.5) +
  facet_wrap(~ om_scenario, scales = "free") +
  scale_colour_manual(values = hcr_pal) +
  labs(x = "Change in median catch (%)", y = "Change in catch AAV (%)", colour = NULL) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

save_fig("tradeoff_medians", p_trade_med, 16, 6)

# Table ----------------------------------------------------------------------

tbl <- delta %>%
  mutate(across(c(perfect, assessment, change), ~ signif(.x, 4))) %>%
  arrange(om_scenario, metric, hcr) %>%
  select(om_scenario, hcr, metric, perfect, assessment, change, unit)

write.csv(tbl, here("outputs", "asmt_vs_noasmt_performance.csv"), row.names = FALSE)

cat("\n--- median change from adding an assessment (across HCRs) ---\n")
print(delta %>%
        group_by(om_scenario, metric, unit) %>%
        summarize(median_change = round(median(change), 2), .groups = "drop") %>%
        pivot_wider(names_from = om_scenario, values_from = median_change) %>%
        as.data.frame())
cat("\nwrote figs/", PREF, "*.png and outputs/asmt_vs_noasmt_performance.csv\n", sep = "")
