# Purpose: Make figures
# Creator: Matthew LH. Cheng

# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
source(here("R", "functions", "mse_viz_utils.R"))

save_fig <- function(name, plot, width, height) {
  ggsave(here("figs", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300)
}

# labels and level orders used across the figures
hcr_labels  <- c("Constant F50%", "F40%, α = 0", "F40%, α = 0.25", "F40%, α = 0.5",
              "F50%, α = 0", "F50%, α = 0.25", "F50%, α = 0.5")
metric_labels <- c(med_catch = "Median catch (t)", med_ssb = "Median SSB (t)",
               catch_aav = "Catch AAV", p_crash = "P(depletion)",
               p_closure = "P(closure)")


# Performance summary -----------------------------------------------------
percent_change_metrics  <- c("med_ssb", "med_catch", "catch_aav", "term_status")
metric_order <- c("med_ssb", "med_catch", "catch_aav", "p_crash", "p_closure")

performance_delta <- performance_metrics() %>%
  pivot_longer(all_of(c(metric_order, "term_status")),
               names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = mode, values_from = value) %>%
  rename(perfect = `Perfect information`, assessment = `With assessment`) %>%
  mutate(
    change = ifelse(metric %in% percent_change_metrics,
                    100 * (assessment - perfect) / perfect,
                    100 * (assessment - perfect)),
    unit = ifelse(metric %in% percent_change_metrics, "% change", "percentage points")
  )

performance_delta %>%
  mutate(across(c(perfect, assessment, change), ~ signif(.x, 4))) %>%
  arrange(om_scenario, metric, hcr) %>%
  select(om_scenario, hcr, metric, perfect, assessment, change, unit) %>%
  write.csv(here("outputs", "asmt_vs_noasmt_performance.csv"), row.names = FALSE)


# Figure 1: recruitment by OM ---------------------------------------------
recruitment <- build_timeseries("Rec", hcrs = reference_hcr, modes = "noasmt") %>%
  filter(year >= terminal_conditioning_year + 1)   # first projection year with simulated recruitment

recruitment_plot <- ggplot(recruitment, aes(year, median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60", colour = NA, alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ om_scenario, scales = "free_y") +
  scale_x_continuous(breaks = seq(2030, 2080, 10)) +
  labs(x = "Year", y = "Age-2 recruitment (millions)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_fig("fig1_recruitment_scenarios", recruitment_plot, 10, 3.6)


# Figure 2: HCR hcr_shape_curves ----------------------------------------------------
hcr_f_at_status <- function(status, f_reference_point, alpha = 0) {
  f_reference_point * pmin(1, pmax(0, (status - alpha) / (1 - alpha)))
}

# fishing mortality vlaues
f_reference_points     <- c("F40%" = 0.12, "F45%" = 0.10, "F50%" = 0.085,
                "F55%" = 0.072, "F60%" = 0.062)
quantity_labels <- c("Fishing mortality (F)", "Relative catch")

hcr_shape_curves <- expand_grid(status = seq(0, 1.5, 0.005), f_label = names(f_reference_points),
                      closure_threshold = c(0, 0.1, 0.25, 0.5)) %>%
  mutate(fishing_mortality = hcr_f_at_status(status, f_reference_points[f_label], closure_threshold),
         catch = fishing_mortality * status) %>%
  pivot_longer(c(fishing_mortality, catch), names_to = "quantity", values_to = "value") %>%
  mutate(f_label = factor(f_label, levels = names(f_reference_points)),
         closure_threshold = factor(closure_threshold),
         quantity = factor(c(fishing_mortality = quantity_labels[1],
                             catch = quantity_labels[2])[quantity],
                           levels = quantity_labels))

hcr_shapes_plot <- ggplot(hcr_shape_curves, aes(status, value, colour = closure_threshold)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 1.1) +
  scale_colour_viridis_d(option = "plasma", end = 0.85) +
  facet_grid(quantity ~ f_label, scales = "free_y", switch = "y") +
  labs(x = expression(SSB / B[x]), y = NULL, colour = expression(alpha)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background.y = element_blank())

save_fig("hcr_shape_both", hcr_shapes_plot, 16, 8)


# Figure 3 / S6-S7: screening profiles ------------------------------------
grid_metrics <- read.csv(here("outputs", "grid_metrics.csv")) %>%
  mutate(om_scenario = factor(om_scenario, om_scenario_names))
closure_thresholds_shown <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)   # 6 of the 11 values run, for legibility

reference_rule_metrics <- grid_metrics %>% filter(hcr_type == "constant") %>%
  select(om_scenario, all_of(names(metric_labels))) %>%
  pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], metric_labels))

shared_y_limits <- grid_metrics %>% filter(hcr_type == "threshold", alpha %in% closure_thresholds_shown) %>%
  select(all_of(names(metric_labels))) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  bind_rows(reference_rule_metrics %>% mutate(metric = names(metric_labels)[match(metric, metric_labels)])) %>%
  group_by(metric) %>% summarise(lowest = min(value), highest = max(value), .groups = "drop") %>%
  pivot_longer(c(lowest, highest), values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], metric_labels),
         om_scenario = factor("Baseline", om_scenario_names), spr_x = 0.5)

screening_profile_plot <- function(stability) {
  profiles <- grid_metrics %>% filter(hcr_type == "threshold", stability == .env$stability,
                       alpha %in% closure_thresholds_shown) %>%
    select(om_scenario, spr_x, alpha, all_of(names(metric_labels))) %>%
    pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric_labels[metric], metric_labels), alpha = factor(alpha))
  ggplot(profiles, aes(spr_x, value, colour = alpha, group = alpha)) +
    geom_hline(data = reference_rule_metrics, aes(yintercept = value, linetype = "Constant F50% (current rule)"),
               colour = "grey30") +
    geom_blank(data = shared_y_limits, aes(spr_x, value), inherit.aes = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    facet_grid(metric ~ om_scenario, scales = "free_y") +
    scale_colour_viridis_d(option = "plasma", end = 0.9,
                           name = expression(Closure~threshold~(alpha))) +
    scale_linetype_manual(values = c("Constant F50% (current rule)" = 2), name = NULL) +
    scale_x_continuous(breaks = seq(0.3, 0.7, 0.1),
                       labels = function(x) paste0("F", x * 100)) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    labs(x = "SPR target", y = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

save_fig("fig3_screening_profiles",      screening_profile_plot("asymmetric"), 9, 10.5)
save_fig("figS_screening_profiles_sym",  screening_profile_plot("symmetric"),  9, 10.5)
save_fig("figS_screening_profiles_none", screening_profile_plot("none"),       9, 10.5)


# Figure 4: candidate rules, both information conditions ------------------
mode_point_shapes <- c("Perfect information" = 21, "With assessment" = 17)

candidate_metrics <- performance_delta %>%
  filter(metric %in% names(metric_labels)) %>%
  mutate(om_scenario = factor(om_scenario, om_scenario_names),
         metric = factor(metric_labels[metric], metric_labels),
         hcr = factor(hcr, hcr_names))

candidate_metrics_plot <- ggplot(candidate_metrics, aes(x = hcr)) +
  geom_segment(aes(xend = hcr, y = perfect, yend = assessment),
               colour = "grey55", linewidth = 0.5) +
  geom_point(aes(y = perfect, colour = "Perfect information", shape = "Perfect information"),
             size = 2.3, stroke = 0.9, fill = "white") +
  geom_point(aes(y = assessment, colour = "With assessment", shape = "With assessment"),
             size = 2.3) +
  facet_grid(metric ~ om_scenario, scales = "free_y") +
  scale_colour_manual(values = mode_colours, name = NULL) +
  scale_shape_manual(values = mode_point_shapes, name = NULL) +
  scale_x_discrete(labels = hcr_labels) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig("fig4_candidate_metrics", candidate_metrics_plot, 9, 10)


# Figure 5: candidate time series -----------------------------------------
candidate_timeseries <- bind_rows(
  build_timeseries("SSB",   modes = "asmt") %>% mutate(metric = "SSB (t)"),
  build_timeseries("Catch", modes = "asmt") %>% mutate(metric = "Catch (t)")
) %>%
  filter(year >= terminal_conditioning_year) %>%
  mutate(metric = factor(metric, c("SSB (t)", "Catch (t)")))

depletion_threshold    <- crash_threshold(get_run("Baseline", "asmt", reference_hcr))[1]
depletion_threshold_df <- data.frame(metric = factor("SSB (t)", levels(candidate_timeseries$metric)), y = depletion_threshold)

candidate_timeseries_plot <- ggplot(candidate_timeseries, aes(year, median, colour = hcr)) +
  geom_ribbon(data = filter(candidate_timeseries, hcr == reference_hcr), aes(ymin = lower, ymax = upper),
              fill = "grey60", colour = NA, alpha = 0.25) +
  geom_hline(data = depletion_threshold_df, aes(yintercept = y), linetype = 2, colour = "grey30") +
  geom_line(linewidth = 0.7) +
  facet_grid(metric ~ om_scenario, scales = "free_y") +
  scale_colour_manual(values = hcr_colours, labels = hcr_labels, name = NULL) +
  scale_x_continuous(breaks = seq(2030, 2080, 10)) +
  labs(x = "Year", y = NULL) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

save_fig("fig5_candidate_timeseries", candidate_timeseries_plot, 10, 6)


# S2-S3: OM demographics and selectivity ---------------------------------
sim_list   <- readRDS(here("outputs", "base_mse_om.RDS"))
assessment_report   <- readRDS(here("inputs", "actual_report.rds"))
n_conditioning_years <- sim_list$feedback_start_yr
conditioning_years        <- 1975:(1975 + n_conditioning_years - 1)   # NSEI assessment starts in 1975
ages       <- 2:(sim_list$n_ages + 1)        # NSEI models ages 2 - 31+
sex_labels     <- c("Female", "Male")

demographics <- bind_rows(
  # weight-at-age (population / survey and fishery sources)
  expand_grid(sex = 1:2, source = c("Population", "Fishery")) %>%
    mutate(values = map2(sex, source, \(sex_index, source_name) {
      if (source_name == "Population") sim_list$WAA[1, 1, 1, 1, , sex_index, 1]
      else sim_list$WAA_fish[1, 1, 1, 1, , sex_index, 1, 1]
    })) %>%
    unnest_longer(values, values_to = "value") %>%
    mutate(age = rep(ages, 4), quantity = "Weight-at-age (kg)"),
  # maturity-at-age (females only; males set to 0 in the OM)
  tibble(sex = 1, source = "Population", age = ages,
         value = sim_list$MatAA[1, 1, 1, 1, , 1, 1],
         quantity = "Maturity-at-age"),
  # natural mortality: flat over years, ages and sims; season dimension is length 1 here
  expand_grid(sex = 1:2, source = "Population") %>%
    mutate(values = map(sex, \(sex_index) sim_list$natmort[1, 1, 1, 1, , sex_index, 1])) %>%
    unnest_longer(values, values_to = "value") %>%
    mutate(age = rep(ages, 2), quantity = "Natural mortality")
) %>%
  mutate(sex = factor(sex_labels[sex], levels = sex_labels),
         quantity = factor(quantity, levels = c("Weight-at-age (kg)", "Maturity-at-age",
                                                "Natural mortality")))

demographics_plot <- ggplot(demographics, aes(age, value, colour = sex, linetype = source)) +
  geom_line(linewidth = 1.2) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  facet_wrap(~ quantity, scales = "free_y", strip.position = "left") +
  labs(x = "Age", y = NULL, colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

save_fig("om_demographics", demographics_plot, 15, 6.5)

process_labels <- c("Fishery selectivity", "Survey selectivity", "Retention")
selectivity <- expand_grid(sex = 1:2, process = process_labels) %>%
  mutate(values = map2(sex, process, \(sex_index, process_name) {
    switch(process_name,
           "Fishery selectivity" = sim_list$fish_sel[1, 1, n_conditioning_years, 1, , sex_index, 1, 1],
           "Survey selectivity"  = sim_list$srv_sel[1, 1, n_conditioning_years, 1, , sex_index, 1, 1],
           "Retention"           = sim_list$ret_sel[1, 1, n_conditioning_years, 1, , sex_index, 1, 1])
  })) %>%
  unnest_longer(values, values_to = "value") %>%
  mutate(age = rep(ages, 6),
         sex = factor(sex_labels[sex], levels = sex_labels),
         process = factor(process, levels = process_labels))

selectivity_plot <- ggplot(selectivity, aes(age, value, colour = sex)) +
  geom_line(linewidth = 1.2) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  facet_wrap(~ process) +
  labs(x = "Age", y = "Proportion selected / retained", colour = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.background = element_blank())

save_fig("om_selectivity", selectivity_plot, 15, 6.5)



# S4-S5: assessment error -------------------------------------------------
relative_error <- build_timeseries("ssb_re", hcrs = hcr_names, modes = "asmt") %>% filter(!is.na(median))

relative_error_plot <- ggplot(relative_error, aes(year, median, colour = hcr)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_ribbon(data = filter(relative_error, hcr == reference_hcr),
              aes(ymin = lower, ymax = upper, fill = hcr), colour = NA, alpha = 0.20) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ om_scenario) +
  scale_colour_manual(values = hcr_colours) +
  scale_fill_manual(values = hcr_fill_colours) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = "none") +   # the interval is the reference rule's; one legend is enough
  labs(x = "Year", y = "Terminal SSB relative error", colour = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top")

save_fig("asmt_em_ssb_relative_error", relative_error_plot, 15, 8)

# What the rule saw vs what was there, for the reference rule
status_true_vs_perceived <- bind_rows(
  build_timeseries("status",           hcrs = reference_hcr, modes = "asmt") %>% mutate(series = "True (OM)"),
  build_timeseries("status_perceived", hcrs = reference_hcr, modes = "asmt") %>% mutate(series = "Perceived (EM)")
) %>%
  filter(!is.na(median)) %>%
  mutate(series = factor(series, c("True (OM)", "Perceived (EM)")))

perceived_status_plot <- ggplot(status_true_vs_perceived, aes(year, median, colour = series)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = series), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ om_scenario, scales = "free_y") +
  scale_colour_manual(values = unname(mode_colours)) +
  scale_fill_manual(values = unname(mode_colours)) +
  labs(x = "Year", y = expression(SSB / B[x]), colour = NULL, fill = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top")

save_fig("asmt_em_perceived_vs_true_status", perceived_status_plot, 15, 8)


# S8-S11: assessment vs perfect information ------------------------------
plot_mode_timeseries <- function(timeseries, y_label) {
  historical <- timeseries %>% filter(year <= terminal_conditioning_year, !is.na(median)) %>%
    distinct(om_scenario, hcr, year, .keep_all = TRUE)
  projection <- timeseries %>% filter(year >= terminal_conditioning_year, !is.na(median))

  ggplot(mapping = aes(year, median)) +
    geom_vline(xintercept = terminal_conditioning_year, linetype = 2, colour = "grey50") +
    geom_ribbon(data = historical, aes(ymin = lower, ymax = upper), colour = NA, alpha = 0.20) +
    geom_line(data = historical, linewidth = 0.9) +
    geom_ribbon(data = projection, aes(ymin = lower, ymax = upper, fill = mode),
                colour = NA, alpha = 0.15) +
    geom_line(data = projection, aes(colour = mode), linewidth = 0.9) +
    facet_grid(hcr ~ om_scenario, scales = "free_y") +
    scale_colour_manual(values = mode_colours) +
    scale_fill_manual(values = mode_colours) +
    labs(x = "Year", y = y_label, colour = NULL, fill = NULL) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw(base_size = 20) +
    theme(legend.position = "top")
}

mode_timeseries_specs <- list(
  ssb    = list(metric = "SSB",    label = "SSB (t)"),
  cat    = list(metric = "Catch",  label = "Catch (t)"),
  f      = list(metric = "Fmort",  label = "Fishing mortality (F)"),
  status = list(metric = "status", label = expression(SSB / B[x]))
)

for (name in names(mode_timeseries_specs)) {
  spec <- mode_timeseries_specs[[name]]
  mode_plot <- plot_mode_timeseries(build_timeseries(spec$metric), spec$label)
  if (name == "status") mode_plot <- mode_plot + geom_hline(yintercept = 1, linetype = 3)
  save_fig(paste0("cmp_all_hcr_", name, "_timeseries"), mode_plot, 15, 20)
}

better_when_higher <- c(med_ssb = TRUE, med_catch = TRUE, catch_aav = FALSE,
                        p_crash = FALSE, p_closure = FALSE)

performance_delta_long <- performance_delta %>%
  filter(metric %in% metric_order) %>%
  mutate(direction = ifelse(change == 0, "No change",
                     ifelse((change > 0) == better_when_higher[metric],
                            "Better with assessment", "Worse with assessment")),
         metric = factor(metric, levels = metric_order,
                         labels = c("Median SSB\n(% change)", "Median catch\n(% change)",
                                    "Catch AAV\n(% change)",
                                    "P(SSB < 10% hist. max)\n(pp change)",
                                    "P(closure)\n(pp change)")))

performance_change_plot <- ggplot(performance_delta_long, aes(change, hcr, fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_col(width = 0.65, alpha = 0.85) +
  facet_grid(om_scenario ~ metric, scales = "free_x") +
  scale_fill_manual(values = c("Better with assessment" = "#0072B2",
                               "Worse with assessment"  = "#D55E00",
                               "No change"              = "grey70")) +
  scale_y_discrete(limits = rev(hcr_names)) +
  labs(x = "Assessment run relative to perfect information", y = NULL, fill = NULL) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top")

save_fig("cmp_performance_change", performance_change_plot, 17, 10)


# S1: OM conditioning fit -------------------------------------------------
library(SPoRC)
sim_env <- Setup_sim_env(sim_list)
set.seed(123)
for (year in 1:n_conditioning_years) run_annual_cycle(year, 1, sim_env)

source_labels <- c("SPoRC OM", "Stock Assessment")
conditioning_fit <- bind_rows(
  tibble(year = conditioning_years, value = sim_env$SSB[1, 1, 1:n_conditioning_years, 1] / 1e3,   # OM kg -> t
         source = source_labels[1], quantity = "Spawning biomass (t)"),
  tibble(year = conditioning_years, value = assessment_report$tot_spawn_biom[1:n_conditioning_years] / 1e3,
         source = source_labels[2], quantity = "Spawning biomass (t)"),
  tibble(year = conditioning_years, value = sim_env$Rec[1, 1, 1:n_conditioning_years, 1] / 1e6,
         source = source_labels[1], quantity = "Age-2 recruits (millions)"),
  tibble(year = conditioning_years, value = assessment_report$pred_rec / 1e6,
         source = source_labels[2], quantity = "Age-2 recruits (millions)")
) %>%
  mutate(source = factor(source, levels = source_labels),
         quantity = factor(quantity, levels = c("Spawning biomass (t)",
                                                "Age-2 recruits (millions)")))

conditioning_fit_plot <- ggplot(conditioning_fit, aes(year, value, colour = source, linetype = source)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = setNames(c("#0D0887", "#CC4678"), source_labels)) +
  scale_linetype_manual(values = setNames(c(1, 2), source_labels)) +
  facet_wrap(~ quantity, scales = "free_y", strip.position = "left") +
  labs(x = "Year", y = NULL, colour = NULL, linetype = NULL) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

save_fig("om_conditioning_fit", conditioning_fit_plot, 14, 10)

ssb_percent_difference <- 100 * (sim_env$SSB[1, 1, 1:n_conditioning_years, 1] - assessment_report$tot_spawn_biom[1:n_conditioning_years]) /
  assessment_report$tot_spawn_biom[1:n_conditioning_years]
cat("OM vs assessment SSB through", max(conditioning_years), "- max abs diff",
    round(max(abs(ssb_percent_difference)), 2), "percent\n")
