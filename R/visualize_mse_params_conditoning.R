# Purpose: To visualize how the MSE is conditioned and parameterizations
# Creator: Matthew LH. Cheng
# Date 8/14/26


# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)

sim_list <- readRDS(here("outputs", "base_mse_om.RDS"))
rep <- readRDS(here("inputs", "actual_report.rds"))
data <- readRDS(here("inputs", "sablefish_data_19May2026.RDS"))

n_yrs_cond <- sim_list$feedback_start_yr # conditioning period (assessment years)
yrs <- 1975:(1975 + n_yrs_cond - 1) # NSEI assessment starts in 1975
ages <- 2:(sim_list$n_ages + 1) # NSEI models ages 2 - 31+
sex_levels <- c("Female", "Male")


# Demographic rates -------------------------------------------------------

# demographics are time- and sim-invariant, so pull from year 1, sim 1
demo_df <- bind_rows(
  # weight-at-age (population / survey and fishery sources)
  expand_grid(sex = 1:2, source = c("Population", "Fishery")) %>%
    mutate(vals = map2(sex, source, \(s, src) {
      if (src == "Population") sim_list$WAA[1, 1, 1, 1, , s, 1]
      else sim_list$WAA_fish[1, 1, 1, 1, , s, 1, 1]
    })) %>%
    unnest_longer(vals, values_to = "value") %>%
    mutate(age = rep(ages, 4), quantity = "Weight-at-age (kg)"),
  # maturity-at-age (females only; males set to 0 in the OM)
  tibble(sex = 1, source = "Population", age = ages,
         value = sim_list$MatAA[1, 1, 1, 1, , 1, 1],
         quantity = "Maturity-at-age"),
  # natural mortality
  expand_grid(sex = 1:2, age = ages) %>%
    mutate(source = "Population",
           value = map2_dbl(age, sex, \(a, s) sim_list$natmort[1, 1, 1, 1, a - 1, s, 1]), # season 1
           quantity = "Natural mortality")
) %>%
  mutate(sex = factor(sex_levels[sex], levels = sex_levels),
         quantity = factor(quantity, levels = c("Weight-at-age (kg)", "Maturity-at-age", "Natural mortality")))

p_demo <- ggplot(demo_df, aes(age, value, colour = sex, linetype = source)) +
  geom_line(linewidth = 1.2) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  facet_wrap(~ quantity, scales = "free_y", strip.position = "left") +
  labs(x = "Age", y = NULL, colour = NULL, linetype = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

print(p_demo)
ggsave(here("figs", "om_demographics.png"), p_demo, width = 15, height = 6.5, dpi = 300)


# Selectivity and retention -----------------------------------------------

# fishery selectivity varies by time block; survey and retention are constant.
# show the terminal-year (current) values used in the projection period.
sel_df <- expand_grid(sex = 1:2,
                      process = c("Fishery selectivity", "Survey selectivity", "Retention")) %>%
  mutate(vals = map2(sex, process, \(s, p) {
    switch(p,
           "Fishery selectivity" = sim_list$fish_sel[1, 1, n_yrs_cond, 1, , s, 1, 1],
           "Survey selectivity"  = sim_list$srv_sel[1, 1, n_yrs_cond, 1, , s, 1, 1],
           "Retention"           = sim_list$ret_sel[1, 1, n_yrs_cond, 1, , s, 1, 1])
  })) %>%
  unnest_longer(vals, values_to = "value") %>%
  mutate(age = rep(ages, 6),
         sex = factor(sex_levels[sex], levels = sex_levels),
         process = factor(process, levels = c("Fishery selectivity", "Survey selectivity", "Retention")))

p_sel <- ggplot(sel_df, aes(age, value, colour = sex)) +
  geom_line(linewidth = 1.2) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  facet_wrap(~ process) +
  labs(x = "Age", y = "Proportion selected / retained", colour = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.background = element_blank())

print(p_sel)
ggsave(here("figs", "om_selectivity.png"), p_sel, width = 15, height = 6.5, dpi = 300)


# Conditioning fit --------------------------------------------------------

# run the OM through the conditioning period (deterministic here: F, recruitment,
# and demographics are all fixed inputs, so a single simulation suffices) and
# compare OM dynamics against the assessment estimates and observed data
sim_env <- Setup_sim_env(sim_list)
set.seed(123)
for (y in 1:n_yrs_cond) run_annual_cycle(y, 1, sim_env)

srv_yrs <- yrs[data$yrs_srv_cpue + 1] # TMB indices are 0-based
src_levels <- c("SPoRC OM", "Stock Assessment")

cond_df <- bind_rows(
  # spawning biomass (OM in kg -> t)
  tibble(year = yrs, value = sim_env$SSB[1, 1, 1:n_yrs_cond, 1] / 1e3,
         source = src_levels[1], quantity = "Spawning biomass (t)"),
  tibble(year = yrs, value = rep$tot_spawn_biom[1:n_yrs_cond] / 1e3,
         source = src_levels[2], quantity = "Spawning biomass (t)"),
  # recruitment
  tibble(year = yrs, value = sim_env$Rec[1, 1, 1:n_yrs_cond, 1] / 1e6,
         source = src_levels[1], quantity = "Age-2 recruits (millions)"),
  tibble(year = yrs, value = rep$pred_rec / 1e6,
         source = src_levels[2], quantity = "Age-2 recruits (millions)"),
) %>%
  mutate(source = factor(source, levels = src_levels),
         quantity = factor(quantity, levels = c("Spawning biomass (t)",
                                                "Age-2 recruits (millions)")))

p_cond <- ggplot(cond_df, aes(year, value, colour = source, linetype = source)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = setNames(c("#0D0887", "#CC4678"), src_levels)) +
  scale_linetype_manual(values = setNames(c(1, 2), src_levels)) +
  facet_wrap(~ quantity, scales = "free_y", strip.position = "left") +
  labs(x = "Year", y = NULL, colour = NULL, linetype = NULL) +
  coord_cartesian(ylim = c(0,NA)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

print(p_cond)
ggsave(here("figs", "om_conditioning_fit.png"), p_cond, width = 14, height = 10, dpi = 300)

# quick summary of how well the OM reproduces the assessment
ssb_pdiff <- 100 * (sim_env$SSB[1, 1, 1:n_yrs_cond, 1] - rep$tot_spawn_biom[1:n_yrs_cond]) / rep$tot_spawn_biom[1:n_yrs_cond]
ssb_pdiff
