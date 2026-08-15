# Purpose: To visualize how HCRs work in the MSE
# Creator: Matthew LH. Cheng
# Date 8/14/26


# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)

hcr_fn <- function(x, frp, brp, local_alpha = 0) {
  stock_status <- x / brp
  if (stock_status >= 1) return(frp)
  if (stock_status > local_alpha && stock_status < 1) return(frp * (stock_status - local_alpha) / (1 - local_alpha))
  return(0)
}

# work in relative terms: brp = 1 so x IS stock status (B / B40)
# illustrative F rates: higher SPR% -> lower fishing mortality
f_vals <- c("F40%" = 0.12, "F45%" = 0.10, "F50%" = 0.085,
            "F55%" = 0.072, "F60%" = 0.062)
alpha_vals <- c(0, 0.1, 0.25, 0.5)

# F and relative catch (F x B, B proportional to status) in long format
quantity_levels <- c("Fishing mortality (F)", "Relative catch")
add_catch <- function(df) {
  df %>%
    mutate(catch = F * status) %>%
    pivot_longer(c(F, catch), names_to = "quantity", values_to = "value") %>%
    mutate(quantity = factor(ifelse(quantity == "F", quantity_levels[1], quantity_levels[2]),
                             levels = quantity_levels))
}


# Effect of increasing Fx% ------------------------------------------------

# threshold rule alongside a constant F rule (F50% applied at all stock statuses)
rule_levels <- c("Threshold", "Constant F")
fx_df <- expand_grid(status = seq(0, 1.5, 0.005), fx = names(f_vals),
                     rule = rule_levels) %>%
  filter(rule == rule_levels[1] | fx == "F50%") %>%
  mutate(F = ifelse(rule == "Threshold",
                    map2_dbl(status, f_vals[fx], \(s, f) hcr_fn(s, f, brp = 1)),
                    f_vals[fx]),
         fx = factor(fx, levels = names(f_vals)),
         rule = factor(rule, levels = rule_levels)) %>%
  add_catch()

p_fx <- ggplot(fx_df, aes(status, value, colour = fx, linetype = rule)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 1.2) +
  scale_colour_viridis_d(option = "plasma", end = 0.85) +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = expression(SSB / B[x]), y = NULL, colour = NULL, linetype = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

print(p_fx)
ggsave(here("figs", "hcr_shape_fx.png"), p_fx, width = 13, height = 10, dpi = 300)


# F50% threshold HCR alone ------------------------------------------------

f50_df <- tibble(status = seq(0, 1.5, 0.005)) %>%
  mutate(F = map_dbl(status, \(s) hcr_fn(s, 0.085, brp = 1))) %>%
  add_catch()

p_f50 <- ggplot(f50_df, aes(status, value)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 1.2, colour = "#0D0887") +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = expression(SSB / B[x]), y = NULL) +
  theme_bw(base_size = 20) +
  theme(strip.placement = "outside", strip.background = element_blank())

print(p_f50)
ggsave(here("figs", "hcr_shape_f50.png"), p_f50, width = 13, height = 10, dpi = 300)


# F50% threshold vs constant F --------------------------------------------

f50_rule_df <- expand_grid(status = seq(0, 1.5, 0.005), rule = rule_levels) %>%
  mutate(F = ifelse(rule == "Threshold",
                    map_dbl(status, \(s) hcr_fn(s, 0.085, brp = 1)),
                    0.085),
         rule = factor(rule, levels = rule_levels)) %>%
  add_catch()

p_f50_rule <- ggplot(f50_rule_df, aes(status, value, linetype = rule)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 1.2, colour = "#0D0887") +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = expression(SSB / B[x]), y = NULL, linetype = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

print(p_f50_rule)
ggsave(here("figs", "hcr_shape_f50_rules.png"), p_f50_rule, width = 13, height = 10, dpi = 300)


# Effect of increasing alpha ----------------------------------------------

# F plateau fixed at F50; curves pivot from (1, F50) to an x-intercept at alpha
alpha_df <- expand_grid(status = seq(0, 1.5, 0.005), a = alpha_vals) %>%
  mutate(F = map2_dbl(status, a, \(s, a) hcr_fn(s, 0.085, brp = 1, local_alpha = a)),
         a = factor(a)) %>%
  add_catch()

p_alpha <- ggplot(alpha_df, aes(status, value, colour = a)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 1.2) +
  scale_colour_viridis_d(option = "plasma", end = 0.85) +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = expression(SSB / B[x]), y = NULL, colour = expression(alpha)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

print(p_alpha)
ggsave(here("figs", "hcr_shape_alpha.png"), p_alpha, width = 13, height = 10, dpi = 300)


# Both: Fx% x alpha -------------------------------------------------------

both_df <- expand_grid(status = seq(0, 1.5, 0.005),
                       fx = names(f_vals),
                       a = alpha_vals) %>%
  mutate(F = pmap_dbl(list(status, f_vals[fx], a),
                      \(s, f, a) hcr_fn(s, f, brp = 1, local_alpha = a)),
         fx = factor(fx, levels = names(f_vals)),
         a = factor(a)) %>%
  add_catch()

p_both <- ggplot(both_df, aes(status, value, colour = a)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 1.1) +
  scale_colour_viridis_d(option = "plasma", end = 0.85) +
  facet_grid(quantity ~ fx, scales = "free_y", switch = "y") +
  labs(x = expression(SSB / B[x]), y = NULL, colour = expression(alpha)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background.y = element_blank())

print(p_both)
ggsave(here("figs", "hcr_shape_both.png"), p_both, width = 16, height = 8, dpi = 300)
