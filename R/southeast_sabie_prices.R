# Derive sablefish prices by grade for 2024 values

library(tidyverse)
library(mgcv)

data <- readRDS(here("inputs", 'sablefish_data_19May2026.RDS'))

# Grade-specific nominal prices from Stephen Rhoads
grade_nominal <- tribble(
  ~grade,     ~`2019`, ~`2020`, ~`2021`, ~`2022`, ~`2023`, ~`2024`,
  "#1 1/2",     1.59,    0.41,    1.21,    1.01,    0.47,    0.32,
  "#1 2/3",     2.00,    1.31,    2.10,    2.08,    1.21,    1.15,
  "#1 3/4",     2.95,    2.05,    2.52,    2.79,    1.79,    1.65,
  "#1 4/5",     3.95,    2.45,    3.01,    4.19,    2.53,    2.45,
  "#1 5/7",     6.40,    3.70,    4.03,    6.66,    5.68,    3.92,
  "#1 7+",      7.64,    5.79,    5.93,    7.80,    6.90,    6.41,
  "Grand Total", 4.22,   2.74,    3.06,    3.94,    2.92,    2.52
)

#  Aggregate nominal price (1984-2024) from ADFG prices https://www.adfg.alaska.gov/index.cfm?adfg=CommercialByFisheryGroundfish.groundfish_exvessel_byarea
# Inflation adjustment from https://fred.stlouisfed.org/series/CPIAUCNS
price_dat <- tibble(
  year = 1984:2024,
  price_nominal = c(
    0.62, 1.02, 1.06, 1.32, 1.62, 1.48, 1.33, 1.91, 1.94, 1.65,
    2.74, 3.11, 3.40, 3.82, 2.53, 3.17, 3.69, 3.26, 3.40, 3.67,
    3.22, 3.50, 3.96, 4.02, 4.69, 4.98, 5.88, 7.90, 6.51, 4.47,
    5.63, 6.03, 6.69, 7.50, 5.56, 4.50, 2.79, 3.25, 3.64, 2.58,
    2.89
  ),
  cpi = c(
    103.9, 107.6, 109.6, 113.6, 118.3, 124.0, 130.7, 136.2, 140.3, 144.5,
    148.2, 152.4, 156.9, 160.5, 163.0, 166.6, 172.2, 177.1, 179.9, 184.0,
    188.9, 195.3, 201.6, 207.3, 215.3, 214.5, 218.1, 224.9, 229.6, 233.0,
    236.7, 237.0, 240.0, 245.1, 251.1, 255.7, 258.8, 271.0, 292.7, 304.7,
    313.7
  )
)

cpi_2024 <- 313.7

# Relativize grade data
grade_long <- grade_nominal %>%
  filter(grade != 'Grand Total') %>%
  mutate(grade = str_remove(grade, "#1 ")) %>%
  pivot_longer(-grade, names_to = "year", values_to = "price_nominal") %>%
  mutate(year = as.integer(year)) %>%
  group_by(year) %>%
  mutate(rel_price = price_nominal / max(price_nominal)) %>%
  group_by(grade) %>%
  summarize(rel_price = mean(rel_price))

# Look at relative price grade info
plot(grade_long$rel_price, type = 'l')

# Get inflation adjusted price / lb
price_dat$price_inf_adj <- price_dat$price_nominal * (cpi_2024 / price_dat$cpi)
price_dat$rel_price_inf_adj <- price_dat$price_inf_adj / max(price_dat$price_inf_adj)

# Plot metrics
plot(price_dat$price_inf_adj, price_dat$price_nominal)

# Join grades by year to NAA, spawning biomass, total biomass
pop_dat <- data.frame(year = 1984:2024,
                      rel_price_inf_adj = price_dat$rel_price_inf_adj,
                      ssb = rep$tot_spawn_biom[which(1975:2025 %in% c(1984:2024))],
                      tot_biom = rep$tot_biom[which(1975:2025 %in% c(1984:2024))],
                      n = apply(rep$N, 1, sum)[which(1975:2025 %in% c(1984:2024))])

# Plot
ggplot(pop_dat, aes(x = ssb, y = rel_price_inf_adj, label = year )) +
  geom_text() +
  geom_smooth() +
  labs(x = 'Spawning Biomass', y = 'Relative Price Adjusted by Inflation')

ggplot(pop_dat, aes(x = tot_biom, y = rel_price_inf_adj, label = year )) +
  geom_text() +
  geom_smooth() +
  labs(x = 'Total Biomass', y = 'Relative Price Adjusted by Inflation')

ggplot(pop_dat, aes(x = n, y = rel_price_inf_adj, label = year )) +
  geom_text() +
  geom_smooth() +
  labs(x = 'Total N', y = 'Relative Price Adjusted by Inflation')

# Fit model
pop_dat$rel_price_inf_adj[pop_dat$rel_price_inf_adj == 1] <- 0.99
mod <- gam(rel_price_inf_adj ~ s(n), data = pop_dat, family = betar)
summary(mod)

# Convert
waa_male_lb   <- data$data_fsh_waa[1,,1] * 2.20462
waa_female_lb <- data$data_fsh_waa[1,,2] * 2.20462

# Grade weight boundaries (round weight, lbs)
grade_breaks <- c(0, 2, 3, 5, 7, 8, 10, Inf)
grade_labels <- c("No grade", "1/2", "2/3", "3/4", "4/5", "5/7", "7+")

ages <- 2:31
grade_male   <- cut(waa_male_lb, breaks = grade_breaks, labels = grade_labels, right = FALSE)
grade_female <- cut(waa_female_lb, breaks = grade_breaks, labels = grade_labels, right = FALSE)

grade_age_df <- data.frame(
  age = ages,
  wt_male_lb = round(waa_male_lb, 2),
  grade_male = grade_male,
  wt_female_lb = round(waa_female_lb, 2),
  grade_female = grade_female
)

grade_age_df <- grade_age_df %>%
  left_join(grade_long, by = c("grade_male" = "grade")) %>%
  rename(rel_price_male = rel_price) %>%
  left_join(grade_long, by = c("grade_female" = "grade")) %>%
  rename(rel_price_female = rel_price)

saveRDS(mod, here("outputs", 'price_model.RDS'))
write.csv(grade_age_df, here("outputs", 'rel_grade_price.csv'))
