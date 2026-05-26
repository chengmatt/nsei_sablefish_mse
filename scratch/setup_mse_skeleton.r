# Purpose: To setup MSE skeleton and validate MSE
# Creator: Matthew LH. Cheng
# Date: 5/20/26

# Setup -------------------------------------------------------------------
set.seed(123) # set seed

# Load in libraries
library(here)
library(SPoRC)
library(furrr)
library(future)
library(progressr)

# Source functions
source(here("R", 'functions', 'utils.R'))
source(here("R", 'functions', 'setup_em.R'))

# Scenario Grid -----------------------------------------------------------
# Factor levels
spr_x_vals    <- seq(0.3, 0.7, by = 0.1)
alpha_vals    <- c(0, 0.25, 0.5)
brp_types     <- c("none", "threshold")
dyn_b0_vals   <- c(FALSE, TRUE)
stability_opts <- c("none", "symmetric", "asymmetric", 'oneway')

# Stability Constraints
stab_pars <- list(
  none       = list(max_increase = Inf, max_decrease = Inf),
  symmetric  = list(max_increase = 0.15, max_decrease = 0.15),
  asymmetric = list(max_increase = 0.15, max_decrease = 0.25), # fast down, slow up
  oneway = list(max_increase = 0.15, max_decrease = Inf) # one way constraint
)

# Build full grid
scenario_grid_full <- expand.grid(
  spr_x     = spr_x_vals,
  alpha     = alpha_vals,
  brp_type  = brp_types,
  dyn_b0    = dyn_b0_vals,
  stability = stability_opts,
  stringsAsFactors = FALSE
)

# When brp_type == "none": constant F at F_SPR_x, no HCR threshold
# alpha is irrelevant (force to 0)
# dyn_b0 is irrelevant (force to FALSE)
scenario_grid <- scenario_grid_full
idx_none <- scenario_grid$brp_type == "none"
scenario_grid$alpha[idx_none]  <- 0
scenario_grid$dyn_b0[idx_none] <- FALSE
scenario_grid <- dplyr::distinct(scenario_grid)

# Add scenario ID
scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))

run_single_sim <- function(sim, sim_list, scenario, use_true_values) {

  library(SPoRC)
  source(here::here("R", "functions", "utils.R"))
  source(here::here("R", "functions", "setup_em.R"))

  # Unpack scenario
  spr_x       <- scenario$spr_x
  alpha       <- scenario$alpha
  brp_type    <- scenario$brp_type
  use_dyn_b0  <- scenario$dyn_b0
  stab_type   <- scenario$stability

  # Stability parameters
  stab <- stab_pars[[stab_type]]

  # Select HCR function
  if (brp_type == "none") {
    hcr_fn <- function(x, frp, brp, alpha = 0) return(frp)
  } else {
    hcr_fn <- function(x, frp, brp, local_alpha = alpha) {
      stock_status <- x / brp
      if (stock_status >= 1) return(frp)
      if (stock_status > local_alpha && stock_status < 1) return(frp * (stock_status - local_alpha) / (1 - local_alpha))
      return(0)
    }
  }

  # Reference point options
  reference_points_opt <- list(
    n_avg_yrs    = 1,
    SPR_x        = spr_x,
    calc_rec_st_yr = 1,
    rec_age      = 0,
    type         = "single_region",
    what         = "SPR"
  )

  # Projection options
  proj_opt <- list(
    n_proj_yrs     = 2,
    HCR_function   = hcr_fn,
    recruitment_opt = "mean_rec",
    fmort_opt      = "HCR",
    bh_rec_opt     = NULL
  )

  # Create sim environment
  sim_env <- Setup_sim_env(sim_list)

  n_yrs    <- sim_env$n_yrs
  asmt_yrs <- sim_env$feedback_start_yr:n_yrs
  n_asmt   <- length(asmt_yrs)

  # Storage
  ssb_em_list <- vector("list", n_asmt)
  pdHess_vec  <- rep(NA, n_asmt)
  grad_vec    <- rep(NA, n_asmt)
  bx_vec      <- rep(NA, n_asmt)

  for (y in 1:n_yrs) {

    # Run Annual Cycle
    run_annual_cycle(y, sim, sim_env)

    if (y >= sim_env$feedback_start_yr) {

      asmt_idx  <- y - sim_env$feedback_start_yr + 1
      local_y   <- y
      local_sim <- sim

      # Run Assessment
      if (!use_true_values) {
        asmt_result <- tryCatch({

          # Setup data lists and model
          tmp_list  <- setup_em(sim_env, local_y, local_sim)
          asmt_data <- tmp_list$data

          obj <- fit_model(
            asmt_data, tmp_list$par, tmp_list$map, NULL,
            newton_loops = 0, silent = TRUE
          )
          sd_rep <- sdreport(obj)

          list(
            success = TRUE,
            ssb_em  = obj$rep$SSB[1:local_y],
            pdHess  = sd_rep$pdHess,
            grad    = max(abs(sd_rep$gradient.fixed))
          )
        }, error = function(e) {
          message(sprintf("Sim %d, year %d failed: %s", local_sim, local_y, e$message))
          list(success = FALSE, ssb_em = NULL, pdHess = NA, grad = NA)
        })
      } else {
        asmt_result <- list(success = FALSE, ssb_em = NULL, pdHess = NA, grad = NA)
      }

      # Get Reference Points
      # input 1 so get_closed_loop_reference_points still optimizes and get_proj_catch works
      tmp_sim_env <- sim_env
      if (sim_env$Fmort[,y,,,sim] == 0) tmp_sim_env$Fmort[,y,,,sim] <- 1
      tmp_obj <- if (use_true_values) NULL else obj
      if (!use_true_values && !is.null(tmp_obj) && tmp_obj$rep$Fmort[,y,,] == 0) {
        tmp_obj$rep$Fmort[] <- 1
      }

      reference_points <- get_closed_loop_reference_points(
        use_true_values    = use_true_values,
        sim_env            = tmp_sim_env,
        asmt_data          = if (use_true_values) NULL else asmt_data,
        asmt_rep           = tmp_obj$rep,
        y                  = local_y,
        sim                = local_sim,
        reference_points_opt = reference_points_opt,
        n_proj_yrs         = proj_opt$n_proj_yrs
      )

      # Get Catch Advice
      tmp_catch <- get_proj_catch(
        tmp_obj, if (use_true_values) NULL else asmt_data,
        proj_opt, reference_points,
        sim_env, local_y, local_sim,
        use_true_values,
        # dyn_bx follows spr_x when dynamic B0 is used
        use_dyn_b0 = use_dyn_b0,
        dyn_bx     = if (use_dyn_b0) spr_x else spr_x
      )

      # Apply Stability Constraint (only when catch from the base HCR is not 0)
      if(tmp_catch != 0) {
        constrained_catch <- apply_stability(
          new_catch    = as.vector(tmp_catch),
          prev_catch   = sim_env$TrueCatch[1,y-1,1,1,sim], # use previous years catch
          max_increase = stab$max_increase,
          max_decrease = stab$max_decrease
        )
      }

      # Catch to F
      if (local_y < n_yrs) {
        if(constrained_catch != 0) catch_to_f(constrained_catch, sim_env, local_y, local_sim)
        else sim_env$Fmort[,y+1,,,sim] <- 0
      }

      # Storage
      ssb_em_list[[asmt_idx]] <- asmt_result$ssb_em
      pdHess_vec[asmt_idx]    <- asmt_result$pdHess
      grad_vec[asmt_idx]      <- asmt_result$grad
      bx_vec[asmt_idx]        <- reference_points$b_ref_pt[,,1]

    } # end feedback
  } # end y

  # Return Lists
  list(
    scenario_id = scenario$scenario_id,
    om = list(
      SSB   = sim_env$SSB[,,,sim],
      Catch = sim_env$TrueCatch[,,,,sim],
      CAA   = sim_env$CAA[,,,,,,,sim],
      CAL   = sim_env$CAL[,,,,,,,sim],
      NAA   = sim_env$NAA[,,,,,,sim],
      Rec   = sim_env$Rec[,,,sim],
      Fmort = sim_env$Fmort[,,,,sim]
    ),
    ssb_em      = ssb_em_list,
    pdHess      = pdHess_vec,
    grad        = grad_vec,
    bx          = bx_vec
  )
}


# Run Base MSE -----------------------------------------------------------------
options(future.globals.maxSize = 5 * 1024^3)
plan(multisession, workers = availableCores() - 3)

# Read in MSE OM
sim_list <- readRDS(here("outputs", 'base_mse_om.RDS'))

# Run each scenario
all_results <- list()

for (sc in seq_len(nrow(scenario_grid))) {

  scenario <- scenario_grid[sc, ]
  cat(sprintf(
    "\n Scenario %d/%d: SPR=%.1f, alpha=%.2f, BRP=%s, dynB0=%s, stab=%sn",
    sc, nrow(scenario_grid),
    scenario$spr_x, scenario$alpha, scenario$brp_type,
    scenario$dyn_b0, scenario$stability
  ))

  start_time <- Sys.time()

  results <- future_map(
    1:sim_list$n_sims,
    ~ run_single_sim(.x, sim_list, scenario, use_true_values = TRUE),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Completed in %.1f min\n", elapsed))

  all_results[[sc]] <- list(
    scenario = scenario,
    results  = results
  )
}

plan(sequential)

# Save
saveRDS(all_results, here("scratch", "base_results.RDS"))

# Run Regime MSE -----------------------------------------------------------------
options(future.globals.maxSize = 5 * 1024^3)
plan(multisession, workers = availableCores() - 3)

# Read in MSE OM
sim_list <- readRDS(here("outputs", 'bh_regime_mse_om.RDS'))

# Run each scenario
all_results <- list()

for (sc in seq_len(nrow(scenario_grid))) {

  scenario <- scenario_grid[sc, ]
  cat(sprintf(
    "\n===== Scenario %d/%d: SPR=%.1f, alpha=%.2f, BRP=%s, dynB0=%s, stab=%s =====\n",
    sc, nrow(scenario_grid),
    scenario$spr_x, scenario$alpha, scenario$brp_type,
    scenario$dyn_b0, scenario$stability
  ))

  start_time <- Sys.time()

  results <- future_map(
    1:sim_list$n_sims,
    ~ run_single_sim(.x, sim_list, scenario, use_true_values = TRUE),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Completed in %.1f min\n", elapsed))

  all_results[[sc]] <- list(
    scenario = scenario,
    results  = results
  )
}

plan(sequential)

# Save
saveRDS(all_results, here("scratch", "regime_results.RDS"))

# Run Crash MSE -----------------------------------------------------------------
options(future.globals.maxSize = 5 * 1024^3)
plan(multisession, workers = availableCores() - 3)

# Read in MSE OM
sim_list <- readRDS(here("outputs", 'bh_crash_mse_om.RDS'))

# Run each scenario
all_results <- list()

for (sc in seq_len(nrow(scenario_grid))) {

  scenario <- scenario_grid[sc, ]
  cat(sprintf(
    "\n===== Scenario %d/%d: SPR=%.1f, alpha=%.2f, BRP=%s, dynB0=%s, stab=%s =====\n",
    sc, nrow(scenario_grid),
    scenario$spr_x, scenario$alpha, scenario$brp_type,
    scenario$dyn_b0, scenario$stability
  ))

  start_time <- Sys.time()

  results <- future_map(
    1:sim_list$n_sims,
    ~ run_single_sim(.x, sim_list, scenario, use_true_values = TRUE),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Completed in %.1f min\n", elapsed))

  all_results[[sc]] <- list(
    scenario = scenario,
    results  = results
  )
}

plan(sequential)

# Save
saveRDS(all_results, here("scratch", "crash_results.RDS"))


# Visualize ---------------------------------------------------------------
library(vegan)
library(ggplot2)
library(here)

# Read in pricing stuff
price_mod <- readRDS(here("outputs", 'price_model.RDS'))
grade_df <- read.csv(here("outputs", 'rel_grade_price.csv'))

om_scenarios <- list(
  list(name = "Baseline",
       results = readRDS(here('scratch', 'base_results.RDS')),
       sim_list = readRDS(here("outputs", 'base_mse_om.RDS'))),
  list(name = "Regime",
       results = readRDS(here('scratch', 'regime_results.RDS')),
       sim_list = readRDS(here("outputs", 'bh_regime_mse_om.RDS'))),
  list(name = "Crash",
       results = readRDS(here('scratch', 'crash_results.RDS')),
       sim_list = readRDS(here("outputs", 'bh_crash_mse_om.RDS')))
)

results_list <- lapply(om_scenarios, function(om) {

  all_results <- om$results
  sim_list <- om$sim_list

  perf <- do.call(rbind, lapply(1:length(all_results), function(i) {
    i = 1
    sc <- all_results[[i]]$scenario
    ssb_mat <- sapply(all_results[[i]]$results, function(x) x$om$SSB)
    catch_mat <- sapply(all_results[[i]]$results, function(x) x$om$Catch)
    catch_mat[catch_mat == 0] <- 1e-3 # add 0 to guard against aav calcs
    naa_array <- abind::abind(lapply(all_results[[i]]$results, function(x) x$om$NAA), along = 4)
    caa_array <- abind::abind(lapply(all_results[[i]]$results, function(x) x$om$CAA), along = 4)
    fb <- sim_list$feedback_start_yr
    proj_rows <- fb:nrow(ssb_mat)
    ssb_threshold <- max(all_results[[i]]$results[[1]]$om$SSB[seq_len(fb - 1)]) * 0.1


    plot(apply(ssb_mat[proj_rows, ], 1, median), ylim = c(0,1e6))
    plot(apply(catch_mat[proj_rows, ], 1, median), ylim = c(0,1e5))


    data.frame(
      scenario_id = sc$scenario_id, spr_x = sc$spr_x, alpha = sc$alpha,
      brp_type = sc$brp_type, dyn_b0 = sc$dyn_b0, stability = sc$stability,
      med_ssb = median(ssb_mat[proj_rows, ]),
      mean_age = mean(apply(naa_array[proj_rows,,,,drop = F], c(1, 4), function(mat) {
        n_total <- sum(mat)
        if (n_total == 0) return(NA)
        sum(mat * 2:31) / n_total
      })),
      pielou_j = mean(apply(naa_array[proj_rows,,,,drop = F], c(1, 4), function(mat) {
        n_total <- sum(mat)
        if (n_total == 0) return(NA)
        p <- mat / n_total
        p <- p[p > 0]
        H <- -sum(p * log(p))
        H / log(length(p))
      })),
      econ_value = mean(sapply(1:ncol(ssb_mat), function(j) {
        proj_vals <- sapply(proj_rows, function(y) {
          n_total <- sum(naa_array[y,,,j]) # get total N
          rel_price_base <- as.numeric(predict(price_mod, newdata = data.frame(n = n_total), type = "response")) # predict price based on N
          caa_female <- caa_array[y,,1,j] # get female catch-at-age
          caa_male   <- caa_array[y,,2,j] # get male catch-at-age
          # multiply, grade relative pricing, and predicted relative pricing
          value_male   <- sum(caa_male * grade_age_df$rel_price_male * rel_price_base, na.rm = TRUE)
          value_female <- sum(caa_female * grade_age_df$rel_price_female * rel_price_base, na.rm = TRUE)
          value_male + value_female
        })
        median(proj_vals)
      })),
      p_crash = mean(ssb_mat[proj_rows, ] < ssb_threshold),
      p_zero_cat = mean(catch_mat[proj_rows, ] == 0),
      med_catch = median(catch_mat[proj_rows, ]),
      catch_aav = median(sapply(1:ncol(catch_mat), function(j) {
        cc <- catch_mat[proj_rows, j]
        mean(abs(diff(cc)) / head(cc, -1), na.rm = TRUE)
      }))
    )
  }))

  # Include p_crash for non-baseline scenarios
  if (om$name == "Crash") {
    use_metrics <- c("med_ssb", "med_catch", "catch_aav", "p_crash",
                     "mean_age", "pielou_j", 'econ_value')
  } else {
    use_metrics <- c("med_ssb", "med_catch", "catch_aav", "mean_age", "pielou_j", 'econ_value')
  }

  perf_scaled <- scale(as.matrix(perf[, use_metrics]))
  set.seed(42)
  nmds <- metaMDS(dist(perf_scaled), k = 2, trymax = 200, autotransform = FALSE)
  ef <- envfit(nmds, perf_scaled, permutations = 999)

  vec_df <- as.data.frame(ef$vectors$arrows * sqrt(ef$vectors$r))
  vec_df$label <- rownames(vec_df)
  arrow_scale <- 0.6 * max(abs(nmds$points))
  vec_df$NMDS1 <- vec_df$NMDS1 * arrow_scale
  vec_df$NMDS2 <- vec_df$NMDS2 * arrow_scale
  vec_df$om_scenario <- om$name

  df <- data.frame(NMDS1 = nmds$points[,1], NMDS2 = nmds$points[,2],
                   perf, om_scenario = om$name)
  dists <- sqrt(df$NMDS1^2 + df$NMDS2^2)
  df$label <- ifelse(dists > quantile(dists, 0.90), df$scenario_id, NA)

  list(points = df, vectors = vec_df)
})

nmds_df <- do.call(rbind, lapply(results_list, `[[`, "points"))
vec_df  <- do.call(rbind, lapply(results_list, `[[`, "vectors"))

nmds_df$om_scenario <- factor(nmds_df$om_scenario, levels = c("Baseline", "Regime", "Crash"))
vec_df$om_scenario  <- factor(vec_df$om_scenario, levels = c("Baseline", "Regime", "Crash"))

nmds_df$brp_type  <- factor(nmds_df$brp_type, levels = c("none", "threshold"),
                            labels = c("Constant F", "Threshold HCR"))
nmds_df$dyn_b0    <- factor(nmds_df$dyn_b0, levels = c(FALSE, TRUE),
                            labels = c("Static B0", "Dynamic B0"))
nmds_df$stability <- factor(nmds_df$stability,
                            levels = c("none", "symmetric", "asymmetric", 'oneway'),
                            labels = c("No Constraint", "Symmetric", "Asymmetric", "OneWay"))

pos <- position_jitter(width = 0.15, height = 0.15, seed = 42)
library(dplyr)

hull_df <- nmds_df %>%
  group_by(om_scenario, alpha) %>%
  slice(chull(NMDS1, NMDS2))

pos <- position_jitter(width = 0.15, height = 0.15, seed = 42)

ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  # geom_polygon(data = hull_df, aes(linetype = factor(alpha), group = alpha),
               # fill = NA, colour = "grey50", linewidth = 0.85) +
  geom_point(aes(colour = stability, shape = dyn_b0, size = spr_x),
             alpha = 0.5, position = pos) +
  geom_text(aes(label = label), size = 2.5, vjust = -1, colour = "grey30",
            na.rm = TRUE, position = pos) +
  geom_segment(data = vec_df,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "grey40", linewidth = 0.5, inherit.aes = FALSE) +
  ggrepel::geom_label_repel(data = vec_df,
                            aes(x = NMDS1, y = NMDS2, label = label),
                            size = 3, colour = "grey30", fontface = "italic",
                            fill = "white", label.size = 0.2,
                            nudge_x = vec_df$NMDS1 * 0.3,
                            nudge_y = vec_df$NMDS2 * 0.3,
                            segment.color = "grey70", segment.size = 0.3,
                            box.padding = 0.4, point.padding = 0.2,
                            min.segment.length = 0, inherit.aes = FALSE) +
  facet_wrap(~ om_scenario) +
  scale_shape_manual(values = c("Static B0" = 16, "Dynamic B0" = 17)) +
  scale_size_continuous(range = c(1.5, 5), breaks = c(0.3, 0.5, 0.7)) +
  scale_linetype_discrete() +
  labs(colour = "Stability Type", shape = "Reference Point",
       size = expression(SPR[x]), linetype = expression(alpha)) +
  theme_bw(base_size = 14) +
  coord_cartesian(ylim = c(-7.5, 7.5), xlim = c(-7.5, 7.5))
