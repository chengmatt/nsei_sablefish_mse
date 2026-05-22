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
    HCR_alpha      = alpha,
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
  catch_vec   <- rep(NA, n_asmt)  # prescribed catch (post-stability)

  # Run MSE
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
      reference_points <- get_closed_loop_reference_points(
        use_true_values    = use_true_values,
        sim_env            = sim_env,
        asmt_data          = if (use_true_values) NULL else asmt_data,
        asmt_rep           = if (use_true_values) NULL else obj$rep,
        y                  = local_y,
        sim                = local_sim,
        reference_points_opt = reference_points_opt,
        n_proj_yrs         = proj_opt$n_proj_yrs
      )

      # Get Catch Advice
      tmp_catch <- get_proj_catch(
        obj, if (use_true_values) NULL else asmt_data,
        proj_opt, reference_points,
        sim_env, local_y, local_sim,
        use_true_values,
        # dyn_bx follows spr_x when dynamic B0 is used
        use_dyn_b0 = use_dyn_b0,
        dyn_bx     = if (use_dyn_b0) spr_x else spr_x
      )

      # Apply Stability Constraint
      constrained_catch <- apply_stability(
        new_catch    = as.vector(tmp_catch),
        prev_catch   = sim_env$TrueCatch[1,y-1,1,1,sim], # use previous years catch
        max_increase = stab$max_increase,
        max_decrease = stab$max_decrease
      )

      # Catch to F (Provides Small F (~1e-9 if catch = 0))
      if (local_y < n_yrs) catch_to_f(constrained_catch, sim_env, local_y, local_sim)

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

# Load all three OM scenarios
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

# Extract performance metrics across all OM scenarios
perf <- do.call(rbind, lapply(om_scenarios, function(om) {
  all_results <- om$results
  sim_list <- om$sim_list

  do.call(rbind, lapply(1:length(all_results), function(i) {
    sc <- all_results[[i]]$scenario
    ssb_mat <- sapply(all_results[[i]]$results, function(x) x$om$SSB)
    catch_mat <- sapply(all_results[[i]]$results, function(x) x$om$Catch)

    fb <- sim_list$feedback_start_yr
    proj_rows <- fb:nrow(ssb_mat)

    data.frame(
      om_scenario = om$name,
      scenario_id = sc$scenario_id,
      spr_x       = sc$spr_x,
      alpha       = sc$alpha,
      brp_type    = sc$brp_type,
      dyn_b0      = sc$dyn_b0,
      stability   = sc$stability,
      med_ssb     = median(ssb_mat[proj_rows, ]),
      p_crash     = mean(apply(ssb_mat[proj_rows, , drop=FALSE], 2, min) <
                           min(all_results[[1]]$results[[1]]$om$SSB[1:51]) * 0.1),
      ssb_cv      = median(apply(ssb_mat[proj_rows, , drop=FALSE], 2, sd) /
                             apply(ssb_mat[proj_rows, , drop=FALSE], 2, mean)),
      med_catch   = median(catch_mat[proj_rows, ]),
      catch_cv    = median(apply(catch_mat[proj_rows, , drop=FALSE], 2, sd) /
                             apply(catch_mat[proj_rows, , drop=FALSE], 2, mean)),
      catch_aav   = median(sapply(1:ncol(catch_mat), function(j) {
        cc <- catch_mat[proj_rows, j]
        mean(abs(diff(cc)) / head(cc, -1), na.rm = TRUE)
      }))
    )
  }))
}))

# NMDS on combined set
metric_cols <- c("med_ssb", "med_catch", "catch_aav")
perf_scaled <- scale(as.matrix(perf[, metric_cols]))
d <- dist(perf_scaled, method = "euclidean")

set.seed(42)
nmds <- metaMDS(d, k = 2, trymax = 200, autotransform = FALSE)
cat("Stress:", nmds$stress, "\n")

nmds_df <- data.frame(
  NMDS1 = nmds$points[, 1],
  NMDS2 = nmds$points[, 2],
  perf
)

nmds_df$om_scenario <- factor(nmds_df$om_scenario,
                              levels = c("Baseline", "Regime", "Crash"))
nmds_df$brp_type  <- factor(nmds_df$brp_type, levels = c("none", "threshold"),
                            labels = c("Constant F", "Threshold HCR"))
nmds_df$dyn_b0    <- factor(nmds_df$dyn_b0, levels = c(FALSE, TRUE),
                            labels = c("Static B0", "Dynamic B0"))
nmds_df$stability <- factor(nmds_df$stability,
                            levels = c("none", "symmetric", "asymmetric"),
                            labels = c("No Constraint", "Symmetric", "Asymmetric"))

# Envfit vectors
ef <- envfit(nmds, perf_scaled, permutations = 999)
vec_df <- as.data.frame(ef$vectors$arrows * sqrt(ef$vectors$r))
vec_df$label <- rownames(vec_df)
arrow_scale <- 0.6 * max(abs(nmds_df[, 1:2]))
vec_df$NMDS1 <- vec_df$NMDS1 * arrow_scale
vec_df$NMDS2 <- vec_df$NMDS2 * arrow_scale

# Outlier labels
dists <- sqrt(nmds_df$NMDS1^2 + nmds_df$NMDS2^2)
nmds_df$label <- ifelse(dists > quantile(dists, 0.90),
                        nmds_df$scenario_id, NA)

ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = brp_type, shape = dyn_b0, size = spr_x),
             alpha = 0.8) +
  geom_text(aes(label = label), size = 2.5, vjust = -1, colour = "grey30",
            na.rm = TRUE) +
  geom_segment(data = vec_df,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "grey40", linewidth = 0.5,
               inherit.aes = FALSE) +
  geom_text(data = vec_df,
            aes(x = NMDS1 * 1.12, y = NMDS2 * 1.12, label = label),
            size = 3, colour = "grey30", fontface = "italic",
            inherit.aes = FALSE) +
  ggh4x::facet_grid2(om_scenario ~ stability, scales = 'free', independent = 'y') +
  scale_colour_manual(values = c("Constant F" = "steelblue",
                                 "Threshold HCR" = "firebrick")) +
  scale_shape_manual(values = c("Static B0" = 16, "Dynamic B0" = 17)) +
  scale_size_continuous(range = c(1.5, 5), breaks = c(0.3, 0.5, 0.7)) +
  labs(colour = "BRP Type", shape = "Reference Point",
       size = expression(SPR[x])) +
  theme_bw(base_size = 14)
  # coord_equal()
