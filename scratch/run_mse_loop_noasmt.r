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
spr_x_vals    <- seq(0.3, 0.7, by = 0.05)
alpha_vals    <- seq(0, 0.5, 0.05)
dyn_b0_vals   <- FALSE
stability_opts <- c("none", "symmetric", "asymmetric")

# Stability Constraints
stab_pars <- list(
  none       = list(max_increase = Inf, max_decrease = Inf),
  symmetric  = list(max_increase = 0.15, max_decrease = 0.15),
  asymmetric = list(max_increase = 0.15, max_decrease = Inf) # fast down, slow up
)

# Threshold rule -- full fine grid
scenario_grid <- expand.grid(
  spr_x     = spr_x_vals,
  alpha     = alpha_vals,
  dyn_b0    = dyn_b0_vals,
  stability = stability_opts,
  stringsAsFactors = FALSE
)
scenario_grid$hcr_type <- "threshold"

# Constant F (no ramp): F50 w/ asymmetric stability only (current management analogue)
scenario_grid <- rbind(
  scenario_grid,
  data.frame(spr_x = 0.5, alpha = NA, dyn_b0 = FALSE,
             stability = "asymmetric", hcr_type = "constant",
             stringsAsFactors = FALSE)
)

# Add scenario ID
scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))

run_single_sim <- function(sim, sim_list, scenario, use_true_values) {

  library(SPoRC)
  source(here::here("R", "functions", "utils.R"))
  source(here::here("R", "functions", "setup_em.R"))

  # Unpack scenario
  spr_x       <- scenario$spr_x
  alpha       <- scenario$alpha
  use_dyn_b0  <- scenario$dyn_b0
  stab_type   <- scenario$stability

  # Stability parameters
  stab <- stab_pars[[stab_type]]

  if (identical(scenario$hcr_type, "constant")) {
    # Constant F HCR: applies frp regardless of stock status (no ramp, no closure)
    hcr_fn <- function(x, frp, brp, local_alpha = NA) frp
  } else {
    # Threshold (hockey-stick) HCR: ramps F linearly to 0 as stock_status falls
    # from 1 down to alpha, closes the fishery below alpha
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
      # closure years (F = 0) need patched F values: the ref-point optimizer needs a
      # nonzero F pattern (input 1), and the advice projection needs a nonzero F for
      # its fleet split but should treat the closed current year as ~unfished
      # (input 1e-4; input 1 would crash the projected stock and delay reopening).
      # both patches go on isolated copies -- sim_env is an environment, so plain
      # assignment aliases it and the patch would end up in the stored OM outputs
      tmp_sim_env <- sim_env
      proj_sim_env <- sim_env
      if (sim_env$Fmort[,y,,,sim] == 0) {
        tmp_sim_env <- list2env(as.list(sim_env))
        tmp_sim_env$Fmort[,y,,,sim] <- 1
        proj_sim_env <- list2env(as.list(sim_env))
        proj_sim_env$Fmort[,y,,,sim] <- 1e-4
      }
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
        proj_sim_env, local_y, local_sim,
        use_true_values,
        # dyn_bx follows spr_x when dynamic B0 is used
        use_dyn_b0 = use_dyn_b0,
        dyn_bx     = if (use_dyn_b0) spr_x else spr_x
      )

      # Apply Stability Constraint (only when catch from the base HCR is not 0;
      # HCR closures bypass the constraint and shut the fishery immediately)
      if(tmp_catch != 0) {
        constrained_catch <- apply_stability(
          new_catch    = as.vector(tmp_catch),
          prev_catch   = sim_env$TrueCatch[1,y,1,1,sim], # use previous years catch
          max_increase = stab$max_increase,
          max_decrease = stab$max_decrease
        )
      } else {
        constrained_catch <- 0 # closure: do not carry over last year's constrained catch
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
    "\n Scenario %d/%d: SPR=%.1f, alpha=%.2f, dynB0=%s, stab=%s, hcr=%s\n",
    sc, nrow(scenario_grid),
    scenario$spr_x, scenario$alpha,
    scenario$dyn_b0, scenario$stability, scenario$hcr_type
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
    "\n===== Scenario %d/%d: SPR=%.1f, alpha=%.2f, dynB0=%s, stab=%s, hcr=%s =====\n",
    sc, nrow(scenario_grid),
    scenario$spr_x, scenario$alpha,
    scenario$dyn_b0, scenario$stability, scenario$hcr_type
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
    "\n===== Scenario %d/%d: SPR=%.1f, alpha=%.2f, dynB0=%s, stab=%s, hcr=%s =====\n",
    sc, nrow(scenario_grid),
    scenario$spr_x, scenario$alpha,
    scenario$dyn_b0, scenario$stability, scenario$hcr_type
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

