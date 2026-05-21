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

# Read in MSE OM
sim_list <- readRDS(here("outputs", 'base_mse_om.RDS'))

# Reference Point Options -------------------------------------------------
reference_points_opt <- list(
  n_avg_yrs = 1, # number of years to average over demographic rates
  SPR_x = 0.5, # SPR rates
  calc_rec_st_yr = 1, # year to use recruitment for bx% calculations
  rec_age = 0, # recruitment age
  type = 'single_region', # single region reference points
  what = "SPR" # using spr calculations
)

# Projection Options ------------------------------------------------------
proj_opt <- list(
  n_proj_yrs = 2,
  HCR_function = function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp # define stock status
    # If stock status is > 1
    if(stock_status >= 1) f <- frp
    # If stock status is between brp and alpha
    if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    # If stock status is less than alpha
    if(stock_status < alpha) f <- 0
    return(f)
  },
  recruitment_opt = 'mean_rec',
  fmort_opt = 'HCR',
  bh_rec_opt = NULL
)


# Helper Functions -----------------------------------------------------------------
run_single_sim <- function(sim, sim_list, reference_points_opt, proj_opt) {

  # Ensure worker has all dependencies
  library(SPoRC)
  source(here::here("R", "functions", "utils.R"))
  source(here::here("R", "functions", "setup_em.R"))

  # Each worker creates its own sim_env from the serializable sim_list
  sim_env <- Setup_sim_env(sim_list)

  # setup dimensions for storage
  n_yrs <- sim_env$n_yrs
  asmt_yrs <- sim_env$feedback_start_yr:n_yrs
  n_asmt_yrs <- length(asmt_yrs)

  # Storage for this sim
  ssb_em_list <- vector("list", n_asmt_yrs)
  conv_vec    <- rep(NA, n_asmt_yrs)

  # Run MSE
  for (y in 1:sim_env$n_yrs) {

    # Run Annual Cycle (anything prior to feedback = conditioning period)
    run_annual_cycle(y, sim, sim_env)

    # Start Assessment / Projection Period
    if (y >= sim_env$feedback_start_yr) {

      # Index for when assessment starts
      asmt_idx <- y - sim_env$feedback_start_yr + 1

      # Get Assessment Results
      local_y <- y
      local_sim <- sim
      asmt_result <- tryCatch({

        # Setup EM
        tmp_list  <- setup_em(sim_env, local_y, local_sim)
        asmt_data <- tmp_list$data

        # Fit EM
        obj <- fit_model(
          asmt_data,
          tmp_list$par,
          tmp_list$map,
          NULL,
          newton_loops = 0,
          silent = TRUE
        )

        # Get sdreport
        sd_rep <- sdreport(obj)

        # Derive Reference Points from Model
        reference_points <- get_closed_loop_reference_points(
          use_true_values = FALSE,
          sim_env = sim_env,
          asmt_data = asmt_data,
          asmt_rep = obj$rep,
          y = local_y,
          sim = local_sim,
          reference_points_opt = reference_points_opt,
          n_proj_yrs = proj_opt$n_proj_yrs
        )

        # Run Projection to Get Catch Advice
        tmp_catch <- get_proj_catch(obj, asmt_data, proj_opt, reference_points, sim_env, local_y, local_sim)

        # Convert Prescribed Catch to True Fishing Mortality
        if (local_y < sim_env$n_yrs) catch_to_f(tmp_catch, sim_env, local_y, local_sim)

        list(
          success = TRUE,
          ssb_em = obj$rep$SSB[1:local_y],
          conv = post_optim_sanity_checks(sd_rep, obj$rep, gradient_tol = 0.01)
        )

      }, error = function(e) {
        message(sprintf("Sim %d, year %d failed: %s", local_sim, local_y, e$message))
        list(success = FALSE, ssb_em = NULL, conv = NA)
      })

      # Store EM Results
      ssb_em_list[[asmt_idx]] <- asmt_result$ssb_em
      conv_vec[asmt_idx]      <- asmt_result$conv

    } # end feedback
  } # end y

  # Return lists
  list(
    om = list(
      SSB   = sim_env$SSB[,,,sim],
      Catch = sim_env$TrueCatch[,,,,sim],
      NAA   = sim_env$NAA[,,,,,,sim],
      Fmort = sim_env$Fmort[,,,,sim]
    ),
    ssb_em      = ssb_em_list,
    convergence = conv_vec,
    failed      = FALSE,
    failure_year = NA
  )

}


# Run in Parrallel --------------------------------------------------------
options(future.globals.maxSize = 1.5 * 1024^3)
plan(multisession, workers = availableCores() - 5)

results <- future_map(
  1:sim_list$n_sims,
  ~ run_single_sim(.x, sim_list, reference_points_opt, proj_opt),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

# save results
saveRDS(results, here("scratch", "om_em_bias_test.RDS"))

# Process Results ---------------------------------------------------------
results <- readRDS(here("scratch", "om_em_bias_test.RDS"))
par(mfrow = c(3, 2), mar = c(4, 5, 3, 1))

n_sims <- length(results)
n_yrs <- sim_list$n_yrs
feedback_start <- sim_list$feedback_start_yr
asmt_yrs <- feedback_start:n_yrs
n_asmt <- length(asmt_yrs)
yr_seq <- 1975:(1975 + n_yrs - 1)

om_mat <- sapply(results, function(x) x$om$SSB)
om_med <- apply(om_mat, 1, median)

for (j in seq_along(asmt_yrs)) {
  em_mat_j <- sapply(results, function(x) x$ssb_em[[j]])
  em_med <- apply(em_mat_j, 1, median)
  em_lo99 <- apply(em_mat_j, 1, quantile, 0.005)
  em_hi99 <- apply(em_mat_j, 1, quantile, 0.995)
  em_lo95 <- apply(em_mat_j, 1, quantile, 0.025)
  em_hi95 <- apply(em_mat_j, 1, quantile, 0.975)
  em_lo90 <- apply(em_mat_j, 1, quantile, 0.05)
  em_hi90 <- apply(em_mat_j, 1, quantile, 0.95)
  yrs_j   <- yr_seq[1:asmt_yrs[j]]

  plot(yr_seq, om_med, type = "n",
       ylim = range(c(em_lo99, em_hi99, om_med)),
       xlab = "Year", ylab = "SSB",
       main = paste("Asmt Year", yrs_j[length(yrs_j)]),
       cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.8)
  polygon(c(yrs_j, rev(yrs_j)), c(em_lo99, rev(em_hi99)),
          col = adjustcolor("steelblue", 0.15), border = NA)
  polygon(c(yrs_j, rev(yrs_j)), c(em_lo95, rev(em_hi95)),
          col = adjustcolor("steelblue", 0.25), border = NA)
  polygon(c(yrs_j, rev(yrs_j)), c(em_lo90, rev(em_hi90)),
          col = adjustcolor("steelblue", 0.35), border = NA)
  lines(yrs_j, em_med, col = "steelblue", lwd = 3)
  lines(yr_seq, om_med, col = "black", lwd = 3, lty = 2)
  abline(v = yr_seq[feedback_start], lty = 3, col = "red", lwd = 2)
  legend("topright", legend = c("OM", "EM median", "90% CI", "95% CI", "99% CI", "Feedback"),
         col = c("black", "steelblue",
                 adjustcolor("steelblue", 0.35),
                 adjustcolor("steelblue", 0.25),
                 adjustcolor("steelblue", 0.15), "red"),
         lwd = c(3, 3, 8, 8, 8, 2), lty = c(2, 1, 1, 1, 1, 3),
         cex = 1.1, bty = "n")
}

par(mfrow = c(1, 1))
