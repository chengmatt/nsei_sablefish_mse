# Purpose: To condition a basic OM for use in MSE
# Creator: Matthew LH. Cheng
# Date 5/18/26

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
source(here("R", 'functions', 'utils.R'))
source(here("R", 'functions', 'setup_em.R'))

# Read in model
out <- readRDS(here("inputs",  "v23_3f_3s_2016_NEW_FINAL.RDS"))
rep <- readRDS(here("inputs",  "actual_report.rds"))
data <- readRDS(here("inputs", 'sablefish_data_19May2026.RDS'))

### OM Dimensions -----------------------------------------------------------
set.seed(123)
n_sims <- 1
n_yrs <- length(rep$Fmort)
closed_loop_yrs <- 1
n_ages <- 30
n_lens <- 30
n_sexes <- 2
n_fish_fleets <- 1
n_srv_fleets <- 1
n_regions <- 1
n_seas <- 1
n_pop <- 1
feedback_start_yr <- n_yrs

### Setup Model Dimensions --------------------------------------------------
sim_list <- Setup_Sim_Dim(n_sims = n_sims,
                          n_yrs = n_yrs + closed_loop_yrs,
                          n_regions = n_regions,
                          n_ages = n_ages,
                          n_lens = n_lens,
                          n_sexes = n_sexes,
                          n_fish_fleets = n_fish_fleets,
                          n_srv_fleets = n_srv_fleets,
                          feedback_start_yr = n_yrs,
                          run_feedback = TRUE
)

# Setup Simulation Containers ---------------------------------------------
sim_list <- Setup_Sim_Containers(sim_list) # set up simulation containers to use

# Setup Fishery Processes -------------------------------------------------

# Fill in fishing mortality
Fmort <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
Fmort[,1:n_yrs,,] <- rep$Fmort

# Fill in fishery selectivity
fish_sel_input <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
fish_sel_input[1,1,1:n_yrs,1,,1,1] <- rep$fsh_slx[,,2] # fill in conditioning period females
fish_sel_input[1,1,1:n_yrs,1,,2,1] <- rep$fsh_slx[,,1] # fill in conditioning period males
for(y in (n_yrs + 1):dim(fish_sel_input)[3]) {
  fish_sel_input[1, 1, y, 1, , 1, 1] <- rep$fsh_slx[n_yrs, , 2]  # projection period females
  fish_sel_input[1, 1, y, 1, , 2, 1] <- rep$fsh_slx[n_yrs, , 1]  # projection period males
}
fish_sel_input <- replicate(n = sim_list$n_sims, fish_sel_input)

# Fill in retention selectivity
retention_m <- round(((rep$Z[1,,1] - 0.1) / rep$F[1,,1] - 0.16) / (1 - 0.16), 10) # derive retention from model object b/c data file retention seems incorrect
retention_f <- round(((rep$Z[1,,2] - 0.1) / rep$F[1,,2] - 0.16) / (1 - 0.16), 10) # derive retention from model object b/c data file retention seems incorrect
ret_sel_input <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
for(i in 1:(n_yrs + closed_loop_yrs)) {
  ret_sel_input[1,1,i,1,,1,1] <- retention_f # fill in conditioning period (females)
  ret_sel_input[1,1,i,1,,2,1] <- retention_m # fill in conditioning period (males)
}
ret_sel_input <- replicate(n = sim_list$n_sims, ret_sel_input)

# Fill in fishery q
fish_q_input = array(NA, dim = c(n_regions, n_yrs + closed_loop_yrs, n_fish_fleets, n_sims))
fish_q_input[,1:20,,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'fsh_logq'][1])
fish_q_input[,21:47,,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'fsh_logq'][2])
fish_q_input[,-c(1:47),,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'fsh_logq'][3])

# Fill in fishery observation error (using the mean)
ObsFishIdx_SE <- array(mean(data$sigma_fsh_cpue),  dim = c(n_regions, n_yrs + closed_loop_yrs, n_seas, n_fish_fleets))

# Fill in ISS for fishery ages (using the mean)
ISS_FishAgeComps <- array(round(mean(data$effn_fsh_age)), dim = c(n_regions, n_yrs  + closed_loop_yrs, n_seas, n_sexes, n_fish_fleets, n_sims))
ISS_FishLenComps <- array(round(mean(data$effn_fsh_len)) * 0.5, dim = c(n_regions, n_yrs  + closed_loop_yrs, n_seas, n_sexes, n_fish_fleets, n_sims))

# setup fishery simulation processes
sim_list <- Setup_Sim_Fishing(
  sim_list = sim_list, # update simulate list
  ln_sigmaC = array(log(0.05), dim = c(n_regions, n_yrs + closed_loop_yrs, n_seas, n_fish_fleets)),
  Fmort_input = SPoRC:::extend_years(replicate(n = n_sims, Fmort), n_years = closed_loop_yrs, 2, fill = 'zeros'),
  dmr_input = array(0.16, dim = c(n_regions, n_yrs + closed_loop_yrs, n_seas, n_fish_fleets, n_sims)),
  fish_sel_input = fish_sel_input,
  ret_sel_input = ret_sel_input,
  fish_q_input = fish_q_input,
  ObsFishIdx_SE = ObsFishIdx_SE,
  ISS_FishAgeComps = ISS_FishAgeComps,
  ISS_FishLenComps = ISS_FishLenComps,
  fish_idx_type = 'biom',
  catch_units = 'biom',
  discard_units = 'biom'
)

# Setup Survey Processes --------------------------------------------------

# Fill in survey selectivity
srv_sel_input <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))

# Conditioning based on selectivity estimates in the middle period (normal looking and estimated in assessment)
for(y in 1:sim_list$n_yrs) {
  srv_sel_input[1, 1, y, 1, , 1, 1] <- rep$srv_slx[30, , 2]  # projection period females
  srv_sel_input[1, 1, y, 1, , 2, 1] <- rep$srv_slx[30, , 1]  # projection period males
}
srv_sel_input <- replicate(n = sim_list$n_sims, srv_sel_input)

# Fill in survey q
srv_q_input = array(NA, dim = c(n_regions, n_yrs + closed_loop_yrs, n_srv_fleets, n_sims))
srv_q_input[,1:25,,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'srv_logq'][1])
srv_q_input[,26:42,,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'srv_logq'][2])
srv_q_input[,-c(1:42),,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'srv_logq'][3])

# Adjust
srv_q_input[] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'srv_logq'][2])

# Fill in survey observation error (using the mean)
ObsSrvIdx_SE <- array(mean(data$sigma_srv_cpue),  dim = c(n_regions, n_yrs  + closed_loop_yrs, n_seas, n_srv_fleets))

# Fill in ISS for survey ages (using the mean)
ISS_SrvAgeComps <- array(round(mean(data$effn_srv_age)), dim = c(n_regions, n_yrs  + closed_loop_yrs, n_seas, n_sexes, n_srv_fleets, n_sims))
ISS_SrvLenComps <- array(round(mean(data$effn_srv_len)) * 0.5, dim = c(n_regions, n_yrs  + closed_loop_yrs, n_seas, n_sexes, n_srv_fleets, n_sims))

sim_list <- Setup_Sim_Survey(
  sim_list = sim_list,
  srv_sel_input = srv_sel_input,
  srv_q_input = srv_q_input,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  srv_idx_type = 'abd',
  comp_srvage_like = "Multinomial",
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  ISS_SrvLenComps = ISS_SrvLenComps,
  t_srv = array(data$srv_month, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets))
)

# Setup Biologicals -------------------------------------------------------
sim_list <- Setup_Sim_Biologicals(
  sim_list = sim_list,
  natmort_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_ages, n_sexes, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) {
      tmp[1,1,i,,1,] <- 0.1  # females
      tmp[1,1,i,,2,] <- 0.1  # males
    }
    tmp
  },
  WAA_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) {
      tmp[1,1,i,1,,1,] <- replicate(n_sims, data$data_srv_waa[1,,2]) # females
      tmp[1,1,i,1,,2,] <- replicate(n_sims, data$data_srv_waa[1,,1]) # males
    }
    tmp
  },
  WAA_fish_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_fish_fleets, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) {
      tmp[1,1,i,1,,1,1,] <- replicate(n_sims, data$data_fsh_waa[1,,2]) # females
      tmp[1,1,i,1,,2,1,] <- replicate(n_sims, data$data_fsh_waa[1,,1]) # males
    }
    tmp
  },
  WAA_srv_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_srv_fleets, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) {
      tmp[1,1,i,1,,1,1,] <- replicate(n_sims, data$data_srv_waa[1,,2]) # females
      tmp[1,1,i,1,,2,1,] <- replicate(n_sims, data$data_srv_waa[1,,1]) # males
    }
    tmp
  },
  MatAA_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) {
      tmp[1,1,i,1,,1,] <- replicate(n_sims, data$prop_mature[1,]) # females
      tmp[1,1,i,1,,2,] <- 0  # males
    }
    tmp
  },
  AgeingError_input = {
    tmp <- array(NA, dim = c(n_yrs + closed_loop_yrs, n_ages, n_ages, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) tmp[i,,,] <- replicate(n_sims, (data$ageing_error / rowSums(data$ageing_error)))
    tmp
  },
  SizeAgeTrans_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_seas, n_lens, n_ages, n_sexes, n_sims))
    for(i in 1:(n_yrs + closed_loop_yrs)) {
      tmp[1,1,i,1,,,1,] <- replicate(n_sims, t(data$agelen_key_fsh[,,2])) # females
      tmp[1,1,i,1,,,2,] <- replicate(n_sims, t(data$agelen_key_fsh[,,1])) # males
    }
    tmp
  }
)

# Setup Recruitment -------------------------------------------------------
sim_list <- Setup_Sim_Rec(
  sim_list = sim_list,
  init_age_strc = "scalar_no_move",
  R0_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_sims))
    tmp[,,-1,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'log_rbar']) # R0
    tmp[,,1,] <- exp(out$rep$par.fixed[names(out$rep$par.fixed) == 'log_rinit']) # initial R0
    tmp
  },
  sexratio_input = array(0.5, dim = c(n_pop, n_regions, n_yrs + closed_loop_yrs, n_sexes, n_sims)),
  ln_sigmaR = array(log(0), dim = c(2, n_pop, n_regions)),
  Rec_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_yrs, n_sims))
    for(i in 1:n_sims) tmp[1,1,,i] <- rep$pred_rec
    tmp
  },
  ln_InitDevs_input = {
    tmp <- array(NA, dim = c(n_pop, n_regions, n_ages - 1, n_sims))
    for(i in 1:n_sims) tmp[1,1,,i] <- c(out$rep$par.fixed[names(out$rep$par.fixed) == 'log_rinit_devs'], 0) # input 0 for plus group
    tmp
  },
  recruitment_opt = 'resample_from_input',
  t_spawn = data$spawn_month
)

# Setup Tagging -----------------------------------------------------------
sim_list <- Setup_Sim_Tagging(
  sim_list = sim_list,
  use_conv_fish_tagging = 0 # not used
)

# Setup Movement ----------------------------------------------------------
sim_list$Movement <- array(1, dim = c(n_pop, n_regions, n_regions, n_yrs + closed_loop_yrs, n_seas, n_ages, n_sexes, n_sims))
sim_list$sgl_seas_spawning_movement <- array(1, dim = c(n_pop, n_regions, n_regions, n_yrs + closed_loop_yrs, n_ages, n_sexes, n_sims))

# Save conitioned OM
saveRDS(sim_list, here("outputs", "base_mse_om.RDS"))
