# Purpose: Function to setup simplified estimation model for NSEI sablefish using SPoRC
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 5/20/26
setup_em <- function(sim_env, y, sim) {

  # Extract simulation data for current year and replicate
  sim_data <- simulation_data_to_SPoRC(sim_env, y, sim)

  # Model dimensions
  input_list <- Setup_Mod_Dim(years = 1:y, # vector of years
                              ages = 1:sim_env$n_ages, # vector of ages
                              lens = 1:sim_env$n_lens, # number of lengths
                              n_regions = sim_env$n_regions, # number of regions
                              n_sexes = sim_env$n_sexes, # number of sexes
                              n_fish_fleets = sim_env$n_fish_fleets, # number of fishery fleet
                              n_srv_fleets = sim_env$n_srv_fleets, # number of survey fleets
                              n_pop = sim_env$n_pop, # number of populations
                              verbose = F
  )

  # Recruitment setup
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    ln_sigmaR = array(log(c(1.2,1.2)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)), # 2 values for early and late sigma
    # Starting values for early and late sigmaR
    rec_model = "mean_rec",
    sigmaR_spec = "fix", # fix early sigmaR and late sigmaR
    init_age_strc = 1, # geometric series to derive initial age structure
    ln_global_R0 = log(13), # starting value for mean_rec
    t_spawn = sim_env$t_spawn # spawn timing
  )

  # Biological setup
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    # Data inputs
    WAA = sim_data$WAA,
    WAA_fish = sim_data$WAA_fish,
    MatAA = sim_data$MatAA,
    # Model options
    fit_lengths = 0,
    SizeAgeTrans = sim_data$SizeAgeTrans,
    AgeingError = sim_data$AgeingError,
    M_spec = "fix",     # fixing natural mortality
    Fixed_natmort = array(0.1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                       length(input_list$data$ages), input_list$data$n_sexes))
  )

  # Movement and tagging
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  # Fishery catch & fishing mortality
  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    # Data inputs
    ObsCatch = sim_data$ObsCatch,
    UseCatch = sim_data$UseCatch,
    # Model options
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    dmr_mean_spec = 'fix',
    # Fixing sigma C and F
    ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets)),
    logit_dmr_mean = array(qlogis(0.16), dim = c(input_list$data$n_regions,
                                                 input_list$data$n_seas, input_list$data$n_fish_fleets)) # fixing dmr at 0.16
  )

  # Survey selectivity and catchability
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    # Data inputs
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = sim_data$UseFishIdx,
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ISS_FishLenComps = sim_data$ISS_FishLenComps,

    # Model options
    fish_idx_type = c("biom"),
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("Multinomial"),
    FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1")
  )

  # Survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    # Data inputs
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    # Model options
    srv_idx_type = c("abd"),
    SrvAgeComps_LikeType = c("Multinomial"),
    SrvLenComps_LikeType = c("Multinomial"),
    SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1")
  )


  # Fishery selectivity and catchability
  ret_sel_array <- array(NA, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                     length(input_list$data$years), input_list$data$n_seas,
                                     length(input_list$data$ages), input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets))
  ret_sel_array[] <- sim_env$ret_sel[,,1:y,,,,,sim]

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    # Model options
    Use_fish_selex_prior = 1,
    fish_selex_prior = expand.grid(region = 1, par = 1:2, block = 1:3, sex = 1:2, fleet = 1, mu = 2, sd = 2),
    fish_sel_model = c("logist1_Fleet_1"), # fishery selex model
    fish_fixed_sel_pars_spec = c("est_all"), # whether to estiamte all fixed effects for fishery selectivity
    fish_sel_blocks = c(
      'Block_1_Year_1-20_Fleet_1',
      'Block_2_Year_21-47_Fleet_1',
      'Block_3_Year_48-terminal_Fleet_1'
    ),
    fish_q_spec = c("est_all"), # whether to estiamte all fixed effects for fishery catchability
    fish_q_blocks = c(
      'Block_1_Year_1-20_Fleet_1',
      'Block_2_Year_21-47_Fleet_1',
      'Block_3_Year_48-terminal_Fleet_1'
    ),
    use_fixed_ret_sel = 1,
    ret_sel_input = ret_sel_array # fixed retention curve
  )

  # Survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,

    # Model options
    Use_srv_selex_prior = 1, # selex priors
    srv_selex_prior = expand.grid(region = 1, par = 1:2, block = 1, sex = 1:2, fleet = 1, mu = 2, sd = 2),
    srv_sel_model = c("logist1_Fleet_1"), # survey selectivity form
    srv_fixed_sel_pars_spec = c("est_all"), # whether to estimate all fixed effects for survey selectivity
    srv_q_spec = c("est_all")  # whether to estiamte all fixed effects for survey catchability
  )

  # Data weighting
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets))
  )

  # Starting values for selectivity parameters
  input_list$par$fish_fixed_sel_pars[] <- log(c(c(4.82, 0.49, 4.34, 1.76, 4.86, 2.48, 8.27, 0.49, 5.49, 0.90, 6.84, 0.48))) # starting values for fishery
  input_list$par$srv_fixed_sel_pars[] <- log(c(5.51, 1.349, 7.52, 0.551)) # starting values

  # Starting values for catchability
  input_list$par$ln_srv_q[] <- -18
  input_list$par$ln_fish_q[] <- -18

  # Data Stuff
  # Set 0 for historical periods
  input_list$data$UseFishIdx[,1:4,,] <- 0
  input_list$data$UseFishAgeComps[,1:26,,] <- 0
  input_list$data$UseSrvAgeComps[,1:21,,] <- input_list$data$UseSrvIdx[,1:21,,] <- 0

  # If fishery catches are 0, don't use fishery age comps and indices
  zero_catch <- apply(sim_data$ObsCatch, 2, function(x) all(x == 0))  # logical vector over years
  if(any(zero_catch)) {
    input_list$data$UseFishAgeComps[, zero_catch, , ] <- 0
    input_list$data$UseFishIdx[, zero_catch, , ] <- 0
  }

  return(input_list)
}
