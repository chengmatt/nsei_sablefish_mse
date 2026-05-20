# Purpose: Utility functions for NSEI Sablefish MSE
# Creator: Matthew LH. Cheng
# Date 5/20/26

catch_to_f <- function(catch, sim_env) {

  # get F
  tmp_f <- catch_to_F_singlefleet(
    f_guess = 0.05, # guess for fishing mortality rate
    catch = catch[,1, 1,, 1], # catch values to use
    NAA = sim_env$NAA[1, 1, y+1,,, , sim], # numbers at age in simulation (truth)
    WAA = sim_env$WAA_fish[1, 1, y+1,, , , 1, sim], # weight-at-age in simulation (truth)
    natmort  = sim_env$natmort[1, 1, y+1, , , sim], # natural mortality in simulation (truth)
    fish_sel = sim_env$fish_sel[1,1, y+1,1, , , 1, sim], # fishery selectivity in simulation (truth)
    dmr = sim_env$dmr[1,y+1,1,1,sim], # dmr rate in simulation (truth)
    ret_sel = sim_env$ret_sel[1,1, y+1,1, , , 1, sim], # retention selectivity in simulation (truth)
    n.iter = 30
  )

  # input F into simulation environment
  sim_env$Fmort[,y+1,,,sim] <- array(tmp_f, dim = c(sim_env$n_regions, sim_env$n_seas, sim_env$n_fish_fleets)) # assign bisection values back into simulation
}

get_proj_catch <- function(obj, asmt_data, proj_opt, reference_points, sim_env) {

  # Get inputs for projection
  tmp_terminal_NAA <- array(obj$rep$NAA[,,y,,,], dim = c(asmt_data$n_pop, asmt_data$n_regions, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes)) # terminal numbers at age
  tmp_terminal_NAA0 <- array(obj$rep$NAA0[,,y,,,], dim = c(asmt_data$n_pop, asmt_data$n_regions, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes)) # terminal unfished numbers at age
  tmp_WAA <- array(rep(asmt_data$WAA[,,y,,,], each = proj_opt$n_proj_yrs), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes)) # weight at age
  tmp_WAA_fish <- array(rep(asmt_data$WAA_fish[,,y,,,,], each = proj_opt$n_proj_yrs), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes, asmt_data$n_fish_fleets)) # weight at age fishery
  tmp_MatAA <- array(rep(asmt_data$MatAA[,,y,,,], each = proj_opt$n_proj_yrs), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes)) # maturity at age
  tmp_fish_sel <- array(rep(obj$rep$fish_sel[,,y,,,,], each = proj_opt$n_proj_yrs), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes, asmt_data$n_fish_fleets)) # total selectivity
  tmp_ret_sel <- array(rep(obj$rep$ret_sel[,,y,,,,], each = proj_opt$n_proj_yrs), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes, asmt_data$n_fish_fleets)) # retained selectivity
  tmp_terminal_F <- array(obj$rep$Fmort[,y,,], dim = c(asmt_data$n_regions, asmt_data$n_seas, asmt_data$n_fish_fleets)) # terminal fishing mortality
  tmp_terminal_dmr <- array(obj$rep$dmr[,y,,], dim = c(asmt_data$n_regions, asmt_data$n_seas, asmt_data$n_fish_fleets)) # terminal discard mortality rate
  tmp_natmort <- array(rep(obj$rep$natmort[,,y,,], each = proj_opt$n_proj_yrs), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, length(asmt_data$ages), asmt_data$n_sexes)) # natural mortality
  tmp_recruitment <- array(obj$rep$Rec[,,1:y], dim = c(asmt_data$n_pop, asmt_data$n_regions, length(1:y))) # recruitment to use for projections
  tmp_sexratio <- array(replicate(n = proj_opt$n_proj_yrs, obj$rep$sexratio[,,y,]), dim = c(asmt_data$n_pop, asmt_data$n_regions, proj_opt$n_proj_yrs, asmt_data$n_sexes)) # recruitment sex ratio
  tmp_Movement <- array(dim = c(asmt_data$n_pop, asmt_data$n_regions, asmt_data$n_regions, proj_opt$n_proj_yrs,asmt_data$n_seas, length(asmt_data$ages), asmt_data$n_sexes))
  for(proj_yr in 1:proj_opt$n_proj_yrs) tmp_Movement[,,,proj_yr,,,] <- obj$rep$Movement[,,,y,,,] # Movement projections

  # Do projection to get TAC
  proj <- Do_Population_Projection(
    n_proj_yrs = proj_opt$n_proj_yrs,
    n_regions = sim_env$n_regions,
    n_ages = sim_env$n_ages,
    n_sexes = sim_env$n_sexes,
    n_pop = sim_env$n_pop,
    sexratio = tmp_sexratio,
    n_fish_fleets = sim_env$n_fish_fleets,
    do_recruits_move = sim_env$do_recruits_move,
    recruitment = tmp_recruitment,
    terminal_NAA = tmp_terminal_NAA,
    terminal_NAA0 = tmp_terminal_NAA0,
    terminal_F = tmp_terminal_F,
    dmr = tmp_terminal_dmr,
    natmort = tmp_natmort,
    WAA = tmp_WAA,
    WAA_fish = tmp_WAA_fish,
    MatAA = tmp_MatAA,
    fish_sel = tmp_fish_sel,
    ret_sel = tmp_ret_sel,
    Movement = tmp_Movement,
    f_ref_pt = reference_points$f_ref_pt,
    b_ref_pt = reference_points$b_ref_pt,
    HCR_function = proj_opt$HCR_function,
    recruitment_opt = proj_opt$recruitment_opt,
    fmort_opt = proj_opt$fmort_opt,
    t_spawn = sim_env$t_spawn,
    bh_rec_opt = proj_opt$bh_rec_opt
  )

  # Get Prescribed Catch
  catch <- proj$proj_Catch[,,-1,,,drop = FALSE] # get catch advice from projected year

  return(catch)
}
