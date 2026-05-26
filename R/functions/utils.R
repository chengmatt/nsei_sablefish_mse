# Purpose: Utility functions for NSEI Sablefish MSE
# Creator: Matthew LH. Cheng
# Date 5/20/26

catch_to_f <- function(catch, sim_env, y, sim) {

  # get F
  tmp_f <- catch_to_F_singlefleet(
    f_guess = 0.05, # guess for fishing mortality rate
    catch = catch, # catch values to use
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

get_proj_catch <- function(obj, asmt_data, proj_opt, reference_points, sim_env, y, sim, use_true_values,
                           use_dyn_b0 = FALSE, dyn_bx = NULL) {

  if(use_true_values) { # use true simulation values
    # Get inputs for projection
    tmp_terminal_NAA <- array(sim_env$NAA[,,y,,,,sim], dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)) # terminal numbers at age
    tmp_terminal_NAA0 <- array(sim_env$NAA0[,,y,,,,sim], dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)) # terminal unfished numbers at age
    tmp_WAA <- array(rep(sim_env$WAA[,,y,,,,sim], each = proj_opt$n_proj_yrs), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)) # weight at age
    tmp_WAA_fish <- array(rep(sim_env$WAA_fish[,,y,,,,,sim], each = proj_opt$n_proj_yrs), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)) # weight at age fishery
    tmp_MatAA <- array(rep(sim_env$MatAA[,,y,,,,sim], each = proj_opt$n_proj_yrs), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)) # maturity at age
    tmp_fish_sel <- array(rep(sim_env$fish_sel[,,y,,,,,sim], each = proj_opt$n_proj_yrs), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)) # total selectivity
    tmp_ret_sel <- array(rep(sim_env$ret_sel[,,y,,,,,sim], each = proj_opt$n_proj_yrs), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)) # retained selectivity
    tmp_terminal_F <- array(sim_env$Fmort[,y,,,sim], dim = c(sim_env$n_regions, sim_env$n_seas, sim_env$n_fish_fleets)) # terminal fishing mortality
    tmp_terminal_dmr <- array(sim_env$dmr[,y,,,sim], dim = c(sim_env$n_regions, sim_env$n_seas, sim_env$n_fish_fleets)) # terminal discard mortality rate
    tmp_natmort <- array(rep(sim_env$natmort[,,y,,,sim], each = proj_opt$n_proj_yrs), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_ages, sim_env$n_sexes)) # natural mortality
    tmp_recruitment <- array(sim_env$Rec[,,1:y,sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y))) # recruitment to use for projections
    tmp_sexratio <- array(replicate(n = proj_opt$n_proj_yrs, sim_env$sexratio[,,y,,sim]), dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs, sim_env$n_sexes)) # recruitment sex ratio
    tmp_Movement <- array(dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, proj_opt$n_proj_yrs,sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
    for(proj_yr in 1:proj_opt$n_proj_yrs) tmp_Movement[,,,proj_yr,,,] <- sim_env$Movement[,,,y,,,,sim] # Movement projections
    if(use_dyn_b0) b_ref <- array(sim_env$Dynamic_SSB0[1,1,y,sim] * dyn_bx, dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs))
    else b_ref <- reference_points$b_ref_pt
  } else {
    # use assessment values
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
    if(use_dyn_b0) b_ref <- array(obj$rep$Dynamic_SSB0[1,1,y] * dyn_bx, dim = c(sim_env$n_pop, sim_env$n_regions, proj_opt$n_proj_yrs))
    else b_ref <- reference_points$b_ref_pt
  }

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
    b_ref_pt = b_ref,
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

#' Apply stability constraint to catch advice
#'
#' @param new_catch Proposed catch from projection
#' @param prev_catch Previous year's realized catch
#' @param max_increase Maximum proportional increase (e.g., 0.15 for 15%)
#' @param max_decrease Maximum proportional decrease (e.g., 0.15 for 15%)
#' @return Constrained catch
apply_stability <- function(new_catch, prev_catch, max_increase = Inf, max_decrease = Inf) {
  if (is.na(prev_catch) || prev_catch <= 0) return(new_catch)
  ratio <- new_catch / prev_catch
  # Cap increase
  if (ratio > (1 + max_increase)) {
    new_catch <- prev_catch * (1 + max_increase)
  }
  # Cap decrease
  if (ratio < (1 - max_decrease)) {
    new_catch <- prev_catch * (1 - max_decrease)
  }
  return(max(new_catch, 0))
}
