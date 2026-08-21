# Purpose: Condense the raw MSE result files (assessment and no-assessment) into a
#          compact cache used by the plotting scripts. The no-assessment files are
#          ~7 GB each, so they are read one at a time and dropped before the next.
# Creator: Matthew LH. Cheng

library(here)
term_cond_yr <- 2026  # terminal conditioning (assessment) year
RETRO_LAG    <- 10    # peels retained per assessment for retrospective diagnostics

# HCRs shared by both run sets. The assessment run uses a reduced grid
# (spr_x in {0.4, 0.5}, alpha in {0, 0.25, 0.5}, asymmetric only), so this is the
# full set of rules that can be compared across modes. Reference rule first.
hcr_set <- list(
  "Constant F50" = list(spr_x = 0.5, alpha = NA,   stability = "asymmetric", hcr_type = "constant"),
  "F40, a=0"     = list(spr_x = 0.4, alpha = 0,    stability = "asymmetric", hcr_type = "threshold"),
  "F40, a=0.25"  = list(spr_x = 0.4, alpha = 0.25, stability = "asymmetric", hcr_type = "threshold"),
  "F40, a=0.5"   = list(spr_x = 0.4, alpha = 0.5,  stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0"     = list(spr_x = 0.5, alpha = 0,    stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.25"  = list(spr_x = 0.5, alpha = 0.25, stability = "asymmetric", hcr_type = "threshold"),
  "F50, a=0.5"   = list(spr_x = 0.5, alpha = 0.5,  stability = "asymmetric", hcr_type = "threshold")
)

# the grids are built with seq(), so alpha/spr_x carry floating-point noise --
# match on tolerance rather than identity
match_hcr <- function(scenario, h) {
  sc_type <- if (is.null(scenario$hcr_type)) "threshold" else scenario$hcr_type
  if (sc_type != h$hcr_type) return(FALSE)
  if (as.character(scenario$stability) != h$stability) return(FALSE)
  if (!isTRUE(all.equal(scenario$spr_x, h$spr_x))) return(FALSE)
  # constant F has no ramp, so alpha is not part of its identity
  if (h$hcr_type == "constant") return(TRUE)
  isTRUE(all.equal(scenario$alpha, h$alpha))
}

# replicate x year matrix helper: rows = year, cols = replicate
rep_mat <- function(res, f) sapply(res, function(x) as.vector(f(x)))

# per-replicate catch AAV over the projection period: sum|dC| / sum C
catch_aav <- function(catch_mat, proj_rows) {
  apply(catch_mat[proj_rows, , drop = FALSE], 2, function(cc) {
    if (sum(head(cc, -1)) == 0) return(NA_real_)
    sum(abs(diff(cc))) / sum(head(cc, -1))
  })
}

# Pull the compact per-HCR record for one scenario entry
summarize_hcr <- function(s, fb, mode) {
  res <- s$results
  out <- list(
    scenario_id = s$scenario$scenario_id,
    SSB   = rep_mat(res, function(x) x$om$SSB),
    Catch = rep_mat(res, function(x) x$om$Catch),
    Fmort = rep_mat(res, function(x) x$om$Fmort),
    Rec   = rep_mat(res, function(x) x$om$Rec),
    bx    = rep_mat(res, function(x) x$bx)
  )
  if (mode == "asmt") {
    # SSB the assessment reported for its own terminal year -- this is the
    # quantity the HCR actually saw, so it drives the assessment-error signal
    out$ssb_em_term <- rep_mat(res, function(x)
      sapply(x$ssb_em, function(v) if (is.null(v)) NA_real_ else v[length(v)]))
    # full SSB series from the final assessment, for retrospective comparison
    n_yrs <- nrow(out$SSB)
    out$ssb_em_final <- sapply(res, function(x) {
      v <- x$ssb_em[[length(x$ssb_em)]]
      if (is.null(v)) rep(NA_real_, n_yrs) else c(v, rep(NA_real_, n_yrs - length(v)))
    })
    out$pdHess <- rep_mat(res, function(x) x$pdHess)
    out$grad   <- rep_mat(res, function(x) x$grad)

    # Retrospective structure. Each closed-loop assessment is a peel relative to
    # every later one, so keep each assessment's estimate of its own terminal
    # year and of the RETRO_LAG years before it -- enough to form Mohn's rho at
    # any reference year without carrying the full estimated series.
    # ssb_em_tail[assessment, lag + 1, replicate] = the SSB that the assessment
    # run in year y reported for year y - lag.
    n_asmt <- length(res[[1]]$ssb_em)
    out$ssb_em_tail <- array(NA_real_, c(n_asmt, RETRO_LAG + 1, length(res)))
    for (j in seq_along(res)) {
      em <- res[[j]]$ssb_em
      for (i in seq_len(n_asmt)) {
        v <- em[[i]]
        if (is.null(v)) next
        idx <- length(v) - (0:RETRO_LAG)
        ok  <- idx >= 1
        out$ssb_em_tail[i, ok, j] <- v[idx[ok]]
      }
    }
    # full estimated series for one replicate, for the retrospective fan plot
    out$ssb_em_example <- res[[1]]$ssb_em
  }
  out
}

# Full-grid per-scenario summary (no-assessment runs only -- the assessment grid
# is too coarse for the heatmaps). Keeps yearly medians so the existing grid
# figures can be rebuilt exactly without re-reading the raw files.
summarize_grid <- function(s, fb) {
  sc    <- s$scenario
  ssb   <- rep_mat(s$results, function(x) x$om$SSB)
  catch <- rep_mat(s$results, function(x) x$om$Catch)
  proj  <- fb:nrow(ssb)
  ssb_threshold <- max(ssb[seq_len(fb - 1), 1]) * 0.1
  data.frame(
    scenario_id = sc$scenario_id, spr_x = sc$spr_x, alpha = sc$alpha,
    dyn_b0 = sc$dyn_b0, stability = sc$stability,
    hcr_type = if (is.null(sc$hcr_type)) "threshold" else sc$hcr_type,
    med_ssb   = median(apply(ssb,   1, median)[proj]),
    med_catch = median(apply(catch, 1, median)[proj]),
    catch_aav = median(catch_aav(catch, proj), na.rm = TRUE),
    p_crash   = mean(ssb[proj, ] < ssb_threshold),
    stringsAsFactors = FALSE
  )
}

# Extract one results file -------------------------------------------------

extract_file <- function(path, om_name, mode, fb, want_grid) {
  cat(sprintf("[%s / %s] reading %s ...\n", om_name, mode, basename(path)))
  t0 <- Sys.time()
  all_results <- readRDS(path)
  cat(sprintf("  loaded %d scenarios in %.1f min\n", length(all_results),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # per-HCR compact records
  hcrs <- lapply(names(hcr_set), function(lab) {
    h  <- hcr_set[[lab]]
    hit <- which(sapply(all_results, function(s) match_hcr(s$scenario, h)))
    if (!length(hit)) {
      warning(sprintf("HCR '%s' not found in %s", lab, basename(path)))
      return(NULL)
    }
    summarize_hcr(all_results[[hit[1]]], fb, mode)
  })
  names(hcrs) <- names(hcr_set)

  grid <- if (want_grid) {
    do.call(rbind, lapply(all_results, summarize_grid, fb = fb))
  } else NULL

  rm(all_results); gc(verbose = FALSE)
  cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(om = om_name, mode = mode, fb = fb, hcrs = hcrs, grid = grid)
}

# Run ----------------------------------------------------------------------

om_files <- list(
  Baseline = "base_mse_om.RDS",
  Regime   = "bh_regime_mse_om.RDS",
  Crash    = "bh_crash_mse_om.RDS"
)
res_files <- list(
  Baseline = c(asmt = "base_results_asmt.RDS",   noasmt = "base_results_noasmt.RDS"),
  Regime   = c(asmt = "regime_results_asmt.RDS", noasmt = "regime_results_noasmt.RDS"),
  Crash    = c(asmt = "crash_results_asmt.RDS",  noasmt = "crash_results_noasmt.RDS")
)

cache <- list()
for (om_name in names(res_files)) {
  fb <- readRDS(here("outputs", om_files[[om_name]]))$feedback_start_yr
  for (mode in c("asmt", "noasmt")) {
    key <- paste(om_name, mode, sep = "_")
    cache[[key]] <- extract_file(here("outputs", res_files[[om_name]][[mode]]),
                                 om_name, mode, fb, want_grid = (mode == "noasmt"))
  }
}

cache$meta <- list(term_cond_yr = term_cond_yr, hcr_set = hcr_set,
                   retro_lag = RETRO_LAG)
saveRDS(cache, here("outputs", "mse_summaries.RDS"))
cat("\nwrote outputs/mse_summaries.RDS\n")
