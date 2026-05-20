# Purpose: To explore model outputs from the NSEI sablefish model
# Creator: Matthew LH. Cheng
# Date: 5/18/26


# Setup -------------------------------------------------------------------

library(here)
library(TMB)

# Read in model
out <- readRDS(here("inputs", "v23_3f_3s_2016_NEW_FINAL.RDS"))
rep <- readRDS(here("inputs", "actual_report.rds"))

# Compile model
compile(here("src", "scaa_mod_v23_Aaron2.cpp"))
dyn.load(TMB::dynlib(here("src", "scaa_mod_v23_Aaron2")))

# Visualize ---------------------------------------------------------------

# spawning biomass
plot(rowSums(rep$spawn_biom), type = 'l')

# fmort
plot(rep$Fmort, type = 'l')

# landed catch
plot(rep$pred_landed, type = 'l')

# fishery selex (3 time blocks)
plot(rep$fsh_slx[1,,1], type = 'l')
lines(rep$fsh_slx[30,,1], type = 'l')
lines(rep$fsh_slx[50,,1], type = 'l')

# recruitment
plot(rep$pred_rec, type = 'l')

# survey selex (3 time blocks)
plot(rep$srv_slx[1,,1], type = 'l')
lines(rep$srv_slx[30,,1], type = 'l')
lines(rep$srv_slx[50,,1], type = 'l')

