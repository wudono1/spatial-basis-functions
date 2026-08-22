################################################################################
# Master function: run all R scripts 
################################################################################
setwd("~/Documents/GitHub/corgi-examples/data-examples")

#-------------------------------------------------------------------------------
# Simulate fine-resolution covariate data
# The data will be used for both simulation studies 1 and 2
#-------------------------------------------------------------------------------
source("simulate_covariate_data.R")


#-------------------------------------------------------------------------------
# Statistical downscaling
# The downscaled data will be used for both simulation studies 1 and 2
#-------------------------------------------------------------------------------
source("sim1-downscaling/downscale_gpbau.R")


#-------------------------------------------------------------------------------
# Parameter-recovery check for the basis-function GPB model
# Verifies that the MCMC recovers the known simulation truth
# (psi = range_vals[i], sigmasq = sill_vals[j], gamma0 = 0) before any of the
# downstream results are trusted.
#-------------------------------------------------------------------------------
source("sim1-downscaling/downscale_param_recovery.R")


#-------------------------------------------------------------------------------
# Simulation study 1
#-------------------------------------------------------------------------------
# Produce results in Section 4.1
source("sim1-downscaling/downscale_comp.R")


#-------------------------------------------------------------------------------
# Simulation study 2 (the downscaling step is done in 'downscale_gpbau.R')
#-------------------------------------------------------------------------------
# Simulate response data
source("sim2-sdm/simulate_pois_data.R")

# Fit Poisson GLMs
source("sim2-sdm/pois_replicate_fit.R")

# Evaluate models
source("sim2-sdm/pois_eval.R")

# Produce results in Section 4.2
source("sim2-sdm/comp_figures.R")


#-------------------------------------------------------------------------------
# Case study
#-------------------------------------------------------------------------------
# Prepare data
source("real-bh-moss/survey_locs.R")
source("real-bh-moss/prepare_moss_data.R")

# Downscale soil temperature
source("real-bh-moss/downscale_soil_temp.R")

# Fit a Bernoulli Berkson-GLM
source("real-bh-moss/sdm_moss.R")

# Bayesian implementation of the Bernoulli Berkson-GLM
source("real-bh-moss/brms_sdm_moss.R")
source("real-bh-moss/nimble_sdm_moss.R")