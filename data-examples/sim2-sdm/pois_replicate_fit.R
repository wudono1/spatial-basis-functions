################################################################################
# Section 4.2 Simulation study: Model fitting and uncertainty propagation
# Fit GLMs to simulated data
################################################################################
rm(list = ls())

inpath1 <- "../outputs/sim1-downscaling/"
inpath2 <- "../outputs/sim2-pois/"
outpath <- "../outputs/sim2-pois/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

set.seed(42)


#-------------------------------------------------------------------------------
# Load covariate data and response data
#-------------------------------------------------------------------------------
# load true covariate
load(paste0(inpath1, "covar_dat.rdata"))
nphi <- dim(fine_xx)[2]
nsill <- dim(fine_xx)[3]

# load downscaled covariate 
load(paste0(inpath1, "sim_gpbau_pred.RData"))

# load response data
load(paste0(inpath2, "pois_rep_data.rdata"))
# nrep <- length(rep_data)
nrep <- min(3, length(rep_data))

#-------------------------------------------------------------------------------
# source model-fit function
#-------------------------------------------------------------------------------
source("sim2-sdm/source_functions.R")


#-------------------------------------------------------------------------------
# Parallel computation settings
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(doFuture)
})

ncores <- parallel::detectCores() - 1
nrange <- length(range_vals)
nsill <- length(sill_vals)

# registerDoFuture()
plan(multisession, workers = min(ncores, nrange * nsill))


#-------------------------------------------------------------------------------
# Fit the model to normal data
#-------------------------------------------------------------------------------
pois_rep_out <- vector("list", length = nrep)

cat("-------------------------------------------------------------------------\n")
cat("fitting GLMs to count data\n")
cat(paste0(nrep, " replications\n\n"))

for (r in 1:nrep) {
  
  cat(paste0("start the ", r, "-th replication.\n"))
  cat(paste0("fitting...\n"))
  
  fine_train_zz <- rep_data[[r]]$train$fine_train_zz
  fine_train_xx <- rep_data[[r]]$train$fine_train_xx
  train_idx <- rep_data[[r]]$train_idx
  
  pois_runtime <- system.time(
    pois_out_list_r <- foreach(i = 1:nphi, .options.future = list(seed = TRUE)) %:%
      foreach(j = 1:nsill, .options.future = list(seed = TRUE)) %dofuture% 
      {
        fitGLMPois(fine_train_zz[, i, j], 
                   fine_train_xx[, i, j], 
                   xx_pred[[i]][[j]][train_idx, ],
                   BB = 300)
      }
  )
  
  pois_rep_out[[r]] <- pois_out_list_r
  
  cat(paste0("done the ", r, "-th replication.\n"))
  cat(paste0("running time: ", round(pois_runtime[3]/60), " minutes\n"))
  cat("-------------------------------------------------------------------------\n")
  
}

save(pois_rep_out, file = paste0(outpath, "pois_rep_out.rdata"))


                       