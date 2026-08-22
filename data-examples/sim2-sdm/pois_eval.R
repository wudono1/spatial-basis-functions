################################################################################
# Section 4.2 Simulation study: Model fitting and uncertainty propagation
# Calculate results using the fitted GLMs
# compute bias and 95% coverage on downscaled covariates B1
# also compute RMSE and 95% predictive coverage across GLM predicted response
################################################################################
rm(list = ls())

inpath1 <- "../outputs/sim1-downscaling/"
inpath2 <- "../outputs/sim2-pois/"
outpath <- "../outputs/sim2-pois/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
load(paste0(inpath1, "covar_dat.rdata"))
load(paste0(inpath1, "sim_gpbau_pred.RData"))
load(paste0(inpath2, "pois_rep_data.rdata"))
load(paste0(inpath2, "pois_rep_out.rdata"))

source("directory/to/source_functions.R")

nphi <- dim(fine_xx)[2]
nsill <- dim(fine_xx)[3]
# nrep <- length(rep_data)
nrep <- length(pois_rep_out)


#-------------------------------------------------------------------------------
# Parallel computation settings
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(doFuture)
})

ncores <- parallel::detectCores() - 1
nrange <- length(range_vals)
nsill <- length(sill_vals)

registerDoFuture()
plan(multisession, workers = min(ncores, nrange * nsill))


#-------------------------------------------------------------------------------
# Evaluate models
#-------------------------------------------------------------------------------
pois_mods_metrics <- vector("list", length = nrep)

# conputes B0 + B1x (no uncertainty)
hfunc1 <- function(x, be) cbind(1, x) %*% be

# computes B0 + B1x + d (with undertainty)
hfunc2 <- function(x, be, delta) cbind(1, x) %*% be + delta

cat("-------------------------------------------------------------------------\n")
cat("evaluating models\n")
cat(paste0(nrep, " replications\n\n"))

for (r in 1:nrep) {
  
  cat(paste0("start the ", r, "-th replication.\n"))
  cat(paste0("evaluating...\n"))
  
  fine_test_xx <- rep_data[[r]]$test$fine_test_xx
  fine_test_yy <- rep_data[[r]]$test$fine_test_yy
  fine_test_zz <- rep_data[[r]]$test$fine_test_zz
  
  train_idx <- rep_data[[r]]$train_idx
  true_be <- rep_data[[r]]$true_be
  
  pois_rep_out_r <- pois_rep_out[[r]]
  
  mod_eval_runtime <- system.time(
    pois_mods_metrics_r <- foreach(i = 1:nphi, .options.future = list(seed = TRUE)) %:%
      foreach(j = 1:nsill, .options.future = list(seed = TRUE)) %dofuture% 
      {
        evalMod(pois_rep_out_r[[i]][[j]], true_be,
                fine_test_zz[, i, j],
                log(fine_test_yy[, i, j]),
                fine_test_xx[, i, j],
                xx_pred[[i]][[j]][-train_idx, ],
                hfunc1,
                hfunc2)
      }
  )
  
  pois_mods_metrics[[r]] <- pois_mods_metrics_r
  
  cat(paste0("done the ", r, "-th replication.\n"))
  cat(paste0("running time: ", round(mod_eval_runtime[3]/60), " minutes\n"))
  cat("-------------------------------------------------------------------------\n")
  
}

save(pois_mods_metrics, file = paste0(outpath, "pois_eval_out.rdata"))


