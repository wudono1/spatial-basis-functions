################################################################################
# Statistical downscaling using GPB
################################################################################
rm(list = ls())

inpath  <- "../outputs/sim1-downscaling/"
outpath <- "../outputs/sim1-downscaling/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


library(ggplot2)
library(spds)

#-------------------------------------------------------------------------------
# Load coarse-resolution covariate data from simulate_covariate_data.r
#-------------------------------------------------------------------------------
load(paste0(inpath, "covar_dat.rdata"))

#-------------------------------------------------------------------------------
# Parallel computation settings
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(doFuture)
})

registerDoFuture()

ncores <- parallel::detectCores() - 1
nrange <- length(range_vals)
nsill <- length(sill_vals)

plan(multisession, workers = min(ncores, nrange*nsill))

#-------------------------------------------------------------------------------
# GPB model settings
#-------------------------------------------------------------------------------
Nd <- nrow(fine_grid)
Md <- nrow(coar_grid)

# modelling covariance between BAUs -- OLD (Eq. 1) METHOD ONLY.
# The basis-function path never touches this: fitGPBAU_new and
# predGPBAU2_new/predGPBAU3_new build K_B from Phi + centroid_distmat instead.
# Left commented out because it is an Nd x Nd matrix (~72 MB at Nd = 3000) that
# would otherwise be copied to every multisession worker for nothing.
# Uncomment (and add it back to mod_param) only if reverting to gpbau().
# ⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄
# bau_distmat <- as.matrix(dist(fine_grid[, c("lon", "lat")]))

# aggregate matrix W_B st. X_B = W_B * X_S (MEE paper 2.2.2.)
aggmat <- agg(fine_grid, coar_grid)

# Set up
prior <- list(
  #### WHY ARE THESE PARAMS CHOSEN?????
  
  
  # bounds u, v, uniform prior on param ψ (Supp S4)
  phi_unif = c(0.1/3, 1),
  
  # params c(u, v) of inverse gamma on param σ^2 (Supp S4)
  sigmasq_invgamma = c(3, 2)
)

# SD of Gaussian random walk (step 1 of MCMC algo in Supp S4)
# SD of distribution glogit(ψ′) ~ N( glogit(ψ⁽ᵗ⁻¹⁾), τ² ), where
# glogit(y) = log((y − u_ψ)/(v_ψ − y))
tuning <- list(se_phi = 2.5)


# list of fixed inputs (aggregation matrix, distance matrix, etc)
mod_param <- list(aggmat = aggmat, Phi = Phi,
                  centroid_distmat = centroid_distmat)

# settings for MCMC algorithm in S$ of supporting information doc
mcmc_settings <- list(niter = 8000, nreport = 1000)

# coar_covar: R_B, M-D x 1 matrix of 1's
# specified in Supp S4 as R̃ =[1_ND, R] where R is N_D x p matrix of covariates
# no BAU-level regressors, so R_B is just matrix of ones
coar_covar <- as.matrix(rep(1, Md))


#-------------------------------------------------------------------------------
# First check to see if r >= MD (unable to compute if not)
# r = number of basis functions, M_D = dimensions of X_B
# need to invert the M_D by M_D matrix K_B(Ψ) = Ψ_B * K_Ψ * Ψ_B^T
# Ψ_β = W_B * ɸ (M_d x r)
# if r < M_d, then K_B (M_D x M_D) is not full rank, and thus not invertible
# Note: r >= M_d does not mean K_B must be full rank (necessary but not sufficient)
#-------------------------------------------------------------------------------

r_basis  <- ncol(Phi)
PsiB_chk <- aggmat %*% Phi  # Ψ_β = W_B * ɸ (M_d x r)
cat(sprintf("Nd = %d, Md = %d, r = %d, rank(W_B Phi) = %d\n",
            Nd, Md, r_basis, qr(PsiB_chk)$rank))
stopifnot(                             # stop if any of 4 conditions are false
  nrow(Phi) == Nd,                     # ɸ: 1 row per BAU (to multiply W_B * ɸ)
  nrow(centroid_distmat) == r_basis,   # K_Ψ is size (r x r)
  r_basis >= Md,                       # r_Basis >= M_d
  qr(PsiB_chk)$rank == Md              # check K_b is full rank (invertibility)
)


#-------------------------------------------------------------------------------
# Fit GPB models (Step 2 in the manuscript)
#-------------------------------------------------------------------------------
set.seed(42)

# Grid search to find initial values of range, variance for MCMC model fitting 
# Evaluating log likelihood over a grid of range and variance values
eps <- 1e-2 * diff(prior$phi_unif)
phi_vals <- seq(prior$phi_unif[1] + eps, prior$phi_unif[2] - eps, length = 50) 
n_sigma_vals <- 5

# for each dataset ({len(range_vals) x len(sill_vals)} unique datasets), 
# find (Ψ(0), σ^2(0)), the starting value pairs w/ 
# max log-likelihood (using findinit())
findint_runtime <- system.time(
  init_list <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture% {
      # findinit -> findinit_new: grid-searches the BASIS-FUNCTION likelihood so
      # the starting values match the model we actually fit below. Also much
      # faster: it builds an r x r kernel per psi instead of an Nd x Nd one.
      findinit_new(coar_xx[, i, j], mod_param, phi_vals, n_sigma_vals)
    }
)

# rearrange the found (Ψ(0), σ^2(0)) pairs into talbes
phi_init <- array(NA, dim = c(length(range_vals), length(sill_vals)))
sigmasq_init <- array(NA, dim = c(length(range_vals), length(sill_vals)))

for (i in seq_along(range_vals)) {
  for (j in seq_along(sill_vals)) {
    phi_init[i, j] <- init_list[[i]][[j]]$phi_init
    sigmasq_init[i, j] <- init_list[[i]][[j]]$sigmasq_init
  }
}

# Fit GPB
cat("----------------------------------------------------\n")
cat("Fitting basis-function GPB models\n")


gpb_runtime <- system.time(
  out_list <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture%
    {
      # coar_covar <- as.matrix(rep(1, Md)) from previous line
      # interpretation: modeling covariate field as const mean with spatial GP residual
      #
      # gpbau -> gpbau_new: dispatches to fitGPBAU_new, which builds
      # K_B = Psi_B K_psi Psi_B' from Phi + centroid_distmat instead of
      # W_B exp(-d/phi) W_B'. Arguments are otherwise IDENTICAL -- the wrapper
      # takes the same prior/tuning/starting keys, and reads the extra basis
      # inputs out of mod_param. To fit the OLD model instead, swap this one
      # name back to gpbau() (mod_param is a superset, so it still works).
      gpbau_new(coarse_resp = coar_xx[,i,j],
                coarse_covar = coar_covar,
                prior = prior,
                starting = list(phi = phi_init[i, j],
                                sigmasq = sigmasq_init[i, j],
                                regcoef = matrix(c(mean(coar_xx[,i,j])))),
                tuning = tuning,
                mod_param = mod_param,
                mcmc_settings = mcmc_settings,
                verbose = FALSE)
    }
)

cat("Done!\n")
# cat(paste0("Running time: ", round(gpb_runtime[3]/60, 3), " minutes.\n"))
cat("----------------------------------------------------\n")


# Thinning the MCMC chain if needed
for (i in seq_along(range_vals)) {
  for (j in seq_along(sill_vals)) {
    out <- out_list[[i]][[j]]
    out <- thinning(out, nburn = 3000, nthin = 10)
    out_list[[i]][[j]] <- out
  }
}


#-------------------------------------------------------------------------------
# Save the fitted objects for parameter comparison
#-------------------------------------------------------------------------------
# out_list now holds objects of class "gpbau_new". Each element's post_sams has
# chains named phi (= the range psi), sigmasq, and beta (= gamma-tilde).
# The data-generating truth for element [i,j] is psi = range_vals[i] and
# sigmasq = sill_vals[j], so these chains can be checked for coverage directly.
save(out_list, range_vals, sill_vals, gpb_runtime, findint_runtime,
     file = paste0(outpath, "sim_gpbau_new_fit.RData"))





#-------------------------------------------------------------------------------
# Posterior prediction (Step 3 in the manuscript)
# does not produce Berkson-error GLM yet
#-------------------------------------------------------------------------------
# RE-ENABLED for the basis-function model. predict() dispatches on class
# "gpbau_new" -> predict.gpbau_new() -> predGPBAU3_new() (independent, default)
# or predGPBAU2_new() (joint). Those build C_S = sigmasq * Phi K_psi Phi' and
# condition in alpha-space, so no Nd x Nd matrix is ever formed.
#
# The OUTPUT CONTRACT is unchanged: xx_pred[[i]][[j]] is an Nd x K matrix, so
# downscale_comp.R and the sim2-sdm scripts need no edits.
set.seed(42)
fine_covar <- as.matrix(rep(1, Nd))

# generating posterior DISTRIBUITON of the field
# getting many draws of covariate vector Xs to create a distribution
xx_pred <- vector("list", length = nrange)

runtime <- system.time(
  xx_pred <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture%
    {
      predict(out_list[[i]][[j]], fine_covar, nreport = 100)
    }
)

save(xx_pred, file = paste0(outpath, "sim_gpbau_pred.RData"))


#-------------------------------------------------------------------------------
# Correctness check: with delta(S) = 0 the constraint W_B X_s = X_B holds
# EXACTLY for every draw (because Psi_B L = I in the alpha-space correction).
# Any appreciable deviation means the conditioning step is wrong.
# NOTE: this is exact only for type = "joint"; the independent predictor draws
# each location from its marginal, so it does NOT reproduce the constraint.
#-------------------------------------------------------------------------------
chk <- predict(out_list[[1]][[1]], fine_covar, type = "joint", nreport = 1e6)
cat(sprintf("joint-draw aggregation error (should be ~1e-12): %.3e\n",
            max(abs(aggmat %*% chk[, 1] - coar_xx[, 1, 1]))))




