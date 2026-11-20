################################################################################
# Refit the 4 checked datasets with a 10x longer chain
# Rebuilds the data, prior, tuning, and starting values exactly as in
# downscale_gpbau.R, so the ONLY things that change are the chain length and
# which datasets are fit.
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

# datasets to check, as (i, j) = (range index, sill index)
# (4,3), (8,5): non-covering fits;  (1,1), (8,4): boundary-start fits that covered
check_sets_range <- c(4, 8, 1, 8)
check_sets_sill <- c(3, 5, 1, 4)

Nd <- nrow(fine_grid)
Md <- nrow(coar_grid)

# aggregate matrix W_B st. X_B = W_B * X_S (MEE paper 2.2.2.)
aggmat <- agg(fine_grid, coar_grid)

# Setup
prior <- list(
  
  # bounds u, v, uniform prior on param ψ (Supp S4)
  psi_unif = c(0.1/3, 1),
  
  # params c(u, v) of inverse gamma on param σ^2 (Supp S4)
  sigmasq_invgamma = c(3, 2)
)

# SD of Gaussian random walk (step 1 of MCMC algo in Supp S4/5)
# distr. glogit(ψ′) ~ N( glogit(ψ⁽ᵗ⁻¹⁾), τ² )
tuning <- list(se_psi = 2.5)


# list of fixed inputs (aggregation matrix, distance matrix, etc)
mod_param <- list(aggmat = aggmat, Phi = Phi,
                  centroid_distmat = centroid_distmat)

# settings for MCMC algorithm in S4 of supporting information doc
mcmc_settings <- list(niter = 80000, nreport = 10000)

# From Supp S4: coar_covar is R̃ =[1_ND, R] (covariates)
# matrix of ones in this case
coar_covar <- as.matrix(rep(1, Md))


#-------------------------------------------------------------------------------
# First check to see if r >= MD (unable to compute if not)
# r = number of basis functions, M_D = dimensions of X_B
# need to invert the M_D by M_D matrix K_B(Ψ) = ɸ_B * K_Ψ * ɸ_B^T
# ɸ_β = W_B * ɸ (M_d x r)
# if r < M_d, then K_B (M_D x M_D) is not full rank, and thus not invertible
# Note: r >= M_d does not mean K_B must be full rank (necessary but not sufficient)

# check if r > M_D
#-------------------------------------------------------------------------------

r_basis  <- ncol(Phi)
PsiB_chk <- aggmat %*% Phi  # ɸ_β = W_B * ɸ (M_d x r)
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
eps <- 1e-2 * diff(prior$psi_unif)
psi_vals <- seq(prior$psi_unif[1] + eps, prior$psi_unif[2] - eps, length = 50) 
n_sigma_vals <- 5

# for each dataset ({len(range_vals) x len(sill_vals)} unique datasets), 
# find (Ψ(0), σ^2(0)), the starting value pairs w/ 
# max log-likelihood (using findinit_new())
findint_runtime <- system.time(
  init_list <- foreach(i = seq_along(range_vals), .options.future = list(packages = "spds", seed = TRUE)) %:%
    foreach(j = seq_along(sill_vals), .options.future = list(packages = "spds", seed = TRUE)) %dofuture% {
      
      # grid search basis function likelihood so starting vals match the 
      # model we actually fit below
      findinit_new(coar_xx[, i, j], mod_param, psi_vals, n_sigma_vals)
    }
)

# rearrange the found (Ψ(0), σ^2(0)) pairs into tables
psi_init <- array(NA, dim = c(length(range_vals), length(sill_vals)))
sigmasq_init <- array(NA, dim = c(length(range_vals), length(sill_vals)))

for (i in seq_along(range_vals)) {
  for (j in seq_along(sill_vals)) {
    psi_init[i, j] <- init_list[[i]][[j]]$psi_init
    sigmasq_init[i, j] <- init_list[[i]][[j]]$sigmasq_init
  }
}

# Fit GPB
cat("----------------------------------------------------\n")
cat("Fitting basis-function GPB models\n")


gpb_runtime <- system.time(
  out_list <- foreach(k = seq_along(check_sets_range), .options.future = list(packages = "spds", seed = TRUE)) %dofuture%
    {
      # coar_covar <- as.matrix(rep(1, Md)) from previous line
      # interpretation: modeling covariate field as const mean with spatial GP residual
      
      spds:::gpbau_new_VARCHANGE(coarse_resp = coar_xx[,check_sets_range[k],check_sets_sill[k]],
                                 coarse_covar = coar_covar,
                                 prior = prior,
                                 starting = list(psi = psi_init[check_sets_range[k],check_sets_sill[k]],
                                                 sigmasq = sigmasq_init[check_sets_range[k],check_sets_sill[k]],
                                                 regcoef = matrix(c(mean(coar_xx[,check_sets_range[k],check_sets_sill[k]])))),
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
for (k in seq_along(check_sets_range)) {
  out <- out_list[[k]]
  out <- thinning(out, nburn = 3000, nthin = 10)
  out_list[[k]] <- out

}

# out_list is FLAT: out_list[[k]] is the fit for dataset
# (check_sets_range[k], check_sets_sill[k]), so save those index vectors too
save(out_list, check_sets_range, check_sets_sill, range_vals, sill_vals,
     gpb_runtime, findint_runtime,
     file = paste0(outpath, "sim_gpbau_new_refit.RData"))








# inpath  <- "../outputs/sim1-downscaling/"
# outpath <- "../outputs/sim1-downscaling/"
# 
# library(ggplot2)
# 
# load(paste0(inpath, "sim_gpbau_new_refit.RData"))
# # gives: out_list, check_sets_range, check_sets_sill, range_vals, sill_vals,
# #        gpb_runtime, findint_runtime

# write recovery outputs to their own folder so the 40-dataset results in
# outputs/sim1-downscaling/ are not overwritten
recov_outpath <- paste0(outpath, "refit-long/")
dir.create(recov_outpath, recursive = TRUE, showWarnings = FALSE)

nrange <- length(check_sets_range)   # number of refit datasets (4)

#-------------------------------------------------------------------------------
# Build a long-format recovery table: one row per (i, j, parameter)
#-------------------------------------------------------------------------------
res <- list()
k <- 0

for (i in seq_len(nrange)) {
    
  obj <- out_list[[i]]
  
  # thinning() stores draws in post_sams_thin; fall back if not thinned
  ps <- if (!is.null(obj$post_sams_thin)) obj$post_sams_thin else obj$post_sams
  
  # NOTE: the list key is still "phi" for backward compatibility, but the
  # quantity it holds is the range parameter psi.
  chains <- list(psi     = as.numeric(ps$phi),
                 sigmasq = as.numeric(ps$sigmasq),
                 gamma0  = as.numeric(ps$beta))
  
  
  # true values
  # check_sets_range/sill hold INDICES into range_vals/sill_vals
  truths <- list(psi     = range_vals[check_sets_range[i]],
                 sigmasq = sill_vals[check_sets_sill[i]],
                 gamma0  = 0)
  
  for (par in names(chains)) {
    
    chain <- chains[[par]]
    truth <- truths[[par]]
    qq <- quantile(chain, c(0.025, 0.975))
    
    # effective sample size if coda is available (MH chains are autocorrelated)
    ess <- if (requireNamespace("coda", quietly = TRUE)) {
      as.numeric(coda::effectiveSize(chain))
    } else NA_real_
    
    k <- k + 1
    res[[k]] <- data.frame(
      i = check_sets_range[i], j = check_sets_sill[i],
      parameter   = par,
      truth       = truth,
      post_mean   = mean(chain),
      post_median = median(chain),
      lower95     = as.numeric(qq[1]),
      upper95     = as.numeric(qq[2]),
      covered     = truth >= qq[1] && truth <= qq[2],
      bias        = mean(chain) - truth,
      rel_bias    = if (truth != 0) (mean(chain) - truth) / truth else NA_real_,
      ci_width    = as.numeric(qq[2] - qq[1]),
      ess         = ess,
      accept_rate = mean(obj$mh_data),
      n_draws     = length(chain),
      row.names   = NULL
    )
  }
}

recov <- do.call(rbind, res)

#-------------------------------------------------------------------------------
# Headline numbers
#-------------------------------------------------------------------------------
# NOTE: only 4 hand-picked datasets (2 chosen BECAUSE they missed), so this is a
# description of these fits, not a coverage estimate
cat("\n=== 95% credible interval coverage (4 selected datasets; descriptive only) ===\n")
print(aggregate(covered ~ parameter, data = recov, FUN = mean))

cat("\n=== mean relative bias of posterior mean ===\n")
print(aggregate(rel_bias ~ parameter, data = recov,
                FUN = function(x) mean(x, na.rm = TRUE)))

cat("\n=== MH acceptance rate (want roughly 0.20-0.40) ===\n")
acc <- unique(recov[, c("i", "j", "accept_rate")])
print(summary(acc$accept_rate))

cat("\n=== effective sample size ===\n")
print(aggregate(ess ~ parameter, data = recov,
                FUN = function(x) mean(x, na.rm = TRUE)))

# any fits where BOTH psi and sigmasq missed
miss <- recov[!recov$covered & recov$parameter %in% c("psi", "sigmasq"), ]
if (nrow(miss) > 0) {
  cat("\n=== non-covering fits ===\n")
  print(miss[, c("i", "j", "parameter", "truth", "post_mean",
                 "lower95", "upper95")])
}

#-------------------------------------------------------------------------------
# Caterpillar plots: posterior interval vs truth, one panel per parameter
#-------------------------------------------------------------------------------
recov$label <- paste0("r", recov$i, "-s", recov$j)

p_cat <- ggplot(subset(recov, parameter %in% c("psi", "sigmasq")),
                aes(x = label)) +
  geom_linerange(aes(ymin = lower95, ymax = upper95, colour = covered)) +
  geom_point(aes(y = post_mean, colour = covered), size = 1) +
  geom_point(aes(y = truth), shape = 4, size = 2) +   # x = truth
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("TRUE" = "grey30", "FALSE" = "red")) +
  labs(x = "range-sill combination", y = "value",
       title = "Posterior mean and 95% CI vs truth (x = truth)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

png(paste0(recov_outpath, "recovery_intervals.png"),
    width = 1000, height = 700, pointsize = 16)
print(p_cat)
dev.off()

#-------------------------------------------------------------------------------
# Estimated vs true scatter (should sit on the 1:1 line)
#-------------------------------------------------------------------------------
p_scat <- ggplot(subset(recov, parameter %in% c("psi", "sigmasq")),
                 aes(x = truth, y = post_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), alpha = 0.4) +
  geom_point() +
  facet_wrap(~ parameter, scales = "free") +
  labs(x = "true value", y = "posterior mean",
       title = "Parameter recovery (dashed = 1:1)") +
  theme_bw()

png(paste0(recov_outpath, "recovery_scatter.png"),
    width = 1000, height = 500, pointsize = 16)
print(p_scat)
dev.off()

#-------------------------------------------------------------------------------
# Trace plots for all refit chains -- convergence sanity check
#-------------------------------------------------------------------------------
png(paste0(recov_outpath, "recovery_traces.png"),
    width = 1000, height = 250 * length(out_list), pointsize = 16)
par(mfrow = c(length(out_list), 2))   # one row per dataset: psi | sigmasq
for (k in seq_along(out_list)) {
  psi_true <- range_vals[check_sets_range[k]]
  s2_true  <- sill_vals[check_sets_sill[k]]
  ps <- out_list[[k]]$post_sams_thin
  plot(as.numeric(ps$phi), type = "l",
       main = paste0("psi trace: range=", psi_true, ", sill=", s2_true),
       ylab = "psi", xlab = "draw")
  abline(h = psi_true, col = "red", lwd = 2)
  plot(as.numeric(ps$sigmasq), type = "l",
       main = paste0("sigmasq trace: range=", psi_true, ", sill=", s2_true),
       ylab = "sigmasq", xlab = "draw")
  abline(h = s2_true, col = "red", lwd = 2)
}
dev.off()

write.csv(recov, paste0(recov_outpath, "recovery_table.csv"), row.names = FALSE)