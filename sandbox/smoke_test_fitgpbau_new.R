################################################################################
# Smoke test for fitGPBAU_new() — basis-function GPB MCMC
# ------------------------------------------------------------------------------
# Purpose: check that the NEW covariance construction KB = PsiB * Kpsi * PsiB'
#          runs and stays positive-definite BEFORE wiring it into the real
#          experiment. This is a throwaway sandbox script — it is deliberately
#          OUTSIDE spds/ so it is not treated as part of the package.
#
# Before running, build/load the package from the repo root (corgi-examples/):
#     Rcpp::compileAttributes("spds")
#     devtools::load_all("spds")
#     exists("fitGPBAU_new")   # should be TRUE
#
# Then either source this file:
#     source("sandbox/smoke_test_fitgpbau_new.R")
# or step through it interactively.
#
# NOTE: this calls fitGPBAU_new() DIRECTLY (bypassing the gpbau() wrapper), so
#       the list keys are the raw C++-level names: u_phi, v_phi, regcoef, se_phi,
#       etc. — NOT the friendlier gpbau() keys (phi_unif, ...).
################################################################################

## ------------------------------------------------------------------ ##
## Knobs you can vary between runs (see the use-cases at the bottom).  ##
## ------------------------------------------------------------------ ##
Md    <- 6      # coarse blocks (= nn, number of coarse observations)
Nd    <- 24     # fine BAUs (must be a multiple of Md for this simple aggregator)
r     <- 3     # basis functions. r >= Md -> KB positive-definite; r < Md -> singular
niter <- 50     # MCMC iterations (small for a smoke test)
seed  <- 1      # RNG seed for reproducibility

run_smoke_test <- function(Md = 6, Nd = 24, r = 8, niter = 50, seed = 1,
                           se_phi = 1.0, psi_start = 0.5) {

  set.seed(seed)

  stopifnot(Nd %% Md == 0)          # simple block-average needs equal-size blocks
  per <- Nd / Md                    # BAUs per block

  ## aggregation matrix W_B : Md x Nd  (each row block-averages its BAUs)
  aggmat <- matrix(0, Md, Nd)
  for (j in 1:Md) aggmat[j, ((j - 1) * per + 1):(j * per)] <- 1 / per

  ## basis functions Phi : Nd x r  (random; full column rank is enough to test)
  Phi <- matrix(rnorm(Nd * r), Nd, r)

  ## centroid distances : r x r  (r knots placed on a line, 0..1)
  centroid_distmat <- as.matrix(dist(seq(0, 1, length.out = r)))

  ## dummy bau_distmat : Nd x Nd  (still READ by the C++, but unused — keep present)
  bau_distmat <- as.matrix(dist(seq(0, 1, length.out = Nd)))

  ## coarse response + design (intercept only -> p+1 = 1)
  xB <- rnorm(Md)
  RB <- matrix(1, Md, 1)

  ## the four input lists (raw C++ keys)
  prior <- list(u_phi = 0.01, v_phi = 2, u_sigmasq = 2, v_sigmasq = 1)
  starting <- list(regcoef = matrix(mean(xB)), sigmasq = 1, phi = psi_start)
  tuning <- list(se_phi = se_phi)
  mod_param <- list(bau_distmat = bau_distmat, aggmat = aggmat,
                    Phi = Phi, centroid_distmat = centroid_distmat)
  mcmc_settings <- list(niter = niter, nreport = max(1, niter %/% 5))

  ## quick sanity check of the covariance BEFORE calling C++:
  ## build KB in R and report its smallest eigenvalue (<= ~0 => not pos-def)
  PsiB <- aggmat %*% Phi
  Kpsi <- exp(-centroid_distmat / psi_start)
  KB <- PsiB %*% Kpsi %*% t(PsiB)
  min_eig <- min(eigen(KB, symmetric = TRUE, only.values = TRUE)$values)
  cat(sprintf("KB is %d x %d, rank(PsiB) = %d, min eigenvalue = %.3e\n",
              nrow(KB), ncol(KB), qr(PsiB)$rank, min_eig))
  if (min_eig <= 1e-8)
    cat("  -> KB is (numerically) singular; expect fitGPBAU_new to fail.\n")

  ## run the sampler
  out <- fitGPBAU_new(xB, RB, prior, starting, tuning, mod_param, mcmc_settings)
  out
}

## ------------------------------------------------------------------ ##
## Default run: r = 8 >= Md = 6, should succeed.                       ##
## ------------------------------------------------------------------ ##
out <- run_smoke_test(Md = Md, Nd = Nd, r = r, niter = niter, seed = seed)

cat("\n--- posterior sample structure ---\n")
str(out$post_sams)
cat("\nacceptance rate:", mean(out$mh_data), "\n")
