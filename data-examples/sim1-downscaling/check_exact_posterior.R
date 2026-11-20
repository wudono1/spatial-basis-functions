################################################################################
# Exact-posterior check for Step 2 (basis-function GPB)
#
# Idea: psi's posterior (with gamma and sigmasq integrated out) is a formula in
# one variable, so we can compute it directly on a grid. That gives an "answer
# key" to compare against the MCMC draws. If the histogram of draws matches the
# grid density, the sampler is drawing from the right distribution.
################################################################################
rm(list = ls())

# same folder layout as the other sim1 scripts (run from data-examples/)
inpath  <- "../outputs/sim1-downscaling/"
outpath <- "../outputs/sim1-downscaling/exact-posterior-check/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# fitted chains from downscale_gpbau.R
# gives: out_list, range_vals, sill_vals, gpb_runtime, findint_runtime
load(paste0(inpath, "sim_gpbau_new_fit.RData"))

# datasets to check, as (i, j) = (range index, sill index)
# (4,3), (8,5): non-covering fits;  (1,1), (8,4): boundary-start fits that covered
check_sets <- list(c(4, 3), c(8, 5), c(1, 1), c(8, 4))


#-------------------------------------------------------------------------------
# building log posterior of psi, up to a constant:
#   log p(psi | X_B) = -1/2 log|K_B| - 1/2 log|R' K_B^-1 R|
#                      - (a + (n-p)/2) log(b + S2/2)
# S4.2
#-------------------------------------------------------------------------------
log_post_psi <- function(psi, xB, RB, PhiB, D, a, b) {
  
  # coarse correlation matrix K_B = Phi_B K_psi Phi_B'
  KB <- PhiB %*% exp(-D / psi) %*% t(PhiB)
  
  # GLS fit of the mean: gam_hat = (R' K_B^-1 R)^-1 R' K_B^-1 x_B
  KBinv_RB <- solve(KB, RB)                           # solve x = K_B^-1 * R
  KBinv_xB <- solve(KB, xB)                           # solve y = K_B^-1 * X_B
  Q        <- crossprod(RB, KBinv_RB)                 # R' K_B^-1 R
  gam_hat  <- solve(Q, crossprod(RB, KBinv_xB))
  
  # generalized residual sum of squares S2 = e' K_B^-1 e
  ee <- xB - RB %*% gam_hat
  S2 <- as.numeric(crossprod(ee, solve(KB, ee)))
  
  # assemble the log posterior
  n <- length(xB); p <- ncol(RB)
  lp <- -0.5 * as.numeric(determinant(KB, logarithm = TRUE)$modulus) -
    0.5 * as.numeric(determinant(Q,  logarithm = TRUE)$modulus) -
    (a + (n - p) / 2) * log(b + S2 / 2)
  
  c(lp = lp, S2 = S2)
}

# log density of InvGamma(shape, scale), used for the σ^2 marginal
ldinvgamma <- function(x, shape, scale) {
  shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
}


results <- list()

for (ij in check_sets) {
  i <- ij[1]; j <- ij[2]
  obj <- out_list[[i]][[j]]
  
  #-----------------------------------------------------------------------------
  # Pull the exact data, fixed inputs, and prior that this fit used
  #-----------------------------------------------------------------------------
  xB   <- obj$coarse_data$coarse_resp
  RB   <- as.matrix(obj$coarse_data$coarse_covar)
  mp   <- obj$mod_spec$mod_param
  PhiB <- mp$aggmat %*% mp$Phi                 # Phi_B = W_B Phi
  D    <- mp$centroid_distmat                  # distances between basis centers
  u <- obj$mod_spec$prior$psi_unif[1]          # psi ~ Uniform(u, v)
  v <- obj$mod_spec$prior$psi_unif[2]
  a <- obj$mod_spec$prior$sigmasq_invgamma[1]  # sigmasq ~ InvGamma(a, b)
  b <- obj$mod_spec$prior$sigmasq_invgamma[2]
  n <- length(xB); p <- ncol(RB)
  
  #-----------------------------------------------------------------------------
  # Exact psi posterior: evaluate on a fine grid, then normalize
  #-----------------------------------------------------------------------------
  # midpoint grid of G values strictly inside (u, v); each cell has width dpsi
  G        <- 1000
  dpsi     <- (v - u) / G
  psi_grid <- u + (seq_len(G) - 0.5) * dpsi
  
  # log posterior (and S2) at every grid point -> 2 x G matrix
  ev <- sapply(psi_grid, log_post_psi, xB = xB, RB = RB, PhiB = PhiB, D = D,
               a = a, b = b)
  lp <- ev["lp", ]; S2 <- ev["S2", ]
  
  # exponentiate (subtract max first to avoid underflow), then scale so the
  # density integrates to 1 over the grid
  w        <- exp(lp - max(lp))
  dens_psi <- w / sum(w * dpsi)
  cdf_psi  <- cumsum(dens_psi * dpsi)          # for exact quantiles
  
  #-----------------------------------------------------------------------------
  # MCMC draws to be tested (thinned chains if available)
  #-----------------------------------------------------------------------------
  ps <- if (!is.null(obj$post_sams_thin)) obj$post_sams_thin else obj$post_sams
  psi_draws <- as.numeric(ps$phi)              # key "phi" holds psi
  s2_draws  <- as.numeric(ps$sigmasq)
  
  #-----------------------------------------------------------------------------
  # Exact sigmasq posterior:
  #   given psi, sigmasq ~ InvGamma(a + (n-p)/2, b + S2(psi)/2), so the marginal
  #   is the average of those InvGammas weighted by p(psi | X_B) on the grid
  #-----------------------------------------------------------------------------
  s2_grid <- seq(0.5 * min(s2_draws), 1.5 * max(s2_draws), length.out = 1000)
  ds2     <- diff(s2_grid)[1]
  dens_s2 <- sapply(s2_grid, function(s)
    sum(exp(ldinvgamma(s, a + (n - p) / 2, b + S2 / 2)) * dens_psi * dpsi))
  cdf_s2  <- cumsum(dens_s2 * ds2)
  
  #-----------------------------------------------------------------------------
  # Numerical comparison: exact vs MCMC quantiles (2.5%, 50%, 97.5%)
  #-----------------------------------------------------------------------------
  probs <- c(0.025, 0.5, 0.975)
  
  # exact quantile = first grid point where the CDF reaches the target prob
  q_psi_exact <- sapply(probs, function(pp) psi_grid[which(cdf_psi >= pp)[1]])
  q_s2_exact  <- sapply(probs, function(pp) s2_grid[which(cdf_s2 >= pp)[1]])
  
  results[[length(results) + 1]] <- data.frame(
    i = i, j = j,
    parameter = rep(c("psi", "sigmasq"), each = 3),
    prob      = rep(probs, 2),
    exact     = c(q_psi_exact, q_s2_exact),
    mcmc      = c(quantile(psi_draws, probs), quantile(s2_draws, probs)),
    truth     = rep(c(range_vals[i], sill_vals[j]), each = 3)
  )
  
  #-----------------------------------------------------------------------------
  # Visual comparison: MCMC histogram (grey) vs exact density (red),
  # dashed line = true value used to simulate the data
  #-----------------------------------------------------------------------------
  png(paste0(outpath, sprintf("exact_vs_mcmc_i%d_j%d.png", i, j)),
      width = 1000, height = 450, pointsize = 16)
  par(mfrow = c(1, 2))
  
  hist(psi_draws, breaks = 40, freq = FALSE, col = "grey85", border = "white",
       xlim = c(u, v), xlab = "psi",
       main = sprintf("psi  (i=%d, j=%d)", i, j))
  lines(psi_grid, dens_psi, col = "red", lwd = 2)
  abline(v = range_vals[i], lty = 2)
  
  hist(s2_draws, breaks = 40, freq = FALSE, col = "grey85", border = "white",
       xlab = "sigmasq", main = sprintf("sigmasq  (i=%d, j=%d)", i, j))
  lines(s2_grid, dens_s2, col = "red", lwd = 2)
  abline(v = sill_vals[j], lty = 2)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Print and save the quantile table
#-------------------------------------------------------------------------------
results <- do.call(rbind, results)
results$exact <- round(results$exact, 3)
results$mcmc  <- round(results$mcmc, 3)
print(results, row.names = FALSE)

write.csv(results, paste0(outpath, "exact_vs_mcmc_quantiles.csv"), row.names = FALSE)