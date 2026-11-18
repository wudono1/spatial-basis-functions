#' @title Aggregation matrix
#'
#' @description This function generates the aggregation matrix considered in the manuscript.
#'
#' @param fine_grid an \eqn{N_D \times 5} data frame; need more details...
#' @param coar_grid an \eqn{M_D \times 2} matrix of coarse-resolution coordinates
#'
#' @return An \eqn{M_D\times N_D} aggregation matrix.
#'   
#' @export
#' 
agg <- function(fine_grid, coar_grid) {
  
  Nd <- nrow(fine_grid)
  Md <- nrow(coar_grid)
  
  aggmat <- array(0, dim = c(Md, Nd))
  
  for (j in 1:Md) {
    j_idx <- which(fine_grid$par_label == j)
    aggmat[j, j_idx] <- 1 / length(j_idx)
  }
  
  aggmat
  
}

#' @title Thinniing of MCMC samples
#'
#' @description This function provides thinning of MCMC samples.
#'
#' @param object an object of class gpbau.
#' @param nburn an integer: number of burn-in samples.
#' @param nthin an integer: thinning degree. 
#'
#' @return An object of class \code{gpbau}
#'   
#' @export
#' 
thinning <- function(object, nburn, nthin) {

    niter <- object$mcmc_settings$niter
    idx <- seq(nburn + 1, niter, by = nthin)
    
    post_sams <- object$post_sams
    post_sams$beta <- t(post_sams$beta)
    post_sams_thin <- lapply(post_sams, function(x) x[idx, ])

    post_sams_thin$beta <- as.matrix(post_sams_thin$beta)

    object$post_sams_thin <- post_sams_thin

    object

}


#' @title Find initial values of GPB parameters
#'
#' @description This function finds initial values of GPB range and variance parameters
#'
#' @param xxb a vector of coarse-resolution data of length \eqn{M_D}.
#' @param mod_param a list of parameters needed for model fitting. Currently the function requires:
#'                  \code{aggmat}: an \eqn{M_D \times N_D} matrix that aggregates fine-resolution data to coarse-resolution data;
#'                  \code{bau_distmat}: an \eqn{N_D \times N_D} matrix about pairwise-locaiton distances.
#' @param phi_vals a vector of gridded values of GPB range parameter \code{phi}.
#' @param n_sigmasq an integer to create a grid of sigmasq.
#'
#' @return The return object is a list that comprises a matrix of log-likelihood values
#'         and a vector of gridded values of sigmasq
#'   
#' @export
#' 

# runs once per MCMC run to get starting state t(0) of MCMC chain. 
# Proposal at t(1) is centered at t(0) state
# note: phi = psi in this function
findinit <- function(xxb, mod_param, phi_vals, n_sigmasq) {
  
  # xxb = X_B, covar = R_B, aggmat = W_B, distmat = N_D x N_D distance matrix between BAU's
  data <- list(resp = xxb, covar = as.matrix(rep(1, length(xxb))),
               aggmat = mod_param$aggmat, distmat = mod_param$bau_distmat)
  
  # empty list of log-likelihood and σ^2 values to be filled
  ll_vals <- array(NA, dim = c(length(phi_vals), n_sigmasq))
  sigmasq_vals <- array(NA, dim = c(length(phi_vals), n_sigmasq))
  
  # Calculate GPB's log likelihood over a grid 
  
  # computing one starting value pair (ψ⁽⁰⁾, σ²⁽⁰⁾) for each true range, sill pair
  
  # outer loop: one pass per ψ candidate (pre-computed from phi_vals) to make
  # amount of candidate σ^2 values determined by n_sigmasq (5 in downscale_gpbau.R)
  for (k in seq_along(phi_vals)) {
    
    ww <- data$aggmat   # aggregation matrix W st. X_B = W * X_S
    
    # build BAU-level correlation matrix for sigma, aggregate to block-level
    bau_covmat <- exp(-data$distmat / phi_vals[k])  # K_ψ = exp(-dist_bau / Ψ), N_D x N_D corr matrix
    coar_covmat <- ww %*% bau_covmat %*% t(ww)      # K_B(Ψ) = W_B * K_Ψ * W_B^T
    
    # choosing the σ^2 models to try
    ratio <- diag(coar_covmat)[1]
    
    # _lb, _ub: lower and upper bounds of σ^2 values to try
    # sample variance: var(X_B) / [K_B(Ψ)]
    sigmasq_lb <- var(xxb) / ratio
    sigmasq_ub <- var(xxb) / ratio * 3  # upper bounds is 3x lower bound
    sigmasq_seq <- seq(sigmasq_lb, sigmasq_ub, length = n_sigmasq)
    sigmasq_vals[k, ] <- sigmasq_seq
    
    # inner loop: one pass per σ^2_{k, l} value to find log likelihood 
    # under (ψ⁽⁰⁾, σ²⁽⁰⁾) pair
    for (l in seq_along(sigmasq_seq)) {
      param <- c(phi_vals[k], sigmasq_seq[l], mean(xxb))
      ll_vals[k, l] <- gpbauLoglik(param, data)  # calc log-likelihood at (ψ_k, σ²_l)
      # cat(paste0("done the ", k, "-", l, "th eval for the ", i,"-", j, " iter.\n"))
    }
  }
  
  # find the (ψ⁽⁰⁾, σ²⁽⁰⁾) pair in ll_vals with largest log-likelihood
  idx <- as.numeric(which(ll_vals == max(ll_vals), arr.ind = TRUE))

  # output entire matrix of candidate value pairs (ψ⁽⁰⁾, σ²⁽⁰⁾) as well
  # as pair with highest log-likelihood
  list(ll_vals = ll_vals, 
       sigmasq_vals = sigmasq_vals,
       phi_init = phi_vals[idx[1]],
       sigmasq_init =  sigmasq_vals[idx[1],idx[2]])
  
}
  
# function to calculate log likelihood
gpbauLoglik <- function(param, data) {
  
  phi <- param[1]
  sigmasq <- param[2]
  be <- as.matrix(param[-(1:2)])
  
  yy <- data$resp
  rr <- data$covar
  ww <- data$aggmat
  distmat <- data$distmat
  
  bau_covmat <- sigmasq * exp(-distmat / phi)
  coar_covmat <- ww %*% bau_covmat %*% t(ww)
  ee <- yy - rr %*% be

  logdet <- as.numeric(determinant(coar_covmat)[[1]])
  ll <- -.5 * t(ee) %*% solve(coar_covmat, ee) - .5 * logdet  # eqn S4.1 in Supporting Info

  ll

}


#' @title Find initial values of basis-function GPB parameters
#'
#' @description Basis-function version of \code{\link{findinit}}. Identical grid search,
#'              but the coarse covariance is built as
#'              \eqn{K_B(\psi) = \Phi_B K_\psi \Phi_B^\top} with \eqn{\Phi_B = W_B \phi(S)^\top},
#'              instead of \eqn{W_B \exp(-d/\psi) W_B^\top} over BAU distances.
#'
#' @param xxb a vector of coarse-resolution data of length \eqn{M_D}.
#' @param mod_param a list requiring \code{aggmat}, \code{Phi} (\eqn{N_D \times r}),
#'                  and \code{centroid_distmat} (\eqn{r \times r}).
#' @param phi_vals a vector of gridded values of the range parameter psi.
#' @param n_sigmasq an integer to create a grid of sigmasq.
#'
#' @return A list with \code{ll_vals}, \code{sigmasq_vals}, \code{phi_init}, \code{sigmasq_init}.
#'
#' @export
#'
findinit_new <- function(xxb, mod_param, psi_vals, n_sigmasq) {
  
  
  # aggmat: aggregate matrix W_B
  PhiB <- mod_param$aggmat %*% mod_param$Phi
  
  
  
  # centroid_distmat: distance matrix b/w BF's
  centroid_distmat <- mod_param$centroid_distmat
  
  data <- list(resp = xxb, covar = as.matrix(rep(1, length(xxb))),
               PhiB = PhiB, centroid_distmat = centroid_distmat)
  
  ll_vals <- array(NA, dim = c(length(psi_vals), n_sigmasq))
  sigmasq_vals <- array(NA, dim = c(length(psi_vals), n_sigmasq))
  
  # Calculate the basis-function GPB log likelihood over a grid
  for (k in seq_along(psi_vals)) {
    
    # K_Ψ = exp(-𝄁s_i-s_i'𝄁/ Ψ); s_i is centroid
    Kpsi <- exp(-centroid_distmat / psi_vals[k])
    
    # K_B = ɸ_B * K_Ψ * ɸ_B^T = W_B * ɸ_S * K_Ψ * ɸ_S^T * W_B^T
    coar_covmat <- PhiB %*% Kpsi %*% t(PhiB)
    
    
    # _lb, _ub: lower and upper bounds of σ^2 values to try
    # sample variance: var(X_B) / [K_B(Ψ)]
    ratio <- diag(coar_covmat)[1]
    sigmasq_lb <- var(xxb) / ratio
    sigmasq_ub <- var(xxb) / ratio * 3
    sigmasq_seq <- seq(sigmasq_lb, sigmasq_ub, length = n_sigmasq)
    sigmasq_vals[k, ] <- sigmasq_seq
    
    # inner loop: one pass per σ^2_{k, l} value to find log likelihood 
    # under (ψ⁽⁰⁾, σ²⁽⁰⁾) pair
    for (l in seq_along(sigmasq_seq)) {
      param <- c(psi_vals[k], sigmasq_seq[l], mean(xxb))
      ll_vals[k, l] <- gpbauLoglik_new(param, data)
    }
  }
  
  # find the (ψ⁽⁰⁾, σ²⁽⁰⁾) pair in ll_vals with largest log-likelihood
  idx <- as.numeric(which(ll_vals == max(ll_vals), arr.ind = TRUE))
  
  # output entire matrix of candidate value pairs (ψ⁽⁰⁾, σ²⁽⁰⁾) as well
  # as pair with highest log-likelihood
  list(ll_vals = ll_vals,
       sigmasq_vals = sigmasq_vals,
       psi_init = psi_vals[idx[1]],
       sigmasq_init = sigmasq_vals[idx[1], idx[2]])
  
}

# log likelihood under the basis-function GPB model
# C_B = Phi_B * (sigmasq * exp(-d_centroid / psi)) * Phi_B'
gpbauLoglik_new <- function(param, data) {
  
  psi <- param[1]
  sigmasq <- param[2]
  be <- as.matrix(param[-(1:2)])
  
  yy <- data$resp
  rr <- data$covar
  PhiB <- data$PhiB
  
  Kpsi <- sigmasq * exp(-data$centroid_distmat / psi)   # r x r, over basis centroids
  coar_covmat <- PhiB %*% Kpsi %*% t(PhiB)              # Md x Md
  ee <- yy - rr %*% be
  
  logdet <- as.numeric(determinant(coar_covmat)[[1]])
  ll <- -.5 * t(ee) %*% solve(coar_covmat, ee) - .5 * logdet
  
  ll
  
}



# NEW FUNCTIONS WITH VARAIBLE CHANGE

#' @title Find initial values of basis-function GPB parameters
#'
#' @description Basis-function version of \code{\link{findinit}}. Identical grid search,
#'              but the coarse covariance is built as
#'              \eqn{K_B(\psi) = \Phi_B K_\psi \Phi_B^\top} with \eqn{\Phi_B = W_B \phi(S)^\top},
#'              instead of \eqn{W_B \exp(-d/\psi) W_B^\top} over BAU distances.
#'
#' @param xxb a vector of coarse-resolution data of length \eqn{M_D}.
#' @param mod_param a list requiring \code{aggmat}, \code{Phi} (\eqn{N_D \times r}),
#'                  and \code{centroid_distmat} (\eqn{r \times r}).
#' @param psi_vals a vector of gridded values of the range parameter psi.
#' @param n_sigmasq an integer to create a grid of sigmasq.
#'
#' @return A list with \code{ll_vals}, \code{sigmasq_vals}, \code{psi_init}, \code{sigmasq_init}.
#'
#' @export
#'
findinit_new_VARCHANGE <- function(xxb, mod_param, psi_vals, n_sigmasq) {
  
  
  # aggmat: aggregate matrix W_B
  PhiB <- mod_param$aggmat %*% mod_param$Phi
  
  
  
  # centroid_distmat: distance matrix b/w BF's
  centroid_distmat <- mod_param$centroid_distmat
  
  data <- list(resp = xxb, covar = as.matrix(rep(1, length(xxb))),
               PhiB = PhiB, centroid_distmat = centroid_distmat)
  
  ll_vals <- array(NA, dim = c(length(psi_vals), n_sigmasq))
  sigmasq_vals <- array(NA, dim = c(length(psi_vals), n_sigmasq))
  
  # Calculate the basis-function GPB log likelihood over a grid
  for (k in seq_along(psi_vals)) {
    
    # K_Ψ = exp(-𝄁c_i-c_i'𝄁/ Ψ); c_i is centroid
    Kpsi <- exp(-centroid_distmat / psi_vals[k])
    
    # K_B = ɸ_B * K_Ψ * ɸ_B^T = W_B * ɸ_S * K_Ψ * ɸ_S^T * W_B^T
    coar_covmat <- PhiB %*% Kpsi %*% t(PhiB)
    
    
    # _lb, _ub: lower and upper bounds of σ^2 values to try
    # sample variance: var(X_B) / |K_B(Ψ)|
    ratio <- diag(coar_covmat)[1]
    sigmasq_lb <- var(xxb) / ratio
    sigmasq_ub <- var(xxb) / ratio * 3
    sigmasq_seq <- seq(sigmasq_lb, sigmasq_ub, length = n_sigmasq)
    sigmasq_vals[k, ] <- sigmasq_seq
    
    # inner loop: one pass per σ^2_{k, l} value to find log likelihood 
    # under (ψ⁽⁰⁾, σ²⁽⁰⁾) pair
    for (l in seq_along(sigmasq_seq)) {
      param <- c(psi_vals[k], sigmasq_seq[l], mean(xxb))
      ll_vals[k, l] <- gpbauLoglik_new(param, data)
    }
  }
  
  # find the (ψ⁽⁰⁾, σ²⁽⁰⁾) pair in ll_vals with largest log-likelihood
  idx <- as.numeric(which(ll_vals == max(ll_vals), arr.ind = TRUE))
  
  # output entire matrix of candidate value pairs (ψ⁽⁰⁾, σ²⁽⁰⁾) as well
  # as pair with highest log-likelihood
  list(ll_vals = ll_vals,
       sigmasq_vals = sigmasq_vals,
       psi_init = psi_vals[idx[1]],
       sigmasq_init = sigmasq_vals[idx[1], idx[2]])
  
}

# log likelihood under the basis-function GPB model
# K_B = Phi_B * (sigmasq * exp(-d_centroid / psi)) * Phi_B'
gpbauLoglik_new_VARCHANGE <- function(param, data) {
  
  psi <- param[1]
  sigmasq <- param[2]
  be <- as.matrix(param[-(1:2)])
  
  yy <- data$resp
  rr <- data$covar
  PhiB <- data$PhiB
  
  Kpsi <- sigmasq * exp(-data$centroid_distmat / psi)   # r x r, over basis centroids
  coar_covmat <- PhiB %*% Kpsi %*% t(PhiB)              # Md x Md
  ee <- yy - rr %*% be
  
  logdet <- as.numeric(determinant(coar_covmat)[[1]])
  ll <- -.5 * t(ee) %*% solve(coar_covmat, ee) - .5 * logdet
  
  ll
  
}