#include <iostream>
#include <chrono>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <boost/math/special_functions/bessel.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;


// ============================================================================
// fitGPBAU : fit the GPB (Gaussian Process over Basic Areal Units) model by MCMC
// ----------------------------------------------------------------------------
// MODEL (current, Eq. 1 of the MEE paper):
//   Fine BAU field :  x = X_fine * beta + xi,   xi ~ N(0, sigmasq * C(phi))
//   Coarse data    :  y = W * x                 (W = aggmat, averages BAUs to blocks)
//   => y ~ N(X*beta, sigmasq * Sigma(phi)),  Sigma(phi) = W C(phi) W'
//   where C(phi) is the exponential correlation exp(-d/phi) over BAU distances.
//
// The sampler MARGINALISES beta and sigmasq analytically (collapsed likelihood),
// so the only parameter explored by Metropolis-Hastings is the range phi.
// beta and sigmasq are then drawn from their closed-form conditionals (Gibbs).
//
// Priors:  phi ~ Uniform(u_phi, v_phi);  sigmasq ~ InvGamma(u_sigmasq, v_sigmasq);
//          beta ~ flat (improper).
// ============================================================================
// [[Rcpp::export]]
List fitGPBAU(const arma::colvec& y,        // coarse-resolution response, length M_D (= nn)
              const arma::mat& X,           // coarse covariate/design matrix, M_D x p
              const List& prior,            // hyperparameters of the priors (see below)
              const List& starting,         // starting values for phi, sigmasq, beta
              const List& tuning,           // Metropolis proposal spread (se_phi)
              const List& mod_param,         // fixed model inputs: bau_distmat, aggmat
              const List& mcmc_settings) {  // niter, nreport
  
  // ---- Prior hyperparameters -------------------------------------------------
  double u_phi = prior["u_phi"];            // lower bound of Uniform prior on range phi
  double v_phi = prior["v_phi"];            // upper bound of Uniform prior on range phi
  double u_sigmasq = prior["u_sigmasq"];    // shape of InvGamma prior on variance sigmasq
  double v_sigmasq = prior["v_sigmasq"];    // scale of InvGamma prior on variance sigmasq
  
  // ---- Fixed model inputs & MCMC controls ------------------------------------
  arma::mat bau_distmat = mod_param["bau_distmat"]; // N_D x N_D pairwise distances between fine BAUs
  arma::mat aggmat = mod_param["aggmat"];           // W: M_D x N_D averaging matrix (fine -> coarse)
  int niter = mcmc_settings["niter"];               // total number of MCMC iterations to run
  int nreport = mcmc_settings["nreport"];           // print a progress report every nreport iters
  
  // ---- Starting values (current state of the chain) --------------------------
  arma::colvec bb = starting["regcoef"];    // beta: regression coefficients, length p
  double sigmasq = starting["sigmasq"];     // sigmasq: partial sill / marginal variance of the GP
  double phi = starting["phi"];             // phi: range parameter of exp correlation
  double logit_phi = logit(phi, u_phi, v_phi); // phi mapped to the whole real line (so the RW proposal
  // stays inside (u_phi, v_phi) automatically)
  
  // ---- Tuning ----------------------------------------------------------------
  double se_phi = tuning["se_phi"];         // std dev of the Gaussian random-walk proposal for logit_phi
  
  // ---- Dimensions ------------------------------------------------------------
  int nn = y.n_rows;                        // M_D: number of coarse-resolution observations
  int pp = X.n_cols;                        // p: number of covariates (columns of X)
  
  // ---- Scratch scalars -------------------------------------------------------
  // "cur_" = quantities at the CURRENT accepted phi; "cand_" = at the PROPOSED (candidate) phi.
  double blk_log_det, cur_loglik, cur_logpost, quad_log_det;
  double cand_logit_phi, cand_phi, cand_blk_log_det, cand_quad_log_det, cand_loglik, cand_logpost;
  double cand_aa, aa;
  //   blk_log_det   = log|Sigma|            (log-determinant of the coarse covariance)
  //   quad_log_det  = log|X' Sigma^{-1} X|  (log-det of the "information" matrix for beta)
  //   aa            = e' Sigma^{-1} e       (generalised residual sum of squares)
  //   cur_loglik/cur_logpost = collapsed log-likelihood / log-posterior at current phi
  
  // ---- Scratch vectors/matrices ----------------------------------------------
  arma::colvec bb_hat, proj_y, ee, cand_bb_hat, cand_proj_y, cand_ee;
  //   bb_hat  = GLS estimate of beta:  (X' Sigma^{-1} X)^{-1} X' Sigma^{-1} y
  //   proj_y  = X' Sigma^{-1} y        (the "projected" response)
  //   ee      = y - X*bb_hat           (GLS residuals)
  arma::mat bau_corrmat, blk_corrmat, dat_corrmat_inv, quad_corrmat;
  arma::mat cand_bau_corrmat, cand_blk_corrmat, cand_dat_corrmat_inv, cand_quad_corrmat;
  arma::mat bb_corr;
  //   bau_corrmat  = C(phi)         : N_D x N_D correlation over fine BAUs (no sigmasq)
  //   blk_corrmat  = Sigma = W C W' : M_D x M_D correlation over coarse blocks
  //   quad_corrmat = X' Sigma^{-1} X: p x p, the GLS information matrix for beta
  //   bb_corr      = (X' Sigma^{-1} X)^{-1}: posterior covariance of beta (up to sigmasq)
  //   (dat_corrmat_inv is declared but unused.)
  
  // ---- Evaluate the collapsed likelihood at the STARTING phi -----------------
  bau_corrmat = expCorr(bau_distmat, phi);                 // C(phi) over BAUs
  blk_corrmat = aggmat * bau_corrmat * aggmat.t();         // Sigma = W C W' over blocks
  blk_log_det = arma::log_det_sympd(blk_corrmat);          // log|Sigma|  (sympd = symmetric pos-def)
  quad_corrmat = X.t() * arma::solve(blk_corrmat, X, arma::solve_opts::likely_sympd); // X' Sigma^{-1} X
  quad_log_det = arma::log_det_sympd(quad_corrmat);        // log|X' Sigma^{-1} X|
  
  proj_y = X.t() * arma::solve(blk_corrmat, y, arma::solve_opts::likely_sympd);       // X' Sigma^{-1} y
  bb_hat = arma::solve(quad_corrmat, proj_y);              // GLS estimate beta_hat
  ee = y - X * bb_hat;                                     // residuals
  aa = as_scalar(ee.t() * arma::solve(blk_corrmat, ee, arma::solve_opts::likely_sympd)); // e' Sigma^{-1} e
  
  // Collapsed log-likelihood of phi (beta & sigmasq integrated out):
  //   -1/2 log|Sigma| - 1/2 log|X'Sigma^{-1}X| - (u_sigmasq + (n-p)/2) * log(v_sigmasq + aa/2)
  cur_loglik = -.5 * blk_log_det - .5* quad_log_det;
  cur_loglik += -(.5 * (nn - pp) + u_sigmasq) * log(.5 * aa + v_sigmasq);
  // Add the log-prior for phi ON THE LOGIT SCALE. The two log terms are the Jacobian of the
  // logit transform, log((phi-u_phi)(v_phi-phi)), which turns the Uniform(u_phi,v_phi) prior
  // into the correct target density for the unconstrained random walk.
  cur_logpost = cur_loglik + log(phi - u_phi) + log(v_phi - phi);
  
  // ---- Storage for the posterior draws ---------------------------------------
  arma::colvec phi_save(niter), sigmasq_save(niter), mh_data(niter); // mh_data = 1 if accepted, else 0
  arma::mat bb_save(pp, niter);                                      // one column of beta per iteration
  
  // ---- Bookkeeping / timing counters -----------------------------------------
  int icount = 0;          // iterations since last progress report
  int batch_accept = 0;    // acceptances since last report (for interval acceptance rate)
  int accept = 0;          // total acceptances so far (for overall acceptance rate)
  double block_time = 0.0; // wall-clock time accumulated since last report
  double time_diff, ert;   // per-iter elapsed time; ert = estimated remaining time
  
  // ===================== MCMC loop ===========================================
  for (int iter = 0; iter < niter; ++iter) {
    
    auto start = std::chrono::steady_clock::now();  // start timer for this iteration
    
    // ---- (1) Metropolis-Hastings update of phi -------------------------------
    // Propose a new phi by a Gaussian random walk on the logit scale, then map back.
    cand_logit_phi= R::rnorm(logit_phi, se_phi);           // proposal on unconstrained scale
    cand_phi = logitInv(cand_logit_phi, u_phi, v_phi);     // back to (u_phi, v_phi)
    
    // Recompute the collapsed likelihood at the candidate phi (same steps as above).
    cand_bau_corrmat = expCorr(bau_distmat, cand_phi);
    cand_blk_corrmat = aggmat * cand_bau_corrmat * aggmat.t();
    cand_blk_log_det = arma::log_det_sympd(cand_blk_corrmat);
    cand_quad_corrmat = X.t() * arma::solve(cand_blk_corrmat, X, arma::solve_opts::likely_sympd);
    cand_quad_log_det = arma::log_det_sympd(cand_quad_corrmat);
    
    cand_proj_y = X.t() * arma::solve(cand_blk_corrmat, y, arma::solve_opts::likely_sympd);
    cand_bb_hat = arma::solve(cand_quad_corrmat, cand_proj_y);
    cand_ee = y - X * cand_bb_hat;
    // NOTE: uses `ee` (current residuals) rather than `cand_ee`. Looks like a bug — the candidate
    // quadratic form should be cand_ee' * cand_Sigma^{-1} * cand_ee. Flagging, not changing.
    cand_aa = as_scalar(ee.t() * arma::solve(cand_blk_corrmat, ee, arma::solve_opts::likely_sympd));
    
    cand_loglik = -.5 * cand_blk_log_det - .5* cand_quad_log_det;
    cand_loglik += -(.5 * (nn - pp) + u_sigmasq) * log(.5 * cand_aa + v_sigmasq);
    cand_logpost = cand_loglik + log(cand_phi - u_phi) + log(v_phi - cand_phi);
    
    // Accept with probability min(1, exp(cand_logpost - cur_logpost)).
    // Comparing log(U) <= (cand - cur) is the standard log-scale acceptance test.
    if (log(arma::randu()) <= (cand_logpost - cur_logpost)) {
      
      // Accepted: promote all candidate quantities to be the new current state.
      phi = cand_phi;
      logit_phi = cand_logit_phi;
      blk_log_det = cand_blk_log_det;
      quad_log_det = cand_quad_log_det;
      cur_logpost = cand_logpost;
      aa = cand_aa;
      batch_accept++;
      accept++;
      mh_data(iter) = 1;                                   // record acceptance
      
      // ---- (2) Gibbs update of sigmasq (only when phi moved) -----------------
      // Full conditional: sigmasq ~ InvGamma(u_sigmasq + (n-p)/2, v_sigmasq + aa/2).
      // Drawn as 1 / Gamma(shape, rate); R::rgamma takes (shape, SCALE)=1/rate.
      sigmasq = 1.0 / R::rgamma(u_sigmasq + .5 * (nn - pp), 1.0 / (v_sigmasq + .5 *aa));
      
      // ---- (3) Gibbs update of beta (only when phi moved) --------------------
      // Full conditional: beta ~ N( beta_hat, sigmasq * (X' Sigma^{-1} X)^{-1} ).
      bb_corr = arma::inv_sympd(cand_quad_corrmat);        // (X' Sigma^{-1} X)^{-1}
      bb = arma::mvnrnd(cand_bb_hat, sigmasq* bb_corr);    // one multivariate-normal draw
      
    }
    // If rejected: phi, sigmasq, beta all keep their previous values (carried forward below).
    
    auto end = std::chrono::steady_clock::now();
    time_diff = std::chrono::duration<double>(end - start).count(); // seconds for this iter
    block_time += time_diff;
    
    // ---- Store this iteration's draws ----------------------------------------
    phi_save(iter) = phi;
    sigmasq_save(iter) = sigmasq;
    bb_save.col(iter) = bb;
    
    // ---- Periodic progress report --------------------------------------------
    icount++;
    if (icount == nreport) {
      Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
      ert = (block_time / nreport) * (niter - iter);       // rough estimated remaining time
      Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
      Rprintf("Interval acceptance rate: %3.2f%%\n", 100.0 * batch_accept / nreport);
      Rprintf("Overall acceptance rate: %3.2f%%\n\n", 100.0 * accept / (iter + 1));
      icount = 0;
      batch_accept = 0;
      block_time = 0.0;
    }
    
  }
  // =================== end MCMC loop =========================================
  
  // Bundle the posterior draws (keys "phi","sigmasq","beta" are read later by predGPBAU2/3).
  List post_sams = List::create(Named("phi") = phi_save,
                                Named("sigmasq") = sigmasq_save,
                                Named("beta") = bb_save);
  
  return List::create(Named("post_sams") = post_sams,   // the chains
                      Named("mh_data") = mh_data);        // per-iter accept/reject flags
}



// GPBAU using new basis functions method
// [[Rcpp::export]]
List fitGPBAU_new(const arma::colvec& xB,       // coarse-resolution response X_B (paper: X_B), length M_D (= nn)
              const arma::mat& RB,          // coarse covariate/design matrix R_B (paper: R_B), M_D x (p+1)
              const List& prior,            // hyperparameters of the priors (see below)
              const List& starting,         // starting values for psi, sigmasq, gam
              const List& tuning,           // Metropolis proposal spread (se_psi)
              const List& mod_param,         // fixed model inputs: bau_distmat, aggmat
              const List& mcmc_settings) {  // niter, nreport

  // ---- Prior hyperparameters -------------------------------------------------
  // (list KEYS stay "u_phi"/"v_phi" so the R wrapper is unchanged; locals renamed to psi.)
  double u_psi = prior["u_phi"];            // lower bound of Uniform prior on range psi
  double v_psi = prior["v_phi"];            // upper bound of Uniform prior on range psi
  double u_sigmasq = prior["u_sigmasq"];    // shape of InvGamma prior on variance sigmasq
  double v_sigmasq = prior["v_sigmasq"];    // scale of InvGamma prior on variance sigmasq

  // ---- Fixed model inputs & MCMC controls ------------------------------------
  // arma::mat bau_distmat = mod_param["bau_distmat"]; // N_D x N_D pairwise distances between fine BAUs
  arma::mat aggmat = mod_param["aggmat"];           // W_B: M_D x N_D averaging matrix (fine -> coarse)
  int niter = mcmc_settings["niter"];               // total number of MCMC iterations to run
  int nreport = mcmc_settings["nreport"];           // print a progress report every nreport iters

  // ---- Starting values (current state of the chain) --------------------------
  arma::colvec gam = starting["regcoef"];   // gam = gamma-tilde: (intercept gamma0, covariate coefs gamma), length p+1
  double sigmasq = starting["sigmasq"];     // sigmasq: variance parameter sigma^2
  double psi = starting["phi"];             // psi: range parameter of exp correlation (paper's psi)
  double logit_psi = logit(psi, u_psi, v_psi); // psi mapped to the whole real line (so the RW proposal
  // stays inside (u_psi, v_psi) automatically)

  // ---- Tuning ----------------------------------------------------------------
  double se_psi = tuning["se_phi"];         // std dev of the Gaussian random-walk proposal for logit_psi

  // ---- Dimensions ------------------------------------------------------------
  int nn = xB.n_rows;                       // M_D: number of coarse-resolution observations
  int pp = RB.n_cols;                       // p+1: number of columns of R_B (covariates incl. intercept)

  // ---- Scratch scalars -------------------------------------------------------
  // "cur_" = quantities at the CURRENT accepted psi; "cand_" = at the PROPOSED (candidate) psi.
  double KB_log_det, cur_loglik, cur_logpost, quad_log_det;
  double cand_logit_psi, cand_psi, cand_KB_log_det, cand_quad_log_det, cand_loglik, cand_logpost;
  double cand_S2, S2;
  //   KB_log_det    = log|K_B|              (log-determinant of the coarse correlation)
  //   quad_log_det  = log|R_B' K_B^{-1} R_B|(log-det of the "information" matrix for gam)
  //   S2            = e' K_B^{-1} e         (generalised residual sum of squares, paper's S^2(psi))
  //   cur_loglik/cur_logpost = collapsed log-likelihood / log-posterior at current psi

  // ---- Scratch vectors/matrices ----------------------------------------------
  arma::colvec gam_hat, proj_xB, ee, cand_gam_hat, cand_proj_xB, cand_ee;
  //   gam_hat = GLS estimate of gam:  (R_B' K_B^{-1} R_B)^{-1} R_B' K_B^{-1} xB
  //   proj_xB = R_B' K_B^{-1} xB      (the "projected" response)
  //   ee      = xB - RB*gam_hat       (GLS residuals)
  arma::mat Kpsi, KB, dat_corrmat_inv, quad_corrmat;
  arma::mat cand_Kpsi, cand_KB, cand_dat_corrmat_inv, cand_quad_corrmat;
  arma::mat gam_corr;
  //   Kpsi   = K_psi = expCorr(...)     : correlation from exp(-d/psi), no sigmasq
  //   KB     = K_B = W_B K_psi W_B'     : M_D x M_D correlation over coarse blocks
  //   quad_corrmat = R_B' K_B^{-1} R_B  : (p+1)x(p+1), the GLS information matrix for gam
  //   gam_corr = (R_B' K_B^{-1} R_B)^{-1}: posterior covariance of gam (up to sigmasq)
  //   (dat_corrmat_inv is declared but unused.)
  
  
  // NEW CHANGES START HERE
  
  // IMPORTANT: specify new input variable for phi (basis function models)
  arma::mat Phi = mod_param["Phi"];                             // φ(S)ᵀ : N_D x r, basis fns at BAU centroids
  arma::mat centroid_distmat = mod_param["centroid_distmat"];   // r x r pairwise distance matrix among centroids
  arma::mat PsiB = aggmat * Phi;                               // Ψ_B = W_B φ(S)ᵀ : M_D x r.
                                                                // Used to simplify calculation of some posterior hyperparameters
  

  // ---- Evaluate the collapsed likelihood at the STARTING psi -----------------
  Kpsi = expCorr(centroid_distmat, psi);                        // K_psi over centroids
  KB = PsiB * Kpsi * PsiB.t();                         // K_B = PsiB K_psi PsiB' over blocks
  KB_log_det = arma::log_det_sympd(KB);                    // log|K_B|  (sympd = symmetric pos-def)
  quad_corrmat = RB.t() * arma::solve(KB, RB, arma::solve_opts::likely_sympd); // R_B' K_B^{-1} R_B
  quad_log_det = arma::log_det_sympd(quad_corrmat);        // log|R_B' K_B^{-1} R_B|

  // fitting regression with generalized least squares
  // code assumes flat prior on V, eg. V_mu^-1 = 0, so implementation is slightly different from og derivation
  // to fully stick to derivation, must change quad_corrmat
  
  // building mu_hat_gamma_hat (aka gam_hat)
  proj_xB = RB.t() * arma::solve(KB, xB, arma::solve_opts::likely_sympd);       // R_B' K_B^{-1} xB
  gam_hat = arma::solve(quad_corrmat, proj_xB);            // GLS estimate gam_hat
  
  // building S2
  ee = xB - RB * gam_hat;                                  // residuals
  S2 = as_scalar(ee.t() * arma::solve(KB, ee, arma::solve_opts::likely_sympd)); // e' K_B^{-1} e
  

  // Collapsed log-likelihood of psi (gam & sigmasq integrated out):
  //   -1/2 log|K_B| - 1/2 log|R_B'K_B^{-1}R_B| - (u_sigmasq + (n-p)/2) * log(v_sigmasq + S2/2)
  cur_loglik = -.5 * KB_log_det - .5* quad_log_det;
  cur_loglik += -(.5 * (nn - pp) + u_sigmasq) * log(.5 * S2 + v_sigmasq);
  // Add the log-prior for psi ON THE LOGIT SCALE. The two log terms are the Jacobian of the
  // logit transform, log((psi-u_psi)(v_psi-psi)), which turns the Uniform(u_psi,v_psi) prior
  // into the correct target density for the unconstrained random walk.
  cur_logpost = cur_loglik + log(psi - u_psi) + log(v_psi - psi);

  // Fail loudly rather than freezing. A non-finite starting log-posterior
  // (e.g. psi exactly on a prior bound -> log(0) = -Inf) makes every
  // acceptance ratio NaN, so the chain never moves and sigmasq/gam stay
  // pinned at their starting values for all niter iterations -- which looks
  // like a converged chain with a zero-width credible interval.
  if (!std::isfinite(cur_logpost)) {
    Rcpp::stop("fitGPBAU_new: initial log-posterior is not finite (psi may be on a prior bound).");
  }


  // NEW CHANGES END HERE 7/22/26 9PM


  // ---- Storage for the posterior draws ---------------------------------------
  arma::colvec psi_save(niter), sigmasq_save(niter), mh_data(niter); // mh_data = 1 if accepted, else 0
  arma::mat gam_save(pp, niter);                                     // one column of gam per iteration

  // ---- Bookkeeping / timing counters -----------------------------------------
  int icount = 0;          // iterations since last progress report
  int batch_accept = 0;    // acceptances since last report (for interval acceptance rate)
  int accept = 0;          // total acceptances so far (for overall acceptance rate)
  double block_time = 0.0; // wall-clock time accumulated since last report
  double time_diff, ert;   // per-iter elapsed time; ert = estimated remaining time

  // ===================== MCMC loop ===========================================
  for (int iter = 0; iter < niter; ++iter) {

    auto start = std::chrono::steady_clock::now();  // start timer for this iteration

    // ---- (1) Metropolis-Hastings update of psi -------------------------------
    // Propose a new psi by a Gaussian random walk on the logit scale, then map back.
    cand_logit_psi= R::rnorm(logit_psi, se_psi);           // proposal on unconstrained scale
    cand_psi = logitInv(cand_logit_psi, u_psi, v_psi);     // back to (u_psi, v_psi)

    // Recompute the collapsed likelihood at the candidate psi (same steps as above).
    cand_Kpsi = expCorr(centroid_distmat, cand_psi);
    cand_KB = PsiB * cand_Kpsi * PsiB.t();
    cand_KB_log_det = arma::log_det_sympd(cand_KB);
    cand_quad_corrmat = RB.t() * arma::solve(cand_KB, RB, arma::solve_opts::likely_sympd);
    cand_quad_log_det = arma::log_det_sympd(cand_quad_corrmat);
    

    cand_proj_xB = RB.t() * arma::solve(cand_KB, xB, arma::solve_opts::likely_sympd);
    cand_gam_hat = arma::solve(cand_quad_corrmat, cand_proj_xB);
    cand_ee = xB - RB * cand_gam_hat;
    // NOTE: changed 'ee' to 'cand_ee'
    // cand_S2 = as_scalar(ee.t() * arma::solve(cand_KB, ee, arma::solve_opts::likely_sympd));
    cand_S2 = as_scalar(cand_ee.t() * arma::solve(cand_KB, cand_ee, arma::solve_opts::likely_sympd));

    cand_loglik = -.5 * cand_KB_log_det - .5* cand_quad_log_det;
    cand_loglik += -(.5 * (nn - pp) + u_sigmasq) * log(.5 * cand_S2 + v_sigmasq);
    cand_logpost = cand_loglik + log(cand_psi - u_psi) + log(v_psi - cand_psi);

    // Accept with probability min(1, exp(cand_logpost - cur_logpost)).
    // Comparing log(U) <= (cand - cur) is the standard log-scale acceptance test.
    if (log(arma::randu()) <= (cand_logpost - cur_logpost)) {

      // Accepted: promote all candidate quantities to be the new current state.
      psi = cand_psi;
      logit_psi = cand_logit_psi;
      KB_log_det = cand_KB_log_det;
      quad_log_det = cand_quad_log_det;
      cur_logpost = cand_logpost;
      S2 = cand_S2;
      batch_accept++;
      accept++;
      mh_data(iter) = 1;                                   // record acceptance

      // ---- (2) Gibbs update of sigmasq (only when psi moved) -----------------
      // Full conditional: sigmasq ~ InvGamma(u_sigmasq + (n-p)/2, v_sigmasq + S2/2).
      // Drawn as 1 / Gamma(shape, rate); R::rgamma takes (shape, SCALE)=1/rate.
      sigmasq = 1.0 / R::rgamma(u_sigmasq + .5 * (nn - pp), 1.0 / (v_sigmasq + .5 *S2));

      // ---- (3) Gibbs update of gam (only when psi moved) ---------------------
      // Full conditional: gam ~ N( gam_hat, sigmasq * (R_B' K_B^{-1} R_B)^{-1} ).
      gam_corr = arma::inv_sympd(cand_quad_corrmat);       // (R_B' K_B^{-1} R_B)^{-1}
      gam = arma::mvnrnd(cand_gam_hat, sigmasq* gam_corr); // one multivariate-normal draw

    }
    // If rejected: psi, sigmasq, gam all keep their previous values (carried forward below).

    auto end = std::chrono::steady_clock::now();
    time_diff = std::chrono::duration<double>(end - start).count(); // seconds for this iter
    block_time += time_diff;

    // ---- Store this iteration's draws ----------------------------------------
    psi_save(iter) = psi;
    sigmasq_save(iter) = sigmasq;
    gam_save.col(iter) = gam;

    // ---- Periodic progress report --------------------------------------------
    icount++;
    if (icount == nreport) {
      Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
      ert = (block_time / nreport) * (niter - iter);       // rough estimated remaining time
      Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
      Rprintf("Interval acceptance rate: %3.2f%%\n", 100.0 * batch_accept / nreport);
      Rprintf("Overall acceptance rate: %3.2f%%\n\n", 100.0 * accept / (iter + 1));
      icount = 0;
      batch_accept = 0;
      block_time = 0.0;
    }

  }
  // =================== end MCMC loop =========================================

  // Bundle the posterior draws. Output KEYS stay "phi","sigmasq","beta" so the R wrapper and
  // predGPBAU2/3 keep working unchanged (only the C++ variable names changed to psi_save/gam_save).
  List post_sams = List::create(Named("phi") = psi_save,
                                Named("sigmasq") = sigmasq_save,
                                Named("beta") = gam_save);

  return List::create(Named("post_sams") = post_sams,   // the chains
                      Named("mh_data") = mh_data);        // per-iter accept/reject flags
}




// BAUGP (joint) prediction
// [[Rcpp::export]]
arma::mat predGPBAU2(const arma::colvec& y,
                     const arma::mat& X,
                     const List& post_sams,
                     const List& mod_param,
                     const int& nreport){
   
   arma::colvec phi_save = post_sams["phi"];
   arma::colvec sigmasq_save = post_sams["sigmasq"];
   arma::mat bb_save = post_sams["beta"];
   
   arma::mat bau_distmat = mod_param["bau_distmat"];
   arma::mat aggmat = mod_param["aggmat"];
   
   double phi, sigmasq;
   int niter = phi_save.n_rows;
   int nn = bau_distmat.n_rows;
   
   arma::mat bau_corrmat, blk_corrmat, bau_blk_corr, cond_var;
   arma::colvec cond_mean, bb, fine_mean, coarse_mean;
   arma::mat latent_save(nn, niter);
   arma::mat blk_corrmat_inv, bau_blk_coef;
   
   int icount = 0;
   double block_time = 0.0;
   double time_diff, ert;
   
   for (int iter = 0; iter < niter; ++iter){
     
     auto start = std::chrono::steady_clock::now();
     
     phi = phi_save(iter);
     sigmasq = sigmasq_save(iter);
     bb = bb_save.row(iter).t();
     
     fine_mean = X * bb;
     coarse_mean = aggmat * fine_mean;
     
     bau_corrmat = sigmasq * expCorr(bau_distmat, phi);
     bau_blk_corr = bau_corrmat * aggmat.t();
     blk_corrmat = aggmat * bau_blk_corr;
     
     // Prediction
     cond_mean = fine_mean + bau_blk_corr * solve(blk_corrmat, (y - coarse_mean), arma::solve_opts::likely_sympd);
     cond_var = bau_corrmat - bau_blk_corr * solve(blk_corrmat, bau_blk_corr.t(), arma::solve_opts::likely_sympd);
     latent_save.col(iter) = arma::mvnrnd(cond_mean, cond_var);
     
     auto end = std::chrono::steady_clock::now();
     time_diff = std::chrono::duration<double>(end - start).count();
     block_time += time_diff;
     
     icount++;
     if (icount == nreport) {
       Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
       ert = (block_time / nreport) * (niter - iter);
       Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
       icount = 0;
       block_time = 0.0;
     }
     
   }
   
   return latent_save;
   
 }


// BAUGP (independent) prediction
// [[Rcpp::export]]
arma::mat predGPBAU3(const arma::colvec& y,
                      const arma::mat& X,
                      const List& post_sams,
                      const List& mod_param,
                      const int& nreport){

   arma::colvec phi_save = post_sams["phi"];
   arma::colvec sigmasq_save = post_sams["sigmasq"];
   arma::mat bb_save = post_sams["beta"];

   arma::mat bau_distmat = mod_param["bau_distmat"];
   arma::mat aggmat = mod_param["aggmat"];

   double phi, sigmasq;
   int niter = phi_save.n_rows;
   int nn = bau_distmat.n_rows;

   arma::mat bau_corrmat, blk_corrmat, bau_blk_corr;
   arma::colvec bb, fine_mean, coarse_mean;
   arma::mat latent_save(nn, niter);
   double cond_mean, cond_var;
   arma::mat blk_corrmat_inv;
   arma::rowvec bau_blk_coef;
   arma::rowvec ss;


   int icount = 0;
   double block_time = 0.0;
   double time_diff, ert;

   for (int iter = 0; iter < niter; ++iter){

     auto start = std::chrono::steady_clock::now();

     phi = phi_save(iter);
     sigmasq = sigmasq_save(iter);
     bb = bb_save.row(iter).t();

     fine_mean = X * bb;
     coarse_mean = aggmat * fine_mean;

     bau_corrmat = sigmasq * expCorr(bau_distmat, phi);
     bau_blk_corr = bau_corrmat * aggmat.t();
     blk_corrmat = aggmat * bau_blk_corr;

     blk_corrmat_inv = arma::inv_sympd(blk_corrmat);
     
     for (int i = 0; i < nn; ++i) {
       
       ss = bau_blk_corr.row(i);
       bau_blk_coef = ss * blk_corrmat_inv;
       cond_mean = fine_mean(i) + arma::as_scalar(bau_blk_coef * (y - coarse_mean));
       cond_var = sigmasq - arma::as_scalar(bau_blk_coef * ss.t());
       latent_save(i, iter) = R::rnorm(cond_mean, sqrt(cond_var));

     }

     auto end = std::chrono::steady_clock::now();
     time_diff = std::chrono::duration<double>(end - start).count();
     block_time += time_diff;

     icount++;
     if (icount == nreport) {
       Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
       ert = (block_time / nreport) * (niter - iter);
       Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
       icount = 0;
       block_time = 0.0;
     }

   }

   return latent_save;

 }




// ============================================================================
// BAUGP (joint) prediction -- BASIS-FUNCTION version
// ----------------------------------------------------------------------------
// MODEL:  X_s = R~ gamma + Phi alpha,   alpha ~ N(0, sigmasq * K_psi)
//         X_B = W_B X_s = R_B gamma + Psi_B alpha    (EXACT: delta(S) = 0)
//
// WHY THIS IS A RESTRUCTURE, NOT A FORMULA SWAP:
// predGPBAU2 built C_S = sigmasq * expCorr(bau_distmat, phi) as an Nd x Nd
// matrix, then drew from N(cond_mean, cond_var) via mvnrnd. Under the basis
// model C_S = sigmasq * Phi K_psi Phi' is still Nd x Nd but has rank only r,
// so the conditional covariance is SINGULAR and mvnrnd cannot factorise it.
// (Note the ORIGINAL was already mildly degenerate: rank(Sigma_1|2) = Nd - Md
// = 2970 of 3000, which usually squeaked through Cholesky on round-off. Here
// the conditional rank drops to r - Md = 90 of 3000 -- no such luck.)
//
// Instead we condition in ALPHA-space using conditional simulation by kriging
// (Journel 1974; Chiles & Delfiner, Geostatistics sec. 7.3). Because
// delta(S) = 0, observing X_B is a LINEAR CONSTRAINT on the r-dimensional
// alpha, not a noisy observation:
//   1. draw alpha_sim ~ N(0, sigmasq * K_psi)                (r-dim, K_psi is PD)
//   2. alpha_cond = alpha_sim + L (c - Psi_B alpha_sim),
//        L = V Psi_B'(Psi_B V Psi_B')^{-1} = K_psi Psi_B' K_B^{-1}   (sigmasq cancels)
//   3. X_s = R~ gamma + Phi alpha_cond
// alpha_cond is exactly N(Lc, V - L Psi_B V), i.e. the same conditional law
// Algorithm 1 specifies -- but the singular covariance is realised IMPLICITLY
// as the image of a full-rank draw under the projection (I - L Psi_B). The
// only matrices factorised are K_psi (full rank r) and K_B (full rank Md when
// r >= Md), so the rank deficiency never has to be handled directly.
// Cost is O(r^3 + Md^3 + Nd*r) instead of O(Nd^3).
//
// INVARIANT: since Psi_B L = I, every draw satisfies W_B X_s = X_B exactly.
// Use this as a correctness check: max(abs(aggmat %*% x_pred - coar_xx)) ~ 1e-12.
// ============================================================================
// ARGUMENTS:
//   xB      coarse response X_B, length M_D
//   Rtilde  FINE-scale design matrix (N_D x p).
//           NOTE: this is R-tilde, NOT R_B. fitGPBAU_new takes the COARSE
//           design R_B (M_D x p); here the caller passes fine_covar, and the
//           two are related by R_B = W_B * R-tilde. Different objects.
//
// Rcpp::compileAttributes CAVEAT: nothing may sit between the export attribute
// and the function signature -- it scans that whole span as the signature and
// will splice comment text into the generated wrapper. Keep notes ABOVE the
// attribute, and keep the parameter list free of multi-line comments.
// [[Rcpp::export]]
arma::mat predGPBAU2_new(const arma::colvec& xB,
                         const arma::mat& Rtilde,
                         const List& post_sams,
                         const List& mod_param,
                         const int& nreport){

   // unchanged: keys stay "phi"/"sigmasq"/"beta" ("phi" holds psi, "beta" holds gamma-tilde)
   arma::colvec psi_save = post_sams["phi"];
   arma::colvec sigmasq_save = post_sams["sigmasq"];
   arma::mat gam_save = post_sams["beta"];

   // CHANGED: basis inputs replace bau_distmat
   arma::mat Phi = mod_param["Phi"];                           // Nd x r
   arma::mat centroid_distmat = mod_param["centroid_distmat"]; // r x r
   arma::mat aggmat = mod_param["aggmat"];                     // Md x Nd

   // NEW: Psi_B = W_B phi(S)' does not depend on psi -> compute ONCE
   arma::mat PsiB = aggmat * Phi;                              // Md x r

   double psi, sigmasq;
   int niter = psi_save.n_rows;
   int nn = Phi.n_rows;        // CHANGED: Nd now comes from Phi, not bau_distmat
   int rr = Phi.n_cols;        // NEW: r, the basis dimension

   // CHANGED: no bau_corrmat / bau_blk_corr / cond_var. Nothing of size
   // Nd x Nd is ever allocated (the old version needed ~72 MB per draw).
   arma::mat Kpsi, KB;
   arma::colvec gam, fine_mean, coarse_mean, alpha_sim, alpha_cond, resid;
   arma::mat latent_save(nn, niter);

   int icount = 0;
   double block_time = 0.0;
   double time_diff, ert;

   for (int iter = 0; iter < niter; ++iter){

     auto start = std::chrono::steady_clock::now();

     psi = psi_save(iter);
     sigmasq = sigmasq_save(iter);
     gam = gam_save.row(iter).t();

     fine_mean = Rtilde * gam;            // R~ gamma  : Nd x 1  (unchanged)
     coarse_mean = aggmat * fine_mean;    // R_B gamma : Md x 1  (unchanged)

     // CHANGED: r x r kernel over basis centroids, then the Psi_B sandwich.
     // KB is kept on the CORRELATION scale because sigmasq cancels in L.
     Kpsi = expCorr(centroid_distmat, psi);          // r x r
     KB = PsiB * Kpsi * PsiB.t();                    // Md x Md

     // --- conditional simulation in alpha-space (replaces the kriging block) ---
     // 1. unconditional draw of the basis coefficients
     alpha_sim = arma::mvnrnd(arma::zeros<arma::colvec>(rr), sigmasq * Kpsi);

     // 2. discrepancy between this draw's implied blocks and the observed data
     resid = (xB - coarse_mean) - PsiB * alpha_sim;   // Md x 1

     // 3. project onto the constraint Psi_B alpha = X_B - R_B gamma
     alpha_cond = alpha_sim +
       Kpsi * PsiB.t() * arma::solve(KB, resid, arma::solve_opts::likely_sympd);

     // 4. lift the r-dim coefficients up to the Nd fine-resolution locations
     latent_save.col(iter) = fine_mean + Phi * alpha_cond;

     auto end = std::chrono::steady_clock::now();
     time_diff = std::chrono::duration<double>(end - start).count();
     block_time += time_diff;

     icount++;
     if (icount == nreport) {
       Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
       ert = (block_time / nreport) * (niter - iter);
       Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
       icount = 0;
       block_time = 0.0;
     }

   }

   return latent_save;

 }


// ============================================================================
// BAUGP (independent) prediction -- BASIS-FUNCTION version
// ----------------------------------------------------------------------------
// Much lighter edit than the joint version: the per-location loop survives
// intact, because only MARGINAL conditional moments are needed and both
// ingredients come straight from Phi without ever forming C_S:
//     C_S[i,i]   = sigmasq * phi_i' K_psi phi_i        (marginal variance)
//     C_S[i,.]W' = sigmasq * phi_i' K_psi Psi_B'       (cross-covariance row)
// Draws are correct pointwise but carry no cross-location correlation within
// a draw -- use predGPBAU2_new when functionals of the whole field are needed.
// ============================================================================
// ARGUMENTS:
//   xB      coarse response X_B, length M_D
//   Rtilde  FINE-scale design matrix (N_D x p) -- R-tilde, not R_B
// [[Rcpp::export]]
arma::mat predGPBAU3_new(const arma::colvec& xB,
                         const arma::mat& Rtilde,
                         const List& post_sams,
                         const List& mod_param,
                         const int& nreport){

   arma::colvec psi_save = post_sams["phi"];
   arma::colvec sigmasq_save = post_sams["sigmasq"];
   arma::mat gam_save = post_sams["beta"];

   // CHANGED: basis inputs replace bau_distmat
   arma::mat Phi = mod_param["Phi"];
   arma::mat centroid_distmat = mod_param["centroid_distmat"];
   arma::mat aggmat = mod_param["aggmat"];

   arma::mat PsiB = aggmat * Phi;      // NEW: Md x r, psi-independent

   double psi, sigmasq;
   int niter = psi_save.n_rows;
   int nn = Phi.n_rows;                // CHANGED: Nd from Phi

   arma::mat Kpsi, KB, KB_inv, PhiKpsi;
   arma::colvec gam, fine_mean, coarse_mean;
   arma::mat latent_save(nn, niter);
   double cond_mean, cond_var, marg_var;
   arma::rowvec ss, coef;

   int icount = 0;
   double block_time = 0.0;
   double time_diff, ert;

   for (int iter = 0; iter < niter; ++iter){

     auto start = std::chrono::steady_clock::now();

     psi = psi_save(iter);
     sigmasq = sigmasq_save(iter);
     gam = gam_save.row(iter).t();

     fine_mean = Rtilde * gam;
     coarse_mean = aggmat * fine_mean;

     // CHANGED: correlation-scale quantities. sigmasq cancels in the mean
     // correction and is applied explicitly to the variance below.
     Kpsi = expCorr(centroid_distmat, psi);     // r x r
     KB = PsiB * Kpsi * PsiB.t();               // Md x Md
     KB_inv = arma::inv_sympd(KB);

     // NEW: hoist Phi * K_psi out of the location loop -- one Nd x r multiply
     // per iteration instead of rebuilding per location.
     PhiKpsi = Phi * Kpsi;                      // Nd x r

     for (int i = 0; i < nn; ++i) {

       // CHANGED: cross-covariance row built from basis quantities rather
       // than read off a row of an explicit Nd x Nd matrix.
       ss = PhiKpsi.row(i) * PsiB.t();                     // 1 x Md
       coef = ss * KB_inv;                                 // 1 x Md

       cond_mean = fine_mean(i) + arma::as_scalar(coef * (xB - coarse_mean));

       // IMPORTANT DIFFERENCE: the old code used `sigmasq` directly as the
       // marginal variance because expCorr has 1 on its diagonal. Under the
       // basis model the marginal variance is sigmasq * phi_i' K_psi phi_i,
       // which is NOT 1 and varies by location (it depends on how much basis
       // support sits at s_i). Using sigmasq here would give uniform, wrong
       // predictive variances -- visible only as mis-calibrated coverage.
       marg_var = arma::as_scalar(PhiKpsi.row(i) * Phi.row(i).t());
       cond_var = sigmasq * (marg_var - arma::as_scalar(coef * ss.t()));

       // NEW: guard against tiny negative values from round-off; the
       // conditional variance is non-negative in exact arithmetic.
       if (cond_var < 0.0) cond_var = 0.0;

       latent_save(i, iter) = R::rnorm(cond_mean, sqrt(cond_var));

     }

     auto end = std::chrono::steady_clock::now();
     time_diff = std::chrono::duration<double>(end - start).count();
     block_time += time_diff;

     icount++;
     if (icount == nreport) {
       Rprintf("Sampled: %i of %i\n", (iter + 1), niter);
       ert = (block_time / nreport) * (niter - iter);
       Rprintf("Estimated remaining time: %3.2f minutes\n", ert / 60.0);
       icount = 0;
       block_time = 0.0;
     }

   }

   return latent_save;

 }
