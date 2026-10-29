#' @title Function for COS-based statistical downscaling
#' 
#' @description This function fits a Gaussian process over basic areal units (GPB) to
#'              downscale coarse-resolution data via Markov Chain Monte Carlo simulations.
#'              Section 2.2.2. in MEE. Outputs posterior sample from MEE 2.2.2. and S4 MCMC algorithm
#' 
#' @param coarse_resp a vector of coarse-resoultion data of length \eqn{M_D}.
#' @param coarse_covar an \eqn{M_D \times p} matrix of covariates for the downscaling variable 
#'                     of interest. If an intercept is desired, the first column of the matrix should be a vector of ones.
#' @param prior a list of priors. Each element of the list corresponds to a combination of an unknown parameter 
#'               and its prior distribution; for example, "\code{phi_unif}" and "\code{sigmasq_invgamma}".
#' @param starting a list of starting values. Each element of the list corresponds to an unknown parameter.
#' @param tuning a list of tuning parameters. Each element of the list corresponds to an unknown parameter.
#'               The value portion of each element defines the variance of the random walk Metropolis sampler 
#'               proposal distribution for the unknown parameter.
#' @param mod_param a list of parameters needed for model fitting. Currently the function requires:
#'                  \code{aggmat}: an \eqn{M_D \times N_D} matrix that aggregates fine-resolution data to coarse-resolution data;
#'                  \code{bau_distmat}: an \eqn{N_D \times N_D} matrix about pairwise-locaiton distances.
#' @param mcmc_settings a list of Markov chain Monte Carlo (MCMC) simulation parameters. 
#'        \code{n_iter}: number of iterations; 
#'        \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If \code{verbose = TRUE};
#'        \code{n_report} defines the interval to report MCMC progress and Metropolis acceptance rates.
#' @param verbose logical; if true, model specification, MCMC progress, and Metropolis acceptance rates
#'                are printed to the screen.
#' 
#' @return
#' An object of class \code{gpbau}. The return object is a list that comprises:
#' 
#' \tabular{ll}{
#' 
#' \code{post_samples} \tab a list of posterior samples. Each element of the list corresponds to an unknown parameter. \cr
#' \tab \cr
#' \code{coarse_data} \tab a list that consists of input data 
#'         (\code{response}, \code{covars}, and \code{coords}). \cr
#' \tab \cr
#' \code{mod_spec}: \tab a list that consists of model specifications 
#'         (\code{mod_param}: parameters needed to fit models; 
#'          \code{priors}: priors for unknown parameters). \cr
#' \tab \cr
#' \code{mcmc_params}: \tab a list that consists of unknown parameters' starting values and Metropolis sampler proposal distribution variances
#'         used in the MCMC simulation (
#'         \code{starting}: starting values;
#'         \code{tuning}: Metropolis sampler proposal distribution variances). \cr
#' \tab \cr
#' \code{mh_data}: \tab a list that consists of each unknown parameter's Metropolis 
#'                      sampler acceptance and rejection indicators.\cr
#' \tab \cr
#' \code{runtime}: \tab running time for MCMC simulation calculated using \code{system.time}. \cr
#' 
#'}
#'
#' @export
#'
gpbau <- function(coarse_resp,
                  coarse_covar,
                  prior,
                  starting = NULL,
                  tuning,
                  mod_param,
                  mcmc_settings,
                  verbose = TRUE) {


    xxb <- coarse_resp
    rrb <- coarse_covar

    #----------------------------------------------------
    # Compute aggregate matrix
    #----------------------------------------------------
    # To be added in...Now the aggregate matrix must be an element of mod_param
    aggmat <- mod_param$aggmat
    Nd <- ncol(aggmat)
    Md <- nrow(aggmat)

    #----------------------------------------------------
    # Check priors
    #----------------------------------------------------
    prior2 <- vector("list")

    prior2$u_phi <- prior$phi_unif[1]
    prior2$v_phi <- prior$phi_unif[2]
    prior2$u_sigmasq <- prior$sigmasq_invgamma[1]
    prior2$v_sigmasq <- prior$sigmasq_invgamma[2]

    #----------------------------------------------------
    # Check tuning parameters
    #----------------------------------------------------
    se_phi <- tuning$se_phi
    if (is.null(se_phi)) stop("error: MH tuning parameter for phi must be specified.")
    if (se_phi <= 0) stop("error: MH tuning parameter for phi must be positive.")    

    #----------------------------------------------------
    # Check starting values
    #----------------------------------------------------
    if (is.null(starting)) {
      
      starting <- list(
        phi = mean(prior$phi_unif), 
        sigmasq = var(xxb), 
        regcoef = matrix(c(mean(xxb)))
        )
    }

    #----------------------------------------------------
    # verbose
    #----------------------------------------------------
    if (verbose) {
        cat("--------------------------------------------------------------\n")
        cat(" << Spatial statistical downscaling using GPB >>\n")
        cat("--------------------------------------------------------------\n")
        cat(paste0("Number of coarse-resolution data: ", Md, "\n"))
        cat(paste0("Number of fine-resolution BAUs: ", Nd, "\n"))
    }
    
    #----------------------------------------------------
    # MCMC
    #----------------------------------------------------
    cat("----------------------------------------\n")
    cat("Running MCMC...\n");
    cat(paste0("Number of iterations: ", mcmc_settings$niter, "\n"))
    cat("Report progress every ", mcmc_settings$nreport, " iterations\n");
    cat("----------------------------------------\n")
    runtime <- system.time(
        mcmc_out <- fitGPBAU(xxb, rrb, prior2, starting, tuning, mod_param, mcmc_settings)
    )
    
    #----------------------------------------------------
    # Output
    #----------------------------------------------------
    out <- list(
        post_sams = mcmc_out$post_sams,
        coarse_data = list(coarse_resp = coarse_resp, 
                           coarse_covar = coarse_covar),
        mod_spec = list(mod_param = mod_param,
                        prior = prior),
        mcmc_params = list(starting = starting, 
                           tuning = tuning),
        mcmc_settings = mcmc_settings,
        mh_data = mcmc_out$mh_data,
        runtime = runtime
    )

    class(out) <- "gpbau"
    
    out

}

#' @title Function for obtaining downscaling samples
#'
#' @description This function produces downscaling samples of the coarse-resolution
#' variable using GPB via posterior predictive distribution.
#'
#' @param object an object of class gpbau.
#' @param fine_covar an \eqn{N_D \times p} matrix of covariates at the fine resolution.
#'                   Each column of the matrix corresponds to the column of the matrix in 
#'                   \code{coarse_covar} when it was used to fit GPB.
#' @param type A quoted keyword that specifies the type of prediction.
#'             Supported keywrods are \code{"independent"} and \code{"joint"}.
#' @param verbose logical; if true, progress of the prediction is printed to the screen.
#' @param nreport If \code{verbose = TRUE}, \code{n_report} defines the interval to report
#'                the progress of prediction.
#' @param ... additional arguments. No additional arguments are supported currently.
#' 
#' @return
#' The return object is an \eqn{N_D \times K} matrix, where \eqn{K} is the number 
#' of posterior predictive samples.
#' 
#' @exportS3Method 
predict.gpbau <- function(object, 
                          fine_covar,
                          type = "independent",
                          verbose = TRUE,
                          nreport = NULL,
                          ...) {

    coarse_resp <- object$coarse_data$coarse_resp

    if (!is.null(object$post_sams_thin)){
        post_sams <- object$post_sams_thin
    } else {
        post_sams <- object$post_sams
        post_sams$beta <- t(post_sams$beta)
    }

    mod_param <- object$mod_spec$mod_param
    niter <- length(post_sams$sigmasq)
    
    if (is.null(nreport)) {
        nreport <- object$mcmc_settings$nreport
        if (is.null(nreport)) {
            nreport <- niter / 10
        }
    }

    if (type == "independent") {
        pred_out <- predGPBAU3(coarse_resp, fine_covar, post_sams, mod_param, nreport)
    } else if (type == "joint") {
        pred_out <- predGPBAU2(coarse_resp, fine_covar, post_sams, mod_param, nreport)
    } else {
        stop("error: the type must be specified: independent or joint.")
    }

    pred_out

}


#' @title COS-based statistical downscaling using a basis-function GPB model
#'
#' @description Basis-function version of \code{\link{gpbau}}. Fits the GPB model in which
#'              the spatial process is represented by \eqn{r} basis functions with coefficients
#'              \eqn{\alpha \sim N(0, \sigma^2 K_\psi)}, so that the block-level covariance is
#'              \eqn{K_B(\psi) = W_B \phi(S)^\top K_\psi \phi(S) W_B^\top}.
#'
#' @param coarse_resp a vector of coarse-resolution data of length \eqn{M_D}.
#' @param coarse_covar an \eqn{M_D \times p} matrix of covariates.
#' @param prior a list of priors: \code{phi_unif}, \code{sigmasq_invgamma}.
#' @param starting a list of starting values.
#' @param tuning a list with \code{se_phi}.
#' @param mod_param a list requiring \code{aggmat} (\eqn{M_D \times N_D}),
#'                  \code{Phi} (\eqn{N_D \times r} basis functions at BAU centroids), and
#'                  \code{centroid_distmat} (\eqn{r \times r} distances among basis centroids).
#' @param mcmc_settings a list with \code{niter}, \code{nreport}.
#' @param verbose logical; print model specification and MCMC progress.
#'
#' @return An object of class \code{gpbau_new}.
#'
#' @export
#'
gpbau_new <- function(coarse_resp,
                      coarse_covar,
                      prior,
                      starting = NULL,
                      tuning,
                      mod_param,
                      mcmc_settings,
                      verbose = TRUE) {

    xxb <- coarse_resp
    rrb <- coarse_covar

    #----------------------------------------------------
    # Compute aggregate matrix
    #----------------------------------------------------
    aggmat <- mod_param$aggmat
    Nd <- ncol(aggmat)
    Md <- nrow(aggmat)

    #----------------------------------------------------
    # >>> NEW: Check basis-function inputs <<<
    #
    # Under the delta(S) = 0 assumption, K_B = Psi_B K_psi Psi_B' has
    # rank at most min(Md, r). It must have rank Md to be invertible, so
    # r >= Md is a WELL-POSEDNESS condition, not a tuning preference.
    # These checks run once, before the (expensive, often parallel) MCMC.
    #----------------------------------------------------
    Phi <- mod_param$Phi
    centroid_distmat <- mod_param$centroid_distmat

    # (a) presence
    if (is.null(Phi))
        stop("error: mod_param$Phi (Nd x r basis functions) must be specified.")
    if (is.null(centroid_distmat))
        stop("error: mod_param$centroid_distmat (r x r) must be specified.")

    # (b) REMOVED: mod_param$bau_distmat is no longer read by fitGPBAU_new()
    #     or by predGPBAU2_new()/predGPBAU3_new(). Callers running only the
    #     basis-function method can drop it from mod_param entirely, which
    #     saves shipping an Nd x Nd matrix (~72 MB at Nd = 3000) to every
    #     parallel worker.

    Phi <- as.matrix(Phi)
    centroid_distmat <- as.matrix(centroid_distmat)
    rr <- ncol(Phi)

    # (c) conformability
    if (nrow(Phi) != Nd)
        stop(paste0("error: nrow(Phi) = ", nrow(Phi),
                    " must equal ncol(aggmat) = ", Nd, "."))
    if (nrow(centroid_distmat) != rr || ncol(centroid_distmat) != rr)
        stop(paste0("error: centroid_distmat must be ", rr, " x ", rr,
                    " to match ncol(Phi)."))

    # (d) rank condition: necessary
    if (rr < Md)
        stop(paste0("error: number of basis functions r = ", rr,
                    " is less than the number of blocks Md = ", Md,
                    ". With delta(S) = 0 this makes K_B singular. Increase r."))

    # (e) rank condition: sufficient. r >= Md does not guarantee Psi_B attains
    #     row rank Md -- smooth bases can collapse after block-averaging.
    PsiB <- aggmat %*% Phi
    psib_rank <- qr(PsiB)$rank
    if (psib_rank < Md)
        stop(paste0("error: rank(W_B %*% Phi) = ", psib_rank,
                    " < Md = ", Md, ". Basis is rank-deficient after aggregation; ",
                    "widen knot spacing or use a more localised basis."))

    #----------------------------------------------------
    # Check priors
    #----------------------------------------------------
    prior2 <- vector("list")

    prior2$u_phi <- prior$phi_unif[1]
    prior2$v_phi <- prior$phi_unif[2]
    prior2$u_sigmasq <- prior$sigmasq_invgamma[1]
    prior2$v_sigmasq <- prior$sigmasq_invgamma[2]

    #----------------------------------------------------
    # Check tuning parameters
    #----------------------------------------------------
    se_phi <- tuning$se_phi
    if (is.null(se_phi)) stop("error: MH tuning parameter for phi must be specified.")
    if (se_phi <= 0) stop("error: MH tuning parameter for phi must be positive.")

    #----------------------------------------------------
    # Check starting values
    #----------------------------------------------------
    if (is.null(starting)) {

      starting <- list(
        phi = mean(prior$phi_unif),
        sigmasq = var(xxb),
        regcoef = matrix(c(mean(xxb)))
        )
    }

    #----------------------------------------------------
    # >>> NEW: keep the starting range strictly inside the prior <<<
    #
    # logit(psi, u, v) = log((psi - u)/(v - psi)) is -Inf/+Inf at the bounds, so
    # cur_logpost becomes -Inf, every MH ratio evaluates to NaN, and the chain
    # never moves. Because sigmasq and gam are only updated on acceptance, they
    # freeze at their starting values too -- producing zero-width posterior
    # "intervals" rather than an error. findinit_new() can return a boundary
    # value whenever the profile likelihood is maximised at an endpoint.
    #----------------------------------------------------
    psi_lo <- prior2$u_phi
    psi_hi <- prior2$v_phi
    inset  <- 0.01 * (psi_hi - psi_lo)          # 1% of the prior range
    if (starting$phi <= psi_lo + inset || starting$phi >= psi_hi - inset) {
        warning(paste0("starting phi (", signif(starting$phi, 4),
                       ") is at or near a prior bound; moving it inside. ",
                       "This usually signals a weakly identified range parameter."))
        starting$phi <- min(max(starting$phi, psi_lo + inset), psi_hi - inset)
    }

    #----------------------------------------------------
    # >>> NEW: conditioning check at the starting psi <<<
    #
    # K_B can be full rank yet numerically unusable. Gaussian RBFs are
    # notoriously collinear, and block-averaging makes it worse. This is a
    # WARNING, not an error -- conditioning varies with psi over the chain.
    #----------------------------------------------------
    KB_start <- PsiB %*% exp(-centroid_distmat / starting$phi) %*% t(PsiB)
    eig <- eigen(KB_start, symmetric = TRUE, only.values = TRUE)$values
    kappa_KB <- max(eig) / min(eig)
    if (min(eig) <= 0 || kappa_KB > 1e12)
        warning(paste0("K_B is poorly conditioned at the starting phi (condition number ",
                       format(kappa_KB, digits = 3),
                       "). Cholesky may fail or lose accuracy during MCMC."))

    #----------------------------------------------------
    # verbose
    #----------------------------------------------------
    if (verbose) {
        cat("--------------------------------------------------------------\n")
        cat(" << Spatial statistical downscaling using basis-function GPB >>\n")
        cat("--------------------------------------------------------------\n")
        cat(paste0("Number of coarse-resolution data: ", Md, "\n"))
        cat(paste0("Number of fine-resolution BAUs: ", Nd, "\n"))
        cat(paste0("Number of basis functions: ", rr, "\n"))
        cat(paste0("rank(W_B %*% Phi): ", psib_rank, " (needs ", Md, ")\n"))
    }

    #----------------------------------------------------
    # MCMC
    #----------------------------------------------------
    cat("----------------------------------------\n")
    cat("Running MCMC...\n");
    cat(paste0("Number of iterations: ", mcmc_settings$niter, "\n"))
    cat("Report progress every ", mcmc_settings$nreport, " iterations\n");
    cat("----------------------------------------\n")
    runtime <- system.time(
        mcmc_out <- fitGPBAU_new(xxb, rrb, prior2, starting, tuning, mod_param, mcmc_settings)
    )

    #----------------------------------------------------
    # Output
    #----------------------------------------------------
    out <- list(
        post_sams = mcmc_out$post_sams,
        coarse_data = list(coarse_resp = coarse_resp,
                           coarse_covar = coarse_covar),
        mod_spec = list(mod_param = mod_param,
                        prior = prior),
        mcmc_params = list(starting = starting,
                           tuning = tuning),
        mcmc_settings = mcmc_settings,
        mh_data = mcmc_out$mh_data,
        runtime = runtime
    )

    class(out) <- "gpbau_new"

    out

}


#' @title Downscaling samples from a basis-function GPB fit
#'
#' @description Basis-function version of \code{\link{predict.gpbau}}. Produces posterior
#'              predictive samples of the fine-resolution field via conditional simulation
#'              from \eqn{[X_s \mid X_B]}, using
#'              \eqn{C_S = \sigma^2 \phi(S)^\top K_\psi \phi(S)}.
#'
#' @param object an object of class \code{gpbau_new}.
#' @param fine_covar an \eqn{N_D \times p} matrix of covariates at the fine resolution.
#' @param type a quoted keyword: \code{"independent"} (pointwise) or \code{"joint"}
#'             (spatially coherent draws). Use \code{"joint"} whenever downstream
#'             calculations involve functionals of the whole field.
#' @param verbose logical; if true, progress is printed to the screen.
#' @param nreport interval at which to report progress.
#' @param ... additional arguments (unused).
#'
#' @return An \eqn{N_D \times K} matrix of predictive samples.
#'
#' @exportS3Method
predict.gpbau_new <- function(object,
                              fine_covar,
                              type = "independent",
                              verbose = TRUE,
                              nreport = NULL,
                              ...) {

    # This method exists for two reasons, neither of which is a bare passthrough:
    #
    #  (1) S3 DISPATCH. predict() is generic; R looks up predict.<class>. Objects
    #      from gpbau_new() carry class "gpbau_new", so without this method
    #      predict() errors -- which is deliberate, since predGPBAU2/3 would
    #      otherwise apply the OLD covariance to basis-function draws and return
    #      silently wrong fields.
    #
    #  (2) ADAPTER WORK. The C++ takes bare matrices; this unpacks them from the
    #      fitted object -- picking thinned vs unthinned draws, fixing the beta
    #      orientation, and defaulting nreport.

    coarse_resp <- object$coarse_data$coarse_resp

    # thinning() already transposes beta to niter x p; the unthinned chains are
    # stored as p x niter by the C++, so they need transposing here. predGPBAU*_new
    # reads bb_save.row(iter), i.e. it expects niter x p.
    if (!is.null(object$post_sams_thin)){
        post_sams <- object$post_sams_thin
    } else {
        post_sams <- object$post_sams
        post_sams$beta <- t(post_sams$beta)
    }

    mod_param <- object$mod_spec$mod_param
    niter <- length(post_sams$sigmasq)

    if (is.null(nreport)) {
        nreport <- object$mcmc_settings$nreport
        if (is.null(nreport)) {
            nreport <- niter / 10
        }
    }

    # CHANGED vs predict.gpbau(): routes to the basis-function predictors.
    if (type == "independent") {
        pred_out <- predGPBAU3_new(coarse_resp, fine_covar, post_sams, mod_param, nreport)
    } else if (type == "joint") {
        pred_out <- predGPBAU2_new(coarse_resp, fine_covar, post_sams, mod_param, nreport)
    } else {
        stop("error: the type must be specified: independent or joint.")
    }

    pred_out

}

