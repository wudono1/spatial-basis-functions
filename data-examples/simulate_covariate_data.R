################################################################################
# Simulate covariate data
################################################################################
rm(list = ls())

outpath <- "../outputs/sim1-downscaling/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(MASS)
library(dplyr)
library(fields)

set.seed(42)

#-------------------------------------------------------------------------------
# Generate true fine-resolution covariate fields
#-------------------------------------------------------------------------------

# current constraint for MCMC simulation to work:
# r ( x_seq_knots * y_seq_knots) > MD (KK_x * kk_y)

#setting up grid of covariance predictors

# 8 range values, 5 sill values
# each (range, sill) pair is a different GP, so 40 GP's total
range_vals <- seq(0.1, 0.8, by = 0.1)
sill_vals <- seq(1, 5, by = 1)

# build fine grid axes (50 points in x, 60 points in y)
x_seq <- seq(0, 1, length = 50)
y_seq <- seq(0, 1, length = 60)
fine_grid <- as.matrix(expand.grid(x_seq, y_seq))


x_seq_knots <- seq(0, 1, length = 10)
y_seq_knots <- seq(0, 1, length = 12)

fine_grid_knots <- as.matrix(expand.grid(x_seq_knots, y_seq_knots))
num_knots <- nrow(fine_grid_knots)

knot_spacing = 0.1  # approximate

# keep scalar for now for simplicity
length_scale <- knot_spacing

# can do the array version later
#length_scale <- seq(0.5, 2.0, by = 0.5) * knot_spacing


colnames(fine_grid) <- c("lon", "lat")
nn <- nrow(fine_grid)

# Simulate from GPB models
# fine-resolution covariate field
fine_xx <- array(NA, dim = c(nn, length(range_vals), length(sill_vals)))

centroid_distmat <- rdist(fine_grid_knots, fine_grid_knots)

# calculating distance between each point on grid and each BF
# dimensions (nn, num_knots)
D <- rdist(fine_grid, fine_grid_knots)

# radial basis function calculation
# will multiple with alpha, which is calculated in the nested loop
Phi <- exp(-D^2 / (2 * length_scale^2))


# psi     = range_vals[i] : correlation range of the basis coefficients
# sigmasq = sill_vals[j]  : marginal variance of the basis coefficients
#
# alpha ~ N(0, sigmasq * K_psi), with K_psi = exp(-d / psi) over the knot
# centroids (manuscript Eq. 8). This is deliberately the same construction as
# expCorr() in the sampler, so the simulation and the MCMC build an identical
# K_psi from the same centroid_distmat.
mu_alpha <- rep(0, num_knots)

for (i in seq_along(range_vals)) {
  for (j in seq_along(sill_vals)) {

    # cov(alpha) = sigmasq * K_psi
    K_psi <- sill_vals[j] * exp(-centroid_distmat / range_vals[i])

    alpha <- mvrnorm(n = 1, mu = mu_alpha, Sigma = K_psi)

    # multiply Phi and alpha to get the complete field
    xx <- Phi %*% alpha

    # store the data
    fine_xx[, i, j] <- xx
  }
}

# old code: the pre-basis-function version. It simulated the field directly
# as a GP over the fine grid (manuscript Eq. 1) rather than as Phi %*% alpha
# (Eq. 8). Kept for reference; needs library(SpatialExtremes) if revived.


# for (i in seq(range_vals)) {
#   
#   for (j in seq(sill_vals)) {
#     
#     # Spatial Extremes::rgp- simulates Gaussian random field
#     # rgp(n, coord, cov.mod = "powexp", mean = 0, nugget = 0, sill = 1, range
#     #     = 1, smooth = 1, grid = FALSE, control = list())
#     
#     # sill = variance of the process
#     # range = how quickly spatial correlation decreases with distance
#     
#     # smooth = 1 means standard exponential covariance
#     # nugget = 0 means no non-spatial-related variance
#     # t() is transpose function
#     
#     # need to replace this with basis function models
#     xx <- t(rgp(1, fine_grid, cov.mod = "powexp", nugget = 0, smooth = 1,
#                 sill = sill_vals[j], range = range_vals[i]))
#     fine_xx[, i, j] <- xx
#     
#     dat <- data.frame(fine_grid, xx = xx)
#     
#   }
# }

#-------------------------------------------------------------------------------
# Aggregate to coarse-resolution data
#-------------------------------------------------------------------------------
kk_x <- 5
kk_y <- 6
mm <- kk_x * kk_y

fine_grid <- as.data.frame(fine_grid)
fine_grid$x_cut <- as.numeric(cut(fine_grid[, "lon"], kk_x))
fine_grid$y_cut <- as.numeric(cut(fine_grid[, "lat"], kk_y))
grid_index <- expand.grid(x_cut = 1:kk_x, y_cut = 1:kk_y)
grid_index$par_label <- as.numeric(row.names(grid_index))
fine_grid <- left_join(fine_grid, grid_index)

coar_lon <- aggregate(fine_grid$lon, by = list(fine_grid$par_label), mean)$x
coar_lat <- aggregate(fine_grid$lat, by = list(fine_grid$par_label), mean)$x
coar_grid <- data.frame(lon = coar_lon, lat = coar_lat)

# coarse resolution covariate field (aggregated from fine grained covariates)
# in 4.1, pretend that we only observe coarse-resolution covariates
coar_xx <- array(NA, dim = c(mm, length(range_vals), length(sill_vals)))

for (i in seq_along(range_vals)) {
  
  for (j in seq_along(sill_vals)) {
    
    xx <- fine_xx[, i, j]
    xxb <- aggregate(xx, by = list(fine_grid$par_label), mean)$x
    coar_xx[, i, j] <- xxb
    
    coar_dat <- data.frame(coar_grid, xxb = xxb)
    
  }
  
}


#-------------------------------------------------------------------------------
# Output covariate data
#-------------------------------------------------------------------------------

# original code
# save(fine_grid, coar_grid, fine_xx, coar_xx, range_vals, sill_vals,
#      file = paste0(outpath, "covar_dat.rdata"))


# new code
save(fine_grid, coar_grid, range_vals, sill_vals, length_scale, fine_grid_knots,
     knot_spacing, fine_xx, coar_xx, Phi, centroid_distmat,
     file = paste0(outpath, "covar_dat.rdata"))
