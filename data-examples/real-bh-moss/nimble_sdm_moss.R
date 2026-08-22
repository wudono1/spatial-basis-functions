################################################################################
# Section 4.3 Case study
# Bayesian implementation of the Bernoulli Berkson-GLM using the 'nimble' package
################################################################################
rm(list = ls())

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

inpath1 <- "../outputs/real-bh-moss/"
inpath2 <- "../data/"
outpath <- "../outputs/real-bh-moss/nimble/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

library(terra)
library(RColorBrewer)
library(dplyr)
library(raster)
library(sp)
library(sf)


#-------------------------------------------------------------------------------
# Import data
#-------------------------------------------------------------------------------
# downscaled MC samples of soil temperature
soil_temp <- terra::rast(paste0(inpath1, "downscaled_soil_temp.tif"))
bh_crds <- as.matrix(crds(soil_temp))
bh_map <- sf::st_read(paste0(inpath2, "bh/BH_mask.shp"))

# moss occurrence data
load(paste0(inpath1, "survey_moss.rdata"))
moss_dat <- data.frame(eastings = as.integer(survey_moss$Eastings),
                       northings = as.integer(survey_moss$Northings),
                       moss = survey_moss$moss)


#-------------------------------------------------------------------------------
# Plot mean and sd of the downscaled MC samples
#-------------------------------------------------------------------------------
map_trunc <- c(564600, 597500, 2639400, 2660500)

soil_temp_bh <- terra::mask(soil_temp$soil_temp_mean, bh_map)
soil_temp_bh <- terra::crop(soil_temp_bh, map_trunc)

getPalette = colorRampPalette(brewer.pal(11, "RdBu"))
rdbu <- rev(getPalette(256))

png(paste0(outpath, "south_bh_soil_temp_mean.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(soil_temp_bh, col = rdbu, xlab = "", ylab = "", plg = list(cex = 1.5), 
     xlim = c(map_trunc[1], map_trunc[2] + 100),
     mar = c(4.1, 0, 1.9, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 2.9, cex = 1.8)
points(moss_dat[which(moss_dat$moss == 1), ], lwd = 2, pch = 2)
dev.off()

soil_temp_sd_bh <- terra::mask(soil_temp$soil_temp_sd, bh_map)
soil_temp_sd_bh <- terra::crop(soil_temp_sd_bh, map_trunc)

getPalette = colorRampPalette(brewer.pal(3, "YlGn"))
YlGn <- getPalette(256)

png(paste0(outpath, "south_bh_soil_temp_sd.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(soil_temp_sd_bh, col = YlGn, xlab = "", ylab = "", plg = list(cex = 1.5),
     xlim = c(map_trunc[1], map_trunc[2] + 100),
     mar = c(4.1, 0, 1.9, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 2.9, cex = 1.8)
dev.off()


#-------------------------------------------------------------------------------
# Stage 2: Fit GLM-Berkson using nimble
#-------------------------------------------------------------------------------
### step(i): Downscaled MC samples of soil temperature
covar_dat <- data.frame(eastings = as.integer(as.integer(bh_crds[, 1])),
                        northings = as.integer(as.integer(bh_crds[, 2])),
                        soil_temp_mean = terra::values(soil_temp$soil_temp_mean),
                        terra::values(soil_temp[[-(1:5)]]))
covar_dat$northings <- covar_dat$northings + 1
bh_dat <- left_join(moss_dat, covar_dat, by = c("eastings", "northings"))

# downscaled MC samples and predictive mean 
xx_ppd <- bh_dat[, -(1:4)]
xx_ppd_mu <- rowMeans(xx_ppd)

### step(ii): beta tilde
mean_out <- glm(moss ~ soil_temp_mean, family = binomial, data = bh_dat)
ppd_mu_be <- coef(mean_out)

### steps (iii, iv, v): 
KK <- ncol(xx_ppd) # number of MC sample
nn <- nrow(bh_dat) # number of survey sites
nburn <- 20000      # number of burn-in samples
niter <- 30000      # number of MCMC samples
nthin <- 5          # thinning degree

bks_be_sam_pool <- vector("list", length = KK)
bks_data_list <- vector("list", length = KK)

# The typical steps are nimbleModel -> configureMCMC -> buildMCMC -> compile model and MCMC -> runMCMC.  
# define the model
library(nimble)
code <- nimbleCode({
  beta0 ~ dflat()
  beta1 ~ dflat()
  for (i in 1:nn) {
    logit(p[i]) <- beta0 + beta1 * xx[i] + eta[i]
    zz[i] ~ dbern(p[i])
  }
})

# constants, data, and initial values
constants <- list(nn = nn)
for (k in 1:KK) {
  xx_k <- xx_ppd[, k]
  delta_k <- ppd_mu_be[2] * (xx_k - xx_ppd_mu)
  bks_data_list[[k]] <- list(zz = bh_dat$moss, xx = xx_ppd_mu, eta = delta_k)     
}
inits <- list(beta0 = 0, beta1 = 0)

# build and compile the model using an initial date set
data <- bks_data_list[[1]]
bks_nim <- nimbleModel(code = code, constants = constants, data = data, inits = inits)
c_bks_nim <- compileNimble(bks_nim)

# configure and compile mcmc
bks_mcmc_config <- configureMCMC(bks_nim)
bks_mcmc_config$addSampler(target = c("beta0", "beta1"), type = "RW_block", control = list(adaptive = TRUE))
bks_mcmc <- buildMCMC(bks_mcmc_config)
c_bks_mcmc <- compileNimble(bks_mcmc, project = bks_nim, resetFunctions = TRUE)

# draw samples in mcmc with the first MC sample
set.seed(42)
post_sams <- runMCMC(c_bks_mcmc, niter = niter, nburnin = nburn, thin = nthin)
bks_be_sam_pool[[1]] <- t(post_sams)

# draw samples with the rest of MC samples
for (k in 2:KK) {
  c_bks_nim$setData(eta = bks_data_list[[k]]$eta)
  post_sams <- runMCMC(c_bks_mcmc, niter = niter, nburnin = nburn, thin = nthin)
  bks_be_sam_pool[[k]] <- t(post_sams)
  cat(paste0("done the ", k, "-th replicate.\n"))
}

# Point and interval estimates of beta
bks_be_sam_pool <- simplify2array(bks_be_sam_pool)
rowMeans(bks_be_sam_pool)
apply(bks_be_sam_pool, 1, function(x) quantile(x, probs = c(0.025, 0.975)))


#-------------------------------------------------------------------------------
# Predict moss presence probabilities
#-------------------------------------------------------------------------------
N_d <- nrow(covar_dat)
nsam <- (niter - nburn) / nthin
bks_pred_sam_pool <- array(NA, dim = c(N_d, nsam, KK))

for (k in 1:KK) {
  
  xx_k <- covar_dat[, 3 + k]
  delta_k <- ppd_mu_be[2] * (xx_k - covar_dat$soil_temp_mean)
  
  for (b in 1:nsam) {
    yy_kb <- plogis(cbind(1, covar_dat$soil_temp_mean) %*% bks_be_sam_pool[, b, k] + delta_k)
    bks_pred_sam_pool[, b, k] <- yy_kb
  }
  
  cat(paste0("done the ", k, "-th fit.\n"))
  
}

# mean
bks_pred <- rowMeans(bks_pred_sam_pool) 
# sd
bks_sd <- apply(bks_pred_sam_pool, 1, sd) 
# 95% pred int
bks_pred_qq <- apply(bks_pred_sam_pool, 1, function(x) quantile(x, probs = c(0.025, 0.975))) 

moss_prob_dat <- data.frame(covar_dat$eastings,
                            covar_dat$northings,
                            prob_mean = bks_pred,
                            prob_sd = bks_sd,
                            prob_qq025 = bks_pred_qq[1, ],
                            prob_qq975 = bks_pred_qq[2, ])

moss_prob_ras <- rast(moss_prob_dat, crs = "+proj=utm +zone=47 +south +datum=WGS84 +units=m +no_defs")
moss_prob_ras <- terra::mask(moss_prob_ras, bh_map)
moss_prob_ras <- terra::crop(moss_prob_ras, map_trunc)


#-------------------------------------------------------------------------------
# Plot moss presence probabilities and uncertainties
#-------------------------------------------------------------------------------
getPalette = colorRampPalette(brewer.pal(3, "Oranges"))
orang <- getPalette(256)

prob_range <- range(c(terra::values(moss_prob_ras$prob_mean), 
                      terra::values(moss_prob_ras$prob_qq025), 
                      terra::values(moss_prob_ras$prob_qq975)), na.rm = TRUE)

png(paste0(outpath, "nimble_south_bh_moss_prob_mean.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(moss_prob_ras$prob_mean, col = orang, range = prob_range, plg = list(cex = 1.5),
     xlab = "", ylab = "", xlim = c(map_trunc[1], map_trunc[2] + 100),
     mar = c(4.1, 0, 1.9, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 2.9, cex = 1.8)
dev.off()

png(paste0(outpath, "nimble_south_bh_moss_prob_025.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(moss_prob_ras$prob_qq025, col = orang, range = prob_range, plg = list(cex = 1.5),
     xlab = "", ylab = "", xlim = c(map_trunc[1], map_trunc[2] + 100),
     mar = c(4.1, 0, 1.9, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 2.9, cex = 1.8)
dev.off()

png(paste0(outpath, "nimble_south_bh_moss_prob_975.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(moss_prob_ras$prob_qq975, col = orang, range = prob_range, plg = list(cex = 1.5),
     xlab = "", ylab = "", xlim = c(map_trunc[1], map_trunc[2] + 100),
     mar = c(4.1, 0, 1.9, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 2.9, cex = 1.8)
dev.off()

getPalette = colorRampPalette(brewer.pal(3, "BuPu"))
bupu <- getPalette(256)

png(paste0(outpath, "nimble_south_bh_moss_prob_sd.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(moss_prob_ras$prob_sd, col = bupu, plg = list(cex = 1.5), 
     xlab = "", ylab = "", xlim = c(map_trunc[1], map_trunc[2] + 100),
     mar = c(4.1, 0, 1.9, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 2.9, cex = 1.8)
dev.off()


