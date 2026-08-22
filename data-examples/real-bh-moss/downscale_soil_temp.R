################################################################################
# Section 4.3 Case study
# Downscale soil temperature
################################################################################
rm(list = ls())

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

inpath  <- "../data/"
outpath <- "../outputs/real-bh-moss/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

library(terra)
library(RColorBrewer)
library(dplyr)
library(raster)
library(sp)
library(sf)
library(spds)


#-------------------------------------------------------------------------------
# Import and prepare data
#-------------------------------------------------------------------------------
# Bunger Hills extent
bh_ext_utm <- c(548554.453407128, 623118.560071669, 2622783.03674264, 2686148.32354778)
bh_ext_utm_trunc <- c(555500, 605500, 2623500, 2670500)

# Create BAU centroids (note: below is not the only way to do so)
# Will align coordinates of different data with the BAU centroids
seq_eastings <- seq(bh_ext_utm_trunc[1], bh_ext_utm_trunc[2], by = 1000)
seq_northings <- seq(bh_ext_utm_trunc[3], bh_ext_utm_trunc[4], by = 1000)
bau_grid <- as.matrix(expand.grid(seq_eastings, seq_northings))

# Coarse-resolution soil temperature (to be downscaled)
coar_temp <- rast(paste0(inpath, "era5-land/soil_temp_raw_mean.tif"))
names(coar_temp) <- "temperature"

res_vec <- c(5, 11)
coar_temp_utm <- terra::project(coar_temp, 
                                y = "+proj=utm +zone=47 +south +datum=WGS84 +units=m +no_defs", 
                                res = res_vec * 1000)
coar_temp_bh_orig <- terra::crop(coar_temp_utm, bh_ext_utm)
new_crds <- t(apply(crds(coar_temp_bh_orig, na.rm = FALSE), 1, function(x) x + c(383.3, 418)))
coar_temp_bh <- rast(data.frame(new_crds, values(coar_temp_bh_orig)),
                     crs = "+proj=utm +zone=47 +south +datum=WGS84 +units=m +no_defs")
coar_temp_bh <- terra::crop(coar_temp_bh, bh_ext_utm_trunc)
coar_temp_bh <- coar_temp_bh$temperature - 273.15

# BAU-level covariates (landscape data)
pred_var <- c("REMA_1km_dem_filled.tif",
              "slope_1km.tif",
              "aspect_1km.tif")
varpath <- paste0(inpath, "dem")
vardir <- paste(varpath, pred_var, sep = "/")
pred_stack <- rast(vardir)
names(pred_stack) <- c("Elevation", "Slope", "Aspect")

pred_stack_utm <- terra::project(pred_stack, y = crs(coar_temp_bh), res = 1000)
pred_stack_bh_fine <- terra::crop(pred_stack_utm, coar_temp_bh)
new_crds <- t(apply(crds(pred_stack_bh_fine, na.rm = FALSE), 1, function(x) x + c(-77.5, -215)))
pred_stack_bh_fine <- rast(data.frame(new_crds, values(pred_stack_bh_fine)),
                     crs = "+proj=utm +zone=47 +south +datum=WGS84 +units=m +no_defs")

pred_stack_bh_coarse <- terra::aggregate(pred_stack_bh_fine, fact = c(11, 5))
pred_stack_bh_fine_full <- as.data.frame(cbind(crds(pred_stack_bh_fine, na.rm = FALSE), values(pred_stack_bh_fine)))
pred_stack_bh_fine_full <- pred_stack_bh_fine_full %>% arrange(x, y)

#-------------------------------------------------------------------------------
# Disaggregate blocks to obtain BAU centroids
#-------------------------------------------------------------------------------
# Data frame of the coarse-resolution soil temperature
coarse_grid_locs <- crds(coar_temp_bh, na.rm = FALSE)
coarse_grid_vals <- values(coar_temp_bh)
coarse_dat_full <- as.data.frame(cbind(coarse_grid_locs, coarse_grid_vals))
coarse_dat_full <- coarse_dat_full %>% arrange(x, y)
coarse_dat_full <- arrange(coarse_dat_full, coarse_dat_full[, 3])

# Disaggregate blocks to to get fine resolution coordinates
fine_temp <- disagg(coar_temp_bh, fact = rev(res_vec))

# Data frame of the fine-resolution data
# (note that the fine-resolution values are coarse-resolution values;
# we haven't downscaled the data yet)
fine_grid_locs <- crds(fine_temp, na.rm = FALSE)
fine_grid_vals <- values(fine_temp)
fine_dat_full <- as.data.frame(cbind(fine_grid_locs, fine_grid_vals))
fine_dat_full <- fine_dat_full %>% arrange(x, y)
fine_dat_full <- cbind(fine_dat_full, pred_stack_bh_fine_full[, -(1:2)])
fine_dat_full <- arrange(fine_dat_full, fine_dat_full[, 3])
colnames(fine_dat_full)[1:2] <- c("eastings", "northings")
  
coarse_dat <- na.omit(coarse_dat_full)
fine_dat <- na.omit(fine_dat_full)

# Aggregation matrix
agg_mat <- array(0, dim = c(nrow(coarse_dat), nrow(fine_dat)))
ncell <- prod(res_vec)
for (j in 1:nrow(coarse_dat)) {
  agg_mat[j, (ncell * (j - 1) + 1):(ncell * j)] <- 1 / ncell
}

# Plot coarse-resolution data and BAUs overlaid with sites.
load(paste0(outpath, "survey_moss.rdata"))
bh_moss <- as.data.frame(survey_moss[, c("Eastings", "Northings")])
getPalette = colorRampPalette(brewer.pal(11, "RdBu"))
rdbu <- rev(getPalette(256))
png(paste0(outpath, "coar_temp_bau_sites.png"), width=1000, height=600, pointsize=20, bg = "transparent")
plot(coar_temp_bh, col = rdbu, xlab = "", ylab = "", plg = list(cex = 1.5),
     mar = c(4.1, 0, 1.1, 1.6), pax = list(cex.axis = 2, mgp = c(3, 1, 0)))
abline(v = c(seq_eastings[1]-500, seq_eastings + 500), col = 'grey')
abline(h = c(seq_northings[1]-500, seq_northings + 500), col = 'grey')
abline(v = seq(seq_eastings[1]-500, tail(seq_eastings, 1) + 500, by = 5000), 
       col = 'mediumseagreen', lwd = 2)
abline(h = seq(seq_northings[1]-500, tail(seq_northings, 1) + 500, by = 11000), 
       col = 'mediumseagreen', lwd = 2)
points(bh_moss, col = 'black', pch = 19, cex = 0.6)
mtext("Eastings", side = 1, line = 3.8, cex = 1.8)
mtext("Northings", side = 2, line = 0, cex = 1.8)
dev.off()

#-------------------------------------------------------------------------------
# Fit the GPB model
#-------------------------------------------------------------------------------
grid_loc <- fine_dat[, c("eastings", "northings")]
grid_distmat <- as.matrix(dist(grid_loc))
dist_range <- range(grid_distmat[grid_distmat!=0])

zz <- coarse_dat[, 3]
coar_xx <- agg_mat %*% as.matrix(cbind(1, fine_dat[,4:6]))

mod_param <- list(bau_distmat = grid_distmat, aggmat = agg_mat)
prior <- list(phi_unif = c(dist_range[1], v_phi = dist_range[2] / 3), sigmasq_invgamma = c(3, 2))
starting <- list(phi = 3000, sigmasq = var(zz), regcoef = matrix(c(mean(zz), rep(0, ncol(coar_xx) - 1))))
tuning <- list(se_phi = 3)

mcmc_settings <- list(niter = 8000, nreport = 500)
  
set.seed(42)
runtime <- system.time(
  out <- gpbau(coarse_resp = zz,
                     coarse_covar = coar_xx,
                     prior = prior,
                     starting = starting,
                     tuning = tuning,
                     mod_param = mod_param, 
                     mcmc_settings = mcmc_settings)
)
 
# Thinning of the MCMC chain 
out <- thinning(out, nburn = 3000, nthin = 10)


#-------------------------------------------------------------------------------
# Prediction 
#-------------------------------------------------------------------------------
fine_covar <- as.matrix(cbind(1, fine_dat[,4:6]))

cat("prediction starts....")
runtime <- system.time(pred_temp <- predict(out, fine_covar, type = "joint", nreport = 100))
cat("Done!\n")
cat(paste0("Running time: ", round(runtime[3]/60, 3), " minutes.\n"))

soil_temp_mean <- rowMeans(pred_temp)
soil_temp_qq <- t(apply(pred_temp, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))))
colnames(soil_temp_qq) <- c("qq025", "qq050", "qq975")
soil_temp_sd <- apply(pred_temp, 1, sd)

soil_temp_dat <- data.frame(Eastings = fine_dat_full$eastings, 
                            Northings = fine_dat_full$northings,
                            soil_temp_mean = soil_temp_mean,
                            soil_temp_sd = soil_temp_sd,
                            soil_temp_qq = soil_temp_qq,
                            soil_temp = pred_temp)

soil_temp_ras <- rast(soil_temp_dat, type = "xyz", crs = crs(coar_temp_bh))

terra::writeRaster(soil_temp_ras,
                   filename = paste0(outpath, "downscaled_soil_temp.tif"),
                   overwrite = TRUE)

