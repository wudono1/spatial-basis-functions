################################################################################
# Section 4.1 Simulation study: Statistical downscaling with uncertainty
# Statistical downscaling and comparison with TPRS
################################################################################
rm(list = ls())

inpath  <- "../outputs/sim1-downscaling/" 
outpath <- "../outputs/sim1-downscaling/downscaling-comp/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

library(mgcv)

# from downscale_gpbau.R (fine grained posterior samples)
load(paste0(inpath, "sim_gpbau_pred.RData"))

# from simulate_covariate_data.R (contains TRUE fine-grained covariate field)
load(paste0(inpath, "covar_dat.rdata"))

set.seed(42)

#-------------------------------------------------------------------------------
# GPB
#-------------------------------------------------------------------------------

# building 95% confidence interval
gpbau_pred <- xx_pred[[1]][[1]]
gpbau_pred_qq <- apply(gpbau_pred, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

xx <- fine_xx[, 1, 1]
low_ind <- xx >= gpbau_pred_qq[1,]
hig_ind <- xx <= gpbau_pred_qq[2,]
gpbau_pred_ci <- mean(low_ind * hig_ind)


# plotting coverage on map
gpbau_field <- data.frame(lon = fine_grid$lon, lat = fine_grid$lat, gpbau = rowMeans(gpbau_pred))
plot_save <- ggplot(gpbau_field, aes(x = lon, y = lat, fill = gpbau)) +
  geom_raster(interpolate = F) + theme_bw() +
  scale_fill_gradient2(high = "red", low = "blue", mid = "white", midpoint = 0, limits = range(xx),
                       space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  labs(fill = "", x = "Easting", y = "Northing") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "gpbau_pred_mean.png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()

gpbau_field$ci_ind <- hig_ind * low_ind
plot_save <- ggplot(gpbau_field, aes(x = lon, y = lat, fill = as.factor(ci_ind))) +
  geom_raster(interpolate = F) + theme_bw() +
  labs(fill = "", x = "Easting", y = "Northing") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "gpbau_pred_95cover.png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()

#-------------------------------------------------------------------------------
# TPRS (thin-plate splines)
#-------------------------------------------------------------------------------
# cross validation to select k
Md <- nrow(coar_xx)
k_start <- 4
cover_rate <- array(NA, dim = Md - k_start)
rmse <- array(NA, dim = Md - k_start)
coar_dat <- data.frame(lon = coar_grid$lon, lat = coar_grid$lat, xxb = coar_xx[, 1, 1])
for (k in k_start:Md) {
  tps <- gam(xxb ~ s(lon, lat, bs = "tp", k = k), data = coar_dat)
  tps_pred <- predict(tps, fine_grid[,1:2], se.fit = TRUE)
  tps_pred_qq025 <- tps_pred$fit - 2 * tps_pred$se.fit
  tps_pred_qq975 <- tps_pred$fit + 2 * tps_pred$se.fit
  low_ind <- xx >= tps_pred_qq025
  hig_ind <- xx <= tps_pred_qq975
  cover_ind <- low_ind * hig_ind
  cover_rate[k - k_start + 1] <- mean(cover_ind)
}
k <- which.max(cover_rate) + k_start - 1

# fit the model with the best k
tps <- gam(xxb ~ s(lon, lat, bs = "tp", k = k),  data = coar_dat)
tps_pred <- predict(tps, fine_grid[,1:2], se.fit = TRUE)

tps_field <- as.data.frame(cbind(fine_grid[,1:2], tps_pred$fit, tps_pred$se.fit))
colnames(tps_field) <- c("lon", "lat", "tps_pred", "tps_pred_se")

tps_field$tps_pred_qq025 <- tps_pred$fit - 2 * tps_pred$se.fit
tps_field$tps_pred_qq975 <- tps_pred$fit + 2 * tps_pred$se.fit

low_ind <- xx >= tps_field$tps_pred_qq025
hig_ind <- xx <= tps_field$tps_pred_qq975
tps_field$cover_ind <- low_ind * hig_ind
mean(tps_field$cover_ind)
  
plot_save <- ggplot(as.data.frame(tps_field), aes(x = lon, y = lat, fill = tps_pred)) +
  geom_raster(interpolate = F) + theme_bw() +
  scale_fill_gradient2(high = "red", low = "blue", mid = "white", midpoint = 0, limits = range(xx),
                       space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  labs(fill = "", x = "Easting", y = "Northing") + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "tprs_pred.png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()


plot_save <- ggplot(as.data.frame(tps_field), 
                    aes(x = lon, y = lat, fill = as.factor(cover_ind))) +
  geom_raster(interpolate = F) + theme_bw() +
  labs(fill = "", x = "Easting", y = "Northing") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.key.height= unit(2, 'cm'),
        axis.text.x = element_text(margin = unit(c(0.2, 0, 0.2, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.2, 0, 0.2), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"))
png(paste0(outpath, "tprs_pred_95cover.png"), width=600, height=500, pointsize=20, bg = "transparent")
print(plot_save)
dev.off()    
