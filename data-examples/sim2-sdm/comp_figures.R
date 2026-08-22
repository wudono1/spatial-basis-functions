################################################################################
# Section 4.2 Simulation study: Model fitting and uncertainty propagation
# Generate figures (Fig. 5 and Fig. 6 in the main text)
################################################################################
rm(list = ls())

inpath1 <- "../outputs/sim1-downscaling/"
inpath2 <- "../outputs/sim2-pois/"
outpath <- "../outputs/sim2-pois/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
load(paste0(inpath1, "covar_dat.rdata"))
load(paste0(inpath2, "pois_eval_out.rdata"))

nrep <- length(pois_mods_metrics)
nphi <- length(pois_mods_metrics[[1]])
nsill <- length(pois_mods_metrics[[1]][[1]])

#-------------------------------------------------------------------------------
# Bias
#-------------------------------------------------------------------------------
orig_be_bias_all <- array(NA, dim = c(nphi, nsill, nrep))
plugin_be_bias_all <- array(NA, dim = c(nphi, nsill, nrep))
ppd_be_bias_all <- array(NA, dim = c(nphi, nsill, nrep))
bks_be_bias_all <- array(NA, dim = c(nphi, nsill, nrep))

for (r in 1:nrep) {
  for (i in 1:nphi) {
    for (j in 1:nsill) {
      orig_be_bias_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$orig$orig_be_bias
      plugin_be_bias_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$plugin$plugin_be_bias
      ppd_be_bias_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$ppd$ppd_be_bias
      bks_be_bias_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$bks$bks_be_bias
    }
  }
}

orig_be_bias <- apply(orig_be_bias_all, c(1,2), mean)
plugin_be_bias <- apply(plugin_be_bias_all, c(1,2), mean)
ppd_be_bias <- apply(ppd_be_bias_all, c(1,2), mean)
bks_be_bias <- apply(bks_be_bias_all, c(1,2), mean)

orig_be_bias_by_phi <- rowMeans(orig_be_bias)
plugin_be_bias_by_phi <- rowMeans(plugin_be_bias)
ppd_be_bias_by_phi <- rowMeans(ppd_be_bias)
bks_be_bias_by_phi <- rowMeans(bks_be_bias)

orig_be_bias_by_sill <- colMeans(orig_be_bias)
plugin_be_bias_by_sill <- colMeans(plugin_be_bias)
ppd_be_bias_by_sill <- colMeans(ppd_be_bias)
bks_be_bias_by_sill <- colMeans(bks_be_bias)

be_bias_total <- c(orig_be_bias_by_phi, orig_be_bias_by_sill,
                   plugin_be_bias_by_phi, plugin_be_bias_by_sill,
                   ppd_be_bias_by_phi, ppd_be_bias_by_sill,
                   bks_be_bias_by_phi, bks_be_bias_by_sill)

png(paste0(outpath, "be_bias_by_phi.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 5.1, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(range_vals), 
     ylim = range(be_bias_total), cex.axis = 1.2,
     ylab = expression(paste("Bias of ", hat(beta)[1])),
     xlab = expression(psi[0]), cex.lab = 1.4)
axis(side = 1, at = range_vals, labels = range_vals, cex.axis = 1.2)
lines(range_vals, orig_be_bias_by_phi, lwd = 5, lty = 1, col = 'black')
lines(range_vals, bks_be_bias_by_phi, lwd = 5, lty = 4, col = 'blue')
lines(range_vals, plugin_be_bias_by_phi, lwd = 5, lty = 2, col = 'green')
lines(range_vals, ppd_be_bias_by_phi, lwd = 5, lty = 3, col = 'red')
dev.off()

png(paste0(outpath, "be_bias_by_sill.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 5.1, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(sill_vals), 
     ylim = range(be_bias_total), cex.axis = 1.2,
     ylab = expression(paste("Bias of ", hat(beta)[1])),
     xlab = expression(sigma[0]^2), cex.lab = 1.4)
axis(side = 1, at = sill_vals, labels = sill_vals, cex.axis = 1.2)
lines(sill_vals, orig_be_bias_by_sill, lwd = 5, lty = 1, col = 'black')
lines(sill_vals, bks_be_bias_by_sill, lwd = 5, lty = 4, col = 'blue')
lines(sill_vals, plugin_be_bias_by_sill, lwd = 5, lty = 2, col = 'green')
lines(sill_vals, ppd_be_bias_by_sill, lwd = 5, lty = 3, col = 'red')
dev.off()


#-------------------------------------------------------------------------------
# 95% CP for beta_1
#-------------------------------------------------------------------------------
orig_be_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
plugin_be_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
ppd_be_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
bks_be_cp_all <- array(NA, dim = c(nphi, nsill, nrep))

for (r in 1:nrep) {
  for (i in 1:nphi) {
    for (j in 1:nsill) {
      orig_be_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$orig$orig_be_ci_ind
      plugin_be_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$plugin$plugin_be_ci_ind
      ppd_be_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$ppd$ppd_be_ci_ind
      bks_be_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$bks$bks_be_ci_ind
    }
  }
}


orig_be_cp <- apply(orig_be_cp_all, c(1,2), mean)
plugin_be_cp <- apply(plugin_be_cp_all, c(1,2), mean)
ppd_be_cp <- apply(ppd_be_cp_all, c(1,2), mean)
bks_be_cp <- apply(bks_be_cp_all, c(1,2), mean)

orig_be_cp_by_phi <- rowMeans(orig_be_cp)
plugin_be_cp_by_phi <- rowMeans(plugin_be_cp)
ppd_be_cp_by_phi <- rowMeans(ppd_be_cp)
bks_be_cp_by_phi <- rowMeans(bks_be_cp)

orig_be_cp_by_sill <- colMeans(orig_be_cp)
plugin_be_cp_by_sill <- colMeans(plugin_be_cp)
ppd_be_cp_by_sill <- colMeans(ppd_be_cp)
bks_be_cp_by_sill <- colMeans(bks_be_cp)

be_cp_total <- c(orig_be_cp_by_phi, orig_be_cp_by_sill,
                   plugin_be_cp_by_phi, plugin_be_cp_by_sill,
                   ppd_be_cp_by_phi, ppd_be_cp_by_sill,
                   bks_be_cp_by_phi, bks_be_cp_by_sill)

png(paste0(outpath, "be_cp_by_phi.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 5.1, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(range_vals), ylim = range(be_cp_total), cex.axis = 1.2,
     ylab = expression(paste("95% CP for ", beta[1])),
     xlab = expression(psi[0]), cex.lab = 1.4)
axis(side = 1, at = range_vals, labels = range_vals, cex.axis = 1.2)
lines(range_vals, orig_be_cp_by_phi, lwd = 5, lty = 1, col = 'black')
lines(range_vals, bks_be_cp_by_phi, lwd = 5, lty = 4, col = 'blue')
lines(range_vals, plugin_be_cp_by_phi, lwd = 5, lty = 2, col = 'green')
lines(range_vals, ppd_be_cp_by_phi, lwd = 5, lty = 3, col = 'red')
dev.off()

png(paste0(outpath, "be_cp_by_sill.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 5.1, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(sill_vals), ylim = range(be_cp_total), cex.axis = 1.2,
     ylab = expression(paste("95% CP for ", beta[1])),
     xlab = expression(sigma[0]^2), cex.lab = 1.4)
axis(side = 1, at = sill_vals, labels = sill_vals, cex.axis = 1.2)
lines(sill_vals, orig_be_cp_by_sill, lwd = 5, lty = 1, col = 'black')
lines(sill_vals, bks_be_cp_by_sill, lwd = 5, lty = 4, col = 'blue')
lines(sill_vals, plugin_be_cp_by_sill, lwd = 5, lty = 2, col = 'green')
lines(sill_vals, ppd_be_cp_by_sill, lwd = 5, lty = 3, col = 'red')
dev.off()


#-------------------------------------------------------------------------------
# RMSPE for linear predictors
#-------------------------------------------------------------------------------
orig_pred_rmse_all <- array(NA, dim = c(nphi, nsill, nrep))
plugin_pred_rmse_all <- array(NA, dim = c(nphi, nsill, nrep))
ppd_pred_rmse_all <- array(NA, dim = c(nphi, nsill, nrep))
bks_pred_rmse_all <- array(NA, dim = c(nphi, nsill, nrep))
tprs_pred_rmse_all <- array(NA, dim = c(nphi, nsill, nrep))

for (r in 1:nrep) {
  for (i in 1:nphi) {
    for (j in 1:nsill) {
      orig_pred_rmse_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$orig$orig_pred_rmse
      plugin_pred_rmse_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$plugin$plugin_pred_rmse
      ppd_pred_rmse_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$ppd$ppd_pred_rmse
      bks_pred_rmse_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$bks$bks_pred_rmse
    }
  }
}


orig_pred_rmse <- apply(orig_pred_rmse_all, c(1,2), mean)
plugin_pred_rmse <- apply(plugin_pred_rmse_all, c(1,2), mean)
ppd_pred_rmse <- apply(ppd_pred_rmse_all, c(1,2), mean)
bks_pred_rmse <- apply(bks_pred_rmse_all, c(1,2), mean)

orig_pred_rmse_by_phi <- rowMeans(orig_pred_rmse)
plugin_pred_rmse_by_phi <- rowMeans(plugin_pred_rmse)
ppd_pred_rmse_by_phi <- rowMeans(ppd_pred_rmse)
bks_pred_rmse_by_phi <- rowMeans(bks_pred_rmse)

orig_pred_rmse_by_sill <- colMeans(orig_pred_rmse)
plugin_pred_rmse_by_sill <- colMeans(plugin_pred_rmse)
ppd_pred_rmse_by_sill <- colMeans(ppd_pred_rmse)
bks_pred_rmse_by_sill <- colMeans(bks_pred_rmse)

pred_rmse_total <- c(orig_pred_rmse_by_phi, orig_pred_rmse_by_sill,
                   plugin_pred_rmse_by_phi, plugin_pred_rmse_by_sill,
                   ppd_pred_rmse_by_phi, ppd_pred_rmse_by_sill,
                   bks_pred_rmse_by_phi, bks_pred_rmse_by_sill)

png(paste0(outpath, "pred_rmse_by_phi.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 4.5, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(range_vals), ylim = range(pred_rmse_total), cex.axis = 1.2,
     ylab = "RMSPE for linear predictor",
     xlab = expression(psi[0]), cex.lab = 1.4)
axis(side = 1, at = range_vals, labels = range_vals, cex.axis = 1.2)
lines(range_vals, orig_pred_rmse_by_phi, lwd = 5, lty = 1, col = 'black')
lines(range_vals, bks_pred_rmse_by_phi, lwd = 5, lty = 4, col = 'blue')
lines(range_vals, plugin_pred_rmse_by_phi, lwd = 5, lty = 2, col = 'green')
lines(range_vals, ppd_pred_rmse_by_phi, lwd = 5, lty = 3, col = 'red')
dev.off()

png(paste0(outpath, "pred_rmse_by_sill.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 4.5, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(sill_vals), ylim = range(pred_rmse_total), cex.axis = 1.2,
     ylab = "RMSPE for linear predictor",
     xlab = expression(sigma[0]^2), cex.lab = 1.4)
axis(side = 1, at = sill_vals, labels = sill_vals, cex.axis = 1.2)
lines(sill_vals, orig_pred_rmse_by_sill, lwd = 5, lty = 1, col = 'black')
lines(sill_vals, bks_pred_rmse_by_sill, lwd = 5, lty = 4, col = 'blue')
lines(sill_vals, plugin_pred_rmse_by_sill, lwd = 5, lty = 2, col = 'green')
lines(sill_vals, ppd_pred_rmse_by_sill, lwd = 5, lty = 3, col = 'red')
dev.off()


#-------------------------------------------------------------------------------
# 95% CP for linear predictors
#-------------------------------------------------------------------------------
orig_pred_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
plugin_pred_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
ppd_pred_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
bks_pred_cp_all <- array(NA, dim = c(nphi, nsill, nrep))
tprs_pred_cp_all <- array(NA, dim = c(nphi, nsill, nrep))

for (r in 1:nrep) {
  for (i in 1:nphi) {
    for (j in 1:nsill) {
      orig_pred_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$orig$orig_pred_ci_cover
      plugin_pred_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$plugin$plugin_pred_ci_cover
      ppd_pred_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$ppd$ppd_pred_ci_cover
      bks_pred_cp_all[i, j, r] <- pois_mods_metrics[[r]][[i]][[j]]$bks$bks_pred_ci_cover
    }
  }
}

orig_pred_cp <- apply(orig_pred_cp_all, c(1,2), mean)
plugin_pred_cp <- apply(plugin_pred_cp_all, c(1,2), mean)
ppd_pred_cp <- apply(ppd_pred_cp_all, c(1,2), mean)
bks_pred_cp <- apply(bks_pred_cp_all, c(1,2), mean)

orig_pred_cp_by_phi <- rowMeans(orig_pred_cp)
plugin_pred_cp_by_phi <- rowMeans(plugin_pred_cp)
ppd_pred_cp_by_phi <- rowMeans(ppd_pred_cp)
bks_pred_cp_by_phi <- rowMeans(bks_pred_cp)

orig_pred_cp_by_sill <- colMeans(orig_pred_cp)
plugin_pred_cp_by_sill <- colMeans(plugin_pred_cp)
ppd_pred_cp_by_sill <- colMeans(ppd_pred_cp)
bks_pred_cp_by_sill <- colMeans(bks_pred_cp)

pred_cp_total <- c(orig_pred_cp_by_phi, orig_pred_cp_by_sill,
                 plugin_pred_cp_by_phi, plugin_pred_cp_by_sill,
                 ppd_pred_cp_by_phi, ppd_pred_cp_by_sill,
                 bks_pred_cp_by_phi, bks_pred_cp_by_sill)

png(paste0(outpath, "pred_cp_by_phi.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 4.5, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(range_vals), ylim = range(pred_cp_total), cex.axis = 1.2,
     ylab = "95% CP for linear predictor",
     xlab = expression(psi[0]), cex.lab = 1.4)
axis(side = 1, at = range_vals, labels = range_vals, cex.axis = 1.2)
lines(range_vals, orig_pred_cp_by_phi, lwd = 5, lty = 1, col = 'black')
lines(range_vals, bks_pred_cp_by_phi, lwd = 5, lty = 4, col = 'blue')
lines(range_vals, plugin_pred_cp_by_phi, lwd = 5, lty = 2, col = 'green')
lines(range_vals, ppd_pred_cp_by_phi, lwd = 5, lty = 3, col = 'red')
dev.off()

png(paste0(outpath, "pred_cp_by_sill.png"), width=800, height=500, pointsize=20, bg = "transparent")
par(mar = c(4.1, 4.5, 1, 1))
plot(NULL, lwd = 5, type = 'l', xaxt = 'n',
     xlim = range(sill_vals), ylim = range(pred_cp_total), cex.axis = 1.2,
     ylab = "95% CP for linear predictor",
     xlab = expression(sigma[0]^2), cex.lab = 1.4)
axis(side = 1, at = sill_vals, labels = sill_vals, cex.axis = 1.2)
lines(sill_vals, orig_pred_cp_by_sill, lwd = 5, lty = 1, col = 'black')
lines(sill_vals, bks_pred_cp_by_sill, lwd = 5, lty = 4, col = 'blue')
lines(sill_vals, plugin_pred_cp_by_sill, lwd = 5, lty = 2, col = 'green')
lines(sill_vals, ppd_pred_cp_by_sill, lwd = 5, lty = 3, col = 'red')
dev.off()


tbl_est_cp <- c(mean(orig_be_cp), mean(bks_be_cp), mean(plugin_be_cp), mean(ppd_be_cp))
tbl_pred_cp <- c(mean(orig_pred_cp), mean(bks_pred_cp), mean(plugin_pred_cp), mean(ppd_pred_cp))
tbl_cp <- rbind(tbl_est_cp, tbl_pred_cp)
colnames(tbl_cp) <- c("GLM-oracle", "GLM-Berkson", "GLM-plugin", "GLM-ensemble")
rownames(tbl_cp) <- c("CP for $\\beta_1$", "CP for linear predictor")
knitr::kable(tbl_cp, format = "latex", digits = 3, escape = FALSE)
