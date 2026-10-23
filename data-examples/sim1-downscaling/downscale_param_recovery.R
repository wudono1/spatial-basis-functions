################################################################################
# Step D: parameter recovery check for the basis-function GPB model
# ------------------------------------------------------------------------------
# The simulation truth is known exactly. In simulate_covariate_data.R the field
# is generated as  xx = Phi %*% alpha,  alpha ~ rgp(sill = sill_vals[j],
# range = range_vals[i], cov.mod = "powexp", smooth = 1, nugget = 0), i.e.
# cov(alpha) = sill * exp(-d / range). The fitted model says cov(alpha) =
# sigmasq * K_psi with K_psi = exp(-d / psi). So:
#     psi     truth = range_vals[i]
#     sigmasq truth = sill_vals[j]
#     gamma0  truth = 0            (alpha has mean zero; no mean term generated)
################################################################################
rm(list = ls())

inpath  <- "../outputs/sim1-downscaling/"
outpath <- "../outputs/sim1-downscaling/"

library(ggplot2)

load(paste0(inpath, "sim_gpbau_new_fit.RData"))
# gives: out_list, range_vals, sill_vals, gpb_runtime, findint_runtime

nrange <- length(range_vals)
nsill  <- length(sill_vals)

#-------------------------------------------------------------------------------
# Build a long-format recovery table: one row per (i, j, parameter)
#-------------------------------------------------------------------------------
res <- list()
k <- 0

for (i in seq_len(nrange)) {
  for (j in seq_len(nsill)) {
    
    obj <- out_list[[i]][[j]]
    
    # thinning() stores draws in post_sams_thin; fall back if not thinned
    ps <- if (!is.null(obj$post_sams_thin)) obj$post_sams_thin else obj$post_sams
    
    # NOTE: the list key is still "phi" for backward compatibility, but the
    # quantity it holds is the range parameter psi.
    chains <- list(psi     = as.numeric(ps$phi),
                   sigmasq = as.numeric(ps$sigmasq),
                   gamma0  = as.numeric(ps$beta))
    
    
    # true values
    truths <- list(psi     = range_vals[i],
                   sigmasq = sill_vals[j],
                   gamma0  = 0)
    
    for (par in names(chains)) {
      
      chain <- chains[[par]]
      truth <- truths[[par]]
      qq <- quantile(chain, c(0.025, 0.975))
      
      # effective sample size if coda is available (MH chains are autocorrelated)
      ess <- if (requireNamespace("coda", quietly = TRUE)) {
        as.numeric(coda::effectiveSize(chain))
      } else NA_real_
      
      k <- k + 1
      res[[k]] <- data.frame(
        i = i, j = j,
        parameter   = par,
        truth       = truth,
        post_mean   = mean(chain),
        post_median = median(chain),
        lower95     = as.numeric(qq[1]),
        upper95     = as.numeric(qq[2]),
        covered     = truth >= qq[1] && truth <= qq[2],
        bias        = mean(chain) - truth,
        rel_bias    = if (truth != 0) (mean(chain) - truth) / truth else NA_real_,
        ci_width    = as.numeric(qq[2] - qq[1]),
        ess         = ess,
        accept_rate = mean(obj$mh_data),
        n_draws     = length(chain),
        row.names   = NULL
      )
    }
  }
}

recov <- do.call(rbind, res)

#-------------------------------------------------------------------------------
# Headline numbers
#-------------------------------------------------------------------------------
cat("\n=== 95% credible interval coverage (target: ~0.95) ===\n")
print(aggregate(covered ~ parameter, data = recov, FUN = mean))

cat("\n=== mean relative bias of posterior mean ===\n")
print(aggregate(rel_bias ~ parameter, data = recov,
                FUN = function(x) mean(x, na.rm = TRUE)))

cat("\n=== MH acceptance rate (want roughly 0.20-0.40) ===\n")
acc <- unique(recov[, c("i", "j", "accept_rate")])
print(summary(acc$accept_rate))

cat("\n=== effective sample size ===\n")
print(aggregate(ess ~ parameter, data = recov,
                FUN = function(x) mean(x, na.rm = TRUE)))

# any fits where BOTH psi and sigmasq missed
miss <- recov[!recov$covered & recov$parameter %in% c("psi", "sigmasq"), ]
if (nrow(miss) > 0) {
  cat("\n=== non-covering fits ===\n")
  print(miss[, c("i", "j", "parameter", "truth", "post_mean",
                 "lower95", "upper95")])
}

#-------------------------------------------------------------------------------
# Caterpillar plots: posterior interval vs truth, one panel per parameter
#-------------------------------------------------------------------------------
recov$label <- paste0("r", recov$i, "-s", recov$j)

p_cat <- ggplot(subset(recov, parameter %in% c("psi", "sigmasq")),
                aes(x = label)) +
  geom_linerange(aes(ymin = lower95, ymax = upper95, colour = covered)) +
  geom_point(aes(y = post_mean, colour = covered), size = 1) +
  geom_point(aes(y = truth), shape = 4, size = 2) +   # x = truth
  facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("TRUE" = "grey30", "FALSE" = "red")) +
  labs(x = "range-sill combination", y = "value",
       title = "Posterior mean and 95% CI vs truth (x = truth)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))

png(paste0(outpath, "recovery_intervals.png"),
    width = 1000, height = 700, pointsize = 16)
print(p_cat)
dev.off()

#-------------------------------------------------------------------------------
# Estimated vs true scatter (should sit on the 1:1 line)
#-------------------------------------------------------------------------------
p_scat <- ggplot(subset(recov, parameter %in% c("psi", "sigmasq")),
                 aes(x = truth, y = post_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), alpha = 0.4) +
  geom_point() +
  facet_wrap(~ parameter, scales = "free") +
  labs(x = "true value", y = "posterior mean",
       title = "Parameter recovery (dashed = 1:1)") +
  theme_bw()

png(paste0(outpath, "recovery_scatter.png"),
    width = 1000, height = 500, pointsize = 16)
print(p_scat)
dev.off()

#-------------------------------------------------------------------------------
# Trace plots for a few chains -- convergence sanity check
#-------------------------------------------------------------------------------
png(paste0(outpath, "recovery_traces.png"),
    width = 1000, height = 700, pointsize = 16)
par(mfrow = c(3, 2))
for (i in c(1, ceiling(nrange / 2), nrange)) {
  j <- ceiling(nsill / 2)
  ps <- out_list[[i]][[j]]$post_sams_thin
  plot(as.numeric(ps$phi), type = "l",
       main = paste0("psi trace: range=", range_vals[i], ", sill=", sill_vals[j]),
       ylab = "psi", xlab = "draw")
  abline(h = range_vals[i], col = "red", lwd = 2)
  plot(as.numeric(ps$sigmasq), type = "l",
       main = paste0("sigmasq trace"), ylab = "sigmasq", xlab = "draw")
  abline(h = sill_vals[j], col = "red", lwd = 2)
}
dev.off()

write.csv(recov, paste0(outpath, "recovery_table.csv"), row.names = FALSE)