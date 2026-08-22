################################################################################
# Section 4.2 Simulation study: Model fitting and uncertainty propagation
# Simulate data using Poisson GLMs
################################################################################
rm(list = ls())

inpath  <- "../outputs/sim1-downscaling/"
outpath <- "../outputs/sim2-pois/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

set.seed(42)


#-------------------------------------------------------------------------------
# Load covariate data
#-------------------------------------------------------------------------------
load(paste0(inpath, "covar_dat.rdata"))
nphi <- dim(fine_xx)[2]
nsill <- dim(fine_xx)[3]


#-------------------------------------------------------------------------------
# Generate poisson data
# Split the data into training and testing set
#-------------------------------------------------------------------------------
true_be <- c(1, 0.5)
nn <- nrow(fine_xx)

nrep <- 100
rep_data <- vector("list", length = nrep)

train_idx <- sample(1:nn, nn / 5)

for (r in 1:nrep) {
  
  fine_yy <- array(NA, dim = c(nn, nphi, nsill))
  fine_zz <- array(NA, dim = c(nn, nphi, nsill))
  
  for (i in 1:nphi) {
    
    for (j in 1:nsill) {
      
      xx <- fine_xx[, i, j]
      yy <- exp(cbind(1, xx) %*% true_be)
      zz <- rpois(nn, yy)
      
      fine_yy[, i, j] <- yy
      fine_zz[, i, j] <- zz
      
    }
    
  }
  
  all <- list(fine_xx = fine_xx, fine_yy = fine_yy, fine_zz = fine_zz)
  train <- list(fine_train_xx = fine_xx[train_idx, , ],
                fine_train_yy = fine_yy[train_idx, , ], 
                fine_train_zz = fine_zz[train_idx, , ])
  test <- list(fine_test_xx = fine_xx[-train_idx, , ], 
               fine_test_yy = fine_yy[-train_idx, , ],
               fine_test_zz = fine_zz[-train_idx, , ])
  
  rep_data[[r]] <- list(all = all, train = train, test = test,
                        train_idx = train_idx, 
                        true_be = true_be)
  
}

save(rep_data, file = paste0(outpath, "pois_rep_data.rdata"))
