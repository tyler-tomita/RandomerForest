# compute bias and variance of mean difference classifier on Trunk

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/R/Data/Trunk_bias_variance/"
library(rerf)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/BiasVariance.R"))

# initialize arrays
ns <- c(10L, 100L, 1000L, 10000L)
nTrials <- 100L
testError <- list()
B <- list()
V <- list()
SE <- list()
VE <- list()
BE <- list()

classifiers <- "md"

for (m in classifiers) {
  
  testError[[m]] <- matrix(as.double(rep(NA, length(ns)*nTrials)), length(ns), nTrials)
  B[[m]] <- rep(as.double(NA), length(ns))
  V[[m]] <- rep(as.double(NA), length(ns))
  SE[[m]] <- rep(as.double(NA), length(ns))
  VE[[m]] <- rep(as.double(NA), length(ns))
  BE[[m]] <- rep(as.double(NA), length(ns))
  
  # read in test data
  D <- read.table(paste0(dataPath, "Test/Trunk_bias_variance_test.csv"), sep = ",")
  p <- ncol(D) - 1L
  Xtest <- as.matrix(D[, 1:p])
  Ytest <- as.integer(D[, p + 1L])
  ntest <- length(Ytest) 
  posteriors <- as.matrix(read.table(paste0(dataPath, "Test/Trunk_bias_variance_test_posteriors.csv"), sep = ","))
  
  # loop over number of train observations
  for (i in 1:length(ns)) {
    ntrain <- ns[i]
    print(paste("n = ", ntrain, sep = ""))
    
    Yhats <- matrix(0L, ntest, nTrials)
    
    # loop over trials
    for (trial in 1:nTrials) {
      print(paste("trial = ",trial, sep = ""))
      
      # read in the nth trial of training data
      D <- read.table(paste0(dataPath, "Train/Trunk_bias_variance_train_n", ntrain, "_trial", trial, ".csv"), sep = ",")
      Xtrain <- as.matrix(D[, 1:p])
      Ytrain <- as.integer(D[, p + 1L])
      
      # learn parameters
      muhat <- sapply(0:1, function(k) apply(Xtrain[Ytrain == k, ], 2, mean))
      delta <- apply(muhat, 1, diff)
      mubar <- apply(muhat, 1, mean)
      
      # make predictions
      Yhats[t(t(Xtest) - mubar)%*%delta > 0, trial] <- 1L
      testError[[m]][i, trial] <- mean(Yhats[, trial] != Ytest)
    }
    bv.decomp <- BiasVariance(Yhats, posteriors)
    B[[m]][i] <- bv.decomp$B
    V[[m]][i] <- bv.decomp$V
    SE[[m]][i] <- bv.decomp$SE
    VE[[m]][i] <- bv.decomp$VE
    BE[[m]][i] <- bv.decomp$BE
    save(B, V, SE, VE, BE, testError, file = paste(rerfPath, "RandomerForest/R/Results/2017.12.10/Trunk_bias_variance_2017_12_10_md.RData", sep = ""))
  }
}
