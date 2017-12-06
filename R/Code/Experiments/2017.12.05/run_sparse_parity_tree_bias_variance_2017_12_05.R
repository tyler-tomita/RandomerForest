# compare bias and variance of trees (i.e. not the forest) on sparse parity

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/Data/Sparse_parity_bias_variance/"
library(rerf)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/BiasVariance.R"))
load("~/RandomerForest/R/Results/2017.10.31/Sparse_parity_bias_variance_2017_10_31.RData")

# initialize arrays
ns <- c(500L, 1000L, 3000L, 5000L)
nTrials <- 100L
classifiers <- c("rf", "rerf", "frc")
testError <- list()
B <- list()
V <- list()
SE <- list()
VE <- list()
BE <- list()

nTrees <- 200L
num.cores <- 40L

for (m in classifiers) {
  
  testError[[m]] <- matrix(as.double(NA), nTrials*nTrees, length(ns))
  B[[m]] <- rep(as.double(NA), length(ns))
  V[[m]] <- rep(as.double(NA), length(ns))
  SE[[m]] <- rep(as.double(NA), length(ns))
  VE[[m]] <- rep(as.double(NA), length(ns))
  BE[[m]] <- rep(as.double(NA), length(ns))
  
  # read in test data
  D <- read.table(paste0(dataPath, "Test/Sparse_parity_bias_variance_test.csv"), sep = ",")
  p <- ncol(D) - 1L
  Xtest <- as.matrix(D[, 1:p])
  Ytest <- as.integer(D[, p + 1L])
  ntest <- length(Ytest) 
  
  if (m == "rf" || m == "rr-rf" || m == "rr-rfr") {
    random.matrix <- "rf"
  } else if (m == "rerf" || m == "rerfr") {
    random.matrix <- "binary"
  } else if (m == "rerfc") {
    random.matrix <- "continuous"
  } else if (m == "rerfp" || m == "rerfpr") {
    random.matrix <- "poisson"
  } else if (m == "frc" || m == "frank") {
    random.matrix <- "frc"
  } else if (m == "frcn") {
    random.matrix <- "frcn"
  }
  
  if (m == "rr-rf" || m == "rr-rfr") {
    rotate <- T
  } else {
    rotate <- F
  }
  
  # loop over number of train observations
  for (i in 1:length(ns)) {
    ntrain <- ns[i]
    print(paste("n = ", ntrain, sep = ""))
    
    Yhats <- matrix(0L, ntest, nTrials*nTrees)
    # loop over trials
    for (trial in 1:nTrials) {
      print(paste("trial = ",trial, sep = ""))
      
      d <- params[[m]][[i]]$d[(bestIdx[[m]][i, trial] - 1L)%%length(params[[m]][[i]]$d) + 1L]
      sparsity <- params[[m]][[i]]$sparsity[floor((bestIdx[[m]][i, trial] - 1L)/length(params[[m]][[i]]$d)) + 1L]
      
      # read in the nth trial of training data
      D <- read.table(paste0(dataPath, "Train/Sparse_parity_bias_variance_train_n", ntrain, "_trial", trial, ".csv"), sep = ",")
      Xtrain <- as.matrix(D[, 1:p])
      Ytrain <- as.integer(D[, p + 1L])
      
      # evaluate model
      forest <- RerF(Xtrain, Ytrain, trees = nTrees, mat.options = list(p, d, random.matrix, sparsity), rank.transform = params[[m]][[i]]$rank.transform,
                     min.parent = params[[m]][[i]]$min.parent, max.depth = params[[m]][[i]]$max.depth, bagging = 0, store.oob = F,
                     store.impurity = F, replacement = F, stratify = F, num.cores = num.cores,
                     seed = trial, cat.map = NULL, rotate = params[[m]][[i]]$rotate)
      
      Yhats[, ((trial - 1L)*nTrees + 1L):(trial*nTrees)] <- Predict(Xtest, forest, num.cores = num.cores, Xtrain = Xtrain, aggregate.output = F)
      
      testError[[m]][((trial - 1L)*nTrees + 1L):(trial*nTrees), i] <- 1 - apply(Yhats[, ((trial - 1L)*nTrees + 1L):(trial*nTrees)] == Ytest, 2, mean)
    }
    bv.decomp <- BiasVariance(Yhats, cbind(1 - Ytest, Ytest))
    B[[m]][i] <- bv.decomp$B
    V[[m]][i] <- bv.decomp$V
    SE[[m]][i] <- bv.decomp$SE
    VE[[m]][i] <- bv.decomp$VE
    BE[[m]][i] <- bv.decomp$BE
    save(B, V, SE, VE, BE, testError, file = paste0(rerfPath, "RandomerForest/R/Results/2017.12.05/Sparse_parity_tree_bias_variance_2017_12_05.RData"))
  }
}
