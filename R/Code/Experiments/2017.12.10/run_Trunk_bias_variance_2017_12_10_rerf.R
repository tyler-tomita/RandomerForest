# compare bias and variance on Trunk dataset

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/work/tyler/"
dataPath <- "~/work/tyler/Data/Trunk_bias_variance/"
library(rerf)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/BiasVariance.R"))

# initialize arrays
ns <- c(10L, 100L, 1000L, 10000L)
nTrials <- 100L
classifiers <- "rerf"
testError <- list()
OOBError <- list()
OOBAUC <- list()
trainTime <- list()
OOBTime <- list()
testTime <- list()
treeStrength <- list()
treeCorr <- list()
numNodes <- list()
bestIdx <- list()
params <- list()
B <- list()
V <- list()
SE <- list()
VE <- list()
BE <- list()

nTrees <- 500L
min.parent <- 2L
max.depth <- 0L
num.cores <- 24L
seed <- 12102017L

for (m in classifiers) {
  
  params[[m]] <- vector(mode = "list", length = length(ns))
  testError[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                          c(length(ns), nTrials, 25L))
  OOBError[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                         c(length(ns), nTrials, 25L))
  OOBAUC[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                       c(length(ns), nTrials, 25L))
  trainTime[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                          c(length(ns), nTrials, 25L))
  OOBTime[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                        c(length(ns), nTrials, 25L))
  testTime[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                         c(length(ns), nTrials, 25L))
  treeStrength[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                             c(length(ns), nTrials, 25L))
  treeCorr[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                         c(length(ns), nTrials, 25L))
  numNodes[[m]] <- array(as.double(rep(NA, length(ns)*25L*nTrials)),
                         c(length(ns), nTrials, 25L))
  bestIdx[[m]] <- matrix(0L, length(ns), nTrials)
  
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
  
  if (m == "rf" || m == "rr-rf" || m == "rr-rfr") {
    random.matrix <- "rf"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- 1/p # this parameter doesn't actually matter for RF
  } else if (m == "rerf" || m == "rerfr") {
    random.matrix <- "binary"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (1:min(p-1, 5))/p
  } else if (m == "rerfc") {
    random.matrix <- "continuous"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (1:min(p-1, 5))/p
  } else if (m == "rerfp" || m == "rerfpr") {
    random.matrix <- "poisson"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (1:min(ceiling(p/2), 4))
  } else if (m == "frc" || m == "frank") {
    random.matrix <- "frc"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (2:min(p, 5))
  } else if (m == "frcn") {
    random.matrix <- "frcn"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (2:min(p, 5))
  }
  
  if (m == "rr-rf" || m == "rr-rfr") {
    rotate <- T
  } else {
    rotate <- F
  }
  
  if (m == "rerfr" || m == "rerfpr" || m == "rerfcr" || m == "frank" || m == "rr-rfr") {
    rank.transform <- T
  } else {
    rank.transform <- F
  }
  
  # loop over number of train observations
  for (i in 1:length(ns)) {
    ntrain <- ns[i]
    print(paste("n = ", ntrain, sep = ""))
    
    params[[m]][[i]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                             rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores, seed = seed)
    Yhats <- matrix(0L, ntest, nTrials)
    # loop over trials
    for (trial in 1:nTrials) {
      print(paste("trial = ",trial, sep = ""))
      
      # read in the nth trial of training data
      D <- read.table(paste0(dataPath, "Train/Trunk_bias_variance_train_n", ntrain, "_trial", trial, ".csv"), sep = ",")
      Xtrain <- as.matrix(D[, 1:p])
      Ytrain <- as.integer(D[, p + 1L])
      
      # evaluate models
      res <- RerFEval(Xtrain, Ytrain, Xtest, Ytest, params[[m]][[i]], store.predictions = T)
      
      testError[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$testError
      OOBError[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$oobError
      OOBAUC[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$oobAUC
      trainTime[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$trainTime
      OOBTime[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$oobTime
      testTime[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$testTime
      treeStrength[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$treeStrength
      treeCorr[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$treeCorrelation
      numNodes[[m]][i, trial, seq.int(length(mtrys)*length(sparsity))] <- res$numNodes
      bestIdx[[m]][i, trial] <- res$best.idx
      Yhats[, trial] <- res$Yhat
    }
    bv.decomp <- BiasVariance(Yhats, posteriors)
    B[[m]][i] <- bv.decomp$B
    V[[m]][i] <- bv.decomp$V
    SE[[m]][i] <- bv.decomp$SE
    VE[[m]][i] <- bv.decomp$VE
    BE[[m]][i] <- bv.decomp$BE
    save(B, V, SE, VE, BE, testError, OOBError, OOBAUC, trainTime, OOBTime, testTime, treeStrength, treeCorr, numNodes, bestIdx, params,
         file = paste(rerfPath, "RandomerForest/R/Results/2017.12.10/Trunk_bias_variance_2017_12_10_rerf.RData", sep = ""))
  }
}
