# compare bias and variance on Gaussian datasets

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/work/tyler/Data/Gaussian/dat/"
library(rerf)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/BiasVariance.R"))

# initialize arrays
ps <- c(2L, 4L)
ns <- c(5L, 10L, 100L, 1000L)
nTrials <- 100L
classifiers <- c("rf", "rerf", "rerfc", "frc", "rr-rf")
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
seed <- 10152017L

for (m in classifiers) {
  
  params[[m]] <- vector(mode = "list", length = length(ps)*length(ns))
  testError[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                          c(length(ns), length(ps), nTrials, 25L))
  OOBError[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                         c(length(ns), length(ps), nTrials, 25L))
  OOBAUC[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                       c(length(ns), length(ps), nTrials, 25L))
  trainTime[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                          c(length(ns), length(ps), nTrials, 25L))
  OOBTime[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                        c(length(ns), length(ps), nTrials, 25L))
  testTime[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                         c(length(ns), length(ps), nTrials, 25L))
  treeStrength[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                             c(length(ns), length(ps), nTrials, 25L))
  treeCorr[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                         c(length(ns), length(ps), nTrials, 25L))
  numNodes[[m]] <- array(as.double(rep(NA, length(ps)*length(ns)*25L*nTrials)),
                         c(length(ns), length(ps), nTrials, 25L))
  bestIdx[[m]] <- array(vector(mode="integer", length=length(ps)*length(ns)*nTrials),
                        c(length(ns), length(ps), nTrials))
  B[[m]] <- matrix(as.double(NA), length(ns), length(ps))
  V[[m]] <- matrix(as.double(NA), length(ns), length(ps))
  SE[[m]] <- matrix(as.double(NA), length(ns), length(ps))
  VE[[m]] <- matrix(as.double(NA), length(ns), length(ps))
  BE[[m]] <- matrix(as.double(NA), length(ns), length(ps))
  
  # loop over number of dimensions
  for (j in 1:length(ps)) {
    p <- ps[j]
    print(paste("p = ", p, sep = ""))
    
    # read in test data for p dimensions
    D <- read.table(paste0(dataPath, "Test/Gaussian_oblique_test_set_p", p, ".dat"), sep = ",")
    Xtest <- as.matrix(D[, 1:ncol(D)-1L])
    Ytest <- as.integer(D[, ncol(D)]) + 1L
    labels <- sort(unique(Ytest))
    ntest = nrow(Xtest)
    
    posteriors <- as.matrix(read.table(paste0(dataPath, "Test/Gaussian_oblique_test_set_posteriors_p", p, ".dat"), sep = ","))
    
    
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
        mtrys <- c(1:p, p^2)
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      }
      sparsity <- (1:min(p-1, 5))/p
    } else if (m == "rerfc") {
      random.matrix <- "continuous"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      }
      sparsity <- (1:min(p-1, 5))/p
    } else if (m == "rerfp" || m == "rerfpr") {
      random.matrix <- "poisson"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      }
      sparsity <- (1:min(ceiling(p/2), 4))
    } else if (m == "frc" || m == "frank") {
      random.matrix <- "frc"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      }
      sparsity <- (2:min(p, 5))
    } else if (m == "frcn") {
      random.matrix <- "frcn"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
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
      
      params[[m]][[(j-1L)*length(ns) + i]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                                                   rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores, seed = seed)
      Yhats <- matrix(0L, ntest, nTrials)
      # loop over trials
      for (trial in 1:nTrials) {
        print(paste("trial = ",trial, sep = ""))
        
        # read in the nth trial of training data for p dimensions 
        D <- read.table(paste0(dataPath, "Train/Gaussian_oblique_train_set_n", ntrain, "_p", p, "_trial", trial, ".dat"), sep = ",")
        Xtrain <- as.matrix(D[, 1:ncol(D) - 1L])
        Ytrain <- as.integer(D[, ncol(D)]) + 1L
        
        # evaluate models
        res <- RerFEval(Xtrain, Ytrain, Xtest, Ytest, params[[m]][[(j-1L)*length(ns) + i]], store.predictions = T)
        
        testError[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$testError
        OOBError[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$oobError
        OOBAUC[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$oobAUC
        trainTime[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$trainTime
        OOBTime[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$oobTime
        testTime[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$testTime
        treeStrength[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$treeStrength
        treeCorr[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$treeCorrelation
        numNodes[[m]][i, j, trial, seq.int(length(mtrys)*length(sparsity))] <- res$numNodes
        bestIdx[[m]][i, j, trial] <- res$best.idx
        Yhats[, trial] <- res$Yhat
      }
      bv.decomp <- BiasVariance(Yhats, posteriors)
      B[[m]][i, j] <- bv.decomp$B
      V[[m]][i, j] <- bv.decomp$V
      SE[[m]][i, j] <- bv.decomp$SE
      VE[[m]][i, j] <- bv.decomp$VE
      BE[[m]][i, j] <- bv.decomp$BE
      save(testError, OOBError, OOBAUC, trainTime, OOBTime, testTime, treeStrength, treeCorr, numNodes, bestIdx, params,
           file = paste(rerfPath, "RandomerForest/R/Results/2017.10.15/Gaussian_oblique_2017_10_15.RData", sep = ""))
    }
  }
}
