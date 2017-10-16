# evaluate RerF on Sparse Parity

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/Data/Sparse_parity/"
library(rerf)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))

# initialize arrays
ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))
nTrials <- 10L
classifiers <- c("rf", "rerf", "rerfc", "rerfp", "frc", "rr-rf")
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

nTrees <- 500L
min.parent <- 2L
max.depth <- 0L
num.cores <- 40L
seed <- 09272017L

for (m in classifiers) {

  params[[m]] <- vector(mode = "list", length = length(ps)*length(ns[[1L]]))
  testError[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                           c(length(ns[[1]]), length(ps), nTrials, 25L))
  OOBError[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                          c(length(ns[[1]]), length(ps), nTrials, 25L))
  OOBAUC[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                        c(length(ns[[1]]), length(ps), nTrials, 25L))
  trainTime[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                           c(length(ns[[1]]), length(ps), nTrials, 25L))
  OOBTime[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                         c(length(ns[[1]]), length(ps), nTrials, 25L))
  testTime[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                          c(length(ns[[1]]), length(ps), nTrials, 25L))
  treeStrength[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                              c(length(ns[[1]]), length(ps), nTrials, 25L))
  treeCorr[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                          c(length(ns[[1]]), length(ps), nTrials, 25L))
  numNodes[[m]] <- array(as.double(rep(NA, length(ps)*length(ns[[1]])*25L*nTrials)),
                          c(length(ns[[1]]), length(ps), nTrials, 25L))
  bestIdx[[m]] <- array(vector(mode="integer", length=length(ps)*length(ns[[1]])*nTrials),
                         c(length(ns[[1]]), length(ps), nTrials))
  
  # loop over number of dimensions
  for (j in 1:length(ps)) {
    p <- ps[j]
    print(paste("p = ", p, sep = ""))
    
    # read in test data for p dimensions
    D <- read.table(paste0(dataPath, "Test/Sparse_parity_test_set_p", p, ".dat"))
    Xtest <- as.matrix(D[, 1:ncol(D)-1L])
    Ytest <- as.integer(D[, ncol(D)]) + 1L
    labels <- sort(unique(Ytest))
    ntest = nrow(Xtest)
    
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
      } else if (p >= 5 && p <= 100) {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
      }
      sparsity <- (1:min(p-1, 5))/p
    } else if (m == "rerfc") {
      random.matrix <- "continuous"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else if (p >= 5 && p <= 100) {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
      }
      sparsity <- (1:min(p, 4))/p
    } else if (m == "rerfp" || m == "rerfpr") {
      random.matrix <- "poisson"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else if (p >= 5 && p <= 100) {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
      }
      sparsity <- (1:min(ceiling(p/2), 4))
    } else if (m == "frc" || m == "frank") {
      random.matrix <- "frc"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else if (p >= 5 && p <= 100) {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
      }
      sparsity <- (2:min(p, 5))
    } else if (m == "frcn") {
      random.matrix <- "frcn"
      if (p < 5) {
        mtrys <- c(1:p, p^2)
      } else if (p >= 5 && p <= 100) {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
      } else {
        mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
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
    for (i in 1:length(ns[[j]])) {
      ntrain <- ns[[j]][i]
      print(paste("n = ", ntrain, sep = ""))
      
      params[[m]][[(j-1L)*length(ns[[j]]) + i]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                           rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores, seed = seed)
      
      # loop over trials
      for (trial in 1:nTrials) {
        print(paste("trial = ",trial, sep = ""))
        
        # read in the nth trial of training data for p dimensions 
        D <- read.table(paste0(dataPath, "Train/Sparse_parity_train_set_n", ntrain, "_p", p, "_trial", trial, ".dat"))
        Xtrain <- as.matrix(D[, 1:ncol(D) - 1L])
        Ytrain <- as.integer(D[, ncol(D)]) + 1L
        ntrain <- nrow(Xtrain)
        
        # evaluate models
        res <- RerFEval(Xtrain, Ytrain, Xtest, Ytest, params[[m]][[(j-1L)*length(ns[[j]]) + i]])
        
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
        
        save(testError, OOBError, OOBAUC, trainTime, OOBTime, testTime, treeStrength, treeCorr, numNodes, bestIdx, params, file = paste(rerfPath, "RandomerForest/R/Results/2017.09.27/Sparse_parity_2017_09_27.RData", sep = ""))
      }
    }
  }
}
