# evaluate classifiers on image stripes

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/R/Data/Image_stripes/"
library(rerf)
library(AUC)
library(dummies)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))

#classifiers <- c("rf", "rerf", "rerfr", "rerfp", "rerfpr", "frc", "frank", "rr-rf", "rr-rfr")
classifiers <- c("rf", "rerf", "strerf", "control")
nCl <- length(classifiers)

nTrees <- 500L
min.parent <- 2L
max.depth <- 0L
num.cores <- 1L
seed <- 01022018L

testError <- list()
testAUC <- list()
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

ns <- c(10L, 20L, 50L)
p <- 20L*20L
num.trials <- 10L


testError <- vector(mode = "list", length = nCl)
names(testError) <- classifiers
testAUC <- vector(mode = "list", length = nCl)
names(testAUC) <- classifiers
OOBError <- vector(mode = "list", length = nCl)
names(OOBError) <- classifiers
OOBAUC <- vector(mode = "list", length = nCl)
names(OOBAUC) <- classifiers
trainTime <- vector(mode = "list", length = nCl)
names(trainTime) <- classifiers
OOBTime <- vector(mode = "list", length = nCl)
names(OOBTime) <- classifiers
testTime <- vector(mode = "list", length = nCl)
names(testTime) <- classifiers
treeStrength <- vector(mode = "list", length = nCl)
names(treeStrength) <- classifiers
treeCorr <- vector(mode = "list", length = nCl)
names(treeCorr) <- classifiers
numNodes <- vector(mode = "list", length = nCl)
names(numNodes) <- classifiers
bestIdx <- vector(mode = "list", length = nCl)
names(bestIdx) <- classifiers
params <- vector(mode = "list", length = nCl)
names(params) <- classifiers

D <- as.matrix(read.table(file = paste0(dataPath, "Test/Image_stripes_test_set.csv"), header = F, sep = ",", quote = "", row.names = NULL))
Xtest <- D[, 1:p]
Ytest <- as.integer(D[, p + 1L])

for (m in classifiers) {
  print(m)
  
  if (m == "rf") {
    random.matrix <- "rf"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    sparsity <- 1/p # this parameter doesn't actually matter for RF
  } else if (m == "rerf") {
    random.matrix <- "binary"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    sparsity <- (1:min(p-1, 5))/p
  } else if (m == "strerf") {
    random.matrix <- "image-patch"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    ih <- 20L
    iw <- 20L
    patch.min <- 1L
    patch.max <- 20L
    sparsity <- 1/p
  } else if (m == "control") {
    random.matrix <- "image-control"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    ih <- 20L
    iw <- 20L
    patch.min <- 1L
    patch.max <- 20L
    sparsity <- 1/p
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
  
  catMap <- NULL
  
  if ((m == "strerf") || (m == "control")) {
    params[[m]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                        rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores,
                        seed = seed, cat.map = catMap, iw = iw, ih = ih, patch.min = patch.min, patch.max = patch.max)
  } else {
    params[[m]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                        rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores,
                        seed = seed, cat.map = catMap)
  }
  
  testError[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                      nrow = num.trials, ncol = length(ns))
  testAUC[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                    nrow = num.trials, ncol = length(ns))
  OOBError[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                     nrow = num.trials, ncol = length(ns))
  OOBAUC[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                   nrow = num.trials, ncol = length(ns))
  trainTime[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                      nrow = num.trials, ncol = length(ns))
  OOBTime[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                    nrow = num.trials, ncol = length(ns))
  testTime[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                     nrow = num.trials, ncol = length(ns))
  treeStrength[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                         nrow = num.trials, ncol = length(ns))
  treeCorr[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                     nrow = num.trials, ncol = length(ns))
  numNodes[[m]] <- matrix(as.double(rep(NA, num.trials*length(ns))),
                                     nrow = num.trials, ncol = length(ns))
  bestIdx[[m]] <- matrix(as.integer(rep(NA, num.trials*length(ns))),
                         nrow = num.trials, ncol = length(ns))
  
  # loop over number of training samples
  for (i in seq_along(ns)) {
    ntrain <- ns[i]
    # loop over trials
    for (trial in seq.int(num.trials)) {
      print(paste0("trial ", trial))
      
      D <- as.matrix(read.table(file = paste0(dataPath, "Train/Image_stripes_train_set_n", ntrain, "_trial", trial, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
      Xtrain <- D[, 1:p]
      Ytrain <- as.integer(D[, p + 1L])
      
      # evaluate models
      res <- RerFEval(Xtrain, Ytrain, Xtest, Ytest, params[[m]])
      
      testError[[m]][trial, i] <- res$testError[res$best.idx]
      testAUC[[m]][trial, i] <- res$testAUC[res$best.idx]
      OOBError[[m]][trial, i] <- res$oobError[res$best.idx]
      OOBAUC[[m]][trial, i] <- res$oobAUC[res$best.idx]
      trainTime[[m]][trial, i] <- res$trainTime[res$best.idx]
      OOBTime[[m]][trial, i] <- res$oobTime[res$best.idx]
      testTime[[m]][trial, i] <- res$testTime[res$best.idx]
      treeStrength[[m]][trial, i] <- res$treeStrength[res$best.idx]
      treeCorr[[m]][trial, i] <- res$treeCorrelation[res$best.idx]
      numNodes[[m]][trial, i] <- res$numNodes[res$best.idx]
      bestIdx[[m]][trial, i] <- res$best.idx
      
      save(testError, testAUC, OOBError, OOBAUC, trainTime, OOBTime, testTime, treeStrength, treeCorr, numNodes, bestIdx, params, file = paste0(rerfPath, "RandomerForest/R/Results/2018.01.02/Image_stripes_2018_01_02.RData"))
    }
  }
}
