# evaluate classifiers on mnist

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/Data/big/processed/"
library(rerf)
library(AUC)
library(dummies)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/GetCatMap.R"))
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/GetFolds.R"))

classifiers <- c("rf", "rerf", "strerf")
nCl <- length(classifiers)

nTrees <- 500L
min.parent <- 2L
max.depth <- 0L
num.cores <- 40L
seed <- 12132017L

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

dataSet <- "mnist"
Xtrain <- as.matrix(read.table(paste0(dataPath, dataSet, ".train.csv"), header = F, sep = ",", quote = "", row.names = NULL))
catMap <- NULL
p <- ncol(Xtrain) - 1L
Ytrain <- as.integer(Xtrain[, p + 1L]) + 1L
Xtrain <- Xtrain[, -(p + 1L)]

Xtest <- as.matrix(read.table(paste0(dataPath, dataSet, ".test.csv"), header = F, sep = ",", quote = "", row.names = NULL))
Ytest <- as.integer(Xtest[, p + 1L]) + 1L
Xtest <- Xtest[, -(p + 1L)]

testError[[dataSet]] <- vector(mode = "list", length = nCl)
names(testError[[dataSet]]) <- classifiers
testAUC[[dataSet]] <- vector(mode = "list", length = nCl)
names(testAUC[[dataSet]]) <- classifiers
OOBError[[dataSet]] <- vector(mode = "list", length = nCl)
names(OOBError[[dataSet]]) <- classifiers
OOBAUC[[dataSet]] <- vector(mode = "list", length = nCl)
names(OOBAUC[[dataSet]]) <- classifiers
trainTime[[dataSet]] <- vector(mode = "list", length = nCl)
names(trainTime[[dataSet]]) <- classifiers
OOBTime[[dataSet]] <- vector(mode = "list", length = nCl)
names(OOBTime[[dataSet]]) <- classifiers
testTime[[dataSet]] <- vector(mode = "list", length = nCl)
names(testTime[[dataSet]]) <- classifiers
treeStrength[[dataSet]] <- vector(mode = "list", length = nCl)
names(treeStrength[[dataSet]]) <- classifiers
treeCorr[[dataSet]] <- vector(mode = "list", length = nCl)
names(treeCorr[[dataSet]]) <- classifiers
numNodes[[dataSet]] <- vector(mode = "list", length = nCl)
names(numNodes[[dataSet]]) <- classifiers
bestIdx[[dataSet]] <- vector(mode = "list", length = nCl)
names(bestIdx[[dataSet]]) <- classifiers
params[[dataSet]] <- vector(mode = "list", length = nCl)
names(params[[dataSet]]) <- classifiers

for (m in classifiers) {
  
  if (m == "rf") {
    random.matrix <- "rf"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    sparsity <- 1/p # this parameter doesn't actually matter for RF
    iw <- NULL
    ih <- NULL
    patch.min <- NULL
    patch.max <- NULL
  } else if (m == "rerf") {
    random.matrix <- "binary"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    sparsity <- 25/p
    iw <- NULL
    ih <- NULL
    patch.min <- NULL
    patch.max <- NULL
  } else if (m == "strerf") {
    random.matrix <- "image-patch"
    mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    sparsity <- 0
    iw <- sqrt(p)
    ih <- sqrt(p)
    patch.min <- 3L
    patch.max <- 7L
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
  
  params[[dataSet]][[m]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                                 rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores,
                                 seed = seed, cat.map = catMap, iw, ih, patch.min, patch.max)
  
  # testError[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # testAUC[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # OOBError[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # OOBAUC[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # trainTime[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # OOBTime[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # testTime[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # treeStrength[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # treeCorr[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))
  # numNodes[[dataSet]][[m]] <- as.double(rep(NA, length(sparsity)*length(mtrys)))

  # evaluate models
  res <- RerFEval(Xtrain, Ytrain, Xtest, Ytest, params[[dataSet]][[m]])
  
  testError[[dataSet]][[m]] <- res$testError
  testAUC[[dataSet]][[m]] <- res$testAUC
  OOBError[[dataSet]][[m]] <- res$oobError
  OOBAUC[[dataSet]][[m]] <- res$oobAUC
  trainTime[[dataSet]][[m]] <- res$trainTime
  OOBTime[[dataSet]][[m]] <- res$oobTime
  testTime[[dataSet]][[m]] <- res$testTime
  treeStrength[[dataSet]][[m]] <- res$treeStrength
  treeCorr[[dataSet]][[m]] <- res$treeCorrelation
  numNodes[[dataSet]][[m]] <- res$numNodes
  bestIdx[[dataSet]][[m]] <- res$best.idx
  
  save(testError, testAUC, OOBError, OOBAUC, trainTime, OOBTime, testTime, treeStrength, treeCorr, numNodes, bestIdx, params, file = paste0(rerfPath, "RandomerForest/R/Results/2017.12.13/", dataSet, "_2017_12_13.RData"))
}
