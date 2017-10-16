# evaluate RerF on Sparse Parity

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/R/Data/Trunk/dat/Raw/"
library(rerf)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/RerFEval.R"))

# initialize arrays
classifiers <- c("rf", "rerf")

nTrees <- 500L
min.parent <- 2L
max.depth <- 0L
num.cores <- 4L
seed <- 10152017L
ntrain <- 1000L
p <- 100L
trial <- 1L
params <- vector("list", length(classifiers))
names(params) <- classifiers
feature.imp <- vector("list", length(classifiers))
names(feature.imp) <- classifiers

for (m in classifiers) {
    
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
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (1:min(p-1, 5))/p
  } else if (m == "rerfc") {
    random.matrix <- "continuous"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (1:min(p-1, 5))/p
  } else if (m == "rerfp" || m == "rerfpr") {
    random.matrix <- "poisson"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (1:min(ceiling(p/2), 4))
  } else if (m == "frc" || m == "frank") {
    random.matrix <- "frc"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- (2:min(p, 5))
  } else if (m == "frcn") {
    random.matrix <- "frcn"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
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
  
  params[[m]] <- list(trees = nTrees, random.matrix = random.matrix, d = mtrys, sparsity = sparsity, rotate = rotate,
                      rank.transform = rank.transform, min.parent = min.parent, max.depth = max.depth, num.cores = num.cores, seed = seed)
  
  # read in the nth trial of training data for p dimensions 
  D <- read.table(paste0(dataPath, "Train/Trunk_train_set_n", ntrain, "_p", p, "_trial", trial, ".dat"))
  Xtrain <- as.matrix(D[, 1:ncol(D) - 1L])
  Ytrain <- as.integer(D[, ncol(D)]) + 1L
  
  D <- read.table(paste0(dataPath, "Test/Trunk_test_set_p", p, ".dat"))
  Xtest <- as.matrix(D[, 1:ncol(D)-1L])
  Ytest <- as.integer(D[, ncol(D)]) + 1L

  # evaluate models
  res <- RerFEval(Xtrain, Ytrain, Xtest, Ytest, params[[m]])
  best.d <- params[[m]]$d[((res$best.idx - 1L) %% length(params[[m]]$d)) + 1L]
  best.sparsity <- params[[m]]$sparsity[floor((res$best.idx - 1L)/length(params[[m]]$d)) + 1L]
  
  forest <- RerF(X = Xtrain, Y = Ytrain, trees = params[[m]]$trees, min.parent = params[[m]]$min.parent, max.depth = params[[m]]$max.depth,
                 replacement = T, stratify = T, mat.options = list(p, best.d, params[[m]]$random.matrix, best.sparsity), store.impurity = T,
                 num.cores = num.cores, seed = params[[m]]$seed)
  
  feature.imp[[m]] <- FeatureImportance(forest, num.cores = num.cores)
        
  save(feature.imp, params, file = paste(rerfPath, "RandomerForest/R/Results/2017.10.15/Trunk_feature_importance_2017_10_15.RData", sep = ""))
}
