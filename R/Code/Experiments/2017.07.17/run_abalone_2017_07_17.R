rm(list = ls())
rerfPath <- "~/work/tyler/"
inPath <- "~/work/tyler/Data/uci/processed/"
outPath <- paste(rerfPath, "RerF/Results/2017.07.23/", sep = "")
dataSet <- "abalone"
require(parallel)
require(mgcv)
require(dummies)
require(AUC)
source(paste(rerfPath, "RerF/Code/Classifiers/rfr_function.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Classifiers/rerf_eval.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Classifiers/randmat.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Utils/rank.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Utils/strcorr.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Utils/num_nodes.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Utils/mode.R", sep = ""))

# load data
Xtrain <- as.matrix(read.table(paste(inPath, dataSet, ".train.csv", sep = ""), header = F, sep = ",", quote = "", row.names = NULL))
ntrain <- nrow(Xtrain)
p <- ncol(Xtrain) - 1L
Ytrain <- Xtrain[, p + 1L] + 1L
Xtrain <- Xtrain[, -(p + 1L)]
Xtest <- as.matrix(read.table(paste(inPath, dataSet, ".test.csv", sep = ""), header = F, sep = ",", quote = "", row.names = NULL))
ntest <- nrow(Xtest)
Ytest <- Xtest[, p + 1L] + 1L
Xtest <- Xtest[, -(p + 1L)]
Yunique <- unique(Ytest)
nClasses <- length(Yunique)

# specify general parameters
# classifiers <- c("rf", "rerf", "rerfr", "rerfp", "rerfpr", "rerfc", "rerfcr", "frc", "frank", "rr-rf", "rr-rfr")
classifiers <- "rerf"
nTrees <- 1000L
NumCores <- 24L
seed <- 123L
store.predictions <- F

# loop over classifiers
for (m in 1:length(classifiers)) {
  cl <- classifiers[m]
  print(cl)
  
  # specify classifier-specific parameters
  if (cl == "rf" || cl == "rr-rf" || cl == "rr-rfr") {
    randomMatrix <- "rf"
    if (p < 5) {
      mtrys <- 1:p
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1))
    }
    sparsity <- 1/p # this parameter doesn't actually matter for RF
  } else if (cl == "rerf" || cl == "rerfr") {
    randomMatrix <- "binary"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else if (p >= 5 && p <= 100) {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
    }
    sparsity <- (1:min(p, 5))/p
  } else if (cl == "rerfp" || cl == "rerfpr") {
    randomMatrix <- "poisson"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else if (p >= 5 && p <= 100) {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
    }
    sparsity <- (1:min(round(p/2 + 0.1), 5))
  } else if (cl == "rerfc" || cl == "rerfcr") {
    randomMatrix <- "continuous"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else if (p >= 5 && p <= 100) {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
    }
    sparsity <- (1:min(p, 5))/p
  } else if (cl == "frc" || cl == "frank") {
    randomMatrix <- "frc"
    if (p < 5) {
      mtrys <- c(1:p, p^2)
    } else if (p >= 5 && p <= 100) {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 2))
    } else {
      mtrys <- ceiling(p^c(1/4, 1/2, 3/4, 1, 1.5))
    }
    sparsity <- (2:min(p, 5))
  }
  
  if (cl == "rr-rf" || cl == "rr-rfr") {
    rotate <- T
  } else {
    rotate <- F
  }
  
  if (cl == "rerfr" || cl == "rerfpr" || cl == "rerfcr" || cl == "frank" || cl == "rr-rfr") {
    rank.transform <- T
  } else {
    rank.transform <- F
  }
  
  params <- list(trees = nTrees, randomMatrix = randomMatrix, d = mtrys, sparsity = sparsity, rotate = rotate, rank.transform = rank.transform, MinParent = 2L, NumCores = NumCores, seed = seed)
  
  # evaluate classifier
  mdl.eval <- rerf_eval(Xtrain, Ytrain, Xtest, Ytest, params, store.predictions)
  
  # save results
  save(mdl.eval, file = paste(outPath, dataSet, ".", cl, ".results.RData", sep = ""))
}
