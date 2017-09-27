# evaluate R-RerF on Sparse Parity

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/R/Data/Sparse_parity/dat/Raw/"
source(paste(rerfPath, "R-RerF/Code/Classifiers/rfr_function.R", sep = ""))
source(paste(rerfPath, "R-RerF/Code/Classifiers/randmat.R", sep = ""))
source(paste(rerfPath, "R-RerF/Code/Classifiers/rerf_train.R", sep = ""))
source(paste(rerfPath, "R-RerF/Code/Classifiers/select_model.R", sep = ""))
source(paste(rerfPath, "R-RerF/Code/Classifiers/classprob.R", sep = ""))

# initialize arrays
ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))
nTrials <- 10
testError <- array(vector(mode="double", length=length(ps)*length(ns[[1]])*nTrials), c(length(ns[[1]]),length(ps),nTrials))
trainTime <- array(vector(mode="double", length=length(ps)*length(ns[[1]])*nTrials), c(length(ns[[1]]),length(ps),nTrials))
best.idx <- array(vector(mode="double", length=length(ps)*length(ns[[1]])*nTrials), c(length(ns[[1]]),length(ps),nTrials))

params <- list()
params$randomMatrix <- "binary"
params$MinParent <- 2
params$MaxDepth <- "inf"
params$bagging <- 0.2
params$COOB <- T
params$rescale <- "none"
params$NumCores <- 12

# loop over number of dimensions
for (j in length(ps)) {
  p <- ps[j]
  print(paste("p = ", p, sep = ""))
  
  # read in test data for p dimensions
  D <- read.table(paste(dataPath, "Test/Sparse_parity_test_set_p", p, ".dat", sep = ""))
  Xtest <- data.matrix(D[,1:ncol(D)-1])
  Ytest <- data.matrix(D[,ncol(D)]) + 1
  labels <- sort(unique(Ytest))
  ntest = nrow(Xtest)
  
  if (p <= 5) {
    params$d <- c(1:p, p^c(2, 3))
  } else {
    params$d <- ceiling(p^c(1/4, 1/2, 3/4, 1:3))
  }
  params$sparsity <- 1/p
  
  # loop over number of train observations
  for (i in 2:length(ns[[j]])) {
    ntrain <- ns[[j]][i]
    print(paste("n = ", ntrain, sep = ""))
    
    if (ntrain <= 1000) {
      params$trees <- 1000
    } else {
      params$trees <- 500
    }
    
    # loop over trials
    for (trial in 1:nTrials) {
      print(paste("trial = ",trial, sep = ""))
      
      # read in the nth trial of training data for p dimensions 
      D <- read.table(paste(dataPath, "Train/Sparse_parity_train_set_n", ntrain, "_p", p, "_trial", trial, ".dat", sep = ""))
      Xtrain <- data.matrix(D[,1:ncol(D)-1])
      Ytrain <- data.matrix(D[,ncol(D)]) + 1
      ntrain <- nrow(Xtrain)
      
      # train forests
      forests <- rerf_train(Xtrain, Ytrain, params)
      
      # select best forest
      forests.summary <- select_model(forests$forests, Xtrain, Ytrain)
      best.idx[i,j,trial] <- forests.summary$best.idx
      trainTime[i,j,trial] <- forests$trainTime[best.idx[i,j,trial]]
      
      # compute test error
      testError[i,j,trial] <- error_rate(Xtest, Ytest, forests$forests[[forests.summary$best.idx]])
      save(testError, trainTime, best.idx, params, file = paste(rerfPath, "R-RerF/Results/2017.06.13/Sparse_parity/Sparse_parity_rerf_2017_06_13.RData", sep = ""))
    }
  }
}
